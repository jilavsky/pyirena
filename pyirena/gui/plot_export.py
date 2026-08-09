"""
Shared export behaviour for every pyqtgraph plot in the GUI.

Igor Pro, Origin and Matlab users expect that any graph can be copied to the
clipboard (Edit→Copy), saved as an image, and written out as text.  pyIrena had
eight parallel implementations of "add export to a plot" — differing in menu
wording, in what they captured (the plot item vs the whole window), and in the
default folder — so every new export feature had to be added eight times and
some panels silently missed out.  This module is the single implementation.

Panels call one function::

    from pyirena.gui.plot_export import attach_plot_export

    attach_plot_export(plot, self, default_name='unified_fit_iq', window=self)

which appends this block to the plot's right-click (ViewBox) menu:

===========================  ==================================================
Copy graph to clipboard      The plot as an image — paste into PowerPoint,
                             Word, email, an electronic notebook.
Save graph as image…         PNG (default), JPEG or SVG, chosen by the file
                             filter or the extension typed.  PNG is the default
                             because JPEG artifacts badly on line plots + text.
Save whole window as image…  All stacked panels (I(Q) + residuals + …), only
                             when the caller passes ``window``.
Save curve data as CSV…      The plotted curves as text, for Excel / Origin /
                             Matlab / pandas.
Save as Igor Pro ITX…        Igor Pro Text waves with display commands.
===========================  ==================================================

Curve collection is shared by the CSV and ITX exporters
(:func:`collect_plot_curves`): it returns the **linear (pre-log) values** —
``getData()`` returns log10-transformed values when the item has log mode
active, which would double-log in Igor — skips error-bar segments, and picks up
the ``_itx_export`` dicts that step-mode bar charts stash on themselves.

All save dialogs open in the same remembered folder (see
:func:`export_folder`), so the graph you save lands next to the data you loaded
rather than in your home directory.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, List, NamedTuple, Optional

import numpy as np
import pyqtgraph as pg

from pyirena.gui._qt import QApplication, QFileDialog, QMessageBox
from pyirena.gui.table_utils import rows_to_csv_text

log = logging.getLogger(__name__)

__all__ = [
    "PlotCurve",
    "attach_plot_export",
    "collect_plot_curves",
    "copy_plot_to_clipboard",
    "curves_to_csv_text",
    "export_folder",
    "remember_export_folder",
    "save_curves_as_csv",
    "save_itx_from_plot",
    "save_plot_image",
    "save_widget_image",
    "tag_curve_uncertainty",
]

# Width in pixels used for every raster export (image file and clipboard).
EXPORT_WIDTH = 1600

IMAGE_FILTER = (
    "PNG image (*.png);;JPEG image (*.jpg *.jpeg);;"
    "SVG vector image (*.svg);;All files (*)"
)
# Whole-window capture is a widget grab, so no vector option.
WINDOW_IMAGE_FILTER = "PNG image (*.png);;JPEG image (*.jpg *.jpeg);;All files (*)"


# ──────────────────────────────────────────────────────────────────────────
#  Remembered folder (A8)
# ──────────────────────────────────────────────────────────────────────────

# Session-level memory.  Persisted to the state file on a best-effort basis so
# it also survives a restart; an unwritable state file must never break an
# export, hence the broad guard.
_LAST_EXPORT_FOLDER: Optional[str] = None


def _saved_folders() -> tuple[Optional[str], Optional[str]]:
    """Read the remembered export folder and data folder from the state file.

    Isolated in its own function so tests can stub it out — a test that reads
    the real state file passes or fails depending on where the developer last
    saved a graph.
    """
    try:
        from pyirena.state.state_manager import StateManager

        state = StateManager()
        return (
            (state.get("exports") or {}).get("last_folder"),
            (state.get("data_selector") or {}).get("last_folder"),
        )
    except Exception:
        log.debug("suppressed exception reading state for export folder", exc_info=True)
        return None, None


def export_folder(hint: Any = None) -> str:
    """Return the folder an export dialog should open in.

    Resolution order, first hit wins:

    1. the folder of the last successful export in this session — where the
       user *chose* to put the previous graph.  This beats the panel's data
       folder: someone who saved the last figure to the Desktop expects the
       next one to start there, and a panel that always reset to its data
       folder would look like the memory was broken;
    2. ``exports/last_folder`` from the saved state (the same thing, one
       session earlier);
    3. *hint* — the panel's own data folder, used before anything has been
       exported, so the very first save lands next to the data;
    4. ``data_selector/last_folder`` from the saved state;
    5. the user's home directory.

    Args:
        hint: Optional folder (str/Path) the caller considers most relevant.

    Returns:
        An existing directory path as a string.
    """
    saved_export, saved_data = _saved_folders()
    candidates = [_LAST_EXPORT_FOLDER, saved_export, hint, saved_data]

    for candidate in candidates:
        if not candidate:
            continue
        path = Path(str(candidate))
        if path.is_file():
            path = path.parent
        if path.is_dir():
            return str(path)
    return str(Path.home())


def remember_export_folder(file_path: str) -> None:
    """Remember where the user last saved, for the next export dialog."""
    global _LAST_EXPORT_FOLDER
    if not file_path:
        return
    folder = str(Path(file_path).expanduser().parent)
    _LAST_EXPORT_FOLDER = folder
    try:
        from pyirena.state.state_manager import StateManager

        state = StateManager()
        state.set("exports", "last_folder", folder)
        state.save()
    except Exception:
        log.debug("suppressed exception saving export folder", exc_info=True)


def _default_path(default_name: str, ext: str, folder: Any = None) -> str:
    stem = re.sub(r"[^\w\-.]", "_", str(default_name)).strip("_") or "pyirena_graph"
    if stem.lower().endswith(ext):
        stem = stem[: -len(ext)]
    return str(Path(export_folder(folder)) / f"{stem}{ext}")


# ──────────────────────────────────────────────────────────────────────────
#  Curve collection (shared by CSV and ITX)
# ──────────────────────────────────────────────────────────────────────────

class PlotCurve(NamedTuple):
    """One named curve taken off a plot, in linear (physical) units."""

    name: str
    x: np.ndarray
    y: np.ndarray
    color: str
    dy: Optional[np.ndarray]
    #: How the curve is drawn on screen: "markers", "line" or "both".  Igor
    #: draws every imported wave as a line by default, so a scatter of data
    #: points would arrive as a zig-zag line unless the ITX says otherwise.
    style: str = "line"


def _item_style(item) -> str:
    """Report whether an item is drawn as markers, a line, or both."""
    opts = getattr(item, "opts", None) or {}
    if isinstance(item, pg.ScatterPlotItem):
        return "markers"
    has_symbol = opts.get("symbol") is not None
    pen = opts.get("pen")
    has_line = pen is not None
    if has_symbol and has_line:
        return "both"
    if has_symbol:
        return "markers"
    return "line"


def _get_item_color_hex(item) -> str:
    """Extract the primary color of a PlotDataItem as an #rrggbb hex string."""
    try:
        pen = item.opts.get("pen")
        if pen is not None:
            return pg.mkPen(pen).color().name()
    except Exception:
        log.debug("suppressed exception", exc_info=True)
    try:
        brush = item.opts.get("symbolBrush")
        if brush is not None:
            return pg.mkBrush(brush).color().name()
    except Exception:
        log.debug("suppressed exception", exc_info=True)
    return "#000000"


def tag_curve_uncertainty(item, dy) -> None:
    """Record the uncertainties belonging to a plotted curve.

    Error bars are drawn as NaN-separated line segments, which no exporter can
    turn back into a ``dY`` column, so the panel has to say which array goes
    with which curve.  ``plot_iq_data`` does this for the plots it builds; call
    this from any panel that draws its own data scatter, right after plotting::

        item = self.main_plot.plot(q, I, symbol='o', name=label)
        tag_curve_uncertainty(item, dI)

    Without it the CSV and ITX exports quietly lose the uncertainties even
    though the error bars are visible on screen.
    """
    if item is not None and dy is not None:
        item._itx_dI = np.asarray(dy, dtype=float)


def _item_xy(item):
    """Return (x, y) in linear units for any pyqtgraph data item, or None.

    ``PlotDataItem.getOriginalDataset()`` gives pre-log-transform values;
    ``getData()`` would return log10 values on a log plot and double-log in
    Igor.  ``ScatterPlotItem`` and bare ``PlotCurveItem`` have no such
    distinction (pyIrena only uses them on linear plots) — but they are also
    not ``PlotDataItem``, so a naive scan misses them entirely, which is why
    the Simple Fits linearization plot and the collected-values scatter used to
    report "no named data curves".
    """
    if isinstance(item, pg.PlotDataItem):
        return item.getOriginalDataset()
    if isinstance(item, (pg.ScatterPlotItem, pg.PlotCurveItem)):
        try:
            return item.getData()
        except Exception:
            log.debug("suppressed exception reading item data", exc_info=True)
    return None, None


def _item_name(item) -> str:
    """Best available name for a data item ('' when it has none)."""
    try:
        if isinstance(item, pg.PlotDataItem):
            return item.name() or ""
    except Exception:
        log.debug("suppressed exception reading item name", exc_info=True)
    opts = getattr(item, "opts", None)
    if isinstance(opts, dict):
        return opts.get("name") or ""
    return ""


def _fallback_name(plot: pg.PlotItem) -> str:
    """Name to use for curves the panel never labelled — the Y axis label."""
    label = getattr(plot.getAxis("left"), "labelText", "") or ""
    label = re.sub(r"\s+", " ", label).strip()
    return label or "Curve"


def collect_plot_curves(plot: pg.PlotItem) -> List[PlotCurve]:
    """Return every data curve on *plot*, in linear (physical) units.

    Named curves are preferred.  If the plot has none — residual panels, the
    Simple Fits linearization, the Contrast (Δρ)² plot and the collected-values
    scatter all draw curves without a legend name — the unnamed ones are
    exported instead, labelled from the Y-axis. Refusing to export a plot the
    user can plainly see was the old behaviour and simply looked broken.

    Always skipped: curves that are more than 30 % NaN, which is how error-bar
    segments are drawn.  Items carrying an explicit ``_itx_export`` dict (step
    mode bar charts, which are not tracked in ``dataItems``) are picked up too,
    and their stored arrays take precedence because they are already linear and
    correctly shaped.

    Uncertainties recorded with :func:`tag_curve_uncertainty` (or by
    ``plot_iq_data``) come back in :attr:`PlotCurve.dy`.
    """
    scan_items = list(plot.listDataItems())
    for it in getattr(plot, "items", []):
        if it not in scan_items and (
            getattr(it, "_itx_export", None) is not None
            or isinstance(it, (pg.ScatterPlotItem, pg.PlotCurveItem))
        ):
            scan_items.append(it)

    named: List[PlotCurve] = []
    unnamed: List[PlotCurve] = []

    for item in scan_items:
        # Explicit export data takes precedence (linear units, correct shape).
        exp = getattr(item, "_itx_export", None)
        if exp is not None:
            name = exp.get("name") or ""
            x_data = np.asarray(exp.get("x"), dtype=float)
            y_data = np.asarray(exp.get("y"), dtype=float)
            if x_data.size >= 2 and y_data.size >= 2:
                mask = np.isfinite(x_data) & np.isfinite(y_data)
                if mask.sum() >= 2:
                    curve = PlotCurve(name, x_data[mask], y_data[mask],
                                      _get_item_color_hex(item), None,
                                      _item_style(item))
                    (named if name else unnamed).append(curve)
            continue

        x_data, y_data = _item_xy(item)
        if x_data is None or y_data is None or len(x_data) < 2:
            continue
        x_data = np.asarray(x_data, dtype=float)
        y_data = np.asarray(y_data, dtype=float)
        # Skip error-bar segments (NaN-separated lines: >30 % NaN values).
        if np.sum(np.isnan(y_data)) > len(y_data) * 0.3:
            continue
        mask = np.isfinite(x_data) & np.isfinite(y_data)
        if mask.sum() < 2:
            continue

        raw_dy = getattr(item, "_itx_dI", None)
        dy_arr: Optional[np.ndarray] = None
        if raw_dy is not None:
            raw_dy = np.asarray(raw_dy, dtype=float)
            if raw_dy.shape == x_data.shape:
                dy_arr = raw_dy[mask]
            elif raw_dy.shape == x_data[mask].shape:
                dy_arr = raw_dy
            if dy_arr is not None and not np.any(np.isfinite(dy_arr) & (dy_arr > 0)):
                dy_arr = None      # all zeros/NaN — no uncertainty to export

        name = _item_name(item)
        curve = PlotCurve(name, x_data[mask], y_data[mask],
                          _get_item_color_hex(item), dy_arr, _item_style(item))
        (named if name else unnamed).append(curve)

    if named:
        return named

    # Nothing was labelled — export what is on screen anyway.
    base = _fallback_name(plot)
    if len(unnamed) == 1:
        return [unnamed[0]._replace(name=base)]
    return [c._replace(name=f"{base} {i + 1}") for i, c in enumerate(unnamed)]


# ──────────────────────────────────────────────────────────────────────────
#  Clipboard and image files (A3)
# ──────────────────────────────────────────────────────────────────────────

def copy_plot_to_clipboard(plot: pg.PlotItem, parent=None,
                           width: int = EXPORT_WIDTH) -> bool:
    """Put *plot* on the clipboard as an image.

    Pastes straight into PowerPoint, Word, email or an electronic notebook —
    the single most-requested graph action for Igor users, who reach for
    Edit→Copy without thinking about it.

    Returns:
        True on success; on failure a warning box is shown and False returned.
    """
    from pyqtgraph.exporters import ImageExporter

    try:
        exporter = ImageExporter(plot)
        exporter.parameters()["width"] = width
        image = exporter.export(toBytes=True)
        clipboard = QApplication.clipboard()
        if clipboard is None:              # headless / no window system
            return False
        clipboard.setImage(image)
        return True
    except Exception as exc:
        log.debug("clipboard export failed", exc_info=True)
        QMessageBox.warning(parent, "Copy failed",
                            f"Could not copy the graph:\n{exc}")
        return False


def save_plot_image(
    plot: pg.PlotItem,
    parent=None,
    default_name: str = "pyirena_graph",
    folder: Any = None,
    width: int = EXPORT_WIDTH,
) -> Optional[str]:
    """Save *plot* as PNG, JPEG or SVG, chosen in the file dialog.

    PNG is the default: JPEG compression artifacts are clearly visible on the
    thin lines and small text of a log-log SAS plot.  SVG uses pyqtgraph's
    vector exporter (best-effort — complex scenes can differ slightly).

    Returns:
        The path written, or None if cancelled or failed.
    """
    path, selected = QFileDialog.getSaveFileName(
        parent, "Save graph as image",
        _default_path(default_name, ".png", folder),
        IMAGE_FILTER,
    )
    if not path:
        return None

    suffix = Path(path).suffix.lower()
    if suffix not in (".png", ".jpg", ".jpeg", ".svg"):
        # Honour the filter the user picked when they typed no extension.
        suffix = ".svg" if "svg" in (selected or "").lower() else (
            ".jpg" if "jpeg" in (selected or "").lower() else ".png"
        )
        path += suffix

    try:
        if suffix == ".svg":
            from pyqtgraph.exporters import SVGExporter

            SVGExporter(plot).export(path)
        else:
            from pyqtgraph.exporters import ImageExporter

            exporter = ImageExporter(plot)
            exporter.parameters()["width"] = width
            exporter.export(path)
    except Exception as exc:
        QMessageBox.warning(parent, "Export failed",
                            f"Could not save image:\n{exc}")
        return None

    remember_export_folder(path)
    return path


def save_widget_image(
    widget,
    parent=None,
    default_name: str = "pyirena_window",
    folder: Any = None,
) -> Optional[str]:
    """Save a whole window (all stacked plots) as a PNG or JPEG screenshot."""
    path, selected = QFileDialog.getSaveFileName(
        parent or widget, "Save whole window as image",
        _default_path(default_name, ".png", folder),
        WINDOW_IMAGE_FILTER,
    )
    if not path:
        return None

    suffix = Path(path).suffix.lower()
    if suffix not in (".png", ".jpg", ".jpeg"):
        suffix = ".jpg" if "jpeg" in (selected or "").lower() else ".png"
        path += suffix

    fmt = "JPEG" if suffix in (".jpg", ".jpeg") else "PNG"
    try:
        ok = widget.grab().save(path, fmt, 95)
    except Exception as exc:
        QMessageBox.warning(parent or widget, "Export failed",
                            f"Could not save image:\n{exc}")
        return None
    if not ok:
        QMessageBox.warning(parent or widget, "Export failed",
                            f"Could not write {path}")
        return None

    remember_export_folder(path)
    return path


# ──────────────────────────────────────────────────────────────────────────
#  Curve data as text (A4)
# ──────────────────────────────────────────────────────────────────────────

def curves_to_csv_text(curves: List[PlotCurve]) -> str:
    """Render plotted curves as CSV.

    Each curve contributes its own ``X (label)`` / ``Y (label)`` pair — and
    ``dY (label)`` where uncertainties exist — because curves on one plot
    rarely share a Q grid (data and model usually do; a fitted background or a
    distribution does not).  Shorter columns are padded with blanks, which is
    what Excel, Origin and pandas all expect.
    """
    headers: List[str] = []
    columns: List[np.ndarray] = []
    for curve in curves:
        headers.append(f"X ({curve.name})")
        columns.append(curve.x)
        headers.append(f"Y ({curve.name})")
        columns.append(curve.y)
        if curve.dy is not None:
            headers.append(f"dY ({curve.name})")
            columns.append(curve.dy)

    n_rows = max((len(c) for c in columns), default=0)
    rows = [
        [(col[r] if r < len(col) else None) for col in columns]
        for r in range(n_rows)
    ]
    return rows_to_csv_text(headers, rows)


def save_curves_as_csv(
    plot: pg.PlotItem,
    parent=None,
    default_name: str = "pyirena_curves",
    folder: Any = None,
) -> Optional[str]:
    """Write the plotted curves to a CSV file chosen by the user."""
    curves = collect_plot_curves(plot)
    if not curves:
        QMessageBox.warning(parent, "No data",
                            "No named data curves found to export.")
        return None

    path, _ = QFileDialog.getSaveFileName(
        parent, "Save curve data as CSV",
        _default_path(default_name, ".csv", folder),
        "CSV files (*.csv);;Text files (*.txt);;All files (*)",
    )
    if not path:
        return None
    if Path(path).suffix.lower() not in (".csv", ".txt"):
        path += ".csv"

    try:
        with open(path, "w", newline="", encoding="utf-8") as fh:
            fh.write(curves_to_csv_text(curves))
    except OSError as exc:
        QMessageBox.warning(parent, "Export failed",
                            f"Could not save file:\n{exc}")
        return None

    remember_export_folder(path)
    return path


# ──────────────────────────────────────────────────────────────────────────
#  Igor Pro ITX  (moved verbatim from sas_plot.py)
# ──────────────────────────────────────────────────────────────────────────

def _itx_folder_cmds(technique: Optional[str],
                     sample: Optional[str]) -> tuple[list[str], list[str]]:
    """Build the Igor ``X`` commands that route loaded waves into a folder.

    Returns ``(open_lines, close_lines)`` such that WAVES blocks emitted
    between them land in ``root:<technique>:<sample>`` instead of ``root:``.
    This lets users import many ITX files into one Igor experiment without
    wave-name collisions — each sample gets its own data folder.

    Returns ``([], [])`` when *technique* is falsy (top-level export, the
    historical behaviour).  The *technique* segment is sanitised to a strict
    Igor object name; the *sample* segment uses a quoted liberal name so it
    stays human-readable (Igor 8+ / Igor 10 long names).
    """
    if not technique:
        return [], []
    tech = re.sub(r'[^A-Za-z0-9_]', '_', str(technique))
    if tech and tech[0].isdigit():
        tech = 'f_' + tech
    tech = tech[:31]
    if not tech:
        return [], []
    open_lines = [f'X NewDataFolder/O/S root:{tech}']
    if sample:
        s = str(sample).strip()
        # Drop a single trailing file extension for a cleaner folder name.
        s = re.sub(r'\.(h5|hdf5|hdf|dat|txt|csv|itx|pxp|xml|nxs)$', '',
                   s, flags=re.IGNORECASE)
        # Single quote is the liberal-name delimiter — replace any in the name.
        s = s.replace("'", '_').strip()[:200]
        if s:
            open_lines.append(f"X NewDataFolder/O/S root:{tech}:'{s}'")
    return open_lines, ['X SetDataFolder root:']


def save_itx_from_plot(
    plot: pg.PlotItem,
    parent,
    default_name: str = 'pyirena_graph',
    x_label: str | None = None,
    y_label: str | None = None,
    title: str | None = None,
    technique: str | None = None,
    sample: str | None = None,
    folder: Any = None,
) -> None:
    """Export named data curves from *plot* as an Igor Pro Text (.itx) file.

    Iterates all named curves on *plot* (scatter data and model curves, via
    :func:`collect_plot_curves`) and writes them as Igor Pro waves with
    display, log-axis, color, label, and legend commands.  Error-bar segments
    (NaN-separated lines without a name) are automatically skipped.

    Auto-extracts axis labels and title from *plot* when the corresponding
    parameters are ``None``.

    When *technique* (and optionally *sample*) is supplied — or discoverable
    as ``parent._itx_technique`` / ``parent._itx_sample_label`` — the waves are
    routed into ``root:<technique>:<sample>`` so multiple ITX imports into one
    Igor experiment do not collide on wave names.
    """
    # Fall back to attributes stashed on the parent window by the tool panels.
    if technique is None:
        technique = getattr(parent, '_itx_technique', None)
    if sample is None:
        sample = getattr(parent, '_itx_sample_label', None)

    # Auto-extract labels / title from the plot_item if not provided.
    if x_label is None:
        x_label = getattr(plot.getAxis('bottom'), 'labelText', '') or ''
    if y_label is None:
        y_label = getattr(plot.getAxis('left'), 'labelText', '') or ''
    if title is None:
        try:
            title = plot.titleLabel.text or ''
        except Exception:
            title = ''

    named_items = collect_plot_curves(plot)
    if not named_items:
        QMessageBox.warning(parent, 'No data',
                            'No named data curves found to export.')
        return

    # Build a default filename from technique + sample so multiple exports
    # don't collide in the same directory (e.g. "Sizes_AeroGel_500C.itx").
    if technique or sample:
        parts = [p for p in (technique, sample) if p]
        stem = re.sub(r'[^\w\-.]', '_', '_'.join(parts)).strip('_') or default_name
        # Drop duplicate/trailing extension that may already be in sample name.
        stem = re.sub(r'\.(h5|hdf5|hdf|dat|txt|csv|itx|pxp|xml|nxs)$', '',
                      stem, flags=re.IGNORECASE)
    else:
        stem = default_name
    filepath, _ = QFileDialog.getSaveFileName(
        parent, 'Save as Igor Pro ITX',
        _default_path(stem, '.itx', folder),
        'Igor Pro Text (*.itx);;All files (*)',
    )
    if not filepath:
        return
    if not filepath.lower().endswith('.itx'):
        filepath += '.itx'

    log_x = bool(getattr(plot.getAxis('bottom'), 'logMode', False))
    log_y = bool(getattr(plot.getAxis('left'),   'logMode', False))

    def _safe_name(prefix: str, label: str, suffix: str = '') -> str:
        """Build an Igor wave name of at most 31 characters.

        The numeric suffix is what keeps waves from different curves apart, so
        the *label* is truncated to make room for it — truncating the whole
        name would give two long-labelled curves the same wave name and Igor
        would silently overwrite the first.
        """
        clean = re.sub(r'[^A-Za-z0-9_]', '_', label)
        budget = 31 - len(prefix) - len(suffix)
        name = f'{prefix}{clean[:max(budget, 1)]}{suffix}'
        if name and name[0].isdigit():
            name = 'w_' + name[:29]
        return name or 'wave'

    def _hex_to_igor(h: str) -> tuple[int, int, int]:
        h = h.lstrip('#')
        if len(h) == 6:
            r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
        else:
            r = g = b = 0
        return r * 257, g * 257, b * 257

    folder_open, folder_close = _itx_folder_cmds(technique, sample)

    lines = ['IGOR'] + folder_open
    # (xn, yn, en_or_None, label, color, style)
    wave_info: list[tuple[str, str, str | None, str, str, str]] = []
    n = len(named_items)

    for i, curve in enumerate(named_items):
        lbl, x_arr, y_arr, color, dI_arr, style = curve
        suffix = f'_{i + 1:02d}' if n > 1 else ''
        xn = _safe_name('X_', lbl, suffix)
        yn = _safe_name('Y_', lbl, suffix)
        en = _safe_name('Yerr_', lbl, suffix) if dI_arr is not None else None
        wave_info.append((xn, yn, en, lbl, color, style))

        lines += [f'WAVES/D  {xn}', 'BEGIN']
        lines += [f'  {v:.10g}' for v in x_arr]
        lines += ['END', f'WAVES/D  {yn}', 'BEGIN']
        lines += [f'  {v:.10g}' for v in y_arr]
        lines.append('END')

        if en is not None and dI_arr is not None:
            lines += [f'WAVES/D  {en}', 'BEGIN']
            lines += [f'  {v:.10g}' for v in dI_arr]
            lines.append('END')

    lines.append('')
    for j, (xn, yn, en, lbl, _, _) in enumerate(wave_info):
        if j == 0:
            lines.append(f'X Display {yn} vs {xn} as "{lbl}"')
        else:
            lines.append(f'X AppendToGraph {yn} vs {xn}')

    if log_x:
        lines.append('X ModifyGraph log(bottom)=1')
    if log_y:
        lines.append('X ModifyGraph log(left)=1')

    for _, yn, en, _, color, style in wave_info:
        r, g, b = _hex_to_igor(color)
        lines.append(f'X ModifyGraph rgb({yn})=({r},{g},{b})')
        # Igor draws every imported wave as a line.  Reproduce how the curve
        # looks in pyIrena, otherwise a scatter of data points arrives as a
        # zig-zag line and is indistinguishable from a model curve.
        #   mode 3 = markers only, 4 = lines + markers, 0 = lines
        #   marker 19 = filled circle (Igor's default data marker)
        if style == 'markers':
            lines.append(f'X ModifyGraph mode({yn})=3,marker({yn})=19,'
                         f'msize({yn})=2')
        elif style == 'both':
            lines.append(f'X ModifyGraph mode({yn})=4,marker({yn})=19,'
                         f'msize({yn})=2')
        else:
            lines.append(f'X ModifyGraph mode({yn})=0,lsize({yn})=1.5')
        if en is not None:
            lines.append(f'X ErrorBars {yn} Y,wave=({en},{en})')

    if x_label:
        lines.append(f'X Label bottom "{x_label}"')
    if y_label:
        lines.append(f'X Label left "{y_label}"')
    if title:
        lines.append(f'X TextBox/C/N=title0/A=MC/X=0/Y=5 "{title}"')

    legend_parts = [f'\\\\s({yn}) {lbl}' for _, yn, en, lbl, _, _ in wave_info]
    if legend_parts:
        legend_text = '\\r'.join(legend_parts)
        lines.append(f'X Legend/C/N=text0 "{legend_text}"')

    # Restore the current data folder to root: after routing waves into a
    # sample subfolder (no-op when no folder routing was requested).
    lines += folder_close

    try:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines) + '\n')
    except Exception as exc:
        QMessageBox.warning(parent, 'Export Failed',
                            f'Could not save file:\n{exc}')
        return

    remember_export_folder(filepath)


# ──────────────────────────────────────────────────────────────────────────
#  The one entry point panels call
# ──────────────────────────────────────────────────────────────────────────

def attach_plot_export(
    plot: pg.PlotItem,
    parent=None,
    default_name: str = "pyirena_graph",
    window=None,
    folder: Any = None,
    x_label: str = "",
    y_label: str = "",
    title: str = "",
    technique: Optional[str] = None,
    sample: Optional[str] = None,
    itx: bool = True,
) -> pg.PlotItem:
    """Append the standard export block to a plot's right-click menu.

    Args:
        plot: The pyqtgraph ``PlotItem`` to equip.
        parent: Widget used as the parent of file dialogs and message boxes.
        default_name: File stem offered in the save dialogs.
        window: Widget to capture for "Save whole window as image…".  Pass the
            window holding several stacked plots; omit for a single plot.
        folder: Folder the dialogs should open in — pass the panel's data
            folder when it has one.  May be a str, Path, or a zero-argument
            callable evaluated at click time (data folders change as the user
            loads files).
        x_label, y_label, title: Overrides for the ITX export; by default they
            are read off the plot.
        technique, sample: Igor data-folder routing for the ITX export.
        itx: Set False where the window already has a richer, purpose-built ITX
            writer (the collected-values window exports labelled X/Y/error
            waves), so the menu does not offer two entries doing the same job.

    Returns:
        The same plot, so the call can be chained.
    """
    vb = plot.getViewBox()
    if vb is None or getattr(vb, "menu", None) is None:
        return plot                        # e.g. a plot with menus disabled

    def _folder():
        return folder() if callable(folder) else folder

    menu = vb.menu
    menu.addSeparator()

    act_copy = menu.addAction("Copy graph to clipboard")
    act_copy.setToolTip("Copy the graph as an image, ready to paste into a "
                        "document, presentation or email.")
    act_copy.triggered.connect(
        lambda checked=False: copy_plot_to_clipboard(plot, parent)
    )

    act_img = menu.addAction("Save graph as image…")
    act_img.setToolTip("Save this graph as PNG (recommended), JPEG or SVG.")
    act_img.triggered.connect(
        lambda checked=False: save_plot_image(plot, parent, default_name, _folder())
    )

    if window is not None:
        act_win = menu.addAction("Save whole window as image…")
        act_win.setToolTip("Save every panel of this window (data, residuals, "
                           "distribution) as one image.")
        act_win.triggered.connect(
            lambda checked=False: save_widget_image(
                window, parent, f"{default_name}_window", _folder()
            )
        )

    act_csv = menu.addAction("Save curve data as CSV…")
    act_csv.setToolTip("Save the plotted curves as text for Excel, Origin, "
                       "Matlab or Python.")
    act_csv.triggered.connect(
        lambda checked=False: save_curves_as_csv(plot, parent, default_name, _folder())
    )

    if itx:
        act_itx = menu.addAction("Save as Igor Pro ITX…")
        act_itx.setToolTip("Save the plotted curves as Igor Pro Text waves.")
        act_itx.triggered.connect(
            lambda checked=False: save_itx_from_plot(
                plot, parent, default_name,
                x_label or None, y_label or None, title or None,
                technique, sample, _folder(),
            )
        )

    return plot
