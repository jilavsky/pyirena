"""
Reusable SAS I(Q) log-log plot utilities for pyIrena.

This module provides a standard set of helpers for creating and populating
pyqtgraph plot items that display I(Q) vs Q data.  Using these helpers
ensures a consistent look-and-feel (colors, error bars, cursors, grid,
labels) and shared features (right-click "Save as JPEG") across ALL
pyIrena panels.

Design notes
------------
* All public helpers work with pyqtgraph ``PlotItem`` objects.
* Data values passed to ``plot_iq_data`` / ``plot_iq_model`` must be in
  **physical (linear) units**.  ``make_sas_plot`` calls
  ``setLogMode(x=True, y=True)`` so pyqtgraph applies the log10 transform
  internally — you never pre-transform the data.
* Cursor positions and ``TextItem`` positions go into the **ViewBox
  coordinate system**, which equals log10(data) when setLogMode is active.
  Use ``get_cursor_q_range()`` to read back cursor positions in linear units.

Quick-start example
-------------------
::

    import pyqtgraph as pg
    from pyirena.gui.sas_plot import (
        make_sas_plot, plot_iq_data, plot_iq_model,
        make_cursors, get_cursor_q_range,
    )

    gl = pg.GraphicsLayoutWidget()
    plot = make_sas_plot(
        gl, row=0, col=0,
        x_label='Q (Å⁻¹)', y_label='I (cm⁻¹)',
        parent_widget=my_widget,
        jpeg_default_name='my_graph',
    )
    plot_iq_data(plot, q, I, dI)
    cursor_a, cursor_b = make_cursors(plot, q.min(), q.max())
    q_min, q_max = get_cursor_q_range(cursor_a, cursor_b)
"""

from __future__ import annotations

import numpy as np
from pathlib import Path
import pyqtgraph as pg

try:
    from PySide6.QtWidgets import QFileDialog, QMessageBox
    from PySide6.QtCore import Qt
except ImportError:
    from PyQt6.QtWidgets import QFileDialog, QMessageBox
    from PyQt6.QtCore import Qt


# ===========================================================================
# Style constants — change here to update all tools simultaneously
# ===========================================================================

class SASPlotStyle:
    """Central definition of visual style for all pyIrena SAS plots.

    Edit these class-level constants to change the appearance of every plot
    that uses the ``sas_plot`` helpers.
    """

    # ── Data scatter ──────────────────────────────────────────────────────
    DATA_BRUSH       = pg.mkBrush(44, 62, 80, 200)    # #2c3e50 dark slate
    DATA_SIZE        = 5                               # symbol diameter (px)

    # ── Error bars ────────────────────────────────────────────────────────
    ERROR_PEN        = pg.mkPen((127, 140, 141, 180), width=1)  # #7f8c8d grey
    ERROR_CAP_FRAC   = 0.05    # half-cap width as fraction of Q in log space

    # ── Model / fit line ─────────────────────────────────────────────────
    FIT_PEN          = pg.mkPen(color=(200, 30, 30), width=2)   # red

    # ── Residuals scatter ─────────────────────────────────────────────────
    RESID_BRUSH      = pg.mkBrush(127, 140, 141, 200)           # grey
    RESID_SIZE       = 4

    # ── Cursors ───────────────────────────────────────────────────────────
    CURSOR_A_PEN     = pg.mkPen('#e74c3c', width=2)             # red
    CURSOR_B_PEN     = pg.mkPen('#3498db', width=2)             # blue
    CURSOR_A_COLOR   = '#e74c3c'
    CURSOR_B_COLOR   = '#3498db'

    # ── Annotation text ───────────────────────────────────────────────────
    ANNOT_COLOR      = (40, 40, 40)    # near-black

    # ── Legend ────────────────────────────────────────────────────────────
    LEGEND_TEXT_COLOR = 'k'     # black — visible on the white background used by all tools

    # ── Grid ──────────────────────────────────────────────────────────────
    GRID_ALPHA       = 0.3


# ===========================================================================
# RadiusAxisItem — top axis showing R = π/Q
# ===========================================================================

class RadiusAxisItem(pg.AxisItem):
    """Top axis for I(Q) log-log plots showing feature radius R = π/Q [Å].

    The ViewBox x-coordinates are log10(Q).  This axis selects nice round
    radius values (1, 2, 5, 10, 20, 50, 100, … Å), converts them to
    log10(Q) positions, and labels them.  Labels decrease from left to
    right because R ∝ 1/Q.
    """

    _NICE_MULTIPLIERS = (1, 2, 5)
    _MAX_TICKS = 8  # keep the top axis readable

    def tickValues(self, minVal, maxVal, size):
        import math
        PI = math.pi
        r_min = PI / (10.0 ** maxVal)
        r_max = PI / (10.0 ** minVal)
        if r_min <= 0 or r_max <= 0 or r_min >= r_max:
            return []
        decade_lo = math.floor(math.log10(r_min)) - 1
        decade_hi = math.ceil(math.log10(r_max)) + 1
        positions = []
        for exp in range(int(decade_lo), int(decade_hi) + 1):
            for mult in self._NICE_MULTIPLIERS:
                r = mult * (10.0 ** exp)
                if r_min <= r <= r_max:
                    positions.append(math.log10(PI / r))
        if not positions:
            return []
        # Thin to at most _MAX_TICKS evenly-spaced entries
        if len(positions) > self._MAX_TICKS:
            step = max(1, len(positions) // self._MAX_TICKS)
            positions = positions[::step]
        return [(1.0, positions)]

    def tickStrings(self, values, scale, spacing):
        import math
        PI = math.pi
        strings = []
        for v in values:
            R = PI / (10.0 ** v)
            if R >= 10:
                strings.append(f'{R:.0f}')
            elif R >= 1:
                strings.append(f'{R:.1f}')
            else:
                strings.append(f'{R:.2f}')
        return strings


# ===========================================================================
# _SafeInfiniteLine — prevents PySide6 segfaults on macOS ARM
# ===========================================================================

class _SafeInfiniteLine(pg.InfiniteLine):
    """pg.InfiniteLine subclass that prevents PySide6/Shiboken segfaults.

    On PySide6 6.x running on macOS ARM, if a Python exception escapes a
    ``QGraphicsItem`` virtual-method override (``mouseMoveEvent``,
    ``mouseDragEvent``), Shiboken calls ``PepException_GetArgs()`` on an
    already-invalid pointer and segfaults.  Wrapping those overrides in
    ``try/except`` ensures that no exception can reach the C++→Python
    boundary.
    """

    def mouseMoveEvent(self, ev):
        try:
            super().mouseMoveEvent(ev)
        except Exception:
            ev.ignore()

    def mouseDragEvent(self, ev):
        try:
            super().mouseDragEvent(ev)
        except Exception:
            ev.ignore()


# ---------------------------------------------------------------------------
# Defensive monkey-patch for pyqtgraph GraphicsScene.sendDragEvent
# ---------------------------------------------------------------------------
# When a plot is rebuilt (e.g. via GraphicsLayoutWidget.clear()), Qt's
# parent-child cleanup destroys the C++ side of cursor (InfiniteLine) items.
# pyqtgraph's GraphicsScene caches ``lastHoverEvent`` which can still hold
# Python references to those now-deleted items.  When the user later starts
# a drag, ``sendDragEvent`` reads ``acceptedItem.scene()`` on the dead item
# and Qt raises ``RuntimeError: Internal C++ object already deleted``.
#
# The patch wraps the unsafe call in try/except so the dead reference is
# cleared and the drag falls back to the "find nearby item" path.
def _install_sendDragEvent_safeguard() -> None:
    import pyqtgraph.GraphicsScene.GraphicsScene as _GSModule
    Scene = _GSModule.GraphicsScene
    if getattr(Scene, '_pyirena_safe_sendDragEvent', False):
        return
    _orig = Scene.sendDragEvent

    def _safe_sendDragEvent(self, ev, init=False, final=False):
        try:
            if init and self.dragItem is None and self.lastHoverEvent is not None:
                # Probe whether any cached drag-target item has been deleted.
                # Iterating the dict is enough to invalidate stale references
                # before the original method dereferences them.
                for btn, item in list(self.lastHoverEvent.dragItems().items()):
                    try:
                        _ = item.scene()
                    except RuntimeError:
                        # Underlying C++ object is gone — drop the cache so the
                        # original method does not try to use it.
                        self.lastHoverEvent = None
                        break
        except Exception:
            self.lastHoverEvent = None
        try:
            return _orig(self, ev, init=init, final=final)
        except RuntimeError:
            # Last-ditch defense: any lingering deleted-object access aborts
            # the drag instead of crashing the GUI.
            self.dragItem = None
            self.lastHoverEvent = None
            return None

    Scene.sendDragEvent = _safe_sendDragEvent
    Scene._pyirena_safe_sendDragEvent = True


_install_sendDragEvent_safeguard()


# ===========================================================================
# Plot factory
# ===========================================================================

def make_sas_plot(
    graphics_layout: pg.GraphicsLayoutWidget,
    row: int,
    col: int,
    x_label: str = 'Q  (Å⁻¹)',
    y_label: str = 'I',
    title: str | None = None,
    x_link=None,
    log_x: bool = True,
    log_y: bool = True,
    parent_widget=None,
    jpeg_default_name: str = 'pyirena_graph',
) -> pg.PlotItem:
    """Create a log-log (or semi-log) plot and configure it with pyIrena style.

    Parameters
    ----------
    graphics_layout : pg.GraphicsLayoutWidget
        The parent graphics layout.
    row, col : int
        Grid position within *graphics_layout*.
    x_label, y_label : str
        Axis labels (support Unicode, e.g. ``'Q  (Å⁻¹)'``).
    title : str or None
        Optional title shown above the plot.
    x_link : PlotItem or None
        If given, the x-axis is linked to this plot (shared pan/zoom).
    log_x, log_y : bool
        Whether to apply log10 to x / y axes.  Pass ``True`` for both for
        an I(Q) vs Q plot.  Pass ``log_x=True, log_y=False`` for residuals.
    parent_widget : QWidget or None
        Parent for file dialogs triggered by the "Save as JPEG" action.
        Pass ``None`` to suppress the JPEG export menu entry.
    jpeg_default_name : str
        Default file stem for the JPEG export dialog.

    Returns
    -------
    pg.PlotItem
        Configured plot item.  Pass **physical (linear)** values to all
        pyqtgraph items added to this plot — pyqtgraph applies the log10
        transform automatically when setLogMode is active.
    """
    axis_items = {}
    if log_x:
        axis_items['top'] = RadiusAxisItem(orientation='top')
    plot = graphics_layout.addPlot(row=row, col=col, axisItems=axis_items)
    plot.setLogMode(x=log_x, y=log_y)
    plot.setLabel('left',   y_label)
    plot.setLabel('bottom', x_label)
    if log_x and 'top' in axis_items:
        plot.getAxis('top').setLabel('R = π/Q  (Å)', **{'color': '#888', 'font-size': '9pt'})
    plot.showGrid(x=True, y=True, alpha=SASPlotStyle.GRID_ALPHA)
    plot.getAxis('left').enableAutoSIPrefix(False)
    plot.getAxis('bottom').enableAutoSIPrefix(False)
    if title:
        plot.setTitle(title)
    if x_link is not None:
        plot.setXLink(x_link)
    if parent_widget is not None:
        _add_jpeg_export(plot, parent_widget, jpeg_default_name,
                         x_label=x_label, y_label=y_label, title=title or '')
    return plot


# ===========================================================================
# JPEG / ITX export
# ===========================================================================

def _add_jpeg_export(
    plot: pg.PlotItem,
    parent_widget,
    default_name: str,
    x_label: str = '',
    y_label: str = '',
    title: str = '',
):
    """Append 'Save as JPEG' and 'Save as Igor Pro ITX' entries to the ViewBox menu."""
    vb = plot.getViewBox()
    vb.menu.addSeparator()
    jpeg_action = vb.menu.addAction("Save graph as JPEG…")
    jpeg_action.triggered.connect(
        lambda checked=False, p=plot, pw=parent_widget, n=default_name:
            _save_plot_as_jpeg(p, pw, n)
    )
    itx_action = vb.menu.addAction("Save as Igor Pro ITX…")
    itx_action.triggered.connect(
        lambda checked=False, p=plot, pw=parent_widget, n=default_name,
               xl=x_label, yl=y_label, t=title:
            save_itx_from_plot(p, pw, n, xl, yl, t)
    )


def _save_plot_as_jpeg(plot: pg.PlotItem, parent, default_name: str):
    """Export *plot* to a JPEG file chosen via a save dialog."""
    from pyqtgraph.exporters import ImageExporter
    default_path = str(Path.home() / f'{default_name}.jpg')
    file_path, _ = QFileDialog.getSaveFileName(
        parent, 'Save Graph as JPEG', default_path,
        'JPEG Images (*.jpg *.jpeg);;All Files (*)',
    )
    if not file_path:
        return
    try:
        exporter = ImageExporter(plot)
        exporter.parameters()['width'] = 1600
        exporter.export(file_path)
    except Exception as exc:
        QMessageBox.warning(parent, 'Export Failed',
                            f'Could not save image:\n{exc}')


def _get_item_color_hex(item: pg.PlotDataItem) -> str:
    """Extract the primary color of a PlotDataItem as an #rrggbb hex string."""
    try:
        pen = item.opts.get('pen')
        if pen is not None:
            return pg.mkPen(pen).color().name()
    except Exception:
        pass
    try:
        brush = item.opts.get('symbolBrush')
        if brush is not None:
            return pg.mkBrush(brush).color().name()
    except Exception:
        pass
    return '#000000'


def save_itx_from_plot(
    plot: pg.PlotItem,
    parent,
    default_name: str = 'pyirena_graph',
    x_label: str | None = None,
    y_label: str | None = None,
    title: str | None = None,
) -> None:
    """Export named data curves from *plot* as an Igor Pro Text (.itx) file.

    Iterates all named ``PlotDataItem`` objects in *plot* (scatter data and
    model curves) and writes them as Igor Pro waves with display, log-axis,
    color, label, and legend commands.  Error-bar segments (NaN-separated
    lines without a name) are automatically skipped.

    Auto-extracts axis labels and title from *plot* when the corresponding
    parameters are ``None``.
    """
    import re

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

    # Collect named, non-NaN-heavy PlotDataItems.
    named_items: list[tuple[str, np.ndarray, np.ndarray, str]] = []
    for item in plot.listDataItems():
        if not isinstance(item, pg.PlotDataItem):
            continue
        name = item.name() or ''
        if not name:
            continue
        # Use getOriginalDataset() to get the pre-log-transform linear values.
        # getData() returns log10-transformed values when the item has logMode
        # active, which would cause a double-log in Igor Pro when combined with
        # the ModifyGraph log=1 commands we also emit.
        x_data, y_data = item.getOriginalDataset()
        if x_data is None or y_data is None or len(x_data) < 2:
            continue
        # Skip error-bar segments (NaN-separated lines: >30 % NaN values)
        if np.sum(np.isnan(y_data)) > len(y_data) * 0.3:
            continue
        mask = np.isfinite(x_data) & np.isfinite(y_data)
        if mask.sum() < 2:
            continue
        named_items.append((name, x_data[mask], y_data[mask],
                             _get_item_color_hex(item)))

    if not named_items:
        QMessageBox.warning(parent, 'No data',
                            'No named data curves found to export.')
        return

    default_path = str(Path.home() / f'{default_name}.itx')
    filepath, _ = QFileDialog.getSaveFileName(
        parent, 'Save as Igor Pro ITX', default_path,
        'Igor Pro Text (*.itx);;All files (*)',
    )
    if not filepath:
        return
    if not filepath.lower().endswith('.itx'):
        filepath += '.itx'

    log_x = bool(getattr(plot.getAxis('bottom'), 'logMode', False))
    log_y = bool(getattr(plot.getAxis('left'),   'logMode', False))

    def _safe_name(s: str) -> str:
        n = re.sub(r'[^A-Za-z0-9_]', '_', s)
        if n and n[0].isdigit():
            n = 'w_' + n
        return n[:31] or 'wave'

    def _hex_to_igor(h: str) -> tuple[int, int, int]:
        h = h.lstrip('#')
        if len(h) == 6:
            r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
        else:
            r = g = b = 0
        return r * 257, g * 257, b * 257

    lines = ['IGOR']
    wave_info: list[tuple[str, str, str, str]] = []   # (xn, yn, label, color)
    n = len(named_items)

    for i, (lbl, x_arr, y_arr, color) in enumerate(named_items):
        suffix = f'_{i + 1:02d}' if n > 1 else ''
        xn = _safe_name(f'X_{lbl}{suffix}')
        yn = _safe_name(f'Y_{lbl}{suffix}')
        wave_info.append((xn, yn, lbl, color))

        lines += [f'WAVES/D  {xn}', 'BEGIN']
        lines += [f'  {v:.10g}' for v in x_arr]
        lines += ['END', f'WAVES/D  {yn}', 'BEGIN']
        lines += [f'  {v:.10g}' for v in y_arr]
        lines.append('END')

    lines.append('')
    for j, (xn, yn, lbl, _) in enumerate(wave_info):
        if j == 0:
            lines.append(f'X Display {yn} vs {xn} as "{lbl}"')
        else:
            lines.append(f'X AppendToGraph {yn} vs {xn}')

    if log_x:
        lines.append('X ModifyGraph log(bottom)=1')
    if log_y:
        lines.append('X ModifyGraph log(left)=1')

    for _, yn, _, color in wave_info:
        r, g, b = _hex_to_igor(color)
        lines.append(f'X ModifyGraph rgb({yn})=({r},{g},{b})')

    if x_label:
        lines.append(f'X Label bottom "{x_label}"')
    if y_label:
        lines.append(f'X Label left "{y_label}"')
    if title:
        lines.append(f'X TextBox/C/N=title0/A=MC/X=0/Y=5 "{title}"')

    legend_parts = [f'\\\\s({yn}) {lbl}' for _, yn, lbl, _ in wave_info]
    if legend_parts:
        legend_text = '\\r'.join(legend_parts)
        lines.append(f'X Legend/C/N=text0 "{legend_text}"')

    try:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines) + '\n')
    except Exception as exc:
        QMessageBox.warning(parent, 'Export Failed',
                            f'Could not save file:\n{exc}')


# ===========================================================================
# Data plotting
# ===========================================================================

def plot_iq_data(
    plot: pg.PlotItem,
    q: np.ndarray,
    I: np.ndarray,
    dI: np.ndarray | None = None,
    label: str = 'Data',
) -> tuple:
    """Add a data scatter and optional error bars to *plot*.

    The plot must have been created with log mode active (via
    ``make_sas_plot``).  Pass **physical (linear)** values for *q* and *I* —
    pyqtgraph applies the log10 transform automatically.

    Parameters
    ----------
    plot : pg.PlotItem
        Target plot (log-log or log-x).
    q, I : array
        Scattering vector and intensity.  Only finite, positive values are
        plotted.
    dI : array or None
        Measurement uncertainties.  When provided, error bars are drawn as
        NaN-separated line segments (vertical bar + top/bottom end caps).
    label : str
        Legend label for the data scatter item.

    Returns
    -------
    (scatter_item, error_item)
        *error_item* is ``None`` when *dI* is ``None``.
    """
    q = np.asarray(q, dtype=float)
    I = np.asarray(I, dtype=float)
    mask = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
    q_, I_ = q[mask], I[mask]

    # Use plot.plot(symbol='o') — PlotDataItem — rather than ScatterPlotItem.
    # PlotDataItem pre-transforms data to log10 when setLogMode is active, giving
    # a correct bounding rect for auto-range.  Standalone ScatterPlotItem does NOT
    # pre-transform data; pyqtgraph's auto-range then uses raw linear Q values as
    # ViewBox (log10-space) coordinates, placing scatter entirely outside the visible
    # window and making it invisible.
    scatter = plot.plot(
        q_, I_,
        pen=None,
        symbol='o', symbolSize=SASPlotStyle.DATA_SIZE,
        symbolBrush=SASPlotStyle.DATA_BRUSH, symbolPen=pg.mkPen(None),
        name=label,
    )

    error_item = None
    if dI is not None:
        dI_ = np.asarray(dI, dtype=float)
        if dI_.shape == q.shape:
            dI_ = dI_[mask]
        elif dI_.shape == q_.shape:
            pass   # already masked
        else:
            dI_ = np.zeros_like(q_)
        # Global cap: prevent error bar tops from exceeding 3 decades above the 99th
        # percentile of the data. Without this, a few extreme outlier data points
        # (e.g. uncalibrated WAXS or Bragg peaks with very large I) drive the error
        # bar caps to astronomically large values (10^38+) visible on extreme zoom-out.
        valid_I = I_[I_ > 0]
        if len(valid_I) >= 5:
            y_global_max = 10.0 ** (float(np.percentile(np.log10(valid_I), 99)) + 3)
        else:
            y_global_max = None
        error_item = _draw_error_bars(plot, q_, I_, dI_, y_global_max=y_global_max)

    # Override auto-range with percentile-based bounds from the data only.
    # Without this, error bar segments are included in pyqtgraph's bounds
    # computation and, under log mode, the y-axis scales to the full
    # double-precision range (~1e-300 … 1e300).
    set_robust_y_range(plot, I_)

    # Constrain x zoom to nearest full decade beyond the data Q range, plus 1 decade slack.
    # e.g. data 0.003–0.8 → hard limits 0.0001–10 (log10 space: -4 to 1).
    # Also explicitly set the visible x range to the data Q span (log10 values = ViewBox
    # coordinates in log mode), matching the pattern used in GraphWindow (data_selector.py).
    valid_q = q_[q_ > 0]
    if len(valid_q) >= 2:
        q_lo = int(np.floor(np.log10(float(valid_q.min())))) - 1
        q_hi = int(np.ceil(np.log10(float(valid_q.max())))) + 1
        plot.getViewBox().setLimits(xMin=q_lo, xMax=q_hi)
        plot.setXRange(
            np.log10(float(valid_q.min())),
            np.log10(float(valid_q.max())),
            padding=0.05,
        )

    return scatter, error_item


def set_robust_y_range(plot: pg.PlotItem, I: np.ndarray) -> None:
    """Set the Y axis range from percentile-based bounds of the actual data.

    CRITICAL — call this after adding error bars.  Without it, pyqtgraph's
    auto-range algorithm includes the error bar line segments in the bounds
    calculation.  When ``setLogMode(y=True)`` is active this can cause the
    axis to scale to the full double-precision range (~1e-300 … 1e300).

    The fix mirrors ``sizes_panel._set_robust_y_range()``: compute log10
    percentile bounds from the *intensity* array only, then call
    ``setYRange(lo, hi, padding=0)``.  In pyqtgraph log mode the ViewBox
    coordinate system is log10 space, so log10 values are passed directly.

    Parameters
    ----------
    plot : pg.PlotItem
        Target plot with ``setLogMode(y=True)`` active.
    I : array
        Linear intensity values (only finite, positive values are used).
    """
    valid = (np.asarray(I) > 0) & np.isfinite(I)
    if np.sum(valid) < 3:
        return
    log_i = np.log10(np.asarray(I)[valid])
    lo = np.percentile(log_i, 2) - 0.5    # 2nd percentile minus half-decade
    hi = np.percentile(log_i, 99) + 0.5   # 99th percentile plus half-decade (excludes extreme outliers)
    plot.setYRange(lo, hi, padding=0)
    # Prevent zooming to extreme y values (e.g. from cosmic rays or uncalibrated data).
    # Allow 3 extra decades of zoom room beyond the data range before hitting a hard limit.
    plot.getViewBox().setLimits(yMin=lo - 3, yMax=hi + 3)


def _draw_error_bars(
    plot: pg.PlotItem,
    q: np.ndarray,
    I: np.ndarray,
    dI: np.ndarray,
    y_global_max: float | None = None,
) -> pg.PlotDataItem | None:
    """Draw I(Q) error bars as NaN-separated line segments.

    Each bar consists of:
    - a vertical line segment from ``I − dI`` to ``I + dI``
    - a top horizontal cap of width ±5 % of Q in log space
    - a bottom horizontal cap of the same width

    Negative lower bounds are clipped to 0.1 % of I to remain positive for
    log-y plots.

    IMPORTANT: Python lists (not numpy arrays) are passed to ``plot.plot()``.
    This matches the pattern from sizes_panel.py and avoids triggering an
    edge case in pyqtgraph's bounds computation with NaN-containing arrays
    under log mode.
    """
    cap = SASPlotStyle.ERROR_CAP_FRAC
    x_lines: list[float] = []
    y_lines: list[float] = []

    for qi, Ii, dIi in zip(q, I, dI):
        if not (np.isfinite(dIi) and dIi > 0):
            continue
        y_top = min(Ii + dIi, Ii * 1000)    # clip to 3 decades above local I
        if y_global_max is not None:
            y_top = min(y_top, y_global_max) # also cap at global 99th-percentile + 3 decades
        y_bot = max(Ii - dIi, Ii * 0.001)   # stay positive for log scale

        # Vertical bar
        x_lines.extend([qi, qi, np.nan])
        y_lines.extend([y_bot, y_top, np.nan])

        # End caps (symmetric in log space)
        xl, xr = qi / (1.0 + cap), qi * (1.0 + cap)
        x_lines.extend([xl, xr, np.nan])
        y_lines.extend([y_top, y_top, np.nan])
        x_lines.extend([xl, xr, np.nan])
        y_lines.extend([y_bot, y_bot, np.nan])

    if not x_lines:
        return None

    # Pass Python lists — NOT np.array() — to match sizes_panel pattern.
    return plot.plot(
        x_lines, y_lines,
        pen=SASPlotStyle.ERROR_PEN,
        connect='finite',
    )


def plot_iq_model(
    plot: pg.PlotItem,
    q: np.ndarray,
    I_model: np.ndarray,
    label: str = 'Model',
) -> pg.PlotDataItem | None:
    """Overlay a model/fit curve on *plot*.

    Parameters
    ----------
    plot : pg.PlotItem
        Target plot (log-log or log-x; must have setLogMode active).
    q, I_model : array
        Q values and corresponding model intensities in **linear** units.
    label : str
        Legend label.

    Returns
    -------
    pg.PlotDataItem or None
        The created item, or ``None`` if fewer than 2 valid points exist.
    """
    q = np.asarray(q, dtype=float)
    I_model = np.asarray(I_model, dtype=float)
    mask = np.isfinite(q) & np.isfinite(I_model) & (q > 0) & (I_model > 0)
    if mask.sum() < 2:
        return None
    return plot.plot(q[mask], I_model[mask],
                     pen=SASPlotStyle.FIT_PEN, name=label)


# ===========================================================================
# Cursor helpers
# ===========================================================================

def make_cursors(
    plot: pg.PlotItem,
    q_min: float,
    q_max: float,
) -> tuple[_SafeInfiniteLine, _SafeInfiniteLine]:
    """Create two movable vertical cursors on *plot*.

    Cursor positions are in **log10 space** because pyqtgraph ``InfiniteLine``
    positions are in ViewBox coordinates, which equal ``log10(data)`` when
    ``setLogMode(x=True)`` is active.

    Parameters
    ----------
    plot : pg.PlotItem
        Target plot.  Must have ``setLogMode(x=True)`` active.
    q_min, q_max : float
        Initial Q range boundaries in **linear** (physical) units.
        Cursors are placed 10 % inside each edge.

    Returns
    -------
    (cursor_a, cursor_b)
        Both are :class:`_SafeInfiniteLine` instances placed in log10 space.
        Use :func:`get_cursor_q_range` to read the current Q range.
    """
    log_min = np.log10(max(float(q_min), 1e-20))
    log_max = np.log10(max(float(q_max), 1e-20))
    span    = log_max - log_min

    cursor_a = _SafeInfiniteLine(
        pos=log_min + 0.1 * span, angle=90, movable=True,
        pen=SASPlotStyle.CURSOR_A_PEN,
        label='A',
        labelOpts={'position': 0.05, 'color': SASPlotStyle.CURSOR_A_COLOR},
    )
    cursor_b = _SafeInfiniteLine(
        pos=log_max - 0.1 * span, angle=90, movable=True,
        pen=SASPlotStyle.CURSOR_B_PEN,
        label='B',
        labelOpts={'position': 0.10, 'color': SASPlotStyle.CURSOR_B_COLOR},
    )
    plot.addItem(cursor_a, ignoreBounds=True)
    plot.addItem(cursor_b, ignoreBounds=True)
    return cursor_a, cursor_b


def get_cursor_q_range(
    cursor_a: _SafeInfiniteLine | None,
    cursor_b: _SafeInfiniteLine | None,
) -> tuple[float | None, float | None]:
    """Return ``(q_min, q_max)`` in **linear** units from cursor positions.

    Parameters
    ----------
    cursor_a, cursor_b : _SafeInfiniteLine or None
        The two cursor lines (as returned by :func:`make_cursors`).
        Returns ``(None, None)`` if either cursor is ``None``.

    Returns
    -------
    (q_min, q_max) : (float, float) or (None, None)
        Q range in linear (physical) units, with ``q_min ≤ q_max``.
    """
    if cursor_a is None or cursor_b is None:
        return None, None
    a = 10.0 ** cursor_a.getPos()[0]
    b = 10.0 ** cursor_b.getPos()[0]
    return (min(a, b), max(a, b))


def set_cursor_q_range(
    plot: pg.PlotItem,
    cursor_a: _SafeInfiniteLine | None,
    cursor_b: _SafeInfiniteLine | None,
    q_min: float,
    q_max: float,
) -> tuple[_SafeInfiniteLine, _SafeInfiniteLine]:
    """Position cursors at *q_min* and *q_max* (linear units).

    If *cursor_a* or *cursor_b* are ``None``, they are created via
    :func:`make_cursors` first.  The cursors (existing or new) are returned.

    Parameters
    ----------
    plot : pg.PlotItem
        Target plot.  Only used when cursors need to be created.
    cursor_a, cursor_b : _SafeInfiniteLine or None
    q_min, q_max : float
        Desired Q range in **linear** units.

    Returns
    -------
    (cursor_a, cursor_b)
    """
    if q_min is not None and q_max is not None and q_min > 0 and q_max > 0:
        if cursor_a is None or cursor_b is None:
            cursor_a, cursor_b = make_cursors(plot, q_min, q_max)
        else:
            cursor_a.setPos(np.log10(q_min))
            cursor_b.setPos(np.log10(q_max))
    return cursor_a, cursor_b


# ===========================================================================
# Annotation helpers
# ===========================================================================

def add_plot_annotation(
    plot: pg.PlotItem,
    text: str,
    corner: str = 'lower_left',
) -> pg.TextItem:
    """Place a text annotation in a corner of *plot*'s visible area.

    Works correctly for plots with ``setLogMode`` active: positions are
    computed from ``viewRange()`` (which returns log10 values) and
    ``TextItem.setPos()`` is called with those log10 values directly (they
    map into the ViewBox coordinate system, which is log10 space).

    Parameters
    ----------
    plot : pg.PlotItem
        Target plot.
    text : str
        Multi-line annotation text.
    corner : {'lower_left', 'upper_left', 'lower_right', 'upper_right'}
        Which corner to anchor the annotation to.

    Returns
    -------
    pg.TextItem
        The added item.  Keep a reference so you can remove it later via
        ``plot.removeItem(item)``.
    """
    vr = plot.viewRange()   # log10 space when setLogMode is active
    dx = vr[0][1] - vr[0][0]
    dy = vr[1][1] - vr[1][0]
    margin_x = 0.02 * dx
    margin_y = 0.03 * dy

    if corner == 'lower_left':
        x = vr[0][0] + margin_x
        y = vr[1][0] + margin_y    # BOTTOM of visible range
        anchor = (0, 1)             # bottom-left of text → text extends upward
    elif corner == 'upper_left':
        x = vr[0][0] + margin_x
        y = vr[1][1] - margin_y    # TOP of visible range
        anchor = (0, 0)             # top-left of text → text extends downward
    elif corner == 'lower_right':
        x = vr[0][1] - margin_x
        y = vr[1][0] + margin_y
        anchor = (1, 1)
    else:  # upper_right
        x = vr[0][1] - margin_x
        y = vr[1][1] - margin_y
        anchor = (1, 0)

    item = pg.TextItem(text=text, color=SASPlotStyle.ANNOT_COLOR, anchor=anchor)
    item.setPos(x, y)    # ViewBox (log10) coordinates
    plot.addItem(item)
    return item
