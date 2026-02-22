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

    # ── Grid ──────────────────────────────────────────────────────────────
    GRID_ALPHA       = 0.3


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
    plot = graphics_layout.addPlot(row=row, col=col)
    plot.setLogMode(x=log_x, y=log_y)
    plot.setLabel('left',   y_label)
    plot.setLabel('bottom', x_label)
    plot.showGrid(x=True, y=True, alpha=SASPlotStyle.GRID_ALPHA)
    plot.getAxis('left').enableAutoSIPrefix(False)
    plot.getAxis('bottom').enableAutoSIPrefix(False)
    if title:
        plot.setTitle(title)
    if x_link is not None:
        plot.setXLink(x_link)
    if parent_widget is not None:
        _add_jpeg_export(plot, parent_widget, jpeg_default_name)
    return plot


# ===========================================================================
# JPEG export
# ===========================================================================

def _add_jpeg_export(
    plot: pg.PlotItem,
    parent_widget,
    default_name: str,
):
    """Append a 'Save graph as JPEG…' entry to the ViewBox right-click menu."""
    vb = plot.getViewBox()
    vb.menu.addSeparator()
    action = vb.menu.addAction("Save graph as JPEG…")
    action.triggered.connect(
        lambda checked=False, p=plot, pw=parent_widget, n=default_name:
            _save_plot_as_jpeg(p, pw, n)
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

    scatter = pg.ScatterPlotItem(
        x=q_, y=I_,
        pen=None,
        brush=SASPlotStyle.DATA_BRUSH,
        size=SASPlotStyle.DATA_SIZE,
        name=label,
    )
    plot.addItem(scatter)

    error_item = None
    if dI is not None:
        dI_ = np.asarray(dI, dtype=float)
        if dI_.shape == q.shape:
            dI_ = dI_[mask]
        elif dI_.shape == q_.shape:
            pass   # already masked
        else:
            dI_ = np.zeros_like(q_)
        error_item = _draw_error_bars(plot, q_, I_, dI_)

    # Override auto-range with percentile-based bounds from the data only.
    # Without this, error bar segments are included in pyqtgraph's bounds
    # computation and, under log mode, the y-axis scales to the full
    # double-precision range (~1e-300 … 1e300).
    set_robust_y_range(plot, I_)

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
    hi = np.percentile(log_i, 100) + 0.5  # max plus half-decade
    plot.setYRange(lo, hi, padding=0)


def _draw_error_bars(
    plot: pg.PlotItem,
    q: np.ndarray,
    I: np.ndarray,
    dI: np.ndarray,
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
        y_top = Ii + dIi
        y_bot = max(Ii - dIi, Ii * 0.001)    # stay positive for log scale

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
