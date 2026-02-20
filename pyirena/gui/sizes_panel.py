"""
Sizes Distribution GUI panel for pyIrena.

Provides SizesFitGraphWindow (3-panel pyqtgraph plot: I(Q), residuals, P(r))
and SizesFitPanel (controls + graph) for interactive particle size distribution
fitting from SAS data using the SizesDistribution core class.
"""

try:
    from PySide6.QtWidgets import (
        QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
        QLabel, QLineEdit, QComboBox, QCheckBox, QSpinBox, QSplitter,
        QMessageBox, QScrollArea, QGroupBox, QSizePolicy, QFrame, QTextEdit,
        QTabWidget, QFileDialog,
    )
    from PySide6.QtCore import Qt, Signal
    from PySide6.QtGui import QDoubleValidator, QIntValidator, QBrush, QColor
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
            QLabel, QLineEdit, QComboBox, QCheckBox, QSpinBox, QSplitter,
            QMessageBox, QScrollArea, QGroupBox, QSizePolicy, QFrame, QTextEdit,
            QTabWidget, QFileDialog,
        )
        from PyQt6.QtCore import Qt, Signal
        from PyQt6.QtGui import QDoubleValidator, QIntValidator, QBrush, QColor
    except ImportError:
        raise ImportError("Neither PySide6 nor PyQt6 found. Install with: pip install PySide6")

import numpy as np
from pathlib import Path

import pyqtgraph as pg

from pyirena.core.sizes import SizesDistribution
from pyirena.state.state_manager import StateManager


# ─────────────────────────────────────────────────────────────────────────────
# ScrubbableLineEdit  (mouse-wheel-enabled numeric input, same as unified_fit)
# ─────────────────────────────────────────────────────────────────────────────

class ScrubbableLineEdit(QLineEdit):
    """A QLineEdit that responds to mouse wheel to increment/decrement values.

    Parameters
    ----------
    step_factor : float
        Fraction of the current magnitude used as step size per wheel tick.
        Default 0.1 (10% of magnitude).  Use a smaller value such as 0.02
        for fields that need fine-grained control (e.g. error_scale).
    """

    def __init__(self, text="", parent=None, step_factor=0.1):
        super().__init__(text, parent)
        self.setFocusPolicy(Qt.FocusPolicy.WheelFocus)
        self._step_factor = step_factor

    def wheelEvent(self, event):
        try:
            text = self.text().strip()
            if not text:
                return
            value = float(text)
            delta = event.angleDelta().y()
            if delta == 0:
                return
            direction = 1 if delta > 0 else -1
            if abs(value) < 1e-300:
                value = direction * 1e-10
            else:
                import math
                magnitude = 10 ** math.floor(math.log10(abs(value)))
                step = magnitude * self._step_factor
                value += direction * step
            # Format nicely
            if abs(value) < 0.01 or abs(value) >= 10000:
                self.setText(f"{value:.3e}")
            else:
                self.setText(f"{value:.4g}")
            self.editingFinished.emit()
            event.accept()   # stop the enclosing QScrollArea from also scrolling
        except (ValueError, OverflowError):
            pass


# ─────────────────────────────────────────────────────────────────────────────
# LogDecadeAxis — log-scale axis that shows only decade (10^n) labels
# ─────────────────────────────────────────────────────────────────────────────

class LogDecadeAxis(pg.AxisItem):
    """
    AxisItem for log-mode plots that:
    - Labels only decade positions (1e-3, 0.01, 0.1, 1, 10, 100, 1e3 …)
    - Shows minor ticks at ×2 and ×5 per decade (no labels)
    """

    def tickValues(self, minVal, maxVal, size):
        if not self.logMode:
            return super().tickValues(minVal, maxVal, size)
        import math
        lo = math.floor(minVal - 0.001)
        hi = math.ceil(maxVal + 0.001)
        # Major ticks at every decade
        major = [float(v) for v in range(lo, hi + 1)
                 if minVal - 0.01 <= float(v) <= maxVal + 0.01]
        # Minor ticks at ×2 (log10≈0.301) and ×5 (log10≈0.699) per decade
        minor = []
        for decade in range(lo, hi + 1):
            for offset in (0.301, 0.699):
                v = float(decade) + offset
                if minVal - 0.01 <= v <= maxVal + 0.01:
                    minor.append(v)
        return [(1.0, major), (0.301, minor)]

    def tickStrings(self, values, scale, spacing):
        if not self.logMode:
            return super().tickStrings(values, scale, spacing)
        strings = []
        for v in values:
            if abs(v - round(v)) < 0.05:          # decade tick → show label
                pwr = int(round(v))
                val = 10.0 ** pwr
                if 0.001 <= abs(val) <= 999:
                    strings.append(f'{val:g}')    # 0.001 … 999
                else:
                    strings.append(f'1e{pwr}')   # 1e-4, 1e3, …
            else:
                strings.append('')               # minor tick → no label
        return strings


# ─────────────────────────────────────────────────────────────────────────────
# SizesFitGraphWindow
# ─────────────────────────────────────────────────────────────────────────────

class SizesFitGraphWindow(QWidget):
    """
    pyqtgraph window with three stacked plots:
      - Row 0 (50%): I(Q) vs Q, log-log, with two draggable cursors
      - Row 1 (12%): Residuals vs Q, log-x linear-y
      - Row 2 (38%): Size distribution P(r) vs r, linear-linear
    Plus a status message label at the bottom.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.data_folder = None
        self._data_item = None
        self._error_item = None
        self._fit_item = None
        self._corrected_item = None    # I(Q) minus complex background
        self._complex_bg_item = None   # Complex background model curve
        self._resid_item = None
        self._dist_item = None
        self._dist_unc_item = None     # Distribution uncertainty error bars
        self._cursor_a = None
        self._cursor_b = None
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(2)

        # ── pyqtgraph widget ─────────────────────────────────────────────────
        self.graphics_layout = pg.GraphicsLayoutWidget()
        self.graphics_layout.setBackground('w')

        # Create three plot areas with custom decade-only log axes
        self.main_plot = self.graphics_layout.addPlot(
            row=0, col=0,
            axisItems={
                'left':   LogDecadeAxis(orientation='left'),
                'bottom': LogDecadeAxis(orientation='bottom'),
            }
        )
        self.residuals_plot = self.graphics_layout.addPlot(
            row=1, col=0,
            axisItems={
                'bottom': LogDecadeAxis(orientation='bottom'),
            }
        )
        self.distribution_plot = self.graphics_layout.addPlot(
            row=2, col=0,
            axisItems={
                'bottom': LogDecadeAxis(orientation='bottom'),
            }
        )

        # Height ratios  5 : 1 : 4  ≈ 50% : 10% : 40%
        gl = self.graphics_layout.ci.layout
        gl.setRowStretchFactor(0, 5)
        gl.setRowStretchFactor(1, 1)
        gl.setRowStretchFactor(2, 4)

        # Link residuals X-axis to main plot so zooming/panning stays in sync
        self.residuals_plot.setXLink(self.main_plot)

        # ── Main I(Q) plot  (log-log) ────────────────────────────────────────
        self.main_plot.setLogMode(x=True, y=True)
        # Use units IN the label string — avoids pyqtgraph's SI-prefix auto-scaling
        self.main_plot.setLabel('left',   'I  (cm⁻¹)')
        self.main_plot.setLabel('bottom', 'Q  (Å⁻¹)')
        self.main_plot.showGrid(x=True, y=True, alpha=0.3)
        self.main_plot.addLegend(offset=(10, 10), labelTextSize='14pt', labelTextColor='k')
        _style_axes(self.main_plot)
        self.main_plot.getAxis('left').setWidth(65)

        # ── Residuals plot  (log-x, linear-y) ───────────────────────────────
        self.residuals_plot.setLogMode(x=True, y=False)
        self.residuals_plot.setLabel('left',   'Residuals  (σ)')
        self.residuals_plot.setLabel('bottom', 'Q  (Å⁻¹)')
        self.residuals_plot.showGrid(x=True, y=True, alpha=0.3)
        self._zero_line = pg.InfiniteLine(
            pos=0, angle=0,
            pen=pg.mkPen('k', width=1, style=Qt.PenStyle.DashLine)
        )
        self.residuals_plot.addItem(self._zero_line)
        _style_axes(self.residuals_plot)
        self.residuals_plot.getAxis('left').setWidth(65)

        # ── Distribution plot  (log-x, linear-y — r in Å) ───────────────────
        self.distribution_plot.setLogMode(x=True, y=False)
        self.distribution_plot.setLabel('left',   'P(r)  (a.u.)')
        self.distribution_plot.setLabel('bottom', 'r  (Å)')
        self.distribution_plot.showGrid(x=True, y=True, alpha=0.3)
        self.distribution_plot.addLegend(offset=(10, 10), labelTextSize='14pt', labelTextColor='k')
        _style_axes(self.distribution_plot)
        self.distribution_plot.getAxis('left').setWidth(65)

        layout.addWidget(self.graphics_layout)

        # ── "Save graph as JPEG" in right-click context menu of each plot ────
        for _plot in (self.main_plot, self.distribution_plot):
            _vb = _plot.getViewBox()
            _vb.menu.addSeparator()
            _action = _vb.menu.addAction("Save graph as JPEG…")
            _action.triggered.connect(
                lambda checked=False, p=_plot: self._save_plot_as_jpeg(p)
            )

        # ── Status message (QTextEdit so the user can select/copy text) ─────────
        self.status_message = QTextEdit("")
        self.status_message.setReadOnly(True)
        self.status_message.setMaximumHeight(60)
        self.status_message.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.status_message.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.status_message.setStyleSheet("""
            QTextEdit {
                color: #2c3e50;
                background-color: #f8f9fa;
                border: 1px solid #dee2e6;
                border-radius: 3px;
                padding: 4px 8px;
                font-size: 11px;
            }
        """)
        layout.addWidget(self.status_message)

        self.setLayout(layout)

    # ── JPEG export ──────────────────────────────────────────────────────────

    def _save_plot_as_jpeg(self, plot):
        """Export a plot panel to a JPEG file chosen via file dialog."""
        from pyqtgraph.exporters import ImageExporter
        default_dir = self.data_folder or str(Path.home())
        default_name = str(Path(default_dir) / "sizes_graph.jpg")
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Graph as JPEG",
            default_name,
            "JPEG Images (*.jpg *.jpeg);;All Files (*)",
        )
        if not file_path:
            return
        try:
            exporter = ImageExporter(plot)
            exporter.parameters()['width'] = 1600
            exporter.export(file_path)
        except Exception as exc:
            QMessageBox.warning(self, "Export Failed", f"Could not save image:\n{exc}")

    # ── Cursor management ────────────────────────────────────────────────────

    def _ensure_cursors(self, q_array=None):
        """Add draggable vertical cursor lines to main plot if not yet done."""
        if self._cursor_a is not None:
            return
        # Default positions: 10th and 90th percentile of Q
        if q_array is not None and len(q_array) > 2:
            lo = float(np.log10(q_array[len(q_array) // 10]))
            hi = float(np.log10(q_array[9 * len(q_array) // 10]))
        else:
            lo, hi = -3.0, -0.5   # log10 values in log mode

        self._cursor_a = pg.InfiniteLine(
            pos=lo, angle=90, movable=True,
            pen=pg.mkPen('#e74c3c', width=2),
            label='A', labelOpts={'color': '#e74c3c', 'position': 0.92}
        )
        self._cursor_b = pg.InfiniteLine(
            pos=hi, angle=90, movable=True,
            pen=pg.mkPen('#3498db', width=2),
            label='B', labelOpts={'color': '#3498db', 'position': 0.85}
        )
        self.main_plot.addItem(self._cursor_a)
        self.main_plot.addItem(self._cursor_b)

    def get_cursor_range(self):
        """Return (q_min, q_max) from cursors, or None if not initialized.
        Values are in linear Q space (cursors live in log10 space)."""
        if self._cursor_a is None or self._cursor_b is None:
            return None
        a = float(self._cursor_a.value())
        b = float(self._cursor_b.value())
        q_min = 10 ** min(a, b)
        q_max = 10 ** max(a, b)
        return q_min, q_max

    def set_cursor_range(self, q_min: float, q_max: float):
        """Set cursor positions from linear Q values.
        Must be called after data has been loaded (cursors must exist)."""
        if self._cursor_a is not None and q_min > 0 and q_max > 0:
            lo = np.log10(min(q_min, q_max))
            hi = np.log10(max(q_min, q_max))
            self._cursor_a.setValue(lo)
            self._cursor_b.setValue(hi)

    # ── Plot methods ─────────────────────────────────────────────────────────

    def init_plots(self):
        """Clear all plots completely, then restore persistent overlays."""
        # Clear legends before clearing items (prevents duplicate legend entries)
        for plot in (self.main_plot, self.distribution_plot):
            if getattr(plot, 'legend', None) is not None:
                plot.legend.clear()

        # Remove all plot items (fastest + most reliable approach)
        self.main_plot.clear()
        self.residuals_plot.clear()
        self.distribution_plot.clear()

        # Restore zero reference line on residuals
        self._zero_line = pg.InfiniteLine(
            pos=0, angle=0,
            pen=pg.mkPen('k', width=1, style=Qt.PenStyle.DashLine)
        )
        self.residuals_plot.addItem(self._zero_line)

        # Restore cursors (they are persistent objects; clear() removed them)
        if self._cursor_a is not None:
            self.main_plot.addItem(self._cursor_a)
        if self._cursor_b is not None:
            self.main_plot.addItem(self._cursor_b)

        # Reset tracked item references
        self._data_item = None
        self._error_item = None
        self._fit_item = None
        self._corrected_item = None
        self._complex_bg_item = None
        self._resid_item = None
        self._dist_item = None
        self._dist_unc_item = None

    def plot_data(self, q, intensity, error=None, label='Data'):
        """Plot I(Q) data on the main plot with optional error bars.

        Pass linear Q and intensity values — pyqtgraph applies log10 internally
        when the plot is in log mode.  Robust Y-range is set automatically to
        avoid a single bad low-intensity point compressing the view.
        """
        q_ = np.asarray(q, dtype=float)
        I_ = np.asarray(intensity, dtype=float)

        # Error bars drawn first so data symbols render on top.
        # Use the same NaN-separated line-segment approach as unified_fit.py:
        # pass raw (linear) data coordinates to a plot() call — pyqtgraph
        # applies log10 internally, so the vertical bars are correct in log-log.
        if error is not None:
            err_ = np.asarray(error, dtype=float)
            cap_width_log = 0.05   # 5% multiplicative cap half-width in log space

            x_lines = []
            y_lines = []
            for i in range(len(q_)):
                if q_[i] <= 0 or I_[i] <= 0 or err_[i] <= 0:
                    continue
                y_top    = I_[i] + err_[i]
                y_bottom = max(I_[i] - err_[i], I_[i] * 0.001)

                # Vertical bar
                x_lines.extend([q_[i], q_[i], np.nan])
                y_lines.extend([y_bottom, y_top, np.nan])

                # Horizontal caps (symmetric in log space)
                x_left  = q_[i] / (1 + cap_width_log)
                x_right = q_[i] * (1 + cap_width_log)
                x_lines.extend([x_left, x_right, np.nan])
                y_lines.extend([y_top, y_top, np.nan])
                x_lines.extend([x_left, x_right, np.nan])
                y_lines.extend([y_bottom, y_bottom, np.nan])

            if x_lines:
                self._error_item = self.main_plot.plot(
                    x_lines, y_lines,
                    pen=pg.mkPen('#7f8c8d', width=1),
                    connect='finite',
                )

        # Data symbols (pyqtgraph handles log10 transform automatically)
        self._data_item = self.main_plot.plot(
            q_, I_,
            pen=None, symbol='o', symbolSize=4,
            symbolPen=pg.mkPen('#2c3e50', width=1),
            symbolBrush=pg.mkBrush('#2c3e50'),
            name=label,
        )
        self._ensure_cursors(q_)
        self._set_robust_y_range(I_)

    def _set_robust_y_range(self, intensity):
        """Set Y range of main plot to a percentile-based range so that
        isolated bad low-intensity points don't compress the whole view."""
        valid = (np.asarray(intensity) > 0) & np.isfinite(intensity)
        if np.sum(valid) < 3:
            return
        log_i = np.log10(np.asarray(intensity)[valid])
        lo = np.percentile(log_i, 2) - 0.5    # 2nd percentile minus half-decade
        hi = np.percentile(log_i, 100) + 0.5  # max plus half-decade
        # In pyqtgraph log mode the ViewBox coordinate space is log10(data),
        # so we pass log10 values directly to setYRange.
        self.main_plot.setYRange(lo, hi, padding=0)

    def update_error_bars(self, q, intensity, error):
        """Redraw only the error bars without touching data or fit lines.

        Call this when the Error scale changes so the user can immediately see
        what the scaled uncertainties look like before running a new fit.
        """
        # Remove old error bar item
        if self._error_item is not None:
            try:
                self.main_plot.removeItem(self._error_item)
            except Exception:
                pass
            self._error_item = None

        if error is None:
            return

        q_   = np.asarray(q,         dtype=float)
        I_   = np.asarray(intensity,  dtype=float)
        err_ = np.asarray(error,      dtype=float)
        cap_width_log = 0.05

        x_lines, y_lines = [], []
        for i in range(len(q_)):
            if q_[i] <= 0 or I_[i] <= 0 or err_[i] <= 0:
                continue
            y_top    = I_[i] + err_[i]
            y_bottom = max(I_[i] - err_[i], I_[i] * 0.001)
            x_lines.extend([q_[i], q_[i], np.nan])
            y_lines.extend([y_bottom, y_top, np.nan])
            x_left  = q_[i] / (1 + cap_width_log)
            x_right = q_[i] * (1 + cap_width_log)
            x_lines.extend([x_left, x_right, np.nan])
            y_lines.extend([y_top, y_top, np.nan])
            x_lines.extend([x_left, x_right, np.nan])
            y_lines.extend([y_bottom, y_bottom, np.nan])

        if x_lines:
            self._error_item = self.main_plot.plot(
                x_lines, y_lines,
                pen=pg.mkPen('#7f8c8d', width=1),
                connect='finite',
            )

    def plot_fit(self, q, intensity_model, label='Fitted Model'):
        """Plot model I(Q) on main plot."""
        if self._fit_item is not None:
            try:
                self.main_plot.removeItem(self._fit_item)
            except Exception:
                pass
        self._fit_item = self.main_plot.plot(
            q, intensity_model,
            pen=pg.mkPen('#e74c3c', width=4),
            name=label
        )

    def plot_residuals(self, q, residuals):
        """Plot residuals on residuals plot."""
        if self._resid_item is not None:
            try:
                self.residuals_plot.removeItem(self._resid_item)
            except Exception:
                pass
        self._resid_item = self.residuals_plot.plot(
            q, residuals,
            pen=None, symbol='o', symbolSize=3,
            symbolPen=pg.mkPen('#7f8c8d', width=1),
            symbolBrush=pg.mkBrush('#7f8c8d')
        )

    def plot_distribution(self, r, distribution, label='P(r)'):
        """Plot size distribution as a bar chart with green slanted-line fill.

        ``stepMode=True`` (boolean) is used instead of the string ``'center'``
        for compatibility with all pyqtgraph versions: it requires the x array
        to have N+1 elements (bin edges) while y has N elements (bar heights).
        Bin edges are computed as geometric midpoints in log space so that
        log-spaced bins appear as equal-width bars on the log-x axis.
        """
        if self._dist_item is not None:
            try:
                self.distribution_plot.removeItem(self._dist_item)
            except Exception:
                pass
            if (hasattr(self.distribution_plot, 'legend')
                    and self.distribution_plot.legend is not None):
                try:
                    self.distribution_plot.legend.removeItem(self._dist_item)
                except Exception:
                    pass

        r = np.asarray(r, dtype=float).ravel()
        y = np.asarray(distribution, dtype=float).ravel()

        # Build N+1 bin edges in log10 space from N bin centres.
        # PlotCurveItem added via addItem() does NOT inherit the plot's log mode:
        # setLogMode() is only applied to items present when PlotItem.setLogMode()
        # is called, and PlotCurveItem doesn't even have that method in most
        # pyqtgraph versions.  Instead we pass x data already in log10 space,
        # which is exactly what the ViewBox expects when in log-x mode.
        if len(r) >= 2:
            log_r = np.log10(np.maximum(r, 1e-300))
            log_edges = np.empty(len(r) + 1)
            log_edges[1:-1] = 0.5 * (log_r[:-1] + log_r[1:])
            log_edges[0]  = log_r[0]  - 0.5 * (log_r[1]  - log_r[0])
            log_edges[-1] = log_r[-1] + 0.5 * (log_r[-1] - log_r[-2])
        else:
            log_r = np.log10(np.maximum(r, 1e-300))
            log_edges = np.array([log_r[0] - np.log10(1.2), log_r[0] + np.log10(1.2)])

        green_brush = QBrush(QColor(0, 160, 0, 230), Qt.BrushStyle.BDiagPattern)
        self._dist_item = pg.PlotCurveItem(
            log_edges,           # already in log10 space — ViewBox renders directly
            y,
            stepMode=True,       # boolean True: compatible with all pyqtgraph versions
            fillLevel=0,
            brush=green_brush,
            pen=pg.mkPen('#006400', width=1),
        )
        self.distribution_plot.addItem(self._dist_item)

        # ── Axis ranges ──────────────────────────────────────────────────────
        # addItem() does not trigger autoRange; set ranges explicitly.
        # distribution_plot is in log-x mode so setXRange takes log10 values.
        x_lo = float(log_edges[0])
        x_hi = float(log_edges[-1])
        self.distribution_plot.setXRange(x_lo, x_hi, padding=0.02)

        # Y: bars extend from 0 to max(P(r)); show a small margin above.
        y_max = float(np.max(y)) if len(y) > 0 else 1.0
        if y_max <= 0:
            y_max = 1.0
        self.distribution_plot.setYRange(0, y_max * 1.15, padding=0)

        if (hasattr(self.distribution_plot, 'legend')
                and self.distribution_plot.legend is not None):
            self.distribution_plot.legend.addItem(self._dist_item, label)

    def plot_distribution_uncertainty(self, r, dist_mean, dist_std):
        """Overlay ±1σ error bars on the distribution plot.

        Uses NaN-separated line segments.  The distribution plot is in log-x
        mode and ``.plot()`` items get log10 applied automatically by pyqtgraph,
        so we pass LINEAR r values here.
        """
        if self._dist_unc_item is not None:
            try:
                self.distribution_plot.removeItem(self._dist_unc_item)
            except Exception:
                pass
            self._dist_unc_item = None

        r = np.asarray(r, dtype=float).ravel()
        mean = np.asarray(dist_mean, dtype=float).ravel()
        std = np.asarray(dist_std, dtype=float).ravel()

        if len(r) < 2:
            return

        # Cap half-width in linear r space: ~5% of neighbouring spacing
        dr_half = np.zeros(len(r))
        dr_half[0]  = 0.05 * (r[1] - r[0])
        dr_half[-1] = 0.05 * (r[-1] - r[-2])
        dr_half[1:-1] = 0.05 * 0.5 * (r[2:] - r[:-2])

        x_lines, y_lines = [], []
        for i in range(len(r)):
            if std[i] <= 0 or not np.isfinite(mean[i]) or r[i] <= 0:
                continue
            y_top = mean[i] + std[i]
            y_bot = max(0.0, mean[i] - std[i])
            # Vertical bar
            x_lines.extend([r[i], r[i], np.nan])
            y_lines.extend([y_bot, y_top, np.nan])
            # Horizontal caps
            x_lines.extend([r[i] - dr_half[i], r[i] + dr_half[i], np.nan])
            y_lines.extend([y_top, y_top, np.nan])
            x_lines.extend([r[i] - dr_half[i], r[i] + dr_half[i], np.nan])
            y_lines.extend([y_bot, y_bot, np.nan])

        if x_lines:
            self._dist_unc_item = self.distribution_plot.plot(
                x_lines, y_lines,
                pen=pg.mkPen('#2c7a2c', width=1),
                connect='finite',
            )

    def plot_corrected_data(self, q, I_corrected, label='I−bg'):
        """Plot background-corrected data (I − complex background) as blue triangles.

        Replaces any existing corrected-data item.
        """
        if self._corrected_item is not None:
            try:
                self.main_plot.removeItem(self._corrected_item)
            except Exception:
                pass
            if (hasattr(self.main_plot, 'legend')
                    and self.main_plot.legend is not None):
                try:
                    self.main_plot.legend.removeItem(self._corrected_item)
                except Exception:
                    pass
        q_ = np.asarray(q, dtype=float)
        I_ = np.asarray(I_corrected, dtype=float)
        self._corrected_item = self.main_plot.plot(
            q_, I_,
            pen=None, symbol='t', symbolSize=5,
            symbolPen=pg.mkPen('#2980b9', width=1),
            symbolBrush=pg.mkBrush('#2980b9'),
            name=label,
        )

    def plot_complex_background(self, q, bg_values, label='Complex bg'):
        """Plot the complex background model as a dashed grey line.

        Replaces any existing background-model item.
        """
        if self._complex_bg_item is not None:
            try:
                self.main_plot.removeItem(self._complex_bg_item)
            except Exception:
                pass
            if (hasattr(self.main_plot, 'legend')
                    and self.main_plot.legend is not None):
                try:
                    self.main_plot.legend.removeItem(self._complex_bg_item)
                except Exception:
                    pass
        q_ = np.asarray(q, dtype=float)
        bg_ = np.asarray(bg_values, dtype=float)
        # Only plot positive background values (log-log plot requires y > 0)
        valid = (bg_ > 0) & np.isfinite(bg_) & (q_ > 0)
        if np.any(valid):
            self._complex_bg_item = self.main_plot.plot(
                q_[valid], bg_[valid],
                pen=pg.mkPen('#95a5a6', width=1.5,
                             style=Qt.PenStyle.DashLine),
                name=label,
            )

    # ── Status helpers ───────────────────────────────────────────────────────

    def show_message(self, text, color='#2c3e50'):
        print(text)   # always echo to terminal so messages can be copy-pasted there too
        self.status_message.setPlainText(text)
        self.status_message.setStyleSheet(f"""
            QTextEdit {{
                color: {color};
                background-color: #f8f9fa;
                border: 1px solid #dee2e6;
                border-radius: 3px;
                padding: 4px 8px;
                font-size: 11px;
            }}
        """)

    def show_error_message(self, text):
        self.show_message(text, '#c0392b')

    def show_success_message(self, text):
        self.show_message(text, '#27ae60')


def _style_axes(plot_item):
    """Apply white background / black axes styling and disable SI prefix scaling."""
    for side in ('left', 'bottom'):
        ax = plot_item.getAxis(side)
        ax.setPen(pg.mkPen('k'))
        ax.setTextPen(pg.mkPen('k'))
        ax.enableAutoSIPrefix(False)   # prevents '1 GCm⁻¹' style mangling


# ─────────────────────────────────────────────────────────────────────────────
# SizesFitPanel
# ─────────────────────────────────────────────────────────────────────────────

class SizesFitPanel(QWidget):
    """
    Main Sizes Distribution panel for pyIrena.

    Left panel: all model controls (shape, grid, method, results).
    Right panel: SizesFitGraphWindow (3-plot layout).
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.graph_window = None
        self.data = None
        self.fit_result = None
        self._last_distribution = None      # Tuple (r, P) from last fit/compute
        self._param_backup = None           # For revert functionality
        self._pending_cursor_q_range = None # Q range restored after data load

        self.state_manager = StateManager()
        self.sizes = SizesDistribution()

        self.init_ui()
        self.load_state()

    # ── UI construction ──────────────────────────────────────────────────────

    def init_ui(self):
        self.setWindowTitle("pyIrena - Size Distribution")

        main_splitter = QSplitter(Qt.Orientation.Horizontal)

        left_panel = self._create_control_panel()
        main_splitter.addWidget(left_panel)

        self.graph_window = SizesFitGraphWindow()
        main_splitter.addWidget(self.graph_window)

        main_splitter.setSizes([400, 800])
        main_splitter.setStretchFactor(0, 1)
        main_splitter.setStretchFactor(1, 2)
        self.main_splitter = main_splitter

        root_layout = QVBoxLayout()
        root_layout.setContentsMargins(0, 0, 0, 0)
        root_layout.addWidget(main_splitter)
        self.setLayout(root_layout)

        self.setMinimumSize(1200, 960)
        self.resize(1200, 960)

    def _create_control_panel(self) -> QWidget:
        """Build the left control panel (tabbed) inside a scroll area."""
        # Outer container with fixed width
        outer = QWidget()
        outer.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Preferred)
        outer.setMinimumWidth(420)
        outer.setMaximumWidth(420)

        outer_layout = QVBoxLayout(outer)
        outer_layout.setContentsMargins(0, 0, 0, 0)
        outer_layout.setSpacing(0)

        # Scroll area so controls don't get clipped on short screens
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setFrameShape(QFrame.Shape.NoFrame)

        inner = QWidget()
        layout = QVBoxLayout(inner)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)

        # ── Title ────────────────────────────────────────────────────────────
        title = QLabel("Sizes Distribution Input")
        title.setStyleSheet("""
            QLabel {
                font-size: 14px; font-weight: bold;
                color: #2c3e50;
                background-color: #ecf0f1;
                padding: 8px;
                border: 1px solid #bdc3c7;
            }
        """)
        layout.addWidget(title)

        # ── Tab widget ───────────────────────────────────────────────────────
        tabs = QTabWidget()
        tabs.setDocumentMode(True)

        # ════════════════════ SIZES TAB ══════════════════════════════════════
        sizes_tab = QWidget()
        sizes_layout = QVBoxLayout(sizes_tab)
        sizes_layout.setContentsMargins(4, 6, 4, 4)
        sizes_layout.setSpacing(6)

        # ── Shape group ───────────────────────────────────────────────────
        shape_box = QGroupBox("Particle Shape")
        shape_layout = QVBoxLayout(shape_box)
        shape_layout.setSpacing(4)

        row = QHBoxLayout()
        row.addWidget(QLabel("Particle type:"))
        self.shape_combo = QComboBox()
        self.shape_combo.addItems(["sphere", "spheroid"])
        self.shape_combo.currentTextChanged.connect(self._on_shape_changed)
        row.addWidget(self.shape_combo)
        row.addStretch()
        shape_layout.addLayout(row)

        row2 = QHBoxLayout()
        self.aspect_label = QLabel("Aspect ratio:")
        row2.addWidget(self.aspect_label)
        self.aspect_ratio_edit = ScrubbableLineEdit("1.5")
        self.aspect_ratio_edit.setValidator(QDoubleValidator(0.01, 1000.0, 6))
        self.aspect_ratio_edit.setMaximumWidth(90)
        self.aspect_ratio_edit.editingFinished.connect(self._on_param_changed)
        row2.addWidget(self.aspect_ratio_edit)
        row2.addStretch()
        shape_layout.addLayout(row2)

        row3 = QHBoxLayout()
        row3.addWidget(QLabel("Contrast (Δρ²):"))
        self.contrast_edit = ScrubbableLineEdit("1.0")
        self.contrast_edit.setValidator(QDoubleValidator(0.0, 1e10, 6))
        self.contrast_edit.setMaximumWidth(90)
        self.contrast_edit.editingFinished.connect(self._on_param_changed)
        row3.addWidget(self.contrast_edit)
        row3.addWidget(QLabel("×10²⁰ cm⁻⁴"))
        row3.addStretch()
        shape_layout.addLayout(row3)

        sizes_layout.addWidget(shape_box)

        # ── Q range display (from cursors) ────────────────────────────────
        q_box = QGroupBox("Q Range (cursors)")
        q_layout = QVBoxLayout(q_box)
        q_layout.setSpacing(3)

        q_row = QHBoxLayout()
        q_row.addWidget(QLabel("Q min:"))
        self.qmin_display = QLineEdit("—")
        self.qmin_display.setReadOnly(True)
        self.qmin_display.setMaximumWidth(90)
        self.qmin_display.setStyleSheet("background-color: #ecf0f1; color: #7f8c8d;")
        q_row.addWidget(self.qmin_display)
        q_row.addWidget(QLabel("  Q max:"))
        self.qmax_display = QLineEdit("—")
        self.qmax_display.setReadOnly(True)
        self.qmax_display.setMaximumWidth(90)
        self.qmax_display.setStyleSheet("background-color: #ecf0f1; color: #7f8c8d;")
        q_row.addWidget(self.qmax_display)
        q_row.addWidget(QLabel("Å⁻¹"))
        q_row.addStretch()
        q_layout.addLayout(q_row)

        hint = QLabel("Drag cursors on I(Q) graph to set Q range")
        hint.setStyleSheet("color: #7f8c8d; font-size: 10px;")
        q_layout.addWidget(hint)

        sizes_layout.addWidget(q_box)

        # ── Size grid ─────────────────────────────────────────────────────
        grid_box = QGroupBox("Size Grid")
        grid_layout = QVBoxLayout(grid_box)
        grid_layout.setSpacing(4)

        for label_text, attr, default in [
            ("r min [Å]:", "rmin_edit", "10"),
            ("r max [Å]:", "rmax_edit", "1000"),
        ]:
            row = QHBoxLayout()
            row.addWidget(QLabel(label_text))
            edit = ScrubbableLineEdit(default)
            edit.setValidator(QDoubleValidator(0.01, 1e7, 6))
            edit.setMaximumWidth(90)
            edit.editingFinished.connect(self._on_param_changed)
            setattr(self, attr, edit)
            row.addWidget(edit)
            row.addWidget(QLabel("Å"))
            btn_div = QPushButton("÷10")
            btn_mul = QPushButton("×10")
            for btn in (btn_div, btn_mul):
                btn.setMaximumWidth(36)
                btn.setMaximumHeight(22)
                btn.setStyleSheet("font-size:10px; padding:1px 3px;")
            btn_div.clicked.connect(lambda _=False, e=edit: self._scale_edit(e, 0.1))
            btn_mul.clicked.connect(lambda _=False, e=edit: self._scale_edit(e, 10.0))
            row.addWidget(btn_div)
            row.addWidget(btn_mul)
            row.addStretch()
            grid_layout.addLayout(row)

        bins_row = QHBoxLayout()
        bins_row.addWidget(QLabel("Number of bins:"))
        self.nbins_spin = QSpinBox()
        self.nbins_spin.setMinimum(5)
        self.nbins_spin.setMaximum(500)
        self.nbins_spin.setValue(200)
        self.nbins_spin.setMaximumWidth(70)
        self.nbins_spin.valueChanged.connect(self._on_param_changed)
        bins_row.addWidget(self.nbins_spin)
        self.log_spacing_check = QCheckBox("Logarithmic spacing")
        self.log_spacing_check.setChecked(True)
        self.log_spacing_check.stateChanged.connect(self._on_param_changed)
        bins_row.addWidget(self.log_spacing_check)
        bins_row.addStretch()
        grid_layout.addLayout(bins_row)

        sizes_layout.addWidget(grid_box)

        # ── Error scaling ─────────────────────────────────────────────────
        err_box = QGroupBox("Error Scaling")
        err_box_layout = QVBoxLayout(err_box)
        err_row = QHBoxLayout()
        err_row.addWidget(QLabel("Error scale:"))
        self.error_scale_edit = ScrubbableLineEdit("1.0", step_factor=0.02)
        self.error_scale_edit.setValidator(QDoubleValidator(0.001, 1000.0, 6))
        self.error_scale_edit.setMaximumWidth(100)
        self.error_scale_edit.setToolTip(
            "Multiply measurement errors by this factor before fitting.\n"
            "Increase (>1) if MaxEnt cannot find a solution within the\n"
            "current error bars."
        )
        self.error_scale_edit.editingFinished.connect(self._on_param_changed)
        self.error_scale_edit.editingFinished.connect(self._update_error_display)
        err_row.addWidget(self.error_scale_edit)
        err_row.addWidget(QLabel("(default 1.0)"))
        err_row.addStretch()
        err_box_layout.addLayout(err_row)
        sizes_layout.addWidget(err_box)

        # ── Method & parameters ───────────────────────────────────────────
        method_box = QGroupBox("Fitting Method")
        method_layout = QVBoxLayout(method_box)
        method_layout.setSpacing(5)

        meth_row = QHBoxLayout()
        meth_row.addWidget(QLabel("Method:"))
        self.method_combo = QComboBox()
        self.method_combo.addItems(["regularization", "maxent", "tnnls", "mcsas"])
        self.method_combo.currentTextChanged.connect(self._on_method_changed)
        meth_row.addWidget(self.method_combo)
        meth_row.addStretch()
        method_layout.addLayout(meth_row)

        # MaxEnt params — Sky bg and Max iter (Stability is hardcoded at 0.01)
        self.maxent_group = QWidget()
        maxent_layout = QVBoxLayout(self.maxent_group)
        maxent_layout.setContentsMargins(0, 0, 0, 0)
        maxent_layout.setSpacing(2)

        r1 = QHBoxLayout()
        r1.addWidget(QLabel("Sky bg:"))
        self.maxent_sky_edit = ScrubbableLineEdit("1e-06")
        self.maxent_sky_edit.setValidator(QDoubleValidator(0.0, 1e10, 6))
        self.maxent_sky_edit.setMaximumWidth(80)
        self.maxent_sky_edit.editingFinished.connect(self._on_param_changed)
        r1.addWidget(self.maxent_sky_edit)
        r1.addWidget(QLabel("  Max iter:"))
        self.maxent_maxiter_spin = QSpinBox()
        self.maxent_maxiter_spin.setMinimum(10)
        self.maxent_maxiter_spin.setMaximum(10000)
        self.maxent_maxiter_spin.setValue(1000)
        self.maxent_maxiter_spin.setMaximumWidth(70)
        self.maxent_maxiter_spin.valueChanged.connect(self._on_param_changed)
        r1.addWidget(self.maxent_maxiter_spin)
        r1.addStretch()
        maxent_layout.addLayout(r1)
        method_layout.addWidget(self.maxent_group)

        # Regularization params — χ² and Min ratio on one row
        self.reg_group = QWidget()
        reg_layout = QVBoxLayout(self.reg_group)
        reg_layout.setContentsMargins(0, 0, 0, 0)
        reg_layout.setSpacing(2)

        r1 = QHBoxLayout()
        r1.addWidget(QLabel("χ² (evalue):"))
        self.reg_evalue_edit = ScrubbableLineEdit("1.0")
        self.reg_evalue_edit.setValidator(QDoubleValidator(0.0, 1e10, 6))
        self.reg_evalue_edit.setMaximumWidth(70)
        self.reg_evalue_edit.editingFinished.connect(self._on_param_changed)
        r1.addWidget(self.reg_evalue_edit)
        r1.addWidget(QLabel("  Min ratio:"))
        self.reg_minratio_edit = ScrubbableLineEdit("1e-04")
        self.reg_minratio_edit.setValidator(QDoubleValidator(0.0, 1e10, 6))
        self.reg_minratio_edit.setMaximumWidth(70)
        self.reg_minratio_edit.editingFinished.connect(self._on_param_changed)
        r1.addWidget(self.reg_minratio_edit)
        r1.addStretch()
        reg_layout.addLayout(r1)
        method_layout.addWidget(self.reg_group)

        # TNNLS params — Approach and Max iter on one row
        self.tnnls_group = QWidget()
        tnnls_layout = QVBoxLayout(self.tnnls_group)
        tnnls_layout.setContentsMargins(0, 0, 0, 0)
        tnnls_layout.setSpacing(2)

        r1 = QHBoxLayout()
        r1.addWidget(QLabel("Approach:"))
        self.tnnls_approach_edit = ScrubbableLineEdit("0.95")
        self.tnnls_approach_edit.setValidator(QDoubleValidator(0.01, 1.0, 6))
        self.tnnls_approach_edit.setMaximumWidth(70)
        self.tnnls_approach_edit.editingFinished.connect(self._on_param_changed)
        r1.addWidget(self.tnnls_approach_edit)
        r1.addWidget(QLabel("  Max iter:"))
        self.tnnls_maxiter_spin = QSpinBox()
        self.tnnls_maxiter_spin.setMinimum(10)
        self.tnnls_maxiter_spin.setMaximum(100000)
        self.tnnls_maxiter_spin.setValue(1000)
        self.tnnls_maxiter_spin.setMaximumWidth(70)
        self.tnnls_maxiter_spin.valueChanged.connect(self._on_param_changed)
        r1.addWidget(self.tnnls_maxiter_spin)
        r1.addStretch()
        tnnls_layout.addLayout(r1)
        method_layout.addWidget(self.tnnls_group)

        # McSAS params — Convergence, Max iter
        # (N contributions = N size bins, shared with the Size Grid above)
        self.mcsas_group = QWidget()
        mcsas_layout = QVBoxLayout(self.mcsas_group)
        mcsas_layout.setContentsMargins(0, 0, 0, 0)
        mcsas_layout.setSpacing(2)

        r2 = QHBoxLayout()
        r2.addWidget(QLabel("Convergence:"))
        self.mcsas_conv_edit = ScrubbableLineEdit("1.0")
        self.mcsas_conv_edit.setValidator(QDoubleValidator(0.01, 100.0, 6))
        self.mcsas_conv_edit.setMaximumWidth(65)
        self.mcsas_conv_edit.editingFinished.connect(self._on_param_changed)
        r2.addWidget(self.mcsas_conv_edit)
        r2.addWidget(QLabel("  Max iter:"))
        self.mcsas_maxiter_spin = QSpinBox()
        self.mcsas_maxiter_spin.setMinimum(1000)
        self.mcsas_maxiter_spin.setMaximum(10000000)
        self.mcsas_maxiter_spin.setSingleStep(10000)
        self.mcsas_maxiter_spin.setValue(100000)
        self.mcsas_maxiter_spin.setMaximumWidth(80)
        self.mcsas_maxiter_spin.valueChanged.connect(self._on_param_changed)
        r2.addWidget(self.mcsas_maxiter_spin)
        r2.addStretch()
        mcsas_layout.addLayout(r2)
        method_layout.addWidget(self.mcsas_group)

        sizes_layout.addWidget(method_box)
        sizes_layout.addStretch()

        tabs.addTab(sizes_tab, "Sizes")

        # ════════════════════ BACKGROUND TAB ═════════════════════════════════
        bg_tab = QWidget()
        bg_tab_layout = QVBoxLayout(bg_tab)
        bg_tab_layout.setContentsMargins(4, 6, 4, 4)
        bg_tab_layout.setSpacing(8)

        info_lbl = QLabel(
            "Corrected data  I − (B·q⁻ᴾ + background)  is shown in blue on "
            "the I(Q) graph."
        )
        info_lbl.setWordWrap(True)
        info_lbl.setStyleSheet("color: #2980b9; font-size: 10px;")
        bg_tab_layout.addWidget(info_lbl)

        # ── Power-law group ───────────────────────────────────────────────
        pl_box = QGroupBox("Power-Law Term  B·q⁻ᴾ")
        pl_layout = QVBoxLayout(pl_box)
        pl_layout.setSpacing(4)

        # B row
        b_row = QHBoxLayout()
        b_row.addWidget(QLabel("B:"))
        self.power_law_B_edit = ScrubbableLineEdit("0.0")
        self.power_law_B_edit.setValidator(QDoubleValidator(0.0, 1e20, 6))
        self.power_law_B_edit.setMaximumWidth(100)
        self.power_law_B_edit.setToolTip("Amplitude of B·q⁻ᴾ power-law background term")
        self.power_law_B_edit.editingFinished.connect(self._on_param_changed)
        b_row.addWidget(self.power_law_B_edit)
        self.fit_B_check = QCheckBox("Fit B?")
        b_row.addWidget(self.fit_B_check)
        b_row.addStretch()
        pl_layout.addLayout(b_row)

        # P row
        p_row = QHBoxLayout()
        p_row.addWidget(QLabel("P:"))
        self.power_law_P_edit = ScrubbableLineEdit("4.0")
        self.power_law_P_edit.setValidator(QDoubleValidator(0.1, 12.0, 6))
        self.power_law_P_edit.setMaximumWidth(100)
        self.power_law_P_edit.setToolTip("Exponent of B·q⁻ᴾ power-law background term")
        self.power_law_P_edit.editingFinished.connect(self._on_param_changed)
        p_row.addWidget(self.power_law_P_edit)
        self.fit_P_check = QCheckBox("Fit P?")
        p_row.addWidget(self.fit_P_check)
        p_row.addStretch()
        pl_layout.addLayout(p_row)

        # Q range for power-law fit
        pl_q_lbl = QLabel("Q range for fit:")
        pl_q_lbl.setStyleSheet("font-size: 10px; color: #555;")
        pl_layout.addWidget(pl_q_lbl)

        pl_qmin_row = QHBoxLayout()
        pl_qmin_row.addWidget(QLabel("  Q min:"))
        self.pl_q_min_edit = ScrubbableLineEdit("")
        self.pl_q_min_edit.setPlaceholderText("cursor")
        self.pl_q_min_edit.setValidator(QDoubleValidator(0.0, 100.0, 8))
        self.pl_q_min_edit.setMaximumWidth(90)
        pl_qmin_row.addWidget(self.pl_q_min_edit)
        pl_qmin_row.addWidget(QLabel("Å⁻¹"))
        pl_qmin_row.addStretch()
        pl_layout.addLayout(pl_qmin_row)

        pl_qmax_row = QHBoxLayout()
        pl_qmax_row.addWidget(QLabel("  Q max:"))
        self.pl_q_max_edit = ScrubbableLineEdit("")
        self.pl_q_max_edit.setPlaceholderText("cursor")
        self.pl_q_max_edit.setValidator(QDoubleValidator(0.0, 100.0, 8))
        self.pl_q_max_edit.setMaximumWidth(90)
        pl_qmax_row.addWidget(self.pl_q_max_edit)
        pl_qmax_row.addWidget(QLabel("Å⁻¹"))
        pl_qmax_row.addStretch()
        pl_layout.addLayout(pl_qmax_row)

        pl_cursor_row = QHBoxLayout()
        pl_cursor_btn = QPushButton("Set Q from cursors")
        pl_cursor_btn.setMaximumHeight(22)
        pl_cursor_btn.setStyleSheet("font-size:10px; padding:1px 4px;")
        pl_cursor_btn.clicked.connect(self._set_pl_q_from_cursors)
        pl_cursor_row.addWidget(pl_cursor_btn)
        pl_cursor_row.addStretch()
        pl_layout.addLayout(pl_cursor_row)

        self.fit_pl_button = QPushButton("Fit P/B")
        self.fit_pl_button.setMinimumHeight(26)
        self.fit_pl_button.setStyleSheet("""
            QPushButton { background-color: #8e44ad; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #7d3c98; }
        """)
        self.fit_pl_button.clicked.connect(self.fit_power_law_action)
        pl_layout.addWidget(self.fit_pl_button)

        bg_tab_layout.addWidget(pl_box)

        # ── Flat background group ─────────────────────────────────────────
        flat_box = QGroupBox("Flat Background (constant)")
        flat_layout = QVBoxLayout(flat_box)
        flat_layout.setSpacing(4)

        bg_val_row = QHBoxLayout()
        bg_val_row.addWidget(QLabel("Background:"))
        self.background_edit = ScrubbableLineEdit("0.0")
        self.background_edit.setValidator(QDoubleValidator(-1e10, 1e10, 6))
        self.background_edit.setMaximumWidth(100)
        self.background_edit.editingFinished.connect(self._on_param_changed)
        bg_val_row.addWidget(self.background_edit)
        bg_val_row.addWidget(QLabel("cm⁻¹"))
        bg_val_row.addStretch()
        flat_layout.addLayout(bg_val_row)

        bg_q_lbl = QLabel("Q range for fit:")
        bg_q_lbl.setStyleSheet("font-size: 10px; color: #555;")
        flat_layout.addWidget(bg_q_lbl)

        bg_qmin_row = QHBoxLayout()
        bg_qmin_row.addWidget(QLabel("  Q min:"))
        self.bg_q_min_edit = ScrubbableLineEdit("")
        self.bg_q_min_edit.setPlaceholderText("cursor")
        self.bg_q_min_edit.setValidator(QDoubleValidator(0.0, 100.0, 8))
        self.bg_q_min_edit.setMaximumWidth(90)
        bg_qmin_row.addWidget(self.bg_q_min_edit)
        bg_qmin_row.addWidget(QLabel("Å⁻¹"))
        bg_qmin_row.addStretch()
        flat_layout.addLayout(bg_qmin_row)

        bg_qmax_row = QHBoxLayout()
        bg_qmax_row.addWidget(QLabel("  Q max:"))
        self.bg_q_max_edit = ScrubbableLineEdit("")
        self.bg_q_max_edit.setPlaceholderText("cursor")
        self.bg_q_max_edit.setValidator(QDoubleValidator(0.0, 100.0, 8))
        self.bg_q_max_edit.setMaximumWidth(90)
        bg_qmax_row.addWidget(self.bg_q_max_edit)
        bg_qmax_row.addWidget(QLabel("Å⁻¹"))
        bg_qmax_row.addStretch()
        flat_layout.addLayout(bg_qmax_row)

        bg_cursor_row = QHBoxLayout()
        bg_cursor_btn = QPushButton("Set Q from cursors")
        bg_cursor_btn.setMaximumHeight(22)
        bg_cursor_btn.setStyleSheet("font-size:10px; padding:1px 4px;")
        bg_cursor_btn.clicked.connect(self._set_bg_q_from_cursors)
        bg_cursor_row.addWidget(bg_cursor_btn)
        bg_cursor_row.addStretch()
        flat_layout.addLayout(bg_cursor_row)

        self.fit_bg_button = QPushButton("Fit Background")
        self.fit_bg_button.setMinimumHeight(26)
        self.fit_bg_button.setStyleSheet("""
            QPushButton { background-color: #8e44ad; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #7d3c98; }
        """)
        self.fit_bg_button.clicked.connect(self.fit_background_action)
        flat_layout.addWidget(self.fit_bg_button)

        bg_tab_layout.addWidget(flat_box)
        bg_tab_layout.addStretch()

        tabs.addTab(bg_tab, "Background")

        layout.addWidget(tabs)

        # ── Action buttons (always visible outside tabs) ──────────────────────
        btn_row1 = QHBoxLayout()

        self.graph_button = QPushButton("Graph Model")
        self.graph_button.setMinimumHeight(28)
        self.graph_button.setStyleSheet("""
            QPushButton { background-color: #52c77a; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #3eb56a; }
        """)
        self.graph_button.clicked.connect(self.compute_model)
        btn_row1.addWidget(self.graph_button)

        self.fit_button = QPushButton("Fit Sizes")
        self.fit_button.setMinimumHeight(28)
        self.fit_button.setStyleSheet("""
            QPushButton { background-color: #27ae60; color: white; font-weight: bold; font-size: 12px; }
            QPushButton:hover { background-color: #1e8449; }
        """)
        self.fit_button.clicked.connect(self.run_fit)
        btn_row1.addWidget(self.fit_button)

        self.fit_all_button = QPushButton("Fit All")
        self.fit_all_button.setMinimumHeight(28)
        self.fit_all_button.setStyleSheet("""
            QPushButton { background-color: #2980b9; color: white; font-weight: bold; font-size: 12px; }
            QPushButton:hover { background-color: #1f618d; }
        """)
        self.fit_all_button.clicked.connect(self.fit_all_action)
        btn_row1.addWidget(self.fit_all_button)

        layout.addLayout(btn_row1)

        btn_row_unc = QHBoxLayout()
        btn_row_unc.addWidget(QLabel("Passes:"))
        self.unc_n_runs_spin = QSpinBox()
        self.unc_n_runs_spin.setMinimum(1)
        self.unc_n_runs_spin.setMaximum(200)
        self.unc_n_runs_spin.setValue(10)
        self.unc_n_runs_spin.setMaximumWidth(52)
        self.unc_n_runs_spin.setToolTip("Number of noise-perturbed fits for uncertainty estimation")
        btn_row_unc.addWidget(self.unc_n_runs_spin)
        self.uncertainty_button = QPushButton("Calc. Uncertainty (MC)")
        self.uncertainty_button.setMinimumHeight(26)
        self.uncertainty_button.setStyleSheet("""
            QPushButton { background-color: #16a085; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #1abc9c; }
        """)
        self.uncertainty_button.clicked.connect(self.calculate_uncertainty)
        btn_row_unc.addWidget(self.uncertainty_button)
        layout.addLayout(btn_row_unc)

        btn_row2 = QHBoxLayout()

        self.revert_button = QPushButton("Revert")
        self.revert_button.setMinimumHeight(26)
        self.revert_button.setStyleSheet("""
            QPushButton { background-color: #e67e22; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #f39c12; }
        """)
        self.revert_button.clicked.connect(self._revert)
        btn_row2.addWidget(self.revert_button)

        self.reset_button = QPushButton("Reset to Defaults")
        self.reset_button.setMinimumHeight(26)
        self.reset_button.setStyleSheet("""
            QPushButton { background-color: #e67e22; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #d35400; }
        """)
        self.reset_button.clicked.connect(self.reset_to_defaults)
        btn_row2.addWidget(self.reset_button)

        layout.addLayout(btn_row2)

        # ── Results section ──────────────────────────────────────────────────
        results_header = QLabel("Results")
        results_header.setStyleSheet("""
            QLabel { font-weight: bold; color: #3498db; font-size: 12px; margin-top: 4px; }
        """)
        layout.addWidget(results_header)

        results_box = QGroupBox("")
        res_layout = QVBoxLayout(results_box)
        res_layout.setSpacing(2)

        # Row 1: χ² and Volume fraction
        chi_vf_row = QHBoxLayout()
        chi_vf_row.addWidget(QLabel("χ²:"))
        self.result_chi2 = QLineEdit("—")
        self.result_chi2.setReadOnly(True)
        self.result_chi2.setStyleSheet(
            "background-color: #ecf0f1; color: #2c3e50; font-weight: bold;")
        chi_vf_row.addWidget(self.result_chi2)
        self.result_chi2_target = QLineEdit("")
        self.result_chi2_target.setReadOnly(True)
        self.result_chi2_target.setMaximumWidth(78)
        self.result_chi2_target.setStyleSheet(
            "background-color: transparent; color: #999999; border: none; font-size: 10px;")
        chi_vf_row.addWidget(self.result_chi2_target)
        chi_vf_row.addWidget(QLabel("Vf:"))
        self.result_vf = QLineEdit("—")
        self.result_vf.setReadOnly(True)
        self.result_vf.setStyleSheet(
            "background-color: #ecf0f1; color: #2c3e50; font-weight: bold;")
        chi_vf_row.addWidget(self.result_vf)
        res_layout.addLayout(chi_vf_row)

        # Row 2: Rg and Peak r  (expand to fill; no addStretch — uncertainties need room)
        rg_pr_row = QHBoxLayout()
        rg_pr_row.addWidget(QLabel("Rg:"))
        self.result_rg = QLineEdit("—")
        self.result_rg.setReadOnly(True)
        self.result_rg.setStyleSheet(
            "background-color: #ecf0f1; color: #2c3e50; font-weight: bold;")
        rg_pr_row.addWidget(self.result_rg)
        rg_pr_row.addWidget(QLabel("Å  Peak r:"))
        self.result_peak_r = QLineEdit("—")
        self.result_peak_r.setReadOnly(True)
        self.result_peak_r.setStyleSheet(
            "background-color: #ecf0f1; color: #2c3e50; font-weight: bold;")
        rg_pr_row.addWidget(self.result_peak_r)
        rg_pr_row.addWidget(QLabel("Å"))
        res_layout.addLayout(rg_pr_row)

        layout.addWidget(results_box)

        # ── Storage buttons ──────────────────────────────────────────────────
        store_row = QHBoxLayout()

        self.save_state_button = QPushButton("Save State")
        self.save_state_button.setMinimumHeight(26)
        self.save_state_button.setStyleSheet("""
            QPushButton { background-color: #3498db; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #2980b9; }
        """)
        self.save_state_button.clicked.connect(self.save_state)
        store_row.addWidget(self.save_state_button)

        self.store_file_button = QPushButton("Store in File")
        self.store_file_button.setMinimumHeight(26)
        self.store_file_button.setStyleSheet("background-color: lightgreen;")
        self.store_file_button.clicked.connect(self.store_results_to_file)
        store_row.addWidget(self.store_file_button)

        layout.addLayout(store_row)

        # ── Export / Import Parameters buttons ───────────────────────────────
        params_row = QHBoxLayout()

        self.export_params_button = QPushButton("Export Parameters")
        self.export_params_button.setMinimumHeight(26)
        self.export_params_button.setStyleSheet("background-color: lightgreen;")
        self.export_params_button.setToolTip(
            "Export current Sizes parameters to a pyIrena JSON config file."
        )
        self.export_params_button.clicked.connect(self.export_parameters)
        params_row.addWidget(self.export_params_button)

        self.import_params_button = QPushButton("Import Parameters")
        self.import_params_button.setMinimumHeight(26)
        self.import_params_button.setStyleSheet("background-color: lightgreen;")
        self.import_params_button.setToolTip(
            "Load Sizes parameters from a pyIrena JSON config file."
        )
        self.import_params_button.clicked.connect(self.import_parameters)
        params_row.addWidget(self.import_params_button)

        layout.addLayout(params_row)

        # ── Status label (left panel) ────────────────────────────────────────
        self.status_label = QLabel("Ready — load data to begin")
        self.status_label.setStyleSheet("""
            QLabel { color: #7f8c8d; padding: 5px; border-top: 1px solid #bdc3c7; font-size: 10px; }
        """)
        layout.addWidget(self.status_label)

        layout.addStretch()

        scroll.setWidget(inner)
        outer_layout.addWidget(scroll)

        # Initialise visibility
        self._on_method_changed(self.method_combo.currentText())
        self._on_shape_changed(self.shape_combo.currentText())

        return outer

    # ── Slot helpers ─────────────────────────────────────────────────────────

    def _on_shape_changed(self, shape):
        show_ar = (shape == "spheroid")
        self.aspect_label.setVisible(show_ar)
        self.aspect_ratio_edit.setVisible(show_ar)

    def _on_method_changed(self, method):
        self.maxent_group.setVisible(method == "maxent")
        self.reg_group.setVisible(method == "regularization")
        self.tnnls_group.setVisible(method == "tnnls")
        self.mcsas_group.setVisible(method == "mcsas")

    def _on_param_changed(self):
        """Called when any parameter field changes."""
        pass  # Could add auto-update here

    def _scale_edit(self, edit: ScrubbableLineEdit, factor: float):
        """Multiply the numeric value in *edit* by *factor* (for ÷10 / ×10 buttons)."""
        try:
            val = float(edit.text()) * factor
            edit.setText(f"{val:.3e}" if (abs(val) < 0.01 or abs(val) >= 10000) else f"{val:.4g}")
            edit.editingFinished.emit()
        except ValueError:
            pass

    def _update_error_display(self):
        """Redraw error bars on I(Q) plot immediately when Error scale changes.

        This lets the user see the effect of the scaling without running a fit.
        The data points and any existing fit line are left untouched.
        """
        if self.data is None or self.graph_window is None:
            return
        raw_err = self.data.get('Error')
        if raw_err is None:
            return
        try:
            error_scale = float(self.error_scale_edit.text() or 1.0)
        except ValueError:
            return
        disp_err = raw_err * error_scale if error_scale != 1.0 else raw_err
        self.graph_window.update_error_bars(
            self.data['Q'], self.data['Intensity'], disp_err
        )

    # ── Data loading ─────────────────────────────────────────────────────────

    def set_data(self, q, intensity, error=None, label='Data', filepath=None, is_nxcansas=False):
        """Set the SAS data to be fitted."""
        self.data = {
            'Q': q,
            'Intensity': intensity,
            'Error': error,
            'label': label,
            'filepath': filepath,
            'is_nxcansas': is_nxcansas,
        }
        self.fit_result = None
        self._last_distribution = None

        if self.graph_window:
            if filepath:
                self.graph_window.data_folder = str(Path(filepath).parent)
            self.graph_window.init_plots()
            self.graph_window.plot_data(q, intensity, error, label)
            # Restore saved cursor Q range (cursors exist after plot_data)
            if self._pending_cursor_q_range is not None:
                self.graph_window.set_cursor_range(*self._pending_cursor_q_range)
                self._pending_cursor_q_range = None

        self.status_label.setText(f"Loaded: {label} ({len(q)} points)")

    # ── Parameter collection ─────────────────────────────────────────────────

    def _collect_params(self) -> SizesDistribution:
        """Read all GUI fields and return a configured SizesDistribution."""
        s = SizesDistribution()
        s.r_min = float(self.rmin_edit.text() or 10.0)
        s.r_max = float(self.rmax_edit.text() or 1000.0)
        s.n_bins = self.nbins_spin.value()
        s.log_spacing = self.log_spacing_check.isChecked()
        s.shape = self.shape_combo.currentText()
        s.contrast = float(self.contrast_edit.text() or 1.0)
        if s.shape == "spheroid":
            s.shape_params = {'aspect_ratio': float(self.aspect_ratio_edit.text() or 1.5)}
        else:
            s.shape_params = {}
        s.background = float(self.background_edit.text() or 0.0)
        s.error_scale = float(self.error_scale_edit.text() or 1.0)
        s.power_law_B = float(self.power_law_B_edit.text() or 0.0)
        s.power_law_P = float(self.power_law_P_edit.text() or 4.0)
        s.method = self.method_combo.currentText()
        s.maxent_sky_background = float(self.maxent_sky_edit.text() or 1e-6)
        s.maxent_stability = 0.01   # hardcoded; see docs/sizes_methods.md
        s.maxent_max_iter = self.maxent_maxiter_spin.value()
        s.regularization_evalue = float(self.reg_evalue_edit.text() or 1.0)
        s.regularization_min_ratio = float(self.reg_minratio_edit.text() or 1e-4)
        s.tnnls_approach_param = float(self.tnnls_approach_edit.text() or 0.95)
        s.tnnls_max_iter = self.tnnls_maxiter_spin.value()
        s.mcsas_n_repetitions = 1  # main fit always uses a single MC run
        s.mcsas_convergence = float(self.mcsas_conv_edit.text() or 1.0)
        s.mcsas_max_iter = self.mcsas_maxiter_spin.value()
        return s

    def _get_bg_fit_q_range(self, q_min_edit, q_max_edit) -> tuple:
        """Return (q_min, q_max) for a background fit using the edit fields,
        falling back to the cursor range if the fields are empty."""
        raw_min = q_min_edit.text().strip()
        raw_max = q_max_edit.text().strip()
        if raw_min and raw_max:
            try:
                return float(raw_min), float(raw_max)
            except ValueError:
                pass
        # Fall back to cursor range
        cr = self.graph_window.get_cursor_range() if self.graph_window else None
        if cr is not None:
            return cr
        # Last resort: full Q range of loaded data
        if self.data is not None:
            q = self.data['Q']
            return float(q.min()), float(q.max())
        return 0.0, 1.0

    def _get_q_filtered_data(self):
        """Return Q-range filtered (q, intensity, error) using cursor positions."""
        q = self.data['Q']
        intensity = self.data['Intensity']
        error = self.data.get('Error')

        cursor_range = self.graph_window.get_cursor_range()
        if cursor_range is not None:
            q_min, q_max = cursor_range
            # Update Q range display
            self.qmin_display.setText(f"{q_min:.4e}")
            self.qmax_display.setText(f"{q_max:.4e}")
            mask = (q >= q_min) & (q <= q_max)
            q = q[mask]
            intensity = intensity[mask]
            error = error[mask] if error is not None else None

        return q, intensity, error

    # ── Compute model (no fitting) ───────────────────────────────────────────

    def compute_model(self):
        """Compute model I(Q) using current parameters (no fitting)."""
        if self.data is None:
            self.graph_window.show_error_message("No data loaded.")
            return

        try:
            s = self._collect_params()
            q, intensity, error = self._get_q_filtered_data()

            # Build G matrix and compute forward model with a flat distribution
            from pyirena.core.form_factors import build_g_matrix, make_r_grid
            r_grid = make_r_grid(s.r_min, s.r_max, s.n_bins, s.log_spacing)
            G = build_g_matrix(q, r_grid, s.shape, s.contrast, **s.shape_params)

            from pyirena.core.form_factors import bin_widths as _bin_widths
            dw = _bin_widths(r_grid)

            # P(r) = x_raw / bin_widths, so to recover x_raw: multiply back
            if self._last_distribution is not None:
                r_prev, p_prev = self._last_distribution
                if len(r_prev) == len(r_grid):
                    x_raw = p_prev * dw
                else:
                    x_raw = np.ones(len(r_grid)) / len(r_grid)
            else:
                x_raw = np.ones(len(r_grid)) / len(r_grid)

            complex_bg = s.compute_complex_background(q)
            I_model = G @ x_raw + complex_bg

            # Plot — scale display errors to match what was used for fitting
            raw_err = self.data.get('Error')
            disp_err = (raw_err * float(s.error_scale)
                        if raw_err is not None and s.error_scale != 1.0
                        else raw_err)
            self.graph_window.init_plots()
            self.graph_window.plot_data(
                self.data['Q'], self.data['Intensity'], disp_err, self.data['label']
            )

            # Show complex background and corrected data when background is non-trivial
            q_full = self.data['Q']
            bg_full = s.compute_complex_background(q_full)
            if s.power_law_B != 0.0 or s.background != 0.0:
                self.graph_window.plot_complex_background(q_full, bg_full, 'Complex bg')
                I_corrected = self.data['Intensity'] - bg_full
                cursor_range = self.graph_window.get_cursor_range()
                if cursor_range is not None:
                    cq_min, cq_max = cursor_range
                    fit_mask = ((q_full >= cq_min) & (q_full <= cq_max)
                                & (I_corrected > 0))
                else:
                    fit_mask = I_corrected > 0
                if np.any(fit_mask):
                    self.graph_window.plot_corrected_data(
                        q_full[fit_mask], I_corrected[fit_mask], 'I−bg'
                    )

            # Plot model on top
            self.graph_window.plot_fit(q, I_model, 'Model')

            if error is not None and np.any(error > 0):
                residuals = (intensity - I_model) / error
            else:
                residuals = (intensity - I_model) / np.maximum(intensity, 1e-40)
            self.graph_window.plot_residuals(q, residuals)

            if self._last_distribution is not None:
                r_prev, p_prev = self._last_distribution
                self.graph_window.plot_distribution(r_prev, p_prev)

            self.status_label.setText("Model computed (not fitted)")

        except Exception as exc:
            self.graph_window.show_error_message(f"Error computing model: {exc}")
            import traceback; traceback.print_exc()

    # ── Fit ──────────────────────────────────────────────────────────────────

    def run_fit(self):
        """Run the size distribution fit with current parameters."""
        if self.data is None:
            self.graph_window.show_error_message("No data loaded.")
            return

        # Backup current results for revert
        self._param_backup = self._collect_params()

        try:
            s = self._collect_params()
            q, intensity, error = self._get_q_filtered_data()

            if len(q) < 5:
                self.graph_window.show_error_message(
                    "Too few data points in Q range (need ≥ 5)."
                )
                return

            self.status_label.setText(f"Fitting {len(q)} points with {s.method}…")

            # McSAS: single-run for the main Fit button (no error bars).
            # Uncertainty comes from "Calculate Uncertainty" like other methods.
            if s.method == 'mcsas':
                s.mcsas_n_repetitions = 1

            result = s.fit(q, intensity, error)

            if not result.get('success', False):
                self.graph_window.show_error_message(
                    f"Fit did not converge: {result.get('message', '')}"
                )
                self.status_label.setText("Fit failed")
                return

            self.fit_result = result
            self.sizes = s  # Store fitted object

            r_grid = result['r_grid']
            distribution = result['distribution']
            self._last_distribution = (r_grid, distribution)

            # ── Update results display ────────────────────────────────────────
            chi2 = result.get('chi_squared', float('nan'))
            vf = result.get('volume_fraction', float('nan'))
            rg = result.get('rg', float('nan'))
            # chi² target = number of data points used in the fit (good fit → chi² ≈ target)
            n_data = result.get('n_data', len(q))

            peak_idx = int(np.argmax(distribution)) if distribution is not None else 0
            peak_r = float(r_grid[peak_idx]) if r_grid is not None else float('nan')

            self.result_chi2.setText(_fmt(chi2))
            self.result_chi2_target.setText(f"(target: {n_data})")
            self.result_vf.setText(_fmt(vf, sig=4))
            self.result_rg.setText(_fmt(rg))
            self.result_peak_r.setText(_fmt(peak_r))

            # ── Plots ─────────────────────────────────────────────────────────
            # model_intensity from fit() is G@x_raw (background already subtracted
            # from data before fitting).  Add the complex background back so the
            # model curve sits on top of the raw measured data in the graph.
            I_model_bg_subtracted = result['model_intensity']
            complex_bg_q = s.compute_complex_background(q)
            I_model_display = I_model_bg_subtracted + complex_bg_q
            residuals = result.get('residuals', None)

            # Re-plot on full Q range — scale display errors to match fitting
            raw_err = self.data.get('Error')
            disp_err = (raw_err * float(s.error_scale)
                        if raw_err is not None and s.error_scale != 1.0
                        else raw_err)
            self.graph_window.init_plots()
            self.graph_window.plot_data(
                self.data['Q'], self.data['Intensity'],
                disp_err, self.data['label']
            )

            # Show complex background model and corrected data
            # (corrected data limited to the cursor Q range used for the fit)
            q_full = self.data['Q']
            bg_full = s.compute_complex_background(q_full)
            if s.power_law_B != 0.0 or s.background != 0.0:
                self.graph_window.plot_complex_background(q_full, bg_full, 'Complex bg')
                I_corrected = self.data['Intensity'] - bg_full
                # Limit to cursor Q range (the range actually used in the fit)
                cursor_range = self.graph_window.get_cursor_range()
                if cursor_range is not None:
                    cq_min, cq_max = cursor_range
                    fit_mask = ((q_full >= cq_min) & (q_full <= cq_max)
                                & (I_corrected > 0))
                else:
                    fit_mask = I_corrected > 0
                if np.any(fit_mask):
                    self.graph_window.plot_corrected_data(
                        q_full[fit_mask], I_corrected[fit_mask], 'I−bg'
                    )

            # Plot fit on top so it is visible over other items
            self.graph_window.plot_fit(q, I_model_display, 'Fitted Model')

            if residuals is not None:
                self.graph_window.plot_residuals(q, residuals)
            self.graph_window.plot_distribution(r_grid, distribution)
            dist_std = result.get('distribution_std')
            if dist_std is not None:
                self.graph_window.plot_distribution_uncertainty(r_grid, distribution, dist_std)

            # ── Status message ────────────────────────────────────────────────
            n_iter = result.get('n_iterations', '?')
            msg = (
                f"Fit complete | method: {s.method} | "
                f"χ²: {chi2:.4g} | "
                f"Vf: {vf:.4g} | "
                f"Rg: {rg:.4g} Å | "
                f"peak r: {peak_r:.4g} Å | "
                f"iterations: {n_iter}"
            )
            self.graph_window.show_success_message(msg)
            self.status_label.setText(f"Fit done: χ²={chi2:.4g}, Vf={vf:.4g}")

        except Exception as exc:
            self.graph_window.show_error_message(f"Error during fitting: {exc}")
            self.status_label.setText("Fit failed")
            import traceback; traceback.print_exc()

    # ── Background fit actions ────────────────────────────────────────────────

    def _replot_background_preview(self):
        """Replot the complex background model and corrected data on the main graph
        without running a full size-distribution fit."""
        if self.data is None or self.graph_window is None:
            return
        s = self._collect_params()
        q_full = self.data['Q']
        bg_full = s.compute_complex_background(q_full)
        # Remove old background/corrected items and redraw
        if self.graph_window._complex_bg_item is not None:
            try:
                self.graph_window.main_plot.removeItem(self.graph_window._complex_bg_item)
            except Exception:
                pass
            self.graph_window._complex_bg_item = None
        if self.graph_window._corrected_item is not None:
            try:
                self.graph_window.main_plot.removeItem(self.graph_window._corrected_item)
            except Exception:
                pass
            self.graph_window._corrected_item = None
        if s.power_law_B != 0.0 or s.background != 0.0:
            self.graph_window.plot_complex_background(q_full, bg_full, 'Complex bg')
            I_corrected = self.data['Intensity'] - bg_full
            cursor_range = self.graph_window.get_cursor_range()
            if cursor_range is not None:
                cq_min, cq_max = cursor_range
                fit_mask = ((q_full >= cq_min) & (q_full <= cq_max)
                            & (I_corrected > 0))
            else:
                fit_mask = I_corrected > 0
            if np.any(fit_mask):
                self.graph_window.plot_corrected_data(
                    q_full[fit_mask], I_corrected[fit_mask], 'I−bg'
                )

    def fit_power_law_action(self):
        """Fit B·q^(-P) to data using the Q range in the Background tab."""
        if self.data is None:
            self.graph_window.show_error_message("No data loaded.")
            return
        fit_B = self.fit_B_check.isChecked()
        fit_P = self.fit_P_check.isChecked()
        if not fit_B and not fit_P:
            self.graph_window.show_error_message(
                "Check at least one of 'Fit B?' or 'Fit P?' to fit the power law."
            )
            return
        try:
            s = self._collect_params()
            q_min, q_max = self._get_bg_fit_q_range(self.pl_q_min_edit, self.pl_q_max_edit)
            result = s.fit_power_law(
                self.data['Q'], self.data['Intensity'],
                q_min, q_max, fit_B=fit_B, fit_P=fit_P
            )
            if result['success']:
                self.power_law_B_edit.setText(f"{result['B']:.6g}")
                self.power_law_P_edit.setText(f"{result['P']:.6g}")
                self.graph_window.show_success_message(result['message'])
                self.status_label.setText(f"P/B fit: B={result['B']:.4g}, P={result['P']:.4g}")
                self._replot_background_preview()
            else:
                self.graph_window.show_error_message(result['message'])
        except Exception as exc:
            self.graph_window.show_error_message(f"Error fitting P/B: {exc}")
            import traceback; traceback.print_exc()

    def fit_background_action(self):
        """Fit flat background by averaging I − B·q^(-P) in the Q range."""
        if self.data is None:
            self.graph_window.show_error_message("No data loaded.")
            return
        try:
            s = self._collect_params()
            q_min, q_max = self._get_bg_fit_q_range(self.bg_q_min_edit, self.bg_q_max_edit)
            result = s.fit_background_term(
                self.data['Q'], self.data['Intensity'], q_min, q_max
            )
            if result['success']:
                self.background_edit.setText(f"{result['background']:.6g}")
                self.graph_window.show_success_message(result['message'])
                self.status_label.setText(f"Background fit: {result['background']:.4g}")
                self._replot_background_preview()
            else:
                self.graph_window.show_error_message(result['message'])
        except Exception as exc:
            self.graph_window.show_error_message(f"Error fitting background: {exc}")
            import traceback; traceback.print_exc()

    def fit_all_action(self):
        """Fit P/B (if checked), then background (if non-empty range), then sizes."""
        if self.data is None:
            self.graph_window.show_error_message("No data loaded.")
            return
        fit_B = self.fit_B_check.isChecked()
        fit_P = self.fit_P_check.isChecked()
        # Step 1: Fit power law if at least one parameter is selected
        if fit_B or fit_P:
            self.fit_power_law_action()
            # Check for failure
            if "Error" in self.status_label.text() or "failed" in self.status_label.text().lower():
                return
        # Step 2: Fit flat background if a Q range is specified
        bg_qmin = self.bg_q_min_edit.text().strip()
        bg_qmax = self.bg_q_max_edit.text().strip()
        if bg_qmin and bg_qmax:
            self.fit_background_action()
            if "Error" in self.status_label.text() or "failed" in self.status_label.text().lower():
                return
        # Step 3: Fit size distribution
        self.run_fit()

    # ── Cursor range helpers (Background tab) ────────────────────────────────

    def _set_pl_q_from_cursors(self):
        """Copy current cursor Q positions into the power-law Q range fields."""
        if self.graph_window is None:
            return
        cr = self.graph_window.get_cursor_range()
        if cr is None:
            self.graph_window.show_error_message("Load data first to initialise cursors.")
            return
        q_min, q_max = cr
        self.pl_q_min_edit.setText(f"{q_min:.6g}")
        self.pl_q_max_edit.setText(f"{q_max:.6g}")

    def _set_bg_q_from_cursors(self):
        """Copy current cursor Q positions into the flat-background Q range fields."""
        if self.graph_window is None:
            return
        cr = self.graph_window.get_cursor_range()
        if cr is None:
            self.graph_window.show_error_message("Load data first to initialise cursors.")
            return
        q_min, q_max = cr
        self.bg_q_min_edit.setText(f"{q_min:.6g}")
        self.bg_q_max_edit.setText(f"{q_max:.6g}")

    # ── Monte-Carlo uncertainty estimation ───────────────────────────────────

    def calculate_uncertainty(self):
        """Run N noise-perturbed fits to estimate size distribution uncertainties.

        Each data point is perturbed by  ΔI = σᵢ · N(0,1)  where σᵢ is the
        measurement error.  The spread of the resulting size distributions gives
        per-bin statistical uncertainties; Rg, Vf, and peak r are also propagated.

        The number of perturbed fits is controlled by the Passes spinner (default 10).
        For McSAS each perturbed fit uses a single internal MC run.
        """
        if self.data is None:
            self.graph_window.show_error_message("No data loaded.")
            return
        if self.fit_result is None:
            self.graph_window.show_error_message(
                "Run a fit first before calculating uncertainty."
            )
            return

        s = self._collect_params()

        # Number of perturbed fits — controlled by the Passes spinner.
        # For McSAS, ensure each perturbed fit uses a single internal MC run.
        N_runs = self.unc_n_runs_spin.value()
        if s.method == 'mcsas':
            s.mcsas_n_repetitions = 1  # one MC run per perturbed fit

        q, intensity, error = self._get_q_filtered_data()

        if len(q) < 5:
            self.graph_window.show_error_message(
                "Too few data points in Q range (need ≥ 5)."
            )
            return

        # If no experimental errors, generate synthetic ones for the MC runs
        if error is None:
            error = intensity * 0.05
            error = np.maximum(error, 1e-30)
        err_mc = error * float(s.error_scale) if s.error_scale != 1.0 else error
        err_mc = np.maximum(err_mc, 1e-30)

        self.status_label.setText(f"Running {N_runs} MC fits for uncertainty…")

        distributions, rg_vals, vf_vals, peak_r_vals = [], [], [], []
        r_grid_ref = None

        for i in range(N_runs):
            # Perturb data with Gaussian noise proportional to the scaled errors
            I_perturbed = intensity + err_mc * np.random.randn(len(intensity))
            try:
                result = s.fit(q, I_perturbed, error)
            except Exception:
                continue
            if not result.get('success', False):
                continue

            dist = result.get('distribution')
            rg   = result.get('rg', np.nan)
            vf   = result.get('volume_fraction', np.nan)
            r_g  = result.get('r_grid')

            if dist is None or r_g is None:
                continue

            if r_grid_ref is None:
                r_grid_ref = r_g

            distributions.append(dist)
            rg_vals.append(rg)
            vf_vals.append(vf)
            peak_idx = int(np.argmax(dist))
            peak_r_vals.append(float(r_g[peak_idx]))

        n_ok = len(distributions)
        if n_ok == 0:
            self.graph_window.show_error_message(
                f"All {N_runs} MC fits failed — cannot estimate uncertainty."
            )
            self.status_label.setText("MC uncertainty: all fits failed")
            return

        # Compute per-bin statistics
        dist_array = np.array(distributions)          # (n_ok, n_bins)
        dist_mean  = np.mean(dist_array, axis=0)
        dist_std   = np.std(dist_array, axis=0, ddof=min(1, n_ok - 1))

        # Scalar statistics
        rg_mean  = float(np.nanmean(rg_vals))
        rg_std   = float(np.nanstd(rg_vals,  ddof=min(1, n_ok - 1)))
        vf_mean  = float(np.nanmean(vf_vals))
        vf_std   = float(np.nanstd(vf_vals,  ddof=min(1, n_ok - 1)))
        pr_mean  = float(np.nanmean(peak_r_vals))
        pr_std   = float(np.nanstd(peak_r_vals, ddof=min(1, n_ok - 1)))

        # Update result fields with ± notation (2 sig figs on uncertainty)
        self.result_rg.setText(_fmt_unc(rg_mean, rg_std))
        self.result_vf.setText(_fmt_unc(vf_mean, vf_std))
        self.result_peak_r.setText(_fmt_unc(pr_mean, pr_std))

        # Re-plot the mean distribution with uncertainty error bars
        self.graph_window.plot_distribution(r_grid_ref, dist_mean)
        self.graph_window.plot_distribution_uncertainty(r_grid_ref, dist_mean, dist_std)

        msg = (
            f"MC uncertainty ({n_ok}/{N_runs} fits OK) | "
            f"Rg: {rg_mean:.4g} ± {rg_std:.4g} Å | "
            f"Vf: {vf_mean:.4g} ± {vf_std:.4g} | "
            f"peak r: {pr_mean:.4g} ± {pr_std:.4g} Å"
        )
        self.graph_window.show_success_message(msg)
        self.status_label.setText(
            f"MC done ({n_ok}/{N_runs}): Rg={rg_mean:.4g}±{rg_std:.4g} Å"
        )

    # ── Revert / Reset ───────────────────────────────────────────────────────

    def _revert(self):
        if self._param_backup is None:
            self.graph_window.show_error_message("No previous state to revert to.")
            return
        # Revert by re-applying the backed-up SizesDistribution params
        s = self._param_backup
        self._apply_sizes_to_gui(s)
        self.graph_window.show_message("Reverted to pre-fit parameters.")

    def reset_to_defaults(self):
        s = SizesDistribution()  # Default values
        self._apply_sizes_to_gui(s)
        self.graph_window.show_message("Reset to default parameters.")

    def _apply_sizes_to_gui(self, s: SizesDistribution):
        self.shape_combo.setCurrentText(s.shape)
        ar = s.shape_params.get('aspect_ratio', 1.5)
        self.aspect_ratio_edit.setText(str(ar))
        self.contrast_edit.setText(str(s.contrast))
        self.rmin_edit.setText(str(s.r_min))
        self.rmax_edit.setText(str(s.r_max))
        self.nbins_spin.setValue(s.n_bins)
        self.log_spacing_check.setChecked(s.log_spacing)
        self.background_edit.setText(str(s.background))
        self.power_law_B_edit.setText(str(s.power_law_B))
        self.power_law_P_edit.setText(str(s.power_law_P))
        self.method_combo.setCurrentText(s.method)
        self.maxent_sky_edit.setText(str(s.maxent_sky_background))
        self.maxent_maxiter_spin.setValue(s.maxent_max_iter)
        self.reg_evalue_edit.setText(str(s.regularization_evalue))
        self.reg_minratio_edit.setText(str(s.regularization_min_ratio))
        self.tnnls_approach_edit.setText(str(s.tnnls_approach_param))
        self.tnnls_maxiter_spin.setValue(s.tnnls_max_iter)
        self.error_scale_edit.setText(str(s.error_scale))

    # ── State save / load ────────────────────────────────────────────────────

    def _get_current_state(self) -> dict:
        """Serialize current GUI state to a flat dict."""
        s = self._collect_params()
        state = {
            'r_min': s.r_min,
            'r_max': s.r_max,
            'n_bins': s.n_bins,
            'log_spacing': s.log_spacing,
            'shape': s.shape,
            'contrast': s.contrast,
            'aspect_ratio': s.shape_params.get('aspect_ratio', 1.5),
            'background': s.background,
            'method': s.method,
            'maxent_sky_background': s.maxent_sky_background,
            'maxent_stability': s.maxent_stability,
            'maxent_max_iter': s.maxent_max_iter,
            'regularization_evalue': s.regularization_evalue,
            'regularization_min_ratio': s.regularization_min_ratio,
            'tnnls_approach_param': s.tnnls_approach_param,
            'tnnls_max_iter': s.tnnls_max_iter,
            'mcsas_convergence': s.mcsas_convergence,
            'unc_n_runs': self.unc_n_runs_spin.value(),
            'mcsas_max_iter': s.mcsas_max_iter,
            'error_scale': s.error_scale,
            'power_law_B': s.power_law_B,
            'power_law_P': s.power_law_P,
            'power_law_q_min': self.pl_q_min_edit.text().strip() or None,
            'power_law_q_max': self.pl_q_max_edit.text().strip() or None,
            'background_q_min': self.bg_q_min_edit.text().strip() or None,
            'background_q_max': self.bg_q_max_edit.text().strip() or None,
        }
        # Save cursor Q range so it survives across sessions
        if self.graph_window is not None:
            cr = self.graph_window.get_cursor_range()
            if cr is not None:
                state['cursor_q_min'] = cr[0]
                state['cursor_q_max'] = cr[1]
        return state

    def _apply_state(self, state: dict):
        # Cursor Q range: store for application once data (and thus cursors) are loaded
        q_min = state.get('cursor_q_min')
        q_max = state.get('cursor_q_max')
        if q_min is not None and q_max is not None:
            try:
                self._pending_cursor_q_range = (float(q_min), float(q_max))
            except (TypeError, ValueError):
                self._pending_cursor_q_range = None

        self.rmin_edit.setText(str(state.get('r_min', 10.0)))
        self.rmax_edit.setText(str(state.get('r_max', 1000.0)))
        self.nbins_spin.setValue(int(state.get('n_bins', 200)))
        self.log_spacing_check.setChecked(bool(state.get('log_spacing', True)))
        self.shape_combo.setCurrentText(state.get('shape', 'sphere'))
        self.contrast_edit.setText(str(state.get('contrast', 1.0)))
        self.aspect_ratio_edit.setText(str(state.get('aspect_ratio', 1.5)))
        self.background_edit.setText(str(state.get('background', 0.0)))
        self.method_combo.setCurrentText(state.get('method', 'regularization'))
        self.maxent_sky_edit.setText(str(state.get('maxent_sky_background', 1e-6)))
        self.maxent_maxiter_spin.setValue(int(state.get('maxent_max_iter', 1000)))
        self.reg_evalue_edit.setText(str(state.get('regularization_evalue', 1.0)))
        self.reg_minratio_edit.setText(str(state.get('regularization_min_ratio', 1e-4)))
        self.tnnls_approach_edit.setText(str(state.get('tnnls_approach_param', 0.95)))
        self.tnnls_maxiter_spin.setValue(int(state.get('tnnls_max_iter', 1000)))
        self.mcsas_conv_edit.setText(str(state.get('mcsas_convergence', 1.0)))
        self.mcsas_maxiter_spin.setValue(int(state.get('mcsas_max_iter', 100000)))
        self.unc_n_runs_spin.setValue(int(state.get('unc_n_runs', 10)))
        self.error_scale_edit.setText(str(state.get('error_scale', 1.0)))
        self.power_law_B_edit.setText(str(state.get('power_law_B', 0.0)))
        self.power_law_P_edit.setText(str(state.get('power_law_P', 4.0)))
        pl_q_min = state.get('power_law_q_min')
        pl_q_max = state.get('power_law_q_max')
        self.pl_q_min_edit.setText(str(pl_q_min) if pl_q_min is not None else '')
        self.pl_q_max_edit.setText(str(pl_q_max) if pl_q_max is not None else '')
        bg_q_min = state.get('background_q_min')
        bg_q_max = state.get('background_q_max')
        self.bg_q_min_edit.setText(str(bg_q_min) if bg_q_min is not None else '')
        self.bg_q_max_edit.setText(str(bg_q_max) if bg_q_max is not None else '')

    def load_state(self):
        state = self.state_manager.get('sizes')
        if state:
            try:
                self._apply_state(state)
            except Exception as exc:
                print(f"Warning: could not restore sizes state: {exc}")

    def save_state(self):
        state = self._get_current_state()
        self.state_manager.update('sizes', state)
        if self.state_manager.save():
            QMessageBox.information(self, "State Saved", "State saved successfully.")
            self.status_label.setText("State saved")
        else:
            QMessageBox.warning(self, "Save Failed", "Failed to save state.")

    def _get_data_folder(self) -> str:
        """Return the folder of the currently loaded data file, or cwd."""
        if self.data is not None:
            fp = self.data.get('filepath')
            if fp:
                return str(Path(fp).parent)
        return str(Path.cwd())

    def export_parameters(self):
        """Export current Sizes parameters to a pyIrena JSON config file."""
        import json
        import datetime
        from pyirena import __version__ as _version

        default_dir = self._get_data_folder()
        default_path = str(Path(default_dir) / "pyirena_config.json")

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Sizes Parameters",
            default_path,
            "pyIrena Config (*.json);;All Files (*)"
        )
        if not file_path:
            return

        file_path = Path(file_path)

        # Load existing config or start fresh
        config = {}
        if file_path.exists():
            try:
                with open(file_path, 'r') as f:
                    config = json.load(f)
            except Exception:
                config = {}

            if '_pyirena_config' not in config:
                QMessageBox.warning(
                    self,
                    "Not a pyIrena File",
                    f"The selected file is not a pyIrena configuration file:\n{file_path}\n\n"
                    "Choose a different file or enter a new filename."
                )
                return

            if 'sizes' in config:
                reply = QMessageBox.question(
                    self,
                    "Overwrite Sizes Parameters?",
                    f"File already contains Sizes parameters:\n{file_path}\n\n"
                    "Overwrite the existing Sizes group?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                )
                if reply != QMessageBox.StandardButton.Yes:
                    return

        now = datetime.datetime.now().isoformat(timespec='seconds')
        if '_pyirena_config' not in config:
            config['_pyirena_config'] = {
                'file_type': 'pyIrena Configuration File',
                'version': _version,
                'created': now,
            }
        config['_pyirena_config']['modified'] = now
        config['_pyirena_config']['written_by'] = f"pyIrena {_version}"

        state = self._get_current_state()
        self.state_manager.update('sizes', state)
        config['sizes'] = self.state_manager.get('sizes')

        try:
            with open(file_path, 'w') as f:
                json.dump(config, f, indent=2)
        except Exception as e:
            QMessageBox.warning(self, "Export Failed", f"Could not write file:\n{e}")
            return

        msg = f"Sizes parameters saved to: {file_path.name}"
        self.status_label.setText(msg)
        if self.graph_window:
            self.graph_window.show_success_message(msg)

    def import_parameters(self):
        """Import Sizes parameters from a pyIrena JSON config file."""
        import json

        default_dir = self._get_data_folder()

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Import Sizes Parameters",
            default_dir,
            "pyIrena Config (*.json);;All Files (*)"
        )
        if not file_path:
            return

        file_path = Path(file_path)

        try:
            with open(file_path, 'r') as f:
                config = json.load(f)
        except Exception as e:
            QMessageBox.warning(self, "Import Failed", f"Could not read file:\n{e}")
            return

        if '_pyirena_config' not in config:
            QMessageBox.warning(
                self,
                "Not a pyIrena File",
                f"The selected file is not a pyIrena configuration file:\n{file_path}"
            )
            return

        if 'sizes' not in config:
            QMessageBox.warning(
                self,
                "No Sizes Parameters",
                f"The file does not contain a Sizes parameter group:\n{file_path}"
            )
            return

        sizes_state = config['sizes']
        self.state_manager.update('sizes', sizes_state)
        self._apply_state(sizes_state)

        written_by = config['_pyirena_config'].get('written_by', 'unknown version')
        msg = f"Sizes parameters loaded from: {file_path.name}  (written by {written_by})"
        self.status_label.setText(f"Loaded config: {file_path.name}")
        if self.graph_window:
            self.graph_window.show_success_message(msg)

    # ── Store results to file ────────────────────────────────────────────────

    def store_results_to_file(self):
        """Save size distribution results to NXcanSAS HDF5 file."""
        from pyirena.io.nxcansas_sizes import save_sizes_results

        if self.data is None:
            self.graph_window.show_error_message("No data loaded.")
            return
        if self.fit_result is None:
            self.graph_window.show_error_message(
                "No fit results yet — run a fit first."
            )
            return

        try:
            source_path = self.data.get('filepath')
            is_nxcansas = self.data.get('is_nxcansas', False)

            if source_path is None:
                self.graph_window.show_error_message(
                    "Cannot determine output file — load data from a file first."
                )
                return

            from pyirena.io.nxcansas_unified import get_output_filepath
            output_path = get_output_filepath(Path(source_path), is_nxcansas)

            q = self.data['Q']
            i_err_raw = self.data.get('Error')

            # Apply cursor Q range if set
            cursor_range = self.graph_window.get_cursor_range()
            if cursor_range is not None:
                q_min, q_max = cursor_range
                mask = (q >= q_min) & (q <= q_max)
                q_fit = q[mask]
                i_data = self.data['Intensity'][mask]
                i_err = i_err_raw[mask] if i_err_raw is not None else None
            else:
                q_fit = q
                i_data = self.data['Intensity']
                i_err = i_err_raw

            result = self.fit_result
            params = self._get_current_state()
            # Fit results
            params['chi_squared'] = result.get('chi_squared')
            params['volume_fraction'] = result.get('volume_fraction')
            params['rg'] = result.get('rg')
            params['n_iterations'] = result.get('n_iterations')

            save_sizes_results(
                filepath=output_path,
                q=q_fit,
                intensity_data=i_data,
                intensity_model=result['model_intensity'],
                residuals=result.get('residuals', np.zeros_like(q_fit)),
                r_grid=result['r_grid'],
                distribution=result['distribution'],
                params=params,
                intensity_error=i_err,
                distribution_std=result.get('distribution_std'),
            )

            self.graph_window.show_success_message(
                f"Results stored to: {output_path.name}"
            )
            self.status_label.setText(f"Stored to {output_path.name}")

        except Exception as exc:
            self.graph_window.show_error_message(f"Error storing results: {exc}")
            import traceback; traceback.print_exc()


# ─────────────────────────────────────────────────────────────────────────────
# Utility
# ─────────────────────────────────────────────────────────────────────────────

def _fmt(value, sig=4):
    """Format a scalar for display (3-4 sig figs)."""
    try:
        v = float(value)
        if not np.isfinite(v):
            return "—"
        if v == 0:
            return "0"
        if abs(v) < 0.01 or abs(v) >= 10000:
            return f"{v:.{sig-1}e}"
        return f"{v:.{sig}g}"
    except Exception:
        return "—"


def _fmt_unc(value, uncertainty):
    """Format 'value ± uncertainty' with 2 significant figures on the uncertainty.

    The value is rounded to the same decimal place as the last significant digit
    of the uncertainty.  Falls back to _fmt for each if uncertainty is not finite
    or positive.

    Examples
    --------
    _fmt_unc(123.456, 1.23)   → "123 ± 1.2"
    _fmt_unc(0.01234, 0.00123) → "0.0123 ± 0.0012"
    _fmt_unc(1.23e-5, 1.1e-6) → "1.23e-05 ± 1.1e-06"
    """
    import math
    try:
        v = float(value)
        u = float(uncertainty)
        if not (np.isfinite(v) and np.isfinite(u) and u > 0):
            return f"{_fmt(v)} ± {_fmt(u)}"
        mag_u = math.floor(math.log10(abs(u)))       # exponent of leading digit of u
        dp = max(0, 1 - mag_u)                        # decimal places for 2 sig figs
        # Use plain fixed-point when values are in a human-friendly range
        use_exp = abs(v) >= 1e5 or (abs(v) < 1e-3 and v != 0)
        if use_exp:
            # Scale both to same power of 10 as the value for consistency
            exp = math.floor(math.log10(abs(v))) if v != 0 else 0
            scale = 10 ** exp
            v_s = v / scale
            u_s = u / scale
            dp_s = max(0, 1 - (mag_u - exp))
            return f"{v_s:.{dp_s}f}e{exp:+03d} ± {u_s:.{dp_s}f}e{exp:+03d}"
        return f"{v:.{dp}f} ± {u:.{dp}f}"
    except Exception:
        return f"{_fmt(value)} ± {_fmt(uncertainty)}"
