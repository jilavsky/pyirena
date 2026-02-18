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
        QMessageBox, QScrollArea, QGroupBox, QSizePolicy, QFrame,
    )
    from PySide6.QtCore import Qt, Signal
    from PySide6.QtGui import QDoubleValidator, QIntValidator
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
            QLabel, QLineEdit, QComboBox, QCheckBox, QSpinBox, QSplitter,
            QMessageBox, QScrollArea, QGroupBox, QSizePolicy, QFrame,
        )
        from PyQt6.QtCore import Qt, Signal
        from PyQt6.QtGui import QDoubleValidator, QIntValidator
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
    """A QLineEdit that responds to mouse wheel to increment/decrement values."""

    def __init__(self, text="", parent=None):
        super().__init__(text, parent)
        self.setFocusPolicy(Qt.FocusPolicy.WheelFocus)

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
                step = magnitude * 0.1
                value += direction * step
            # Format nicely
            if abs(value) < 0.01 or abs(value) >= 10000:
                self.setText(f"{value:.3e}")
            else:
                self.setText(f"{value:.4g}")
            self.editingFinished.emit()
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
        self._fit_item = None
        self._resid_item = None
        self._dist_item = None
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
        self.main_plot.addLegend(offset=(10, 10))
        _style_axes(self.main_plot)

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

        # ── Distribution plot  (log-x, linear-y — r in Å) ───────────────────
        self.distribution_plot.setLogMode(x=True, y=False)
        self.distribution_plot.setLabel('left',   'P(r)  (a.u.)')
        self.distribution_plot.setLabel('bottom', 'r  (Å)')
        self.distribution_plot.showGrid(x=True, y=True, alpha=0.3)
        self.distribution_plot.addLegend(offset=(10, 10))
        _style_axes(self.distribution_plot)

        layout.addWidget(self.graphics_layout)

        # ── Status message ───────────────────────────────────────────────────
        self.status_message = QLabel("")
        self.status_message.setMaximumHeight(60)
        self.status_message.setWordWrap(True)
        self.status_message.setStyleSheet("""
            QLabel {
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
        self._fit_item = None
        self._resid_item = None
        self._dist_item = None

    def plot_data(self, q, intensity, error=None, label='Data'):
        """Plot I(Q) data on the main plot.

        Pass linear Q and intensity values — pyqtgraph applies log10 internally
        when the plot is in log mode.  Robust Y-range is set automatically to
        avoid a single bad low-intensity point compressing the view.
        """
        # pyqtgraph handles the log10 transform; just pass raw values
        self._data_item = self.main_plot.plot(
            q, intensity,
            pen=None, symbol='o', symbolSize=4,
            symbolPen=pg.mkPen('#2c3e50', width=1),
            symbolBrush=pg.mkBrush('#2c3e50'),
            name=label,
        )
        self._ensure_cursors(q)
        self._set_robust_y_range(intensity)

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

    def plot_fit(self, q, intensity_model, label='Fitted Model'):
        """Plot model I(Q) on main plot."""
        if self._fit_item is not None:
            try:
                self.main_plot.removeItem(self._fit_item)
            except Exception:
                pass
        self._fit_item = self.main_plot.plot(
            q, intensity_model,
            pen=pg.mkPen('#e74c3c', width=2),
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
        """Plot size distribution on distribution plot."""
        if self._dist_item is not None:
            try:
                self.distribution_plot.removeItem(self._dist_item)
            except Exception:
                pass
        self._dist_item = self.distribution_plot.plot(
            r, distribution,
            pen=pg.mkPen('#2980b9', width=2),
            symbol='o', symbolSize=4,
            symbolPen=pg.mkPen('#2980b9', width=1),
            symbolBrush=pg.mkBrush('#2980b9'),
            name=label
        )

    # ── Status helpers ───────────────────────────────────────────────────────

    def show_message(self, text, color='#2c3e50'):
        self.status_message.setText(text)
        self.status_message.setStyleSheet(f"""
            QLabel {{
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
        self._last_distribution = None   # Tuple (r, P) from last fit/compute
        self._param_backup = None        # For revert functionality

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
        """Build the left control panel inside a scroll area."""
        # Outer container with fixed width
        outer = QWidget()
        outer.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Preferred)
        outer.setMinimumWidth(400)
        outer.setMaximumWidth(400)

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

        # ── Shape group ──────────────────────────────────────────────────────
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

        layout.addWidget(shape_box)

        # ── Q range display (from cursors) ───────────────────────────────────
        q_box = QGroupBox("Q Range (cursors)")
        q_layout = QVBoxLayout(q_box)
        q_layout.setSpacing(3)

        qmin_row = QHBoxLayout()
        qmin_row.addWidget(QLabel("Q min:"))
        self.qmin_display = QLineEdit("—")
        self.qmin_display.setReadOnly(True)
        self.qmin_display.setMaximumWidth(100)
        self.qmin_display.setStyleSheet("background-color: #ecf0f1; color: #7f8c8d;")
        qmin_row.addWidget(self.qmin_display)
        qmin_row.addWidget(QLabel("Å⁻¹"))
        qmin_row.addStretch()
        q_layout.addLayout(qmin_row)

        qmax_row = QHBoxLayout()
        qmax_row.addWidget(QLabel("Q max:"))
        self.qmax_display = QLineEdit("—")
        self.qmax_display.setReadOnly(True)
        self.qmax_display.setMaximumWidth(100)
        self.qmax_display.setStyleSheet("background-color: #ecf0f1; color: #7f8c8d;")
        qmax_row.addWidget(self.qmax_display)
        qmax_row.addWidget(QLabel("Å⁻¹"))
        qmax_row.addStretch()
        q_layout.addLayout(qmax_row)

        hint = QLabel("Drag cursors on I(Q) graph to set Q range")
        hint.setStyleSheet("color: #7f8c8d; font-size: 10px;")
        q_layout.addWidget(hint)

        layout.addWidget(q_box)

        # ── Size grid ────────────────────────────────────────────────────────
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
            row.addStretch()
            grid_layout.addLayout(row)

        bins_row = QHBoxLayout()
        bins_row.addWidget(QLabel("Number of bins:"))
        self.nbins_spin = QSpinBox()
        self.nbins_spin.setMinimum(5)
        self.nbins_spin.setMaximum(300)
        self.nbins_spin.setValue(50)
        self.nbins_spin.setMaximumWidth(70)
        self.nbins_spin.valueChanged.connect(self._on_param_changed)
        bins_row.addWidget(self.nbins_spin)
        bins_row.addStretch()
        grid_layout.addLayout(bins_row)

        log_row = QHBoxLayout()
        self.log_spacing_check = QCheckBox("Logarithmic spacing")
        self.log_spacing_check.stateChanged.connect(self._on_param_changed)
        log_row.addWidget(self.log_spacing_check)
        log_row.addStretch()
        grid_layout.addLayout(log_row)

        layout.addWidget(grid_box)

        # ── Background ───────────────────────────────────────────────────────
        bg_box = QGroupBox("Background")
        bg_layout = QHBoxLayout(bg_box)
        bg_layout.addWidget(QLabel("Background:"))
        self.background_edit = ScrubbableLineEdit("0.0")
        self.background_edit.setValidator(QDoubleValidator(-1e10, 1e10, 6))
        self.background_edit.setMaximumWidth(100)
        self.background_edit.editingFinished.connect(self._on_param_changed)
        bg_layout.addWidget(self.background_edit)
        bg_layout.addWidget(QLabel("cm⁻¹"))
        bg_layout.addStretch()
        layout.addWidget(bg_box)

        # ── Method & parameters ──────────────────────────────────────────────
        method_box = QGroupBox("Fitting Method")
        method_layout = QVBoxLayout(method_box)
        method_layout.setSpacing(5)

        meth_row = QHBoxLayout()
        meth_row.addWidget(QLabel("Method:"))
        self.method_combo = QComboBox()
        self.method_combo.addItems(["regularization", "maxent", "tnnls"])
        self.method_combo.currentTextChanged.connect(self._on_method_changed)
        meth_row.addWidget(self.method_combo)
        meth_row.addStretch()
        method_layout.addLayout(meth_row)

        # ── MaxEnt params ────────────────────────────────────────────────────
        self.maxent_group = QWidget()
        maxent_layout = QVBoxLayout(self.maxent_group)
        maxent_layout.setContentsMargins(0, 0, 0, 0)
        maxent_layout.setSpacing(3)

        for label_text, attr, default in [
            ("Sky background:", "maxent_sky_edit", "1e-06"),
            ("Stability:", "maxent_stab_edit", "0.01"),
        ]:
            r = QHBoxLayout()
            r.addWidget(QLabel(label_text))
            e = ScrubbableLineEdit(default)
            e.setValidator(QDoubleValidator(0.0, 1e10, 6))
            e.setMaximumWidth(90)
            e.editingFinished.connect(self._on_param_changed)
            setattr(self, attr, e)
            r.addWidget(e)
            r.addStretch()
            maxent_layout.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Max iterations:"))
        self.maxent_maxiter_spin = QSpinBox()
        self.maxent_maxiter_spin.setMinimum(10)
        self.maxent_maxiter_spin.setMaximum(10000)
        self.maxent_maxiter_spin.setValue(1000)
        self.maxent_maxiter_spin.setMaximumWidth(80)
        self.maxent_maxiter_spin.valueChanged.connect(self._on_param_changed)
        r.addWidget(self.maxent_maxiter_spin)
        r.addStretch()
        maxent_layout.addLayout(r)
        method_layout.addWidget(self.maxent_group)

        # ── Regularization params ────────────────────────────────────────────
        self.reg_group = QWidget()
        reg_layout = QVBoxLayout(self.reg_group)
        reg_layout.setContentsMargins(0, 0, 0, 0)
        reg_layout.setSpacing(3)

        for label_text, attr, default in [
            ("χ² tolerance (evalue):", "reg_evalue_edit", "1.0"),
            ("Min ratio:", "reg_minratio_edit", "1e-04"),
        ]:
            r = QHBoxLayout()
            r.addWidget(QLabel(label_text))
            e = ScrubbableLineEdit(default)
            e.setValidator(QDoubleValidator(0.0, 1e10, 6))
            e.setMaximumWidth(90)
            e.editingFinished.connect(self._on_param_changed)
            setattr(self, attr, e)
            r.addWidget(e)
            r.addStretch()
            reg_layout.addLayout(r)

        method_layout.addWidget(self.reg_group)

        # ── TNNLS params ─────────────────────────────────────────────────────
        self.tnnls_group = QWidget()
        tnnls_layout = QVBoxLayout(self.tnnls_group)
        tnnls_layout.setContentsMargins(0, 0, 0, 0)
        tnnls_layout.setSpacing(3)

        r = QHBoxLayout()
        r.addWidget(QLabel("Approach parameter:"))
        self.tnnls_approach_edit = ScrubbableLineEdit("0.95")
        self.tnnls_approach_edit.setValidator(QDoubleValidator(0.01, 1.0, 6))
        self.tnnls_approach_edit.setMaximumWidth(90)
        self.tnnls_approach_edit.editingFinished.connect(self._on_param_changed)
        r.addWidget(self.tnnls_approach_edit)
        r.addStretch()
        tnnls_layout.addLayout(r)

        r2 = QHBoxLayout()
        r2.addWidget(QLabel("Max iterations:"))
        self.tnnls_maxiter_spin = QSpinBox()
        self.tnnls_maxiter_spin.setMinimum(10)
        self.tnnls_maxiter_spin.setMaximum(100000)
        self.tnnls_maxiter_spin.setValue(1000)
        self.tnnls_maxiter_spin.setMaximumWidth(80)
        self.tnnls_maxiter_spin.valueChanged.connect(self._on_param_changed)
        r2.addWidget(self.tnnls_maxiter_spin)
        r2.addStretch()
        tnnls_layout.addLayout(r2)

        method_layout.addWidget(self.tnnls_group)

        layout.addWidget(method_box)

        # ── Action buttons ───────────────────────────────────────────────────
        btn_row1 = QHBoxLayout()

        self.graph_button = QPushButton("Graph Model")
        self.graph_button.setMinimumHeight(28)
        self.graph_button.setStyleSheet("""
            QPushButton { background-color: #52c77a; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #3eb56a; }
        """)
        self.graph_button.clicked.connect(self.compute_model)
        btn_row1.addWidget(self.graph_button)
        btn_row1.addStretch()

        self.fit_button = QPushButton("Fit")
        self.fit_button.setMinimumHeight(28)
        self.fit_button.setStyleSheet("""
            QPushButton { background-color: #27ae60; color: white; font-weight: bold; font-size: 13px; }
            QPushButton:hover { background-color: #1e8449; }
        """)
        self.fit_button.clicked.connect(self.run_fit)
        btn_row1.addWidget(self.fit_button)

        self.revert_button = QPushButton("Revert")
        self.revert_button.setMinimumHeight(28)
        self.revert_button.setStyleSheet("""
            QPushButton { background-color: #e67e22; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #f39c12; }
        """)
        self.revert_button.clicked.connect(self._revert)
        btn_row1.addWidget(self.revert_button)

        layout.addLayout(btn_row1)

        self.reset_button = QPushButton("Reset to Defaults")
        self.reset_button.setMinimumHeight(26)
        self.reset_button.setStyleSheet("""
            QPushButton { background-color: #e67e22; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #d35400; }
        """)
        self.reset_button.clicked.connect(self.reset_to_defaults)
        layout.addWidget(self.reset_button)

        # ── Results section ──────────────────────────────────────────────────
        results_header = QLabel("Results")
        results_header.setStyleSheet("""
            QLabel { font-weight: bold; color: #3498db; font-size: 12px; margin-top: 4px; }
        """)
        layout.addWidget(results_header)

        results_box = QGroupBox("")
        res_layout = QVBoxLayout(results_box)
        res_layout.setSpacing(3)

        for label_text, attr, units in [
            ("χ²:", "result_chi2", ""),
            ("Volume fraction:", "result_vf", ""),
            ("Rg:", "result_rg", "Å"),
            ("Peak r:", "result_peak_r", "Å"),
        ]:
            r = QHBoxLayout()
            r.addWidget(QLabel(label_text))
            field = QLineEdit("—")
            field.setReadOnly(True)
            field.setMaximumWidth(100)
            field.setStyleSheet("background-color: #ecf0f1; color: #2c3e50; font-weight: bold;")
            setattr(self, attr, field)
            r.addWidget(field)
            if units:
                r.addWidget(QLabel(units))
            r.addStretch()
            res_layout.addLayout(r)

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

        # ── Status label (left panel) ────────────────────────────────────────
        self.status_label = QLabel("Ready — load data to begin")
        self.status_label.setStyleSheet("""
            QLabel { color: #7f8c8d; padding: 5px; border-top: 1px solid #bdc3c7; font-size: 10px; }
        """)
        layout.addWidget(self.status_label)

        layout.addStretch()

        scroll.setWidget(inner)
        outer_layout.addWidget(scroll)

        # Hide method groups; show the active one
        self._on_method_changed(self.method_combo.currentText())
        # Hide aspect ratio by default (sphere selected)
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

    def _on_param_changed(self):
        """Called when any parameter field changes."""
        pass  # Could add auto-update here

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
        s.method = self.method_combo.currentText()
        s.maxent_sky_background = float(self.maxent_sky_edit.text() or 1e-6)
        s.maxent_stability = float(self.maxent_stab_edit.text() or 0.01)
        s.maxent_max_iter = self.maxent_maxiter_spin.value()
        s.regularization_evalue = float(self.reg_evalue_edit.text() or 1.0)
        s.regularization_min_ratio = float(self.reg_minratio_edit.text() or 1e-4)
        s.tnnls_approach_param = float(self.tnnls_approach_edit.text() or 0.95)
        s.tnnls_max_iter = self.tnnls_maxiter_spin.value()
        return s

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

            I_model = G @ x_raw + s.background

            # Plot
            self.graph_window.init_plots()
            self.graph_window.plot_data(
                self.data['Q'], self.data['Intensity'], self.data.get('Error'), self.data['label']
            )
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

            peak_idx = int(np.argmax(distribution)) if distribution is not None else 0
            peak_r = float(r_grid[peak_idx]) if r_grid is not None else float('nan')

            self.result_chi2.setText(_fmt(chi2))
            self.result_vf.setText(_fmt(vf, sig=4))
            self.result_rg.setText(_fmt(rg))
            self.result_peak_r.setText(_fmt(peak_r))

            # ── Plots ─────────────────────────────────────────────────────────
            I_model = result['model_intensity']
            residuals = result.get('residuals', None)

            # Re-plot on full Q range but with fitted model only on fitted range
            self.graph_window.init_plots()
            self.graph_window.plot_data(
                self.data['Q'], self.data['Intensity'],
                self.data.get('Error'), self.data['label']
            )
            self.graph_window.plot_fit(q, I_model, 'Fitted Model')
            if residuals is not None:
                self.graph_window.plot_residuals(q, residuals)
            self.graph_window.plot_distribution(r_grid, distribution)

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
        self.method_combo.setCurrentText(s.method)
        self.maxent_sky_edit.setText(str(s.maxent_sky_background))
        self.maxent_stab_edit.setText(str(s.maxent_stability))
        self.maxent_maxiter_spin.setValue(s.maxent_max_iter)
        self.reg_evalue_edit.setText(str(s.regularization_evalue))
        self.reg_minratio_edit.setText(str(s.regularization_min_ratio))
        self.tnnls_approach_edit.setText(str(s.tnnls_approach_param))
        self.tnnls_maxiter_spin.setValue(s.tnnls_max_iter)
        self.qpower_edit.setText(str(s.q_power))

    # ── State save / load ────────────────────────────────────────────────────

    def _get_current_state(self) -> dict:
        """Serialize current GUI state to a flat dict."""
        s = self._collect_params()
        return {
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
        }

    def _apply_state(self, state: dict):
        self.rmin_edit.setText(str(state.get('r_min', 10.0)))
        self.rmax_edit.setText(str(state.get('r_max', 1000.0)))
        self.nbins_spin.setValue(int(state.get('n_bins', 50)))
        self.log_spacing_check.setChecked(bool(state.get('log_spacing', False)))
        self.shape_combo.setCurrentText(state.get('shape', 'sphere'))
        self.contrast_edit.setText(str(state.get('contrast', 1.0)))
        self.aspect_ratio_edit.setText(str(state.get('aspect_ratio', 1.5)))
        self.background_edit.setText(str(state.get('background', 0.0)))
        self.method_combo.setCurrentText(state.get('method', 'regularization'))
        self.maxent_sky_edit.setText(str(state.get('maxent_sky_background', 1e-6)))
        self.maxent_stab_edit.setText(str(state.get('maxent_stability', 0.01)))
        self.maxent_maxiter_spin.setValue(int(state.get('maxent_max_iter', 1000)))
        self.reg_evalue_edit.setText(str(state.get('regularization_evalue', 1.0)))
        self.reg_minratio_edit.setText(str(state.get('regularization_min_ratio', 1e-4)))
        self.tnnls_approach_edit.setText(str(state.get('tnnls_approach_param', 0.95)))
        self.tnnls_maxiter_spin.setValue(int(state.get('tnnls_max_iter', 1000)))

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
            # Get the Q range used for fitting
            cursor_range = self.graph_window.get_cursor_range()
            if cursor_range is not None:
                q_min, q_max = cursor_range
                mask = (q >= q_min) & (q <= q_max)
                q_fit = q[mask]
                i_data = self.data['Intensity'][mask]
            else:
                q_fit = q
                i_data = self.data['Intensity']

            result = self.fit_result
            params = self._get_current_state()
            params['chi_squared'] = result.get('chi_squared')
            params['volume_fraction'] = result.get('volume_fraction')
            params['rg'] = result.get('rg')

            save_sizes_results(
                filepath=output_path,
                q=q_fit,
                intensity_data=i_data,
                intensity_model=result['model_intensity'],
                residuals=result.get('residuals', np.zeros_like(q_fit)),
                r_grid=result['r_grid'],
                distribution=result['distribution'],
                params=params,
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
