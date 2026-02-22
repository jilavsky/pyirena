"""
Simple Fits GUI panel for pyIrena.

Provides ``SimpleFitsGraphWindow`` (three-panel pyqtgraph display: I(Q),
residuals, linearization) and ``SimpleFitsPanel`` (controls + graph) for
interactive single-model fitting of SAS data.
"""

try:
    from PySide6.QtWidgets import (
        QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
        QPushButton, QLabel, QLineEdit, QComboBox, QCheckBox, QSpinBox,
        QSplitter, QMessageBox, QScrollArea, QGroupBox, QSizePolicy, QFrame,
    )
    from PySide6.QtCore import Qt, Signal
    from PySide6.QtGui import QDoubleValidator
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
            QPushButton, QLabel, QLineEdit, QComboBox, QCheckBox, QSpinBox,
            QSplitter, QMessageBox, QScrollArea, QGroupBox, QSizePolicy, QFrame,
        )
        from PyQt6.QtCore import Qt, Signal
        from PyQt6.QtGui import QDoubleValidator
    except ImportError:
        raise ImportError("Neither PySide6 nor PyQt6 found.  Install with: pip install PySide6")

import numpy as np
from pathlib import Path

import pyqtgraph as pg

from pyirena.core.simple_fits import SimpleFitModel, MODEL_NAMES, MODEL_REGISTRY
from pyirena.state.state_manager import StateManager

# Re-use the ScrubbableLineEdit and LogDecadeAxis from sizes_panel to keep
# a consistent look and feel.
from pyirena.gui.sizes_panel import ScrubbableLineEdit, LogDecadeAxis


# ===========================================================================
# SimpleFitsGraphWindow
# ===========================================================================

class SimpleFitsGraphWindow(QWidget):
    """
    Three-panel pyqtgraph display:
      - Row 0 (50%): log-log I(Q) vs Q with data scatter + model line + cursors
      - Row 1 (10%): residuals (I−model)/err vs Q  (log-x, linear-y)
      - Row 2 (40%): linearization (Guinier/Porod) or disabled placeholder
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._data_item  = None
        self._fit_item   = None
        self._resid_item = None
        self._lin_data_item = None
        self._lin_fit_item  = None
        self._cursor_a   = None
        self._cursor_b   = None
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(2)
        self.setLayout(layout)

        self.graphics_layout = pg.GraphicsLayoutWidget()
        self.graphics_layout.setBackground('w')

        # ── I(Q) main plot ────────────────────────────────────────────────────
        self.main_plot = self.graphics_layout.addPlot(
            row=0, col=0,
            axisItems={
                'left':   LogDecadeAxis(orientation='left'),
                'bottom': LogDecadeAxis(orientation='bottom'),
            },
        )
        self.main_plot.setLabel('left',   'I')
        self.main_plot.setLabel('bottom', 'Q  (Å⁻¹)')
        self.main_plot.showGrid(x=True, y=True, alpha=0.3)
        self.main_plot.getAxis('left').enableAutoSIPrefix(False)
        self.main_plot.getAxis('bottom').enableAutoSIPrefix(False)

        # ── Residuals plot ────────────────────────────────────────────────────
        self.residuals_plot = self.graphics_layout.addPlot(
            row=1, col=0,
            axisItems={'bottom': LogDecadeAxis(orientation='bottom')},
        )
        self.residuals_plot.setLabel('left',   '(I−fit)/err')
        self.residuals_plot.setLabel('bottom', 'Q  (Å⁻¹)')
        self.residuals_plot.showGrid(x=True, y=True, alpha=0.3)
        self.residuals_plot.setXLink(self.main_plot)
        self.residuals_plot.getAxis('left').enableAutoSIPrefix(False)
        self.residuals_plot.getAxis('bottom').enableAutoSIPrefix(False)
        # Zero line
        zero = pg.InfiniteLine(pos=0, angle=0, pen=pg.mkPen('k', width=1, style=Qt.PenStyle.DashLine))
        self.residuals_plot.addItem(zero)

        # ── Linearization plot ─────────────────────────────────────────────────
        self.lin_plot = self.graphics_layout.addPlot(
            row=2, col=0,
        )
        self.lin_plot.setLabel('left',   'Y')
        self.lin_plot.setLabel('bottom', 'X')
        self.lin_plot.showGrid(x=True, y=True, alpha=0.3)
        self.lin_plot.getAxis('left').enableAutoSIPrefix(False)
        self.lin_plot.getAxis('bottom').enableAutoSIPrefix(False)
        self._lin_title = self.lin_plot.setTitle('Linearization', color='#555555', size='10pt')

        # ── Height ratios 5:1:4 ───────────────────────────────────────────────
        ci = self.graphics_layout.ci
        ci.layout.setRowStretchFactor(0, 5)
        ci.layout.setRowStretchFactor(1, 1)
        ci.layout.setRowStretchFactor(2, 4)

        layout.addWidget(self.graphics_layout)

        # ── Status label ──────────────────────────────────────────────────────
        self.status_label = QLabel('')
        self.status_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.status_label.setStyleSheet('font-size: 11px; color: #444;')
        layout.addWidget(self.status_label)

    # ── Data plotting ─────────────────────────────────────────────────────────

    def plot_data(self, q: np.ndarray, I: np.ndarray, dI=None, label: str = 'Data'):
        """Plot the raw I(Q) data as scatter + optional error bars."""
        self.main_plot.clear()
        self._fit_item = None

        # Restore cursors removed by clear() — always re-add if they exist
        if self._cursor_a is not None:
            self.main_plot.addItem(self._cursor_a)
        if self._cursor_b is not None:
            self.main_plot.addItem(self._cursor_b)

        mask = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
        q, I = q[mask], I[mask]

        if dI is not None:
            dI = np.asarray(dI, dtype=float)[mask]
            dI = np.maximum(dI, 1e-30 * I)

        scatter = pg.ScatterPlotItem(
            x=np.log10(q), y=np.log10(I),
            pen=None, brush=pg.mkBrush(80, 80, 80, 180), size=5,
        )
        self.main_plot.addItem(scatter)
        self._data_item = scatter
        # Only create cursors on the first data load; afterwards keep user positions
        self._ensure_cursors(np.log10(q.min()), np.log10(q.max()))

    def plot_fit(self, q_fit: np.ndarray, I_model: np.ndarray):
        """Overlay the model fit line on the main plot."""
        if self._fit_item is not None:
            self.main_plot.removeItem(self._fit_item)

        mask = np.isfinite(q_fit) & np.isfinite(I_model) & (q_fit > 0) & (I_model > 0)
        if mask.sum() < 2:
            return
        pen = pg.mkPen(color=(200, 30, 30), width=2)
        fit_line = pg.PlotDataItem(
            x=np.log10(q_fit[mask]), y=np.log10(I_model[mask]),
            pen=pen,
        )
        self.main_plot.addItem(fit_line)
        self._fit_item = fit_line

    def plot_residuals(self, q: np.ndarray, residuals: np.ndarray):
        """Plot (I−model)/err residuals."""
        self.residuals_plot.clear()
        zero = pg.InfiniteLine(pos=0, angle=0,
                               pen=pg.mkPen('k', width=1,
                                            style=Qt.PenStyle.DashLine))
        self.residuals_plot.addItem(zero)
        self._resid_item = None

        mask = np.isfinite(q) & np.isfinite(residuals) & (q > 0)
        if mask.sum() < 1:
            return
        scatter = pg.ScatterPlotItem(
            x=np.log10(q[mask]), y=residuals[mask],
            pen=None, brush=pg.mkBrush(80, 80, 80, 180), size=4,
        )
        self.residuals_plot.addItem(scatter)
        self._resid_item = scatter

    def plot_linearization(self, lin_result):
        """Plot linearized data + fit.  lin_result is from SimpleFitModel.linearize()."""
        self.lin_plot.clear()
        self._lin_data_item = None
        self._lin_fit_item  = None

        if lin_result is None:
            self.lin_plot.setTitle('Linearization — not available for this model',
                                   color='#888888', size='9pt')
            return

        X = lin_result['X']
        Y = lin_result['Y']
        X_fit = lin_result['X_fit']
        Y_fit = lin_result['Y_fit']

        self.lin_plot.setTitle(
            f"{lin_result['y_label']} vs {lin_result['x_label']}   "
            f"slope={lin_result['slope']:.4g}  intercept={lin_result['intercept']:.4g}  "
            f"R²={lin_result['r_squared']:.4f}",
            color='#333333', size='9pt',
        )
        self.lin_plot.setLabel('left',   lin_result['y_label'])
        self.lin_plot.setLabel('bottom', lin_result['x_label'])

        mask = np.isfinite(X) & np.isfinite(Y)
        data_item = pg.ScatterPlotItem(
            x=X[mask], y=Y[mask],
            pen=None, brush=pg.mkBrush(80, 80, 80, 180), size=5,
        )
        self.lin_plot.addItem(data_item)
        self._lin_data_item = data_item

        mask2 = np.isfinite(X_fit) & np.isfinite(Y_fit)
        fit_item = pg.PlotDataItem(
            x=X_fit[mask2], y=Y_fit[mask2],
            pen=pg.mkPen(color=(200, 30, 30), width=2),
        )
        self.lin_plot.addItem(fit_item)
        self._lin_fit_item = fit_item

    # ── Cursors ───────────────────────────────────────────────────────────────

    def _ensure_cursors(self, log_min: float, log_max: float):
        """Create cursors on the first call and place them within the data range.

        On subsequent calls the user's dragged positions are preserved; the
        cursors have already been re-added to the plot by plot_data().
        """
        if self._cursor_a is not None:
            return   # Already exist; positions kept as user left them

        span = log_max - log_min
        pos_a = log_min + 0.1 * span
        pos_b = log_max - 0.1 * span

        self._cursor_a = pg.InfiniteLine(
            pos=pos_a, angle=90, movable=True,
            pen=pg.mkPen(color=(200, 50, 50), width=2),
            label='A', labelOpts={'position': 0.05, 'color': (200, 50, 50)},
        )
        self.main_plot.addItem(self._cursor_a)

        self._cursor_b = pg.InfiniteLine(
            pos=pos_b, angle=90, movable=True,
            pen=pg.mkPen(color=(50, 100, 200), width=2),
            label='B', labelOpts={'position': 0.10, 'color': (50, 100, 200)},
        )
        self.main_plot.addItem(self._cursor_b)

    def get_cursor_range(self):
        """Return (q_min, q_max) in linear space from cursor positions."""
        if self._cursor_a is None or self._cursor_b is None:
            return None, None
        a = 10.0 ** self._cursor_a.getPos()[0]
        b = 10.0 ** self._cursor_b.getPos()[0]
        return (min(a, b), max(a, b))

    def set_cursor_range(self, q_min: float, q_max: float):
        """Position cursors at the given Q values (linear space)."""
        if q_min is not None and q_max is not None and q_min > 0 and q_max > 0:
            self._ensure_cursors(np.log10(q_min), np.log10(q_max))
            self._cursor_a.setPos(np.log10(q_min))
            self._cursor_b.setPos(np.log10(q_max))


# ===========================================================================
# SimpleFitsPanel
# ===========================================================================

class SimpleFitsPanel(QWidget):
    """
    Simple Fits interactive panel.

    Combines a control panel (model selection, parameters, action buttons)
    with a ``SimpleFitsGraphWindow`` in a horizontal QSplitter.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Simple Fits')
        self.resize(1200, 700)

        self.state_manager = StateManager()
        self.model = SimpleFitModel()
        self.data: dict | None = None
        self.fit_result: dict | None = None

        # Dynamic widget references (rebuilt when model changes)
        self._param_value_edits:   dict[str, QLineEdit] = {}
        self._param_lo_edits:      dict[str, QLineEdit] = {}
        self._param_hi_edits:      dict[str, QLineEdit] = {}
        self._param_unc_labels:    dict[str, QLabel]    = {}
        self._param_grid_widget:   QWidget | None = None

        self.init_ui()
        self.load_state()

    # ── UI construction ───────────────────────────────────────────────────────

    def init_ui(self):
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(4, 4, 4, 4)
        self.setLayout(main_layout)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self._create_control_panel())
        self.graph_window = SimpleFitsGraphWindow()
        splitter.addWidget(self.graph_window)
        splitter.setSizes([380, 820])
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        main_layout.addWidget(splitter)

        # Status label at bottom
        self.status_label = QLabel('No data loaded.')
        self.status_label.setStyleSheet('font-size: 11px; color: #555;')
        main_layout.addWidget(self.status_label)

    def _create_control_panel(self) -> QWidget:
        panel = QWidget()
        panel.setMinimumWidth(350)
        panel.setMaximumWidth(420)
        layout = QVBoxLayout()
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)
        panel.setLayout(layout)

        # ── Model selector ────────────────────────────────────────────────────
        model_row = QHBoxLayout()
        model_row.addWidget(QLabel('Model:'))
        self.model_combo = QComboBox()
        self.model_combo.addItems(MODEL_NAMES)
        self.model_combo.setCurrentText(self.model.model)
        self.model_combo.currentTextChanged.connect(self._on_model_changed)
        model_row.addWidget(self.model_combo, 1)
        layout.addLayout(model_row)

        # ── Q range ───────────────────────────────────────────────────────────
        q_box = QGroupBox('Q range for fit')
        q_layout = QGridLayout()
        q_layout.setContentsMargins(6, 4, 6, 4)
        q_box.setLayout(q_layout)

        q_layout.addWidget(QLabel('Q min:'), 0, 0)
        self.q_min_edit = QLineEdit('')
        self.q_min_edit.setPlaceholderText('all data')
        self.q_min_edit.setValidator(QDoubleValidator(0.0, 1e10, 6))
        self.q_min_edit.setMaximumWidth(80)
        q_layout.addWidget(self.q_min_edit, 0, 1)
        q_layout.addWidget(QLabel('Å⁻¹'), 0, 2)

        q_layout.addWidget(QLabel('Q max:'), 1, 0)
        self.q_max_edit = QLineEdit('')
        self.q_max_edit.setPlaceholderText('all data')
        self.q_max_edit.setValidator(QDoubleValidator(0.0, 1e10, 6))
        self.q_max_edit.setMaximumWidth(80)
        q_layout.addWidget(self.q_max_edit, 1, 1)
        q_layout.addWidget(QLabel('Å⁻¹'), 1, 2)

        self.cursors_btn = QPushButton('Set from cursors')
        self.cursors_btn.setFixedHeight(22)
        self.cursors_btn.clicked.connect(self._set_q_from_cursors)
        q_layout.addWidget(self.cursors_btn, 2, 0, 1, 3)
        layout.addWidget(q_box)

        # ── Background ────────────────────────────────────────────────────────
        self.complex_bg_check = QCheckBox('Complex background  (A·Q⁻ⁿ + flat)')
        self.complex_bg_check.setChecked(False)
        self.complex_bg_check.stateChanged.connect(self._on_complex_bg_changed)
        layout.addWidget(self.complex_bg_check)

        # ── Parameters ────────────────────────────────────────────────────────
        self.params_box = QGroupBox('Parameters')
        self.params_box_layout = QVBoxLayout()
        self.params_box_layout.setContentsMargins(4, 4, 4, 4)
        self.params_box.setLayout(self.params_box_layout)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setMinimumHeight(200)
        self.params_box_layout.addWidget(scroll)

        self._params_scroll_area = scroll
        layout.addWidget(self.params_box)

        # Build initial parameter widgets
        self._build_param_widgets()

        # ── Fit results ───────────────────────────────────────────────────────
        results_box = QGroupBox('Fit results')
        results_layout = QGridLayout()
        results_layout.setContentsMargins(6, 4, 6, 4)
        results_box.setLayout(results_layout)

        results_layout.addWidget(QLabel('χ²:'), 0, 0)
        self.chi2_label = QLabel('—')
        results_layout.addWidget(self.chi2_label, 0, 1)
        results_layout.addWidget(QLabel('Reduced χ²:'), 1, 0)
        self.rchi2_label = QLabel('—')
        results_layout.addWidget(self.rchi2_label, 1, 1)
        layout.addWidget(results_box)

        # ── Uncertainty options ───────────────────────────────────────────────
        unc_row = QHBoxLayout()
        unc_row.addWidget(QLabel('MC runs:'))
        self.n_runs_spin = QSpinBox()
        self.n_runs_spin.setRange(5, 500)
        self.n_runs_spin.setValue(50)
        self.n_runs_spin.setMaximumWidth(70)
        unc_row.addWidget(self.n_runs_spin)
        unc_row.addStretch()
        layout.addLayout(unc_row)

        # ── Buttons ───────────────────────────────────────────────────────────
        self.fit_btn = QPushButton('Fit')
        self.fit_btn.setMinimumHeight(28)
        self.fit_btn.setStyleSheet('font-weight: bold;')
        self.fit_btn.clicked.connect(self._run_fit)
        layout.addWidget(self.fit_btn)

        self.unc_btn = QPushButton('Calculate Uncertainty (MC)')
        self.unc_btn.setMinimumHeight(26)
        self.unc_btn.clicked.connect(self._calculate_uncertainty)
        layout.addWidget(self.unc_btn)

        self.store_btn = QPushButton('Store in File')
        self.store_btn.setMinimumHeight(26)
        self.store_btn.clicked.connect(self._store_results)
        layout.addWidget(self.store_btn)

        layout.addStretch()
        return panel

    # ── Parameter widget management ───────────────────────────────────────────

    def _build_param_widgets(self):
        """Destroy and rebuild the parameter grid for the current model."""
        if self._param_grid_widget is not None:
            self._params_scroll_area.takeWidget()
            self._param_grid_widget.deleteLater()

        self._param_value_edits.clear()
        self._param_lo_edits.clear()
        self._param_hi_edits.clear()
        self._param_unc_labels.clear()

        container = QWidget()
        grid = QGridLayout()
        grid.setContentsMargins(2, 2, 2, 2)
        grid.setSpacing(3)
        container.setLayout(grid)

        # Header row
        for col, text in enumerate(['Parameter', 'Value', 'lo', 'hi', '± std']):
            lbl = QLabel(text)
            lbl.setStyleSheet('font-weight: bold; font-size: 10px;')
            lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            grid.addWidget(lbl, 0, col)

        entry = MODEL_REGISTRY[self.model.model]
        param_specs = list(entry['params'])
        if self.model.use_complex_bg and entry['complex_bg']:
            from pyirena.core.simple_fits import _BG_PARAMS
            param_specs += _BG_PARAMS

        for row_i, (name, default, lo_def, hi_def) in enumerate(param_specs):
            row = row_i + 1
            val  = self.model.params.get(name, default)
            lo_v, hi_v = self.model.limits.get(name, (lo_def, hi_def))

            # Separator between model params and BG params
            if (self.model.use_complex_bg and entry['complex_bg']
                    and row_i == len(entry['params'])):
                sep = QFrame()
                sep.setFrameShape(QFrame.Shape.HLine)
                sep.setFrameShadow(QFrame.Shadow.Sunken)
                grid.addWidget(sep, row, 0, 1, 5)
                row_i += 1
                row += 1
                row = row_i + 1  # keep synced

            # Name label
            name_lbl = QLabel(name + ':')
            name_lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            name_lbl.setStyleSheet('font-size: 11px;')
            grid.addWidget(name_lbl, row, 0)

            # Value edit
            val_edit = ScrubbableLineEdit(f'{val:.6g}')
            val_edit.setValidator(QDoubleValidator(-1e30, 1e30, 10))
            val_edit.setMaximumWidth(85)
            val_edit.setMinimumWidth(70)
            val_edit.editingFinished.connect(self._on_param_edited)
            grid.addWidget(val_edit, row, 1)
            self._param_value_edits[name] = val_edit

            # Lo edit
            lo_edit = QLineEdit('' if lo_v is None else f'{lo_v:.4g}')
            lo_edit.setPlaceholderText('−∞')
            lo_edit.setValidator(QDoubleValidator(-1e30, 1e30, 8))
            lo_edit.setMaximumWidth(70)
            lo_edit.setMinimumWidth(55)
            grid.addWidget(lo_edit, row, 2)
            self._param_lo_edits[name] = lo_edit

            # Hi edit
            hi_edit = QLineEdit('' if hi_v is None else f'{hi_v:.4g}')
            hi_edit.setPlaceholderText('+∞')
            hi_edit.setValidator(QDoubleValidator(-1e30, 1e30, 8))
            hi_edit.setMaximumWidth(70)
            hi_edit.setMinimumWidth(55)
            grid.addWidget(hi_edit, row, 3)
            self._param_hi_edits[name] = hi_edit

            # Uncertainty label
            unc_lbl = QLabel('')
            unc_lbl.setStyleSheet('font-size: 10px; color: #666;')
            unc_lbl.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            grid.addWidget(unc_lbl, row, 4)
            self._param_unc_labels[name] = unc_lbl

        grid.setColumnStretch(4, 1)
        self._param_grid_widget = container
        self._params_scroll_area.setWidget(container)

    # ── Slots ─────────────────────────────────────────────────────────────────

    def _on_model_changed(self, model_name: str):
        """Switch model and rebuild parameter widgets."""
        self.save_state()
        self.model.set_model(model_name)
        self.model.use_complex_bg = self.complex_bg_check.isChecked()
        if self.model.use_complex_bg and MODEL_REGISTRY[model_name]['complex_bg']:
            from pyirena.core.simple_fits import _BG_PARAMS
            for name, default, lo, hi in _BG_PARAMS:
                self.model.params.setdefault(name, default)
                self.model.limits.setdefault(name, (lo, hi))
        self._build_param_widgets()
        self.fit_result = None
        self.chi2_label.setText('—')
        self.rchi2_label.setText('—')

    def _on_complex_bg_changed(self, state):
        """Toggle complex background and rebuild param widgets."""
        enabled = bool(state)
        self.model.use_complex_bg = enabled
        if enabled and MODEL_REGISTRY[self.model.model]['complex_bg']:
            from pyirena.core.simple_fits import _BG_PARAMS
            for name, default, lo, hi in _BG_PARAMS:
                self.model.params.setdefault(name, default)
                self.model.limits.setdefault(name, (lo, hi))
        self._build_param_widgets()

    def _on_param_edited(self):
        """Read the edited value back into self.model.params immediately."""
        for name, edit in self._param_value_edits.items():
            txt = edit.text().strip()
            if txt:
                try:
                    self.model.params[name] = float(txt)
                except ValueError:
                    pass

    # ── Data loading ──────────────────────────────────────────────────────────

    def set_data(self, q, intensity, error=None, label='Data',
                 filepath=None, is_nxcansas=False):
        """Load SAS data into the panel and plot it."""
        self.data = {
            'Q': np.asarray(q, dtype=float),
            'Intensity': np.asarray(intensity, dtype=float),
            'Error': np.asarray(error, dtype=float) if error is not None else None,
            'label': label,
            'filepath': filepath,
            'is_nxcansas': is_nxcansas,
        }
        self.fit_result = None
        self.chi2_label.setText('—')
        self.rchi2_label.setText('—')

        self.graph_window.plot_data(
            self.data['Q'],
            self.data['Intensity'],
            self.data['Error'],
            label=label,
        )
        self.setWindowTitle(f'Simple Fits — {label}')
        self.status_label.setText(f"Loaded: {label} ({len(q)} points)")

    # ── Q range helpers ───────────────────────────────────────────────────────

    def _set_q_from_cursors(self):
        """Copy cursor Q positions into the Q min/max fields."""
        q_min, q_max = self.graph_window.get_cursor_range()
        if q_min is not None:
            self.q_min_edit.setText(f'{q_min:.6g}')
        if q_max is not None:
            self.q_max_edit.setText(f'{q_max:.6g}')

    def _get_q_range(self):
        """Return (q_min, q_max) or (None, None) from the Q range fields."""
        def _parse(edit):
            txt = edit.text().strip()
            try:
                return float(txt) if txt else None
            except ValueError:
                return None
        return _parse(self.q_min_edit), _parse(self.q_max_edit)

    def _get_filtered_data(self):
        """Return (q, I, dI) arrays filtered to the active Q range."""
        if self.data is None:
            return None, None, None
        q  = self.data['Q']
        I  = self.data['Intensity']
        dI = self.data['Error']
        q_min, q_max = self._get_q_range()
        mask = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
        if q_min is not None:
            mask &= (q >= q_min)
        if q_max is not None:
            mask &= (q <= q_max)
        if dI is not None:
            dI_f = dI[mask]
        else:
            dI_f = None
        return q[mask], I[mask], dI_f

    # ── Collect model from widgets ─────────────────────────────────────────────

    def _collect_model(self) -> SimpleFitModel:
        """Read widget values into self.model and return it."""
        for name, edit in self._param_value_edits.items():
            txt = edit.text().strip()
            if txt:
                try:
                    self.model.params[name] = float(txt)
                except ValueError:
                    pass
        for name, edit in self._param_lo_edits.items():
            txt = edit.text().strip()
            lo = float(txt) if txt else None
            _, hi = self.model.limits.get(name, (None, None))
            self.model.limits[name] = (lo, hi)
        for name, edit in self._param_hi_edits.items():
            txt = edit.text().strip()
            hi = float(txt) if txt else None
            lo, _ = self.model.limits.get(name, (None, None))
            self.model.limits[name] = (lo, hi)
        return self.model

    def _apply_result_to_widgets(self, result: dict):
        """Write fitted parameter values and ± uncertainties back to widgets."""
        params = result.get('params', {})
        stds   = result.get('params_std', {})
        for name, edit in self._param_value_edits.items():
            if name in params:
                edit.setText(f'{params[name]:.6g}')
        for name, lbl in self._param_unc_labels.items():
            std = stds.get(name)
            if std is not None and np.isfinite(std):
                lbl.setText(f'±{std:.3g}')
            else:
                lbl.setText('')

    # ── Fit ───────────────────────────────────────────────────────────────────

    def _run_fit(self):
        if self.data is None:
            QMessageBox.warning(self, 'No data', 'Load data first.')
            return

        model = self._collect_model()
        q, I, dI = self._get_filtered_data()
        if len(q) < 2:
            QMessageBox.warning(self, 'Too few points',
                                'Not enough data points in the selected Q range.')
            return

        result = model.fit(q, I, dI)
        if not result['success']:
            QMessageBox.critical(self, 'Fit failed',
                                 f"Fit failed:\n{result.get('error', 'Unknown error')}")
            self.status_label.setText('Fit FAILED.')
            return

        self.fit_result = result

        # Update chi² display
        chi2 = result.get('chi2', float('nan'))
        rchi2 = result.get('reduced_chi2', float('nan'))
        self.chi2_label.setText(f'{chi2:.4g}')
        self.rchi2_label.setText(f'{rchi2:.4g}')

        # Write fitted values back to widgets
        self._apply_result_to_widgets(result)

        # Update derived quantities if any
        derived = result.get('derived', {})
        if derived:
            derived_txt = '   '.join(f'{k}={v:.4g}' for k, v in derived.items())
            self.status_label.setText(f'Fit OK | χ²_red={rchi2:.3g} | {derived_txt}')
        else:
            self.status_label.setText(f'Fit OK | Reduced χ² = {rchi2:.4g}')

        # Update all plots
        self._update_plots(result)
        self.save_state()

    # ── Uncertainty ───────────────────────────────────────────────────────────

    def _calculate_uncertainty(self):
        if self.data is None:
            QMessageBox.warning(self, 'No data', 'Load data first.')
            return
        if self.fit_result is None:
            QMessageBox.warning(self, 'No fit', 'Run Fit first.')
            return

        model = self._collect_model()
        q, I, dI = self._get_filtered_data()
        n_runs = self.n_runs_spin.value()
        dI_safe = dI if dI is not None else np.maximum(I * 0.05, 1e-30)

        mc_params: dict[str, list] = {k: [] for k in self.fit_result['params']}
        n_ok = 0
        for _ in range(n_runs):
            I_pert = I + dI_safe * np.random.randn(len(I))
            mc_res = model.fit(q, I_pert, dI)
            if mc_res.get('success'):
                for k in mc_params:
                    mc_params[k].append(mc_res['params'].get(k, float('nan')))
                n_ok += 1

        if n_ok < 2:
            QMessageBox.warning(self, 'Uncertainty',
                                f'Too few successful MC runs ({n_ok}/{n_runs}).')
            return

        stds = {k: float(np.std(v, ddof=1)) if len(v) > 1 else float('nan')
                for k, v in mc_params.items()}
        self.fit_result['params_std'] = stds
        self._apply_result_to_widgets(self.fit_result)
        self.status_label.setText(
            f'MC uncertainty done ({n_ok}/{n_runs} successful runs).')

    # ── Store in file ─────────────────────────────────────────────────────────

    def _store_results(self):
        if self.data is None:
            QMessageBox.warning(self, 'No data', 'Load data first.')
            return
        if self.fit_result is None:
            QMessageBox.warning(self, 'No fit', 'Run Fit first.')
            return
        if not self.data.get('is_nxcansas') or not self.data.get('filepath'):
            QMessageBox.warning(self, 'Text file',
                                'Results can only be stored in HDF5 (NXcanSAS) files.\n'
                                'The current data was loaded from a text file.')
            return

        from pyirena.io.nxcansas_simple_fits import save_simple_fit_results
        filepath = Path(self.data['filepath'])
        try:
            q, I, dI = self._get_filtered_data()
            save_simple_fit_results(
                filepath=filepath,
                result=self.fit_result,
                model_obj=self.model,
                intensity_data=I,
                intensity_error=dI,
            )
            self.status_label.setText(f'Results saved to {filepath.name}')
        except Exception as exc:
            QMessageBox.critical(self, 'Save error', f'Could not save results:\n{exc}')

    # ── Plot update ───────────────────────────────────────────────────────────

    def _update_plots(self, result: dict):
        """Refresh all three graph panels after a fit."""
        q_fit = result.get('q')
        I_model = result.get('I_model')
        residuals = result.get('residuals')

        if q_fit is not None and I_model is not None:
            self.graph_window.plot_fit(q_fit, I_model)

        if q_fit is not None and residuals is not None:
            self.graph_window.plot_residuals(q_fit, residuals)

        # Linearization (use full data range, not just fit range)
        q_all  = self.data['Q'] if self.data else q_fit
        I_all  = self.data['Intensity'] if self.data else I_model
        dI_all = self.data.get('Error') if self.data else None
        if q_all is not None and I_all is not None:
            lin = self.model.linearize(q_all, I_all, dI_all)
        else:
            lin = None
        self.graph_window.plot_linearization(lin)

    # ── State persistence ─────────────────────────────────────────────────────

    def load_state(self):
        """Restore panel state from StateManager."""
        state = self.state_manager.get('simple_fits') or {}
        model_name = state.get('model', 'Guinier')
        if model_name not in MODEL_NAMES:
            model_name = 'Guinier'

        # Restore model object
        self.model = SimpleFitModel()
        self.model.model = model_name
        self.model.use_complex_bg = bool(state.get('use_complex_bg', False))
        self.model.n_mc_runs = int(state.get('n_mc_runs', 50))

        saved_params = state.get('params', {})
        saved_limits = state.get('param_limits', {})
        if saved_params:
            self.model.params.update(saved_params)
        if saved_limits:
            self.model.limits.update(
                {k: tuple(v) for k, v in saved_limits.items()}
            )

        # Update UI controls (signals blocked to avoid cascading rebuilds)
        self.model_combo.blockSignals(True)
        self.model_combo.setCurrentText(model_name)
        self.model_combo.blockSignals(False)
        self.complex_bg_check.blockSignals(True)
        self.complex_bg_check.setChecked(self.model.use_complex_bg)
        self.complex_bg_check.blockSignals(False)
        self.n_runs_spin.setValue(self.model.n_mc_runs)

        q_min = state.get('q_min')
        q_max = state.get('q_max')
        if q_min is not None:
            self.q_min_edit.setText(str(q_min))
        if q_max is not None:
            self.q_max_edit.setText(str(q_max))

        self._build_param_widgets()

    def save_state(self):
        """Persist current panel state via StateManager."""
        q_min_txt = self.q_min_edit.text().strip()
        q_max_txt = self.q_max_edit.text().strip()
        state = {
            'model': self.model.model,
            'q_min': float(q_min_txt) if q_min_txt else None,
            'q_max': float(q_max_txt) if q_max_txt else None,
            'use_complex_bg': self.model.use_complex_bg,
            'params': dict(self.model.params),
            'param_limits': {k: list(v) for k, v in self.model.limits.items()},
            'n_mc_runs': self.n_runs_spin.value(),
        }
        self.state_manager.update('simple_fits', state)
        self.state_manager.save()

    def closeEvent(self, event):
        self.save_state()
        super().closeEvent(event)
