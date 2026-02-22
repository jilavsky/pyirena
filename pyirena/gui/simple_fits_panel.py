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
        QFileDialog,
    )
    from PySide6.QtCore import Qt, Signal
    from PySide6.QtGui import QDoubleValidator
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
            QPushButton, QLabel, QLineEdit, QComboBox, QCheckBox, QSpinBox,
            QSplitter, QMessageBox, QScrollArea, QGroupBox, QSizePolicy, QFrame,
            QFileDialog,
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
from pyirena.gui.sizes_panel import ScrubbableLineEdit
from pyirena.gui.sas_plot import (
    make_sas_plot, plot_iq_data, plot_iq_model,
    make_cursors, get_cursor_q_range, set_cursor_q_range,
    add_plot_annotation, _SafeInfiniteLine, SASPlotStyle,
)


# ===========================================================================
# SimpleFitsGraphWindow
# ===========================================================================

class SimpleFitsGraphWindow(QWidget):
    """
    Three-panel pyqtgraph display using the standard pyIrena SAS plot style:
      - Row 0 (50%): log-log I(Q) vs Q — data scatter + error bars + model line + cursors
      - Row 1 (10%): residuals (I−model)/err vs Q  (log-x, linear-y)
      - Row 2 (40%): linearization (Guinier/Porod) or 'not available' label

    All data passed to ``plot_data`` / ``plot_fit`` must be in physical
    (linear) units.  ``setLogMode(x=True, y=True)`` is applied to the top two
    plots; pyqtgraph handles the log10 transform internally.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._data_item      = None   # ScatterPlotItem (data points)
        self._error_item     = None   # PlotDataItem (error bars)
        self._fit_item       = None   # PlotDataItem (model curve)
        self._resid_item     = None   # ScatterPlotItem (residuals)
        self._lin_data_item  = None
        self._lin_fit_item   = None
        self._cursor_a       = None   # _SafeInfiniteLine in log10 space
        self._cursor_b       = None
        self._annotation_items: list = []
        self._parent_ref     = parent
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(2)
        self.setLayout(layout)

        self.graphics_layout = pg.GraphicsLayoutWidget()
        self.graphics_layout.setBackground('w')

        # ── I(Q) main plot: log-log, cursors, JPEG export ─────────────────────
        self.main_plot = make_sas_plot(
            self.graphics_layout, row=0, col=0,
            x_label='Q  (Å⁻¹)', y_label='I',
            log_x=True, log_y=True,
            parent_widget=self._parent_ref,
            jpeg_default_name='simple_fits_IQ',
        )

        # ── Residuals plot: log-x, linear-y, JPEG export ──────────────────────
        self.residuals_plot = make_sas_plot(
            self.graphics_layout, row=1, col=0,
            x_label='Q  (Å⁻¹)', y_label='(I−fit)/σ',
            log_x=True, log_y=False,
            x_link=self.main_plot,
            parent_widget=self._parent_ref,
            jpeg_default_name='simple_fits_residuals',
        )
        # Zero line for residuals
        zero = pg.InfiniteLine(
            pos=0, angle=0,
            pen=pg.mkPen('k', width=1, style=Qt.PenStyle.DashLine),
        )
        self.residuals_plot.addItem(zero)
        self._resid_zero = zero

        # ── Linearization plot: fully linear, JPEG export ─────────────────────
        self.lin_plot = self.graphics_layout.addPlot(row=2, col=0)
        self.lin_plot.setLabel('left',   'Y')
        self.lin_plot.setLabel('bottom', 'X')
        self.lin_plot.showGrid(x=True, y=True, alpha=SASPlotStyle.GRID_ALPHA)
        self.lin_plot.getAxis('left').enableAutoSIPrefix(False)
        self.lin_plot.getAxis('bottom').enableAutoSIPrefix(False)
        self.lin_plot.setTitle('Linearization', color='#555555', size='10pt')
        # JPEG export for linearization plot
        _vb_lin = self.lin_plot.getViewBox()
        _vb_lin.menu.addSeparator()
        _lin_action = _vb_lin.menu.addAction("Save graph as JPEG…")
        _lin_action.triggered.connect(
            lambda checked=False: self._save_lin_plot_as_jpeg()
        )

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

    def _save_lin_plot_as_jpeg(self):
        from pyqtgraph.exporters import ImageExporter
        file_path, _ = QFileDialog.getSaveFileName(
            self._parent_ref, 'Save Graph as JPEG',
            str(Path.home() / 'simple_fits_linearization.jpg'),
            'JPEG Images (*.jpg *.jpeg);;All Files (*)',
        )
        if not file_path:
            return
        try:
            exporter = ImageExporter(self.lin_plot)
            exporter.parameters()['width'] = 1600
            exporter.export(file_path)
        except Exception as exc:
            QMessageBox.warning(self._parent_ref, 'Export Failed',
                                f'Could not save image:\n{exc}')

    # ── Data plotting ─────────────────────────────────────────────────────────

    def plot_data(self, q: np.ndarray, I: np.ndarray, dI=None, label: str = 'Data'):
        """Plot I(Q) data scatter + error bars.

        Data must be in **physical (linear) units** — pyqtgraph applies log10
        internally via ``setLogMode``.  Existing data, error, and fit items
        are removed before re-plotting; cursors and annotations are kept.
        """
        # Remove old data / error / fit items without clearing the whole plot
        for attr in ('_data_item', '_error_item', '_fit_item'):
            item = getattr(self, attr, None)
            if item is not None:
                try:
                    self.main_plot.removeItem(item)
                except Exception:
                    pass
            setattr(self, attr, None)

        # Plot new data with error bars
        scatter, error = plot_iq_data(self.main_plot, q, I, dI, label=label)
        self._data_item  = scatter
        self._error_item = error

        # Create cursors on the very first data load only
        if self._cursor_a is None:
            mask = np.isfinite(q) & (q > 0) & np.isfinite(I) & (I > 0)
            if mask.any():
                self._cursor_a, self._cursor_b = make_cursors(
                    self.main_plot, q[mask].min(), q[mask].max()
                )

    def plot_fit(self, q_fit: np.ndarray, I_model: np.ndarray):
        """Overlay the model fit / forward-model curve on the main I(Q) plot.

        Data must be in **physical (linear) units**.
        """
        if self._fit_item is not None:
            try:
                self.main_plot.removeItem(self._fit_item)
            except Exception:
                pass
            self._fit_item = None

        item = plot_iq_model(self.main_plot, q_fit, I_model)
        self._fit_item = item

    def clear_fit(self):
        """Remove the model curve and residuals (called on model/param change)."""
        if self._fit_item is not None:
            try:
                self.main_plot.removeItem(self._fit_item)
            except Exception:
                pass
            self._fit_item = None
        # Clear residuals plot back to just the zero line
        self.residuals_plot.clear()
        zero = pg.InfiniteLine(
            pos=0, angle=0,
            pen=pg.mkPen('k', width=1, style=Qt.PenStyle.DashLine),
        )
        self.residuals_plot.addItem(zero)
        self._resid_zero = zero
        self._resid_item = None
        # Clear linearization
        self.lin_plot.clear()
        self.lin_plot.setTitle('Linearization', color='#555555', size='10pt')

    def plot_residuals(self, q: np.ndarray, residuals: np.ndarray):
        """Plot (I−model)/err residuals.  *q* in linear units; *residuals* linear."""
        self.residuals_plot.clear()
        zero = pg.InfiniteLine(
            pos=0, angle=0,
            pen=pg.mkPen('k', width=1, style=Qt.PenStyle.DashLine),
        )
        self.residuals_plot.addItem(zero)
        self._resid_item = None

        mask = np.isfinite(q) & np.isfinite(residuals) & (q > 0)
        if mask.sum() < 1:
            return
        scatter = pg.ScatterPlotItem(
            x=q[mask], y=residuals[mask],
            pen=None,
            brush=SASPlotStyle.RESID_BRUSH,
            size=SASPlotStyle.RESID_SIZE,
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
            pen=None, brush=SASPlotStyle.DATA_BRUSH, size=5,
        )
        self.lin_plot.addItem(data_item)
        self._lin_data_item = data_item

        mask2 = np.isfinite(X_fit) & np.isfinite(Y_fit)
        fit_item = pg.PlotDataItem(
            x=X_fit[mask2], y=Y_fit[mask2],
            pen=SASPlotStyle.FIT_PEN,
        )
        self.lin_plot.addItem(fit_item)
        self._lin_fit_item = fit_item

    # ── Cursors ───────────────────────────────────────────────────────────────

    def get_cursor_range(self):
        """Return (q_min, q_max) in linear units from cursor positions."""
        return get_cursor_q_range(self._cursor_a, self._cursor_b)

    def set_cursor_range(self, q_min: float, q_max: float):
        """Position cursors at the given Q values (linear units)."""
        self._cursor_a, self._cursor_b = set_cursor_q_range(
            self.main_plot, self._cursor_a, self._cursor_b, q_min, q_max
        )

    # ── Result annotations ────────────────────────────────────────────────────

    def clear_result_annotations(self):
        """Remove all fit-result text annotations from the main plot."""
        for item in list(self._annotation_items):
            try:
                self.main_plot.removeItem(item)
            except Exception:
                pass
        self._annotation_items = []

    def add_result_annotation(self, text: str):
        """Add a fit-result annotation in the lower-left of the main I(Q) plot."""
        item = add_plot_annotation(self.main_plot, text, corner='lower_left')
        self._annotation_items.append(item)


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
        self._param_backup: dict | None = None   # snapshot taken before each fit

        # Dynamic widget references (rebuilt when model changes)
        self._param_value_edits:   dict[str, ScrubbableLineEdit] = {}
        self._param_lo_edits:      dict[str, QLineEdit] = {}
        self._param_hi_edits:      dict[str, QLineEdit] = {}
        self._param_unc_labels:    dict[str, QLabel]    = {}
        self._param_fit_checks:    dict[str, QCheckBox] = {}   # "Fit?" per param
        self._param_grid_widget:   QWidget | None = None
        self._lo_header_lbl:       QLabel | None = None
        self._hi_header_lbl:       QLabel | None = None
        self._fit_col_header_lbl:  QLabel | None = None

        self.init_ui()
        self.load_state()

    # ── UI construction ───────────────────────────────────────────────────────

    def init_ui(self):
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(4, 4, 4, 4)
        self.setLayout(main_layout)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self._create_control_panel())
        self.graph_window = SimpleFitsGraphWindow(parent=self)
        splitter.addWidget(self.graph_window)
        splitter.setSizes([420, 780])
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        main_layout.addWidget(splitter)

        # Status label at bottom
        self.status_label = QLabel('No data loaded.')
        self.status_label.setStyleSheet('font-size: 11px; color: #555;')
        main_layout.addWidget(self.status_label)

    def _create_control_panel(self) -> QWidget:
        panel = QWidget()
        panel.setMinimumWidth(380)
        panel.setMaximumWidth(460)
        layout = QVBoxLayout()
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(5)
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

        # ── Q range (cursor-driven, read-only display) ────────────────────────
        q_box = QGroupBox('Q range for fit')
        q_layout = QGridLayout()
        q_layout.setContentsMargins(6, 4, 6, 4)
        q_box.setLayout(q_layout)

        cursor_hint = QLabel('Drag cursors on I(Q) graph to set Q range')
        cursor_hint.setStyleSheet('font-size: 10px; color: #7f8c8d; font-style: italic;')
        cursor_hint.setWordWrap(True)
        q_layout.addWidget(cursor_hint, 0, 0, 1, 3)

        q_layout.addWidget(QLabel('Q min:'), 1, 0)
        self.q_min_display = QLineEdit()
        self.q_min_display.setReadOnly(True)
        self.q_min_display.setPlaceholderText('(cursor A)')
        self.q_min_display.setStyleSheet(
            'background-color: #ecf0f1; color: #7f8c8d;'
        )
        self.q_min_display.setMaximumWidth(90)
        q_layout.addWidget(self.q_min_display, 1, 1)
        q_layout.addWidget(QLabel('Å⁻¹'), 1, 2)

        q_layout.addWidget(QLabel('Q max:'), 2, 0)
        self.q_max_display = QLineEdit()
        self.q_max_display.setReadOnly(True)
        self.q_max_display.setPlaceholderText('(cursor B)')
        self.q_max_display.setStyleSheet(
            'background-color: #ecf0f1; color: #7f8c8d;'
        )
        self.q_max_display.setMaximumWidth(90)
        q_layout.addWidget(self.q_max_display, 2, 1)
        q_layout.addWidget(QLabel('Å⁻¹'), 2, 2)
        layout.addWidget(q_box)

        # ── Global fitting options ─────────────────────────────────────────────
        options_row = QHBoxLayout()
        self.no_limits_check = QCheckBox('No limits?')
        self.no_limits_check.setToolTip(
            'When checked, ignore all parameter bounds and fit unconstrained.\n'
            'The lo/hi input fields are hidden.'
        )
        self.no_limits_check.stateChanged.connect(self._on_no_limits_changed)
        options_row.addWidget(self.no_limits_check)

        self.complex_bg_check = QCheckBox('Complex background  (A·Q⁻ⁿ + flat)')
        self.complex_bg_check.setChecked(False)
        self.complex_bg_check.stateChanged.connect(self._on_complex_bg_changed)
        options_row.addWidget(self.complex_bg_check)
        layout.addLayout(options_row)

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

        # ── Primary action buttons: [Graph model] [Fit] ───────────────────────
        btn_row1 = QHBoxLayout()

        self.graph_model_btn = QPushButton('Graph model')
        self.graph_model_btn.setMinimumHeight(28)
        self.graph_model_btn.setStyleSheet("""
            QPushButton { background-color: #52c77a; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #3eb56a; }
        """)
        self.graph_model_btn.setToolTip('Compute and display the current model curve without fitting.')
        self.graph_model_btn.clicked.connect(self._graph_model)
        btn_row1.addWidget(self.graph_model_btn)

        self.fit_btn = QPushButton('Fit')
        self.fit_btn.setMinimumHeight(28)
        self.fit_btn.setStyleSheet("""
            QPushButton { background-color: #27ae60; color: white; font-weight: bold; font-size: 13px; }
            QPushButton:hover { background-color: #1e8449; }
        """)
        self.fit_btn.clicked.connect(self._run_fit)
        btn_row1.addWidget(self.fit_btn)
        layout.addLayout(btn_row1)

        # ── Uncertainty: [Passes: [10]] [Calc. Uncertainty (MC)] ─────────────
        btn_row_unc = QHBoxLayout()
        btn_row_unc.addWidget(QLabel('Passes:'))
        self.n_runs_spin = QSpinBox()
        self.n_runs_spin.setRange(1, 500)
        self.n_runs_spin.setValue(10)
        self.n_runs_spin.setMaximumWidth(55)
        self.n_runs_spin.setToolTip('Number of noise-perturbed fits for uncertainty estimation.')
        btn_row_unc.addWidget(self.n_runs_spin)

        self.unc_btn = QPushButton('Calc. Uncertainty (MC)')
        self.unc_btn.setMinimumHeight(26)
        self.unc_btn.setStyleSheet("""
            QPushButton { background-color: #16a085; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #1abc9c; }
        """)
        self.unc_btn.clicked.connect(self._calculate_uncertainty)
        btn_row_unc.addWidget(self.unc_btn)
        layout.addLayout(btn_row_unc)

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

        # ── Output / display buttons ──────────────────────────────────────────
        # Row 1: Results to graph | Revert back
        row_out1 = QHBoxLayout()
        self.results_to_graph_btn = QPushButton('Results to graph')
        self.results_to_graph_btn.setMinimumHeight(26)
        self.results_to_graph_btn.setStyleSheet("""
            QPushButton { background-color: #81c784; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #66bb6a; }
        """)
        self.results_to_graph_btn.setToolTip(
            'Annotate the I(Q) plot with the current fitted parameter values.'
        )
        self.results_to_graph_btn.clicked.connect(self._results_to_graph)
        row_out1.addWidget(self.results_to_graph_btn)

        self.revert_btn = QPushButton('Revert back')
        self.revert_btn.setMinimumHeight(26)
        self.revert_btn.setStyleSheet("""
            QPushButton { background-color: #e67e22; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #f39c12; }
        """)
        self.revert_btn.setToolTip(
            'Restore the parameter values that were in place before the last fit.'
        )
        self.revert_btn.clicked.connect(self._revert_to_backup)
        row_out1.addWidget(self.revert_btn)
        layout.addLayout(row_out1)

        # Row 2: Save State | Store in File
        row_out2 = QHBoxLayout()
        self.save_state_btn = QPushButton('Save state')
        self.save_state_btn.setMinimumHeight(26)
        self.save_state_btn.setStyleSheet("""
            QPushButton { background-color: #3498db; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #2980b9; }
        """)
        self.save_state_btn.setToolTip(
            'Save current model choice and parameter values to the state file.'
        )
        self.save_state_btn.clicked.connect(self._save_state_explicit)
        row_out2.addWidget(self.save_state_btn)

        self.store_btn = QPushButton('Store in File')
        self.store_btn.setMinimumHeight(26)
        self.store_btn.setStyleSheet('background-color: lightgreen;')
        self.store_btn.setToolTip('Save fit results to the HDF5 (NXcanSAS) file.')
        self.store_btn.clicked.connect(self._store_results)
        row_out2.addWidget(self.store_btn)
        layout.addLayout(row_out2)

        # Row 3: Export | Import parameters
        row_out3 = QHBoxLayout()
        self.export_btn = QPushButton('Export parameters')
        self.export_btn.setMinimumHeight(26)
        self.export_btn.setStyleSheet('background-color: lightgreen;')
        self.export_btn.setToolTip(
            'Export current parameters to a pyIrena JSON configuration file.'
        )
        self.export_btn.clicked.connect(self._export_parameters)
        row_out3.addWidget(self.export_btn)

        self.import_btn = QPushButton('Import parameters')
        self.import_btn.setMinimumHeight(26)
        self.import_btn.setStyleSheet('background-color: lightgreen;')
        self.import_btn.setToolTip(
            'Import parameters from a pyIrena JSON configuration file.'
        )
        self.import_btn.clicked.connect(self._import_parameters)
        row_out3.addWidget(self.import_btn)
        layout.addLayout(row_out3)

        # Row 4: Reset to defaults (full width)
        self.reset_btn = QPushButton('Reset to defaults')
        self.reset_btn.setMinimumHeight(26)
        self.reset_btn.setStyleSheet("""
            QPushButton { background-color: #e67e22; color: white; font-weight: bold; }
            QPushButton:hover { background-color: #f39c12; }
        """)
        self.reset_btn.setToolTip(
            'Reset all parameters for the current model to their registry defaults.'
        )
        self.reset_btn.clicked.connect(self._reset_to_defaults)
        layout.addWidget(self.reset_btn)

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
        self._param_fit_checks.clear()
        self._lo_header_lbl = None
        self._hi_header_lbl = None
        self._fit_col_header_lbl = None

        container = QWidget()
        grid = QGridLayout()
        grid.setContentsMargins(2, 2, 2, 2)
        grid.setSpacing(3)
        container.setLayout(grid)

        no_limits = self.no_limits_check.isChecked() if hasattr(self, 'no_limits_check') else False

        # Header row: [Fit? | Parameter | Value | lo | hi | ±std]
        headers = [('Fit?', 0), ('Parameter', 1), ('Value', 2),
                   ('lo', 3), ('hi', 4), ('± std', 5)]
        for text, col in headers:
            lbl = QLabel(text)
            lbl.setStyleSheet('font-weight: bold; font-size: 10px;')
            lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            lbl.setVisible(not (no_limits and text in ('lo', 'hi')))
            grid.addWidget(lbl, 0, col)
            if text == 'lo':
                self._lo_header_lbl = lbl
            elif text == 'hi':
                self._hi_header_lbl = lbl
            elif text == 'Fit?':
                self._fit_col_header_lbl = lbl

        entry = MODEL_REGISTRY[self.model.model]
        n_model_params = len(entry['params'])
        use_bg = self.model.use_complex_bg and entry['complex_bg']

        param_specs = list(entry['params'])
        if use_bg:
            from pyirena.core.simple_fits import _BG_PARAMS
            param_specs = param_specs + list(_BG_PARAMS)

        # Saved "Fit?" states from previous build or state
        saved_fixed = getattr(self, '_saved_param_fixed', {})

        row = 1  # header is at row 0
        for i, (name, default, lo_def, hi_def) in enumerate(param_specs):
            # Insert separator before the first BG param
            if use_bg and i == n_model_params:
                sep = QFrame()
                sep.setFrameShape(QFrame.Shape.HLine)
                sep.setFrameShadow(QFrame.Shadow.Sunken)
                grid.addWidget(sep, row, 0, 1, 6)
                row += 1

            val   = self.model.params.get(name, default)
            lo_v, hi_v = self.model.limits.get(name, (lo_def, hi_def))

            # Col 0: Fit? checkbox (default = True = fitted)
            fit_chk = QCheckBox()
            fit_chk.setChecked(not saved_fixed.get(name, False))
            fit_chk.setToolTip(f'Fit {name}?  Uncheck to hold fixed during fitting.')
            grid.addWidget(fit_chk, row, 0, alignment=Qt.AlignmentFlag.AlignCenter)
            self._param_fit_checks[name] = fit_chk

            # Col 1: Name label
            name_lbl = QLabel(name + ':')
            name_lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            name_lbl.setStyleSheet('font-size: 11px;')
            grid.addWidget(name_lbl, row, 1)

            # Col 2: Value edit (scrubable)
            val_edit = ScrubbableLineEdit(f'{val:.6g}')
            val_edit.setValidator(QDoubleValidator(-1e30, 1e30, 10))
            val_edit.setMaximumWidth(85)
            val_edit.setMinimumWidth(70)
            val_edit.editingFinished.connect(self._auto_graph_model)
            grid.addWidget(val_edit, row, 2)
            self._param_value_edits[name] = val_edit

            # Col 3: lo bound
            lo_edit = QLineEdit('' if lo_v is None else f'{lo_v:.4g}')
            lo_edit.setPlaceholderText('−∞')
            lo_edit.setValidator(QDoubleValidator(-1e30, 1e30, 8))
            lo_edit.setMaximumWidth(65)
            lo_edit.setMinimumWidth(50)
            lo_edit.setVisible(not no_limits)
            grid.addWidget(lo_edit, row, 3)
            self._param_lo_edits[name] = lo_edit

            # Col 4: hi bound
            hi_edit = QLineEdit('' if hi_v is None else f'{hi_v:.4g}')
            hi_edit.setPlaceholderText('+∞')
            hi_edit.setValidator(QDoubleValidator(-1e30, 1e30, 8))
            hi_edit.setMaximumWidth(65)
            hi_edit.setMinimumWidth(50)
            hi_edit.setVisible(not no_limits)
            grid.addWidget(hi_edit, row, 4)
            self._param_hi_edits[name] = hi_edit

            # Col 5: ±std label
            unc_lbl = QLabel('')
            unc_lbl.setStyleSheet('font-size: 10px; color: #666;')
            unc_lbl.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            grid.addWidget(unc_lbl, row, 5)
            self._param_unc_labels[name] = unc_lbl

            row += 1

        grid.setColumnStretch(5, 1)
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
        self._saved_param_fixed = {}   # reset "Fit?" state for new model
        self._build_param_widgets()
        self.fit_result = None
        self.chi2_label.setText('—')
        self.rchi2_label.setText('—')
        # Clear stale model curve and annotations from previous model
        self.graph_window.clear_fit()
        self.graph_window.clear_result_annotations()

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
        self.graph_window.clear_fit()
        self.graph_window.clear_result_annotations()

    def _on_no_limits_changed(self, state):
        """Show/hide the lo/hi bound columns when 'No limits?' is toggled."""
        no_limits = bool(state)
        if self._lo_header_lbl:
            self._lo_header_lbl.setVisible(not no_limits)
        if self._hi_header_lbl:
            self._hi_header_lbl.setVisible(not no_limits)
        for edit in self._param_lo_edits.values():
            edit.setVisible(not no_limits)
        for edit in self._param_hi_edits.values():
            edit.setVisible(not no_limits)

    def _on_param_edited(self):
        """Read edited values back into self.model.params."""
        for name, edit in self._param_value_edits.items():
            txt = edit.text().strip()
            if txt:
                try:
                    self.model.params[name] = float(txt)
                except ValueError:
                    pass

    def _auto_graph_model(self):
        """Sync parameter values and auto-redisplay model curve after edit."""
        self._on_param_edited()
        # Stale annotations are now wrong — clear them before re-graphing
        self.graph_window.clear_result_annotations()
        if self.data is not None:
            self._graph_model()

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
        self._update_q_display()

    # ── Q range helpers ───────────────────────────────────────────────────────

    def _update_q_display(self):
        """Read cursor positions and update the read-only Q range display fields."""
        q_min, q_max = self.graph_window.get_cursor_range()
        if q_min is not None:
            self.q_min_display.setText(f'{q_min:.6g}')
        if q_max is not None:
            self.q_max_display.setText(f'{q_max:.6g}')

    def _get_q_range(self):
        """Return (q_min, q_max) from cursor positions, updating the display."""
        q_min, q_max = self.graph_window.get_cursor_range()
        # Update display fields
        if q_min is not None:
            self.q_min_display.setText(f'{q_min:.6g}')
        if q_max is not None:
            self.q_max_display.setText(f'{q_max:.6g}')
        return q_min, q_max

    def _get_filtered_data(self):
        """Return (q, I, dI) arrays filtered to the active Q range (from cursors)."""
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
        dI_f = dI[mask] if dI is not None else None
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

    def _collect_fixed_params(self) -> dict:
        """Return dict of param_name → current_value for unchecked 'Fit?' params."""
        fixed = {}
        for name, chk in self._param_fit_checks.items():
            if not chk.isChecked():
                val_edit = self._param_value_edits.get(name)
                if val_edit:
                    try:
                        fixed[name] = float(val_edit.text())
                    except ValueError:
                        fixed[name] = self.model.params.get(name, 0.0)
                else:
                    fixed[name] = self.model.params.get(name, 0.0)
        return fixed

    def _apply_result_to_widgets(self, result: dict):
        """Write fitted parameter values and ± uncertainties back to widgets."""
        params = result.get('params', {})
        stds   = result.get('params_std', {})
        for name, edit in self._param_value_edits.items():
            if name in params:
                edit.setText(f'{params[name]:.6g}')
        for name, lbl in self._param_unc_labels.items():
            std = stds.get(name)
            if std is not None and np.isfinite(std) and std > 0:
                lbl.setText(f'±{std:.3g}')
            else:
                lbl.setText('')

    # ── Graph model (no fitting) ──────────────────────────────────────────────

    def _graph_model(self):
        """Compute and display the current model curve without fitting."""
        if self.data is None:
            return

        self._collect_model()   # sync widget values → model.params

        q_all = self.data['Q']
        I_all = self.data['Intensity']
        mask_all = np.isfinite(q_all) & np.isfinite(I_all) & (q_all > 0) & (I_all > 0)
        q_plot = q_all[mask_all]

        # Restrict to cursor Q range if available
        q_min, q_max = self.graph_window.get_cursor_range()
        if q_min is not None:
            q_plot = q_plot[q_plot >= q_min]
        if q_max is not None:
            q_plot = q_plot[q_plot <= q_max]

        if len(q_plot) < 2:
            return

        try:
            I_model = self.model.compute(q_plot)
        except Exception as exc:
            self.status_label.setText(f'Model error: {exc}')
            return

        valid = np.isfinite(I_model) & (I_model > 0)
        if valid.sum() < 2:
            self.status_label.setText('Model returned no valid (positive) values.')
            return

        self.graph_window.plot_fit(q_plot[valid], I_model[valid])
        self.status_label.setText(f'Model: {self.model.model}  (not fitted)')

    # ── Fit ───────────────────────────────────────────────────────────────────

    def _run_fit(self):
        if self.data is None:
            QMessageBox.warning(self, 'No data', 'Load data first.')
            return

        # Read cursor positions into display fields and get Q range
        self._get_q_range()

        # Snapshot parameters before fitting so "Revert back" can restore them
        self._param_backup = {
            'model': self.model.model,
            'use_complex_bg': self.model.use_complex_bg,
            'params': dict(self.model.params),
            'limits': {k: tuple(v) for k, v in self.model.limits.items()},
        }

        model = self._collect_model()
        fixed_params = self._collect_fixed_params()
        no_limits = self.no_limits_check.isChecked()

        q, I, dI = self._get_filtered_data()
        if q is None or len(q) < 2:
            QMessageBox.warning(self, 'Too few points',
                                'Not enough data points in the selected Q range.')
            return

        result = model.fit(q, I, dI,
                           fixed_params=fixed_params if fixed_params else None,
                           no_limits=no_limits)
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
        fixed_params = self._collect_fixed_params()
        no_limits = self.no_limits_check.isChecked()
        q, I, dI = self._get_filtered_data()
        n_runs = self.n_runs_spin.value()
        dI_safe = dI if dI is not None else np.maximum(I * 0.05, 1e-30)

        mc_params: dict[str, list] = {k: [] for k in self.fit_result['params']}
        n_ok = 0
        for _ in range(n_runs):
            I_pert = I + dI_safe * np.random.randn(len(I))
            mc_res = model.fit(q, I_pert, dI,
                               fixed_params=fixed_params if fixed_params else None,
                               no_limits=no_limits)
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

    # ── Helper ────────────────────────────────────────────────────────────────

    def _get_data_folder(self) -> str:
        """Return the folder of the currently loaded file, or home dir."""
        if self.data and self.data.get('filepath'):
            return str(Path(self.data['filepath']).parent)
        return str(Path.home())

    # ── Revert / reset ────────────────────────────────────────────────────────

    def _revert_to_backup(self):
        """Restore parameters to the snapshot taken before the last fit."""
        if self._param_backup is None:
            QMessageBox.warning(self, 'No backup',
                                'No parameter backup available — run a fit first.')
            return
        b = self._param_backup
        if b['model'] != self.model.model:
            self.model_combo.setCurrentText(b['model'])  # triggers _on_model_changed
        self.model.params.update(b['params'])
        self.model.limits.update(b['limits'])
        self._build_param_widgets()
        self.status_label.setText('Reverted to pre-fit parameters.')

    def _reset_to_defaults(self):
        """Reset the active model's parameters to registry defaults."""
        reply = QMessageBox.question(
            self, 'Reset to defaults',
            f'Reset all {self.model.model} parameters to their default values?',
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if reply == QMessageBox.StandardButton.Yes:
            self.model._reset_to_defaults()
            self._saved_param_fixed = {}
            self._build_param_widgets()
            self.status_label.setText('Parameters reset to defaults.')

    # ── Save state (explicit button) ─────────────────────────────────────────

    def _save_state_explicit(self):
        """Save state and show a confirmation dialog."""
        self.save_state()
        QMessageBox.information(self, 'State saved',
                                'Current model and parameters saved successfully.')

    # ── Export / import parameters ────────────────────────────────────────────

    def _export_parameters(self):
        """Export current parameters to a pyIrena JSON configuration file."""
        import json, datetime
        try:
            from pyirena import __version__ as _version
        except Exception:
            _version = 'unknown'

        default_path = str(Path(self._get_data_folder()) / 'pyirena_config.json')
        file_path, _ = QFileDialog.getSaveFileName(
            self, 'Export pyIrena Configuration', default_path,
            'pyIrena Config (*.json);;All Files (*)',
        )
        if not file_path:
            return
        file_path = Path(file_path)

        config: dict = {}
        if file_path.exists():
            try:
                with open(file_path, 'r') as f:
                    config = json.load(f)
            except Exception:
                config = {}
            if '_pyirena_config' not in config:
                QMessageBox.warning(
                    self, 'Not a pyIrena file',
                    f'The selected file is not a pyIrena configuration file:\n{file_path}\n\n'
                    'Choose a different file or enter a new filename.',
                )
                return
            if 'simple_fits' in config:
                reply = QMessageBox.question(
                    self, 'Overwrite Simple Fits parameters?',
                    f'File already contains Simple Fits parameters:\n{file_path}\n\n'
                    'Overwrite them?',
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
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
        config['_pyirena_config']['written_by'] = f'pyIrena {_version}'

        self.save_state()
        config['simple_fits'] = self.state_manager.get('simple_fits')

        try:
            with open(file_path, 'w') as f:
                json.dump(config, f, indent=2)
        except Exception as exc:
            QMessageBox.warning(self, 'Export failed', f'Could not write file:\n{exc}')
            return

        self.status_label.setText(f'Parameters exported to: {file_path.name}')

    def _import_parameters(self):
        """Import parameters from a pyIrena JSON configuration file."""
        import json

        file_path, _ = QFileDialog.getOpenFileName(
            self, 'Import pyIrena Configuration', self._get_data_folder(),
            'pyIrena Config (*.json);;All Files (*)',
        )
        if not file_path:
            return
        file_path = Path(file_path)

        try:
            with open(file_path, 'r') as f:
                config = json.load(f)
        except Exception as exc:
            QMessageBox.warning(self, 'Import failed', f'Could not read file:\n{exc}')
            return

        if '_pyirena_config' not in config:
            QMessageBox.warning(
                self, 'Not a pyIrena file',
                f'The selected file is not a pyIrena configuration file:\n{file_path}',
            )
            return
        if 'simple_fits' not in config:
            QMessageBox.warning(
                self, 'No Simple Fits parameters',
                f'The file does not contain a Simple Fits parameter group:\n{file_path}',
            )
            return

        sf_state = config['simple_fits']
        self.state_manager.update('simple_fits', sf_state)
        self.load_state()
        written_by = config['_pyirena_config'].get('written_by', 'unknown version')
        self.status_label.setText(
            f'Parameters imported from: {file_path.name}  (written by {written_by})')

    # ── Results to graph ──────────────────────────────────────────────────────

    def _results_to_graph(self):
        """Annotate the I(Q) plot with the current fitted parameter values."""
        if self.fit_result is None or not self.fit_result.get('success'):
            QMessageBox.warning(self, 'No fit results',
                                'Run a fit first to generate results.')
            return

        params  = self.fit_result.get('params', {})
        stds    = self.fit_result.get('params_std', {})
        derived = self.fit_result.get('derived', {})
        chi2    = self.fit_result.get('reduced_chi2')

        lines = [f'Model: {self.model.model}']
        if chi2 is not None:
            lines.append(f'χ²_red = {chi2:.4g}')
        for name, val in params.items():
            std = stds.get(name)
            if std is not None and np.isfinite(std) and std > 0:
                lines.append(f'{name} = {val:.4g} ± {std:.3g}')
            else:
                lines.append(f'{name} = {val:.4g}')
        for name, val in derived.items():
            lines.append(f'{name} = {val:.4g}  [derived]')

        self.graph_window.clear_result_annotations()
        self.graph_window.add_result_annotation('\n'.join(lines))
        self.status_label.setText('Results annotated on graph.')

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
        self.model.n_mc_runs = int(state.get('n_mc_runs', 10))

        saved_params = state.get('params', {})
        saved_limits = state.get('param_limits', {})
        if saved_params:
            self.model.params.update(saved_params)
        if saved_limits:
            self.model.limits.update(
                {k: tuple(v) for k, v in saved_limits.items()}
            )

        # Restore "Fit?" (fixed) states — stored as {name: True if fixed}
        self._saved_param_fixed = state.get('param_fixed', {})

        # Update UI controls (signals blocked to avoid cascading rebuilds)
        self.model_combo.blockSignals(True)
        self.model_combo.setCurrentText(model_name)
        self.model_combo.blockSignals(False)

        self.complex_bg_check.blockSignals(True)
        self.complex_bg_check.setChecked(self.model.use_complex_bg)
        self.complex_bg_check.blockSignals(False)

        no_limits = bool(state.get('no_limits', False))
        self.no_limits_check.blockSignals(True)
        self.no_limits_check.setChecked(no_limits)
        self.no_limits_check.blockSignals(False)

        self.n_runs_spin.setValue(self.model.n_mc_runs)

        # Restore cursor positions if saved
        q_min = state.get('q_min')
        q_max = state.get('q_max')
        if q_min is not None and q_max is not None:
            try:
                self.graph_window.set_cursor_range(float(q_min), float(q_max))
                self.q_min_display.setText(f'{float(q_min):.6g}')
                self.q_max_display.setText(f'{float(q_max):.6g}')
            except Exception:
                pass

        self._build_param_widgets()

    def save_state(self):
        """Persist current panel state via StateManager."""
        # Read cursor Q range (and update display)
        q_min, q_max = self.graph_window.get_cursor_range()
        if q_min is not None:
            self.q_min_display.setText(f'{q_min:.6g}')
        if q_max is not None:
            self.q_max_display.setText(f'{q_max:.6g}')

        # Collect which params are fixed (Fit? unchecked)
        param_fixed = {
            name: not chk.isChecked()
            for name, chk in self._param_fit_checks.items()
        }

        state = {
            'model': self.model.model,
            'q_min': q_min,
            'q_max': q_max,
            'use_complex_bg': self.model.use_complex_bg,
            'no_limits': self.no_limits_check.isChecked(),
            'params': dict(self.model.params),
            'param_limits': {k: list(v) for k, v in self.model.limits.items()},
            'param_fixed': param_fixed,
            'n_mc_runs': self.n_runs_spin.value(),
        }
        self.state_manager.update('simple_fits', state)
        self.state_manager.save()

    def closeEvent(self, event):
        self.save_state()
        super().closeEvent(event)
