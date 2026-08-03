"""
Simple Fits GUI panel for pyIrena.

Provides ``SimpleFitsGraphWindow`` (three-panel pyqtgraph display: I(Q),
residuals, linearization) and ``SimpleFitsPanel`` (controls + graph) for
interactive single-model fitting of SAS data.
"""

from __future__ import annotations
import logging

log = logging.getLogger(__name__)


from pyirena.gui._qt import (
    QCheckBox, QComboBox, QDesktopServices, QDoubleValidator, QFileDialog, QFrame, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPushButton, QScrollArea, QSizePolicy, QSpinBox, QSplitter, QUrl, QVBoxLayout, QWidget, Qt,
)

import numpy as np
from pathlib import Path

import pyqtgraph as pg

from pyirena.core.simple_fits import (
    SimpleFitModel, MODEL_NAMES, MODEL_REGISTRY, _resolve_model_name,
)
from pyirena.gui.data_loading import DataFileLoaderRow
from pyirena.gui.slit_smearing_ui import SlitSmearingMixin
from pyirena.state.state_manager import StateManager
from pyirena.gui.sizes_panel import ScrubbableLineEdit
from pyirena.gui.sas_plot import (
    make_sas_plot, plot_iq_data, plot_iq_model,
    make_cursors, get_cursor_q_range, set_cursor_q_range,
    add_plot_annotation, SASPlotStyle,
)

# Friendly display labels for the complex-background parameters.  The dict
# *keys* stay ``BG_*`` (so persistence, HDF5 I/O and ``linearize()`` are
# unaffected); only the label shown in the parameter grid uses these symbols,
# matching the Unified Fit convention (B = prefactor, P = exponent).
_BG_DISPLAY = {'BG_B': 'B', 'BG_P': 'P', 'BG_flat': 'flat'}


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
        # Invariant display: corrected data (main vb) + running integral
        # on a second ViewBox tied to the right axis (Igor-style overlay).
        self._inv_corr_item  = None   # PlotDataItem, bg-corrected data (green)
        self._inv_int_item   = None   # PlotDataItem, running integral (right axis)
        self._right_vb       = None   # pg.ViewBox for the right axis
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
            x_label='Q  (Å⁻¹)', y_label="Residuals r' (rescaled)",
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
        self._lin_row_stretch = 4      # remembered so hide/show can restore it
        self._lin_visible = True
        self._resid_row_stretch = 1    # remembered so hide/show can restore it
        self._resid_visible = True

        # Stretch=1 so the graphics widget expands vertically when the
        # window grows (default is 0 = stay at preferred size).  Other tools
        # (sizes_panel, unified_fit) get this for free because their
        # graph windows hold a QTabWidget, which has Expanding policy by
        # default; pg.GraphicsLayoutWidget needs the explicit factor.
        self.graphics_layout.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        layout.addWidget(self.graphics_layout, stretch=1)

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
                    log.debug("suppressed exception", exc_info=True)
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
                log.debug("suppressed exception", exc_info=True)
            self._fit_item = None

        item = plot_iq_model(self.main_plot, q_fit, I_model)
        self._fit_item = item

    # ── Invariant overlay (corrected data + running integral, right axis) ────

    def _ensure_right_axis(self):
        """Create (once) a second ViewBox linked to the right axis and make
        the axis visible and labelled.

        The main plot is in log-log mode, so its ViewBox coordinates are
        log10 units; items added to the right ViewBox must therefore be
        given log10(x), log10(y) data, and the right AxisItem is switched
        to log tick rendering.

        The show/label part runs on EVERY call (not only on first creation):
        switching to another model hides the right axis via
        ``set_right_axis_visible(False)``, and returning to the Invariant
        must bring both the axis and its label back.
        """
        if self._right_vb is None:
            vb = pg.ViewBox()
            self._right_vb = vb
            self.main_plot.scene().addItem(vb)
            self.main_plot.getAxis('right').linkToView(vb)
            vb.setXLink(self.main_plot.vb)
            self.main_plot.vb.sigResized.connect(self._sync_right_vb)
        # (Re)show + (re)label every time — hideAxis() may have been called.
        self.main_plot.showAxis('right')
        axis = self.main_plot.getAxis('right')
        try:
            axis.setLogMode(True)
        except Exception:
            log.debug("right axis setLogMode failed", exc_info=True)
        axis.setLabel('∫q²I dq  (cm⁻⁴)')
        self._sync_right_vb()

    def _sync_right_vb(self):
        if self._right_vb is not None:
            self._right_vb.setGeometry(self.main_plot.vb.sceneBoundingRect())

    def plot_invariant(self, q_int, running_integral, q=None, I_corrected=None):
        """Overlay the invariant running integral (right axis, log) and the
        background-corrected data (green, main axis).

        All inputs in physical (linear) units; log10 applied here because the
        right ViewBox bypasses pyqtgraph's automatic log transform.
        """
        self.clear_invariant()

        # Background-corrected data on the main (left) axis
        if q is not None and I_corrected is not None:
            q_arr = np.asarray(q, float)
            Ic = np.asarray(I_corrected, float)
            m = np.isfinite(q_arr) & np.isfinite(Ic) & (q_arr > 0) & (Ic > 0)
            if m.sum() > 1:
                self._inv_corr_item = self.main_plot.plot(
                    q_arr[m], Ic[m],
                    pen=pg.mkPen(0, 150, 0, width=1),
                    name='I − bg',
                )

        # Running integral on the right axis
        qi = np.asarray(q_int, float)
        ri = np.asarray(running_integral, float)
        m = np.isfinite(qi) & np.isfinite(ri) & (qi > 0) & (ri > 0)
        if m.sum() > 1:
            self._ensure_right_axis()
            item = pg.PlotDataItem(
                np.log10(qi[m]), np.log10(ri[m]),
                pen=pg.mkPen(200, 0, 0, width=2),
            )
            self._right_vb.addItem(item)
            self._inv_int_item = item
            self._right_vb.setYRange(
                float(np.log10(ri[m]).min()), float(np.log10(ri[m]).max()),
                padding=0.05,
            )

    def clear_invariant(self):
        """Remove invariant overlay items (corrected data + running integral)."""
        if self._inv_corr_item is not None:
            try:
                self.main_plot.removeItem(self._inv_corr_item)
            except Exception:
                log.debug("suppressed exception", exc_info=True)
            self._inv_corr_item = None
        if self._inv_int_item is not None and self._right_vb is not None:
            try:
                self._right_vb.removeItem(self._inv_int_item)
            except Exception:
                log.debug("suppressed exception", exc_info=True)
            self._inv_int_item = None

    def set_residuals_visible(self, visible: bool):
        """Show or hide the residuals panel.

        Calculation models (Invariant) have no residuals; hiding the row lets
        the I(Q) plot use the space.  Mirrors set_linearization_visible().
        """
        visible = bool(visible)
        if visible == self._resid_visible:
            return
        self._resid_visible = visible
        ci = self.graphics_layout.ci
        if visible:
            self.residuals_plot.show()
            ci.layout.setRowStretchFactor(1, self._resid_row_stretch)
        else:
            self.residuals_plot.hide()
            ci.layout.setRowStretchFactor(1, 0)

    def set_right_axis_visible(self, visible: bool):
        """Show/hide the right (invariant integral) axis."""
        if visible:
            self.main_plot.showAxis('right')
        else:
            try:
                self.main_plot.hideAxis('right')
            except Exception:
                log.debug("suppressed exception", exc_info=True)

    def clear_fit(self):
        """Remove the model curve and residuals (called on model/param change)."""
        if self._fit_item is not None:
            try:
                self.main_plot.removeItem(self._fit_item)
            except Exception:
                log.debug("suppressed exception", exc_info=True)
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
        # Clear invariant overlay
        self.clear_invariant()

    def plot_residuals(self, q: np.ndarray, residuals: np.ndarray):
        """Plot (I−model)/err residuals.  *q* in linear units; *residuals* linear.

        Uses ``residuals_plot.plot()`` (returns PlotDataItem) rather than a
        standalone ScatterPlotItem so that pyqtgraph applies the log-x transform
        correctly via its internal dataItems bookkeeping.
        """
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
        item = self.residuals_plot.plot(
            q[mask], residuals[mask],
            pen=None,
            symbol='o',
            symbolSize=SASPlotStyle.RESID_SIZE,
            symbolPen=None,
            symbolBrush=SASPlotStyle.RESID_BRUSH,
        )
        self._resid_item = item

    def plot_linearization(self, lin_result, q_min=None, q_max=None):
        """Plot linearized data + fit.  lin_result is from SimpleFitModel.linearize().

        Parameters
        ----------
        lin_result : dict or None
            Output of ``SimpleFitModel.linearize()``.
        q_min, q_max : float or None
            Q-range boundaries (linear units) from the cursor positions.
            Points whose original Q lies within [q_min, q_max] are plotted
            in the standard dark colour; points outside that range are plotted
            in light grey so the fitting region is immediately visible.
        """
        self.lin_plot.clear()
        self._lin_data_item = None
        self._lin_fit_item  = None

        if lin_result is None:
            self.lin_plot.setTitle('Linearization — not available for this model',
                                   color='#888888', size='9pt')
            return

        X     = lin_result['X']
        Y     = lin_result['Y']
        X_fit = lin_result['X_fit']
        Y_fit = lin_result['Y_fit']
        q_raw = lin_result.get('q')   # original Q values from linearize()

        _title = (
            f"{lin_result['y_label']} vs {lin_result['x_label']}   "
            f"slope={lin_result['slope']:.4g}  intercept={lin_result['intercept']:.4g}  "
            f"R²={lin_result['r_squared']:.4f}"
        )
        if lin_result.get('best_effort_smeared'):
            _title = "⚠ ideal-space (best effort) — data are slit smeared\n" + _title
        self.lin_plot.setTitle(
            _title,
            color='#333333', size='9pt',
        )
        self.lin_plot.setLabel('left',   lin_result['y_label'])
        self.lin_plot.setLabel('bottom', lin_result['x_label'])

        good = np.isfinite(X) & np.isfinite(Y)

        # Determine which data points lie inside the fitted Q range
        if q_raw is not None and (q_min is not None or q_max is not None):
            in_range = np.ones(len(q_raw), dtype=bool)
            if q_min is not None:
                in_range &= (q_raw >= q_min)
            if q_max is not None:
                in_range &= (q_raw <= q_max)
        else:
            in_range = np.ones(len(X), dtype=bool)

        # ── Out-of-range points: light grey, smaller ──────────────────────────
        out_mask = good & ~in_range
        if out_mask.any():
            out_item = pg.ScatterPlotItem(
                x=X[out_mask], y=Y[out_mask],
                pen=None,
                brush=pg.mkBrush(180, 180, 180, 140),
                size=4,
            )
            self.lin_plot.addItem(out_item)

        # ── In-range points: dark, standard size ──────────────────────────────
        in_mask = good & in_range
        if in_mask.any():
            data_item = pg.ScatterPlotItem(
                x=X[in_mask], y=Y[in_mask],
                pen=None, brush=SASPlotStyle.DATA_BRUSH, size=5,
            )
            self.lin_plot.addItem(data_item)
            self._lin_data_item = data_item

        # ── Model line ────────────────────────────────────────────────────────
        mask2 = np.isfinite(X_fit) & np.isfinite(Y_fit)
        fit_item = pg.PlotDataItem(
            x=X_fit[mask2], y=Y_fit[mask2],
            pen=SASPlotStyle.FIT_PEN,
        )
        self.lin_plot.addItem(fit_item)
        self._lin_fit_item = fit_item

        # ── Bound both axes to the cursor Q range (converted to X units) ───────
        # The transform is monotone in Q, so the in-range X span is simply the
        # min/max of the X values whose original Q lies within [q_min, q_max].
        # A robust Y window is taken from the in-range Y values only so tiny
        # high-Q intensities (huge negative ln values) can't dominate the scale.
        self._apply_lin_ranges(X, Y, in_mask)

    def _apply_lin_ranges(self, X, Y, in_mask):
        """Fix the linearization plot to the in-range data window.

        *in_mask* selects the points inside the cursor Q range (and finite).
        Falls back to all finite points when the cursor range excludes
        everything.  Disables auto-range so the view stays put.
        """
        sel = in_mask
        if not np.any(sel):
            sel = np.isfinite(X) & np.isfinite(Y)
        if not np.any(sel):
            return

        Xs = X[sel]
        Ys = Y[sel]

        x_lo, x_hi = float(np.min(Xs)), float(np.max(Xs))
        if x_hi <= x_lo:
            x_hi = x_lo + max(abs(x_lo) * 1e-3, 1e-30)

        # Robust linear Y window (1st–99th percentile) + 8% padding.
        y_lo, y_hi = np.nanpercentile(Ys, [1.0, 99.0])
        span = y_hi - y_lo
        if span <= 0:
            span = max(abs(y_hi), 1.0) * 0.1
        pad = span * 0.08
        y_lo -= pad
        y_hi += pad

        self.lin_plot.enableAutoRange(enable=False)
        self.lin_plot.setXRange(x_lo, x_hi, padding=0.02)
        self.lin_plot.setYRange(y_lo, y_hi, padding=0)

    def set_linearization_visible(self, visible: bool):
        """Show or hide the bottom linearization panel.

        When hidden, the row's stretch factor is set to 0 so the I(Q) and
        residuals plots expand to fill the freed vertical space.
        """
        visible = bool(visible)
        if visible == self._lin_visible:
            return
        self._lin_visible = visible
        ci = self.graphics_layout.ci
        if visible:
            self.lin_plot.show()
            ci.layout.setRowStretchFactor(2, self._lin_row_stretch)
        else:
            self.lin_plot.hide()
            ci.layout.setRowStretchFactor(2, 0)

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
                log.debug("suppressed exception", exc_info=True)
        self._annotation_items = []

    def add_result_annotation(self, text: str):
        """Add a fit-result annotation in the lower-left of the main I(Q) plot."""
        item = add_plot_annotation(self.main_plot, text, corner='lower_left')
        self._annotation_items.append(item)


# ===========================================================================
# SimpleFitsPanel
# ===========================================================================

class SimpleFitsPanel(SlitSmearingMixin, QWidget):
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
        # stretch=1 so the splitter (graph + controls) absorbs all extra
        # vertical space when the window grows; the status label below
        # stays at its preferred (small) height.  Without this, the layout
        # leaves the splitter at preferred size on resize.
        main_layout.addWidget(splitter, stretch=1)

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

        # ── Title + Help ──────────────────────────────────────────────────────
        title_row = QHBoxLayout()
        title_lbl = QLabel("Simple Fits")
        title_lbl.setStyleSheet("font-size: 14px; font-weight: bold; color: #2c3e50;")
        title_row.addWidget(title_lbl)
        title_row.addStretch()
        _help_btn = QPushButton("? Help")
        _help_btn.setFixedSize(60, 22)
        _help_btn.setStyleSheet(
            "QPushButton{background:#c0392b;color:white;font-size:11px;border-radius:3px;}"
            "QPushButton:hover{background:#e74c3c;}"
        )
        _help_btn.setToolTip("Open online documentation in your browser")
        _help_btn.clicked.connect(
            lambda: QDesktopServices.openUrl(QUrl(
                "https://github.com/jilavsky/pyirena/blob/main/docs/simple_fits_gui.md"
            ))
        )
        title_row.addWidget(_help_btn)
        layout.addLayout(title_row)

        # ── Data file loader ──────────────────────────────────────────────────
        self.data_loader = DataFileLoaderRow(state_manager=self.state_manager)
        self.data_loader.data_loaded.connect(self._on_loader_data_loaded)
        layout.addWidget(self.data_loader)
        self._build_slit_row(layout)

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

        self.complex_bg_check = QCheckBox('Complex background  (B·Q⁻ᴾ + flat)')
        self.complex_bg_check.setChecked(False)
        self.complex_bg_check.stateChanged.connect(self._on_complex_bg_changed)
        options_row.addWidget(self.complex_bg_check)
        layout.addLayout(options_row)

        # ── Invariant-only options (hidden unless model is a calculation) ─────
        inv_col = QVBoxLayout()
        inv_col.setContentsMargins(0, 0, 0, 0)
        inv_col.setSpacing(2)

        inv_row1 = QHBoxLayout()
        self.porod_tail_check = QCheckBox('Extend by Porod tail (Kp/Qmax)')
        self.porod_tail_check.setToolTip(
            'Extend the invariant integral beyond QmaxUsed by the analytic\n'
            'Porod tail ∫Kp·q⁻²dq = Kp/Qmax.  Kp is estimated as the median\n'
            'of (I−bg)·q⁴ over the last half-decade before QmaxUsed.\n'
            'Compensates for the high-Q truncation of the measurement\n'
            '(Igor Irena does not apply this correction).'
        )
        self.porod_tail_check.stateChanged.connect(self._on_porod_tail_changed)
        inv_row1.addWidget(self.porod_tail_check)
        inv_row1.addStretch()
        inv_col.addLayout(inv_row1)

        inv_row2 = QHBoxLayout()
        self.bg_refit_check = QCheckBox('Refit background from saved ranges')
        self.bg_refit_check.setToolTip(
            'Before every Calculate — in the GUI and in scripted/batch runs —\n'
            're-determine the complex background from the Q ranges remembered\n'
            'when you last used "Fit B/P btwn cursors" / "Fit Flat btwn\n'
            'cursors" (power-law first, then flat).  This makes scripted\n'
            'invariant results robust: the background is refit per file\n'
            'instead of assuming the exported B/P/flat values apply.\n'
            'BG_P is refit only if its "Fit?" box was checked when the\n'
            'range was recorded.  Requires "Complex background".'
        )
        self.bg_refit_check.stateChanged.connect(self._on_bg_refit_changed)
        inv_row2.addWidget(self.bg_refit_check)
        inv_row2.addStretch()
        inv_col.addLayout(inv_row2)

        self._bg_prefit_ranges_lbl = QLabel('')
        self._bg_prefit_ranges_lbl.setStyleSheet(
            'font-size: 10px; color: #7f8c8d; font-style: italic;')
        self._bg_prefit_ranges_lbl.setWordWrap(True)
        inv_col.addWidget(self._bg_prefit_ranges_lbl)

        self._invariant_options_row = QWidget()
        self._invariant_options_row.setLayout(inv_col)
        self._invariant_options_row.setVisible(False)
        layout.addWidget(self._invariant_options_row)

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

        # ── Complex-background prefit helpers (shown only when bg is active) ───
        # Small helper buttons that prefit the power-law (B/P) or flat term over
        # the cursor-selected Q window, so users don't have to guess B by hand.
        # Mirrors the "Fit … btwn cursors" pattern from the Unified Fit tool.
        self._bg_prefit_row = QWidget()
        _bg_prefit_layout = QHBoxLayout(self._bg_prefit_row)
        _bg_prefit_layout.setContentsMargins(0, 0, 0, 0)
        _bg_prefit_layout.setSpacing(4)
        _bg_helper_style = (
            'QPushButton { font-size: 10px; padding: 1px 6px; background-color: #ecf0f1; }'
            'QPushButton:hover { background-color: #dfe4e6; }'
        )

        self.prefit_bp_btn = QPushButton('Fit B/P btwn cursors')
        self.prefit_bp_btn.setMaximumHeight(22)
        self.prefit_bp_btn.setStyleSheet(_bg_helper_style)
        self.prefit_bp_btn.setToolTip(
            'Prefit the background power-law B·Q⁻ᴾ to the data between the two '
            'cursors.\nIf P’s "Fit?" box is unchecked, only B is fit at the '
            'current (model-guided) P.\nResults are written as starting values '
            'for the full Fit.'
        )
        self.prefit_bp_btn.clicked.connect(self._prefit_bg_power_law)
        _bg_prefit_layout.addWidget(self.prefit_bp_btn)

        self.prefit_flat_btn = QPushButton('Fit Flat btwn cursors')
        self.prefit_flat_btn.setMaximumHeight(22)
        self.prefit_flat_btn.setStyleSheet(_bg_helper_style)
        self.prefit_flat_btn.setToolTip(
            'Prefit the flat background term to the data between the two cursors '
            '(median of I − B·Q⁻ᴾ over the range, using the current B and P).'
        )
        self.prefit_flat_btn.clicked.connect(self._prefit_bg_flat)
        _bg_prefit_layout.addWidget(self.prefit_flat_btn)
        _bg_prefit_layout.addStretch()

        self._bg_prefit_row.setVisible(False)
        layout.addWidget(self._bg_prefit_row)

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
        self.fit_btn.setToolTip(
            "Fit the selected model to the loaded data using the current Q range.\n"
            "Checked parameters are varied; unchecked parameters are held fixed."
        )
        self.fit_btn.clicked.connect(self._run_fit)
        btn_row1.addWidget(self.fit_btn)
        layout.addLayout(btn_row1)

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

        # Row 2: Store in File
        # (State is saved automatically on close — see closeEvent — so an
        # explicit "Save state" button is no longer needed.)
        row_out2 = QHBoxLayout()

        self.store_btn = QPushButton('Store in File')
        self.store_btn.setMinimumHeight(26)
        self.store_btn.setStyleSheet('background-color: lightgreen;')
        self.store_btn.setToolTip(
            'Save fit results to the HDF5 (NXcanSAS) file.\n'
            'The full GUI setup is embedded so "Load Setup from File…" can\n'
            'later restore every control.'
        )
        self.store_btn.clicked.connect(self._store_results)
        row_out2.addWidget(self.store_btn)

        self.load_setup_btn = QPushButton('Load Setup from File…')
        self.load_setup_btn.setMinimumHeight(26)
        self.load_setup_btn.setStyleSheet('background-color: #ffe082;')
        self.load_setup_btn.setToolTip(
            'Restore every Simple Fits control (model, parameter values,\n'
            'bounds, fit flags, q-range, …) from a NXcanSAS file previously\n'
            'saved by pyirena or by the pyirena-ai agent.'
        )
        self.load_setup_btn.clicked.connect(self._load_setup_from_file)
        row_out2.addWidget(self.load_setup_btn)

        layout.addLayout(row_out2)

        # Row 3: Export | Import parameters
        row_out3 = QHBoxLayout()
        self.export_btn = QPushButton('Save params to JSON')
        self.export_btn.setMinimumHeight(26)
        self.export_btn.setStyleSheet('background-color: lightgreen;')
        self.export_btn.setToolTip(
            'Save current fit parameters to a pyIrena JSON file.\n'
            'Use "Load params from JSON" to restore them later.'
        )
        self.export_btn.clicked.connect(self._export_parameters)
        row_out3.addWidget(self.export_btn)

        self.import_btn = QPushButton('Load params from JSON')
        self.import_btn.setMinimumHeight(26)
        self.import_btn.setStyleSheet('background-color: lightgreen;')
        self.import_btn.setToolTip(
            'Load fit parameters from a previously saved pyIrena JSON file.\n'
            'Use "Save params to JSON" to create a compatible file.'
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

        # Row 5: Passes + Calc. Uncertainty (MC) — last row in the output section,
        # matches the shared template used by Unified Fit / Modeling / Sizes.
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
        self.unc_btn.setToolTip(
            "Estimate parameter uncertainties by repeating the fit on noise-perturbed data.\n"
            "Set 'Passes' to control how many Monte Carlo replicates are used."
        )
        self.unc_btn.clicked.connect(self._calculate_uncertainty)
        btn_row_unc.addWidget(self.unc_btn)
        layout.addLayout(btn_row_unc)

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
        is_calc = bool(entry.get('calculation', False))
        # Calculation models have no least-squares step: the "Fit?" column is
        # hidden for model params (values are taken as entered).  The BG_*
        # rows keep their checkbox — it controls whether the background
        # prefit (buttons + saved-range replay) refits that term (currently
        # BG_P: fit both B and P vs. hold P and fit B only).
        if is_calc and self._fit_col_header_lbl is not None:
            self._fit_col_header_lbl.setVisible(
                self.model.use_complex_bg and entry['complex_bg'])

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
            fit_chk.setVisible((not is_calc) or name.startswith('BG_'))
            grid.addWidget(fit_chk, row, 0, alignment=Qt.AlignmentFlag.AlignCenter)
            self._param_fit_checks[name] = fit_chk

            # Col 1: Name label (friendly symbol for BG_* params)
            name_lbl = QLabel(_BG_DISPLAY.get(name, name) + ':')
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
        model_name = _resolve_model_name(model_name)
        self.model.set_model(model_name)
        self.model.use_complex_bg = self.complex_bg_check.isChecked()
        if self.model.use_complex_bg and MODEL_REGISTRY[model_name]['complex_bg']:
            from pyirena.core.simple_fits import _BG_PARAMS
            for name, default, lo, hi in _BG_PARAMS:
                self.model.params.setdefault(name, default)
                self.model.limits.setdefault(name, (lo, hi))
        self._saved_param_fixed = {}   # reset "Fit?" state for new model
        self._build_param_widgets()
        self._update_bg_prefit_visibility()
        self._update_invariant_ui()
        self.fit_result = None
        self.chi2_label.setText('—')
        self.rchi2_label.setText('—')
        # Clear stale model curve and annotations from previous model
        self.graph_window.clear_fit()
        self.graph_window.clear_result_annotations()
        # Refresh the bottom linearization panel for the new model (or hide it).
        # clear_fit() re-titled the empty lin plot; _refresh_linearization draws
        # the correct transform (using current param widget values) or hides it.
        self._refresh_linearization()

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
        self._update_bg_prefit_visibility()
        self._update_invariant_ui()
        self.graph_window.clear_fit()
        self.graph_window.clear_result_annotations()

    def _on_porod_tail_changed(self, state):
        """Sync the Porod-tail option into the model object."""
        self.model.invariant_porod_tail = bool(state)

    def _on_bg_refit_changed(self, state):
        """Sync the background-refit-replay option into the model object."""
        bp = dict(self.model.bg_prefit or {})
        bp['enabled'] = bool(state)
        self.model.bg_prefit = bp
        self._refresh_bg_prefit_label()

    def _refresh_bg_prefit_label(self):
        """Show which background windows are remembered for prefit replay."""
        if not hasattr(self, '_bg_prefit_ranges_lbl'):
            return
        bp = self.model.bg_prefit or {}
        parts = []
        pl = bp.get('power_law') or {}
        if pl.get('use'):
            p_txt = 'B,P' if pl.get('fit_P', True) else 'B only'
            parts.append(f"power-law ({p_txt}): {pl['q_min']:.4g}–{pl['q_max']:.4g} Å⁻¹")
        fl = bp.get('flat') or {}
        if fl.get('use'):
            parts.append(f"flat: {fl['q_min']:.4g}–{fl['q_max']:.4g} Å⁻¹")
        if parts:
            self._bg_prefit_ranges_lbl.setText('Saved ranges — ' + ' · '.join(parts))
        else:
            self._bg_prefit_ranges_lbl.setText(
                'No saved ranges yet — use the "Fit B/P / Fit Flat btwn '
                'cursors" buttons to record them.')

    def _update_invariant_ui(self):
        """Adjust controls for calculation models (Invariant) vs fit models."""
        is_calc = self.model.is_calculation
        self._invariant_options_row.setVisible(is_calc)
        if is_calc:
            # Refit-replay only makes sense with the complex background on
            has_bg = self.model.use_complex_bg
            self.bg_refit_check.setEnabled(has_bg)
            if not has_bg:
                self.bg_refit_check.setToolTip(
                    'Enable "Complex background" first — the refit replay\n'
                    're-determines B·Q⁻ᴾ + flat from the saved Q ranges.')
            self._refresh_bg_prefit_label()
            self.fit_btn.setText('Calculate Invariant')
            self.fit_btn.setToolTip(
                'Calculate the invariant Q* = ∫q²·(I−bg) dq over the cursor Q range.\n'
                'The invariant is read where the running integral flattens\n'
                '(QmaxUsed).  With contrast and absolute intensities, the\n'
                'volume fraction φ is derived from Q* = 2π²·Δρ²·φ(1−φ).'
            )
            self.graph_model_btn.setToolTip(
                'Display the current background curve (B·Q⁻ᴾ + flat) — the\n'
                'curve subtracted from the data before integration.'
            )
            self.unc_btn.setEnabled(False)
            self.unc_btn.setToolTip(
                'Monte Carlo uncertainty is not applicable to the Invariant\n'
                'calculation (no least-squares fit).'
            )
        else:
            self.fit_btn.setText('Fit')
            self.fit_btn.setToolTip(
                "Fit the selected model to the loaded data using the current Q range.\n"
                "Checked parameters are varied; unchecked parameters are held fixed."
            )
            self.graph_model_btn.setToolTip(
                'Compute and display the current model curve without fitting.')
            self.unc_btn.setEnabled(True)
            self.unc_btn.setToolTip(
                "Estimate parameter uncertainties by repeating the fit on noise-perturbed data.\n"
                "Set 'Passes' to control how many Monte Carlo replicates are used."
            )
        # Residuals do not apply to calculation models — free the space
        self.graph_window.set_residuals_visible(not is_calc)
        # Hide the right (integral) axis when leaving the Invariant model;
        # plot_invariant() re-shows it (with label) on the next Calculate.
        self.graph_window.set_right_axis_visible(False if not is_calc else
                                                 self.graph_window._inv_int_item is not None)

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
                    log.debug("suppressed exception", exc_info=True)

    def _auto_graph_model(self):
        """Sync parameter values and auto-redisplay model curve after edit."""
        self._on_param_edited()
        # Stale annotations are now wrong — clear them before re-graphing
        self.graph_window.clear_result_annotations()
        if self.data is not None:
            self._graph_model()

    # ── Data loading ──────────────────────────────────────────────────────────

    def _on_loader_data_loaded(self, data, hdf5_path: str, display_name: str):
        """Slot wired to DataFileLoaderRow.data_loaded — calls set_data."""
        self.set_data(
            np.asarray(data['Q'],        dtype=float),
            np.asarray(data['Intensity'], dtype=float),
            data.get('Error'),
            label=display_name,
            filepath=hdf5_path,
            is_nxcansas=True,
            slit_length=float(data.get('slit_length', 0.0) or 0.0),
            is_slit_smeared=bool(data.get('is_slit_smeared', False)),
        )

    def _reload_data_with_smearing(self, prefer_slit_smeared):
        """Reload the current file's desmeared or slit-smeared dataset."""
        fp = self.data.get('filepath') if self.data else None
        if not fp:
            return
        try:
            from pyirena.io.hdf5 import readGenericNXcanSAS
            from pathlib import Path as _P
            d = readGenericNXcanSAS(str(_P(fp).parent), _P(fp).name,
                                    prefer_slit_smeared=prefer_slit_smeared)
            if d is None:
                return
            self.set_data(
                np.asarray(d['Q'], dtype=float), np.asarray(d['Intensity'], dtype=float),
                d.get('Error'), label=self.data.get('label', 'Data'),
                filepath=fp, is_nxcansas=True,
                slit_length=float(d.get('slit_length', 0.0) or 0.0),
                is_slit_smeared=bool(d.get('is_slit_smeared', False)),
            )
        except Exception as exc:
            self.status_label.setText(f"Could not reload data: {exc}")

    def _replot_after_slit_change(self):
        if self.data is not None:
            self.graph_window.plot_data(
                self.data['Q'], self.data['Intensity'],
                self.data.get('Error'), label=self.data.get('label', 'Data'))

    def set_data(self, q, intensity, error=None, label='Data',
                 filepath=None, is_nxcansas=False, slit_length=0.0, is_slit_smeared=False):
        """Load SAS data into the panel and plot it."""
        self.data = {
            'Q': np.asarray(q, dtype=float),
            'Intensity': np.asarray(intensity, dtype=float),
            'Error': np.asarray(error, dtype=float) if error is not None else None,
            'label': label,
            'filepath': filepath,
            'is_nxcansas': is_nxcansas,
            'slit_length': slit_length,
            'is_slit_smeared': is_slit_smeared,
        }
        self._refresh_slit_ui_from_data(slit_length, is_slit_smeared, filepath)
        self.fit_result = None
        self.chi2_label.setText('—')
        self.rchi2_label.setText('—')

        if hasattr(self, 'data_loader'):
            self.data_loader.set_filename(label)

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

    # ── Complex-background prefit ──────────────────────────────────────────────

    def _update_bg_prefit_visibility(self):
        """Show the B/P + Flat prefit buttons only when complex bg is active."""
        if not hasattr(self, '_bg_prefit_row'):
            return
        visible = (
            self.model.use_complex_bg
            and MODEL_REGISTRY.get(self.model.model, {}).get('complex_bg', False)
        )
        self._bg_prefit_row.setVisible(visible)

    def _prefit_bg_power_law(self):
        """Prefit the background power-law B·Q⁻ᴾ over the cursor Q range.

        Respects the P "Fit?" checkbox: if P is unchecked, only B is estimated
        at the current (fixed / model-guided) P; otherwise both B and P are fit
        by a log-log least-squares line.  Writes results into the BG_B / BG_P
        widgets as starting values (does NOT run the full fit).
        """
        from pyirena.core.saxs_morph import fit_power_law_bg, fit_power_law_bg_fixed_p

        q, I, _ = self._get_filtered_data()
        if q is None or len(q) < 2:
            self.status_label.setText('Prefit B/P: need ≥2 points between the cursors.')
            return
        q_min, q_max = float(q.min()), float(q.max())

        fit_p_chk = self._param_fit_checks.get('BG_P')
        fit_p = fit_p_chk.isChecked() if fit_p_chk is not None else True

        # Smear the prefit model when the data are slit smeared, so B/P are
        # ideal-space (the main fit smears the background too — F1).
        sl = self.current_slit_length() if self.slit_active() else 0.0

        if fit_p:
            B, P = fit_power_law_bg(q, I, q_min, q_max, slit_length=sl)
            self._set_param_value('BG_P', P)
        else:
            # Hold P at the user's current value; fit only the prefactor B.
            try:
                P = float(self._param_value_edits['BG_P'].text())
            except (KeyError, ValueError):
                P = self.model.params.get('BG_P', 4.0)
            B = fit_power_law_bg_fixed_p(q, I, q_min, q_max, P, slit_length=sl)

        self._set_param_value('BG_B', B)
        # Remember the window so scripted runs can replay this prefit
        # per file (used by the Invariant; harmless for fit models).
        bp = dict(self.model.bg_prefit or {})
        bp['power_law'] = {'use': True, 'q_min': q_min, 'q_max': q_max,
                           'fit_P': bool(fit_p)}
        self.model.bg_prefit = bp
        self._refresh_bg_prefit_label()
        self._collect_model()
        self._auto_graph_model()
        self._refresh_linearization()
        self.status_label.setText(f'Prefit background: B={B:.4g}  P={P:.4g}')

    def _prefit_bg_flat(self):
        """Prefit the flat background term over the cursor Q range.

        Uses the current B and P to subtract the power-law first, then takes the
        median residual as the flat estimate.  Writes BG_flat as a starting value.
        """
        from pyirena.core.saxs_morph import fit_flat_bg

        q, I, _ = self._get_filtered_data()
        if q is None or len(q) < 1:
            self.status_label.setText('Prefit flat: need ≥1 point between the cursors.')
            return
        q_min, q_max = float(q.min()), float(q.max())

        try:
            B = float(self._param_value_edits['BG_B'].text())
        except (KeyError, ValueError):
            B = self.model.params.get('BG_B', 0.0)
        try:
            P = float(self._param_value_edits['BG_P'].text())
        except (KeyError, ValueError):
            P = self.model.params.get('BG_P', 4.0)

        sl = self.current_slit_length() if self.slit_active() else 0.0
        flat = fit_flat_bg(q, I, q_min, q_max, power_law_B=B, power_law_P=P,
                           slit_length=sl)
        self._set_param_value('BG_flat', flat)
        # Remember the window so scripted runs can replay this prefit per file
        bp = dict(self.model.bg_prefit or {})
        bp['flat'] = {'use': True, 'q_min': q_min, 'q_max': q_max}
        self.model.bg_prefit = bp
        self._refresh_bg_prefit_label()
        self._collect_model()
        self._auto_graph_model()
        self._refresh_linearization()
        self.status_label.setText(f'Prefit background: flat={flat:.4g}')

    def _set_param_value(self, name: str, value: float):
        """Write a value into a parameter's edit field and model.params."""
        edit = self._param_value_edits.get(name)
        if edit is not None:
            edit.setText(f'{value:.6g}')
        self.model.params[name] = float(value)

    # ── Collect model from widgets ─────────────────────────────────────────────

    def _collect_model(self) -> SimpleFitModel:
        """Read widget values into self.model and return it."""
        for name, edit in self._param_value_edits.items():
            txt = edit.text().strip()
            if txt:
                try:
                    self.model.params[name] = float(txt)
                except ValueError:
                    log.debug("suppressed exception", exc_info=True)
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
        # Slit smearing (analytic models smeared; Invariant disabled on SMR data).
        self.model.use_slit_smearing = self.slit_active()
        self.model.slit_length = self.current_slit_length()
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
        self._refresh_linearization()
        self.status_label.setText(f'Model: {self.model.model}  (not fitted)')

    # ── Fit ───────────────────────────────────────────────────────────────────

    def _run_fit(self):
        if self.data is None:
            QMessageBox.warning(self, 'No data', 'Load data first.')
            return

        # Push slit-smearing state onto the model before any path (the Invariant
        # calculation path reads it to refuse smeared data with a clear message).
        self.model.use_slit_smearing = self.slit_active()
        self.model.slit_length = self.current_slit_length()

        # Calculation models (Invariant) take a separate path — no least squares
        if self.model.is_calculation:
            self._run_invariant()
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
        from pyirena.gui.fmt_utils import eng_fmt
        chi2 = result.get('chi2', float('nan'))
        rchi2 = result.get('reduced_chi2', float('nan'))
        self.chi2_label.setText(eng_fmt(chi2))
        self.rchi2_label.setText(eng_fmt(rchi2))

        # Write fitted values back to widgets
        self._apply_result_to_widgets(result)

        # Robust fit-quality suffix (uniform across all fit tools)
        from pyirena.gui.quality_display import compute_quality_display
        _q_suffix = ''
        I_model = result.get('I_model')
        if I_model is not None:
            n_free = len(result.get('params', {}) or {})
            _, _, _q_suffix, _ = compute_quality_display(
                q, I, I_model, dI, n_params=max(1, n_free))

        # Update derived quantities if any
        derived = result.get('derived', {})
        if derived:
            derived_txt = '   '.join(f'{k}={eng_fmt(v)}' for k, v in derived.items())
            self.status_label.setText(f'Fit OK | χ²_red={eng_fmt(rchi2, sig=3)} | {derived_txt}{_q_suffix}')
        else:
            self.status_label.setText(f'Fit OK | Reduced χ² = {eng_fmt(rchi2)}{_q_suffix}')

        # Update all plots
        self._update_plots(result)
        self.save_state()

    # ── Invariant calculation ─────────────────────────────────────────────────

    def _run_invariant(self):
        """Run the Invariant calculation over the cursor Q range and plot it."""
        self._get_q_range()
        model = self._collect_model()
        model.invariant_porod_tail = self.porod_tail_check.isChecked()

        # Background prefit replay: refit BG terms over the saved Q windows
        # using the FULL data (windows usually lie outside the integration
        # range), exactly as a scripted run will.  Refit values are written
        # back into the widgets so the user sees what was used.
        prefit_note = ''
        if (model.bg_prefit or {}).get('enabled'):
            # Honor the current BG_P "Fit?" checkbox at replay time
            bp = dict(model.bg_prefit)
            pl = dict(bp.get('power_law') or {})
            if pl:
                chk = self._param_fit_checks.get('BG_P')
                if chk is not None:
                    pl['fit_P'] = chk.isChecked()
                bp['power_law'] = pl
                model.bg_prefit = bp
            applied = model.prefit_background(
                self.data['Q'], self.data['Intensity'])
            for name in ('BG_B', 'BG_P', 'BG_flat'):
                if name in applied:
                    self._set_param_value(name, applied[name])
            if applied.get('warning'):
                prefit_note = f"  ⚠ prefit: {applied['warning']}"
            elif applied:
                prefit_note = '  (background refit from saved ranges)'

        q, I, dI = self._get_filtered_data()
        if q is None or len(q) < 5:
            QMessageBox.warning(self, 'Too few points',
                                'Not enough data points in the selected Q range '
                                '(need ≥ 5).')
            return

        result = model.fit(q, I, dI)   # delegates to the calculation path
        if not result['success']:
            QMessageBox.critical(self, 'Calculation failed',
                                 f"Invariant failed:\n{result.get('error', 'Unknown error')}")
            self.status_label.setText('Invariant calculation FAILED.')
            return

        self.fit_result = result
        self.chi2_label.setText('—')
        self.rchi2_label.setText('—')

        # Background curve (what was subtracted) via the standard model-curve slot
        I_model = result.get('I_model')
        if I_model is not None and np.any(np.asarray(I_model) > 0):
            valid = np.isfinite(I_model) & (I_model > 0)
            if valid.sum() > 1:
                self.graph_window.plot_fit(q[valid], I_model[valid])

        # Corrected data + running integral overlay (Igor-style, right axis)
        extra = result.get('extra_arrays') or {}
        if extra.get('Q_integral') is not None:
            self.graph_window.plot_invariant(
                extra['Q_integral'], extra['running_integral'],
                q=q, I_corrected=extra.get('I_corrected'),
            )

        # Hide the linearization panel (not applicable)
        self._refresh_linearization()

        # Status line with the key numbers
        from pyirena.gui.fmt_utils import eng_fmt
        d = result.get('derived', {})
        inv  = d.get('Invariant', float('nan'))
        phi  = d.get('VolumeFraction', float('nan'))
        qmax = d.get('QmaxUsed', float('nan'))
        msg = (f'Invariant Q* = {eng_fmt(inv)} cm⁻⁴ | '
               f'φ = {phi:.4g} | QmaxUsed = {qmax:.4g} Å⁻¹')
        if self.porod_tail_check.isChecked():
            tail = d.get('PorodTail', 0.0)
            if inv > 0:
                msg += f' | Porod tail = {100.0 * tail / inv:.2g}%'
        warning = result.get('warning', '')
        if warning:
            msg += f'  ⚠ {warning}'
        msg += prefit_note
        self.status_label.setText(msg)

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
            # Snapshot the full GUI state for round-trip restore via
            # "Load Setup from File…".
            try:
                setup_state = self._collect_state()
            except Exception:
                setup_state = None
            save_simple_fit_results(
                filepath=filepath,
                result=self.fit_result,
                model_obj=self.model,
                intensity_data=I,
                intensity_error=dI,
                setup_state=setup_state,
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
            # Rescaled residuals (robust, MAD-based) — uniform across all fit tools.
            from pyirena.core.fit_metrics import rescale_residuals
            r_prime, _ = rescale_residuals(residuals)
            self.graph_window.plot_residuals(q_fit, r_prime)

        # Linearization (bottom panel) — hidden for models with no linearization
        self._refresh_linearization()

    def _refresh_linearization(self):
        """Refresh (or hide) the bottom linearization panel for the current model.

        Called after Fit, after "Graph model", and on model change so the panel
        never gets stuck on a previous model's transform.  Models whose registry
        entry has ``linearization is None`` have the panel hidden entirely.
        """
        has_lin = MODEL_REGISTRY[self.model.model]['linearization'] is not None
        self.graph_window.set_linearization_visible(has_lin)
        if not has_lin:
            self.graph_window.plot_linearization(None)
            return

        if self.data is None:
            self.graph_window.plot_linearization(None)
            return

        # Sync widget values into model.params so the analytic model line matches.
        self._collect_model()

        # Use full data range so out-of-range points are visible in grey alongside
        # the dark in-range points; axes are then bounded to the cursor Q range.
        q_all  = self.data['Q']
        I_all  = self.data['Intensity']
        dI_all = self.data.get('Error')
        q_min, q_max = self.graph_window.get_cursor_range()
        lin = self.model.linearize(q_all, I_all, dI_all, q_min=q_min, q_max=q_max)
        self.graph_window.plot_linearization(lin, q_min=q_min, q_max=q_max)

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

    # ── Export / import parameters ────────────────────────────────────────────

    def _export_parameters(self):
        """Export current parameters to a pyIrena JSON configuration file."""
        import json
        import datetime
        try:
            from pyirena import __version__ as _version
        except Exception:
            _version = 'unknown'

        default_path = str(Path(self._get_data_folder()) / 'pyirena_config.json')
        try:
            _save_opts = QFileDialog.Option.DontConfirmOverwrite | QFileDialog.Option.DontUseNativeDialog
        except AttributeError:
            _save_opts = QFileDialog.DontConfirmOverwrite | QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(
            self, 'Export pyIrena Configuration', default_path,
            'pyIrena Config (*.json);;All Files (*)',
            options=_save_opts,
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
                    self, 'Update Simple Fits Section?',
                    f'This file already has a Simple Fits section. Only that section will be '
                    f'updated — all other tool settings in this file are preserved.\n\n'
                    f'{file_path}',
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.Yes,
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
        # Accept legacy (pre-1.1.0) model-name spellings, e.g. the
        # "Treubner-Strey" typo, so old saved state / result files still
        # restore the correct model instead of silently falling back.
        model_name = _resolve_model_name(state.get('model', 'Guinier'))
        if model_name not in MODEL_NAMES:
            model_name = 'Guinier'

        # Restore model object
        self.model = SimpleFitModel()
        self.model.model = model_name
        self.model.use_complex_bg = bool(state.get('use_complex_bg', False))
        self.model.n_mc_runs = int(state.get('n_mc_runs', 10))
        self.model.invariant_porod_tail = bool(
            state.get('invariant_porod_tail', False))
        self.model.bg_prefit = dict(state.get('bg_prefit') or {})

        saved_params = state.get('params', {})
        saved_limits = state.get('param_limits', {})
        if saved_params:
            self.model.params.update(saved_params)
        if saved_limits:
            self.model.limits.update(
                {k: tuple(v) for k, v in saved_limits.items()}
            )

        # Remove stale params that are no longer in the registry for this model
        # (e.g. Background was removed from Sphere/Spheroid)
        _registry_param_names = {
            name for name, *_ in MODEL_REGISTRY.get(model_name, {}).get('params', [])
        }
        _bg_keys = {'BG_B', 'BG_P', 'BG_flat'}
        for stale_key in list(self.model.params.keys()):
            if stale_key not in _registry_param_names and stale_key not in _bg_keys:
                del self.model.params[stale_key]
        for stale_key in list(self.model.limits.keys()):
            if stale_key not in _registry_param_names and stale_key not in _bg_keys:
                del self.model.limits[stale_key]

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

        self.porod_tail_check.blockSignals(True)
        self.porod_tail_check.setChecked(self.model.invariant_porod_tail)
        self.porod_tail_check.blockSignals(False)

        self.bg_refit_check.blockSignals(True)
        self.bg_refit_check.setChecked(
            bool((self.model.bg_prefit or {}).get('enabled', False)))
        self.bg_refit_check.blockSignals(False)

        # Restore cursor positions if saved
        q_min = state.get('q_min')
        q_max = state.get('q_max')
        if q_min is not None and q_max is not None:
            try:
                self.graph_window.set_cursor_range(float(q_min), float(q_max))
                self.q_min_display.setText(f'{float(q_min):.6g}')
                self.q_max_display.setText(f'{float(q_max):.6g}')
            except Exception:
                log.debug("suppressed exception", exc_info=True)

        self._build_param_widgets()
        self._update_bg_prefit_visibility()
        self._update_invariant_ui()

    def _collect_state(self) -> dict:
        """Return a snapshot of the current panel state.

        Same shape that ``save_state()`` writes to StateManager and that
        ``load_state()`` re-applies — so this dict can be embedded in an
        NXcanSAS file (``setup_state`` kwarg) and later restored verbatim
        via "Load Setup from File…".
        """
        q_min, q_max = self.graph_window.get_cursor_range()
        param_fixed = {
            name: not chk.isChecked()
            for name, chk in self._param_fit_checks.items()
        }
        return {
            'model': self.model.model,
            'q_min': q_min,
            'q_max': q_max,
            'use_complex_bg': self.model.use_complex_bg,
            'no_limits': self.no_limits_check.isChecked(),
            'params': dict(self.model.params),
            'param_limits': {k: list(v) for k, v in self.model.limits.items()},
            'param_fixed': param_fixed,
            'n_mc_runs': self.n_runs_spin.value(),
            'invariant_porod_tail': self.porod_tail_check.isChecked(),
            'bg_prefit': dict(self.model.bg_prefit or {}),
        }

    def save_state(self):
        """Persist current panel state via StateManager."""
        # Refresh cursor Q-range display
        q_min, q_max = self.graph_window.get_cursor_range()
        if q_min is not None:
            self.q_min_display.setText(f'{q_min:.6g}')
        if q_max is not None:
            self.q_max_display.setText(f'{q_max:.6g}')

        state = self._collect_state()
        self.state_manager.update('simple_fits', state)
        self.state_manager.save()

    def _load_setup_from_file(self):
        """Restore the full Simple Fits setup from a NXcanSAS file.

        Reads the ``_pyirena_config`` attribute embedded by
        :func:`pyirena.io.nxcansas_simple_fits.save_simple_fit_results` (or
        by the pyirena-ai agent) and re-applies it through
        :meth:`load_state` so every control matches what was stored.
        """
        from pyirena.gui.setup_loader import prompt_and_load_setup

        if self.data and self.data.get('filepath'):
            default_folder = str(Path(self.data['filepath']).parent)
            suggested = str(self.data['filepath'])
        else:
            default_folder = self._get_data_folder()
            suggested = None

        def _apply(state: dict) -> None:
            # SimpleFitsPanel's load_state reads from StateManager rather than
            # taking a dict.  Push the restored state into the manager first,
            # then trigger a reload so every widget refreshes consistently.
            self.state_manager.update('simple_fits', state)
            self.load_state()

        prompt_and_load_setup(
            parent=self,
            tool="simple_fits",
            default_folder=default_folder,
            apply_state=_apply,
            on_status=lambda msg: self.status_label.setText(msg),
            suggested_path=suggested,
        )

    def closeEvent(self, event):
        self.save_state()
        super().closeEvent(event)
