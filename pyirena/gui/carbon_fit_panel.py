"""
Carbon model GUI panel for pyIrena.

Provides ``CarbonFitGraphWindow`` (full-range I(Q) with the three model
components overlaid, residuals, and an optional linear-Q WAXS zoom) and
``CarbonFitPanel`` (tabbed controls + graph) for interactive full-range
SAXS+WAXS fitting of disordered carbonaceous materials.

The panel is deliberately thin: every control writes straight into a
:class:`~pyirena.core.carbon_fit.CarbonFitModel` attribute and the model does
all the maths.  Because every fittable quantity in that model follows the same
``X`` / ``fit_X`` / ``X_limits`` convention, the whole parameter UI is built by
one :class:`_ParamRow` factory rather than by four hand-written grids — which
is also why adding a parameter to the core needs no change here beyond naming
it in a tab's row list.

See ``docs/carbon_fit_gui.md`` for the user-facing documentation.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

import numpy as np
import pyqtgraph as pg

from pyirena.core.carbon_fit import (
    CARBON_SAXS_MODES,
    CARBON_WAXS_ENVELOPES,
    CarbonFitAborted,
    CarbonFitModel,
    CarbonWaxsPeak,
)
from pyirena.gui._qt import (
    QCheckBox,
    QComboBox,
    QDesktopServices,
    QFileDialog,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QSplitter,
    Qt,
    QTableWidget,
    QTabWidget,
    QUrl,
    QVBoxLayout,
    QWidget,
)
from pyirena.gui.data_loading import DataFileLoaderRow
from pyirena.gui.plot_export import attach_plot_export
from pyirena.gui.q_range_ui import QRangeFields
from pyirena.gui.report_buttons import make_report_buttons
from pyirena.gui.sas_plot import (
    SASPlotStyle,
    get_cursor_q_range,
    make_cursors,
    make_sas_plot,
    plot_iq_data,
    plot_iq_model,
    set_cursor_q_range,
)
from pyirena.gui.sizes_panel import ScrubbableLineEdit
from pyirena.gui.table_utils import (
    attach_table_copy,
    make_numeric_item,
    populating,
    save_rows_as_csv,
)
from pyirena.gui.theme import (
    ACCENT_GREEN,
    SOFT_AMBER,
    SOFT_GREEN,
    accent_button_css,
    readonly_field_css,
    soft_button_css,
)
from pyirena.gui.window_state import install_window_state
from pyirena.state.state_manager import StateManager

log = logging.getLogger(__name__)

_DOC_URL = "https://github.com/jilavsky/pyirena/blob/main/docs/carbon_fit_gui.md"

#: Component curve colours.  Chosen to match the source paper's own figures
#: (Saurel et al. 2019, Figs. 5–7) so a user holding the paper recognises the
#: plot: total fit red, grain Porod blue, micropores orange, diffraction green.
_COMPONENT_PENS = {
    'porod': (31, 119, 180),
    'mp': (255, 127, 14),
    'waxs': (44, 160, 44),
}
_COMPONENT_LABELS = {
    'porod': 'Grain Porod',
    'mp': 'Micropores',
    'waxs': 'Diffraction',
}


# ===========================================================================
# Graph window
# ===========================================================================

class CarbonFitGraphWindow(QWidget):
    """Full-range I(Q) with components, residuals, and an optional WAXS zoom.

    The whole point of the tool is seeing all three Q regions on one curve, so
    the main plot is a single log-log panel spanning the entire measured range
    — commonly five decades — with the total fit and the three components
    overlaid.  That inevitably squeezes the diffraction peaks into the last
    half-decade, which is why there is a second, optional panel showing the
    WAXS region on a linear Q axis: it is off by default and the panel's
    "Show WAXS zoom" box turns it on when the peaks need real attention.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._parent_ref = parent
        self._data_item = None
        self._error_item = None
        self._fit_item = None
        self._component_items: dict = {}
        self._resid_item = None
        self._cursor_a = None
        self._cursor_b = None
        self._zoom_items: list = []
        self._zoom_visible = False
        self._last = {}            # cached curves, for redrawing the zoom panel
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(2)
        self.setLayout(layout)

        self.graphics_layout = pg.GraphicsLayoutWidget()
        self.graphics_layout.setBackground('w')

        self.main_plot = make_sas_plot(
            self.graphics_layout, row=0, col=0,
            x_label='Q  (Å⁻¹)', y_label='I  (cm⁻¹)',
            log_x=True, log_y=True,
            parent_widget=self._parent_ref,
            jpeg_default_name='carbon_fit_IQ',
        )
        self.main_plot.addLegend(offset=(-10, 10), labelTextSize='8pt')

        self.residuals_plot = make_sas_plot(
            self.graphics_layout, row=1, col=0,
            x_label='Q  (Å⁻¹)', y_label='Residuals',
            log_x=True, log_y=False,
            x_link=self.main_plot,
            parent_widget=self._parent_ref,
            jpeg_default_name='carbon_fit_residuals',
        )
        self._resid_zero = pg.InfiniteLine(
            pos=0, angle=0,
            pen=pg.mkPen('k', width=1, style=Qt.PenStyle.DashLine))
        self.residuals_plot.addItem(self._resid_zero)

        # WAXS zoom: linear in Q *and* in I, because a diffraction peak's
        # shape — which is the thing the Voigt widths are fitted to — is
        # unreadable on a log-log plot.
        self.zoom_plot = self.graphics_layout.addPlot(row=2, col=0)
        self.zoom_plot.setLabel('bottom', 'Q  (Å⁻¹)')
        self.zoom_plot.setLabel('left', 'I  (cm⁻¹)')
        self.zoom_plot.showGrid(x=True, y=True, alpha=SASPlotStyle.GRID_ALPHA)
        self.zoom_plot.getAxis('left').enableAutoSIPrefix(False)
        self.zoom_plot.getAxis('bottom').enableAutoSIPrefix(False)
        self.zoom_plot.setTitle('WAXS region (linear)', color='#555555', size='9pt')
        self.zoom_plot.addLegend(offset=(-10, 10), labelTextSize='8pt')
        attach_plot_export(self.zoom_plot, self._parent_ref, 'carbon_fit_waxs',
                           window=self.graphics_layout)

        ci = self.graphics_layout.ci
        ci.layout.setRowStretchFactor(0, 6)
        ci.layout.setRowStretchFactor(1, 1)
        ci.layout.setRowStretchFactor(2, 0)
        self.zoom_plot.hide()

        self.graphics_layout.setSizePolicy(QSizePolicy.Policy.Expanding,
                                           QSizePolicy.Policy.Expanding)
        layout.addWidget(self.graphics_layout, stretch=1)

        self.status_label = QLabel('')
        self.status_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.status_label.setStyleSheet('font-size: 11px; color: #444;')
        layout.addWidget(self.status_label)

    # ── Zoom panel ──────────────────────────────────────────────────────────

    def set_zoom_visible(self, visible: bool):
        """Show or hide the linear-Q WAXS panel."""
        visible = bool(visible)
        if visible == self._zoom_visible:
            return
        self._zoom_visible = visible
        ci = self.graphics_layout.ci
        if visible:
            self.zoom_plot.show()
            ci.layout.setRowStretchFactor(2, 4)
            self._redraw_zoom()
        else:
            self.zoom_plot.hide()
            ci.layout.setRowStretchFactor(2, 0)

    def set_zoom_range(self, q_lo: float, q_hi: float):
        """Set the WAXS panel's Q window (linear units)."""
        if q_hi > q_lo > 0:
            self.zoom_plot.setXRange(float(q_lo), float(q_hi), padding=0.02)
            self._redraw_zoom()

    def _redraw_zoom(self):
        """Repaint the linear WAXS panel from the cached curves."""
        for item in self._zoom_items:
            try:
                self.zoom_plot.removeItem(item)
            except Exception:
                log.debug("suppressed exception", exc_info=True)
        self._zoom_items = []
        if not self._zoom_visible or not self._last:
            return

        q = self._last.get('q')
        if q is None:
            return
        data = self._last.get('I_data')
        if data is not None:
            self._zoom_items.append(self.zoom_plot.plot(
                q, data, pen=None, symbol='o', symbolSize=4,
                symbolPen=None, symbolBrush=(60, 60, 60), name='Data'))
        total = self._last.get('total')
        if total is not None:
            self._zoom_items.append(self.zoom_plot.plot(
                q, total, pen=pg.mkPen((200, 30, 30), width=2), name='Total fit'))
        waxs = self._last.get('waxs')
        if waxs is not None:
            self._zoom_items.append(self.zoom_plot.plot(
                q, waxs, pen=pg.mkPen(_COMPONENT_PENS['waxs'], width=1,
                                      style=Qt.PenStyle.DashLine),
                name=_COMPONENT_LABELS['waxs']))

    # ── Main plotting ───────────────────────────────────────────────────────

    def plot_data(self, q, I, dI=None, label='Data'):
        """Draw the measured curve, creating the cursors on first load."""
        for attr in ('_data_item', '_error_item'):
            item = getattr(self, attr, None)
            if item is not None:
                try:
                    self.main_plot.removeItem(item)
                except Exception:
                    log.debug("suppressed exception", exc_info=True)
            setattr(self, attr, None)

        self._data_item, self._error_item = plot_iq_data(
            self.main_plot, q, I, dI, label=label)
        self._last['q'] = np.asarray(q, dtype=float)
        self._last['I_data'] = np.asarray(I, dtype=float)

        if self._cursor_a is None:
            mask = np.isfinite(q) & (q > 0) & np.isfinite(I) & (I > 0)
            if mask.any():
                lo, hi = float(q[mask].min()), float(q[mask].max())
                self._cursor_a, self._cursor_b = make_cursors(self.main_plot, lo, hi)
                # make_cursors insets by 10 % of the log span, which is right
                # for a single-region fit and wrong here: 10 % of five decades
                # is half a decade, and the highest half-decade is where the
                # (100) reflection lives.  This tool fits the whole range by
                # default, so the cursors start on the data's own edges.
                self.set_cursor_range(lo, hi)
        self._redraw_zoom()

    def plot_model(self, q, components: dict, show_components: bool = True):
        """Draw the total fit plus, optionally, the three component curves.

        Args:
            q: Q grid the components were evaluated on [Å⁻¹].
            components: The dict :meth:`CarbonFitModel.evaluate_components`
                returns.
            show_components: Overlay ``porod``/``mp``/``waxs`` as dashed
                curves.  Seeing which region each part of the data belongs to
                is most of what makes a full-range fit debuggable, so this is
                on by default.
        """
        self.clear_model()
        q = np.asarray(q, dtype=float)
        total = np.asarray(components['total'], dtype=float)
        self._fit_item = plot_iq_model(self.main_plot, q, total, label='Total fit')
        self._last['q'] = q
        self._last['total'] = total
        self._last['waxs'] = np.asarray(components['waxs'], dtype=float)

        if show_components:
            for key in ('porod', 'mp', 'waxs'):
                curve = np.asarray(components[key], dtype=float)
                mask = np.isfinite(curve) & (curve > 0) & np.isfinite(q) & (q > 0)
                if mask.sum() < 2:
                    continue
                item = self.main_plot.plot(
                    q[mask], curve[mask],
                    pen=pg.mkPen(_COMPONENT_PENS[key], width=1,
                                 style=Qt.PenStyle.DashLine),
                    name=_COMPONENT_LABELS[key])
                self._component_items[key] = item
        self._redraw_zoom()

    def plot_residuals(self, q, residuals):
        """Draw weighted residuals against Q."""
        self.residuals_plot.clear()
        self._resid_zero = pg.InfiniteLine(
            pos=0, angle=0,
            pen=pg.mkPen('k', width=1, style=Qt.PenStyle.DashLine))
        self.residuals_plot.addItem(self._resid_zero)
        self._resid_item = None

        q = np.asarray(q, dtype=float)
        residuals = np.asarray(residuals, dtype=float)
        mask = np.isfinite(q) & np.isfinite(residuals) & (q > 0)
        if mask.sum() < 1:
            return
        self._resid_item = self.residuals_plot.plot(
            q[mask], residuals[mask], pen=None, symbol='o',
            symbolSize=SASPlotStyle.RESID_SIZE, symbolPen=None,
            symbolBrush=SASPlotStyle.RESID_BRUSH, name='Residuals')

    def clear_model(self):
        """Remove the fit and component curves, leaving the data and cursors."""
        if self._fit_item is not None:
            try:
                self.main_plot.removeItem(self._fit_item)
            except Exception:
                log.debug("suppressed exception", exc_info=True)
            self._fit_item = None
        for item in self._component_items.values():
            try:
                self.main_plot.removeItem(item)
            except Exception:
                log.debug("suppressed exception", exc_info=True)
        self._component_items = {}

    # ── Cursors ─────────────────────────────────────────────────────────────

    def get_cursor_range(self):
        """(q_min, q_max) in linear units, or (None, None) before first load."""
        return get_cursor_q_range(self._cursor_a, self._cursor_b)

    def set_cursor_range(self, q_min: float, q_max: float):
        """Move the cursors; they stay the single source of truth for the range."""
        self._cursor_a, self._cursor_b = set_cursor_q_range(
            self.main_plot, self._cursor_a, self._cursor_b, q_min, q_max)


# ===========================================================================
# One parameter, one row
# ===========================================================================

class _ParamRow:
    """A ``Fit? | label | value | lo | hi | ±std`` row bound to a model attribute.

    Every fittable quantity in :mod:`pyirena.core.carbon_fit` is declared as
    the triple ``X`` / ``fit_X`` / ``X_limits``, so one widget factory can
    drive all of them.  That is what keeps this panel thin: a new core
    parameter needs its name added to a tab's row list and nothing else — no
    new handler, no new state key (the whole model serialises itself), no new
    line in the fit code.

    The row owns no value of its own.  ``sync()`` pushes the widgets into the
    model, ``refresh()`` pulls the model into the widgets, and the model is
    the only place a number is ever stored.
    """

    def __init__(self, owner, attr: str, label: str, unit: str = '',
                 on_change=None, tooltip: str = ''):
        self.owner = owner
        self.attr = attr
        self.on_change = on_change

        self.fit_check = QCheckBox()
        self.fit_check.setToolTip('Refine this parameter during the fit')
        self.fit_check.setChecked(bool(getattr(owner, f'fit_{attr}', False)))
        self.fit_check.stateChanged.connect(self._changed)

        text = f'{label}  ({unit})' if unit else label
        self.label = QLabel(text)
        self.label.setStyleSheet('font-size: 11px;')
        if tooltip:
            self.label.setToolTip(tooltip)

        self.value_edit = ScrubbableLineEdit(_fmt_value(getattr(owner, attr)))
        self.value_edit.setMaximumWidth(95)
        self._base_tooltip = ((tooltip + '\n\n' if tooltip else '')
                              + 'Scroll over the field to nudge the value.')
        self.value_edit.setToolTip(self._base_tooltip)
        self.value_edit.editingFinished.connect(self._changed)

        lo, hi = getattr(owner, f'{attr}_limits', (-np.inf, np.inf))
        self.lo_edit = QLineEdit(_fmt_value(lo))
        self.hi_edit = QLineEdit(_fmt_value(hi))
        for edit in (self.lo_edit, self.hi_edit):
            edit.setMaximumWidth(70)
            edit.setStyleSheet('font-size: 10px;')
            edit.editingFinished.connect(self._changed)

        self.std_label = QLabel('—')
        self.std_label.setStyleSheet('font-size: 10px; color: #666;')
        self.std_label.setMinimumWidth(60)

    def add_to(self, grid: QGridLayout, row: int):
        """Place the six widgets on one row of ``grid``."""
        grid.addWidget(self.fit_check, row, 0)
        grid.addWidget(self.label, row, 1)
        grid.addWidget(self.value_edit, row, 2)
        grid.addWidget(self.lo_edit, row, 3)
        grid.addWidget(self.hi_edit, row, 4)
        grid.addWidget(self.std_label, row, 5)

    def _changed(self, *_):
        self.sync()
        if self.on_change is not None:
            self.on_change()

    def sync(self) -> None:
        """Widgets → model.  Unparseable text is reverted, not silently dropped."""
        setattr(self.owner, f'fit_{self.attr}', bool(self.fit_check.isChecked()))
        current = float(getattr(self.owner, self.attr))
        setattr(self.owner, self.attr, _parse(self.value_edit.text(), current))
        lo_now, hi_now = getattr(self.owner, f'{self.attr}_limits', (-np.inf, np.inf))
        setattr(self.owner, f'{self.attr}_limits',
                (_parse(self.lo_edit.text(), lo_now),
                 _parse(self.hi_edit.text(), hi_now)))
        self.refresh(keep_std=True)

    def refresh(self, keep_std: bool = False) -> None:
        """Model → widgets, without re-triggering the change handlers."""
        for widget, text in (
            (self.value_edit, _fmt_value(getattr(self.owner, self.attr))),
            (self.lo_edit, _fmt_value(
                getattr(self.owner, f'{self.attr}_limits', (0, 0))[0])),
            (self.hi_edit, _fmt_value(
                getattr(self.owner, f'{self.attr}_limits', (0, 0))[1])),
        ):
            widget.blockSignals(True)
            widget.setText(text)
            widget.blockSignals(False)
        self.fit_check.blockSignals(True)
        self.fit_check.setChecked(bool(getattr(self.owner, f'fit_{self.attr}', False)))
        self.fit_check.blockSignals(False)
        if not keep_std:
            self.std_label.setText('—')

    def show_uncertainty(self, std) -> None:
        """Display a 1-σ uncertainty next to the value."""
        try:
            value = float(std)
        except (TypeError, ValueError):
            value = float('nan')
        self.std_label.setText('—' if not np.isfinite(value)
                               else f'± {_fmt_value(value, 3)}')

    def mark_pinned(self, pinned: bool, side: str = '') -> None:
        """Flag a value that is sitting on its own bound.

        The number in the field is then dictated by the limit rather than by
        the data, and it is the one case where a parameter can look fitted and
        be meaningless, so it gets a visible mark rather than only a line in
        the summary.
        """
        if pinned:
            self.value_edit.setStyleSheet(
                'background-color: #fdebd0; color: #7e5109;')
            self.value_edit.setToolTip(
                f'Pinned at its {side} limit — this value is set by the bound, '
                f'not by the data. Widen the bound or untick Fit?.')
        else:
            self.value_edit.setStyleSheet('')
            self.value_edit.setToolTip(self._base_tooltip)

    def set_enabled(self, enabled: bool) -> None:
        """Grey the row out without removing it, so the layout does not jump."""
        for widget in (self.fit_check, self.label, self.value_edit,
                       self.lo_edit, self.hi_edit, self.std_label):
            widget.setEnabled(bool(enabled))


def _fmt_value(value, sig: int = 6) -> str:
    """Compact, round-trippable text for a float shown in an editable field."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return ''
    if not np.isfinite(v):
        return 'inf' if v > 0 else ('-inf' if v < 0 else 'nan')
    if v == 0:
        return '0'
    if 1e-3 <= abs(v) < 1e6:
        return f'{v:.{sig}g}'
    return f'{v:.{sig - 2}e}'


def _parse(text: str, fallback: float) -> float:
    """Float from a field, falling back to the model's current value."""
    try:
        v = float(str(text).strip())
    except (TypeError, ValueError):
        return float(fallback)
    return v if np.isfinite(v) or abs(v) == np.inf else float(fallback)


def _grid_header(grid: QGridLayout) -> None:
    """Column headings for a parameter grid."""
    for col, text in ((0, 'Fit?'), (1, 'Parameter'), (2, 'Value'),
                      (3, 'lo'), (4, 'hi'), (5, '± std')):
        lbl = QLabel(text)
        lbl.setStyleSheet('font-weight: bold; font-size: 10px;')
        lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        grid.addWidget(lbl, 0, col)


# ===========================================================================
# The panel
# ===========================================================================

class CarbonFitPanel(QWidget):
    """Interactive Carbon model panel — four control tabs beside one graph.

    Tabs follow the model's own sections, which are the three Q regions plus
    the material: **Background** (grain Porod + roughness), **SAXS region**
    (micropores), **WAXS region** (diffraction peaks), **Material**
    (composition → density → contrast), and a read-only **Results** tab with
    every derived quantity.

    Everything is fitted together.  There is deliberately no per-region fit
    button: the contrast that scales the background and the micropore term is
    computed from the WAXS peak positions and the porosity, so fitting the
    regions in sequence converges somewhere the simultaneous fit does not.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Carbon model')
        self.resize(1300, 780)

        self.state_manager = StateManager()
        self.model = CarbonFitModel()
        self.data: dict | None = None
        self.fit_result = None
        self._rows: dict[str, _ParamRow] = {}      # key → row, for uncertainties
        self._peak_rows: list = []                 # one dict per peak widget block
        self._updating = False                     # guard against feedback loops
        self._stop_requested = False
        self._last_progress = 0.0

        self.init_ui()
        self.load_state()
        self._refresh_all()

    # ── UI construction ─────────────────────────────────────────────────────

    def init_ui(self):
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(4, 4, 4, 4)
        self.setLayout(main_layout)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self._create_control_panel())
        self.graph_window = CarbonFitGraphWindow(parent=self)
        splitter.addWidget(self.graph_window)
        splitter.setSizes([470, 830])
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        main_layout.addWidget(splitter, stretch=1)

        self.status_label = QLabel('No data loaded.')
        self.status_label.setStyleSheet('font-size: 11px; color: #555;')
        main_layout.addWidget(self.status_label)

        install_window_state(self, 'carbon_fit_panel', splitters={'main': splitter})

    def _create_control_panel(self) -> QWidget:
        panel = QWidget()
        panel.setMinimumWidth(430)
        panel.setMaximumWidth(560)
        layout = QVBoxLayout()
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(5)
        panel.setLayout(layout)

        # Title + help
        title_row = QHBoxLayout()
        title = QLabel('Carbon model')
        title.setStyleSheet('font-size: 14px; font-weight: bold; color: #2c3e50;')
        title_row.addWidget(title)
        title_row.addStretch()
        help_btn = QPushButton('? Help')
        help_btn.setFixedSize(60, 22)
        help_btn.setStyleSheet(
            'QPushButton{background:#c0392b;color:white;font-size:11px;border-radius:3px;}'
            'QPushButton:hover{background:#e74c3c;}')
        help_btn.setToolTip('Open online documentation in your browser')
        help_btn.clicked.connect(
            lambda: QDesktopServices.openUrl(QUrl(_DOC_URL)))
        title_row.addWidget(help_btn)
        layout.addLayout(title_row)

        self.data_loader = DataFileLoaderRow(state_manager=self.state_manager)
        self.data_loader.data_loaded.connect(self._on_loader_data_loaded)
        layout.addWidget(self.data_loader)

        # Tabs — one per model section
        self.tabs = QTabWidget()
        self.tabs.addTab(_scrolled(self._build_background_tab()), 'Background')
        self.tabs.addTab(_scrolled(self._build_saxs_tab()), 'SAXS region')
        self.tabs.addTab(_scrolled(self._build_waxs_tab()), 'WAXS region')
        self.tabs.addTab(_scrolled(self._build_material_tab()), 'Material')
        self.tabs.addTab(self._build_results_tab(), 'Results')
        layout.addWidget(self.tabs, stretch=1)

        layout.addWidget(self._build_q_range_box())
        layout.addLayout(self._build_action_rows())
        return panel

    # ── Tab: Background ─────────────────────────────────────────────────────

    def _build_background_tab(self) -> QWidget:
        bg = self.model.background
        page = QWidget()
        v = QVBoxLayout()
        v.setContentsMargins(6, 6, 6, 6)
        page.setLayout(v)

        v.addWidget(_note(
            'Scattering from the outer surface of the powder grains — '
            'I = 2π(Δρ)²·[S_macro·Q⁻ⁿ + S_rough·f_rough(Q, R_rough)] '
            '(Saurel et al. eq. 3), plus a flat instrumental background. '
            'Dominates the lowest Q.'))

        self.bg_enabled = _check('Include grain Porod scattering', bg.enabled,
                                 self._on_structure_changed)
        v.addWidget(self.bg_enabled)

        grid = QGridLayout()
        grid.setSpacing(3)
        _grid_header(grid)
        self._add_rows(grid, 1, [
            ('background.S_macro', bg, 'S_macro', 'cm²/cm³',
             'Specific surface area of the grain surface. Only a surface area '
             'while the Porod exponent is exactly 4.'),
            ('background.porod_exponent', bg, 'porod_exponent', '',
             'Porod exponent n. Fixed at 4 (sharp interface) by default; free '
             'it only when the low-Q slope clearly is not −4.'),
        ])
        v.addLayout(grid)

        self.bg_roughness = _check(
            'Add nanoscale surface roughness', bg.use_roughness,
            self._on_structure_changed)
        self.bg_roughness.setToolTip(
            'Adds a second Porod term that switches on above Q ≈ 1/R_rough. '
            'Its surface area adds to S_macro at high Q (eq. 4).')
        v.addWidget(self.bg_roughness)

        self.bg_rough_box = QGroupBox('Surface roughness')
        rg = QGridLayout()
        rg.setSpacing(3)
        _grid_header(rg)
        self._add_rows(rg, 1, [
            ('background.S_rough', bg, 'S_rough', 'cm²/cm³',
             'Extra surface area contributed by the roughness.'),
            ('background.R_rough', bg, 'R_rough', 'Å',
             'Characteristic roughness length. The roughness term is invisible '
             'below Q ≈ 1/R_rough.'),
        ])
        self.bg_rough_box.setLayout(rg)
        v.addWidget(self.bg_rough_box)

        fg = QGridLayout()
        fg.setSpacing(3)
        _grid_header(fg)
        self._add_rows(fg, 1, [
            ('background.flat_background', bg, 'flat_background', 'cm⁻¹',
             'Flat instrumental background added to the whole model.'),
        ])
        v.addLayout(fg)
        v.addStretch()
        return page

    # ── Tab: SAXS region ────────────────────────────────────────────────────

    def _build_saxs_tab(self) -> QWidget:
        s = self.model.saxs
        page = QWidget()
        v = QVBoxLayout()
        v.setContentsMargins(6, 6, 6, 6)
        page.setLayout(v)

        v.addWidget(_note(
            'Micropore scattering, dominant in the middle of the range. The '
            'two models are alternative descriptions of the same pores and '
            'are never summed — pick the one that fits the sample.'))

        self.saxs_enabled = _check('Include micropore scattering', s.enabled,
                                   self._on_structure_changed)
        v.addWidget(self.saxs_enabled)

        mode_row = QHBoxLayout()
        mode_row.addWidget(QLabel('Model:'))
        self.saxs_mode_combo = QComboBox()
        for key, label in CARBON_SAXS_MODES.items():
            self.saxs_mode_combo.addItem(label, key)
        self.saxs_mode_combo.currentIndexChanged.connect(self._on_structure_changed)
        mode_row.addWidget(self.saxs_mode_combo, 1)
        v.addLayout(mode_row)

        # ── fractal + globule branch ──
        self.saxs_fractal_box = QGroupBox('Pores + fractal aggregation')
        fv = QVBoxLayout()
        grid = QGridLayout()
        grid.setSpacing(3)
        _grid_header(grid)
        self._add_rows(grid, 1, [
            ('saxs.phi', s, 'phi', '',
             'Micropore volume fraction φ. Also drives the sample density, '
             'and through it the grain contrast.'),
            ('saxs.pore_radius', s, 'pore_radius', 'Å',
             'Pore radius r. Micropores in carbons are typically 3–10 Å.'),
            ('saxs.globule_k', s, 'globule_k', '',
             'Globule shape factor k: 1 for monodisperse spheres, larger for '
             'pores of less well-defined shape.'),
        ])
        fv.addLayout(grid)

        self.saxs_use_fractal = _check(
            'Pores aggregate into a mass fractal', s.use_fractal,
            self._on_structure_changed)
        self.saxs_use_fractal.setToolTip(
            'Multiplies by the Teixeira structure factor. Leave off for '
            'dilute, uncorrelated pores.')
        fv.addWidget(self.saxs_use_fractal)

        self.saxs_fractal_params = QWidget()
        fg = QGridLayout()
        fg.setContentsMargins(0, 0, 0, 0)
        fg.setSpacing(3)
        _grid_header(fg)
        self._add_rows(fg, 1, [
            ('saxs.fractal_D', s, 'fractal_D', '',
             'Mass fractal dimension D of the aggregate.'),
            ('saxs.fractal_sigma', s, 'fractal_sigma', 'Å',
             'Cutoff length Σ above which the fractal ordering is lost.'),
        ])
        self.saxs_fractal_params.setLayout(fg)
        fv.addWidget(self.saxs_fractal_params)
        self.saxs_fractal_box.setLayout(fv)
        v.addWidget(self.saxs_fractal_box)

        # ── Teubner-Strey branch ──
        self.saxs_ts_box = QGroupBox('Teubner-Strey (two-phase)')
        tg = QGridLayout()
        tg.setSpacing(3)
        _grid_header(tg)
        self._add_rows(tg, 1, [
            ('saxs.ts_I0', s, 'ts_I0', 'cm⁻¹', 'Prefactor I₀ of I = I₀/(1 + C₁Q² + C₂Q⁴).'),
            ('saxs.ts_C1', s, 'ts_C1', 'Å²',
             'C₁. Negative values produce the correlation peak; ξ, d and the '
             'amphiphilicity f_a are derived from C₁ and C₂.'),
            ('saxs.ts_C2', s, 'ts_C2', 'Å⁴', 'C₂; also sets the micropore surface area.'),
        ])
        self.saxs_ts_box.setLayout(tg)
        v.addWidget(self.saxs_ts_box)
        v.addStretch()
        return page

    # ── Tab: WAXS region ────────────────────────────────────────────────────

    def _build_waxs_tab(self) -> QWidget:
        w = self.model.waxs
        page = QWidget()
        v = QVBoxLayout()
        v.setContentsMargins(6, 6, 6, 6)
        page.setLayout(v)

        v.addWidget(_note(
            'Turbostratic stacking: true Voigt peaks (Gaussian = crystallite '
            'size, Lorentzian = layer curvature), a shared Debye-Waller '
            'factor, and the 1/Q² powder-orientation average.'))

        self.waxs_enabled = _check('Include diffraction', w.enabled,
                                   self._on_structure_changed)
        v.addWidget(self.waxs_enabled)

        grid = QGridLayout()
        grid.setSpacing(3)
        _grid_header(grid)
        self._add_rows(grid, 1, [
            ('waxs.delta_z2', w, 'delta_z2', 'Å²',
             'Mean-square interlayer-spacing fluctuation ⟨δz²⟩. Shared by every '
             'peak, because it is a property of the material — that is what '
             'makes it identifiable from the (004)/(002) intensity ratio.'),
        ])
        v.addLayout(grid)

        self.waxs_orientation = _check(
            'Apply 1/Q² orientation average', w.use_orientation_factor,
            self._on_structure_changed)
        self.waxs_orientation.setToolTip(
            'Powder average of a locally one-dimensional stacking correlation '
            '(Annex 3 eq. A3.25). Leave on unless you know why not: a fit '
            'without it looks fine and reports a wrong amplitude.')
        v.addWidget(self.waxs_orientation)

        env_row = QHBoxLayout()
        env_row.addWidget(QLabel('Layer shape:'))
        self.waxs_env_combo = QComboBox()
        for key, label in CARBON_WAXS_ENVELOPES.items():
            self.waxs_env_combo.addItem(label, key)
        self.waxs_env_combo.currentIndexChanged.connect(self._on_structure_changed)
        env_row.addWidget(self.waxs_env_combo, 1)
        v.addLayout(env_row)

        self.waxs_crumple_box = QGroupBox('Crumpled-layer geometry')
        cg_outer = QVBoxLayout()
        cg_outer.addWidget(_note(
            'For layers that bend rather than forming flat nanocrystallites '
            '(eq. 16–17). The same geometry may be what the SAXS region sees, '
            'so each parameter can be linked to its SAXS counterpart instead '
            'of being fitted twice.'))
        cg = QGridLayout()
        cg.setSpacing(3)
        _grid_header(cg)
        self._add_rows(cg, 1, [
            ('waxs.R_layer', w, 'R_layer', 'Å',
             'Radius over which a layer stays flat before it bends.'),
            ('waxs.fractal_D', w, 'fractal_D', '', 'Crumpling fractal dimension D.'),
            ('waxs.fractal_sigma', w, 'fractal_sigma', 'Å', 'Crumpling cutoff Σ.'),
        ])
        cg_outer.addLayout(cg)

        self.link_R = _check('Link R to the SAXS pore radius', w.link_R_to_pore,
                             self._on_structure_changed)
        self.link_D = _check('Link D to the SAXS fractal dimension',
                             w.link_D_to_saxs, self._on_structure_changed)
        self.link_sigma = _check('Link Σ to the SAXS fractal cutoff',
                                 w.link_sigma_to_saxs, self._on_structure_changed)
        for cb in (self.link_R, self.link_D, self.link_sigma):
            cb.setToolTip('A linked parameter takes the SAXS region\'s value '
                          'and leaves the fit vector, so the optimiser never '
                          'sees the same degree of freedom twice.')
            cg_outer.addWidget(cb)
        self.waxs_crumple_box.setLayout(cg_outer)
        v.addWidget(self.waxs_crumple_box)

        # ── Peak list ──
        peak_header = QHBoxLayout()
        peak_header.addWidget(QLabel('Diffraction peaks'))
        peak_header.addStretch()
        add_btn = QPushButton('+ Add peak')
        add_btn.setMaximumWidth(90)
        add_btn.setStyleSheet(soft_button_css(SOFT_GREEN))
        add_btn.setToolTip('Add another reflection to the model')
        add_btn.clicked.connect(self._on_add_peak)
        peak_header.addWidget(add_btn)
        v.addLayout(peak_header)

        self.peaks_container = QWidget()
        self.peaks_layout = QVBoxLayout()
        self.peaks_layout.setContentsMargins(0, 0, 0, 0)
        self.peaks_layout.setSpacing(4)
        self.peaks_container.setLayout(self.peaks_layout)
        v.addWidget(self.peaks_container)
        self._rebuild_peak_widgets()

        v.addStretch()
        return page

    # ── Tab: Material ───────────────────────────────────────────────────────

    def _build_material_tab(self) -> QWidget:
        mat = self.model.material
        page = QWidget()
        v = QVBoxLayout()
        v.setContentsMargins(6, 6, 6, 6)
        page.setLayout(v)

        v.addWidget(_note(
            'The fitted peak positions give the lattice spacings, the spacings '
            'give the structural density, the porosity gives the sample '
            'density, and each density gives an SLD and a contrast. Nothing '
            'is re-typed between tabs — the chain re-runs on every fit '
            'iteration.'))

        form = QGridLayout()
        form.setSpacing(4)
        r = 0
        form.addWidget(QLabel('Chemical formula:'), r, 0)
        self.formula_edit = QLineEdit(mat.formula)
        self.formula_edit.setToolTip(
            'Element composition of the solid, e.g. "C", or "C0.95N0.05" for '
            'an N-doped carbon. Parsed by the Scattering Contrast engine.')
        self.formula_edit.editingFinished.connect(self._on_structure_changed)
        form.addWidget(self.formula_edit, r, 1, 1, 2)

        r += 1
        form.addWidget(QLabel('Structural density:'), r, 0)
        self.rho_mode_combo = QComboBox()
        self.rho_mode_combo.addItem('From fitted peaks', 'from_peaks')
        self.rho_mode_combo.addItem('Manual', 'manual')
        self.rho_mode_combo.currentIndexChanged.connect(self._on_structure_changed)
        form.addWidget(self.rho_mode_combo, r, 1)
        self.rho_manual_edit = ScrubbableLineEdit(_fmt_value(mat.rho_struc_manual))
        self.rho_manual_edit.setMaximumWidth(80)
        self.rho_manual_edit.setToolTip('ρ_struc in g/cm³, used when the mode '
                                        'is Manual or no (002)/(100) pair exists.')
        self.rho_manual_edit.editingFinished.connect(self._on_structure_changed)
        form.addWidget(self.rho_manual_edit, r, 2)

        r += 1
        form.addWidget(QLabel('Porosity for density:'), r, 0)
        self.phi_mode_combo = QComboBox()
        self.phi_mode_combo.addItem('From SAXS region', 'auto')
        self.phi_mode_combo.addItem('Manual', 'manual')
        self.phi_mode_combo.setToolTip(
            'The Teubner-Strey branch derives φ from I₀, which itself depends '
            'on the contrast — taking it live would close a loop around the '
            'fit, so that branch needs a manual estimate here.')
        self.phi_mode_combo.currentIndexChanged.connect(self._on_structure_changed)
        form.addWidget(self.phi_mode_combo, r, 1)
        self.phi_manual_edit = ScrubbableLineEdit(_fmt_value(mat.porosity_manual))
        self.phi_manual_edit.setMaximumWidth(80)
        self.phi_manual_edit.editingFinished.connect(self._on_structure_changed)
        form.addWidget(self.phi_manual_edit, r, 2)

        r += 1
        form.addWidget(QLabel('Contrast:'), r, 0)
        self.contrast_mode_combo = QComboBox()
        self.contrast_mode_combo.addItem('From formula + density', 'auto')
        self.contrast_mode_combo.addItem('Manual', 'manual')
        self.contrast_mode_combo.currentIndexChanged.connect(self._on_structure_changed)
        form.addWidget(self.contrast_mode_combo, r, 1, 1, 2)

        r += 1
        form.addWidget(QLabel('  Grain (Δρ)²  (10²⁰ cm⁻⁴):'), r, 0)
        self.contrast_porod_edit = ScrubbableLineEdit(
            _fmt_value(mat.contrast_porod_manual))
        self.contrast_porod_edit.setMaximumWidth(90)
        self.contrast_porod_edit.editingFinished.connect(self._on_structure_changed)
        form.addWidget(self.contrast_porod_edit, r, 1)

        r += 1
        form.addWidget(QLabel('  Pore (Δρ)²  (10²⁰ cm⁻⁴):'), r, 0)
        self.contrast_mp_edit = ScrubbableLineEdit(
            _fmt_value(mat.contrast_micropore_manual))
        self.contrast_mp_edit.setMaximumWidth(90)
        self.contrast_mp_edit.editingFinished.connect(self._on_structure_changed)
        form.addWidget(self.contrast_mp_edit, r, 1)
        v.addLayout(form)

        ref_box = QGroupBox('Peak labels and graphite reference')
        rg = QGridLayout()
        rg.setSpacing(4)
        self.d002_label_edit = QLineEdit(mat.d002_label)
        self.d100_label_edit = QLineEdit(mat.d100_label)
        self.rho_graphite_edit = ScrubbableLineEdit(_fmt_value(mat.rho_graphite))
        self.d002_graphite_edit = ScrubbableLineEdit(_fmt_value(mat.d002_graphite))
        self.d100_graphite_edit = ScrubbableLineEdit(_fmt_value(mat.d100_graphite))
        self.d100_graphite_edit.setToolTip(
            'The graphite (100) d-spacing, a·√3/2 = 2.1315 Å. The source paper '
            'prints 0.246 nm, which is the a lattice parameter, not the '
            'd-spacing — using it here would scale the density wrongly.')
        for i, (label, widget) in enumerate((
            ('(002) peak label:', self.d002_label_edit),
            ('(100) peak label:', self.d100_label_edit),
            ('ρ graphite  (g/cm³):', self.rho_graphite_edit),
            ('d₀₀₂ graphite  (Å):', self.d002_graphite_edit),
            ('d₁₀₀ graphite  (Å):', self.d100_graphite_edit),
        )):
            widget.setMaximumWidth(90)
            widget.editingFinished.connect(self._on_structure_changed)
            rg.addWidget(QLabel(label), i, 0)
            rg.addWidget(widget, i, 1)
        ref_box.setLayout(rg)
        v.addWidget(ref_box)

        chain_box = QGroupBox('Computed now')
        cv = QVBoxLayout()
        self.material_readout = QLabel('')
        self.material_readout.setStyleSheet(readonly_field_css()
                                            + 'font-family: monospace; font-size: 10px;')
        self.material_readout.setWordWrap(True)
        self.material_readout.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse)
        cv.addWidget(self.material_readout)
        chain_box.setLayout(cv)
        v.addWidget(chain_box)
        v.addStretch()
        return page

    # ── Tab: Results ────────────────────────────────────────────────────────

    def _build_results_tab(self) -> QWidget:
        page = QWidget()
        v = QVBoxLayout()
        v.setContentsMargins(6, 6, 6, 6)
        page.setLayout(v)

        v.addWidget(_note(
            'Everything the fitted parameters imply, in the units a carbon '
            'paper quotes. Updates live as the model changes.'))

        self.derived_table = QTableWidget(0, 3)
        self.derived_table.setHorizontalHeaderLabels(['Quantity', 'Value', 'Units'])
        self.derived_table.verticalHeader().setVisible(False)
        self.derived_table.setAlternatingRowColors(True)
        attach_table_copy(self.derived_table, on_save_csv=self._save_derived_csv)
        v.addWidget(self.derived_table, stretch=1)

        self.fit_summary = QLabel('Not fitted yet.')
        self.fit_summary.setStyleSheet(readonly_field_css() + 'font-size: 10px;')
        self.fit_summary.setWordWrap(True)
        v.addWidget(self.fit_summary)
        return page

    # ── Q range and action buttons ──────────────────────────────────────────

    def _build_q_range_box(self) -> QWidget:
        box = QGroupBox('Q range for fit')
        v = QVBoxLayout()
        v.setContentsMargins(6, 4, 6, 4)
        box.setLayout(v)

        self.q_range_fields = QRangeFields(
            get_range=lambda: self.graph_window.get_cursor_range(),
            set_range=lambda lo, hi: self.graph_window.set_cursor_range(lo, hi),
            get_data_range=self._data_q_range,
        )
        self.q_range_fields.message.connect(self._set_status)
        self.q_range_fields.range_changed.connect(self._on_q_range_typed)
        v.addWidget(self.q_range_fields)

        opt_row = QHBoxLayout()
        self.show_components_check = _check('Show components', True,
                                            self._on_display_option_changed)
        self.show_components_check.setToolTip(
            'Overlay the grain Porod, micropore and diffraction curves. '
            'Seeing which region owns which part of the data is most of what '
            'makes a full-range fit debuggable.')
        opt_row.addWidget(self.show_components_check)

        self.zoom_check = _check('WAXS zoom', False, self._on_display_option_changed)
        self.zoom_check.setToolTip(
            'Add a linear-Q panel below the residuals. The full range is five '
            'decades wide, which crushes the diffraction peaks.')
        opt_row.addWidget(self.zoom_check)

        self.auto_update_check = _check('Auto-update', True, None)
        self.auto_update_check.setToolTip(
            'Redraw the model whenever a control changes. Turn off for very '
            'large datasets.')
        opt_row.addWidget(self.auto_update_check)
        v.addLayout(opt_row)
        return box

    def _build_action_rows(self) -> QVBoxLayout:
        col = QVBoxLayout()
        col.setSpacing(4)

        row1 = QHBoxLayout()
        self.graph_btn = QPushButton('Graph model')
        self.graph_btn.setMinimumHeight(28)
        self.graph_btn.setToolTip('Evaluate the model at the current parameters '
                                  'without fitting.')
        self.graph_btn.clicked.connect(self._on_graph_model)
        row1.addWidget(self.graph_btn)

        self.fit_btn = QPushButton('Fit all')
        self.fit_btn.setMinimumHeight(28)
        self.fit_btn.setStyleSheet(accent_button_css(ACCENT_GREEN, '#1e8449'))
        self.fit_btn.setToolTip(
            'Refine every ticked parameter across the whole Q range at once.\n'
            'The three regions share the contrast, so they cannot be fitted '
            'separately.')
        self.fit_btn.clicked.connect(self._on_fit)
        row1.addWidget(self.fit_btn)

        self.stop_btn = QPushButton('Stop')
        self.stop_btn.setMinimumHeight(28)
        self.stop_btn.setMaximumWidth(60)
        self.stop_btn.setEnabled(False)
        self.stop_btn.setToolTip('Abort the running fit and restore the '
                                 'parameters it started from.')
        self.stop_btn.clicked.connect(self._on_stop)
        row1.addWidget(self.stop_btn)
        col.addLayout(row1)

        row2 = QHBoxLayout()
        row2.addWidget(QLabel('MC passes:'))
        self.mc_spin = QSpinBox()
        self.mc_spin.setRange(2, 500)
        self.mc_spin.setValue(10)
        self.mc_spin.setMaximumWidth(60)
        self.mc_spin.setToolTip('Noise-perturbed refits used for the '
                                'uncertainty estimate.')
        row2.addWidget(self.mc_spin)

        self.mc_btn = QPushButton('Calc. Uncertainty (MC)')
        self.mc_btn.setMinimumHeight(26)
        self.mc_btn.setStyleSheet(accent_button_css('#16a085', '#1abc9c'))
        self.mc_btn.setToolTip(
            'Re-fit noise-perturbed copies of the data. More honest than the '
            'covariance estimate when parameters are correlated — and in this '
            'model they always are, because contrast couples the regions.')
        self.mc_btn.clicked.connect(self._on_monte_carlo)
        row2.addWidget(self.mc_btn)
        col.addLayout(row2)

        row3 = QHBoxLayout()
        self.save_btn = QPushButton('Store in File')
        self.save_btn.setMinimumHeight(26)
        self.save_btn.setStyleSheet(soft_button_css(SOFT_GREEN))
        self.save_btn.setToolTip('Write the fit into entry/carbon_fit_results '
                                 'in the loaded HDF5 file.')
        self.save_btn.clicked.connect(self._on_save_to_file)
        row3.addWidget(self.save_btn)

        self.load_setup_btn = QPushButton('Load Setup from File…')
        self.load_setup_btn.setMinimumHeight(26)
        self.load_setup_btn.setStyleSheet(soft_button_css(SOFT_AMBER))
        self.load_setup_btn.setToolTip('Restore every control from a saved '
                                       'result file.')
        self.load_setup_btn.clicked.connect(self._on_load_setup)
        row3.addWidget(self.load_setup_btn)
        col.addLayout(row3)

        row4 = QHBoxLayout()
        self.export_btn = QPushButton('Save params to JSON')
        self.export_btn.setMinimumHeight(26)
        self.export_btn.setStyleSheet(soft_button_css(SOFT_GREEN))
        self.export_btn.setToolTip('Write a carbon_fit section into a pyIrena '
                                   'config file for batch or scripted runs.')
        self.export_btn.clicked.connect(self._on_export_json)
        row4.addWidget(self.export_btn)

        self.import_btn = QPushButton('Load params from JSON')
        self.import_btn.setMinimumHeight(26)
        self.import_btn.setStyleSheet(soft_button_css(SOFT_GREEN))
        self.import_btn.setToolTip('Read the carbon_fit section back out of a '
                                   'pyIrena config file.')
        self.import_btn.clicked.connect(self._on_import_json)
        row4.addWidget(self.import_btn)
        col.addLayout(row4)

        col.addLayout(make_report_buttons(
            self,
            self.results_for_report,
            tool_key='carbon_fit_results',
            default_stem='carbon_fit',
            file_path_provider=lambda: (self.data or {}).get('filepath', ''),
            data_info_provider=self._data_info_for_report,
            status_setter=self._set_status,
            folder_provider=self._data_folder,
            image_widget_provider=lambda: getattr(
                self.graph_window, 'graphics_layout', None),
        ))

        self.reset_btn = QPushButton('Reset to defaults')
        self.reset_btn.setMinimumHeight(26)
        self.reset_btn.setStyleSheet(accent_button_css('#e67e22', '#f39c12'))
        self.reset_btn.setToolTip('Discard the current model and start from '
                                  'the built-in carbon defaults.')
        self.reset_btn.clicked.connect(self._on_reset)
        col.addWidget(self.reset_btn)
        return col

    # ── Row plumbing ────────────────────────────────────────────────────────

    def _add_rows(self, grid: QGridLayout, start_row: int, specs) -> None:
        """Build and place one :class:`_ParamRow` per spec, remembering each."""
        for offset, (key, owner, attr, unit, tooltip) in enumerate(specs):
            row = _ParamRow(owner, attr, attr.replace('_', ' '), unit,
                            on_change=self._on_param_changed, tooltip=tooltip)
            row.add_to(grid, start_row + offset)
            self._rows[key] = row

    def _rebuild_peak_widgets(self) -> None:
        """Rebuild the per-peak control blocks after the peak list changes.

        One group box per peak rather than a table: a peak carries four
        parameters that each want a Fit? box and bounds, which is exactly what
        :class:`_ParamRow` already draws, and a 24-column table in a 470-pixel
        panel would not be usable.
        """
        while self.peaks_layout.count():
            item = self.peaks_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        for key in [k for k in self._rows if k.startswith('peak.')]:
            del self._rows[key]
        self._peak_rows = []

        for index, peak in enumerate(self.model.peaks):
            self.peaks_layout.addWidget(self._build_peak_box(index, peak))

    def _build_peak_box(self, index: int, peak: CarbonWaxsPeak) -> QWidget:
        box = QGroupBox()
        v = QVBoxLayout()
        v.setContentsMargins(5, 4, 5, 4)
        v.setSpacing(3)

        head = QHBoxLayout()
        enabled = _check('', peak.enabled, None)
        enabled.setToolTip('Include this reflection in the model')
        enabled.stateChanged.connect(
            lambda _s, p=peak: self._on_peak_enabled(p, enabled))
        head.addWidget(enabled)

        head.addWidget(QLabel('Label:'))
        label_edit = QLineEdit(peak.label)
        label_edit.setMaximumWidth(70)
        label_edit.setToolTip(
            'Miller index. The Material tab looks peaks up by this label, so '
            'the (002) and (100) names must match what it expects.')
        label_edit.editingFinished.connect(
            lambda p=peak, e=label_edit: self._on_peak_label(p, e))
        head.addWidget(label_edit)

        self.d_labels = getattr(self, 'd_labels', {})
        d_label = QLabel('')
        d_label.setStyleSheet('font-size: 10px; color: #666;')
        d_label.setToolTip('Bragg spacing d = 2π/Q₀')
        head.addWidget(d_label)
        head.addStretch()

        remove = QPushButton('✕')
        remove.setFixedSize(22, 20)
        remove.setToolTip('Remove this peak')
        remove.clicked.connect(lambda _c, i=index: self._on_remove_peak(i))
        head.addWidget(remove)
        v.addLayout(head)

        grid = QGridLayout()
        grid.setSpacing(3)
        _grid_header(grid)
        specs = [
            ('Q0', 'Å⁻¹', 'Peak centre. d = 2π/Q₀; the Material tab reads the '
                          'spacings from here.'),
            ('K', '', 'Peak amplitude, before the 1/Q² and Debye-Waller factors.'),
            ('FWHM_G', 'Å⁻¹', 'Gaussian component — finite crystallite size. '
                              'The coherence length is 2π·0.9/FWHM_G.'),
            ('FWHM_L', 'Å⁻¹', 'Lorentzian component — layer bending and '
                              'curvature (distortions of the second kind).'),
        ]
        for offset, (attr, unit, tooltip) in enumerate(specs):
            row = _ParamRow(peak, attr, attr, unit,
                            on_change=self._on_param_changed, tooltip=tooltip)
            row.add_to(grid, 1 + offset)
            self._rows[f'peak.{peak.label}.{attr}'] = row
        v.addLayout(grid)
        box.setLayout(v)

        self._peak_rows.append({'peak': peak, 'box': box, 'enabled': enabled,
                                'label_edit': label_edit, 'd_label': d_label})
        return box

    # ── Handlers ────────────────────────────────────────────────────────────

    def _on_param_changed(self) -> None:
        """A value, bound or Fit? box changed — redraw if auto-update is on."""
        if self._updating:
            return
        self._update_material_readout()
        self._update_derived_table()
        self._update_peak_labels()
        if self.auto_update_check.isChecked():
            self._graph_model(quiet=True)

    def _on_structure_changed(self) -> None:
        """A control that changes the model's *shape* changed.

        Mode combos, enable boxes and links do not just alter a number, they
        alter which parameters exist — so the visibility pass has to run before
        anything is redrawn.
        """
        if self._updating:
            return
        self._pull_structure_widgets()
        self._apply_visibility()
        self._on_param_changed()

    def _on_display_option_changed(self) -> None:
        self.graph_window.set_zoom_visible(self.zoom_check.isChecked())
        if self.zoom_check.isChecked():
            self._set_zoom_to_peaks()
        if not self._updating:
            self._graph_model(quiet=True)

    def _pull_structure_widgets(self) -> None:
        """Copy the non-parameter controls into the model."""
        m = self.model
        m.background.enabled = self.bg_enabled.isChecked()
        m.background.use_roughness = self.bg_roughness.isChecked()

        m.saxs.enabled = self.saxs_enabled.isChecked()
        m.saxs.mode = self.saxs_mode_combo.currentData() or 'fractal'
        m.saxs.use_fractal = self.saxs_use_fractal.isChecked()

        m.waxs.enabled = self.waxs_enabled.isChecked()
        m.waxs.use_orientation_factor = self.waxs_orientation.isChecked()
        m.waxs.envelope = self.waxs_env_combo.currentData() or 'none'
        m.waxs.link_R_to_pore = self.link_R.isChecked()
        m.waxs.link_D_to_saxs = self.link_D.isChecked()
        m.waxs.link_sigma_to_saxs = self.link_sigma.isChecked()

        mat = m.material
        mat.formula = self.formula_edit.text().strip() or 'C'
        mat.rho_struc_mode = self.rho_mode_combo.currentData() or 'from_peaks'
        mat.rho_struc_manual = _parse(self.rho_manual_edit.text(),
                                      mat.rho_struc_manual)
        mat.porosity_mode = self.phi_mode_combo.currentData() or 'auto'
        mat.porosity_manual = _parse(self.phi_manual_edit.text(),
                                     mat.porosity_manual)
        mat.contrast_mode = self.contrast_mode_combo.currentData() or 'auto'
        mat.contrast_porod_manual = _parse(self.contrast_porod_edit.text(),
                                           mat.contrast_porod_manual)
        mat.contrast_micropore_manual = _parse(self.contrast_mp_edit.text(),
                                               mat.contrast_micropore_manual)
        mat.d002_label = self.d002_label_edit.text().strip() or '002'
        mat.d100_label = self.d100_label_edit.text().strip() or '100'
        mat.rho_graphite = _parse(self.rho_graphite_edit.text(), mat.rho_graphite)
        mat.d002_graphite = _parse(self.d002_graphite_edit.text(), mat.d002_graphite)
        mat.d100_graphite = _parse(self.d100_graphite_edit.text(), mat.d100_graphite)

    def _apply_visibility(self) -> None:
        """Show only the controls the current model shape actually uses."""
        m = self.model
        fractal = m.saxs.mode == 'fractal'
        self.saxs_fractal_box.setVisible(fractal)
        self.saxs_ts_box.setVisible(not fractal)
        self.saxs_fractal_params.setVisible(fractal and m.saxs.use_fractal)
        self.bg_rough_box.setVisible(m.background.use_roughness)
        self.waxs_crumple_box.setVisible(m.waxs.envelope == 'crumpled')

        for widget, enabled in ((self.bg_rough_box, m.background.enabled),
                                (self.saxs_fractal_box, m.saxs.enabled),
                                (self.saxs_ts_box, m.saxs.enabled),
                                (self.waxs_crumple_box, m.waxs.enabled),
                                (self.peaks_container, m.waxs.enabled)):
            widget.setEnabled(bool(enabled))

        # A linked geometry parameter is not a free parameter; grey it out so
        # the panel says the same thing the fit vector does.
        for key, linked in (('waxs.R_layer', m.waxs.link_R_to_pore),
                            ('waxs.fractal_D', m.waxs.link_D_to_saxs),
                            ('waxs.fractal_sigma', m.waxs.link_sigma_to_saxs)):
            row = self._rows.get(key)
            if row is not None:
                row.set_enabled(not linked)

        manual_contrast = m.material.contrast_mode == 'manual'
        self.contrast_porod_edit.setEnabled(manual_contrast)
        self.contrast_mp_edit.setEnabled(manual_contrast)
        self.rho_manual_edit.setEnabled(m.material.rho_struc_mode == 'manual'
                                        or not m.waxs.enabled)
        self.phi_manual_edit.setEnabled(m.material.porosity_mode == 'manual'
                                        or m.saxs.mode != 'fractal')

    def _on_peak_enabled(self, peak: CarbonWaxsPeak, check: QCheckBox) -> None:
        peak.enabled = bool(check.isChecked())
        self._on_param_changed()

    def _on_peak_label(self, peak: CarbonWaxsPeak, edit: QLineEdit) -> None:
        """Renaming a peak re-keys its rows, since results are keyed by label."""
        new = edit.text().strip()
        if not new or new == peak.label:
            return
        peak.label = new
        self._rebuild_peak_widgets()
        self._on_param_changed()

    def _on_add_peak(self) -> None:
        self.model.add_peak()
        self._rebuild_peak_widgets()
        self._apply_visibility()
        self._on_param_changed()

    def _on_remove_peak(self, index: int) -> None:
        self.model.remove_peak(index)
        self._rebuild_peak_widgets()
        self._apply_visibility()
        self._on_param_changed()

    def _update_peak_labels(self) -> None:
        """Refresh each peak block's live d-spacing readout."""
        for entry in self._peak_rows:
            peak = entry['peak']
            d = 2.0 * np.pi / peak.Q0 if peak.Q0 > 0 else float('nan')
            entry['d_label'].setText(f'd = {d:.4g} Å' if np.isfinite(d) else '')
            box = entry['box']
            box.setTitle(f'Peak ({peak.label})')

    # ── Data ────────────────────────────────────────────────────────────────

    def _on_loader_data_loaded(self, data, hdf5_path: str, display_name: str):
        self.set_data(np.asarray(data['Q'], dtype=float),
                      np.asarray(data['Intensity'], dtype=float),
                      data.get('Error'), label=display_name, filepath=hdf5_path)

    def set_data(self, q, I, dI=None, label='Data', filepath=''):
        """Load a measured curve and draw it, then draw the current model.

        Args:
            q: Scattering vector [Å⁻¹].
            I: Intensity [cm⁻¹].
            dI: 1-σ uncertainties, optional.
            label: Legend label.
            filepath: HDF5 file results will be written back into.
        """
        q = np.asarray(q, dtype=float)
        I = np.asarray(I, dtype=float)
        dI = np.asarray(dI, dtype=float) if dI is not None else None
        self.data = {'Q': q, 'Intensity': I, 'Error': dI,
                     'label': label, 'filepath': filepath}
        self.fit_result = None

        self.graph_window.plot_data(q, I, dI, label=label)
        self.q_range_fields.refresh()
        self._set_zoom_to_peaks()
        self._set_status(f'Loaded {label} — {q.size} points, '
                         f'Q {np.nanmin(q):.4g} to {np.nanmax(q):.4g} Å⁻¹.')
        self._graph_model(quiet=True)

    def _data_q_range(self):
        if not self.data:
            return None
        q = self.data.get('Q')
        if q is None or len(q) == 0:
            return None
        return float(np.nanmin(q)), float(np.nanmax(q))

    def _data_folder(self):
        fp = (self.data or {}).get('filepath', '')
        return str(Path(fp).parent) if fp else None

    def _data_info_for_report(self):
        if not self.data:
            return None
        return {'Q': self.data['Q'], 'I': self.data['Intensity'],
                'I_error': self.data.get('Error')}

    def _on_q_range_typed(self, q_min: float, q_max: float):
        self._graph_model(quiet=True)

    def _set_status(self, text: str):
        label = getattr(self, 'status_label', None)
        if label is not None:
            label.setText(str(text))

    def _set_zoom_to_peaks(self):
        """Point the linear WAXS panel at the enabled peaks."""
        centres = [p.Q0 for p in self.model.peaks if p.enabled and p.Q0 > 0]
        if not centres:
            return
        widths = [max(p.FWHM_G, p.FWHM_L, 0.05)
                  for p in self.model.peaks if p.enabled]
        lo = max(min(centres) - 4.0 * max(widths), 1e-3)
        hi = max(centres) + 4.0 * max(widths)
        self.graph_window.set_zoom_range(lo, hi)

    # ── Model evaluation and fitting ────────────────────────────────────────

    def _current_q(self):
        """Q grid to evaluate on: the loaded data, or a synthetic spread."""
        if self.data is not None:
            return self.data['Q']
        return np.logspace(-4, 0.7, 600)

    def _on_graph_model(self):
        self._graph_model(quiet=False)

    def _graph_model(self, quiet: bool = True) -> None:
        """Evaluate and draw the model without fitting."""
        q = self._current_q()
        try:
            components = self.model.evaluate_components(q)
        except Exception as exc:
            self._set_status(f'Model could not be evaluated: {exc}')
            log.debug('carbon_fit: evaluate failed', exc_info=True)
            return
        self.graph_window.plot_model(
            q, components, show_components=self.show_components_check.isChecked())
        self._update_material_readout()
        self._update_derived_table()
        if not quiet:
            self._set_status('Model graphed at the current parameters.')

    def _sync_q_range_into_model(self) -> None:
        """Copy the cursor range into the model's own q_min/q_max.

        Only when the cursors exist: before any data is loaded there are none,
        and zeroing the range then would throw away a limit that had just been
        restored from a saved setup.
        """
        q_min, q_max = self.graph_window.get_cursor_range()
        if q_min and q_max:
            self.model.q_min = float(q_min)
            self.model.q_max = float(q_max)

    def _on_stop(self):
        self._stop_requested = True
        self._set_status('Stopping…')

    def _fit_progress(self, iteration: int, chi2: float) -> None:
        """Progress callback: keeps the Stop button responsive.

        The fit is fast enough (well under a second for a few hundred points)
        that a worker thread would cost more in complexity than it saves, so
        the event loop is pumped here instead — the same approach WAXS Peak
        Fit takes.

        Throttled on the clock rather than on the evaluation count.  A
        badly-conditioned model can run tens of thousands of evaluations, and
        a repaint every tenth one made the fit several times slower than the
        arithmetic it was reporting on.  Ten updates a second is as much as
        anyone can read.
        """
        now = time.monotonic()
        if now - self._last_progress < 0.1:
            return
        self._last_progress = now
        self._set_status(f'Fitting… {iteration} evaluations, χ² = {chi2:.6g}')
        pg.QtWidgets.QApplication.processEvents()
        if self._stop_requested:
            raise CarbonFitAborted('stopped by the user')

    def _on_fit(self):
        self._run_fit(n_mc_runs=0)

    def _on_monte_carlo(self):
        self._run_fit(n_mc_runs=int(self.mc_spin.value()))

    def _run_fit(self, n_mc_runs: int) -> None:
        if self.data is None:
            QMessageBox.information(self, 'Carbon model', 'Load a data file first.')
            return
        self._sync_q_range_into_model()
        self.model.n_mc_runs = int(n_mc_runs)

        self._stop_requested = False
        self._last_progress = 0.0
        self.stop_btn.setEnabled(True)
        self.fit_btn.setEnabled(False)
        self.mc_btn.setEnabled(False)
        try:
            result = self.model.fit(self.data['Q'], self.data['Intensity'],
                                    self.data.get('Error'),
                                    progress=self._fit_progress)
        except ValueError as exc:
            self._set_status(str(exc))
            QMessageBox.warning(self, 'Carbon model', str(exc))
            return
        except Exception as exc:
            self._set_status(f'Fit failed: {exc}')
            log.error('carbon_fit: fit failed', exc_info=True)
            QMessageBox.critical(self, 'Carbon model', f'Fit failed:\n{exc}')
            return
        finally:
            self.stop_btn.setEnabled(False)
            self.fit_btn.setEnabled(True)
            self.mc_btn.setEnabled(True)
            self.model.n_mc_runs = 0

        self.fit_result = result
        self._refresh_all()
        self.graph_window.plot_residuals(result.q, result.residuals)
        self.graph_window.plot_model(
            result.q, self.model.evaluate_components(result.q),
            show_components=self.show_components_check.isChecked())
        for key, row in self._rows.items():
            if key in result.errors:
                row.show_uncertainty(result.errors[key])
        pinned = dict(self.model.pinned_fitted_parameters())
        for key, row in self._rows.items():
            row.mark_pinned(key in pinned, pinned.get(key, ''))

        verdict = 'converged' if result.success else 'did not converge'
        status = (f'Fit {verdict}: reduced χ² = {result.reduced_chi_squared:.5g} '
                  f'over {result.n_points} points, {result.n_params} free '
                  f'parameters.')
        if result.warnings:
            # A pinned parameter looks exactly like one that is not wired up —
            # it uses the whole evaluation budget and never moves — so it has
            # to be said on the status line, not left in the Results tab.
            status += '  ⚠ ' + result.warnings[0].split(' — ')[0]
        self._set_status(status)
        self._update_fit_summary()

    # ── Readouts ────────────────────────────────────────────────────────────

    def _refresh_all(self) -> None:
        """Push the whole model back into the widgets (after a fit or a load)."""
        self._updating = True
        try:
            self._refresh_structure_widgets()
            for row in self._rows.values():
                row.refresh(keep_std=True)
            self._update_peak_labels()
            self._apply_visibility()
        finally:
            self._updating = False
        self._update_material_readout()
        self._update_derived_table()

    def _refresh_structure_widgets(self) -> None:
        """Model → the non-parameter controls, without firing handlers."""
        m = self.model
        pairs = [
            (self.bg_enabled, m.background.enabled),
            (self.bg_roughness, m.background.use_roughness),
            (self.saxs_enabled, m.saxs.enabled),
            (self.saxs_use_fractal, m.saxs.use_fractal),
            (self.waxs_enabled, m.waxs.enabled),
            (self.waxs_orientation, m.waxs.use_orientation_factor),
            (self.link_R, m.waxs.link_R_to_pore),
            (self.link_D, m.waxs.link_D_to_saxs),
            (self.link_sigma, m.waxs.link_sigma_to_saxs),
        ]
        for widget, value in pairs:
            widget.blockSignals(True)
            widget.setChecked(bool(value))
            widget.blockSignals(False)

        for combo, value in ((self.saxs_mode_combo, m.saxs.mode),
                             (self.waxs_env_combo, m.waxs.envelope),
                             (self.rho_mode_combo, m.material.rho_struc_mode),
                             (self.phi_mode_combo, m.material.porosity_mode),
                             (self.contrast_mode_combo, m.material.contrast_mode)):
            index = combo.findData(value)
            combo.blockSignals(True)
            if index >= 0:
                combo.setCurrentIndex(index)
            combo.blockSignals(False)

        mat = m.material
        for widget, text in (
            (self.formula_edit, mat.formula),
            (self.rho_manual_edit, _fmt_value(mat.rho_struc_manual)),
            (self.phi_manual_edit, _fmt_value(mat.porosity_manual)),
            (self.contrast_porod_edit, _fmt_value(mat.contrast_porod_manual)),
            (self.contrast_mp_edit, _fmt_value(mat.contrast_micropore_manual)),
            (self.d002_label_edit, mat.d002_label),
            (self.d100_label_edit, mat.d100_label),
            (self.rho_graphite_edit, _fmt_value(mat.rho_graphite)),
            (self.d002_graphite_edit, _fmt_value(mat.d002_graphite)),
            (self.d100_graphite_edit, _fmt_value(mat.d100_graphite)),
        ):
            widget.blockSignals(True)
            widget.setText(str(text))
            widget.blockSignals(False)

    def _update_material_readout(self) -> None:
        """Show the live composition → density → SLD → contrast chain."""
        try:
            mat = self.model.resolve_material()
        except Exception:
            self.material_readout.setText('Contrast chain unavailable.')
            return
        lines = [
            f"d002        {mat['d002']:>10.4g} Å      "
            f"d100  {mat['d100']:>10.4g} Å",
            f"ρ_struc     {mat['rho_struc']:>10.4g} g/cm³   "
            f"φ     {mat['porosity']:>10.4g}",
            f"ρ_sample    {mat['rho_sample']:>10.4g} g/cm³",
            f"SLD matrix  {mat['sld_struc']:>10.4g} ×10¹⁰ cm⁻²",
            f"SLD grain   {mat['sld_sample']:>10.4g} ×10¹⁰ cm⁻²",
            f"(Δρ)² grain {mat['contrast_porod']:>10.4g} ×10²⁰ cm⁻⁴",
            f"(Δρ)² pore  {mat['contrast_micropore']:>10.4g} ×10²⁰ cm⁻⁴",
        ]
        self.material_readout.setText('\n'.join(lines))

    #: Derived keys in display order, with their labels and units.  Keys the
    #: current mode does not define come back as NaN and are skipped, so the
    #: table shows what the model actually says rather than a wall of N/A.
    _DERIVED_DISPLAY = (
        ('S_part_m2_g', 'Particle surface area (BET-comparable)', 'm²/g'),
        ('S_macro_m2_g', '  of which macroscopic', 'm²/g'),
        ('S_rough_m2_g', '  of which roughness', 'm²/g'),
        ('S_part', 'Particle surface area', 'cm²/cm³'),
        ('mp_phi', 'Micropore volume fraction φ', ''),
        ('mp_radius', 'Pore radius r', 'Å'),
        ('mp_I0', 'Micropore I₀', 'cm⁻¹'),
        ('S_mp_m2_g', 'Micropore surface area', 'm²/g'),
        ('S_mp', 'Micropore surface area', 'cm²/cm³'),
        ('ts_xi', 'Correlation length ξ', 'Å'),
        ('ts_d', 'Repeat distance d', 'Å'),
        ('ts_fa', 'Amphiphilicity f_a', ''),
        ('w_pore', 'Average pore width w_P', 'Å'),
        ('w_carbon', 'Average wall width w_C', 'Å'),
        ('ts_r_spheroid', 'Pore radius (spheroid limit)', 'Å'),
        ('d002', 'Interlayer spacing d₀₀₂', 'Å'),
        ('d100', 'In-plane spacing d₁₀₀', 'Å'),
        ('L_c', 'Stack height L_c', 'Å'),
        ('N_layers', 'Layers per stack', ''),
        ('L_a', 'Layer extent L_a', 'Å'),
        ('delta_z2', 'Stacking disorder ⟨δz²⟩', 'Å²'),
        ('crumple_D', 'Crumpling fractal dimension D', ''),
        ('crumple_sigma', 'Crumpling cutoff Σ', 'Å'),
        ('crumple_R', 'Layer transition radius R', 'Å'),
        ('rho_struc', 'Structural density ρ_struc', 'g/cm³'),
        ('rho_sample', 'Sample density ρ_sample', 'g/cm³'),
        ('sld_struc', 'SLD, carbon matrix', '10¹⁰ cm⁻²'),
        ('sld_sample', 'SLD, grain', '10¹⁰ cm⁻²'),
        ('contrast_porod', 'Contrast, grain vs. vacuum', '10²⁰ cm⁻⁴'),
        ('contrast_micropore', 'Contrast, pore vs. matrix', '10²⁰ cm⁻⁴'),
    )

    def _derived_rows(self):
        """(label, value, units) triples for the Results table and its CSV."""
        try:
            derived = self.model.compute_derived()
        except Exception:
            log.debug('carbon_fit: derived quantities unavailable', exc_info=True)
            return []
        rows = []
        for key, label, unit in self._DERIVED_DISPLAY:
            value = derived.get(key)
            if value is None or not np.isfinite(value):
                continue
            rows.append((label, float(value), unit))
        for peak in self.model.peaks:
            if not peak.enabled:
                continue
            for suffix, label, unit in (
                ('d', 'spacing d', 'Å'),
                ('FWHM', 'total FWHM', 'Å⁻¹'),
                ('L', 'coherence length', 'Å'),
                ('height', 'peak height', 'cm⁻¹'),
            ):
                value = derived.get(f'peak_{peak.label}_{suffix}')
                if value is None or not np.isfinite(value):
                    continue
                rows.append((f'({peak.label}) {label}', float(value), unit))
        return rows

    def _update_derived_table(self) -> None:
        rows = self._derived_rows()
        table = self.derived_table
        with populating(table):
            table.setRowCount(len(rows))
            for r, (label, value, unit) in enumerate(rows):
                table.setItem(r, 0, _text_item(label))
                table.setItem(r, 1, make_numeric_item(value, fmt='{:.5g}'))
                table.setItem(r, 2, _text_item(unit))
        table.resizeColumnsToContents()

    def _save_derived_csv(self) -> None:
        rows = [[label, f'{value:.8g}', unit]
                for label, value, unit in self._derived_rows()]
        folder = self._data_folder() or ''
        save_rows_as_csv(self, ['Quantity', 'Value', 'Units'], rows,
                         default_path=str(Path(folder) / 'carbon_fit_derived.csv')
                         if folder else 'carbon_fit_derived.csv',
                         title='Save derived quantities as CSV')

    def _update_fit_summary(self) -> None:
        result = self.fit_result
        if result is None:
            self.fit_summary.setText('Not fitted yet.')
            return
        quality = result.quality or {}
        bits = [
            f"χ² = {result.chi_squared:.6g}",
            f"reduced χ² = {result.reduced_chi_squared:.6g}",
            f"{result.n_points} points, {result.n_params} free parameters",
            f"{result.n_iterations} model evaluations",
        ]
        verdict = quality.get('verdict') or quality.get('summary')
        if verdict:
            bits.append(str(verdict))
        bits.append(result.message)
        text = ' · '.join(str(b) for b in bits if b)
        for warning in result.warnings:
            text += f'\n⚠ {warning}'
        self.fit_summary.setText(text)

    # ── Persistence ─────────────────────────────────────────────────────────

    def _collect_state(self) -> dict:
        """The panel's full state, ready for StateManager, JSON or HDF5.

        The physics lives under ``model`` as the core object's own
        ``to_dict()``, so a new core field is persisted, scriptable and
        batch-runnable without a line here.  Only genuinely panel-level
        settings are listed alongside it.
        """
        self._pull_structure_widgets()
        for row in self._rows.values():
            row.sync()
        self._sync_q_range_into_model()
        return {
            'schema_version': 1,
            'model': self.model.to_dict(),
            'q_min': self.model.q_min or None,
            'q_max': self.model.q_max or None,
            'auto_update': bool(self.auto_update_check.isChecked()),
            'show_components': bool(self.show_components_check.isChecked()),
            'waxs_zoom_visible': bool(self.zoom_check.isChecked()),
            'active_tab': int(self.tabs.currentIndex()),
            'last_folder': self._data_folder() or '',
        }

    def _apply_state(self, state: dict) -> None:
        """Apply a state dict to the model and the widgets."""
        state = state or {}
        self.model = CarbonFitModel.from_dict(state.get('model') or {})
        if state.get('q_min') is not None:
            self.model.q_min = float(state['q_min'])
        if state.get('q_max') is not None:
            self.model.q_max = float(state['q_max'])

        self._updating = True
        try:
            for widget, key, default in (
                (self.auto_update_check, 'auto_update', True),
                (self.show_components_check, 'show_components', True),
                (self.zoom_check, 'waxs_zoom_visible', False),
            ):
                widget.blockSignals(True)
                widget.setChecked(bool(state.get(key, default)))
                widget.blockSignals(False)
            index = int(state.get('active_tab', 0) or 0)
            if 0 <= index < self.tabs.count():
                self.tabs.setCurrentIndex(index)
            # The peak list may be a different length, so the per-peak widget
            # blocks have to be thrown away and rebuilt against the new model.
            self._rebuild_peak_widgets()
            self._rebind_rows()
        finally:
            self._updating = False

        self.graph_window.set_zoom_visible(self.zoom_check.isChecked())
        self._refresh_all()
        if self.model.q_min and self.model.q_max:
            try:
                self.graph_window.set_cursor_range(self.model.q_min,
                                                   self.model.q_max)
                self.q_range_fields.refresh()
            except Exception:
                log.debug('carbon_fit: could not restore cursors', exc_info=True)
        self._graph_model(quiet=True)

    def _rebind_rows(self) -> None:
        """Point the section rows at the new model object after a state load.

        ``_ParamRow`` holds a reference to the dataclass it edits, and
        ``_apply_state`` replaces the whole model — so every non-peak row has
        to be re-pointed.  The peak rows are rebuilt from scratch instead,
        because the peak list itself may have changed length.
        """
        owners = {'background': self.model.background, 'saxs': self.model.saxs,
                  'waxs': self.model.waxs}
        for key, row in self._rows.items():
            section = key.split('.', 1)[0]
            if section in owners:
                row.owner = owners[section]

    def save_state(self) -> None:
        """Write the panel state into ``StateManager`` and persist it."""
        try:
            self.state_manager.set('carbon_fit', None, self._collect_state())
            self.state_manager.save()
        except Exception:
            log.warning('carbon_fit: could not save panel state', exc_info=True)

    def load_state(self) -> None:
        """Restore the panel state from ``StateManager``."""
        try:
            self._apply_state(self.state_manager.get('carbon_fit') or {})
        except Exception:
            log.warning('carbon_fit: could not restore panel state', exc_info=True)

    # Public aliases, matching the names the setup/control layers quote.
    def get_current_state(self) -> dict:
        """Public alias of :meth:`_collect_state`, used by the api layer."""
        return self._collect_state()

    def apply_state(self, state: dict) -> None:
        """Public alias of :meth:`_apply_state`, used by the api layer."""
        self._apply_state(state)

    def closeEvent(self, event):                     # noqa: N802 (Qt override)
        self.save_state()
        super().closeEvent(event)

    # ── File and config actions ─────────────────────────────────────────────

    def _on_save_to_file(self) -> None:
        if self.fit_result is None:
            QMessageBox.information(self, 'Carbon model',
                                    'Fit the data before storing results.')
            return
        path = (self.data or {}).get('filepath', '')
        if not path:
            QMessageBox.information(self, 'Carbon model',
                                    'No HDF5 file is associated with this data.')
            return
        try:
            from pyirena.io.nxcansas_carbon_fit import save_carbon_fit_results
            save_carbon_fit_results(Path(path), self.fit_result, self.model,
                                    setup_state=self._collect_state())
        except Exception as exc:
            log.error('carbon_fit: save failed', exc_info=True)
            QMessageBox.critical(self, 'Carbon model', f'Could not save:\n{exc}')
            return
        self._set_status(f'Stored in {Path(path).name}:entry/carbon_fit_results')

    def _on_load_setup(self) -> None:
        from pyirena.gui.setup_loader import prompt_and_load_setup

        prompt_and_load_setup(
            self, 'carbon_fit', self._data_folder() or '',
            apply_state=self._apply_state, on_status=self._set_status,
            suggested_path=(self.data or {}).get('filepath', '') or None)

    def _on_export_json(self) -> None:
        """Append a ``carbon_fit`` section to a shared pyIrena config file."""
        import json
        from datetime import datetime

        folder = self._data_folder() or ''
        path, _ = QFileDialog.getSaveFileName(
            self, 'Save parameters to JSON',
            str(Path(folder) / 'pyirena_config.json') if folder
            else 'pyirena_config.json',
            'JSON files (*.json);;All files (*)')
        if not path:
            return
        config = {}
        if Path(path).exists():
            try:
                config = json.loads(Path(path).read_text())
            except Exception:
                log.debug('carbon_fit: existing config unreadable', exc_info=True)
                config = {}
        config.setdefault('_pyirena_config', {})
        config['_pyirena_config'].update(
            {'tool': 'carbon_fit', 'saved_at': datetime.now().isoformat()})
        config['carbon_fit'] = self._collect_state()
        try:
            Path(path).write_text(json.dumps(config, indent=2))
        except Exception as exc:
            QMessageBox.critical(self, 'Carbon model', f'Could not write:\n{exc}')
            return
        self._set_status(f'Parameters written to {Path(path).name}')

    def _on_import_json(self) -> None:
        import json

        folder = self._data_folder() or ''
        path, _ = QFileDialog.getOpenFileName(
            self, 'Load parameters from JSON', folder,
            'JSON files (*.json);;All files (*)')
        if not path:
            return
        try:
            config = json.loads(Path(path).read_text())
        except Exception as exc:
            QMessageBox.critical(self, 'Carbon model', f'Could not read:\n{exc}')
            return
        section = config.get('carbon_fit')
        if section is None:
            QMessageBox.warning(self, 'Carbon model',
                                f'{Path(path).name} has no carbon_fit section.')
            return
        self._apply_state(section)
        self._set_status(f'Parameters loaded from {Path(path).name}')

    def _on_reset(self) -> None:
        if QMessageBox.question(
                self, 'Carbon model',
                'Discard the current model and start from the defaults?'
        ) != QMessageBox.StandardButton.Yes:
            return
        self.fit_result = None
        self._apply_state({})
        self._set_status('Model reset to defaults.')

    # ── Report ──────────────────────────────────────────────────────────────

    def results_for_report(self):
        """The fit in the dict shape ``load_carbon_fit_results()`` produces.

        The saved key names, not the model's internal ones, so the panel's
        report, the Data Selector's report and ``export_fit_report`` cannot
        drift apart — :mod:`pyirena.core.reporting` renders all three.
        """
        if self.fit_result is None:
            return None
        result = self.fit_result
        return {
            'timestamp': result.timestamp,
            'success': result.success,
            'message': result.message,
            'formula': self.model.material.formula,
            'saxs_mode': self.model.saxs.mode,
            'waxs_envelope': self.model.waxs.envelope,
            'n_peaks': sum(1 for p in self.model.peaks if p.enabled),
            'chi_squared': result.chi_squared,
            'reduced_chi_squared': result.reduced_chi_squared,
            'n_points': result.n_points,
            'n_params': result.n_params,
            'q_min': self.model.q_min or None,
            'q_max': self.model.q_max or None,
            'params': result.params,
            'params_std': result.errors,
            'derived': result.derived,
            'fit_quality': result.quality,
            'warnings': list(result.warnings),
        }


# ===========================================================================
# Small widget helpers
# ===========================================================================

def _check(text: str, checked: bool, handler) -> QCheckBox:
    """A checkbox wired to ``handler`` (which may be None)."""
    box = QCheckBox(text)
    box.setChecked(bool(checked))
    if handler is not None:
        box.stateChanged.connect(lambda _s: handler())
    return box


def _note(text: str) -> QLabel:
    """A small wrapped explanatory paragraph above a group of controls."""
    label = QLabel(text)
    label.setWordWrap(True)
    label.setStyleSheet('font-size: 10px; color: #555; padding: 2px;')
    return label


def _text_item(text: str):
    """A plain, non-editable table cell."""
    from pyirena.gui.table_utils import make_text_item

    item = make_text_item(text)
    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
    return item


def _scrolled(page: QWidget) -> QScrollArea:
    """Wrap a tab page so a long parameter list scrolls instead of clipping."""
    area = QScrollArea()
    area.setWidgetResizable(True)
    area.setFrameShape(QScrollArea.Shape.NoFrame)
    area.setWidget(page)
    return area
