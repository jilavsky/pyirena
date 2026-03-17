"""
Modeling panel GUI for pyIrena.

Provides a combined left-panel (controls) + right-panel (graphs) window for
fitting multi-population size-distribution models to SAS data.

Up to 10 populations (P1 … P10) can be simultaneously active.  In Step 1
every population is a parametric size distribution; future steps add Unified
Fit levels and diffraction peaks.

Entry points
------------
ModelingPanel  — main QWidget (used as a tab in the Data Selector).
main()         — standalone CLI entry: python -m pyirena.gui.modeling_panel <file>
"""

from __future__ import annotations

import os
import sys
from copy import deepcopy
from pathlib import Path
from typing import Optional

import numpy as np

try:
    from PySide6.QtWidgets import (
        QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
        QPushButton, QLabel, QLineEdit, QCheckBox, QSpinBox, QTabWidget,
        QGroupBox, QMessageBox, QSplitter, QFileDialog, QComboBox,
        QScrollArea, QFrame, QSizePolicy, QDoubleSpinBox,
    )
    from PySide6.QtCore import Qt, Signal, QTimer
    from PySide6.QtGui import QFont, QDoubleValidator
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
            QPushButton, QLabel, QLineEdit, QCheckBox, QSpinBox, QTabWidget,
            QGroupBox, QMessageBox, QSplitter, QFileDialog, QComboBox,
            QScrollArea, QFrame, QSizePolicy, QDoubleSpinBox,
        )
        from PyQt6.QtCore import Qt, pyqtSignal as Signal, QTimer
        from PyQt6.QtGui import QFont, QDoubleValidator
    except ImportError:
        from PyQt5.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
            QPushButton, QLabel, QLineEdit, QCheckBox, QSpinBox, QTabWidget,
            QGroupBox, QMessageBox, QSplitter, QFileDialog, QComboBox,
            QScrollArea, QFrame, QSizePolicy, QDoubleSpinBox,
        )
        from PyQt5.QtCore import Qt, pyqtSignal as Signal, QTimer
        from PyQt5.QtGui import QFont, QDoubleValidator

import pyqtgraph as pg

from pyirena.core.distributions import DIST_PARAM_NAMES, DIST_LABELS, DIST_DEFAULTS
from pyirena.core.modeling import (
    ModelingEngine, ModelingConfig, SizeDistPopulation, ModelingResult,
)
from pyirena.gui.sas_plot import (
    make_sas_plot, plot_iq_data, set_robust_y_range,
)
from pyirena.gui.unified_fit import ScrubbableLineEdit, _SafeInfiniteLine
from pyirena.io.nxcansas_modeling import save_modeling_results, load_modeling_results
from pyirena.state import StateManager

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────

N_POPULATIONS = 10

# Colors for each population's curve (P1 … P10)
POP_COLORS = [
    '#2980b9',  # P1  blue
    '#e74c3c',  # P2  red
    '#27ae60',  # P3  green
    '#8e44ad',  # P4  purple
    '#e67e22',  # P5  orange
    '#16a085',  # P6  teal
    '#c0392b',  # P7  dark red
    '#2471a3',  # P8  dark blue
    '#1e8449',  # P9  dark green
    '#6c3483',  # P10 dark purple
]

# Structure-factor display labels and their sf_params keys
SF_LABELS = {
    'none':          ('None (dilute)', []),
    'interferences': ('Interferences', ['eta', 'pack']),
    'hard_sphere':   ('Hard Sphere (P-Y)', ['radius', 'volume_fraction']),
}
SF_PARAM_LABELS = {
    'eta': 'Correlation dist. η [Å]',
    'pack': 'Packing factor',
    'radius': 'HS radius [Å]',
    'volume_fraction': 'Volume fraction',
}

# Form-factor display labels and their extra params
FF_LABELS = {
    'sphere':   ('Sphere', []),
    'spheroid': ('Spheroid', ['aspect_ratio']),
}
FF_PARAM_LABELS = {
    'aspect_ratio': 'Aspect ratio',
}

# ──────────────────────────────────────────────────────────────────────────────
# Small helper widgets
# ──────────────────────────────────────────────────────────────────────────────

def _sep(orientation='h') -> QFrame:
    """Thin horizontal or vertical separator line."""
    f = QFrame()
    if orientation == 'h':
        f.setFrameShape(QFrame.Shape.HLine)
    else:
        f.setFrameShape(QFrame.Shape.VLine)
    f.setFrameShadow(QFrame.Shadow.Sunken)
    f.setStyleSheet('color: #cccccc;')
    return f


def _fmt(v: float) -> str:
    """Format float for display in parameter fields."""
    if v == 0.0:
        return '0'
    if abs(v) < 0.01 or abs(v) >= 1e5:
        return f'{v:.3e}'
    if abs(v) < 10:
        return f'{v:.4g}'
    return f'{v:.4g}'


def _parse(text: str, default: float = 0.0) -> float:
    try:
        return float(text)
    except (ValueError, TypeError):
        return default


# ──────────────────────────────────────────────────────────────────────────────
# Population tab widget
# ──────────────────────────────────────────────────────────────────────────────

class PopulationTab(QWidget):
    """Widget for one population's parameters, shown inside a QTabWidget tab."""

    changed = Signal()   # emitted whenever any parameter changes

    def __init__(self, pop_index: int, parent=None):
        super().__init__(parent)
        self.pop_index = pop_index
        self._building = True   # suppress signals during construction

        self._row_widgets: dict = {}  # {key: (label, value, fit_cb, lo, hi)}
        self._sf_rows: dict = {}
        self._ff_rows: dict = {}

        self._build_ui()
        self._building = False

    # ── UI construction ──────────────────────────────────────────────────────

    def _build_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.setSpacing(4)

        # ── Enable toggle ──────────────────────────────────────────────────
        self.use_cb = QCheckBox(f'Use population P{self.pop_index + 1}')
        self.use_cb.setChecked(self.pop_index == 0)
        self.use_cb.setStyleSheet('font-weight: bold;')
        self.use_cb.stateChanged.connect(self._on_use_changed)
        outer.addWidget(self.use_cb)

        # Container that gets enabled/disabled as a group
        self._content = QWidget()
        self._content.setEnabled(self.pop_index == 0)
        c = QVBoxLayout(self._content)
        c.setContentsMargins(0, 0, 0, 0)
        c.setSpacing(4)

        # ── Distribution type ──────────────────────────────────────────────
        row = QHBoxLayout()
        row.addWidget(QLabel('Distribution:'))
        self.dist_combo = QComboBox()
        for key, label in DIST_LABELS.items():
            self.dist_combo.addItem(label, key)
        self.dist_combo.currentIndexChanged.connect(self._on_dist_changed)
        row.addWidget(self.dist_combo)
        row.addStretch()
        c.addLayout(row)

        # ── Parameter table ────────────────────────────────────────────────
        self._param_group = QGroupBox('Distribution Parameters')
        self._param_grid = QGridLayout(self._param_group)
        self._param_grid.setSpacing(2)
        self._build_param_header(self._param_grid, 0)
        c.addWidget(self._param_group)

        # ── Distribution type: Volume / Number ─────────────────────────────
        dist_type_row = QHBoxLayout()
        self.vol_dist_rb = QCheckBox('Volume distribution (default)')
        self.vol_dist_rb.setChecked(True)
        self.vol_dist_rb.stateChanged.connect(self._emit_changed)
        dist_type_row.addWidget(self.vol_dist_rb)
        self.num_dist_rb = QCheckBox('Number distribution')
        self.num_dist_rb.setChecked(False)
        self.num_dist_rb.stateChanged.connect(self._on_num_dist_changed)
        dist_type_row.addWidget(self.num_dist_rb)
        dist_type_row.addStretch()
        c.addLayout(dist_type_row)

        # ── N bins ─────────────────────────────────────────────────────────
        nbins_row = QHBoxLayout()
        nbins_row.addWidget(QLabel('# Bins:'))
        self.nbins_spin = QSpinBox()
        self.nbins_spin.setRange(20, 2000)
        self.nbins_spin.setValue(200)
        self.nbins_spin.setFixedWidth(70)
        self.nbins_spin.valueChanged.connect(self._emit_changed)
        nbins_row.addWidget(self.nbins_spin)
        nbins_row.addStretch()
        c.addLayout(nbins_row)

        c.addWidget(_sep())

        # ── Form factor ────────────────────────────────────────────────────
        ff_group = QGroupBox('Form Factor')
        ff_layout = QVBoxLayout(ff_group)
        ff_layout.setSpacing(2)

        ff_row = QHBoxLayout()
        ff_row.addWidget(QLabel('Shape:'))
        self.ff_combo = QComboBox()
        for key, (label, _) in FF_LABELS.items():
            self.ff_combo.addItem(label, key)
        self.ff_combo.currentIndexChanged.connect(self._on_ff_changed)
        ff_row.addWidget(self.ff_combo)
        ff_row.addStretch()
        ff_layout.addLayout(ff_row)

        self._ff_grid = QGridLayout()
        self._ff_grid.setSpacing(2)
        ff_layout.addLayout(self._ff_grid)
        c.addWidget(ff_group)

        # ── Structure factor ───────────────────────────────────────────────
        sf_group = QGroupBox('Structure Factor')
        sf_layout = QVBoxLayout(sf_group)
        sf_layout.setSpacing(2)

        sf_row = QHBoxLayout()
        sf_row.addWidget(QLabel('Type:'))
        self.sf_combo = QComboBox()
        for key, (label, _) in SF_LABELS.items():
            self.sf_combo.addItem(label, key)
        self.sf_combo.currentIndexChanged.connect(self._on_sf_changed)
        sf_row.addWidget(self.sf_combo)
        sf_row.addStretch()
        sf_layout.addLayout(sf_row)

        self._sf_grid = QGridLayout()
        self._sf_grid.setSpacing(2)
        sf_layout.addLayout(self._sf_grid)
        c.addWidget(sf_group)

        c.addWidget(_sep())

        # ── Contrast / Scale / Volume fraction ─────────────────────────────
        phys_group = QGroupBox('Physical Parameters')
        phys_lay = QGridLayout(phys_group)
        phys_lay.setSpacing(3)

        phys_lay.addWidget(QLabel('Contrast (Δρ)² [10²⁰ cm⁻⁴]:'), 0, 0)
        self.contrast_edit = ScrubbableLineEdit()
        self.contrast_edit.setText('1.0')
        self.contrast_edit.setFixedWidth(90)
        self.contrast_edit.editingFinished.connect(self._emit_changed)
        phys_lay.addWidget(self.contrast_edit, 0, 1)
        self.contrast_fit_cb = QCheckBox('Fit')
        self.contrast_fit_cb.stateChanged.connect(self._emit_changed)
        phys_lay.addWidget(self.contrast_fit_cb, 0, 2)

        phys_lay.addWidget(QLabel('Scale [= Vf(1−Vf)]:'), 1, 0)
        self.scale_edit = ScrubbableLineEdit()
        self.scale_edit.setText('0.001')
        self.scale_edit.setFixedWidth(90)
        self.scale_edit.editingFinished.connect(self._on_scale_changed)
        phys_lay.addWidget(self.scale_edit, 1, 1)
        self.scale_fit_cb = QCheckBox('Fit')
        self.scale_fit_cb.setChecked(True)
        self.scale_fit_cb.stateChanged.connect(self._emit_changed)
        phys_lay.addWidget(self.scale_fit_cb, 1, 2)

        phys_lay.addWidget(QLabel('Volume fraction (Vf):'), 2, 0)
        self.vf_label = QLabel('0.001')
        self.vf_label.setStyleSheet(
            'color: #555; background: #f5f5f5; padding: 2px 4px; border-radius: 2px;'
        )
        phys_lay.addWidget(self.vf_label, 2, 1)

        c.addWidget(phys_group)
        c.addStretch()

        outer.addWidget(self._content)

        # ── Build initial parameter rows ───────────────────────────────────
        self._rebuild_dist_params()
        self._rebuild_ff_params()
        self._rebuild_sf_params()

    def _build_param_header(self, grid: QGridLayout, start_row: int):
        hdr_style = 'font-weight: bold; font-size: 10px;'
        for col, text in enumerate(['Parameter', 'Value', 'Fit?', 'Min', 'Max']):
            lbl = QLabel(text)
            lbl.setStyleSheet(hdr_style)
            grid.addWidget(lbl, start_row, col)

    def _add_param_row(self, grid: QGridLayout, row: int, key: str, label: str,
                       value: float, fit: bool, lo: float, hi: float,
                       store: dict) -> dict:
        """Add a parameter row to *grid* and register widgets in *store*."""
        lbl = QLabel(label)
        val_edit = ScrubbableLineEdit(_fmt(value))
        val_edit.setFixedWidth(90)
        val_edit.editingFinished.connect(self._emit_changed)

        fit_cb = QCheckBox()
        fit_cb.setChecked(fit)
        fit_cb.stateChanged.connect(self._emit_changed)

        lo_edit = ScrubbableLineEdit(_fmt(lo))
        lo_edit.setFixedWidth(80)
        lo_edit.editingFinished.connect(self._emit_changed)

        hi_edit = ScrubbableLineEdit(_fmt(hi))
        hi_edit.setFixedWidth(80)
        hi_edit.editingFinished.connect(self._emit_changed)

        grid.addWidget(lbl,      row, 0)
        grid.addWidget(val_edit, row, 1)
        grid.addWidget(fit_cb,   row, 2)
        grid.addWidget(lo_edit,  row, 3)
        grid.addWidget(hi_edit,  row, 4)

        store[key] = (lbl, val_edit, fit_cb, lo_edit, hi_edit)
        return store

    # ── Dynamic parameter rebuilding ─────────────────────────────────────────

    def _clear_grid(self, grid: QGridLayout):
        while grid.count():
            item = grid.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

    def _rebuild_dist_params(self):
        """Rebuild distribution parameter rows for the current dist_type."""
        self._clear_grid(self._param_grid)
        self._row_widgets.clear()
        dist_type = self.dist_combo.currentData() or 'lognormal'
        defaults = DIST_DEFAULTS.get(dist_type, {})
        param_names = DIST_PARAM_NAMES.get(dist_type, [])

        self._build_param_header(self._param_grid, 0)
        for row_i, pname in enumerate(param_names, start=1):
            dval = defaults.get(pname, 1.0)
            lo, hi = self._default_dist_limits(pname, dval)
            self._add_param_row(
                self._param_grid, row_i, pname, pname.replace('_', ' ').title(),
                dval, False, lo, hi, self._row_widgets,
            )

    def _default_dist_limits(self, pname: str, val: float) -> tuple:
        default_limits = {
            'mean_size':  (1.0, 1e6),
            'width':      (0.01, 1e5),
            'min_size':   (0.0, 1e6),
            'sdeviation': (0.01, 5.0),
            'location':   (1.0, 1e6),
            'parameter':  (2.01, 2.99),
        }
        return default_limits.get(pname, (0.0, 1e10))

    def _rebuild_ff_params(self):
        self._clear_grid(self._ff_grid)
        self._ff_rows.clear()
        ff_key = self.ff_combo.currentData() or 'sphere'
        _, extra_keys = FF_LABELS.get(ff_key, ('', []))
        for row_i, pname in enumerate(extra_keys):
            label = FF_PARAM_LABELS.get(pname, pname)
            self._add_param_row(
                self._ff_grid, row_i, pname, label,
                1.0, False, 0.01, 100.0, self._ff_rows,
            )

    def _rebuild_sf_params(self):
        self._clear_grid(self._sf_grid)
        self._sf_rows.clear()
        sf_key = self.sf_combo.currentData() or 'none'
        _, sf_param_keys = SF_LABELS.get(sf_key, ('', []))
        defaults = {'eta': 100.0, 'pack': 0.1, 'radius': 50.0, 'volume_fraction': 0.1}
        limits = {
            'eta':             (1.0, 1e6),
            'pack':            (0.0, 16.0),
            'radius':          (1.0, 1e6),
            'volume_fraction': (0.0, 0.74),
        }
        for row_i, pname in enumerate(sf_param_keys):
            label = SF_PARAM_LABELS.get(pname, pname)
            dval = defaults.get(pname, 1.0)
            lo, hi = limits.get(pname, (0.0, 1e10))
            self._add_param_row(
                self._sf_grid, row_i, pname, label,
                dval, False, lo, hi, self._sf_rows,
            )

    # ── Signal handlers ──────────────────────────────────────────────────────

    def _on_use_changed(self, state):
        enabled = (state == Qt.CheckState.Checked.value
                   if hasattr(Qt.CheckState, 'Checked')
                   else state == 2)
        self._content.setEnabled(enabled)
        self._emit_changed()

    def _on_dist_changed(self, _=None):
        if self._building:
            return
        self._rebuild_dist_params()
        self._emit_changed()

    def _on_ff_changed(self, _=None):
        if self._building:
            return
        self._rebuild_ff_params()
        self._emit_changed()

    def _on_sf_changed(self, _=None):
        if self._building:
            return
        self._rebuild_sf_params()
        self._emit_changed()

    def _on_num_dist_changed(self, state):
        """Num-dist and vol-dist are mutually exclusive checkboxes."""
        if self._building:
            return
        checked = (state == Qt.CheckState.Checked.value
                   if hasattr(Qt.CheckState, 'Checked')
                   else state == 2)
        self._building = True
        self.vol_dist_rb.setChecked(not checked)
        self._building = False
        self._emit_changed()

    def _on_scale_changed(self):
        scale = _parse(self.scale_edit.text(), 0.001)
        scale = max(0.0, min(1.0, scale))
        # Compute Vf from scale = Vf*(1-Vf) → quadratic
        disc = max(1.0 - 4.0 * scale, 0.0)
        vf = 0.5 * (1.0 - disc ** 0.5)
        self.vf_label.setText(f'{vf:.5f}')
        self._emit_changed()

    def _emit_changed(self, *_):
        if not self._building:
            self.changed.emit()

    # ── No-limits mode ───────────────────────────────────────────────────────

    def set_no_limits(self, no_limits: bool):
        """Show/hide min/max columns based on global 'No limits' checkbox."""
        for store in (self._row_widgets, self._ff_rows, self._sf_rows):
            for key, (lbl, val_edit, fit_cb, lo_edit, hi_edit) in store.items():
                lo_edit.setVisible(not no_limits)
                hi_edit.setVisible(not no_limits)

    # ── Read/write population state ──────────────────────────────────────────

    def to_population(self) -> SizeDistPopulation:
        """Convert GUI state to a SizeDistPopulation dataclass."""
        pop = SizeDistPopulation()
        pop.enabled = self.use_cb.isChecked()
        pop.dist_type = self.dist_combo.currentData() or 'lognormal'

        # Distribution params
        pop.dist_params = {}
        pop.dist_params_fit = {}
        pop.dist_params_limits = {}
        for key, (lbl, val_edit, fit_cb, lo_edit, hi_edit) in self._row_widgets.items():
            pop.dist_params[key] = _parse(val_edit.text(), 1.0)
            pop.dist_params_fit[key] = fit_cb.isChecked()
            lo = _parse(lo_edit.text(), 0.0)
            hi = _parse(hi_edit.text(), 1e10)
            pop.dist_params_limits[key] = (lo, hi)

        # Form factor
        pop.form_factor = self.ff_combo.currentData() or 'sphere'
        pop.ff_params = {}
        pop.ff_params_fit = {}
        pop.ff_params_limits = {}
        for key, (lbl, val_edit, fit_cb, lo_edit, hi_edit) in self._ff_rows.items():
            pop.ff_params[key] = _parse(val_edit.text(), 1.0)
            pop.ff_params_fit[key] = fit_cb.isChecked()
            lo = _parse(lo_edit.text(), 0.0)
            hi = _parse(hi_edit.text(), 1e10)
            pop.ff_params_limits[key] = (lo, hi)

        # Structure factor
        pop.structure_factor = self.sf_combo.currentData() or 'none'
        pop.sf_params = {}
        pop.sf_params_fit = {}
        pop.sf_params_limits = {}
        for key, (lbl, val_edit, fit_cb, lo_edit, hi_edit) in self._sf_rows.items():
            pop.sf_params[key] = _parse(val_edit.text(), 1.0)
            pop.sf_params_fit[key] = fit_cb.isChecked()
            lo = _parse(lo_edit.text(), 0.0)
            hi = _parse(hi_edit.text(), 1e10)
            pop.sf_params_limits[key] = (lo, hi)

        pop.contrast = _parse(self.contrast_edit.text(), 1.0)
        pop.fit_contrast = self.contrast_fit_cb.isChecked()
        pop.scale = _parse(self.scale_edit.text(), 0.001)
        pop.fit_scale = self.scale_fit_cb.isChecked()
        pop.use_number_dist = self.num_dist_rb.isChecked()
        pop.n_bins = self.nbins_spin.value()
        return pop

    def from_population(self, pop: SizeDistPopulation):
        """Load a SizeDistPopulation into the GUI widgets."""
        self._building = True
        try:
            self.use_cb.setChecked(pop.enabled)
            self._content.setEnabled(pop.enabled)

            # Distribution type
            for i in range(self.dist_combo.count()):
                if self.dist_combo.itemData(i) == pop.dist_type:
                    self.dist_combo.setCurrentIndex(i)
                    break
            self._rebuild_dist_params()

            for key, (lbl, val_edit, fit_cb, lo_edit, hi_edit) in self._row_widgets.items():
                val_edit.setText(_fmt(pop.dist_params.get(key, 1.0)))
                fit_cb.setChecked(pop.dist_params_fit.get(key, False))
                lim = pop.dist_params_limits.get(key, (0.0, 1e10))
                lo_edit.setText(_fmt(lim[0]))
                hi_edit.setText(_fmt(lim[1]))

            # Form factor
            for i in range(self.ff_combo.count()):
                if self.ff_combo.itemData(i) == pop.form_factor:
                    self.ff_combo.setCurrentIndex(i)
                    break
            self._rebuild_ff_params()
            for key, (lbl, val_edit, fit_cb, lo_edit, hi_edit) in self._ff_rows.items():
                val_edit.setText(_fmt(pop.ff_params.get(key, 1.0)))
                fit_cb.setChecked(pop.ff_params_fit.get(key, False))
                lim = pop.ff_params_limits.get(key, (0.0, 1e10))
                lo_edit.setText(_fmt(lim[0]))
                hi_edit.setText(_fmt(lim[1]))

            # Structure factor
            for i in range(self.sf_combo.count()):
                if self.sf_combo.itemData(i) == pop.structure_factor:
                    self.sf_combo.setCurrentIndex(i)
                    break
            self._rebuild_sf_params()
            for key, (lbl, val_edit, fit_cb, lo_edit, hi_edit) in self._sf_rows.items():
                val_edit.setText(_fmt(pop.sf_params.get(key, 1.0)))
                fit_cb.setChecked(pop.sf_params_fit.get(key, False))
                lim = pop.sf_params_limits.get(key, (0.0, 1e10))
                lo_edit.setText(_fmt(lim[0]))
                hi_edit.setText(_fmt(lim[1]))

            self.contrast_edit.setText(_fmt(pop.contrast))
            self.contrast_fit_cb.setChecked(pop.fit_contrast)
            self.scale_edit.setText(_fmt(pop.scale))
            self.scale_fit_cb.setChecked(pop.fit_scale)
            self.vol_dist_rb.setChecked(not pop.use_number_dist)
            self.num_dist_rb.setChecked(pop.use_number_dist)
            self.nbins_spin.setValue(pop.n_bins)
            self._on_scale_changed()
        finally:
            self._building = False


# ──────────────────────────────────────────────────────────────────────────────
# Graph window
# ──────────────────────────────────────────────────────────────────────────────

class ModelingGraphWindow(QWidget):
    """Right-side graph area: I(Q) log-log + Residuals + Size distribution."""

    # Emitted when user drags a Q cursor
    cursor_moved = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.q_data = None
        self.data_folder = str(Path.cwd())

        self._cursor_left_line: Optional[_SafeInfiniteLine] = None
        self._cursor_right_line: Optional[_SafeInfiniteLine] = None
        self._cursor_updating = False

        self._pop_items: dict = {}     # {pop_index: PlotDataItem} for I(Q)
        self._total_item = None        # total model PlotDataItem
        self._dist_items: dict = {}    # {pop_index: PlotDataItem} for size dist

        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.gl = pg.GraphicsLayoutWidget()
        self.gl.setBackground('w')
        layout.addWidget(self.gl)

        # Status bar
        self.status_lbl = QLabel('')
        self.status_lbl.setWordWrap(True)
        self.status_lbl.setMaximumHeight(50)
        self.status_lbl.setStyleSheet(
            'padding: 6px; border: 1px solid #ccc; border-radius: 4px; font-size: 10pt;'
        )
        layout.addWidget(self.status_lbl)

        self._init_plots()

    def _init_plots(self):
        self.gl.clear()

        # ── I(Q) main plot ────────────────────────────────────────────────
        self.iq_plot = self.gl.addPlot(row=0, col=0)
        self.iq_plot.setLabel('bottom', 'Q (Å⁻¹)', **{'color': 'k', 'font-size': '11pt'})
        self.iq_plot.setLabel('left', 'Intensity (cm⁻¹)', **{'color': 'k', 'font-size': '11pt'})
        self.iq_plot.setLogMode(x=True, y=True)
        self.iq_plot.showGrid(x=True, y=True, alpha=0.3)
        self.iq_plot.setTitle('Modeling — I(Q)', size='12pt', color='k')
        self.iq_plot.addLegend(offset=(-10, 10), labelTextSize='9pt')
        self._style_axes(self.iq_plot)

        # ── Residuals plot ────────────────────────────────────────────────
        self.resid_plot = self.gl.addPlot(row=1, col=0)
        self.resid_plot.setLabel('bottom', 'Q (Å⁻¹)', **{'color': 'k', 'font-size': '10pt'})
        self.resid_plot.setLabel('left', 'Residuals', **{'color': 'k', 'font-size': '10pt'})
        self.resid_plot.setLogMode(x=True, y=False)
        self.resid_plot.showGrid(x=True, y=True, alpha=0.3)
        self._style_axes(self.resid_plot)
        self.resid_plot.addLine(y=0, pen=pg.mkPen('k', style=Qt.PenStyle.DashLine))

        # ── Size distribution plot ────────────────────────────────────────
        self.dist_plot = self.gl.addPlot(row=2, col=0)
        self.dist_plot.setLabel('bottom', 'Radius (Å)', **{'color': 'k', 'font-size': '10pt'})
        self.dist_plot.setLabel('left', 'Distribution', **{'color': 'k', 'font-size': '10pt'})
        self.dist_plot.showGrid(x=True, y=True, alpha=0.3)
        self.dist_plot.setTitle('Size Distribution', size='10pt', color='k')
        self.dist_plot.addLegend(offset=(-10, 10), labelTextSize='9pt')
        self._style_axes(self.dist_plot)

        # Height ratios: 5 : 2 : 3
        self.gl.ci.layout.setRowStretchFactor(0, 5)
        self.gl.ci.layout.setRowStretchFactor(1, 2)
        self.gl.ci.layout.setRowStretchFactor(2, 3)

        # Link X axes so Q range is synchronised
        self.resid_plot.setXLink(self.iq_plot)

        # JPEG export right-click
        vb = self.iq_plot.getViewBox()
        vb.menu.addSeparator()
        act = vb.menu.addAction('Save I(Q) graph as JPEG…')
        act.triggered.connect(self._save_iq_jpeg)
        vb2 = self.dist_plot.getViewBox()
        vb2.menu.addSeparator()
        act2 = vb2.menu.addAction('Save size distribution as JPEG…')
        act2.triggered.connect(self._save_dist_jpeg)

    @staticmethod
    def _style_axes(plot):
        for side in ('left', 'bottom', 'top', 'right'):
            ax = plot.getAxis(side)
            ax.setPen('k')
            ax.setTextPen('k')
            ax.enableAutoSIPrefix(False)
        plot.showAxis('top')
        plot.showAxis('right')
        plot.getAxis('top').setStyle(showValues=False)
        plot.getAxis('right').setStyle(showValues=False)

    # ── Cursors ──────────────────────────────────────────────────────────────

    def add_cursors(self, q_min_log: float, q_max_log: float):
        """Add Q-range cursors to the I(Q) plot."""
        if self._cursor_left_line is not None:
            try:
                self.iq_plot.removeItem(self._cursor_left_line)
                self.iq_plot.removeItem(self._cursor_right_line)
            except Exception:
                pass

        self._cursor_left_line = _SafeInfiniteLine(
            pos=q_min_log, angle=90, movable=True,
            pen=pg.mkPen('#e74c3c', width=2),
            label='Qmin',
        )
        self._cursor_right_line = _SafeInfiniteLine(
            pos=q_max_log, angle=90, movable=True,
            pen=pg.mkPen('#2980b9', width=2),
            label='Qmax',
        )
        self._cursor_left_line.sigPositionChanged.connect(self._on_cursor_moved)
        self._cursor_right_line.sigPositionChanged.connect(self._on_cursor_moved)
        self.iq_plot.addItem(self._cursor_left_line)
        self.iq_plot.addItem(self._cursor_right_line)

    def _on_cursor_moved(self, _=None):
        if not self._cursor_updating:
            self._cursor_updating = True
            self.cursor_moved.emit()
            self._cursor_updating = False

    def get_q_range(self) -> tuple[float, float]:
        """Return current cursor Q range (linear Å⁻¹)."""
        if self._cursor_left_line is None:
            return (0.001, 1.0)
        lo = 10 ** self._cursor_left_line.getPos()[0]
        hi = 10 ** self._cursor_right_line.getPos()[0]
        return (min(lo, hi), max(lo, hi))

    def set_q_range(self, q_min: float, q_max: float):
        """Set cursor positions from linear Q values."""
        if self._cursor_left_line is not None and q_min > 0 and q_max > 0:
            self._cursor_left_line.setPos(np.log10(q_min))
            self._cursor_right_line.setPos(np.log10(q_max))

    # ── Data plotting ─────────────────────────────────────────────────────────

    def plot_data(self, q, I, dI=None):
        """Plot experimental SAS data."""
        self.q_data = q
        plot_iq_data(self.iq_plot, q, I, dI, label='Data')

    def plot_model(self, result: ModelingResult):
        """Plot total model + per-population curves + distributions + residuals."""
        # Remove old model items
        self._clear_model_items()

        q = result.model_q
        I_total = result.model_I

        # ── Total model ────────────────────────────────────────────────────
        self._total_item = self.iq_plot.plot(
            q, I_total,
            pen=pg.mkPen('#000000', width=2),
            name='Total model',
        )

        # ── Per-population curves ──────────────────────────────────────────
        for k, pi in enumerate(result.pop_indices):
            color = POP_COLORS[pi % len(POP_COLORS)]
            I_pop = result.pop_model_I[k]
            self._pop_items[pi] = self.iq_plot.plot(
                q, I_pop,
                pen=pg.mkPen(color, width=1.5, style=Qt.PenStyle.DashLine),
                name=f'P{pi+1}',
            )

        # ── Residuals ──────────────────────────────────────────────────────
        if self.q_data is not None and len(self.q_data) == len(q):
            # attempt residuals from cached data — may be out of sync
            pass  # residuals are computed in ModelingPanel.run_fit and passed explicitly

        # ── Size distributions ─────────────────────────────────────────────
        self.dist_plot.clear()
        for k, pi in enumerate(result.pop_indices):
            color = POP_COLORS[pi % len(POP_COLORS)]
            rg = result.radius_grids[k]
            vd = result.volume_dists[k]
            nd = result.number_dists[k]
            self._dist_items[pi] = self.dist_plot.plot(
                rg, vd,
                pen=pg.mkPen(color, width=2),
                name=f'P{pi+1} vol',
            )
            self.dist_plot.plot(
                rg, nd,
                pen=pg.mkPen(color, width=1.5, style=Qt.PenStyle.DotLine),
                name=f'P{pi+1} num',
            )

    def plot_residuals(self, q: np.ndarray, residuals: np.ndarray):
        """Plot normalised residuals on the residuals panel."""
        self.resid_plot.clear()
        self.resid_plot.addLine(y=0, pen=pg.mkPen('k', style=Qt.PenStyle.DashLine))
        self.resid_plot.plot(
            q, residuals,
            pen=None, symbol='o', symbolSize=4,
            symbolBrush='#2980b9', symbolPen=None,
        )

    def _clear_model_items(self):
        """Remove all previously plotted model curves."""
        if self._total_item is not None:
            try:
                self.iq_plot.removeItem(self._total_item)
            except Exception:
                pass
            self._total_item = None

        for pi, item in list(self._pop_items.items()):
            try:
                self.iq_plot.removeItem(item)
            except Exception:
                pass
        self._pop_items.clear()

    def set_status(self, msg: str, style: str = 'info'):
        colours = {
            'info':    ('', '#333333'),
            'success': ('#d4edda', '#155724'),
            'error':   ('#f8d7da', '#721c24'),
        }
        bg, fg = colours.get(style, ('', '#333333'))
        bg_css = f'background-color: {bg};' if bg else ''
        self.status_lbl.setText(msg)
        self.status_lbl.setStyleSheet(
            f'padding: 6px; border: 1px solid #ccc; border-radius: 4px; '
            f'font-size: 10pt; color: {fg}; {bg_css}'
        )

    # ── JPEG exports ──────────────────────────────────────────────────────────

    def _save_jpeg(self, plot_item, default_name: str):
        from pyqtgraph.exporters import ImageExporter
        path, _ = QFileDialog.getSaveFileName(
            self, 'Save as JPEG',
            str(Path(self.data_folder) / default_name),
            'JPEG (*.jpg *.jpeg);;All Files (*)',
        )
        if not path:
            return
        try:
            exp = ImageExporter(plot_item)
            exp.parameters()['width'] = 1600
            exp.export(path)
        except Exception as e:
            QMessageBox.warning(self, 'Export failed', str(e))

    def _save_iq_jpeg(self):
        self._save_jpeg(self.iq_plot, 'modeling_iq.jpg')

    def _save_dist_jpeg(self):
        self._save_jpeg(self.dist_plot, 'modeling_dist.jpg')


# ──────────────────────────────────────────────────────────────────────────────
# Main panel
# ──────────────────────────────────────────────────────────────────────────────

class ModelingPanel(QWidget):
    """Main Modeling tool panel: left controls + right graph in a QSplitter."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('pyIrena — Modeling')
        self.setMinimumSize(1300, 900)

        self._data_q: Optional[np.ndarray] = None
        self._data_I: Optional[np.ndarray] = None
        self._data_dI: Optional[np.ndarray] = None
        self._file_path: Optional[Path] = None
        self._last_result: Optional[ModelingResult] = None

        self._engine = ModelingEngine()
        self._state = StateManager()
        self._mc_worker = None   # future: QThread for MC

        self._build_ui()
        self._load_state()

    # ── UI construction ──────────────────────────────────────────────────────

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self._build_left_panel())
        self.graph = ModelingGraphWindow()
        splitter.addWidget(self.graph)
        splitter.setSizes([440, 860])
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 2)
        layout.addWidget(splitter)

        self.graph.cursor_moved.connect(self._on_cursor_moved)

    def _build_left_panel(self) -> QWidget:
        panel = QWidget()
        panel.setFixedWidth(440)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        inner = QWidget()
        lay = QVBoxLayout(inner)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(6)

        # ── Title ────────────────────────────────────────────────────────
        title = QLabel('Modeling Input')
        title.setStyleSheet('font-size: 14px; font-weight: bold;')
        lay.addWidget(title)

        # ── Data file ────────────────────────────────────────────────────
        lay.addWidget(_sep())
        file_row = QHBoxLayout()
        file_row.addWidget(QLabel('Data file:'))
        self.file_edit = QLineEdit()
        self.file_edit.setReadOnly(True)
        self.file_edit.setPlaceholderText('(no file selected)')
        file_row.addWidget(self.file_edit)
        btn_open = QPushButton('Open…')
        btn_open.setFixedWidth(60)
        btn_open.clicked.connect(self._open_file)
        file_row.addWidget(btn_open)
        lay.addLayout(file_row)

        # ── Q range display ──────────────────────────────────────────────
        q_row = QHBoxLayout()
        q_row.addWidget(QLabel('Q fit range:'))
        self.qmin_lbl = QLabel('—')
        self.qmax_lbl = QLabel('—')
        q_row.addWidget(QLabel('min'))
        q_row.addWidget(self.qmin_lbl)
        q_row.addWidget(QLabel('max'))
        q_row.addWidget(self.qmax_lbl)
        q_row.addStretch()
        lay.addLayout(q_row)

        # ── Background ───────────────────────────────────────────────────
        bg_row = QHBoxLayout()
        bg_row.addWidget(QLabel('Background:'))
        self.bg_edit = ScrubbableLineEdit('0.0')
        self.bg_edit.setFixedWidth(90)
        bg_row.addWidget(self.bg_edit)
        self.bg_fit_cb = QCheckBox('Fit')
        self.bg_fit_cb.setChecked(True)
        bg_row.addWidget(self.bg_fit_cb)
        bg_row.addStretch()
        lay.addLayout(bg_row)

        # ── No limits checkbox ───────────────────────────────────────────
        self.no_limits_cb = QCheckBox('No limits? (unconstrained fit)')
        self.no_limits_cb.stateChanged.connect(self._on_no_limits_changed)
        lay.addWidget(self.no_limits_cb)

        lay.addWidget(_sep())

        # ── Population tabs ──────────────────────────────────────────────
        self.pop_tabs = QTabWidget()
        self.pop_tabs.setTabPosition(QTabWidget.TabPosition.North)
        self._pop_widgets: list[PopulationTab] = []

        for i in range(N_POPULATIONS):
            pw = PopulationTab(i)
            pw.changed.connect(self._on_pop_changed)

            scroll_i = QScrollArea()
            scroll_i.setWidgetResizable(True)
            scroll_i.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            scroll_i.setWidget(pw)

            self.pop_tabs.addTab(scroll_i, f'P{i+1}')
            self._pop_widgets.append(pw)

        lay.addWidget(self.pop_tabs)
        lay.addWidget(_sep())

        # ── Action buttons ───────────────────────────────────────────────
        btn_row1 = QHBoxLayout()
        self.btn_graph = QPushButton('Graph Model')
        self.btn_graph.setMinimumHeight(34)
        self.btn_graph.setStyleSheet(
            'QPushButton {background: #52c77a; color: white; font-weight: bold;'
            ' border-radius: 4px;}'
            'QPushButton:hover {background: #42b86a;}'
            'QPushButton:disabled {background: #bdc3c7;}'
        )
        self.btn_graph.clicked.connect(self.graph_model)
        self.btn_graph.setEnabled(False)

        self.btn_fit = QPushButton('Fit')
        self.btn_fit.setMinimumHeight(34)
        self.btn_fit.setStyleSheet(
            'QPushButton {background: #27ae60; color: white; font-weight: bold;'
            ' border-radius: 4px;}'
            'QPushButton:hover {background: #229954;}'
            'QPushButton:disabled {background: #bdc3c7;}'
        )
        self.btn_fit.clicked.connect(self.run_fit)
        self.btn_fit.setEnabled(False)

        btn_row1.addWidget(self.btn_graph)
        btn_row1.addWidget(self.btn_fit)
        lay.addLayout(btn_row1)

        btn_row2 = QHBoxLayout()
        btn_row2.addWidget(QLabel('Passes:'))
        self.n_runs_spin = QSpinBox()
        self.n_runs_spin.setRange(1, 500)
        self.n_runs_spin.setValue(10)
        self.n_runs_spin.setFixedWidth(60)
        btn_row2.addWidget(self.n_runs_spin)

        self.btn_mc = QPushButton('Calc. Uncertainty (MC)')
        self.btn_mc.setMinimumHeight(34)
        self.btn_mc.setStyleSheet(
            'QPushButton {background: #16a085; color: white; font-weight: bold;'
            ' border-radius: 4px;}'
            'QPushButton:hover {background: #1abc9c;}'
            'QPushButton:disabled {background: #bdc3c7;}'
        )
        self.btn_mc.clicked.connect(self.calc_uncertainty_mc)
        self.btn_mc.setEnabled(False)
        btn_row2.addWidget(self.btn_mc)
        lay.addLayout(btn_row2)

        lay.addWidget(_sep())

        # ── Results area ─────────────────────────────────────────────────
        res_group = QGroupBox('Fit Results')
        res_lay = QGridLayout(res_group)
        res_lay.addWidget(QLabel('χ²:'), 0, 0)
        self.chi2_lbl = QLabel('—')
        res_lay.addWidget(self.chi2_lbl, 0, 1)
        res_lay.addWidget(QLabel('χ²/dof:'), 0, 2)
        self.rchi2_lbl = QLabel('—')
        res_lay.addWidget(self.rchi2_lbl, 0, 3)
        res_lay.addWidget(QLabel('DOF:'), 1, 0)
        self.dof_lbl = QLabel('—')
        res_lay.addWidget(self.dof_lbl, 1, 1)
        lay.addWidget(res_group)

        # ── Output buttons ───────────────────────────────────────────────
        out_row = QHBoxLayout()
        self.btn_save = QPushButton('Save to HDF5')
        self.btn_save.setMinimumHeight(30)
        self.btn_save.setStyleSheet(
            'QPushButton {background: #2980b9; color: white; font-weight: bold;'
            ' border-radius: 4px;}'
            'QPushButton:hover {background: #2471a3;}'
            'QPushButton:disabled {background: #bdc3c7;}'
        )
        self.btn_save.clicked.connect(self.save_results)
        self.btn_save.setEnabled(False)

        self.btn_export = QPushButton('Export Parameters (JSON)')
        self.btn_export.setMinimumHeight(30)
        self.btn_export.setStyleSheet(
            'QPushButton {background: #8e44ad; color: white; font-weight: bold;'
            ' border-radius: 4px;}'
            'QPushButton:hover {background: #7d3c98;}'
            'QPushButton:disabled {background: #bdc3c7;}'
        )
        self.btn_export.clicked.connect(self.export_json)
        out_row.addWidget(self.btn_save)
        out_row.addWidget(self.btn_export)
        lay.addLayout(out_row)
        lay.addStretch()

        scroll.setWidget(inner)

        panel_lay = QVBoxLayout(panel)
        panel_lay.setContentsMargins(0, 0, 0, 0)
        panel_lay.addWidget(scroll)
        return panel

    # ── State (load/save) ────────────────────────────────────────────────────

    def _load_state(self):
        mod_state = self._state.get('modeling') or {}
        self.bg_edit.setText(_fmt(mod_state.get('background', 0.0)))
        self.bg_fit_cb.setChecked(mod_state.get('fit_background', True))
        nl = mod_state.get('no_limits', False)
        self.no_limits_cb.setChecked(nl)
        self.n_runs_spin.setValue(mod_state.get('n_mc_runs', 10))

        pops = mod_state.get('populations', [])
        for i, pw in enumerate(self._pop_widgets):
            if i < len(pops):
                pop = _pop_from_dict(pops[i])
                pw.from_population(pop)

    def _save_state(self):
        state = {
            'background':    _parse(self.bg_edit.text(), 0.0),
            'fit_background': self.bg_fit_cb.isChecked(),
            'no_limits':     self.no_limits_cb.isChecked(),
            'n_mc_runs':     self.n_runs_spin.value(),
            'q_min':         None,
            'q_max':         None,
            'populations':   [_pop_to_dict(pw.to_population()) for pw in self._pop_widgets],
        }
        if self._data_q is not None:
            q_lo, q_hi = self.graph.get_q_range()
            state['q_min'] = float(q_lo)
            state['q_max'] = float(q_hi)
        self._state.update('modeling', state)
        self._state.save()

    # ── Signal handlers ──────────────────────────────────────────────────────

    def _on_pop_changed(self):
        pass  # Real-time update is expensive; user presses Graph Model instead

    def _on_no_limits_changed(self, state):
        no_lim = (state == Qt.CheckState.Checked.value
                  if hasattr(Qt.CheckState, 'Checked') else state == 2)
        for pw in self._pop_widgets:
            pw.set_no_limits(no_lim)

    def _on_cursor_moved(self):
        q_lo, q_hi = self.graph.get_q_range()
        self.qmin_lbl.setText(f'{q_lo:.4g}')
        self.qmax_lbl.setText(f'{q_hi:.4g}')

    # ── File I/O ─────────────────────────────────────────────────────────────

    def _open_file(self):
        last = (self._state.get('data_selector') or {}).get('last_folder', '')
        path, _ = QFileDialog.getOpenFileName(
            self, 'Open HDF5 file', last,
            'HDF5 files (*.hdf5 *.h5 *.nxs);;All files (*)',
        )
        if path:
            self.load_file(Path(path))

    def set_data(
        self,
        q: np.ndarray,
        I: np.ndarray,
        dI,
        filename: str,
        filepath: str = '',
        is_nxcansas: bool = True,
    ):
        """Accept pre-loaded Q/I/dI arrays (called from Data Selector)."""
        q = np.asarray(q, dtype=float)
        I = np.asarray(I, dtype=float)
        dI = np.asarray(dI, dtype=float) if dI is not None else I * 0.05

        self._file_path = Path(filepath) if filepath else None
        self._data_q = q
        self._data_I = I
        self._data_dI = dI

        self.file_edit.setText(filename)
        if filepath:
            self.graph.data_folder = str(Path(filepath).parent)

        q_min = float(q.min())
        q_max = float(q.max())
        self.graph.plot_data(q, I, dI)
        set_robust_y_range(self.graph.iq_plot, I)
        self.graph.add_cursors(np.log10(q_min), np.log10(q_max))
        self.graph.set_status(
            f'Loaded {filename}: {len(q)} points, '
            f'Q ∈ [{q_min:.4g}, {q_max:.4g}] Å⁻¹', 'success',
        )

        q_lo, q_hi = self.graph.get_q_range()
        self.qmin_lbl.setText(f'{q_lo:.4g}')
        self.qmax_lbl.setText(f'{q_hi:.4g}')

        self.btn_graph.setEnabled(True)
        self.btn_fit.setEnabled(True)

    def load_file(self, file_path: Path):
        """Load SAS data from an HDF5 NXcanSAS file."""
        try:
            from pyirena.io.hdf5 import load_nxcansas
            q, I, dI = load_nxcansas(file_path)
        except Exception as e:
            QMessageBox.critical(self, 'Load error', str(e))
            return

        self.set_data(
            q, I, dI,
            filename=str(file_path.name),
            filepath=str(file_path),
            is_nxcansas=True,
        )

    # ── Build ModelingConfig from GUI ────────────────────────────────────────

    def _build_config(self) -> ModelingConfig:
        populations = [pw.to_population() for pw in self._pop_widgets]
        q_lo, q_hi = self.graph.get_q_range()
        return ModelingConfig(
            populations=populations,
            background=_parse(self.bg_edit.text(), 0.0),
            fit_background=self.bg_fit_cb.isChecked(),
            background_limits=(0.0, 1e10),
            q_min=q_lo,
            q_max=q_hi,
            no_limits=self.no_limits_cb.isChecked(),
            n_mc_runs=self.n_runs_spin.value(),
        )

    def _write_config_back(self, config: ModelingConfig):
        """Write fitted config back into GUI widgets."""
        self.bg_edit.setText(_fmt(config.background))
        for i, pop in enumerate(config.populations):
            if i < len(self._pop_widgets):
                self._pop_widgets[i].from_population(pop)

    # ── Graph Model ──────────────────────────────────────────────────────────

    def graph_model(self):
        """Calculate and display the model without fitting."""
        if self._data_q is None:
            QMessageBox.information(self, 'No data', 'Please open an HDF5 file first.')
            return

        try:
            config = self._build_config()
            active = [p for p in config.populations if p.enabled]
            if not active:
                QMessageBox.information(self, 'No populations',
                                        'Enable at least one population tab.')
                return

            q_range = (config.q_min, config.q_max)
            mask = (self._data_q >= q_range[0]) & (self._data_q <= q_range[1])
            q = self._data_q[mask]
            if len(q) == 0:
                QMessageBox.warning(self, 'Empty range',
                                    'No data points in the selected Q range.')
                return

            self._engine.clear_cache()
            I_total, pop_idx, pop_I, pop_dist = self._engine.total_intensity(
                config, q, use_cache=True
            )

            # Build a mock ModelingResult for plotting
            r_grids, v_dists, n_dists = [], [], []
            derived = []
            for (rg, vd, nd), pi in zip(pop_dist, pop_idx):
                r_grids.append(rg)
                v_dists.append(vd)
                n_dists.append(nd)
                from pyirena.core.modeling import compute_derived
                derived.append(compute_derived(rg, vd, nd, config.populations[pi]))

            from datetime import datetime
            mock = ModelingResult(
                config=config,
                chi_squared=0.0, reduced_chi_squared=0.0, dof=0,
                timestamp=datetime.now().isoformat(timespec='seconds'),
                model_q=q, model_I=I_total,
                pop_indices=pop_idx, pop_model_I=pop_I,
                radius_grids=r_grids, volume_dists=v_dists,
                number_dists=n_dists, derived=derived,
            )

            # Replot data then overlay model
            self.graph.iq_plot.clear()
            self.graph._pop_items.clear()
            self.graph._total_item = None
            if self.graph._cursor_left_line is not None:
                self.graph.iq_plot.addItem(self.graph._cursor_left_line)
                self.graph.iq_plot.addItem(self.graph._cursor_right_line)

            plot_iq_data(self.graph.iq_plot, self._data_q, self._data_I,
                         self._data_dI, label='Data')
            set_robust_y_range(self.graph.iq_plot, self._data_I)
            self.graph.plot_model(mock)
            self.graph.set_status('Model graphed successfully.', 'success')

        except Exception as e:
            import traceback
            self.graph.set_status(f'Error: {e}', 'error')
            traceback.print_exc()

    # ── Fit ──────────────────────────────────────────────────────────────────

    def run_fit(self):
        """Run the fitting routine and update GUI with results."""
        if self._data_q is None:
            QMessageBox.information(self, 'No data', 'Please open an HDF5 file first.')
            return

        try:
            config = self._build_config()
            active = [p for p in config.populations if p.enabled]
            if not active:
                QMessageBox.information(self, 'No populations',
                                        'Enable at least one population tab.')
                return

            self.btn_fit.setEnabled(False)
            self.btn_graph.setEnabled(False)
            self.graph.set_status('Fitting … please wait.', 'info')
            QApplication.processEvents()

            self._engine.clear_cache()
            result = self._engine.fit(
                config, self._data_q, self._data_I, self._data_dI
            )
            self._last_result = result

            # Update GUI with best-fit parameters
            self._write_config_back(result.config)
            self.bg_edit.setText(_fmt(result.config.background))

            # Chi² display
            self.chi2_lbl.setText(f'{result.chi_squared:.4f}')
            self.rchi2_lbl.setText(f'{result.reduced_chi_squared:.4f}')
            self.dof_lbl.setText(str(result.dof))

            # Plot
            self.graph.iq_plot.clear()
            self.graph._pop_items.clear()
            self.graph._total_item = None
            if self.graph._cursor_left_line is not None:
                self.graph.iq_plot.addItem(self.graph._cursor_left_line)
                self.graph.iq_plot.addItem(self.graph._cursor_right_line)

            plot_iq_data(self.graph.iq_plot, self._data_q, self._data_I,
                         self._data_dI, label='Data')
            set_robust_y_range(self.graph.iq_plot, self._data_I)
            self.graph.plot_model(result)

            # Residuals: only over fitted Q range
            mask = ((self._data_q >= result.config.q_min) &
                    (self._data_q <= result.config.q_max))
            q_fit = self._data_q[mask]
            I_fit = self._data_I[mask]
            dI_fit = self._data_dI[mask]
            with np.errstate(invalid='ignore', divide='ignore'):
                resid = (I_fit - result.model_I) / np.maximum(dI_fit, 1e-30)
            self.graph.plot_residuals(q_fit, resid)

            self.btn_mc.setEnabled(True)
            self.btn_save.setEnabled(True)
            self.btn_export.setEnabled(True)

            msg = (f'Fit done.  χ² = {result.chi_squared:.4f},  '
                   f'χ²/dof = {result.reduced_chi_squared:.4f},  '
                   f'DOF = {result.dof}')
            self.graph.set_status(msg, 'success')

        except Exception as e:
            import traceback
            self.graph.set_status(f'Fit error: {e}', 'error')
            traceback.print_exc()
        finally:
            self.btn_fit.setEnabled(True)
            self.btn_graph.setEnabled(True)

        self._save_state()

    # ── MC uncertainty ────────────────────────────────────────────────────────

    def calc_uncertainty_mc(self):
        """Run Monte-Carlo uncertainty estimation in the GUI thread."""
        if self._last_result is None or self._data_q is None:
            QMessageBox.information(self, 'No result', 'Run a fit first.')
            return

        n = self.n_runs_spin.value()
        self.btn_mc.setEnabled(False)
        self.graph.set_status(f'Running {n} MC passes … please wait.', 'info')
        QApplication.processEvents()

        try:
            config = deepcopy(self._last_result.config)
            stds = self._engine.calculate_uncertainty_mc(
                config, self._data_q, self._data_I, self._data_dI, n_runs=n,
            )
            self._last_result.params_std = stds

            # Show results as a simple dialog
            if stds:
                lines = [f'  {k}: ± {v:.4g}' for k, v in sorted(stds.items())]
                QMessageBox.information(
                    self, 'MC Uncertainties',
                    'Parameter standard deviations:\n' + '\n'.join(lines),
                )
            else:
                QMessageBox.information(self, 'MC Uncertainties',
                                        'No fittable parameters or too few successful runs.')
        except Exception as e:
            QMessageBox.critical(self, 'MC error', str(e))
        finally:
            self.btn_mc.setEnabled(True)
            self.graph.set_status('MC uncertainty estimation complete.', 'success')

    # ── Save / Export ─────────────────────────────────────────────────────────

    def save_results(self):
        """Save ModelingResult to the current HDF5 file."""
        if self._last_result is None:
            return
        if self._file_path is None:
            self._file_path = Path(QFileDialog.getSaveFileName(
                self, 'Save HDF5', '', 'HDF5 (*.hdf5 *.h5);;All (*)'
            )[0])
            if not self._file_path.name:
                return
        try:
            save_modeling_results(self._file_path, self._last_result)
            self.graph.set_status(f'Results saved to {self._file_path.name}', 'success')
        except Exception as e:
            QMessageBox.critical(self, 'Save error', str(e))

    def export_json(self):
        """Export current population parameters as a JSON config file."""
        import json
        path, _ = QFileDialog.getSaveFileName(
            self, 'Export Parameters', '', 'JSON (*.json);;All (*)',
        )
        if not path:
            return
        config = self._build_config()
        data = {
            'tool': 'modeling',
            'background': config.background,
            'fit_background': config.fit_background,
            'no_limits': config.no_limits,
            'populations': [_pop_to_dict(p) for p in config.populations],
        }
        try:
            with open(path, 'w') as f:
                json.dump(data, f, indent=2)
            self.graph.set_status(f'Parameters exported to {Path(path).name}', 'success')
        except Exception as e:
            QMessageBox.critical(self, 'Export error', str(e))


# ──────────────────────────────────────────────────────────────────────────────
# Serialisation helpers for populations ↔ state dicts
# ──────────────────────────────────────────────────────────────────────────────

def _pop_to_dict(pop: SizeDistPopulation) -> dict:
    return {
        'enabled':          pop.enabled,
        'dist_type':        pop.dist_type,
        'dist_params':      pop.dist_params,
        'dist_params_fit':  pop.dist_params_fit,
        'dist_params_limits': {k: list(v) for k, v in pop.dist_params_limits.items()},
        'form_factor':      pop.form_factor,
        'ff_params':        pop.ff_params,
        'ff_params_fit':    pop.ff_params_fit,
        'ff_params_limits': {k: list(v) for k, v in pop.ff_params_limits.items()},
        'structure_factor': pop.structure_factor,
        'sf_params':        pop.sf_params,
        'sf_params_fit':    pop.sf_params_fit,
        'sf_params_limits': {k: list(v) for k, v in pop.sf_params_limits.items()},
        'contrast':         pop.contrast,
        'fit_contrast':     pop.fit_contrast,
        'contrast_limits':  list(pop.contrast_limits),
        'scale':            pop.scale,
        'fit_scale':        pop.fit_scale,
        'scale_limits':     list(pop.scale_limits),
        'use_number_dist':  pop.use_number_dist,
        'n_bins':           pop.n_bins,
    }


def _pop_from_dict(d: dict) -> SizeDistPopulation:
    pop = SizeDistPopulation()
    pop.enabled          = d.get('enabled', False)
    pop.dist_type        = d.get('dist_type', 'lognormal')
    pop.dist_params      = d.get('dist_params', dict(DIST_DEFAULTS['lognormal']))
    pop.dist_params_fit  = d.get('dist_params_fit', {})
    pop.dist_params_limits = {k: tuple(v) for k, v in
                               d.get('dist_params_limits', {}).items()}
    pop.form_factor      = d.get('form_factor', 'sphere')
    pop.ff_params        = d.get('ff_params', {})
    pop.ff_params_fit    = d.get('ff_params_fit', {})
    pop.ff_params_limits = {k: tuple(v) for k, v in
                             d.get('ff_params_limits', {}).items()}
    pop.structure_factor = d.get('structure_factor', 'none')
    pop.sf_params        = d.get('sf_params', {})
    pop.sf_params_fit    = d.get('sf_params_fit', {})
    pop.sf_params_limits = {k: tuple(v) for k, v in
                             d.get('sf_params_limits', {}).items()}
    pop.contrast         = float(d.get('contrast', 1.0))
    pop.fit_contrast     = bool(d.get('fit_contrast', False))
    cl = d.get('contrast_limits', [0.0, 1e10])
    pop.contrast_limits  = (cl[0], cl[1])
    pop.scale            = float(d.get('scale', 0.001))
    pop.fit_scale        = bool(d.get('fit_scale', True))
    sl = d.get('scale_limits', [1e-8, 1.0])
    pop.scale_limits     = (sl[0], sl[1])
    pop.use_number_dist  = bool(d.get('use_number_dist', False))
    pop.n_bins           = int(d.get('n_bins', 200))
    return pop


# ──────────────────────────────────────────────────────────────────────────────
# StateManager helper (get/set section)
# ──────────────────────────────────────────────────────────────────────────────
# StateManager does not expose get/set directly; patch via state dict attribute.

def _sm_get(sm: StateManager, section: str, default=None):
    return sm.state.get(section, default if default is not None else {})


def _sm_set(sm: StateManager, section: str, value):
    sm.state[section] = value


# Monkey-patch if not available (StateManager usually has these):
if not hasattr(StateManager, 'get'):
    StateManager.get = lambda self, k, d=None: self.state.get(k, d)
if not hasattr(StateManager, 'set'):
    StateManager.set = lambda self, k, v: self.state.update({k: v})


# ──────────────────────────────────────────────────────────────────────────────
# CLI entry
# ──────────────────────────────────────────────────────────────────────────────

def main():
    """Standalone entry: python -m pyirena.gui.modeling_panel [hdf5_file]"""
    app = QApplication.instance() or QApplication(sys.argv)
    panel = ModelingPanel()
    panel.show()

    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
        if path.exists():
            panel.load_file(path)
        else:
            print(f'File not found: {path}')

    sys.exit(app.exec())


if __name__ == '__main__':
    main()
