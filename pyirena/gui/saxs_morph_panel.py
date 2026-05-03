"""
SAXS Morph panel GUI for pyIrena.

Provides the left-panel (controls) + right-panel (graphs + 3D) window for
generating a 3D voxelgram of a two-phase porous structure from experimental
I(Q) using the Gaussian Random Fields method.

Phase 2 of the implementation:
  - Left panel: Voxel grid box, Two-phase parameters box, Background tabs
    (Power-law + Flat), action buttons.
  - Right panel: I(Q) plot with three traces (data, data-bg, model) and two
    cursors; placeholder area below for the 2D slice + 3D viewer added in
    Phase 3.
  - "Graph Model" button is wired (synchronous evaluation).
  - "Fit"/"Cancel"/"MC uncertainty" buttons are constructed but disabled
    until Phase 4 wires the QThread workers.

Entry points
------------
SaxsMorphPanel  — main QWidget; reusable as a tab in the Data Selector.
main()          — CLI: python -m pyirena.gui.saxs_morph_panel <file>
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
        QScrollArea, QFrame, QSizePolicy,
    )
    from PySide6.QtCore import Qt, Signal
    from PySide6.QtGui import QFont, QDoubleValidator
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
            QPushButton, QLabel, QLineEdit, QCheckBox, QSpinBox, QTabWidget,
            QGroupBox, QMessageBox, QSplitter, QFileDialog, QComboBox,
            QScrollArea, QFrame, QSizePolicy,
        )
        from PyQt6.QtCore import Qt, pyqtSignal as Signal
        from PyQt6.QtGui import QFont, QDoubleValidator
    except ImportError:
        from PyQt5.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
            QPushButton, QLabel, QLineEdit, QCheckBox, QSpinBox, QTabWidget,
            QGroupBox, QMessageBox, QSplitter, QFileDialog, QComboBox,
            QScrollArea, QFrame, QSizePolicy,
        )
        from PyQt5.QtCore import Qt, pyqtSignal as Signal
        from PyQt5.QtGui import QFont, QDoubleValidator

import pyqtgraph as pg

from pyirena.core.saxs_morph import (
    SaxsMorphEngine, SaxsMorphConfig, SaxsMorphResult,
    ALLOWED_VOXEL_SIZES, MAX_FIT_VOXEL_SIZE,
)
from pyirena.io.nxcansas_saxs_morph import (
    save_saxs_morph_results, load_saxs_morph_results,
)
from pyirena.gui.sas_plot import (
    make_sas_plot, plot_iq_data, set_robust_y_range, add_plot_annotation,
    RadiusAxisItem, save_itx_from_plot, SASPlotStyle,
)
from pyirena.gui.unified_fit import _SafeInfiniteLine
from pyirena.gui.saxs_morph_3d import (
    Voxel3DViewer, Slice2DViewer, make_popout_button,
    HAS_PYVISTA, PYVISTA_INSTALL_HINT,
)
from pyirena.state import StateManager


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _sep(orientation='h') -> QFrame:
    f = QFrame()
    f.setFrameShape(QFrame.Shape.HLine if orientation == 'h' else QFrame.Shape.VLine)
    f.setFrameShadow(QFrame.Shadow.Sunken)
    f.setStyleSheet('color: #cccccc;')
    return f


def _fmt(v) -> str:
    if v is None:
        return ''
    try:
        v = float(v)
    except (ValueError, TypeError):
        return str(v)
    if v == 0.0:
        return '0'
    if abs(v) < 0.01 or abs(v) >= 1e5:
        return f'{v:.4g}'
    return f'{v:.4g}'


def _parse(text: str, default: float = 0.0) -> float:
    try:
        return float(text)
    except (ValueError, TypeError):
        return default


def _parse_int_optional(text: str) -> Optional[int]:
    if not text or not text.strip():
        return None
    try:
        return int(text.strip())
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Graph window: I(Q) plot with cursors + placeholder for 3D viewers
# ---------------------------------------------------------------------------

class SaxsMorphGraphWindow(QWidget):
    """Right-side graph area: top I(Q) plot + bottom placeholder (Phase 3)."""

    cursor_moved = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.data_folder = '.'

        self._cursor_left = None
        self._cursor_right = None
        self._cursor_updating = False

        self._data_items = []
        self._corr_item = None
        self._model_item = None
        self._annotation_items: list = []

        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.splitter = QSplitter(Qt.Orientation.Vertical)

        # ── Top: I(Q) plot ────────────────────────────────────────────────
        top = QWidget()
        top_lay = QVBoxLayout(top)
        top_lay.setContentsMargins(0, 0, 0, 0)
        self.gl = pg.GraphicsLayoutWidget()
        self.gl.setBackground('w')
        top_lay.addWidget(self.gl)

        self.iq_plot = make_sas_plot(
            self.gl, row=0, col=0,
            title='SAXS Morph — I(Q)',
            parent_widget=self,
            jpeg_default_name='saxs_morph_iq.jpg',
        )

        self.splitter.addWidget(top)

        # ── Bottom: 2D slice + 3D viewer side by side ────────────────────
        self.viewer_container = QWidget()
        viewer_row = QHBoxLayout(self.viewer_container)
        viewer_row.setContentsMargins(4, 4, 4, 4)
        viewer_row.setSpacing(6)

        # Left half: slice viewer + popout
        slice_box = QWidget()
        slice_col = QVBoxLayout(slice_box)
        slice_col.setContentsMargins(0, 0, 0, 0)
        self.slice_viewer = Slice2DViewer()
        slice_col.addWidget(self.slice_viewer)
        slice_col.addWidget(make_popout_button(self.slice_viewer, '2D slice viewer'))
        viewer_row.addWidget(slice_box, 1)

        # Right half: 3D viewer + popout
        voxel_box = QWidget()
        voxel_col = QVBoxLayout(voxel_box)
        voxel_col.setContentsMargins(0, 0, 0, 0)
        self.voxel3d_viewer = Voxel3DViewer()
        voxel_col.addWidget(self.voxel3d_viewer)
        voxel_col.addWidget(make_popout_button(self.voxel3d_viewer, '3D viewer'))
        viewer_row.addWidget(voxel_box, 1)

        self.splitter.addWidget(self.viewer_container)

        # Top : bottom = 3 : 5 so the 3D viewer gets adequate room.
        self.splitter.setSizes([360, 600])
        layout.addWidget(self.splitter)

        # ── Status bar ────────────────────────────────────────────────────
        self.status_lbl = QLabel('')
        self.status_lbl.setWordWrap(True)
        self.status_lbl.setMaximumHeight(50)
        self.status_lbl.setStyleSheet(
            'padding: 6px; border: 1px solid #ccc; border-radius: 4px; font-size: 10pt;'
        )
        layout.addWidget(self.status_lbl)

    # ── Cursors ──────────────────────────────────────────────────────────

    def add_cursors(self, q_min_log: float, q_max_log: float):
        if self._cursor_left is not None:
            try:
                self.iq_plot.removeItem(self._cursor_left)
                self.iq_plot.removeItem(self._cursor_right)
            except Exception:
                pass

        self._cursor_left = _SafeInfiniteLine(
            pos=q_min_log, angle=90, movable=True,
            pen=pg.mkPen('#e74c3c', width=2),
            label='Qmin',
        )
        self._cursor_right = _SafeInfiniteLine(
            pos=q_max_log, angle=90, movable=True,
            pen=pg.mkPen('#2980b9', width=2),
            label='Qmax',
        )
        self._cursor_left.sigPositionChanged.connect(self._on_cursor_moved)
        self._cursor_right.sigPositionChanged.connect(self._on_cursor_moved)
        self.iq_plot.addItem(self._cursor_left)
        self.iq_plot.addItem(self._cursor_right)

    def _on_cursor_moved(self, _=None):
        if not self._cursor_updating:
            self._cursor_updating = True
            self.cursor_moved.emit()
            self._cursor_updating = False

    def get_q_range(self) -> tuple[float, float]:
        if self._cursor_left is None:
            return (0.001, 1.0)
        lo = 10 ** self._cursor_left.getPos()[0]
        hi = 10 ** self._cursor_right.getPos()[0]
        return (min(lo, hi), max(lo, hi))

    def set_q_range(self, q_min: float, q_max: float):
        if self._cursor_left is not None and q_min > 0 and q_max > 0:
            self._cursor_left.setPos(np.log10(q_min))
            self._cursor_right.setPos(np.log10(q_max))

    # ── Data plotting ────────────────────────────────────────────────────

    def plot_data(self, q, I, dI=None):
        """Plot raw experimental data; clears previous items."""
        self.iq_plot.clear()
        self._data_items.clear()
        self._corr_item = None
        self._model_item = None
        self.clear_annotations()

        scatter, errbar = plot_iq_data(self.iq_plot, q, I, dI, label='Data')
        self._data_items = [scatter, errbar]

    def plot_data_minus_bg(self, q, I_corr):
        """Overlay data - background as a black scatter."""
        if self._corr_item is not None:
            try:
                self.iq_plot.removeItem(self._corr_item)
            except Exception:
                pass
        I_pos = np.where(I_corr > 0, I_corr, np.nan)
        self._corr_item = self.iq_plot.plot(
            q, I_pos, pen=None,
            symbol='t', symbolSize=5,
            symbolBrush='#222', symbolPen=None,
            name='Data − Background',
        )

    def plot_model_iq(self, q, I_model):
        """Overlay the GRF model as a red line."""
        if self._model_item is not None:
            try:
                self.iq_plot.removeItem(self._model_item)
            except Exception:
                pass
        self._model_item = self.iq_plot.plot(
            q, I_model,
            pen=pg.mkPen('#c0392b', width=2),
            name='Model',
        )

    def add_annotation(self, text: str, corner: str = 'lower_left'):
        item = add_plot_annotation(self.iq_plot, text, corner=corner)
        self._annotation_items.append(item)

    def show_voxelgram(self, voxelgram, pitch_A: float):
        """Push a fresh voxelgram into both the 2D slice and 3D viewers."""
        self.slice_viewer.set_voxelgram(voxelgram, pitch_A)
        self.voxel3d_viewer.set_voxelgram(voxelgram, pitch_A)

    def clear_annotations(self):
        for item in list(self._annotation_items):
            try:
                self.iq_plot.removeItem(item)
            except Exception:
                pass
        self._annotation_items.clear()

    def set_status(self, msg: str, style: str = 'info'):
        colours = {
            'info':    ('', '#333333'),
            'working': ('#fff3cd', '#856404'),
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


# ---------------------------------------------------------------------------
# Helper widgets for parameter rows
# ---------------------------------------------------------------------------

class ParamRow(QWidget):
    """One labelled parameter row: name | value | Fit? | lo | hi.

    All edits emit ``changed``.  ``no_limits=True`` hides lo/hi columns.
    """
    changed = Signal()

    def __init__(self, label: str, value: float, fit_flag: bool,
                 limits: tuple, parent=None):
        super().__init__(parent)
        self._building = True
        self._build_ui(label, value, fit_flag, limits)
        self._building = False

    def _build_ui(self, label, value, fit_flag, limits):
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(4)

        self.lbl = QLabel(label)
        self.lbl.setMinimumWidth(120)
        row.addWidget(self.lbl)

        self.val = QLineEdit(_fmt(value))
        self.val.setMaximumWidth(90)
        self.val.setValidator(QDoubleValidator())
        self.val.editingFinished.connect(self._emit_changed)
        row.addWidget(self.val)

        self.fit_cb = QCheckBox('Fit?')
        self.fit_cb.setChecked(bool(fit_flag))
        self.fit_cb.stateChanged.connect(self._emit_changed)
        row.addWidget(self.fit_cb)

        self.lo = QLineEdit(_fmt(limits[0]))
        self.lo.setMaximumWidth(70)
        self.lo.setValidator(QDoubleValidator())
        self.lo.editingFinished.connect(self._emit_changed)
        self.lo_lbl = QLabel('lo:')
        row.addWidget(self.lo_lbl)
        row.addWidget(self.lo)

        self.hi = QLineEdit(_fmt(limits[1]))
        self.hi.setMaximumWidth(70)
        self.hi.setValidator(QDoubleValidator())
        self.hi.editingFinished.connect(self._emit_changed)
        self.hi_lbl = QLabel('hi:')
        row.addWidget(self.hi_lbl)
        row.addWidget(self.hi)

        row.addStretch()

    def _emit_changed(self, *_):
        if not self._building:
            self.changed.emit()

    # API
    def value(self) -> float:
        return _parse(self.val.text(), 0.0)

    def fit_flag(self) -> bool:
        return self.fit_cb.isChecked()

    def limits(self) -> tuple:
        return (_parse(self.lo.text(), 0.0), _parse(self.hi.text(), 1e10))

    def set_value(self, v):
        self._building = True
        self.val.setText(_fmt(v))
        self._building = False

    def set_fit(self, f: bool):
        self._building = True
        self.fit_cb.setChecked(bool(f))
        self._building = False

    def set_limits(self, lo, hi):
        self._building = True
        self.lo.setText(_fmt(lo))
        self.hi.setText(_fmt(hi))
        self._building = False

    def set_no_limits(self, no_lim: bool):
        for w in (self.lo, self.hi, self.lo_lbl, self.hi_lbl):
            w.setVisible(not no_lim)


# ---------------------------------------------------------------------------
# Main panel
# ---------------------------------------------------------------------------

class SaxsMorphPanel(QWidget):
    """SAXS Morph main panel: left controls + right graph window."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('pyIrena — SAXS Morph')
        self.setMinimumSize(1300, 900)

        self._data_q: Optional[np.ndarray] = None
        self._data_I: Optional[np.ndarray] = None
        self._data_dI: Optional[np.ndarray] = None
        self._file_path: Optional[Path] = None
        self._last_result: Optional[SaxsMorphResult] = None

        self._engine = SaxsMorphEngine()
        self._state = StateManager()

        self._building = True
        self._build_ui()
        self._building = False
        self._load_state()

    # ── UI construction ──────────────────────────────────────────────────

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self._build_left_panel())
        self.graph = SaxsMorphGraphWindow()
        splitter.addWidget(self.graph)
        splitter.setSizes([460, 880])
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 2)
        layout.addWidget(splitter)

        self.graph.cursor_moved.connect(self._on_cursor_moved)

    def _build_left_panel(self) -> QWidget:
        panel = QWidget()
        panel.setMinimumWidth(420)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)

        inner = QWidget()
        lay = QVBoxLayout(inner)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(6)

        # ── Title + Help ─────────────────────────────────────────────────
        title_row = QHBoxLayout()
        title = QLabel('SAXS Morph')
        title.setStyleSheet('font-size: 14px; font-weight: bold; color: #6c3483;')
        title_row.addWidget(title)
        title_row.addStretch()
        help_btn = QPushButton('? Help')
        help_btn.setFixedSize(60, 22)
        help_btn.setStyleSheet(
            'background:#c0392b;color:white;font-size:11px;border-radius:3px;')
        help_btn.clicked.connect(self._open_help)
        title_row.addWidget(help_btn)
        lay.addLayout(title_row)

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

        # ── Q range readout (driven by cursors) ──────────────────────────
        q_row = QHBoxLayout()
        q_row.addWidget(QLabel('Q range:'))
        self.qmin_lbl = QLabel('—')
        self.qmin_lbl.setStyleSheet('font-family: monospace;')
        self.qmax_lbl = QLabel('—')
        self.qmax_lbl.setStyleSheet('font-family: monospace;')
        q_row.addWidget(self.qmin_lbl)
        q_row.addWidget(QLabel('to'))
        q_row.addWidget(self.qmax_lbl)
        q_row.addWidget(QLabel('Å⁻¹'))
        q_row.addStretch()
        lay.addLayout(q_row)

        # ── Voxel grid ───────────────────────────────────────────────────
        lay.addWidget(self._build_voxel_grid_box())

        # ── Two-phase parameters ─────────────────────────────────────────
        lay.addWidget(self._build_two_phase_box())

        # ── Background tabs ──────────────────────────────────────────────
        lay.addWidget(self._build_background_tabs())

        # ── No-limits toggle ─────────────────────────────────────────────
        nl_row = QHBoxLayout()
        self.no_limits_cb = QCheckBox('No limits (use Nelder-Mead)')
        self.no_limits_cb.stateChanged.connect(self._on_no_limits_changed)
        nl_row.addWidget(self.no_limits_cb)
        nl_row.addStretch()
        lay.addLayout(nl_row)

        # ── Action row 1: Graph Model + Fit + Cancel ─────────────────────
        lay.addWidget(_sep())
        act1 = QHBoxLayout()
        self.btn_graph = QPushButton('Graph Model')
        self.btn_graph.setStyleSheet(
            'background:#52c77a;color:white;font-weight:bold;'
            'font-size:11pt;padding:6px;border-radius:4px;border:none;')
        self.btn_graph.clicked.connect(self._on_graph_model)
        act1.addWidget(self.btn_graph)

        self.btn_fit = QPushButton('Fit')
        self.btn_fit.setStyleSheet(
            'background:#27ae60;color:white;font-weight:bold;'
            'font-size:11pt;padding:6px;border-radius:4px;border:none;')
        self.btn_fit.setEnabled(False)
        self.btn_fit.setToolTip('Fit (wired in Phase 4)')
        act1.addWidget(self.btn_fit)

        self.btn_cancel = QPushButton('Cancel')
        self.btn_cancel.setEnabled(False)
        act1.addWidget(self.btn_cancel)
        lay.addLayout(act1)

        # ── Action row 2: MC uncertainty ─────────────────────────────────
        act2 = QHBoxLayout()
        act2.addWidget(QLabel('Passes:'))
        self.n_runs_spin = QSpinBox()
        self.n_runs_spin.setRange(1, 500)
        self.n_runs_spin.setValue(10)
        act2.addWidget(self.n_runs_spin)
        self.btn_mc = QPushButton('Calc. Uncertainty (MC)')
        self.btn_mc.setStyleSheet(
            'background:#16a085;color:white;font-weight:bold;'
            'font-size:11pt;padding:6px;border-radius:4px;border:none;')
        self.btn_mc.setEnabled(False)
        self.btn_mc.setToolTip('Monte-Carlo uncertainty (wired in Phase 4)')
        act2.addWidget(self.btn_mc)
        self.btn_revert = QPushButton('Revert')
        self.btn_revert.setEnabled(False)
        self.btn_revert.setToolTip('Restore parameters from before last fit (Phase 4)')
        act2.addWidget(self.btn_revert)
        lay.addLayout(act2)

        # ── Result block ─────────────────────────────────────────────────
        lay.addWidget(_sep())
        self.result_lbl = QLabel('No model computed yet.')
        self.result_lbl.setWordWrap(True)
        self.result_lbl.setStyleSheet(
            'background:#fafafa;border:1px solid #ddd;border-radius:4px;'
            'padding:6px;font-family:monospace;font-size:9pt;')
        lay.addWidget(self.result_lbl)

        # ── Save row ─────────────────────────────────────────────────────
        save_row = QHBoxLayout()
        self.btn_save = QPushButton('Save Result to HDF5…')
        self.btn_save.setStyleSheet(
            'background:#2980b9;color:white;font-weight:bold;'
            'font-size:11pt;padding:6px;border-radius:4px;border:none;')
        self.btn_save.setEnabled(False)
        self.btn_save.clicked.connect(self._save_result)
        save_row.addWidget(self.btn_save)
        lay.addLayout(save_row)

        lay.addStretch()
        scroll.setWidget(inner)

        wrap = QVBoxLayout(panel)
        wrap.setContentsMargins(0, 0, 0, 0)
        wrap.addWidget(scroll)
        return panel

    def _build_voxel_grid_box(self) -> QGroupBox:
        gb = QGroupBox('Voxel grid')
        g = QGridLayout(gb)
        g.setVerticalSpacing(4)

        g.addWidget(QLabel('Cube side (fit):'), 0, 0)
        self.fit_size_combo = QComboBox()
        for n in ALLOWED_VOXEL_SIZES:
            if n <= MAX_FIT_VOXEL_SIZE:
                self.fit_size_combo.addItem(f'{n}³', n)
        self.fit_size_combo.setCurrentIndex(1)  # 128
        self.fit_size_combo.setToolTip(
            f'Voxel cube side used during the fit loop. '
            f'Hard limit {MAX_FIT_VOXEL_SIZE}³ for memory safety.'
        )
        g.addWidget(self.fit_size_combo, 0, 1)

        g.addWidget(QLabel('Cube side (render):'), 1, 0)
        self.render_size_combo = QComboBox()
        for n in ALLOWED_VOXEL_SIZES:
            self.render_size_combo.addItem(f'{n}³', n)
        self.render_size_combo.setCurrentIndex(2)  # 256
        self.render_size_combo.setToolTip(
            'Final voxel cube used for the saved result and 3D viewer. '
            'Larger values cost RAM and disk; 256³ is a good default.'
        )
        g.addWidget(self.render_size_combo, 1, 1)

        g.addWidget(QLabel('Box size [Å]:'), 2, 0)
        self.box_size_edit = QLineEdit('1000')
        self.box_size_edit.setValidator(QDoubleValidator(1.0, 1e9, 4))
        self.box_size_edit.setToolTip(
            'Physical edge length of the cubic simulation box. '
            'Determines voxel pitch = box_size / cube_side.'
        )
        g.addWidget(self.box_size_edit, 2, 1)

        g.addWidget(QLabel('RNG seed:'), 3, 0)
        self.seed_edit = QLineEdit('')
        self.seed_edit.setPlaceholderText('(blank = random)')
        self.seed_edit.setToolTip(
            'Integer for reproducible voxelgrams. Leave blank to draw a '
            'fresh seed each Graph Model / Fit.'
        )
        g.addWidget(self.seed_edit, 3, 1)

        return gb

    def _build_two_phase_box(self) -> QGroupBox:
        gb = QGroupBox('Two-phase parameters')
        v = QVBoxLayout(gb)

        self.phi_row = ParamRow(
            'Volume fraction φ:', 0.30, True, (0.05, 0.95),
        )
        v.addWidget(self.phi_row)

        self.contrast_row = ParamRow(
            'Contrast (Δρ)²:', 1.0, False, (0.0, 1e10),
        )
        v.addWidget(self.contrast_row)

        link_row = QHBoxLayout()
        self.link_cb = QCheckBox('Link contrast to φ via invariant')
        self.link_cb.setChecked(True)
        self.link_cb.setToolTip(
            'When checked, contrast is derived from the data invariant '
            'Q* = 2π² φ(1−φ) Δρ² rather than being a free parameter.'
        )
        self.link_cb.stateChanged.connect(self._on_link_changed)
        link_row.addWidget(self.link_cb)
        link_row.addStretch()
        v.addLayout(link_row)

        return gb

    def _build_background_tabs(self) -> QTabWidget:
        tabs = QTabWidget()

        # ── Power-law tab ────────────────────────────────────────────────
        pl_tab = QWidget()
        pl_lay = QVBoxLayout(pl_tab)
        self.pl_B_row = ParamRow('Amplitude B:', 0.0, False, (0.0, 1e10))
        self.pl_P_row = ParamRow('Exponent P:', 4.0, False, (0.0, 6.0))
        pl_lay.addWidget(self.pl_B_row)
        pl_lay.addWidget(self.pl_P_row)
        pl_help = QLabel('I_pl(Q) = B · Q^(−P)')
        pl_help.setStyleSheet('color:#888;font-size:9pt;')
        pl_lay.addWidget(pl_help)
        pl_lay.addStretch()
        tabs.addTab(pl_tab, 'Power-law')

        # ── Flat tab ─────────────────────────────────────────────────────
        flat_tab = QWidget()
        flat_lay = QVBoxLayout(flat_tab)
        self.bg_row = ParamRow('Background:', 0.0, False, (0.0, 1e10))
        flat_lay.addWidget(self.bg_row)
        flat_help = QLabel('Constant offset added to the model.')
        flat_help.setStyleSheet('color:#888;font-size:9pt;')
        flat_lay.addWidget(flat_help)
        flat_lay.addStretch()
        tabs.addTab(flat_tab, 'Flat')

        return tabs

    # ── State (load/save) ────────────────────────────────────────────────

    def _current_state(self) -> dict:
        st = {
            'voxel_size_fit':    self.fit_size_combo.currentData(),
            'voxel_size_render': self.render_size_combo.currentData(),
            'box_size_A':        _parse(self.box_size_edit.text(), 1000.0),
            'rng_seed':          _parse_int_optional(self.seed_edit.text()),

            'volume_fraction':         self.phi_row.value(),
            'fit_volume_fraction':     self.phi_row.fit_flag(),
            'volume_fraction_limits':  list(self.phi_row.limits()),

            'contrast':         self.contrast_row.value(),
            'fit_contrast':     self.contrast_row.fit_flag(),
            'contrast_limits':  list(self.contrast_row.limits()),
            'link_phi_contrast': self.link_cb.isChecked(),

            'power_law_B':         self.pl_B_row.value(),
            'fit_power_law_B':     self.pl_B_row.fit_flag(),
            'power_law_B_limits':  list(self.pl_B_row.limits()),
            'power_law_P':         self.pl_P_row.value(),
            'fit_power_law_P':     self.pl_P_row.fit_flag(),
            'power_law_P_limits':  list(self.pl_P_row.limits()),

            'background':         self.bg_row.value(),
            'fit_background':     self.bg_row.fit_flag(),
            'background_limits':  list(self.bg_row.limits()),

            'no_limits': self.no_limits_cb.isChecked(),
            'n_mc_runs': self.n_runs_spin.value(),
            'q_min': None,
            'q_max': None,
        }
        if self._data_q is not None:
            q_lo, q_hi = self.graph.get_q_range()
            st['q_min'] = float(q_lo)
            st['q_max'] = float(q_hi)
        return st

    def _load_state(self):
        st = self._state.get('saxs_morph') or {}

        # Voxel grid
        self._select_combo_value(self.fit_size_combo, st.get('voxel_size_fit', 128))
        self._select_combo_value(self.render_size_combo, st.get('voxel_size_render', 256))
        self.box_size_edit.setText(_fmt(st.get('box_size_A', 1000.0)))
        seed = st.get('rng_seed')
        self.seed_edit.setText('' if seed is None else str(seed))

        # Two-phase
        self.phi_row.set_value(st.get('volume_fraction', 0.30))
        self.phi_row.set_fit(st.get('fit_volume_fraction', True))
        self.phi_row.set_limits(*st.get('volume_fraction_limits', [0.05, 0.95]))

        self.contrast_row.set_value(st.get('contrast', 1.0))
        self.contrast_row.set_fit(st.get('fit_contrast', False))
        self.contrast_row.set_limits(*st.get('contrast_limits', [0.0, 1e10]))
        self.link_cb.setChecked(st.get('link_phi_contrast', True))
        self._on_link_changed()

        # Background
        self.pl_B_row.set_value(st.get('power_law_B', 0.0))
        self.pl_B_row.set_fit(st.get('fit_power_law_B', False))
        self.pl_B_row.set_limits(*st.get('power_law_B_limits', [0.0, 1e10]))
        self.pl_P_row.set_value(st.get('power_law_P', 4.0))
        self.pl_P_row.set_fit(st.get('fit_power_law_P', False))
        self.pl_P_row.set_limits(*st.get('power_law_P_limits', [0.0, 6.0]))

        self.bg_row.set_value(st.get('background', 0.0))
        self.bg_row.set_fit(st.get('fit_background', False))
        self.bg_row.set_limits(*st.get('background_limits', [0.0, 1e10]))

        self.no_limits_cb.setChecked(st.get('no_limits', False))
        self.n_runs_spin.setValue(st.get('n_mc_runs', 10))
        self._on_no_limits_changed()

    def _save_state(self):
        self._state.update('saxs_morph', self._current_state())
        self._state.save()

    @staticmethod
    def _select_combo_value(combo: QComboBox, value):
        for i in range(combo.count()):
            if combo.itemData(i) == value:
                combo.setCurrentIndex(i)
                return
        # fallback: keep current selection

    # ── Signal handlers ──────────────────────────────────────────────────

    def _on_link_changed(self, *_):
        linked = self.link_cb.isChecked()
        self.contrast_row.setEnabled(not linked)
        if linked:
            self.contrast_row.set_fit(False)

    def _on_no_limits_changed(self, *_):
        nl = self.no_limits_cb.isChecked()
        for row in (self.phi_row, self.contrast_row,
                    self.pl_B_row, self.pl_P_row, self.bg_row):
            row.set_no_limits(nl)

    def _on_cursor_moved(self):
        q_lo, q_hi = self.graph.get_q_range()
        self.qmin_lbl.setText(f'{q_lo:.4g}')
        self.qmax_lbl.setText(f'{q_hi:.4g}')

    def _open_help(self):
        try:
            from PySide6.QtGui import QDesktopServices
            from PySide6.QtCore import QUrl
        except ImportError:
            from PyQt6.QtGui import QDesktopServices
            from PyQt6.QtCore import QUrl
        QDesktopServices.openUrl(QUrl(
            'https://github.com/jilavsky/pyirena/blob/main/docs/saxs_morph_gui.md'
        ))

    # ── File I/O ─────────────────────────────────────────────────────────

    def _open_file(self):
        last = (self._state.get('data_selector') or {}).get('last_folder', '')
        path, _ = QFileDialog.getOpenFileName(
            self, 'Open HDF5 file', last,
            'HDF5 files (*.hdf5 *.h5 *.nxs);;All files (*)',
        )
        if path:
            self.load_file(Path(path))

    def load_file(self, path: Path):
        """Load Q/I/dI from an NXcanSAS HDF5 file."""
        try:
            from pyirena.io.hdf5 import readGenericNXcanSAS
            data = readGenericNXcanSAS(str(path.parent), path.name)
            q = np.asarray(data['Q'], dtype=float)
            I = np.asarray(data['Intensity'], dtype=float)
            err = data.get('Error')
            dI = (np.asarray(err, dtype=float)
                  if err is not None else np.maximum(I * 0.05, 1e-30))
        except Exception as e:
            QMessageBox.critical(self, 'Open failed', f'Could not load:\n{e}')
            return
        self.set_data(q, I, dI, filename=path.name, filepath=str(path))

    def set_data(self, q, I, dI, filename: str, filepath: str = '',
                 is_nxcansas: bool = True):
        q = np.asarray(q, dtype=float)
        I = np.asarray(I, dtype=float)
        dI = np.asarray(dI, dtype=float) if dI is not None else np.maximum(I * 0.05, 1e-30)

        self._file_path = Path(filepath) if filepath else None
        self._data_q = q
        self._data_I = I
        self._data_dI = dI

        self.file_edit.setText(filename)
        if filepath:
            self.graph.data_folder = str(Path(filepath).parent)

        # Filter to positive Q
        mask = q > 0
        q_pos = q[mask]
        I_pos = I[mask]
        dI_pos = dI[mask]

        self.graph.plot_data(q_pos, I_pos, dI_pos)
        set_robust_y_range(self.graph.iq_plot, I_pos)

        q_min, q_max = float(q_pos.min()), float(q_pos.max())
        self.graph.add_cursors(np.log10(q_min), np.log10(q_max))

        # Restore saved Q range if it falls within data
        st = self._state.get('saxs_morph') or {}
        sq_min = st.get('q_min')
        sq_max = st.get('q_max')
        if sq_min and sq_max and q_min <= sq_min < sq_max <= q_max:
            self.graph.set_q_range(sq_min, sq_max)

        self._on_cursor_moved()

    # ── Action: Graph Model ──────────────────────────────────────────────

    def _make_config(self) -> SaxsMorphConfig:
        st = self._current_state()
        cfg = SaxsMorphConfig(
            q_min=st['q_min'], q_max=st['q_max'],
            voxel_size_fit=int(st['voxel_size_fit']),
            voxel_size_render=int(st['voxel_size_render']),
            box_size_A=float(st['box_size_A']),
            volume_fraction=float(st['volume_fraction']),
            fit_volume_fraction=bool(st['fit_volume_fraction']),
            volume_fraction_limits=tuple(st['volume_fraction_limits']),
            contrast=float(st['contrast']),
            fit_contrast=bool(st['fit_contrast']),
            contrast_limits=tuple(st['contrast_limits']),
            link_phi_contrast=bool(st['link_phi_contrast']),
            power_law_B=float(st['power_law_B']),
            fit_power_law_B=bool(st['fit_power_law_B']),
            power_law_B_limits=tuple(st['power_law_B_limits']),
            power_law_P=float(st['power_law_P']),
            fit_power_law_P=bool(st['fit_power_law_P']),
            power_law_P_limits=tuple(st['power_law_P_limits']),
            background=float(st['background']),
            fit_background=bool(st['fit_background']),
            background_limits=tuple(st['background_limits']),
            no_limits=bool(st['no_limits']),
            n_mc_runs=int(st['n_mc_runs']),
            rng_seed=st['rng_seed'],
        )
        return cfg

    def _on_graph_model(self):
        if self._data_q is None:
            QMessageBox.warning(self, 'No data', 'Open a data file first.')
            return

        try:
            cfg = self._make_config()
        except Exception as e:
            QMessageBox.warning(self, 'Bad parameters', str(e))
            return

        self.graph.set_status(
            f'Computing voxelgram at {cfg.voxel_size_fit}³ … please wait.',
            style='working',
        )
        QApplication.processEvents()

        try:
            result = self._engine.compute_voxelgram(
                cfg, self._data_q, self._data_I, self._data_dI,
            )
        except Exception as e:
            self.graph.set_status(f'Failed: {e}', style='error')
            QMessageBox.critical(self, 'Compute failed', str(e))
            return

        self._last_result = result
        self._update_after_eval(result, label='Graph Model')
        self._save_state()

    def _update_after_eval(self, result: SaxsMorphResult, label: str):
        # Plot data − bg + model
        self.graph.plot_data_minus_bg(result.data_q, result.data_I_corr)
        self.graph.plot_model_iq(result.model_q, result.model_I)
        # Push the voxelgram to the 2D slice + 3D viewers
        self.graph.show_voxelgram(result.voxelgram, result.voxel_pitch_A)
        self.graph.set_status(
            f'{label}: χ² = {result.chi_squared:.4g}, '
            f'φ_actual = {result.phi_actual:.4g}, '
            f'voxel = {result.voxel_size}³',
            style='success',
        )
        # Result block
        cfg = result.config
        text = (
            f"χ²            = {result.chi_squared:.6g}\n"
            f"red. χ²       = {result.reduced_chi_squared:.6g}\n"
            f"dof           = {result.dof}\n"
            f"voxel size    = {result.voxel_size}³\n"
            f"box size      = {result.box_size_A:.4g} Å\n"
            f"voxel pitch   = {result.voxel_pitch_A:.4g} Å\n"
            f"φ (target)    = {cfg.volume_fraction:.4g}\n"
            f"φ (actual)    = {result.phi_actual:.4g}\n"
            f"contrast Δρ²  = {cfg.contrast:.4g}"
            f"  ({'derived' if cfg.link_phi_contrast else 'set'})\n"
            f"power-law B   = {cfg.power_law_B:.4g}\n"
            f"power-law P   = {cfg.power_law_P:.4g}\n"
            f"flat bg       = {cfg.background:.4g}\n"
            f"RNG seed used = {result.rng_seed_used}"
        )
        self.result_lbl.setText(text)
        self.btn_save.setEnabled(True)

    def _save_result(self):
        if self._last_result is None:
            return
        default = (self._file_path.with_suffix('.h5')
                   if self._file_path else Path('saxs_morph_result.h5'))
        path, _ = QFileDialog.getSaveFileName(
            self, 'Save SAXS Morph result', str(default),
            'HDF5 files (*.hdf5 *.h5);;All files (*)',
        )
        if not path:
            return
        try:
            save_saxs_morph_results(Path(path), self._last_result)
            QMessageBox.information(self, 'Saved', f'Result saved to:\n{path}')
        except Exception as e:
            QMessageBox.critical(self, 'Save failed', str(e))


# ---------------------------------------------------------------------------
# CLI entry
# ---------------------------------------------------------------------------

def main():
    app = QApplication.instance() or QApplication(sys.argv)
    panel = SaxsMorphPanel()
    panel.show()
    if len(sys.argv) > 1:
        panel.load_file(Path(sys.argv[1]))
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
