"""
Fractals tool — left panel + right graph window for the mass-fractal
aggregate generator/analyzer.

GUI mirrors `saxs_morph_panel.py`: scrollable left panel + right vertical
splitter (top: log-log I(Q), bottom: horizontal splitter with `Slice2DViewer`
and `Voxel3DViewer`).  Aggregate growth runs in a background QThread queue
so multiple aggregates can be grown while the user inspects already-completed
ones.  Aggregates live in an in-session list; explicit save writes a single
aggregate to a NeXus file.

Public entry points
-------------------
FractalsGraphWindow  : main window — `data_selector` constructs and shows it.
"""

from __future__ import annotations
import logging

log = logging.getLogger(__name__)


import math
from pathlib import Path
from typing import Optional

import numpy as np

from pyirena.gui._qt import (
    QAbstractItemView, QAction, QComboBox, QDoubleSpinBox, QFileDialog, QFrame, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QListWidget, QListWidgetItem, QMainWindow, QMenu, QMessageBox, QPushButton, QScrollArea, QSpinBox, QSplitter, QTabWidget, QVBoxLayout, QWidget, Qt,
)

import pyqtgraph as pg

from pyirena.core.fractals import (
    FractalAggregate, GrowthConfig, OptimizerConfig,
    intensity_unified, voxelize,
)
from pyirena.gui.fractals_workers import GrowthQueueWorker, MCIntensityWorker
from pyirena.gui.sas_plot import (
    make_sas_plot, plot_iq_data, plot_iq_model, set_robust_y_range, SASPlotStyle,
)
from pyirena.gui.saxs_morph_3d import (
    Voxel3DViewer, Slice2DViewer, make_popout_button,
)
from pyirena.io.nxcansas_fractals import (
    save_fractal_aggregate, list_fractal_aggregates, load_fractal_aggregate,
)
from pyirena.io.nxcansas_unified import load_unified_fit_results
from pyirena.state import StateManager


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _sep() -> QFrame:
    f = QFrame()
    f.setFrameShape(QFrame.Shape.HLine)
    f.setFrameShadow(QFrame.Shadow.Sunken)
    f.setStyleSheet("color: #cccccc;")
    return f


def _fmt(v) -> str:
    if v is None:
        return "—"
    try:
        v = float(v)
    except (ValueError, TypeError):
        return str(v)
    if not math.isfinite(v):
        return "—"
    if v == 0.0:
        return "0"
    if abs(v) < 0.01 or abs(v) >= 1e5:
        return f"{v:.4g}"
    return f"{v:.4g}"


_ATTRACTION_OPTIONS = ["Neutral", "Attractive", "Repulsive", "Not allowed"]
_NEAR_DIST_OPTIONS = [
    ("Edge (a, 6 neighbors)", 1),
    ("Face diagonal (a·√2, 18 neighbors)", 2),
    ("Body diagonal (a·√3, 26 neighbors)", 3),
]


# ---------------------------------------------------------------------------
# Right-side graph window
# ---------------------------------------------------------------------------

class FractalsGraphWindow(QMainWindow):
    """Main window: left controls + right graph area."""

    def __init__(self, state_manager: Optional[StateManager] = None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena — Fractals")
        self.setMinimumSize(1300, 900)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)

        self._state = state_manager or StateManager()

        # Workers
        self._growth_worker = GrowthQueueWorker(self)
        self._growth_worker.start()
        self._mc_worker = MCIntensityWorker(self)

        # Build UI
        self.panel = FractalsPanel(self._state, self._growth_worker, self._mc_worker, parent=self)
        self.setCentralWidget(self.panel)

    def closeEvent(self, evt):
        try:
            self._growth_worker.shutdown()
        except Exception:
            log.debug("suppressed exception", exc_info=True)
        try:
            self._mc_worker.shutdown()
        except Exception:
            log.debug("suppressed exception", exc_info=True)
        super().closeEvent(evt)


# ---------------------------------------------------------------------------
# Main panel
# ---------------------------------------------------------------------------

class FractalsPanel(QWidget):
    """Left controls + right graph area."""

    def __init__(self, state_manager: StateManager,
                 growth_worker: GrowthQueueWorker,
                 mc_worker: MCIntensityWorker,
                 parent=None):
        super().__init__(parent)
        self._state = state_manager
        self._growth_worker = growth_worker
        self._mc_worker = mc_worker

        # State
        self._aggregates: list[FractalAggregate] = []   # in-session list
        self._active: Optional[FractalAggregate] = None
        self._jobs: dict[str, dict] = {}                # job_id → {label, list_item}
        self._loaded_q: Optional[np.ndarray] = None
        self._loaded_I: Optional[np.ndarray] = None
        self._loaded_dI: Optional[np.ndarray] = None
        self._loaded_unified_q: Optional[np.ndarray] = None
        self._loaded_unified_I: Optional[np.ndarray] = None
        self._targets: dict = {}      # populated when a Unified-fit pair is detected

        self._building = True
        self._build_ui()
        self._building = False
        self._load_state()
        self._wire_workers()

    # ── UI construction ──────────────────────────────────────────────────

    def _build_ui(self):
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(4, 4, 4, 4)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self._build_left_panel())
        splitter.addWidget(self._build_right_graph_area())
        splitter.setSizes([480, 880])
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        main_layout.addWidget(splitter)

    # ── Left panel (controls) ────────────────────────────────────────────

    def _build_left_panel(self) -> QWidget:
        panel = QWidget()
        panel.setMinimumWidth(440)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)

        inner = QWidget()
        lay = QVBoxLayout(inner)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(6)

        # Title + Help
        title_row = QHBoxLayout()
        title = QLabel("Fractals (Mass Fractal Aggregate Model)")
        title.setStyleSheet("font-size: 14px; font-weight: bold; color: #16a085;")
        title_row.addWidget(title)
        title_row.addStretch()
        help_btn = QPushButton("? Help")
        help_btn.setFixedSize(60, 22)
        help_btn.setStyleSheet(
            "background:#c0392b;color:white;font-size:11px;border-radius:3px;")
        help_btn.clicked.connect(self._open_help)
        title_row.addWidget(help_btn)
        lay.addLayout(title_row)
        lay.addWidget(_sep())

        # Optional NeXus
        lay.addWidget(self._build_nexus_box())

        # Growth parameters
        lay.addWidget(self._build_growth_box())

        # Mode tabs (Grow One / Grow Many / Optimizer)
        lay.addWidget(self._build_mode_tabs())

        # Job queue status
        lay.addWidget(self._build_jobs_box())

        # Stored aggregates list
        lay.addWidget(self._build_aggregates_list())

        # Active aggregate parameters
        lay.addWidget(self._build_params_box())

        # Intensity (Q range + MC)
        lay.addWidget(self._build_intensity_box())

        # Save
        lay.addWidget(self._build_save_box())

        lay.addStretch()
        scroll.setWidget(inner)

        wrap = QVBoxLayout(panel)
        wrap.setContentsMargins(0, 0, 0, 0)
        wrap.addWidget(scroll)
        return panel

    def _build_nexus_box(self) -> QGroupBox:
        gb = QGroupBox("NeXus file with Unified-fit (optional)")
        v = QVBoxLayout(gb)

        file_row = QHBoxLayout()
        self.nexus_path_edit = QLineEdit()
        self.nexus_path_edit.setReadOnly(True)
        self.nexus_path_edit.setPlaceholderText("(no file loaded)")
        file_row.addWidget(self.nexus_path_edit)
        btn = QPushButton("Open NeXus…")
        btn.setFixedWidth(90)
        btn.setToolTip("Load a NeXus file with Unified-fit results.\n"
                       "Detects two consecutive levels where the higher level's\n"
                       "RgCutoff matches the lower level's Rg, and uses them as\n"
                       "comparison targets for grown fractals.")
        btn.clicked.connect(self._on_open_nexus)
        file_row.addWidget(btn)
        v.addLayout(file_row)

        self.targets_label = QLabel("No targets — load a NeXus file or use the optimizer.")
        self.targets_label.setStyleSheet(
            "color:#7f8c8d;font-size:10pt;font-style:italic;padding:4px;")
        self.targets_label.setWordWrap(True)
        v.addWidget(self.targets_label)
        return gb

    def _build_growth_box(self) -> QGroupBox:
        gb = QGroupBox("Growth parameters")
        g = QGridLayout(gb)
        g.setVerticalSpacing(4)

        # Z
        g.addWidget(QLabel("Degree of aggregation Z:"), 0, 0)
        self.z_spin = QSpinBox()
        self.z_spin.setRange(10, 10000)
        self.z_spin.setValue(250)
        self.z_spin.setToolTip("Number of primary particles in the aggregate.\n"
                                "Larger Z → more reliable fractal parameters but slower growth.")
        g.addWidget(self.z_spin, 0, 1)

        # Sticking probability
        g.addWidget(QLabel("Sticking probability [%]:"), 1, 0)
        self.sp_spin = QDoubleSpinBox()
        self.sp_spin.setRange(1.0, 100.0)
        self.sp_spin.setDecimals(1)
        self.sp_spin.setValue(75.0)
        self.sp_spin.setToolTip("Probability (per contact event) that the random-walking particle\n"
                                 "sticks when it touches the existing aggregate.\n"
                                 "Low SP → compact, high df. High SP → loose, low df.")
        g.addWidget(self.sp_spin, 1, 1)

        # Number of test paths
        g.addWidget(QLabel("Number of test paths:"), 2, 0)
        self.ntp_spin = QSpinBox()
        self.ntp_spin.setRange(100, 100000)
        self.ntp_spin.setSingleStep(100)
        self.ntp_spin.setValue(2500)
        self.ntp_spin.setToolTip("Maximum unique paths enumerated per endpoint when computing\n"
                                  "fractal parameters. Affects c (and therefore dmin) statistics.")
        g.addWidget(self.ntp_spin, 2, 1)

        # Rg primary
        g.addWidget(QLabel("Rg primary [Å]:"), 3, 0)
        self.rg_spin = QDoubleSpinBox()
        self.rg_spin.setRange(0.1, 10000.0)
        self.rg_spin.setDecimals(3)
        self.rg_spin.setValue(10.0)
        self.rg_spin.setToolTip("Radius of gyration of the primary sphere [Å].\n"
                                 "Sets the physical scale: primary diameter = 2·√(5/3)·Rg ≈ 2.58 Rg.")
        g.addWidget(self.rg_spin, 3, 1)

        # Allowed neighbor distance
        g.addWidget(QLabel("Allowed neighbor distance:"), 4, 0)
        self.near_combo = QComboBox()
        for label, value in _NEAR_DIST_OPTIONS:
            self.near_combo.addItem(label, value)
        self.near_combo.setCurrentIndex(2)   # body diagonal default
        self.near_combo.setToolTip(
            "Which lattice neighbors count as 'in contact' for sticking:\n"
            "  Edge: only ±1 in x/y/z (6 neighbors)\n"
            "  Face: also include 2D diagonals (18 neighbors)\n"
            "  Body: also include 3D diagonals (26 neighbors)\n"
            "Body-diagonal is the most common Irena setting.")
        g.addWidget(self.near_combo, 4, 1)

        # Multi-particle attraction
        g.addWidget(QLabel("Multi-particle attraction:"), 5, 0)
        self.attr_combo = QComboBox()
        for opt in _ATTRACTION_OPTIONS:
            self.attr_combo.addItem(opt)
        self.attr_combo.setToolTip(
            "How sticking probability changes when ≥2 existing particles are within reach:\n"
            "  Neutral: same SP as for 1 contact\n"
            "  Attractive: SP boosted toward 100% (favors compaction)\n"
            "  Repulsive: SP suppressed (favors openness)\n"
            "  Not allowed: SP = 0 (only chains, no junctions)")
        g.addWidget(self.attr_combo, 5, 1)

        # Random seed
        g.addWidget(QLabel("Random seed (0 = random):"), 6, 0)
        self.seed_spin = QSpinBox()
        self.seed_spin.setRange(0, 2_000_000_000)
        self.seed_spin.setValue(0)
        self.seed_spin.setToolTip("Set non-zero for reproducible aggregates. 0 = fully random.")
        g.addWidget(self.seed_spin, 6, 1)
        return gb

    def _build_mode_tabs(self) -> QTabWidget:
        tabs = QTabWidget()

        # ── Grow One ─────────────────────────────────────────────────────
        one_tab = QWidget()
        one_lay = QVBoxLayout(one_tab)
        self.btn_grow = QPushButton("Grow")
        self.btn_grow.setStyleSheet(
            "background:#27ae60;color:white;font-weight:bold;"
            "font-size:12pt;padding:8px;border-radius:4px;border:none;")
        self.btn_grow.setToolTip("Grow one aggregate with the parameters above.\n"
                                  "Runs in a background thread; the panel stays responsive.")
        self.btn_grow.clicked.connect(self._on_grow_one)
        one_lay.addWidget(self.btn_grow)
        one_lay.addStretch()
        tabs.addTab(one_tab, "Grow One")

        # ── Grow Many ────────────────────────────────────────────────────
        many_tab = QWidget()
        many_lay = QGridLayout(many_tab)
        many_lay.addWidget(QLabel("N aggregates:"), 0, 0)
        self.many_n_spin = QSpinBox()
        self.many_n_spin.setRange(1, 100)
        self.many_n_spin.setValue(5)
        self.many_n_spin.setToolTip("Queue this many growths with the same parameters\n"
                                     "(seed varies per run). Useful for assessing variance.")
        many_lay.addWidget(self.many_n_spin, 0, 1)
        self.btn_grow_many = QPushButton("Grow N")
        self.btn_grow_many.setStyleSheet(
            "background:#27ae60;color:white;font-weight:bold;"
            "font-size:12pt;padding:8px;border-radius:4px;border:none;")
        self.btn_grow_many.clicked.connect(self._on_grow_many)
        many_lay.addWidget(self.btn_grow_many, 1, 0, 1, 2)
        many_lay.setRowStretch(2, 1)
        tabs.addTab(many_tab, "Grow Many")

        # ── Optimizer ────────────────────────────────────────────────────
        opt_tab = QWidget()
        opt_lay = QGridLayout(opt_tab)
        opt_lay.addWidget(QLabel("Target dmin:"), 0, 0)
        self.opt_dmin_spin = QDoubleSpinBox()
        self.opt_dmin_spin.setRange(1.0, 3.0)
        self.opt_dmin_spin.setDecimals(3)
        self.opt_dmin_spin.setValue(2.0)
        self.opt_dmin_spin.setToolTip("Target value of dmin for the bisection search.")
        opt_lay.addWidget(self.opt_dmin_spin, 0, 1)
        opt_lay.addWidget(QLabel("Target c:"), 1, 0)
        self.opt_c_spin = QDoubleSpinBox()
        self.opt_c_spin.setRange(1.0, 3.0)
        self.opt_c_spin.setDecimals(3)
        self.opt_c_spin.setValue(1.2)
        self.opt_c_spin.setToolTip("Target value of c for the bisection search.")
        opt_lay.addWidget(self.opt_c_spin, 1, 1)
        opt_lay.addWidget(QLabel("Tolerance:"), 2, 0)
        self.opt_tol_spin = QDoubleSpinBox()
        self.opt_tol_spin.setRange(0.001, 1.0)
        self.opt_tol_spin.setDecimals(3)
        self.opt_tol_spin.setValue(0.05)
        self.opt_tol_spin.setToolTip("Stop when sqrt(objective) drops below this value.")
        opt_lay.addWidget(self.opt_tol_spin, 2, 1)
        opt_lay.addWidget(QLabel("Max iterations:"), 3, 0)
        self.opt_iter_spin = QSpinBox()
        self.opt_iter_spin.setRange(1, 50)
        self.opt_iter_spin.setValue(10)
        self.opt_iter_spin.setToolTip("Each iteration grows up to 3 aggregates (low/mid/high SP).")
        opt_lay.addWidget(self.opt_iter_spin, 3, 1)
        self.btn_optimize = QPushButton("Find Best Growth")
        self.btn_optimize.setStyleSheet(
            "background:#16a085;color:white;font-weight:bold;"
            "font-size:12pt;padding:8px;border-radius:4px;border:none;")
        self.btn_optimize.clicked.connect(self._on_optimize)
        opt_lay.addWidget(self.btn_optimize, 4, 0, 1, 2)
        opt_lay.setRowStretch(5, 1)
        tabs.addTab(opt_tab, "Optimizer")

        return tabs

    def _build_jobs_box(self) -> QGroupBox:
        gb = QGroupBox("Active jobs")
        v = QVBoxLayout(gb)
        self.jobs_list = QListWidget()
        self.jobs_list.setMaximumHeight(80)
        self.jobs_list.setStyleSheet("font-size:9pt;font-family:monospace;")
        self.jobs_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        v.addWidget(self.jobs_list)
        cancel_row = QHBoxLayout()
        self.btn_cancel_job = QPushButton("Cancel selected job")
        self.btn_cancel_job.setStyleSheet(
            "background:#e74c3c;color:white;font-weight:bold;"
            "padding:3px 8px;border-radius:3px;border:none;")
        self.btn_cancel_job.clicked.connect(self._on_cancel_job)
        cancel_row.addWidget(self.btn_cancel_job)
        cancel_row.addStretch()
        v.addLayout(cancel_row)
        return gb

    def _build_aggregates_list(self) -> QGroupBox:
        gb = QGroupBox("Stored aggregates (session)")
        v = QVBoxLayout(gb)
        self.agg_list = QListWidget()
        self.agg_list.setMinimumHeight(120)
        self.agg_list.setStyleSheet("font-size:9pt;font-family:monospace;")
        self.agg_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.agg_list.itemSelectionChanged.connect(self._on_select_aggregate)
        self.agg_list.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.agg_list.customContextMenuRequested.connect(self._show_agg_context_menu)
        v.addWidget(self.agg_list)
        return gb

    def _build_params_box(self) -> QGroupBox:
        gb = QGroupBox("Active aggregate parameters")
        outer = QVBoxLayout(gb)

        # Target summary that mirrors targets_label at top of panel — duplicated
        # here so the user can compare actual vs target without scrolling.
        self.targets_summary_label = QLabel("Targets: —")
        self.targets_summary_label.setWordWrap(True)
        self.targets_summary_label.setStyleSheet(
            "color:#7f8c8d;font-size:10pt;font-style:italic;padding:2px;")
        outer.addWidget(self.targets_summary_label)

        g = QGridLayout()
        g.setVerticalSpacing(2)
        outer.addLayout(g)
        self._param_labels: dict[str, QLabel] = {}

        rows = [
            ("Z", "Z"),
            ("dmin", "dmin"),
            ("c", "c"),
            ("df", "df"),
            ("R (lattice)", "R"),
            ("p", "p"),
            ("s", "s"),
            ("True sticking [%]", "true_sp"),
            ("Rg primary [Å]", "rg_primary"),
            ("Rg aggregate [Å]", "rg_aggregate"),
            ("Primary diameter [Å]", "primary_d"),
            ("# endpoints", "n_ends"),
        ]
        for i, (label, key) in enumerate(rows):
            r, c = divmod(i, 2)
            g.addWidget(QLabel(label + ":"), r, c * 2)
            lbl = QLabel("—")
            lbl.setStyleSheet("font-family:monospace;font-weight:bold;")
            g.addWidget(lbl, r, c * 2 + 1)
            self._param_labels[key] = lbl
        return gb

    def _build_intensity_box(self) -> QGroupBox:
        gb = QGroupBox("Back-calculated intensity I(Q)")
        g = QGridLayout(gb)
        g.addWidget(QLabel("Q min [Å⁻¹]:"), 0, 0)
        self.qmin_edit = QLineEdit("0.001")
        g.addWidget(self.qmin_edit, 0, 1)
        g.addWidget(QLabel("Q max [Å⁻¹]:"), 0, 2)
        self.qmax_edit = QLineEdit("1.0")
        g.addWidget(self.qmax_edit, 0, 3)
        g.addWidget(QLabel("# points:"), 1, 0)
        self.qn_spin = QSpinBox()
        self.qn_spin.setRange(20, 20000)
        self.qn_spin.setValue(200)
        self.qn_spin.setToolTip(
            "Number of Q-points in the model intensity grid.\n"
            "Both the analytical Unified curve and the Monte-Carlo curve\n"
            "are evaluated on this grid.  Computation is fast even at 20000;\n"
            "a finer grid mainly helps see Porod-region structure."
        )
        g.addWidget(self.qn_spin, 1, 1)

        self.btn_mc = QPushButton("Compute Monte Carlo I(Q)")
        self.btn_mc.setStyleSheet(
            "background:#16a085;color:white;font-weight:bold;"
            "padding:5px 8px;border-radius:3px;border:none;")
        self.btn_mc.setToolTip("Slow PDF-based intensity from the 3D voxelized model.\n"
                                "Uses random pair sampling over solid voxels.\n"
                                "Run this when you want to cross-check the analytical Unified curve.")
        self.btn_mc.clicked.connect(self._on_compute_mc)
        g.addWidget(self.btn_mc, 2, 0, 1, 3)

        self.btn_mc_cancel = QPushButton("Cancel MC")
        self.btn_mc_cancel.setStyleSheet(
            "background:#e74c3c;color:white;font-weight:bold;"
            "padding:5px 8px;border-radius:3px;border:none;")
        self.btn_mc_cancel.clicked.connect(self._on_cancel_mc)
        self.btn_mc_cancel.setVisible(False)
        g.addWidget(self.btn_mc_cancel, 2, 3)

        self.mc_status = QLabel("")
        self.mc_status.setStyleSheet("color:#7f8c8d;font-size:9pt;")
        g.addWidget(self.mc_status, 3, 0, 1, 4)
        return gb

    def _build_save_box(self) -> QGroupBox:
        gb = QGroupBox("Save to NeXus")
        v = QVBoxLayout(gb)
        self.btn_save = QPushButton("Save selected aggregate to NeXus…")
        self.btn_save.setStyleSheet(
            "background:#2980b9;color:white;font-weight:bold;"
            "font-size:11pt;padding:6px;border-radius:4px;border:none;")
        self.btn_save.setToolTip(
            "Append the selected aggregate as a new aggregate_{N} group inside\n"
            "an HDF5/NeXus file.  Existing groups (other aggregates, Unified fit,\n"
            "raw data) are preserved.")
        self.btn_save.clicked.connect(self._on_save)
        v.addWidget(self.btn_save)

        load_row = QHBoxLayout()
        self.btn_load_agg = QPushButton("Load aggregate(s) from NeXus…")
        self.btn_load_agg.setStyleSheet(
            "background:#7f8c8d;color:white;font-weight:bold;"
            "padding:5px 8px;border-radius:3px;border:none;")
        self.btn_load_agg.clicked.connect(self._on_load_aggregates)
        load_row.addWidget(self.btn_load_agg)
        v.addLayout(load_row)
        return gb

    # ── Right side: I(Q) plot + 2D + 3D viewers ──────────────────────────

    def _build_right_graph_area(self) -> QWidget:
        right = QWidget()
        rlay = QVBoxLayout(right)
        rlay.setContentsMargins(0, 0, 0, 0)

        v_split = QSplitter(Qt.Orientation.Vertical)

        # I(Q) plot
        gl = pg.GraphicsLayoutWidget()
        gl.setBackground("w")
        self.iq_plot = make_sas_plot(
            gl, row=0, col=0,
            title="Fractals — I(Q)",
            parent_widget=self,
            jpeg_default_name="fractals_iq.jpg",
        )
        # Legend: keep a reference so we can clear and rebuild it on every
        # _refresh_plot() (otherwise stale entries from prior aggregates
        # accumulate as you click through the list).
        self._legend = self.iq_plot.addLegend(
            offset=(-10, 10),
            labelTextColor=SASPlotStyle.LEGEND_TEXT_COLOR,
        )
        v_split.addWidget(gl)

        # 2D + 3D viewers in horizontal splitter
        bottom = QWidget()
        h_lay = QHBoxLayout(bottom)
        h_lay.setContentsMargins(4, 4, 4, 4)
        h_lay.setSpacing(6)

        slice_box = QWidget()
        slice_col = QVBoxLayout(slice_box)
        slice_col.setContentsMargins(0, 0, 0, 0)
        self.slice_viewer = Slice2DViewer()
        slice_col.addWidget(self.slice_viewer)
        slice_col.addWidget(make_popout_button(self.slice_viewer, "2D slice viewer"))
        h_lay.addWidget(slice_box, 1)

        voxel_box = QWidget()
        voxel_col = QVBoxLayout(voxel_box)
        voxel_col.setContentsMargins(0, 0, 0, 0)
        self.voxel3d_viewer = Voxel3DViewer()
        voxel_col.addWidget(self.voxel3d_viewer)
        voxel_col.addWidget(make_popout_button(self.voxel3d_viewer, "3D viewer"))
        h_lay.addWidget(voxel_box, 1)

        v_split.addWidget(bottom)
        v_split.setSizes([400, 600])
        rlay.addWidget(v_split)
        return right

    # ── Worker wiring ────────────────────────────────────────────────────

    def _wire_workers(self):
        self._growth_worker.job_started.connect(self._on_job_started)
        self._growth_worker.job_progress.connect(self._on_job_progress)
        self._growth_worker.aggregate_complete.connect(self._on_aggregate_complete)
        self._growth_worker.optimizer_complete.connect(self._on_optimizer_complete)
        self._growth_worker.job_failed.connect(self._on_job_failed)

        self._mc_worker.progress.connect(self._on_mc_progress)
        self._mc_worker.finished_ok.connect(self._on_mc_finished)
        self._mc_worker.failed.connect(self._on_mc_failed)

    # ── State load/save ──────────────────────────────────────────────────

    def _current_growth_config(self) -> GrowthConfig:
        return GrowthConfig(
            z=int(self.z_spin.value()),
            sticking_prob=float(self.sp_spin.value()),
            num_test_paths=int(self.ntp_spin.value()),
            rg_primary=float(self.rg_spin.value()),
            allowed_near_dist=int(self.near_combo.currentData()),
            attraction=str(self.attr_combo.currentText()),
            seed=int(self.seed_spin.value()),
        )

    def _current_optimizer_config(self) -> OptimizerConfig:
        gc = self._current_growth_config()
        return OptimizerConfig(
            target_dmin=float(self.opt_dmin_spin.value()),
            target_c=float(self.opt_c_spin.value()),
            tolerance=float(self.opt_tol_spin.value()),
            max_iter=int(self.opt_iter_spin.value()),
            z=gc.z, num_test_paths=gc.num_test_paths,
            rg_primary=gc.rg_primary, allowed_near_dist=gc.allowed_near_dist,
            attraction=gc.attraction, seed=gc.seed,
        )

    def _load_state(self):
        st = self._state.get("fractals") or {}
        gp = st.get("last_growth_params") or {}
        if "z" in gp:
            self.z_spin.setValue(int(gp["z"]))
        if "sticking_prob" in gp:
            self.sp_spin.setValue(float(gp["sticking_prob"]))
        if "num_test_paths" in gp:
            self.ntp_spin.setValue(int(gp["num_test_paths"]))
        if "rg_primary" in gp:
            self.rg_spin.setValue(float(gp["rg_primary"]))
        if "allowed_near_dist" in gp:
            for i in range(self.near_combo.count()):
                if self.near_combo.itemData(i) == int(gp["allowed_near_dist"]):
                    self.near_combo.setCurrentIndex(i)
                    break
        if "attraction" in gp:
            idx = self.attr_combo.findText(str(gp["attraction"]))
            if idx >= 0:
                self.attr_combo.setCurrentIndex(idx)
        if "seed" in gp:
            self.seed_spin.setValue(int(gp["seed"]))

        qr = st.get("q_range") or {}
        if "q_min" in qr:
            self.qmin_edit.setText(_fmt(qr["q_min"]))
        if "q_max" in qr:
            self.qmax_edit.setText(_fmt(qr["q_max"]))
        if "n_points" in qr:
            self.qn_spin.setValue(int(qr["n_points"]))

        opt = st.get("optimizer") or {}
        if "target_dmin" in opt:
            self.opt_dmin_spin.setValue(float(opt["target_dmin"]))
        if "target_c" in opt:
            self.opt_c_spin.setValue(float(opt["target_c"]))
        if "tolerance" in opt:
            self.opt_tol_spin.setValue(float(opt["tolerance"]))
        if "max_iter" in opt:
            self.opt_iter_spin.setValue(int(opt["max_iter"]))

        self.many_n_spin.setValue(int(st.get("grow_many_n", 5)))

        # Auto-restore last NeXus file: if it still exists on disk, load
        # everything (data + Unified-fit + targets) so the user does not
        # have to re-pick the file every session.  If the file has been
        # moved or deleted, silently clear the path so the field is empty.
        last_path = st.get("last_loaded_nexus_path", "")
        if last_path:
            nexus_path = Path(last_path)
            if nexus_path.exists():
                self.nexus_path_edit.setText(last_path)
                try:
                    self._load_unified_from_nexus(nexus_path)
                except Exception as exc:
                    log.warning("[fractals] auto-load of %s failed: %s", last_path, exc)
            else:
                # File no longer exists — drop the stale path
                self.nexus_path_edit.clear()

    def _save_state(self):
        st = {
            "schema_version": 1,
            "last_growth_params": {
                "z": int(self.z_spin.value()),
                "sticking_prob": float(self.sp_spin.value()),
                "num_test_paths": int(self.ntp_spin.value()),
                "rg_primary": float(self.rg_spin.value()),
                "allowed_near_dist": int(self.near_combo.currentData()),
                "attraction": str(self.attr_combo.currentText()),
                "seed": int(self.seed_spin.value()),
            },
            "q_range": {
                "q_min": float(self.qmin_edit.text() or "0.001"),
                "q_max": float(self.qmax_edit.text() or "1.0"),
                "n_points": int(self.qn_spin.value()),
            },
            "optimizer": {
                "target_dmin": float(self.opt_dmin_spin.value()),
                "target_c": float(self.opt_c_spin.value()),
                "tolerance": float(self.opt_tol_spin.value()),
                "max_iter": int(self.opt_iter_spin.value()),
            },
            "grow_many_n": int(self.many_n_spin.value()),
            "last_loaded_nexus_path": str(self.nexus_path_edit.text()),
        }
        self._state.update("fractals", st)
        self._state.save()

    # ── Help / NeXus loading ─────────────────────────────────────────────

    def _open_help(self):
        try:
            from PySide6.QtGui import QDesktopServices
            from PySide6.QtCore import QUrl
        except ImportError:
            from PyQt6.QtGui import QDesktopServices
            from PyQt6.QtCore import QUrl
        QDesktopServices.openUrl(QUrl(
            "https://github.com/jilavsky/pyirena/blob/main/docs/fractals_gui.md"
        ))

    def _on_open_nexus(self):
        last = self.nexus_path_edit.text() or str(Path.home())
        path, _ = QFileDialog.getOpenFileName(
            self, "Open NeXus file with Unified-fit results", last,
            "HDF5 files (*.h5 *.hdf5 *.nxs);;All files (*)",
        )
        if not path:
            return
        self.nexus_path_edit.setText(path)
        self._load_unified_from_nexus(Path(path))
        self._save_state()

    def _load_unified_from_nexus(self, filepath: Path):
        # Try to load Unified fit
        try:
            uf = load_unified_fit_results(filepath)
        except Exception as exc:
            QMessageBox.warning(
                self, "No Unified fit",
                f"Could not load Unified-fit results from\n{filepath}:\n{exc}"
            )
            self._loaded_unified_q = None
            self._loaded_unified_I = None
            self._targets = {}
            self._update_targets_label()
            self._refresh_plot()
            return

        self._loaded_unified_q = np.asarray(uf.get("Q"), dtype=float)
        self._loaded_unified_I = np.asarray(uf.get("intensity_model"), dtype=float)

        # Try to load primary I(Q) from same file (NXcanSAS)
        try:
            from pyirena.io.hdf5 import readGenericNXcanSAS
            data = readGenericNXcanSAS(str(filepath.parent), filepath.name)
            self._loaded_q = np.asarray(data.get("Q"), dtype=float)
            self._loaded_I = np.asarray(data.get("Intensity"), dtype=float)
            err = data.get("Error")
            self._loaded_dI = np.asarray(err, dtype=float) if err is not None else None
        except Exception:
            # Use the Unified-fit's stored experimental data as a fallback
            self._loaded_q = np.asarray(uf.get("Q"), dtype=float)
            self._loaded_I = np.asarray(uf.get("intensity_data"), dtype=float)
            err = uf.get("intensity_error")
            self._loaded_dI = np.asarray(err, dtype=float) if err is not None else None

        # Detect a fractal-pair (consecutive levels with RgCutoff_high ≈ Rg_low)
        levels = uf.get("levels") or []
        self._targets = self._derive_targets_from_levels(levels)
        self._update_targets_label()

        # Auto-set Q range fields from data
        if self._loaded_q is not None and self._loaded_q.size:
            qmin = float(np.min(self._loaded_q[self._loaded_q > 0]))
            qmax = float(np.max(self._loaded_q[self._loaded_q > 0]))
            self.qmin_edit.setText(_fmt(qmin))
            self.qmax_edit.setText(_fmt(qmax))

        # If a target Rg is detected, push it to the Rg-primary spinbox as a hint
        if "rg_primary" in self._targets:
            self.rg_spin.setValue(float(self._targets["rg_primary"]))

        self._refresh_plot()

    @staticmethod
    def _derive_targets_from_levels(levels: list[dict]) -> dict:
        """Find the first consecutive (low, high) level pair where the high
        level's RgCutoff ≈ the low level's Rg (within ±25 %), then derive
        the target fractal parameters.

        Igor (IR3A_*) formulas:
          df    = P2
          dmin  = B2 · Rg2^P2 / ( Γ(P2/2) · G2 )
          c     = P2 / dmin
          Z     = G2 / G1 + 1            (degree of aggregation)
          fBr   = 1 − (G2/G1)^(1/(c-1))  (branched fraction, informational)
          fM    = 1 − (G2/G1)^(1/(dmin-1))  (mass fraction, informational)

        Returns a dict with keys:
          level_low, level_high, rg_primary, rg_aggregate,
          df, z_target, dmin_target, c_target, fBr, fM
        Empty dict if no qualifying pair.
        """
        if not levels or len(levels) < 2:
            return {}
        for i in range(len(levels) - 1):
            lo = levels[i]
            hi = levels[i + 1]
            rg_lo = float(lo.get("Rg", 0.0))
            rg_hi = float(hi.get("Rg", 0.0))
            cutoff_hi = float(hi.get("RgCutoff", 0.0))
            if rg_lo <= 0 or rg_hi <= 0 or cutoff_hi <= 0:
                continue
            if abs(cutoff_hi - rg_lo) / rg_lo > 0.25:
                continue
            G1 = float(lo.get("G", 0.0))
            G2 = float(hi.get("G", 0.0))
            B2 = float(hi.get("B", 0.0))
            P2 = float(hi.get("P", 0.0))
            if not (G1 > 0 and G2 > 0 and B2 > 0 and P2 > 0):
                continue
            try:
                gamma_term = math.exp(math.lgamma(P2 / 2.0))
                dmin = (B2 * (rg_hi ** P2)) / (gamma_term * G2)
                if not math.isfinite(dmin) or dmin <= 0:
                    continue
                c = P2 / dmin
                z_target = G2 / G1 + 1.0
                df = P2
                # Auxiliary fractions (informational; NaN-safe if exponent invalid)
                ratio = G2 / G1
                try:
                    fBr = 1.0 - ratio ** (1.0 / (c - 1.0)) if c != 1.0 else float("nan")
                except (OverflowError, ValueError):
                    fBr = float("nan")
                try:
                    fM = 1.0 - ratio ** (1.0 / (dmin - 1.0)) if dmin != 1.0 else float("nan")
                except (OverflowError, ValueError):
                    fM = float("nan")
                return {
                    "level_low": i + 1, "level_high": i + 2,
                    "rg_primary": rg_lo, "rg_aggregate": rg_hi,
                    "df": df, "z_target": z_target,
                    "dmin_target": dmin, "c_target": c,
                    "fBr": fBr, "fM": fM,
                }
            except (OverflowError, ValueError, ZeroDivisionError):
                continue
        return {}

    def _update_targets_label(self):
        if not self._targets:
            empty = ("No targets — load a NeXus file with Unified-fit where two "
                     "consecutive levels have RgCutoff_high ≈ Rg_low, or set "
                     "manual targets in the Optimizer tab.")
            self.targets_label.setText(empty)
            self.targets_label.setStyleSheet(
                "color:#7f8c8d;font-size:10pt;font-style:italic;padding:4px;")
            if hasattr(self, 'targets_summary_label'):
                self.targets_summary_label.setText("Targets: —")
                self.targets_summary_label.setStyleSheet(
                    "color:#7f8c8d;font-size:10pt;font-style:italic;padding:2px;")
            return
        t = self._targets
        # Sanity check on c (connectivity dimension): physically c MUST be
        # ≥ 1 for any branched / fractal aggregate (chains have c=1, more
        # connected branches have c>1).  Igor's code carries the same
        # warning — when the Unified-fit numbers give c<1 the data simply
        # cannot be modelled by a mass-fractal aggregate (Z, dmin, df all
        # become non-physical too).  Show in RED so the user notices.
        c = t.get('c_target')
        c_invalid = (
            isinstance(c, (int, float))
            and math.isfinite(float(c))
            and float(c) < 1.0
        )
        # Top label (full detail) — visible at top of panel
        if c_invalid:
            msg = (
                f"<b>⚠ Detected fractal pair (levels {t['level_low']} + {t['level_high']}) "
                f"is NOT a valid fractal:</b><br>"
                f"connectivity dimension c = {_fmt(c)} &lt; 1.<br>"
                f"This data cannot be modelled by a mass-fractal aggregate. "
                f"Check the Unified-fit input or pick a different level pair.<br>"
                f"<small>(c=1: linear chain; c≈1.1–1.4: branched fractal.  "
                f"c&lt;1 implies non-physical Z, dmin, df.)</small>"
            )
            self.targets_label.setText(msg)
            self.targets_label.setStyleSheet(
                "color:#c0392b;font-size:10pt;padding:6px;"
                "background-color:#fdecea;border:1px solid #c0392b;border-radius:4px;"
            )
        else:
            msg = (
                f"<b>Detected fractal pair: levels {t['level_low']} + {t['level_high']}</b><br>"
                f"<b>Z = {_fmt(t.get('z_target'))}</b>, "
                f"<b>dmin = {_fmt(t.get('dmin_target'))}</b>, "
                f"<b>c = {_fmt(t.get('c_target'))}</b>, "
                f"<b>df = {_fmt(t.get('df'))}</b><br>"
                f"Rg primary = {_fmt(t['rg_primary'])} Å, "
                f"Rg aggregate = {_fmt(t['rg_aggregate'])} Å"
            )
            self.targets_label.setText(msg)
            self.targets_label.setStyleSheet(
                "color:#16a085;font-size:10pt;padding:4px;"
            )
        # Compact summary near "Active aggregate parameters" widget
        if hasattr(self, 'targets_summary_label'):
            if c_invalid:
                short = (
                    f"<b>⚠ Invalid fractal: c = {_fmt(c)} &lt; 1 — "
                    f"this Unified-fit data is not modellable as a "
                    f"mass-fractal aggregate.</b>"
                )
                self.targets_summary_label.setText(short)
                self.targets_summary_label.setStyleSheet(
                    "color:#c0392b;font-size:10pt;padding:2px;"
                )
            else:
                short = (
                    f"<b>Targets — Z={_fmt(t.get('z_target'))}, "
                    f"dmin={_fmt(t.get('dmin_target'))}, "
                    f"c={_fmt(t.get('c_target'))}, "
                    f"df={_fmt(t.get('df'))}, "
                    f"Rg_p={_fmt(t['rg_primary'])} Å, "
                    f"Rg_a={_fmt(t['rg_aggregate'])} Å</b>"
                )
                self.targets_summary_label.setText(short)
                self.targets_summary_label.setStyleSheet(
                    "color:#16a085;font-size:10pt;padding:2px;"
                )

    # ── Action handlers ──────────────────────────────────────────────────

    def _on_grow_one(self):
        cfg = self._current_growth_config()
        self._growth_worker.enqueue_grow(cfg, label=f"Grow Z={cfg.z} SP={cfg.sticking_prob:.0f}%")
        self._save_state()

    def _on_grow_many(self):
        cfg = self._current_growth_config()
        n = int(self.many_n_spin.value())
        self._growth_worker.enqueue_grow_many(cfg, n)
        self._save_state()

    def _on_optimize(self):
        opt_cfg = self._current_optimizer_config()
        self._growth_worker.enqueue_optimize(opt_cfg)
        # Also push the targets into the comparison panel
        self._targets = {
            **self._targets,
            "dmin_target": float(opt_cfg.target_dmin),
            "c_target": float(opt_cfg.target_c),
        }
        self._save_state()

    def _on_cancel_job(self):
        items = self.jobs_list.selectedItems()
        if not items:
            return
        job_id = items[0].data(Qt.ItemDataRole.UserRole)
        self._growth_worker.request_cancel(job_id)

    # ── Worker callbacks ─────────────────────────────────────────────────

    def _on_job_started(self, job_id: str, label: str):
        item = QListWidgetItem(f"⏳ [{job_id}] {label} — 0%")
        item.setData(Qt.ItemDataRole.UserRole, job_id)
        self.jobs_list.addItem(item)
        self._jobs[job_id] = {"label": label, "list_item": item}

    def _on_job_progress(self, job_id: str, percent: int, message: str):
        info = self._jobs.get(job_id)
        if not info:
            return
        info["list_item"].setText(f"⏳ [{job_id}] {info['label']} — {percent}% ({message})")

    def _on_aggregate_complete(self, job_id: str, agg: FractalAggregate):
        info = self._jobs.get(job_id, {})
        # Auto-compute Unified intensity
        try:
            q = self._make_q_grid()
            agg.q = q
            agg.i_unified = intensity_unified(agg.params, q)
        except Exception:
            log.debug("suppressed exception", exc_info=True)
        self._aggregates.append(agg)
        self._refresh_aggregates_list()
        # Auto-select the newest one
        self.agg_list.setCurrentRow(self.agg_list.count() - 1)

        if info.get("list_item") is not None:
            info["list_item"].setText(
                f"✓ [{job_id}] {info.get('label', '')} — done"
            )
            # Remove from jobs after a short visible moment
            self._remove_job_after_delay(job_id, 1500)

    def _remove_job_after_delay(self, job_id: str, ms: int = 1500):
        try:
            from PySide6.QtCore import QTimer
        except ImportError:
            from PyQt6.QtCore import QTimer
        QTimer.singleShot(ms, lambda: self._remove_job(job_id))

    def _remove_job(self, job_id: str):
        info = self._jobs.pop(job_id, None)
        if not info:
            return
        item = info.get("list_item")
        if item is not None:
            row = self.jobs_list.row(item)
            if row >= 0:
                self.jobs_list.takeItem(row)

    def _on_optimizer_complete(self, job_id: str, best: FractalAggregate, all_aggs: list):
        # Mark the best one with an 'best' suffix in label
        best.label = (best.label or "") + " ★ BEST"
        # The per-iteration aggregates have already been appended via aggregate_complete
        info = self._jobs.get(job_id, {})
        if info.get("list_item") is not None:
            info["list_item"].setText(
                f"✓ [{job_id}] Optimizer done — best dmin={_fmt(best.params.dmin)}, c={_fmt(best.params.c)}"
            )
            self._remove_job_after_delay(job_id, 2500)
        # Highlight the best aggregate in the list
        for i, a in enumerate(self._aggregates):
            if a is best:
                self.agg_list.setCurrentRow(i)
                break
        self._refresh_aggregates_list()

    def _on_job_failed(self, job_id: str, message: str):
        info = self._jobs.get(job_id, {})
        if info.get("list_item") is not None:
            info["list_item"].setText(f"✗ [{job_id}] FAILED — {message}")
        QMessageBox.critical(self, "Job failed",
                              f"Job {job_id} failed:\n{message}")
        self._remove_job_after_delay(job_id, 4000)

    def _on_mc_progress(self, percent: int):
        self.mc_status.setText(f"MC running — {percent}%")

    def _on_mc_finished(self, agg_uuid: str, q: np.ndarray, I_mc: np.ndarray):
        self.btn_mc_cancel.setVisible(False)
        self.mc_status.setText("MC done.")
        # Find the aggregate by uuid and store result.
        # Critical: keep `a.q` and `a.i_unified` in sync.  If the user changed
        # the # points spinbox between Grow and MC, the new MC q-grid has a
        # different length than the stored i_unified array, and any later
        # plot mask `(a.q > 0) & (a.i_unified > 0)` would raise a numpy
        # broadcast error and silently leave the curve off the plot.
        for a in self._aggregates:
            if a.uuid == agg_uuid:
                a.q = q
                a.i_montecarlo = I_mc
                # Re-evaluate the analytical Unified intensity on the same
                # grid so the two arrays always have matching lengths.
                try:
                    a.i_unified = intensity_unified(a.params, q)
                except Exception:
                    a.i_unified = None
                break
        if self._active is not None and self._active.uuid == agg_uuid:
            self._refresh_plot()

    def _on_mc_failed(self, agg_uuid: str, message: str):
        self.btn_mc_cancel.setVisible(False)
        self.mc_status.setText(f"MC failed: {message}")
        QMessageBox.warning(self, "MC failed", message)

    def _on_compute_mc(self):
        if self._active is None:
            QMessageBox.information(self, "No aggregate",
                                     "Select a grown aggregate first.")
            return
        if self._mc_worker.isRunning():
            QMessageBox.information(self, "MC busy",
                                     "A Monte-Carlo computation is already running.")
            return
        try:
            q = self._make_q_grid()
        except Exception as exc:
            QMessageBox.warning(self, "Bad Q range", str(exc))
            return
        # Synchronize the aggregate's q-grid + analytical Unified intensity
        # to the same array we'll send to the MC worker.  This guarantees
        # that when MC returns, all three arrays (q, i_unified, i_mc) have
        # matching shapes — fixes the bug where changing # points caused
        # all curves to vanish from the plot.
        self._active.q = q
        try:
            self._active.i_unified = intensity_unified(self._active.params, q)
        except Exception:
            self._active.i_unified = None
        self._mc_worker.configure(self._active, q)
        self.btn_mc_cancel.setVisible(True)
        self.mc_status.setText("MC starting…")
        self._mc_worker.start()

    def _on_cancel_mc(self):
        self._mc_worker.request_cancel()

    # ── Aggregate selection / display ────────────────────────────────────

    def _refresh_aggregates_list(self):
        # Preserve current selection by uuid
        sel_uuid = self._active.uuid if self._active else None
        self.agg_list.blockSignals(True)
        self.agg_list.clear()
        for i, a in enumerate(self._aggregates):
            params = a.params
            line = (f"[{i+1:02d}] Z={params.z}  df={_fmt(params.df)}  "
                    f"c={_fmt(params.c)}  dmin={_fmt(params.dmin)}  "
                    f"— {a.label}")
            it = QListWidgetItem(line)
            it.setData(Qt.ItemDataRole.UserRole, i)
            if "BEST" in (a.label or ""):
                fnt = it.font(); fnt.setBold(True); it.setFont(fnt)
            self.agg_list.addItem(it)
            if sel_uuid and a.uuid == sel_uuid:
                self.agg_list.setCurrentRow(i)
        self.agg_list.blockSignals(False)

    def _on_select_aggregate(self):
        items = self.agg_list.selectedItems()
        if not items:
            self._active = None
            return
        idx = items[0].data(Qt.ItemDataRole.UserRole)
        if 0 <= idx < len(self._aggregates):
            self._active = self._aggregates[idx]
            self._update_params_display()
            self._refresh_voxel_views()
            self._refresh_plot()

    def _update_params_display(self):
        p = self._active.params if self._active else None
        if p is None:
            for lbl in self._param_labels.values():
                lbl.setText("—")
                lbl.setStyleSheet("font-family:monospace;font-weight:bold;")
            return
        targets = self._targets
        def _color(actual: float, target_key: str):
            t = targets.get(target_key)
            if t is None or not math.isfinite(actual):
                return ""
            try:
                rel = abs(actual - float(t)) / max(abs(float(t)), 1e-12)
            except Exception:
                return ""
            if rel < 0.10:
                return "color:#27ae60;"
            if rel < 0.25:
                return "color:#e67e22;"
            return "color:#c0392b;"

        def _with_target(actual_text: str, target_key: str) -> str:
            t = targets.get(target_key)
            if t is None or not isinstance(t, (int, float)) or not math.isfinite(float(t)):
                return actual_text
            return f"{actual_text}   (target: {_fmt(t)})"

        mapping = [
            ("Z",            _with_target(str(p.z), "z_target"),                ""),
            ("dmin",         _with_target(_fmt(p.dmin), "dmin_target"),         _color(p.dmin, "dmin_target")),
            ("c",            _with_target(_fmt(p.c), "c_target"),               _color(p.c, "c_target")),
            ("df",           _with_target(_fmt(p.df), "df"),                    _color(p.df, "df")),
            ("R",            _fmt(p.R_dimensionless),                           ""),
            ("p",            _fmt(p.p),                                         ""),
            ("s",            _fmt(p.s),                                         ""),
            ("true_sp",      _fmt(p.true_sticking_prob),                        ""),
            ("rg_primary",   _with_target(_fmt(p.rg_primary), "rg_primary"),    _color(p.rg_primary, "rg_primary")),
            ("rg_aggregate", _with_target(_fmt(p.rg_aggregate), "rg_aggregate"), _color(p.rg_aggregate, "rg_aggregate")),
            ("primary_d",    _fmt(p.primary_diameter),                          ""),
            ("n_ends",       str(p.num_endpoints),                              ""),
        ]
        for key, text, css in mapping:
            lbl = self._param_labels.get(key)
            if lbl is not None:
                lbl.setText(text)
                lbl.setStyleSheet("font-family:monospace;font-weight:bold;" + css)

    def _refresh_voxel_views(self):
        if self._active is None:
            return
        try:
            # DISPLAY uses Irena-style "fat spheres": oversample=10 with
            # sphere_voxel_radius=10 so that the sphere RADIUS equals the
            # full lattice spacing D (= primary_diameter).  This guarantees
            # that all neighbor types (edge at d=1, face-diag at d=√2,
            # body-diag at d=√3) overlap visually — producing the familiar
            # connected-blob fractal look.  The sphere is physically 2× too
            # large, but that is acceptable for a qualitative 3D visualis-
            # ation.  The MC scattering calculation (intensity_montecarlo)
            # uses the physically-correct geometry independently
            # (oversample=20, sphere_voxel_radius=10 → sphere radius = R).
            voxelgram, _pitch_lattice = voxelize(
                self._active.positions,
                oversample=10,
                sphere_voxel_radius=10,
            )
            # pitch_A = D / oversample = primary_diameter / 10
            pitch_A = self._active.params.primary_diameter / 10.0
            self.slice_viewer.set_voxelgram(voxelgram, pitch_A)
            self.voxel3d_viewer.set_voxelgram(voxelgram, pitch_A)
        except Exception as exc:
            log.warning("[fractals] voxelization failed: %s", exc)

    # ── I(Q) plotting ────────────────────────────────────────────────────

    def _make_q_grid(self) -> np.ndarray:
        try:
            qmin = float(self.qmin_edit.text())
            qmax = float(self.qmax_edit.text())
        except ValueError:
            raise ValueError("Q min and Q max must be numeric.")
        if qmin <= 0 or qmax <= qmin:
            raise ValueError("Need 0 < Q min < Q max.")
        n = int(self.qn_spin.value())
        return np.logspace(np.log10(qmin), np.log10(qmax), n)

    def _refresh_plot(self):
        # Drop ALL plot items + clear the legend so the next batch of
        # entries replaces (rather than appends to) any previous ones.
        self.iq_plot.clear()
        if self._legend is not None:
            try:
                self._legend.clear()
            except Exception:
                log.debug("suppressed exception", exc_info=True)

        items_for_y: list[float] = []

        # ── Loaded experimental data ─────────────────────────────────────
        if self._loaded_q is not None and self._loaded_I is not None:
            mask = (self._loaded_q > 0) & (self._loaded_I > 0)
            q_d = self._loaded_q[mask]; I_d = self._loaded_I[mask]
            dI = self._loaded_dI[mask] if self._loaded_dI is not None else None
            if q_d.size > 0:
                plot_iq_data(self.iq_plot, q_d, I_d, dI, label="Data")
                items_for_y.extend(I_d.tolist())

        # ── Loaded Unified-fit model ─────────────────────────────────────
        if self._loaded_unified_q is not None and self._loaded_unified_I is not None:
            qm = self._loaded_unified_q
            Im = self._loaded_unified_I
            if qm.shape == Im.shape:
                mask = (qm > 0) & (Im > 0)
                if np.any(mask):
                    plot_iq_model(self.iq_plot, qm[mask], Im[mask],
                                  label="Unified fit (loaded)")

        # ── Active aggregate ─────────────────────────────────────────────
        if self._active is not None:
            self._ensure_active_q_and_unified()
            q_active = self._active.q
            i_uni = self._active.i_unified
            i_mc = self._active.i_montecarlo

            # Defensive shape check — if any earlier code path left arrays
            # of mismatched length, recompute i_unified once more so the
            # mask operation below cannot raise.
            if (q_active is not None and i_uni is not None
                    and q_active.shape != i_uni.shape):
                try:
                    i_uni = intensity_unified(self._active.params, q_active)
                    self._active.i_unified = i_uni
                except Exception:
                    i_uni = None
            if (q_active is not None and i_mc is not None
                    and q_active.shape != i_mc.shape):
                # MC array doesn't match: drop it so we don't crash.
                i_mc = None
                self._active.i_montecarlo = None

            # Aggregate Unified-from-fractals (green) — always shown when available
            if q_active is not None and i_uni is not None:
                I_show = self._invariant_rescale(q_active, i_uni)
                mask = (q_active > 0) & np.isfinite(I_show) & (I_show > 0)
                if np.any(mask):
                    self.iq_plot.plot(
                        q_active[mask], I_show[mask],
                        pen=pg.mkPen("#27ae60", width=2),
                        name="Aggregate Unified (analytical)",
                    )
                    items_for_y.extend(I_show[mask].tolist())

            # Aggregate Monte-Carlo (blue dashed) — only if computed.
            # MC is NaN above the voxel-grid Nyquist Q, so np.isfinite is
            # required (a bare > 0 check would still let NaNs through some
            # downstream paths and confuse set_robust_y_range).
            if q_active is not None and i_mc is not None:
                I_mc_show = self._invariant_rescale(q_active, i_mc)
                mask = (q_active > 0) & np.isfinite(I_mc_show) & (I_mc_show > 0)
                if np.any(mask):
                    self.iq_plot.plot(
                        q_active[mask], I_mc_show[mask],
                        pen=pg.mkPen("#2980b9", width=2,
                                     style=Qt.PenStyle.DashLine),
                        name="Aggregate Monte-Carlo",
                    )
                    items_for_y.extend(I_mc_show[mask].tolist())

        if items_for_y:
            set_robust_y_range(self.iq_plot, np.asarray(items_for_y))

    def _ensure_active_q_and_unified(self):
        """Make sure the active aggregate has a q-grid and matching i_unified.

        Lazily generates them on first plot or after `# points` changes when
        the user re-selects an older aggregate that has no MC computed.
        """
        if self._active is None:
            return
        a = self._active
        try:
            q_target = self._make_q_grid()
        except Exception:
            return
        # If we have nothing yet, or the lengths don't match the desired grid
        # AND no MC has locked in a specific grid, regenerate from the GUI grid.
        needs_new = (
            a.q is None or a.i_unified is None
            or (a.q.shape != q_target.shape and a.i_montecarlo is None)
        )
        if needs_new:
            a.q = q_target
            a.i_unified = intensity_unified(a.params, q_target)

    def _invariant_rescale(self, q: np.ndarray, I: np.ndarray) -> np.ndarray:
        """Scale a model curve to the data using the *fractal-regime* invariant.

        Matches Igor's `IR3A_Calculate1DIntensity`: integrate over
        Q ∈ [0.5π/Rg_aggregate, 1.5π/Rg_primary].  This window covers the
        Q range where the aggregate's mass-fractal scattering dominates,
        and deliberately EXCLUDES:
          * very low Q where a sample-level power-law from larger structures
            inflates the data integral
          * very high Q where flat instrumental background dominates the data
        Without this window the model would always sit too high relative to
        the data when either contamination is present.

        Falls back to the visible-Q overlap when no usable Rg pair exists.
        Returns I unchanged when no data is loaded.
        """
        if self._loaded_q is None or self._loaded_I is None:
            return I
        try:
            q_d = self._loaded_q
            I_d = self._loaded_I
            # Choose the integration window
            qlo, qhi = None, None
            if self._active is not None:
                rg_p = self._active.params.rg_primary
                rg_a = self._active.params.rg_aggregate
                if (math.isfinite(rg_p) and math.isfinite(rg_a)
                        and rg_p > 0 and rg_a > rg_p):
                    qlo = 0.5 * math.pi / rg_a
                    qhi = 1.5 * math.pi / rg_p
            if qlo is None or qhi is None:
                # Fall back to visible-Q overlap
                qlo = max(float(np.min(q_d[q_d > 0])),
                          float(np.min(q[q > 0])))
                qhi = min(float(np.max(q_d[q_d > 0])),
                          float(np.max(q[q > 0])))
            # Always clamp to the union of both Q ranges
            qlo = max(qlo, float(np.min(q_d[q_d > 0])),
                      float(np.min(q[q > 0])))
            qhi = min(qhi, float(np.max(q_d[q_d > 0])),
                      float(np.max(q[q > 0])))
            if qhi <= qlo:
                return I
            mask_d = (q_d >= qlo) & (q_d <= qhi) & (I_d > 0)
            mask_m = (q >= qlo) & (q <= qhi) & (I > 0)
            if not (np.any(mask_d) and np.any(mask_m)):
                return I
            int_d = float(np.trapezoid(I_d[mask_d], q_d[mask_d]))
            int_m = float(np.trapezoid(I[mask_m], q[mask_m]))
            if int_m <= 0:
                return I
            return I * (int_d / int_m)
        except Exception:
            return I

    # ── Save / Load aggregates ───────────────────────────────────────────

    def _on_save(self):
        if self._active is None:
            QMessageBox.information(self, "No aggregate",
                                     "Select an aggregate to save first.")
            return
        last = self.nexus_path_edit.text() or str(Path.home())
        path, _ = QFileDialog.getSaveFileName(
            self, "Save aggregate to NeXus file", last,
            "HDF5 files (*.h5 *.hdf5);;All files (*)",
        )
        if not path:
            return
        # Ensure .h5 extension
        if not path.lower().endswith((".h5", ".hdf5", ".nxs")):
            path = path + ".h5"
        try:
            gpath = save_fractal_aggregate(Path(path), self._active)
        except Exception as exc:
            QMessageBox.critical(self, "Save failed", str(exc))
            return
        QMessageBox.information(self, "Saved",
                                 f"Aggregate written to:\n{path}\n\nGroup: {gpath}")

    def _on_load_aggregates(self):
        last = self.nexus_path_edit.text() or str(Path.home())
        path, _ = QFileDialog.getOpenFileName(
            self, "Load aggregate(s) from NeXus", last,
            "HDF5 files (*.h5 *.hdf5);;All files (*)",
        )
        if not path:
            return
        try:
            entries = list_fractal_aggregates(Path(path))
        except Exception as exc:
            QMessageBox.critical(self, "List failed", str(exc))
            return
        if not entries:
            QMessageBox.information(self, "No aggregates",
                                     "No fractal aggregates found in this file.")
            return
        # Load all of them
        loaded = 0
        for e in entries:
            try:
                agg = load_fractal_aggregate(Path(path), e["group_path"])
                agg.label = (agg.label or e["name"]) + " (loaded)"
                self._aggregates.append(agg)
                loaded += 1
            except Exception as exc:
                log.warning("[fractals] failed to load %s: %s", e['group_path'], exc)
        self._refresh_aggregates_list()
        QMessageBox.information(self, "Loaded",
                                 f"Loaded {loaded} aggregate(s) from {path}.")

    # ── Aggregates list context menu ─────────────────────────────────────

    def _show_agg_context_menu(self, pos):
        item = self.agg_list.itemAt(pos)
        if item is None:
            return
        idx = item.data(Qt.ItemDataRole.UserRole)
        if not (0 <= idx < len(self._aggregates)):
            return
        menu = QMenu(self)
        a_save = QAction("Save to NeXus…", self)
        a_save.triggered.connect(self._on_save)
        menu.addAction(a_save)
        a_mc = QAction("Compute Monte Carlo I(Q)", self)
        a_mc.triggered.connect(self._on_compute_mc)
        menu.addAction(a_mc)
        a_remove = QAction("Remove from list", self)
        a_remove.triggered.connect(lambda: self._remove_aggregate(idx))
        menu.addAction(a_remove)
        menu.exec(self.agg_list.mapToGlobal(pos))

    def _remove_aggregate(self, idx: int):
        if not (0 <= idx < len(self._aggregates)):
            return
        a = self._aggregates.pop(idx)
        if self._active is a:
            self._active = None
        self._refresh_aggregates_list()
        if self._active is None:
            self.iq_plot.clear()
            self.slice_viewer.clear()
            self.voxel3d_viewer.clear()
            self._update_params_display()
