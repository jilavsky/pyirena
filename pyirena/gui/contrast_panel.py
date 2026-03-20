"""
contrast_panel.py — Scattering Contrast Calculator GUI.

Computes X-ray and neutron scattering length densities (SLDs) and contrast for
two user-defined compounds using:
  * periodictable  — neutron b_c, molecular weights, element data
  * xraydb         — Chantler anomalous scattering factors (f1) and absorption

Entry points
------------
* ``ContrastPanel`` — QWidget launched from Data Selector (Support Tools group).
* ``main()``        — ``pyirena-contrast`` CLI command.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

try:
    from PySide6.QtWidgets import (
        QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
        QPushButton, QLabel, QLineEdit, QComboBox, QDoubleSpinBox, QSpinBox,
        QTableWidget, QTableWidgetItem, QHeaderView, QGroupBox,
        QScrollArea, QMessageBox, QFileDialog, QSplitter,
        QAbstractItemView, QCheckBox, QInputDialog, QMenu,
    )
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QColor, QFont, QAction
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
            QPushButton, QLabel, QLineEdit, QComboBox, QDoubleSpinBox, QSpinBox,
            QTableWidget, QTableWidgetItem, QHeaderView, QGroupBox,
            QScrollArea, QMessageBox, QFileDialog, QSplitter,
            QAbstractItemView, QCheckBox, QInputDialog, QMenu,
        )
        from PyQt6.QtCore import Qt
        from PyQt6.QtGui import QColor, QFont, QAction
    except ImportError:
        raise ImportError(
            "Neither PySide6 nor PyQt6 found.  Install with: pip install PySide6"
        )

import pyqtgraph as pg
from pyqtgraph.exporters import ImageExporter

from pyirena.core.scattering_contrast import (
    VACUUM, CompoundProperties, ContrastResult,
    compute_compound, compute_contrast, compute_contrast_anomalous,
    compute_anomalous_scan, get_isotopes_for_element, parse_formula,
)
from pyirena.io.contrast_io import (
    DEFAULT_LIBRARY_PATH,
    list_compounds_in_library,
    save_compound_to_library,
    load_compound_from_library,
    delete_compound_from_library,
    export_results_csv,
    export_scan_csv,
    save_scan_to_hdf5,
)

# ─── Button styles ─────────────────────────────────────────────────────────────
_S_BLUE = (
    "QPushButton{background:#2980b9;color:white;font-weight:bold;"
    "border-radius:4px;padding:5px 10px;}"
    "QPushButton:hover{background:#2471a3;}"
    "QPushButton:disabled{background:#95a5a6;}"
)
_S_PURPLE = (
    "QPushButton{background:#8e44ad;color:white;font-weight:bold;"
    "border-radius:4px;padding:5px 10px;}"
    "QPushButton:hover{background:#7d3c98;}"
    "QPushButton:disabled{background:#95a5a6;}"
)
_S_TEAL = (
    "QPushButton{background:#16a085;color:white;font-weight:bold;"
    "border-radius:4px;padding:5px 8px;}"
    "QPushButton:hover{background:#138d75;}"
    "QPushButton:disabled{background:#95a5a6;}"
)
_S_SMALL = (
    "QPushButton{background:#7f8c8d;color:white;font-size:10px;"
    "border-radius:3px;padding:2px 6px;}"
    "QPushButton:hover{background:#95a5a6;}"
)
_S_DANGER = (
    "QPushButton{background:#c0392b;color:white;font-size:10px;"
    "border-radius:3px;padding:2px 6px;}"
    "QPushButton:hover{background:#e74c3c;}"
)
_S_GREEN_SM = (
    "QPushButton{background:#27ae60;color:white;font-size:10px;"
    "border-radius:3px;padding:2px 6px;}"
    "QPushButton:hover{background:#2ecc71;}"
)
_S_FILE_SM = (
    "QPushButton{background:#2471a3;color:white;font-size:10px;"
    "border-radius:3px;padding:2px 6px;}"
    "QPushButton:hover{background:#1a5276;}"
)

# Composition modes: (display label, internal key)
_MODES: List[tuple] = [
    ("Atomic formula  (e.g. H2O, Fe3O4)", "atomic_ratio"),
    ("Wt-fractions of elements  (e.g. Au0.35Ag0.65)", "weight_fraction_elements"),
    ("Wt-fractions of compounds  (e.g. Y2O3:0.10 ZrO2:0.90)", "weight_fraction_compounds"),
]
_MODE_KEYS: List[str] = [m[1] for m in _MODES]

# GroupBox title colours
_C_COMP1  = "#1a5276"   # dark blue
_C_COMP2  = "#7b241c"   # dark red
_C_PARAMS = "#1e8449"   # dark green
_C_SCAN   = "#6c3483"   # dark purple

# Results table row definitions  (row_type, label, units, comp1_key, comp2_key)
_T_ROWS: List[tuple] = [
    ("hdr", "── Molecular Properties ──", "", None, None),
    ("cpd", "Molecular weight", "g/mol", "mol_weight", "mol_weight"),
    ("cpd", "Weight / formula unit", "g", "weight_1mol", "weight_1mol"),
    ("cpd", "Formula units / cm³", "cm⁻³", "n_mol_per_cm3", "n_mol_per_cm3"),
    ("cpd", "Electrons / formula unit", "", "n_electrons_per_mol", "n_electrons_per_mol"),
    ("cpd", "Electrons / cm³", "cm⁻³", "n_electrons_per_cm3", "n_electrons_per_cm3"),
    ("cpd", "Volume / formula unit", "cm³", "volume_1mol", "volume_1mol"),
    ("hdr", "── X-ray (free electron) ──", "", None, None),
    ("cpd", "X-ray SLD", "10¹⁰ cm⁻²", "xray_sld", "xray_sld"),
    ("cpd", "X-ray SLD / gram", "10¹⁰ cm/g", "xray_sld_per_gram", "xray_sld_per_gram"),
    ("hdr", "── Neutron ──", "", None, None),
    ("cpd", "Total neutron b", "cm", "neutron_total_b", "neutron_total_b"),
    ("cpd", "Neutron SLD", "10¹⁰ cm⁻²", "neutron_sld", "neutron_sld"),
    ("cpd", "Neutron SLD / gram", "10¹⁰ cm/g", "neutron_sld_per_gram", "neutron_sld_per_gram"),
    ("hdr", "── Contrast ──", "", None, None),
    ("ctr", "X-ray contrast (Δρ)²", "10²⁰ cm⁻⁴", "xray_contrast", None),
    ("ctr", "Neutron contrast (Δρ)²", "10²⁰ cm⁻⁴", "neutron_contrast", None),
    ("ctr", "X-ray / Neutron contrast ratio", "", "ratio_xn", None),
    ("hdr", "── Anomalous X-ray (Chantler tables) ──", "", None, None),
    ("ano", "X-ray SLD (anomalous)", "10¹⁰ cm⁻²", "xray_sld_anom_1", "xray_sld_anom_2"),
    ("ano", "Linear absorption μ", "cm⁻¹", "mu_1", "mu_2"),
    ("ano", "Transmission (compound)", "", "transmission_1", "transmission_2"),
    ("ano", "Sample transmission", "", "transmission_sample", None),
    ("ano", "X-ray contrast (anomalous) (Δρ)²", "10²⁰ cm⁻⁴", "xray_contrast_anom", None),
]


# ─── Module helpers ────────────────────────────────────────────────────────────

def _fmt(val: Any, sig: int = 4) -> str:
    if val is None:
        return "—"
    try:
        v = float(val)
        if not np.isfinite(v):
            return str(v)
        if v == 0:
            return "0"
        mag = abs(v)
        if 0.001 <= mag < 1e5:
            return f"{v:.{sig}g}"
        return f"{v:.{sig - 1}e}"
    except Exception:
        return str(val)


def _grp_style(color: str) -> str:
    """Return a QGroupBox stylesheet that colours the title and border."""
    return (
        f"QGroupBox{{font-weight:bold;border:2px solid {color};"
        f"border-radius:5px;margin-top:10px;padding-top:2px;}}"
        f"QGroupBox::title{{subcontrol-origin:margin;"
        f"subcontrol-position:top left;padding:0 6px;color:{color};font-weight:bold;}}"
    )


def _add_jpeg_action(plot_item, parent: QWidget, default_name: str = "graph") -> None:
    action = QAction("Save graph as JPEG…", parent)

    def _go() -> None:
        path, _ = QFileDialog.getSaveFileName(
            parent,
            "Save graph as JPEG",
            str(Path.home() / f"{default_name}.jpg"),
            "JPEG Images (*.jpg *.jpeg)",
        )
        if not path:
            return
        try:
            exp = ImageExporter(plot_item)
            exp.parameters()["width"] = 1600
            exp.export(path)
        except Exception as exc:
            QMessageBox.warning(parent, "Export failed", str(exc))

    action.triggered.connect(_go)
    plot_item.getViewBox().menu.addAction(action)


def _ro_item(text: str) -> QTableWidgetItem:
    item = QTableWidgetItem(text)
    item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
    return item


def _hdr_item(text: str) -> QTableWidgetItem:
    item = QTableWidgetItem(text)
    item.setFlags(Qt.ItemFlag.ItemIsEnabled)
    f = item.font()
    f.setBold(True)
    item.setFont(f)
    item.setBackground(QColor("#dfe6e9"))
    return item


def _configure_white_plot(ax: pg.PlotItem) -> None:
    """Configure axis pens and labels for a white-background pyqtgraph plot."""
    dark = pg.mkPen("#2c3e50", width=1)
    dark_text = pg.mkPen("#2c3e50")
    for side in ("left", "bottom", "right", "top"):
        axis = ax.getAxis(side)
        axis.setPen(dark)
        axis.setTextPen(dark_text)


# ─── Energy Scan Graph Window ─────────────────────────────────────────────────

class ContrastGraphWindow(QWidget):
    """Floating window with three energy-scan plots and a crosshair cursor."""

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle("Scattering Contrast — Energy Scan")
        self.resize(820, 720)
        self._proxies: list = []          # keep SignalProxy refs alive
        self._vlines: list = []
        self._hlines: list = []
        self._scan_data: Optional[Dict] = None
        self._build_ui()

    def _build_ui(self) -> None:
        lay = QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)

        # Toolbar row
        tb = QHBoxLayout()
        self._log_y_cb = QCheckBox("Log Y axis (contrast plot)")
        self._log_y_cb.setChecked(False)
        self._log_y_cb.toggled.connect(self._toggle_log_y)
        tb.addWidget(self._log_y_cb)
        tb.addStretch()
        lay.addLayout(tb)

        self.gl = pg.GraphicsLayoutWidget()
        self.gl.setBackground("w")        # white background
        lay.addWidget(self.gl)

        # ── Plot 1: X-ray contrast (Δρ)² ──
        self._ax_ctr = self.gl.addPlot(row=0, col=0)
        self._ax_ctr.setLabel("left", "(Δρ)² X-ray", units="10²⁰ cm⁻⁴")
        self._ax_ctr.setLabel("bottom", "X-ray energy", units="keV")
        self._ax_ctr.showGrid(x=True, y=True, alpha=0.25)
        _configure_white_plot(self._ax_ctr)
        _add_jpeg_action(self._ax_ctr, self, "xray_contrast_vs_energy")

        # ── Plot 2: Linear absorption ──
        self._ax_mu = self.gl.addPlot(row=1, col=0)
        self._ax_mu.setLabel("left", "Linear absorption μ", units="cm⁻¹")
        self._ax_mu.setLabel("bottom", "X-ray energy", units="keV")
        self._ax_mu.showGrid(x=True, y=True, alpha=0.25)
        _configure_white_plot(self._ax_mu)
        _add_jpeg_action(self._ax_mu, self, "absorption_vs_energy")

        # ── Plot 3: Transmission ──
        self._ax_tr = self.gl.addPlot(row=2, col=0)
        self._ax_tr.setLabel("left", "Transmission")
        self._ax_tr.setLabel("bottom", "X-ray energy", units="keV")
        self._ax_tr.showGrid(x=True, y=True, alpha=0.25)
        self._ax_tr.setYRange(0, 1.05, padding=0)
        _configure_white_plot(self._ax_tr)
        _add_jpeg_action(self._ax_tr, self, "transmission_vs_energy")

        # Link x-axes
        self._ax_mu.setXLink(self._ax_ctr)
        self._ax_tr.setXLink(self._ax_ctr)

        # Crosshair cursor lines (one pair per plot)
        dash = pg.mkPen("#888888", width=1, style=Qt.PenStyle.DashLine)
        self._axes_list = [self._ax_ctr, self._ax_mu, self._ax_tr]
        for ax in self._axes_list:
            vl = pg.InfiniteLine(angle=90, movable=False, pen=dash)
            hl = pg.InfiniteLine(angle=0, movable=False, pen=dash)
            hl.setVisible(False)
            ax.addItem(vl, ignoreBounds=True)
            ax.addItem(hl, ignoreBounds=True)
            self._vlines.append(vl)
            self._hlines.append(hl)

        # Set up mouse-move proxies for each plot
        for i, ax in enumerate(self._axes_list):
            proxy = pg.SignalProxy(
                ax.scene().sigMouseMoved,
                rateLimit=60,
                slot=lambda ev, _ax=ax, _i=i: self._mouse_moved(ev, _ax, _i),
            )
            self._proxies.append(proxy)

        # Status / cursor readout row
        bot = QHBoxLayout()
        self._scan_lbl = QLabel("No scan data.")
        self._scan_lbl.setStyleSheet("color:#7f8c8d; font-size:11px;")
        bot.addWidget(self._scan_lbl, stretch=1)
        self._cursor_lbl = QLabel("")
        self._cursor_lbl.setStyleSheet(
            "color:#2c3e50; font-size:11px; font-weight:bold; padding:0 4px;"
        )
        bot.addWidget(self._cursor_lbl)
        lay.addLayout(bot)

    # ------------------------------------------------------------------
    def _mouse_moved(self, event: tuple, ax: pg.PlotItem, plot_idx: int) -> None:
        pos = event[0]
        if not ax.sceneBoundingRect().contains(pos):
            return
        mp = ax.getViewBox().mapSceneToView(pos)
        x, y = mp.x(), mp.y()

        # Move all vertical lines to the same energy
        for vl in self._vlines:
            vl.setPos(x)

        # Show horizontal line only in the active plot
        for i, hl in enumerate(self._hlines):
            hl.setVisible(i == plot_idx)
            if i == plot_idx:
                hl.setPos(y)

        # Update readout label
        unit_map = {0: "10²⁰ cm⁻⁴", 1: "cm⁻¹", 2: ""}
        lbl_map  = {0: "(Δρ)²", 1: "μ", 2: "T"}
        units = unit_map.get(plot_idx, "")
        lbl   = lbl_map.get(plot_idx, "Y")
        u_str = f" {units}" if units else ""
        self._cursor_lbl.setText(f"E: {x:.4f} keV  │  {lbl}: {y:.4g}{u_str}")

    # ------------------------------------------------------------------
    def _toggle_log_y(self, checked: bool) -> None:
        self._ax_ctr.setLogMode(y=checked)
        self._ax_ctr.autoRange()

    # ------------------------------------------------------------------
    def update_scan(
        self,
        scan_data: Dict[str, np.ndarray],
        comp1_name: str,
        comp2_name: str,
    ) -> None:
        self._scan_data = scan_data
        E = scan_data.get("energy")
        if E is None or len(E) == 0:
            return

        self._ax_ctr.clear()
        self._ax_mu.clear()
        self._ax_tr.clear()
        # Re-add crosshair items (clear() removes them)
        dash = pg.mkPen("#888888", width=1, style=Qt.PenStyle.DashLine)
        for i, ax in enumerate(self._axes_list):
            ax.addItem(self._vlines[i], ignoreBounds=True)
            ax.addItem(self._hlines[i], ignoreBounds=True)
        for hl in self._hlines:
            hl.setVisible(False)

        self._ax_ctr.setLogMode(y=self._log_y_cb.isChecked())

        c1 = comp1_name or "Compound 1"
        c2 = comp2_name or "Compound 2"

        # Contrast
        contrast = scan_data.get("xray_contrast_anom", np.zeros(len(E)))
        self._ax_ctr.plot(E, contrast, pen=pg.mkPen("#2980b9", width=2))
        self._ax_ctr.autoRange()

        # Absorption
        for key, label, color in [
            ("mu_1", c1, "#e74c3c"),
            ("mu_2", c2, "#2980b9"),
        ]:
            arr = scan_data.get(key)
            if arr is not None:
                self._ax_mu.plot(E, arr, pen=pg.mkPen(color, width=1.5), name=label)
        self._ax_mu.addLegend(offset=(5, 5))
        self._ax_mu.autoRange()

        # Transmission
        for key, label, color in [
            ("transmission_1", c1, "#e74c3c"),
            ("transmission_2", c2, "#2980b9"),
            ("transmission_sample", "Sample", "#27ae60"),
        ]:
            arr = scan_data.get(key)
            if arr is not None:
                self._ax_tr.plot(E, arr, pen=pg.mkPen(color, width=1.5), name=label)
        self._ax_tr.addLegend(offset=(5, 5))
        self._ax_tr.setYRange(0, 1.05, padding=0)

        self._scan_lbl.setText(
            f"Energy: {E[0]:.2f}–{E[-1]:.2f} keV, {len(E)} points | {c1} vs {c2}"
            "   (move mouse over plot for cursor readout)"
        )


# ─── Main contrast panel ───────────────────────────────────────────────────────

class ContrastPanel(QWidget):
    """
    Scattering Contrast Calculator.

    Left panel: two compound editors + calculation parameters.
    Right panel: results table (with right-click Copy).
    Separate ContrastGraphWindow: energy-scan plots with crosshair.
    """

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        state_manager=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Scattering Contrast Calculator")
        self.resize(1120, 740)

        self._state_mgr = state_manager
        self._scan_data: Optional[Dict[str, np.ndarray]] = None
        self._comp1: Optional[CompoundProperties] = None
        self._comp2: Optional[CompoundProperties] = None
        self._contrast: Optional[ContrastResult] = None
        self._graph_win: Optional[ContrastGraphWindow] = None
        self._c1_refs: Dict[str, Any] = {}
        self._c2_refs: Dict[str, Any] = {}

        self._build_ui()
        self._refresh_library(1)
        self._refresh_library(2)
        self._restore_state()

    # ══════════════════════════════════════════════════════════════════
    #  UI construction
    # ══════════════════════════════════════════════════════════════════

    def _build_ui(self) -> None:
        root = QHBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(6)

        # ── Left scroll panel ──────────────────────────────────────────
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMinimumWidth(390)
        scroll.setMaximumWidth(450)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        left_w = QWidget()
        ll = QVBoxLayout(left_w)
        ll.setContentsMargins(4, 4, 4, 4)
        ll.setSpacing(6)

        self._build_compound_group(1, "Compound 1 (Phase)", "H2O", 1.0, ll)
        self._build_compound_group(2, "Compound 2 (Matrix / Solvent)", "", 0.0, ll)
        ll.addWidget(self._make_params_group())
        ll.addWidget(self._make_scan_range_group())
        ll.addStretch()
        ll.addWidget(self._make_action_buttons())

        scroll.setWidget(left_w)
        root.addWidget(scroll)

        # ── Right panel ────────────────────────────────────────────────
        right_w = QWidget()
        rl = QVBoxLayout(right_w)
        rl.setContentsMargins(0, 0, 0, 0)
        rl.setSpacing(4)

        self._status_lbl = QLabel("Enter compound definitions and click Calculate.")
        self._status_lbl.setStyleSheet("color:#7f8c8d; padding:3px; font-style:italic;")
        rl.addWidget(self._status_lbl)

        self._tbl = self._make_results_table()
        rl.addWidget(self._tbl, stretch=1)
        rl.addWidget(self._make_export_buttons())

        root.addWidget(right_w, stretch=1)

    # ------------------------------------------------------------------
    def _build_compound_group(
        self, n: int, title: str, formula: str, density: float,
        parent_lay: QVBoxLayout,
    ) -> None:
        refs = self._c1_refs if n == 1 else self._c2_refs
        color = _C_COMP1 if n == 1 else _C_COMP2

        grp = QGroupBox(title)
        grp.setStyleSheet(_grp_style(color))
        lay = QVBoxLayout(grp)
        lay.setSpacing(4)

        # Name
        row = QHBoxLayout()
        row.addWidget(QLabel("Name:"))
        name_le = QLineEdit()
        name_le.setPlaceholderText("e.g. Water")
        refs["name"] = name_le
        row.addWidget(name_le)
        lay.addLayout(row)

        # Formula + Parse
        row = QHBoxLayout()
        row.addWidget(QLabel("Formula:"))
        formula_le = QLineEdit(formula)
        formula_le.setPlaceholderText("e.g. H2O")
        refs["formula"] = formula_le
        row.addWidget(formula_le)
        pb = QPushButton("Parse ▸")
        pb.setStyleSheet(_S_SMALL)
        pb.setFixedWidth(58)
        pb.setToolTip("Parse formula and update the isotope selection table below.")
        # NOTE: use *_ to absorb the 'checked' bool that clicked() emits
        pb.clicked.connect(lambda *_, _n=n: self._update_isotope_table(_n))
        row.addWidget(pb)
        lay.addLayout(row)

        # Composition mode
        row = QHBoxLayout()
        row.addWidget(QLabel("Mode:"))
        mode_cb = QComboBox()
        for lbl, _ in _MODES:
            mode_cb.addItem(lbl)
        refs["mode"] = mode_cb
        row.addWidget(mode_cb)
        lay.addLayout(row)

        # Density
        row = QHBoxLayout()
        row.addWidget(QLabel("Density:"))
        dspin = QDoubleSpinBox()
        dspin.setRange(0.0, 99.0)
        dspin.setDecimals(4)
        dspin.setSuffix("  g/cm³")
        dspin.setValue(density)
        refs["density"] = dspin
        row.addWidget(dspin)
        vac = QPushButton("Set Vacuum")
        vac.setStyleSheet(_S_SMALL)
        vac.setToolTip("Set this compound to vacuum (empty formula, density = 0).")
        vac.clicked.connect(lambda *_, _n=n: self._set_vacuum(_n))
        row.addWidget(vac)
        lay.addLayout(row)

        # Isotope table
        iso_grp = QGroupBox("Neutron Isotope Selection")
        iso_grp.setStyleSheet("QGroupBox{font-size:10px;}")
        iso_lay = QVBoxLayout(iso_grp)
        iso_lay.setContentsMargins(4, 8, 4, 4)

        hint = QLabel(
            "Click Parse ▸ or Calculate to populate.\n"
            "'natural (default)' uses natural-abundance b_c."
        )
        hint.setStyleSheet("font-size:10px; color:#7f8c8d;")
        hint.setWordWrap(True)
        iso_lay.addWidget(hint)

        iso_tbl = QTableWidget(0, 2)
        iso_tbl.setHorizontalHeaderLabels(["Element", "Isotope"])
        iso_tbl.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.ResizeToContents)
        iso_tbl.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.ResizeMode.Stretch)
        iso_tbl.verticalHeader().setVisible(False)
        iso_tbl.setMaximumHeight(110)
        iso_tbl.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        iso_tbl.setStyleSheet("font-size:11px;")
        refs["iso_table"] = iso_tbl
        iso_lay.addWidget(iso_tbl)
        lay.addWidget(iso_grp)

        # Library controls
        lib_grp = QGroupBox("Saved Compounds Library")
        lib_grp.setStyleSheet("QGroupBox{font-size:10px;}")
        lib_lay = QVBoxLayout(lib_grp)
        lib_lay.setContentsMargins(4, 8, 4, 4)
        lib_lay.setSpacing(3)

        lc = QComboBox()
        lc.setPlaceholderText("─ select compound ─")
        lc.setToolTip(f"Library: {DEFAULT_LIBRARY_PATH}")
        refs["lib_combo"] = lc

        # Row 1: combo + action buttons
        row1 = QHBoxLayout()
        row1.setSpacing(3)
        row1.addWidget(lc, stretch=1)
        for sym, tip, meth, style, w in [
            ("↻",    "Refresh list",        "_refresh_library",       _S_SMALL,    28),
            ("Load", "Load into this slot", "_load_from_library",     _S_SMALL,    40),
            ("Save", "Save to library",     "_save_to_library",       _S_GREEN_SM, 40),
            ("Del",  "Delete from library", "_delete_from_library",   _S_DANGER,   38),
        ]:
            b = QPushButton(sym)
            b.setStyleSheet(style)
            b.setFixedWidth(w)
            b.setToolTip(tip)
            # CRITICAL fix: use *_ to absorb the bool 'checked' emitted by clicked()
            _meth = meth
            b.clicked.connect(
                lambda *_, _n=n, _m=_meth: getattr(self, _m)(_n)
            )
            row1.addWidget(b)
        lib_lay.addLayout(row1)

        # Row 2: file-based export / import
        row2 = QHBoxLayout()
        row2.setSpacing(3)
        exp_f = QPushButton("Export to File…")
        exp_f.setStyleSheet(_S_FILE_SM)
        exp_f.setToolTip(
            "Save this compound definition to a portable HDF5 file\n"
            "that can be shared with other users / computers."
        )
        exp_f.clicked.connect(lambda *_, _n=n: self._export_compound_to_file(_n))
        row2.addWidget(exp_f, stretch=1)

        imp_f = QPushButton("Import from File…")
        imp_f.setStyleSheet(_S_FILE_SM)
        imp_f.setToolTip(
            "Load a compound from a shared HDF5 file into this slot.\n"
            "Optionally also adds it to the local library."
        )
        imp_f.clicked.connect(lambda *_, _n=n: self._import_compound_from_file(_n))
        row2.addWidget(imp_f, stretch=1)
        lib_lay.addLayout(row2)

        lay.addWidget(lib_grp)
        parent_lay.addWidget(grp)

    # ------------------------------------------------------------------
    def _make_params_group(self) -> QGroupBox:
        grp = QGroupBox("Calculation Parameters")
        grp.setStyleSheet(_grp_style(_C_PARAMS))
        g = QGridLayout(grp)
        g.setColumnStretch(1, 1)

        g.addWidget(QLabel("X-ray energy:"), 0, 0)
        self._energy_spin = QDoubleSpinBox()
        self._energy_spin.setRange(0.5, 200.0)
        self._energy_spin.setDecimals(3)
        self._energy_spin.setValue(12.0)
        self._energy_spin.setSuffix("  keV")
        g.addWidget(self._energy_spin, 0, 1)

        g.addWidget(QLabel("Sample thickness:"), 1, 0)
        self._thick_spin = QDoubleSpinBox()
        self._thick_spin.setRange(0.0, 10000.0)
        self._thick_spin.setDecimals(3)
        self._thick_spin.setValue(1.0)
        self._thick_spin.setSuffix("  mm")
        g.addWidget(self._thick_spin, 1, 1)

        g.addWidget(QLabel("Vol. frac. compound 1:"), 2, 0)
        self._vf_spin = QDoubleSpinBox()
        self._vf_spin.setRange(0.0, 1.0)
        self._vf_spin.setDecimals(4)
        self._vf_spin.setValue(0.01)
        self._vf_spin.setToolTip(
            "Volume fraction of Compound 1 in the mixed sample.\n"
            "Used to calculate the combined sample transmission."
        )
        g.addWidget(self._vf_spin, 2, 1)

        return grp

    # ------------------------------------------------------------------
    def _make_scan_range_group(self) -> QGroupBox:
        grp = QGroupBox("Energy Scan Range")
        grp.setStyleSheet(_grp_style(_C_SCAN))
        g = QGridLayout(grp)
        g.setColumnStretch(1, 1)

        g.addWidget(QLabel("Start energy:"), 0, 0)
        self._e_start = QDoubleSpinBox()
        self._e_start.setRange(0.5, 200.0)
        self._e_start.setDecimals(2)
        self._e_start.setValue(5.0)
        self._e_start.setSuffix("  keV")
        g.addWidget(self._e_start, 0, 1)

        g.addWidget(QLabel("End energy:"), 1, 0)
        self._e_end = QDoubleSpinBox()
        self._e_end.setRange(0.5, 200.0)
        self._e_end.setDecimals(2)
        self._e_end.setValue(25.0)
        self._e_end.setSuffix("  keV")
        g.addWidget(self._e_end, 1, 1)

        g.addWidget(QLabel("Number of points:"), 2, 0)
        self._n_pts = QSpinBox()
        self._n_pts.setRange(10, 5000)
        self._n_pts.setValue(500)
        g.addWidget(self._n_pts, 2, 1)

        return grp

    # ------------------------------------------------------------------
    def _make_action_buttons(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(5)

        calc_btn = QPushButton("Calculate")
        calc_btn.setStyleSheet(_S_BLUE)
        calc_btn.setMinimumHeight(36)
        calc_btn.setToolTip(
            "Compute SLDs, free-electron and anomalous contrast at the\n"
            "specified X-ray energy and fill the results table."
        )
        calc_btn.clicked.connect(self._on_calculate)
        lay.addWidget(calc_btn)

        scan_btn = QPushButton("Energy Scan")
        scan_btn.setStyleSheet(_S_PURPLE)
        scan_btn.setMinimumHeight(34)
        scan_btn.setToolTip(
            "Compute anomalous X-ray contrast and transmission over the\n"
            "energy scan range and open the graph window."
        )
        scan_btn.clicked.connect(self._on_energy_scan)
        lay.addWidget(scan_btn)

        return w

    # ------------------------------------------------------------------
    def _make_results_table(self) -> QTableWidget:
        tbl = QTableWidget(len(_T_ROWS), 4)
        tbl.setHorizontalHeaderLabels(
            ["Property", "Units", "Compound 1", "Compound 2"]
        )
        tbl.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        tbl.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.ResizeMode.ResizeToContents)
        tbl.horizontalHeader().setSectionResizeMode(
            2, QHeaderView.ResizeMode.ResizeToContents)
        tbl.horizontalHeader().setSectionResizeMode(
            3, QHeaderView.ResizeMode.ResizeToContents)
        tbl.verticalHeader().setVisible(False)
        tbl.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        tbl.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        tbl.setAlternatingRowColors(True)

        # Right-click Copy context menu
        tbl.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        tbl.customContextMenuRequested.connect(self._table_context_menu)

        for row, (rtype, label, units, _k1, _k2) in enumerate(_T_ROWS):
            if rtype == "hdr":
                tbl.setItem(row, 0, _hdr_item(label))
                for col in range(1, 4):
                    it = QTableWidgetItem("")
                    it.setFlags(Qt.ItemFlag.ItemIsEnabled)
                    it.setBackground(QColor("#dfe6e9"))
                    tbl.setItem(row, col, it)
                tbl.setRowHeight(row, 22)
            else:
                tbl.setItem(row, 0, _ro_item(label))
                tbl.setItem(row, 1, _ro_item(units))
                tbl.setItem(row, 2, _ro_item(""))
                tbl.setItem(row, 3, _ro_item(""))

        return tbl

    # ------------------------------------------------------------------
    def _make_export_buttons(self) -> QWidget:
        w = QWidget()
        lay = QHBoxLayout(w)
        lay.setContentsMargins(0, 2, 0, 2)
        lay.setSpacing(5)

        for lbl, tip, slot in [
            ("Export Results CSV",
             "Save the results table to CSV.  Ctrl+C also copies selected cells.",
             self._export_csv),
            ("Export Scan CSV",
             "Save the energy-scan arrays to CSV.",
             self._export_scan_csv),
            ("Save Scan HDF5",
             "Save the energy-scan arrays plus compound definitions to HDF5.",
             self._save_scan_hdf5),
            ("Show Graphs",
             "Open the energy-scan graph window.",
             self._show_graph_window),
        ]:
            btn = QPushButton(lbl)
            btn.setStyleSheet(_S_TEAL if lbl != "Show Graphs" else _S_PURPLE)
            btn.setToolTip(tip)
            btn.clicked.connect(slot)
            lay.addWidget(btn)

        return w

    # ══════════════════════════════════════════════════════════════════
    #  Table context menu (right-click Copy)
    # ══════════════════════════════════════════════════════════════════

    def _table_context_menu(self, pos) -> None:
        menu = QMenu(self._tbl)
        copy_act = QAction("Copy", self._tbl)
        copy_act.setShortcut("Ctrl+C")
        copy_act.triggered.connect(self._copy_table_selection)
        menu.addAction(copy_act)
        menu.exec(self._tbl.viewport().mapToGlobal(pos))

    def _copy_table_selection(self) -> None:
        items = self._tbl.selectedItems()
        if not items:
            return
        # Collect by row, then column, then join
        rows: Dict[int, Dict[int, str]] = {}
        for it in items:
            rows.setdefault(it.row(), {})[it.column()] = it.text()
        text = "\n".join(
            "\t".join(rows[r].get(c, "") for c in sorted(rows[r]))
            for r in sorted(rows)
        )
        QApplication.clipboard().setText(text)

    # ══════════════════════════════════════════════════════════════════
    #  Isotope table
    # ══════════════════════════════════════════════════════════════════

    def _refs(self, n: int) -> Dict[str, Any]:
        return self._c1_refs if n == 1 else self._c2_refs

    def _update_isotope_table(self, n: int) -> None:
        refs = self._refs(n)
        formula = refs["formula"].text().strip()
        mode = _MODE_KEYS[refs["mode"].currentIndex()]

        # Pending overrides (from state restore) take priority; otherwise read table
        pending = refs.pop("_pending_iso", None)
        saved = pending if pending is not None else self._read_isotope_overrides(n)

        iso_tbl: QTableWidget = refs["iso_table"]
        iso_tbl.setRowCount(0)

        if not formula or formula.lower() == "vacuum":
            return

        try:
            element_counts = parse_formula(formula, mode)
            elements = sorted(element_counts.keys())
        except Exception as exc:
            self._set_status(f"Formula parse error: {exc}", error=True)
            return

        iso_tbl.setRowCount(len(elements))
        for row, sym in enumerate(elements):
            el_item = QTableWidgetItem(sym)
            el_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
            iso_tbl.setItem(row, 0, el_item)

            combo = QComboBox()
            combo.addItem("natural (default)")
            for label, bc in get_isotopes_for_element(sym):
                if label != "natural":
                    combo.addItem(f"{label}  [b_c = {bc:.3f} fm]")

            if sym in saved:
                target = saved[sym]
                for i in range(combo.count()):
                    if combo.itemText(i).split()[0] == target:
                        combo.setCurrentIndex(i)
                        break

            iso_tbl.setCellWidget(row, 1, combo)

        iso_tbl.resizeRowsToContents()

    def _read_isotope_overrides(self, n: int) -> Dict[str, str]:
        tbl: QTableWidget = self._refs(n)["iso_table"]
        overrides: Dict[str, str] = {}
        for row in range(tbl.rowCount()):
            sym_item = tbl.item(row, 0)
            if not sym_item:
                continue
            sym = sym_item.text()
            combo: Optional[QComboBox] = tbl.cellWidget(row, 1)
            if combo is None:
                continue
            text = combo.currentText().strip()
            if text and not text.startswith("natural"):
                iso_num = text.split()[0]
                if iso_num.isdigit():
                    overrides[sym] = iso_num
        return overrides

    # ══════════════════════════════════════════════════════════════════
    #  Compound input accessors
    # ══════════════════════════════════════════════════════════════════

    def _get_compound_spec(self, n: int) -> Dict[str, Any]:
        refs = self._refs(n)
        return {
            "name": refs["name"].text().strip() or f"Compound {n}",
            "formula_str": refs["formula"].text().strip(),
            "composition_mode": _MODE_KEYS[refs["mode"].currentIndex()],
            "density": refs["density"].value(),
            "isotope_overrides": self._read_isotope_overrides(n),
        }

    def _set_vacuum(self, n: int) -> None:
        refs = self._refs(n)
        refs["formula"].setText("")
        refs["density"].setValue(0.0)
        refs["iso_table"].setRowCount(0)
        refs["name"].setText("vacuum")

    def _compute_one(self, n: int) -> Optional[CompoundProperties]:
        spec = self._get_compound_spec(n)
        if not spec["formula_str"] or spec["density"] == 0.0:
            return VACUUM
        try:
            return compute_compound(
                formula_str=spec["formula_str"],
                density=spec["density"],
                mode=spec["composition_mode"],
                isotope_overrides=spec["isotope_overrides"],
                name=spec["name"],
            )
        except Exception as exc:
            self._set_status(f"Compound {n} error: {exc}", error=True)
            return None

    # ══════════════════════════════════════════════════════════════════
    #  Calculation
    # ══════════════════════════════════════════════════════════════════

    def _on_calculate(self) -> None:
        for n in (1, 2):
            if self._refs(n)["iso_table"].rowCount() == 0:
                self._update_isotope_table(n)

        comp1 = self._compute_one(1)
        comp2 = self._compute_one(2)
        if comp1 is None or comp2 is None:
            return

        energy_keV = self._energy_spin.value()
        thickness_mm = self._thick_spin.value()
        vol_frac = self._vf_spin.value()

        try:
            contrast = compute_contrast_anomalous(
                comp1, comp2,
                energy_keV=energy_keV,
                thickness_mm=thickness_mm,
                vol_frac_comp1=vol_frac,
            )
        except Exception as exc:
            self._set_status(f"Calculation error: {exc}", error=True)
            return

        self._comp1 = comp1
        self._comp2 = comp2
        self._contrast = contrast
        self._fill_results_table(comp1, comp2, contrast)

        self._set_status(
            f"Calculated at {energy_keV:.3f} keV  |  "
            f"X-ray contrast (anom): {_fmt(contrast.xray_contrast_anom)} × 10²⁰ cm⁻⁴  |  "
            f"Neutron contrast: {_fmt(contrast.neutron_contrast)} × 10²⁰ cm⁻⁴"
        )
        self._save_state()

    def _on_energy_scan(self) -> None:
        for n in (1, 2):
            if self._refs(n)["iso_table"].rowCount() == 0:
                self._update_isotope_table(n)

        comp1 = self._compute_one(1)
        comp2 = self._compute_one(2)
        if comp1 is None or comp2 is None:
            return

        e_start = self._e_start.value()
        e_end = self._e_end.value()
        if e_end <= e_start:
            QMessageBox.warning(
                self, "Invalid range", "End energy must be > start energy."
            )
            return

        thickness_mm = self._thick_spin.value()
        vol_frac = self._vf_spin.value()
        n_pts = self._n_pts.value()

        self._set_status(
            f"Computing energy scan {e_start:.2f}–{e_end:.2f} keV, {n_pts} points…"
        )
        QApplication.processEvents()

        try:
            scan_data = compute_anomalous_scan(
                comp1, comp2,
                e_start_keV=e_start,
                e_end_keV=e_end,
                n_points=n_pts,
                thickness_mm=thickness_mm,
                vol_frac_comp1=vol_frac,
            )
        except Exception as exc:
            self._set_status(f"Energy scan error: {exc}", error=True)
            return

        self._scan_data = scan_data
        self._comp1 = comp1
        self._comp2 = comp2

        try:
            contrast = compute_contrast_anomalous(
                comp1, comp2,
                energy_keV=self._energy_spin.value(),
                thickness_mm=thickness_mm,
                vol_frac_comp1=vol_frac,
            )
            self._contrast = contrast
            self._fill_results_table(comp1, comp2, contrast)
        except Exception:
            pass

        self._show_graph_window()
        self._graph_win.update_scan(scan_data, comp1.name, comp2.name)

        self._set_status(
            f"Energy scan complete: {e_start:.2f}–{e_end:.2f} keV, {n_pts} points."
        )
        self._save_state()

    # ══════════════════════════════════════════════════════════════════
    #  Results table
    # ══════════════════════════════════════════════════════════════════

    def _fill_results_table(
        self,
        comp1: CompoundProperties,
        comp2: CompoundProperties,
        contrast: ContrastResult,
    ) -> None:
        c1d = comp1.__dict__
        c2d = comp2.__dict__
        ctd = contrast.__dict__

        self._tbl.setHorizontalHeaderLabels([
            "Property", "Units",
            comp1.name or "Compound 1",
            comp2.name or "Compound 2",
        ])

        for row, (rtype, _lbl, _units, key1, key2) in enumerate(_T_ROWS):
            if rtype == "hdr":
                continue
            if rtype == "cpd":
                v1 = _fmt(c1d.get(key1))
                v2 = _fmt(c2d.get(key2)) if key2 else "—"
            elif rtype == "ctr":
                v1 = _fmt(ctd.get(key1))
                v2 = "—"
            elif rtype == "ano":
                v1 = _fmt(ctd.get(key1)) if key1 else "—"
                v2 = _fmt(ctd.get(key2)) if key2 else "—"
            else:
                v1 = v2 = "—"

            for col, val in [(2, v1), (3, v2)]:
                it = self._tbl.item(row, col)
                if it:
                    it.setText(val)

        self._tbl.resizeColumnsToContents()

    # ══════════════════════════════════════════════════════════════════
    #  Library operations (local library)
    # ══════════════════════════════════════════════════════════════════

    def _refresh_library(self, n: int) -> None:
        combo: QComboBox = self._refs(n)["lib_combo"]
        names = list_compounds_in_library()
        current = combo.currentText()
        combo.blockSignals(True)
        combo.clear()
        combo.addItems(names)
        idx = combo.findText(current)
        if idx >= 0:
            combo.setCurrentIndex(idx)
        combo.blockSignals(False)

    def _save_to_library(self, n: int) -> None:
        spec = self._get_compound_spec(n)
        if not spec["formula_str"]:
            QMessageBox.warning(self, "Nothing to save", "Formula is empty.")
            return
        name, ok = QInputDialog.getText(
            self, "Save Compound",
            "Compound name in library:",
            text=spec["name"],
        )
        if not ok or not name.strip():
            return
        try:
            saved = save_compound_to_library(spec, name=name.strip())
            self._refresh_library(1)
            self._refresh_library(2)
            self._set_status(f"Compound '{saved}' saved to library.")
        except Exception as exc:
            QMessageBox.critical(self, "Save failed", str(exc))

    def _load_from_library(self, n: int) -> None:
        combo: QComboBox = self._refs(n)["lib_combo"]
        name = combo.currentText()
        if not name:
            QMessageBox.information(self, "Load", "No compound selected in list.")
            return
        try:
            d = load_compound_from_library(name)
        except Exception as exc:
            QMessageBox.critical(self, "Load failed", str(exc))
            return
        self._apply_compound_dict(n, d)
        self._set_status(f"Loaded '{name}' into Compound {n}.")

    def _delete_from_library(self, n: int) -> None:
        combo: QComboBox = self._refs(n)["lib_combo"]
        name = combo.currentText()
        if not name:
            return
        ans = QMessageBox.question(
            self, "Delete compound",
            f"Delete '{name}' from the library?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if ans != QMessageBox.StandardButton.Yes:
            return
        try:
            delete_compound_from_library(name)
            self._refresh_library(1)
            self._refresh_library(2)
            self._set_status(f"Deleted '{name}' from library.")
        except Exception as exc:
            QMessageBox.critical(self, "Delete failed", str(exc))

    # ══════════════════════════════════════════════════════════════════
    #  File-based export / import (portable sharing)
    # ══════════════════════════════════════════════════════════════════

    def _export_compound_to_file(self, n: int) -> None:
        spec = self._get_compound_spec(n)
        if not spec["formula_str"]:
            QMessageBox.warning(self, "Nothing to export", "Formula is empty.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Compound to File",
            str(Path.home() / f"{spec['name']}_compounds.h5"),
            "HDF5 Files (*.h5 *.hdf5);;All Files (*)",
        )
        if not path:
            return
        try:
            saved = save_compound_to_library(
                spec, library_path=Path(path), name=spec["name"]
            )
            self._set_status(f"Exported '{saved}' to {Path(path).name}.")
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))

    def _import_compound_from_file(self, n: int) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Import Compound from File",
            str(Path.home()),
            "HDF5 Files (*.h5 *.hdf5);;All Files (*)",
        )
        if not path:
            return

        names = list_compounds_in_library(Path(path))
        if not names:
            QMessageBox.information(
                self, "No compounds",
                f"No compounds found in {Path(path).name}."
            )
            return

        name, ok = QInputDialog.getItem(
            self, "Import Compound",
            f"Select compound from {Path(path).name}:",
            names, 0, False,
        )
        if not ok:
            return

        try:
            d = load_compound_from_library(name, library_path=Path(path))
        except Exception as exc:
            QMessageBox.critical(self, "Import failed", str(exc))
            return

        self._apply_compound_dict(n, d)
        self._set_status(f"Imported '{name}' from {Path(path).name} into Compound {n}.")

        # Offer to also save to local library
        ans = QMessageBox.question(
            self, "Add to local library?",
            f"Also add '{name}' to the local library (~/.pyirena/)?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if ans == QMessageBox.StandardButton.Yes:
            try:
                save_compound_to_library(d, name=name)
                self._refresh_library(1)
                self._refresh_library(2)
            except Exception as exc:
                QMessageBox.warning(self, "Library save failed", str(exc))

    # ------------------------------------------------------------------
    def _apply_compound_dict(self, n: int, d: Dict[str, Any]) -> None:
        """Populate compound n's widgets from a compound dict (from library or file)."""
        refs = self._refs(n)
        refs["name"].setText(d.get("name", ""))
        refs["formula"].setText(d.get("formula_str", ""))
        mode_key = d.get("composition_mode", "atomic_ratio")
        refs["mode"].setCurrentIndex(
            _MODE_KEYS.index(mode_key) if mode_key in _MODE_KEYS else 0
        )
        refs["density"].setValue(float(d.get("density", 0.0)))
        refs["_pending_iso"] = d.get("isotope_overrides", {})
        self._update_isotope_table(n)

    # ══════════════════════════════════════════════════════════════════
    #  Export (results CSV, scan CSV/HDF5)
    # ══════════════════════════════════════════════════════════════════

    def _export_csv(self) -> None:
        if self._comp1 is None or self._contrast is None:
            QMessageBox.information(self, "No data", "Run Calculate first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Results CSV",
            str(Path.home() / "contrast_results.csv"),
            "CSV Files (*.csv)",
        )
        if not path:
            return
        try:
            results = {
                "comp1": self._comp1.__dict__,
                "comp2": self._comp2.__dict__,
                "contrast": self._contrast.__dict__,
                "anomalous": {"energy_keV": self._energy_spin.value()},
            }
            export_results_csv(results, Path(path))
            self._set_status(f"Results exported to {Path(path).name}")
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))

    def _export_scan_csv(self) -> None:
        if self._scan_data is None:
            QMessageBox.information(self, "No scan data", "Run Energy Scan first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Scan CSV",
            str(Path.home() / "contrast_scan.csv"),
            "CSV Files (*.csv)",
        )
        if not path:
            return
        try:
            export_scan_csv(self._scan_data, Path(path))
            self._set_status(f"Scan CSV exported to {Path(path).name}")
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))

    def _save_scan_hdf5(self) -> None:
        if self._scan_data is None:
            QMessageBox.information(self, "No scan data", "Run Energy Scan first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Scan HDF5",
            str(Path.home() / "contrast_scan.h5"),
            "HDF5 Files (*.h5 *.hdf5)",
        )
        if not path:
            return
        try:
            meta = {
                "comp1_name": self._comp1.name if self._comp1 else "",
                "comp2_name": self._comp2.name if self._comp2 else "",
                "thickness_mm": self._thick_spin.value(),
                "vol_frac_comp1": self._vf_spin.value(),
            }
            save_scan_to_hdf5(self._scan_data, Path(path), metadata=meta)

            # Write compound definitions and all computed properties
            import h5py
            with h5py.File(path, "a") as f:
                grp = f["contrast_scan"]
                _prop_fields = [
                    "mol_weight", "weight_1mol", "n_mol_per_cm3",
                    "n_electrons_per_mol", "n_electrons_per_cm3", "volume_1mol",
                    "xray_sld", "xray_sld_per_gram",
                    "neutron_total_b", "neutron_sld", "neutron_sld_per_gram",
                ]
                for label, comp in [("compound_1", self._comp1),
                                     ("compound_2", self._comp2)]:
                    if comp is None:
                        continue
                    cg = grp.require_group(label)
                    cg.attrs["name"]             = comp.name
                    cg.attrs["formula_str"]       = comp.formula_str
                    cg.attrs["composition_mode"]  = comp.composition_mode
                    cg.attrs["density_g_cm3"]     = comp.density
                    for field in _prop_fields:
                        val = getattr(comp, field, None)
                        if val is not None:
                            cg.attrs[field] = float(val)

            self._set_status(f"Scan + compound data saved to {Path(path).name}")
        except Exception as exc:
            QMessageBox.critical(self, "Save failed", str(exc))

    # ══════════════════════════════════════════════════════════════════
    #  Graph window
    # ══════════════════════════════════════════════════════════════════

    def _show_graph_window(self) -> None:
        if self._graph_win is None:
            self._graph_win = ContrastGraphWindow()
        self._graph_win.show()
        self._graph_win.raise_()
        self._graph_win.activateWindow()

    # ══════════════════════════════════════════════════════════════════
    #  State save / restore
    # ══════════════════════════════════════════════════════════════════

    def _restore_state(self) -> None:
        if self._state_mgr is None:
            return
        st = self._state_mgr.state.get("contrast", {})
        if not st:
            return
        self._apply_compound_state(1, st.get("comp1", {}))
        self._apply_compound_state(2, st.get("comp2", {}))
        self._energy_spin.setValue(float(st.get("energy_keV", 12.0)))
        self._thick_spin.setValue(float(st.get("thickness_mm", 1.0)))
        self._vf_spin.setValue(float(st.get("vol_frac_comp1", 0.01)))
        self._e_start.setValue(float(st.get("e_start_keV", 5.0)))
        self._e_end.setValue(float(st.get("e_end_keV", 25.0)))
        self._n_pts.setValue(int(st.get("n_points", 500)))

    def _apply_compound_state(self, n: int, cd: Dict[str, Any]) -> None:
        if not cd:
            return
        refs = self._refs(n)
        refs["name"].setText(cd.get("name", ""))
        refs["formula"].setText(cd.get("formula_str", ""))
        mode_key = cd.get("composition_mode", "atomic_ratio")
        refs["mode"].setCurrentIndex(
            _MODE_KEYS.index(mode_key) if mode_key in _MODE_KEYS else 0
        )
        refs["density"].setValue(float(cd.get("density", 0.0)))
        refs["_pending_iso"] = cd.get("isotope_overrides", {})
        self._update_isotope_table(n)

    def _save_state(self) -> None:
        if self._state_mgr is None:
            return
        self._state_mgr.state["contrast"] = {
            "comp1": self._get_compound_spec(1),
            "comp2": self._get_compound_spec(2),
            "energy_keV": self._energy_spin.value(),
            "thickness_mm": self._thick_spin.value(),
            "vol_frac_comp1": self._vf_spin.value(),
            "e_start_keV": self._e_start.value(),
            "e_end_keV": self._e_end.value(),
            "n_points": self._n_pts.value(),
        }
        self._state_mgr.save()

    def closeEvent(self, event) -> None:
        self._save_state()
        if self._graph_win is not None:
            self._graph_win.close()
        super().closeEvent(event)

    # ══════════════════════════════════════════════════════════════════
    #  Status bar
    # ══════════════════════════════════════════════════════════════════

    def _set_status(self, msg: str, error: bool = False) -> None:
        self._status_lbl.setText(msg)
        color = "#c0392b" if error else "#2c3e50"
        self._status_lbl.setStyleSheet(
            f"color:{color}; padding:3px;"
            + (" font-weight:bold;" if error else "")
        )


# ─── CLI entry point ───────────────────────────────────────────────────────────

def main() -> None:
    """Entry point for ``pyirena-contrast`` CLI command."""
    from pyirena.state.state_manager import StateManager

    app = QApplication.instance() or QApplication(sys.argv)
    state_manager = StateManager()
    win = ContrastPanel(state_manager=state_manager)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
