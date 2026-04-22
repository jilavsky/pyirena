"""
Diffraction Lines tab for the WAXS Peak Fit panel.

Lets the user import CIF files, computes their theoretical powder diffraction
stick patterns via Dans_Diffraction, and emits a list of visible patterns to
the WAXS graph window for overlay on the experimental I(Q) curve.
"""
from __future__ import annotations

import webbrowser
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from PySide6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QCheckBox,
        QDoubleSpinBox, QFileDialog, QMessageBox, QScrollArea, QFrame,
        QGroupBox, QColorDialog, QMenu, QSizePolicy,
    )
    from PySide6.QtCore import Qt, Signal
    from PySide6.QtGui import QColor
except ImportError:
    from PyQt6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QCheckBox,
        QDoubleSpinBox, QFileDialog, QMessageBox, QScrollArea, QFrame,
        QGroupBox, QColorDialog, QMenu, QSizePolicy,
    )
    from PyQt6.QtCore import Qt, pyqtSignal as Signal
    from PyQt6.QtGui import QColor

from pyirena.core.diffraction_lines import DiffractionPattern, compute_pattern


# ── colour palette: visually distinct phases ──────────────────────────────
_CIF_COLORS = [
    "#2980b9", "#e74c3c", "#27ae60", "#8e44ad", "#e67e22",
    "#16a085", "#c0392b", "#2c3e50", "#d35400", "#7f8c8d",
]


def _next_color(used_colors: List[str]) -> str:
    """Pick the next unused color, recycling if all are taken."""
    for c in _CIF_COLORS:
        if c not in used_colors:
            return c
    return _CIF_COLORS[len(used_colors) % len(_CIF_COLORS)]


# AMCSD search page (free CIFs from American Mineralogist)
AMCSD_URL = "https://rruff.geo.arizona.edu/AMS/amcsd.php"
# COD (Crystallography Open Database) search page
COD_URL = "https://www.crystallography.net/cod/search.html"


# =========================================================================
# Per-CIF row widget
# =========================================================================

class _CifRowWidget(QFrame):
    """One CIF entry in the list: visibility, color, name, scale, hkl, delete."""

    changed = Signal()      # any control changed → recompute/redraw
    delete_requested = Signal()

    def __init__(self, entry: Dict, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.setFrameShape(QFrame.Shape.StyledPanel)
        self.entry = entry  # mutated in place when controls change
        self._build_ui()
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self._on_context_menu)

    def _build_ui(self) -> None:
        row = QHBoxLayout(self)
        row.setContentsMargins(4, 2, 4, 2)
        row.setSpacing(4)

        self._visible_chk = QCheckBox()
        self._visible_chk.setChecked(bool(self.entry.get("visible", True)))
        self._visible_chk.setToolTip("Show / hide this phase")
        self._visible_chk.stateChanged.connect(self._on_visible)
        row.addWidget(self._visible_chk)

        self._color_btn = QPushButton()
        self._color_btn.setFixedSize(20, 20)
        self._color_btn.setToolTip("Click to change colour")
        self._color_btn.clicked.connect(self._on_color)
        self._update_color_swatch()
        row.addWidget(self._color_btn)

        self._name_lbl = QLabel(self.entry.get("name", "?"))
        self._name_lbl.setToolTip(self.entry.get("path", ""))
        self._name_lbl.setMinimumWidth(110)
        self._name_lbl.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred
        )
        row.addWidget(self._name_lbl, 1)

        self._scale_spin = QDoubleSpinBox()
        self._scale_spin.setRange(0.01, 10.0)
        self._scale_spin.setDecimals(2)
        self._scale_spin.setSingleStep(0.1)
        self._scale_spin.setValue(float(self.entry.get("scale", 1.0)))
        self._scale_spin.setFixedWidth(60)
        self._scale_spin.setToolTip("Manual scale (× auto-scale to data peak)")
        self._scale_spin.valueChanged.connect(self._on_scale)
        row.addWidget(self._scale_spin)

        self._hkl_chk = QCheckBox("hkl")
        self._hkl_chk.setChecked(bool(self.entry.get("show_hkl", False)))
        self._hkl_chk.setToolTip("Show Miller-index labels above sticks")
        self._hkl_chk.stateChanged.connect(self._on_hkl)
        row.addWidget(self._hkl_chk)

        self._del_btn = QPushButton("×")
        self._del_btn.setFixedSize(22, 22)
        self._del_btn.setToolTip("Remove this CIF")
        self._del_btn.setStyleSheet(
            "QPushButton{background:#e74c3c;color:white;font-weight:bold;border-radius:3px;}"
            "QPushButton:hover{background:#c0392b;}"
        )
        self._del_btn.clicked.connect(self.delete_requested.emit)
        row.addWidget(self._del_btn)

    def _update_color_swatch(self) -> None:
        col = self.entry.get("color", "#2980b9")
        self._color_btn.setStyleSheet(
            f"QPushButton{{background:{col};border:1px solid #555;border-radius:3px;}}"
        )

    # ── slots ────────────────────────────────────────────────────────────
    def _on_visible(self) -> None:
        self.entry["visible"] = self._visible_chk.isChecked()
        self.changed.emit()

    def _on_scale(self, v: float) -> None:
        self.entry["scale"] = float(v)
        self.changed.emit()

    def _on_hkl(self) -> None:
        self.entry["show_hkl"] = self._hkl_chk.isChecked()
        self.changed.emit()

    def _on_color(self) -> None:
        cur = QColor(self.entry.get("color", "#2980b9"))
        col = QColorDialog.getColor(cur, self, "Pick phase colour")
        if col.isValid():
            self.entry["color"] = col.name()
            self._update_color_swatch()
            self.changed.emit()

    def _on_context_menu(self, pos) -> None:
        menu = QMenu(self)
        act = menu.addAction("Delete")
        chosen = menu.exec(self.mapToGlobal(pos))
        if chosen is act:
            self.delete_requested.emit()


# =========================================================================
# Main panel
# =========================================================================

class DiffractionLinesPanel(QWidget):
    """Tab content: wavelength control, CIF list, import/database buttons."""

    # Emitted whenever the visible patterns or their visual state changes.
    # Payload: list of dicts {pattern: DiffractionPattern, color, scale, show_hkl}
    patterns_changed = Signal(list)

    def __init__(self, state_mgr, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._state_mgr = state_mgr
        self._entries: List[Dict] = []
        self._row_widgets: List[_CifRowWidget] = []
        # Cache: (resolved_path, wavelength_a) -> DiffractionPattern
        self._pattern_cache: Dict[Tuple[str, float], DiffractionPattern] = {}
        # Q-range for pattern computation, populated when data is loaded
        self._q_min: float = 0.1
        self._q_max: float = 10.0

        self._build_ui()
        self._apply_state(self._state_mgr.get("diffraction_lines", default={}))

    # ─────────────────────────────────────────────────────────────────────
    # UI
    # ─────────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.setSpacing(6)

        # ── Wavelength row ───────────────────────────────────────────────
        wl_row = QHBoxLayout()
        wl_row.addWidget(QLabel("λ (Å):"))
        self._wl_spin = QDoubleSpinBox()
        self._wl_spin.setRange(0.01, 20.0)
        self._wl_spin.setDecimals(4)
        self._wl_spin.setSingleStep(0.001)
        self._wl_spin.setValue(1.5406)
        self._wl_spin.setFixedWidth(90)
        self._wl_spin.setToolTip("X-ray wavelength used to compute peak Q-positions")
        self._wl_spin.valueChanged.connect(self._on_wavelength_changed)
        wl_row.addWidget(self._wl_spin)

        self._wl_auto_chk = QCheckBox("Auto from file")
        self._wl_auto_chk.setChecked(True)
        self._wl_auto_chk.setToolTip(
            "When checked, wavelength is taken from the loaded NXcanSAS file's "
            "/entry/instrument/wavelength when available."
        )
        self._wl_auto_chk.stateChanged.connect(self._on_wl_auto_toggled)
        wl_row.addWidget(self._wl_auto_chk)
        wl_row.addStretch()
        outer.addLayout(wl_row)

        # ── CIF list (scrollable) ────────────────────────────────────────
        list_box = QGroupBox("CIF files")
        list_layout = QVBoxLayout(list_box)
        list_layout.setSpacing(2)
        list_layout.setContentsMargins(4, 4, 4, 4)

        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setMinimumHeight(160)
        self._rows_container = QWidget()
        self._rows_vbox = QVBoxLayout(self._rows_container)
        self._rows_vbox.setSpacing(2)
        self._rows_vbox.setContentsMargins(2, 2, 2, 2)
        self._rows_vbox.addStretch()
        self._scroll.setWidget(self._rows_container)
        list_layout.addWidget(self._scroll)

        self._empty_lbl = QLabel(
            "No CIFs loaded. Use Import CIF… or download one from AMCSD / COD."
        )
        self._empty_lbl.setStyleSheet("color:#7f8c8d;font-style:italic;padding:4px;")
        self._empty_lbl.setWordWrap(True)
        list_layout.addWidget(self._empty_lbl)

        outer.addWidget(list_box)

        # ── Action buttons ───────────────────────────────────────────────
        btn_row = QHBoxLayout()
        self._import_btn = QPushButton("Import CIF…")
        self._import_btn.setStyleSheet(
            "QPushButton{background:#3498db;color:white;font-weight:bold;}"
        )
        self._import_btn.clicked.connect(self._import_cif)
        btn_row.addWidget(self._import_btn)

        self._amcsd_btn = QPushButton("AMCSD…")
        self._amcsd_btn.setToolTip(
            f"Open the American Mineralogist Crystal Structure Database in your browser ({AMCSD_URL})"
        )
        self._amcsd_btn.clicked.connect(lambda: webbrowser.open(AMCSD_URL))
        btn_row.addWidget(self._amcsd_btn)

        self._cod_btn = QPushButton("COD…")
        self._cod_btn.setToolTip(
            f"Open the Crystallography Open Database in your browser ({COD_URL})"
        )
        self._cod_btn.clicked.connect(lambda: webbrowser.open(COD_URL))
        btn_row.addWidget(self._cod_btn)
        outer.addLayout(btn_row)

        # ── Reset button ─────────────────────────────────────────────────
        reset_row = QHBoxLayout()
        self._reset_btn = QPushButton("Clear all CIFs")
        self._reset_btn.setStyleSheet(
            "QPushButton{background:#e67e22;color:white;}"
        )
        self._reset_btn.clicked.connect(self._clear_all_with_confirm)
        reset_row.addWidget(self._reset_btn)
        reset_row.addStretch()
        outer.addLayout(reset_row)

        outer.addStretch()

    # ─────────────────────────────────────────────────────────────────────
    # Public API used by WAXSPeakFitPanel
    # ─────────────────────────────────────────────────────────────────────

    def set_data_q_range(self, q_min: float, q_max: float) -> None:
        """Inform panel of the loaded data's Q-range (used to clip patterns)."""
        # Slightly extend range so reflections at the edge are included
        new_min = max(0.01, float(q_min) * 0.95)
        new_max = float(q_max) * 1.05
        if new_min == self._q_min and new_max == self._q_max:
            return
        self._q_min, self._q_max = new_min, new_max
        # Q-range invalidates cached patterns (since we filtered by Q-range)
        self._pattern_cache.clear()
        self._recompute_and_emit()

    def set_wavelength_from_data(self, wavelength_a: Optional[float]) -> None:
        """Apply wavelength loaded from data file when 'Auto' is enabled.

        Pass None when no wavelength was found in the file.
        """
        if wavelength_a is None or wavelength_a <= 0:
            self._wl_auto_chk.setEnabled(False)
            return
        self._wl_auto_chk.setEnabled(True)
        if self._wl_auto_chk.isChecked():
            # Block signals to avoid double-emit; we recompute explicitly
            self._wl_spin.blockSignals(True)
            self._wl_spin.setValue(float(wavelength_a))
            self._wl_spin.blockSignals(False)
            self._pattern_cache.clear()
            self._recompute_and_emit()

    def clear_all(self) -> None:
        """Remove all CIF entries (used by WAXS panel reset-to-defaults)."""
        self._entries.clear()
        for w in self._row_widgets:
            w.setParent(None)
            w.deleteLater()
        self._row_widgets.clear()
        self._pattern_cache.clear()
        self._update_empty_label()
        self._save_state()
        self.patterns_changed.emit([])

    def collect_state(self) -> Dict:
        """Return current state dict for persistence."""
        return {
            "schema_version": 1,
            "wavelength_a": float(self._wl_spin.value()),
            "wavelength_auto": bool(self._wl_auto_chk.isChecked()),
            "last_folder": self._state_mgr.get("diffraction_lines", default={}).get(
                "last_folder"
            ),
            "cif_files": [
                {k: e[k] for k in ("path", "name", "color", "visible", "scale", "show_hkl")}
                for e in self._entries
            ],
        }

    # ─────────────────────────────────────────────────────────────────────
    # State load / save
    # ─────────────────────────────────────────────────────────────────────

    def _apply_state(self, state: Dict) -> None:
        if not state:
            return
        wl = state.get("wavelength_a")
        if wl is not None:
            self._wl_spin.blockSignals(True)
            self._wl_spin.setValue(float(wl))
            self._wl_spin.blockSignals(False)
        if "wavelength_auto" in state:
            self._wl_auto_chk.blockSignals(True)
            self._wl_auto_chk.setChecked(bool(state["wavelength_auto"]))
            self._wl_auto_chk.blockSignals(False)

        for cif_entry in state.get("cif_files", []) or []:
            path = cif_entry.get("path")
            if not path or not Path(path).is_file():
                continue   # silently drop entries whose CIFs vanished
            self._add_entry(
                path=path,
                name=cif_entry.get("name", Path(path).stem),
                color=cif_entry.get("color", _next_color([e["color"] for e in self._entries])),
                visible=bool(cif_entry.get("visible", True)),
                scale=float(cif_entry.get("scale", 1.0)),
                show_hkl=bool(cif_entry.get("show_hkl", False)),
                emit=False,
            )
        self._update_empty_label()
        # Initial emit so graph shows persisted patterns at startup
        self._recompute_and_emit()

    def _save_state(self) -> None:
        self._state_mgr.update("diffraction_lines", self.collect_state())
        # Don't call save() here — main panel saves on its own cadence

    # ─────────────────────────────────────────────────────────────────────
    # Adding / removing entries
    # ─────────────────────────────────────────────────────────────────────

    def _add_entry(
        self,
        path: str,
        name: Optional[str] = None,
        color: Optional[str] = None,
        visible: bool = True,
        scale: float = 1.0,
        show_hkl: bool = False,
        emit: bool = True,
    ) -> None:
        used_colors = [e["color"] for e in self._entries]
        entry = {
            "path": str(path),
            "name": name or Path(path).stem,
            "color": color or _next_color(used_colors),
            "visible": visible,
            "scale": scale,
            "show_hkl": show_hkl,
        }
        # If name not provided, try to use formula from CIF
        if name is None:
            try:
                p = self._compute_cached(entry["path"], float(self._wl_spin.value()))
                if p is not None and p.name:
                    entry["name"] = p.name
            except Exception:
                pass

        self._entries.append(entry)
        row = _CifRowWidget(entry)
        row.changed.connect(self._on_row_changed)
        row.delete_requested.connect(lambda r=row: self._remove_row(r))
        self._row_widgets.append(row)
        # Insert before the trailing stretch
        self._rows_vbox.insertWidget(self._rows_vbox.count() - 1, row)
        self._update_empty_label()
        if emit:
            self._save_state()
            self._recompute_and_emit()

    def _remove_row(self, row: _CifRowWidget) -> None:
        try:
            i = self._row_widgets.index(row)
        except ValueError:
            return
        self._row_widgets.pop(i)
        self._entries.pop(i)
        row.setParent(None)
        row.deleteLater()
        self._update_empty_label()
        self._save_state()
        self._recompute_and_emit()

    def _update_empty_label(self) -> None:
        self._empty_lbl.setVisible(len(self._entries) == 0)

    def _clear_all_with_confirm(self) -> None:
        if not self._entries:
            return
        ans = QMessageBox.question(
            self, "Clear all CIFs",
            f"Remove all {len(self._entries)} CIF file(s) from the list?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if ans == QMessageBox.StandardButton.Yes:
            self.clear_all()

    # ─────────────────────────────────────────────────────────────────────
    # Slots
    # ─────────────────────────────────────────────────────────────────────

    def _on_row_changed(self) -> None:
        self._save_state()
        self._recompute_and_emit()

    def _on_wavelength_changed(self, _v: float) -> None:
        # Manual change → disable auto so file load doesn't overwrite it
        if self._wl_auto_chk.isChecked():
            self._wl_auto_chk.blockSignals(True)
            self._wl_auto_chk.setChecked(False)
            self._wl_auto_chk.blockSignals(False)
        self._pattern_cache.clear()
        self._save_state()
        self._recompute_and_emit()

    def _on_wl_auto_toggled(self) -> None:
        self._save_state()
        # If user just re-enabled auto, the next set_data() call will overwrite

    def _import_cif(self) -> None:
        st = self._state_mgr.get("diffraction_lines", default={})
        last_folder = st.get("last_folder") or str(Path.home())
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Import CIF file(s)", last_folder,
            "Crystallographic Information File (*.cif);;All files (*)",
        )
        if not paths:
            return
        # Remember folder
        st["last_folder"] = str(Path(paths[0]).parent)
        self._state_mgr.update("diffraction_lines", st)

        added = 0
        for p in paths:
            try:
                # Pre-compute to validate; failures get reported
                self._compute_cached(p, float(self._wl_spin.value()))
            except Exception as exc:
                QMessageBox.warning(
                    self, "CIF import failed",
                    f"Could not parse {Path(p).name}:\n{exc}",
                )
                continue
            self._add_entry(p, emit=False)
            added += 1
        if added:
            self._save_state()
            self._recompute_and_emit()

    # ─────────────────────────────────────────────────────────────────────
    # Pattern computation + emission
    # ─────────────────────────────────────────────────────────────────────

    def _compute_cached(
        self, cif_path: str, wavelength_a: float
    ) -> Optional[DiffractionPattern]:
        key = (str(cif_path), float(wavelength_a))
        cached = self._pattern_cache.get(key)
        if cached is not None:
            return cached
        pat = compute_pattern(
            cif_path,
            wavelength_a=wavelength_a,
            q_min=self._q_min,
            q_max=self._q_max,
        )
        self._pattern_cache[key] = pat
        return pat

    def _recompute_and_emit(self) -> None:
        wl = float(self._wl_spin.value())
        out: List[Dict] = []
        for entry in self._entries:
            if not entry.get("visible", True):
                continue
            try:
                pat = self._compute_cached(entry["path"], wl)
            except Exception:
                # Silently skip — CIF may have moved or wavelength yields no peaks
                continue
            if pat is None:
                continue
            out.append({
                "pattern": pat,
                "color": entry["color"],
                "scale": float(entry.get("scale", 1.0)),
                "show_hkl": bool(entry.get("show_hkl", False)),
                "name": entry.get("name", pat.name),
            })
        self.patterns_changed.emit(out)
