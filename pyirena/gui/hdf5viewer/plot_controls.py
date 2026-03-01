"""
PlotControlsPanel — tabbed control panel for the HDF5 Viewer.

Tab 1 — "1D Graph":
  • Pyirena presets: checkboxes for known data types
  • Custom data: X/Y/Yerr/Xerr slots filled by the HDF5 browser
  • Source: first-selected file or all selected files
  • Buttons: New Graph, Add to active graph

Tab 2 — "Collect Values":
  • Choose what to collect (Unified Fit / WAXS / Simple Fit / Custom path)
  • Choose X axis (file order, sort key, HDF5 metadata path)
  • Button: Collect from all selected files → opens CollectWindow
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

try:
    from PySide6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QTabWidget, QGroupBox,
        QCheckBox, QLabel, QLineEdit, QPushButton, QRadioButton,
        QButtonGroup, QComboBox, QSpinBox, QGridLayout, QFrame,
        QMessageBox, QSizePolicy,
    )
    from PySide6.QtCore import Qt, Signal
except ImportError:
    from PyQt6.QtWidgets import (  # type: ignore[no-redef]
        QWidget, QVBoxLayout, QHBoxLayout, QTabWidget, QGroupBox,
        QCheckBox, QLabel, QLineEdit, QPushButton, QRadioButton,
        QButtonGroup, QComboBox, QSpinBox, QGridLayout, QFrame,
        QMessageBox, QSizePolicy,
    )
    from PyQt6.QtCore import Qt, pyqtSignal as Signal  # type: ignore[no-redef]

from . import pyirena_readers as _readers


def _label(text: str, bold: bool = False) -> QLabel:
    lbl = QLabel(text)
    if bold:
        lbl.setStyleSheet("font-weight:bold;")
    return lbl


class PlotControlsPanel(QWidget):
    """
    Tabbed panel providing plot and collect-values controls.

    Signals (all handled by HDF5ViewerWindow)
    -----------------------------------------
    new_graph_requested(list[dict])
        Open a fresh GraphWindow and add these curves to it.
        Each dict: {label, x, y, yerr, xerr, suggest_log_x, suggest_log_y}

    add_to_active_graph_requested(list[dict])
        Add curves to the currently active GraphWindow.

    collect_requested(spec, x_spec, label_spec)
        Collect 0D values from all selected files.
        spec: dict describing what to collect (passed to pyirena_readers.collect_value)
        x_spec: dict describing the X axis
        label_spec: dict with y_label, x_label, title

    status_message(str)
        A status message for the main window status bar.
    """

    new_graph_requested        = Signal(list)
    add_to_active_graph_requested = Signal(list)
    collect_requested          = Signal(dict, dict, dict)
    status_message             = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._selected_files: list[str] = []

        # Slots for custom X/Y/Yerr/Xerr dataset paths (from HDF5 browser)
        self._slot_x    = ""
        self._slot_y    = ""
        self._slot_yerr = ""
        self._slot_xerr = ""

        self._build_ui()

    # ── UI construction ────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        self._tabs = QTabWidget()
        layout.addWidget(self._tabs)

        self._tabs.addTab(self._build_graph_tab(), "1D Graph")
        self._tabs.addTab(self._build_collect_tab(), "Collect Values")

    def _build_graph_tab(self) -> QWidget:
        w = QWidget()
        vl = QVBoxLayout(w)
        vl.setContentsMargins(4, 4, 4, 4)
        vl.setSpacing(6)

        # ── Pyirena presets ───────────────────────────────────────────────
        presets_box = QGroupBox("Pyirena presets")
        pb = QGridLayout(presets_box)
        pb.setHorizontalSpacing(8)
        pb.setVerticalSpacing(3)

        self._cb_nxcansas    = QCheckBox("NXcanSAS  (I vs Q)")
        self._cb_unified     = QCheckBox("Unified Fit model")
        self._cb_sizes_iq    = QCheckBox("Size Dist. model (I vs Q)")
        self._cb_sizes_pr    = QCheckBox("Size Dist.  (P(r) vs r)")
        self._cb_waxs        = QCheckBox("WAXS model")
        self._cb_simple      = QCheckBox("Simple Fit model")

        self._cb_nxcansas.setChecked(True)
        self._cb_unified.setChecked(True)

        pb.addWidget(self._cb_nxcansas, 0, 0)
        pb.addWidget(self._cb_unified,  0, 1)
        pb.addWidget(self._cb_sizes_iq, 1, 0)
        pb.addWidget(self._cb_sizes_pr, 1, 1)
        pb.addWidget(self._cb_waxs,     2, 0)
        pb.addWidget(self._cb_simple,   2, 1)
        vl.addWidget(presets_box)

        # ── Custom data from HDF5 browser ─────────────────────────────────
        custom_box = QGroupBox("Custom data  (from HDF5 browser right-click)")
        cg = QGridLayout(custom_box)
        cg.setHorizontalSpacing(6)
        cg.setVerticalSpacing(3)

        self._slot_labels = {}
        self._slot_edits  = {}
        for row, (role, label) in enumerate([
            ("x",    "X:"),
            ("y",    "Y:"),
            ("yerr", "Y error:"),
            ("xerr", "X error:"),
        ]):
            lbl = QLabel(label)
            cg.addWidget(lbl, row, 0)
            edit = QLabel("(none)")
            edit.setStyleSheet(
                "font-size:9pt; color:#555; background:#f0f0f0; "
                "padding:1px 4px; border:1px solid #ccc; border-radius:2px;"
            )
            edit.setWordWrap(False)
            edit.setMaximumWidth(240)
            cg.addWidget(edit, row, 1)
            self._slot_edits[role] = edit

        clear_btn = QPushButton("Clear")
        clear_btn.setFixedWidth(55)
        clear_btn.clicked.connect(self._clear_slots)
        cg.addWidget(clear_btn, 4, 1, Qt.AlignmentFlag.AlignRight)
        vl.addWidget(custom_box)

        # ── Source selection ──────────────────────────────────────────────
        src_box = QGroupBox("Source")
        sb = QHBoxLayout(src_box)
        self._rb_first = QRadioButton("First selected file")
        self._rb_all   = QRadioButton("All selected files")
        self._rb_all.setChecked(True)
        sb.addWidget(self._rb_first)
        sb.addWidget(self._rb_all)
        vl.addWidget(src_box)

        # ── Buttons ───────────────────────────────────────────────────────
        btn_row = QHBoxLayout()
        self._new_graph_btn = QPushButton("New Graph")
        self._new_graph_btn.setStyleSheet(
            "background:#2980b9; color:white; font-weight:bold;"
        )
        self._new_graph_btn.setMinimumHeight(30)
        self._new_graph_btn.clicked.connect(self._on_new_graph)

        self._add_btn = QPushButton("Add to active graph")
        self._add_btn.setStyleSheet(
            "background:#27ae60; color:white; font-weight:bold;"
        )
        self._add_btn.setMinimumHeight(30)
        self._add_btn.clicked.connect(self._on_add_to_active)

        btn_row.addWidget(self._new_graph_btn)
        btn_row.addWidget(self._add_btn)
        vl.addLayout(btn_row)
        vl.addStretch(1)
        return w

    def _build_collect_tab(self) -> QWidget:
        w = QWidget()
        vl = QVBoxLayout(w)
        vl.setContentsMargins(4, 4, 4, 4)
        vl.setSpacing(6)

        # ── What to collect ───────────────────────────────────────────────
        what_box = QGroupBox("What to collect")
        wg = QGridLayout(what_box)
        wg.setHorizontalSpacing(6)
        wg.setVerticalSpacing(4)

        wg.addWidget(_label("Type:"), 0, 0)
        self._collect_type = QComboBox()
        self._collect_type.addItems([
            "Unified Fit",
            "Size Distribution",
            "WAXS Peak Fit",
            "Simple Fits",
            "Custom HDF5 path",
        ])
        self._collect_type.currentIndexChanged.connect(self._update_collect_ui)
        wg.addWidget(self._collect_type, 0, 1, 1, 2)

        self._collect_item_lbl = _label("Item:")
        wg.addWidget(self._collect_item_lbl, 1, 0)
        self._collect_item = QComboBox()
        self._collect_item.setMinimumWidth(160)
        wg.addWidget(self._collect_item, 1, 1)

        self._collect_index_lbl = _label("Level/Peak:")
        wg.addWidget(self._collect_index_lbl, 2, 0)
        self._collect_index = QSpinBox()
        self._collect_index.setRange(1, 20)
        self._collect_index.setValue(1)
        self._collect_index.setFixedWidth(60)
        wg.addWidget(self._collect_index, 2, 1)

        wg.addWidget(_label("Custom path:"), 3, 0)
        self._collect_path = QLineEdit()
        self._collect_path.setPlaceholderText("/entry/...  or @attr")
        wg.addWidget(self._collect_path, 3, 1, 1, 2)

        self._collect_hint = QLabel("← right-click a dataset in HDF5 browser")
        self._collect_hint.setStyleSheet("font-size:8pt; color:#888; font-style:italic;")
        wg.addWidget(self._collect_hint, 4, 0, 1, 3)

        vl.addWidget(what_box)

        # ── X axis ───────────────────────────────────────────────────────
        x_box = QGroupBox("X axis")
        xg = QVBoxLayout(x_box)
        xg.setSpacing(3)

        self._xrb_order   = QRadioButton("File order  (1, 2, 3…)")
        self._xrb_sortkey = QRadioButton("Filename sort key:")
        self._xrb_path    = QRadioButton("HDF5 metadata path:")
        self._xrb_order.setChecked(True)

        self._x_sortkey_combo = QComboBox()
        self._x_sortkey_combo.addItems([
            "Temperature (°C)", "Time (min)", "Order number", "Pressure (PSI)",
        ])
        self._x_sortkey_combo.setEnabled(False)

        self._x_meta_path = QLineEdit()
        self._x_meta_path.setPlaceholderText("/entry/...  or @attr")
        self._x_meta_path.setEnabled(False)

        self._xrb_sortkey.toggled.connect(self._x_sortkey_combo.setEnabled)
        self._xrb_path.toggled.connect(self._x_meta_path.setEnabled)

        row1 = QHBoxLayout()
        row1.addWidget(self._xrb_sortkey)
        row1.addWidget(self._x_sortkey_combo)

        row2 = QHBoxLayout()
        row2.addWidget(self._xrb_path)
        row2.addWidget(self._x_meta_path, 1)

        xg.addWidget(self._xrb_order)
        xg.addLayout(row1)
        xg.addLayout(row2)
        vl.addWidget(x_box)

        # ── Collect button ────────────────────────────────────────────────
        collect_btn = QPushButton("Collect from all selected files")
        collect_btn.setStyleSheet(
            "background:#8e44ad; color:white; font-weight:bold;"
        )
        collect_btn.setMinimumHeight(30)
        collect_btn.clicked.connect(self._on_collect)
        vl.addWidget(collect_btn)
        vl.addStretch(1)

        # Initialise item lists
        self._update_collect_ui(0)
        return w

    # ── Public API ─────────────────────────────────────────────────────────

    def set_selected_files(self, paths: list[str]) -> None:
        self._selected_files = list(paths)

    def set_dataset_role(self, role: str, hdf5_path: str) -> None:
        """Called when the HDF5 browser emits add_dataset_requested."""
        if role == "x":
            self._slot_x = hdf5_path
        elif role == "y":
            self._slot_y = hdf5_path
        elif role == "yerr":
            self._slot_yerr = hdf5_path
        elif role == "xerr":
            self._slot_xerr = hdf5_path

        short = hdf5_path.split("/")[-1] if hdf5_path else "(none)"
        edit = self._slot_edits.get(role)
        if edit:
            edit.setText(short)
            edit.setToolTip(hdf5_path)

    def set_collect_custom_path(self, hdf5_path: str) -> None:
        """Called when the HDF5 browser emits collect_value_requested."""
        self._tabs.setCurrentIndex(1)
        self._collect_type.setCurrentText("Custom HDF5 path")
        self._collect_path.setText(hdf5_path)
        self._update_collect_ui(4)

    def set_x_axis_path(self, hdf5_path: str) -> None:
        """Called when the HDF5 browser emits set_x_axis_path_requested."""
        self._tabs.setCurrentIndex(1)
        self._xrb_path.setChecked(True)
        self._x_meta_path.setText(hdf5_path)
        self._x_meta_path.setEnabled(True)

    # ── Slot management ────────────────────────────────────────────────────

    def _clear_slots(self) -> None:
        self._slot_x = self._slot_y = self._slot_yerr = self._slot_xerr = ""
        for edit in self._slot_edits.values():
            edit.setText("(none)")
            edit.setToolTip("")

    # ── Build curve dicts ─────────────────────────────────────────────────

    def _get_target_files(self) -> list[str]:
        if not self._selected_files:
            return []
        if self._rb_first.isChecked():
            return [self._selected_files[0]]
        return list(self._selected_files)

    def _build_curves(self) -> list[dict]:
        """Build a list of curve dicts for all target files and selected presets."""
        files = self._get_target_files()
        if not files:
            self.status_message.emit("No files selected.")
            return []

        curves = []
        errors = []

        for filepath in files:
            stem = Path(filepath).stem

            # NXcanSAS
            if self._cb_nxcansas.isChecked():
                result = _readers.read_nxcansas(filepath)
                if result:
                    curves.append({
                        "label": stem,
                        "x": result["Q"], "y": result["I"],
                        "yerr": result.get("dI"),
                        "xerr": None,
                        "suggest_log_x": True,
                        "suggest_log_y": True,
                    })
                else:
                    errors.append(f"No NXcanSAS in {stem}")

            # Unified Fit model
            if self._cb_unified.isChecked():
                result = _readers.read_unified_fit(filepath)
                if result:
                    curves.append({
                        "label": f"{stem}  UF model",
                        "x": result["Q"], "y": result["I_model"],
                        "yerr": None, "xerr": None,
                        "suggest_log_x": True,
                        "suggest_log_y": True,
                    })
                else:
                    errors.append(f"No Unified Fit in {stem}")

            # Size Distribution I(Q) model
            if self._cb_sizes_iq.isChecked():
                result = _readers.read_sizes(filepath)
                if result:
                    curves.append({
                        "label": f"{stem}  Sizes I(Q)",
                        "x": result["Q"], "y": result["I_model"],
                        "yerr": None, "xerr": None,
                        "suggest_log_x": True,
                        "suggest_log_y": True,
                    })
                else:
                    errors.append(f"No Size Dist. in {stem}")

            # Size Distribution P(r) — goes to a separate graph (different axes)
            if self._cb_sizes_pr.isChecked():
                result = _readers.read_sizes(filepath)
                if result:
                    yerr = result.get("distribution_std")
                    curves.append({
                        "label": f"{stem}  P(r)",
                        "x": result["r"], "y": result["distribution"],
                        "yerr": yerr, "xerr": None,
                        "suggest_log_x": False,
                        "suggest_log_y": False,
                        "separate_graph": True,   # P(r) vs r lives on its own axes
                    })
                else:
                    errors.append(f"No Size Dist. in {stem}")

            # WAXS model
            if self._cb_waxs.isChecked():
                result = _readers.read_waxs(filepath)
                if result and len(result["Q"]) > 0:
                    curves.append({
                        "label": f"{stem}  WAXS",
                        "x": result["Q"], "y": result["I_fit"],
                        "yerr": None, "xerr": None,
                        "suggest_log_x": False,
                        "suggest_log_y": False,
                    })
                else:
                    errors.append(f"No WAXS results in {stem}")

            # Simple Fit model
            if self._cb_simple.isChecked():
                result = _readers.read_simple_fit(filepath)
                if result:
                    curves.append({
                        "label": f"{stem}  {result.get('model_name', 'Simple Fit')}",
                        "x": result["Q"], "y": result["I_model"],
                        "yerr": None, "xerr": None,
                        "suggest_log_x": True,
                        "suggest_log_y": True,
                    })
                else:
                    errors.append(f"No Simple Fit in {stem}")

            # Custom X/Y
            if self._slot_y:
                y_arr = _readers.read_dataset(filepath, self._slot_y)
                x_arr = None
                if self._slot_x:
                    x_arr = _readers.read_dataset(filepath, self._slot_x)
                yerr_arr = None
                if self._slot_yerr:
                    yerr_arr = _readers.read_dataset(filepath, self._slot_yerr)

                if y_arr is not None:
                    if x_arr is None:
                        x_arr = np.arange(len(y_arr), dtype=float)
                    x_lbl = self._slot_x.split("/")[-1] if self._slot_x else "index"
                    y_lbl = self._slot_y.split("/")[-1]
                    curves.append({
                        "label": f"{stem}  {y_lbl}",
                        "x": x_arr, "y": y_arr,
                        "yerr": yerr_arr, "xerr": None,
                        "suggest_log_x": False,
                        "suggest_log_y": False,
                    })
                else:
                    errors.append(f"Could not read Y data from {stem}")

        if errors and not curves:
            QMessageBox.warning(self, "No data found",
                                "Could not read any data:\n" + "\n".join(errors))
        elif errors:
            self.status_message.emit(
                f"Plotted {len(curves)} curve(s). Warnings: " + "; ".join(errors)
            )
        else:
            self.status_message.emit(f"Plotted {len(curves)} curve(s).")

        return curves

    # ── Button handlers ────────────────────────────────────────────────────

    @staticmethod
    def _split_curves(curves: list[dict]) -> tuple[list[dict], list[dict]]:
        """Split curves into (regular, separate_graph) groups."""
        regular = [c for c in curves if not c.get("separate_graph")]
        separate = [c for c in curves if c.get("separate_graph")]
        return regular, separate

    def _on_new_graph(self) -> None:
        curves = self._build_curves()
        if not curves:
            return
        regular, separate = self._split_curves(curves)
        if regular:
            self.new_graph_requested.emit(regular)
        if separate:
            self.new_graph_requested.emit(separate)

    def _on_add_to_active(self) -> None:
        curves = self._build_curves()
        if not curves:
            return
        regular, separate = self._split_curves(curves)
        if regular:
            self.add_to_active_graph_requested.emit(regular)
        if separate:
            # P(r) always opens its own window — different axes
            self.new_graph_requested.emit(separate)

    def _on_collect(self) -> None:
        if not self._selected_files:
            QMessageBox.warning(self, "No files", "Select files in the file tree first.")
            return

        # Build value spec
        type_text = self._collect_type.currentText()
        item_text = self._collect_item.currentText()
        idx = self._collect_index.value()

        type_map = {
            "Unified Fit":        "unified_fit",
            "Size Distribution":  "sizes",
            "WAXS Peak Fit":      "waxs",
            "Simple Fits":        "simple_fit",
            "Custom HDF5 path":   "custom",
        }
        type_key = type_map.get(type_text, "custom")

        if type_key == "custom":
            spec = {"type": "custom", "path": self._collect_path.text().strip()}
        elif type_key == "unified_fit":
            spec = {"type": "unified_fit", "item": item_text, "level": idx}
        elif type_key == "sizes":
            spec = {"type": "sizes", "item": item_text}
        elif type_key == "waxs":
            spec = {"type": "waxs", "item": item_text, "peak": idx - 1}
        elif type_key == "simple_fit":
            spec = {"type": "simple_fit", "item": "param", "param_name": item_text}
        else:
            spec = {"type": "custom", "path": ""}

        # Build X spec
        if self._xrb_order.isChecked():
            x_spec = {"type": "order"}
        elif self._xrb_sortkey.isChecked():
            sortkey_map = {
                "Temperature (°C)": 2,
                "Time (min)":       4,
                "Order number":     6,
                "Pressure (PSI)":   8,
            }
            x_spec = {
                "type": "sortkey",
                "sort_index": sortkey_map.get(self._x_sortkey_combo.currentText(), 6),
            }
        else:
            x_spec = {"type": "path", "path": self._x_meta_path.text().strip()}

        # Label spec
        label_spec = {
            "y_label": f"{type_text}: {item_text}",
            "x_label": self._x_sortkey_combo.currentText() if self._xrb_sortkey.isChecked()
                       else ("File order" if self._xrb_order.isChecked() else "Metadata"),
            "title":   f"{type_text} — {item_text}",
        }

        self.collect_requested.emit(spec, x_spec, label_spec)

    # ── Collect tab dynamic UI ─────────────────────────────────────────────

    def _update_collect_ui(self, type_index: int) -> None:
        """Update item combo and visibility based on selected type."""
        self._collect_item.clear()
        type_text = self._collect_type.currentText()

        is_custom = (type_text == "Custom HDF5 path")
        needs_index = type_text in ("Unified Fit", "WAXS Peak Fit")

        # Custom path row visibility
        self._collect_path.setVisible(is_custom)
        self._collect_hint.setVisible(is_custom)

        # Item shown for all non-custom types; Level/Peak only for indexed types
        self._collect_item_lbl.setVisible(not is_custom)
        self._collect_item.setVisible(not is_custom)
        self._collect_index_lbl.setVisible(needs_index)
        self._collect_index.setVisible(needs_index)

        if type_text == "Unified Fit":
            items = ["Rg", "G", "B", "P", "ETA", "PACK",
                     "Rg_err", "G_err", "B_err", "P_err",
                     "chi2", "background"]
            self._collect_item.addItems(items)
            self._collect_index.setPrefix("Level ")
            self._collect_index.setRange(1, 5)
            self._collect_index.setEnabled(True)

        elif type_text == "Size Distribution":
            items = ["chi_squared", "volume_fraction", "rg"]
            self._collect_item.addItems(items)
            self._collect_index.setEnabled(False)

        elif type_text == "WAXS Peak Fit":
            items = ["Q0", "A", "FWHM", "Q0_err", "A_err", "FWHM_err", "chi2"]
            self._collect_item.addItems(items)
            self._collect_index.setPrefix("Peak ")
            self._collect_index.setRange(1, 20)
            self._collect_index.setEnabled(True)

        elif type_text == "Simple Fits":
            self._collect_item.addItems(["Rg", "I0", "Rg_err", "I0_err", "chi2"])
            self._collect_index.setEnabled(False)

        # Custom: item/index are hidden so no need to configure them
