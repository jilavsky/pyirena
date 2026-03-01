"""
HDF5BrowserWidget — lazy tree browser for HDF5 file contents.

Shows groups and datasets as a collapsible tree.  Children are loaded
from h5py only when a node is expanded (lazy loading) to keep large
files responsive.

Context menu on dataset nodes:  Add as X / Y / Y error / X error axis.
Context menu on group nodes:  Plot NXcanSAS / Unified Fit / … if known type.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

try:
    from PySide6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QTreeWidget, QTreeWidgetItem,
        QLabel, QLineEdit, QMenu, QSizePolicy, QAbstractItemView,
    )
    from PySide6.QtCore import Qt, Signal
    from PySide6.QtGui import QFont, QAction
except ImportError:
    from PyQt6.QtWidgets import (  # type: ignore[no-redef]
        QWidget, QVBoxLayout, QHBoxLayout, QTreeWidget, QTreeWidgetItem,
        QLabel, QLineEdit, QMenu, QSizePolicy, QAbstractItemView,
    )
    from PyQt6.QtCore import Qt, pyqtSignal as Signal  # type: ignore[no-redef]
    from PyQt6.QtGui import QFont, QAction             # type: ignore[no-redef]

import h5py
import numpy as np

# Qt user-data roles
_HDF5_PATH_ROLE = Qt.ItemDataRole.UserRole        # full hdf5 path string
_IS_GROUP_ROLE  = Qt.ItemDataRole.UserRole + 1    # bool: True = group
_LOADED_ROLE    = Qt.ItemDataRole.UserRole + 2    # bool: children loaded?

# Known pyirena result group name suffixes → type key
_KNOWN_GROUPS = {
    "unified_fit_results":  "unified_fit",
    "sizes_results":        "sizes",
    "waxs_peakfit_results": "waxs",
    "simple_fit_results":   "simple_fit",
}


def _classify_group(h5_group: h5py.Group) -> str | None:
    """Return 'nxcansas' if the group looks like NXcanSAS data entry."""
    try:
        css = h5_group.attrs.get("canSAS_class", b"")
        if isinstance(css, bytes):
            css = css.decode("utf-8", errors="ignore")
        if css == "SASentry":
            return "nxcansas"
    except Exception:
        pass
    return None


class HDF5BrowserWidget(QWidget):
    """
    Widget that shows the internal tree of an HDF5 file.

    Signals
    -------
    add_dataset_requested(str role, str hdf5_path)
        User chose "Add as X/Y/Yerr/Xerr" from context menu.
        role ∈ {"x", "y", "yerr", "xerr"}.

    plot_known_type_requested(str type_key, str group_path)
        User chose "Plot NXcanSAS" / "Plot Unified Fit model" / … from context
        menu.  type_key ∈ {"nxcansas", "unified_fit", "sizes", "waxs",
        "simple_fit"}.

    collect_value_requested(str hdf5_path)
        User chose "Collect value across selected files".
    """

    add_dataset_requested   = Signal(str, str)
    plot_known_type_requested = Signal(str, str)
    collect_value_requested = Signal(str)
    set_x_axis_path_requested = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._filepath: str | None = None
        self._h5file:   h5py.File | None = None
        self._build_ui()

    # ── UI construction ────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(3)

        # Panel subtitle
        subtitle = QLabel("Data Select")
        subtitle.setStyleSheet(
            "font-size:10pt; font-weight:bold; color:#2c3e50;"
            "padding:3px 4px; background:#dfe6e9; border-bottom:1px solid #b2bec3;"
        )
        layout.addWidget(subtitle)

        # Prominent filename
        self._file_label = QLabel("(no file selected)")
        self._file_label.setWordWrap(True)
        self._file_label.setStyleSheet(
            "font-size:11pt; font-weight:bold; color:#1a1a1a;"
            "padding:3px 4px;"
        )
        self._file_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        layout.addWidget(self._file_label)

        self._tree = QTreeWidget()
        self._tree.setHeaderLabels(["Name", "Info"])
        self._tree.setColumnWidth(0, 200)
        self._tree.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self._tree.setIndentation(14)
        self._tree.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self._tree.customContextMenuRequested.connect(self._on_context_menu)
        self._tree.itemExpanded.connect(self._on_item_expanded)
        self._tree.currentItemChanged.connect(self._on_current_item_changed)
        layout.addWidget(self._tree, 1)

        # Scalar value display
        self._value_edit = QLineEdit()
        self._value_edit.setReadOnly(True)
        self._value_edit.setPlaceholderText("Select a scalar dataset or attribute to see its value")
        self._value_edit.setStyleSheet("font-size:9pt; color:#333;")
        layout.addWidget(self._value_edit)

    # ── Public API ─────────────────────────────────────────────────────────

    def load_file(self, filepath: str | None) -> None:
        """Show the structure of *filepath*.  Pass None to clear."""
        self._close_h5()
        self._tree.clear()
        self._filepath = filepath
        self._value_edit.clear()

        if not filepath or not os.path.isfile(filepath):
            self._file_label.setText("(no file selected)")
            return

        self._file_label.setText(os.path.basename(filepath))
        self._file_label.setToolTip(filepath)

        try:
            self._h5file = h5py.File(filepath, "r")
        except Exception as exc:
            self._file_label.setText(f"Error opening file: {exc}")
            return

        # Populate root level
        self._populate_group(
            parent=self._tree.invisibleRootItem(),
            h5_node=self._h5file,
            h5_path="/",
        )

    def clear(self) -> None:
        self._close_h5()
        self._tree.clear()
        self._file_label.setText("(no file selected)")
        self._value_edit.clear()

    # ── Tree population ────────────────────────────────────────────────────

    def _populate_group(
        self,
        parent: QTreeWidgetItem,
        h5_node: h5py.Group | h5py.File,
        h5_path: str,
    ) -> None:
        """Add immediate children of *h5_node* to *parent*."""
        if self._h5file is None:
            return

        try:
            names = list(h5_node.keys())
        except Exception:
            return

        for name in sorted(names, key=str.lower):
            child_path = f"{h5_path.rstrip('/')}/{name}"
            try:
                child = h5_node[name]
            except Exception:
                continue

            is_group = isinstance(child, (h5py.Group, h5py.File))
            if is_group:
                info = f"Group  ({len(child)} items)"
                item = QTreeWidgetItem(parent, [name, info])
                item.setData(0, _HDF5_PATH_ROLE, child_path)
                item.setData(0, _IS_GROUP_ROLE, True)
                item.setData(0, _LOADED_ROLE, False)
                # Force triangle even before children are loaded
                if len(child) > 0:
                    item.setChildIndicatorPolicy(
                        QTreeWidgetItem.ChildIndicatorPolicy.ShowIndicator
                    )
                font = QFont()
                font.setBold(True)
                item.setFont(0, font)
            else:
                # Dataset
                ds: h5py.Dataset = child
                shape_str = "×".join(str(d) for d in ds.shape) if ds.shape else "scalar"
                dtype_str = str(ds.dtype)
                # Read scalar string values inline
                val_str = ""
                if ds.shape == () or ds.ndim == 0:
                    try:
                        v = ds[()]
                        if isinstance(v, (bytes, np.bytes_)):
                            v = v.decode("utf-8", errors="replace")
                        val_str = f" = {v}"
                    except Exception:
                        pass
                info = f"{shape_str}  [{dtype_str}]{val_str}"
                item = QTreeWidgetItem(parent, [name, info])
                item.setData(0, _HDF5_PATH_ROLE, child_path)
                item.setData(0, _IS_GROUP_ROLE, False)
                item.setData(0, _LOADED_ROLE, True)

            # Show HDF5 attributes as italic children (for key attrs)
            self._add_key_attrs(item, child, child_path)

    def _add_key_attrs(
        self,
        parent: QTreeWidgetItem,
        h5_node: h5py.HLObject,
        h5_path: str,
    ) -> None:
        """Add a few important attributes as sub-items for discoverability."""
        interesting = {
            "NX_class", "canSAS_class", "signal", "axes", "I_axes",
            "analysis_type", "program", "units", "n_peaks",
        }
        try:
            for key in h5_node.attrs:
                if key not in interesting:
                    continue
                val = h5_node.attrs[key]
                if isinstance(val, (bytes, np.bytes_)):
                    val = val.decode("utf-8", errors="replace")
                attr_item = QTreeWidgetItem(parent, [f"@{key}", str(val)])
                attr_item.setData(0, _HDF5_PATH_ROLE, None)  # not selectable as data
                attr_item.setData(0, _IS_GROUP_ROLE, False)
                font = QFont()
                font.setItalic(True)
                attr_item.setFont(0, font)
                attr_item.setForeground(0, self._tree.palette().color(
                    self._tree.foregroundRole()
                ))
                attr_item.setForeground(1, self._tree.palette().color(
                    self._tree.foregroundRole()
                ))
                attr_item.setFlags(attr_item.flags() & ~Qt.ItemFlag.ItemIsSelectable)
        except Exception:
            pass

    def _on_item_expanded(self, item: QTreeWidgetItem) -> None:
        """Lazily load children of a group item on first expand."""
        if not item.data(0, _IS_GROUP_ROLE):
            return
        if item.data(0, _LOADED_ROLE):
            return
        item.setData(0, _LOADED_ROLE, True)

        h5_path = item.data(0, _HDF5_PATH_ROLE)
        if not h5_path or self._h5file is None:
            return
        try:
            h5_node = self._h5file[h5_path]
        except Exception:
            return

        self._populate_group(item, h5_node, h5_path)

    def _on_current_item_changed(self, current: QTreeWidgetItem | None, _previous) -> None:
        """Show scalar value of the selected item in the value display field."""
        if current is None:
            self._value_edit.clear()
            return
        info = current.text(1)   # "Info" column
        if " = " in info:
            # Scalar dataset: info looks like "scalar  [float64] = 1.23"
            value_str = info.split(" = ", 1)[1]
            name = current.text(0)
            self._value_edit.setText(f"{name} = {value_str}")
        elif current.text(0).startswith("@"):
            # Attribute item: info column is the value
            name = current.text(0)
            self._value_edit.setText(f"{name} = {info}")
        else:
            self._value_edit.clear()

    # ── Context menu ───────────────────────────────────────────────────────

    def _on_context_menu(self, pos) -> None:
        item = self._tree.itemAt(pos)
        if item is None:
            return

        h5_path = item.data(0, _HDF5_PATH_ROLE)
        if h5_path is None:
            return   # attribute sub-item — no menu

        is_group = item.data(0, _IS_GROUP_ROLE)
        menu = QMenu(self)

        if is_group:
            self._build_group_menu(menu, item, h5_path)
        else:
            self._build_dataset_menu(menu, item, h5_path)

        if not menu.isEmpty():
            menu.exec(self._tree.viewport().mapToGlobal(pos))

    def _build_dataset_menu(self, menu: QMenu, item: QTreeWidgetItem, h5_path: str) -> None:
        for role, label in [
            ("y",    "Add as Y axis"),
            ("x",    "Add as X axis"),
            ("yerr", "Add as Y error"),
            ("xerr", "Add as X error"),
        ]:
            act = QAction(label, menu)
            act.triggered.connect(
                lambda checked=False, r=role, p=h5_path:
                self.add_dataset_requested.emit(r, p)
            )
            menu.addAction(act)

        menu.addSeparator()
        act = QAction("Collect value across selected files", menu)
        act.triggered.connect(
            lambda checked=False, p=h5_path: self.collect_value_requested.emit(p)
        )
        menu.addAction(act)

        act_xpath = QAction("Set as X-axis metadata path", menu)
        act_xpath.triggered.connect(
            lambda checked=False, p=h5_path: self.set_x_axis_path_requested.emit(p)
        )
        menu.addAction(act_xpath)

    def _build_group_menu(self, menu: QMenu, item: QTreeWidgetItem, h5_path: str) -> None:
        if self._h5file is None:
            return

        # Check if it is an NXcanSAS entry
        try:
            grp = self._h5file[h5_path]
            nxcansas_type = _classify_group(grp)
        except Exception:
            nxcansas_type = None

        if nxcansas_type:
            act = QAction("Plot NXcanSAS data  (I vs Q)", menu)
            act.triggered.connect(
                lambda checked=False, p=h5_path:
                self.plot_known_type_requested.emit("nxcansas", p)
            )
            menu.addAction(act)

        # Check for known pyirena result groups
        group_name = h5_path.rstrip("/").split("/")[-1]
        type_key = _KNOWN_GROUPS.get(group_name)
        if type_key:
            label_map = {
                "unified_fit": "Plot Unified Fit model",
                "sizes":       "Plot Size Distribution model",
                "waxs":        "Plot WAXS Peak Fit model",
                "simple_fit":  "Plot Simple Fit model",
            }
            label = label_map.get(type_key, f"Plot {type_key}")
            act = QAction(label, menu)
            act.triggered.connect(
                lambda checked=False, tk=type_key, p=h5_path:
                self.plot_known_type_requested.emit(tk, p)
            )
            menu.addAction(act)

    # ── Cleanup ────────────────────────────────────────────────────────────

    def _close_h5(self) -> None:
        if self._h5file is not None:
            try:
                self._h5file.close()
            except Exception:
                pass
            self._h5file = None

    def closeEvent(self, event) -> None:
        self._close_h5()
        super().closeEvent(event)

    def __del__(self) -> None:
        self._close_h5()
