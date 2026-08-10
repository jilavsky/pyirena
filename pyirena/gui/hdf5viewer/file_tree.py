"""
FileTreeWidget — collapsible file-tree for browsing HDF5 files in a folder
and its subfolders.

Key design:
- Top-level items: HDF5 files in the root folder + subfolder nodes.
- Subfolder nodes expand lazily (children populated on itemExpanded).
- Each file item stores the full absolute path as Qt.UserRole data so that
  identically-named files in different subfolders are unambiguous.
- Multi-selection is enabled across the whole tree.
- Sorting is driven by the same sort-key functions as the Data Selector.
"""

from __future__ import annotations

import logging

log = logging.getLogger(__name__)


import os
from pathlib import Path
from typing import Callable

from pyirena.core.file_sorting import SORT_LABELS, SORT_TOOLTIP, sort_names
from pyirena.gui._qt import (
    QBrush,
    QColor,
    QComboBox,
    QFileDialog,
    QFont,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSizePolicy,
    Qt,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
    Signal,
)
from pyirena.gui.file_drop import drop_hint, enable_file_drop, first_folder
from pyirena.gui.file_filter import FILTER_PLACEHOLDER, FILTER_TOOLTIP, make_file_matcher

# ── Extensions treated as HDF5 files ───────────────────────────────────────
HDF5_EXTENSIONS = {".h5", ".hdf5", ".hdf", ".nxs", ".h5xp"}

# Sort modes come from the one shared implementation; see
# pyirena.core.file_sorting (they used to be copied into every browser).

# Qt user-data role used to store the absolute file path on file items
_PATH_ROLE = Qt.ItemDataRole.UserRole


class FileTreeWidget(QWidget):
    """
    A folder browser that lists HDF5 files and supports collapsible subfolders.

    Signals
    -------
    selection_changed(list[str])
        Emitted whenever the tree selection changes.
        Carries the absolute paths of all currently selected *file* items.
    folder_changed(str)
        Emitted when the user picks a new root folder.
    """

    selection_changed = Signal(list)
    folder_changed = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._root_folder: str = ""
        self._sort_index: int = 6   # Order number ↑ (same default as data_selector)
        self._filter_text: str = ""
        self._building: bool = False  # guard against recursive expand signals

        self._build_ui()

    # ── UI construction ────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(3)

        # Panel subtitle
        subtitle = QLabel("File Select")
        subtitle.setStyleSheet(
            "font-size:10pt; font-weight:bold; color:#2c3e50;"
            "padding:3px 4px; background:#dfe6e9; border-bottom:1px solid #b2bec3;"
        )
        layout.addWidget(subtitle)

        # Folder selection row
        folder_row = QHBoxLayout()
        self._folder_btn = QPushButton("Select Folder")
        self._folder_btn.setMaximumWidth(100)
        self._folder_btn.clicked.connect(self._select_folder)
        self._folder_label = QLabel("(no folder)")
        self._folder_label.setWordWrap(False)
        self._folder_label.setSizePolicy(
            QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred
        )
        folder_row.addWidget(self._folder_btn)
        folder_row.addWidget(self._folder_label, 1)
        layout.addLayout(folder_row)

        # Sort row
        self._sort_combo = QComboBox()
        self._sort_combo.addItems(SORT_LABELS)
        self._sort_combo.setToolTip(SORT_TOOLTIP)
        self._sort_combo.setCurrentIndex(self._sort_index)
        self._sort_combo.currentIndexChanged.connect(self._on_sort_changed)
        layout.addWidget(self._sort_combo)

        # Filter row (below sort)
        self._filter_edit = QLineEdit()
        self._filter_edit.setPlaceholderText(FILTER_PLACEHOLDER)
        self._filter_edit.setToolTip(FILTER_TOOLTIP)
        self._filter_edit.textChanged.connect(self._on_filter_changed)
        layout.addWidget(self._filter_edit)

        # Count label
        self._count_label = QLabel("")
        self._count_label.setStyleSheet("font-size:9pt; color:#777;")
        layout.addWidget(self._count_label)

        # Tree
        self._tree = QTreeWidget()
        self._tree.setHeaderHidden(True)
        self._tree.setSelectionMode(
            QTreeWidget.SelectionMode.ExtendedSelection
        )
        self._tree.setIndentation(14)
        self._tree.itemExpanded.connect(self._on_item_expanded)
        self._tree.itemSelectionChanged.connect(self._on_selection_changed)
        self._tree.setToolTip(drop_hint("a folder or HDF5 files"))
        layout.addWidget(self._tree, 1)

        # Dropping a folder browses it; dropping files browses their folder and
        # selects them.  Installed on the whole widget so the drop target is the
        # entire left pane, not just the tree rows.
        enable_file_drop(self, self.open_dropped_paths, extensions=tuple(HDF5_EXTENSIONS))

    # ── Public API ─────────────────────────────────────────────────────────

    @property
    def root_folder(self) -> str:
        return self._root_folder

    def set_folder(self, folder: str) -> None:
        """Load a new root folder into the tree."""
        if not folder or not os.path.isdir(folder):
            return
        self._root_folder = folder
        self._folder_label.setText(os.path.basename(folder) or folder)
        self._folder_label.setToolTip(folder)
        self._refresh_tree()
        self.folder_changed.emit(folder)

    def open_dropped_paths(self, paths: list[str]) -> None:
        """Browse what was dropped: switch to its folder and select the files.

        Dropping a folder is the common case (a measurement directory), and
        ``collect_dropped_paths`` has already expanded it to the HDF5 files
        inside, so both gestures land here the same way.
        """
        if not paths:
            return

        folder = first_folder(paths)
        if folder and folder != self._root_folder:
            self.set_folder(folder)

        wanted = {os.path.basename(p) for p in paths
                  if os.path.dirname(p) == self._root_folder}
        if not wanted:
            return

        self._tree.clearSelection()
        matched = 0
        root = self._tree.invisibleRootItem()
        for i in range(root.childCount()):
            item = root.child(i)
            if item.data(0, _PATH_ROLE) and item.text(0) in wanted:
                item.setSelected(True)
                matched += 1
                if matched == 1:
                    self._tree.scrollToItem(item)
        if not matched and self._filter_text:
            # The dropped files exist but the filter hides them; clearing it is
            # less confusing than a drop that appears to do nothing.
            log.info("dropped files hidden by filter %r — clearing it", self._filter_text)
            self.set_filter("")
            self.open_dropped_paths(paths)

    def get_selected_paths(self) -> list[str]:
        """Return absolute paths of all currently selected file items."""
        paths = []
        for item in self._tree.selectedItems():
            path = item.data(0, _PATH_ROLE)
            if path:
                paths.append(str(path))
        return paths

    def set_sort_index(self, idx: int) -> None:
        self._sort_combo.setCurrentIndex(idx)

    def set_filter(self, text: str) -> None:
        self._filter_edit.setText(text)

    # ── Folder selection ───────────────────────────────────────────────────

    def _select_folder(self) -> None:
        start = self._root_folder or os.path.expanduser("~")
        folder = QFileDialog.getExistingDirectory(
            self, "Select Folder", start,
            QFileDialog.Option.ShowDirsOnly,
        )
        if folder:
            self.set_folder(folder)

    # ── Tree building ──────────────────────────────────────────────────────

    def _refresh_tree(self) -> None:
        """Rebuild the tree from the current root folder."""
        self._building = True
        try:
            self._tree.clear()
            if not self._root_folder:
                return
            self._populate_folder(
                parent=self._tree.invisibleRootItem(),
                folder=self._root_folder,
            )
            self._apply_filter()
            self._update_count()
        finally:
            self._building = False

    def _populate_folder(
        self,
        parent: QTreeWidgetItem,
        folder: str,
    ) -> None:
        """
        Add file items and collapsed subfolder nodes to *parent*.
        Subfolders are only included if they contain at least one HDF5 file
        (recursively).
        """
        try:
            entries = os.listdir(folder)
        except PermissionError:
            return

        # Separate files and subdirs
        files = []
        subdirs = []
        for name in entries:
            full = os.path.join(folder, name)
            if os.path.isfile(full) and Path(name).suffix.lower() in HDF5_EXTENSIONS:
                files.append(name)
            elif os.path.isdir(full) and not name.startswith("."):
                # Only include if it (recursively) contains HDF5 files
                if self._dir_has_hdf5(full):
                    subdirs.append(name)

        # Sort files by current sort key
        files = self._sort_names(files)

        # Add subfolder nodes first (collapsed, at top)
        for dirname in sorted(subdirs, key=str.lower):
            dir_item = QTreeWidgetItem(parent, [dirname])
            dir_item.setData(0, _PATH_ROLE, None)   # no file path
            # Bold italic for folder nodes
            font = QFont()
            font.setBold(True)
            font.setItalic(True)
            dir_item.setFont(0, font)
            dir_item.setForeground(0, QBrush(QColor("#2980b9")))
            # Force the expand triangle even before children are loaded
            dir_item.setChildIndicatorPolicy(
                QTreeWidgetItem.ChildIndicatorPolicy.ShowIndicator
            )
            # Store the absolute path so we can populate on expand
            dir_item.setData(0, Qt.ItemDataRole.UserRole + 1, os.path.join(folder, dirname))
            dir_item.setData(0, Qt.ItemDataRole.UserRole + 2, False)  # not yet loaded

        # Add file items
        for name in files:
            full_path = os.path.join(folder, name)
            file_item = QTreeWidgetItem(parent, [name])
            file_item.setData(0, _PATH_ROLE, full_path)
            file_item.setToolTip(0, full_path)

    def _on_item_expanded(self, item: QTreeWidgetItem) -> None:
        """Lazily populate a subfolder node when the user expands it."""
        if self._building:
            return
        already_loaded = item.data(0, Qt.ItemDataRole.UserRole + 2)
        if already_loaded:
            return
        folder = item.data(0, Qt.ItemDataRole.UserRole + 1)
        if not folder:
            return
        # Mark as loaded before populating (prevent re-entry)
        item.setData(0, Qt.ItemDataRole.UserRole + 2, True)
        self._populate_folder(item, folder)
        self._apply_filter()
        self._update_count()

    def _dir_has_hdf5(self, path: str) -> bool:
        """Return True if *path* or any descendant contains an HDF5 file."""
        try:
            for entry in os.scandir(path):
                if entry.is_file() and Path(entry.name).suffix.lower() in HDF5_EXTENSIONS:
                    return True
                if entry.is_dir() and not entry.name.startswith("."):
                    if self._dir_has_hdf5(entry.path):
                        return True
        except PermissionError:
            log.debug("suppressed exception", exc_info=True)
        return False

    # ── Sort / filter ──────────────────────────────────────────────────────

    def _sort_names(self, names: list[str]) -> list[str]:
        return sort_names(names, self._sort_index)

    def _on_sort_changed(self, idx: int) -> None:
        self._sort_index = idx
        self._refresh_tree()

    def _on_filter_changed(self, text: str) -> None:
        self._filter_text = text.strip()
        self._apply_filter()
        self._update_count()

    def _apply_filter(self) -> None:
        """Show/hide file items based on the filter text."""
        # Compile once here rather than per filename inside the recursion.
        matcher = make_file_matcher(self._filter_text)
        root = self._tree.invisibleRootItem()
        self._filter_node(root, matcher)

    def _filter_node(self, parent: QTreeWidgetItem, flt: Callable[[str], bool]) -> bool:
        """Recursively filter items; return True if any child is visible."""
        any_visible = False
        for i in range(parent.childCount()):
            child = parent.child(i)
            path = child.data(0, _PATH_ROLE)
            if path is None:
                # Folder node
                loaded = child.data(0, Qt.ItemDataRole.UserRole + 2)
                if not loaded:
                    # Not expanded yet — assume it has content, keep visible
                    child.setHidden(False)
                    any_visible = True
                else:
                    child_visible = self._filter_node(child, flt)
                    child.setHidden(not child_visible)
                    any_visible |= child_visible
            else:
                # File item
                visible = flt(child.text(0))
                child.setHidden(not visible)
                any_visible |= visible
        return any_visible

    def _update_count(self) -> None:
        visible = self._count_visible(self._tree.invisibleRootItem())
        self._count_label.setText(f"{visible} file(s)")

    def _count_visible(self, parent: QTreeWidgetItem) -> int:
        count = 0
        for i in range(parent.childCount()):
            child = parent.child(i)
            if child.isHidden():
                continue
            if child.data(0, _PATH_ROLE):
                count += 1
            else:
                count += self._count_visible(child)
        return count

    # ── Selection ──────────────────────────────────────────────────────────

    def _on_selection_changed(self) -> None:
        paths = self.get_selected_paths()
        self.selection_changed.emit(paths)
