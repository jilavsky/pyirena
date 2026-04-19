"""
data_manipulation_panel.py — GUI panel and CLI entry point for the Data Manipulation tool.

Provides scaling, trimming, rebinning, averaging, subtracting, and dividing
of SAS datasets.

Entry points
------------
* ``DataManipulationPanel`` — standalone QWidget; launched from the Data Selector hub.
* ``main()`` — ``pyirena-datamanip`` CLI command.
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path
from typing import Callable, Optional, List

import numpy as np

try:
    from PySide6.QtWidgets import (
        QApplication, QWidget, QVBoxLayout, QHBoxLayout,
        QPushButton, QLabel, QLineEdit, QComboBox, QCheckBox,
        QListWidget, QMessageBox, QGroupBox, QFrame, QFileDialog,
        QAbstractItemView, QSizePolicy, QListWidgetItem,
        QTabWidget, QSpinBox, QMenu,
    )
    from PySide6.QtCore import Qt, QUrl, QTimer
    from PySide6.QtGui import QDesktopServices, QDoubleValidator, QAction, QPixmap, QIcon, QColor
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout,
            QPushButton, QLabel, QLineEdit, QComboBox, QCheckBox,
            QListWidget, QMessageBox, QGroupBox, QFrame, QFileDialog,
            QAbstractItemView, QSizePolicy, QListWidgetItem,
            QTabWidget, QSpinBox, QMenu,
        )
        from PyQt6.QtCore import Qt, QUrl, QTimer
        from PyQt6.QtGui import QDesktopServices, QDoubleValidator, QAction, QPixmap, QIcon, QColor
    except ImportError:
        raise ImportError(
            "Neither PySide6 nor PyQt6 found. Install with: pip install PySide6"
        )

import pyqtgraph as pg

from pyirena.core.data_manipulation import (
    DataManipulation, ManipResult,
    ScaleConfig, TrimConfig, RebinConfig, SubtractConfig, DivideConfig,
)
from pyirena.state.state_manager import StateManager
from pyirena.gui.sas_plot import (
    make_sas_plot, make_cursors, get_cursor_q_range,
    set_robust_y_range, _SafeInfiniteLine, SASPlotStyle,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_FILE_TYPES = ["HDF5 Nexus", "HDF5 Generic", "Text (.dat/.txt)"]
_FILE_TYPE_EXTS = {
    "HDF5 Nexus":       ('.h5', '.hdf5', '.hdf'),
    "HDF5 Generic":     ('.h5', '.hdf5', '.hdf'),
    "Text (.dat/.txt)": ('.dat', '.txt'),
}

# Sort helpers — same as hdf5viewer/file_tree.py
def _sort_key_name(name: str) -> str:
    return name.lower()

def _sort_key_temperature(name: str) -> float:
    m = re.search(r'_(-?\d+(?:\.\d+)?)C(?=_|\.|$)', name, re.IGNORECASE)
    return float(m.group(1)) if m else float('inf')

def _sort_key_time(name: str) -> float:
    m = re.search(r'_(\d+(?:\.\d+)?)min(?=_|\.|$)', name, re.IGNORECASE)
    return float(m.group(1)) if m else float('inf')

def _sort_key_order(name: str) -> float:
    m = re.search(r'_(\d+)(?:_merged|_scaled|_trimmed|_rebinned|_avg|_sub|_div)?(?:\.[^.]+)?$', name)
    return float(m.group(1)) if m else float('inf')

def _sort_key_pressure(name: str) -> float:
    m = re.search(r'_(\d+(?:\.\d+)?)PSI(?=_|\.|$)', name, re.IGNORECASE)
    return float(m.group(1)) if m else float('inf')

_SORT_LABELS = [
    "Filename A\u2192Z", "Filename Z\u2192A",
    "Temperature \u2191", "Temperature \u2193",
    "Time \u2191", "Time \u2193",
    "Order number \u2191", "Order number \u2193",
    "Pressure \u2191", "Pressure \u2193",
]
_SORT_KEYS: list = [
    _sort_key_name, _sort_key_name,
    _sort_key_temperature, _sort_key_temperature,
    _sort_key_time, _sort_key_time,
    _sort_key_order, _sort_key_order,
    _sort_key_pressure, _sort_key_pressure,
]

# Distinct colors for multi-dataset plots (average tab)
_DATASET_COLORS = [
    '#2980b9', '#e74c3c', '#27ae60', '#8e44ad', '#e67e22',
    '#16a085', '#c0392b', '#2c3e50', '#f39c12', '#1abc9c',
    '#d35400', '#7f8c8d', '#9b59b6', '#34495e', '#e91e63',
]

_RESULT_PEN = pg.mkPen('#27ae60', width=3)   # green result line (thick, drawn on top)

_BTN_GREEN = ("QPushButton { background: #27ae60; color: white; font-weight: bold; "
              "border-radius: 4px; padding: 6px 10px; }"
              "QPushButton:hover { background: #2ecc71; }"
              "QPushButton:disabled { background: #95a5a6; }")
_BTN_BLUE  = ("QPushButton { background: #2980b9; color: white; font-weight: bold; "
              "border-radius: 4px; padding: 6px 10px; }"
              "QPushButton:hover { background: #3498db; }"
              "QPushButton:disabled { background: #95a5a6; }")
_BTN_GREY  = ("QPushButton { background: #7f8c8d; color: white; font-weight: bold; "
              "border-radius: 4px; padding: 4px 8px; }"
              "QPushButton:hover { background: #95a5a6; }")
_RDONLY_STYLE = "background: #ecf0f1; color: #2c3e50; border: 1px solid #bdc3c7;"

# Tab indices
_TAB_SCALE = 0
_TAB_TRIM = 1
_TAB_REBIN = 2
_TAB_AVERAGE = 3
_TAB_SUBTRACT = 4
_TAB_DIVIDE = 5


def _hline() -> QFrame:
    sep = QFrame()
    sep.setFrameShape(QFrame.Shape.HLine)
    sep.setFrameShadow(QFrame.Shadow.Sunken)
    return sep


# ===========================================================================
# _ManipFileBrowser — file browser with sort
# ===========================================================================

class _ManipFileBrowser(QWidget):
    """Vertical file-browser column with sort and filter."""

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.current_folder: Optional[str] = None
        self._all_files: List[str] = []
        self.folder_changed_callback: Optional[Callable] = None
        self._context_menu_builder: Optional[Callable] = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        hdr = QLabel("Data Files")
        hdr.setStyleSheet("font-weight: bold; color: #2c3e50;")
        layout.addWidget(hdr)

        # Folder selector
        fld_row = QHBoxLayout()
        self.folder_btn = QPushButton("Select Folder\u2026")
        self.folder_btn.setMinimumHeight(28)
        self.folder_btn.clicked.connect(self._select_folder)
        fld_row.addWidget(self.folder_btn)
        layout.addLayout(fld_row)

        self.folder_label = QLabel("(no folder selected)")
        self.folder_label.setStyleSheet("color: #7f8c8d; font-size: 11px;")
        self.folder_label.setWordWrap(True)
        layout.addWidget(self.folder_label)

        # File type
        type_row = QHBoxLayout()
        type_row.addWidget(QLabel("Type:"))
        self.type_combo = QComboBox()
        self.type_combo.addItems(_FILE_TYPES)
        self.type_combo.currentTextChanged.connect(self._refresh_file_list)
        type_row.addWidget(self.type_combo, stretch=1)
        layout.addLayout(type_row)

        # Sort
        sort_row = QHBoxLayout()
        sort_row.addWidget(QLabel("Sort:"))
        self.sort_combo = QComboBox()
        self.sort_combo.addItems(_SORT_LABELS)
        self.sort_combo.setCurrentIndex(6)  # Order number up
        self.sort_combo.currentIndexChanged.connect(self._refresh_file_list)
        sort_row.addWidget(self.sort_combo, stretch=1)
        layout.addLayout(sort_row)

        # Filter
        filt_row = QHBoxLayout()
        filt_row.addWidget(QLabel("Filter:"))
        self.filter_edit = QLineEdit()
        self.filter_edit.setPlaceholderText("text filter\u2026")
        self.filter_edit.textChanged.connect(self._apply_filter)
        filt_row.addWidget(self.filter_edit, stretch=1)
        layout.addLayout(filt_row)

        # File count
        self.count_label = QLabel("")
        self.count_label.setStyleSheet("color: #7f8c8d; font-size: 10px;")
        layout.addWidget(self.count_label)

        # File list
        self.file_list = QListWidget()
        self.file_list.setMinimumWidth(180)
        self.file_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.file_list.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.file_list.customContextMenuRequested.connect(self._show_context_menu)
        layout.addWidget(self.file_list, stretch=1)

        self.setMinimumWidth(200)
        self.setMaximumWidth(260)

    # -- Public helpers ---------------------------------------------------

    def set_folder(self, folder: str) -> None:
        self.current_folder = folder
        self.folder_label.setText(os.path.basename(folder))
        self.folder_label.setToolTip(folder)
        self._refresh_file_list()

    def get_file_type(self) -> str:
        return self.type_combo.currentText()

    def get_selected_filenames(self) -> List[str]:
        selected = {id(item): item for item in self.file_list.selectedItems()}
        return [
            self.file_list.item(i).text()
            for i in range(self.file_list.count())
            if id(self.file_list.item(i)) in selected
        ]

    def get_all_visible_filenames(self) -> List[str]:
        return [self.file_list.item(i).text() for i in range(self.file_list.count())]

    def set_context_menu_builder(self, builder: Optional[Callable]) -> None:
        """Set callback ``builder(menu, item_text, item_row)`` for right-click."""
        self._context_menu_builder = builder

    # -- Private helpers --------------------------------------------------

    def _select_folder(self) -> None:
        start = self.current_folder or str(Path.home())
        folder = QFileDialog.getExistingDirectory(self, "Select Data Folder", start)
        if folder and folder != self.current_folder:
            self.set_folder(folder)
            if self.folder_changed_callback is not None:
                self.folder_changed_callback(folder)

    def _refresh_file_list(self) -> None:
        if not self.current_folder or not os.path.isdir(self.current_folder):
            return
        exts = _FILE_TYPE_EXTS[self.type_combo.currentText()]
        try:
            files = [
                f for f in os.listdir(self.current_folder)
                if os.path.isfile(os.path.join(self.current_folder, f))
                and Path(f).suffix.lower() in exts
            ]
        except PermissionError:
            files = []

        # Sort
        idx = self.sort_combo.currentIndex()
        key_fn = _SORT_KEYS[idx]
        reverse = (idx % 2 == 1)
        files.sort(key=key_fn, reverse=reverse)

        self._all_files = files
        self._apply_filter()

    def _apply_filter(self) -> None:
        text = self.filter_edit.text().lower()
        self.file_list.clear()
        for f in self._all_files:
            if text in f.lower():
                self.file_list.addItem(f)
        n = self.file_list.count()
        self.count_label.setText(f"{n} file(s)")

    def _show_context_menu(self, pos) -> None:
        item = self.file_list.itemAt(pos)
        if item is None or self._context_menu_builder is None:
            return
        menu = QMenu(self)
        row = self.file_list.row(item)
        self._context_menu_builder(menu, item.text(), row)
        if menu.actions():
            menu.exec(self.file_list.mapToGlobal(pos))


# ===========================================================================
# DataManipulationGraphWindow — pyqtgraph plot widget
# ===========================================================================

def _color_icon(color_str: str, size: int = 12) -> QIcon:
    """Create a small square icon filled with the given color."""
    pix = QPixmap(size, size)
    pix.fill(QColor(color_str))
    return QIcon(pix)


class DataManipulationGraphWindow(QWidget):
    """Single-panel I(Q) log-log plot for the Data Manipulation tool."""

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._show_errorbars: bool = True

        # Plot items: list of (scatter_item, color_str, label_str)
        self._data_items: list = []
        self._result_item = None
        self._cursor_a: Optional[_SafeInfiniteLine] = None
        self._cursor_b: Optional[_SafeInfiniteLine] = None

        # Average-mode dataset management (set by parent panel)
        self._avg_datasets: List[str] = []       # ordered filenames
        self._avg_remove_callback = None         # fn(filename) called on remove
        self._avg_remove_after_callback = None   # fn(filename) called on remove-all-after

        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._gl = pg.GraphicsLayoutWidget()
        self._gl.setBackground('w')
        layout.addWidget(self._gl)

        self._plot = make_sas_plot(
            self._gl, row=0, col=0,
            x_label='Q  (\u00c5\u207b\u00b9)',
            y_label='I  (cm\u207b\u00b9)',
            log_x=True, log_y=True,
            parent_widget=self,
            jpeg_default_name='data_manipulation',
        )
        self._eb_action = self._plot.getViewBox().menu.addAction("Hide Error Bars")
        self._eb_action.triggered.connect(self._toggle_error_bars)

        # Intercept the ViewBox context menu to inject Average dataset actions
        vb = self._plot.getViewBox()
        vb.menu.aboutToShow.connect(self._inject_avg_menu)

    def _inject_avg_menu(self) -> None:
        """Add 'Remove dataset' / 'Remove all after' submenus when in Average mode."""
        menu = self._plot.getViewBox().menu

        # Clean up previously injected actions
        for act in getattr(self, '_injected_actions', []):
            menu.removeAction(act)
        self._injected_actions = []

        if not self._avg_datasets:
            return

        sep = menu.addSeparator()
        self._injected_actions.append(sep)

        # "Remove dataset" submenu
        remove_menu = QMenu("Remove dataset", menu)
        for i, fname in enumerate(self._avg_datasets):
            color = _DATASET_COLORS[i % len(_DATASET_COLORS)]
            act = remove_menu.addAction(_color_icon(color), fname)
            act.triggered.connect(
                lambda _checked=False, fn=fname: self._on_avg_remove(fn)
            )
        remove_action = menu.addMenu(remove_menu)
        self._injected_actions.append(remove_action)

        # "Remove all after" submenu
        remove_after_menu = QMenu("Remove all after\u2026", menu)
        for i, fname in enumerate(self._avg_datasets):
            color = _DATASET_COLORS[i % len(_DATASET_COLORS)]
            act = remove_after_menu.addAction(_color_icon(color), fname)
            act.triggered.connect(
                lambda _checked=False, fn=fname: self._on_avg_remove_after(fn)
            )
        remove_after_action = menu.addMenu(remove_after_menu)
        self._injected_actions.append(remove_after_action)

    def _on_avg_remove(self, filename: str) -> None:
        if self._avg_remove_callback:
            self._avg_remove_callback(filename)

    def _on_avg_remove_after(self, filename: str) -> None:
        if self._avg_remove_after_callback:
            self._avg_remove_after_callback(filename)

    def set_avg_datasets(self, filenames: List[str],
                         remove_cb=None, remove_after_cb=None) -> None:
        """Set the list of average-mode datasets for graph right-click menu."""
        self._avg_datasets = list(filenames)
        self._avg_remove_callback = remove_cb
        self._avg_remove_after_callback = remove_after_cb

    def clear_avg_datasets(self) -> None:
        self._avg_datasets = []
        self._avg_remove_callback = None
        self._avg_remove_after_callback = None

    # -- Public API -------------------------------------------------------

    def clear_all(self) -> None:
        """Remove all data and result items (preserves cursors)."""
        for scatter, _color, _label in self._data_items:
            if scatter is not None:
                self._plot.removeItem(scatter)
        self._data_items.clear()
        if self._result_item is not None:
            self._plot.removeItem(self._result_item)
            self._result_item = None

    def plot_data(
        self, q: np.ndarray, I: np.ndarray, dI: Optional[np.ndarray] = None,
        color: str = '#2980b9', label: str = 'Data',
    ) -> None:
        """Add a dataset scatter to the plot with custom color."""
        q = np.asarray(q, dtype=float)
        I = np.asarray(I, dtype=float)
        mask = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
        q_, I_ = q[mask], I[mask]

        brush = pg.mkBrush(color)
        scatter = self._plot.plot(
            q_, I_, pen=None,
            symbol='o', symbolSize=SASPlotStyle.DATA_SIZE,
            symbolBrush=brush, symbolPen=pg.mkPen(None),
            name=label,
        )
        self._data_items.append((scatter, color, label))

    def plot_result(self, q: np.ndarray, I: np.ndarray) -> None:
        """Overlay the manipulation result as a thick green line on top."""
        if self._result_item is not None:
            self._plot.removeItem(self._result_item)
        valid = (q > 0) & (I > 0) & np.isfinite(q) & np.isfinite(I)
        if valid.sum() < 2:
            self._result_item = None
            return
        self._result_item = self._plot.plot(
            q[valid], I[valid], pen=_RESULT_PEN, name='Result',
        )
        self._result_item.setZValue(10)  # draw on top of data scatter (Z=0)

    def set_y_range_from_data(self, *I_arrays: np.ndarray) -> None:
        """Set robust Y range from one or more intensity arrays."""
        combined = np.concatenate([arr for arr in I_arrays if arr is not None and len(arr) > 0])
        set_robust_y_range(self._plot, combined)

    def ensure_cursors(self, q_min: float, q_max: float,
                       on_moved=None) -> None:
        """Create or reposition cursors at the given Q values.

        Parameters
        ----------
        on_moved : callable or None
            If provided, connected to ``sigPositionChanged`` on both cursors
            (only on first creation to avoid duplicate connections).
        """
        if self._cursor_a is None or self._cursor_b is None:
            self._cursor_a, self._cursor_b = make_cursors(self._plot, q_min, q_max)
            if on_moved is not None:
                self._cursor_a.sigPositionChanged.connect(on_moved)
                self._cursor_b.sigPositionChanged.connect(on_moved)
        else:
            self._cursor_a.setValue(np.log10(q_min))
            self._cursor_b.setValue(np.log10(q_max))

    def remove_cursors(self) -> None:
        if self._cursor_a is not None:
            self._plot.removeItem(self._cursor_a)
            self._cursor_a = None
        if self._cursor_b is not None:
            self._plot.removeItem(self._cursor_b)
            self._cursor_b = None

    def get_cursor_q_range(self):
        return get_cursor_q_range(self._cursor_a, self._cursor_b)

    # -- Private ----------------------------------------------------------

    def _toggle_error_bars(self) -> None:
        self._show_errorbars = not self._show_errorbars
        self._eb_action.setText(
            "Hide Error Bars" if self._show_errorbars else "Show Error Bars"
        )


# ===========================================================================
# _ScalableLineEdit — QLineEdit with mouse-wheel adjustment
# ===========================================================================

class _ScalableLineEdit(QLineEdit):
    """QLineEdit that increments/decrements its float value on mouse wheel."""

    def __init__(self, text: str = "", step: float = 0.001, parent=None):
        super().__init__(text, parent)
        self._step = step
        self.setValidator(QDoubleValidator(-1e10, 1e10, 6))

    def wheelEvent(self, event):
        try:
            val = float(self.text())
        except (ValueError, TypeError):
            return
        delta = event.angleDelta().y()
        if delta > 0:
            val += self._step
        elif delta < 0:
            val -= self._step
        self.setText(f"{val:.6g}")
        self.editingFinished.emit()
        event.accept()


# ===========================================================================
# DataManipulationPanel — main window
# ===========================================================================

class DataManipulationPanel(QWidget):
    """Standalone Data Manipulation window.

    Provides a file browser, six operation tabs, a pyqtgraph plot, and
    save/export controls.
    """

    def __init__(
        self,
        state_manager: Optional[StateManager] = None,
        parent: Optional[QWidget] = None,
    ):
        super().__init__(parent)
        self.setWindowTitle("pyIrena \u2014 Data Manipulation")
        self.resize(1350, 720)

        self._sm = state_manager or StateManager()
        self._engine = DataManipulation()

        # Loaded data cache: filename -> data dict
        self._loaded: dict = {}
        self._last_result: Optional[ManipResult] = None
        self._buffer_file: Optional[str] = None      # for Subtract tab
        self._denominator_file: Optional[str] = None  # for Divide tab
        self._out_folder: Optional[str] = None

        # Throttle timer for auto-recalculate (500 ms)
        self._auto_timer = QTimer(self)
        self._auto_timer.setSingleShot(True)
        self._auto_timer.setInterval(500)
        self._auto_timer.timeout.connect(self._auto_update)

        self._build_ui()
        self._connect_auto_signals()
        self.load_state()

    # ================================================================== #
    #  UI construction                                                     #
    # ================================================================== #

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(6, 6, 6, 6)
        main_layout.setSpacing(6)

        # -- Top content row -------------------------------------------------
        top_row = QHBoxLayout()
        top_row.setSpacing(6)

        # File browser (left)
        self._fb = _ManipFileBrowser()
        self._fb.file_list.itemSelectionChanged.connect(self._on_selection_changed)
        self._fb.file_list.itemDoubleClicked.connect(self._on_file_double_clicked)
        self._fb.folder_changed_callback = self._on_folder_changed
        self._fb.set_context_menu_builder(self._build_file_context_menu)
        top_row.addWidget(self._fb)

        # Center: title + tabs + graph
        center_col = QVBoxLayout()
        center_col.setSpacing(4)

        # Title row with Help button
        title_row = QHBoxLayout()
        title_row.setContentsMargins(0, 0, 0, 0)
        title_row.setSpacing(6)
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
                "https://github.com/jilavsky/pyirena/blob/main/docs/data_manipulation_gui.md"
            ))
        )
        title_row.addWidget(_help_btn)
        center_col.addLayout(title_row)

        self._tabs = QTabWidget()
        self._tabs.setMaximumHeight(200)
        self._build_scale_tab()
        self._build_trim_tab()
        self._build_rebin_tab()
        self._build_average_tab()
        self._build_subtract_tab()
        self._build_divide_tab()
        self._tabs.currentChanged.connect(self._on_tab_changed)
        center_col.addWidget(self._tabs)

        self._graph = DataManipulationGraphWindow()
        center_col.addWidget(self._graph, stretch=1)
        top_row.addLayout(center_col, stretch=1)

        # Save panel (right)
        top_row.addWidget(self._build_save_panel())

        main_layout.addLayout(top_row, stretch=1)

        # -- Status bar ------------------------------------------------------
        self._status = QLabel("Ready \u2014 Select a folder and data files to begin.")
        self._status.setStyleSheet("color: #555; font-size: 11px;")
        main_layout.addWidget(self._status)

    # ------------------------------------------------------------------ #
    #  Scale tab                                                           #
    # ------------------------------------------------------------------ #

    def _build_scale_tab(self) -> None:
        tab = QWidget()
        layout = QHBoxLayout(tab)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(12)

        # Scale I
        layout.addWidget(QLabel("Scale I:"))
        self._scale_I_edit = _ScalableLineEdit("1.0", step=0.01)
        self._scale_I_edit.setMaximumWidth(100)
        layout.addWidget(self._scale_I_edit)

        # Background
        layout.addWidget(QLabel("Background:"))
        self._scale_bg_edit = _ScalableLineEdit("0.0", step=0.001)
        self._scale_bg_edit.setMaximumWidth(100)
        layout.addWidget(self._scale_bg_edit)

        # Scale uncertainty
        layout.addWidget(QLabel("Scale uncert.:"))
        self._scale_unc_edit = QLineEdit("")
        self._scale_unc_edit.setPlaceholderText("= Scale I")
        self._scale_unc_edit.setValidator(QDoubleValidator(-1e10, 1e10, 6))
        self._scale_unc_edit.setMaximumWidth(100)
        layout.addWidget(self._scale_unc_edit)

        layout.addStretch()
        self._tabs.addTab(tab, "Scale")

    # ------------------------------------------------------------------ #
    #  Trim tab                                                            #
    # ------------------------------------------------------------------ #

    def _build_trim_tab(self) -> None:
        tab = QWidget()
        layout = QHBoxLayout(tab)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(12)

        layout.addWidget(QLabel("Q min:"))
        self._trim_qmin = QLineEdit("")
        self._trim_qmin.setReadOnly(True)
        self._trim_qmin.setStyleSheet(_RDONLY_STYLE)
        self._trim_qmin.setMaximumWidth(100)
        layout.addWidget(self._trim_qmin)

        layout.addWidget(QLabel("Q max:"))
        self._trim_qmax = QLineEdit("")
        self._trim_qmax.setReadOnly(True)
        self._trim_qmax.setStyleSheet(_RDONLY_STYLE)
        self._trim_qmax.setMaximumWidth(100)
        layout.addWidget(self._trim_qmax)

        layout.addWidget(QLabel("\u00c5\u207b\u00b9   (drag cursors on plot)"))

        layout.addStretch()
        self._tabs.addTab(tab, "Trim")

    # ------------------------------------------------------------------ #
    #  Rebin tab                                                           #
    # ------------------------------------------------------------------ #

    def _build_rebin_tab(self) -> None:
        tab = QWidget()
        layout = QHBoxLayout(tab)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(12)

        layout.addWidget(QLabel("Grid:"))
        self._rebin_mode = QComboBox()
        self._rebin_mode.addItems(["Log-spaced", "Linear", "From reference file"])
        self._rebin_mode.currentIndexChanged.connect(self._on_rebin_mode_changed)
        layout.addWidget(self._rebin_mode)

        layout.addWidget(QLabel("Points:"))
        self._rebin_npts = QSpinBox()
        self._rebin_npts.setRange(10, 10000)
        self._rebin_npts.setValue(200)
        self._rebin_npts.setMaximumWidth(80)
        layout.addWidget(self._rebin_npts)

        layout.addWidget(QLabel("Q min:"))
        self._rebin_qmin = QLineEdit("")
        self._rebin_qmin.setPlaceholderText("auto")
        self._rebin_qmin.setValidator(QDoubleValidator(0, 100, 8))
        self._rebin_qmin.setMaximumWidth(80)
        layout.addWidget(self._rebin_qmin)

        layout.addWidget(QLabel("Q max:"))
        self._rebin_qmax = QLineEdit("")
        self._rebin_qmax.setPlaceholderText("auto")
        self._rebin_qmax.setValidator(QDoubleValidator(0, 100, 8))
        self._rebin_qmax.setMaximumWidth(80)
        layout.addWidget(self._rebin_qmax)

        self._rebin_ref_btn = QPushButton("Load Reference\u2026")
        self._rebin_ref_btn.setMaximumHeight(26)
        self._rebin_ref_btn.setEnabled(False)
        self._rebin_ref_btn.setToolTip(
            "Select 'From reference file' in the Grid combo first.\n"
            "Then click here to load any data file — its Q values\n"
            "will be used as the new Q grid for rebinning."
        )
        self._rebin_ref_btn.clicked.connect(self._load_rebin_reference)
        layout.addWidget(self._rebin_ref_btn)
        self._rebin_ref_q: Optional[np.ndarray] = None
        self._rebin_ref_label = QLabel("")
        self._rebin_ref_label.setStyleSheet("color: #555; font-size: 10px;")
        layout.addWidget(self._rebin_ref_label)

        layout.addStretch()
        self._tabs.addTab(tab, "Rebin")

    # ------------------------------------------------------------------ #
    #  Average tab                                                         #
    # ------------------------------------------------------------------ #

    def _build_average_tab(self) -> None:
        tab = QWidget()
        layout = QHBoxLayout(tab)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(12)

        info = QLabel(
            "Select multiple files to average.\n"
            "Right-click on graph to remove bad datasets."
        )
        info.setStyleSheet("color: #555;")
        layout.addWidget(info)

        self._avg_count_label = QLabel("0 files selected")
        self._avg_count_label.setStyleSheet("font-weight: bold; color: #2c3e50;")
        layout.addWidget(self._avg_count_label)

        layout.addStretch()
        self._tabs.addTab(tab, "Average")

    # ------------------------------------------------------------------ #
    #  Subtract tab                                                        #
    # ------------------------------------------------------------------ #

    def _build_subtract_tab(self) -> None:
        tab = QWidget()
        outer = QVBoxLayout(tab)
        outer.setContentsMargins(8, 4, 8, 4)
        outer.setSpacing(4)

        # Row 1: dataset names
        row1 = QHBoxLayout()
        row1.setSpacing(16)
        self._sub_sample_label = QLabel("Sample: (none)")
        self._sub_sample_label.setStyleSheet("color: #2980b9;")
        row1.addWidget(self._sub_sample_label)
        self._sub_buffer_label = QLabel("Buffer: (none)")
        self._sub_buffer_label.setStyleSheet("color: #e74c3c; font-weight: bold;")
        row1.addWidget(self._sub_buffer_label)
        row1.addStretch()
        outer.addLayout(row1)

        # Row 2: controls
        row2 = QHBoxLayout()
        row2.setSpacing(8)
        row2.addWidget(QLabel("Buffer scale:"))
        self._sub_scale_edit = _ScalableLineEdit("1.0", step=0.001)
        self._sub_scale_edit.setMaximumWidth(100)
        row2.addWidget(self._sub_scale_edit)

        self._sub_auto_chk = QCheckBox("Auto-scale")
        self._sub_auto_chk.stateChanged.connect(self._on_sub_auto_changed)
        row2.addWidget(self._sub_auto_chk)

        self._sub_qmin_label = QLabel("Q min:")
        self._sub_qmin_label.setVisible(False)
        row2.addWidget(self._sub_qmin_label)
        self._sub_qmin = QLineEdit("")
        self._sub_qmin.setReadOnly(True)
        self._sub_qmin.setStyleSheet(_RDONLY_STYLE)
        self._sub_qmin.setMaximumWidth(80)
        self._sub_qmin.setVisible(False)
        row2.addWidget(self._sub_qmin)

        self._sub_qmax_label = QLabel("Q max:")
        self._sub_qmax_label.setVisible(False)
        row2.addWidget(self._sub_qmax_label)
        self._sub_qmax = QLineEdit("")
        self._sub_qmax.setReadOnly(True)
        self._sub_qmax.setStyleSheet(_RDONLY_STYLE)
        self._sub_qmax.setMaximumWidth(80)
        self._sub_qmax.setVisible(False)
        row2.addWidget(self._sub_qmax)

        row2.addStretch()
        outer.addLayout(row2)

        self._tabs.addTab(tab, "Subtract")

    # ------------------------------------------------------------------ #
    #  Divide tab                                                          #
    # ------------------------------------------------------------------ #

    def _build_divide_tab(self) -> None:
        tab = QWidget()
        outer = QVBoxLayout(tab)
        outer.setContentsMargins(8, 4, 8, 4)
        outer.setSpacing(4)

        # Row 1: dataset names
        row1 = QHBoxLayout()
        row1.setSpacing(16)
        self._div_num_label = QLabel("Numerator: (none)")
        self._div_num_label.setStyleSheet("color: #2980b9;")
        row1.addWidget(self._div_num_label)
        self._div_den_label = QLabel("Denominator: (none)")
        self._div_den_label.setStyleSheet("color: #e74c3c; font-weight: bold;")
        row1.addWidget(self._div_den_label)
        row1.addStretch()
        outer.addLayout(row1)

        # Row 2: controls
        row2 = QHBoxLayout()
        row2.setSpacing(8)
        row2.addWidget(QLabel("Denom. scale:"))
        self._div_scale_edit = QLineEdit("1.0")
        self._div_scale_edit.setValidator(QDoubleValidator(-1e10, 1e10, 6))
        self._div_scale_edit.setMaximumWidth(100)
        row2.addWidget(self._div_scale_edit)

        row2.addWidget(QLabel("Denom. bg:"))
        self._div_bg_edit = QLineEdit("0.0")
        self._div_bg_edit.setValidator(QDoubleValidator(-1e10, 1e10, 6))
        self._div_bg_edit.setMaximumWidth(100)
        row2.addWidget(self._div_bg_edit)

        row2.addStretch()
        outer.addLayout(row2)

        self._tabs.addTab(tab, "Divide")

    # ------------------------------------------------------------------ #
    #  Save panel (right column)                                           #
    # ------------------------------------------------------------------ #

    def _build_save_panel(self) -> QWidget:
        panel = QWidget()
        panel.setFixedWidth(180)
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)

        layout.addWidget(QLabel("Output folder:"))
        self._out_folder_label = QLabel("(not set)")
        self._out_folder_label.setWordWrap(True)
        self._out_folder_label.setStyleSheet("color: #2c3e50; font-size: 11px;")
        layout.addWidget(self._out_folder_label)

        create_btn = QPushButton("Create Default Folder")
        create_btn.setMinimumHeight(28)
        create_btn.setStyleSheet(_BTN_GREY)
        create_btn.setToolTip(
            "Create a new folder named '<data folder>_manip' next to the data folder."
        )
        create_btn.clicked.connect(self._create_default_output_folder)
        layout.addWidget(create_btn)

        select_btn = QPushButton("Select Existing Folder\u2026")
        select_btn.setMinimumHeight(28)
        select_btn.setStyleSheet(_BTN_GREY)
        select_btn.clicked.connect(self._select_output_folder)
        layout.addWidget(select_btn)

        layout.addWidget(_hline())

        self._save_btn = QPushButton("Apply && Save")
        self._save_btn.setMinimumHeight(34)
        self._save_btn.setStyleSheet(_BTN_GREEN)
        self._save_btn.setEnabled(False)
        self._save_btn.clicked.connect(self._apply_and_save)
        layout.addWidget(self._save_btn)

        layout.addWidget(_hline())

        self._batch_btn = QPushButton("Batch All Selected")
        self._batch_btn.setMinimumHeight(34)
        self._batch_btn.setStyleSheet(_BTN_BLUE)
        self._batch_btn.setEnabled(False)
        self._batch_btn.clicked.connect(self._batch_run)
        layout.addWidget(self._batch_btn)

        layout.addStretch()
        return panel

    # ================================================================== #
    #  Public helpers                                                      #
    # ================================================================== #

    def set_folder(self, path: str) -> None:
        """Pre-populate the file browser folder (called from Data Selector)."""
        self._fb.set_folder(path)

    # ================================================================== #
    #  Data loading                                                        #
    # ================================================================== #

    def _load_file(self, filepath: str) -> Optional[dict]:
        """Load a SAS data file. Returns dict with Q, Intensity, Error, dQ, etc."""
        from pyirena.io.hdf5 import readGenericNXcanSAS, readSimpleHDF5, readTextFile
        fp = Path(filepath)
        file_type = self._fb.get_file_type()
        try:
            if file_type == "Text (.dat/.txt)":
                data = readTextFile(str(fp.parent), fp.name)
                data['is_nxcansas'] = False
            elif file_type == "HDF5 Generic":
                data = readSimpleHDF5(str(fp.parent), fp.name)
                data['is_nxcansas'] = False
            else:  # HDF5 Nexus
                data = readGenericNXcanSAS(str(fp.parent), fp.name)
                data['is_nxcansas'] = True
            data['filepath'] = filepath
            return data
        except Exception as exc:
            QMessageBox.warning(self, "Load Error", f"Could not read {fp.name}:\n{exc}")
            return None

    def _get_or_load(self, filename: str) -> Optional[dict]:
        """Get cached data or load from disk."""
        folder = self._fb.current_folder
        if not folder:
            return None
        filepath = os.path.join(folder, filename)
        if filepath not in self._loaded:
            data = self._load_file(filepath)
            if data is None:
                return None
            self._loaded[filepath] = data
        return self._loaded[filepath]

    # ================================================================== #
    #  Auto-recalculate                                                    #
    # ================================================================== #

    def _connect_auto_signals(self) -> None:
        """Wire parameter-change signals to throttled auto-update."""
        # Scale tab
        self._scale_I_edit.editingFinished.connect(self._schedule_auto_update)
        self._scale_bg_edit.editingFinished.connect(self._schedule_auto_update)
        self._scale_unc_edit.editingFinished.connect(self._schedule_auto_update)

        # Rebin tab
        self._rebin_mode.currentIndexChanged.connect(self._schedule_auto_update)
        self._rebin_npts.valueChanged.connect(self._schedule_auto_update)
        self._rebin_qmin.editingFinished.connect(self._schedule_auto_update)
        self._rebin_qmax.editingFinished.connect(self._schedule_auto_update)

        # Subtract tab — editingFinished also fires from wheelEvent
        self._sub_scale_edit.editingFinished.connect(self._schedule_auto_update)
        self._sub_auto_chk.stateChanged.connect(self._schedule_auto_update)

        # Divide tab
        self._div_scale_edit.editingFinished.connect(self._schedule_auto_update)
        self._div_bg_edit.editingFinished.connect(self._schedule_auto_update)

    def _schedule_auto_update(self) -> None:
        """Restart the throttle timer (500 ms single-shot)."""
        self._auto_timer.start()

    def _auto_update(self) -> None:
        """Fire the appropriate preview for the active tab."""
        selected = self._fb.get_selected_filenames()
        if not selected:
            return

        tab = self._tabs.currentIndex()
        if tab == _TAB_SCALE:
            self._preview_scale()
        elif tab == _TAB_TRIM:
            self._preview_trim()
        elif tab == _TAB_REBIN:
            self._preview_rebin()
        elif tab == _TAB_AVERAGE:
            self._preview_average(auto_triggered=True)
        elif tab == _TAB_SUBTRACT:
            self._preview_subtract()
        elif tab == _TAB_DIVIDE:
            self._preview_divide()

    # ================================================================== #
    #  Event handlers                                                      #
    # ================================================================== #

    def _on_folder_changed(self, folder: str) -> None:
        self._loaded.clear()
        self._buffer_file = None
        self._denominator_file = None
        self._last_result = None
        self._graph.clear_all()
        self._graph.remove_cursors()
        self._update_sub_labels()
        self._update_div_labels()
        self._check_enable_buttons()

    def _on_selection_changed(self) -> None:
        self._check_enable_buttons()
        tab = self._tabs.currentIndex()
        if tab == _TAB_AVERAGE:
            n = len(self._fb.get_selected_filenames())
            self._avg_count_label.setText(f"{n} file(s) selected")
            self._schedule_auto_update()
        elif tab == _TAB_SUBTRACT:
            self._update_sub_labels()
        elif tab == _TAB_DIVIDE:
            self._update_div_labels()

    def _on_file_double_clicked(self, item: QListWidgetItem) -> None:
        """Double-click: plot the file (tab-aware) and auto-calculate."""
        tab = self._tabs.currentIndex()

        if tab == _TAB_SUBTRACT:
            self._update_sub_labels()
            self._plot_subtract_inputs()
            self._status.setText(f"Sample: {item.text()}")
            if self._sub_auto_chk.isChecked():
                data = self._get_or_load(item.text())
                if data is not None:
                    self._ensure_subtract_cursors(data)
            self._schedule_auto_update()
            return

        if tab == _TAB_DIVIDE:
            self._update_div_labels()
            self._plot_divide_inputs()
            self._status.setText(f"Numerator: {item.text()}")
            self._schedule_auto_update()
            return

        # Default: plot the file and trigger auto-calculate
        data = self._get_or_load(item.text())
        if data is None:
            return
        self._graph.clear_all()
        self._graph.plot_data(data['Q'], data['Intensity'], data.get('Error'))
        self._graph.set_y_range_from_data(data['Intensity'])
        self._status.setText(f"Plotted: {item.text()}")

        if tab == _TAB_TRIM:
            self._ensure_trim_cursors(data)

        self._schedule_auto_update()

    def _on_tab_changed(self, index: int) -> None:
        """Handle tab change — clear stale graph and manage cursors."""
        # Clear graph and result when switching tabs
        self._graph.clear_all()
        self._graph.remove_cursors()
        self._last_result = None

        # Clear average dataset tracking when leaving Average tab
        if index != _TAB_AVERAGE:
            self._graph.clear_avg_datasets()

        self._update_cursor_display()
        self._check_enable_buttons()

    def _on_rebin_mode_changed(self, index: int) -> None:
        is_ref = (index == 2)  # "From reference file"
        self._rebin_ref_btn.setEnabled(is_ref)
        self._rebin_npts.setEnabled(not is_ref)
        self._rebin_qmin.setEnabled(not is_ref)
        self._rebin_qmax.setEnabled(not is_ref)

    def _on_sub_auto_changed(self, state) -> None:
        auto = self._sub_auto_chk.isChecked()
        for w in (self._sub_qmin_label, self._sub_qmin,
                  self._sub_qmax_label, self._sub_qmax):
            w.setVisible(auto)
        if auto:
            # Show cursors if data is loaded
            selected = self._fb.get_selected_filenames()
            if selected:
                data = self._get_or_load(selected[0])
                if data is not None:
                    self._ensure_subtract_cursors(data)
        else:
            if self._tabs.currentIndex() == _TAB_SUBTRACT:
                self._graph.remove_cursors()

    # ================================================================== #
    #  Context menu                                                        #
    # ================================================================== #

    def _build_file_context_menu(self, menu: QMenu, filename: str, row: int) -> None:
        tab = self._tabs.currentIndex()

        # NOTE: QAction.triggered emits a bool `checked` argument.  Lambdas
        # that use a default-value parameter (e.g. ``lambda fn=val:``) would
        # receive that bool as the first positional arg, shadowing the default.
        # Use an explicit ``_checked`` parameter to absorb it.

        # Average tab context menu
        if tab == _TAB_AVERAGE:
            act_remove = menu.addAction("Remove from selection")
            act_remove.triggered.connect(
                lambda _checked=False, _r=row: self._fb.file_list.item(_r).setSelected(False)
            )
            act_remove_after = menu.addAction("Remove all after this")
            act_remove_after.triggered.connect(
                lambda _checked=False, r=row: self._deselect_after(r)
            )

        # Subtract tab context menu
        elif tab == _TAB_SUBTRACT:
            act_buffer = menu.addAction("Set as buffer")
            act_buffer.triggered.connect(
                lambda _checked=False, fn=filename: self._set_buffer(fn)
            )

        # Divide tab context menu
        elif tab == _TAB_DIVIDE:
            act_denom = menu.addAction("Set as denominator")
            act_denom.triggered.connect(
                lambda _checked=False, fn=filename: self._set_denominator(fn)
            )

    def _deselect_after(self, row: int) -> None:
        for i in range(row + 1, self._fb.file_list.count()):
            self._fb.file_list.item(i).setSelected(False)
        self._on_selection_changed()

    def _set_buffer(self, filename: str) -> None:
        self._buffer_file = filename
        self._update_sub_labels()
        self._plot_subtract_inputs()

    def _set_denominator(self, filename: str) -> None:
        self._denominator_file = filename
        self._update_div_labels()
        self._plot_divide_inputs()

    def _update_sub_labels(self) -> None:
        buf = self._buffer_file or "(none \u2014 right-click to set)"
        self._sub_buffer_label.setText(f"Buffer: {buf}")
        # Sample = first selected that is not buffer
        selected = self._fb.get_selected_filenames()
        samples = [f for f in selected if f != self._buffer_file]
        if samples:
            self._sub_sample_label.setText(f"Sample: {samples[0]}")
        else:
            self._sub_sample_label.setText("Sample: (none)")

    def _update_div_labels(self) -> None:
        den = self._denominator_file or "(none \u2014 right-click to set)"
        self._div_den_label.setText(f"Denominator: {den}")
        selected = self._fb.get_selected_filenames()
        nums = [f for f in selected if f != self._denominator_file]
        if nums:
            self._div_num_label.setText(f"Numerator: {nums[0]}")
        else:
            self._div_num_label.setText("Numerator: (none)")

    # ================================================================== #
    #  Auto-plot helpers (Subtract / Divide)                               #
    # ================================================================== #

    def _plot_subtract_inputs(self) -> None:
        """Plot sample + buffer datasets when both are identified."""
        selected = self._fb.get_selected_filenames()
        samples = [f for f in selected if f != self._buffer_file]
        self._graph.clear_all()
        all_I = []
        if samples:
            data = self._get_or_load(samples[0])
            if data is not None:
                self._graph.plot_data(data['Q'], data['Intensity'],
                                      data.get('Error'), color='#2980b9', label='Sample')
                all_I.append(data['Intensity'])
        if self._buffer_file:
            data = self._get_or_load(self._buffer_file)
            if data is not None:
                self._graph.plot_data(data['Q'], data['Intensity'],
                                      data.get('Error'), color='#e74c3c', label='Buffer')
                all_I.append(data['Intensity'])
        if all_I:
            self._graph.set_y_range_from_data(*all_I)

    def _plot_divide_inputs(self) -> None:
        """Plot numerator + denominator datasets when both are identified."""
        selected = self._fb.get_selected_filenames()
        nums = [f for f in selected if f != self._denominator_file]
        self._graph.clear_all()
        all_I = []
        if nums:
            data = self._get_or_load(nums[0])
            if data is not None:
                self._graph.plot_data(data['Q'], data['Intensity'],
                                      data.get('Error'), color='#2980b9', label='Numerator')
                all_I.append(data['Intensity'])
        if self._denominator_file:
            data = self._get_or_load(self._denominator_file)
            if data is not None:
                self._graph.plot_data(data['Q'], data['Intensity'],
                                      data.get('Error'), color='#e74c3c', label='Denominator')
                all_I.append(data['Intensity'])
        if all_I:
            self._graph.set_y_range_from_data(*all_I)

    # ================================================================== #
    #  Cursor helpers                                                      #
    # ================================================================== #

    def _on_cursor_moved(self, _line=None) -> None:
        """Called when a cursor is dragged — update display and recalculate."""
        self._update_cursor_display()
        self._schedule_auto_update()

    def _ensure_trim_cursors(self, data: dict) -> None:
        q = data['Q']
        q_lo = float(q.min()) * 1.05
        q_hi = float(q.max()) * 0.95
        self._graph.ensure_cursors(q_lo, q_hi, on_moved=self._on_cursor_moved)
        self._update_cursor_display()

    def _ensure_subtract_cursors(self, data: dict) -> None:
        q = data['Q']
        # Place cursors at high-Q end (last 20%)
        q_lo = float(np.percentile(q, 80))
        q_hi = float(q.max()) * 0.95
        self._graph.ensure_cursors(q_lo, q_hi, on_moved=self._on_cursor_moved)
        self._update_cursor_display()

    def _update_cursor_display(self) -> None:
        q_min, q_max = self._graph.get_cursor_q_range()
        tab = self._tabs.currentIndex()
        if tab == _TAB_TRIM:
            self._trim_qmin.setText(f"{q_min:.6g}" if q_min is not None else "")
            self._trim_qmax.setText(f"{q_max:.6g}" if q_max is not None else "")
        elif tab == _TAB_SUBTRACT and self._sub_auto_chk.isChecked():
            self._sub_qmin.setText(f"{q_min:.6g}" if q_min is not None else "")
            self._sub_qmax.setText(f"{q_max:.6g}" if q_max is not None else "")

    # ================================================================== #
    #  Preview operations                                                  #
    # ================================================================== #

    def _preview_scale(self) -> None:
        selected = self._fb.get_selected_filenames()
        if not selected:
            self._status.setText("Select at least one file to preview.")
            return
        data = self._get_or_load(selected[0])
        if data is None:
            return

        try:
            s = float(self._scale_I_edit.text() or "1")
        except ValueError:
            s = 1.0
        try:
            bg = float(self._scale_bg_edit.text() or "0")
        except ValueError:
            bg = 0.0
        unc_text = self._scale_unc_edit.text().strip()
        s_unc = float(unc_text) if unc_text else None

        config = ScaleConfig(scale_I=s, background=bg, scale_uncertainty=s_unc)
        result = self._engine.scale(
            data['Q'], data['Intensity'], data.get('Error', data['Intensity'] * 0.05),
            data.get('dQ'), config,
        )
        self._show_preview(data, result)

    def _preview_trim(self) -> None:
        selected = self._fb.get_selected_filenames()
        if not selected:
            self._status.setText("Select at least one file to preview.")
            return
        data = self._get_or_load(selected[0])
        if data is None:
            return

        self._update_cursor_display()
        q_min, q_max = self._graph.get_cursor_q_range()
        if q_min is None or q_max is None:
            self._ensure_trim_cursors(data)
            q_min, q_max = self._graph.get_cursor_q_range()
        if q_min is None or q_max is None:
            self._status.setText("Could not determine cursor positions.")
            return

        config = TrimConfig(q_min=q_min, q_max=q_max)
        result = self._engine.trim(
            data['Q'], data['Intensity'], data.get('Error', data['Intensity'] * 0.05),
            data.get('dQ'), config,
        )
        self._show_preview(data, result)

    def _preview_rebin(self) -> None:
        selected = self._fb.get_selected_filenames()
        if not selected:
            self._status.setText("Select at least one file to preview.")
            return
        data = self._get_or_load(selected[0])
        if data is None:
            return

        mode_idx = self._rebin_mode.currentIndex()
        mode = ['log', 'linear', 'reference'][mode_idx]
        n_pts = self._rebin_npts.value()
        qmin_text = self._rebin_qmin.text().strip()
        qmax_text = self._rebin_qmax.text().strip()
        q_min = float(qmin_text) if qmin_text else None
        q_max = float(qmax_text) if qmax_text else None

        config = RebinConfig(
            mode=mode, n_points=n_pts, q_min=q_min, q_max=q_max,
            reference_q=self._rebin_ref_q,
        )
        try:
            result = self._engine.rebin(
                data['Q'], data['Intensity'], data.get('Error', data['Intensity'] * 0.05),
                data.get('dQ'), config,
            )
        except ValueError as e:
            self._status.setText(f"Error: {e}")
            return
        self._show_preview(data, result)

    def _preview_average(self, auto_triggered: bool = False) -> None:
        selected = self._fb.get_selected_filenames()
        if len(selected) < 1:
            if not auto_triggered:
                self._status.setText("Select at least 2 files to average.")
            return

        datasets = []
        loaded_names = []
        self._graph.clear_all()
        all_I = []
        for i, fname in enumerate(selected):
            data = self._get_or_load(fname)
            if data is None:
                continue
            color = _DATASET_COLORS[i % len(_DATASET_COLORS)]
            self._graph.plot_data(data['Q'], data['Intensity'], data.get('Error'), color=color, label=fname)
            q = data['Q']
            I = data['Intensity']
            dI = data.get('Error', I * 0.05)
            dQ = data.get('dQ')
            datasets.append((q, I, dI, dQ))
            loaded_names.append(fname)
            all_I.append(I)

        # Register datasets for graph right-click menu
        self._graph.set_avg_datasets(
            loaded_names,
            remove_cb=self._avg_remove_dataset,
            remove_after_cb=self._avg_remove_after_dataset,
        )

        if len(datasets) < 2:
            if not auto_triggered:
                self._status.setText("Select at least 2 files to average.")
            self._last_result = None
            if all_I:
                self._graph.set_y_range_from_data(*all_I)
            return

        result = self._engine.average(datasets, reference_index=0)
        self._last_result = result
        self._graph.plot_result(result.q, result.I)
        if all_I:
            self._graph.set_y_range_from_data(*all_I)
        self._save_btn.setEnabled(self._out_folder is not None)
        self._status.setText(f"Average of {len(datasets)} datasets previewed.")

    def _avg_remove_dataset(self, filename: str) -> None:
        """Remove a single dataset from the Average selection and re-plot."""
        # Deselect in file list
        for i in range(self._fb.file_list.count()):
            item = self._fb.file_list.item(i)
            if item.text() == filename and item.isSelected():
                item.setSelected(False)
                break
        self._on_selection_changed()
        self._preview_average(auto_triggered=True)

    def _avg_remove_after_dataset(self, filename: str) -> None:
        """Remove all datasets after the given one and re-plot."""
        selected = self._fb.get_selected_filenames()
        try:
            idx = selected.index(filename)
        except ValueError:
            return
        # Deselect all files after this one in the selected order
        to_deselect = set(selected[idx + 1:])
        for i in range(self._fb.file_list.count()):
            item = self._fb.file_list.item(i)
            if item.text() in to_deselect:
                item.setSelected(False)
        self._on_selection_changed()
        self._preview_average(auto_triggered=True)

    def _preview_subtract(self) -> None:
        selected = self._fb.get_selected_filenames()
        samples = [f for f in selected if f != self._buffer_file]
        if not samples:
            self._status.setText("Select a sample file.")
            return
        if not self._buffer_file:
            self._status.setText("Right-click a file to set it as buffer.")
            return

        sample_data = self._get_or_load(samples[0])
        buffer_data = self._get_or_load(self._buffer_file)
        if sample_data is None or buffer_data is None:
            return

        try:
            scale = float(self._sub_scale_edit.text() or "1")
        except ValueError:
            scale = 1.0

        auto = self._sub_auto_chk.isChecked()
        auto_qmin = auto_qmax = None
        if auto:
            self._update_cursor_display()
            q_min, q_max = self._graph.get_cursor_q_range()
            auto_qmin = q_min
            auto_qmax = q_max

        config = SubtractConfig(
            buffer_scale=scale, auto_scale=auto,
            auto_q_min=auto_qmin, auto_q_max=auto_qmax,
        )
        result = self._engine.subtract(
            sample_data['Q'], sample_data['Intensity'],
            sample_data.get('Error', sample_data['Intensity'] * 0.05),
            sample_data.get('dQ'),
            buffer_data['Q'], buffer_data['Intensity'],
            buffer_data.get('Error', buffer_data['Intensity'] * 0.05),
            config,
        )
        # Update scale display if auto-scaled
        if auto:
            self._sub_scale_edit.setText(f"{result.metadata.get('buffer_scale', scale):.6g}")

        self._graph.clear_all()
        self._graph.plot_data(sample_data['Q'], sample_data['Intensity'],
                              sample_data.get('Error'), color='#2980b9', label='Sample')
        self._graph.plot_data(buffer_data['Q'], buffer_data['Intensity'],
                              buffer_data.get('Error'), color='#e74c3c', label='Buffer')
        self._last_result = result
        self._graph.plot_result(result.q, result.I)
        self._graph.set_y_range_from_data(sample_data['Intensity'], buffer_data['Intensity'])
        self._save_btn.setEnabled(self._out_folder is not None)
        self._status.setText(f"Subtraction previewed (scale={result.metadata.get('buffer_scale', scale):.4g}).")

    def _preview_divide(self) -> None:
        selected = self._fb.get_selected_filenames()
        nums = [f for f in selected if f != self._denominator_file]
        if not nums:
            self._status.setText("Select a numerator file.")
            return
        if not self._denominator_file:
            self._status.setText("Right-click a file to set it as denominator.")
            return

        num_data = self._get_or_load(nums[0])
        den_data = self._get_or_load(self._denominator_file)
        if num_data is None or den_data is None:
            return

        try:
            scale = float(self._div_scale_edit.text() or "1")
        except ValueError:
            scale = 1.0
        try:
            bg = float(self._div_bg_edit.text() or "0")
        except ValueError:
            bg = 0.0

        config = DivideConfig(denominator_scale=scale, denominator_background=bg)
        result = self._engine.divide(
            num_data['Q'], num_data['Intensity'],
            num_data.get('Error', num_data['Intensity'] * 0.05),
            num_data.get('dQ'),
            den_data['Q'], den_data['Intensity'],
            den_data.get('Error', den_data['Intensity'] * 0.05),
            config,
        )

        self._graph.clear_all()
        self._graph.plot_data(num_data['Q'], num_data['Intensity'],
                              num_data.get('Error'), color='#2980b9', label='Numerator')
        self._graph.plot_data(den_data['Q'], den_data['Intensity'],
                              den_data.get('Error'), color='#e74c3c', label='Denominator')
        self._last_result = result
        self._graph.plot_result(result.q, result.I)
        self._graph.set_y_range_from_data(num_data['Intensity'], den_data['Intensity'])
        self._save_btn.setEnabled(self._out_folder is not None)
        self._status.setText("Division previewed.")

    # ================================================================== #
    #  Common preview helper                                               #
    # ================================================================== #

    def _show_preview(self, data: dict, result: ManipResult) -> None:
        """Show input data + result on the plot."""
        self._graph.clear_all()
        self._graph.plot_data(data['Q'], data['Intensity'], data.get('Error'))
        self._last_result = result
        self._graph.plot_result(result.q, result.I)
        self._graph.set_y_range_from_data(data['Intensity'])
        self._save_btn.setEnabled(self._out_folder is not None)
        n_in = len(data['Q'])
        n_out = len(result.q)
        self._status.setText(
            f"{result.operation} preview: {n_in} \u2192 {n_out} points."
        )

    # ================================================================== #
    #  Reference file loading (Rebin tab)                                  #
    # ================================================================== #

    def _load_rebin_reference(self) -> None:
        start = self._fb.current_folder or str(Path.home())
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Reference Q Grid", start,
            "HDF5/Text (*.h5 *.hdf5 *.hdf *.dat *.txt);;All files (*)",
        )
        if not path:
            return
        data = self._load_reference_file(path)
        if data is None:
            return
        self._rebin_ref_q = data['Q']
        self._rebin_ref_label.setText(f"({len(data['Q'])} pts from {Path(path).name})")
        self._status.setText(f"Reference Q grid loaded: {len(data['Q'])} points.")

    def _load_reference_file(self, filepath: str) -> Optional[dict]:
        """Load a reference data file, auto-detecting type by extension."""
        from pyirena.io.hdf5 import readGenericNXcanSAS, readSimpleHDF5, readTextFile
        fp = Path(filepath)
        ext = fp.suffix.lower()
        try:
            if ext in ('.dat', '.txt'):
                data = readTextFile(str(fp.parent), fp.name)
            elif ext in ('.h5', '.hdf5', '.hdf'):
                data = readGenericNXcanSAS(str(fp.parent), fp.name)
            else:
                data = readGenericNXcanSAS(str(fp.parent), fp.name)
            return data
        except Exception as exc:
            QMessageBox.warning(self, "Load Error",
                                f"Could not read reference file {fp.name}:\n{exc}")
            return None

    # ================================================================== #
    #  Save / Batch                                                        #
    # ================================================================== #

    def _apply_and_save(self) -> None:
        """Save the last previewed result to disk."""
        if self._last_result is None:
            self._status.setText("Nothing to save \u2014 preview an operation first.")
            return
        if not self._out_folder:
            self._status.setText("Set an output folder first.")
            return

        tab = self._tabs.currentIndex()
        selected = self._fb.get_selected_filenames()

        if tab == _TAB_AVERAGE:
            # For average, use first selected file as source reference
            if not selected:
                return
            source_file = selected[0]
        elif tab == _TAB_SUBTRACT:
            samples = [f for f in selected if f != self._buffer_file]
            source_file = samples[0] if samples else None
        elif tab == _TAB_DIVIDE:
            nums = [f for f in selected if f != self._denominator_file]
            source_file = nums[0] if nums else None
        else:
            source_file = selected[0] if selected else None

        if not source_file:
            self._status.setText("No source file identified.")
            return

        self._save_result(source_file, self._last_result)

    def _save_result(self, source_filename: str, result: ManipResult) -> Optional[Path]:
        """Save a ManipResult to HDF5."""
        from pyirena.io.nxcansas_data_manipulation import save_manipulated_data

        folder = self._fb.current_folder
        if not folder:
            return None
        source_path = Path(os.path.join(folder, source_filename))

        # Determine if source is NXcanSAS
        data = self._get_or_load(source_filename)
        is_nxcansas = data.get('is_nxcansas', False) if data else False

        try:
            out_path = save_manipulated_data(
                output_folder=Path(self._out_folder),
                source_path=source_path,
                source_is_nxcansas=is_nxcansas,
                q=result.q,
                I=result.I,
                dI=result.dI,
                dQ=result.dQ,
                operation=result.operation,
                provenance=result.metadata,
            )
            self._status.setText(f"Saved: {out_path.name}")
            return out_path
        except Exception as exc:
            QMessageBox.warning(self, "Save Error", f"Error saving:\n{exc}")
            self._status.setText(f"Save error: {exc}")
            return None

    def _batch_run(self) -> None:
        """Apply current operation to all selected files and save."""
        if not self._out_folder:
            self._status.setText("Set an output folder first.")
            return

        tab = self._tabs.currentIndex()
        selected = self._fb.get_selected_filenames()
        if not selected:
            self._status.setText("No files selected.")
            return

        saved = 0
        errors = 0

        if tab == _TAB_AVERAGE:
            # Average is a single operation on all selected files
            self._preview_average()
            if self._last_result is not None:
                out = self._save_result(selected[0], self._last_result)
                saved = 1 if out else 0
        elif tab in (_TAB_SUBTRACT, _TAB_DIVIDE):
            # Process each sample file (excluding buffer/denominator)
            if tab == _TAB_SUBTRACT:
                exclude = self._buffer_file
                if not exclude:
                    self._status.setText("Set a buffer first.")
                    return
            else:
                exclude = self._denominator_file
                if not exclude:
                    self._status.setText("Set a denominator first.")
                    return

            samples = [f for f in selected if f != exclude]
            for fname in samples:
                data = self._get_or_load(fname)
                if data is None:
                    errors += 1
                    continue
                result = self._run_operation(tab, fname)
                if result is not None:
                    out = self._save_result(fname, result)
                    if out:
                        saved += 1
                    else:
                        errors += 1
                else:
                    errors += 1
        else:
            # Scale, Trim, Rebin: process each selected file independently
            for fname in selected:
                result = self._run_operation(tab, fname)
                if result is not None:
                    out = self._save_result(fname, result)
                    if out:
                        saved += 1
                    else:
                        errors += 1
                else:
                    errors += 1

        msg = f"Batch complete: {saved} saved"
        if errors:
            msg += f", {errors} error(s)"
        self._status.setText(msg)

    def _run_operation(self, tab: int, filename: str) -> Optional[ManipResult]:
        """Run the current tab's operation on a single file."""
        data = self._get_or_load(filename)
        if data is None:
            return None

        q = data['Q']
        I = data['Intensity']
        dI = data.get('Error', I * 0.05)
        dQ = data.get('dQ')

        try:
            if tab == _TAB_SCALE:
                s = float(self._scale_I_edit.text() or "1")
                bg = float(self._scale_bg_edit.text() or "0")
                unc_text = self._scale_unc_edit.text().strip()
                s_unc = float(unc_text) if unc_text else None
                return self._engine.scale(q, I, dI, dQ,
                                          ScaleConfig(s, bg, s_unc))
            elif tab == _TAB_TRIM:
                self._update_cursor_display()
                q_min, q_max = self._graph.get_cursor_q_range()
                if q_min is None or q_max is None:
                    return None
                return self._engine.trim(q, I, dI, dQ,
                                         TrimConfig(q_min, q_max))
            elif tab == _TAB_REBIN:
                mode_idx = self._rebin_mode.currentIndex()
                mode = ['log', 'linear', 'reference'][mode_idx]
                n_pts = self._rebin_npts.value()
                qmin_text = self._rebin_qmin.text().strip()
                qmax_text = self._rebin_qmax.text().strip()
                return self._engine.rebin(q, I, dI, dQ, RebinConfig(
                    mode=mode, n_points=n_pts,
                    q_min=float(qmin_text) if qmin_text else None,
                    q_max=float(qmax_text) if qmax_text else None,
                    reference_q=self._rebin_ref_q,
                ))
            elif tab == _TAB_SUBTRACT:
                buf_data = self._get_or_load(self._buffer_file)
                if buf_data is None:
                    return None
                scale = float(self._sub_scale_edit.text() or "1")
                auto = self._sub_auto_chk.isChecked()
                auto_qmin = auto_qmax = None
                if auto:
                    q_min, q_max = self._graph.get_cursor_q_range()
                    auto_qmin, auto_qmax = q_min, q_max
                return self._engine.subtract(
                    q, I, dI, dQ,
                    buf_data['Q'], buf_data['Intensity'],
                    buf_data.get('Error', buf_data['Intensity'] * 0.05),
                    SubtractConfig(scale, auto, auto_qmin, auto_qmax),
                )
            elif tab == _TAB_DIVIDE:
                den_data = self._get_or_load(self._denominator_file)
                if den_data is None:
                    return None
                s = float(self._div_scale_edit.text() or "1")
                bg = float(self._div_bg_edit.text() or "0")
                return self._engine.divide(
                    q, I, dI, dQ,
                    den_data['Q'], den_data['Intensity'],
                    den_data.get('Error', den_data['Intensity'] * 0.05),
                    DivideConfig(s, bg),
                )
        except Exception as exc:
            self._status.setText(f"Error processing {filename}: {exc}")
            return None
        return None

    # ================================================================== #
    #  Output folder helpers                                               #
    # ================================================================== #

    def _create_default_output_folder(self) -> None:
        folder = self._fb.current_folder
        if not folder:
            self._status.setText("Select a data folder first.")
            return
        out = folder.rstrip('/\\') + '_manip'
        os.makedirs(out, exist_ok=True)
        self._out_folder = out
        self._out_folder_label.setText(os.path.basename(out))
        self._out_folder_label.setToolTip(out)
        self._check_enable_buttons()
        self._status.setText(f"Output folder: {out}")

    def _select_output_folder(self) -> None:
        start = self._fb.current_folder or str(Path.home())
        folder = QFileDialog.getExistingDirectory(self, "Select Output Folder", start)
        if folder:
            self._out_folder = folder
            self._out_folder_label.setText(os.path.basename(folder))
            self._out_folder_label.setToolTip(folder)
            self._check_enable_buttons()

    # ================================================================== #
    #  Button enable logic                                                 #
    # ================================================================== #

    def _check_enable_buttons(self) -> None:
        has_output = self._out_folder is not None
        has_result = self._last_result is not None
        has_selection = len(self._fb.get_selected_filenames()) > 0
        self._save_btn.setEnabled(has_output and has_result)
        self._batch_btn.setEnabled(has_output and has_selection)

    # ================================================================== #
    #  State persistence                                                   #
    # ================================================================== #

    def load_state(self) -> None:
        s = self._sm.get('data_manipulation') or {}
        folder = s.get('folder')
        if folder and os.path.isdir(folder):
            ft = s.get('file_type', 'HDF5 Nexus')
            self._fb.type_combo.setCurrentText(ft)
            sort_idx = s.get('sort_index', 6)
            self._fb.sort_combo.setCurrentIndex(int(sort_idx))
            self._fb.filter_edit.setText(s.get('filter', ''))
            self._fb.set_folder(folder)

        out = s.get('output_folder')
        if out and os.path.isdir(out):
            self._out_folder = out
            self._out_folder_label.setText(os.path.basename(out))
            self._out_folder_label.setToolTip(out)

        tab = s.get('active_tab', 0)
        self._tabs.setCurrentIndex(int(tab))

        # Scale tab
        self._scale_I_edit.setText(str(s.get('scale_I', 1.0)))
        self._scale_bg_edit.setText(str(s.get('scale_background', 0.0)))
        su = s.get('scale_uncertainty')
        self._scale_unc_edit.setText(str(su) if su is not None else "")

        # Rebin tab
        mode = s.get('rebin_mode', 'log')
        idx = {'log': 0, 'linear': 1, 'reference': 2}.get(mode, 0)
        self._rebin_mode.setCurrentIndex(idx)
        self._rebin_npts.setValue(int(s.get('rebin_n_points', 200)))
        rqmin = s.get('rebin_q_min')
        rqmax = s.get('rebin_q_max')
        self._rebin_qmin.setText(str(rqmin) if rqmin is not None else "")
        self._rebin_qmax.setText(str(rqmax) if rqmax is not None else "")

        # Subtract tab
        self._sub_scale_edit.setText(str(s.get('subtract_buffer_scale', 1.0)))
        self._sub_auto_chk.setChecked(bool(s.get('subtract_auto_scale', False)))

        # Divide tab
        self._div_scale_edit.setText(str(s.get('divide_denominator_scale', 1.0)))
        self._div_bg_edit.setText(str(s.get('divide_denominator_background', 0.0)))

    def save_state(self) -> None:
        su_text = self._scale_unc_edit.text().strip()
        rqmin_text = self._rebin_qmin.text().strip()
        rqmax_text = self._rebin_qmax.text().strip()

        self._sm.update('data_manipulation', {
            'folder': self._fb.current_folder,
            'file_type': self._fb.get_file_type(),
            'sort_index': self._fb.sort_combo.currentIndex(),
            'filter': self._fb.filter_edit.text(),
            'output_folder': self._out_folder,
            'active_tab': self._tabs.currentIndex(),
            # Scale
            'scale_I': float(self._scale_I_edit.text() or 1),
            'scale_background': float(self._scale_bg_edit.text() or 0),
            'scale_uncertainty': float(su_text) if su_text else None,
            # Rebin
            'rebin_mode': ['log', 'linear', 'reference'][self._rebin_mode.currentIndex()],
            'rebin_n_points': self._rebin_npts.value(),
            'rebin_q_min': float(rqmin_text) if rqmin_text else None,
            'rebin_q_max': float(rqmax_text) if rqmax_text else None,
            # Subtract
            'subtract_buffer_scale': float(self._sub_scale_edit.text() or 1),
            'subtract_auto_scale': self._sub_auto_chk.isChecked(),
            # Divide
            'divide_denominator_scale': float(self._div_scale_edit.text() or 1),
            'divide_denominator_background': float(self._div_bg_edit.text() or 0),
        })
        self._sm.save()

    def closeEvent(self, event) -> None:
        self.save_state()
        super().closeEvent(event)


# ===========================================================================
# Standalone entry point
# ===========================================================================

def main() -> None:
    """Entry point: ``pyirena-datamanip``."""
    import argparse

    parser = argparse.ArgumentParser(
        prog='pyirena-datamanip',
        description='Data Manipulation tool for SAS datasets.',
    )
    parser.add_argument('--folder', metavar='DIR',
                        help='Pre-populate the data folder.')
    args = parser.parse_args()

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
        app.setStyle('Fusion')

    window = DataManipulationPanel()
    if args.folder:
        window.set_folder(args.folder)
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
