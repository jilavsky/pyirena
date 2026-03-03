"""
data_merge_panel.py — GUI panel and CLI entry point for the Data Merge tool.

Merges two SAS datasets (lower-Q DS1 + higher-Q DS2) by optimising a scale
factor, optional Q-shift, and constant background in a user-selected overlap
region.

Entry points
------------
* ``DataMergePanel`` — standalone QWidget; launched from the Data Selector hub.
* ``main()`` — ``pyirena-datamerge`` CLI command.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Optional, List, Tuple

import numpy as np

try:
    from PySide6.QtWidgets import (
        QApplication, QWidget, QVBoxLayout, QHBoxLayout,
        QPushButton, QLabel, QLineEdit, QComboBox, QCheckBox,
        QListWidget, QMessageBox, QGroupBox, QFrame, QFileDialog,
        QAbstractItemView, QSizePolicy, QListWidgetItem,
    )
    from PySide6.QtCore import Qt, QUrl
    from PySide6.QtGui import QDesktopServices
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout,
            QPushButton, QLabel, QLineEdit, QComboBox, QCheckBox,
            QListWidget, QMessageBox, QGroupBox, QFrame, QFileDialog,
            QAbstractItemView, QSizePolicy, QListWidgetItem,
        )
        from PyQt6.QtCore import Qt, QUrl
        from PyQt6.QtGui import QDesktopServices
    except ImportError:
        raise ImportError(
            "Neither PySide6 nor PyQt6 found. Install with: pip install PySide6"
        )

import pyqtgraph as pg

from pyirena.core.data_merge import DataMerge, MergeConfig, MergeResult
from pyirena.state.state_manager import StateManager
from pyirena.gui.sas_plot import (
    make_sas_plot, plot_iq_data, plot_iq_model,
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

# Colours specific to this tool
_DS1_BRUSH  = pg.mkBrush('#2980b9')   # blue — absolute-scale data
_DS2_BRUSH  = pg.mkBrush('#e74c3c')   # red  — data to be scaled
_DS1_PEN    = pg.mkPen('#2980b9', width=1)
_DS2_PEN    = pg.mkPen('#e74c3c', width=1)
_MERGED_PEN = pg.mkPen('#27ae60', width=2)   # dark green

_BTN_GREEN  = ("QPushButton { background: #27ae60; color: white; font-weight: bold; "
               "border-radius: 4px; padding: 6px 10px; }"
               "QPushButton:hover { background: #2ecc71; }"
               "QPushButton:disabled { background: #95a5a6; }")
_BTN_BLUE   = ("QPushButton { background: #2980b9; color: white; font-weight: bold; "
               "border-radius: 4px; padding: 6px 10px; }"
               "QPushButton:hover { background: #3498db; }"
               "QPushButton:disabled { background: #95a5a6; }")
_BTN_GREY   = ("QPushButton { background: #7f8c8d; color: white; font-weight: bold; "
               "border-radius: 4px; padding: 4px 8px; }"
               "QPushButton:hover { background: #95a5a6; }")
_RDONLY_STYLE = "background: #ecf0f1; color: #2c3e50; border: 1px solid #bdc3c7;"


# ===========================================================================
# _DatasetSelectorWidget — reusable single-dataset file browser
# ===========================================================================

class _DatasetSelectorWidget(QWidget):
    """A vertical file-browser column for one dataset.

    Attributes
    ----------
    dataset_number : int
        1 or 2 (displayed in the header).
    """

    def __init__(self, dataset_number: int, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.dataset_number = dataset_number
        self.current_folder: Optional[str] = None
        self._all_files: List[str] = []
        # Optional callback invoked when the user clicks "Select Folder" and a new
        # folder is chosen.  Set by DataMergePanel after construction.
        self.folder_changed_callback = None

        self._build_ui()

    # ------------------------------------------------------------------ #
    #  UI construction                                                     #
    # ------------------------------------------------------------------ #

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        dn = self.dataset_number
        label_text = (
            f"Dataset {dn} — lower Q (absolute scale)"
            if dn == 1 else
            f"Dataset {dn} — higher Q"
        )
        hdr = QLabel(label_text)
        hdr.setStyleSheet("font-weight: bold; color: #2c3e50;")
        layout.addWidget(hdr)

        # Folder selector
        fld_row = QHBoxLayout()
        self.folder_btn = QPushButton("Select Folder…")
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

        # Filter
        filt_row = QHBoxLayout()
        filt_row.addWidget(QLabel("Filter:"))
        self.filter_edit = QLineEdit()
        self.filter_edit.setPlaceholderText("text filter…")
        self.filter_edit.textChanged.connect(self._apply_filter)
        filt_row.addWidget(self.filter_edit, stretch=1)
        layout.addLayout(filt_row)

        # File list — ExtendedSelection allows Ctrl/Shift multi-select.
        # Selected pairs (first DS1 with first DS2, etc.) are used for batch run.
        self.file_list = QListWidget()
        self.file_list.setMinimumWidth(180)
        self.file_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        layout.addWidget(self.file_list, stretch=1)

        self.setMinimumWidth(200)
        self.setMaximumWidth(260)

    # ------------------------------------------------------------------ #
    #  Public helpers                                                      #
    # ------------------------------------------------------------------ #

    def set_folder(self, folder: str) -> None:
        """Programmatically set the folder and refresh the file list."""
        self.current_folder = folder
        self.folder_label.setText(os.path.basename(folder))
        self.folder_label.setToolTip(folder)
        self._refresh_file_list()

    def get_file_type(self) -> str:
        return self.type_combo.currentText()

    def get_selected_filename(self) -> Optional[str]:
        items = self.file_list.selectedItems()
        return items[0].text() if items else None

    def get_selected_filenames(self) -> List[str]:
        """Return all selected filenames in list order."""
        selected = {id(item): item for item in self.file_list.selectedItems()}
        return [
            self.file_list.item(i).text()
            for i in range(self.file_list.count())
            if id(self.file_list.item(i)) in selected
        ]

    def set_displayed_files(self, names: List[str]) -> None:
        """Replace the list contents (used by match mode)."""
        self.file_list.clear()
        for n in names:
            self.file_list.addItem(n)

    def get_all_files(self) -> List[str]:
        """Return the unfiltered file list for the current folder/type."""
        return list(self._all_files)

    # ------------------------------------------------------------------ #
    #  Private helpers                                                     #
    # ------------------------------------------------------------------ #

    def _select_folder(self) -> None:
        start = self.current_folder or str(Path.home())
        folder = QFileDialog.getExistingDirectory(
            self, f"Select Dataset {self.dataset_number} Folder", start
        )
        if folder and folder != self.current_folder:
            self.set_folder(folder)
            if self.folder_changed_callback is not None:
                self.folder_changed_callback(folder)

    def _refresh_file_list(self) -> None:
        if not self.current_folder or not os.path.isdir(self.current_folder):
            return
        exts = _FILE_TYPE_EXTS[self.type_combo.currentText()]
        try:
            files = sorted(
                f for f in os.listdir(self.current_folder)
                if os.path.isfile(os.path.join(self.current_folder, f))
                and Path(f).suffix.lower() in exts
            )
        except PermissionError:
            files = []
        self._all_files = files
        self._apply_filter()

    def _apply_filter(self) -> None:
        text = self.filter_edit.text().lower()
        self.file_list.clear()
        for f in self._all_files:
            if text in f.lower():
                self.file_list.addItem(f)


# ===========================================================================
# DataMergeGraphWindow — pyqtgraph plot widget
# ===========================================================================

class DataMergeGraphWindow(QWidget):
    """Single-panel I(Q) plot for the Data Merge tool.

    Supports SAXS (log-log) and WAXS/diffraction (lin-lin) modes.
    The two movable cursors A and B define the overlap/optimisation region.
    """

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._log_mode: bool = True   # True = log-log (SAXS), False = lin-lin (WAXS)
        self._show_errorbars: bool = True

        # Plot items (None = not yet plotted)
        self._ds1_scatter = None
        self._ds1_err     = None
        self._ds2_scatter = None
        self._ds2_err     = None
        self._merged_item = None
        self._cursor_a: Optional[_SafeInfiniteLine] = None
        self._cursor_b: Optional[_SafeInfiniteLine] = None

        # Stash raw data for re-plotting after mode switch
        self._q1: Optional[np.ndarray] = None
        self._I1: Optional[np.ndarray] = None
        self._dI1: Optional[np.ndarray] = None
        self._q2: Optional[np.ndarray] = None
        self._I2: Optional[np.ndarray] = None
        self._dI2: Optional[np.ndarray] = None

        self._build_ui()

    # ------------------------------------------------------------------ #
    #  UI construction                                                     #
    # ------------------------------------------------------------------ #

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._gl = pg.GraphicsLayoutWidget()
        self._gl.setBackground('w')
        layout.addWidget(self._gl)

        self._plot = make_sas_plot(
            self._gl, row=0, col=0,
            x_label='Q  (Å⁻¹)',
            y_label='I  (cm⁻¹)',
            log_x=self._log_mode,
            log_y=self._log_mode,
            parent_widget=self,
            jpeg_default_name='data_merge',
        )
        self._eb_action = self._plot.getViewBox().menu.addAction("Hide Error Bars")
        self._eb_action.triggered.connect(self._toggle_error_bars)

    # ------------------------------------------------------------------ #
    #  Public API                                                          #
    # ------------------------------------------------------------------ #

    def set_mode(self, saxs: bool) -> None:
        """Switch between SAXS (log-log) and WAXS (lin-lin) display mode.

        Rebuilds the plot and re-plots any loaded data.
        """
        if saxs == self._log_mode:
            return
        # Save cursor Q positions (physical units) before destroying them
        q_cur_min, q_cur_max = self.get_overlap_range()

        self._log_mode = saxs
        # Rebuild the PlotItem with the new mode
        self._gl.clear()
        self._ds1_scatter = self._ds1_err = None
        self._ds2_scatter = self._ds2_err = None
        self._merged_item = None
        self._cursor_a = self._cursor_b = None

        self._plot = make_sas_plot(
            self._gl, row=0, col=0,
            x_label='Q  (Å⁻¹)',
            y_label='I  (cm⁻¹)',
            log_x=self._log_mode,
            log_y=self._log_mode,
            parent_widget=self,
            jpeg_default_name='data_merge',
        )
        eb_label = "Hide Error Bars" if self._show_errorbars else "Show Error Bars"
        self._eb_action = self._plot.getViewBox().menu.addAction(eb_label)
        self._eb_action.triggered.connect(self._toggle_error_bars)

        # Re-plot loaded datasets
        if self._q1 is not None:
            self.plot_ds1(self._q1, self._I1, self._dI1)
        if self._q2 is not None:
            self.plot_ds2(self._q2, self._I2, self._dI2)

        # Restore cursors at their previous physical Q positions
        if q_cur_min is not None and q_cur_max is not None:
            self.init_cursors(q_cur_min, q_cur_max)

    def plot_ds1(
        self, q: np.ndarray, I: np.ndarray, dI: Optional[np.ndarray] = None
    ) -> None:
        """Plot DS1 scatter (blue). Replaces any previous DS1 plot."""
        self._q1, self._I1, self._dI1 = q, I, dI
        self._remove_ds1()

        # Use plot.plot(symbol='o') — PlotDataItem — rather than ScatterPlotItem.
        # PlotDataItem explicitly applies log10 to data when setLogMode is active,
        # giving a correct bounding rect for auto-range in log-log mode.
        # Standalone ScatterPlotItem does not pre-transform data; pyqtgraph's
        # auto-range then uses linear Q values as ViewBox (log10-space) coordinates,
        # pushing the actual scatter points entirely out of the visible window.
        mask = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
        self._ds1_scatter = self._plot.plot(
            q[mask], I[mask],
            pen=None,
            symbol='o', symbolSize=SASPlotStyle.DATA_SIZE,
            symbolBrush=_DS1_BRUSH, symbolPen=pg.mkPen(None),
        )

        if dI is not None:
            err = self._make_error_bars(q, I, dI)
            err.setVisible(self._show_errorbars)
            self._ds1_err = err

        self._update_y_range()

    def plot_ds2(
        self, q: np.ndarray, I: np.ndarray, dI: Optional[np.ndarray] = None
    ) -> None:
        """Plot DS2 scatter (red). Replaces any previous DS2 plot."""
        self._q2, self._I2, self._dI2 = q, I, dI
        self._remove_ds2()

        mask = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
        self._ds2_scatter = self._plot.plot(
            q[mask], I[mask],
            pen=None,
            symbol='o', symbolSize=SASPlotStyle.DATA_SIZE,
            symbolBrush=_DS2_BRUSH, symbolPen=pg.mkPen(None),
        )

        if dI is not None:
            err = self._make_error_bars(q, I, dI)
            err.setVisible(self._show_errorbars)
            self._ds2_err = err

        self._update_y_range()

    def plot_merged(
        self, q: np.ndarray, I: np.ndarray, dI: Optional[np.ndarray] = None
    ) -> None:
        """Overlay merged result (dark green line)."""
        self._remove_merged()
        item = pg.PlotDataItem(x=q, y=I, pen=_MERGED_PEN)
        self._plot.addItem(item)
        self._merged_item = item

    def clear_merged(self) -> None:
        self._remove_merged()

    def init_cursors(self, q_min: float, q_max: float) -> None:
        """Create or reposition the two overlap cursors."""
        if self._cursor_a is None:
            self._cursor_a = self._make_cursor(q_min, SASPlotStyle.CURSOR_A_PEN, 'A')
            self._cursor_b = self._make_cursor(q_max, SASPlotStyle.CURSOR_B_PEN, 'B')
            self._plot.addItem(self._cursor_a)
            self._plot.addItem(self._cursor_b)
        else:
            self._set_cursor_pos(self._cursor_a, q_min)
            self._set_cursor_pos(self._cursor_b, q_max)

    def get_overlap_range(self) -> Tuple[Optional[float], Optional[float]]:
        """Return (q_min, q_max) in physical units from cursor positions."""
        if self._cursor_a is None or self._cursor_b is None:
            return None, None
        pos_a = self._cursor_a.getPos()[0]
        pos_b = self._cursor_b.getPos()[0]
        qa = 10 ** pos_a if self._log_mode else pos_a
        qb = 10 ** pos_b if self._log_mode else pos_b
        q_min, q_max = (qa, qb) if qa < qb else (qb, qa)
        return q_min, q_max

    def get_cursor_positions_raw(self) -> Tuple[float, float]:
        """Return raw ViewBox positions of both cursors (for state save)."""
        if self._cursor_a is None:
            return 0.0, 0.0
        return self._cursor_a.getPos()[0], self._cursor_b.getPos()[0]

    # ------------------------------------------------------------------ #
    #  Private helpers                                                     #
    # ------------------------------------------------------------------ #

    def _make_cursor(
        self, q: float, pen, label: str
    ) -> _SafeInfiniteLine:
        pos = np.log10(q) if self._log_mode else q
        c = _SafeInfiniteLine(
            pos=pos, angle=90, movable=True, pen=pen,
            label=label, labelOpts={'position': 0.95, 'color': pen.color()},
        )
        return c

    def _set_cursor_pos(self, cursor: _SafeInfiniteLine, q: float) -> None:
        pos = np.log10(q) if self._log_mode else q
        cursor.setPos(pos)

    def _make_error_bars(
        self, q: np.ndarray, I: np.ndarray, dI: np.ndarray
    ) -> pg.PlotDataItem:
        """Build NaN-separated error-bar line segments."""
        valid = (q > 0) & (I > 0) & np.isfinite(q) & np.isfinite(I) & np.isfinite(dI)
        q_v, I_v, dI_v = q[valid], I[valid], dI[valid]

        xs: list = []
        ys: list = []
        cap = SASPlotStyle.ERROR_CAP_FRAC

        for qi, Ii, dIi in zip(q_v, I_v, dI_v):
            I_lo = max(Ii - dIi, Ii * 0.001)
            I_hi = Ii + dIi
            # Vertical bar
            xs += [qi, qi, float('nan')]
            ys += [I_lo, I_hi, float('nan')]
            if self._log_mode:
                q_lo_cap = qi * (10 ** (-cap))
                q_hi_cap = qi * (10 ** cap)
            else:
                q_lo_cap = qi * (1 - cap)
                q_hi_cap = qi * (1 + cap)
            # Top cap
            xs += [q_lo_cap, q_hi_cap, float('nan')]
            ys += [I_hi, I_hi, float('nan')]
            # Bottom cap
            xs += [q_lo_cap, q_hi_cap, float('nan')]
            ys += [I_lo, I_lo, float('nan')]

        item = self._plot.plot(
            xs, ys,
            pen=SASPlotStyle.ERROR_PEN,
            connect='finite',
        )
        return item

    def _update_y_range(self) -> None:
        """Combine I from all loaded datasets and set robust Y (and X) range."""
        parts = []
        if self._I1 is not None:
            parts.append(self._I1)
        if self._I2 is not None:
            parts.append(self._I2)
        if not parts:
            return
        I_all = np.concatenate(parts)
        if self._log_mode:
            set_robust_y_range(self._plot, I_all)
        else:
            # Lin-lin auto-range
            valid = np.isfinite(I_all) & (I_all > 0)
            if valid.any():
                lo = float(np.percentile(I_all[valid], 1))
                hi = float(np.percentile(I_all[valid], 99))
                pad = (hi - lo) * 0.1
                self._plot.setYRange(lo - pad, hi + pad, padding=0)

        # Constrain x-axis so auto-range does not snap to linear Q values.
        # Without this, NaN-separated error bar segments cause pyqtgraph to
        # auto-range x to the full linear Q span (e.g. 0–20 Å⁻¹) instead of
        # the correct log10 range, hiding the low-Q data.
        q_parts = []
        if self._q1 is not None:
            q_parts.append(self._q1)
        if self._q2 is not None:
            q_parts.append(self._q2)
        if not q_parts:
            return
        q_all = np.concatenate(q_parts)
        valid_q = q_all[(q_all > 0) & np.isfinite(q_all)]
        if len(valid_q) < 2:
            return
        if self._log_mode:
            q_lo = int(np.floor(np.log10(float(valid_q.min())))) - 1
            q_hi = int(np.ceil(np.log10(float(valid_q.max())))) + 1
            # setLimits constrains pan/zoom in ViewBox (log10) coordinates.
            self._plot.getViewBox().setLimits(xMin=q_lo, xMax=q_hi)
            # Explicitly set the visible x range to the actual data Q span.
            # Pass log10(physical Q) — these are ViewBox coordinates in log mode.
            # This matches the pattern used in GraphWindow (data_selector.py).
            self._plot.setXRange(
                np.log10(float(valid_q.min())),
                np.log10(float(valid_q.max())),
                padding=0.05,
            )
        else:
            lo = float(valid_q.min())
            hi = float(valid_q.max())
            pad = (hi - lo) * 0.05
            self._plot.setXRange(lo, hi + pad, padding=0)

    def _toggle_error_bars(self) -> None:
        """Show/hide error bar items (callable from right-click menu)."""
        self._show_errorbars = not self._show_errorbars
        label = "Show Error Bars" if not self._show_errorbars else "Hide Error Bars"
        self._eb_action.setText(label)
        for item in (self._ds1_err, self._ds2_err):
            if item is not None:
                item.setVisible(self._show_errorbars)

    def _remove_ds1(self) -> None:
        for item in (self._ds1_scatter, self._ds1_err):
            if item is not None:
                self._plot.removeItem(item)
        self._ds1_scatter = self._ds1_err = None

    def _remove_ds2(self) -> None:
        for item in (self._ds2_scatter, self._ds2_err):
            if item is not None:
                self._plot.removeItem(item)
        self._ds2_scatter = self._ds2_err = None

    def _remove_merged(self) -> None:
        if self._merged_item is not None:
            self._plot.removeItem(self._merged_item)
            self._merged_item = None


# ===========================================================================
# DataMergePanel — main window
# ===========================================================================

class DataMergePanel(QWidget):
    """Standalone Data Merge window.

    Provides two file selectors (DS1 and DS2), a pyqtgraph plot with cursors,
    merge controls, and save/export controls.  Can be launched from the Data
    Selector hub or run standalone via the ``pyirena-datamerge`` CLI.
    """

    def __init__(
        self,
        state_manager: Optional[StateManager] = None,
        parent: Optional[QWidget] = None,
    ):
        super().__init__(parent)
        self.setWindowTitle("pyIrena — Data Merge")
        self.resize(1350, 720)

        self._sm = state_manager or StateManager()
        self._engine = DataMerge()

        self._data1: Optional[dict] = None
        self._data2: Optional[dict] = None
        self._last_result: Optional[MergeResult] = None
        self._last_config: Optional[MergeConfig] = None
        self._last_q_merged: Optional[np.ndarray] = None
        self._last_I_merged: Optional[np.ndarray] = None
        self._last_dI_merged: Optional[np.ndarray] = None
        self._last_dQ_merged: Optional[np.ndarray] = None

        self._build_ui()
        self.load_state()

    # ================================================================== #
    #  UI construction                                                     #
    # ================================================================== #

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(6, 6, 6, 6)
        main_layout.setSpacing(6)

        # ── Top content row ──────────────────────────────────────────────
        top_row = QHBoxLayout()
        top_row.setSpacing(6)

        # DS1 and DS2 columns
        self._ds1 = _DatasetSelectorWidget(1)
        self._ds2 = _DatasetSelectorWidget(2)
        self._ds1.file_list.itemDoubleClicked.connect(self._on_ds1_double_clicked)
        self._ds2.file_list.itemDoubleClicked.connect(self._on_ds2_double_clicked)
        # Re-check Batch button enable state whenever selection changes
        self._ds1.file_list.itemSelectionChanged.connect(self._check_enable_batch_btn)
        self._ds2.file_list.itemSelectionChanged.connect(self._check_enable_batch_btn)
        # When the user picks a new folder via the folder button: clear output folder
        # (DS1 only) and uncheck "Match files" (both datasets).
        self._ds1.folder_changed_callback = self._on_ds1_folder_changed
        self._ds2.folder_changed_callback = self._on_ds2_folder_changed

        top_row.addWidget(self._ds1)
        top_row.addWidget(self._ds2)

        # Center: controls + graph
        center_col = QVBoxLayout()
        center_col.setSpacing(4)
        center_col.addLayout(self._build_controls())
        center_col.addWidget(self._build_q_range_box())
        self._graph = DataMergeGraphWindow()
        center_col.addWidget(self._graph, stretch=1)
        top_row.addLayout(center_col, stretch=1)

        # Save panel (right column)
        top_row.addWidget(self._build_save_panel())

        main_layout.addLayout(top_row, stretch=1)

        # ── Bottom row: match checkbox + status ──────────────────────────
        bottom_row = QHBoxLayout()
        self._match_chk = QCheckBox("Match files  "
                                    "(prefix before '_'  +  last number must match)")
        self._match_chk.stateChanged.connect(self._on_match_changed)
        bottom_row.addWidget(self._match_chk)
        bottom_row.addStretch()
        main_layout.addLayout(bottom_row)

        self._status = QLabel("Ready — load a DS1 and DS2 file to begin.")
        self._status.setStyleSheet("color: #7f8c8d; padding: 2px;")
        main_layout.addWidget(self._status)

    def _build_controls(self) -> QVBoxLayout:
        vbox = QVBoxLayout()
        vbox.setSpacing(4)

        # ── Row 1: Mode / Scale / Q-shift / Background ───────────────────
        row1 = QHBoxLayout()
        row1.setSpacing(8)

        row1.addWidget(QLabel("Mode:"))
        self._mode_combo = QComboBox()
        self._mode_combo.addItems(["SAXS (log-log)", "WAXS / diffraction (lin-lin)"])
        self._mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        row1.addWidget(self._mode_combo)

        row1.addWidget(_vline())

        # ── Scale ────────────────────────────────────────────────────────
        row1.addWidget(QLabel("Scale:"))
        self._scale_ds_combo = QComboBox()
        self._scale_ds_combo.addItems(["DS1", "DS2"])
        self._scale_ds_combo.setCurrentIndex(1)   # default: scale DS2
        row1.addWidget(self._scale_ds_combo)

        self._fit_scale_chk = QCheckBox("Fit")
        self._fit_scale_chk.setChecked(True)
        self._fit_scale_chk.setToolTip("Fit scale factor during optimisation")
        self._fit_scale_chk.toggled.connect(self._on_fit_scale_toggled)
        row1.addWidget(self._fit_scale_chk)

        self._scale_result = QLineEdit("1.0000")
        self._scale_result.setReadOnly(True)
        self._scale_result.setFixedWidth(75)
        self._scale_result.setStyleSheet(_RDONLY_STYLE)
        self._scale_result.setToolTip(
            "Optimised scale factor (read-only when Fit is checked).\n"
            "Uncheck Fit to enter a fixed scale value manually."
        )
        row1.addWidget(self._scale_result)

        row1.addWidget(_vline())

        # ── Q shift ──────────────────────────────────────────────────────
        row1.addWidget(QLabel("Q shift:"))
        self._qshift_combo = QComboBox()
        self._qshift_combo.addItems(["None", "DS1", "DS2"])
        row1.addWidget(self._qshift_combo)

        self._fit_qshift_chk = QCheckBox("Fit")
        self._fit_qshift_chk.setChecked(False)
        self._fit_qshift_chk.setToolTip("Fit Q-shift during optimisation")
        self._fit_qshift_chk.toggled.connect(self._on_fit_qshift_toggled)
        row1.addWidget(self._fit_qshift_chk)

        self._qshift_result = QLineEdit("0.000000")
        self._qshift_result.setReadOnly(False)   # editable by default (Fit unchecked)
        self._qshift_result.setFixedWidth(80)
        self._qshift_result.setToolTip(
            "Optimised Q shift (Å⁻¹).\n"
            "Editable when Fit is unchecked — enter a known Q offset (e.g. 0)."
        )
        row1.addWidget(self._qshift_result)

        row1.addWidget(_vline())

        # ── Background ────────────────────────────────────────────────────
        row1.addWidget(QLabel("BG (DS1):"))
        self._bg_result = QLineEdit("0.000000")
        self._bg_result.setReadOnly(True)
        self._bg_result.setFixedWidth(80)
        self._bg_result.setStyleSheet(_RDONLY_STYLE)
        self._bg_result.setToolTip("Optimised constant background subtracted from DS1")
        row1.addWidget(self._bg_result)

        row1.addStretch()
        _help_btn = QPushButton("? Help")
        _help_btn.setFixedSize(60, 22)
        _help_btn.setStyleSheet(
            "QPushButton{background:#c0392b;color:white;font-size:11px;border-radius:3px;}"
            "QPushButton:hover{background:#e74c3c;}"
        )
        _help_btn.setToolTip("Open online documentation in your browser")
        _help_btn.clicked.connect(
            lambda: QDesktopServices.openUrl(QUrl(
                "https://github.com/jilavsky/pyirena/blob/main/docs/data_merge_gui.md"
            ))
        )
        row1.addWidget(_help_btn)
        vbox.addLayout(row1)

        # ── Row 2: Split option / Optimize button / Method ───────────────
        row2 = QHBoxLayout()
        row2.setSpacing(8)

        self._split_chk = QCheckBox("Split at left cursor")
        self._split_chk.setChecked(False)
        self._split_chk.setToolTip(
            "If checked, the final merged curve takes DS1 only below the left cursor "
            "and DS2 from the left cursor onward.  Default (unchecked): include all data "
            "from both datasets (interleaved in the overlap region)."
        )
        row2.addWidget(self._split_chk)

        row2.addWidget(_vline())

        self._optimize_btn = QPushButton("Optimize Merge")
        self._optimize_btn.setMinimumHeight(30)
        self._optimize_btn.setStyleSheet(_BTN_GREEN)
        self._optimize_btn.setEnabled(False)
        self._optimize_btn.clicked.connect(self._optimize_merge)
        row2.addWidget(self._optimize_btn)

        # Method combo (for future extensibility)
        self._method_combo = QComboBox()
        for key, label in DataMerge.METHODS.items():
            self._method_combo.addItem(label, key)
        self._method_combo.setToolTip("Optimisation method")
        row2.addWidget(self._method_combo)

        row2.addStretch()
        vbox.addLayout(row2)

        return vbox

    def _build_q_range_box(self) -> QGroupBox:
        box = QGroupBox("Overlap Q range")
        row = QHBoxLayout(box)
        row.setSpacing(8)

        row.addWidget(QLabel("Q min (Å⁻¹):"))
        self._q_min_edit = QLineEdit()
        self._q_min_edit.setReadOnly(True)
        self._q_min_edit.setFixedWidth(90)
        self._q_min_edit.setStyleSheet(_RDONLY_STYLE)
        row.addWidget(self._q_min_edit)

        row.addWidget(QLabel("Q max (Å⁻¹):"))
        self._q_max_edit = QLineEdit()
        self._q_max_edit.setReadOnly(True)
        self._q_max_edit.setFixedWidth(90)
        self._q_max_edit.setStyleSheet(_RDONLY_STYLE)
        row.addWidget(self._q_max_edit)

        hint = QLabel("← Drag cursors A / B on the graph to set the overlap range")
        hint.setStyleSheet("color: #7f8c8d; font-style: italic; font-size: 11px;")
        row.addWidget(hint)
        row.addStretch()

        # Update displays whenever cursors move
        # (connected lazily in _update_cursor_display() called from a timer or on data load)
        return box

    def _build_save_panel(self) -> QWidget:
        panel = QWidget()
        panel.setFixedWidth(210)
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)

        layout.addWidget(QLabel("Output folder:"))
        self._out_folder_label = QLabel("(not set)")
        self._out_folder_label.setWordWrap(True)
        self._out_folder_label.setStyleSheet("color: #2c3e50; font-size: 11px;")
        layout.addWidget(self._out_folder_label)

        create_folder_btn = QPushButton("Create Default Folder")
        create_folder_btn.setMinimumHeight(28)
        create_folder_btn.setStyleSheet(_BTN_GREY)
        create_folder_btn.setToolTip(
            "Create a new folder named '<DS1 folder>_merged' next to the DS1 folder "
            "and use it as the output folder."
        )
        create_folder_btn.clicked.connect(self._create_default_output_folder)
        layout.addWidget(create_folder_btn)

        select_folder_btn = QPushButton("Select Existing Folder…")
        select_folder_btn.setMinimumHeight(28)
        select_folder_btn.setStyleSheet(_BTN_GREY)
        select_folder_btn.setToolTip("Choose an existing folder for saving merged files.")
        select_folder_btn.clicked.connect(self._select_existing_output_folder)
        layout.addWidget(select_folder_btn)

        layout.addWidget(_hline())

        self._save_btn = QPushButton("Save Merged Data")
        self._save_btn.setMinimumHeight(34)
        self._save_btn.setStyleSheet(_BTN_GREEN)
        self._save_btn.setEnabled(False)
        self._save_btn.clicked.connect(self._save_merged)
        layout.addWidget(self._save_btn)

        layout.addWidget(_hline())

        self._batch_btn = QPushButton("Batch Run")
        self._batch_btn.setMinimumHeight(34)
        self._batch_btn.setStyleSheet(_BTN_BLUE)
        self._batch_btn.setEnabled(False)
        self._batch_btn.clicked.connect(self._batch_merge)
        layout.addWidget(self._batch_btn)

        layout.addWidget(_hline())

        save_cfg_btn = QPushButton("Save JSON Config…")
        save_cfg_btn.setMinimumHeight(28)
        save_cfg_btn.setStyleSheet(_BTN_GREY)
        save_cfg_btn.clicked.connect(self._save_json_config)
        layout.addWidget(save_cfg_btn)

        load_cfg_btn = QPushButton("Load JSON Config…")
        load_cfg_btn.setMinimumHeight(28)
        load_cfg_btn.setStyleSheet(_BTN_GREY)
        load_cfg_btn.clicked.connect(self._load_json_config)
        layout.addWidget(load_cfg_btn)

        layout.addStretch()
        return panel

    # ================================================================== #
    #  Public helpers (called by DataSelectorPanel.launch_data_merge())   #
    # ================================================================== #

    def set_folder(self, dataset: int, path: str) -> None:
        """Pre-populate a dataset folder (e.g. from the Data Selector hub)."""
        if dataset == 1:
            self._ds1.set_folder(path)
        else:
            self._ds2.set_folder(path)

    # ================================================================== #
    #  Event handlers                                                      #
    # ================================================================== #

    def _on_ds1_double_clicked(self, item: QListWidgetItem) -> None:
        fname = item.text()
        folder = self._ds1.current_folder
        if not folder:
            return
        data = self._load_file(os.path.join(folder, fname), self._ds1.get_file_type())
        if data is None:
            return
        self._data1 = data
        self._graph.plot_ds1(data['Q'], data['Intensity'], data.get('Error'))
        self._init_cursors_from_data()
        self._update_cursor_display()
        self._check_enable_buttons()
        self._status.setText(f"DS1 loaded: {fname}")

    def _on_ds2_double_clicked(self, item: QListWidgetItem) -> None:
        fname = item.text()
        folder = self._ds2.current_folder
        if not folder:
            return
        data = self._load_file(os.path.join(folder, fname), self._ds2.get_file_type())
        if data is None:
            return
        self._data2 = data
        self._graph.plot_ds2(data['Q'], data['Intensity'], data.get('Error'))
        self._init_cursors_from_data()
        self._update_cursor_display()
        self._check_enable_buttons()
        self._status.setText(f"DS2 loaded: {fname}")

    def _on_fit_scale_toggled(self, checked: bool) -> None:
        """When 'Fit scale' is unchecked, let the user type a fixed scale value."""
        self._scale_result.setReadOnly(checked)
        self._scale_result.setStyleSheet(_RDONLY_STYLE if checked else "")

    def _on_fit_qshift_toggled(self, checked: bool) -> None:
        """When 'Fit Q-shift' is unchecked, let the user type a fixed Q shift."""
        self._qshift_result.setReadOnly(checked)
        self._qshift_result.setStyleSheet(_RDONLY_STYLE if checked else "")

    def _on_mode_changed(self, _idx: int) -> None:
        saxs = self._mode_combo.currentIndex() == 0
        self._graph.set_mode(saxs)
        self._update_cursor_display()
        self.save_state()

    def _on_match_changed(self, _state) -> None:
        if self._match_chk.isChecked():
            self._apply_match_filter()
        else:
            # Restore full lists
            self._ds1._refresh_file_list()
            self._ds2._refresh_file_list()
        self._check_enable_batch_btn()

    # ================================================================== #
    #  Core operations                                                     #
    # ================================================================== #

    def _optimize_merge(self) -> None:
        if self._data1 is None or self._data2 is None:
            QMessageBox.warning(self, "Missing Data", "Load both DS1 and DS2 first.")
            return

        q_min, q_max = self._graph.get_overlap_range()
        if q_min is None:
            QMessageBox.warning(self, "No Cursors", "Drag the cursors to set the overlap region.")
            return
        if q_min >= q_max:
            QMessageBox.warning(
                self, "Invalid Q Range",
                f"Left cursor ({q_min:.4g}) must be to the left of the right cursor ({q_max:.4g})."
            )
            return

        # Validate that overlap region has sufficient data
        q1, I1 = self._data1['Q'], self._data1['Intensity']
        q2, I2 = self._data2['Q'], self._data2['Intensity']
        n1 = int(np.sum((q1 >= q_min) & (q1 <= q_max) & (I1 > 0)))
        n2 = int(np.sum((q2 >= q_min) & (q2 <= q_max) & (I2 > 0)))
        if n1 < 2:
            QMessageBox.warning(
                self, "Insufficient DS1 Data",
                f"Only {n1} DS1 point(s) in the overlap Q range [{q_min:.4g}, {q_max:.4g}]. "
                "Move the cursors to a region with more DS1 data."
            )
            return
        if n2 < 1:
            QMessageBox.warning(
                self, "Insufficient DS2 Data",
                f"No DS2 points in the overlap Q range [{q_min:.4g}, {q_max:.4g}]. "
                "Move the cursors to a region with DS2 data."
            )
            return

        # Build MergeConfig from UI
        qshift_map = {"None": 0, "DS1": 1, "DS2": 2}
        try:
            fixed_scale = float(self._scale_result.text())
        except ValueError:
            fixed_scale = 1.0
        try:
            fixed_qshift = float(self._qshift_result.text())
        except ValueError:
            fixed_qshift = 0.0
        config = MergeConfig(
            q_overlap_min=q_min,
            q_overlap_max=q_max,
            fit_scale=self._fit_scale_chk.isChecked(),
            scale_dataset=self._scale_ds_combo.currentIndex() + 1,   # 0→1, 1→2
            fixed_scale_value=fixed_scale,
            fit_qshift=self._fit_qshift_chk.isChecked(),
            fixed_qshift_value=fixed_qshift,
            qshift_dataset=qshift_map[self._qshift_combo.currentText()],
            method=self._method_combo.currentData() or 'interpolation',
            split_at_left_cursor=self._split_chk.isChecked(),
        )

        _dI1 = self._data1.get('Error')
        dI1 = _dI1 if _dI1 is not None else I1 * 0.05
        _dI2 = self._data2.get('Error')
        dI2 = _dI2 if _dI2 is not None else I2 * 0.05

        self._status.setText("Optimising… please wait.")
        QApplication.processEvents()

        result = self._engine.optimize(q1, I1, dI1, q2, I2, dI2, config)
        self._last_result = result
        self._last_config = config

        # Update result displays
        self._scale_result.setText(f"{result.scale:.5g}")
        self._qshift_result.setText(f"{result.q_shift:.6g}")
        self._bg_result.setText(f"{result.background:.6g}")

        status_txt = (
            f"Optimise {'OK' if result.success else 'FAILED'} — "
            f"scale={result.scale:.4g}  BG={result.background:.4g}  "
            f"Q-shift={result.q_shift:.5g}  "
            f"χ²={result.chi_squared:.4g}  n_pts={result.n_overlap_points}"
        )
        self._status.setText(status_txt)

        # Compute merged arrays and plot preview
        dQ1 = self._data1.get('dQ')
        dQ2 = self._data2.get('dQ')
        q_m, I_m, dI_m, dQ_m = self._engine.merge(
            q1, I1, dI1, dQ1, q2, I2, dI2, dQ2, result, config
        )
        self._last_q_merged = q_m
        self._last_I_merged = I_m
        self._last_dI_merged = dI_m
        self._last_dQ_merged = dQ_m

        self._graph.plot_merged(q_m, I_m)
        self._save_btn.setEnabled(self._out_folder is not None)
        self.save_state()

    def _save_merged(self) -> None:
        if self._last_q_merged is None:
            QMessageBox.warning(self, "No Result", "Run 'Optimize Merge' first.")
            return
        if not self._out_folder:
            QMessageBox.warning(self, "No Output Folder", "Select an output folder first.")
            return

        from pyirena.io.nxcansas_data_merge import save_merged_data

        ds1_path = Path(self._data1['filepath'])
        ds2_path = Path(self._data2['filepath']) if self._data2 else None

        merge_result_dict = self._build_provenance_dict()
        try:
            out = save_merged_data(
                output_folder=Path(self._out_folder),
                ds1_path=ds1_path,
                ds1_is_nxcansas=self._data1.get('is_nxcansas', False),
                q=self._last_q_merged,
                I=self._last_I_merged,
                dI=self._last_dI_merged,
                dQ=self._last_dQ_merged,
                merge_result_dict=merge_result_dict,
                ds2_path=ds2_path,
            )
            self._status.setText(f"Saved → {out}")
        except Exception as exc:
            QMessageBox.critical(self, "Save Error", str(exc))
            self._status.setText(f"Save error: {exc}")

    def _batch_merge(self) -> None:
        """Process file pairs sequentially.

        Pair source priority:
        1. Selected files — if both lists have selections, zip them in order
           (first selected DS1 with first selected DS2, etc.).
        2. Matched pairs — fall back to the 'Match files' auto-match.
        """
        sel1 = self._ds1.get_selected_filenames()
        sel2 = self._ds2.get_selected_filenames()
        if sel1 and sel2:
            pairs = list(zip(sel1, sel2))
        else:
            pairs = self._get_matched_pairs()
        if not pairs:
            QMessageBox.warning(self, "No Pairs",
                                "Select files in both lists, or enable 'Match files'.")
            return
        if not self._out_folder:
            QMessageBox.warning(self, "No Output Folder", "Select an output folder first.")
            return

        from pyirena.io.nxcansas_data_merge import save_merged_data

        qshift_map = {"None": 0, "DS1": 1, "DS2": 2}
        q_min, q_max = self._graph.get_overlap_range()
        if q_min is None or q_min >= q_max:
            QMessageBox.warning(self, "Q Range",
                                "Set the overlap Q range with cursors before batch run.")
            return

        n_ok = 0
        n_fail = 0
        for f1, f2 in pairs:
            path1 = os.path.join(self._ds1.current_folder, f1)
            path2 = os.path.join(self._ds2.current_folder, f2)
            d1 = self._load_file(path1, self._ds1.get_file_type())
            d2 = self._load_file(path2, self._ds2.get_file_type())
            if d1 is None or d2 is None:
                print(f"[batch] Skipping {f1} / {f2} — load error.")
                n_fail += 1
                continue

            q1, I1 = d1['Q'], d1['Intensity']
            q2, I2 = d2['Q'], d2['Intensity']
            _dI1 = d1.get('Error')
            dI1 = _dI1 if _dI1 is not None else I1 * 0.05
            _dI2 = d2.get('Error')
            dI2 = _dI2 if _dI2 is not None else I2 * 0.05

            try:
                fixed_scale = float(self._scale_result.text())
            except ValueError:
                fixed_scale = 1.0
            try:
                fixed_qshift = float(self._qshift_result.text())
            except ValueError:
                fixed_qshift = 0.0
            config = MergeConfig(
                q_overlap_min=q_min, q_overlap_max=q_max,
                fit_scale=self._fit_scale_chk.isChecked(),
                scale_dataset=self._scale_ds_combo.currentIndex() + 1,
                fixed_scale_value=fixed_scale,
                fit_qshift=self._fit_qshift_chk.isChecked(),
                fixed_qshift_value=fixed_qshift,
                qshift_dataset=qshift_map[self._qshift_combo.currentText()],
                method=self._method_combo.currentData() or 'interpolation',
                split_at_left_cursor=self._split_chk.isChecked(),
            )
            result = self._engine.optimize(q1, I1, dI1, q2, I2, dI2, config)
            q_m, I_m, dI_m, dQ_m = self._engine.merge(
                q1, I1, dI1, d1.get('dQ'), q2, I2, dI2, d2.get('dQ'), result, config
            )
            prov = {
                'scale': result.scale, 'q_shift': result.q_shift,
                'background': result.background, 'chi_squared': result.chi_squared,
                'n_overlap_points': result.n_overlap_points,
                'q_overlap_min': config.q_overlap_min, 'q_overlap_max': config.q_overlap_max,
                'scale_dataset': config.scale_dataset, 'fit_scale': config.fit_scale,
                'qshift_dataset': config.qshift_dataset, 'fit_qshift': config.fit_qshift,
                'split_at_left_cursor': config.split_at_left_cursor,
            }
            try:
                out = save_merged_data(
                    output_folder=Path(self._out_folder),
                    ds1_path=Path(path1),
                    ds1_is_nxcansas=d1.get('is_nxcansas', False),
                    q=q_m, I=I_m, dI=dI_m, dQ=dQ_m,
                    merge_result_dict=prov,
                    ds2_path=Path(path2),
                )
                status = "OK" if result.success else "FAILED"
                print(f"[batch] {status}  {f1} + {f2} → {out.name}  "
                      f"scale={result.scale:.4g}  BG={result.background:.4g}")
                n_ok += 1
            except Exception as exc:
                print(f"[batch] ERROR {f1} + {f2}: {exc}")
                n_fail += 1

            QApplication.processEvents()

        msg = f"Batch done: {n_ok} OK, {n_fail} failed."
        self._status.setText(msg)
        QMessageBox.information(self, "Batch Complete", msg)

    # ================================================================== #
    #  File matching                                                       #
    # ================================================================== #

    def _apply_match_filter(self) -> None:
        files1 = self._ds1.get_all_files()
        files2 = self._ds2.get_all_files()
        pairs = self._engine.match_files(files1, files2)
        matched1 = [p[0] for p in pairs]
        matched2 = [p[1] for p in pairs]
        self._ds1.set_displayed_files(matched1)
        self._ds2.set_displayed_files(matched2)
        self._status.setText(
            f"Match mode: {len(pairs)} paired file(s) shown. "
            "Lines correspond to matched pairs."
        )

    def _get_matched_pairs(self) -> List[Tuple[str, str]]:
        if not self._match_chk.isChecked():
            return []
        files1 = self._ds1.get_all_files()
        files2 = self._ds2.get_all_files()
        return self._engine.match_files(files1, files2)

    # ================================================================== #
    #  JSON config save/load                                              #
    # ================================================================== #

    def _save_json_config(self) -> None:
        start_dir = str(Path(self._ds1.current_folder or Path.home()).parent)
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Merge Config", os.path.join(start_dir, "merge_config.json"),
            "JSON files (*.json)"
        )
        if not path:
            return
        cfg = self._current_config_dict()
        try:
            with open(path, 'w') as fh:
                json.dump(cfg, fh, indent=2)
            self._status.setText(f"Config saved → {os.path.basename(path)}")
        except Exception as exc:
            QMessageBox.critical(self, "Save Error", str(exc))

    def _load_json_config(self, path: Optional[str] = None) -> None:
        if path is None:
            start_dir = str(Path(self._ds1.current_folder or Path.home()).parent)
            path, _ = QFileDialog.getOpenFileName(
                self, "Load Merge Config", start_dir, "JSON files (*.json)"
            )
        if not path:
            return
        try:
            with open(path) as fh:
                cfg = json.load(fh)
        except Exception as exc:
            QMessageBox.critical(self, "Load Error", str(exc))
            return
        self._apply_config_dict(cfg)
        self._status.setText(f"Config loaded from {os.path.basename(path)}")

    # ================================================================== #
    #  State persistence                                                   #
    # ================================================================== #

    def load_state(self) -> None:
        s = self._sm.get('data_merge') or {}

        folder1 = s.get('folder1')
        folder2 = s.get('folder2')
        if folder1 and os.path.isdir(folder1):
            self._ds1.type_combo.setCurrentText(s.get('file_type1', 'HDF5 Nexus'))
            self._ds1.set_folder(folder1)
        if folder2 and os.path.isdir(folder2):
            self._ds2.type_combo.setCurrentText(s.get('file_type2', 'HDF5 Nexus'))
            self._ds2.set_folder(folder2)

        out = s.get('output_folder')
        if out and os.path.isdir(out):
            self._out_folder = out
            self._out_folder_label.setText(os.path.basename(out))
            self._out_folder_label.setToolTip(out)
        else:
            self._out_folder: Optional[str] = None

        self._fit_scale_chk.setChecked(bool(s.get('fit_scale', True)))
        idx = s.get('scale_dataset', 2)
        self._scale_ds_combo.setCurrentIndex(int(idx) - 1)
        qshift_ds = int(s.get('qshift_dataset', 0))
        self._qshift_combo.setCurrentIndex(qshift_ds)
        self._fit_qshift_chk.setChecked(bool(s.get('fit_qshift', False)))
        self._split_chk.setChecked(bool(s.get('split_at_left_cursor', False)))

        mode = s.get('plot_mode', 'saxs')
        self._mode_combo.setCurrentIndex(0 if mode == 'saxs' else 1)

        self._match_chk.setChecked(bool(s.get('match_files', False)))

        self._ds1.filter_edit.setText(s.get('filter1', ''))
        self._ds2.filter_edit.setText(s.get('filter2', ''))

    def save_state(self) -> None:
        q_min, q_max = self._graph.get_overlap_range()
        self._sm.update('data_merge', {
            'folder1': self._ds1.current_folder,
            'folder2': self._ds2.current_folder,
            'file_type1': self._ds1.get_file_type(),
            'file_type2': self._ds2.get_file_type(),
            'filter1': self._ds1.filter_edit.text(),
            'filter2': self._ds2.filter_edit.text(),
            'output_folder': self._out_folder,
            'q_overlap_min': q_min,
            'q_overlap_max': q_max,
            'fit_scale': self._fit_scale_chk.isChecked(),
            'scale_dataset': self._scale_ds_combo.currentIndex() + 1,
            'qshift_dataset': self._qshift_combo.currentIndex(),
            'fit_qshift': self._fit_qshift_chk.isChecked(),
            'split_at_left_cursor': self._split_chk.isChecked(),
            'plot_mode': 'saxs' if self._mode_combo.currentIndex() == 0 else 'waxs',
            'match_files': self._match_chk.isChecked(),
        })
        self._sm.save()

    def closeEvent(self, event) -> None:
        self.save_state()
        super().closeEvent(event)

    # ================================================================== #
    #  Private helpers                                                     #
    # ================================================================== #

    def _load_file(self, filepath: str, file_type: str) -> Optional[dict]:
        """Load a SAS data file according to the selected file type."""
        from pyirena.io.hdf5 import readGenericNXcanSAS, readSimpleHDF5, readTextFile
        fp = Path(filepath)
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
            QMessageBox.warning(self, "Load Error",
                                f"Could not read {fp.name}:\n{exc}")
            return None

    def _init_cursors_from_data(self) -> None:
        """Place cursors at the intersection of DS1 and DS2 Q ranges."""
        if self._data1 is None or self._data2 is None:
            return
        q1, q2 = self._data1['Q'], self._data2['Q']
        q_ov_lo = max(float(q1.min()), float(q2.min()))
        q_ov_hi = min(float(q1.max()), float(q2.max()))
        if q_ov_lo >= q_ov_hi:
            # No real overlap — place cursors 20% inside each dataset's range
            q_ov_lo = float(q1.max()) * 0.8
            q_ov_hi = float(q2.min()) * 1.2
        # Trim 10% on each side to avoid edge artefacts
        span = q_ov_hi - q_ov_lo
        q_min = q_ov_lo + 0.1 * span
        q_max = q_ov_hi - 0.1 * span

        # Check if state has a saved Q range
        s = self._sm.get('data_merge') or {}
        saved_min = s.get('q_overlap_min')
        saved_max = s.get('q_overlap_max')
        if saved_min is not None and saved_max is not None:
            q_min, q_max = float(saved_min), float(saved_max)

        # Don't reinitialise if cursors are already placed (user adjusted them)
        if self._graph._cursor_a is None:
            self._graph.init_cursors(q_min, q_max)
        self._update_cursor_display()

        # Connect cursor movement to display update
        if self._graph._cursor_a is not None:
            try:
                self._graph._cursor_a.sigPositionChanged.connect(
                    lambda: self._update_cursor_display()
                )
                self._graph._cursor_b.sigPositionChanged.connect(
                    lambda: self._update_cursor_display()
                )
            except Exception:
                pass

    def _update_cursor_display(self) -> None:
        q_min, q_max = self._graph.get_overlap_range()
        if q_min is not None:
            self._q_min_edit.setText(f"{q_min:.5g}")
            self._q_max_edit.setText(f"{q_max:.5g}")

    def _on_ds1_folder_changed(self, _folder: str) -> None:
        """Called when the user selects a new DS1 folder via the button."""
        # Clear output folder — the old path is no longer relevant for the new dataset.
        self._out_folder = None
        self._out_folder_label.setText("(not set)")
        self._out_folder_label.setToolTip("")
        # Uncheck Match — existing pairs are invalid after a folder change.
        if self._match_chk.isChecked():
            self._match_chk.setChecked(False)

    def _on_ds2_folder_changed(self, _folder: str) -> None:
        """Called when the user selects a new DS2 folder via the button."""
        # Uncheck Match — existing pairs are invalid after a folder change.
        if self._match_chk.isChecked():
            self._match_chk.setChecked(False)

    def _create_default_output_folder(self) -> None:
        """Create '<DS1 folder name>_merged' next to the DS1 folder and set it."""
        folder1 = self._ds1.current_folder
        if not folder1:
            QMessageBox.warning(self, "No DS1 Folder",
                                "Select a DS1 folder first.")
            return
        p = Path(folder1)
        default = p.parent / f"{p.name}_merged"
        try:
            default.mkdir(parents=True, exist_ok=True)
        except Exception as exc:
            QMessageBox.critical(self, "Create Folder Error",
                                 f"Could not create {default}:\n{exc}")
            return
        self._out_folder = str(default)
        self._out_folder_label.setText(default.name)
        self._out_folder_label.setToolTip(str(default))
        self._save_btn.setEnabled(self._last_q_merged is not None)
        self._status.setText(f"Output folder: {default.name}")
        self.save_state()

    def _select_existing_output_folder(self) -> None:
        """Open a folder-selection dialog to choose an existing output folder."""
        if self._out_folder:
            start = self._out_folder
        elif self._ds1.current_folder:
            start = str(Path(self._ds1.current_folder).parent)
        else:
            start = str(Path.home())
        folder = QFileDialog.getExistingDirectory(
            self, "Select Output Folder", start
        )
        if folder:
            self._out_folder = folder
            self._out_folder_label.setText(os.path.basename(folder))
            self._out_folder_label.setToolTip(folder)
            self._save_btn.setEnabled(self._last_q_merged is not None)
            self.save_state()

    def _check_enable_buttons(self) -> None:
        both_loaded = self._data1 is not None and self._data2 is not None
        self._optimize_btn.setEnabled(both_loaded)
        self._check_enable_batch_btn()

    def _check_enable_batch_btn(self) -> None:
        has_folders = (
            self._ds1.current_folder is not None
            and self._ds2.current_folder is not None
        )
        has_match = self._match_chk.isChecked()
        has_selection = (
            len(self._ds1.file_list.selectedItems()) > 0
            and len(self._ds2.file_list.selectedItems()) > 0
        )
        self._batch_btn.setEnabled(has_folders and (has_match or has_selection))

    def _build_provenance_dict(self) -> dict:
        r = self._last_result
        c = self._last_config
        if r is None or c is None:
            return {}
        return {
            'scale': r.scale, 'q_shift': r.q_shift, 'background': r.background,
            'chi_squared': r.chi_squared, 'n_overlap_points': r.n_overlap_points,
            'q_overlap_min': c.q_overlap_min, 'q_overlap_max': c.q_overlap_max,
            'scale_dataset': c.scale_dataset, 'fit_scale': c.fit_scale,
            'qshift_dataset': c.qshift_dataset, 'fit_qshift': c.fit_qshift,
            'split_at_left_cursor': c.split_at_left_cursor,
        }

    def _current_config_dict(self) -> dict:
        q_min, q_max = self._graph.get_overlap_range()
        qshift_map = {"None": 0, "DS1": 1, "DS2": 2}
        return {
            "version": "1.0",
            "q_overlap_min": q_min,
            "q_overlap_max": q_max,
            "fit_scale": self._fit_scale_chk.isChecked(),
            "scale_dataset": self._scale_ds_combo.currentIndex() + 1,
            "fit_qshift": self._fit_qshift_chk.isChecked(),
            "qshift_dataset": qshift_map[self._qshift_combo.currentText()],
            "split_at_left_cursor": self._split_chk.isChecked(),
        }

    def _apply_config_dict(self, cfg: dict) -> None:
        qds_map = {0: "None", 1: "DS1", 2: "DS2"}
        self._fit_scale_chk.setChecked(bool(cfg.get('fit_scale', True)))
        sds = int(cfg.get('scale_dataset', 2))
        self._scale_ds_combo.setCurrentIndex(sds - 1)
        self._fit_qshift_chk.setChecked(bool(cfg.get('fit_qshift', False)))
        qds = int(cfg.get('qshift_dataset', 0))
        self._qshift_combo.setCurrentText(qds_map.get(qds, "None"))
        self._split_chk.setChecked(bool(cfg.get('split_at_left_cursor', False)))

        q_min = cfg.get('q_overlap_min')
        q_max = cfg.get('q_overlap_max')
        if q_min is not None and q_max is not None:
            self._graph.init_cursors(float(q_min), float(q_max))
            self._update_cursor_display()


# ---------------------------------------------------------------------------
# Utility widget helpers
# ---------------------------------------------------------------------------

def _vline() -> QFrame:
    """Vertical separator line for use in QHBoxLayout."""
    sep = QFrame()
    sep.setFrameShape(QFrame.Shape.VLine)
    sep.setFrameShadow(QFrame.Shadow.Sunken)
    return sep


def _hline() -> QFrame:
    """Horizontal separator line for use in QVBoxLayout."""
    sep = QFrame()
    sep.setFrameShape(QFrame.Shape.HLine)
    sep.setFrameShadow(QFrame.Shadow.Sunken)
    return sep


# ===========================================================================
# CLI entry point
# ===========================================================================

def main() -> None:
    """Entry point: ``pyirena-datamerge``.

    With no arguments: launches the full GUI.
    With ``--file1`` and ``--file2``: runs a headless merge (no Qt needed).
    With ``--folder1`` and/or ``--folder2``: launches GUI with pre-filled folders.
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog='pyirena-datamerge',
        description='Merge two SAS datasets (USAXS+SAXS, SAXS+WAXS, etc.).',
    )
    parser.add_argument('--folder1', metavar='DIR',
                        help='DS1 folder (lower Q).  Pre-populates GUI.')
    parser.add_argument('--folder2', metavar='DIR',
                        help='DS2 folder (higher Q).  Pre-populates GUI.')
    parser.add_argument('--file1', metavar='FILE',
                        help='DS1 single file — triggers headless merge mode.')
    parser.add_argument('--file2', metavar='FILE',
                        help='DS2 single file — triggers headless merge mode.')
    parser.add_argument('--config', metavar='JSON',
                        help='JSON config file with merge parameters.')
    parser.add_argument('--match', action='store_true',
                        help='Enable file-matching mode in GUI.')
    parser.add_argument('--output', metavar='DIR',
                        help='Output folder for headless mode.')

    args = parser.parse_args()

    if args.file1 and args.file2:
        # Headless mode — no Qt
        from pyirena.batch import merge_data
        result = merge_data(
            args.file1, args.file2,
            config_file=args.config,
            output_folder=args.output,
        )
        sys.exit(0 if (result and result.get('success')) else 1)

    # GUI mode
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
        app.setStyle('Fusion')

    window = DataMergePanel()
    if args.folder1:
        window.set_folder(1, args.folder1)
    if args.folder2:
        window.set_folder(2, args.folder2)
    if args.match:
        window._match_chk.setChecked(True)

    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
