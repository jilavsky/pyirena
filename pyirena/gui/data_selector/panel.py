"""
pyirena.gui.data_selector.panel — the DataSelectorPanel main widget and main() entry point.

Split from the original monolithic data_selector.py (no behavior change).
"""
import logging

log = logging.getLogger(__name__)


import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pyqtgraph as pg

from pyirena.gui.data_selector._qt import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QPushButton, QListWidget, QLabel, QLineEdit, QFileDialog, QComboBox, QAbstractItemView, QMessageBox, QMenuBar, QMenu, QFrame, QScrollArea, QDialog, QDialogButtonBox, QGroupBox, QCheckBox, Qt, QDir, QUrl, QAction, QDesktopServices,
)
from pyirena.io.hdf5 import readGenericNXcanSAS
from pyirena.io.text_import import ensure_nxcansas_sibling
from pyirena.gui.data_loading import (
    prompt_dataset_choice as _prompt_dataset_choice_fn,
    read_nxcansas_with_picker as _read_nxcansas_with_picker_fn,
    load_data_file as _load_data_file_fn,
)
from pyirena.io.nxcansas_unified import load_unified_fit_results
from pyirena.gui.unified_fit import UnifiedFitPanel
from pyirena.gui.sizes_panel import SizesFitPanel
from pyirena.state import StateManager

from pyirena.gui.data_selector.config_dialogs import ConfigManagerDialog, DataSelectorConfigDialog
from pyirena.gui.data_selector.igor_import import _IgorImportDialog
from pyirena.gui.data_selector.plot_utils import _gen_colors, _legend_indices
from pyirena.gui.data_selector.report import _build_report
from pyirena.gui.data_selector.results_windows import GraphWindow, SimpleFitResultsWindow, SizeDistResultsWindow, TabulateResultsWindow, UnifiedFitResultsWindow, WAXSPeakFitResultsWindow
from pyirena.gui.data_selector.sorting import _SORT_KEYS
from pyirena.gui.data_selector.workers import BatchWorker


class DataSelectorPanel(QWidget):
    """
    Main data selector panel for pyIrena.

    Provides file browsing, filtering, and data visualization capabilities.
    """

    # File extensions for different data types
    HDF5_EXTENSIONS = ['.hdf', '.h5', '.hdf5']
    TEXT_EXTENSIONS = ['.txt', '.dat']

    def __init__(self):
        super().__init__()
        self.current_folder = None
        self.last_folder = None  # Remember last selected folder
        self.graph_window = None
        self.unified_fit_results_window = None  # Graph of stored fit results
        self.size_dist_results_window = None   # Graph of stored size dist results
        self.tabulate_results_window = None    # Tabulated results window
        self.unified_fit_window = None  # Unified fit panel
        self.modeling_window = None    # Modeling panel
        self.sizes_fit_window = None   # Size distribution panel
        self.simple_fits_window = None         # Simple Fits panel
        self.simple_fits_results_window = None # Graph of stored simple fit results
        self.waxs_peakfit_window = None        # WAXS Peak Fit panel
        self.waxs_peakfit_results_window = None  # Graph of stored WAXS peak-fit results
        self.hdf5_viewer_window = None         # Data Explorer window
        self.data_merge_window = None          # Data Merge panel
        self.data_manip_window = None          # Data Manipulation panel
        self.contrast_window = None            # Scattering Contrast Calculator
        self.fractals_window = None            # Fractals (mass fractal aggregate) panel
        self.saxs_morph_window = None          # SAXS Morph (3D voxelgram) panel
        # Standalone read-only voxel viewers spawned by Create Graph with
        # the Fractals / 3D-saxsMorph checkboxes.  We keep references so
        # Python does not GC them; each viewer removes itself on close.
        self._standalone_voxel_viewers: list = []
        self._batch_worker = None      # Batch fitting thread

        # Initialize state manager
        self.state_manager = StateManager()

        # Load last used folder from state
        self.load_last_folder()

        self.init_ui()

        # Restore saved sort selection and re-sort the already-populated list
        saved_sort = int(self.state_manager.get('data_selector', 'sort_index', 0) or 0)
        self.sort_combo.blockSignals(True)
        self.sort_combo.setCurrentIndex(saved_sort)
        self.sort_combo.blockSignals(False)
        self.sort_file_list()   # apply restored sort order to the initial file list

    def init_ui(self):
        """Initialize the user interface."""
        self.setWindowTitle("pyIrena - Data Selector")

        # Main layout
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Menu bar
        menu_bar = self.create_menu_bar()
        main_layout.addWidget(menu_bar)

        # Content layout (with margins)
        content_layout = QVBoxLayout()
        content_layout.setContentsMargins(20, 20, 20, 20)
        content_layout.setSpacing(15)

        # Title row — centred title with Help button at top-right
        title_row = QHBoxLayout()
        title_row.setContentsMargins(0, 0, 0, 0)
        title_row.setSpacing(0)
        # Invisible left spacer balances the help button so the title stays centred
        _title_spacer = QWidget()
        _title_spacer.setFixedWidth(70)
        title_row.addWidget(_title_spacer)
        title_label = QLabel("pyIrena")
        title_label.setStyleSheet("""
            QLabel {
                font-size: 18px;
                font-weight: bold;
                color: #2c3e50;
                padding: 4px;
            }
        """)
        title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title_row.addWidget(title_label, 1)
        _help_btn = QPushButton("? Help")
        _help_btn.setFixedSize(60, 22)
        _help_btn.setStyleSheet(
            "QPushButton{background:#c0392b;color:white;font-size:11px;border-radius:3px;}"
            "QPushButton:hover{background:#e74c3c;}"
        )
        _help_btn.setToolTip("Open online documentation in your browser")
        _help_btn.clicked.connect(
            lambda: QDesktopServices.openUrl(QUrl(
                "https://github.com/jilavsky/pyirena/blob/main/docs/GUI_README.md"
            ))
        )
        title_row.addWidget(_help_btn)
        title_row.addSpacing(5)
        content_layout.addLayout(title_row)

        # Folder selection section
        folder_layout = QHBoxLayout()
        self.folder_button = QPushButton("Select Data Folder")
        self.folder_button.setMinimumHeight(28)
        self.folder_button.setToolTip("Browse to a folder containing NXcanSAS/HDF5 data files.")
        self.folder_button.clicked.connect(self.select_folder)
        folder_layout.addWidget(self.folder_button)

        self.refresh_button = QPushButton("Refresh")
        self.refresh_button.setMinimumHeight(28)
        self.refresh_button.setMaximumWidth(100)
        self.refresh_button.setToolTip(
            "Re-scan the current folder and update the file list.\n"
            "Use after adding or removing files from the folder."
        )
        self.refresh_button.clicked.connect(self.refresh_file_list)
        self.refresh_button.setEnabled(False)
        folder_layout.addWidget(self.refresh_button)

        self.folder_label = QLabel("No folder selected")
        self.folder_label.setStyleSheet("color: #7f8c8d; font-style: italic;")
        folder_layout.addWidget(self.folder_label)
        folder_layout.addStretch()

        content_layout.addLayout(folder_layout)

        # File type + Sort selection (combined row)
        type_layout = QHBoxLayout()
        type_layout.addWidget(QLabel("File Type:"))

        self.file_type_combo = QComboBox()
        self.file_type_combo.addItem("HDF5 Files (.hdf, .h5, .hdf5)", "hdf5")
        self.file_type_combo.addItem("Text Files (.txt, .dat)", "text")
        self.file_type_combo.addItem("All Supported Files", "all")
        self.file_type_combo.currentIndexChanged.connect(self.refresh_file_list)
        self.file_type_combo.setMaximumWidth(140)
        type_layout.addWidget(self.file_type_combo)

        type_layout.addSpacing(12)
        type_layout.addWidget(QLabel("Sort:"))
        self.sort_combo = QComboBox()
        self.sort_combo.addItems([
            "Filename  A→Z",
            "Filename  Z→A",
            "Temperature  ↑",
            "Temperature  ↓",
            "Time  ↑",
            "Time  ↓",
            "Order number  ↑",
            "Order number  ↓",
            "Pressure  ↑",
            "Pressure  ↓",
        ])
        self.sort_combo.setToolTip(
            "Sort order for the file list.\n"
            "Patterns recognised in filenames:\n"
            "  Temperature : _25C\n"
            "  Time        : _50min\n"
            "  Order number: _354  (last underscore-number before extension)\n"
            "  Pressure    : _35PSI"
        )
        self.sort_combo.currentIndexChanged.connect(self._on_sort_changed)
        self.sort_combo.setMaximumWidth(140)
        type_layout.addWidget(self.sort_combo)
        type_layout.addStretch()

        # Content area (listbox + graph button)
        file_area_layout = QHBoxLayout()

        # Left side: file list section
        left_layout = QVBoxLayout()

        # File Type + Sort row lives at the top of the left column so the right
        # column's tool palette can extend up to the same vertical level —
        # gives the right column more usable height before scrolling kicks in
        # on high-DPI / scaled displays.
        left_layout.addLayout(type_layout)

        # File filter input
        filter_layout = QHBoxLayout()
        filter_layout.addWidget(QLabel("Filter:"))
        self.filter_input = QLineEdit()
        self.filter_input.setPlaceholderText("Enter text to filter files...")
        self.filter_input.setMaximumWidth(350)
        self.filter_input.textChanged.connect(self.filter_files)
        filter_layout.addWidget(self.filter_input)
        filter_layout.addStretch()
        left_layout.addLayout(filter_layout)

        # File list
        self.file_list = QListWidget()
        self.file_list.setMinimumWidth(300)  # ~22 characters wide
        # Cap the list's preferred width so long filenames don't bloat the
        # left column and squeeze the right column. Users can scroll
        # horizontally or hover for the full name; the right column needs
        # the space more.
        self.file_list.setMaximumWidth(400)
        self.file_list.setMinimumHeight(400)  # Show at least 15 items
        self.file_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.file_list.itemDoubleClicked.connect(self.plot_selected_files)
        self.file_list.itemSelectionChanged.connect(self.update_plot_button_state)
        left_layout.addWidget(self.file_list)

        # Configure + Manage Config — small buttons below the file list
        configure_row = QHBoxLayout()
        configure_row.setSpacing(4)
        self.configure_button = QPushButton("Configure...")
        self.configure_button.setMaximumWidth(110)
        self.configure_button.setMinimumHeight(24)
        self.configure_button.setToolTip("Configure data loading options")
        self.configure_button.clicked.connect(self.open_configure_dialog)
        configure_row.addWidget(self.configure_button)
        _mc_style = (
            "QPushButton { background:#7f8c8d; color:white; font-size:10px; "
            "border-radius:3px; padding:2px 6px; }"
            "QPushButton:hover { background:#95a5a6; }"
        )
        self.manage_config_button = QPushButton("Manage Config…")
        self.manage_config_button.setMinimumHeight(24)
        self.manage_config_button.setStyleSheet(_mc_style)
        self.manage_config_button.setToolTip(
            "Inspect and remove tool sections from pyirena_config.json.\n"
            "Use this to prevent unwanted tools from running during batch fitting."
        )
        self.manage_config_button.clicked.connect(self.open_config_manager)
        configure_row.addWidget(self.manage_config_button)
        self.locate_logs_button = QPushButton("Locate Logs…")
        self.locate_logs_button.setMinimumHeight(24)
        self.locate_logs_button.setStyleSheet(_mc_style)
        self.locate_logs_button.setToolTip(
            "Open the pyirena log folder (~/.pyirena/logs).\n"
            "Attach the newest gui.log when reporting a problem."
        )
        self.locate_logs_button.clicked.connect(self.locate_logs)
        configure_row.addWidget(self.locate_logs_button)
        configure_row.addStretch()
        left_layout.addLayout(configure_row)

        # Left column takes only the width its widest capped child needs
        # (~400px from the file list). All extra horizontal space goes to
        # the right column so the two columns sit flush against each other
        # — no gap between the file list and the right-column scroll area.
        file_area_layout.addLayout(left_layout, stretch=0)

        # Right side: action buttons — starts at same level as Filter row
        right_layout = QVBoxLayout()
        right_layout.setSpacing(6)

        # ── Batch script status display ────────────────────────────────────
        self.batch_status_label = QLabel("")
        self.batch_status_label.setWordWrap(True)
        self.batch_status_label.setMinimumHeight(23)
        self.batch_status_label.setMaximumHeight(56)
        self.batch_status_label.setVisible(False)
        self.batch_status_label.setStyleSheet(
            "QLabel { padding: 5px 8px; border-radius: 4px; font-size: 12px; }"
        )
        right_layout.addWidget(self.batch_status_label)

        # ── GROUP 1: View & Export ─────────────────────────────────────────────
        grp_view = QGroupBox("")
        grp_view_lay = QVBoxLayout(grp_view)
        grp_view_lay.setSpacing(4)

        # Checkboxes — control what graph / reports show
        _cb_lbl = QLabel("Show in graph / reports:")
        _cb_lbl.setStyleSheet("font-weight:bold; color:#2c3e50; font-size:11px;")
        grp_view_lay.addWidget(_cb_lbl)

        # Checkboxes split across two rows (4+4) so they're not crowded —
        # this also leaves vertical room while keeping the right column
        # narrow enough for the 60%-width button grid below to fit.
        cb_row1 = QHBoxLayout()
        cb_row1.setSpacing(8)
        cb_row2 = QHBoxLayout()
        cb_row2.setSpacing(8)

        self.data_checkbox = QCheckBox("Data")
        self.data_checkbox.setChecked(True)
        self.data_checkbox.setToolTip("Plot experimental data for selected files")
        cb_row1.addWidget(self.data_checkbox)

        self.unified_fit_result_checkbox = QCheckBox("Unified Fit")
        self.unified_fit_result_checkbox.setChecked(False)
        self.unified_fit_result_checkbox.setToolTip(
            "Plot stored Unified Fit results (data + model + residuals).\n"
            "Only HDF5 files are checked; files without fit results are skipped."
        )
        cb_row1.addWidget(self.unified_fit_result_checkbox)

        self.size_dist_checkbox = QCheckBox("Size Dist.")
        self.size_dist_checkbox.setChecked(False)
        self.size_dist_checkbox.setToolTip(
            "Open Size Distribution panel with selected data (Create Graph),\n"
            "or include stored sizes results in report (Create Report).\n"
            "Only HDF5 files with stored sizes results are used for reports."
        )
        cb_row1.addWidget(self.size_dist_checkbox)

        self.simple_fits_checkbox = QCheckBox("Simple Fits")
        self.simple_fits_checkbox.setChecked(False)
        self.simple_fits_checkbox.setToolTip(
            "Plot stored Simple Fits results (data + model + residuals).\n"
            "Only HDF5 files with stored simple fit results are used."
        )
        cb_row1.addWidget(self.simple_fits_checkbox)
        cb_row1.addStretch()

        self.waxs_peakfit_checkbox = QCheckBox("WAXS Peaks")
        self.waxs_peakfit_checkbox.setChecked(False)
        self.waxs_peakfit_checkbox.setToolTip(
            "Plot stored WAXS peak-fit results (data + model + residuals).\n"
            "Only HDF5 files with stored WAXS peak-fit results are used."
        )
        cb_row2.addWidget(self.waxs_peakfit_checkbox)

        self.modeling_checkbox = QCheckBox("Modeling")
        self.modeling_checkbox.setChecked(False)
        self.modeling_checkbox.setToolTip(
            "Plot stored Modeling results (total model I(Q)).\n"
            "Only HDF5 files with stored Modeling results are used."
        )
        cb_row2.addWidget(self.modeling_checkbox)

        self.saxs_morph_checkbox = QCheckBox("3D saxsMorph")
        self.saxs_morph_checkbox.setChecked(False)
        self.saxs_morph_checkbox.setToolTip(
            "Open 3D saxsMorph viewer for the first selected file with\n"
            "stored saxs_morph results (2D slice + 3D voxelgram)."
        )
        cb_row2.addWidget(self.saxs_morph_checkbox)

        self.fractals_checkbox = QCheckBox("Fractals")
        self.fractals_checkbox.setChecked(False)
        self.fractals_checkbox.setToolTip(
            "Open Fractals viewer for the first selected file with stored\n"
            "fractal aggregates (2D slice + 3D voxelgram + I(Q))."
        )
        cb_row2.addWidget(self.fractals_checkbox)
        cb_row2.addStretch()

        grp_view_lay.addLayout(cb_row1)
        grp_view_lay.addLayout(cb_row2)

        # 4 output buttons — smaller, 2×2 grid
        _btn_sm_blue = (
            "QPushButton { background:#3498db; color:white; font-size:11px; "
            "font-weight:bold; border-radius:4px; padding:3px 6px; }"
            "QPushButton:hover { background:#2980b9; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        _btn_sm_purple = (
            "QPushButton { background:#8e44ad; color:white; font-size:11px; "
            "font-weight:bold; border-radius:4px; padding:3px 6px; }"
            "QPushButton:hover { background:#7d3c98; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        _btn_sm_teal = (
            "QPushButton { background:#16a085; color:white; font-size:11px; "
            "font-weight:bold; border-radius:4px; padding:3px 6px; }"
            "QPushButton:hover { background:#138d75; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        _btn_sm_orange = (
            "QPushButton { background:#2980b9; color:white; font-size:11px; "
            "font-weight:bold; border-radius:4px; padding:3px 6px; }"
            "QPushButton:hover { background:#2471a3; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )

        self.plot_button = QPushButton("Create Graph")
        self.plot_button.setMinimumHeight(23)
        self.plot_button.setStyleSheet(_btn_sm_blue)
        self.plot_button.setToolTip(
            "Plot I(Q) for all selected files in a new graph window.\n"
            "Select files in the list above, then click this button."
        )
        self.plot_button.clicked.connect(self.plot_selected_files)
        self.plot_button.setEnabled(False)

        self.report_button = QPushButton("Create Report")
        self.report_button.setMinimumHeight(23)
        self.report_button.setStyleSheet(_btn_sm_purple)
        self.report_button.setToolTip(
            "Generate a Markdown report (.md) summarising the Unified Fit\n"
            "results for each selected HDF5 file.  Files without stored\n"
            "fit results are skipped."
        )
        self.report_button.clicked.connect(self.create_report)
        self.report_button.setEnabled(False)

        self.tabulate_button = QPushButton("Tabulate Results")
        self.tabulate_button.setMinimumHeight(23)
        self.tabulate_button.setStyleSheet(_btn_sm_teal)
        self.tabulate_button.setToolTip(
            "Build a table of fit results for selected files and display it.\n"
            "Results included depend on the checked checkboxes.\n"
            "You can save the table as a CSV file (Excel, Igor Pro, Origin, etc.)."
        )
        self.tabulate_button.clicked.connect(self.tabulate_results)
        self.tabulate_button.setEnabled(False)

        self.export_ascii_button = QPushButton("Export to ASCII")
        self.export_ascii_button.setMinimumHeight(23)
        self.export_ascii_button.setStyleSheet(_btn_sm_orange)
        self.export_ascii_button.setToolTip(
            "Write space-separated .dat files into an 'ascii_export' subfolder\n"
            "next to each selected HDF5 file.  Three columns: Q  I  dI.\n"
            "USAXS files: only desmeared data is exported (slit-smeared\n"
            "variants are skipped).\n"
            "If 'Also write model curves' is enabled in Configure, an extra\n"
            "4-column file is added per enabled fit-result checkbox\n"
            "(_unif, _simp, _mod, _sd, _waxs).\n"
            "Defaults (Configure...): space delimiter, 7 sig figs,\n"
            "header on, model curves on."
        )
        self.export_ascii_button.clicked.connect(self.export_to_ascii)
        self.export_ascii_button.setEnabled(False)

        # Data Explorer — row 2 col 0 of the same grid (half-width); distinct
        # red colour signals it is a different kind of tool from the four above.
        _de_style = (
            "QPushButton { background:#c0392b; color:white; font-size:11px; "
            "font-weight:bold; border-radius:4px; padding:3px 6px; }"
            "QPushButton:hover { background:#e74c3c; }"
        )
        self.hdf5_viewer_button = QPushButton("Data Explorer")
        self.hdf5_viewer_button.setMinimumHeight(23)
        self.hdf5_viewer_button.setStyleSheet(_de_style)
        self.hdf5_viewer_button.setToolTip(
            "Open the Data Explorer to browse, plot, and export data\n"
            "from HDF5 files — including export to Igor Pro h5xp format."
        )
        self.hdf5_viewer_button.clicked.connect(self.launch_hdf5_viewer)

        _out_grid = QGridLayout()
        _out_grid.setHorizontalSpacing(4)
        _out_grid.setVerticalSpacing(4)
        _out_grid.setColumnStretch(0, 1)
        _out_grid.setColumnStretch(1, 1)
        _out_grid.addWidget(self.plot_button,          0, 0)
        _out_grid.addWidget(self.report_button,        0, 1)
        _out_grid.addWidget(self.tabulate_button,      1, 0)
        _out_grid.addWidget(self.export_ascii_button,  1, 1)
        _out_grid.addWidget(self.hdf5_viewer_button,   2, 0)
        grp_view_lay.addLayout(_out_grid)

        right_layout.addWidget(grp_view)

        right_layout.addSpacing(4)

        # ── GROUP 2: Analysis Tools ────────────────────────────────────────────
        grp_analysis = QGroupBox("Analysis Tools")
        analysis_grid = QGridLayout(grp_analysis)
        analysis_grid.setHorizontalSpacing(4)
        analysis_grid.setVerticalSpacing(4)
        analysis_grid.setColumnStretch(0, 1)
        analysis_grid.setColumnStretch(1, 1)

        _uf_gui_style = (
            "QPushButton { background:#27ae60; color:white; font-size:12px; "
            "font-weight:bold; border-radius:4px; padding:6px 8px; }"
            "QPushButton:hover { background:#229954; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        _uf_script_style = (
            "QPushButton { background:#1e8449; color:white; font-size:12px; "
            "font-weight:bold; border-radius:4px; padding:6px 8px; }"
            "QPushButton:hover { background:#196f3d; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        self.unified_fit_button = QPushButton("Unified Fit (GUI)")
        self.unified_fit_button.setMinimumHeight(23)
        self.unified_fit_button.setStyleSheet(_uf_gui_style)
        self.unified_fit_button.setToolTip("Open Unified Fit panel for the first selected file.")
        self.unified_fit_button.clicked.connect(self.launch_unified_fit)
        self.unified_fit_button.setEnabled(False)

        self.unified_script_button = QPushButton("Unified Fit (script)")
        self.unified_script_button.setMinimumHeight(23)
        self.unified_script_button.setStyleSheet(_uf_script_style)
        self.unified_script_button.setToolTip(
            "Batch-fit all selected files with Unified Fit using pyirena_config.json.\n"
            "Results are saved into each file's NXcanSAS record."
        )
        self.unified_script_button.clicked.connect(self.run_unified_script)
        self.unified_script_button.setEnabled(False)

        _sz_gui_style = (
            "QPushButton { background:#2980b9; color:white; font-size:12px; "
            "font-weight:bold; border-radius:4px; padding:6px 8px; }"
            "QPushButton:hover { background:#2471a3; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        _sz_script_style = (
            "QPushButton { background:#1f618d; color:white; font-size:12px; "
            "font-weight:bold; border-radius:4px; padding:6px 8px; }"
            "QPushButton:hover { background:#1a5276; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        self.sizes_fit_button = QPushButton("Size Distribution (GUI)")
        self.sizes_fit_button.setMinimumHeight(23)
        self.sizes_fit_button.setStyleSheet(_sz_gui_style)
        self.sizes_fit_button.setToolTip("Open Size Distribution panel for the first selected file.")
        self.sizes_fit_button.clicked.connect(self.launch_sizes_fit)
        self.sizes_fit_button.setEnabled(False)

        self.sizes_script_button = QPushButton("Size Distribution (script)")
        self.sizes_script_button.setMinimumHeight(23)
        self.sizes_script_button.setStyleSheet(_sz_script_style)
        self.sizes_script_button.setToolTip(
            "Batch-fit all selected files with Size Distribution using pyirena_config.json.\n"
            "Results are saved into each file's NXcanSAS record."
        )
        self.sizes_script_button.clicked.connect(self.run_sizes_script)
        self.sizes_script_button.setEnabled(False)

        _mod_gui_style = (
            "QPushButton { background:#e67e22; color:white; font-size:12px; "
            "font-weight:bold; border-radius:4px; padding:4px; border:none; }"
            "QPushButton:hover { background:#ca6f1e; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        _mod_script_style = (
            "QPushButton { background:#d35400; color:white; font-size:12px; "
            "font-weight:bold; border-radius:4px; padding:4px; border:none; }"
            "QPushButton:hover { background:#ba4a00; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        self.modeling_button = QPushButton("Modeling (GUI)")
        self.modeling_button.setMinimumHeight(23)
        self.modeling_button.setStyleSheet(_mod_gui_style)
        self.modeling_button.setToolTip("Open Modeling panel for the first selected file.")
        self.modeling_button.clicked.connect(self.launch_modeling)
        self.modeling_button.setEnabled(False)

        self.modeling_script_button = QPushButton("Modeling (script)")
        self.modeling_script_button.setMinimumHeight(23)
        self.modeling_script_button.setStyleSheet(_mod_script_style)
        self.modeling_script_button.setToolTip(
            "Batch-fit all selected files with Modeling using pyirena_config.json.\n"
            "Results are saved into each file's NXcanSAS record."
        )
        self.modeling_script_button.clicked.connect(self.run_modeling_script)
        self.modeling_script_button.setEnabled(False)

        _sf_gui_style = (
            "QPushButton { background:#27ae60; color:white; font-size:12px; "
            "font-weight:bold; border-radius:4px; padding:4px; border:none; }"
            "QPushButton:hover { background:#1e8449; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        _sf_script_style = (
            "QPushButton { background:#1e8449; color:white; font-size:12px; "
            "font-weight:bold; border-radius:4px; padding:4px; border:none; }"
            "QPushButton:hover { background:#196f3d; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        self.simple_fits_button = QPushButton("Simple Fits (GUI)")
        self.simple_fits_button.setMinimumHeight(23)
        self.simple_fits_button.setStyleSheet(_sf_gui_style)
        self.simple_fits_button.setToolTip("Open Simple Fits panel for the first selected file.")
        self.simple_fits_button.clicked.connect(self.launch_simple_fits)
        self.simple_fits_button.setEnabled(False)

        self.simple_fits_script_button = QPushButton("Simple Fits (script)")
        self.simple_fits_script_button.setMinimumHeight(23)
        self.simple_fits_script_button.setStyleSheet(_sf_script_style)
        self.simple_fits_script_button.setToolTip(
            "Batch-fit all selected files with Simple Fits using pyirena_config.json.\n"
            "Results are saved into each file's NXcanSAS record."
        )
        self.simple_fits_script_button.clicked.connect(self.run_simple_fits_script)
        self.simple_fits_script_button.setEnabled(False)

        _waxs_gui_style = (
            "QPushButton { background:#2980b9; color:white; font-size:12px; "
            "font-weight:bold; border-radius:4px; padding:4px; border:none; }"
            "QPushButton:hover { background:#1f618d; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        _waxs_script_style = (
            "QPushButton { background:#1f618d; color:white; font-size:12px; "
            "font-weight:bold; border-radius:4px; padding:4px; border:none; }"
            "QPushButton:hover { background:#1a5276; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        self.waxs_peakfit_button = QPushButton("WAXS Peaks (GUI)")
        self.waxs_peakfit_button.setMinimumHeight(23)
        self.waxs_peakfit_button.setStyleSheet(_waxs_gui_style)
        self.waxs_peakfit_button.setToolTip("Open WAXS Peak Fit panel for the first selected file.")
        self.waxs_peakfit_button.clicked.connect(self.launch_waxs_peakfit)
        self.waxs_peakfit_button.setEnabled(False)

        self.waxs_peakfit_script_button = QPushButton("WAXS Peaks (script)")
        self.waxs_peakfit_script_button.setMinimumHeight(23)
        self.waxs_peakfit_script_button.setStyleSheet(_waxs_script_style)
        self.waxs_peakfit_script_button.setToolTip(
            "Batch-fit all selected files with WAXS Peak Fit using pyirena_config.json.\n"
            "Results are saved into each file's NXcanSAS record."
        )
        self.waxs_peakfit_script_button.clicked.connect(self.run_waxs_peakfit_script)
        self.waxs_peakfit_script_button.setEnabled(False)

        _sm_gui_style = (
            "QPushButton { background:#8e44ad; color:white; font-size:12px; "
            "font-weight:bold; border-radius:4px; padding:4px; border:none; }"
            "QPushButton:hover { background:#7d3c98; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        _sm_script_style = (
            "QPushButton { background:#6c3483; color:white; font-size:12px; "
            "font-weight:bold; border-radius:4px; padding:4px; border:none; }"
            "QPushButton:hover { background:#5b2c6f; }"
            "QPushButton:disabled { background:#bdc3c7; }"
        )
        self.saxs_morph_button = QPushButton("3D saxsMorph (GUI)")
        self.saxs_morph_button.setMinimumHeight(23)
        self.saxs_morph_button.setStyleSheet(_sm_gui_style)
        self.saxs_morph_button.setToolTip(
            "Open SAXS Morph (3D voxelgram) panel for the first selected file."
        )
        self.saxs_morph_button.clicked.connect(self.launch_saxs_morph)
        self.saxs_morph_button.setEnabled(False)

        self.saxs_morph_script_button = QPushButton("3D saxsMorph (script)")
        self.saxs_morph_script_button.setMinimumHeight(23)
        self.saxs_morph_script_button.setStyleSheet(_sm_script_style)
        self.saxs_morph_script_button.setToolTip(
            "Batch-fit all selected files with SAXS Morph using pyirena_config.json.\n"
            "Compressed voxelgrams + parameters are saved into each NXcanSAS file."
        )
        self.saxs_morph_script_button.clicked.connect(self.run_saxs_morph_script)
        self.saxs_morph_script_button.setEnabled(False)

        analysis_grid.addWidget(self.unified_fit_button,         0, 0)
        analysis_grid.addWidget(self.unified_script_button,      0, 1)
        analysis_grid.addWidget(self.modeling_button,            1, 0)
        analysis_grid.addWidget(self.modeling_script_button,     1, 1)
        analysis_grid.addWidget(self.sizes_fit_button,           2, 0)
        analysis_grid.addWidget(self.sizes_script_button,        2, 1)
        analysis_grid.addWidget(self.simple_fits_button,         3, 0)
        analysis_grid.addWidget(self.simple_fits_script_button,  3, 1)
        analysis_grid.addWidget(self.waxs_peakfit_button,        4, 0)
        analysis_grid.addWidget(self.waxs_peakfit_script_button, 4, 1)
        analysis_grid.addWidget(self.saxs_morph_button,          5, 0)
        analysis_grid.addWidget(self.saxs_morph_script_button,   5, 1)
        right_layout.addWidget(grp_analysis)

        # ── GROUP 3: Data Processing & Reference ──────────────────────────────
        grp_proc = QGroupBox("Data Processing & Reference")
        proc_grid = QGridLayout(grp_proc)
        proc_grid.setHorizontalSpacing(4)
        proc_grid.setVerticalSpacing(4)
        proc_grid.setColumnStretch(0, 1)
        proc_grid.setColumnStretch(1, 1)

        _utility_style = (
            "QPushButton { background:#16a085; color:white; "
            "font-weight:bold; border-radius:4px; padding:4px 8px; }"
            "QPushButton:hover { background:#1abc9c; }"
            "QPushButton:disabled { background:#95a5a6; }"
        )
        self.data_merge_button = QPushButton("Data Merge")
        self.data_merge_button.setMinimumHeight(23)
        self.data_merge_button.setStyleSheet(_utility_style)
        self.data_merge_button.setToolTip(
            "Open the Data Merge tool to combine USAXS and SAXS/WAXS datasets."
        )
        self.data_merge_button.clicked.connect(self.launch_data_merge)

        self.data_manip_button = QPushButton("Data Manipulation")
        self.data_manip_button.setMinimumHeight(23)
        self.data_manip_button.setStyleSheet(_utility_style)
        self.data_manip_button.setToolTip(
            "Open the Data Manipulation tool for scaling, trimming,\n"
            "rebinning, averaging, subtracting, and dividing datasets."
        )
        self.data_manip_button.clicked.connect(self.launch_data_manipulation)

        _ref_style = (
            "QPushButton { background:#16a085; color:white; "
            "font-weight:bold; border-radius:4px; padding:4px 8px; }"
            "QPushButton:hover { background:#138d75; }"
            "QPushButton:disabled { background:#95a5a6; }"
        )
        self.contrast_button = QPushButton("Scattering Contrast Calculator")
        self.contrast_button.setMinimumHeight(23)
        self.contrast_button.setStyleSheet(_ref_style)
        self.contrast_button.setToolTip(
            "Open the Scattering Contrast Calculator.\n"
            "Computes X-ray and neutron SLDs and contrast for two compounds\n"
            "(free-electron and anomalous Chantler-corrected X-ray values)."
        )
        self.contrast_button.clicked.connect(self.launch_contrast)

        self.fractals_button = QPushButton("Fractals")
        self.fractals_button.setMinimumHeight(23)
        self.fractals_button.setStyleSheet(_ref_style)
        self.fractals_button.setToolTip(
            "Open the Fractals tool: grow random mass-fractal aggregates by\n"
            "Monte-Carlo on a simple cubic lattice, compute their fractal\n"
            "parameters (Z, dmin, c, df, Rg primary, Rg aggregate), and\n"
            "back-calculate I(Q).  Optionally compares against Unified-fit\n"
            "results from a NeXus file."
        )
        self.fractals_button.clicked.connect(self.launch_fractals)

        _igor_import_style = (
            "QPushButton { background:#8e44ad; color:white; "
            "font-weight:bold; border-radius:4px; padding:4px 8px; }"
            "QPushButton:hover { background:#6c3483; }"
            "QPushButton:disabled { background:#95a5a6; }"
        )
        self.igor_import_button = QPushButton("Import Igor Experiment…")
        self.igor_import_button.setMinimumHeight(23)
        self.igor_import_button.setStyleSheet(_igor_import_style)
        self.igor_import_button.setToolTip(
            "Open an Igor Pro packed experiment (.pxp or .h5xp) and export\n"
            "each reduced USAXS/SAXS/WAXS sample into a NeXus (.h5) file.\n"
            "Use this to bring legacy Igor data into pyIrena for analysis."
        )
        self.igor_import_button.clicked.connect(self.launch_igor_import)

        proc_grid.addWidget(self.data_merge_button, 0, 0)
        proc_grid.addWidget(self.data_manip_button, 0, 1)
        proc_grid.addWidget(self.contrast_button,   1, 0)
        proc_grid.addWidget(self.fractals_button,   1, 1)
        proc_grid.addWidget(self.igor_import_button, 2, 0, 1, 2)  # span both cols
        right_layout.addWidget(grp_proc)

        right_layout.addStretch()

        # Wrap the right column in a scroll area so users on small or DPI-scaled
        # displays (Windows 125–150% scaling) see a vertical scrollbar instead
        # of squashed buttons. On large displays the scrollbar never appears
        # and the layout looks identical to before.
        right_container = QWidget()
        right_container.setLayout(right_layout)
        right_scroll = QScrollArea()
        right_scroll.setWidget(right_container)
        right_scroll.setWidgetResizable(True)
        right_scroll.setFrameShape(QFrame.Shape.NoFrame)
        right_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        right_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        # Match the original right-column natural width so the two-column
        # button grid never gets horizontally compressed (and no horizontal
        # scrollbar ever appears).
        right_scroll.setMinimumWidth(450)
        file_area_layout.addWidget(right_scroll, stretch=1)

        content_layout.addLayout(file_area_layout)

        # Status bar
        self.status_label = QLabel("Ready - Select a folder to begin")
        self.status_label.setStyleSheet("""
            QLabel {
                color: #7f8c8d;
                padding: 5px;
                border-top: 1px solid #bdc3c7;
            }
        """)
        content_layout.addWidget(self.status_label)

        # Add content layout to main layout
        main_layout.addLayout(content_layout)

        self.setLayout(main_layout)

        # Set minimum window size. Lowered from 600 to 500 now that the right
        # column is in a QScrollArea — small screens / scaled displays no longer
        # squash the buttons. File list still has setMinimumHeight(400) which
        # dominates anyway.
        self.setMinimumSize(900, 500)

        # If we have a restored folder, display it and list files
        if self.current_folder:
            self.folder_label.setText(self.current_folder)
            self.folder_label.setStyleSheet("color: #2c3e50;")
            self.refresh_button.setEnabled(True)
            self.refresh_file_list()
            self.status_label.setText(f"Restored folder: {os.path.basename(self.current_folder)}")

    def create_menu_bar(self) -> QMenuBar:
        """Create the menu bar with Models menu."""
        menu_bar = QMenuBar()

        # Models menu
        models_menu = QMenu("&Models", self)

        # Unified Fit action
        unified_fit_action = QAction("&Unified Fit", self)
        unified_fit_action.setStatusTip("Open Unified Fit model panel")
        unified_fit_action.triggered.connect(self.launch_unified_fit)
        models_menu.addAction(unified_fit_action)

        # Size Distribution action
        sizes_fit_action = QAction("&Size Distribution", self)
        sizes_fit_action.setStatusTip("Open Size Distribution fitting panel")
        sizes_fit_action.triggered.connect(self.launch_sizes_fit)
        models_menu.addAction(sizes_fit_action)

        # Simple Fits action
        simple_fits_action = QAction("Simple &Fits", self)
        simple_fits_action.setStatusTip("Open Simple Fits panel")
        simple_fits_action.triggered.connect(self.launch_simple_fits)
        models_menu.addAction(simple_fits_action)

        # Modeling action
        modeling_action = QAction("&Modeling", self)
        modeling_action.setStatusTip("Open Modeling panel")
        modeling_action.triggered.connect(self.launch_modeling)
        models_menu.addAction(modeling_action)

        # WAXS Peak Fit action
        waxs_peakfit_action = QAction("&WAXS Peak Fit", self)
        waxs_peakfit_action.setStatusTip("Open WAXS Peak Fit panel")
        waxs_peakfit_action.triggered.connect(self.launch_waxs_peakfit)
        models_menu.addAction(waxs_peakfit_action)

        # SAXS Morph action (3D voxelgram)
        saxs_morph_action = QAction("S&AXS Morph (3D)", self)
        saxs_morph_action.setStatusTip("Open SAXS Morph (3D voxelgram) panel")
        saxs_morph_action.triggered.connect(self.launch_saxs_morph)
        models_menu.addAction(saxs_morph_action)

        menu_bar.addMenu(models_menu)

        # Tools menu
        tools_menu = QMenu("&Tools", self)
        data_merge_action = QAction("&Data Merge", self)
        data_merge_action.setStatusTip("Open Data Merge tool for USAXS+SAXS/WAXS merging")
        data_merge_action.triggered.connect(self.launch_data_merge)
        tools_menu.addAction(data_merge_action)
        data_manip_action = QAction("Data &Manipulation", self)
        data_manip_action.setStatusTip("Open Data Manipulation tool for scaling, trimming, averaging, etc.")
        data_manip_action.triggered.connect(self.launch_data_manipulation)
        tools_menu.addAction(data_manip_action)
        tools_menu.addSeparator()
        hdf5_viewer_action = QAction("&Data Explorer", self)
        hdf5_viewer_action.setStatusTip("Open Data Explorer (browse, plot, and export HDF5 data)")
        hdf5_viewer_action.triggered.connect(self.launch_hdf5_viewer)
        tools_menu.addAction(hdf5_viewer_action)
        tools_menu.addSeparator()
        contrast_action = QAction("&Scattering Contrast", self)
        contrast_action.setStatusTip("Open Scattering Contrast Calculator")
        contrast_action.triggered.connect(self.launch_contrast)
        tools_menu.addAction(contrast_action)
        menu_bar.addMenu(tools_menu)

        # Help menu
        help_menu = QMenu("&Help", self)

        about_action = QAction("&About pyIrena", self)
        about_action.setStatusTip("About pyIrena")
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

        menu_bar.addMenu(help_menu)

        return menu_bar

    def select_folder(self):
        """Open folder selection dialog."""
        # Use last folder if available, otherwise home directory
        start_dir = self.last_folder if self.last_folder else QDir.homePath()

        folder = QFileDialog.getExistingDirectory(
            self,
            "Select Data Folder",
            start_dir,
            QFileDialog.Option.ShowDirsOnly
        )

        if folder:
            self.current_folder = folder
            self.last_folder = folder  # Remember for next time
            self.save_last_folder(folder)  # Save to state
            self.folder_label.setText(folder)
            self.folder_label.setStyleSheet("color: #2c3e50;")
            self.refresh_button.setEnabled(True)
            self.refresh_file_list()
            self.status_label.setText(f"Folder: {os.path.basename(folder)}")

    def get_file_extensions(self) -> List[str]:
        """Get current file extensions based on selection."""
        file_type = self.file_type_combo.currentData()

        if file_type == "hdf5":
            return self.HDF5_EXTENSIONS
        elif file_type == "text":
            return self.TEXT_EXTENSIONS
        else:  # all
            return self.HDF5_EXTENSIONS + self.TEXT_EXTENSIONS

    def refresh_file_list(self):
        """Refresh the file list based on current folder and file type."""
        self.file_list.clear()

        if not self.current_folder:
            return

        extensions = self.get_file_extensions()

        try:
            files = []
            for file in os.listdir(self.current_folder):
                file_path = os.path.join(self.current_folder, file)
                if os.path.isfile(file_path):
                    _, ext = os.path.splitext(file)
                    if ext.lower() in extensions:
                        files.append(file)

            self.file_list.addItems(files)
            self.status_label.setText(f"Found {len(files)} files")
            self.sort_file_list()   # apply current sort-combo order
            self.update_plot_button_state()

        except Exception as e:
            self.status_label.setText(f"Error reading folder: {e}")

    def filter_files(self, filter_text: str):
        """Filter the file list based on search text."""
        if not filter_text:
            # Show all items
            for i in range(self.file_list.count()):
                self.file_list.item(i).setHidden(False)
            return

        # Use grep-like filtering (case-insensitive substring match)
        filter_text = filter_text.lower()
        visible_count = 0

        for i in range(self.file_list.count()):
            item = self.file_list.item(i)
            matches = filter_text in item.text().lower()
            item.setHidden(not matches)
            if matches:
                visible_count += 1

        self.status_label.setText(f"Showing {visible_count} of {self.file_list.count()} files")

    def _on_sort_changed(self, index: int):
        """Save sort selection to state, then apply it to the file list."""
        self.state_manager.set('data_selector', 'sort_index', index)
        self.state_manager.save()
        self.sort_file_list()

    def sort_file_list(self):
        """
        Re-order all items in the file list according to the sort combo selection.
        The current text filter is re-applied afterwards so hidden items stay hidden.
        """
        n = self.file_list.count()
        if n == 0:
            return

        items = [self.file_list.item(i).text() for i in range(n)]

        idx     = self.sort_combo.currentIndex()
        reverse = bool(idx % 2)                 # odd indices → descending
        key_fn  = _SORT_KEYS[min(idx, len(_SORT_KEYS) - 1)]

        items.sort(key=key_fn, reverse=reverse)

        self.file_list.clear()
        self.file_list.addItems(items)

        # Re-apply whatever text filter is active
        self.filter_files(self.filter_input.text())

    def update_plot_button_state(self):
        """Enable or disable buttons based on file selection."""
        has_selection = len(self.file_list.selectedItems()) > 0
        self.plot_button.setEnabled(has_selection)
        self.report_button.setEnabled(has_selection)
        self.tabulate_button.setEnabled(has_selection)
        self.export_ascii_button.setEnabled(has_selection)
        self.unified_fit_button.setEnabled(has_selection)
        self.unified_script_button.setEnabled(has_selection)
        self.modeling_button.setEnabled(has_selection)
        self.modeling_script_button.setEnabled(has_selection)
        self.sizes_fit_button.setEnabled(has_selection)
        self.sizes_script_button.setEnabled(has_selection)
        self.simple_fits_button.setEnabled(has_selection)
        self.simple_fits_script_button.setEnabled(has_selection)
        self.waxs_peakfit_button.setEnabled(has_selection)
        self.waxs_peakfit_script_button.setEnabled(has_selection)
        self.saxs_morph_button.setEnabled(has_selection)
        self.saxs_morph_script_button.setEnabled(has_selection)

    def plot_selected_files(self):
        """Plot the selected files according to the Data / Unified Fit checkboxes."""
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(
                self,
                "No Selection",
                "Please select one or more files to plot."
            )
            return

        file_paths = [
            os.path.join(self.current_folder, item.text())
            for item in selected_items
        ]

        error_fraction    = self.state_manager.get('data_selector', 'error_fraction', 0.05)
        max_legend_items  = int(self.state_manager.get('data_selector', 'max_legend_items', 12))
        plotted = []

        # ── Experimental data ──────────────────────────────────────────────
        if self.data_checkbox.isChecked():
            if self.graph_window is None:
                self.graph_window = GraphWindow()
            try:
                self.graph_window.plot_data(
                    file_paths,
                    error_fraction=error_fraction,
                    max_legend_items=max_legend_items,
                )
                plotted.append("data")
            except Exception as e:
                QMessageBox.critical(
                    self, "Plot Error", f"Error creating data plot:\n{str(e)}"
                )

        # ── Unified Fit results ────────────────────────────────────────────
        if self.unified_fit_result_checkbox.isChecked():
            if self.unified_fit_results_window is None:
                self.unified_fit_results_window = UnifiedFitResultsWindow()
            try:
                self.unified_fit_results_window.plot_results(
                    file_paths,
                    max_legend_items=max_legend_items,
                )
                plotted.append("Unified Fit results")
            except Exception as e:
                QMessageBox.critical(
                    self, "Plot Error", f"Error creating Unified Fit plot:\n{str(e)}"
                )

        # ── Size Distribution results ───────────────────────────────────────
        if self.size_dist_checkbox.isChecked():
            if self.size_dist_results_window is None:
                self.size_dist_results_window = SizeDistResultsWindow()
            try:
                self.size_dist_results_window.plot_results(
                    file_paths,
                    max_legend_items=max_legend_items,
                )
                plotted.append("Size Distribution results")
            except Exception as e:
                QMessageBox.critical(
                    self, "Plot Error", f"Error creating Size Distribution plot:\n{str(e)}"
                )

        # ── Simple Fits results ────────────────────────────────────────────
        if self.simple_fits_checkbox.isChecked():
            if self.simple_fits_results_window is None:
                self.simple_fits_results_window = SimpleFitResultsWindow()
            try:
                self.simple_fits_results_window.plot_results(
                    file_paths,
                    max_legend_items=max_legend_items,
                )
                plotted.append("Simple Fits results")
            except Exception as e:
                QMessageBox.critical(
                    self, "Plot Error", f"Error creating Simple Fits plot:\n{str(e)}"
                )

        # ── WAXS Peak Fit results ──────────────────────────────────────────
        if self.waxs_peakfit_checkbox.isChecked():
            if self.waxs_peakfit_results_window is None:
                self.waxs_peakfit_results_window = WAXSPeakFitResultsWindow()
            try:
                self.waxs_peakfit_results_window.plot_results(
                    file_paths,
                    max_legend_items=max_legend_items,
                )
                plotted.append("WAXS Peak Fit results")
            except Exception as e:
                QMessageBox.critical(
                    self, "Plot Error", f"Error creating WAXS Peak Fit plot:\n{str(e)}"
                )

        # ── Modeling results ───────────────────────────────────────────────
        if self.modeling_checkbox.isChecked():
            from pyirena.gui.hdf5viewer.pyirena_readers import read_modeling
            if self.graph_window is None:
                self.graph_window = GraphWindow()
            colors = _gen_colors(len(file_paths))
            mod_plotted = 0
            legend_idx = _legend_indices(len(file_paths), max_legend_items)
            for i, fp in enumerate(file_paths):
                mod = read_modeling(fp)
                if mod is not None and len(mod["Q"]) > 0:
                    q = np.asarray(mod["Q"], dtype=float)
                    I = np.asarray(mod["I_model"], dtype=float)
                    mask = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
                    if not mask.any():
                        continue
                    color = colors[i]
                    name = f'{Path(fp).stem} model' if i in legend_idx else None
                    self.graph_window.plot.plot(
                        q[mask], I[mask],
                        pen=pg.mkPen(color, width=2, style=Qt.PenStyle.DashLine),
                        name=name,
                    )
                    mod_plotted += 1
            if mod_plotted:
                self.graph_window.show()
                plotted.append(f"Modeling ({mod_plotted} file(s))")

        # ── 3D saxsMorph viewer ─────────────────────────────────────────────
        if self.saxs_morph_checkbox.isChecked():
            opened = self._open_saxs_morph_viewer(file_paths)
            if opened:
                plotted.append(f"3D saxsMorph ({opened} file(s))")

        # ── Fractals viewer ─────────────────────────────────────────────────
        if self.fractals_checkbox.isChecked():
            opened = self._open_fractals_viewer(file_paths)
            if opened:
                plotted.append(f"Fractals ({opened} file(s))")

        if plotted:
            self.status_label.setText(
                f"Plotted {len(file_paths)} file(s): {', '.join(plotted)}"
            )
        else:
            self.status_label.setText(
                "Nothing to plot — check 'Data', 'Unified Fit', 'Size Dist.', "
                "'Simple Fits', 'WAXS Peaks' or 'Modeling' checkbox"
            )

    def create_report(self):
        """
        Generate a Markdown report for each selected file.

        Content depends on which checkboxes are ticked:
          - 'Data'        → basic data summary (Q range, I range, N points)
          - 'Unified Fit' → fit quality metrics and level parameters (HDF5 only)
          - 'Size Dist.'  → size distribution results (HDF5 only)
        Multiple checkboxes can be active simultaneously.
        """
        selected_items = self.file_list.selectedItems()
        if not selected_items:
            QMessageBox.warning(self, "No Selection", "Please select files to report.")
            return

        show_data         = self.data_checkbox.isChecked()
        show_fit          = self.unified_fit_result_checkbox.isChecked()
        show_sizes        = self.size_dist_checkbox.isChecked()
        show_simple_fits  = self.simple_fits_checkbox.isChecked()
        show_waxs_peakfit = self.waxs_peakfit_checkbox.isChecked()
        show_modeling     = self.modeling_checkbox.isChecked()
        show_saxs_morph   = self.saxs_morph_checkbox.isChecked()

        if not any([show_data, show_fit, show_sizes, show_simple_fits,
                    show_waxs_peakfit, show_modeling, show_saxs_morph]):
            self.status_label.setText(
                "Nothing to report — check 'Data', 'Unified Fit', 'Size Dist.', "
                "'Simple Fits', 'WAXS Peaks', 'Modeling' or '3D saxsMorph' checkbox"
            )
            return

        file_paths = [
            os.path.join(self.current_folder, item.text())
            for item in selected_items
        ]

        error_fraction = self.state_manager.get('data_selector', 'error_fraction', 0.05)
        saved, skipped = [], []

        for file_path in file_paths:
            _, ext = os.path.splitext(file_path)
            is_hdf = ext.lower() in ['.h5', '.hdf5', '.hdf']

            data_info            = None
            fit_results          = None
            sizes_results        = None
            simple_fit_results   = None
            waxs_peakfit_results = None
            modeling_results     = None
            saxs_morph_results   = None

            # ── Load raw data ──────────────────────────────────────────────
            if show_data:
                try:
                    dir_path, filename = os.path.split(file_path)
                    if ext.lower() in ['.txt', '.dat']:
                        h5_path = ensure_nxcansas_sibling(
                            Path(file_path), error_fraction=error_fraction)
                        raw = readGenericNXcanSAS(str(h5_path.parent), h5_path.name)
                    else:
                        raw = readGenericNXcanSAS(dir_path, filename)

                    if raw is not None:
                        data_info = {
                            'Q':       raw['Q'],
                            'I':       raw['Intensity'],
                            'I_error': raw.get('Error'),
                        }
                except Exception:
                    log.debug("suppressed exception", exc_info=True)   # data load failure is non-fatal; section simply omitted

            # ── Load Unified Fit results (HDF5 only) ───────────────────────
            if show_fit and is_hdf:
                try:
                    fit_results = load_unified_fit_results(Path(file_path))
                except Exception:
                    log.debug("suppressed exception", exc_info=True)   # no fit group or unreadable — section omitted silently

            # ── Load Size Distribution results (HDF5 only) ─────────────────
            if show_sizes and is_hdf:
                try:
                    from pyirena.io.nxcansas_sizes import load_sizes_results
                    sizes_results = load_sizes_results(Path(file_path))
                except Exception:
                    log.debug("suppressed exception", exc_info=True)   # no sizes group or unreadable — section omitted silently

            # ── Load Simple Fits results (HDF5 only) ───────────────────────
            if show_simple_fits and is_hdf:
                try:
                    from pyirena.io.nxcansas_simple_fits import load_simple_fit_results
                    simple_fit_results = load_simple_fit_results(Path(file_path))
                except Exception:
                    log.debug("suppressed exception", exc_info=True)   # no simple_fit_results group — section omitted silently

            # ── Load WAXS Peak Fit results (HDF5 only) ─────────────────────
            if show_waxs_peakfit and is_hdf:
                try:
                    from pyirena.io.nxcansas_waxs_peakfit import load_waxs_peakfit_results
                    waxs_peakfit_results = load_waxs_peakfit_results(Path(file_path))
                except Exception:
                    log.debug("suppressed exception", exc_info=True)   # no waxs_peakfit_results group — section omitted silently

            # ── Load Modeling results (HDF5 only) ──────────────────────────
            if show_modeling and is_hdf:
                try:
                    from pyirena.io.nxcansas_modeling import load_modeling_results
                    modeling_results = load_modeling_results(Path(file_path))
                except Exception:
                    log.debug("suppressed exception", exc_info=True)   # no modeling_results group — section omitted silently

            # ── Load SAXS Morph results (HDF5 only) ────────────────────────
            if show_saxs_morph and is_hdf:
                try:
                    from pyirena.io.nxcansas_saxs_morph import load_saxs_morph_results
                    saxs_morph_results = load_saxs_morph_results(Path(file_path))
                except Exception:
                    log.debug("suppressed exception", exc_info=True)   # no saxs_morph_results group — section omitted silently

            # Nothing to write for this file?
            if (data_info is None and fit_results is None and sizes_results is None
                    and simple_fit_results is None and waxs_peakfit_results is None
                    and modeling_results is None and saxs_morph_results is None):
                skipped.append(os.path.basename(file_path))
                continue

            md = _build_report(
                file_path,
                data_info=data_info,
                fit_results=fit_results,
                sizes_results=sizes_results,
                simple_fit_results=simple_fit_results,
                waxs_peakfit_results=waxs_peakfit_results,
                modeling_results=modeling_results,
                saxs_morph_results=saxs_morph_results,
            )
            out_path = Path(file_path).parent / (Path(file_path).stem + '_report.md')
            out_path.write_text(md, encoding='utf-8')
            saved.append(out_path.name)

            # Open in system default application for .md files
            try:
                if sys.platform == 'darwin':
                    subprocess.run(['open', str(out_path)], check=False)
                elif sys.platform == 'win32':
                    os.startfile(str(out_path))
                else:
                    subprocess.run(['xdg-open', str(out_path)], check=False)
            except Exception:
                log.debug("suppressed exception", exc_info=True)   # Opening is best-effort; don't fail the whole save

        if saved:
            msg = f"Report(s) saved: {', '.join(saved)}"
            if skipped:
                msg += f"  (skipped: {', '.join(skipped)})"
            self.status_label.setText(msg)
        else:
            self.status_label.setText(
                "No reportable content found — no reports generated"
            )

    # ─────────────────────────────────────────────────────────────────────────
    def tabulate_results(self):
        """
        Collect fit results from selected HDF5 files and display them in a
        spreadsheet-like table.  The active checkboxes control which result
        sets are included:
          • 'Unified Fit' checkbox → Unified Fit scalars + per-level parameters
          • 'Size Dist.'  checkbox → Size Distribution scalars

        The table can be saved as a CSV file from within the result window.
        """
        selected_items = self.file_list.selectedItems()
        if not selected_items:
            QMessageBox.warning(self, "No Selection", "Please select files to tabulate.")
            return

        show_unified      = self.unified_fit_result_checkbox.isChecked()
        show_sizes        = self.size_dist_checkbox.isChecked()
        show_simple_fits  = self.simple_fits_checkbox.isChecked()
        show_waxs_peakfit = self.waxs_peakfit_checkbox.isChecked()
        show_modeling     = self.modeling_checkbox.isChecked()
        show_saxs_morph   = self.saxs_morph_checkbox.isChecked()

        if not any([show_unified, show_sizes, show_simple_fits,
                    show_waxs_peakfit, show_modeling, show_saxs_morph]):
            self.status_label.setText(
                "Nothing to tabulate — check 'Unified Fit', 'Size Dist.', 'Simple Fits', "
                "'WAXS Peaks', 'Modeling' or '3D saxsMorph' checkbox"
            )
            return

        selected_set = {item.text() for item in selected_items}
        file_paths = [
            os.path.join(self.current_folder, self.file_list.item(i).text())
            for i in range(self.file_list.count())
            if self.file_list.item(i).text() in selected_set
        ]

        from pyirena.io.nxcansas_unified import load_unified_fit_results
        from pyirena.io.nxcansas_sizes import load_sizes_results

        # ── First pass: load all data, determine max unified-fit levels ──────
        loaded = []   # list of (fname, uf, sz, sf, wp, mod, sm)
        max_levels   = 0
        max_wp_peaks = 0
        max_mod_pops = 0
        # Collect all SF param names seen across files (order preserved via dict)
        sf_param_names_seen: dict = {}
        sf_derived_names_seen: dict = {}

        for fp in file_paths:
            _, ext = os.path.splitext(fp)
            is_hdf = ext.lower() in ('.h5', '.hdf5', '.hdf')
            fname  = os.path.basename(fp)

            uf  = None
            sz  = None
            sf  = None
            wp  = None
            mod = None
            sm  = None

            if is_hdf:
                if show_unified:
                    try:
                        uf = load_unified_fit_results(Path(fp))
                        max_levels = max(max_levels, len(uf.get('levels', [])))
                    except Exception:
                        log.debug("suppressed exception", exc_info=True)
                if show_sizes:
                    try:
                        sz = load_sizes_results(Path(fp))
                    except Exception:
                        log.debug("suppressed exception", exc_info=True)
                if show_simple_fits:
                    try:
                        from pyirena.io.nxcansas_simple_fits import load_simple_fit_results
                        sf = load_simple_fit_results(Path(fp))
                        for pname in sf.get('params', {}):
                            sf_param_names_seen[pname] = None
                        for dname in sf.get('derived', {}):
                            sf_derived_names_seen[dname] = None
                    except Exception:
                        log.debug("suppressed exception", exc_info=True)
                if show_waxs_peakfit:
                    try:
                        from pyirena.io.nxcansas_waxs_peakfit import load_waxs_peakfit_results
                        wp = load_waxs_peakfit_results(Path(fp))
                        max_wp_peaks = max(max_wp_peaks, int(wp.get('n_peaks', 0)))
                    except Exception:
                        log.debug("suppressed exception", exc_info=True)
                if show_modeling:
                    try:
                        from pyirena.io.nxcansas_modeling import load_modeling_results
                        mod = load_modeling_results(Path(fp))
                        n_enabled = sum(
                            1 for p in mod.get('populations', []) if p.get('enabled', True)
                        )
                        max_mod_pops = max(max_mod_pops, n_enabled)
                    except Exception:
                        log.debug("suppressed exception", exc_info=True)
                if show_saxs_morph:
                    try:
                        from pyirena.io.nxcansas_saxs_morph import load_saxs_morph_results
                        sm = load_saxs_morph_results(Path(fp))
                    except Exception:
                        log.debug("suppressed exception", exc_info=True)

            loaded.append((fname, uf, sz, sf, wp, mod, sm))

        # ── Build column headers ──────────────────────────────────────────────
        headers = ['filename']

        if show_unified:
            headers += ['UF_chi2', 'UF_background', 'UF_background_err', 'UF_num_levels']
            _uf_level_params = ['G', 'G_err', 'Rg', 'Rg_err',
                                 'B', 'B_err', 'P', 'P_err',
                                 'RgCutoff', 'ETA', 'ETA_err',
                                 'PACK', 'PACK_err', 'correlated',
                                 'Sv', 'Invariant']
            for lvl in range(1, max_levels + 1):
                for p in _uf_level_params:
                    headers.append(f'L{lvl}_{p}')

        if show_sizes:
            headers += [
                'SD_method', 'SD_shape', 'SD_chi2', 'SD_volume_fraction',
                'SD_rg', 'SD_n_iterations', 'SD_q_power',
                'SD_contrast', 'SD_aspect_ratio',
                'SD_background', 'SD_error_scale',
                'SD_power_law_B', 'SD_power_law_P',
                'SD_r_min', 'SD_r_max', 'SD_n_bins', 'SD_log_spacing',
                'SD_cursor_q_min', 'SD_cursor_q_max',
            ]

        _sf_param_names   = list(sf_param_names_seen.keys())
        _sf_derived_names = list(sf_derived_names_seen.keys())
        if show_simple_fits:
            headers += ['SF_model', 'SF_chi2', 'SF_reduced_chi2', 'SF_dof',
                        'SF_q_min', 'SF_q_max', 'SF_use_complex_bg']
            for pname in _sf_param_names:
                headers.append(f'SF_{pname}')
                headers.append(f'SF_{pname}_std')
            for dname in _sf_derived_names:
                headers.append(f'SF_derived_{dname}')

        if show_waxs_peakfit:
            headers += ['WP_n_peaks', 'WP_bg_shape', 'WP_chi2', 'WP_reduced_chi2',
                        'WP_dof', 'WP_q_min', 'WP_q_max']
            for pk in range(1, max_wp_peaks + 1):
                headers += [
                    f'WP_peak{pk}_shape',
                    f'WP_peak{pk}_Q0', f'WP_peak{pk}_Q0_std',
                    f'WP_peak{pk}_A',  f'WP_peak{pk}_A_std',
                    f'WP_peak{pk}_FWHM', f'WP_peak{pk}_FWHM_std',
                    f'WP_peak{pk}_eta', f'WP_peak{pk}_eta_std',
                    f'WP_peak{pk}_area', f'WP_peak{pk}_area_std',
                ]

        # Distribution / form-factor / structure-factor parameter sets vary by
        # model (sphere vs core-shell vs core-shell-shell vs hard-sphere SF, …),
        # so collect the union of keys actually present across all loaded
        # modeling results and emit a column for each. This keeps the export
        # complete as new form/structure factors are added, instead of relying
        # on a hand-maintained column list that silently drops new parameters.
        _mod_dist_keys, _mod_ff_keys, _mod_sf_keys = [], [], []
        if show_modeling:
            _dk, _fk, _sk = set(), set(), set()
            for _t in loaded:
                _mod_res = _t[5]
                if not _mod_res:
                    continue
                for _p in _mod_res.get('populations', []):
                    if not _p.get('enabled', True):
                        continue
                    if _p.get('pop_type', 'size_dist') != 'size_dist':
                        continue
                    _dk.update((_p.get('dist_params') or {}).keys())
                    _fk.update((_p.get('ff_params') or {}).keys())
                    _sk.update((_p.get('sf_params') or {}).keys())
            _mod_dist_keys = sorted(_dk)
            _mod_ff_keys = sorted(_fk)
            _mod_sf_keys = sorted(_sk)

        _mod_pop_cols = (
            ['type', 'label',
             # size_dist
             'dist_type', 'scale', 'contrast', 'form_factor']
            + [f'dist_{k}' for k in _mod_dist_keys]
            + [f'ff_{k}' for k in _mod_ff_keys]
            + ['structure_factor']
            + [f'sf_{k}' for k in _mod_sf_keys]
            + ['vol_fraction', 'mean_r', 'r_total_mean',
               # unified_level
               'G', 'Rg', 'B', 'P', 'RgCO', 'ETA', 'PACK',
               # diffraction_peak
               'peak_type', 'position', 'amplitude', 'width', 'eta_voigt']
        )
        if show_modeling:
            headers += ['MOD_chi2', 'MOD_background', 'MOD_q_min', 'MOD_q_max', 'MOD_n_pops']
            for k in range(1, max_mod_pops + 1):
                for col in _mod_pop_cols:
                    headers.append(f'MOD_P{k}_{col}')

        if show_saxs_morph:
            # Fit-related scalars + Porod-derived S/V + topology metrics
            # of the minority phase (open/closed porosity, percolation,
            # Euler number, pore-size percentiles).
            headers += [
                'SM_chi2', 'SM_phi', 'SM_contrast',
                'SM_rg_A', 'SM_S_per_V_inv_A',
                'SM_voxel_size', 'SM_box_size_A', 'SM_pitch_A',
                'SM_minority_phi',
                'SM_n_clusters',
                'SM_open_porosity_pct', 'SM_closed_porosity_pct',
                'SM_perc_x', 'SM_perc_y', 'SM_perc_z',
                'SM_euler',
                'SM_pore_q25_A', 'SM_pore_median_A', 'SM_pore_q75_A',
            ]

        # ── Build rows ────────────────────────────────────────────────────────
        def _fmt(v):
            """Format a value for CSV/table: round floats to 6 sig-figs."""
            if v is None:
                return None
            if isinstance(v, float):
                import math
                if math.isnan(v) or math.isinf(v):
                    return None
                return float(f'{v:.6g}')
            return v

        _wp_n_fixed_cols = 7  # WP_n_peaks … WP_q_max
        _wp_peak_cols    = 11  # shape, Q0, Q0_std, A, A_std, FWHM, FWHM_std,
                               # eta, eta_std, area, area_std

        rows = []
        for fname, uf, sz, sf, wp, mod, sm in loaded:
            row = [fname]

            if show_unified:
                if uf is not None:
                    row += [
                        _fmt(uf.get('chi_squared')),
                        _fmt(uf.get('background')),
                        _fmt(uf.get('background_err')),
                        uf.get('num_levels'),
                    ]
                    levels = uf.get('levels', [])
                    for lvl_idx in range(max_levels):
                        lv = levels[lvl_idx] if lvl_idx < len(levels) else {}
                        row += [
                            _fmt(lv.get('G')),
                            _fmt(lv.get('G_err')),
                            _fmt(lv.get('Rg')),
                            _fmt(lv.get('Rg_err')),
                            _fmt(lv.get('B')),
                            _fmt(lv.get('B_err')),
                            _fmt(lv.get('P')),
                            _fmt(lv.get('P_err')),
                            _fmt(lv.get('RgCutoff')),
                            _fmt(lv.get('ETA')),
                            _fmt(lv.get('ETA_err')),
                            _fmt(lv.get('PACK')),
                            _fmt(lv.get('PACK_err')),
                            lv.get('correlated'),
                            _fmt(lv.get('Sv')),
                            _fmt(lv.get('Invariant')),
                        ]
                else:
                    row += [None] * (4 + max_levels * 16)

            if show_sizes:
                if sz is not None:
                    row += [
                        sz.get('method'),
                        sz.get('shape'),
                        _fmt(sz.get('chi_squared')),
                        _fmt(sz.get('volume_fraction')),
                        _fmt(sz.get('rg')),
                        sz.get('n_iterations'),
                        _fmt(sz.get('q_power')),
                        _fmt(sz.get('contrast')),
                        _fmt(sz.get('aspect_ratio')),
                        _fmt(sz.get('background')),
                        _fmt(sz.get('error_scale')),
                        _fmt(sz.get('power_law_B')),
                        _fmt(sz.get('power_law_P')),
                        _fmt(sz.get('r_min')),
                        _fmt(sz.get('r_max')),
                        sz.get('n_bins'),
                        sz.get('log_spacing'),
                        _fmt(sz.get('cursor_q_min')),
                        _fmt(sz.get('cursor_q_max')),
                    ]
                else:
                    row += [None] * 19

            if show_simple_fits:
                if sf is not None:
                    sf_params_d  = sf.get('params', {})
                    sf_std_d     = sf.get('params_std', {})
                    sf_derived_d = sf.get('derived', {})
                    row += [
                        sf.get('model'),
                        _fmt(sf.get('chi_squared')),
                        _fmt(sf.get('reduced_chi_squared')),
                        sf.get('dof'),
                        _fmt(sf.get('q_min')),
                        _fmt(sf.get('q_max')),
                        sf.get('use_complex_bg'),
                    ]
                    for pname in _sf_param_names:
                        row.append(_fmt(sf_params_d.get(pname)))
                        row.append(_fmt(sf_std_d.get(pname)))
                    for dname in _sf_derived_names:
                        row.append(_fmt(sf_derived_d.get(dname)))
                else:
                    n_sf_cols = 7 + 2 * len(_sf_param_names) + len(_sf_derived_names)
                    row += [None] * n_sf_cols

            if show_waxs_peakfit:
                if wp is not None:
                    row += [
                        int(wp.get('n_peaks', 0)),
                        wp.get('bg_shape'),
                        _fmt(wp.get('chi_squared')),
                        _fmt(wp.get('reduced_chi_squared')),
                        int(wp.get('dof', 0)),
                        _fmt(wp.get('q_min')),
                        _fmt(wp.get('q_max')),
                    ]
                    peaks_list = wp.get('peaks', [])
                    peaks_std  = wp.get('peaks_std', [])
                    for pk_idx in range(max_wp_peaks):
                        if pk_idx < len(peaks_list):
                            pk   = peaks_list[pk_idx]
                            pstd = peaks_std[pk_idx] if pk_idx < len(peaks_std) else {}
                            # Area: prefer stored scalar; recompute for older files
                            area_v = pk.get('area')
                            area_s = pk.get('area_std')
                            if area_v is None:
                                from pyirena.core.waxs_peakfit import (
                                    peak_area, peak_area_std,
                                )
                                area_v = peak_area(pk.get('shape', 'Gauss'), pk)
                                area_s = peak_area_std(pk.get('shape', 'Gauss'),
                                                       pk, pstd)
                            row += [
                                pk.get('shape'),
                                _fmt(pk.get('Q0',   {}).get('value')),
                                _fmt(pstd.get('Q0')),
                                _fmt(pk.get('A',    {}).get('value')),
                                _fmt(pstd.get('A')),
                                _fmt(pk.get('FWHM', {}).get('value')),
                                _fmt(pstd.get('FWHM')),
                                _fmt(pk.get('eta',  {}).get('value') if 'eta' in pk else None),
                                _fmt(pstd.get('eta')),
                                _fmt(area_v),
                                _fmt(area_s),
                            ]
                        else:
                            row += [None] * _wp_peak_cols
                else:
                    row += [None] * (_wp_n_fixed_cols + max_wp_peaks * _wp_peak_cols)

            if show_modeling:
                if mod is not None:
                    pops = [p for p in mod.get('populations', [])
                            if p.get('enabled', True)]
                    row += [
                        _fmt(mod.get('chi_squared')),
                        _fmt(mod.get('background')),
                        _fmt(mod.get('q_min')),
                        _fmt(mod.get('q_max')),
                        len(pops),
                    ]
                    for k in range(max_mod_pops):
                        if k < len(pops):
                            pop = pops[k]
                            pt      = pop.get('pop_type', 'size_dist')
                            derived = pop.get('derived', {})
                            # All params are flat keys on pop dict
                            _ffp = pop.get('ff_params', {}) or {}
                            _sfp = pop.get('sf_params', {}) or {}
                            _dp = pop.get('dist_params', {}) or {}
                            row += [
                                pt,
                                pop.get('label', ''),
                                # size_dist
                                pop.get('dist_type'),
                                _fmt(pop.get('scale')),
                                _fmt(pop.get('contrast')),
                                pop.get('form_factor'),
                            ]
                            row += [_fmt(_dp.get(k)) for k in _mod_dist_keys]
                            row += [_fmt(_ffp.get(k)) for k in _mod_ff_keys]
                            row += [pop.get('structure_factor')]
                            row += [_fmt(_sfp.get(k)) for k in _mod_sf_keys]
                            row += [
                                _fmt(derived.get('volume_fraction')),
                                _fmt(derived.get('vol_mean_r')),
                                _fmt(derived.get('r_total_mean')),
                                # unified_level
                                _fmt(pop.get('G')),
                                _fmt(pop.get('Rg')),
                                _fmt(pop.get('B')),
                                _fmt(pop.get('P')),
                                _fmt(pop.get('RgCO')),
                                _fmt(pop.get('ETA')),
                                _fmt(pop.get('PACK')),
                                # diffraction_peak
                                pop.get('peak_type'),
                                _fmt(pop.get('position')),
                                _fmt(pop.get('amplitude')),
                                _fmt(pop.get('width')),
                                _fmt(pop.get('eta_voigt')),
                            ]
                        else:
                            row += [None] * len(_mod_pop_cols)
                else:
                    row += [None] * (5 + max_mod_pops * len(_mod_pop_cols))

            if show_saxs_morph:
                if sm is not None:
                    mm = sm.get('morphology_metrics')
                    # Convert booleans to "yes"/"no" strings for CSV-friendly output;
                    # falsy values (None when no metrics) → empty.
                    def _yn(v):
                        if v is None:
                            return None
                        return 'yes' if bool(v) else 'no'
                    row += [
                        _fmt(sm.get('chi_squared')),
                        _fmt(sm.get('volume_fraction')),
                        _fmt(sm.get('contrast')),
                        _fmt(sm.get('rg_A')),
                        # Specific surface area was added in the recent
                        # Porod-extension feature; pull from the file's
                        # morphology metrics group if not present in legacy
                        # results.
                        _fmt(sm.get('specific_surface_area_inv_A')),
                        sm.get('voxel_size'),
                        _fmt(sm.get('box_size_A')),
                        _fmt(sm.get('voxel_pitch_A')),
                        _fmt(getattr(mm, 'minority_volume_fraction', None)) if mm else None,
                        getattr(mm, 'n_clusters', None) if mm else None,
                        _fmt(getattr(mm, 'open_porosity_fraction', None) * 100
                              if (mm and mm.open_porosity_fraction is not None) else None),
                        _fmt(getattr(mm, 'closed_porosity_fraction', None) * 100
                              if (mm and mm.closed_porosity_fraction is not None) else None),
                        _yn(getattr(mm, 'percolating_x', None) if mm else None),
                        _yn(getattr(mm, 'percolating_y', None) if mm else None),
                        _yn(getattr(mm, 'percolating_z', None) if mm else None),
                        getattr(mm, 'euler_number', None) if mm else None,
                        _fmt(getattr(mm, 'pore_size_q25_A', None) if mm else None),
                        _fmt(getattr(mm, 'pore_size_median_A', None) if mm else None,),
                        _fmt(getattr(mm, 'pore_size_q75_A', None) if mm else None),
                    ]
                else:
                    row += [None] * 19

            rows.append(row)

        # ── Display ───────────────────────────────────────────────────────────
        default_path = os.path.join(
            self.current_folder or '',
            'pyIrena_TableOfResults.csv',
        )

        if self.tabulate_results_window is None:
            self.tabulate_results_window = TabulateResultsWindow()

        self.tabulate_results_window.set_data(headers, rows, default_path)
        self.status_label.setText(
            f"Tabulated {len(rows)} file(s) — use 'Save as CSV' in the results window."
        )

    # ─────────────────────────────────────────────────────────────────────────
    def export_to_ascii(self):
        """
        Export selected HDF5 files to plain-text .dat files in an
        'ascii_export' subfolder next to each source file.

        Checkbox semantics
        ------------------
        - 'Data' checkbox: gates the primary {stem}.dat file (Q, I, dI).
          Uncheck to export only model curves without re-writing data.
        - Each fit-result checkbox (Unified Fit, Simple Fits, Sizes,
          WAXS Peaks, Modeling) gates its own model .dat files when the
          'Also write model curves' option is enabled in Configure.

        File naming
        -----------
        - Primary:        {stem}.dat                (Q, I, dI)
        - Unified Fit:    {stem}_unif.dat           (Q, I_model, I_data, dI)
        - Simple Fits:    {stem}_simp.dat           (Q, I_model, I_data, dI)
        - WAXS Peaks:     {stem}_waxs.dat           (Q, I_fit, I_data, dI)
        - Size Dist:      {stem}_sdQI.dat           (Q, I_model, I_data, dI)
                          {stem}_sdSD.dat           (r, vol_dist, num_dist [, std])
        - Modeling:       {stem}_modQI.dat          (Q, I_total, I_data, dI)
                          {stem}_modP1.dat, _modP2.dat, …
                              size_dist:        (r, vol_dist, num_dist)
                              diffraction_peak: (Q, I_peak)
                              unified_level:    (no file written)

        For USAXS files only the desmeared sasdata group is exported;
        files with only a slit-smeared (_SMR) variant are silently skipped
        — that condition usually means the reduction never produced
        desmeared data, which is a problem with the file itself.
        """
        selected_items = self.file_list.selectedItems()
        if not selected_items:
            QMessageBox.warning(self, "No Selection",
                                "Please select files to export.")
            return

        sm = self.state_manager
        delim          = sm.get('data_selector', 'ascii_delimiter', ' ')
        precision      = int(sm.get('data_selector', 'ascii_precision', 7))
        include_header = bool(sm.get('data_selector', 'ascii_include_header', True))
        include_models = bool(sm.get('data_selector', 'ascii_include_models', True))
        err_frac       = float(sm.get('data_selector', 'error_fraction', 0.05))

        # Data checkbox gates the primary .dat file; model checkboxes gate
        # their respective model .dat files.
        include_data = self.data_checkbox.isChecked()
        model_flags = {
            "unif": self.unified_fit_result_checkbox.isChecked(),
            "simp": self.simple_fits_checkbox.isChecked(),
            "mod":  self.modeling_checkbox.isChecked(),
            "sd":   self.size_dist_checkbox.isChecked(),
            "waxs": self.waxs_peakfit_checkbox.isChecked(),
        }

        if not include_data and not any(model_flags.values()):
            self.status_label.setText(
                "Nothing to export — check 'Data' or a fit-result checkbox."
            )
            return

        file_paths = [
            Path(self.current_folder) / item.text()
            for item in selected_items
        ]

        from pyirena.io.ascii_export import export_dataset_to_ascii

        n_data        = 0
        n_models      = 0
        n_silent_skip = 0
        errors        = []
        out_dir_first = None

        for fp in file_paths:
            # Auto-convert text files to NXcanSAS before ASCII export
            if fp.suffix.lower() in ('.txt', '.dat'):
                try:
                    fp = ensure_nxcansas_sibling(fp, error_fraction=err_frac)
                except Exception as exc:
                    errors.append((fp.name, f'conversion failed: {exc}'))
                    continue
            if fp.suffix.lower() not in ('.h5', '.hdf5', '.hdf'):
                errors.append((fp.name, 'not an HDF5 file'))
                continue
            try:
                out_dir = fp.parent / 'ascii_export'
                manifest = export_dataset_to_ascii(
                    fp, out_dir,
                    delimiter=delim,
                    precision=precision,
                    include_header=include_header,
                    include_data=include_data,
                    include_models=include_models,
                    model_flags=model_flags,
                    error_fraction=err_frac,
                )
                if manifest.get('silently_skipped'):
                    n_silent_skip += 1
                    continue
                if manifest['data'] is not None:
                    n_data += 1
                if manifest['models']:
                    n_models += len(manifest['models'])
                if (manifest['data'] is not None or manifest['models']) \
                        and out_dir_first is None:
                    out_dir_first = out_dir
            except Exception as e:
                errors.append((fp.name, str(e)))

        parts = []
        if include_data:
            parts.append(f"{n_data} data file(s)")
        parts.append(f"{n_models} model file(s)")
        msg = "Exported " + ", ".join(parts)
        if out_dir_first:
            msg += f" → {out_dir_first}"
        if n_silent_skip:
            msg += f".  {n_silent_skip} file(s) silently skipped (no desmeared data)"
        if errors:
            preview = '; '.join(f"{n} ({r})" for n, r in errors[:3])
            msg += f".  Errors {len(errors)}: {preview}"
            if len(errors) > 3:
                msg += f" (+{len(errors) - 3} more)"
        self.status_label.setText(msg)

    def _prompt_dataset_choice(self, filename, datasets):
        """Thin wrapper — delegates to the shared function in data_loading."""
        return _prompt_dataset_choice_fn(self, filename, datasets)

    def _read_nxcansas_with_picker(self, path, filename):
        """Thin wrapper — delegates to the shared function in data_loading."""
        return _read_nxcansas_with_picker_fn(
            self, path, filename,
            status_cb=lambda msg: self.status_label.setText(msg),
        )

    def _load_data_for_tool(self, file_path: str):
        """Thin wrapper — delegates to the shared function in data_loading."""
        error_fraction = self.state_manager.get('data_selector', 'error_fraction', 0.05)
        return _load_data_file_fn(self, file_path, error_fraction=error_fraction)

    def launch_unified_fit(self):
        """Launch the Unified Fit model panel with selected data."""
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(
                self, "No Selection",
                "Please select one or more files to analyze with Unified Fit."
            )
            return

        file_path = os.path.join(self.current_folder, selected_items[0].text())
        res = self._load_data_for_tool(file_path)
        if res is None:
            return
        data, hdf5_path, display_name = res

        try:
            if self.unified_fit_window is None:
                self.unified_fit_window = UnifiedFitPanel()

            self.unified_fit_window.set_data(
                data['Q'], data['Intensity'], data.get('Error'),
                display_name, filepath=hdf5_path, is_nxcansas=True,
            )
            self.unified_fit_window.show()
            self.unified_fit_window.raise_()
            self.unified_fit_window.activateWindow()
            self.status_label.setText(f"Opened Unified Fit for {display_name}")

        except Exception as e:
            QMessageBox.critical(self, "Error",
                                 f"Error loading data for Unified Fit:\n{e}")
            self.status_label.setText(f"Error: {e}")

    def launch_sizes_fit(self):
        """Launch the Size Distribution fitting panel with selected data."""
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(
                self, "No Selection",
                "Please select one or more files to analyze with Size Distribution."
            )
            return

        file_path = os.path.join(self.current_folder, selected_items[0].text())
        res = self._load_data_for_tool(file_path)
        if res is None:
            return
        data, hdf5_path, display_name = res

        try:
            if self.sizes_fit_window is None:
                self.sizes_fit_window = SizesFitPanel()

            self.sizes_fit_window.set_data(
                data['Q'], data['Intensity'], data.get('Error'),
                display_name, filepath=hdf5_path, is_nxcansas=True,
            )
            self.sizes_fit_window.show()
            self.sizes_fit_window.raise_()
            self.sizes_fit_window.activateWindow()
            self.status_label.setText(f"Opened Size Distribution for {display_name}")

        except Exception as e:
            QMessageBox.critical(self, "Error",
                                 f"Error loading data for Size Distribution:\n{e}")
            self.status_label.setText(f"Error: {e}")

    def launch_modeling(self):
        """Launch the Modeling panel with the first selected file."""
        from pyirena.gui.modeling_panel import ModelingPanel
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(
                self, "No Selection",
                "Please select a file to open in Modeling.",
            )
            return

        file_path = os.path.join(self.current_folder, selected_items[0].text())
        res = self._load_data_for_tool(file_path)
        if res is None:
            return
        data, hdf5_path, display_name = res

        try:
            if self.modeling_window is None:
                self.modeling_window = ModelingPanel()

            self.modeling_window.set_data(
                data['Q'], data['Intensity'], data.get('Error'),
                display_name, filepath=hdf5_path, is_nxcansas=True,
            )
            self.modeling_window.show()
            self.modeling_window.raise_()
            self.modeling_window.activateWindow()
            self.status_label.setText(f"Opened Modeling for {display_name}")

        except Exception as e:
            QMessageBox.critical(self, "Error",
                                 f"Error loading data for Modeling:\n{e}")
            self.status_label.setText(f"Error: {e}")

    # ── Batch script fitting ───────────────────────────────────────────────

    def run_unified_script(self):
        """Batch-fit all selected files with Unified Fit."""
        self._run_batch_fit('unified')

    def run_sizes_script(self):
        """Batch-fit all selected files with Size Distribution."""
        self._run_batch_fit('sizes')

    def run_modeling_script(self):
        """Batch-fit all selected files with Modeling."""
        self._run_batch_fit('modeling')

    def launch_simple_fits(self):
        """Launch the Simple Fits panel with the first selected file."""
        from pyirena.gui.simple_fits_panel import SimpleFitsPanel
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(
                self, "No Selection",
                "Please select a file to open in Simple Fits.",
            )
            return

        file_path = os.path.join(self.current_folder, selected_items[0].text())
        res = self._load_data_for_tool(file_path)
        if res is None:
            return
        data, hdf5_path, display_name = res

        try:
            if self.simple_fits_window is None:
                self.simple_fits_window = SimpleFitsPanel()

            self.simple_fits_window.set_data(
                data['Q'], data['Intensity'], data.get('Error'),
                display_name, filepath=hdf5_path, is_nxcansas=True,
            )
            self.simple_fits_window.show()
            self.simple_fits_window.raise_()
            self.simple_fits_window.activateWindow()
            self.status_label.setText(f"Opened Simple Fits for {display_name}")

        except Exception as e:
            QMessageBox.critical(self, "Error",
                                 f"Error loading data for Simple Fits:\n{e}")
            self.status_label.setText(f"Error: {e}")

    def run_simple_fits_script(self):
        """Batch-fit all selected files with Simple Fits."""
        self._run_batch_fit('simple_fits')

    def launch_waxs_peakfit(self):
        """Open the WAXS Peak Fit panel with the first selected file."""
        from pyirena.gui.waxs_peakfit_panel import WAXSPeakFitPanel

        selected_items = self.file_list.selectedItems()
        if not selected_items:
            QMessageBox.warning(
                self, "No Selection",
                "Please select a file to open in WAXS Peak Fit.",
            )
            return

        file_path = os.path.join(self.current_folder, selected_items[0].text())
        res = self._load_data_for_tool(file_path)
        if res is None:
            return
        data, hdf5_path, display_name = res

        try:
            if self.waxs_peakfit_window is None:
                self.waxs_peakfit_window = WAXSPeakFitPanel()

            self.waxs_peakfit_window.set_data(
                data['Q'], data['Intensity'], data.get('Error'),
                label=display_name, filepath=hdf5_path, is_nxcansas=True,
            )
            self.waxs_peakfit_window.show()
            self.waxs_peakfit_window.raise_()
            self.waxs_peakfit_window.activateWindow()
            self.status_label.setText(f"Opened WAXS Peak Fit for {display_name}")

        except Exception as e:
            QMessageBox.critical(self, "Error",
                                 f"Error loading data for WAXS Peak Fit:\n{e}")
            self.status_label.setText(f"Error: {e}")

    def run_waxs_peakfit_script(self):
        """Batch-fit all selected files with WAXS Peak Fit."""
        self._run_batch_fit('waxs_peakfit')

    def launch_saxs_morph(self):
        """Open the SAXS Morph (3D voxelgram) panel with the first selected file."""
        from pyirena.gui.saxs_morph_panel import SaxsMorphPanel

        selected_items = self.file_list.selectedItems()
        if not selected_items:
            QMessageBox.warning(
                self, "No Selection",
                "Please select a file to open in SAXS Morph.",
            )
            return

        file_path = os.path.join(self.current_folder, selected_items[0].text())
        res = self._load_data_for_tool(file_path)
        if res is None:
            return
        data, hdf5_path, display_name = res

        try:
            if self.saxs_morph_window is None:
                self.saxs_morph_window = SaxsMorphPanel()
                # SaxsMorphPanel sets WA_DeleteOnClose so the embedded
                # PyVista QtInteractor is fully destroyed on close.  When
                # that happens, clear our cached reference so the next
                # launch creates a fresh panel with a live VTK window.
                self.saxs_morph_window.destroyed.connect(
                    self._on_saxs_morph_destroyed
                )

            self.saxs_morph_window.set_data(
                data['Q'], data['Intensity'], data.get('Error'),
                filename=display_name, filepath=hdf5_path, is_nxcansas=True,
            )
            self.saxs_morph_window.show()
            self.saxs_morph_window.raise_()
            self.saxs_morph_window.activateWindow()
            self.status_label.setText(f"Opened SAXS Morph for {display_name}")

        except Exception as e:
            QMessageBox.critical(self, "Error",
                                 f"Error loading data for SAXS Morph:\n{e}")
            self.status_label.setText(f"Error: {e}")

    def _on_saxs_morph_destroyed(self, _obj=None):
        """Clear the cached SAXS Morph panel reference once Qt has
        destroyed the underlying widget (triggered by WA_DeleteOnClose).
        """
        self.saxs_morph_window = None

    def run_saxs_morph_script(self):
        """Batch-fit all selected files with SAXS Morph."""
        self._run_batch_fit('saxs_morph')

    def launch_data_merge(self):
        """Open the Data Merge tool, pre-populated with the current folder as DS1."""
        from pyirena.gui.data_merge_panel import DataMergePanel

        if self.data_merge_window is None:
            self.data_merge_window = DataMergePanel(
                state_manager=self.state_manager,
            )
            # Pre-populate DS1 only when there is no saved state — if the user
            # has previously set up folders they are restored from state instead.
            saved_dm = self.state_manager.get('data_merge') or {}
            if not saved_dm.get('folder1') and self.current_folder:
                self.data_merge_window.set_folder(1, self.current_folder)

        self.data_merge_window.show()
        self.data_merge_window.raise_()
        self.data_merge_window.activateWindow()
        self.status_label.setText("Data Merge tool opened.")

    def launch_data_manipulation(self):
        """Open the Data Manipulation tool."""
        from pyirena.gui.data_manipulation_panel import DataManipulationPanel

        if self.data_manip_window is None:
            self.data_manip_window = DataManipulationPanel(
                state_manager=self.state_manager,
            )

        if self.current_folder:
            self.data_manip_window.set_folder(self.current_folder)

        self.data_manip_window.show()
        self.data_manip_window.raise_()
        self.data_manip_window.activateWindow()
        self.status_label.setText("Data Manipulation tool opened.")

    def launch_contrast(self):
        """Open the Scattering Contrast Calculator."""
        from pyirena.gui.contrast_panel import ContrastPanel

        if self.contrast_window is None:
            self.contrast_window = ContrastPanel(
                state_manager=self.state_manager,
            )

        self.contrast_window.show()
        self.contrast_window.raise_()
        self.contrast_window.activateWindow()
        self.status_label.setText("Scattering Contrast Calculator opened.")

    def launch_fractals(self):
        """Open the Fractals (mass fractal aggregate) visualization tool."""
        from pyirena.gui.fractals_panel import FractalsGraphWindow

        if self.fractals_window is None:
            self.fractals_window = FractalsGraphWindow(
                state_manager=self.state_manager,
            )
            # WA_DeleteOnClose is set on the window; clear our cached reference
            # once Qt destroys the panel so the next launch creates a fresh one
            # with a live VTK render window.
            self.fractals_window.destroyed.connect(self._on_fractals_destroyed)

        self.fractals_window.show()
        self.fractals_window.raise_()
        self.fractals_window.activateWindow()
        self.status_label.setText("Fractals tool opened.")

    def _on_fractals_destroyed(self, _obj=None):
        """Clear the cached Fractals window reference once Qt has destroyed
        the underlying widget (triggered by WA_DeleteOnClose).
        """
        self.fractals_window = None

    # ── Igor packed experiment import ──────────────────────────────────────
    def launch_igor_import(self):
        """Pick an Igor experiment (.pxp or .h5xp), extract its
        USAXS/SAXS/WAXS data to NeXus files, and offer to load the
        output folder as the current data folder.
        """
        from pyirena.batch import igor_to_nexus

        start_dir = self.last_folder if self.last_folder else QDir.homePath()
        pxp_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open Igor Experiment",
            start_dir,
            "Igor Experiment (*.pxp *.h5xp);;"
            "Igor Packed Experiment (*.pxp);;"
            "Igor HDF5 Packed Experiment (*.h5xp);;"
            "All Files (*)",
        )
        if not pxp_path:
            return

        pxp_path = Path(pxp_path)
        # Show the import options dialog.
        dlg = _IgorImportDialog(pxp_path, parent=self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return

        opts = dlg.options()
        # Run extraction synchronously — these files are small (~100s of KB
        # each) and the parse is fast even for 16 MB pxp inputs. If users
        # report 100+ MB experiments hanging the UI, move this to a
        # QThread; for now keep it simple.
        try:
            self.status_label.setText(f"Importing {pxp_path.name}…")
            QApplication.processEvents()
            result = igor_to_nexus(
                igor_file=str(pxp_path),
                output_folder=opts['output_folder'],
                techniques=opts['techniques'],
                overwrite=opts['overwrite'],
            )
        except Exception as exc:
            QMessageBox.critical(
                self, "Import error",
                f"Failed to import {pxp_path.name}:\n\n{exc}",
            )
            self.status_label.setText("Import failed.")
            return

        if result is None:
            QMessageBox.warning(
                self, "Import returned no result",
                f"No data was extracted from {pxp_path.name}.",
            )
            return

        # Summarise
        msg_lines = [
            f"Imported {pxp_path.name}",
            "",
            f"Output folder: {result['output_folder']}",
            f"Files written: {result['n_written']}",
            f"Skipped:       {result['n_skipped']}  (folders with no recognised wave triple)",
            f"Errors:        {result['n_errors']}",
        ]
        if result.get('n_unparseable_records'):
            msg_lines.append(
                f"Note: {result['n_unparseable_records']} wave record(s) in the "
                f"pxp could not be parsed (skipped)."
            )

        # Igor 8/10 long-name records (folders/waves > 31 chars) —
        # igor2 cannot decode these, so samples are silently missing
        # from the output. Surface this prominently so the user knows
        # to re-save the experiment as .h5xp.
        n_longname = result.get('n_igor8_longname_markers') or 0
        if n_longname > 0:
            msg_lines.append("")
            msg_lines.append(
                f"⚠ WARNING: this .pxp contains {n_longname} Igor-8 long-name "
                f"record(s) (folders or waves with names > 31 chars)."
            )
            msg_lines.append(
                "  igor2 cannot decode these, so SOME SAMPLES ARE MISSING "
                "from the output above."
            )
            msg_lines.append(
                "  Workaround: in Igor Pro, save the experiment as .h5xp "
                "instead and re-import here."
            )

        # Per-technique tally
        tech_count: Dict[str, int] = {}
        for fr in result['files']:
            if fr['status'] == 'ok':
                tech_count[fr['technique']] = tech_count.get(fr['technique'], 0) + 1
        if tech_count:
            msg_lines.append("")
            msg_lines.append("By technique:")
            for t in sorted(tech_count):
                msg_lines.append(f"  {t}: {tech_count[t]}")

        self.status_label.setText(
            f"Imported {result['n_written']} files to {Path(result['output_folder']).name}"
        )

        # Offer to load the output folder
        msg_box = QMessageBox(self)
        # Warning icon when long-name markers were seen, since the
        # user almost certainly has missing samples.
        if n_longname > 0:
            msg_box.setIcon(QMessageBox.Icon.Warning)
            msg_box.setWindowTitle("Import complete — some samples missing")
        else:
            msg_box.setIcon(QMessageBox.Icon.Information)
            msg_box.setWindowTitle("Import complete")
        msg_box.setText("\n".join(msg_lines))
        if result['n_written'] > 0:
            load_btn = msg_box.addButton(
                "Load output folder",
                QMessageBox.ButtonRole.AcceptRole,
            )
            tech_subdirs = sorted(tech_count)
            if len(tech_subdirs) == 1:
                load_btn.setText(f"Load {tech_subdirs[0]} folder")
        msg_box.addButton(QMessageBox.StandardButton.Close)
        msg_box.exec()

        if (msg_box.clickedButton() is not None
                and msg_box.clickedButton().text().startswith("Load")):
            # If exactly one technique was written, jump directly into it
            # so the file list shows files immediately. Otherwise open the
            # parent output folder and let the user pick.
            target = Path(result['output_folder'])
            if len(tech_count) == 1:
                only_tech = next(iter(tech_count))
                tech_dir = target / only_tech
                if tech_dir.is_dir():
                    target = tech_dir
            self._load_folder(str(target))

    def _load_folder(self, folder: str) -> None:
        """Set *folder* as the current data folder and refresh the file list.

        Same effect as a successful ``select_folder`` dialog, but accepts a
        path string directly. Used by import workflows.
        """
        self.current_folder = folder
        self.last_folder = folder
        self.save_last_folder(folder)
        self.folder_label.setText(folder)
        self.folder_label.setStyleSheet("color: #2c3e50;")
        self.refresh_button.setEnabled(True)
        self.refresh_file_list()
        self.status_label.setText(f"Folder: {os.path.basename(folder)}")

    # ── Visualization for "Create Graph" with 3D-tool checkboxes ─────────
    #
    # When a user checks "Fractals" or "3D saxsMorph" and clicks
    # **Create Graph**, we do NOT open the full SAXS-Morph / Fractals
    # tool (which would force a recompute).  We just load the saved 3D
    # voxelgram from the NeXus file and display it in a standalone
    # `VoxelViewerWindow` (the same `Slice2DViewer` + `Voxel3DViewer`
    # pair the parent tools use, but read-only).  Tabulate / Report /
    # Export-ASCII don't apply to either tool, so they ignore these
    # checkboxes.

    def _open_saxs_morph_viewer(self, file_paths) -> int:
        """Open a standalone 2D + 3D viewer pre-loaded with stored saxsMorph
        voxelgrams from the selected files.  Returns the number of files
        whose results were displayed."""
        from pyirena.io.nxcansas_saxs_morph import load_saxs_morph_results
        from pyirena.gui.saxs_morph_3d import VoxelViewerWindow

        items = []
        for fp in file_paths:
            try:
                res = load_saxs_morph_results(Path(fp))
            except Exception:
                continue
            if not res or 'voxelgram' not in res:
                continue
            vox = res['voxelgram']
            pitch = float(res.get('voxel_pitch_A', 1.0)) or 1.0
            label = (f'{Path(fp).name}  '
                     f'(saxsMorph, {vox.shape[0]}³, pitch={pitch:.3g} Å)')
            items.append({'voxelgram': vox, 'pitch_A': pitch, 'label': label})

        if not items:
            QMessageBox.information(
                self, "No saxsMorph results",
                "None of the selected files contain stored 3D saxsMorph results.",
            )
            return 0

        win = VoxelViewerWindow(items, title='3D saxsMorph — stored results',
                                 parent=self)
        # Keep a reference so Python doesn't GC it; clear on close.
        self._standalone_voxel_viewers.append(win)
        win.destroyed.connect(
            lambda _o=None, w=win: self._standalone_voxel_viewers.remove(w)
            if w in self._standalone_voxel_viewers else None
        )
        win.show()
        win.raise_()
        win.activateWindow()
        return len(items)

    def _open_fractals_viewer(self, file_paths) -> int:
        """Open a standalone 2D + 3D viewer pre-loaded with stored fractal
        aggregates from the selected files.  Each saved aggregate is
        re-voxelised on the fly with the same display geometry the
        Fractals tool uses (`oversample=10, sphere_voxel_radius=10` —
        the chunky-sphere render that makes the aggregate look connected).
        Returns the number of aggregates displayed."""
        from pyirena.io.nxcansas_fractals import (
            list_fractal_aggregates, load_fractal_aggregate,
        )
        from pyirena.core.fractals import voxelize
        from pyirena.gui.saxs_morph_3d import VoxelViewerWindow

        items = []
        for fp in file_paths:
            try:
                entries = list_fractal_aggregates(Path(fp))
            except Exception:
                continue
            for e in entries:
                try:
                    agg = load_fractal_aggregate(Path(fp), e['group_path'])
                except Exception:
                    continue
                try:
                    voxelgram, _pitch_lattice = voxelize(
                        agg.positions,
                        oversample=10, sphere_voxel_radius=10,
                    )
                except Exception:
                    continue
                pitch_A = float(agg.params.primary_diameter) / 10.0
                label = (f'{Path(fp).name} : {e["name"]}  '
                         f'(Z={agg.params.z}, df={agg.params.df:.2f})')
                items.append({'voxelgram': voxelgram,
                              'pitch_A': pitch_A,
                              'label': label})

        if not items:
            QMessageBox.information(
                self, "No Fractals results",
                "None of the selected files contain stored Fractals aggregates.",
            )
            return 0

        win = VoxelViewerWindow(items, title='Fractals — stored aggregates',
                                 parent=self)
        self._standalone_voxel_viewers.append(win)
        win.destroyed.connect(
            lambda _o=None, w=win: self._standalone_voxel_viewers.remove(w)
            if w in self._standalone_voxel_viewers else None
        )
        win.show()
        win.raise_()
        win.activateWindow()
        return len(items)

    def launch_hdf5_viewer(self):
        """Open the Data Explorer for the current folder."""
        from pyirena.gui.hdf5viewer import HDF5ViewerWindow

        if self.hdf5_viewer_window is None:
            self.hdf5_viewer_window = HDF5ViewerWindow(
                initial_folder=self.current_folder or None,
                state_manager=self.state_manager,
            )

        self.hdf5_viewer_window.show()
        self.hdf5_viewer_window.raise_()
        self.hdf5_viewer_window.activateWindow()
        self.status_label.setText("Data Explorer opened.")

    def open_config_manager(self):
        """Open the Config Manager dialog for the current folder's pyirena_config.json."""
        config_file = self._find_config_file()
        if not config_file:
            QMessageBox.information(
                self, "No Config File",
                "No pyirena_config.json found.\n"
                "Export parameters from a fit panel first, then try again.",
            )
            return
        dlg = ConfigManagerDialog(config_file, parent=self)
        dlg.exec()

    def locate_logs(self):
        """Open the pyirena log folder (~/.pyirena/logs) in the file browser."""
        from pyirena.logging_setup import get_log_dir
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(get_log_dir())))

    def _find_config_file(self) -> Optional[str]:
        """
        Return the path to pyirena_config.json.
        Looks first in the current folder; if not found, prompts the user.
        """
        if self.current_folder:
            candidate = os.path.join(self.current_folder, 'pyirena_config.json')
            if os.path.isfile(candidate):
                return candidate
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select pyirena_config.json",
            self.current_folder or str(Path.home()),
            "pyIrena Config (*.json);;All Files (*)",
        )
        return path or None

    def _run_batch_fit(self, tool: str):
        """Start a BatchWorker thread for the given tool ('unified' or 'sizes')."""
        if self._batch_worker is not None and self._batch_worker.isRunning():
            QMessageBox.information(
                self, "Busy", "A batch fit is already in progress — please wait."
            )
            return

        selected_items = self.file_list.selectedItems()
        if not selected_items:
            QMessageBox.warning(
                self, "No Selection", "Please select one or more files to fit."
            )
            return

        config_file = self._find_config_file()
        if not config_file:
            QMessageBox.warning(
                self,
                "No Config File",
                "No pyirena_config.json found.\n"
                "Export parameters from a GUI fit panel first, then try again.",
            )
            return

        file_paths = [
            os.path.join(self.current_folder, item.text())
            for item in selected_items
        ]
        _tool_display = {'unified': "Unified Fit", 'modeling': "Modeling",
                         'sizes': "Size Distribution",
                         'simple_fits': "Simple Fits", 'waxs_peakfit': "WAXS Peak Fit"}
        tool_name = _tool_display.get(tool, tool)
        self._set_batch_status(
            f"⏳  Starting {tool_name} batch on {len(file_paths)} file(s)…", 'working'
        )

        with_unc = bool(self.state_manager.get('data_selector', 'batch_mc_uncertainty', False))
        n_runs   = int(self.state_manager.get('data_selector', 'batch_mc_n_runs', 10))
        self._batch_worker = BatchWorker(
            tool, file_paths, config_file,
            with_uncertainty=with_unc, n_mc_runs=n_runs,
            parent=self,
        )
        self._batch_worker.progress.connect(self._on_batch_progress)
        self._batch_worker.finished.connect(self._on_batch_finished)
        self._batch_worker.start()

    def _set_batch_status(self, text: str, state: str):
        """Update the batch status label colour and text."""
        colour_map = {
            'working': ('#f39c12', 'white'),
            'done':    ('#27ae60', 'white'),
            'partial': ('#e67e22', 'white'),
            'error':   ('#e74c3c', 'white'),
        }
        bg, fg = colour_map.get(state, ('#7f8c8d', 'white'))
        self.batch_status_label.setStyleSheet(
            f"QLabel {{ padding: 5px 8px; border-radius: 4px; font-size: 12px; "
            f"background-color: {bg}; color: {fg}; }}"
        )
        self.batch_status_label.setText(text)
        self.batch_status_label.setVisible(True)

    def _on_batch_progress(self, msg: str):
        self._set_batch_status(f"⏳  {msg}", 'working')

    def _on_batch_finished(self, n_ok: int, n_fail: int, messages: list):
        total = n_ok + n_fail
        if n_fail == 0:
            summary = f"✓  Done: all {total} fit(s) succeeded"
            state = 'done'
        elif n_ok == 0:
            summary = f"✗  Done: all {total} fit(s) failed"
            state = 'error'
        else:
            summary = f"⚠  Done: {n_ok} succeeded, {n_fail} failed (out of {total})"
            state = 'partial'
        self._set_batch_status(summary, state)
        self.status_label.setText(f"Batch complete — {summary.lstrip('✓✗⚠ ')}")

        # Always show a scrollable results dialog
        self._show_batch_results_dialog("Batch Fit — Results", summary, messages)

    def _show_batch_results_dialog(
        self, title: str, summary: str, messages: list
    ) -> None:
        """Show batch results in a scrollable, non-blocking floating window."""
        dlg = QDialog(self, Qt.WindowType.Window)
        dlg.setWindowTitle(title)
        dlg.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)

        layout = QVBoxLayout(dlg)
        layout.setContentsMargins(12, 12, 12, 8)
        layout.setSpacing(8)

        # Summary line at top
        lbl = QLabel(summary)
        lbl.setWordWrap(True)
        layout.addWidget(lbl)

        # Scrollable per-file results list
        list_w = QListWidget()
        list_w.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        for msg in messages:
            list_w.addItem(msg)
        layout.addWidget(list_w, 1)

        # OK / Close button always visible at bottom
        btn_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok)
        btn_box.accepted.connect(dlg.close)
        layout.addWidget(btn_box)

        # Fixed reasonable size — list widget scrolls automatically
        dlg.resize(540, 460)

        dlg.show()
        dlg.raise_()
        dlg.activateWindow()

    def open_configure_dialog(self):
        """Open the extensible configuration dialog for data loading options."""
        current = {
            spec['key']: self.state_manager.get(
                'data_selector', spec['key'], spec.get('default')
            )
            for spec in DataSelectorConfigDialog.FIELD_SPECS
        }
        dialog = DataSelectorConfigDialog(current, parent=self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            values = dialog.get_values()
            for key, value in values.items():
                self.state_manager.set('data_selector', key, value)
            self.state_manager.save()

    def load_last_folder(self):
        """Load the last used folder from state and set it if it exists."""
        last_folder = self.state_manager.get('data_selector', 'last_folder')

        if last_folder and os.path.isdir(last_folder):
            # Folder exists, use it
            self.last_folder = last_folder
            self.current_folder = last_folder
            log.info("Restored last folder: %s", last_folder)
        else:
            # Folder doesn't exist or wasn't saved, use home directory
            self.last_folder = str(Path.home())
            self.current_folder = None
            if last_folder:
                log.info("Last folder no longer exists: %s", last_folder)
                log.info("Starting in home directory: %s", self.last_folder)

    def save_last_folder(self, folder: str):
        """Save the current folder to state."""
        self.state_manager.set('data_selector', 'last_folder', folder)
        self.state_manager.save()

    def show_about(self):
        """Show about dialog."""
        from pyirena import __version__ as _version
        QMessageBox.about(
            self,
            "About pyIrena",
            f"""<h3>pyIrena</h3>
            <p><b>Python tools for small-angle scattering data analysis</b></p>
            <p>Version: {_version}</p>
            <p>pyIrena provides tools for analyzing SAXS/SANS/USAXS data,
            including the Unified Fit model for hierarchical structures.</p>
            <p>Based on Irena SAS package for Igor Pro by Jan Ilavsky</p>
            <p><a href='https://github.com/jilavsky/SAXS_IgorCode'>
            Original Irena Package</a></p>
            """
        )


def main():
    """Main entry point for the data selector GUI."""
    app = QApplication(sys.argv)

    # Set application style
    app.setStyle('Fusion')

    # Create and show the main window
    window = DataSelectorPanel()
    window.show()

    sys.exit(app.exec())
