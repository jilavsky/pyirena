"""
Data Selector GUI for pyIrena.

This module provides a GUI panel for selecting data files and displaying
their content as graphs.
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import List, Optional

try:
    from PySide6.QtWidgets import (
        QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
        QListWidget, QLabel, QLineEdit, QFileDialog, QComboBox,
        QAbstractItemView, QMessageBox, QMenuBar, QMenu,
        QDialog, QFormLayout, QDialogButtonBox, QGroupBox, QCheckBox, QColorDialog,
    )
    from PySide6.QtCore import Qt, QDir
    from PySide6.QtGui import QAction, QDoubleValidator
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
            QListWidget, QLabel, QLineEdit, QFileDialog, QComboBox,
            QAbstractItemView, QMessageBox, QMenuBar, QMenu,
            QDialog, QFormLayout, QDialogButtonBox, QGroupBox, QCheckBox, QColorDialog,
        )
        from PyQt6.QtCore import Qt, QDir
        from PyQt6.QtGui import QAction, QDoubleValidator
    except ImportError:
        raise ImportError(
            "Neither PySide6 nor PyQt6 found. Install with: pip install PySide6"
        )

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')  # Use Qt backend for matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from pyirena.io.hdf5 import readGenericNXcanSAS, readTextFile
from pyirena.io.nxcansas_unified import load_unified_fit_results
from pyirena.gui.unified_fit import UnifiedFitPanel
from pyirena.gui.sizes_panel import SizesFitPanel
from pyirena.state import StateManager


def _build_report(file_path: str,
                  data_info: Optional[dict] = None,
                  fit_results: Optional[dict] = None,
                  sizes_results: Optional[dict] = None) -> str:
    """
    Build a Markdown report string.

    Args:
        file_path:     Absolute path to the source file.
        data_info:     Dict with keys 'Q', 'I', 'I_error' (optional array).
                       Pass None to omit the data section.
        fit_results:   Dict from load_unified_fit_results().
                       Pass None to omit the unified fit section.
        sizes_results: Dict from load_sizes_results().
                       Pass None to omit the size distribution section.

    Returns:
        Multi-line Markdown string ready to be written to a .md file.
    """
    from datetime import datetime as _dt

    filename = os.path.basename(file_path)
    now = _dt.now().strftime('%Y-%m-%d %H:%M:%S')

    L = []

    # ── Header ───────────────────────────────────────────────────────────────
    L += [
        "# pyIrena Report",
        "",
        "| | |",
        "|---|---|",
        f"| **File** | `{filename}` |",
        f"| **Report generated** | {now} |",
    ]
    if fit_results is not None:
        L += [
            f"| **Fit timestamp** | {fit_results.get('timestamp', 'unknown')} |",
            f"| **Program** | {fit_results.get('program', 'pyirena')} |",
        ]
    L.append("")

    # ── Data summary ─────────────────────────────────────────────────────────
    if data_info is not None:
        Q = data_info['Q']
        I = data_info['I']
        I_error = data_info.get('I_error')
        L += [
            "## Data Summary",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Q range | {Q.min():.4g} – {Q.max():.4g} Å⁻¹ |",
            f"| Intensity range | {I.min():.4e} – {I.max():.4e} cm⁻¹ |",
            f"| Data points | {len(Q)} |",
        ]
        if I_error is not None:
            L.append(
                f"| Uncertainty range | {I_error.min():.4e} – {I_error.max():.4e} cm⁻¹ |"
            )
        L.append("")

    # ── Fit quality ──────────────────────────────────────────────────────────
    if fit_results is not None:
        chi2      = fit_results['chi_squared']
        bg        = fit_results['background']
        n_levels  = fit_results['num_levels']
        Q_fit     = fit_results['Q']
        residuals = fit_results['residuals']
        levels    = fit_results.get('levels', [])

        L += [
            "## Fit Quality",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Chi-squared (χ²) | {chi2:.4f} |",
            f"| Number of levels | {n_levels} |",
            f"| Background | {bg:.4e} cm⁻¹ |",
            f"| Q range (fit) | {Q_fit.min():.4g} – {Q_fit.max():.4g} Å⁻¹ |",
            f"| Data points (fit) | {len(Q_fit)} |",
            f"| Residuals mean | {np.mean(residuals):.4f} |",
            f"| Residuals std dev | {np.std(residuals):.4f} |",
            f"| Max \\|residual\\| | {np.max(np.abs(residuals)):.4f} |",
            "",
        ]

        # ── Level parameters ──────────────────────────────────────────────
        for i, level in enumerate(levels):
            lnum = i + 1
            L.append(f"## Level {lnum} Parameters")
            L.append("")

            has_mc = any(f'{p}_err' in level for p in ('G', 'Rg', 'B', 'P', 'ETA', 'PACK'))
            if has_mc:
                L += ["| Parameter | Value | Uncertainty (1σ) |",
                      "|-----------|-------|------------------|"]
            else:
                L += ["| Parameter | Value |",
                      "|-----------|-------|"]

            def _row(label, key, unit='', fmt='.4e'):
                val = level.get(key)
                if val is None:
                    return
                val_str = f"{val:{fmt}}{unit}"
                if has_mc:
                    err = level.get(f'{key}_err', 0.0)
                    err_str = f"± {err:{fmt}}{unit}" if err > 0 else "—"
                    L.append(f"| {label} | {val_str} | {err_str} |")
                else:
                    L.append(f"| {label} | {val_str} |")

            _row('G',  'G')
            _row('Rg', 'Rg', ' Å')
            _row('B',  'B')
            _row('P',  'P',  fmt='.4f')

            rgcut = level.get('RgCutoff', 0.0)
            if isinstance(rgcut, float) and rgcut > 0.01:
                _row('RgCutoff', 'RgCutoff', ' Å')

            if level.get('correlated', False):
                _row('ETA',  'ETA',  ' Å', fmt='.2f')
                _row('PACK', 'PACK', '',   fmt='.4f')

            sv = level.get('Sv')
            if sv is not None and sv != 'N/A' and isinstance(sv, (int, float)):
                _row('Sv', 'Sv', ' m²/cm³')

            inv = level.get('Invariant')
            if inv is not None and inv != 'N/A' and isinstance(inv, (int, float)):
                _row('Invariant', 'Invariant', ' cm⁻⁴')

            L.append("")

    # ── Size Distribution results ─────────────────────────────────────────────
    if sizes_results is not None:
        params    = sizes_results.get('params', {})
        r_grid    = sizes_results.get('r_grid')
        dist      = sizes_results.get('distribution')
        residuals = sizes_results.get('residuals')

        chi2   = params.get('chi_squared',     float('nan'))
        vf     = params.get('volume_fraction',  float('nan'))
        rg     = params.get('rg',              float('nan'))
        method = params.get('method',          'unknown')
        shape  = params.get('shape',           'unknown')
        n_bins = params.get('n_bins',          0)
        r_min  = params.get('r_min',           float('nan'))
        r_max  = params.get('r_max',           float('nan'))

        peak_r = float('nan')
        if dist is not None and r_grid is not None and len(dist) > 0:
            peak_r = float(r_grid[int(np.argmax(dist))])

        L += [
            "## Size Distribution",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Method | {method} |",
            f"| Particle shape | {shape} |",
            f"| Chi-squared (χ²) | {chi2:.4g} |",
            f"| Volume fraction | {vf:.4e} |",
            f"| Rg | {rg:.4g} Å |",
            f"| Peak r | {peak_r:.4g} Å |",
            f"| r min | {r_min:.4g} Å |",
            f"| r max | {r_max:.4g} Å |",
            f"| Bins | {n_bins} |",
        ]
        if residuals is not None:
            L += [
                f"| Residuals mean | {np.mean(residuals):.4f} |",
                f"| Residuals std dev | {np.std(residuals):.4f} |",
            ]
        L.append("")

    L += ["---", "*Generated by pyIrena*", ""]
    return "\n".join(L)


class DataSelectorConfigDialog(QDialog):
    """
    Extensible configuration dialog for the Data Selector.

    Settings are defined as a list of field specifications (FIELD_SPECS).
    Adding a new configurable parameter only requires adding one entry to that list.

    Supported field types
    ---------------------
    'float'  — QLineEdit with QDoubleValidator
    'int'    — QLineEdit with integer validation
    'str'    — plain QLineEdit
    'bool'   — QCheckBox
    'color'  — QPushButton that opens QColorDialog (stores hex color string)
    """

    # ---------------------------------------------------------------------------
    # Field specifications — add new settings here
    # ---------------------------------------------------------------------------
    FIELD_SPECS = [
        {
            'group':    'Text File Options',
            'key':      'error_fraction',
            'label':    'Generated uncertainty fraction',
            'tooltip':  (
                'When a text file has only two columns (Q and I) and no uncertainty\n'
                'column, the uncertainty is generated as:  σ = I × this_value\n'
                'Default: 0.05  (5 % of intensity)'
            ),
            'type':     'float',
            'default':  0.05,
            'min':      0.0,
            'max':      100.0,
            'decimals': 4,
        },
        # -----------------------------------------------------------------------
        # Future settings — just append a dict here, no other code changes needed
        # -----------------------------------------------------------------------
        # {
        #     'group':   'Text File Options',
        #     'key':     'q_units_scale',
        #     'label':   'Q unit scale factor',
        #     'tooltip': 'Multiply Q by this factor on load (1.0 = no change)',
        #     'type':    'float',
        #     'default': 1.0,
        #     'min':     0.0,
        #     'max':     1000.0,
        #     'decimals': 6,
        # },
        # {
        #     'group':   'Display',
        #     'key':     'plot_color',
        #     'label':   'Default plot color',
        #     'tooltip': 'Color used for single-file plots',
        #     'type':    'color',
        #     'default': '#3498db',
        # },
    ]

    def __init__(self, current_values: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Data Selector — Configure")
        self.setMinimumWidth(420)
        self._widgets = {}   # key -> (field_type, widget)
        self._init_ui(current_values)

    def _init_ui(self, current_values: dict):
        outer = QVBoxLayout()
        outer.setSpacing(12)

        # Group fields by their 'group' key
        groups = {}
        for spec in self.FIELD_SPECS:
            g = spec.get('group', 'General')
            groups.setdefault(g, []).append(spec)

        for group_name, specs in groups.items():
            box = QGroupBox(group_name)
            form = QFormLayout()
            form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)

            for spec in specs:
                key       = spec['key']
                label_txt = spec['label']
                ftype     = spec['type']
                value     = current_values.get(key, spec.get('default', ''))
                tooltip   = spec.get('tooltip', '')

                if ftype in ('float', 'int'):
                    widget = QLineEdit(str(value))
                    validator = QDoubleValidator(
                        float(spec.get('min', -1e300)),
                        float(spec.get('max',  1e300)),
                        int(spec.get('decimals', 6)),
                    )
                    widget.setValidator(validator)
                    widget.setMaximumWidth(120)

                elif ftype == 'bool':
                    widget = QCheckBox()
                    widget.setChecked(bool(value))

                elif ftype == 'color':
                    widget = QPushButton()
                    widget._color = str(value)
                    widget.setStyleSheet(f"background-color: {value};")
                    widget.setFixedSize(60, 24)
                    widget.clicked.connect(
                        lambda checked, btn=widget: self._pick_color(btn)
                    )

                else:   # 'str'
                    widget = QLineEdit(str(value))

                if tooltip:
                    widget.setToolTip(tooltip)

                lbl = QLabel(label_txt)
                if tooltip:
                    lbl.setToolTip(tooltip)

                form.addRow(lbl, widget)
                self._widgets[key] = (ftype, widget)

            box.setLayout(form)
            outer.addWidget(box)

        # OK / Cancel buttons
        btn_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        btn_box.accepted.connect(self.accept)
        btn_box.rejected.connect(self.reject)
        outer.addWidget(btn_box)

        self.setLayout(outer)

    def _pick_color(self, btn):
        color = QColorDialog.getColor(parent=self)
        if color.isValid():
            btn._color = color.name()
            btn.setStyleSheet(f"background-color: {color.name()};")

    def get_values(self) -> dict:
        """Return validated values from all widgets keyed by field key."""
        result = {}
        for spec in self.FIELD_SPECS:
            key   = spec['key']
            ftype = spec['type']
            _, widget = self._widgets[key]

            if ftype in ('float', 'int'):
                try:
                    result[key] = float(widget.text()) if ftype == 'float' else int(widget.text())
                except ValueError:
                    result[key] = spec.get('default', 0)
            elif ftype == 'bool':
                result[key] = widget.isChecked()
            elif ftype == 'color':
                result[key] = widget._color
            else:
                result[key] = widget.text()
        return result


class GraphWindow(QWidget):
    """
    Separate window for displaying graphs.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena - Data Viewer")
        self.setGeometry(100, 100, 800, 600)

        # Create matplotlib figure and canvas
        self.figure = Figure(figsize=(8, 6))
        self.canvas = FigureCanvas(self.figure)

        # Layout
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def plot_data(self, file_paths: List[str], error_fraction: float = 0.05):
        """
        Plot data from the selected files.

        Args:
            file_paths: List of file paths to plot
            error_fraction: Fraction used to generate uncertainty for 2-column text files
        """
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        for file_path in file_paths:
            try:
                # Load data based on file extension
                path, filename = os.path.split(file_path)
                _, ext = os.path.splitext(filename)

                if ext.lower() in ['.txt', '.dat']:
                    # Read text file
                    data = readTextFile(path, filename, error_fraction=error_fraction)
                else:
                    # Read HDF5 file
                    data = readGenericNXcanSAS(path, filename)

                if data is None:
                    continue

                q = data['Q']
                intensity = data['Intensity']

                # Plot
                label = os.path.basename(file_path)
                ax.plot(q, intensity, 'o-', label=label, markersize=4, linewidth=1.5)

            except Exception as e:
                print(f"Error loading {file_path}: {e}")
                continue

        # Configure plot
        ax.set_xlabel('Q (Å⁻¹)', fontsize=12)
        ax.set_ylabel('Intensity (cm⁻¹)', fontsize=12)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3, which='both')
        ax.legend(fontsize=10)
        ax.set_title('Small-Angle Scattering Data', fontsize=14)

        self.canvas.draw()
        self.show()


class UnifiedFitResultsWindow(QWidget):
    """
    Separate window for displaying Unified Fit results stored in HDF5 files.

    Shows a two-panel matplotlib figure: data + model fit (top) and
    normalised residuals (bottom).  Files that contain no unified fit group
    are silently skipped.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena - Unified Fit Results")
        self.setGeometry(130, 130, 900, 700)

        self.figure = Figure(figsize=(9, 7))
        self.canvas = FigureCanvas(self.figure)

        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def plot_results(self, file_paths: List[str]):
        """
        Load and plot Unified Fit results from the given file paths.

        Only HDF5 files are considered.  Files without a unified fit group
        are skipped silently.
        """
        self.figure.clear()
        ax_main  = self.figure.add_subplot(211)
        ax_resid = self.figure.add_subplot(212)

        found_any = False

        for file_path in file_paths:
            _, ext = os.path.splitext(file_path)
            if ext.lower() not in ['.h5', '.hdf5', '.hdf']:
                continue   # text files cannot carry unified fit results

            try:
                results = load_unified_fit_results(Path(file_path))
            except Exception:
                continue   # no unified fit group or unreadable — skip silently

            label = os.path.basename(file_path)
            Q         = results['Q']
            I_data    = results['intensity_data']
            I_model   = results['intensity_model']
            residuals = results['residuals']
            I_error   = results.get('intensity_error')
            chi2      = results.get('chi_squared', float('nan'))

            # ── data points ────────────────────────────────────────────────
            if I_error is not None:
                pts = ax_main.errorbar(
                    Q, I_data, yerr=I_error,
                    fmt='o', markersize=3,
                    label=f'{label}  (data)', alpha=0.8,
                    capsize=2,
                )
                color = pts[0].get_color()
            else:
                pts, = ax_main.plot(
                    Q, I_data, 'o', markersize=3,
                    label=f'{label}  (data)', alpha=0.8,
                )
                color = pts.get_color()

            # ── model line ─────────────────────────────────────────────────
            ax_main.plot(
                Q, I_model, '-', linewidth=1.5, color=color,
                label=f'{label}  fit  χ²={chi2:.3f}',
            )

            # ── residuals ──────────────────────────────────────────────────
            ax_resid.plot(
                Q, residuals, 'o-', markersize=3, linewidth=1,
                color=color, label=label, alpha=0.8,
            )

            found_any = True

        if found_any:
            ax_main.set_xscale('log')
            ax_main.set_yscale('log')
            ax_main.set_ylabel('Intensity (cm⁻¹)', fontsize=11)
            ax_main.grid(True, alpha=0.3, which='both')
            ax_main.legend(fontsize=9)
            ax_main.set_title('Unified Fit Results', fontsize=13)

            ax_resid.axhline(y=0, color='k', linestyle='-', linewidth=0.8)
            ax_resid.set_xscale('log')
            ax_resid.set_xlabel('Q (Å⁻¹)', fontsize=11)
            ax_resid.set_ylabel('Residuals (normalised)', fontsize=11)
            ax_resid.grid(True, alpha=0.3)
            ax_resid.legend(fontsize=9)
        else:
            ax_main.text(
                0.5, 0.5,
                'No Unified Fit results found\nin selected files',
                ha='center', va='center',
                transform=ax_main.transAxes,
                fontsize=14, color='gray',
            )
            self.figure.delaxes(ax_resid)

        self.figure.tight_layout()
        self.canvas.draw()
        self.show()


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
        self.unified_fit_window = None  # Unified fit panel
        self.sizes_fit_window = None   # Size distribution panel

        # Initialize state manager
        self.state_manager = StateManager()

        # Load last used folder from state
        self.load_last_folder()

        self.init_ui()

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

        # Title
        title_label = QLabel("pyIrena")
        title_label.setStyleSheet("""
            QLabel {
                font-size: 24px;
                font-weight: bold;
                color: #2c3e50;
                padding: 10px;
            }
        """)
        title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        content_layout.addWidget(title_label)

        # Folder selection section
        folder_layout = QHBoxLayout()
        self.folder_button = QPushButton("Select Data Folder")
        self.folder_button.setMinimumHeight(40)
        self.folder_button.clicked.connect(self.select_folder)
        folder_layout.addWidget(self.folder_button)

        self.refresh_button = QPushButton("Refresh")
        self.refresh_button.setMinimumHeight(40)
        self.refresh_button.setMaximumWidth(100)
        self.refresh_button.clicked.connect(self.refresh_file_list)
        self.refresh_button.setEnabled(False)
        folder_layout.addWidget(self.refresh_button)

        self.folder_label = QLabel("No folder selected")
        self.folder_label.setStyleSheet("color: #7f8c8d; font-style: italic;")
        folder_layout.addWidget(self.folder_label)
        folder_layout.addStretch()

        content_layout.addLayout(folder_layout)

        # File type selection
        type_layout = QHBoxLayout()
        type_layout.addWidget(QLabel("File Type:"))

        self.file_type_combo = QComboBox()
        self.file_type_combo.addItem("HDF5 Files (.hdf, .h5, .hdf5)", "hdf5")
        self.file_type_combo.addItem("Text Files (.txt, .dat)", "text")
        self.file_type_combo.addItem("All Supported Files", "all")
        self.file_type_combo.currentIndexChanged.connect(self.refresh_file_list)
        type_layout.addWidget(self.file_type_combo)
        type_layout.addStretch()

        content_layout.addLayout(type_layout)

        # Content area (listbox + graph button)
        file_area_layout = QHBoxLayout()

        # Left side: file list section
        left_layout = QVBoxLayout()

        # File filter input
        filter_layout = QHBoxLayout()
        filter_layout.addWidget(QLabel("Filter:"))
        self.filter_input = QLineEdit()
        self.filter_input.setPlaceholderText("Enter text to filter files...")
        self.filter_input.textChanged.connect(self.filter_files)
        filter_layout.addWidget(self.filter_input)
        left_layout.addLayout(filter_layout)

        # File list
        self.file_list = QListWidget()
        self.file_list.setMinimumWidth(400)  # At least 30 characters wide
        self.file_list.setMinimumHeight(400)  # Show at least 15 items
        self.file_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.file_list.itemDoubleClicked.connect(self.plot_selected_files)
        self.file_list.itemSelectionChanged.connect(self.update_plot_button_state)
        left_layout.addWidget(self.file_list)

        # Configure button — small, sits below the file list
        configure_row = QHBoxLayout()
        self.configure_button = QPushButton("Configure...")
        self.configure_button.setMaximumWidth(110)
        self.configure_button.setMinimumHeight(24)
        self.configure_button.setToolTip("Configure data loading options")
        self.configure_button.clicked.connect(self.open_configure_dialog)
        configure_row.addWidget(self.configure_button)
        configure_row.addStretch()
        left_layout.addLayout(configure_row)

        file_area_layout.addLayout(left_layout, stretch=2)

        # Right side: action buttons — starts at same level as Filter row
        right_layout = QVBoxLayout()
        right_layout.setSpacing(6)

        # ── Graph content checkboxes ───────────────────────────────────────
        graph_content_label = QLabel("Show in graph:")
        graph_content_label.setStyleSheet("font-weight: bold; color: #2c3e50;")
        right_layout.addWidget(graph_content_label)

        cb_row = QHBoxLayout()
        cb_row.setSpacing(10)
        self.data_checkbox = QCheckBox("Data")
        self.data_checkbox.setChecked(True)
        self.data_checkbox.setToolTip("Plot experimental data for selected files")
        cb_row.addWidget(self.data_checkbox)

        self.unified_fit_result_checkbox = QCheckBox("Unified Fit")
        self.unified_fit_result_checkbox.setChecked(False)
        self.unified_fit_result_checkbox.setToolTip(
            "Plot stored Unified Fit results (data + model + residuals).\n"
            "Only HDF5 files are checked; files without fit results are skipped."
        )
        cb_row.addWidget(self.unified_fit_result_checkbox)

        self.size_dist_checkbox = QCheckBox("Size Dist.")
        self.size_dist_checkbox.setChecked(False)
        self.size_dist_checkbox.setToolTip(
            "Open Size Distribution panel with selected data (Create Graph),\n"
            "or include stored sizes results in report (Create Report).\n"
            "Only HDF5 files with stored sizes results are used for reports."
        )
        cb_row.addWidget(self.size_dist_checkbox)
        cb_row.addStretch()
        right_layout.addLayout(cb_row)

        right_layout.addSpacing(4)

        # ── Create Graph / Create Report side by side ──────────────────────
        _btn_style_blue = """
            QPushButton {
                background-color: #3498db; color: white;
                font-size: 12px; font-weight: bold;
                border-radius: 4px; padding: 5px 8px;
            }
            QPushButton:hover { background-color: #2980b9; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """
        _btn_style_purple = """
            QPushButton {
                background-color: #8e44ad; color: white;
                font-size: 12px; font-weight: bold;
                border-radius: 4px; padding: 5px 8px;
            }
            QPushButton:hover { background-color: #7d3c98; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """

        graph_row = QHBoxLayout()
        graph_row.setSpacing(6)

        self.plot_button = QPushButton("Create Graph")
        self.plot_button.setMinimumHeight(34)
        self.plot_button.setStyleSheet(_btn_style_blue)
        self.plot_button.clicked.connect(self.plot_selected_files)
        self.plot_button.setEnabled(False)
        graph_row.addWidget(self.plot_button)

        self.report_button = QPushButton("Create Report")
        self.report_button.setMinimumHeight(34)
        self.report_button.setStyleSheet(_btn_style_purple)
        self.report_button.setToolTip(
            "Generate a Markdown report (.md) summarising the Unified Fit\n"
            "results for each selected HDF5 file.  Files without stored\n"
            "fit results are skipped."
        )
        self.report_button.clicked.connect(self.create_report)
        self.report_button.setEnabled(False)
        graph_row.addWidget(self.report_button)

        right_layout.addLayout(graph_row)

        right_layout.addSpacing(10)

        # ── Unified Fit model button ───────────────────────────────────────
        self.unified_fit_button = QPushButton("Unified Fit")
        self.unified_fit_button.setMinimumHeight(38)
        self.unified_fit_button.setStyleSheet("""
            QPushButton {
                background-color: #27ae60; color: white;
                font-size: 13px; font-weight: bold;
                border-radius: 4px; padding: 6px 10px;
            }
            QPushButton:hover { background-color: #229954; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """)
        self.unified_fit_button.clicked.connect(self.launch_unified_fit)
        self.unified_fit_button.setEnabled(False)
        right_layout.addWidget(self.unified_fit_button)

        # ── Size Distribution model button ─────────────────────────────────
        self.sizes_fit_button = QPushButton("Size Distribution")
        self.sizes_fit_button.setMinimumHeight(38)
        self.sizes_fit_button.setStyleSheet("""
            QPushButton {
                background-color: #2980b9; color: white;
                font-size: 13px; font-weight: bold;
                border-radius: 4px; padding: 6px 10px;
            }
            QPushButton:hover { background-color: #2471a3; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """)
        self.sizes_fit_button.clicked.connect(self.launch_sizes_fit)
        self.sizes_fit_button.setEnabled(False)
        right_layout.addWidget(self.sizes_fit_button)

        right_layout.addStretch()
        file_area_layout.addLayout(right_layout, stretch=1)

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

        # Set minimum window size (at least twice the listbox width)
        self.setMinimumSize(900, 600)

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

        # Add separator and future models placeholder
        models_menu.addSeparator()
        placeholder_action = QAction("More models coming soon...", self)
        placeholder_action.setEnabled(False)
        models_menu.addAction(placeholder_action)

        menu_bar.addMenu(models_menu)

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

            files.sort()
            self.file_list.addItems(files)
            self.status_label.setText(f"Found {len(files)} files")
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

    def update_plot_button_state(self):
        """Enable or disable buttons based on file selection."""
        has_selection = len(self.file_list.selectedItems()) > 0
        self.plot_button.setEnabled(has_selection)
        self.report_button.setEnabled(has_selection)
        self.unified_fit_button.setEnabled(has_selection)
        self.sizes_fit_button.setEnabled(has_selection)

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

        error_fraction = self.state_manager.get('data_selector', 'error_fraction', 0.05)
        plotted = []

        # ── Experimental data ──────────────────────────────────────────────
        if self.data_checkbox.isChecked():
            if self.graph_window is None:
                self.graph_window = GraphWindow()
            try:
                self.graph_window.plot_data(file_paths, error_fraction=error_fraction)
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
                self.unified_fit_results_window.plot_results(file_paths)
                plotted.append("Unified Fit results")
            except Exception as e:
                QMessageBox.critical(
                    self, "Plot Error", f"Error creating Unified Fit plot:\n{str(e)}"
                )

        # ── Size Distribution — open the Sizes panel with the first file ──
        if self.size_dist_checkbox.isChecked():
            self.launch_sizes_fit()
            plotted.append("Size Distribution")

        if plotted:
            self.status_label.setText(
                f"Plotted {len(file_paths)} file(s): {', '.join(plotted)}"
            )
        else:
            self.status_label.setText(
                "Nothing to plot — check 'Data', 'Unified Fit', or 'Size Dist.' checkbox"
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

        show_data  = self.data_checkbox.isChecked()
        show_fit   = self.unified_fit_result_checkbox.isChecked()
        show_sizes = self.size_dist_checkbox.isChecked()

        if not show_data and not show_fit and not show_sizes:
            self.status_label.setText(
                "Nothing to report — check 'Data', 'Unified Fit', or 'Size Dist.' checkbox"
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

            data_info    = None
            fit_results  = None
            sizes_results = None

            # ── Load raw data ──────────────────────────────────────────────
            if show_data:
                try:
                    dir_path, filename = os.path.split(file_path)
                    if ext.lower() in ['.txt', '.dat']:
                        raw = readTextFile(dir_path, filename, error_fraction=error_fraction)
                    else:
                        raw = readGenericNXcanSAS(dir_path, filename)

                    if raw is not None:
                        data_info = {
                            'Q':       raw['Q'],
                            'I':       raw['Intensity'],
                            'I_error': raw.get('Error'),
                        }
                except Exception:
                    pass   # data load failure is non-fatal; section simply omitted

            # ── Load Unified Fit results (HDF5 only) ───────────────────────
            if show_fit and is_hdf:
                try:
                    fit_results = load_unified_fit_results(Path(file_path))
                except Exception:
                    pass   # no fit group or unreadable — section omitted silently

            # ── Load Size Distribution results (HDF5 only) ─────────────────
            if show_sizes and is_hdf:
                try:
                    from pyirena.io.nxcansas_sizes import load_sizes_results
                    sizes_results = load_sizes_results(Path(file_path))
                except Exception:
                    pass   # no sizes group or unreadable — section omitted silently

            # Nothing to write for this file?
            if data_info is None and fit_results is None and sizes_results is None:
                skipped.append(os.path.basename(file_path))
                continue

            md = _build_report(
                file_path,
                data_info=data_info,
                fit_results=fit_results,
                sizes_results=sizes_results,
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
                pass   # Opening is best-effort; don't fail the whole save

        if saved:
            msg = f"Report(s) saved: {', '.join(saved)}"
            if skipped:
                msg += f"  (skipped: {', '.join(skipped)})"
            self.status_label.setText(msg)
        else:
            self.status_label.setText(
                "No reportable content found — no reports generated"
            )

    def launch_unified_fit(self):
        """Launch the Unified Fit model panel with selected data."""
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(
                self,
                "No Selection",
                "Please select one or more files to analyze with Unified Fit."
            )
            return

        # Get full file paths
        file_paths = []
        for item in selected_items:
            file_path = os.path.join(self.current_folder, item.text())
            file_paths.append(file_path)

        # Load first selected file for fitting
        # (Multiple files can be loaded and fitted separately)
        file_path = file_paths[0]
        path, filename = os.path.split(file_path)
        _, ext = os.path.splitext(filename)

        error_fraction = self.state_manager.get('data_selector', 'error_fraction', 0.05)
        try:
            # Load data based on file extension
            if ext.lower() in ['.txt', '.dat']:
                data = readTextFile(path, filename, error_fraction=error_fraction)
                is_nxcansas = False
            else:
                data = readGenericNXcanSAS(path, filename)
                is_nxcansas = True  # HDF5 files loaded with NXcanSAS reader

            if data is None:
                QMessageBox.critical(
                    self,
                    "Load Error",
                    f"Could not load data from {filename}"
                )
                return

            # Create or show unified fit window
            if self.unified_fit_window is None:
                self.unified_fit_window = UnifiedFitPanel()

            # Set the data with filepath and format information
            self.unified_fit_window.set_data(
                data['Q'],
                data['Intensity'],
                data.get('Error'),
                filename,
                filepath=file_path,  # Pass full path to file
                is_nxcansas=is_nxcansas  # Pass format information
            )

            # Show the window
            self.unified_fit_window.show()
            self.unified_fit_window.raise_()
            self.unified_fit_window.activateWindow()

            self.status_label.setText(f"Opened Unified Fit for {filename}")

        except Exception as e:
            QMessageBox.critical(
                self,
                "Error",
                f"Error loading data for Unified Fit:\n{str(e)}"
            )
            self.status_label.setText(f"Error: {str(e)}")

    def launch_sizes_fit(self):
        """Launch the Size Distribution fitting panel with selected data."""
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(
                self,
                "No Selection",
                "Please select one or more files to analyze with Size Distribution."
            )
            return

        file_path = os.path.join(self.current_folder, selected_items[0].text())
        path, filename = os.path.split(file_path)
        _, ext = os.path.splitext(filename)

        error_fraction = self.state_manager.get('data_selector', 'error_fraction', 0.05)
        try:
            if ext.lower() in ['.txt', '.dat']:
                data = readTextFile(path, filename, error_fraction=error_fraction)
                is_nxcansas = False
            else:
                data = readGenericNXcanSAS(path, filename)
                is_nxcansas = True

            if data is None:
                QMessageBox.critical(
                    self,
                    "Load Error",
                    f"Could not load data from {filename}"
                )
                return

            if self.sizes_fit_window is None:
                self.sizes_fit_window = SizesFitPanel()

            self.sizes_fit_window.set_data(
                data['Q'],
                data['Intensity'],
                data.get('Error'),
                filename,
                filepath=file_path,
                is_nxcansas=is_nxcansas,
            )

            self.sizes_fit_window.show()
            self.sizes_fit_window.raise_()
            self.sizes_fit_window.activateWindow()

            self.status_label.setText(f"Opened Size Distribution for {filename}")

        except Exception as e:
            QMessageBox.critical(
                self,
                "Error",
                f"Error loading data for Size Distribution:\n{str(e)}"
            )
            self.status_label.setText(f"Error: {str(e)}")

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
            print(f"Restored last folder: {last_folder}")
        else:
            # Folder doesn't exist or wasn't saved, use home directory
            self.last_folder = str(Path.home())
            self.current_folder = None
            if last_folder:
                print(f"Last folder no longer exists: {last_folder}")
                print(f"Starting in home directory: {self.last_folder}")

    def save_last_folder(self, folder: str):
        """Save the current folder to state."""
        self.state_manager.set('data_selector', 'last_folder', folder)
        self.state_manager.save()

    def show_about(self):
        """Show about dialog."""
        QMessageBox.about(
            self,
            "About pyIrena",
            """<h3>pyIrena</h3>
            <p><b>Python tools for small-angle scattering data analysis</b></p>
            <p>Version: 0.1.0</p>
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


if __name__ == "__main__":
    main()
