"""
Unified Fit Model GUI for pyIrena.

This module provides a GUI panel for fitting small-angle scattering data
using the Unified Fit model (Beaucage model).
"""

import os
import sys
from typing import List, Optional, Dict
import numpy as np

try:
    from PySide6.QtWidgets import (
        QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
        QLabel, QLineEdit, QCheckBox, QSpinBox, QTabWidget, QGroupBox,
        QGridLayout, QMessageBox, QSplitter
    )
    from PySide6.QtCore import Qt, Signal
    from PySide6.QtGui import QFont, QDoubleValidator
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
            QLabel, QLineEdit, QCheckBox, QSpinBox, QTabWidget, QGroupBox,
            QGridLayout, QMessageBox, QSplitter
        )
        from PyQt6.QtCore import Qt, pyqtSignal as Signal
        from PyQt6.QtGui import QFont, QDoubleValidator
    except ImportError:
        raise ImportError(
            "Neither PySide6 nor PyQt6 found. Install with: pip install PySide6"
        )

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from pyirena.core.unified import UnifiedFitModel, UnifiedLevel


class UnifiedFitGraphWindow(QWidget):
    """
    Graph window for displaying data and unified fit results.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena - Unified Fit")
        self.setGeometry(100, 100, 900, 700)

        # Create matplotlib figure with two subplots
        self.figure = Figure(figsize=(9, 7))
        self.canvas = FigureCanvas(self.figure)

        # Add navigation toolbar
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Cursor attributes
        self.cursor_left = None
        self.cursor_right = None
        self.cursor_left_line = None
        self.cursor_right_line = None
        self.dragging_cursor = None
        self.q_data = None  # Store Q data for cursor snapping

        # Layout
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)

        # Initialize plots
        self.init_plots()

        # Connect mouse events for cursor dragging
        self.canvas.mpl_connect('button_press_event', self.on_mouse_press)
        self.canvas.mpl_connect('button_release_event', self.on_mouse_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)

    def init_plots(self):
        """Initialize the plot axes."""
        self.figure.clear()

        # Main plot (data + fit)
        self.ax_main = self.figure.add_subplot(211)
        self.ax_main.set_xlabel('Q (Å⁻¹)', fontsize=11)
        self.ax_main.set_ylabel('Intensity (cm⁻¹)', fontsize=11)
        self.ax_main.set_xscale('log')
        self.ax_main.set_yscale('log')
        self.ax_main.grid(True, alpha=0.3, which='both')
        self.ax_main.set_title('Unified Fit Model', fontsize=12, fontweight='bold')

        # Residuals plot
        self.ax_residuals = self.figure.add_subplot(212)
        self.ax_residuals.set_xlabel('Q (Å⁻¹)', fontsize=11)
        self.ax_residuals.set_ylabel('Residuals', fontsize=11)
        self.ax_residuals.set_xscale('log')
        self.ax_residuals.axhline(y=0, color='k', linestyle='--', alpha=0.5)
        self.ax_residuals.grid(True, alpha=0.3, which='both')

        self.figure.tight_layout()
        self.canvas.draw()

    def plot_data(self, q, intensity, error=None, label='Data'):
        """Plot experimental data."""
        self.q_data = q  # Store Q data for cursor positioning
        self.ax_main.errorbar(
            q, intensity, yerr=error,
            fmt='o', markersize=4, capsize=2,
            label=label, alpha=0.7
        )
        self.ax_main.legend(fontsize=9)

        # Initialize cursor positions if not already set
        if self.cursor_left is None and len(q) > 0:
            # Set left cursor to 20% of Q range
            q_min, q_max = q.min(), q.max()
            self.cursor_left = q_min * (q_max / q_min) ** 0.2
            self.cursor_right = q_min * (q_max / q_min) ** 0.8

        # Always add cursors if positions exist (to restore after init_plots clears the figure)
        if self.cursor_left is not None and self.cursor_right is not None:
            self.add_cursors()

        self.canvas.draw()

    def plot_fit(self, q, fit, label='Unified Fit'):
        """Plot fit curve."""
        self.ax_main.plot(q, fit, '-', linewidth=2, label=label)
        self.ax_main.legend(fontsize=9)
        self.canvas.draw()

    def plot_residuals(self, q, residuals):
        """Plot fit residuals."""
        self.ax_residuals.clear()
        self.ax_residuals.plot(q, residuals, 'o', markersize=3, alpha=0.6)
        self.ax_residuals.axhline(y=0, color='k', linestyle='--', alpha=0.5)
        self.ax_residuals.set_xlabel('Q (Å⁻¹)', fontsize=11)
        self.ax_residuals.set_ylabel('Residuals', fontsize=11)
        self.ax_residuals.set_xscale('log')
        self.ax_residuals.grid(True, alpha=0.3, which='both')
        self.canvas.draw()

    def add_cursors(self):
        """Add draggable vertical cursor lines to the main plot."""
        # Remove old cursor lines if they exist
        # Use broad exception handling since matplotlib can raise various exceptions
        # when trying to remove artists that were already cleared
        try:
            if self.cursor_left_line is not None:
                self.cursor_left_line.remove()
        except Exception:
            # Line was already removed (e.g., by figure.clear()) or can't be removed
            pass

        try:
            if self.cursor_right_line is not None:
                self.cursor_right_line.remove()
        except Exception:
            # Line was already removed (e.g., by figure.clear()) or can't be removed
            pass

        # Add vertical lines at cursor positions
        self.cursor_left_line = self.ax_main.axvline(
            x=self.cursor_left,
            color='red',
            linestyle='--',
            linewidth=2,
            label='Left Cursor',
            alpha=0.8
        )
        self.cursor_right_line = self.ax_main.axvline(
            x=self.cursor_right,
            color='blue',
            linestyle='--',
            linewidth=2,
            label='Right Cursor',
            alpha=0.8
        )

        self.ax_main.legend(fontsize=9)
        self.canvas.draw()

    def on_mouse_press(self, event):
        """Handle mouse press for cursor dragging."""
        if event.inaxes != self.ax_main:
            return

        if self.cursor_left is None or self.cursor_right is None:
            return

        # Check if click is near a cursor (within 10% of the log distance)
        if event.xdata is not None:
            log_x = np.log10(event.xdata)
            log_left = np.log10(self.cursor_left)
            log_right = np.log10(self.cursor_right)

            # Calculate tolerance based on axis range
            log_xlim = np.log10(self.ax_main.get_xlim())
            tolerance = 0.05 * (log_xlim[1] - log_xlim[0])

            if abs(log_x - log_left) < tolerance:
                self.dragging_cursor = 'left'
            elif abs(log_x - log_right) < tolerance:
                self.dragging_cursor = 'right'

    def on_mouse_release(self, event):
        """Handle mouse release for cursor dragging."""
        self.dragging_cursor = None

    def on_mouse_move(self, event):
        """Handle mouse move for cursor dragging."""
        if self.dragging_cursor is None:
            return

        if event.inaxes != self.ax_main or event.xdata is None:
            return

        # Update cursor position
        new_x = event.xdata

        # Ensure cursors don't cross
        if self.dragging_cursor == 'left':
            if new_x < self.cursor_right:
                self.cursor_left = new_x
                self.cursor_left_line.set_xdata([new_x, new_x])
        elif self.dragging_cursor == 'right':
            if new_x > self.cursor_left:
                self.cursor_right = new_x
                self.cursor_right_line.set_xdata([new_x, new_x])

        self.canvas.draw()

    def get_cursor_range(self):
        """Get the current cursor Q range."""
        if self.cursor_left is not None and self.cursor_right is not None:
            return (self.cursor_left, self.cursor_right)
        return None


class LevelParametersWidget(QWidget):
    """
    Widget for a single level's parameters in the Unified Fit model.
    """

    def __init__(self, level_number: int, parent=None):
        super().__init__(parent)
        self.level_number = level_number
        self.init_ui()

    def init_ui(self):
        """Initialize the user interface for this level."""
        layout = QVBoxLayout()
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)

        # Level header
        header = QLabel(f"Level {self.level_number}")
        header.setStyleSheet("""
            QLabel {
                background-color: #e74c3c;
                color: white;
                font-weight: bold;
                font-size: 12px;
                padding: 5px;
            }
        """)
        layout.addWidget(header)

        # Controls header
        controls_header = QLabel("Controls")
        controls_header.setStyleSheet("""
            QLabel {
                background-color: #e74c3c;
                color: white;
                font-weight: bold;
                font-size: 11px;
                padding: 3px;
            }
        """)
        layout.addWidget(controls_header)

        # Parameters grid
        grid = QGridLayout()
        grid.setSpacing(8)

        # Column headers
        grid.addWidget(QLabel(""), 0, 0)
        grid.addWidget(QLabel("Fit?"), 0, 2, Qt.AlignmentFlag.AlignCenter)
        grid.addWidget(QLabel("Low limit:"), 0, 3, Qt.AlignmentFlag.AlignRight)
        grid.addWidget(QLabel("High Limit:"), 0, 4, Qt.AlignmentFlag.AlignRight)

        # G parameter
        row = 1
        grid.addWidget(QLabel("G"), row, 0)
        self.g_value = QLineEdit("100")
        self.g_value.setValidator(QDoubleValidator())
        self.g_value.setMinimumWidth(120)
        self.g_value.editingFinished.connect(self.on_parameter_changed)
        grid.addWidget(self.g_value, row, 1)
        self.g_fit = QCheckBox()
        grid.addWidget(self.g_fit, row, 2, Qt.AlignmentFlag.AlignCenter)
        self.g_low = QLineEdit()
        self.g_low.setMaximumWidth(80)
        grid.addWidget(self.g_low, row, 3)
        self.g_high = QLineEdit()
        self.g_high.setMaximumWidth(80)
        grid.addWidget(self.g_high, row, 4)

        # Rg parameter
        row = 2
        grid.addWidget(QLabel("Rg"), row, 0)
        self.rg_value = QLineEdit("100")
        self.rg_value.setValidator(QDoubleValidator())
        self.rg_value.editingFinished.connect(self.on_parameter_changed)
        grid.addWidget(self.rg_value, row, 1)
        self.rg_fit = QCheckBox()
        grid.addWidget(self.rg_fit, row, 2, Qt.AlignmentFlag.AlignCenter)
        self.rg_low = QLineEdit()
        self.rg_low.setMaximumWidth(80)
        grid.addWidget(self.rg_low, row, 3)
        self.rg_high = QLineEdit()
        self.rg_high.setMaximumWidth(80)
        grid.addWidget(self.rg_high, row, 4)

        layout.addLayout(grid)

        # Fit Rg/G button
        self.fit_rg_g_button = QPushButton("Fit Rg/G btwn cursors")
        self.fit_rg_g_button.setMinimumHeight(30)
        layout.addWidget(self.fit_rg_g_button)

        # Estimate B checkbox
        self.estimate_b_check = QCheckBox("Estimate B from G/Rg/P?")
        layout.addWidget(self.estimate_b_check)

        # B and P parameters
        grid2 = QGridLayout()
        grid2.setSpacing(8)

        # B parameter
        row = 0
        grid2.addWidget(QLabel("B"), row, 0)
        self.b_value = QLineEdit("0.01")
        self.b_value.setValidator(QDoubleValidator())
        self.b_value.setMinimumWidth(120)
        self.b_value.editingFinished.connect(self.on_parameter_changed)
        grid2.addWidget(self.b_value, row, 1)
        self.b_fit = QCheckBox()
        grid2.addWidget(self.b_fit, row, 2, Qt.AlignmentFlag.AlignCenter)
        self.b_low = QLineEdit()
        self.b_low.setMaximumWidth(80)
        grid2.addWidget(self.b_low, row, 3)
        self.b_high = QLineEdit()
        self.b_high.setMaximumWidth(80)
        grid2.addWidget(self.b_high, row, 4)

        # P parameter
        row = 1
        grid2.addWidget(QLabel("P"), row, 0)
        self.p_value = QLineEdit("4")
        self.p_value.setValidator(QDoubleValidator())
        self.p_value.editingFinished.connect(self.on_parameter_changed)
        grid2.addWidget(self.p_value, row, 1)
        self.p_fit = QCheckBox()
        grid2.addWidget(self.p_fit, row, 2, Qt.AlignmentFlag.AlignCenter)
        self.p_low = QLineEdit()
        self.p_low.setMaximumWidth(80)
        grid2.addWidget(self.p_low, row, 3)
        self.p_high = QLineEdit()
        self.p_high.setMaximumWidth(80)
        grid2.addWidget(self.p_high, row, 4)

        layout.addLayout(grid2)

        # Warning label
        self.warning_label = QLabel("Level may not be physically feasible")
        self.warning_label.setStyleSheet("""
            QLabel {
                color: #e74c3c;
                font-style: italic;
                font-size: 10px;
            }
        """)
        self.warning_label.setVisible(False)
        layout.addWidget(self.warning_label)

        # Fit P/B button
        self.fit_p_b_button = QPushButton("Fit P/B btwn cursors")
        self.fit_p_b_button.setMinimumHeight(30)
        layout.addWidget(self.fit_p_b_button)

        # RgCutoff
        rg_cutoff_layout = QHBoxLayout()
        rg_cutoff_layout.addWidget(QLabel("RgCutoff"))
        self.rg_cutoff = QLineEdit("0")
        self.rg_cutoff.setValidator(QDoubleValidator())
        self.rg_cutoff.setMaximumWidth(100)
        rg_cutoff_layout.addWidget(self.rg_cutoff)
        rg_cutoff_layout.addStretch()
        layout.addLayout(rg_cutoff_layout)

        # pi B/Q display
        pi_bq_layout = QHBoxLayout()
        pi_bq_layout.addWidget(QLabel("pi B/Q [m2/cm3]"))
        self.pi_bq_value = QLineEdit("0")
        self.pi_bq_value.setReadOnly(True)
        self.pi_bq_value.setMaximumWidth(100)
        pi_bq_layout.addWidget(self.pi_bq_value)
        pi_bq_layout.addStretch()
        layout.addLayout(pi_bq_layout)

        # Correlated system checkbox
        self.correlated_check = QCheckBox("Is this correlated system?")
        layout.addWidget(self.correlated_check)

        # Copy/Move/swap button
        self.copy_move_button = QPushButton("Copy/Move/swap level")
        self.copy_move_button.setMinimumHeight(30)
        layout.addWidget(self.copy_move_button)

        layout.addStretch()
        self.setLayout(layout)

    def get_parameters(self) -> Dict:
        """Get all parameters for this level."""
        return {
            'G': float(self.g_value.text() or 0),
            'Rg': float(self.rg_value.text() or 0),
            'B': float(self.b_value.text() or 0),
            'P': float(self.p_value.text() or 0),
            'RgCutoff': float(self.rg_cutoff.text() or 0),
            'fit_G': self.g_fit.isChecked(),
            'fit_Rg': self.rg_fit.isChecked(),
            'fit_B': self.b_fit.isChecked(),
            'fit_P': self.p_fit.isChecked(),
            'estimate_B': self.estimate_b_check.isChecked(),
            'correlated': self.correlated_check.isChecked(),
        }

    def set_parameters(self, params: Dict):
        """Set parameters for this level."""
        if 'G' in params:
            self.g_value.setText(str(params['G']))
        if 'Rg' in params:
            self.rg_value.setText(str(params['Rg']))
        if 'B' in params:
            self.b_value.setText(str(params['B']))
        if 'P' in params:
            self.p_value.setText(str(params['P']))

    def on_parameter_changed(self):
        """Called when a parameter value changes."""
        # Emit signal to parent (UnifiedFitPanel) to trigger auto-update
        parent = self.parent()
        while parent is not None:
            if isinstance(parent, UnifiedFitPanel):
                parent.on_parameter_changed()
                break
            parent = parent.parent()


class UnifiedFitPanel(QWidget):
    """
    Main Unified Fit panel for pyIrena.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.graph_window = None
        self.data = None  # Will store current data
        self.model = UnifiedFitModel()  # Unified fit model
        self.fit_result = None  # Store last fit result
        self.init_ui()

    def init_ui(self):
        """Initialize the user interface."""
        self.setWindowTitle("pyIrena - Unified Fit Model")

        # Main horizontal splitter (panel on left, graph on right)
        main_splitter = QSplitter(Qt.Orientation.Horizontal)

        # Left panel (controls)
        left_panel = self.create_control_panel()
        main_splitter.addWidget(left_panel)

        # Right panel (graph)
        self.graph_window = UnifiedFitGraphWindow()
        main_splitter.addWidget(self.graph_window)

        # Set initial sizes (1:2 ratio)
        main_splitter.setSizes([400, 800])

        # Main layout
        main_layout = QVBoxLayout()
        main_layout.addWidget(main_splitter)
        self.setLayout(main_layout)

        # Set minimum window size
        self.setMinimumSize(1200, 700)

    def create_control_panel(self) -> QWidget:
        """Create the left control panel."""
        panel = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)

        # Title
        title_label = QLabel("Unified model input")
        title_label.setStyleSheet("""
            QLabel {
                font-size: 14px;
                font-weight: bold;
                color: #2c3e50;
                background-color: #ecf0f1;
                padding: 8px;
                border: 1px solid #bdc3c7;
            }
        """)
        layout.addWidget(title_label)

        # Top controls row
        top_controls = QHBoxLayout()

        self.graph_unified_button = QPushButton("Graph Unified")
        self.graph_unified_button.setMinimumHeight(35)
        self.graph_unified_button.setStyleSheet("""
            QPushButton {
                background-color: #27ae60;
                color: white;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #229954;
            }
        """)
        self.graph_unified_button.clicked.connect(self.graph_unified)
        top_controls.addWidget(self.graph_unified_button)

        top_controls.addStretch()

        top_controls.addWidget(QLabel("Number of levels:"))
        self.num_levels_spin = QSpinBox()
        self.num_levels_spin.setMinimum(1)
        self.num_levels_spin.setMaximum(5)
        self.num_levels_spin.setValue(1)
        self.num_levels_spin.setMinimumHeight(30)
        self.num_levels_spin.valueChanged.connect(self.on_num_levels_changed)
        top_controls.addWidget(self.num_levels_spin)

        layout.addLayout(top_controls)

        # Checkboxes row
        checkboxes_layout = QVBoxLayout()
        self.update_auto_check = QCheckBox("Update automatically?")
        checkboxes_layout.addWidget(self.update_auto_check)

        self.display_local_check = QCheckBox("Display local (Porod & Guinier) fits?")
        checkboxes_layout.addWidget(self.display_local_check)

        self.no_limits_check = QCheckBox("No limits?")
        checkboxes_layout.addWidget(self.no_limits_check)

        layout.addLayout(checkboxes_layout)

        # Level tabs
        self.level_tabs = QTabWidget()
        self.level_widgets = []

        # Create widgets for all 5 possible levels
        for i in range(1, 6):
            level_widget = LevelParametersWidget(i)
            self.level_widgets.append(level_widget)
            self.level_tabs.addTab(level_widget, f"{i}. Level")

        # Initially disable unused levels
        self.update_level_tabs()

        layout.addWidget(self.level_tabs)

        # SAS Background
        background_layout = QHBoxLayout()
        background_layout.addWidget(QLabel("SAS Background"))
        self.background_value = QLineEdit("1e-06")
        self.background_value.setValidator(QDoubleValidator())
        self.background_value.setMaximumWidth(100)
        self.background_value.editingFinished.connect(self.on_parameter_changed)
        background_layout.addWidget(self.background_value)
        self.fit_background_check = QCheckBox("Fit Bckg?")
        background_layout.addWidget(self.fit_background_check)
        self.skip_fit_check = QCheckBox("Skip Fit Check?")
        background_layout.addWidget(self.skip_fit_check)
        background_layout.addStretch()
        layout.addLayout(background_layout)

        # Fitting method label
        fit_method_label = QLabel("Fit using least square fitting ?")
        fit_method_label.setStyleSheet("color: #3498db; font-style: italic;")
        layout.addWidget(fit_method_label)

        # Fit buttons row
        fit_buttons = QHBoxLayout()
        self.fit_button = QPushButton("Fit")
        self.fit_button.setMinimumHeight(35)
        self.fit_button.setStyleSheet("""
            QPushButton {
                background-color: #27ae60;
                color: white;
                font-weight: bold;
                font-size: 13px;
            }
            QPushButton:hover {
                background-color: #229954;
            }
        """)
        self.fit_button.clicked.connect(self.run_fit)
        fit_buttons.addWidget(self.fit_button)

        self.revert_button = QPushButton("Revert back")
        self.revert_button.setMinimumHeight(35)
        fit_buttons.addWidget(self.revert_button)

        layout.addLayout(fit_buttons)

        # Additional buttons row
        additional_buttons = QHBoxLayout()
        self.reset_unif_button = QPushButton("reset unif?")
        self.reset_unif_button.setMinimumHeight(30)
        additional_buttons.addWidget(self.reset_unif_button)

        self.fix_limits_button = QPushButton("Fix limits?")
        self.fix_limits_button.setMinimumHeight(30)
        additional_buttons.addWidget(self.fix_limits_button)

        self.store_local_check = QCheckBox("Store local (Porod & Guinier) fits?")
        additional_buttons.addWidget(self.store_local_check)

        layout.addLayout(additional_buttons)

        # Results section header
        results_header = QLabel("Results")
        results_header.setStyleSheet("""
            QLabel {
                font-weight: bold;
                color: #3498db;
                font-size: 12px;
                margin-top: 5px;
            }
        """)
        layout.addWidget(results_header)

        # Results buttons row 1
        results_buttons1 = QHBoxLayout()
        self.store_data_button = QPushButton("Store in Data Folder")
        self.store_data_button.setMinimumHeight(30)
        results_buttons1.addWidget(self.store_data_button)

        self.export_ascii_button = QPushButton("Export ASCII")
        self.export_ascii_button.setMinimumHeight(30)
        results_buttons1.addWidget(self.export_ascii_button)

        self.results_graphs_button = QPushButton("Results to graphs")
        self.results_graphs_button.setMinimumHeight(30)
        results_buttons1.addWidget(self.results_graphs_button)

        layout.addLayout(results_buttons1)

        # Results buttons row 2
        results_buttons2 = QHBoxLayout()
        self.analyze_results_button = QPushButton("Analyze Results")
        self.analyze_results_button.setMinimumHeight(30)
        results_buttons2.addWidget(self.analyze_results_button)

        self.anal_uncertainty_button = QPushButton("Anal. Uncertainty")
        self.anal_uncertainty_button.setMinimumHeight(30)
        results_buttons2.addWidget(self.anal_uncertainty_button)

        self.ext_warnings_check = QCheckBox("Ext. warnings?")
        results_buttons2.addWidget(self.ext_warnings_check)

        layout.addLayout(results_buttons2)

        # Status label
        self.status_label = QLabel("Ready - Load data to begin")
        self.status_label.setStyleSheet("""
            QLabel {
                color: #7f8c8d;
                padding: 5px;
                border-top: 1px solid #bdc3c7;
                font-size: 10px;
            }
        """)
        layout.addWidget(self.status_label)

        panel.setLayout(layout)
        return panel

    def on_num_levels_changed(self, value):
        """Handle change in number of levels."""
        self.update_level_tabs()

    def update_level_tabs(self):
        """Update which level tabs are enabled."""
        num_levels = self.num_levels_spin.value()
        for i in range(5):
            self.level_tabs.setTabEnabled(i, i < num_levels)

    def on_parameter_changed(self):
        """Called when any parameter changes. Auto-update if enabled."""
        if self.update_auto_check.isChecked() and self.data is not None:
            # Automatically recalculate and graph
            self.graph_unified()

    def set_data(self, q, intensity, error=None, label='Data'):
        """Set the data to be fitted."""
        self.data = {
            'Q': q,
            'Intensity': intensity,
            'Error': error,
            'label': label
        }
        # Plot the data
        if self.graph_window:
            self.graph_window.init_plots()
            self.graph_window.plot_data(q, intensity, error, label)

        # Update status
        if hasattr(self, 'status_label'):
            self.status_label.setText(f"Loaded: {label} ({len(q)} points)")

    def graph_unified(self):
        """Graph the unified fit with current parameters (no fitting)."""
        if self.data is None:
            QMessageBox.warning(
                self,
                "No Data",
                "Please load data first from the main panel."
            )
            return

        try:
            # Get current parameters from GUI
            num_levels = self.num_levels_spin.value()
            levels = []

            for i in range(num_levels):
                params = self.level_widgets[i].get_parameters()

                # Create UnifiedLevel object
                level = UnifiedLevel(
                    Rg=params['Rg'],
                    G=params['G'],
                    P=params['P'],
                    B=params['B'],
                    RgCO=params['RgCutoff'],
                    correlations=params['correlated']
                )
                levels.append(level)

            # Get background
            background = float(self.background_value.text() or 0)

            # Update model
            self.model.num_levels = num_levels
            self.model.levels = levels
            self.model.background = background

            # Calculate model intensity for all data
            q_calc = self.data['Q']
            intensity_calc = self.model.calculate_intensity(q_calc)

            # Plot the data and fit
            self.graph_window.init_plots()
            self.graph_window.plot_data(
                self.data['Q'],
                self.data['Intensity'],
                self.data.get('Error'),
                self.data['label']
            )
            self.graph_window.plot_fit(q_calc, intensity_calc, 'Unified Fit')

            # Calculate and plot residuals
            if self.data.get('Error') is not None:
                residuals = (self.data['Intensity'] - intensity_calc) / self.data['Error']
            else:
                residuals = (self.data['Intensity'] - intensity_calc) / self.data['Intensity']

            self.graph_window.plot_residuals(q_calc, residuals)

            self.status_label.setText(f"Calculated unified fit with {num_levels} level(s)")

        except Exception as e:
            QMessageBox.critical(
                self,
                "Calculation Error",
                f"Error calculating unified fit:\n{str(e)}"
            )
            print(f"Error in graph_unified: {e}")
            import traceback
            traceback.print_exc()

    def run_fit(self):
        """Run the unified fit."""
        if self.data is None:
            QMessageBox.warning(
                self,
                "No Data",
                "Please load data first from the main panel."
            )
            return

        try:
            # Get current parameters from GUI
            num_levels = self.num_levels_spin.value()
            levels = []

            for i in range(num_levels):
                params = self.level_widgets[i].get_parameters()

                # Create UnifiedLevel object
                level = UnifiedLevel(
                    Rg=params['Rg'],
                    G=params['G'],
                    P=params['P'],
                    B=params['B'],
                    RgCO=params['RgCutoff'],
                    correlations=params['correlated'],
                    # Set which parameters to fit
                    fit_Rg=params['fit_Rg'],
                    fit_G=params['fit_G'],
                    fit_P=params['fit_P'],
                    fit_B=params['fit_B']
                )
                levels.append(level)

            # Get background
            background = float(self.background_value.text() or 0)
            fit_background = self.fit_background_check.isChecked()

            # Update model
            self.model.num_levels = num_levels
            self.model.levels = levels
            self.model.background = background
            self.model.fit_background = fit_background  # Set fit_background flag

            # Get cursor range and filter data
            cursor_range = self.graph_window.get_cursor_range()
            if cursor_range is not None:
                q_min, q_max = cursor_range
                # Filter data to only include points between cursors
                mask = (self.data['Q'] >= q_min) & (self.data['Q'] <= q_max)
                q_fit = self.data['Q'][mask]
                intensity_fit = self.data['Intensity'][mask]
                error_fit = self.data.get('Error')[mask] if self.data.get('Error') is not None else None
                num_points = len(q_fit)
            else:
                # Use all data if no cursors
                q_fit = self.data['Q']
                intensity_fit = self.data['Intensity']
                error_fit = self.data.get('Error')
                num_points = len(q_fit)

            # Run the fit with filtered data
            self.status_label.setText(f"Fitting {num_levels} level(s) with {num_points} points... please wait")

            result = self.model.fit(
                q_fit,
                intensity_fit,
                error_fit
            )

            self.fit_result = result

            # Update GUI with fitted parameters
            for i in range(num_levels):
                fitted_level = result['levels'][i]  # Dictionary access, not attribute
                self.level_widgets[i].set_parameters({
                    'G': fitted_level.G,
                    'Rg': fitted_level.Rg,
                    'B': fitted_level.B,
                    'P': fitted_level.P
                })

            # Update background
            self.background_value.setText(f"{result['background']:.6e}")

            # Calculate model intensity with fitted parameters (for all data)
            intensity_calc = self.model.calculate_intensity(self.data['Q'])

            # Plot the data and fit
            self.graph_window.init_plots()
            self.graph_window.plot_data(
                self.data['Q'],
                self.data['Intensity'],
                self.data.get('Error'),
                self.data['label']
            )
            self.graph_window.plot_fit(self.data['Q'], intensity_calc, 'Fitted Model')

            # Calculate and plot residuals (for all data)
            if self.data.get('Error') is not None:
                residuals = (self.data['Intensity'] - intensity_calc) / self.data['Error']
            else:
                residuals = (self.data['Intensity'] - intensity_calc) / self.data['Intensity']

            self.graph_window.plot_residuals(self.data['Q'], residuals)

            # Show fit statistics
            chi2 = result.get('chi_squared', 0.0)  # Dictionary access
            reduced_chi2 = result.get('reduced_chi_squared', 0.0)

            # Add cursor range info to status
            if cursor_range is not None:
                self.status_label.setText(
                    f"Fit complete: {num_levels} level(s), {num_points} pts (Q: {q_min:.3e} - {q_max:.3e}), χ² = {chi2:.4f}"
                )
            else:
                self.status_label.setText(
                    f"Fit complete: {num_levels} level(s), χ² = {chi2:.4f}"
                )

            # Include cursor range in success message
            cursor_info = f"\nFit range: Q = {q_min:.3e} to {q_max:.3e}\nPoints used: {num_points}" if cursor_range else ""

            QMessageBox.information(
                self,
                "Fit Complete",
                f"Fit completed successfully!\n\n"
                f"Number of levels: {num_levels}\n"
                f"Chi-squared: {chi2:.4f}\n"
                f"Reduced χ²: {reduced_chi2:.4f}\n"
                f"Parameters updated in GUI.{cursor_info}"
            )

        except Exception as e:
            QMessageBox.critical(
                self,
                "Fit Error",
                f"Error during fitting:\n{str(e)}"
            )
            self.status_label.setText("Fit failed - see error message")
            print(f"Error in run_fit: {e}")
            import traceback
            traceback.print_exc()


def main():
    """Main entry point for the unified fit GUI."""
    app = QApplication(sys.argv)

    # Set application style
    app.setStyle('Fusion')

    # Create and show the main window
    window = UnifiedFitPanel()

    # For testing, add some dummy data
    q = np.logspace(-3, 0, 100)
    intensity = 1e10 * np.exp(-(q * 100)**2 / 3) + 1e6 * q**-4
    error = 0.05 * intensity
    window.set_data(q, intensity, error, "Test Data")

    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
