"""
Data Selector GUI for pyIrena.

This module provides a GUI panel for selecting data files and displaying
their content as graphs.
"""

import os
import sys
from pathlib import Path
from typing import List, Optional

try:
    from PySide6.QtWidgets import (
        QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
        QListWidget, QLabel, QLineEdit, QFileDialog, QComboBox,
        QAbstractItemView, QMessageBox, QMenuBar, QMenu
    )
    from PySide6.QtCore import Qt, QDir
    from PySide6.QtGui import QAction
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
            QListWidget, QLabel, QLineEdit, QFileDialog, QComboBox,
            QAbstractItemView, QMessageBox, QMenuBar, QMenu
        )
        from PyQt6.QtCore import Qt, QDir
        from PyQt6.QtGui import QAction
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
from pyirena.gui.unified_fit import UnifiedFitPanel
from pyirena.state import StateManager


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

    def plot_data(self, file_paths: List[str]):
        """
        Plot data from the selected files.

        Args:
            file_paths: List of file paths to plot
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
                    data = readTextFile(path, filename)
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
        self.unified_fit_window = None  # Unified fit panel

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

        file_area_layout.addLayout(left_layout, stretch=2)

        # Right side: action buttons
        right_layout = QVBoxLayout()
        right_layout.addStretch()

        self.plot_button = QPushButton("Create Graph")
        self.plot_button.setMinimumWidth(150)
        self.plot_button.setMinimumHeight(50)
        self.plot_button.setStyleSheet("""
            QPushButton {
                background-color: #3498db;
                color: white;
                font-size: 14px;
                font-weight: bold;
                border-radius: 5px;
                padding: 10px;
            }
            QPushButton:hover {
                background-color: #2980b9;
            }
            QPushButton:disabled {
                background-color: #bdc3c7;
            }
        """)
        self.plot_button.clicked.connect(self.plot_selected_files)
        self.plot_button.setEnabled(False)
        right_layout.addWidget(self.plot_button)

        right_layout.addSpacing(20)

        # Unified Fit button
        self.unified_fit_button = QPushButton("Unified Fit")
        self.unified_fit_button.setMinimumWidth(150)
        self.unified_fit_button.setMinimumHeight(50)
        self.unified_fit_button.setStyleSheet("""
            QPushButton {
                background-color: #27ae60;
                color: white;
                font-size: 14px;
                font-weight: bold;
                border-radius: 5px;
                padding: 10px;
            }
            QPushButton:hover {
                background-color: #229954;
            }
            QPushButton:disabled {
                background-color: #bdc3c7;
            }
        """)
        self.unified_fit_button.clicked.connect(self.launch_unified_fit)
        self.unified_fit_button.setEnabled(False)
        right_layout.addWidget(self.unified_fit_button)

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
        """Enable or disable the plot and unified fit buttons based on file selection."""
        has_selection = len(self.file_list.selectedItems()) > 0
        self.plot_button.setEnabled(has_selection)
        self.unified_fit_button.setEnabled(has_selection)

    def plot_selected_files(self):
        """Plot the selected files."""
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(
                self,
                "No Selection",
                "Please select one or more files to plot."
            )
            return

        # Get full file paths
        file_paths = []
        for item in selected_items:
            file_path = os.path.join(self.current_folder, item.text())
            file_paths.append(file_path)

        # Create or update graph window
        if self.graph_window is None:
            self.graph_window = GraphWindow()

        try:
            self.graph_window.plot_data(file_paths)
            self.status_label.setText(f"Plotted {len(file_paths)} file(s)")
        except Exception as e:
            QMessageBox.critical(
                self,
                "Plot Error",
                f"Error creating plot:\n{str(e)}"
            )
            self.status_label.setText(f"Error: {str(e)}")

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

        try:
            # Load data based on file extension
            if ext.lower() in ['.txt', '.dat']:
                data = readTextFile(path, filename)
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
