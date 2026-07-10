"""
pyirena.gui.data_selector._qt — single PySide6/PyQt6 import point for the
data_selector package (PySide6 preferred, PyQt6 fallback).
"""

try:
    from PySide6.QtWidgets import (
        QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QPushButton,
        QListWidget, QLabel, QLineEdit, QFileDialog, QComboBox,
        QAbstractItemView, QMessageBox, QMenuBar, QMenu, QFrame, QScrollArea,
        QDialog, QFormLayout, QDialogButtonBox, QGroupBox, QCheckBox, QColorDialog,
        QTableWidget, QTableWidgetItem, QInputDialog,
    )
    from PySide6.QtCore import Qt, QDir, QThread, Signal, QUrl
    from PySide6.QtGui import QAction, QDoubleValidator, QDesktopServices
except ImportError:
    try:
        from PyQt6.QtWidgets import (  # type: ignore[no-redef]
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QPushButton,
        QListWidget, QLabel, QLineEdit, QFileDialog, QComboBox,
        QAbstractItemView, QMessageBox, QMenuBar, QMenu, QFrame, QScrollArea,
        QDialog, QFormLayout, QDialogButtonBox, QGroupBox, QCheckBox, QColorDialog,
        QTableWidget, QTableWidgetItem, QInputDialog,
        )
        from PyQt6.QtCore import Qt, QDir, QThread, pyqtSignal as Signal, QUrl  # type: ignore[no-redef]
        from PyQt6.QtGui import QAction, QDoubleValidator, QDesktopServices  # type: ignore[no-redef]
    except ImportError:
        raise ImportError(
            "Neither PySide6 nor PyQt6 found. Install with: pip install PySide6"
        )

__all__ = [
    "QApplication", "QWidget", "QVBoxLayout", "QHBoxLayout", "QGridLayout",
    "QPushButton", "QListWidget", "QLabel", "QLineEdit", "QFileDialog",
    "QComboBox", "QAbstractItemView", "QMessageBox", "QMenuBar", "QMenu",
    "QFrame", "QScrollArea", "QDialog", "QFormLayout", "QDialogButtonBox",
    "QGroupBox", "QCheckBox", "QColorDialog", "QTableWidget",
    "QTableWidgetItem", "QInputDialog",
    "Qt", "QDir", "QThread", "Signal", "QUrl",
    "QAction", "QDoubleValidator", "QDesktopServices",
]
