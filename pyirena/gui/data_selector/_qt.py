"""
pyirena.gui.data_selector._qt — single PySide6/PyQt6 import point for the
data_selector package (PySide6 preferred, PyQt6 fallback).

Failure diagnosis is delegated to ``pyirena.diagnostics``; see
``pyirena/gui/_qt.py`` for the rationale.
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
    from PySide6.QtGui import QAction, QDoubleValidator, QDesktopServices, QKeySequence, QShortcut
except ImportError as _pyside_error:
    try:
        from PyQt6.QtWidgets import (  # type: ignore[no-redef]
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QPushButton,
        QListWidget, QLabel, QLineEdit, QFileDialog, QComboBox,
        QAbstractItemView, QMessageBox, QMenuBar, QMenu, QFrame, QScrollArea,
        QDialog, QFormLayout, QDialogButtonBox, QGroupBox, QCheckBox, QColorDialog,
        QTableWidget, QTableWidgetItem, QInputDialog,
        )
        from PyQt6.QtCore import Qt, QDir, QThread, pyqtSignal as Signal, QUrl  # type: ignore[no-redef]
        from PyQt6.QtGui import QAction, QDoubleValidator, QDesktopServices, QKeySequence, QShortcut  # type: ignore[no-redef]
    except ImportError:
        from pyirena.diagnostics import format_qt_import_failure

        raise ImportError(
            "pyIrena could not load a Qt binding.\n\n"
            + format_qt_import_failure()
        ) from _pyside_error

__all__ = [
    "QApplication", "QWidget", "QVBoxLayout", "QHBoxLayout", "QGridLayout",
    "QPushButton", "QListWidget", "QLabel", "QLineEdit", "QFileDialog",
    "QComboBox", "QAbstractItemView", "QMessageBox", "QMenuBar", "QMenu",
    "QFrame", "QScrollArea", "QDialog", "QFormLayout", "QDialogButtonBox",
    "QGroupBox", "QCheckBox", "QColorDialog", "QTableWidget",
    "QTableWidgetItem", "QInputDialog",
    "Qt", "QDir", "QThread", "Signal", "QUrl",
    "QAction", "QDoubleValidator", "QDesktopServices", "QKeySequence", "QShortcut",
]
