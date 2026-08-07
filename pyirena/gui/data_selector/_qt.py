"""
pyirena.gui.data_selector._qt — single PySide6/PyQt6 import point for the
data_selector package (PySide6 preferred, PyQt6 fallback).

Failure diagnosis is delegated to ``pyirena.diagnostics``; see
``pyirena/gui/_qt.py`` for the rationale.
"""

try:
    from PySide6.QtCore import QDir, Qt, QThread, QUrl, Signal
    from PySide6.QtGui import QAction, QDesktopServices, QDoubleValidator, QKeySequence, QShortcut
    from PySide6.QtWidgets import (
        QAbstractItemView,
        QApplication,
        QCheckBox,
        QColorDialog,
        QComboBox,
        QDialog,
        QDialogButtonBox,
        QFileDialog,
        QFormLayout,
        QFrame,
        QGridLayout,
        QGroupBox,
        QHBoxLayout,
        QInputDialog,
        QLabel,
        QLineEdit,
        QListWidget,
        QMenu,
        QMenuBar,
        QMessageBox,
        QPushButton,
        QScrollArea,
        QTableWidget,
        QTableWidgetItem,
        QVBoxLayout,
        QWidget,
    )
except ImportError as _pyside_error:
    try:
        from PyQt6.QtCore import QDir, Qt, QThread, QUrl  # type: ignore[no-redef]
        from PyQt6.QtCore import pyqtSignal as Signal
        from PyQt6.QtGui import (  # type: ignore[no-redef]
            QAction,
            QDesktopServices,
            QDoubleValidator,
            QKeySequence,
            QShortcut,
        )
        from PyQt6.QtWidgets import (  # type: ignore[no-redef]
            QAbstractItemView,
            QApplication,
            QCheckBox,
            QColorDialog,
            QComboBox,
            QDialog,
            QDialogButtonBox,
            QFileDialog,
            QFormLayout,
            QFrame,
            QGridLayout,
            QGroupBox,
            QHBoxLayout,
            QInputDialog,
            QLabel,
            QLineEdit,
            QListWidget,
            QMenu,
            QMenuBar,
            QMessageBox,
            QPushButton,
            QScrollArea,
            QTableWidget,
            QTableWidgetItem,
            QVBoxLayout,
            QWidget,
        )
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
