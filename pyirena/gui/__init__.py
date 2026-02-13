"""
GUI module for pyIrena.

This module provides graphical user interfaces for data selection,
analysis, and visualization.

Requires:
    PySide6 or PyQt6 (install with: pip install pyirena[gui])
"""

try:
    # Try PySide6 first
    from PySide6 import QtWidgets, QtCore, QtGui
    QT_BACKEND = "PySide6"
except ImportError:
    try:
        # Fall back to PyQt6
        from PyQt6 import QtWidgets, QtCore, QtGui
        QT_BACKEND = "PyQt6"
    except ImportError:
        QT_BACKEND = None
        QtWidgets = None
        QtCore = None
        QtGui = None

__all__ = ["QT_BACKEND", "QtWidgets", "QtCore", "QtGui"]
