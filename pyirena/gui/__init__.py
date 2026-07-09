"""
GUI module for pyIrena.

This module provides graphical user interfaces for data selection,
analysis, and visualization.

Requires:
    PySide6 or PyQt6 (install with: pip install pyirena[gui])
"""

try:
    # Single import point for the whole gui package (PySide6 → PyQt6).
    from pyirena.gui._qt import QtWidgets, QtCore, QtGui, QT_BINDING as QT_BACKEND
except ImportError:
    QT_BACKEND = None
    QtWidgets = None
    QtCore = None
    QtGui = None

__all__ = ["QT_BACKEND", "QtWidgets", "QtCore", "QtGui"]
