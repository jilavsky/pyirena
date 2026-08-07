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

    QT_IMPORT_ERROR = None
except ImportError as _qt_error:
    # Degrade gracefully (importing pyirena must work without GUI extras) but
    # keep the diagnosed failure: callers that need Qt should re-raise this
    # rather than report a bare "QT_BACKEND is None".
    QT_BACKEND = None
    QtWidgets = None
    QtCore = None
    QtGui = None
    QT_IMPORT_ERROR = _qt_error

__all__ = ["QT_BACKEND", "QT_IMPORT_ERROR", "QtWidgets", "QtCore", "QtGui"]


def require_qt():
    """Raise the diagnosed Qt import failure if no binding is available.

    Use at the top of any entry point that cannot proceed without Qt, so the
    user sees the full diagnosis instead of an ``AttributeError`` on ``None``.

    Raises:
        ImportError: The original, diagnosed failure from ``pyirena.gui._qt``.
    """
    if QT_IMPORT_ERROR is not None:
        raise QT_IMPORT_ERROR
