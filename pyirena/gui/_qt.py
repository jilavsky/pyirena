"""
pyirena.gui._qt — single PySide6/PyQt6 import point for the whole ``pyirena.gui``
package (PySide6 preferred, PyQt6 fallback).

Every GUI module imports the Qt names it needs from here instead of repeating a
``try: from PySide6 … except ImportError: from PyQt6 …`` block.  This removes
~15 duplicated import blocks and the associated ruff F401 noise (the "unused"
names in each fallback branch).

Two import styles are supported:

* flat names — ``from pyirena.gui._qt import QWidget, Qt, Signal``
* submodules — ``from pyirena.gui._qt import QtWidgets, QtCore, QtGui``

``Signal`` is normalised across bindings (``pyqtSignal`` under PyQt6).

When neither binding can be imported the failure is diagnosed by
``pyirena.diagnostics`` rather than reported as "not installed" — an
installed-but-unloadable Qt raises ``ImportError`` too, and telling the user
to reinstall it sends them in a circle.  Run ``pyirena-doctor`` for the full
report.
"""

try:
    from PySide6 import QtCore, QtGui, QtWidgets
    from PySide6.QtCore import Signal
    QT_BINDING = "PySide6"
except ImportError as _pyside_error:  # pragma: no cover - only without PySide6
    try:
        from PyQt6 import QtCore, QtGui, QtWidgets  # type: ignore[no-redef]
        from PyQt6.QtCore import pyqtSignal as Signal  # type: ignore[no-redef]
        QT_BINDING = "PyQt6"
    except ImportError as _pyqt_error:  # pragma: no cover
        # Do NOT report this as "not installed": ImportError also covers an
        # installed-but-unloadable binding (architecture mismatch, missing
        # system libraries, shiboken6 version skew).  pyirena.diagnostics
        # tells the two apart and names the interpreter that failed.
        from pyirena.diagnostics import format_qt_import_failure

        del _pyqt_error  # diagnosed below; PySide6 is the binding we report on
        raise ImportError(
            "pyIrena could not load a Qt binding.\n\n"
            + format_qt_import_failure()
        ) from _pyside_error

# --- QtWidgets ---------------------------------------------------------------
QAbstractItemView = QtWidgets.QAbstractItemView
QApplication = QtWidgets.QApplication
QButtonGroup = QtWidgets.QButtonGroup
QCheckBox = QtWidgets.QCheckBox
QColorDialog = QtWidgets.QColorDialog
QComboBox = QtWidgets.QComboBox
QDialog = QtWidgets.QDialog
QDialogButtonBox = QtWidgets.QDialogButtonBox
QDoubleSpinBox = QtWidgets.QDoubleSpinBox
QFileDialog = QtWidgets.QFileDialog
QFormLayout = QtWidgets.QFormLayout
QFrame = QtWidgets.QFrame
QGridLayout = QtWidgets.QGridLayout
QGroupBox = QtWidgets.QGroupBox
QHBoxLayout = QtWidgets.QHBoxLayout
QHeaderView = QtWidgets.QHeaderView
QInputDialog = QtWidgets.QInputDialog
QLabel = QtWidgets.QLabel
QLineEdit = QtWidgets.QLineEdit
QListWidget = QtWidgets.QListWidget
QListWidgetItem = QtWidgets.QListWidgetItem
QMainWindow = QtWidgets.QMainWindow
QMenu = QtWidgets.QMenu
QMenuBar = QtWidgets.QMenuBar
QMessageBox = QtWidgets.QMessageBox
QPlainTextEdit = QtWidgets.QPlainTextEdit
QProgressBar = QtWidgets.QProgressBar
QPushButton = QtWidgets.QPushButton
QRadioButton = QtWidgets.QRadioButton
QScrollArea = QtWidgets.QScrollArea
QSizePolicy = QtWidgets.QSizePolicy
QSlider = QtWidgets.QSlider
QSpacerItem = QtWidgets.QSpacerItem
QSpinBox = QtWidgets.QSpinBox
QSplitter = QtWidgets.QSplitter
QStatusBar = QtWidgets.QStatusBar
QTabWidget = QtWidgets.QTabWidget
QTableWidget = QtWidgets.QTableWidget
QTableWidgetItem = QtWidgets.QTableWidgetItem
QTextBrowser = QtWidgets.QTextBrowser
QTextEdit = QtWidgets.QTextEdit
QToolBar = QtWidgets.QToolBar
QTreeWidget = QtWidgets.QTreeWidget
QTreeWidgetItem = QtWidgets.QTreeWidgetItem
QVBoxLayout = QtWidgets.QVBoxLayout
QWidget = QtWidgets.QWidget

# --- QtCore ------------------------------------------------------------------
Qt = QtCore.Qt
QBuffer = QtCore.QBuffer
QByteArray = QtCore.QByteArray
QDir = QtCore.QDir
QEvent = QtCore.QEvent
QIODevice = QtCore.QIODevice
QMimeData = QtCore.QMimeData
QPointF = QtCore.QPointF
QRectF = QtCore.QRectF
QThread = QtCore.QThread
QTimer = QtCore.QTimer
QUrl = QtCore.QUrl

# --- QtGui -------------------------------------------------------------------
QAction = QtGui.QAction
QBrush = QtGui.QBrush
QCloseEvent = QtGui.QCloseEvent
QColor = QtGui.QColor
QDesktopServices = QtGui.QDesktopServices
QDoubleValidator = QtGui.QDoubleValidator
QFont = QtGui.QFont
QIcon = QtGui.QIcon
QIntValidator = QtGui.QIntValidator
QKeySequence = QtGui.QKeySequence
QPainterPath = QtGui.QPainterPath
QPainterPathStroker = QtGui.QPainterPathStroker
QPen = QtGui.QPen
QPixmap = QtGui.QPixmap
QShortcut = QtGui.QShortcut
QTransform = QtGui.QTransform

__all__ = [
    "QtWidgets", "QtCore", "QtGui", "Signal", "QT_BINDING",
    # QtWidgets
    "QAbstractItemView", "QApplication", "QButtonGroup", "QCheckBox",
    "QColorDialog", "QComboBox", "QDialog", "QDialogButtonBox", "QDoubleSpinBox",
    "QFileDialog", "QFormLayout", "QFrame", "QGridLayout", "QGroupBox",
    "QHBoxLayout", "QHeaderView", "QInputDialog", "QLabel", "QLineEdit",
    "QListWidget", "QListWidgetItem", "QMainWindow", "QMenu", "QMenuBar",
    "QMessageBox", "QPlainTextEdit", "QProgressBar", "QPushButton",
    "QRadioButton", "QScrollArea", "QSizePolicy", "QSlider", "QSpacerItem",
    "QSpinBox", "QSplitter", "QStatusBar", "QTabWidget", "QTableWidget",
    "QTableWidgetItem", "QTextBrowser", "QTextEdit", "QToolBar", "QTreeWidget",
    "QTreeWidgetItem", "QVBoxLayout", "QWidget",
    # QtCore
    "Qt", "QBuffer", "QByteArray", "QDir", "QEvent", "QIODevice", "QMimeData",
    "QPointF",
    "QRectF", "QThread", "QTimer", "QUrl",
    # QtGui
    "QAction", "QBrush", "QCloseEvent", "QColor", "QDesktopServices",
    "QDoubleValidator", "QFont", "QIcon", "QIntValidator", "QKeySequence",
    "QPainterPath", "QPainterPathStroker", "QPen", "QPixmap", "QShortcut",
    "QTransform",
]
