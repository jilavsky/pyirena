
import os
import sys
from typing import List, Optional, Dict
from pathlib import Path
import numpy as np

try:
    from PySide6.QtWidgets import (
        QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
        QLabel, QLineEdit, QCheckBox, QSpinBox, QTabWidget, QGroupBox,
        QGridLayout, QMessageBox, QSplitter, QFileDialog
    )
    from PySide6.QtCore import Qt, Signal
    from PySide6.QtGui import QFont, QDoubleValidator
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
            QLabel, QLineEdit, QCheckBox, QSpinBox, QTabWidget, QGroupBox,
            QGridLayout, QMessageBox, QSplitter, QFileDialog
        )
        from PyQt6.QtCore import Qt, pyqtSignal as Signal
        from PyQt6.QtGui import QFont, QDoubleValidator
    except ImportError:
        try:
            from PyQt5.QtWidgets import (
                QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
                QLabel, QLineEdit, QCheckBox, QSpinBox, QTabWidget, QGroupBox,
                QGridLayout, QMessageBox, QSplitter, QFileDialog
            )
            from PyQt5.QtCore import Qt, pyqtSignal as Signal
            from PyQt5.QtGui import QFont, QDoubleValidator
        except ImportError:
            raise ImportError(
                "No Qt found. Install with: pip install PySide6 or PyQt6 or PyQt5"
            )

import pyqtgraph as pg

from pyirena.core.unified import UnifiedFitModel, UnifiedLevel
from pyirena.state import StateManager


class DraggableCursor(pg.GraphicsObject):
    """
    Igor-style draggable cursor with vertical line and marker symbol.

    Combines a vertical dashed line with a draggable marker (circle/square/cross)
    at the top, similar to Igor Pro cursors.
    """

    sigPositionChanged = Signal(object)  # Signal emitted when cursor moves

    def __init__(self, pos, symbol='o', color=(255, 0, 0), label=''):
        super().__init__()
        self.pos_x = pos
        self.symbol = symbol
        self.color = color
        self.label_text = label
        self.plot_item = None
        self.marker_size = 10

        # Enable mouse interaction
        self.setAcceptHoverEvents(True)
        self.setFlag(pg.GraphicsObject.ItemIsMovable, False)  # Custom movement
        self.setZValue(1000)  # Draw on top

        # CRITICAL: Tell ViewBox to ignore this item when autoscaling
        # This prevents the cursor from affecting the data range
        self._shape = None
        self._boundingRect = None

    def set_plot_item(self, plot_item):
        """Set the parent plot item for coordinate conversion."""
        self.plot_item = plot_item

    def boundingRect(self):
        """
        Return bounding rectangle in item coordinates.

        NOTE: Returns a large static rectangle to ensure rendering works.
        dataBounds() returning None ensures we don't affect autoscale.
        """
        try:
            from PySide6.QtCore import QRectF
        except ImportError:
            try:
                from PyQt6.QtCore import QRectF
            except ImportError:
                from PyQt5.QtCore import QRectF

        # Return a large rectangle that covers any reasonable plot range
        # This is safe because dataBounds() returns None, so this won't
        # affect autoscaling
        return QRectF(-1000, -1e10, 2000, 2e10)

    def paint(self, painter, option, widget):
        """Paint the cursor."""
        print(f"DEBUG: paint() called for cursor at x={self.pos_x}")  # DEBUG
        if self.plot_item is None:
            print("DEBUG: plot_item is None, skipping paint")  # DEBUG
            return

        try:
            from PySide6.QtCore import QPointF, QRectF
            from PySide6.QtGui import QPen, QBrush, QColor
            from PySide6.QtCore import Qt as QtCore
        except ImportError:
            try:
                from PyQt6.QtCore import QPointF, QRectF
                from PyQt6.QtGui import QPen, QBrush, QColor
                from PyQt6.QtCore import Qt as QtCore
            except ImportError:
                from PyQt5.QtCore import QPointF, QRectF, Qt as QtCore
                from PyQt5.QtGui import QPen, QBrush, QColor

        # Get plot Y range in data coordinates
        view_range = self.plot_item.viewRange()
        y_range = view_range[1]  # (y_min, y_max) in log space if log scale
        x_range = view_range[0]  # (x_min, x_max)

        # For log scale, y_range is already in log space (e.g., 2 to 10 for 10^2 to 10^10)
        # Convert to linear for drawing
        if self.plot_item.ctrl.logYCheck.isChecked():
            y_min = max(10**y_range[0], 1e-300)  # Clamp to avoid underflow
            y_max = min(10**y_range[1], 1e300)   # Clamp to avoid overflow
        else:
            y_min, y_max = y_range

        # Safety checks
        if not (abs(y_min) < 1e100 and abs(y_max) < 1e100):
            return  # Skip drawing if values are too extreme
        if y_max <= y_min:
            return  # Invalid range

        # Draw vertical dashed line spanning the full plot height
        pen = QPen(QColor(*self.color))
        pen.setStyle(QtCore.PenStyle.DashLine)
        pen.setWidth(2)
        painter.setPen(pen)

        # Line from bottom to near top (in data coordinates)
        painter.drawLine(QPointF(0, y_min), QPointF(0, y_max * 0.92))

        # Draw marker at top
        marker_y = y_max * 0.96
        pen.setStyle(QtCore.PenStyle.SolidLine)
        pen.setWidth(2)
        painter.setPen(pen)

        # Calculate marker size - use log-space for better scaling
        x_span = x_range[1] - x_range[0]
        if x_span <= 0:
            return

        # Marker size as fraction of view range
        if self.plot_item.ctrl.logYCheck.isChecked():
            # In log scale, work in log space for Y
            log_y_span = y_range[1] - y_range[0]
            # Convert to data space carefully
            marker_size_y_fraction = 0.03  # 3% of log range
            marker_size_y = 10**(y_range[1] - marker_size_y_fraction * log_y_span / 2)
            marker_size_y = min(marker_size_y, y_max * 0.05)  # Clamp to 5% of max
        else:
            marker_size_y = (y_max - y_min) * 0.02

        # X marker size
        marker_size_x = x_span * 0.015  # 1.5% of X range

        # Additional safety: clamp marker sizes
        marker_size_x = max(min(marker_size_x, x_span * 0.1), x_span * 0.005)
        if self.plot_item.ctrl.logYCheck.isChecked():
            marker_size_y = max(min(marker_size_y, y_max * 0.1), y_max * 0.001)
        else:
            y_span = y_max - y_min
            marker_size_y = max(min(marker_size_y, y_span * 0.1), y_span * 0.001)

        # Draw symbols
        if self.symbol == 'o':  # Circle
            painter.setBrush(QBrush(QtCore.BrushStyle.NoBrush))
            # Use QRectF for ellipse to ensure float precision
            try:
                painter.drawEllipse(QPointF(0, marker_y), float(marker_size_x), float(marker_size_y))
            except (OverflowError, ValueError):
                pass  # Skip if still overflows

        elif self.symbol == 's':  # Square
            painter.setBrush(QBrush(QtCore.BrushStyle.NoBrush))
            try:
                rect = QRectF(
                    float(-marker_size_x),
                    float(marker_y - marker_size_y),
                    float(marker_size_x * 2),
                    float(marker_size_y * 2)
                )
                painter.drawRect(rect)
            except (OverflowError, ValueError):
                pass  # Skip if still overflows

        elif self.symbol == '+':  # Cross
            try:
                painter.drawLine(
                    QPointF(float(-marker_size_x), float(marker_y)),
                    QPointF(float(marker_size_x), float(marker_y))
                )
                painter.drawLine(
                    QPointF(0, float(marker_y - marker_size_y)),
                    QPointF(0, float(marker_y + marker_size_y))
                )
            except (OverflowError, ValueError):
                pass  # Skip if still overflows

    def setPos(self, x, y=0):
        """Override setPos to only allow X movement."""
        self.pos_x = x
        super().setPos(x, 0)  # Y is always 0 in item coordinates

    def value(self):
        """Get current X position."""
        return self.pos_x

    def setValue(self, x):
        """Set X position."""
        self.setPos(x)

    def dataBounds(self, ax, frac=1.0, orthoRange=None):
        """
        Override dataBounds to return None, preventing this item from
        affecting the ViewBox's autoscale range.
        """
        return None

    def mouseDragEvent(self, ev):
        """Handle mouse drag - only horizontal movement."""
        if ev.button() == Qt.MouseButton.LeftButton:
            if ev.isStart():
                self.moving = True
                self.cursorOffset = self.pos() - self.mapToParent(ev.buttonDownPos())
                ev.accept()
            elif ev.isFinish():
                self.moving = False
                ev.accept()
            else:
                if self.moving:
                    # Only allow horizontal movement
                    new_pos = self.mapToParent(ev.pos()) + self.cursorOffset
                    self.setPos(new_pos.x(), 0)  # Lock Y to 0
                    self.pos_x = new_pos.x()
                    self.sigPositionChanged.emit(self)
                    ev.accept()


class UnifiedFitGraphWindow(QWidget):
    """
    Graph window for displaying data and unified fit results using pyqtgraph.

    Much faster than matplotlib version with built-in cursor support.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena - Unified Fit")
        self.setGeometry(100, 100, 900, 700)

        # Create pyqtgraph layout widget
        self.graphics_layout = pg.GraphicsLayoutWidget()

        # Plot widget references (renamed to avoid conflicts with methods)
        self.main_plot = None
        self.residual_plot = None

        # Cursor attributes
        self.cursor_left = None
        self.cursor_right = None
        self.cursor_left_line = None  # pg.InfiniteLine object
        self.cursor_right_line = None  # pg.InfiniteLine object
        self.q_data = None

        # Layout
        layout = QVBoxLayout()
        layout.addWidget(self.graphics_layout)
        self.setLayout(layout)

        # Initialize plots
        self.init_plots()

    def init_plots(self):
        """Initialize the plot widgets."""
        self.graphics_layout.clear()

        # Main plot (data + fit)
        self.main_plot = self.graphics_layout.addPlot(row=0, col=0)
        self.main_plot.setLabel('bottom', 'Q (Å⁻¹)')
        self.main_plot.setLabel('left', 'Intensity (cm⁻¹)')
        self.main_plot.setLogMode(x=True, y=True)
        self.main_plot.showGrid(x=True, y=True, alpha=0.3)
        self.main_plot.setTitle('Unified Fit Model', size='12pt')
        self.main_plot.addLegend()

        # Enable auto-range
        self.main_plot.enableAutoRange()

        # Residuals plot
        self.residual_plot = self.graphics_layout.addPlot(row=1, col=0)
        self.residual_plot.setLabel('bottom', 'Q (Å⁻¹)')
        self.residual_plot.setLabel('left', 'Residuals')
        self.residual_plot.setLogMode(x=True, y=False)
        self.residual_plot.showGrid(x=True, y=True, alpha=0.3)

        # Enable auto-range
        self.residual_plot.enableAutoRange()

        # Add zero line for residuals
        self.residual_plot.addLine(y=0, pen=pg.mkPen('k', style=Qt.PenStyle.DashLine))

    def plot_data(self, q, intensity, error=None, label='Data'):
        """Plot experimental data."""
        self.q_data = q

        # Plot data points
        self.main_plot.plot(
            q, intensity,
            pen=None,
            symbol='o',
            symbolSize=4,
            symbolBrush=(100, 100, 255, 150),
            name=label
        )

        # Add error bars if available
        if error is not None:
            error_bars = pg.ErrorBarItem(
                x=q, y=intensity,
                height=error,
                pen=pg.mkPen((100, 100, 255, 100), width=1)
            )
            # IMPORTANT: ignoreBounds=True prevents error bars from affecting autoscale
            self.main_plot.addItem(error_bars, ignoreBounds=True)

        # ALWAYS reset cursor positions when new data is loaded
        # This ensures cursors are within the data range
        if len(q) > 0:
            q_min, q_max = q.min(), q.max()
            self.cursor_left = q_min * (q_max / q_min) ** 0.2
            self.cursor_right = q_min * (q_max / q_min) ** 0.8

        # Always add cursors if positions exist
        if self.cursor_left is not None and self.cursor_right is not None:
            self.add_cursors()

        # Force autoscale to data range only
        self.main_plot.enableAutoRange()
        self.main_plot.setMouseEnabled(x=True, y=True)  # Ensure both axes are zoomable

    def plot_fit(self, q, fit, label='Unified Fit'):
        """Plot fit curve."""
        self.main_plot.plot(
            q, fit,
            pen=pg.mkPen('r', width=2),
            name=label
        )

    def plot_residuals(self, q, residuals):
        """Plot fit residuals."""
        self.residual_plot.clear()
        self.residual_plot.addLine(y=0, pen=pg.mkPen('k', style=Qt.PenStyle.DashLine))
        self.residual_plot.plot(
            q, residuals,
            pen=None,
            symbol='o',
            symbolSize=3,
            symbolBrush=(100, 100, 255, 150)
        )

    def add_cursors(self):
        """
        Add draggable cursor lines.

        Using pyqtgraph's InfiniteLine for reliability.
        """
        # Remove old cursors if they exist
        if self.cursor_left_line is not None:
            try:
                self.main_plot.removeItem(self.cursor_left_line)
            except:
                pass
        if self.cursor_right_line is not None:
            try:
                self.main_plot.removeItem(self.cursor_right_line)
            except:
                pass

        # IMPORTANT: Since X-axis is in log mode, InfiniteLine expects positions in log space!
        # Convert linear positions to log10
        import numpy as np
        cursor_left_log = np.log10(self.cursor_left) if self.cursor_left > 0 else -10
        cursor_right_log = np.log10(self.cursor_right) if self.cursor_right > 0 else -10

        # Left cursor: RED, dashed
        self.cursor_left_line = pg.InfiniteLine(
            pos=cursor_left_log,  # Position in log space
            angle=90,  # Vertical
            movable=True,
            pen=pg.mkPen(color=(255, 0, 0), width=2, style=Qt.PenStyle.DashLine),
            label='A',
            labelOpts={'position': 0.95, 'color': (200, 0, 0), 'fill': (200, 0, 0, 50)}
        )
        self.cursor_left_line.setZValue(100)  # Draw on top

        # Right cursor: BLUE, dashed
        self.cursor_right_line = pg.InfiniteLine(
            pos=cursor_right_log,  # Position in log space
            angle=90,  # Vertical
            movable=True,
            pen=pg.mkPen(color=(0, 0, 255), width=2, style=Qt.PenStyle.DashLine),
            label='B',
            labelOpts={'position': 0.95, 'color': (0, 0, 200), 'fill': (0, 0, 200, 50)}
        )
        self.cursor_right_line.setZValue(100)  # Draw on top

        # Connect signals for cursor movement
        self.cursor_left_line.sigPositionChanged.connect(self.on_left_cursor_moved)
        self.cursor_right_line.sigPositionChanged.connect(self.on_right_cursor_moved)

        # Add to plot with ignoreBounds to prevent affecting autoscale
        self.main_plot.addItem(self.cursor_left_line, ignoreBounds=True)
        self.main_plot.addItem(self.cursor_right_line, ignoreBounds=True)

    def on_left_cursor_moved(self, line):
        """Handle left cursor movement."""
        import numpy as np
        # InfiniteLine returns position in log space (since plot is in log mode)
        # Convert back to linear for storage
        new_pos_log = line.value()
        new_pos_linear = 10**new_pos_log

        # Ensure left cursor doesn't cross right cursor (compare in linear space)
        if new_pos_linear < self.cursor_right:
            self.cursor_left = new_pos_linear
        else:
            # Reset to valid position (convert back to log)
            line.setValue(np.log10(self.cursor_left))

    def on_right_cursor_moved(self, line):
        """Handle right cursor movement."""
        import numpy as np
        # InfiniteLine returns position in log space (since plot is in log mode)
        # Convert back to linear for storage
        new_pos_log = line.value()
        new_pos_linear = 10**new_pos_log

        # Ensure right cursor doesn't cross left cursor (compare in linear space)
        if new_pos_linear > self.cursor_left:
            self.cursor_right = new_pos_linear
        else:
            # Reset to valid position (convert back to log)
            line.setValue(np.log10(self.cursor_right))

    def get_cursor_range(self):
        """Get the current cursor Q range."""
        if self.cursor_left is not None and self.cursor_right is not None:
            return (self.cursor_left, self.cursor_right)
        return None


class LevelParametersWidget(QWidget):
    """
    Widget for a single level's parameters in the Unified Fit model.
    """

    parameter_changed = Signal()  # Signal emitted when any parameter changes

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
        self.g_value.editingFinished.connect(self.parameter_changed.emit)
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
        self.rg_value.editingFinished.connect(self.parameter_changed.emit)
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
        self.b_value.editingFinished.connect(self.parameter_changed.emit)
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
        self.p_value.editingFinished.connect(self.parameter_changed.emit)
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


"""
UnifiedFitPanel class - main control panel for Unified Fit GUI.
This will be merged into unified_fit_pyqtgraph.py
"""


class UnifiedFitPanel(QWidget):
    """
    Main Unified Fit panel for pyIrena with state management.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.graph_window = None
        self.data = None
        self.model = UnifiedFitModel()
        self.fit_result = None

        # State management
        self.state_manager = StateManager()

        self.init_ui()
        self.load_state()

    def init_ui(self):
        """Initialize the user interface."""
        self.setWindowTitle("pyIrena - Unified Fit Model")

        # Main horizontal splitter
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
            # Connect parameter changed signal
            level_widget.parameter_changed.connect(self.on_parameter_changed)
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

        # Additional buttons row with Reset
        additional_buttons = QHBoxLayout()

        self.reset_button = QPushButton("Reset to Defaults")
        self.reset_button.setMinimumHeight(30)
        self.reset_button.setStyleSheet("""
            QPushButton {
                background-color: #e67e22;
                color: white;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #d35400;
            }
        """)
        self.reset_button.clicked.connect(self.reset_to_defaults)
        additional_buttons.addWidget(self.reset_button)

        self.reset_unif_button = QPushButton("reset unif?")
        self.reset_unif_button.setMinimumHeight(30)
        additional_buttons.addWidget(self.reset_unif_button)

        self.fix_limits_button = QPushButton("Fix limits?")
        self.fix_limits_button.setMinimumHeight(30)
        additional_buttons.addWidget(self.fix_limits_button)

        layout.addLayout(additional_buttons)

        # Store local checkbox
        self.store_local_check = QCheckBox("Store local (Porod & Guinier) fits?")
        layout.addWidget(self.store_local_check)

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

        self.save_state_button = QPushButton("Save State")
        self.save_state_button.setMinimumHeight(30)
        self.save_state_button.setStyleSheet("""
            QPushButton {
                background-color: #3498db;
                color: white;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #2980b9;
            }
        """)
        self.save_state_button.clicked.connect(self.save_state)
        results_buttons1.addWidget(self.save_state_button)

        self.store_data_button = QPushButton("Store in Data Folder")
        self.store_data_button.setMinimumHeight(30)
        results_buttons1.addWidget(self.store_data_button)

        self.export_ascii_button = QPushButton("Export ASCII")
        self.export_ascii_button.setMinimumHeight(30)
        results_buttons1.addWidget(self.export_ascii_button)

        layout.addLayout(results_buttons1)

        # Results buttons row 2
        results_buttons2 = QHBoxLayout()

        self.export_params_button = QPushButton("Export Parameters")
        self.export_params_button.setMinimumHeight(30)
        self.export_params_button.clicked.connect(self.export_parameters)
        results_buttons2.addWidget(self.export_params_button)

        self.import_params_button = QPushButton("Import Parameters")
        self.import_params_button.setMinimumHeight(30)
        self.import_params_button.clicked.connect(self.import_parameters)
        results_buttons2.addWidget(self.import_params_button)

        self.results_graphs_button = QPushButton("Results to graphs")
        self.results_graphs_button.setMinimumHeight(30)
        results_buttons2.addWidget(self.results_graphs_button)

        layout.addLayout(results_buttons2)

        # Results buttons row 3
        results_buttons3 = QHBoxLayout()
        self.analyze_results_button = QPushButton("Analyze Results")
        self.analyze_results_button.setMinimumHeight(30)
        results_buttons3.addWidget(self.analyze_results_button)

        self.anal_uncertainty_button = QPushButton("Anal. Uncertainty")
        self.anal_uncertainty_button.setMinimumHeight(30)
        results_buttons3.addWidget(self.anal_uncertainty_button)

        self.ext_warnings_check = QCheckBox("Ext. warnings?")
        results_buttons3.addWidget(self.ext_warnings_check)

        layout.addLayout(results_buttons3)

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
        """Graph the unified fit with current parameters."""
        if self.data is None:
            QMessageBox.warning(self, "No Data", "Please load data first.")
            return

        try:
            # Get parameters
            num_levels = self.num_levels_spin.value()
            levels = []

            for i in range(num_levels):
                params = self.level_widgets[i].get_parameters()
                level = UnifiedLevel(
                    Rg=params['Rg'],
                    G=params['G'],
                    P=params['P'],
                    B=params['B'],
                    RgCO=params['RgCutoff'],
                    correlations=params['correlated']
                )
                levels.append(level)

            background = float(self.background_value.text() or 0)

            # Update model
            self.model.num_levels = num_levels
            self.model.levels = levels
            self.model.background = background

            # Calculate
            intensity_calc = self.model.calculate_intensity(self.data['Q'])

            # Plot
            self.graph_window.init_plots()
            self.graph_window.plot_data(
                self.data['Q'],
                self.data['Intensity'],
                self.data.get('Error'),
                self.data['label']
            )
            self.graph_window.plot_fit(self.data['Q'], intensity_calc, 'Unified Fit')

            # Residuals
            if self.data.get('Error') is not None:
                residuals = (self.data['Intensity'] - intensity_calc) / self.data['Error']
            else:
                residuals = (self.data['Intensity'] - intensity_calc) / self.data['Intensity']

            self.graph_window.plot_residuals(self.data['Q'], residuals)
            self.status_label.setText(f"Calculated unified fit with {num_levels} level(s)")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error calculating unified fit:\n{str(e)}")
            import traceback
            traceback.print_exc()

    def run_fit(self):
        """Run the unified fit."""
        if self.data is None:
            QMessageBox.warning(self, "No Data", "Please load data first.")
            return

        try:
            # Get parameters
            num_levels = self.num_levels_spin.value()
            levels = []

            for i in range(num_levels):
                params = self.level_widgets[i].get_parameters()
                level = UnifiedLevel(
                    Rg=params['Rg'],
                    G=params['G'],
                    P=params['P'],
                    B=params['B'],
                    RgCO=params['RgCutoff'],
                    correlations=params['correlated'],
                    fit_Rg=params['fit_Rg'],
                    fit_G=params['fit_G'],
                    fit_P=params['fit_P'],
                    fit_B=params['fit_B']
                )
                levels.append(level)

            background = float(self.background_value.text() or 0)
            fit_background = self.fit_background_check.isChecked()

            # Update model
            self.model.num_levels = num_levels
            self.model.levels = levels
            self.model.background = background
            self.model.fit_background = fit_background

            # Get cursor range and filter data
            cursor_range = self.graph_window.get_cursor_range()
            if cursor_range is not None:
                q_min, q_max = cursor_range
                mask = (self.data['Q'] >= q_min) & (self.data['Q'] <= q_max)
                q_fit = self.data['Q'][mask]
                intensity_fit = self.data['Intensity'][mask]
                error_fit = self.data.get('Error')[mask] if self.data.get('Error') is not None else None
                num_points = len(q_fit)
            else:
                q_fit = self.data['Q']
                intensity_fit = self.data['Intensity']
                error_fit = self.data.get('Error')
                num_points = len(q_fit)

            # Run fit
            self.status_label.setText(f"Fitting {num_levels} level(s) with {num_points} points...")

            result = self.model.fit(q_fit, intensity_fit, error_fit)
            self.fit_result = result

            # Update GUI
            for i in range(num_levels):
                fitted_level = result['levels'][i]
                self.level_widgets[i].set_parameters({
                    'G': fitted_level.G,
                    'Rg': fitted_level.Rg,
                    'B': fitted_level.B,
                    'P': fitted_level.P
                })

            self.background_value.setText(f"{result['background']:.6e}")

            # Calculate and plot
            intensity_calc = self.model.calculate_intensity(self.data['Q'])

            self.graph_window.init_plots()
            self.graph_window.plot_data(
                self.data['Q'],
                self.data['Intensity'],
                self.data.get('Error'),
                self.data['label']
            )
            self.graph_window.plot_fit(self.data['Q'], intensity_calc, 'Fitted Model')

            # Residuals
            if self.data.get('Error') is not None:
                residuals = (self.data['Intensity'] - intensity_calc) / self.data['Error']
            else:
                residuals = (self.data['Intensity'] - intensity_calc) / self.data['Intensity']

            self.graph_window.plot_residuals(self.data['Q'], residuals)

            # Show statistics
            chi2 = result.get('chi_squared', 0.0)
            reduced_chi2 = result.get('reduced_chi_squared', 0.0)

            if cursor_range is not None:
                self.status_label.setText(
                    f"Fit complete: {num_levels} level(s), {num_points} pts (Q: {q_min:.3e} - {q_max:.3e}), χ² = {chi2:.4f}"
                )
            else:
                self.status_label.setText(f"Fit complete: {num_levels} level(s), χ² = {chi2:.4f}")

            cursor_info = f"\nFit range: Q = {q_min:.3e} to {q_max:.3e}\nPoints used: {num_points}" if cursor_range else ""

            QMessageBox.information(
                self, "Fit Complete",
                f"Fit completed successfully!\n\n"
                f"Number of levels: {num_levels}\n"
                f"Chi-squared: {chi2:.4f}\n"
                f"Reduced χ²: {reduced_chi2:.4f}\n"
                f"Parameters updated in GUI.{cursor_info}"
            )

        except Exception as e:
            QMessageBox.critical(self, "Fit Error", f"Error during fitting:\n{str(e)}")
            self.status_label.setText("Fit failed")
            import traceback
            traceback.print_exc()

    # STATE MANAGEMENT METHODS

    def get_current_state(self) -> Dict:
        """Get current GUI state for saving."""
        state = {
            'num_levels': self.num_levels_spin.value(),
            'levels': [],
            'background': {
                'value': float(self.background_value.text() or 0),
                'fit': self.fit_background_check.isChecked()
            },
            'cursor_left': self.graph_window.cursor_left,
            'cursor_right': self.graph_window.cursor_right,
            'update_auto': self.update_auto_check.isChecked(),
            'display_local': self.display_local_check.isChecked(),
            'no_limits': self.no_limits_check.isChecked(),
            'skip_fit_check': self.skip_fit_check.isChecked(),
            'store_local': self.store_local_check.isChecked()
        }

        # Get all level parameters
        for i in range(5):
            params = self.level_widgets[i].get_parameters()
            level_state = {
                'level': i + 1,
                'G': {
                    'value': params['G'],
                    'fit': params['fit_G'],
                    'low_limit': None,  # TODO: Add limit fields to GUI
                    'high_limit': None
                },
                'Rg': {
                    'value': params['Rg'],
                    'fit': params['fit_Rg'],
                    'low_limit': None,
                    'high_limit': None
                },
                'B': {
                    'value': params['B'],
                    'fit': params['fit_B'],
                    'low_limit': None,
                    'high_limit': None
                },
                'P': {
                    'value': params['P'],
                    'fit': params['fit_P'],
                    'low_limit': None,
                    'high_limit': None
                },
                'RgCutoff': params['RgCutoff'],
                'correlated': params['correlated'],
                'estimate_B': params['estimate_B']
            }
            state['levels'].append(level_state)

        return state

    def apply_state(self, state: Dict):
        """Apply saved state to GUI."""
        # Set number of levels
        self.num_levels_spin.setValue(state.get('num_levels', 1))

        # Set background
        bg = state.get('background', {})
        self.background_value.setText(str(bg.get('value', 1e-6)))
        self.fit_background_check.setChecked(bg.get('fit', False))

        # Set checkboxes
        self.update_auto_check.setChecked(state.get('update_auto', False))
        self.display_local_check.setChecked(state.get('display_local', False))
        self.no_limits_check.setChecked(state.get('no_limits', False))
        self.skip_fit_check.setChecked(state.get('skip_fit_check', False))
        self.store_local_check.setChecked(state.get('store_local', False))

        # Set cursor positions
        if state.get('cursor_left') is not None:
            self.graph_window.cursor_left = state['cursor_left']
        if state.get('cursor_right') is not None:
            self.graph_window.cursor_right = state['cursor_right']

        # Set level parameters
        levels = state.get('levels', [])
        for i, level_state in enumerate(levels):
            if i < 5:
                level_widget = self.level_widgets[i]

                # Set parameter values
                level_widget.g_value.setText(str(level_state['G']['value']))
                level_widget.g_fit.setChecked(level_state['G']['fit'])

                level_widget.rg_value.setText(str(level_state['Rg']['value']))
                level_widget.rg_fit.setChecked(level_state['Rg']['fit'])

                level_widget.b_value.setText(str(level_state['B']['value']))
                level_widget.b_fit.setChecked(level_state['B']['fit'])

                level_widget.p_value.setText(str(level_state['P']['value']))
                level_widget.p_fit.setChecked(level_state['P']['fit'])

                level_widget.rg_cutoff.setText(str(level_state.get('RgCutoff', 0)))
                level_widget.correlated_check.setChecked(level_state.get('correlated', False))
                level_widget.estimate_b_check.setChecked(level_state.get('estimate_B', False))

    def load_state(self):
        """Load state from state manager."""
        state = self.state_manager.get('unified_fit')
        if state:
            self.apply_state(state)
            print("Loaded saved state")

    def save_state(self):
        """Save current state."""
        state = self.get_current_state()
        self.state_manager.update('unified_fit', state)
        if self.state_manager.save():
            QMessageBox.information(self, "State Saved", "Current state has been saved successfully!")
            self.status_label.setText("State saved")
        else:
            QMessageBox.warning(self, "Save Failed", "Failed to save state")

    def reset_to_defaults(self):
        """Reset all parameters to defaults."""
        reply = QMessageBox.question(
            self,
            "Reset to Defaults",
            "Are you sure you want to reset all parameters to their default values?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
        )

        if reply == QMessageBox.StandardButton.Yes:
            self.state_manager.reset('unified_fit')
            self.load_state()
            QMessageBox.information(self, "Reset Complete", "All parameters reset to defaults")
            self.status_label.setText("Reset to defaults")

    def export_parameters(self):
        """Export current parameters to JSON file."""
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Parameters",
            "",
            "JSON Files (*.json);;All Files (*)"
        )

        if file_path:
            state = self.get_current_state()
            if self.state_manager.export_tool_state('unified_fit', Path(file_path)):
                QMessageBox.information(self, "Export Complete", f"Parameters exported to:\n{file_path}")
            else:
                QMessageBox.warning(self, "Export Failed", "Failed to export parameters")

    def import_parameters(self):
        """Import parameters from JSON file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Import Parameters",
            "",
            "JSON Files (*.json);;All Files (*)"
        )

        if file_path:
            if self.state_manager.import_tool_state('unified_fit', Path(file_path)):
                self.load_state()
                QMessageBox.information(self, "Import Complete", f"Parameters imported from:\n{file_path}")
            else:
                QMessageBox.warning(self, "Import Failed", "Failed to import parameters")

    def closeEvent(self, event):
        """Handle window close - auto-save state."""
        state = self.get_current_state()
        self.state_manager.update('unified_fit', state)
        self.state_manager.save()
        event.accept()


def main():
    """Main entry point."""
    app = QApplication(sys.argv)
    app.setStyle('Fusion')

    window = UnifiedFitPanel()

    # Test data
    q = np.logspace(-3, 0, 100)
    intensity = 1e10 * np.exp(-(q * 100)**2 / 3) + 1e6 * q**-4
    error = 0.05 * intensity
    window.set_data(q, intensity, error, "Test Data")

    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
