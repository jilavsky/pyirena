
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


class ScrubbableLineEdit(QLineEdit):
    """
    Enhanced QLineEdit that allows changing values with mouse wheel.

    Features:
    - Mouse wheel alone: Change by 1% of current value (or by 1 if small)
    - Shift + wheel: Change by 10% of current value (coarse)
    - Ctrl/Cmd + wheel: Change by 0.1% of current value (fine)

    For most parameters (B, Rg, G), step adapts to current value so you can explore large ranges.
    For P parameter, use fixed step based on typical values (3-4) for precision.
    """

    def __init__(self, *args, use_fixed_step=False, fixed_reference=4.0, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.use_fixed_step = use_fixed_step
        self.fixed_reference = fixed_reference

    def wheelEvent(self, event):
        """Handle mouse wheel events to change numerical values."""
        # Only process if widget has focus or mouse is over it
        if not self.hasFocus():
            # Give focus when wheel is used
            self.setFocus()

        # Try to get current value
        try:
            current_value = float(self.text())
        except (ValueError, AttributeError):
            # Not a valid number, don't process
            super().wheelEvent(event)
            return

        # Determine delta based on wheel movement
        # event.angleDelta().y() is typically +120 or -120 for one notch
        delta_notches = event.angleDelta().y() / 120.0

        # Determine step size based on modifier keys
        modifiers = event.modifiers()

        if modifiers & Qt.KeyboardModifier.ShiftModifier:
            # Shift = coarse adjustment (10% of value)
            step_factor = 0.1
        elif modifiers & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier):
            # Ctrl/Cmd = fine adjustment (0.1% of value)
            step_factor = 0.001
        else:
            # No modifier = medium adjustment (1% of value)
            step_factor = 0.01

        # Calculate step size
        # Use fixed reference for P parameter (stays in narrow range)
        # Use current value for B, Rg, G (can span many orders of magnitude)
        if self.use_fixed_step:
            # P parameter: step based on fixed reference (e.g., 4.0)
            reference_value = self.fixed_reference
        else:
            # B, Rg, G: step adapts to current value
            reference_value = current_value

        # Special case: if exactly at zero, use small absolute step to get started
        if abs(current_value) == 0:
            step = 0.01  # Small starting value when at zero
        else:
            # Pure percentage-based step (works from 1e-9 to 1e9)
            step = abs(reference_value) * step_factor

        # Apply delta
        new_value = current_value + (delta_notches * step)

        # Physical parameters (G, Rg, B, P) should never be negative
        # Enforce minimum of 0 (or validator minimum if set)
        minimum_value = 0.0

        # Check validator limits if present
        if self.validator():
            validator = self.validator()
            if hasattr(validator, 'bottom') and hasattr(validator, 'top'):
                # QDoubleValidator has bottom/top
                bottom = validator.bottom()
                top = validator.top()
                minimum_value = max(0.0, bottom)  # Use validator min but not less than 0
                new_value = max(minimum_value, min(top, new_value))
            else:
                # No explicit limits, but still enforce >= 0
                new_value = max(minimum_value, new_value)
        else:
            # No validator, enforce >= 0
            new_value = max(minimum_value, new_value)

        # Update the text field
        # Determine appropriate precision
        if abs(new_value) < 0.001:
            text = f"{new_value:.2e}"
        elif abs(new_value) < 1:
            text = f"{new_value:.6f}".rstrip('0').rstrip('.')
        elif abs(new_value) < 1000:
            text = f"{new_value:.4f}".rstrip('0').rstrip('.')
        else:
            text = f"{new_value:.2e}"

        self.setText(text)
        self.editingFinished.emit()

        # Accept the event so it doesn't propagate
        event.accept()


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

        # Status message area (2-3 lines under graphs)
        self.status_message = QLabel("")
        self.status_message.setWordWrap(True)
        self.status_message.setMaximumHeight(60)
        self.status_message.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.status_message.setStyleSheet("""
            QLabel {
                padding: 8px;
                border: 1px solid #ccc;
                border-radius: 4px;
                font-size: 11pt;
            }
        """)
        layout.addWidget(self.status_message)

        self.setLayout(layout)

        # Initialize plots
        self.init_plots()

    def show_success_message(self, message):
        """Show success message with green background."""
        self.status_message.setText(message)
        self.status_message.setStyleSheet("""
            QLabel {
                background-color: #d4edda;
                color: #155724;
                padding: 8px;
                border: 1px solid #c3e6cb;
                border-radius: 4px;
                font-size: 11pt;
            }
        """)

    def show_error_message(self, message):
        """Show error message with red background."""
        self.status_message.setText(message)
        self.status_message.setStyleSheet("""
            QLabel {
                background-color: #f8d7da;
                color: #721c24;
                padding: 8px;
                border: 1px solid #f5c6cb;
                border-radius: 4px;
                font-size: 11pt;
            }
        """)

    def clear_message(self):
        """Clear the status message."""
        self.status_message.setText("")
        self.status_message.setStyleSheet("""
            QLabel {
                padding: 8px;
                border: 1px solid #ccc;
                border-radius: 4px;
                font-size: 11pt;
            }
        """)

    def init_plots(self):
        """Initialize the plot widgets."""
        self.graphics_layout.clear()

        # Set white background for the graphics layout
        self.graphics_layout.setBackground('w')

        # Main plot (data + fit) - 80% of height
        self.main_plot = self.graphics_layout.addPlot(row=0, col=0)
        self.main_plot.setLabel('bottom', 'Q (Å⁻¹)', **{'color': 'k', 'font-size': '11pt'})
        self.main_plot.setLabel('left', 'Intensity (cm⁻¹)', **{'color': 'k', 'font-size': '11pt'})
        self.main_plot.setLogMode(x=True, y=True)
        self.main_plot.showGrid(x=True, y=True, alpha=0.3)
        self.main_plot.setTitle('Unified Fit Model', size='12pt', color='k')
        self.main_plot.addLegend()

        # Set axis colors to black for visibility on white background
        self.main_plot.getAxis('bottom').setPen('k')
        self.main_plot.getAxis('left').setPen('k')
        self.main_plot.getAxis('bottom').setTextPen('k')
        self.main_plot.getAxis('left').setTextPen('k')

        # Enable auto-range
        self.main_plot.enableAutoRange()

        # Residuals plot - 20% of height
        self.residual_plot = self.graphics_layout.addPlot(row=1, col=0)
        self.residual_plot.setLabel('bottom', 'Q (Å⁻¹)', **{'color': 'k', 'font-size': '11pt'})
        self.residual_plot.setLabel('left', 'Residuals', **{'color': 'k', 'font-size': '11pt'})
        self.residual_plot.setLogMode(x=True, y=False)
        self.residual_plot.showGrid(x=True, y=True, alpha=0.3)

        # Set axis colors to black for visibility on white background
        self.residual_plot.getAxis('bottom').setPen('k')
        self.residual_plot.getAxis('left').setPen('k')
        self.residual_plot.getAxis('bottom').setTextPen('k')
        self.residual_plot.getAxis('left').setTextPen('k')

        # Enable auto-range
        self.residual_plot.enableAutoRange()

        # Add zero line for residuals
        self.residual_plot.addLine(y=0, pen=pg.mkPen('k', style=Qt.PenStyle.DashLine))

        # Set height ratios: main plot gets 4 parts (80%), residuals get 1 part (20%)
        self.graphics_layout.ci.layout.setRowStretchFactor(0, 4)
        self.graphics_layout.ci.layout.setRowStretchFactor(1, 1)

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

        # Add error bars if available - batch all segments for fast rendering
        if error is not None and len(error) > 0:
            # Pre-allocate arrays for all line segments (vertical bars + caps)
            # Each error bar has 3 segments: vertical, top cap, bottom cap
            # Each segment needs 2 points + 1 NaN to disconnect
            n_points = len(q)

            # Create arrays with NaN separators for disconnected line segments
            x_lines = []
            y_lines = []

            cap_width_log = 0.05  # 5% in log space for cap width

            for i in range(n_points):
                # Calculate error bar limits, ensuring we stay positive for log scale
                y_top = intensity[i] + error[i]
                y_bottom = max(intensity[i] - error[i], intensity[i] * 0.001)

                # Vertical bar
                x_lines.extend([q[i], q[i], np.nan])
                y_lines.extend([y_bottom, y_top, np.nan])

                # Cap coordinates
                x_left = q[i] / (1 + cap_width_log)
                x_right = q[i] * (1 + cap_width_log)

                # Top cap
                x_lines.extend([x_left, x_right, np.nan])
                y_lines.extend([y_top, y_top, np.nan])

                # Bottom cap
                x_lines.extend([x_left, x_right, np.nan])
                y_lines.extend([y_bottom, y_bottom, np.nan])

            # Draw all error bars as a single plot item for speed
            error_pen = pg.mkPen((100, 100, 255, 120), width=1)
            self.main_plot.plot(
                x_lines, y_lines,
                pen=error_pen,
                connect='finite'  # Connect all non-NaN points, NaN breaks segments
            )

        # Initialize cursor positions ONLY if not already set (first data load)
        # This preserves user's cursor selections during fitting
        if len(q) > 0 and self.cursor_left is None:
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
        """Plot fit residuals with symmetric Y-axis around 0."""
        self.residual_plot.clear()
        self.residual_plot.addLine(y=0, pen=pg.mkPen('k', style=Qt.PenStyle.DashLine))
        self.residual_plot.plot(
            q, residuals,
            pen=None,
            symbol='o',
            symbolSize=3,
            symbolBrush=(100, 100, 255, 150)
        )

        # Set symmetric Y-axis range around 0
        if len(residuals) > 0:
            max_abs = np.max(np.abs(residuals))
            if max_abs > 0:
                self.residual_plot.setYRange(-max_abs * 1.1, max_abs * 1.1)

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

        # Define colors for each level (Red, Green, Blue, Orange, Purple)
        level_colors = {
            1: '#d32f2f',  # Red
            2: '#388e3c',  # Green
            3: '#1976d2',  # Blue
            4: '#f57c00',  # Orange
            5: '#7b1fa2'   # Purple
        }
        header_color = level_colors.get(self.level_number, '#e74c3c')

        # Level header
        header = QLabel(f"Level {self.level_number}")
        header.setStyleSheet(f"""
            QLabel {{
                background-color: {header_color};
                color: white;
                font-weight: bold;
                font-size: 12px;
                padding: 5px;
            }}
        """)
        layout.addWidget(header)

        # Controls header
        controls_header = QLabel("Controls")
        controls_header.setStyleSheet(f"""
            QLabel {{
                background-color: {header_color};
                color: white;
                font-weight: bold;
                font-size: 11px;
                padding: 3px;
            }}
        """)
        layout.addWidget(controls_header)

        # Parameters grid
        grid = QGridLayout()
        grid.setSpacing(8)

        # Set minimum column widths to prevent layout shift when hiding limit fields
        grid.setColumnMinimumWidth(3, 85)  # Low limit column
        grid.setColumnMinimumWidth(4, 85)  # High limit column

        # Column headers
        grid.addWidget(QLabel(""), 0, 0)
        grid.addWidget(QLabel("Fit?"), 0, 2, Qt.AlignmentFlag.AlignCenter)
        self.low_limit_header = QLabel("Low limit:")
        self.low_limit_header.setAlignment(Qt.AlignmentFlag.AlignRight)
        grid.addWidget(self.low_limit_header, 0, 3)
        self.high_limit_header = QLabel("High Limit:")
        self.high_limit_header.setAlignment(Qt.AlignmentFlag.AlignRight)
        grid.addWidget(self.high_limit_header, 0, 4)

        # G parameter
        row = 1
        grid.addWidget(QLabel("G"), row, 0)
        self.g_value = ScrubbableLineEdit("100")
        self.g_value.setValidator(QDoubleValidator())
        self.g_value.setMinimumWidth(120)
        self.g_value.setMaximumWidth(120)
        self.g_value.editingFinished.connect(self._on_g_changed)
        grid.addWidget(self.g_value, row, 1)
        self.g_fit = QCheckBox()
        grid.addWidget(self.g_fit, row, 2, Qt.AlignmentFlag.AlignCenter)
        self.g_low = QLineEdit("20")
        self.g_low.setMaximumWidth(80)
        self.g_low.setValidator(QDoubleValidator())
        grid.addWidget(self.g_low, row, 3)
        self.g_high = QLineEdit("500")
        self.g_high.setMaximumWidth(80)
        self.g_high.setValidator(QDoubleValidator())
        grid.addWidget(self.g_high, row, 4)

        # Rg parameter
        row = 2
        grid.addWidget(QLabel("Rg"), row, 0)
        self.rg_value = ScrubbableLineEdit("100")
        self.rg_value.setValidator(QDoubleValidator())
        self.rg_value.setMinimumWidth(120)
        self.rg_value.setMaximumWidth(120)
        self.rg_value.editingFinished.connect(self._on_rg_changed)
        grid.addWidget(self.rg_value, row, 1)
        self.rg_fit = QCheckBox()
        grid.addWidget(self.rg_fit, row, 2, Qt.AlignmentFlag.AlignCenter)
        self.rg_low = QLineEdit("20")
        self.rg_low.setMaximumWidth(80)
        self.rg_low.setValidator(QDoubleValidator())
        grid.addWidget(self.rg_low, row, 3)
        self.rg_high = QLineEdit("500")
        self.rg_high.setMaximumWidth(80)
        self.rg_high.setValidator(QDoubleValidator())
        grid.addWidget(self.rg_high, row, 4)

        layout.addLayout(grid)

        # Fit Rg/G button - 35% width, right aligned with high limit
        fit_rg_g_layout = QHBoxLayout()
        fit_rg_g_layout.addStretch()
        self.fit_rg_g_button = QPushButton("Fit Rg/G btwn cursors")
        self.fit_rg_g_button.setMinimumHeight(24)
        # Calculate width to align with high limit field (roughly 35% of total width)
        button_width = 140  # Approximate width for 35%
        self.fit_rg_g_button.setMaximumWidth(button_width)
        fit_rg_g_layout.addWidget(self.fit_rg_g_button)
        layout.addLayout(fit_rg_g_layout)

        # Estimate B checkbox
        estimate_b_layout = QHBoxLayout()
        self.estimate_b_check = QCheckBox("Estimate B from G/Rg/P?")
        self.estimate_b_check.stateChanged.connect(self._on_estimate_b_changed)
        estimate_b_layout.addWidget(self.estimate_b_check)
        estimate_b_layout.addStretch()
        layout.addLayout(estimate_b_layout)

        # B and P parameters
        grid2 = QGridLayout()
        grid2.setSpacing(8)

        # Set minimum column widths to prevent layout shift when hiding limit fields
        grid2.setColumnMinimumWidth(3, 85)  # Low limit column
        grid2.setColumnMinimumWidth(4, 85)  # High limit column

        # B parameter
        row = 0
        grid2.addWidget(QLabel("B"), row, 0)
        self.b_value = ScrubbableLineEdit("0.01")
        self.b_value.setValidator(QDoubleValidator())
        self.b_value.setMinimumWidth(120)
        self.b_value.setMaximumWidth(120)
        self.b_value.editingFinished.connect(self._on_b_changed)
        grid2.addWidget(self.b_value, row, 1)
        self.b_fit = QCheckBox()
        grid2.addWidget(self.b_fit, row, 2, Qt.AlignmentFlag.AlignCenter)
        self.b_low = QLineEdit("0.002")
        self.b_low.setMaximumWidth(80)
        self.b_low.setValidator(QDoubleValidator())
        grid2.addWidget(self.b_low, row, 3)
        self.b_high = QLineEdit("0.05")
        self.b_high.setMaximumWidth(80)
        self.b_high.setValidator(QDoubleValidator())
        grid2.addWidget(self.b_high, row, 4)

        # P parameter
        row = 1
        grid2.addWidget(QLabel("P"), row, 0)
        self.p_value = ScrubbableLineEdit("4", use_fixed_step=True, fixed_reference=4.0)
        self.p_value.setValidator(QDoubleValidator())
        self.p_value.setMinimumWidth(120)
        self.p_value.setMaximumWidth(120)
        self.p_value.editingFinished.connect(self._on_p_changed)
        grid2.addWidget(self.p_value, row, 1)
        self.p_fit = QCheckBox()
        grid2.addWidget(self.p_fit, row, 2, Qt.AlignmentFlag.AlignCenter)
        self.p_low = QLineEdit("2")
        self.p_low.setMaximumWidth(80)
        self.p_low.setValidator(QDoubleValidator())
        grid2.addWidget(self.p_low, row, 3)
        self.p_high = QLineEdit("5")
        self.p_high.setMaximumWidth(80)
        self.p_high.setValidator(QDoubleValidator())
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

        # Fit P/B button - 35% width, right aligned with high limit
        fit_p_b_layout = QHBoxLayout()
        fit_p_b_layout.addStretch()
        self.fit_p_b_button = QPushButton("Fit P/B btwn cursors")
        self.fit_p_b_button.setMinimumHeight(24)
        self.fit_p_b_button.setMaximumWidth(button_width)
        fit_p_b_layout.addWidget(self.fit_p_b_button)
        layout.addLayout(fit_p_b_layout)

        # Link RgCutoff checkbox (only for levels 2-5)
        if self.level_number > 1:
            self.link_rgco_check = QCheckBox("Link RgCutoff")
            self.link_rgco_check.setToolTip(f"Link RgCutoff to Rg of Level {self.level_number - 1}")
            self.link_rgco_check.stateChanged.connect(lambda: self.parameter_changed.emit())
            layout.addWidget(self.link_rgco_check)
        else:
            self.link_rgco_check = None

        # RgCutoff and Correlations checkbox on same row
        rg_cutoff_corr_layout = QHBoxLayout()
        rg_cutoff_corr_layout.addWidget(QLabel("RgCutoff"))
        self.rg_cutoff = ScrubbableLineEdit("0")
        self.rg_cutoff.setValidator(QDoubleValidator())
        self.rg_cutoff.setMaximumWidth(100)
        rg_cutoff_corr_layout.addWidget(self.rg_cutoff)
        rg_cutoff_corr_layout.addSpacing(20)
        self.correlated_check = QCheckBox("Correlations?")
        self.correlated_check.stateChanged.connect(self._on_correlations_changed)
        rg_cutoff_corr_layout.addWidget(self.correlated_check)
        rg_cutoff_corr_layout.addStretch()
        layout.addLayout(rg_cutoff_corr_layout)

        # ETA and PACK parameters (visible only when Correlations is checked)
        self.corr_params_widget = QWidget()
        corr_params_layout = QVBoxLayout()
        corr_params_layout.setContentsMargins(0, 0, 0, 0)
        corr_params_layout.setSpacing(8)

        # Grid for ETA and PACK
        corr_grid = QGridLayout()
        corr_grid.setSpacing(8)

        # Set minimum column widths to prevent layout shift when hiding limit fields
        corr_grid.setColumnMinimumWidth(3, 85)  # Low limit column
        corr_grid.setColumnMinimumWidth(4, 85)  # High limit column

        # ETA parameter
        row = 0
        corr_grid.addWidget(QLabel("ETA"), row, 0)
        self.eta_value = ScrubbableLineEdit("0")
        self.eta_value.setValidator(QDoubleValidator())
        self.eta_value.setMinimumWidth(95)
        self.eta_value.setMaximumWidth(95)
        self.eta_value.editingFinished.connect(self._on_eta_changed)
        corr_grid.addWidget(self.eta_value, row, 1)
        self.eta_fit = QCheckBox()
        corr_grid.addWidget(self.eta_fit, row, 2, Qt.AlignmentFlag.AlignCenter)
        self.eta_low = QLineEdit("0")
        self.eta_low.setMaximumWidth(80)
        self.eta_low.setValidator(QDoubleValidator())
        corr_grid.addWidget(self.eta_low, row, 3)
        self.eta_high = QLineEdit("1")
        self.eta_high.setMaximumWidth(80)
        self.eta_high.setValidator(QDoubleValidator())
        corr_grid.addWidget(self.eta_high, row, 4)

        # PACK parameter
        row = 1
        corr_grid.addWidget(QLabel("PACK"), row, 0)
        self.pack_value = ScrubbableLineEdit("0")
        self.pack_value.setValidator(QDoubleValidator())
        self.pack_value.setMinimumWidth(95)
        self.pack_value.setMaximumWidth(95)
        self.pack_value.editingFinished.connect(self._on_pack_changed)
        corr_grid.addWidget(self.pack_value, row, 1)
        self.pack_fit = QCheckBox()
        corr_grid.addWidget(self.pack_fit, row, 2, Qt.AlignmentFlag.AlignCenter)
        self.pack_low = QLineEdit("0")
        self.pack_low.setMaximumWidth(80)
        self.pack_low.setValidator(QDoubleValidator())
        corr_grid.addWidget(self.pack_low, row, 3)
        self.pack_high = QLineEdit("1")
        self.pack_high.setMaximumWidth(80)
        self.pack_high.setValidator(QDoubleValidator())
        corr_grid.addWidget(self.pack_high, row, 4)

        corr_params_layout.addLayout(corr_grid)
        self.corr_params_widget.setLayout(corr_params_layout)
        self.corr_params_widget.setVisible(False)  # Hidden by default
        layout.addWidget(self.corr_params_widget)

        # Copy/Move/swap button
        self.copy_move_button = QPushButton("Copy/Move/swap level")
        self.copy_move_button.setMinimumHeight(24)  # Reduced by ~20% from 30
        layout.addWidget(self.copy_move_button)

        # Calculated values display - Sv and Invariant
        calc_values_layout = QVBoxLayout()
        calc_values_layout.setSpacing(4)

        # Sv (Surface to Volume ratio)
        sv_layout = QHBoxLayout()
        sv_layout.addWidget(QLabel("Sv [m2/cm3]:"))
        self.sv_value = QLineEdit("0")
        self.sv_value.setReadOnly(True)
        self.sv_value.setMaximumWidth(100)
        sv_layout.addWidget(self.sv_value)
        sv_layout.addStretch()
        calc_values_layout.addLayout(sv_layout)

        # Invariant
        invariant_layout = QHBoxLayout()
        invariant_layout.addWidget(QLabel("Invariant [cm^-4]:"))
        self.invariant_value = QLineEdit("0")
        self.invariant_value.setReadOnly(True)
        self.invariant_value.setMaximumWidth(100)
        invariant_layout.addWidget(self.invariant_value)
        invariant_layout.addStretch()
        calc_values_layout.addLayout(invariant_layout)

        layout.addLayout(calc_values_layout)

        layout.addStretch()
        self.setLayout(layout)

        # Store all limit fields for easy show/hide
        self.limit_fields = [
            self.g_low, self.g_high,
            self.rg_low, self.rg_high,
            self.b_low, self.b_high,
            self.p_low, self.p_high,
            self.eta_low, self.eta_high,
            self.pack_low, self.pack_high
        ]
        self.limit_headers = [self.low_limit_header, self.high_limit_header]

        # Initialize limits based on default values
        self.fix_limits()

    def _format_value(self, value: float) -> str:
        """Format a value to 3 significant digits."""
        if value == 0:
            return "0"
        # Use scientific notation for very small or very large numbers
        if abs(value) < 0.001 or abs(value) >= 10000:
            return f"{value:.2e}"
        else:
            # Use fixed decimal for normal range numbers
            return f"{value:.3g}"

    def toggle_limits_visibility(self, show: bool):
        """Show or hide all limit fields and headers."""
        for field in self.limit_fields:
            # Don't show B limits if estimate_B is checked
            if field in [self.b_low, self.b_high] and self.estimate_b_check.isChecked():
                field.setVisible(False)
            else:
                field.setVisible(show)
        for header in self.limit_headers:
            header.setVisible(show)

    def _on_g_changed(self):
        """Update G limits when G value changes."""
        try:
            value = float(self.g_value.text() or 0)
            # If G is essentially 0, set Rg to 1e10 and uncheck fit boxes
            if value < 1e-10:
                self.rg_value.setText(self._format_value(1e10))
                self.g_fit.setChecked(False)
                self.rg_fit.setChecked(False)
                self.g_value.setText("0")
            elif value > 0:
                self.g_low.setText(self._format_value(value * 0.2))
                self.g_high.setText(self._format_value(value * 5))
                self.g_value.setText(self._format_value(value))

            # Recalculate B if estimate_B is checked
            if self.estimate_b_check.isChecked():
                self._calculate_and_set_b()
        except ValueError:
            pass
        self.parameter_changed.emit()

    def _on_rg_changed(self):
        """Update Rg limits when Rg value changes."""
        try:
            value = float(self.rg_value.text() or 0)
            if value > 0:
                # Clamp to absolute limits [0.1, 1e4]
                low_val = max(0.1, value * 0.2)
                high_val = min(1e4, value * 5)
                self.rg_low.setText(self._format_value(low_val))
                self.rg_high.setText(self._format_value(high_val))
                self.rg_value.setText(self._format_value(value))

            # Recalculate B if estimate_B is checked
            if self.estimate_b_check.isChecked():
                self._calculate_and_set_b()
        except ValueError:
            pass
        self.parameter_changed.emit()

    def _on_b_changed(self):
        """Update B limits when B value changes."""
        try:
            value = float(self.b_value.text() or 0)
            if value > 0:
                self.b_low.setText(self._format_value(value * 0.2))
                self.b_high.setText(self._format_value(value * 5))
                self.b_value.setText(self._format_value(value))
        except ValueError:
            pass
        # Update Sv and Invariant when B changes
        self._update_calculated_values()
        self.parameter_changed.emit()

    def _on_p_changed(self):
        """Update P limits when P value changes."""
        try:
            value = float(self.p_value.text() or 0)
            if value > 0:
                # P has special limits: 0.5-1.5x but clamped to [1, 5]
                low_val = max(1.0, value * 0.5)
                high_val = min(5.0, value * 1.5)
                self.p_low.setText(self._format_value(low_val))
                self.p_high.setText(self._format_value(high_val))
                self.p_value.setText(self._format_value(value))

            # Recalculate B if estimate_B is checked
            if self.estimate_b_check.isChecked():
                self._calculate_and_set_b()
        except ValueError:
            pass
        # Update Sv and Invariant when P changes
        self._update_calculated_values()
        self.parameter_changed.emit()

    def _on_eta_changed(self):
        """Update ETA limits when ETA value changes."""
        try:
            value = float(self.eta_value.text() or 0)
            if value > 0:
                # Same logic as Rg: clamp to absolute limits [0.1, 1e4]
                low_val = max(0.1, value * 0.2)
                high_val = min(1e4, value * 5)
                self.eta_low.setText(self._format_value(low_val))
                self.eta_high.setText(self._format_value(high_val))
                self.eta_value.setText(self._format_value(value))
        except ValueError:
            pass
        self.parameter_changed.emit()

    def _on_pack_changed(self):
        """Update PACK limits when PACK value changes."""
        try:
            value = float(self.pack_value.text() or 0)
            if value > 0:
                # Same fractional logic as ETA/Rg but absolute limits [0, 12]
                low_val = max(0.0, value * 0.2)
                high_val = min(12.0, value * 5)
                self.pack_low.setText(self._format_value(low_val))
                self.pack_high.setText(self._format_value(high_val))
                self.pack_value.setText(self._format_value(value))
        except ValueError:
            pass
        self.parameter_changed.emit()

    def _on_correlations_changed(self, state):
        """Show/hide correlation parameters when checkbox changes."""
        self.corr_params_widget.setVisible(self.correlated_check.isChecked())
        self.parameter_changed.emit()

    def _on_estimate_b_changed(self, state):
        """Handle Estimate B checkbox change."""
        estimate_b = self.estimate_b_check.isChecked()

        if estimate_b:
            # Calculate B from G, Rg, and P
            self._calculate_and_set_b()

            # Uncheck and disable B fit checkbox
            self.b_fit.setChecked(False)
            self.b_fit.setEnabled(False)

            # Hide B limit fields
            self.b_low.setVisible(False)
            self.b_high.setVisible(False)
        else:
            # Re-enable B fit checkbox
            self.b_fit.setEnabled(True)

            # Show B limit fields (unless "No limits?" is checked globally)
            # This will be handled by the parent's toggle_limits_visibility if needed
            self.b_low.setVisible(True)
            self.b_high.setVisible(True)

        self.parameter_changed.emit()

    def _calculate_and_set_b(self):
        """Calculate B from formula: B = G * exp(-P/2) * (3*P/2)^(P/2) * (1/Rg^P)"""
        try:
            import numpy as np

            G = float(self.g_value.text() or 0)
            Rg = float(self.rg_value.text() or 0)
            P = float(self.p_value.text() or 0)

            if G > 0 and Rg > 0 and P > 0:
                # B = G * exp(-P/2) * (3*P/2)^(P/2) * (1/Rg^P)
                B = G * np.exp(-P/2.0) * ((3.0 * P / 2.0) ** (P / 2.0)) * (1.0 / Rg ** P)

                # Set B value
                self.b_value.setText(self._format_value(B))

                # Update B limits based on new value
                self.b_low.setText(self._format_value(B / 5.0))
                self.b_high.setText(self._format_value(B * 5.0))
        except (ValueError, ZeroDivisionError):
            pass

    def _update_calculated_values(self):
        """Trigger calculation update - needs to be called from parent with Q vector."""
        # This will be called by the parent UnifiedFitPanel when it has Q data
        # For now, just emit signal to notify parent
        pass

    def update_porod_surface_and_invariant(self, q_vector):
        """
        Calculate and update Porod surface (Sv) and Invariant for this level.
        Based on IR1A_UpdatePorodSfcandInvariant from Igor code.

        Args:
            q_vector: Q vector from the experimental data
        """
        try:
            import numpy as np
            try:
                from scipy.integrate import simpson
            except ImportError:
                from scipy.integrate import simps as simpson
            from scipy.special import erf, gamma

            # Get parameters
            P = float(self.p_value.text() or 0)
            B = float(self.b_value.text() or 0)
            Rg = float(self.rg_value.text() or 0)
            G = float(self.g_value.text() or 0)
            RgCO = float(self.rg_cutoff.text() or 0)
            ETA = float(self.eta_value.text() or 0)
            PACK = float(self.pack_value.text() or 0)
            correlated = self.correlated_check.isChecked()

            if Rg <= 0:
                self.sv_value.setText("N/A")
                self.invariant_value.setText("N/A")
                return

            # Create Q vector for integration: maxQ = 2*pi / (Rg/10)
            maxQ = 2 * np.pi / (Rg / 10)
            npoints = 2000
            surf_q = np.linspace(0, maxQ, npoints)

            # Calculate intensity using unified model
            surf_int = self._calculate_unified_intensity(surf_q, G, Rg, P, B, RgCO, ETA, PACK, correlated)

            # Handle Q=0 point (use Q=dQ value)
            if surf_int[0] == 0 or np.isnan(surf_int[0]):
                surf_int[0] = surf_int[1]

            # Calculate invariant integrand: I * Q^2
            surf_invariant = surf_int * surf_q**2

            # Integrate to get invariant (using Simpson's rule)
            invariant = simpson(surf_invariant, x=surf_q)

            # Add Porod tail if RgCO < 0.1 (no cutoff)
            if RgCO < 0.1:
                # Porod tail contribution: -B * maxQ^(3-P) / (3-P)
                invariant += -B * maxQ**(3 - abs(P)) / (3 - abs(P))

            # Check if invariant is valid (negative means bad extrapolation)
            if invariant < 0:
                self.invariant_value.setText("N/A")
                self.sv_value.setText("N/A")
                return

            # Calculate surface to volume ratio (Sv = pi*B/Q) BEFORE scaling invariant
            # Only valid when P is close to 4 (Porod regime)
            # At this point invariant is still in A^-3 * cm^-1
            if 3.95 <= P <= 4.05:
                sv = 1e4 * np.pi * B / invariant
                self.sv_value.setText(self._format_value(sv))
            else:
                self.sv_value.setText("N/A")

            # Convert invariant from A^-3 * cm^-1 to cm^-4
            invariant = invariant * 1e24

            # Update invariant display
            self.invariant_value.setText(self._format_value(invariant))

        except Exception as e:
            print(f"Error calculating Sv and Invariant: {e}")
            import traceback
            traceback.print_exc()
            self.sv_value.setText("Error")
            self.invariant_value.setText("Error")

    def _calculate_unified_intensity(self, q, G, Rg, P, B, RgCO, ETA, PACK, correlated):
        """
        Calculate unified model intensity for given Q vector.
        Based on IR1A_SurfToVolCalcInvarVec from Igor code.
        """
        import numpy as np
        from scipy.special import erf, gamma

        # Calculate K value
        K = 1.0 if P > 3 else 1.06

        # Calculate Q* (erf correction)
        qstar = q / (erf(K * q * Rg / np.sqrt(6)))**3

        # Avoid division by zero
        qstar = np.where(qstar == 0, 1e-10, qstar)

        # Calculate unified intensity
        intensity = G * np.exp(-q**2 * Rg**2 / 3) + (B / qstar**P) * np.exp(-RgCO**2 * q**2 / 3)

        # Apply correlation correction if enabled
        if correlated and PACK > 0 and ETA > 0:
            # Calculate sphere amplitude (hard sphere structure factor)
            sphere_amp = self._sphere_amplitude(q, ETA)
            intensity = intensity / (1 + PACK * sphere_amp)

        return intensity

    def _sphere_amplitude(self, q, eta):
        """Calculate sphere amplitude for structure factor correction."""
        import numpy as np

        # Based on hard sphere structure factor
        # This is a simplified version - full implementation in Igor is more complex
        qr = q * eta

        # Avoid division by zero
        qr = np.where(qr == 0, 1e-10, qr)

        # Sphere form factor amplitude
        amplitude = 3 * (np.sin(qr) - qr * np.cos(qr)) / qr**3

        return amplitude

    def fix_limits(self):
        """
        Fix all fitting limits based on current parameter values.
        Called when limits get out of sync with values (e.g., after fitting).
        """
        # Fix G limits
        try:
            g_val = float(self.g_value.text() or 0)
            if g_val < 1e-10:
                # G is essentially 0 - set Rg to 1e10 and uncheck fit boxes
                self.rg_value.setText(self._format_value(1e10))
                self.g_fit.setChecked(False)
                self.rg_fit.setChecked(False)
                self.g_value.setText("0")
            elif g_val > 0:
                self.g_low.setText(self._format_value(g_val * 0.2))
                self.g_high.setText(self._format_value(g_val * 5))
        except ValueError:
            pass

        # Fix Rg limits - absolute limits [0.1, 1e4]
        try:
            rg_val = float(self.rg_value.text() or 0)
            if rg_val > 0:
                low_val = max(0.1, rg_val * 0.2)
                high_val = min(1e4, rg_val * 5)
                self.rg_low.setText(self._format_value(low_val))
                self.rg_high.setText(self._format_value(high_val))
        except ValueError:
            pass

        # Fix B limits
        try:
            b_val = float(self.b_value.text() or 0)
            if b_val > 0:
                self.b_low.setText(self._format_value(b_val * 0.2))
                self.b_high.setText(self._format_value(b_val * 5))
        except ValueError:
            pass

        # Fix P limits - 0.5-1.5x, clamped to [1, 5]
        try:
            p_val = float(self.p_value.text() or 0)
            if p_val > 0:
                low_val = max(1.0, p_val * 0.5)
                high_val = min(5.0, p_val * 1.5)
                self.p_low.setText(self._format_value(low_val))
                self.p_high.setText(self._format_value(high_val))
        except ValueError:
            pass

        # Fix ETA limits - same as Rg: absolute limits [0.1, 1e4]
        try:
            eta_val = float(self.eta_value.text() or 0)
            if eta_val > 0:
                low_val = max(0.1, eta_val * 0.2)
                high_val = min(1e4, eta_val * 5)
                self.eta_low.setText(self._format_value(low_val))
                self.eta_high.setText(self._format_value(high_val))
        except ValueError:
            pass

        # Fix PACK limits - absolute limits [0, 12]
        try:
            pack_val = float(self.pack_value.text() or 0)
            if pack_val > 0:
                low_val = max(0.0, pack_val * 0.2)
                high_val = min(12.0, pack_val * 5)
                self.pack_low.setText(self._format_value(low_val))
                self.pack_high.setText(self._format_value(high_val))
        except ValueError:
            pass

    def get_parameters(self) -> Dict:
        """Get all parameters for this level."""
        return {
            'G': float(self.g_value.text() or 0),
            'Rg': float(self.rg_value.text() or 0),
            'B': float(self.b_value.text() or 0),
            'P': float(self.p_value.text() or 0),
            'RgCutoff': float(self.rg_cutoff.text() or 0),
            'ETA': float(self.eta_value.text() or 0),
            'PACK': float(self.pack_value.text() or 0),
            'fit_G': self.g_fit.isChecked(),
            'fit_Rg': self.rg_fit.isChecked(),
            'fit_B': self.b_fit.isChecked(),
            'fit_P': self.p_fit.isChecked(),
            'fit_ETA': self.eta_fit.isChecked(),
            'fit_PACK': self.pack_fit.isChecked(),
            'estimate_B': self.estimate_b_check.isChecked(),
            'correlated': self.correlated_check.isChecked(),
            'link_rgco': self.link_rgco_check.isChecked() if self.link_rgco_check else False,
            'G_low': float(self.g_low.text() or 0),
            'G_high': float(self.g_high.text() or 0),
            'Rg_low': float(self.rg_low.text() or 0),
            'Rg_high': float(self.rg_high.text() or 0),
            'B_low': float(self.b_low.text() or 0),
            'B_high': float(self.b_high.text() or 0),
            'P_low': float(self.p_low.text() or 0),
            'P_high': float(self.p_high.text() or 0),
            'ETA_low': float(self.eta_low.text() or 0),
            'ETA_high': float(self.eta_high.text() or 0),
            'PACK_low': float(self.pack_low.text() or 0),
            'PACK_high': float(self.pack_high.text() or 0),
        }

    def set_parameters(self, params: Dict):
        """Set parameters for this level."""
        if 'G' in params:
            self.g_value.setText(self._format_value(params['G']))
        if 'Rg' in params:
            self.rg_value.setText(self._format_value(params['Rg']))
        if 'B' in params:
            self.b_value.setText(self._format_value(params['B']))
        if 'P' in params:
            self.p_value.setText(self._format_value(params['P']))
        if 'RgCutoff' in params:
            self.rg_cutoff.setText(self._format_value(params['RgCutoff']))
        if 'G_low' in params:
            self.g_low.setText(self._format_value(params['G_low']))
        if 'G_high' in params:
            self.g_high.setText(self._format_value(params['G_high']))
        if 'Rg_low' in params:
            self.rg_low.setText(self._format_value(params['Rg_low']))
        if 'Rg_high' in params:
            self.rg_high.setText(self._format_value(params['Rg_high']))
        if 'B_low' in params:
            self.b_low.setText(self._format_value(params['B_low']))
        if 'B_high' in params:
            self.b_high.setText(self._format_value(params['B_high']))
        if 'P_low' in params:
            self.p_low.setText(self._format_value(params['P_low']))
        if 'P_high' in params:
            self.p_high.setText(self._format_value(params['P_high']))
        if 'ETA' in params:
            self.eta_value.setText(self._format_value(params['ETA']))
        if 'PACK' in params:
            self.pack_value.setText(self._format_value(params['PACK']))
        if 'ETA_low' in params:
            self.eta_low.setText(self._format_value(params['ETA_low']))
        if 'ETA_high' in params:
            self.eta_high.setText(self._format_value(params['ETA_high']))
        if 'PACK_low' in params:
            self.pack_low.setText(self._format_value(params['PACK_low']))
        if 'PACK_high' in params:
            self.pack_high.setText(self._format_value(params['PACK_high']))
        if 'fit_G' in params:
            self.g_fit.setChecked(params['fit_G'])
        if 'fit_Rg' in params:
            self.rg_fit.setChecked(params['fit_Rg'])
        if 'fit_B' in params:
            self.b_fit.setChecked(params['fit_B'])
        if 'fit_P' in params:
            self.p_fit.setChecked(params['fit_P'])
        if 'fit_ETA' in params:
            self.eta_fit.setChecked(params['fit_ETA'])
        if 'fit_PACK' in params:
            self.pack_fit.setChecked(params['fit_PACK'])
        if 'estimate_B' in params:
            self.estimate_b_check.setChecked(params['estimate_B'])
            # Update B fit checkbox and limits visibility
            if params['estimate_B']:
                self.b_fit.setChecked(False)
                self.b_fit.setEnabled(False)
                self.b_low.setVisible(False)
                self.b_high.setVisible(False)
            else:
                self.b_fit.setEnabled(True)
                self.b_low.setVisible(True)
                self.b_high.setVisible(True)
        if 'correlated' in params:
            self.correlated_check.setChecked(params['correlated'])
            # Update visibility of correlation parameters
            self.corr_params_widget.setVisible(params['correlated'])
        if 'link_rgco' in params and self.link_rgco_check:
            self.link_rgco_check.setChecked(params['link_rgco'])


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

        # Storage for local fit curves
        self.local_fits = {}  # Dict to store local fit curves: {level: {'guinier': (q, I), 'porod': (q, I)}}

        # State management
        self.state_manager = StateManager()

        # Backup for revert functionality
        self.parameter_backup = None

        self.init_ui()
        self.load_state()

    def init_ui(self):
        """Initialize the user interface."""
        self.setWindowTitle("pyIrena - Unified Fit Model")

        # Main horizontal splitter
        main_splitter = QSplitter(Qt.Orientation.Horizontal)

        # Left panel (controls) - maintain 33% width ratio
        left_panel = self.create_control_panel()
        main_splitter.addWidget(left_panel)

        # Right panel (graph) - maintain 67% width ratio
        self.graph_window = UnifiedFitGraphWindow()
        main_splitter.addWidget(self.graph_window)

        # Set initial sizes (1:2 ratio = 33%:67%) and stretch factors to maintain ratio
        main_splitter.setSizes([400, 800])
        main_splitter.setStretchFactor(0, 1)  # Left panel: 1 part (33%)
        main_splitter.setStretchFactor(1, 2)  # Right panel: 2 parts (67%)

        # Store splitter reference to maintain ratio
        self.main_splitter = main_splitter

        # Main layout
        main_layout = QVBoxLayout()
        main_layout.addWidget(main_splitter)
        self.setLayout(main_layout)

        # Set minimum and initial window size (taller for better field visibility)
        self.setMinimumSize(1200, 960)  # Same width, 20% taller (800 * 1.2 = 960)
        self.resize(1200, 960)  # Set initial size

    def format_value_3sig(self, value: float) -> str:
        """Format a value to 3 significant digits for display."""
        if value == 0:
            return "0"
        # Use scientific notation with 2 decimal places (3 sig figs total)
        if abs(value) < 0.01 or abs(value) >= 1000:
            return f"{value:.2e}"
        else:
            # For values in normal range, use 3 significant figures
            return f"{value:.3g}"

    def create_control_panel(self) -> QWidget:
        """Create the left control panel."""
        panel = QWidget()

        # Set size policy to prevent content-driven expansion
        try:
            from PySide6.QtWidgets import QSizePolicy
        except ImportError:
            try:
                from PyQt6.QtWidgets import QSizePolicy
            except ImportError:
                from PyQt5.QtWidgets import QSizePolicy

        size_policy = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Preferred)
        panel.setSizePolicy(size_policy)
        panel.setMinimumWidth(400)
        panel.setMaximumWidth(400)

        layout = QVBoxLayout()
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)  # Reduced from 10 to 6 to save vertical space

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

        # Top controls row - Number of levels and No limits
        top_controls = QHBoxLayout()
        top_controls.addWidget(QLabel("Number of levels:"))
        self.num_levels_spin = QSpinBox()
        self.num_levels_spin.setMinimum(1)
        self.num_levels_spin.setMaximum(5)
        self.num_levels_spin.setValue(1)
        self.num_levels_spin.setMinimumHeight(26)
        self.num_levels_spin.valueChanged.connect(self.on_num_levels_changed)
        top_controls.addWidget(self.num_levels_spin)

        top_controls.addSpacing(20)

        self.no_limits_check = QCheckBox("No limits?")
        self.no_limits_check.stateChanged.connect(self.on_no_limits_changed)
        top_controls.addWidget(self.no_limits_check)
        top_controls.addStretch()
        layout.addLayout(top_controls)

        # Level tabs
        self.level_tabs = QTabWidget()
        self.level_widgets = []

        # Create widgets for all 5 possible levels
        for i in range(1, 6):
            level_widget = LevelParametersWidget(i)
            # Connect parameter changed signal
            level_widget.parameter_changed.connect(self.on_parameter_changed)
            # Connect fit buttons
            level_widget.fit_rg_g_button.clicked.connect(lambda checked, level=i: self.fit_local_guinier(level))
            level_widget.fit_p_b_button.clicked.connect(lambda checked, level=i: self.fit_local_porod(level))
            self.level_widgets.append(level_widget)
            self.level_tabs.addTab(level_widget, f"{i}. Level")

        # Apply stylesheet for tab colors (Level 1: Red, Level 2: Green, etc.)
        tab_stylesheet = """
            QTabBar::tab:nth-child(1) { background-color: #d32f2f; color: white; }
            QTabBar::tab:nth-child(2) { background-color: #388e3c; color: white; }
            QTabBar::tab:nth-child(3) { background-color: #1976d2; color: white; }
            QTabBar::tab:nth-child(4) { background-color: #f57c00; color: white; }
            QTabBar::tab:nth-child(5) { background-color: #7b1fa2; color: white; }
            QTabBar::tab:selected { border: 2px solid #ffd700; font-weight: bold; }
            QTabBar::tab { padding: 6px 10px; }
        """
        self.level_tabs.setStyleSheet(tab_stylesheet)

        # Initially disable unused levels
        self.update_level_tabs()

        # Set reasonable minimum height for tab widget
        self.level_tabs.setMinimumHeight(350)  # Balanced with other layout improvements

        layout.addWidget(self.level_tabs)

        # SAS Background
        background_layout = QHBoxLayout()
        background_layout.addWidget(QLabel("SAS Background"))
        self.background_value = ScrubbableLineEdit("1e-06")
        self.background_value.setValidator(QDoubleValidator())
        self.background_value.setMaximumWidth(100)
        self.background_value.editingFinished.connect(self.on_parameter_changed)
        background_layout.addWidget(self.background_value)
        self.fit_background_check = QCheckBox("Fit Bckg?")
        background_layout.addWidget(self.fit_background_check)
        background_layout.addStretch()
        layout.addLayout(background_layout)

        # Checkboxes row - Update automatically and Display local fits
        checkboxes_layout = QHBoxLayout()
        self.update_auto_check = QCheckBox("Update automatically?")
        checkboxes_layout.addWidget(self.update_auto_check)
        checkboxes_layout.addSpacing(20)
        self.display_local_check = QCheckBox("Display local fits?")
        self.display_local_check.stateChanged.connect(self.on_display_local_changed)
        checkboxes_layout.addWidget(self.display_local_check)
        checkboxes_layout.addStretch()
        layout.addLayout(checkboxes_layout)

        # Fit buttons row - Graph Unified (lighter green), Fit (darker green), Revert back (orange)
        fit_buttons = QHBoxLayout()

        self.graph_unified_button = QPushButton("Graph Unified")
        self.graph_unified_button.setMinimumHeight(28)
        self.graph_unified_button.setMaximumWidth(120)
        self.graph_unified_button.setStyleSheet("""
            QPushButton {
                background-color: #52c77a;
                color: white;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #3eb56a;
            }
        """)
        self.graph_unified_button.clicked.connect(self.graph_unified)
        fit_buttons.addWidget(self.graph_unified_button)

        self.fit_button = QPushButton("Fit")
        self.fit_button.setMinimumHeight(28)
        self.fit_button.setMaximumWidth(120)
        self.fit_button.setStyleSheet("""
            QPushButton {
                background-color: #27ae60;
                color: white;
                font-weight: bold;
                font-size: 13px;
            }
            QPushButton:hover {
                background-color: #1e8449;
            }
        """)
        self.fit_button.clicked.connect(self.run_fit)
        fit_buttons.addWidget(self.fit_button)

        self.revert_button = QPushButton("Revert back")
        self.revert_button.setMinimumHeight(28)
        self.revert_button.setMaximumWidth(120)
        self.revert_button.setStyleSheet("""
            QPushButton {
                background-color: #e67e22;
                color: white;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #f39c12;
            }
        """)
        self.revert_button.clicked.connect(self.revert_to_backup)
        fit_buttons.addWidget(self.revert_button)

        fit_buttons.addStretch()
        layout.addLayout(fit_buttons)

        # Additional buttons row with Reset
        additional_buttons = QHBoxLayout()

        self.reset_button = QPushButton("Reset to Defaults")
        self.reset_button.setMinimumHeight(26)  # Reduced from 30
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
        self.reset_unif_button.setMinimumHeight(26)  # Reduced from 30
        additional_buttons.addWidget(self.reset_unif_button)

        self.fix_limits_button = QPushButton("Fix limits?")
        self.fix_limits_button.setMinimumHeight(26)  # Reduced from 30
        self.fix_limits_button.setStyleSheet("""
            QPushButton {
                background-color: #27ae60;
                color: white;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #2ecc71;
            }
        """)
        self.fix_limits_button.clicked.connect(self.fix_all_limits)
        additional_buttons.addWidget(self.fix_limits_button)

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

        self.save_state_button = QPushButton("Save State")
        self.save_state_button.setMinimumHeight(26)
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

        self.store_data_button = QPushButton("Store in File")
        self.store_data_button.setMinimumHeight(26)
        results_buttons1.addWidget(self.store_data_button)

        self.results_graphs_button = QPushButton("Results to graphs")
        self.results_graphs_button.setMinimumHeight(26)
        results_buttons1.addWidget(self.results_graphs_button)

        layout.addLayout(results_buttons1)

        # Results buttons row 2
        results_buttons2 = QHBoxLayout()

        self.export_params_button = QPushButton("Export Parameters")
        self.export_params_button.setMinimumHeight(26)
        self.export_params_button.clicked.connect(self.export_parameters)
        results_buttons2.addWidget(self.export_params_button)

        self.import_params_button = QPushButton("Import Parameters")
        self.import_params_button.setMinimumHeight(26)
        self.import_params_button.clicked.connect(self.import_parameters)
        results_buttons2.addWidget(self.import_params_button)

        layout.addLayout(results_buttons2)

        # Results buttons row 3
        results_buttons3 = QHBoxLayout()
        self.analyze_results_button = QPushButton("Analyze Results")
        self.analyze_results_button.setMinimumHeight(30)
        results_buttons3.addWidget(self.analyze_results_button)

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

        # Recalculate if "Update automatically" is checked
        if self.update_auto_check.isChecked() and self.data is not None:
            self.graph_unified()

    def update_level_tabs(self):
        """Update which level tabs are enabled."""
        num_levels = self.num_levels_spin.value()
        for i in range(5):
            self.level_tabs.setTabEnabled(i, i < num_levels)

    def on_parameter_changed(self):
        """Called when any parameter changes. Auto-update if enabled."""
        # Sync RgCutoff links whenever parameters change
        self.sync_rgcutoff_links()

        if self.update_auto_check.isChecked() and self.data is not None:
            self.graph_unified()

    def on_no_limits_changed(self, state):
        """Handle change in 'No limits?' checkbox."""
        show_limits = not self.no_limits_check.isChecked()
        # Toggle visibility of limit fields in all level widgets
        for level_widget in self.level_widgets:
            level_widget.toggle_limits_visibility(show_limits)

    def on_display_local_changed(self, state):
        """Handle change in 'Display local fits?' checkbox."""
        if self.display_local_check.isChecked():
            # Show local fits if any exist
            if self.local_fits and self.data is not None:
                # Redraw the graph with local fits
                self.graph_unified()
        else:
            # Hide local fits by redrawing without them
            if self.data is not None:
                self.graph_unified()

    def set_data(self, q, intensity, error=None, label='Data'):
        """Set the data to be fitted."""
        self.data = {
            'Q': q,
            'Intensity': intensity,
            'Error': error,
            'label': label
        }

        # Clear local fits when new data is loaded (they would be invalid for different data)
        self.clear_local_fits()

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
            self.graph_window.show_error_message("No data loaded. Please load data first.")
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
                    ETA=params['ETA'],
                    PACK=params['PACK'],
                    correlations=params['correlated'],
                    Rg_limits=(params['Rg_low'], params['Rg_high']),
                    G_limits=(params['G_low'], params['G_high']),
                    P_limits=(params['P_low'], params['P_high']),
                    B_limits=(params['B_low'], params['B_high']),
                    ETA_limits=(params['ETA_low'], params['ETA_high']),
                    PACK_limits=(params['PACK_low'], params['PACK_high'])
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

            # Plot local fits if checkbox is enabled and there are local fits to display
            if self.display_local_check.isChecked() and self.local_fits:
                self.plot_local_fits()

            # Update Sv and Invariant for all active levels
            for i in range(num_levels):
                self.level_widgets[i].update_porod_surface_and_invariant(self.data['Q'])

            self.status_label.setText(f"Calculated unified fit with {num_levels} level(s)")

        except Exception as e:
            self.graph_window.show_error_message(f"Error calculating unified fit: {str(e)}")
            import traceback
            traceback.print_exc()

    def run_fit(self):
        """Run the unified fit."""
        if self.data is None:
            self.graph_window.show_error_message("No data loaded. Please load data first.")
            return

        # Sync RgCutoff links before fitting
        self.sync_rgcutoff_links()

        # Backup current parameters before fitting
        self.backup_parameters()

        try:
            # Get parameters
            num_levels = self.num_levels_spin.value()
            levels = []

            # Check if limits should be applied
            no_limits = self.no_limits_check.isChecked()

            for i in range(num_levels):
                params = self.level_widgets[i].get_parameters()

                # If "No limits?" is checked, use very wide bounds instead of user-specified limits
                if no_limits:
                    level = UnifiedLevel(
                        Rg=params['Rg'],
                        G=params['G'],
                        P=params['P'],
                        B=params['B'],
                        RgCO=params['RgCutoff'],
                        ETA=params['ETA'],
                        PACK=params['PACK'],
                        correlations=params['correlated'],
                        fit_Rg=params['fit_Rg'],
                        fit_G=params['fit_G'],
                        fit_P=params['fit_P'],
                        fit_B=params['fit_B'],
                        fit_ETA=params['fit_ETA'],
                        fit_PACK=params['fit_PACK'],
                        Rg_limits=(0.1, 1e6),      # Default wide bounds
                        G_limits=(1e-10, 1e10),    # Default wide bounds
                        P_limits=(0.0, 6.0),       # Default wide bounds
                        B_limits=(1e-20, 1e10),    # Default wide bounds
                        ETA_limits=(0.1, 1e6),     # Default wide bounds
                        PACK_limits=(0.0, 16.0)    # Default wide bounds
                    )
                else:
                    level = UnifiedLevel(
                        Rg=params['Rg'],
                        G=params['G'],
                        P=params['P'],
                        B=params['B'],
                        RgCO=params['RgCutoff'],
                        ETA=params['ETA'],
                        PACK=params['PACK'],
                        correlations=params['correlated'],
                        fit_Rg=params['fit_Rg'],
                        fit_G=params['fit_G'],
                        fit_P=params['fit_P'],
                        fit_B=params['fit_B'],
                        fit_ETA=params['fit_ETA'],
                        fit_PACK=params['fit_PACK'],
                        Rg_limits=(params['Rg_low'], params['Rg_high']),
                        G_limits=(params['G_low'], params['G_high']),
                        P_limits=(params['P_low'], params['P_high']),
                        B_limits=(params['B_low'], params['B_high']),
                        ETA_limits=(params['ETA_low'], params['ETA_high']),
                        PACK_limits=(params['PACK_low'], params['PACK_high'])
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
                    'P': fitted_level.P,
                    'ETA': fitted_level.ETA,
                    'PACK': fitted_level.PACK
                })
                # Fix limits after updating fitted values
                self.level_widgets[i].fix_limits()

            self.background_value.setText(self.format_value_3sig(result['background']))

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

            cursor_info = f" | Q: {q_min:.3e} - {q_max:.3e} ({num_points} pts)" if cursor_range else ""

            self.graph_window.show_success_message(
                f"Fit completed successfully! "
                f"Levels: {num_levels} | "
                f"χ²: {chi2:.4f} | "
                f"Reduced χ²: {reduced_chi2:.4f}{cursor_info}"
            )

            # Update Sv and Invariant for all fitted levels
            for i in range(num_levels):
                self.level_widgets[i].update_porod_surface_and_invariant(self.data['Q'])

        except Exception as e:
            self.graph_window.show_error_message(f"Error during fitting: {str(e)}")
            self.status_label.setText("Fit failed")
            import traceback
            traceback.print_exc()

    def fit_local_guinier(self, level):
        """
        Fit Rg and G for a single level using data between cursors.
        This is a local Guinier fit that ignores limits.

        Based on IR1A_FitLocalGuinier from Igor Pro code.
        """
        if self.data is None:
            self.graph_window.show_error_message("No data loaded. Please load data first.")
            return

        # Get the level widget
        level_widget = self.level_widgets[level - 1]
        params = level_widget.get_parameters()

        # Check that at least one parameter is selected for fitting
        if not params['fit_G'] and not params['fit_Rg']:
            self.graph_window.show_error_message(
                f"No fitting parameters selected for Level {level}. "
                "Please check 'Fit?' for G and/or Rg before fitting."
            )
            return

        # Get cursor range
        cursor_range = self.graph_window.get_cursor_range()
        if cursor_range is None:
            self.graph_window.show_error_message(
                "Both cursors must be set on the graph. "
                "Use the cursors to select the Q range for Guinier fit."
            )
            return

        try:
            q_min, q_max = cursor_range

            # Filter data to cursor range
            q = self.data['Q']
            intensity = self.data['Intensity']
            error = self.data.get('Error')

            mask = (q >= q_min) & (q <= q_max)
            q_fit = q[mask]
            intensity_fit = intensity[mask]

            if len(q_fit) < 3:
                self.graph_window.show_error_message(
                    f"Not enough data points between cursors ({len(q_fit)} points). "
                    "Need at least 3 points for fitting."
                )
                return

            # Estimate starting parameters from cursor range
            # LocalRg = 2*pi / Q_avg (from Igor code line 494)
            q_avg = (q_fit[0] + q_fit[-1]) / 2
            local_rg = 2 * np.pi / q_avg

            # LocalG = I_avg (from Igor code line 495)
            local_g = (intensity_fit[0] + intensity_fit[-1]) / 2

            # If not fitting a parameter, use current GUI value
            if not params['fit_G']:
                local_g = params['G']
            if not params['fit_Rg']:
                local_rg = params['Rg']

            # Define Guinier model: I(q) = G * exp(-q^2 * Rg^2 / 3)
            def guinier_model(q, G, Rg):
                return G * np.exp(-q**2 * Rg**2 / 3)

            # Prepare parameters for fitting
            from scipy.optimize import curve_fit

            # Initial guess
            p0 = [local_g, local_rg]

            # Setup which parameters to fit
            # If a parameter is not being fit, we need to fix it
            if not params['fit_G'] and not params['fit_Rg']:
                # Both fixed - already handled above, shouldn't get here
                return
            elif not params['fit_G']:
                # Fix G, fit only Rg
                def model_fixed_g(q, Rg):
                    return guinier_model(q, local_g, Rg)
                p0_fit = [local_rg]
                popt, pcov = curve_fit(model_fixed_g, q_fit, intensity_fit, p0=p0_fit)
                fitted_g = local_g
                fitted_rg = abs(popt[0])
            elif not params['fit_Rg']:
                # Fix Rg, fit only G
                def model_fixed_rg(q, G):
                    return guinier_model(q, G, local_rg)
                p0_fit = [local_g]
                popt, pcov = curve_fit(model_fixed_rg, q_fit, intensity_fit, p0=p0_fit)
                fitted_g = abs(popt[0])
                fitted_rg = local_rg
            else:
                # Fit both G and Rg
                popt, pcov = curve_fit(guinier_model, q_fit, intensity_fit, p0=p0)
                fitted_g = abs(popt[0])
                fitted_rg = abs(popt[1])

            # Update GUI with fitted values
            level_widget.set_parameters({
                'G': fitted_g,
                'Rg': fitted_rg,
                'G_low': fitted_g / 5,
                'G_high': fitted_g * 5,
                'Rg_low': max(0.1, fitted_rg / 5),
                'Rg_high': min(1e4, fitted_rg * 5)
            })

            # Fix limits for this level
            level_widget.fix_limits()

            # Store local fit curve for plotting if display is enabled
            # Calculate Guinier curve over the Q range used for fitting
            guinier_calc = fitted_g * np.exp(-q_fit**2 * fitted_rg**2 / 3)

            # Store in local_fits dictionary
            if level not in self.local_fits:
                self.local_fits[level] = {}
            self.local_fits[level]['guinier'] = (q_fit, guinier_calc)

            # Recalculate and update plot (this will also plot local fits if checkbox is enabled)
            if self.update_auto_check.isChecked():
                self.graph_unified()
            elif self.display_local_check.isChecked():
                # If not auto-updating, manually redraw the graph with local fits
                self.graph_unified()

            # Show success message
            self.status_label.setText(
                f"Local Guinier fit complete for Level {level}: "
                f"G = {fitted_g:.3e}, Rg = {fitted_rg:.3e} "
                f"(Q: {q_min:.3e} - {q_max:.3e}, {len(q_fit)} pts)"
            )

            self.graph_window.show_success_message(
                f"Local Guinier fit completed for Level {level}! "
                f"G = {fitted_g:.4e}, Rg = {fitted_rg:.4e} | "
                f"Q: {q_min:.3e} - {q_max:.3e} ({len(q_fit)} pts) | "
                f"Limits updated automatically."
            )

        except Exception as e:
            self.graph_window.show_error_message(f"Error during local Guinier fit: {str(e)}")
            self.status_label.setText("Fit failed")
            import traceback
            traceback.print_exc()

    def fit_local_porod(self, level):
        """
        Fit P and B for a single level using data between cursors.
        This is a local Porod/power law fit that ignores limits.

        Based on IR1A_FitLocalPorod from Igor Pro code.
        """
        if self.data is None:
            self.graph_window.show_error_message("No data loaded. Please load data first.")
            return

        # Get the level widget
        level_widget = self.level_widgets[level - 1]
        params = level_widget.get_parameters()

        # Check that at least one parameter is selected for fitting
        if not params['fit_B'] and not params['fit_P']:
            self.graph_window.show_error_message(
                f"No fitting parameters selected for Level {level}. "
                "Please check 'Fit?' for B and/or P before fitting."
            )
            return

        # Get cursor range
        cursor_range = self.graph_window.get_cursor_range()
        if cursor_range is None:
            self.graph_window.show_error_message(
                "Both cursors must be set on the graph. "
                "Use the cursors to select the Q range for Porod fit."
            )
            return

        try:
            q_min, q_max = cursor_range

            # Filter data to cursor range
            q = self.data['Q']
            intensity = self.data['Intensity']

            mask = (q >= q_min) & (q <= q_max)
            q_fit = q[mask]
            intensity_fit = intensity[mask]

            if len(q_fit) < 3:
                self.graph_window.show_error_message(
                    f"Not enough data points between cursors ({len(q_fit)} points). "
                    "Need at least 3 points for fitting."
                )
                return

            # Estimate starting parameters from cursor range
            # P (slope) from log-log slope between cursors (Igor line 368)
            # P = abs((log(I_A) - log(I_B)) / (log(Q_B) - log(Q_A)))
            # Using first and last points in range
            local_p = abs(
                (np.log(intensity_fit[0]) - np.log(intensity_fit[-1])) /
                (np.log(q_fit[-1]) - np.log(q_fit[0]))
            )

            # B (prefactor) from I * Q^P at first cursor (Igor line 376)
            local_b = intensity_fit[0] * (q_fit[0] ** local_p)

            # If not fitting a parameter, use current GUI value
            if not params['fit_P']:
                local_p = params['P']
            if not params['fit_B']:
                local_b = params['B']

            # Define power law model: I(q) = B * q^(-P)
            def power_law_model(q, B, P):
                return B * q**(-P)

            # Prepare parameters for fitting
            from scipy.optimize import curve_fit

            # Initial guess
            p0 = [local_b, local_p]

            # Setup which parameters to fit
            if not params['fit_B'] and not params['fit_P']:
                # Both fixed - already handled above, shouldn't get here
                return
            elif not params['fit_B']:
                # Fix B, fit only P
                def model_fixed_b(q, P):
                    return power_law_model(q, local_b, P)
                p0_fit = [local_p]
                popt, pcov = curve_fit(model_fixed_b, q_fit, intensity_fit, p0=p0_fit)
                fitted_b = local_b
                fitted_p = abs(popt[0])
            elif not params['fit_P']:
                # Fix P, fit only B
                def model_fixed_p(q, B):
                    return power_law_model(q, B, local_p)
                p0_fit = [local_b]
                popt, pcov = curve_fit(model_fixed_p, q_fit, intensity_fit, p0=p0_fit)
                fitted_b = abs(popt[0])
                fitted_p = local_p
            else:
                # Fit both B and P
                popt, pcov = curve_fit(power_law_model, q_fit, intensity_fit, p0=p0)
                fitted_b = abs(popt[0])
                fitted_p = abs(popt[1])

            # Set P limits based on Igor code (lines 427-432)
            # P low limit = 1
            # P high limit = 3 (if mass fractal) or 4 (otherwise)
            # For now, we'll use 4 as the high limit (non-mass fractal)
            p_low = 1.0
            p_high = 4.0

            # Update GUI with fitted values
            level_widget.set_parameters({
                'B': fitted_b,
                'P': fitted_p,
                'B_low': fitted_b / 5,
                'B_high': fitted_b * 5,
                'P_low': p_low,
                'P_high': p_high
            })

            # Fix limits for this level
            level_widget.fix_limits()

            # Store local fit curve for plotting if display is enabled
            # Calculate Porod/power law curve over the Q range used for fitting
            porod_calc = fitted_b * q_fit**(-fitted_p)

            # Store in local_fits dictionary
            if level not in self.local_fits:
                self.local_fits[level] = {}
            self.local_fits[level]['porod'] = (q_fit, porod_calc)

            # Recalculate and update plot (this will also plot local fits if checkbox is enabled)
            if self.update_auto_check.isChecked():
                self.graph_unified()
            elif self.display_local_check.isChecked():
                # If not auto-updating, manually redraw the graph with local fits
                self.graph_unified()

            # Show success message
            self.status_label.setText(
                f"Local Porod fit complete for Level {level}: "
                f"B = {fitted_b:.3e}, P = {fitted_p:.3e} "
                f"(Q: {q_min:.3e} - {q_max:.3e}, {len(q_fit)} pts)"
            )

            self.graph_window.show_success_message(
                f"Local Porod fit completed for Level {level}! "
                f"B = {fitted_b:.4e}, P = {fitted_p:.4e} | "
                f"Q: {q_min:.3e} - {q_max:.3e} ({len(q_fit)} pts) | "
                f"Limits updated automatically."
            )

        except Exception as e:
            self.graph_window.show_error_message(f"Error during local Porod fit: {str(e)}")
            self.status_label.setText("Fit failed")
            import traceback
            traceback.print_exc()

    def plot_local_fits(self):
        """
        Plot all stored local fit curves (Guinier and Porod) on the graph.
        Uses green color and different line styles for Guinier (dashed) and Porod (dotted).
        """
        if not self.local_fits:
            return

        # Use green color for all local fits (distinct from blue data points and red unified fit line)
        local_fit_color = (0, 180, 0)  # Green color

        # Plot local fits for each level
        for level, fits in self.local_fits.items():
            if level < 1 or level > 5:
                continue

            # Plot Guinier fit (dashed line)
            if 'guinier' in fits:
                q_data, i_data = fits['guinier']
                self.graph_window.main_plot.plot(
                    q_data, i_data,
                    pen=pg.mkPen(color=local_fit_color, width=2, style=Qt.PenStyle.DashLine),
                    name=f'Level {level} Guinier fit'
                )

            # Plot Porod fit (dotted line)
            if 'porod' in fits:
                q_data, i_data = fits['porod']
                self.graph_window.main_plot.plot(
                    q_data, i_data,
                    pen=pg.mkPen(color=local_fit_color, width=2, style=Qt.PenStyle.DotLine),
                    name=f'Level {level} Porod fit'
                )

    def clear_local_fits(self):
        """Clear all stored local fit curves."""
        self.local_fits = {}

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
            'no_limits': self.no_limits_check.isChecked()
        }

        # Get all level parameters
        for i in range(5):
            params = self.level_widgets[i].get_parameters()
            level_state = {
                'level': i + 1,
                'G': {
                    'value': params['G'],
                    'fit': params['fit_G'],
                    'low_limit': params['G_low'],
                    'high_limit': params['G_high']
                },
                'Rg': {
                    'value': params['Rg'],
                    'fit': params['fit_Rg'],
                    'low_limit': params['Rg_low'],
                    'high_limit': params['Rg_high']
                },
                'B': {
                    'value': params['B'],
                    'fit': params['fit_B'],
                    'low_limit': params['B_low'],
                    'high_limit': params['B_high']
                },
                'P': {
                    'value': params['P'],
                    'fit': params['fit_P'],
                    'low_limit': params['P_low'],
                    'high_limit': params['P_high']
                },
                'ETA': {
                    'value': params['ETA'],
                    'fit': params['fit_ETA'],
                    'low_limit': params['ETA_low'],
                    'high_limit': params['ETA_high']
                },
                'PACK': {
                    'value': params['PACK'],
                    'fit': params['fit_PACK'],
                    'low_limit': params['PACK_low'],
                    'high_limit': params['PACK_high']
                },
                'RgCutoff': params['RgCutoff'],
                'correlated': params['correlated'],
                'estimate_B': params['estimate_B'],
                'link_rgco': params['link_rgco']
            }
            state['levels'].append(level_state)

        return state

    def apply_state(self, state: Dict):
        """Apply saved state to GUI."""
        # Set number of levels
        self.num_levels_spin.setValue(state.get('num_levels', 1))

        # Set background
        bg = state.get('background', {})
        self.background_value.setText(self.format_value_3sig(bg.get('value', 1e-6)))
        self.fit_background_check.setChecked(bg.get('fit', False))

        # Set checkboxes
        self.update_auto_check.setChecked(state.get('update_auto', False))
        self.display_local_check.setChecked(state.get('display_local', False))
        self.no_limits_check.setChecked(state.get('no_limits', False))

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

                # Set parameter values and limits
                level_widget.g_value.setText(str(level_state['G']['value']))
                level_widget.g_fit.setChecked(level_state['G']['fit'])
                if level_state['G'].get('low_limit') is not None:
                    level_widget.g_low.setText(str(level_state['G']['low_limit']))
                if level_state['G'].get('high_limit') is not None:
                    level_widget.g_high.setText(str(level_state['G']['high_limit']))

                level_widget.rg_value.setText(str(level_state['Rg']['value']))
                level_widget.rg_fit.setChecked(level_state['Rg']['fit'])
                if level_state['Rg'].get('low_limit') is not None:
                    level_widget.rg_low.setText(str(level_state['Rg']['low_limit']))
                if level_state['Rg'].get('high_limit') is not None:
                    level_widget.rg_high.setText(str(level_state['Rg']['high_limit']))

                level_widget.b_value.setText(str(level_state['B']['value']))
                level_widget.b_fit.setChecked(level_state['B']['fit'])
                if level_state['B'].get('low_limit') is not None:
                    level_widget.b_low.setText(str(level_state['B']['low_limit']))
                if level_state['B'].get('high_limit') is not None:
                    level_widget.b_high.setText(str(level_state['B']['high_limit']))

                level_widget.p_value.setText(str(level_state['P']['value']))
                level_widget.p_fit.setChecked(level_state['P']['fit'])
                if level_state['P'].get('low_limit') is not None:
                    level_widget.p_low.setText(str(level_state['P']['low_limit']))
                if level_state['P'].get('high_limit') is not None:
                    level_widget.p_high.setText(str(level_state['P']['high_limit']))

                # ETA and PACK parameters
                if 'ETA' in level_state:
                    level_widget.eta_value.setText(str(level_state['ETA']['value']))
                    level_widget.eta_fit.setChecked(level_state['ETA']['fit'])
                    if level_state['ETA'].get('low_limit') is not None:
                        level_widget.eta_low.setText(str(level_state['ETA']['low_limit']))
                    if level_state['ETA'].get('high_limit') is not None:
                        level_widget.eta_high.setText(str(level_state['ETA']['high_limit']))

                if 'PACK' in level_state:
                    level_widget.pack_value.setText(str(level_state['PACK']['value']))
                    level_widget.pack_fit.setChecked(level_state['PACK']['fit'])
                    if level_state['PACK'].get('low_limit') is not None:
                        level_widget.pack_low.setText(str(level_state['PACK']['low_limit']))
                    if level_state['PACK'].get('high_limit') is not None:
                        level_widget.pack_high.setText(str(level_state['PACK']['high_limit']))

                level_widget.rg_cutoff.setText(str(level_state.get('RgCutoff', 0)))
                level_widget.correlated_check.setChecked(level_state.get('correlated', False))
                level_widget.estimate_b_check.setChecked(level_state.get('estimate_B', False))
                if level_widget.link_rgco_check:
                    level_widget.link_rgco_check.setChecked(level_state.get('link_rgco', False))

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

    def fix_all_limits(self):
        """Fix fitting limits for all levels based on current parameter values."""
        num_levels = self.num_levels_spin.value()
        for i in range(num_levels):
            self.level_widgets[i].fix_limits()
        self.status_label.setText(f"Fixed limits for {num_levels} level(s)")

    def sync_rgcutoff_links(self):
        """
        Sync RgCutoff values based on Link RgCutoff checkboxes.
        For each level with Link RgCutoff checked, set RgCutoff to Rg of previous level.
        """
        num_levels = self.num_levels_spin.value()
        for i in range(1, num_levels):  # Start from level 2 (index 1)
            level_widget = self.level_widgets[i]
            # Check if this level has the link checkbox and if it's checked
            if level_widget.link_rgco_check and level_widget.link_rgco_check.isChecked():
                # Get Rg from previous level (i-1)
                prev_rg = float(self.level_widgets[i-1].rg_value.text() or 0)
                # Set current level's RgCutoff to previous level's Rg
                level_widget.rg_cutoff.setText(level_widget._format_value(prev_rg))

    def backup_parameters(self):
        """Backup current parameters before fitting."""
        self.parameter_backup = {
            'num_levels': self.num_levels_spin.value(),
            'background': float(self.background_value.text() or 0),
            'levels': []
        }

        # Backup all level parameters
        for i in range(5):
            params = self.level_widgets[i].get_parameters()
            self.parameter_backup['levels'].append(params.copy())

        self.status_label.setText("Parameters backed up before fitting")

    def revert_to_backup(self):
        """Restore parameters from backup and recalculate the model."""
        if self.parameter_backup is None:
            QMessageBox.warning(
                self,
                "No Backup",
                "No parameter backup available. Run a fit first to create a backup."
            )
            return

        if self.data is None:
            QMessageBox.warning(self, "No Data", "Please load data first.")
            return

        try:
            # Restore number of levels
            self.num_levels_spin.setValue(self.parameter_backup['num_levels'])

            # Restore background
            self.background_value.setText(self.format_value_3sig(self.parameter_backup['background']))

            # Restore all level parameters
            for i in range(5):
                if i < len(self.parameter_backup['levels']):
                    self.level_widgets[i].set_parameters(self.parameter_backup['levels'][i])

            # Recalculate the model with restored parameters
            num_levels = self.parameter_backup['num_levels']
            levels = []

            for i in range(num_levels):
                params = self.parameter_backup['levels'][i]
                level = UnifiedLevel(
                    Rg=params['Rg'],
                    G=params['G'],
                    P=params['P'],
                    B=params['B'],
                    RgCO=params['RgCutoff'],
                    ETA=params['ETA'],
                    PACK=params['PACK'],
                    correlations=params['correlated'],
                    Rg_limits=(params['Rg_low'], params['Rg_high']),
                    G_limits=(params['G_low'], params['G_high']),
                    P_limits=(params['P_low'], params['P_high']),
                    B_limits=(params['B_low'], params['B_high']),
                    ETA_limits=(params['ETA_low'], params['ETA_high']),
                    PACK_limits=(params['PACK_low'], params['PACK_high'])
                )
                levels.append(level)

            # Update model
            self.model.num_levels = num_levels
            self.model.levels = levels
            self.model.background = self.parameter_backup['background']

            # Recalculate
            intensity_calc = self.model.calculate_intensity(self.data['Q'])

            # Re-plot
            self.graph_window.init_plots()
            self.graph_window.plot_data(
                self.data['Q'],
                self.data['Intensity'],
                self.data.get('Error'),
                self.data['label']
            )
            self.graph_window.plot_fit(self.data['Q'], intensity_calc, 'Restored Model')

            # Residuals
            if self.data.get('Error') is not None:
                residuals = (self.data['Intensity'] - intensity_calc) / self.data['Error']
            else:
                residuals = (self.data['Intensity'] - intensity_calc) / self.data['Intensity']

            self.graph_window.plot_residuals(self.data['Q'], residuals)
            self.status_label.setText("Reverted to pre-fit parameters")

        except Exception as e:
            QMessageBox.critical(self, "Revert Error", f"Error reverting to backup:\n{str(e)}")
            import traceback
            traceback.print_exc()

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
