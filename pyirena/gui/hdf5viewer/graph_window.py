"""
GraphWindow — floating plot window with pyqtgraph.

Each GraphWindow is an independent QWidget (not QMainWindow, to avoid
extra taskbar entries).  It contains a compact toolbar and a pyqtgraph
PlotItem.  Multiple GraphWindows can coexist; the most recently focused
one is the "active" graph that receives new curves from PlotControlsPanel.

Key features
------------
- Toolbar: X/Y log toggles, Labels dialog, Legend toggle, Save JPEG, Save CSV
- Right-click (ViewBox menu): set range, curve styles, Save PNG, HDF5, ITX,
  matplotlib figure
- Curve registry: _curves list tracks all plotted data for export
- Active-graph tracking via `became_active` signal
"""

from __future__ import annotations

import itertools
from typing import Any

import numpy as np
import pyqtgraph as pg

try:
    from PySide6.QtWidgets import (
        QWidget, QVBoxLayout, QToolBar, QDialog, QFormLayout, QLineEdit,
        QDialogButtonBox, QLabel, QPushButton, QColorDialog, QComboBox,
        QDoubleSpinBox, QHBoxLayout, QScrollArea, QGridLayout, QGroupBox,
        QMessageBox, QAction, QSizePolicy,
    )
    from PySide6.QtCore import Qt, Signal, QEvent
    from PySide6.QtGui import QColor, QIcon
except ImportError:
    from PyQt6.QtWidgets import (  # type: ignore[no-redef]
        QWidget, QVBoxLayout, QToolBar, QDialog, QFormLayout, QLineEdit,
        QDialogButtonBox, QLabel, QPushButton, QColorDialog, QComboBox,
        QDoubleSpinBox, QHBoxLayout, QScrollArea, QGridLayout, QGroupBox,
        QMessageBox, QAction, QSizePolicy,
    )
    from PyQt6.QtCore import Qt, pyqtSignal as Signal, QEvent  # type: ignore[no-redef]
    from PyQt6.QtGui import QColor, QIcon                       # type: ignore[no-redef]

from . import export as _export

# ── Default color cycle (distinct, colorblind-friendly palette) ────────────
_COLOR_CYCLE = [
    "#2980b9",  # blue
    "#c0392b",  # red
    "#27ae60",  # green
    "#8e44ad",  # purple
    "#e67e22",  # orange
    "#16a085",  # teal
    "#2c3e50",  # dark slate
    "#f39c12",  # amber
    "#1abc9c",  # emerald
    "#e74c3c",  # tomato
]

_color_iter = itertools.cycle(_COLOR_CYCLE)

# Class-level counter for window titles
_graph_count = 0


def _next_color() -> str:
    return next(_color_iter)


class GraphWindow(QWidget):
    """
    A floating plot window with pyqtgraph.

    Signals
    -------
    became_active()  — emitted when this window gains focus; used by the main
                       window to track which graph is currently "active".
    closed()         — emitted on close; for cleanup in the main window.
    """

    became_active = Signal()
    closed        = Signal()

    def __init__(self, parent: QWidget | None = None, title: str | None = None) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        global _graph_count
        _graph_count += 1
        self._graph_number = _graph_count
        self._title = title or f"Graph {_graph_count}"

        self._curves: list[dict[str, Any]] = []
        self._log_x = False
        self._log_y = False
        self._x_label = "X"
        self._y_label = "Y"
        self._legend: pg.LegendItem | None = None

        self.setWindowTitle(self._title)
        self.resize(750, 520)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)

        self._build_ui()

    # ── UI construction ────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(0)

        # Toolbar
        tb = QToolBar()
        tb.setMovable(False)

        self._log_x_btn = QPushButton("X: Lin")
        self._log_x_btn.setCheckable(True)
        self._log_x_btn.setFixedWidth(55)
        self._log_x_btn.setToolTip("Toggle X axis between linear and logarithmic")
        self._log_x_btn.clicked.connect(self._toggle_log_x)
        tb.addWidget(self._log_x_btn)

        self._log_y_btn = QPushButton("Y: Lin")
        self._log_y_btn.setCheckable(True)
        self._log_y_btn.setFixedWidth(55)
        self._log_y_btn.setToolTip("Toggle Y axis between linear and logarithmic")
        self._log_y_btn.clicked.connect(self._toggle_log_y)
        tb.addWidget(self._log_y_btn)

        tb.addSeparator()

        labels_btn = QPushButton("Labels…")
        labels_btn.setFixedWidth(70)
        labels_btn.setToolTip("Edit title and axis labels")
        labels_btn.clicked.connect(self._edit_labels)
        tb.addWidget(labels_btn)

        self._legend_btn = QPushButton("Legend")
        self._legend_btn.setCheckable(True)
        self._legend_btn.setChecked(True)
        self._legend_btn.setFixedWidth(60)
        self._legend_btn.setToolTip("Show/hide legend")
        self._legend_btn.clicked.connect(self._toggle_legend)
        tb.addWidget(self._legend_btn)

        tb.addSeparator()

        jpeg_btn = QPushButton("Save JPEG")
        jpeg_btn.setFixedWidth(80)
        jpeg_btn.clicked.connect(lambda: _export.save_jpeg(self))
        tb.addWidget(jpeg_btn)

        csv_btn = QPushButton("Save CSV")
        csv_btn.setFixedWidth(75)
        csv_btn.clicked.connect(lambda: _export.save_csv(self))
        tb.addWidget(csv_btn)

        layout.addWidget(tb)

        # pyqtgraph plot
        self._glw = pg.GraphicsLayoutWidget()
        self._plot: pg.PlotItem = self._glw.addPlot()
        self._plot.showGrid(x=True, y=True, alpha=0.3)
        self._plot.setLabel("left", self._y_label)
        self._plot.setLabel("bottom", self._x_label)

        # Add right-click menu items
        self._add_viewbox_menu()

        layout.addWidget(self._glw, 1)

        # Status label
        self._status = QLabel("")
        self._status.setStyleSheet("font-size:9pt; color:#555; padding:1px 4px;")
        layout.addWidget(self._status)

    def _add_viewbox_menu(self) -> None:
        vb = self._plot.getViewBox()

        act_xrange = QAction("Set X range…", self)
        act_xrange.triggered.connect(self._set_x_range)
        vb.menu.addAction(act_xrange)

        act_yrange = QAction("Set Y range…", self)
        act_yrange.triggered.connect(self._set_y_range)
        vb.menu.addAction(act_yrange)

        vb.menu.addSeparator()

        act_styles = QAction("Curve styles…", self)
        act_styles.triggered.connect(self._edit_curve_styles)
        vb.menu.addAction(act_styles)

        vb.menu.addSeparator()

        act_png = QAction("Save PNG…", self)
        act_png.triggered.connect(lambda: _export.save_png(self))
        vb.menu.addAction(act_png)

        act_hdf5 = QAction("Save HDF5…", self)
        act_hdf5.triggered.connect(lambda: _export.save_hdf5(self))
        vb.menu.addAction(act_hdf5)

        act_itx = QAction("Save ITX (Igor Pro)…", self)
        act_itx.triggered.connect(lambda: _export.save_itx(self))
        vb.menu.addAction(act_itx)

        act_mpl = QAction("Open as matplotlib figure…", self)
        act_mpl.triggered.connect(lambda: _export.open_matplotlib(self))
        vb.menu.addAction(act_mpl)

    # ── Public API ─────────────────────────────────────────────────────────

    def add_curve(
        self,
        x: np.ndarray,
        y: np.ndarray,
        label: str = "",
        yerr: np.ndarray | None = None,
        xerr: np.ndarray | None = None,
        color: str | None = None,
        suggest_log_x: bool = False,
        suggest_log_y: bool = False,
    ) -> None:
        """
        Add a curve to the graph.

        Parameters
        ----------
        x, y        : data arrays (physical values — pyqtgraph handles log scaling)
        label       : legend label
        yerr        : optional Y error bars
        xerr        : optional X error bars (not yet displayed but stored for export)
        color       : hex color string; auto-assigned from color cycle if None
        suggest_log_x/y : if the axes are currently linear, switch to log
        """
        if color is None:
            color = _next_color()

        style = {"color": color, "width": 1.5, "symbol": None}
        self._curves.append({
            "label": label,
            "x": np.asarray(x, float),
            "y": np.asarray(y, float),
            "yerr": np.asarray(yerr, float) if yerr is not None else None,
            "xerr": np.asarray(xerr, float) if xerr is not None else None,
            "style": style,
            "plot_ref": None,
            "err_ref":  None,
        })

        # Optionally switch to log axes
        if suggest_log_x and not self._log_x:
            self._set_log_x(True)
        if suggest_log_y and not self._log_y:
            self._set_log_y(True)

        self._redraw_curve(self._curves[-1])
        self._update_legend()
        self._status.setText(f"{len(self._curves)} curve(s)")

    def clear_curves(self) -> None:
        self._plot.clear()
        self._curves.clear()
        self._legend = None
        self._status.setText("(cleared)")

    def get_curves(self) -> list[dict]:
        return list(self._curves)

    # Axis mode getters (used by export.open_matplotlib)
    def is_log_x(self) -> bool: return self._log_x
    def is_log_y(self) -> bool: return self._log_y
    def get_title(self) -> str: return self._title
    def get_x_label(self) -> str: return self._x_label
    def get_y_label(self) -> str: return self._y_label

    # ── Drawing ────────────────────────────────────────────────────────────

    def _redraw_curve(self, curve: dict) -> None:
        """Draw a single curve and its error bars onto the PlotItem."""
        x = curve["x"]
        y = curve["y"]
        style = curve["style"]
        color = style.get("color", "#2980b9")
        width = style.get("width", 1.5)
        symbol = style.get("symbol", None)

        pen = pg.mkPen(color=color, width=width)
        if symbol:
            brush = pg.mkBrush(color)
            ref = self._plot.plot(
                x, y, pen=pen, symbol=symbol,
                symbolBrush=brush, symbolSize=6,
                name=curve["label"],
            )
        else:
            ref = self._plot.plot(x, y, pen=pen, name=curve["label"])

        curve["plot_ref"] = ref

        # Y error bars as NaN-separated lines
        if curve.get("yerr") is not None:
            ex, ey = self._build_error_bars(x, y, curve["yerr"])
            err_ref = self._plot.plot(
                ex, ey,
                pen=pg.mkPen(color=color, width=max(1, width - 0.5)),
                connect="finite",
            )
            curve["err_ref"] = err_ref

    def _redraw_all(self) -> None:
        """Clear the PlotItem and redraw all curves from scratch."""
        # Save and restore legend state
        had_legend = self._legend_btn.isChecked()
        self._plot.clear()
        self._legend = None
        for curve in self._curves:
            curve["plot_ref"] = None
            curve["err_ref"]  = None
            self._redraw_curve(curve)
        if had_legend:
            self._update_legend()

    @staticmethod
    def _build_error_bars(x, y, yerr, cap_frac=0.02):
        """Build NaN-separated (x,y) arrays for vertical error bars + caps."""
        xlines, ylines = [], []
        for i in range(len(x)):
            if not (np.isfinite(y[i]) and np.isfinite(yerr[i])):
                continue
            ytop = y[i] + yerr[i]
            ybot = y[i] - yerr[i]
            xlines += [x[i], x[i], np.nan]
            ylines += [ybot, ytop, np.nan]
            # Horizontal caps
            if x[i] > 0:
                xl = x[i] / (1 + cap_frac)
                xr = x[i] * (1 + cap_frac)
                xlines += [xl, xr, np.nan, xl, xr, np.nan]
                ylines += [ytop, ytop, np.nan, ybot, ybot, np.nan]
        return list(xlines), list(ylines)

    def _update_legend(self) -> None:
        if not self._legend_btn.isChecked():
            return
        if self._legend is None:
            self._legend = self._plot.addLegend(offset=(10, 10))
        # The legend auto-updates from named plot items

    # ── Toolbar actions ────────────────────────────────────────────────────

    def _toggle_log_x(self, checked: bool) -> None:
        self._set_log_x(checked)

    def _toggle_log_y(self, checked: bool) -> None:
        self._set_log_y(checked)

    def _set_log_x(self, on: bool) -> None:
        self._log_x = on
        self._log_x_btn.setChecked(on)
        self._log_x_btn.setText("X: Log" if on else "X: Lin")
        self._plot.setLogMode(x=on, y=self._log_y)

    def _set_log_y(self, on: bool) -> None:
        self._log_y = on
        self._log_y_btn.setChecked(on)
        self._log_y_btn.setText("Y: Log" if on else "Y: Lin")
        self._plot.setLogMode(x=self._log_x, y=on)

    def _toggle_legend(self, checked: bool) -> None:
        if checked:
            if self._legend is None:
                self._update_legend()
            elif self._legend.scene():
                self._legend.setVisible(True)
        else:
            if self._legend is not None and self._legend.scene():
                self._legend.setVisible(False)

    def _edit_labels(self) -> None:
        dlg = _LabelsDialog(self._title, self._x_label, self._y_label, parent=self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self._title, self._x_label, self._y_label = dlg.get_values()
            self.setWindowTitle(self._title)
            self._plot.setTitle(self._title)
            self._plot.setLabel("bottom", self._x_label)
            self._plot.setLabel("left",   self._y_label)

    def _set_x_range(self) -> None:
        vr = self._plot.getViewBox().viewRange()
        dlg = _RangeDialog("X range", vr[0][0], vr[0][1], parent=self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            lo, hi = dlg.get_values()
            self._plot.setXRange(lo, hi, padding=0)

    def _set_y_range(self) -> None:
        vr = self._plot.getViewBox().viewRange()
        dlg = _RangeDialog("Y range", vr[1][0], vr[1][1], parent=self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            lo, hi = dlg.get_values()
            self._plot.setYRange(lo, hi, padding=0)

    def _edit_curve_styles(self) -> None:
        if not self._curves:
            QMessageBox.information(self, "No curves", "No curves to style.")
            return
        dlg = _CurveStylesDialog(self._curves, parent=self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self._redraw_all()

    # ── Focus / lifecycle ──────────────────────────────────────────────────

    def changeEvent(self, event: QEvent) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.ActivationChange:
            if self.isActiveWindow():
                self.became_active.emit()

    def closeEvent(self, event) -> None:
        self.closed.emit()
        super().closeEvent(event)


# ── Dialogs ────────────────────────────────────────────────────────────────

class _LabelsDialog(QDialog):
    def __init__(self, title, x_label, y_label, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Edit Labels")
        form = QFormLayout(self)
        self._title_edit = QLineEdit(title)
        self._x_edit     = QLineEdit(x_label)
        self._y_edit     = QLineEdit(y_label)
        form.addRow("Title:", self._title_edit)
        form.addRow("X label:", self._x_edit)
        form.addRow("Y label:", self._y_edit)
        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        form.addRow(btns)

    def get_values(self):
        return (
            self._title_edit.text(),
            self._x_edit.text(),
            self._y_edit.text(),
        )


class _RangeDialog(QDialog):
    def __init__(self, axis_name: str, lo: float, hi: float, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Set {axis_name}")
        form = QFormLayout(self)
        self._lo = QDoubleSpinBox()
        self._hi = QDoubleSpinBox()
        for sp in (self._lo, self._hi):
            sp.setRange(-1e15, 1e15)
            sp.setDecimals(6)
            sp.setSingleStep(0.1)
        self._lo.setValue(lo)
        self._hi.setValue(hi)
        form.addRow("Min:", self._lo)
        form.addRow("Max:", self._hi)
        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        form.addRow(btns)

    def get_values(self):
        return self._lo.value(), self._hi.value()


class _CurveStylesDialog(QDialog):
    """Dialog to edit color, line width, and symbol for each curve."""

    _SYMBOLS = ["none", "o", "s", "t", "d", "+", "x"]

    def __init__(self, curves: list[dict], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Curve Styles")
        self._curves = curves
        self._color_btns = []
        self._width_spins = []
        self._symbol_combos = []

        main = QVBoxLayout(self)
        grid = QGridLayout()
        grid.addWidget(QLabel("Curve"), 0, 0)
        grid.addWidget(QLabel("Color"), 0, 1)
        grid.addWidget(QLabel("Width"), 0, 2)
        grid.addWidget(QLabel("Symbol"), 0, 3)

        for i, curve in enumerate(curves):
            style = curve.get("style", {})
            color = style.get("color", "#2980b9")
            width = style.get("width", 1.5)
            symbol = style.get("symbol") or "none"

            lbl = QLabel(curve["label"][:40])
            grid.addWidget(lbl, i + 1, 0)

            btn = QPushButton()
            btn.setFixedSize(28, 20)
            btn.setStyleSheet(f"background-color:{color};")
            btn.setProperty("curve_idx", i)
            btn.setProperty("color", color)
            btn.clicked.connect(self._pick_color)
            grid.addWidget(btn, i + 1, 1)
            self._color_btns.append(btn)

            sp = QDoubleSpinBox()
            sp.setRange(0.5, 8.0)
            sp.setSingleStep(0.5)
            sp.setDecimals(1)
            sp.setValue(width)
            sp.setFixedWidth(60)
            grid.addWidget(sp, i + 1, 2)
            self._width_spins.append(sp)

            combo = QComboBox()
            combo.addItems(self._SYMBOLS)
            idx = self._SYMBOLS.index(symbol) if symbol in self._SYMBOLS else 0
            combo.setCurrentIndex(idx)
            combo.setFixedWidth(70)
            grid.addWidget(combo, i + 1, 3)
            self._symbol_combos.append(combo)

        main.addLayout(grid)
        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        btns.accepted.connect(self._apply)
        btns.rejected.connect(self.reject)
        main.addWidget(btns)

    def _pick_color(self) -> None:
        btn = self.sender()
        current = QColor(btn.property("color"))
        color = QColorDialog.getColor(current, self)
        if color.isValid():
            hex_color = color.name()
            btn.setProperty("color", hex_color)
            btn.setStyleSheet(f"background-color:{hex_color};")

    def _apply(self) -> None:
        for i, curve in enumerate(self._curves):
            curve["style"]["color"]  = self._color_btns[i].property("color")
            curve["style"]["width"]  = self._width_spins[i].value()
            sym = self._symbol_combos[i].currentText()
            curve["style"]["symbol"] = None if sym == "none" else sym
        self.accept()
