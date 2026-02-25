"""
WAXS Peak-Fitting GUI panel for pyIrena.

Classes
-------
WAXSPeakFitGraphWindow  — standalone graph window (2 linear/linear panels)
WAXSPeakFitPanel        — full GUI: left controls + right graph window
PeakRowWidget           — collapsible widget for one peak's parameters
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pyqtgraph as pg

try:
    from PySide6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QSplitter,
        QLabel, QComboBox, QCheckBox, QPushButton, QLineEdit,
        QDoubleSpinBox, QScrollArea, QGroupBox, QFileDialog,
        QMessageBox, QFrame, QSizePolicy, QSpacerItem,
        QDialog, QDialogButtonBox,
    )
    from PySide6.QtCore import Qt, Signal, QTimer
    from PySide6.QtGui import QFont, QDoubleValidator
except ImportError:
    from PyQt6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QSplitter,
        QLabel, QComboBox, QCheckBox, QPushButton, QLineEdit,
        QDoubleSpinBox, QScrollArea, QGroupBox, QFileDialog,
        QMessageBox, QFrame, QSizePolicy, QSpacerItem,
        QDialog, QDialogButtonBox,
    )
    from PyQt6.QtCore import Qt, pyqtSignal as Signal, QTimer
    from PyQt6.QtGui import QFont, QDoubleValidator

from pyirena.core.waxs_peakfit import (
    PEAK_SHAPES, BG_SHAPES,
    _PEAK_PARAM_NAMES, _BG_NCOEFFS, _PARAM_DEFAULTS,
    bg_param_names, default_bg_params, default_peak_params,
    eval_background, eval_peak, eval_model,
    find_peaks_in_data, WAXSPeakFitModel,
)


# ── colour palette for peaks ──────────────────────────────────────────────
_PEAK_COLORS = [
    '#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6',
    '#1abc9c', '#e67e22', '#34495e', '#d35400', '#27ae60',
    '#2980b9', '#8e44ad', '#16a085', '#c0392b', '#7f8c8d',
]


def _peak_color(i: int) -> str:
    return _PEAK_COLORS[i % len(_PEAK_COLORS)]


# ── shared field styles ───────────────────────────────────────────────────
_FIELD_STYLE = "QLineEdit { max-width: 90px; }"
_SPIN_STYLE  = "QDoubleSpinBox { max-width: 90px; }"

_BTN_GRAPH   = "QPushButton { background:#52c77a;color:white;font-weight:bold; } QPushButton:disabled{background:#bdc3c7;}"
_BTN_FIT     = "QPushButton { background:#27ae60;color:white;font-weight:bold; } QPushButton:disabled{background:#bdc3c7;}"
_BTN_REVERT  = "QPushButton { background:#e67e22;color:white;font-weight:bold; } QPushButton:disabled{background:#bdc3c7;}"
_BTN_FIND    = "QPushButton { background:#3498db;color:white;font-weight:bold; } QPushButton:disabled{background:#bdc3c7;}"
_BTN_RESET   = "QPushButton { background:#e67e22;color:white; } QPushButton:disabled{background:#bdc3c7;}"
_BTN_SAVE    = "QPushButton { background:#3498db;color:white; } QPushButton:disabled{background:#bdc3c7;}"
_BTN_STORE   = "QPushButton { background:#52c77a;color:white; } QPushButton:disabled{background:#bdc3c7;}"
_BTN_RGRAPH  = "QPushButton { background:#81c784;color:white; } QPushButton:disabled{background:#bdc3c7;}"
_BTN_EXPORT  = "QPushButton { background:#a9d18e;color:black; } QPushButton:disabled{background:#bdc3c7;}"
_BTN_ADD     = "QPushButton { background:#95a5a6;color:white; } QPushButton:disabled{background:#bdc3c7;}"
_BTN_REMOVE  = "QPushButton { background:#e74c3c;color:white; } QPushButton:disabled{background:#bdc3c7;}"
_BTN_FIXLIM  = "QPushButton { background:#8e44ad;color:white; } QPushButton:disabled{background:#bdc3c7;}"

# Sensible hi-limit defaults for peak parameters (when user hasn't set them)
_PEAK_PARAM_HI_DEFAULTS = {"A": 1e9, "Q0": 100.0, "FWHM": 10.0, "eta": 1.0}
_PEAK_PARAM_LO_DEFAULTS = {"A": 0.0,  "Q0": 0.0,   "FWHM": 1e-6, "eta": 0.0}

# Wheel-spin min/max for peak param value fields (prevents negative)
_PEAK_PARAM_VAL_MIN = {"A": 0.0, "Q0": 0.0, "FWHM": 0.0, "eta": 0.0}
_PEAK_PARAM_VAL_MAX = {"A": np.inf, "Q0": np.inf, "FWHM": np.inf, "eta": 1.0}


def _connect_limit_check(val_fld, lo_fld, hi_fld):
    """Connect val_fld so that typing outside [lo, hi] auto-expands the limit."""
    def _check():
        v  = val_fld.float_value()
        lo = lo_fld.float_value()
        hi = hi_fld.float_value()
        if v < lo:
            lo_fld.set_float(v)
        elif v > hi:
            hi_fld.set_float(v)
    val_fld.editingFinished.connect(_check)


def _label(text: str, bold: bool = False, size: int = 0) -> QLabel:
    lbl = QLabel(text)
    if bold or size:
        f = lbl.font()
        if bold:
            f.setBold(True)
        if size:
            f.setPointSize(size)
        lbl.setFont(f)
    return lbl


def _hline() -> QFrame:
    line = QFrame()
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Sunken)
    return line


class _ValidatedField(QLineEdit):
    """QLineEdit that accepts float input with optional bounds and mouse-wheel support."""

    def __init__(self, value: float = 0.0,
                 min_val: float = -np.inf, max_val: float = np.inf,
                 wheel_enabled: bool = True,
                 wheel_frac: float = 0.01,
                 parent=None):
        super().__init__(parent)
        self._min_val = min_val
        self._max_val = max_val
        self._wheel_enabled = wheel_enabled
        self._wheel_frac = wheel_frac
        self.setValidator(QDoubleValidator())
        self.setText(f"{value:.6g}")
        self.setFixedWidth(90)

    def float_value(self) -> float:
        try:
            return float(self.text())
        except ValueError:
            return 0.0

    def set_float(self, v: float):
        self.setText(f"{v:g}")

    def wheelEvent(self, ev):
        """Shift = 10× step; normal = 1% of |value| (min 1e-6 absolute)."""
        if not self._wheel_enabled:
            ev.ignore()
            return
        val  = self.float_value()
        mods = ev.modifiers()
        frac = self._wheel_frac * 10.0 if (mods & Qt.KeyboardModifier.ShiftModifier) else self._wheel_frac
        step = max(abs(val) * frac, 1e-6)
        delta = ev.angleDelta().y()
        if delta > 0:
            val += step
        elif delta < 0:
            val -= step
        val = float(np.clip(val, self._min_val, self._max_val))
        self.set_float(val)
        self.editingFinished.emit()
        ev.accept()


# ===========================================================================
# PeakRowWidget — one peak's UI row
# ===========================================================================

class PeakRowWidget(QWidget):
    """Widget displaying controls for a single peak (shape + parameters).

    Emits ``changed`` whenever a value is edited.
    Emits ``remove_requested`` with its own index when the Remove button is clicked.
    """
    changed         = Signal()
    remove_requested = Signal(object)  # self

    # Column indices in the grid
    _COL_NAME  = 0
    _COL_VAL   = 1
    _COL_FIT   = 2
    _COL_LO    = 3
    _COL_HI    = 4

    def __init__(self, index: int, peak_dict: Dict, parent=None):
        super().__init__(parent)
        self._index = index
        self._limit_widgets: List[QWidget] = []  # Lo/Hi fields + header labels
        self._param_rows: Dict[str, dict] = {}   # pname → {val, fit, lo, hi widgets}
        self._setup_ui(peak_dict)

    # ── UI construction ──────────────────────────────────────────────────

    def _setup_ui(self, peak_dict: Dict):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(2, 2, 2, 2)
        outer.setSpacing(2)

        # ── Header row: "Peak N" + shape combo + Remove button ────────────
        hdr = QHBoxLayout()
        hdr.setSpacing(4)
        color = _peak_color(self._index)
        lbl = QLabel(f"Peak {self._index + 1}")
        lbl.setStyleSheet(f"color:{color}; font-weight:bold;")
        hdr.addWidget(lbl)

        self._shape_combo = QComboBox()
        self._shape_combo.addItems(PEAK_SHAPES)
        shape = peak_dict.get("shape", "Gauss")
        self._shape_combo.setCurrentText(shape)
        self._shape_combo.currentTextChanged.connect(self._on_shape_changed)
        hdr.addWidget(self._shape_combo)

        hdr.addStretch()
        remove_btn = QPushButton("Remove")
        remove_btn.setStyleSheet(_BTN_REMOVE)
        remove_btn.setFixedHeight(20)
        remove_btn.setMaximumWidth(65)
        remove_btn.clicked.connect(lambda: self.remove_requested.emit(self))
        hdr.addWidget(remove_btn)
        outer.addLayout(hdr)

        # ── Parameter grid ────────────────────────────────────────────────
        self._grid = QGridLayout()
        self._grid.setHorizontalSpacing(4)
        self._grid.setVerticalSpacing(2)
        self._grid.setColumnMinimumWidth(self._COL_LO, 81)
        self._grid.setColumnMinimumWidth(self._COL_HI, 81)

        # Column headers (order: Name | Value | Fit? | Lo limit | Hi limit)
        for col, txt in [(self._COL_NAME, "Param"),
                         (self._COL_VAL,  "Value"),
                         (self._COL_FIT,  "Fit?"),
                         (self._COL_LO,   "Lo limit"),
                         (self._COL_HI,   "Hi limit")]:
            h = _label(txt, bold=True)
            self._grid.addWidget(h, 0, col)
            if col in (self._COL_LO, self._COL_HI):
                self._limit_widgets.append(h)

        outer.addLayout(self._grid)
        outer.addWidget(_hline())

        # Build rows for current shape
        self._build_param_rows(peak_dict)

    def _build_param_rows(self, peak_dict: Dict):
        """Create parameter rows from scratch for the current shape."""
        shape  = self._shape_combo.currentText()
        pnames = _PEAK_PARAM_NAMES.get(shape, [])

        # Clear old rows (keep header at row 0)
        for row_widgets in self._param_rows.values():
            for w in row_widgets.values():
                if hasattr(w, "setParent"):
                    w.setParent(None)
        self._limit_widgets = [w for w in self._limit_widgets if w.text() in ("Lo limit", "Hi limit", "Param", "Value", "Fit?")]
        self._param_rows = {}

        for row_i, pname in enumerate(pnames, start=1):
            pd = peak_dict.get(pname, {})
            if isinstance(pd, dict):
                val  = float(pd.get("value", _PARAM_DEFAULTS.get(pname, {}).get("value", 0.0)))
                fit  = bool(pd.get("fit", True))
                lo   = pd.get("lo")
                hi   = pd.get("hi")
            else:
                val, fit, lo, hi = float(pd), True, None, None

            name_lbl = _label(pname)
            # Q0 gets a 10× finer wheel step (positions need fine control)
            _wfrac = 0.001 if pname == "Q0" else 0.01
            val_fld  = _ValidatedField(
                val,
                min_val=_PEAK_PARAM_VAL_MIN.get(pname, -np.inf),
                max_val=_PEAK_PARAM_VAL_MAX.get(pname, np.inf),
                wheel_frac=_wfrac,
            )
            val_fld.setFixedWidth(104)
            fit_chk  = QCheckBox()
            fit_chk.setChecked(fit)
            lo_fld   = _ValidatedField(lo if lo is not None else _PEAK_PARAM_LO_DEFAULTS.get(pname, -1e6),
                                       wheel_enabled=False)
            lo_fld.setFixedWidth(81)
            hi_fld   = _ValidatedField(hi if hi is not None else _PEAK_PARAM_HI_DEFAULTS.get(pname, 1e6),
                                       wheel_enabled=False)
            hi_fld.setFixedWidth(81)

            # Auto-expand limits when value is typed outside range
            _connect_limit_check(val_fld, lo_fld, hi_fld)

            # Connect change signals
            val_fld.editingFinished.connect(self.changed)
            fit_chk.stateChanged.connect(self.changed)
            lo_fld.editingFinished.connect(self.changed)
            hi_fld.editingFinished.connect(self.changed)

            self._grid.addWidget(name_lbl, row_i, self._COL_NAME)
            self._grid.addWidget(val_fld,  row_i, self._COL_VAL)
            self._grid.addWidget(fit_chk,  row_i, self._COL_FIT)
            self._grid.addWidget(lo_fld,   row_i, self._COL_LO)
            self._grid.addWidget(hi_fld,   row_i, self._COL_HI)

            self._limit_widgets.extend([lo_fld, hi_fld])
            self._param_rows[pname] = {
                "val": val_fld, "fit": fit_chk, "lo": lo_fld, "hi": hi_fld,
            }

    def _on_shape_changed(self, shape: str):
        """Rebuild parameter rows when shape changes; preserve A, Q0, FWHM."""
        old_dict = self.get_params()
        old_dict["shape"] = shape
        # Keep existing A/Q0/FWHM values
        self._build_param_rows(old_dict)
        self.changed.emit()

    # ── Public API ────────────────────────────────────────────────────────

    @property
    def index(self) -> int:
        return self._index

    def set_index(self, i: int):
        self._index = i
        # Update header label color
        for child in self.children():
            if isinstance(child, QHBoxLayout):
                for j in range(child.count()):
                    w = child.itemAt(j).widget()
                    if isinstance(w, QLabel) and w.text().startswith("Peak "):
                        w.setText(f"Peak {i + 1}")
                        w.setStyleSheet(f"color:{_peak_color(i)};font-weight:bold;")
                        break
                break

    def get_params(self) -> Dict:
        shape  = self._shape_combo.currentText()
        pnames = _PEAK_PARAM_NAMES.get(shape, [])
        d: Dict = {"shape": shape}
        for pname in pnames:
            if pname not in self._param_rows:
                continue
            row = self._param_rows[pname]
            lo_v = row["lo"].float_value()
            hi_v = row["hi"].float_value()
            d[pname] = {
                "value": row["val"].float_value(),
                "fit":   row["fit"].isChecked(),
                "lo":    lo_v,
                "hi":    hi_v,
            }
        return d

    def set_params(self, peak_dict: Dict):
        shape = peak_dict.get("shape", self._shape_combo.currentText())
        self._shape_combo.blockSignals(True)
        self._shape_combo.setCurrentText(shape)
        self._shape_combo.blockSignals(False)
        self._build_param_rows(peak_dict)

    def toggle_limits(self, show: bool):
        """Show/hide Lo and Hi columns."""
        for w in self._limit_widgets:
            w.setVisible(show)


# ===========================================================================
# WAXSPeakFitGraphWindow
# ===========================================================================

class WAXSPeakFitGraphWindow(QWidget):
    """Two-panel LINEAR/LINEAR graph for WAXS peak fitting.

    Top panel: data scatter/line + total model (red) + per-peak overlays
               (coloured) + background (dashed grey).
    Bottom panel: normalised residuals with zero line.

    The ``add_peak_requested`` signal is emitted (with Q position) when the
    user selects "Add Peak at Q = …" from the right-click menu.
    """

    add_peak_requested           = Signal(float, float)  # Q position, I amplitude
    remove_nearest_peak_requested = Signal(float)        # Q position
    q_range_changed              = Signal(float, float)  # qmin, qmax from cursors

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena – WAXS Peak Fit")
        self.setGeometry(120, 120, 900, 700)

        self.gl = pg.GraphicsLayoutWidget()
        self.gl.setBackground('w')

        # ── Top plot: data + model ────────────────────────────────────────
        self.main_plot = self.gl.addPlot(row=0, col=0)
        self.main_plot.setLogMode(False, False)
        self.main_plot.setLabel('left',   'Intensity')
        self.main_plot.setLabel('bottom', 'Q  (Å⁻¹)')
        self.main_plot.showGrid(x=True, y=True, alpha=0.3)
        self._legend = self.main_plot.addLegend(offset=(-10, 10), labelTextSize='13pt')
        self._style_axes(self.main_plot)

        # ── Bottom plot: residuals ────────────────────────────────────────
        self.resid_plot = self.gl.addPlot(row=1, col=0)
        self.resid_plot.setLogMode(False, False)
        self.resid_plot.setLabel('left',   '(I−fit)/σ')
        self.resid_plot.setLabel('bottom', 'Q  (Å⁻¹)')
        self.resid_plot.showGrid(x=True, y=True, alpha=0.3)
        self.resid_plot.setXLink(self.main_plot)
        self._style_axes(self.resid_plot)
        self._zero_line = pg.InfiniteLine(
            pos=0, angle=0,
            pen=pg.mkPen('k', width=1, style=Qt.PenStyle.DashLine),
        )
        self.resid_plot.addItem(self._zero_line)

        # Height ratios 4:1
        ci = self.gl.ci
        ci.layout.setRowStretchFactor(0, 4)
        ci.layout.setRowStretchFactor(1, 1)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.gl)

        # ── Right-click: JPEG export + add-peak + display toggles ─────────
        self._add_jpeg_export(self.main_plot, "waxs_peakfit_main")
        self._add_jpeg_export(self.resid_plot, "waxs_peakfit_resid")
        self._last_right_click_q: Optional[float] = None
        self._setup_right_click_add_peak()
        self._setup_display_toggles()

        # ── Q-range cursors ───────────────────────────────────────────────
        _cursor_pen_lo = pg.mkPen('#3498db', width=2, style=Qt.PenStyle.DashLine)
        _cursor_pen_hi = pg.mkPen('#e74c3c', width=2, style=Qt.PenStyle.DashLine)
        self._cursor_qmin = pg.InfiniteLine(
            pos=0.1, angle=90, movable=True, pen=_cursor_pen_lo,
            label='Q_min', labelOpts={'position': 0.92, 'color': '#3498db',
                                      'fill': (240, 248, 255, 180)},
        )
        self._cursor_qmax = pg.InfiniteLine(
            pos=1.0, angle=90, movable=True, pen=_cursor_pen_hi,
            label='Q_max', labelOpts={'position': 0.92, 'color': '#e74c3c',
                                      'fill': (255, 245, 245, 180)},
        )
        self.main_plot.addItem(self._cursor_qmin)
        self.main_plot.addItem(self._cursor_qmax)
        self._cursor_qmin.sigPositionChanged.connect(self._on_cursor_moved)
        self._cursor_qmax.sigPositionChanged.connect(self._on_cursor_moved)

        # ── Data display state ────────────────────────────────────────────
        self._data_as_lines: bool = False
        self._show_error_bars: bool = True
        self._dI_stored: Optional[np.ndarray] = None  # saved for error bar rebuild

        # ── Data item (added/removed when data is loaded) ─────────────────
        self._data_item:  Optional[pg.PlotDataItem] = None
        self._err_item:   Optional[pg.PlotDataItem] = None

        # ── Persistent overlay items — created once, updated via setData() ─
        # Keeping items permanently in the scene and calling setData() is
        # significantly faster than removing and re-adding PlotDataItems on
        # every model redraw, because it avoids scene-graph layout passes.
        _ex: np.ndarray = np.array([], dtype=float)
        _ey: np.ndarray = np.array([], dtype=float)
        self._model_item: pg.PlotDataItem = self.main_plot.plot(
            _ex, _ey, pen=pg.mkPen('#e74c3c', width=2),
        )
        self._bg_item: pg.PlotDataItem = self.main_plot.plot(
            _ex, _ey,
            pen=pg.mkPen('#7f8c8d', width=1.5, style=Qt.PenStyle.DashLine),
        )
        self._resid_item: pg.PlotDataItem = self.resid_plot.plot(
            _ex, _ey,
            pen=None, symbol='o', symbolSize=3,
            symbolPen=pg.mkPen('#7f8c8d', width=1),
            symbolBrush=pg.mkBrush('#7f8c8d'),
        )
        # Peak item pool — grows on demand; excess slots hidden via setData([])
        self._peak_items:       List[pg.PlotDataItem] = []
        self._peak_label_items: List[pg.TextItem]     = []

    # ── Internal helpers ─────────────────────────────────────────────────

    @staticmethod
    def _style_axes(plot: pg.PlotItem):
        for side in ('left', 'bottom'):
            ax = plot.getAxis(side)
            ax.setPen(pg.mkPen('k'))
            ax.setTextPen(pg.mkPen('k'))
            ax.enableAutoSIPrefix(False)

    def _add_jpeg_export(self, plot: pg.PlotItem, stem: str):
        """Add 'Save as JPEG…' to the ViewBox right-click menu."""
        vb = plot.getViewBox()
        vb.menu.addSeparator()
        act = vb.menu.addAction("Save as JPEG…")
        act.triggered.connect(
            lambda checked=False, p=plot, s=stem: self._save_jpeg(p, s)
        )

    def _save_jpeg(self, plot: pg.PlotItem, stem: str):
        from pyqtgraph.exporters import ImageExporter
        path, _ = QFileDialog.getSaveFileName(
            self, "Save as JPEG",
            str(Path.home() / f"{stem}.jpg"),
            "JPEG Images (*.jpg *.jpeg);;All Files (*)",
        )
        if not path:
            return
        try:
            exp = ImageExporter(plot)
            exp.parameters()['width'] = 1600
            exp.export(path)
        except Exception as exc:
            QMessageBox.warning(self, "Export failed", str(exc))

    # ── Cursor helpers ────────────────────────────────────────────────────

    def _on_cursor_moved(self):
        qmin = float(self._cursor_qmin.value())
        qmax = float(self._cursor_qmax.value())
        self.q_range_changed.emit(min(qmin, qmax), max(qmin, qmax))

    def get_q_range(self) -> tuple:
        """Return (qmin, qmax) from cursor positions."""
        qmin = float(self._cursor_qmin.value())
        qmax = float(self._cursor_qmax.value())
        return min(qmin, qmax), max(qmin, qmax)

    def set_q_range(self, qmin: float, qmax: float):
        """Move cursors to given Q range."""
        self._cursor_qmin.setValue(min(qmin, qmax))
        self._cursor_qmax.setValue(max(qmin, qmax))

    def _setup_right_click_add_peak(self):
        """Add peak-management entries to the ViewBox right-click menu.

        Intercepts ``raiseContextMenu`` (called just before the menu is
        shown) to capture the exact click position and update action labels.
        """
        vb = self.main_plot.getViewBox()
        vb.menu.addSeparator()
        self._add_peak_action    = vb.menu.addAction("Add Peak at Q = ?")
        self._rm_peak_action     = vb.menu.addAction("Remove Nearest Peak")
        self._last_click_pos     = (0.1, 1.0)   # (Q, I) updated on each right-click

        # Triggered: read stored position (set by the patched raiseContextMenu)
        self._add_peak_action.triggered.connect(
            lambda: self.add_peak_requested.emit(*self._last_click_pos)
        )
        self._rm_peak_action.triggered.connect(
            lambda: self.remove_nearest_peak_requested.emit(self._last_click_pos[0])
        )

        # Patch raiseContextMenu so we capture position BEFORE the menu is shown
        _orig_raise = vb.raiseContextMenu

        def _patched_raise(ev):
            pos = vb.mapToView(ev.pos())
            q   = float(pos.x())
            I   = float(pos.y())
            self._last_click_pos = (q, I)
            self._add_peak_action.setText(f"Add Peak at Q = {q:.4f} Å⁻¹")
            _orig_raise(ev)

        vb.raiseContextMenu = _patched_raise

    def _setup_display_toggles(self):
        """Add data-display toggle actions to the main plot right-click menu."""
        vb = self.main_plot.getViewBox()
        vb.menu.addSeparator()
        self._data_mode_action = vb.menu.addAction("Show data as line")
        self._errbar_action    = vb.menu.addAction("Hide error bars")
        self._data_mode_action.triggered.connect(self._toggle_data_mode)
        self._errbar_action.triggered.connect(self._toggle_error_bars)

    def _toggle_data_mode(self):
        """Switch data display between scatter points and connected line."""
        self._data_as_lines = not self._data_as_lines
        if self._data_item is not None:
            if self._data_as_lines:
                self._data_item.setData(
                    pen=pg.mkPen('#2c3e50', width=1.5),
                    symbol=None,
                )
                self._data_mode_action.setText("Show data as points")
            else:
                self._data_item.setData(
                    pen=None, symbol='o', symbolSize=4,
                    symbolPen=pg.mkPen('#2c3e50', width=1),
                    symbolBrush=pg.mkBrush('#2c3e50'),
                )
                self._data_mode_action.setText("Show data as line")

    def _toggle_error_bars(self):
        """Show or hide error bar lines."""
        self._show_error_bars = not self._show_error_bars
        if self._err_item is not None:
            self._err_item.setVisible(self._show_error_bars)
        self._errbar_action.setText(
            "Hide error bars" if self._show_error_bars else "Show error bars"
        )

    # ── Plot API ─────────────────────────────────────────────────────────

    def _clear_overlay_data(self):
        """Clear all overlay curves without removing items from the scene graph."""
        _ex: np.ndarray = np.array([], dtype=float)
        _ey: np.ndarray = np.array([], dtype=float)
        self._model_item.setData(_ex, _ey)
        self._bg_item.setData(_ex, _ey)
        for item in self._peak_items:
            item.setData(_ex, _ey)
        for lbl in self._peak_label_items:
            lbl.setVisible(False)

    @staticmethod
    def _make_errbar_data(q: np.ndarray, I: np.ndarray,
                          dI: np.ndarray) -> tuple:
        """Build NaN-segmented vertical error bar arrays for linear/linear scale."""
        n = len(q)
        xs = np.empty(3 * n)
        ys = np.empty(3 * n)
        xs[0::3] = q
        xs[1::3] = q
        xs[2::3] = np.nan
        ys[0::3] = I - dI
        ys[1::3] = I + dI
        ys[2::3] = np.nan
        return xs, ys

    def plot_data(self, q: np.ndarray, I: np.ndarray, dI: Optional[np.ndarray] = None,
                  label: str = "Data"):
        """Plot raw data (scatter) and optional error bars."""
        if self._data_item is not None:
            self.main_plot.removeItem(self._data_item)
        if self._err_item is not None:
            self.main_plot.removeItem(self._err_item)
            self._err_item = None
        q_, I_ = np.asarray(q, float), np.asarray(I, float)
        mask = np.isfinite(q_) & np.isfinite(I_)
        self._data_item = self.main_plot.plot(
            q_[mask], I_[mask],
            pen=None, symbol='o', symbolSize=4,
            symbolPen=pg.mkPen('#2c3e50', width=1),
            symbolBrush=pg.mkBrush('#2c3e50'),
            name=label,
        )
        # Error bars (thin gray vertical segments, NaN-separated)
        if dI is not None:
            dI_ = np.asarray(dI, float)
            em = mask & np.isfinite(dI_) & (dI_ > 0)
            if em.any():
                self._dI_stored = dI_
                ex, ey = self._make_errbar_data(q_[em], I_[em], dI_[em])
                self._err_item = self.main_plot.plot(
                    ex, ey,
                    pen=pg.mkPen('#7f8c8d', width=0.8),
                    connect='finite',
                )
                self._err_item.setVisible(self._show_error_bars)
        # Style the legend entry: bold, larger, dark
        try:
            for _, lbl_item in self._legend.items:
                f = QFont()
                f.setPointSize(13)
                f.setBold(True)
                f.setWeight(QFont.Weight.Black)
                lbl_item.item.setFont(f)
                lbl_item.item.setDefaultTextColor(pg.mkColor('#1a1a2e'))
        except Exception:
            pass
        # Set x-range and cursor positions to data range
        qv = q_[mask & (q_ > 0)]
        if len(qv) >= 2:
            self.main_plot.setXRange(float(qv.min()), float(qv.max()), padding=0.02)
            self.set_q_range(float(qv.min()), float(qv.max()))
            self.q_range_changed.emit(float(qv.min()), float(qv.max()))

    def plot_model(
        self,
        q: np.ndarray,
        I_model: np.ndarray,
        I_bg: Optional[np.ndarray],
        peak_curves: List[Dict],
    ):
        """Overlay model, background, and per-peak curves via setData() (no add/remove)."""
        _ex: np.ndarray = np.array([], dtype=float)
        _ey: np.ndarray = np.array([], dtype=float)
        q_ = np.asarray(q, float)
        mask = np.isfinite(q_)

        # Total model (red)
        if I_model is not None:
            Im = np.asarray(I_model, float)
            self._model_item.setData(q_[mask], Im[mask])
        else:
            self._model_item.setData(_ex, _ey)

        # Background (dashed grey)
        if I_bg is not None:
            Ib = np.asarray(I_bg, float)
            self._bg_item.setData(q_[mask], Ib[mask])
        else:
            self._bg_item.setData(_ex, _ey)

        # Per-peak curves + number labels — grow pool on demand
        n_peaks = len(peak_curves)
        while len(self._peak_items) < n_peaks:
            idx = len(self._peak_items)
            col = _peak_color(idx)
            item = self.main_plot.plot(_ex, _ey, pen=pg.mkPen(col, width=1.5))
            lbl = pg.TextItem(text='', color=col, anchor=(0.5, 1.2))
            f = QFont()
            f.setBold(True)
            f.setPointSize(20)
            lbl.setFont(f)
            self.main_plot.addItem(lbl)
            self._peak_items.append(item)
            self._peak_label_items.append(lbl)

        for i, pc in enumerate(peak_curves):
            item = self._peak_items[i]
            lbl  = self._peak_label_items[i]
            col  = _peak_color(i)
            item.setPen(pg.mkPen(col, width=1.5))  # colour may shift if peaks reordered
            qp = pc.get("q")
            Ip = pc.get("I")
            if qp is None or Ip is None:
                item.setData(_ex, _ey)
                lbl.setVisible(False)
                continue
            qp_, Ip_ = np.asarray(qp, float), np.asarray(Ip, float)
            pm = np.isfinite(qp_) & np.isfinite(Ip_)
            item.setData(qp_[pm], Ip_[pm])
            if pm.any():
                top_idx = int(np.argmax(Ip_[pm]))
                lbl.setText(str(i + 1))
                lbl.setColor(col)
                lbl.setPos(float(qp_[pm][top_idx]), float(Ip_[pm][top_idx]))
                lbl.setVisible(True)
            else:
                lbl.setVisible(False)

        # Hide excess pool slots
        for i in range(n_peaks, len(self._peak_items)):
            self._peak_items[i].setData(_ex, _ey)
            self._peak_label_items[i].setVisible(False)

    def plot_residuals(self, q: np.ndarray, residuals: np.ndarray):
        """Plot normalised residuals in the bottom panel."""
        q_  = np.asarray(q, float)
        res = np.asarray(residuals, float)
        mask = np.isfinite(q_) & np.isfinite(res)
        self._resid_item.setData(q_[mask], res[mask])

    def clear_all(self):
        """Reset graph without calling plot.clear() — preserves persistent items."""
        if self._data_item is not None:
            self.main_plot.removeItem(self._data_item)
            self._data_item = None
        if self._err_item is not None:
            self.main_plot.removeItem(self._err_item)
            self._err_item = None
        self._dI_stored = None
        self._legend.clear()
        self._clear_overlay_data()
        _ex: np.ndarray = np.array([], dtype=float)
        _ey: np.ndarray = np.array([], dtype=float)
        self._resid_item.setData(_ex, _ey)


# ===========================================================================
# _FixLimitsDialog — fraction-based limit setter
# ===========================================================================

class _FixLimitsDialog(QDialog):
    """Dialog for setting ±fraction limits on all fit parameters."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Fix Limits — Fractions")
        self.setFixedWidth(340)
        layout = QGridLayout(self)
        layout.setHorizontalSpacing(8)
        layout.setVerticalSpacing(4)

        layout.addWidget(_label("Parameter", bold=True), 0, 0)
        layout.addWidget(_label("±Fraction", bold=True), 0, 1)
        layout.addWidget(_label("Example at val=1.0", bold=True), 0, 2)

        _params = [
            ("Q0  (position)",  "Q0",    0.05),
            ("FWHM  (width)",   "FWHM",  0.10),
            ("A  (amplitude)",  "A",     0.50),
            ("Other params",    "other", 0.30),
        ]
        self._spins: Dict[str, QDoubleSpinBox] = {}
        for row, (label, key, default) in enumerate(_params, start=1):
            layout.addWidget(_label(label), row, 0)
            sp = QDoubleSpinBox()
            sp.setRange(0.001, 10.0)
            sp.setDecimals(3)
            sp.setSingleStep(0.01)
            sp.setValue(default)
            sp.setFixedWidth(80)
            layout.addWidget(sp, row, 1)
            example_lbl = QLabel(f"[{1-default:.2f}, {1+default:.2f}]")
            example_lbl.setStyleSheet("color:#777;font-size:9pt;")
            sp.valueChanged.connect(
                lambda v, lbl=example_lbl: lbl.setText(f"[{1-v:.2f}, {1+v:.2f}]")
            )
            layout.addWidget(example_lbl, row, 2)
            self._spins[key] = sp

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons, len(_params) + 1, 0, 1, 3)

    def get_fracs(self) -> Dict[str, float]:
        return {k: sp.value() for k, sp in self._spins.items()}


# ===========================================================================
# WAXSPeakFitPanel — main panel
# ===========================================================================

class WAXSPeakFitPanel(QWidget):
    """Main WAXS Peak-Fitting GUI: left controls + right graph."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena – WAXS Peak Fit")
        self.setMinimumSize(1000, 750)
        self.resize(1200, 900)

        # ── Application state ─────────────────────────────────────────────
        from pyirena.state.state_manager import StateManager
        self._state_mgr = StateManager()

        # ── Working data ──────────────────────────────────────────────────
        self._q:   Optional[np.ndarray] = None
        self._I:   Optional[np.ndarray] = None
        self._dI:  Optional[np.ndarray] = None
        self._filepath: Optional[Path]  = None
        self._is_nxcansas: bool         = False

        self._pre_fit_bg_params: Optional[Dict] = None   # for Revert
        self._pre_fit_peaks:     Optional[List] = None

        # ── Graph window (embedded on right side) ─────────────────────────
        self._graph = WAXSPeakFitGraphWindow()
        self._graph.add_peak_requested.connect(self._add_peak_at_q)
        self._graph.remove_nearest_peak_requested.connect(self._remove_nearest_peak)
        self._graph.q_range_changed.connect(self._on_q_range_changed)

        # ── Debounce timer for auto-recalculation (wheel events) ──────────
        self._debounce_timer = QTimer(self)
        self._debounce_timer.setSingleShot(True)
        self._debounce_timer.setInterval(100)
        self._debounce_timer.timeout.connect(self._graph_model)

        # ── Build layout ──────────────────────────────────────────────────
        splitter = QSplitter(Qt.Orientation.Horizontal, self)
        splitter.setChildrenCollapsible(False)

        left_container = QWidget()
        left_container.setFixedWidth(440)
        left_container.setMaximumWidth(440)
        self._left_layout = QVBoxLayout(left_container)
        self._left_layout.setContentsMargins(4, 4, 4, 4)
        self._left_layout.setSpacing(4)
        self._build_left_panel()

        splitter.addWidget(left_container)
        splitter.addWidget(self._graph)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

        root = QVBoxLayout(self)
        root.setContentsMargins(4, 4, 4, 4)
        root.addWidget(splitter)

        # ── Apply saved state ─────────────────────────────────────────────
        self._apply_state(self._state_mgr.get("waxs_peakfit", default={}))

    # ===========================================================================
    # Left panel construction
    # ===========================================================================

    def _build_left_panel(self):
        ll = self._left_layout

        # ── Title + No limits ─────────────────────────────────────────────
        title_row = QHBoxLayout()
        title_lbl = _label("WAXS Peak Fit", bold=True, size=12)
        title_lbl.setStyleSheet(
            "background:#2c3e50;color:white;padding:4px 6px;border-radius:3px;"
        )
        title_row.addWidget(title_lbl)
        title_row.addStretch()
        self._no_limits_chk = QCheckBox("No limits?")
        self._no_limits_chk.stateChanged.connect(self._on_no_limits_toggled)
        title_row.addWidget(self._no_limits_chk)
        ll.addLayout(title_row)

        # ── Q fit range (cursor positions, read-only display) ─────────────
        qr_row = QHBoxLayout()
        qr_row.addWidget(_label("Fit Q range:"))
        self._qmin_label = QLabel("–")
        self._qmin_label.setStyleSheet("color:#3498db;")
        self._qmax_label = QLabel("–")
        self._qmax_label.setStyleSheet("color:#e74c3c;")
        qr_row.addWidget(self._qmin_label)
        qr_row.addWidget(_label("–"))
        qr_row.addWidget(self._qmax_label)
        qr_row.addWidget(_label("Å⁻¹"))
        qr_row.addStretch()
        ll.addLayout(qr_row)

        # ── Background section ────────────────────────────────────────────
        bg_box = QGroupBox("Background")
        bg_layout = QVBoxLayout(bg_box)
        bg_layout.setSpacing(2)

        bg_shape_row = QHBoxLayout()
        bg_shape_row.addWidget(_label("Shape:"))
        self._bg_combo = QComboBox()
        self._bg_combo.addItems(BG_SHAPES)
        self._bg_combo.currentTextChanged.connect(self._on_bg_shape_changed)
        bg_shape_row.addWidget(self._bg_combo)
        bg_shape_row.addStretch()
        bg_layout.addLayout(bg_shape_row)

        self._bg_grid = QGridLayout()
        self._bg_grid.setHorizontalSpacing(4)
        self._bg_grid.setVerticalSpacing(2)
        bg_layout.addLayout(self._bg_grid)
        self._bg_param_rows: Dict[str, dict] = {}
        self._bg_limit_widgets: List[QWidget] = []
        ll.addWidget(bg_box)

        # ── Peak finding section ──────────────────────────────────────────
        pf_box = QGroupBox("Peak Finding")
        pf_layout = QGridLayout(pf_box)
        pf_layout.setHorizontalSpacing(6)
        pf_layout.setVerticalSpacing(2)

        def _dspin(val, lo, hi, dec=4, step=0.001):
            s = QDoubleSpinBox()
            s.setRange(lo, hi)
            s.setDecimals(dec)
            s.setSingleStep(step)
            s.setValue(val)
            s.setFixedWidth(90)
            return s

        self._pf_prom  = _dspin(0.05,  0.0,  1.0, dec=3, step=0.01)
        self._pf_minfwhm = _dspin(0.001, 0.0, 10.0, dec=4, step=0.0005)
        self._pf_maxfwhm = _dspin(0.5,  0.001,100.0,dec=3, step=0.05)
        self._pf_mindist = _dspin(0.005, 0.0, 10.0, dec=4, step=0.001)
        self._pf_sgwin  = _dspin(15.0,  1.0, 80.0, dec=1, step=1.0)

        rows = [
            ("Min. Prominence (frac.)", self._pf_prom),
            ("Min. FWHM (Å⁻¹)",        self._pf_minfwhm),
            ("Max. FWHM (Å⁻¹)",        self._pf_maxfwhm),
            ("Min. Distance (Å⁻¹)",     self._pf_mindist),
            ("BG smooth window (%)",    self._pf_sgwin),
        ]
        for i, (lbl_txt, widget) in enumerate(rows):
            pf_layout.addWidget(_label(lbl_txt), i, 0)
            pf_layout.addWidget(widget, i, 1)

        self._find_btn = QPushButton("Find Peaks")
        self._find_btn.setStyleSheet(_BTN_FIND)
        self._find_btn.setMinimumHeight(28)
        self._find_btn.clicked.connect(self._find_peaks)
        pf_layout.addWidget(self._find_btn, len(rows), 0, 1, 2)
        ll.addWidget(pf_box)

        # ── Peaks scroll area ─────────────────────────────────────────────
        peaks_outer = QGroupBox("Peaks")
        peaks_outer_layout = QVBoxLayout(peaks_outer)
        peaks_outer_layout.setSpacing(2)

        self._peaks_scroll = QScrollArea()
        self._peaks_scroll.setWidgetResizable(True)
        self._peaks_scroll.setMinimumHeight(180)
        self._peaks_container = QWidget()
        self._peaks_vbox = QVBoxLayout(self._peaks_container)
        self._peaks_vbox.setSpacing(4)
        self._peaks_vbox.setContentsMargins(2, 2, 2, 2)
        self._peaks_vbox.addStretch()
        self._peaks_scroll.setWidget(self._peaks_container)
        peaks_outer_layout.addWidget(self._peaks_scroll)

        peak_btn_row = QHBoxLayout()
        add_btn = QPushButton("Add Peak Manually")
        add_btn.setStyleSheet(_BTN_ADD)
        add_btn.setMinimumHeight(26)
        add_btn.clicked.connect(lambda: self._add_peak_at_q(None, None))
        peak_btn_row.addWidget(add_btn)

        sort_btn = QPushButton("Sort by Q")
        sort_btn.setStyleSheet(_BTN_SAVE)
        sort_btn.setMinimumHeight(26)
        sort_btn.clicked.connect(self._sort_peaks_by_q)
        peak_btn_row.addWidget(sort_btn)
        peaks_outer_layout.addLayout(peak_btn_row)
        ll.addWidget(peaks_outer)

        # ── Fit controls ──────────────────────────────────────────────────
        fit_row = QHBoxLayout()
        self._graph_btn = QPushButton("Graph Model")
        self._graph_btn.setStyleSheet(_BTN_GRAPH)
        self._graph_btn.setMinimumHeight(28)
        self._graph_btn.setMaximumWidth(120)
        self._graph_btn.clicked.connect(self._graph_model)
        fit_row.addWidget(self._graph_btn)
        fit_row.addStretch()

        self._fit_btn = QPushButton("Fit")
        self._fit_btn.setStyleSheet(_BTN_FIT)
        self._fit_btn.setMinimumHeight(28)
        self._fit_btn.setMaximumWidth(90)
        self._fit_btn.clicked.connect(self._run_fit)
        fit_row.addWidget(self._fit_btn)

        self._revert_btn = QPushButton("Revert")
        self._revert_btn.setStyleSheet(_BTN_REVERT)
        self._revert_btn.setMinimumHeight(28)
        self._revert_btn.setMaximumWidth(90)
        self._revert_btn.setEnabled(False)
        self._revert_btn.clicked.connect(self._revert)
        fit_row.addWidget(self._revert_btn)
        ll.addLayout(fit_row)

        # ── Fix Limits + Reset (same row) ─────────────────────────────────
        extra_row = QHBoxLayout()
        fix_limits_btn = QPushButton("Fix Limits")
        fix_limits_btn.setStyleSheet(_BTN_FIXLIM)
        fix_limits_btn.setMinimumHeight(28)
        fix_limits_btn.clicked.connect(self._fix_limits)
        extra_row.addWidget(fix_limits_btn)

        reset_btn = QPushButton("Reset to Defaults")
        reset_btn.setStyleSheet(_BTN_RESET)
        reset_btn.setMinimumHeight(28)
        reset_btn.clicked.connect(self._reset_defaults)
        extra_row.addWidget(reset_btn)
        ll.addLayout(extra_row)

        # ── Results section ───────────────────────────────────────────────
        res_lbl = _label("Results", bold=True)
        res_lbl.setStyleSheet("color:#3498db;margin-top:4px;")
        ll.addWidget(res_lbl)

        row1 = QHBoxLayout()
        for txt, style, slot in [
            ("Save State",        _BTN_SAVE,   self._save_state),
            ("Store in File",     _BTN_STORE,  self._store_in_file),
            ("Results to graphs", _BTN_RGRAPH, self._results_to_graphs),
        ]:
            b = QPushButton(txt)
            b.setStyleSheet(style)
            b.setMinimumHeight(26)
            b.clicked.connect(slot)
            row1.addWidget(b)
        ll.addLayout(row1)

        row2 = QHBoxLayout()
        for txt, style, slot in [
            ("Export Parameters", _BTN_EXPORT, self._export_params),
            ("Import Parameters", _BTN_EXPORT, self._import_params),
        ]:
            b = QPushButton(txt)
            b.setStyleSheet(style)
            b.setMinimumHeight(26)
            b.clicked.connect(slot)
            row2.addWidget(b)
        ll.addLayout(row2)

        # ── Status label ──────────────────────────────────────────────────
        self._status = QLabel("")
        self._status.setStyleSheet(
            "font-size:10pt;color:#555;border-top:1px solid #ccc;padding-top:2px;"
        )
        self._status.setWordWrap(True)
        self._status.setMaximumHeight(60)
        ll.addWidget(self._status)
        ll.addStretch()

        # ── Populate background parameter grid ────────────────────────────
        self._rebuild_bg_grid("Constant")

    # ===========================================================================
    # Background parameter grid
    # ===========================================================================

    def _rebuild_bg_grid(self, bg_shape: str, saved_params: Optional[Dict] = None):
        """Rebuild the background parameter grid for a new shape."""
        # Clear old widgets
        for row_d in self._bg_param_rows.values():
            for w in row_d.values():
                if hasattr(w, "setParent"):
                    w.setParent(None)
        self._bg_param_rows = {}
        self._bg_limit_widgets = []

        # Clear grid
        while self._bg_grid.count():
            item = self._bg_grid.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        # Column headers (order: Param | Value | Fit? | Lo limit | Hi limit)
        for col, txt in [(0, "Param"), (1, "Value"), (2, "Fit?"),
                         (3, "Lo limit"), (4, "Hi limit")]:
            h = _label(txt, bold=True)
            self._bg_grid.addWidget(h, 0, col)
            if col in (3, 4):
                self._bg_limit_widgets.append(h)

        names = bg_param_names(bg_shape)
        defaults = default_bg_params(bg_shape)

        for row_i, name in enumerate(names, start=1):
            pd = (saved_params or {}).get(name, defaults[name])
            val = float(pd.get("value", 0.0))
            fit = bool(pd.get("fit", True))
            lo  = pd.get("lo")
            hi  = pd.get("hi")

            name_lbl = _label(name)
            val_fld  = _ValidatedField(val)
            val_fld.setFixedWidth(104)
            fit_chk  = QCheckBox()
            fit_chk.setChecked(fit)
            lo_fld  = _ValidatedField(lo if lo is not None else -1e6, wheel_enabled=False)
            lo_fld.setFixedWidth(81)
            hi_fld  = _ValidatedField(hi if hi is not None else 1e6, wheel_enabled=False)
            hi_fld.setFixedWidth(81)

            # Auto-expand limits and trigger model redraw on value change
            _connect_limit_check(val_fld, lo_fld, hi_fld)
            val_fld.editingFinished.connect(self._request_model_update)
            fit_chk.stateChanged.connect(self._request_model_update)

            self._bg_grid.addWidget(name_lbl, row_i, 0)
            self._bg_grid.addWidget(val_fld,  row_i, 1)
            self._bg_grid.addWidget(fit_chk,  row_i, 2)
            self._bg_grid.addWidget(lo_fld,   row_i, 3)
            self._bg_grid.addWidget(hi_fld,   row_i, 4)

            self._bg_limit_widgets.extend([lo_fld, hi_fld])
            self._bg_param_rows[name] = {
                "val": val_fld, "fit": fit_chk, "lo": lo_fld, "hi": hi_fld,
            }

        # Apply current no-limits toggle
        show = not self._no_limits_chk.isChecked()
        for w in self._bg_limit_widgets:
            w.setVisible(show)

    # ===========================================================================
    # Peaks list management
    # ===========================================================================

    @property
    def _peak_rows(self) -> List[PeakRowWidget]:
        rows = []
        layout = self._peaks_vbox
        for i in range(layout.count()):
            item = layout.itemAt(i)
            if item and isinstance(item.widget(), PeakRowWidget):
                rows.append(item.widget())
        return rows

    def _add_peak_row(self, peak_dict: Dict):
        idx    = len(self._peak_rows)
        row    = PeakRowWidget(idx, peak_dict)
        row.changed.connect(self._request_model_update)
        row.remove_requested.connect(self._remove_peak_row)

        # Apply no-limits state
        row.toggle_limits(not self._no_limits_chk.isChecked())

        # Insert before the trailing stretch
        layout = self._peaks_vbox
        stretch_idx = layout.count() - 1
        layout.insertWidget(stretch_idx, row)

    def _remove_peak_row_silent(self, row: PeakRowWidget):
        """Remove a peak row without triggering a model redraw."""
        self._peaks_vbox.removeWidget(row)
        row.setParent(None)
        row.deleteLater()
        # Re-number remaining rows
        for i, r in enumerate(self._peak_rows):
            r.set_index(i)

    def _remove_peak_row(self, row: PeakRowWidget):
        self._remove_peak_row_silent(row)
        self._graph_model()

    def _clear_peaks(self):
        for row in list(self._peak_rows):
            self._peaks_vbox.removeWidget(row)
            row.setParent(None)
            row.deleteLater()

    # ===========================================================================
    # Slots: background shape change + no-limits
    # ===========================================================================

    def _on_bg_shape_changed(self, shape: str):
        # Preserve existing values for coefficients that still exist
        old_params = self._get_bg_params()
        self._rebuild_bg_grid(shape, old_params)

    def _on_no_limits_toggled(self, state):
        show = (state == Qt.CheckState.Unchecked.value or state == 0)
        for w in self._bg_limit_widgets:
            w.setVisible(show)
        for row in self._peak_rows:
            row.toggle_limits(show)

    # ===========================================================================
    # Parameter collection helpers
    # ===========================================================================

    def _get_bg_params(self) -> Dict:
        shape = self._bg_combo.currentText()
        result = {}
        for name in bg_param_names(shape):
            if name not in self._bg_param_rows:
                continue
            row = self._bg_param_rows[name]
            result[name] = {
                "value": row["val"].float_value(),
                "fit":   row["fit"].isChecked(),
                "lo":    row["lo"].float_value(),
                "hi":    row["hi"].float_value(),
            }
        return result

    def _get_peaks(self) -> List[Dict]:
        return [r.get_params() for r in self._peak_rows]

    # ===========================================================================
    # Slots: find / add peaks
    # ===========================================================================

    def _find_peaks(self):
        if self._q is None or self._I is None:
            self._set_status("No data loaded.", error=True)
            return
        # Restrict search to cursor Q range
        qmin, qmax = self._graph.get_q_range()
        mask = (self._q >= qmin) & (self._q <= qmax) & np.isfinite(self._I)
        q_search = self._q[mask]
        I_search = self._I[mask]
        if len(q_search) < 10:
            self._set_status(
                f"Too few points in Q range [{qmin:.4g}, {qmax:.4g}].", error=True
            )
            return
        self._set_status("Finding peaks…")
        peaks = find_peaks_in_data(
            q_search, I_search,
            prominence_frac=self._pf_prom.value(),
            min_fwhm=self._pf_minfwhm.value(),
            max_fwhm=self._pf_maxfwhm.value(),
            min_distance=self._pf_mindist.value(),
            sg_window_frac=self._pf_sgwin.value() / 100.0,
        )
        self._clear_peaks()
        for p in peaks:
            self._add_peak_row(p)
        self._set_status(f"Found {len(peaks)} peak(s) in Q=[{qmin:.4g}, {qmax:.4g}].")
        self._graph_model()

    def _add_peak_at_q(self, q0: Optional[float], amplitude: Optional[float] = None):
        """Add a new peak row (optionally centred at q0 with given amplitude)."""
        if q0 is None:
            q0 = float(np.median(self._q)) if self._q is not None else 0.1

        # Determine amplitude A
        if amplitude is not None and np.isfinite(amplitude) and amplitude > 0:
            A = float(amplitude)
        elif self._q is not None and self._I is not None:
            idx = int(np.argmin(np.abs(self._q - q0)))
            A = max(float(self._I[idx]), 1e-30)
        else:
            A = 1.0

        # Estimate FWHM: prefer median of existing peaks, then fraction of Q range
        existing_fwhms = []
        for r in self._peak_rows:
            p = r.get_params()
            v = p.get("FWHM", {})
            if isinstance(v, dict):
                fv = float(v.get("value", 0))
            else:
                fv = float(v)
            if fv > 0:
                existing_fwhms.append(fv)

        if existing_fwhms:
            fwhm = float(np.median(existing_fwhms))
        elif self._q is not None and len(self._q) >= 2:
            fwhm = max(0.005, 0.02 * (float(self._q.max()) - float(self._q.min())))
        else:
            fwhm = 0.01

        pd = default_peak_params("Gauss", Q0=q0, A=A, FWHM=fwhm)
        self._add_peak_row(pd)
        self._graph_model()

    # ===========================================================================
    # Slots: Graph / Fit / Revert / Reset
    # ===========================================================================

    def _compute_model(self):
        """Evaluate model from current GUI state. Returns (I_model, I_bg, peak_curves)."""
        if self._q is None:
            return None, None, []
        q         = self._q
        bg_shape  = self._bg_combo.currentText()
        bg_params = self._get_bg_params()
        peaks     = self._get_peaks()

        bg_names = bg_param_names(bg_shape)
        coeffs   = [float(bg_params[n]["value"]) for n in bg_names if n in bg_params]
        I_bg     = eval_background(q, bg_shape, coeffs)

        I_model = I_bg.copy()
        peak_curves = []
        for peak in peaks:
            shape  = peak["shape"]
            pnames = _PEAK_PARAM_NAMES.get(shape, [])
            pvals  = {pn: float(peak[pn]["value"]) for pn in pnames if pn in peak}
            I_pk   = eval_peak(q, shape, pvals)
            I_model += I_pk
            peak_curves.append({"q": q, "I": I_pk})

        return I_model, I_bg, peak_curves

    def _prune_peaks_to_range(self) -> int:
        """Remove peaks whose Q0 falls outside the cursor Q range. Returns count removed."""
        if self._q is None:
            return 0
        qmin, qmax = self._graph.get_q_range()
        to_remove = []
        for row in self._peak_rows:
            p = row.get_params()
            v = p.get("Q0", {})
            q0 = float(v.get("value", 0) if isinstance(v, dict) else v)
            if q0 < qmin or q0 > qmax:
                to_remove.append(row)
        for row in to_remove:
            self._remove_peak_row_silent(row)
        if to_remove:
            self._set_status(
                f"Removed {len(to_remove)} peak(s) outside Q range "
                f"[{qmin:.4g}, {qmax:.4g}] Å⁻¹."
            )
        return len(to_remove)

    def _request_model_update(self):
        """Debounced model update — used by auto-recalculation on value changes."""
        self._debounce_timer.start()  # restarts the timer if already running

    def _graph_model(self):
        self._debounce_timer.stop()  # cancel any pending debounce
        self._prune_peaks_to_range()
        I_model, I_bg, peak_curves = self._compute_model()
        if I_model is None:
            return  # No data — silently do nothing
        self._graph.plot_model(self._q, I_model, I_bg, peak_curves)
        # Residuals (if dI available) — NaN outside cursor Q range so they don't blow up
        if self._dI is not None and np.any(self._dI > 0):
            s = np.where(self._dI > 0, self._dI, 1.0)
            r = (self._I - I_model) / s
            qmin, qmax = self._graph.get_q_range()
            r = r.copy()
            r[(self._q < qmin) | (self._q > qmax)] = np.nan
            self._graph.plot_residuals(self._q, r)

    def _run_fit(self):
        if self._q is None:
            self._set_status("No data loaded.", error=True)
            return

        # Crop data to cursor Q range
        qmin, qmax = self._graph.get_q_range()
        mask = (self._q >= qmin) & (self._q <= qmax) & np.isfinite(self._I)
        if mask.sum() < 3:
            self._set_status(
                f"Too few data points in Q range [{qmin:.4g}, {qmax:.4g}] Å⁻¹.",
                error=True,
            )
            return
        q_fit  = self._q[mask]
        I_fit  = self._I[mask]
        dI_fit = self._dI[mask] if self._dI is not None else None

        # Remove peaks outside Q range before fitting
        self._prune_peaks_to_range()

        # Save state for revert
        self._pre_fit_bg_params = copy.deepcopy(self._get_bg_params())
        self._pre_fit_peaks     = copy.deepcopy(self._get_peaks())

        bg_shape  = self._bg_combo.currentText()
        no_limits = self._no_limits_chk.isChecked()
        engine    = WAXSPeakFitModel(bg_shape=bg_shape, peaks=[], no_limits=no_limits)

        # Show orange "Fitting…" and force Qt to repaint before the blocking call
        self._set_status("Fitting…", progress=True)
        pg.QtWidgets.QApplication.processEvents()

        try:
            result = engine.fit(
                q_fit, I_fit, dI_fit,
                bg_params=self._get_bg_params(),
                peaks=self._get_peaks(),
            )
        except Exception as exc:
            self._set_status(f"Fit error: {exc}", error=True)
            return

        if not result.get("success", False):
            self._set_status(f"Fit failed: {result.get('message', '')}", error=True)
        else:
            chi2 = result.get("reduced_chi2", float("nan"))
            self._set_status(
                f"Fit converged.  Reduced χ² = {chi2:.4g}  "
                f"(DOF = {result.get('dof', 0)})",
                success=True,
            )

        # Apply fitted parameters back to GUI
        self._apply_fit_result(result)
        self._revert_btn.setEnabled(True)

        # Re-graph
        self._graph_model()

    def _apply_fit_result(self, result: Dict):
        """Write fitted param values back into the GUI fields."""
        bg_params_new = result.get("bg_params", {})
        for name, row in self._bg_param_rows.items():
            if name in bg_params_new:
                row["val"].set_float(float(bg_params_new[name].get("value", 0.0)))

        peaks_new = result.get("peaks", [])
        for i, row_widget in enumerate(self._peak_rows):
            if i < len(peaks_new):
                row_widget.set_params(peaks_new[i])

    def _revert(self):
        if self._pre_fit_bg_params is not None:
            for name, row in self._bg_param_rows.items():
                if name in self._pre_fit_bg_params:
                    row["val"].set_float(
                        float(self._pre_fit_bg_params[name].get("value", 0.0))
                    )
        if self._pre_fit_peaks is not None:
            for i, row_widget in enumerate(self._peak_rows):
                if i < len(self._pre_fit_peaks):
                    row_widget.set_params(self._pre_fit_peaks[i])
        self._revert_btn.setEnabled(False)
        self._graph_model()

    def _reset_defaults(self):
        self._bg_combo.setCurrentText("Constant")
        self._rebuild_bg_grid("Constant")
        self._clear_peaks()
        self._graph.clear_all()
        if self._q is not None and self._I is not None:
            self._graph.plot_data(self._q, self._I, self._dI)
        self._set_status("Reset to defaults.")

    def _on_q_range_changed(self, qmin: float, qmax: float):
        """Update Q-range display labels when cursors move."""
        self._qmin_label.setText(f"{qmin:.5g}")
        self._qmax_label.setText(f"{qmax:.5g}")

    def _fix_limits(self):
        """Show Fix Limits dialog and apply ± fraction limits to all params."""
        dlg = _FixLimitsDialog(self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        fracs = dlg.get_fracs()

        # Determine a data-scale reference for zero-valued params
        if self._I is not None and len(self._I) > 0:
            _data_ref = max(float(np.nanmax(np.abs(self._I))), 1.0)
        else:
            _data_ref = 1.0

        def _apply_limits(val, frac, row_dict):
            if abs(val) < 1e-10:
                # val ≈ 0: use fraction of data range as absolute half-range
                half = _data_ref * frac
                row_dict["lo"].set_float(-half)
                row_dict["hi"].set_float(half)
            else:
                row_dict["lo"].set_float(val * (1.0 - frac))
                row_dict["hi"].set_float(val * (1.0 + frac))

        # Apply to background params
        for name, row in self._bg_param_rows.items():
            _apply_limits(row["val"].float_value(), fracs.get("other", 0.30), row)

        # Apply to peak params
        for peak_row in self._peak_rows:
            for pname, row in peak_row._param_rows.items():
                frac = fracs.get(pname, fracs.get("other", 0.30))
                _apply_limits(row["val"].float_value(), frac, row)

        self._set_status("Limits set to ±fraction of current values.")

    def _sort_peaks_by_q(self):
        """Re-order peak rows by ascending Q0 value."""
        rows = self._peak_rows
        if len(rows) <= 1:
            return
        params_list = [r.get_params() for r in rows]

        def _q0(p):
            v = p.get("Q0", {})
            return float(v.get("value", 0) if isinstance(v, dict) else v)

        sorted_params = sorted(params_list, key=_q0)
        self._clear_peaks()
        for p in sorted_params:
            self._add_peak_row(p)
        self._graph_model()
        self._set_status("Peaks sorted by Q position.")

    def _remove_nearest_peak(self, q: float):
        """Remove the peak whose Q0 is closest to q."""
        rows = self._peak_rows
        if not rows:
            return

        def _q0(r):
            p = r.get_params()
            v = p.get("Q0", {})
            return float(v.get("value", 0) if isinstance(v, dict) else v)

        nearest = min(rows, key=lambda r: abs(_q0(r) - q))
        self._remove_peak_row(nearest)
        self._graph_model()

    # ===========================================================================
    # Slots: Results buttons
    # ===========================================================================

    def _save_state(self):
        state = self._get_current_state()
        self._state_mgr.update("waxs_peakfit", state)
        if self._state_mgr.save():
            self._set_status("State saved.")
        else:
            self._set_status("Failed to save state.", error=True)

    def _store_in_file(self):
        if self._filepath is None:
            QMessageBox.warning(self, "No file", "No data file is open.")
            return
        # Need a fit result to store; run model evaluation first
        I_model, I_bg, _ = self._compute_model()
        if I_model is None:
            self._set_status("Compute model first.", error=True)
            return
        bg_shape  = self._bg_combo.currentText()
        bg_params = self._get_bg_params()
        peaks     = self._get_peaks()

        dI = self._dI
        s  = dI if dI is not None else np.ones_like(self._I)
        s  = np.where(s > 0, s, 1.0)
        residuals = (self._I - I_model) / s

        result = {
            "bg_shape":    bg_shape,
            "bg_params":   bg_params,
            "bg_params_std": {},
            "peaks":       peaks,
            "peaks_std":   [{} for _ in peaks],
            "I_model":     I_model,
            "I_bg":        I_bg,
            "residuals":   residuals,
            "chi2":        float(np.sum((residuals[np.isfinite(residuals)]) ** 2)),
            "dof":         max(1, np.sum(np.isfinite(self._I)) - len(peaks) * 3),
            "reduced_chi2": float(np.nanmean(residuals ** 2)),
        }
        try:
            from pyirena.io.nxcansas_waxs_peakfit import save_waxs_peakfit_results
            q_fit_min, q_fit_max = self._graph.get_q_range()
            save_waxs_peakfit_results(
                self._filepath, result, self._q,
                intensity_data=self._I,
                intensity_error=self._dI,
                q_min=q_fit_min,
                q_max=q_fit_max,
            )
            self._set_status(f"Saved to {self._filepath.name}.")
        except Exception as exc:
            self._set_status(f"Save error: {exc}", error=True)

    def _results_to_graphs(self):
        """Add annotation with key fit results to the graph."""
        bg_shape  = self._bg_combo.currentText()
        bg_params = self._get_bg_params()
        peaks     = self._get_peaks()

        lines = [f"Background: {bg_shape}"]
        for name, pd in bg_params.items():
            lines.append(f"  {name} = {pd['value']:.4g}")
        lines.append(f"Peaks: {len(peaks)}")
        for i, pk in enumerate(peaks):
            q0   = pk.get("Q0",   {}).get("value", 0.0)
            fwhm = pk.get("FWHM", {}).get("value", 0.0)
            A    = pk.get("A",    {}).get("value", 0.0)
            lines.append(f"  P{i+1}: Q0={q0:.4g}  FWHM={fwhm:.4g}  A={A:.4g}")

        text = "\n".join(lines)
        vr = self._graph.main_plot.viewRange()
        dx = vr[0][1] - vr[0][0]
        dy = vr[1][1] - vr[1][0]
        # anchor=(1, 0): top-right corner of the TextItem at the set position,
        # so text is right-aligned and placed below the legend (top-right area)
        ann = pg.TextItem(text=text, color=(40, 40, 40), anchor=(1, 0))
        ann.setPos(vr[0][1] - 0.02 * dx, vr[1][1] - 0.12 * dy)
        self._graph.main_plot.addItem(ann)

    def _export_params(self):
        # Default to the folder where the current data file lives
        default_dir = str(self._filepath.parent) if self._filepath else str(Path.cwd())
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Parameters",
            str(Path(default_dir) / "pyirena_config.json"),
            "JSON (*.json);;All Files (*)",
        )
        if not path:
            return
        try:
            state = self._get_current_state()
            # Load existing config and merge so other tool sections are preserved
            config_path = Path(path)
            config: Dict = {}
            if config_path.exists():
                try:
                    config = json.loads(config_path.read_text(encoding='utf-8'))
                except Exception:
                    pass
            config['waxs_peakfit'] = state
            config_path.write_text(json.dumps(config, indent=2), encoding='utf-8')
            self._set_status(f"Exported to {config_path.name}.")
        except Exception as exc:
            self._set_status(f"Export error: {exc}", error=True)

    def _import_params(self):
        default_dir = str(self._filepath.parent) if self._filepath else str(Path.cwd())
        path, _ = QFileDialog.getOpenFileName(
            self, "Import Parameters", default_dir,
            "JSON (*.json);;All Files (*)",
        )
        if not path:
            return
        try:
            raw = json.loads(Path(path).read_text(encoding='utf-8'))
            # Support both wrapped {"waxs_peakfit": {...}} and unwrapped formats
            state = raw.get('waxs_peakfit', raw)
            self._apply_state(state)
            self._set_status(f"Imported from {Path(path).name}.")
        except Exception as exc:
            self._set_status(f"Import error: {exc}", error=True)

    # ===========================================================================
    # State helpers
    # ===========================================================================

    def _get_current_state(self) -> Dict:
        qmin, qmax = self._graph.get_q_range()
        return {
            "q_min":     qmin,
            "q_max":     qmax,
            "bg_shape":  self._bg_combo.currentText(),
            "bg_params": self._get_bg_params(),
            "peaks":     self._get_peaks(),
            "no_limits": self._no_limits_chk.isChecked(),
            "peak_find": {
                "prominence":    self._pf_prom.value(),
                "min_fwhm":      self._pf_minfwhm.value(),
                "max_fwhm":      self._pf_maxfwhm.value(),
                "min_distance":  self._pf_mindist.value(),
                "sg_window_frac": self._pf_sgwin.value() / 100.0,
            },
        }

    def _apply_state(self, state: Dict):
        if not state:
            return

        # No limits
        self._no_limits_chk.setChecked(bool(state.get("no_limits", False)))

        # Background
        bg_shape = state.get("bg_shape", "Constant")
        self._bg_combo.blockSignals(True)
        self._bg_combo.setCurrentText(bg_shape)
        self._bg_combo.blockSignals(False)
        self._rebuild_bg_grid(bg_shape, state.get("bg_params", {}))

        # Peaks
        self._clear_peaks()
        for pk in state.get("peaks", []):
            self._add_peak_row(pk)

        # Peak-find params
        pf = state.get("peak_find", {})
        if "prominence"    in pf: self._pf_prom.setValue(pf["prominence"])
        if "min_fwhm"      in pf: self._pf_minfwhm.setValue(pf["min_fwhm"])
        if "max_fwhm"      in pf: self._pf_maxfwhm.setValue(pf["max_fwhm"])
        if "min_distance"  in pf: self._pf_mindist.setValue(pf["min_distance"])
        if "sg_window_frac" in pf: self._pf_sgwin.setValue(pf["sg_window_frac"] * 100.0)

        # Q range — restore cursor positions when data is already loaded
        q_min = state.get("q_min")
        q_max = state.get("q_max")
        if q_min is not None and q_max is not None and self._q is not None:
            self._graph.set_q_range(float(q_min), float(q_max))

    # ===========================================================================
    # Public data API
    # ===========================================================================

    def set_data(
        self,
        q: np.ndarray,
        I: np.ndarray,
        dI: Optional[np.ndarray] = None,
        label: str = "Data",
        filepath: Optional[Path] = None,
        is_nxcansas: bool = False,
    ):
        """Load data into the panel and display it."""
        self._q  = np.asarray(q,  float)
        self._I  = np.asarray(I,  float)
        self._dI = np.asarray(dI, float) if dI is not None else None
        self._filepath   = Path(filepath) if filepath else None
        self._is_nxcansas = is_nxcansas

        self._graph.plot_data(self._q, self._I, self._dI, label=label)

        # Load stored results if the file has them
        if is_nxcansas and filepath is not None:
            self._try_load_stored_results(Path(filepath))

    def _try_load_stored_results(self, filepath: Path):
        try:
            from pyirena.io.nxcansas_waxs_peakfit import load_waxs_peakfit_results
            res = load_waxs_peakfit_results(filepath)
            state = {
                "bg_shape":  res.get("bg_shape", "Constant"),
                "bg_params": res.get("bg_params", {}),
                "peaks":     res.get("peaks", []),
                "no_limits": False,
                "q_min":     res.get("q_min"),
                "q_max":     res.get("q_max"),
            }
            self._apply_state(state)
            self._graph_model()
            self._set_status(f"Loaded stored results from {filepath.name}.")
        except KeyError:
            pass   # no stored results — silently ignore
        except Exception as exc:
            self._set_status(f"Could not load stored results: {exc}", error=True)

    # ===========================================================================
    # Status
    # ===========================================================================

    def _set_status(self, msg: str, error: bool = False,
                    progress: bool = False, success: bool = False):
        if error:
            color = "#c0392b"   # red
        elif progress:
            color = "#e67e22"   # orange
        elif success:
            color = "#27ae60"   # green
        else:
            color = "#555"
        self._status.setStyleSheet(
            f"font-size:10pt;color:{color};border-top:1px solid #ccc;padding-top:2px;"
        )
        self._status.setText(msg)
