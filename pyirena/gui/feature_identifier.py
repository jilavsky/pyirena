"""Feature Identifier dialog for the Unified Fit panel.

Opens as a non-modal window from `UnifiedFitPanel`.  Detects plateaus,
power-law regions, and structure-factor peaks in the loaded I(Q) data and
draws them as overlays on the parent's graph window.  Visualisation-only —
never modifies the Unified Fit model or level parameters.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np

try:
    from PySide6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
        QLabel, QPushButton, QComboBox, QDoubleSpinBox, QCheckBox,
        QGroupBox, QTextEdit, QSizePolicy,
    )
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QColor
except ImportError:
    try:
        from PyQt6.QtWidgets import (  # type: ignore[no-redef]
            QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
            QLabel, QPushButton, QComboBox, QDoubleSpinBox, QCheckBox,
            QGroupBox, QTextEdit, QSizePolicy,
        )
        from PyQt6.QtCore import Qt  # type: ignore[no-redef]
        from PyQt6.QtGui import QColor  # type: ignore[no-redef]
    except ImportError:
        from PyQt5.QtWidgets import (  # type: ignore[no-redef]
            QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
            QLabel, QPushButton, QComboBox, QDoubleSpinBox, QCheckBox,
            QGroupBox, QTextEdit, QSizePolicy,
        )
        from PyQt5.QtCore import Qt  # type: ignore[no-redef]
        from PyQt5.QtGui import QColor  # type: ignore[no-redef]

import pyqtgraph as pg

from pyirena.core.feature_detect import (
    FeatureDetectConfig,
    detect_features,
)

if TYPE_CHECKING:
    from pyirena.gui.unified_fit import UnifiedFitPanel


# Colors (matching the rest of the pyirena GUI palette where possible)
_PLATEAU_BRUSH = (39, 174, 96, 60)   # translucent green
_PLATEAU_PEN   = (39, 174, 96, 200)
_PORODBR       = (243, 156, 18, 50)  # translucent orange
_PORODPEN      = (211, 84, 0, 200)
_PEAK_PEN      = (231, 76, 60, 220)  # red
_TEXT_COLOR    = (40, 40, 40)


# Fields that the user can edit in custom mode, with sensible ranges
_FIELDS = [
    ("sigma_smooth",       "σ smooth (dec)",     0.05, 0.50, 0.01),
    ("span_deriv",         "slope span (dec)",   0.10, 1.00, 0.05),
    ("plateau_slope_max",  "plateau |slope| max", 0.05, 2.0,  0.05),
    ("plateau_min_width",  "plateau min width",   0.05, 1.0,  0.05),
    ("power_law_max_slope","power-law slope max",-5.0, -0.1, 0.1),
    ("power_law_min_width","power-law min width", 0.05, 2.0,  0.05),
    ("peak_min_slope_mag", "peak slope mag",      0.05, 2.0,  0.05),
]


class FeatureIdentifierDialog(QWidget):
    """Non-modal window that detects and visualises features on the parent
    `UnifiedFitPanel`'s graph.  Never modifies the panel's model.
    """

    def __init__(self, parent: "UnifiedFitPanel"):
        super().__init__(None)  # top-level, no Qt parent → independent window
        self._panel = parent
        self.setWindowTitle("Feature Identifier")
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.Window)
        self.resize(420, 540)

        # Items we add to the parent's graph; cleared on Clear button or close
        self._markers: list = []

        # Preset combo
        preset_row = QHBoxLayout()
        preset_row.addWidget(QLabel("Preset:"))
        self.preset_combo = QComboBox()
        self.preset_combo.addItems(["Auto", "SAXS", "USAXS", "Custom"])
        self.preset_combo.setCurrentText("Auto")
        self.preset_combo.currentTextChanged.connect(self._on_preset_changed)
        preset_row.addWidget(self.preset_combo)
        preset_row.addStretch()

        # Advanced toggle
        self.show_advanced = QCheckBox("Show advanced params")
        self.show_advanced.setChecked(False)
        self.show_advanced.stateChanged.connect(self._update_advanced_visibility)

        # Advanced parameter grid
        self._advanced_box = QGroupBox("Detection parameters (log-Q decades)")
        adv_layout = QGridLayout()
        adv_layout.setColumnStretch(1, 1)
        self._spin_boxes: dict[str, QDoubleSpinBox] = {}
        for row, (key, label, lo, hi, step) in enumerate(_FIELDS):
            adv_layout.addWidget(QLabel(label), row, 0)
            sb = QDoubleSpinBox()
            sb.setRange(lo, hi)
            sb.setSingleStep(step)
            sb.setDecimals(3)
            sb.setMinimumWidth(90)
            adv_layout.addWidget(sb, row, 1)
            self._spin_boxes[key] = sb
        self._advanced_box.setLayout(adv_layout)
        self._advanced_box.setVisible(False)

        # Detect / Clear buttons
        btn_row = QHBoxLayout()
        self.detect_btn = QPushButton("Detect features")
        self.detect_btn.setStyleSheet(
            "QPushButton{background:#2980b9;color:white;font-weight:bold;padding:6px 10px;}"
            "QPushButton:hover{background:#3498db;}"
        )
        self.detect_btn.setMinimumHeight(28)
        self.detect_btn.clicked.connect(self._on_detect)
        btn_row.addWidget(self.detect_btn)

        self.clear_btn = QPushButton("Clear markers")
        self.clear_btn.setStyleSheet(
            "QPushButton{background:#95a5a6;color:white;padding:6px 10px;}"
            "QPushButton:hover{background:#7f8c8d;}"
        )
        self.clear_btn.setMinimumHeight(28)
        self.clear_btn.clicked.connect(self._clear_markers)
        btn_row.addWidget(self.clear_btn)
        btn_row.addStretch()

        # Summary
        self.summary = QTextEdit()
        self.summary.setReadOnly(True)
        self.summary.setMinimumHeight(180)
        self.summary.setStyleSheet(
            "QTextEdit{background:#fbfbfb;border:1px solid #ddd;font-family:Consolas,monospace;font-size:11px;}"
        )

        # Help text
        help_label = QLabel(
            "<b>Visualisation only</b> — detection never modifies your fit "
            "parameters.  Use the markers as visual guides for placing the "
            "Unified Fit Q cursors and choosing the number of levels."
        )
        help_label.setWordWrap(True)
        help_label.setStyleSheet("QLabel{color:#555;font-size:10px;}")

        # Assemble
        layout = QVBoxLayout()
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)
        layout.addLayout(preset_row)
        layout.addWidget(self.show_advanced)
        layout.addWidget(self._advanced_box)
        layout.addLayout(btn_row)
        layout.addWidget(QLabel("<b>Result summary</b>"))
        layout.addWidget(self.summary, stretch=1)
        layout.addWidget(help_label)
        self.setLayout(layout)

        # Initial population of spin boxes from auto preset (uses current data
        # range if available, else SAXS preset)
        self._on_preset_changed(self.preset_combo.currentText())

    # ----------------------------------------------------------------------
    # Preset handling
    # ----------------------------------------------------------------------

    def _current_cfg(self) -> FeatureDetectConfig:
        """Build a config from the preset combo + current spin-box values."""
        text = self.preset_combo.currentText().lower()
        if text == "auto" and self._has_data():
            q = np.asarray(self._panel.data["Q"], dtype=float)
            cfg = FeatureDetectConfig.auto(q)
        elif text == "saxs":
            cfg = FeatureDetectConfig.saxs_preset()
        elif text == "usaxs":
            cfg = FeatureDetectConfig.usaxs_preset()
        elif text == "custom":
            # Start from defaults, override from spin boxes
            cfg = FeatureDetectConfig()
            for key, sb in self._spin_boxes.items():
                if hasattr(cfg, key):
                    setattr(cfg, key, float(sb.value()))
        else:
            cfg = FeatureDetectConfig.saxs_preset()
        return cfg

    def _on_preset_changed(self, text: str):
        """Repopulate spin-box values from the new preset (unless Custom)."""
        if text.lower() == "custom":
            # Leave current spin values; user edits them
            return
        cfg = self._current_cfg()
        for key, sb in self._spin_boxes.items():
            if hasattr(cfg, key):
                sb.blockSignals(True)
                sb.setValue(float(getattr(cfg, key)))
                sb.blockSignals(False)

    def _update_advanced_visibility(self):
        self._advanced_box.setVisible(self.show_advanced.isChecked())

    # ----------------------------------------------------------------------
    # Data access
    # ----------------------------------------------------------------------

    def _has_data(self) -> bool:
        data = getattr(self._panel, "data", None)
        return bool(data and "Q" in data and "Intensity" in data
                    and data["Q"] is not None and data["Intensity"] is not None)

    # ----------------------------------------------------------------------
    # Detection
    # ----------------------------------------------------------------------

    def _on_detect(self):
        if not self._has_data():
            self.summary.setPlainText("No data loaded.  Load a dataset in "
                                       "the Unified Fit panel first.")
            return
        q = np.asarray(self._panel.data["Q"], dtype=float)
        I = np.asarray(self._panel.data["Intensity"], dtype=float)
        err = self._panel.data.get("Error")
        if err is not None:
            err = np.asarray(err, dtype=float)
        cfg = self._current_cfg()
        # If a preset was selected (not custom), still let the user edit spin
        # boxes to override individual fields by re-reading them
        if self.preset_combo.currentText().lower() != "auto":
            for key, sb in self._spin_boxes.items():
                if hasattr(cfg, key):
                    setattr(cfg, key, float(sb.value()))

        try:
            result = detect_features(q, I, sigma_I=err, config=cfg)
        except Exception as exc:  # pragma: no cover — guard against unexpected
            self.summary.setPlainText(f"Detection failed: {exc}")
            return

        self._clear_markers()
        self._render_markers(result)
        self._render_summary(result)

    # ----------------------------------------------------------------------
    # Rendering
    # ----------------------------------------------------------------------

    def _plot(self):
        """Return the parent's main pyqtgraph PlotItem (log-log I vs Q)."""
        gw = getattr(self._panel, "graph_window", None)
        if gw is None:
            return None
        return getattr(gw, "main_plot", None)

    def _render_markers(self, result):
        plot = self._plot()
        if plot is None:
            return
        # ViewBox is in log10 space (setLogMode(x=True, y=True)).  Region/line
        # positions must be passed as log10(q).

        # Plateaus — green semi-transparent vertical regions
        for i, p in enumerate(result.plateaus, start=1):
            q_lo, q_hi = p["q_min"], p["q_max"]
            region = pg.LinearRegionItem(
                values=(float(np.log10(q_lo)), float(np.log10(q_hi))),
                orientation="vertical",
                brush=_PLATEAU_BRUSH,
                pen=pg.mkPen(_PLATEAU_PEN, width=1),
                movable=False,
            )
            region.setZValue(-10)  # behind the data
            plot.addItem(region)
            self._markers.append(region)
            label = pg.TextItem(f"Plateau {i}", color=_PORODPEN if False else _PLATEAU_PEN,
                                anchor=(0.5, 1.0))
            label.setPos(float(np.log10(p["q_center"])),
                         self._top_y(plot))
            plot.addItem(label)
            self._markers.append(label)

        # Power-law regions — orange semi-transparent regions
        for i, pl in enumerate(result.power_law_regions, start=1):
            q_lo, q_hi = pl["q_min"], pl["q_max"]
            region = pg.LinearRegionItem(
                values=(float(np.log10(q_lo)), float(np.log10(q_hi))),
                orientation="vertical",
                brush=_PORODBR,
                pen=pg.mkPen(_PORODPEN, width=1, style=Qt.PenStyle.DashLine),
                movable=False,
            )
            region.setZValue(-9)
            plot.addItem(region)
            self._markers.append(region)
            mid = float(np.sqrt(q_lo * q_hi))
            label = pg.TextItem(f"P≈{abs(pl['slope']):.1f}",
                                color=_PORODPEN, anchor=(0.5, 0.0))
            label.setPos(float(np.log10(mid)), self._bottom_y(plot))
            plot.addItem(label)
            self._markers.append(label)

        # Peaks — red vertical line
        for i, pk in enumerate(result.peaks, start=1):
            line = pg.InfiniteLine(
                pos=float(np.log10(pk["q_peak"])),
                angle=90,
                pen=pg.mkPen(_PEAK_PEN, width=2),
                movable=False,
            )
            plot.addItem(line)
            self._markers.append(line)
            label = pg.TextItem(f"Peak Q={pk['q_peak']:.3g}",
                                color=_PEAK_PEN, anchor=(0.0, 0.0))
            label.setPos(float(np.log10(pk["q_peak"])), self._top_y(plot))
            plot.addItem(label)
            self._markers.append(label)

    def _top_y(self, plot) -> float:
        """Get a y position near the top of the visible view (log10 space)."""
        try:
            (_, _), (y_lo, y_hi) = plot.viewRange()
            return y_hi - 0.05 * (y_hi - y_lo)
        except Exception:
            return 0.0

    def _bottom_y(self, plot) -> float:
        try:
            (_, _), (y_lo, y_hi) = plot.viewRange()
            return y_lo + 0.05 * (y_hi - y_lo)
        except Exception:
            return 0.0

    def _clear_markers(self):
        plot = self._plot()
        if plot is None:
            self._markers.clear()
            return
        for item in self._markers:
            try:
                plot.removeItem(item)
            except Exception:
                pass
        self._markers.clear()

    # ----------------------------------------------------------------------
    # Summary
    # ----------------------------------------------------------------------

    def _render_summary(self, result):
        lines: list[str] = []
        lines.append(f"Data: {result.n_points} points,  "
                     f"{result.log_decades:.2f} log decades")
        lines.append(f"Preset used: {result.preset_used}")
        lines.append("")
        lines.append(f"Found {len(result.plateaus)} plateau(s), "
                     f"{len(result.peaks)} peak(s), "
                     f"{len(result.power_law_regions)} power-law region(s).")
        lines.append(f"→ Suggested Unified Fit levels: "
                     f"{result.recommended_nlevels}")
        lines.append("")
        if result.plateaus:
            lines.append("Plateaus (Guinier knees):")
            for i, p in enumerate(result.plateaus, 1):
                lines.append(f"  P{i}: Q ∈ [{p['q_min']:.4g}, {p['q_max']:.4g}]  "
                             f"width = {p['width_decades']:.2f} dec")
        if result.power_law_regions:
            lines.append("Power-law regions:")
            for i, pl in enumerate(result.power_law_regions, 1):
                lines.append(f"  PL{i}: Q ∈ [{pl['q_min']:.4g}, {pl['q_max']:.4g}]  "
                             f"slope = {pl['slope']:.2f} ± {pl['slope_std']:.2f}")
        if result.peaks:
            lines.append("Peaks (Guinier + structure factor):")
            for i, pk in enumerate(result.peaks, 1):
                lines.append(f"  Pk{i}: Q_peak = {pk['q_peak']:.4g}  "
                             f"(Q ∈ [{pk['q_low']:.4g}, {pk['q_high']:.4g}])")
        if result.recommended_guinier_windows:
            lines.append("")
            lines.append("Recommended Guinier windows:")
            for i, w in enumerate(result.recommended_guinier_windows, 1):
                lines.append(f"  L{i} ({w['feature_type']}): "
                             f"Guinier Q ∈ [{w['q_min_guinier']:.4g}, "
                             f"{w['q_max_guinier']:.4g}]   "
                             f"Power-law starts at Q ≈ {w['q_min_powerlaw']:.4g}")
        if result.background_q_min is not None:
            lines.append("")
            lines.append(f"High-Q flat background begins at "
                         f"Q ≈ {result.background_q_min:.4g}")
        self.summary.setPlainText("\n".join(lines))

    # ----------------------------------------------------------------------
    # Window lifecycle
    # ----------------------------------------------------------------------

    def closeEvent(self, event):
        """Clear markers from the parent's graph and notify the panel so it
        creates a fresh dialog next time."""
        self._clear_markers()
        try:
            self._panel._feature_dialog = None  # noqa: SLF001
        except Exception:
            pass
        super().closeEvent(event)


__all__ = ["FeatureIdentifierDialog"]
