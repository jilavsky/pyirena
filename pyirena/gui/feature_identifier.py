"""Feature Identifier dialog for the Unified Fit panel.

Opens as a non-modal window from `UnifiedFitPanel`.  Segments the loaded
I(Q) curve into power-law-slope regions and draws them as overlays on the
parent's graph window.  Visualisation only — never modifies the Unified
Fit model or level parameters.
"""
from __future__ import annotations
import logging

log = logging.getLogger(__name__)


from typing import TYPE_CHECKING

import numpy as np

from pyirena.gui._qt import (
    QCheckBox, QDoubleSpinBox, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QPushButton, QTextEdit, QVBoxLayout, QWidget, Qt,
)

import pyqtgraph as pg

from pyirena.core.feature_detect import (
    FeatureDetectConfig,
    detect_features,
)

if TYPE_CHECKING:
    from pyirena.gui.unified_fit import UnifiedFitPanel


# Colours for the four marker kinds
_COLORS = {
    "background":      ((128, 128, 128, 60),  (90, 90, 90, 200)),    # grey
    "guinier_plateau": ((39, 174, 96, 60),    (39, 174, 96, 200)),   # green
    "power_law":       ((243, 156, 18, 60),   (211, 84, 0, 200)),    # orange
}
_KNEE_BRUSH = (231, 76, 60, 50)
_KNEE_PEN   = (192, 57, 43, 220)


# Advanced spin-box fields with sensible ranges.
# v0.8.5: change-point detection has two passes (coarse + tight refinement
# of wide regions); old stability_window / stability_std_max are gone.
_FIELDS = [
    ("change_window_1",          "Pass-1 window (dec)",        0.10, 1.00, 0.05),
    ("change_threshold_1",       "Pass-1 threshold",           0.10, 1.50, 0.05),
    ("change_window_2",          "Pass-2 window (dec)",        0.10, 1.00, 0.05),
    ("change_threshold_2",       "Pass-2 threshold",           0.05, 1.00, 0.05),
    ("wide_region_decades",      "Wide-region threshold (dec)",0.30, 3.00, 0.10),
    ("min_segment_decades",      "Min segment width (dec)",    0.05, 1.00, 0.05),
    ("edge_min_segment_decades", "Edge min width (dec)",       0.02, 0.50, 0.01),
    ("merge_slope_tol",          "Merge slope tolerance",      0.05, 1.50, 0.05),
    ("guinier_knee_min_delta_slope", "Min knee Δslope",        0.05, 2.00, 0.05),
    ("q_max_clip",               "Q max clip (Å⁻¹, 0=off)",   0.0,  1.0,  0.05),
]


class FeatureIdentifierDialog(QWidget):
    """Non-modal window that segments the parent panel's I(Q) into
    power-law-slope regions.  Never modifies the panel's model.
    """

    def __init__(self, parent: "UnifiedFitPanel"):
        super().__init__(None)  # top-level, no Qt parent
        self._panel = parent
        self.setWindowTitle(self._window_title())
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.Window)
        self.resize(440, 600)

        self._markers: list = []  # items added to parent's graph

        # Advanced parameter grid
        self._advanced_box = QGroupBox(
            "Segmentation parameters (log-Q decades unless noted)"
        )
        adv_layout = QGridLayout()
        adv_layout.setColumnStretch(1, 1)
        self._spin_boxes: dict[str, QDoubleSpinBox] = {}
        default_cfg = FeatureDetectConfig()
        for row, (key, label, lo, hi, step) in enumerate(_FIELDS):
            adv_layout.addWidget(QLabel(label), row, 0)
            sb = QDoubleSpinBox()
            sb.setRange(lo, hi)
            sb.setSingleStep(step)
            sb.setDecimals(3)
            sb.setMinimumWidth(90)
            if key == "q_max_clip":
                default = default_cfg.q_max_clip or 0.0
            else:
                default = float(getattr(default_cfg, key))
            sb.setValue(default)
            adv_layout.addWidget(sb, row, 1)
            self._spin_boxes[key] = sb
        self._advanced_box.setLayout(adv_layout)
        self._advanced_box.setVisible(False)

        # Show-advanced toggle
        self.show_advanced = QCheckBox("Show advanced params")
        self.show_advanced.setChecked(False)
        self.show_advanced.stateChanged.connect(
            lambda: self._advanced_box.setVisible(self.show_advanced.isChecked())
        )

        # Restore saved parameters (must happen after spin boxes are created)
        self._restore_params()

        # Detect / Clear buttons
        btn_row = QHBoxLayout()
        self.detect_btn = QPushButton("Detect segments")
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
        self.summary.setMinimumHeight(220)
        self.summary.setStyleSheet(
            "QTextEdit{background:#fbfbfb;border:1px solid #ddd;"
            "font-family:Consolas,monospace;font-size:11px;}"
        )

        # Help label
        help_label = QLabel(self._help_text())
        help_label.setWordWrap(True)
        help_label.setStyleSheet("QLabel{color:#555;font-size:10px;}")

        layout = QVBoxLayout()
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)
        layout.addWidget(self.show_advanced)
        layout.addWidget(self._advanced_box)
        layout.addLayout(btn_row)
        layout.addWidget(QLabel("<b>Result summary</b>"))
        layout.addWidget(self.summary, stretch=1)
        layout.addWidget(help_label)
        self.setLayout(layout)

    # ----------------------------------------------------------------------
    # Overridable presentation hooks (subclasses customise title / help / summary)
    # ----------------------------------------------------------------------

    def _window_title(self) -> str:
        return "Feature Identifier — Power-law segmentation"

    def _help_text(self) -> str:
        return (
            "<b>Visualisation only</b> — segmentation never modifies your fit "
            "parameters.  Segments listed from high-Q to low-Q (Unified Fit "
            "Level 1 = high-Q end).  Use markers as visual guides for placing "
            "cursors and choosing the number of levels."
        )

    def _extra_summary_lines(self, result) -> list:
        """Extra summary lines appended after the standard segmentation report.

        Default: none.  Subclasses (e.g. the Sizes dialog) override this to add
        model-specific recommendations.
        """
        return []

    # ----------------------------------------------------------------------
    # State persistence (via parent panel's StateManager)
    # ----------------------------------------------------------------------

    # Bumped when the algorithm or field set changes incompatibly so that
    # old saved values (which mean different things under different defaults)
    # are not re-applied to a different algorithm.
    _STATE_SCHEMA_VERSION = 2

    def _save_params(self) -> None:
        """Persist current spin-box values + advanced-toggle state."""
        sm = getattr(self._panel, "state_manager", None)
        if sm is None:
            return
        params = {key: float(sb.value()) for key, sb in self._spin_boxes.items()}
        params["_show_advanced"] = bool(self.show_advanced.isChecked())
        params["_schema_version"] = self._STATE_SCHEMA_VERSION
        sm.update("feature_detect", params)
        sm.save()

    def _restore_params(self) -> None:
        """Restore spin-box values saved by a previous session."""
        sm = getattr(self._panel, "state_manager", None)
        if sm is None:
            return
        saved = sm.get("feature_detect")
        if not isinstance(saved, dict):
            return
        # Reject state saved under an older schema (v0.8.4 used different
        # field names; some field names overlap, so just-key-matching would
        # produce a half-old-half-new config).
        if saved.get("_schema_version") != self._STATE_SCHEMA_VERSION:
            return
        for key, sb in self._spin_boxes.items():
            if key in saved:
                try:
                    sb.setValue(float(saved[key]))
                except (TypeError, ValueError):
                    log.debug("suppressed exception", exc_info=True)
        if "_show_advanced" in saved:
            self.show_advanced.setChecked(bool(saved["_show_advanced"]))
            self._advanced_box.setVisible(bool(saved["_show_advanced"]))

    # ----------------------------------------------------------------------
    # Config + data access
    # ----------------------------------------------------------------------

    def _current_cfg(self) -> FeatureDetectConfig:
        cfg = FeatureDetectConfig()
        for key, sb in self._spin_boxes.items():
            val = float(sb.value())
            if key == "q_max_clip":
                cfg.q_max_clip = None if val <= 0.0 else val
            else:
                setattr(cfg, key, val)
        return cfg

    def _has_data(self) -> bool:
        data = getattr(self._panel, "data", None)
        return bool(data and "Q" in data and "Intensity" in data
                    and data["Q"] is not None and data["Intensity"] is not None)

    # ----------------------------------------------------------------------
    # Detection
    # ----------------------------------------------------------------------

    def _on_detect(self):
        if not self._has_data():
            self.summary.setPlainText(
                "No data loaded.  Load a dataset in the Unified Fit panel first."
            )
            return
        q = np.asarray(self._panel.data["Q"], dtype=float)
        I = np.asarray(self._panel.data["Intensity"], dtype=float)
        err = self._panel.data.get("Error")
        if err is not None:
            err = np.asarray(err, dtype=float)
        cfg = self._current_cfg()
        # Slit-smeared data → detection is approximate (E3); pass the slit
        # length so the result carries an advisory shown in the summary.
        slit_length = 0.0
        if bool(self._panel.data.get("is_slit_smeared", False)):
            slit_length = float(self._panel.data.get("slit_length", 0.0) or 0.0)
        try:
            result = detect_features(q, I, sigma_I=err, config=cfg,
                                     slit_length=slit_length)
        except Exception as exc:  # pragma: no cover
            self.summary.setPlainText(f"Detection failed: {exc}")
            return

        self._save_params()
        self._clear_markers()
        self._render_markers(result)
        self._render_summary(result)

    # ----------------------------------------------------------------------
    # Rendering on parent graph
    # ----------------------------------------------------------------------

    def _plot(self):
        gw = getattr(self._panel, "graph_window", None)
        if gw is None:
            return None
        return getattr(gw, "main_plot", None)

    def _render_markers(self, result):
        plot = self._plot()
        if plot is None:
            return
        # Plot uses setLogMode(x=True, y=True) → x positions are log10(q)

        # Segments: translucent vertical regions + label
        for i, seg in enumerate(result.segments, start=1):
            brush, pen = _COLORS.get(seg["kind"], _COLORS["power_law"])
            q_lo, q_hi = seg["q_min"], seg["q_max"]
            region = pg.LinearRegionItem(
                values=(float(np.log10(q_lo)), float(np.log10(q_hi))),
                orientation="vertical",
                brush=brush,
                pen=pg.mkPen(pen, width=1),
                movable=False,
            )
            region.setZValue(-10)
            plot.addItem(region)
            self._markers.append(region)
            # Label
            tag = self._label_for_kind(seg["kind"])
            text = f"{tag} P={seg['P']:.1f}" if seg["kind"] == "power_law" else tag
            label = pg.TextItem(text, color=pen, anchor=(0.5, 1.0))
            log_q_mid = 0.5 * (np.log10(q_lo) + np.log10(q_hi))
            label.setPos(log_q_mid, self._top_y(plot))
            plot.addItem(label)
            self._markers.append(label)

        # Guinier knees: red translucent band between segments
        for i, knee in enumerate(result.guinier_knees, start=1):
            q_lo, q_hi = knee["q_min"], knee["q_max"]
            # Ensure visible width even when the gap is small (less than ~0.01 dec)
            log_lo = float(np.log10(q_lo))
            log_hi = float(np.log10(q_hi))
            if log_hi - log_lo < 0.04:
                pad = (0.04 - (log_hi - log_lo)) / 2
                log_lo -= pad
                log_hi += pad
            region = pg.LinearRegionItem(
                values=(log_lo, log_hi),
                orientation="vertical",
                brush=_KNEE_BRUSH,
                pen=pg.mkPen(_KNEE_PEN, width=1, style=Qt.PenStyle.DashLine),
                movable=False,
            )
            region.setZValue(-8)
            plot.addItem(region)
            self._markers.append(region)
            label = pg.TextItem(
                f"GK ΔP={knee['delta_P']:.1f}",
                color=_KNEE_PEN, anchor=(0.5, 0.0)
            )
            label.setPos(0.5 * (log_lo + log_hi), self._bottom_y(plot))
            plot.addItem(label)
            self._markers.append(label)

    @staticmethod
    def _label_for_kind(kind: str) -> str:
        return {
            "background":      "Background",
            "guinier_plateau": "GP",
            "power_law":       "PLS",
        }.get(kind, kind)

    def _top_y(self, plot) -> float:
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
                log.debug("suppressed exception", exc_info=True)
        self._markers.clear()

    # ----------------------------------------------------------------------
    # Summary text
    # ----------------------------------------------------------------------

    def _render_summary(self, result):
        lines: list[str] = []
        for w in getattr(result, "warnings", []) or []:
            lines.append(f"⚠ {w}")
        if getattr(result, "warnings", None):
            lines.append("")
        lines.append(f"Q analysed: [{result.q_min_analysed:.4g}, "
                     f"{result.q_max_analysed:.4g}]   "
                     f"({result.log_decades:.2f} decades, "
                     f"{result.n_points} points)")
        lines.append("")
        n_pls = sum(1 for s in result.segments if s["kind"] == "power_law")
        n_gp  = sum(1 for s in result.segments if s["kind"] == "guinier_plateau")
        n_bg  = sum(1 for s in result.segments if s["kind"] == "background")
        lines.append(f"Found {len(result.segments)} segment(s): "
                     f"{n_pls} power-law, {n_gp} guinier-plateau, "
                     f"{n_bg} background; {len(result.guinier_knees)} knee(s).")
        lines.append(f"Suggested Unified Fit levels: {result.recommended_nlevels}")
        lines.append("")
        lines.append("Segments (high-Q → low-Q; Unified Fit level 1 = high-Q end):")
        # Display segments in reverse (high-Q first) per user's preferred order
        for idx, seg in enumerate(reversed(result.segments), start=1):
            tag = self._label_for_kind(seg["kind"]).ljust(11)
            extra = (f"P = {seg['P']:.2f} ± {seg['P_std']:.2f}"
                     if seg["kind"] == "power_law" else
                     f"P = {seg['P']:.2f}")
            lines.append(
                f"  L{idx}: {tag} Q ∈ [{seg['q_min']:.4g}, {seg['q_max']:.4g}]  "
                f"({seg['width_decades']:.2f} dec)  {extra}"
            )
        if result.guinier_knees:
            lines.append("")
            lines.append("Guinier knees (between segments):")
            for k in result.guinier_knees:
                lines.append(
                    f"  Q ∈ [{k['q_min']:.4g}, {k['q_max']:.4g}]  "
                    f"P: {k['P_high_q']:.2f} (high-Q) → "
                    f"{k['P_low_q']:.2f} (low-Q)  "
                    f"(ΔP = {k['delta_P']:.2f})"
                )
        if result.recommended_guinier_windows:
            lines.append("")
            lines.append("Suggested Guinier fit windows "
                         "(for fit_local_guinier between cursors):")
            for w in result.recommended_guinier_windows:
                lines.append(
                    f"  Q ∈ [{w['q_min_guinier']:.4g}, "
                    f"{w['q_max_guinier']:.4g}]   "
                    f"power-law starts at Q ≈ {w['q_min_powerlaw']:.4g}"
                )
        if result.background_q_min is not None:
            lines.append("")
            lines.append(f"Background begins at Q ≈ {result.background_q_min:.4g}")
        lines.extend(self._extra_summary_lines(result))
        self.summary.setPlainText("\n".join(lines))

    # ----------------------------------------------------------------------
    # Window lifecycle
    # ----------------------------------------------------------------------

    def closeEvent(self, event):
        self._save_params()
        self._clear_markers()
        try:
            self._panel._feature_dialog = None  # noqa: SLF001
        except Exception:
            log.debug("suppressed exception", exc_info=True)
        super().closeEvent(event)


__all__ = ["FeatureIdentifierDialog"]
