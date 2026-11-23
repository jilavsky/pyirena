"""
Demo: pyIrena button looks, side by side, plus a performance benchmark.

Run from the repo root with the pyirena environment:

    /opt/miniconda3/envs/pyirena/bin/python planning/raised-buttons/demo_raised_buttons.py

Columns
    A  Today        — pyIrena's current theme rule + the inline styles panels use now.
    B  Raised QSS   — the proposal (raised_buttons.py), via btnRole properties.
    C  Fusion only  — no colour stylesheet; colours come from each button's QPalette.

The sliders retune column B live. "Run benchmark" builds 300 buttons per
approach off-screen and times construction+polish and repaint.

This file lives outside the package on purpose (planning/), so it is not
subject to the _qt.py contract test.
"""

from __future__ import annotations

import re
import sys
import time
from pathlib import Path

from pyirena.gui import raised_buttons as rb
from pyirena.gui import theme
from pyirena.gui._qt import (
    QApplication,
    QCheckBox,
    QColor,
    QComboBox,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPalette,
    QPlainTextEdit,
    QPushButton,
    QSlider,
    Qt,
    QVBoxLayout,
    QWidget,
)

try:  # only used by the benchmark's "what not to do" row
    from PySide6.QtWidgets import QGraphicsDropShadowEffect
except ImportError:  # pragma: no cover
    QGraphicsDropShadowEffect = None

# (label, role, size, checkable, enabled) — a realistic slice of a fit panel
BUTTONS = [
    ("Graph Unified", "graph", None, False, True),
    ("Fit", "primary", "large", False, True),
    ("Revert back", "revert", None, False, True),
    ("Calculate uncertainty", "analysis", None, False, True),
    ("Fit power law", "helper", None, False, True),
    ("Report", "tool", None, False, True),
    ("Cancel fit", "stop", None, False, True),
    ("Store in File", "settings", None, False, True),
    ("Load setup", "setup", None, False, True),
    ("Browse…", None, None, False, True),
    ("Fit B/P btwn cursors", None, "small", False, True),
    ("Log X (toggle me)", None, None, True, True),
    ("Fit (disabled)", "primary", None, False, False),
    ("Reset (disabled)", None, None, False, False),
]


# ---------------------------------------------------------------------------
# Today's look, reproduced from pyirena.gui.theme and the panels' inline styles
# ---------------------------------------------------------------------------

_BASE = theme._base_stylesheet()
_m = re.search(r"/\* --- buttons.*?(?=/\* --- selection controls)", _BASE, re.S)
BASE_WITHOUT_BUTTONS = _BASE.replace(_m.group(0), "") if _m else _BASE

# Frozen copy of the pre-2026-11-23 flat look (theme.py has since moved on).
TODAY_BUTTON_BLOCK = """
QPushButton, QToolButton { background-color: #e6eaed; color: #20303c;
    border: 1px solid #b8c0c6; border-radius: 3px; padding: 3px 8px; }
QPushButton:hover, QToolButton:hover { background-color: #d5dbe0; color: #20303c; }
QPushButton:pressed, QToolButton:pressed { background-color: #d6dbdf; color: #20303c; }
QPushButton:checked, QToolButton:checked { background-color: #2f6fb5; color: #ffffff;
    border-color: #2f6fb5; }
QPushButton:disabled, QToolButton:disabled { background-color: #f0f0f0; color: #9aa3ab;
    border-color: #d6dbdf; }
"""
_TODAY_CHIP = ("QPushButton { background-color: #e6eaed; color: #20303c; font-size: 10px;"
               " padding: 1px 6px; border: 1px solid #b8c0c6; border-radius: 3px; }"
               "QPushButton:hover { background-color: #d5dbe0; color: #20303c; }")


def _today_accent(bg, hover, font_size=None):
    size = f"font-size: {font_size};" if font_size else ""
    return (f"QPushButton {{ background-color: {bg}; color: #ffffff; font-weight: bold;"
            f" {size} border: 1px solid {bg}; border-radius: 3px; padding: 3px 8px; }}"
            f"QPushButton:hover {{ background-color: {hover}; border-color: {hover}; }}"
            f"QPushButton:disabled {{ background-color: #95a5a6; color: #f2f4f5; }}")


def _today_soft(bg):
    return (f"QPushButton {{ background-color: {bg}; color: #20303c;"
            f" border: 1px solid #b8c0c6; border-radius: 3px; padding: 3px 8px; }}"
            f"QPushButton:hover {{ border-color: #2f6fb5; }}")


def today_inline(role: str | None, size: str | None) -> str:
    """The inline stylesheet a pyIrena panel attaches today for this role."""
    if role is None:
        return _TODAY_CHIP if size == "small" else ""
    bg, fg, _bold = rb.DEFAULT_ROLES[role]
    if role in ("settings", "setup"):
        return _today_soft(bg)
    font = "12px" if size == "large" else None
    return _today_accent(bg, rb.darken(bg, 0.15), font_size=font)


def make_button(label, role, size, checkable, enabled, mode, parent=None):
    b = QPushButton(label, parent)
    b.setCheckable(checkable)
    b.setEnabled(enabled)
    if mode == "today":
        css = today_inline(role, size)
        if css:
            b.setStyleSheet(css)
    elif mode == "raised-role":
        rb.set_role(b, role, size)
    elif mode == "raised-inline":
        b.setStyleSheet(rb.inline_qss(role or "neutral", size=size or "normal"))
    elif mode == "fusion":
        if role:
            bg, fg, bold = rb.DEFAULT_ROLES[role]
            # Active/Inactive only, so the disabled group keeps the grey.
            # NB: the label colour is overridden by pyIrena's app-wide
            # `QWidget { color: ... }` rule — a real cost of this option.
            pal = b.palette()
            for grp in (QPalette.ColorGroup.Active, QPalette.ColorGroup.Inactive):
                pal.setColor(grp, QPalette.ColorRole.Button, QColor(bg))
                pal.setColor(grp, QPalette.ColorRole.ButtonText, QColor(fg))
            b.setPalette(pal)
            if bold:
                f = b.font()
                f.setBold(True)
                b.setFont(f)
    return b


def column(title: str, mode: str) -> tuple[QGroupBox, QWidget]:
    box = QGroupBox(title)
    lay = QVBoxLayout(box)
    inner = QWidget()
    il = QVBoxLayout(inner)
    il.setContentsMargins(4, 4, 4, 4)
    for spec in BUTTONS:
        il.addWidget(make_button(*spec, mode))
    # Non-button controls: the complaint is that buttons look like these.
    il.addWidget(QLabel("Next to fields and check boxes:"))
    il.addWidget(QLineEdit("0.0123"))
    cb = QComboBox()
    cb.addItems(["Spheroid", "Cylinder"])
    il.addWidget(cb)
    il.addWidget(QCheckBox("Fit?"))
    il.addStretch(1)
    lay.addWidget(inner)
    return box, inner


class Demo(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("pyIrena — button look comparison")

        top = QHBoxLayout()
        a, self.col_a = column("A  Today", "today")
        b, self.col_b = column("B  Raised QSS (proposal)", "raised-role")
        c, self.col_c = column("C  Fusion only, palette tint", "fusion")
        self.col_a.setStyleSheet(TODAY_BUTTON_BLOCK)
        for w in (a, b, c):
            top.addWidget(w)

        # live tuning of column B
        tune = QGridLayout()
        self.bevel = QSlider(Qt.Orientation.Horizontal)
        self.bevel.setRange(0, 45)
        self.bevel.setValue(round(rb.DEFAULT_BEVEL * 100))
        self.radius = QSlider(Qt.Orientation.Horizontal)
        self.radius.setRange(0, 12)
        self.radius.setValue(rb.DEFAULT_RADIUS)
        self.bevel_lbl = QLabel()
        self.radius_lbl = QLabel()
        tune.addWidget(QLabel("Bevel depth (B):"), 0, 0)
        tune.addWidget(self.bevel, 0, 1)
        tune.addWidget(self.bevel_lbl, 0, 2)
        tune.addWidget(QLabel("Corner radius (B):"), 1, 0)
        tune.addWidget(self.radius, 1, 1)
        tune.addWidget(self.radius_lbl, 1, 2)
        self.bevel.valueChanged.connect(self.restyle)
        self.radius.valueChanged.connect(self.restyle)

        self.qss_view = QPlainTextEdit()
        self.qss_view.setReadOnly(True)
        self.qss_view.setMaximumHeight(110)

        bench_row = QHBoxLayout()
        self.bench_btn = QPushButton("Run benchmark (300 buttons per approach)")
        self.bench_btn.setStyleSheet(rb.inline_qss("tool"))
        self.bench_btn.clicked.connect(self.run_benchmark)
        bench_row.addWidget(self.bench_btn)
        bench_row.addStretch(1)
        self.bench_out = QPlainTextEdit()
        self.bench_out.setReadOnly(True)
        self.bench_out.setMaximumHeight(130)
        self.bench_out.setPlaceholderText("Benchmark results appear here.")

        root = QVBoxLayout(self)
        root.addLayout(top)
        root.addLayout(tune)
        root.addWidget(QLabel("Generated QSS for column B (neutral rule shown; ~20 lines per role):"))
        root.addWidget(self.qss_view)
        root.addLayout(bench_row)
        root.addWidget(self.bench_out)
        self.restyle()

    def restyle(self):
        bevel = self.bevel.value() / 100
        radius = self.radius.value()
        self.bevel_lbl.setText(f"{bevel:.2f}")
        self.radius_lbl.setText(f"{radius}px")
        self.col_b.setStyleSheet(rb.stylesheet(bevel=bevel, radius=radius))
        bg, fg, bold = rb.DEFAULT_ROLES["neutral"]
        self.qss_view.setPlainText(
            rb.button_qss(bg, fg, bold=bold, bevel=bevel, radius=radius)
            + rb.size_qss("QPushButton", 3, 10))

    # -- benchmark ---------------------------------------------------------

    def _bench_one(self, mode: str, n: int = 300, shadow: bool = False):
        app = QApplication.instance()
        t0 = time.perf_counter()
        host = QWidget()
        host.setAttribute(Qt.WidgetAttribute.WA_DontShowOnScreen)
        if mode == "today":
            host.setStyleSheet(TODAY_BUTTON_BLOCK)
        elif mode == "raised-role":
            host.setStyleSheet(rb.stylesheet())
        grid = QGridLayout(host)
        for i in range(n):
            spec = BUTTONS[i % len(BUTTONS)]
            btn = make_button(*spec, mode)
            if shadow and QGraphicsDropShadowEffect is not None:
                eff = QGraphicsDropShadowEffect(btn)
                eff.setBlurRadius(6)
                eff.setOffset(1, 2)
                btn.setGraphicsEffect(eff)
            grid.addWidget(btn, i // 15, i % 15)
        host.show()
        app.processEvents()
        build_ms = (time.perf_counter() - t0) * 1e3
        reps = 10
        t1 = time.perf_counter()
        for _ in range(reps):
            host.grab()
        paint_ms = (time.perf_counter() - t1) * 1e3 / reps
        host.hide()
        host.deleteLater()
        app.processEvents()
        return build_ms, paint_ms

    def run_benchmark(self):
        self.bench_btn.setEnabled(False)
        self.bench_out.setPlainText("Running…")
        QApplication.instance().processEvents()
        rows = [
            ("A today (inline per button)", "today", False),
            ("B raised via btnRole (one sheet)", "raised-role", False),
            ("B raised inline per button", "raised-inline", False),
            ("C Fusion palette tint", "fusion", False),
            ("Raised + drop-shadow effect (DON'T)", "raised-role", True),
        ]
        self._bench_one("today", 30)  # warm-up
        lines = [f"{'approach':40s} {'build+polish':>13s} {'repaint all':>12s}"]
        for label, mode, shadow in rows:
            if shadow and QGraphicsDropShadowEffect is None:
                continue
            build, paint = self._bench_one(mode, shadow=shadow)
            lines.append(f"{label:40s} {build:10.1f} ms {paint:9.2f} ms")
        lines.append("(300 buttons; repaint = one full off-screen render, mean of 10)")
        self.bench_out.setPlainText("\n".join(lines))
        self.bench_btn.setEnabled(True)


def main():
    app = QApplication(sys.argv)
    # Same as pyirena.gui.theme.apply_theme, minus the button block, so the
    # three columns can each carry their own button styling.
    app.setStyle("Fusion")
    app.setPalette(theme._palette())
    app.setStyleSheet(BASE_WITHOUT_BUTTONS)
    w = Demo()
    w.resize(900, 900)
    w.show()
    if "--bench" in sys.argv:
        w.run_benchmark()
        print(w.bench_out.toPlainText())
        return 0
    if "--shot" in sys.argv:
        app.processEvents()
        out = Path(sys.argv[sys.argv.index("--shot") + 1])
        w.grab().save(str(out))
        return 0
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
