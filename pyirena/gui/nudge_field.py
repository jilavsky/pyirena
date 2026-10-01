"""
A numeric field you can nudge — mouse wheel, arrow keys, or ▲▼ buttons.

Merging two SAS curves by eye needs small, repeated adjustments to a scale or
a background while watching the overlap: type-a-number-and-press-Enter is the
wrong interaction for that, which is the gap against Igor that GitHub issue
#35 is about.

Steps are a **percentage of the current value**, not an absolute amount, so
one control works for a scale near 1 and a background near 1e-5 without
anybody configuring it. Modifiers follow the convention already used by the
fit panels' scrubbable fields: plain = 1%, Shift = 10% (coarse), Ctrl/Cmd =
0.1% (fine).

Zero is the case a percentage cannot handle — a background usually *starts* at
zero, and 1% of zero is zero, so the field would be dead exactly when the user
reaches for it. :meth:`NudgeField.set_fallback_step` takes an absolute step
for that case; the Data Merge panel sets it from the data it has loaded, which
is the only place that knows what a meaningful background step is.

Related, and deliberately not merged with this
----------------------------------------------
Four panels already grew their own wheel-editable line edit — ``unified_fit``
and ``sizes_panel`` (both named ``ScrubbableLineEdit``, different code),
``data_manipulation_panel._ScalableLineEdit`` and
``waxs_peakfit_panel._ValidatedField``. None has arrow buttons or a zero
fallback. Migrating them is a behaviour change in four fitting panels and is
not part of issue #35; this module is where that consolidation should land
when someone takes it on.
"""
from __future__ import annotations

import logging
import math
from typing import Optional

from pyirena.gui._qt import (
    QApplication,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    Qt,
    QWidget,
    Signal,
)

log = logging.getLogger(__name__)

#: Step as a fraction of the current value, by modifier.
STEP_PLAIN = 0.01     # 1%
STEP_COARSE = 0.10    # 10% — Shift
STEP_FINE = 0.001     # 0.1% — Ctrl/Cmd


class NudgeField(QLineEdit):
    """A line edit holding a float that the wheel and arrow keys can nudge.

    Emits :attr:`nudged` after every change it makes itself, so a panel can
    recompute without also firing on every keystroke of a typed value.
    ``editingFinished`` still behaves as usual for typing.
    """

    #: Emitted after a wheel/arrow/button nudge changes the value.
    nudged = Signal()

    def __init__(self, text: str = "", parent: Optional[QWidget] = None,
                 fallback_step: float = 1e-6, decimals: int = 6):
        super().__init__(text, parent)
        self.setFocusPolicy(Qt.FocusPolicy.WheelFocus)
        self._fallback_step = float(fallback_step)
        self._decimals = int(decimals)

    # ── configuration ────────────────────────────────────────────────────

    def set_fallback_step(self, step: float) -> None:
        """Set the absolute step used when the current value is zero.

        A percentage of zero is zero, so without this a field sitting at 0 —
        which is where a background starts — cannot be nudged at all.
        """
        try:
            step = float(step)
        except (TypeError, ValueError):
            return
        if math.isfinite(step) and step > 0:
            self._fallback_step = step

    def value(self) -> Optional[float]:
        """The current value, or None when the text is not a number."""
        try:
            return float(self.text().strip())
        except (ValueError, AttributeError):
            return None

    # ── the nudge itself ─────────────────────────────────────────────────

    def nudge(self, direction: int, fraction: float = STEP_PLAIN) -> None:
        """Change the value by *fraction* of itself, *direction* being ±1."""
        value = self.value()
        if value is None:
            return
        step = abs(value) * fraction
        if step <= 0:
            # Zero (or denormal): fall back to the absolute step.
            step = self._fallback_step
        new = value + direction * step
        if not math.isfinite(new):
            return
        self.setText(f"{new:.{self._decimals}g}")
        self.nudged.emit()

    @staticmethod
    def _fraction_for(modifiers) -> float:
        if modifiers & Qt.KeyboardModifier.ShiftModifier:
            return STEP_COARSE
        # ControlModifier is ⌘ on macOS, which is what a Mac user reaches for.
        if modifiers & Qt.KeyboardModifier.ControlModifier:
            return STEP_FINE
        return STEP_PLAIN

    # ── Qt events ────────────────────────────────────────────────────────

    def wheelEvent(self, event):              # noqa: N802 (Qt override)
        if self.isReadOnly():
            super().wheelEvent(event)
            return
        delta = event.angleDelta().y()
        if delta == 0:
            super().wheelEvent(event)
            return
        self.nudge(1 if delta > 0 else -1, self._fraction_for(event.modifiers()))
        # Swallow it, or an enclosing scroll area scrolls as well and the view
        # jumps away from the curve the user is watching.
        event.accept()

    def keyPressEvent(self, event):           # noqa: N802 (Qt override)
        key = event.key()
        if not self.isReadOnly() and key in (Qt.Key.Key_Up, Qt.Key.Key_Down):
            self.nudge(1 if key == Qt.Key.Key_Up else -1,
                       self._fraction_for(event.modifiers()))
            event.accept()
            return
        super().keyPressEvent(event)


def make_nudge_buttons(field: NudgeField, width: int = 16) -> QWidget:
    """A stacked ▲ / ▼ pair that nudges *field*, for people who do not scroll.

    Returned as one widget so a caller can drop it into a row next to the
    field with a single ``addWidget``.
    """
    holder = QWidget()
    box = QHBoxLayout(holder)
    box.setContentsMargins(0, 0, 0, 0)
    box.setSpacing(0)

    for glyph, direction, tip in (
        ("▲", +1, "Increase by 1% (Shift: 10%, Ctrl/Cmd: 0.1%)"),
        ("▼", -1, "Decrease by 1% (Shift: 10%, Ctrl/Cmd: 0.1%)"),
    ):
        btn = QPushButton(glyph)
        btn.setFixedSize(width, 22)
        btn.setToolTip(tip + "\nThe mouse wheel and ↑/↓ over the field do the same.")
        btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        # A background with no `color:` inherits the system text colour and
        # the glyph vanishes on a dark desktop — see
        # pyirena/tests/test_gui_theme_contract.py, which caught exactly that.
        btn.setStyleSheet(
            "QPushButton{border:1px solid #bdc3c7;background:#ecf0f1;"
            "color:#2c3e50;font-size:8pt;padding:0px;}"
            "QPushButton:hover{background:#d5dbdb;color:#2c3e50;}"
        )
        # Read the modifiers at click time, so Shift-click is coarse exactly
        # as Shift-wheel is.
        def _on_click(_checked=False, d=direction):
            if field.isReadOnly():
                return
            field.nudge(d, NudgeField._fraction_for(
                QApplication.keyboardModifiers()))

        btn.clicked.connect(_on_click)
        box.addWidget(btn)
    return holder
