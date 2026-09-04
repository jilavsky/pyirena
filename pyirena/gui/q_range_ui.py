"""
pyirena.gui.q_range_ui — one editable "Q range for fit" control for every tool.

Irena (Igor) let the user *type* the fit limits as well as drag the graph
cursors; pyIrena originally showed them as read-only readouts, so a user who
knew the exact Q they wanted had to nudge a cursor until the number matched.
:class:`QRangeFields` restores the Igor behaviour and — because every analysis
tool now embeds the same widget — makes the Q-range control identical in
Simple Fits, Modeling and Unified Fit.

The widget is deliberately graph-agnostic: it talks to the host's plot through
three callables, because the tools do not share a graph base class (Simple Fits
and Size Distribution expose ``get_cursor_range``/``set_cursor_range``, the
Modeling graph exposes ``get_q_range``/``set_q_range``).

Typical use::

    self.q_fields = QRangeFields(
        get_range=self.graph_window.get_cursor_range,
        set_range=self.graph_window.set_cursor_range,
        get_data_range=lambda: (self.data['Q'].min(), self.data['Q'].max()),
    )
    self.q_fields.message.connect(self.status_label.setText)
    self.q_fields.range_changed.connect(self._on_q_range_typed)

Editing rules (applied on *editingFinished*, i.e. Return or focus-out):

* Blank or unparsable input reverts to the current cursor positions.
* ``Qmin > Qmax`` is treated as the user entering the pair the wrong way round
  and is silently swapped, rather than rejected — that is what they meant.
* ``Qmin == Qmax`` is an empty range and is rejected (reverted, with a
  message), because every downstream fit would see zero points.
* Values are clamped to the loaded data's Q range when the host supplies one;
  a cursor outside the data does nothing useful and confuses the autoscale.
* The cursors are the single source of truth: after moving them the fields are
  re-read from the cursors, so what is displayed is always what will be fitted.
"""

from __future__ import annotations

import logging
from typing import Callable, Optional

from pyirena.gui._qt import (
    QDoubleValidator,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    Qt,
    QVBoxLayout,
    QWidget,
    Signal,
)
from pyirena.gui.theme import muted_css

log = logging.getLogger(__name__)

__all__ = ["QRangeFields"]

#: Text shown when there is no data / no cursors yet.
_EMPTY = ""


class QRangeFields(QWidget):
    """Editable ``Q min`` / ``Q max`` pair kept in sync with the plot cursors.

    Args:
        get_range: Returns the current cursor range as ``(q_min, q_max)``.
            May return ``None`` or ``(None, None)`` before data is loaded.
        set_range: Moves the cursors to ``(q_min, q_max)`` (linear Å⁻¹).
        get_data_range: Optional; returns the loaded data's
            ``(q_lo, q_hi)`` used to clamp typed values.  Return ``None``
            when no data is loaded and no clamping should happen.
        unit: Unit label shown after each field.
        hint: Whether to show the "drag cursors or type values" hint line.
        parent: Qt parent.

    Signals:
        range_changed(float, float): The user typed a new, accepted range and
            the cursors have already been moved to it.
        message(str): Human-readable feedback for the host's status label.
    """

    range_changed = Signal(float, float)
    message = Signal(str)

    def __init__(
        self,
        get_range: Callable[[], Optional[tuple]],
        set_range: Callable[[float, float], object],
        get_data_range: Optional[Callable[[], Optional[tuple]]] = None,
        *,
        unit: str = "Å⁻¹",
        hint: bool = True,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._get_range = get_range
        self._set_range = set_range
        self._get_data_range = get_data_range
        # Guards re-entrancy: moving the cursors makes the host call refresh(),
        # which writes the fields, which must not look like a fresh user edit.
        self._updating = False

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(2)

        if hint:
            lbl = QLabel("Drag cursors on the I(Q) graph, or type Q values here")
            lbl.setStyleSheet(muted_css("10px", italic=True))
            lbl.setWordWrap(True)
            outer.addWidget(lbl)

        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(4)

        row.addWidget(QLabel("Q min:"))
        self.q_min_edit = self._make_edit("lower fit limit (cursor A)")
        row.addWidget(self.q_min_edit)
        row.addWidget(QLabel(unit))

        row.addSpacing(8)
        row.addWidget(QLabel("Q max:"))
        self.q_max_edit = self._make_edit("upper fit limit (cursor B)")
        row.addWidget(self.q_max_edit)
        row.addWidget(QLabel(unit))
        row.addStretch()
        outer.addLayout(row)

        self.q_min_edit.editingFinished.connect(self._on_edited)
        self.q_max_edit.editingFinished.connect(self._on_edited)

    # -- construction helpers -------------------------------------------------

    def _make_edit(self, what: str) -> QLineEdit:
        e = QLineEdit()
        e.setMaximumWidth(90)
        e.setMinimumWidth(70)
        e.setAlignment(Qt.AlignmentFlag.AlignRight)
        # Bottom is 0, not -inf: Q is positive and the plots are log-x, so a
        # non-positive limit has no representable cursor position.
        e.setValidator(QDoubleValidator(0.0, 1e30, 12))
        e.setToolTip(
            f"Type the {what}, or drag the cursor on the graph.\n"
            "Values are clamped to the data's Q range; entering them the wrong\n"
            "way round swaps them."
        )
        return e

    # -- public API -----------------------------------------------------------

    def refresh(self) -> None:
        """Re-read the cursor positions into the fields.

        Call this whenever the cursors move (drag, data load, state restore).
        Safe to call while the user is typing: it only writes the fields, and
        Qt keeps the caret where it was for an unchanged string.
        """
        rng = self._current_range()
        self._updating = True
        try:
            q_min, q_max = (rng if rng else (None, None))
            self.q_min_edit.setText(_fmt(q_min))
            self.q_max_edit.setText(_fmt(q_max))
        finally:
            self._updating = False

    def values(self) -> tuple[Optional[float], Optional[float]]:
        """Return the currently displayed ``(q_min, q_max)``, or ``(None, None)``."""
        return _parse(self.q_min_edit.text()), _parse(self.q_max_edit.text())

    def setEnabled(self, enabled: bool) -> None:  # noqa: N802 (Qt naming)
        """Enable/disable both fields together."""
        super().setEnabled(enabled)

    # -- internals ------------------------------------------------------------

    def _current_range(self) -> Optional[tuple]:
        try:
            rng = self._get_range()
        except Exception:
            log.debug("cursor range read failed", exc_info=True)
            return None
        if not rng:
            return None
        q_min, q_max = rng
        if q_min is None or q_max is None:
            return None
        return float(q_min), float(q_max)

    def _data_range(self) -> Optional[tuple]:
        if self._get_data_range is None:
            return None
        try:
            rng = self._get_data_range()
        except Exception:
            log.debug("data range read failed", exc_info=True)
            return None
        if not rng:
            return None
        lo, hi = rng
        if lo is None or hi is None:
            return None
        return float(lo), float(hi)

    def _on_edited(self) -> None:
        """Validate the typed pair and move the cursors to match."""
        if self._updating:
            return

        q_min = _parse(self.q_min_edit.text())
        q_max = _parse(self.q_max_edit.text())

        if q_min is None or q_max is None or q_min <= 0 or q_max <= 0:
            self.refresh()
            self.message.emit("Q range unchanged — both Q min and Q max must be positive numbers.")
            return

        swapped = False
        if q_min > q_max:
            q_min, q_max = q_max, q_min
            swapped = True

        if q_min == q_max:
            self.refresh()
            self.message.emit("Q range unchanged — Q min and Q max must differ.")
            return

        note = ""
        data = self._data_range()
        if data is not None:
            lo, hi = min(data), max(data)
            clamped_min = min(max(q_min, lo), hi)
            clamped_max = min(max(q_max, lo), hi)
            if clamped_min == clamped_max:
                self.refresh()
                self.message.emit(
                    f"Q range unchanged — the requested range lies outside the "
                    f"data ({lo:.4g} … {hi:.4g} Å⁻¹)."
                )
                return
            if (clamped_min, clamped_max) != (q_min, q_max):
                note = f" (clamped to the data range {lo:.4g} … {hi:.4g} Å⁻¹)"
            q_min, q_max = clamped_min, clamped_max

        try:
            self._set_range(q_min, q_max)
        except Exception:
            log.debug("set cursor range failed", exc_info=True)
            self.refresh()
            self.message.emit("Could not move the cursors — load data first.")
            return

        # The cursors are authoritative: read back what they actually took.
        self.refresh()
        actual = self._current_range() or (q_min, q_max)
        if swapped and not note:
            note = " (Q min and Q max were swapped)"
        self.message.emit(
            f"Q range set to {actual[0]:.4g} … {actual[1]:.4g} Å⁻¹{note}"
        )
        self.range_changed.emit(actual[0], actual[1])


# ---------------------------------------------------------------------------
# tiny formatting helpers (module-level so tests can exercise them directly)
# ---------------------------------------------------------------------------

def _fmt(value: Optional[float]) -> str:
    """Format a Q value for display, or the empty string for ``None``."""
    if value is None:
        return _EMPTY
    try:
        return f"{float(value):.6g}"
    except (TypeError, ValueError):
        return _EMPTY


def _parse(text: str) -> Optional[float]:
    """Parse a field's text to float, or ``None`` when it is not a number."""
    text = (text or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None
