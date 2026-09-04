"""
Tests for the shared, editable Q-range control (:mod:`pyirena.gui.q_range_ui`).

Irena let the user type the fit limits; pyIrena showed them read-only, so a
user who wanted Q = 0.01 had to nudge a cursor until the readout matched.
:class:`~pyirena.gui.q_range_ui.QRangeFields` restores typing *and* is the same
widget in Simple Fits, Modeling and Unified Fit, so the three tools cannot
drift apart again.

The behaviours worth pinning are the edge cases the user will actually hit:
reversed entry, a degenerate range, junk text, and values outside the data.
"""

from __future__ import annotations

import numpy as np
import pytest


def _require_qt() -> None:
    try:
        import pyirena.gui._qt  # noqa: F401
    except ImportError:
        pytest.skip("Qt (PySide6/PyQt6) not available", allow_module_level=True)


_require_qt()
pytest.importorskip("pyqtgraph")

from pyirena.gui._qt import QApplication  # noqa: E402
from pyirena.gui.q_range_ui import QRangeFields  # noqa: E402

_APP = None


@pytest.fixture(scope="session")
def qapp():
    global _APP
    _APP = QApplication.instance() or QApplication([])
    yield _APP


class _FakeCursors:
    """Minimal stand-in for a panel's graph window."""

    def __init__(self, lo=0.01, hi=0.5):
        self.lo, self.hi = lo, hi

    def get(self):
        return (self.lo, self.hi)

    def set(self, lo, hi):
        self.lo, self.hi = lo, hi


@pytest.fixture
def fields(qapp):
    cursors = _FakeCursors()
    w = QRangeFields(
        get_range=cursors.get,
        set_range=cursors.set,
        get_data_range=lambda: (1e-3, 1.0),
    )
    w.refresh()
    w._cursors = cursors      # test convenience
    return w


def _type(w, lo, hi):
    w.q_min_edit.setText(str(lo))
    w.q_max_edit.setText(str(hi))
    w._on_edited()


def test_refresh_shows_the_cursor_range(fields):
    assert fields.values() == (0.01, 0.5)


def test_typed_range_moves_the_cursors(fields):
    _type(fields, 0.02, 0.3)
    assert fields._cursors.get() == (0.02, 0.3)
    assert fields.values() == (0.02, 0.3)


def test_reversed_entry_is_swapped_not_rejected(fields):
    """Typing the pair the wrong way round means the range, not an error."""
    _type(fields, 0.3, 0.02)
    lo, hi = fields.values()
    assert lo < hi
    assert (lo, hi) == (0.02, 0.3)


def test_equal_limits_are_rejected(fields):
    before = fields.values()
    _type(fields, 0.05, 0.05)
    assert fields.values() == before, "an empty Q range must not be accepted"
    assert fields._cursors.get() == before


def test_non_numeric_input_reverts(fields):
    before = fields.values()
    fields.q_min_edit.setText("not a number")
    fields._on_edited()
    assert fields.values() == before


def test_non_positive_input_reverts(fields):
    """Q is positive and the plots are log-x: 0 has no cursor position."""
    before = fields.values()
    _type(fields, 0, 0.3)
    assert fields.values() == before


def test_values_are_clamped_to_the_data_range(fields):
    _type(fields, 1e-9, 1e6)
    lo, hi = fields.values()
    assert lo == pytest.approx(1e-3)
    assert hi == pytest.approx(1.0)


def test_range_entirely_outside_the_data_is_rejected(fields):
    before = fields.values()
    _type(fields, 10.0, 100.0)
    assert fields.values() == before


def test_range_changed_signal_carries_the_accepted_values(fields):
    seen = []
    fields.range_changed.connect(lambda lo, hi: seen.append((lo, hi)))
    _type(fields, 0.02, 0.3)
    assert seen == [(0.02, 0.3)]


def test_rejected_edit_emits_no_range_changed(fields):
    seen = []
    fields.range_changed.connect(lambda lo, hi: seen.append((lo, hi)))
    _type(fields, 0.05, 0.05)
    assert seen == []


def test_message_signal_explains_a_rejection(fields):
    msgs = []
    fields.message.connect(msgs.append)
    _type(fields, 0.05, 0.05)
    assert msgs and "must differ" in msgs[-1]


# --- the three panels really share this widget ------------------------------

_PANELS = [
    ("pyirena.gui.simple_fits_panel", "SimpleFitsPanel", "set_data", "label"),
    ("pyirena.gui.modeling_panel", "ModelingPanel", "set_data", "filename"),
    ("pyirena.gui.unified_fit", "UnifiedFitPanel", "set_data", "label"),
]


@pytest.mark.parametrize("module, cls_name, setter, name_kw", _PANELS)
def test_every_fit_panel_exposes_editable_q_range(qapp, module, cls_name,
                                                  setter, name_kw):
    """Simple Fits, Modeling and Unified Fit must behave identically here."""
    import importlib

    panel = getattr(importlib.import_module(module), cls_name)()
    q = np.logspace(-3, 0, 200)
    intensity = 1e3 * q ** -3 + 0.1
    getattr(panel, setter)(q, intensity, 0.02 * intensity, **{name_kw: "t"})

    w = panel.q_range_fields
    assert not w.q_min_edit.isReadOnly(), "Q min must be editable"
    assert not w.q_max_edit.isReadOnly(), "Q max must be editable"

    _type(w, 0.01, 0.1)
    lo, hi = w.values()
    assert lo == pytest.approx(0.01)
    assert hi == pytest.approx(0.1)
    # the cursors, not the fields, are what the fit reads
    cursor_lo, cursor_hi = w._get_range()
    assert cursor_lo == pytest.approx(0.01)
    assert cursor_hi == pytest.approx(0.1)
