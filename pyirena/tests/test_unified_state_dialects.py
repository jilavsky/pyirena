"""The Unified Fit panel reads both config dialects and writes only one.

The panel's setup state is the one thing in pyIrena that users already hold in
two vocabularies: the older *panel* dialect, in every ``state.json`` and every
``_pyirena_config`` attribute written before 1.2, and the *core* dialect
(``UnifiedLevel.to_dict()``) that everything writes from now on
(``planning/config-dialects/`` §4 rule 3).

Reading forever is the promise; these tests are what makes it one. A live
panel is driven through both, because the failure mode is not an exception —
it is a level that comes back at its dataclass defaults and fits to a
plausible wrong answer.
"""

from __future__ import annotations

import pytest


def _require_qt() -> None:
    try:
        import pyirena.gui._qt  # noqa: F401
    except ImportError:
        pytest.skip("Qt (PySide6/PyQt6) not available", allow_module_level=True)


_require_qt()
pytest.importorskip("pyqtgraph")

from pyirena.core.unified import (  # noqa: E402
    UNIFIED_STATE_SCHEMA_VERSION,
    UnifiedLevel,
)
from pyirena.gui._qt import QApplication  # noqa: E402

_APP = None


@pytest.fixture(scope="session")
def qapp():
    global _APP
    _APP = QApplication.instance() or QApplication([])
    yield _APP


@pytest.fixture
def panel(qapp):
    from pyirena.gui.unified_fit import UnifiedFitPanel

    p = UnifiedFitPanel()
    yield p
    p.deleteLater()


# ── A legacy state file still restores every control ─────────────────────

def _legacy_level_state(level_number: int) -> dict:
    """One level exactly as pyIrena wrote it before 1.2."""
    return {
        "level": level_number,
        "G":    {"value": 31300.0, "fit": True,  "low_limit": 6260.0, "high_limit": 156000.0},
        "Rg":   {"value": 286.0,   "fit": False, "low_limit": 57.2,   "high_limit": 1430.0},
        "B":    {"value": 1.13e-4, "fit": True,  "low_limit": 2.26e-5, "high_limit": 5.65e-4},
        "P":    {"value": 4.1,     "fit": False, "low_limit": 2.05,   "high_limit": 5.0},
        "ETA":  {"value": 146.0,   "fit": True,  "low_limit": 29.2,   "high_limit": 730.0},
        "PACK": {"value": 0.314,   "fit": True,  "low_limit": 0.0628, "high_limit": 1.57},
        "RgCutoff":   12.0,
        "correlated": True,
        "estimate_B": False,
        "link_rgco":  False,
    }


LEGACY_STATE = {
    "num_levels": 2,
    "levels": [_legacy_level_state(i + 1) for i in range(5)],
    "background": {"value": 0.017, "fit": False},
    "cursor_left": 0.001,
    "cursor_right": 0.1,
    "update_auto": False,
    "display_local": False,
    "no_limits": False,
}


def test_a_state_file_written_before_1_2_still_restores_every_control(panel):
    """The panel dialect is read forever. This is that promise, executed."""
    panel.apply_state(LEGACY_STATE)

    assert panel.num_levels_spin.value() == 2
    assert float(panel.background_value.text()) == pytest.approx(0.017, rel=1e-3)
    assert panel.fit_background_check.isChecked() is False

    widget = panel.level_widgets[0]
    params = widget.get_parameters()
    assert params["G"] == pytest.approx(31300.0)
    assert params["Rg"] == pytest.approx(286.0)
    assert params["P"] == pytest.approx(4.1)
    # The decisions, not just the values — this is the half that was lost.
    assert params["fit_G"] is True
    assert params["fit_Rg"] is False
    assert params["fit_P"] is False
    assert params["fit_ETA"] is True
    assert params["G_low"] == pytest.approx(6260.0)
    assert params["G_high"] == pytest.approx(156000.0)
    assert params["Rg_low"] == pytest.approx(57.2)
    assert params["correlated"] is True


def test_a_legacy_state_comes_back_out_in_the_core_dialect(panel):
    """Read the old vocabulary, write the new one — the migration, in one line."""
    panel.apply_state(LEGACY_STATE)
    state = panel.get_current_state()

    assert state["schema_version"] == UNIFIED_STATE_SCHEMA_VERSION
    assert isinstance(state["background"], float)
    assert "fit_background" in state

    level = state["levels"][0]
    assert isinstance(level["Rg"], float), "a level parameter must be a bare number now"
    assert level["fit_Rg"] is False
    assert level["Rg_limits"] == [pytest.approx(57.2), pytest.approx(1430.0)]
    assert level["correlations"] is True
    # RgCO is asserted on level 2: the panel forces level 1's cutoff to zero,
    # because there is no previous level for it to cut off against.
    assert state["levels"][1]["RgCO"] == pytest.approx(12.0)
    # The panel's own key names are gone from what is written.
    for retired in ("RgCutoff", "correlated", "estimate_B", "link_rgco", "value"):
        assert retired not in level


# ── The round trip the panel does on every save and launch ───────────────

def test_the_panel_round_trips_its_own_state_unchanged(panel):
    """collect → apply → collect. Anything the panel drops shows up here."""
    panel.apply_state(LEGACY_STATE)
    once = panel.get_current_state()
    panel.apply_state(once)
    twice = panel.get_current_state()

    assert twice == once


def test_the_shipped_defaults_apply_without_widening_the_limit_fields(panel):
    """A default state states no bounds, so the widget keeps the ones it ships.

    The level widget starts at G low "20", B low "0.002", P low "2" — narrow,
    GUI-friendly values chosen for the panel. The model's defaults are
    1e-10, 1e-20 and 0. Applying a bound-less state must not replace one with
    the other, or every user's limit fields change on upgrade.
    """
    from pyirena.state.state_manager import StateManager

    before = panel.level_widgets[0].get_parameters()
    panel.apply_state(StateManager.DEFAULT_STATE["unified_fit"])
    after = panel.level_widgets[0].get_parameters()

    for key in ("G_low", "G_high", "Rg_low", "Rg_high", "B_low", "B_high",
                "P_low", "P_high"):
        assert after[key] == before[key], f"{key} was overwritten by a state that did not state it"


def test_a_core_dialect_state_restores_the_same_controls(panel):
    """The dialect the panel writes has to read back into the same widgets."""
    level = UnifiedLevel(Rg=250.0, G=5000.0, P=3.2, B=1e-5, ETA=40.0, PACK=2.0)
    level.fit_Rg = False
    level.fit_P = True
    level.Rg_limits = (100.0, 400.0)
    level.correlations = True
    level.RgCO = 17.0

    # Level 2, because the panel pins level 1's RgCutoff to zero.
    panel.apply_state({
        "schema_version": UNIFIED_STATE_SCHEMA_VERSION,
        "num_levels": 2,
        "levels": [UnifiedLevel().to_dict(), level.to_dict()],
        "background": 0.02,
        "fit_background": True,
    })

    params = panel.level_widgets[1].get_parameters()
    assert params["Rg"] == pytest.approx(250.0)
    assert params["G"] == pytest.approx(5000.0)
    assert params["fit_Rg"] is False
    assert params["fit_P"] is True
    assert params["Rg_low"] == pytest.approx(100.0)
    assert params["Rg_high"] == pytest.approx(400.0)
    assert params["RgCutoff"] == pytest.approx(17.0)
    assert params["correlated"] is True
    assert panel.fit_background_check.isChecked() is True


# ── The state the agent API embeds is the same shape the panel writes ────

def test_the_control_api_embeds_what_the_panel_can_read(panel):
    """``save_fit`` writes a setup the GUI has to be able to open.

    The two used to be built by separate hand-written mappings — the panel's
    ``_collect_state`` and ``api/control/unified_fit._session_to_gui_state``
    — which is two chances to drop a field and no test that they agree.
    """
    import numpy as np

    from pyirena.api import control as ctrl

    q = np.logspace(-3, -1, 80)
    intensity = 5000.0 * np.exp(-(q * 250.0) ** 2 / 3.0) + 1e-4 * q ** -4 + 0.01
    opened = ctrl.open_dataset_from_data(q=q.tolist(), intensity=intensity.tolist())
    session_id = opened["session_id"]
    try:
        ctrl.select_model(session_id, "unified_fit", nlevels=2)
        from pyirena.api.control.session import get_session
        from pyirena.api.control.unified_fit import _session_to_gui_state

        state = _session_to_gui_state(get_session(session_id))
        panel.apply_state(state)
        assert panel.get_current_state()["levels"][0] == state["levels"][0]
    finally:
        ctrl.close_session(session_id)
