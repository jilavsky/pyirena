""""No limits?" means the same thing everywhere, and an agent can ask for it.

The check box has been in the GUI for four of the fitting tools and reachable
from the agent API for exactly one of them (Simple Fits), so a scientist could
release a stuck fit and an agent driving the same session could not. These
tests pin both halves: that every tool with bounds exposes ``no_limits``, and
that asking for it actually changes the fit rather than being accepted and
ignored — which is the failure mode that looks like success.

What it *means* differs by tool, and deliberately:

* Unified Fit and the Carbon model release the user's bounds back to each
  parameter's declared default. That is what the panels have always done, so
  the numbers are unchanged and the GUI and the API cannot disagree.
* Simple Fits and WAXS Peak Fit solve genuinely unconstrained (±inf).
* Modeling switches to an unconstrained Nelder-Mead solve.

The contract is the same in every case — *ignore the bounds I set* — and the
bounds themselves must survive the fit, because they belong to the setup and
not to one run of it.
"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from pyirena.api import control as ctrl

#: Tools whose fit is bounded, and so must offer the escape hatch. Size
#: Distribution is absent on purpose: it is an inversion, with no
#: per-parameter bounds for "no limits" to mean anything about.
TOOLS_WITH_BOUNDS = ["run_fit", "run_simple_fit", "run_modeling_fit",
                     "run_waxs_fit", "run_carbon_fit"]


@pytest.mark.parametrize("runner", TOOLS_WITH_BOUNDS)
def test_every_bounded_tool_offers_no_limits(runner):
    """Add a bounded tool, and it has to be reachable here too."""
    sig = inspect.signature(getattr(ctrl, runner))
    assert "no_limits" in sig.parameters, (
        f"{runner} has bounds but no way for an agent to release them"
    )
    assert sig.parameters["no_limits"].default is False, (
        f"{runner} must default to honouring the bounds"
    )


@pytest.mark.parametrize("runner", TOOLS_WITH_BOUNDS)
def test_the_mcp_schema_advertises_it(runner):
    """An argument an agent cannot see is an argument it will not use."""
    from pyirena.api.control.schemas import TOOL_SCHEMAS

    schema = next(t for t in TOOL_SCHEMAS if t["name"] == runner)
    props = schema["input_schema"]["properties"]
    assert "no_limits" in props, f"{runner}'s schema does not mention no_limits"
    assert props["no_limits"]["type"] == "boolean"
    assert len(props["no_limits"].get("description", "")) > 40, (
        "an agent needs to be told what it does, not just that it exists"
    )


def _curve():
    q = np.logspace(np.log10(5e-4), np.log10(0.4), 200)
    intensity = 3e4 * np.exp(-(q * 280.0) ** 2 / 3.0) + 1.2e-4 * q ** -3.9 + 0.02
    return q, intensity, intensity * 0.03


def _unified_session():
    q, intensity, error = _curve()
    opened = ctrl.open_dataset_from_data(q=q.tolist(), intensity=intensity.tolist(),
                                         error=error.tolist())
    sid = opened["session_id"]
    ctrl.select_model(sid, "unified_fit", nlevels=1)
    ctrl.set_parameter_value(sid, "Rg_1", 150.0)
    # Deliberately too tight: the real Rg is ~280 Å.
    ctrl.set_parameter_bounds(sid, "Rg_1", 100.0, 200.0)
    return sid


def test_a_bounded_unified_fit_stops_at_the_bound():
    """The baseline the escape hatch exists for."""
    sid = _unified_session()
    try:
        result = ctrl.run_fit(sid, walk_limits=False)
        rg = next(p["value"] for p in result["parameters_updated"]
                  if p["name"] == "Rg_1")
        assert rg == pytest.approx(200.0, rel=1e-6)
        assert "level 1 Rg" in result["pinned_parameters"]
    finally:
        ctrl.close_session(sid)


def test_no_limits_lets_the_unified_fit_past_the_bound():
    sid = _unified_session()
    try:
        bounded = ctrl.run_fit(sid, walk_limits=False)
        ctrl.set_parameter_value(sid, "Rg_1", 150.0)
        released = ctrl.run_fit(sid, walk_limits=False, no_limits=True)

        rg = next(p["value"] for p in released["parameters_updated"]
                  if p["name"] == "Rg_1")
        assert rg > 200.0, "no_limits did not release the bound"
        assert released["chi_squared"] < bounded["chi_squared"], (
            "releasing the bound should not make the fit worse"
        )
        assert released["no_limits"] is True
        # Not "pinned at a limit" — the fit never used those limits. The
        # values landing outside them is the answer, and is reported as such.
        assert released["pinned_parameters"] == []
        assert any("no_limits" in w for w in released["warnings"])
    finally:
        ctrl.close_session(sid)


def test_the_bounds_survive_a_no_limits_fit():
    """They belong to the setup, not to one run of it.

    If a no_limits fit left the wide defaults behind, the next ordinary fit
    would silently be unbounded too — and the saved setup would no longer be
    the one the scientist configured.
    """
    sid = _unified_session()
    try:
        ctrl.run_fit(sid, walk_limits=False, no_limits=True)
        row = next(p for p in ctrl.get_model_parameters(sid)["parameters"]
                   if p["name"] == "Rg_1")
        assert (row["lo"], row["hi"]) == (100.0, 200.0)
    finally:
        ctrl.close_session(sid)


def test_analyze_honours_no_limits_from_a_config():
    """A saved setup with the box ticked has to replay with it ticked."""
    q, intensity, error = _curve()
    data = {"q": q.tolist(), "intensity": intensity.tolist(), "error": error.tolist()}
    level = {
        "G": 3e4, "Rg": 150.0, "B": 1.2e-4, "P": 3.9, "ETA": 10.0, "PACK": 0.0,
        "RgCO": 0.0, "K": 1.0,
        "fit_G": True, "fit_Rg": True, "fit_B": True, "fit_P": False,
        "fit_ETA": False, "fit_PACK": False, "fit_RgCO": False,
        "correlations": False, "mass_fractal": False,
        "link_B": False, "link_RGCO": False,
        "Rg_limits": [100.0, 200.0], "G_limits": [1.0, 1e10],
        "B_limits": [1e-20, 1e10], "P_limits": [0.0, 6.0],
        "ETA_limits": [0.1, 1e6], "PACK_limits": [0.0, 16.0],
        "RgCO_limits": [0.0, 1e6],
    }

    def _run(no_limits):
        section = {"num_levels": 1, "levels": [dict(level)],
                   "cursor_left": 5e-4, "cursor_right": 0.4,
                   "no_limits": no_limits}
        report = ctrl.analyze(data, {"_pyirena_config": {"tool": "unified_fit"},
                                     "unified_fit": section})
        assert "error" not in report, report
        return report["config"]["levels"][0]["Rg"]

    assert _run(True) != _run(False), (
        "the config's no_limits was accepted and ignored"
    )


# ── The core context managers restore what they borrowed ─────────────────

def test_unified_default_bounds_restores_the_users_bounds():
    from pyirena.core.unified import UnifiedFitModel

    model = UnifiedFitModel(num_levels=2)
    model.levels[0].Rg_limits = (100.0, 400.0)
    model.background_limits = (0.0, 1.0)
    with model.default_bounds(True):
        assert model.levels[0].Rg_limits != (100.0, 400.0)
    assert model.levels[0].Rg_limits == (100.0, 400.0)
    assert model.background_limits == (0.0, 1.0)


def test_unified_default_bounds_is_a_no_op_when_off():
    from pyirena.core.unified import UnifiedFitModel

    model = UnifiedFitModel(num_levels=1)
    model.levels[0].Rg_limits = (100.0, 400.0)
    with model.default_bounds(False):
        assert model.levels[0].Rg_limits == (100.0, 400.0)


def test_carbon_default_bounds_restores_and_keeps_the_safety_clamp():
    """Releasing the user's bounds must not release the gradient clamp.

    Outside (1.001, 2.999) the Teixeira structure factor is flat, so the
    finite-difference derivative is exactly zero and the fit burns its whole
    budget without moving — indistinguishable from a parameter that was never
    wired up. "No limits" must not mean "no gradient".
    """
    from pyirena.core.carbon_fit import _SAFE_BOUNDS, CarbonFitModel

    model = CarbonFitModel()
    model.background.S_macro_limits = (1000.0, 2000.0)
    # The natural thing for a user to type, and the one that put D in the
    # dead zone before the clamp existed.
    model.saxs.fractal_D_limits = (1.0, 3.0)
    safe_lo, safe_hi = _SAFE_BOUNDS["fractal_D"]

    with model.default_bounds(True):
        assert model.background.S_macro_limits != (1000.0, 2000.0)
        for ref in model.parameter_refs(active_only=False):
            if ref.attr != "fractal_D":
                continue
            (lo, hi), _ = model.safe_bounds(ref)
            assert lo >= safe_lo and hi <= safe_hi, (
                f"{ref.key} escaped the gradient clamp under no_limits: "
                f"({lo}, {hi})"
            )
    assert model.background.S_macro_limits == (1000.0, 2000.0)
    assert model.saxs.fractal_D_limits == (1.0, 3.0)
