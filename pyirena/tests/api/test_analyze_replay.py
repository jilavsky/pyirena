"""Replaying a saved setup must reproduce the fit, not a plausible neighbour.

``analyze()`` exists so that a setup made once in the GUI can be applied to
every measurement that follows. Each check here is one way that promise was
broken without anything failing: the fit returned ``success``, the report
looked ordinary, and a number was wrong — by a factor of two hundred in the
Invariant case.

They use synthetic curves with a known answer rather than
``testData/Scripting/``, deliberately. The exported fixtures are the right
test of *what the GUI writes*; these are a test of *what happens to a setting
on the way through*, and that needs a case constructed so the loss changes
the answer visibly. Each test states the number both paths must agree on.

Companion to ``tests/test_config_contract.py``, which checks the same
settings survive the round trip, and to ``tests/api/test_analyze.py``, which
checks the real exported files still replay.
"""
from __future__ import annotations

import numpy as np
import pytest

import pyirena.api.control as ctrl
from pyirena.core.simple_fits import SimpleFitModel
from pyirena.core.tool_config import build_setup, resolve_fixed_params


def _payload(q, intensity, error=None, **extra):
    if error is None:
        error = 0.03 * intensity
    return dict(q=np.asarray(q).tolist(),
                intensity=np.asarray(intensity).tolist(),
                error=np.asarray(error).tolist(), **extra)


def _guinier(rg_true=20.0, rg_start=12.0, i0=10.0, n=100):
    """A Guinier curve whose true Rg is nowhere near the starting value.

    The gap is the point: a parameter that was held fixed stays at the start,
    and one that was wrongly freed runs to the truth, so the two outcomes are
    never confusable with fit noise.
    """
    q = np.linspace(0.01, 0.12, n)
    intensity = i0 * np.exp(-q * q * rg_true**2 / 3.0)
    model = SimpleFitModel()
    model.params.update(I0=i0, Rg=rg_start)
    return q, intensity, model


@pytest.fixture(autouse=True)
def _no_sessions_left_behind():
    yield
    assert ctrl.list_open_sessions()["count"] == 0, "analyze leaked a session"


# ── Held parameters survive every way they are written down ──────────────

def test_a_report_fed_back_in_still_holds_what_was_held():
    """The documented replay loop: analyze, then analyze the reply.

    ``model.to_dict()`` has no room for the held-parameter choices, so a
    report that did not state them separately freed on the second pass
    whatever the scientist had pinned on the first.
    """
    q, intensity, model = _guinier()
    data = _payload(q, intensity)

    first = ctrl.analyze(data, {"simple_fits": {**model.to_dict(),
                                                "param_fixed": {"Rg": True}}})
    assert "error" not in first, first
    assert first["config"]["params"]["Rg"] == pytest.approx(12.0)
    assert first["analyze"]["fixed_parameters"] == ["Rg"]

    second = ctrl.analyze(data, first)
    assert "error" not in second, second
    assert second["config"]["params"]["Rg"] == pytest.approx(12.0)
    assert second["analyze"]["fixed_parameters"] == ["Rg"]


def test_the_setup_saved_into_a_result_file_still_holds_what_was_held():
    """``save_simple_fit()`` writes a ``fixed_params`` list; read it as one.

    The writer and the reader disagreed about the spelling — a list of names
    against a dict of checkbox states — so a setup read back out of an
    agent-written ``.h5`` lost the choices it was written to preserve.
    """
    _, _, model = _guinier()
    setup = build_setup({"_pyirena_config": {"tool": "simple_fits"},
                         "state": {**model.to_dict(), "fixed_params": ["Rg"]}})
    assert setup.fixed_params == {"Rg": 12.0}


@pytest.mark.parametrize("spec", [
    {"Rg": True, "I0": False},          # the GUI's checkbox dict
    ["Rg"],                             # the list saved into an .h5
    {"Rg": 12.0},                       # already resolved to {name: value}
])
def test_every_spelling_of_held_reads_the_same(spec):
    params = {"I0": 10.0, "Rg": 12.0}
    assert resolve_fixed_params(spec, params) == {"Rg": 12.0}


def test_a_name_with_no_parameter_behind_it_is_dropped_not_invented():
    assert resolve_fixed_params(["Rg", "Nonsense"], {"Rg": 1.0}) == {"Rg": 1.0}


# ── The wrapper must not hide the execution settings ─────────────────────

@pytest.mark.parametrize("wrap", [
    pytest.param(lambda sec: {"simple_fits": sec}, id="sidecar"),
    pytest.param(lambda sec: {"_pyirena_config": {"tool": "simple_fits"},
                              "state": sec}, id="hdf5"),
    pytest.param(lambda sec: {"tool": "simple_fits", "model": sec},
                 id="minimal"),
])
def test_no_limits_reaches_the_runner_in_every_envelope(wrap):
    """A setting beside the model is only as good as the reader that finds it.

    Model construction unwrapped all four envelopes; the runner dispatch
    still looked one level up, so for three of them ``no_limits`` — and the
    fit method, and the weighting — simply were not there. The fit stopped
    dead at a bound it had been told to ignore, and reported success.
    """
    q, intensity, model = _guinier()
    model.limits["Rg"] = (5.0, 15.0)
    section = {**model.to_dict(), "no_limits": True}

    report = ctrl.analyze(_payload(q, intensity), wrap(section))
    assert "error" not in report, report
    # Without no_limits the fit stops at the 15 Å bound instead of the truth.
    assert report["config"]["params"]["Rg"] == pytest.approx(20.0, rel=1e-4)


def test_carbon_keeps_the_weighting_its_config_restored():
    """Carbon nests ``weighting`` inside its model block, where the runner
    mapping was not looking — so the runner's ``auto`` default overwrote it.

    Weighting is not cosmetic here: this curve rises linearly with unit
    errors, so sigma-equivalent and relative weighting pull the background to
    visibly different places.
    """
    from pyirena.core.carbon_fit import CarbonFitModel

    q = np.linspace(0.01, 0.12, 100)
    intensity = np.linspace(1.0, 10.0, 100)
    error = np.ones(100)

    model = CarbonFitModel()
    model.saxs.enabled = False
    model.waxs.enabled = False
    model.background.S_macro = 0.0
    model.background.fit_S_macro = False
    model.background.fit_porod_exponent = False
    model.background.flat_background = 1.0
    model.background.fit_flat_background = True
    model.weighting = "relative"
    model.n_mc_runs = 0

    expected = model.copy().fit(q, intensity, error)
    report = ctrl.analyze(_payload(q, intensity, error),
                          {"carbon_fit": {"model": model.to_dict(),
                                          "q_min": None, "q_max": None}})
    assert "error" not in report, report
    assert report["config"]["weighting"] == "relative"
    assert report["config"]["background"]["flat_background"] == pytest.approx(
        expected.params["background.flat_background"], rel=1e-6)


# ── Slit smearing is a property of the measurement ───────────────────────

def test_a_slit_smeared_curve_is_fitted_smeared_even_from_a_pinhole_config():
    """The session API smears from the data's own metadata; so must a replay.

    A config exported from a pinhole run says nothing about slit smearing.
    Replayed on a slit-smeared curve it fitted the raw data against an
    unsmeared model, and the recovered I0 was ~10% low with no warning.
    """
    q = np.linspace(0.01, 0.1, 100)
    truth = SimpleFitModel()
    truth.params.update(I0=10.0, Rg=20.0)
    truth.use_slit_smearing = True
    truth.slit_length = 0.05
    data = _payload(q, truth.compute(q), is_slit_smeared=True, slit_length=0.05)

    report = ctrl.analyze(data, {"simple_fits": SimpleFitModel().to_dict()})
    assert "error" not in report, report
    assert report["config"]["use_slit_smearing"] is True
    assert report["config"]["slit_length"] == pytest.approx(0.05)
    assert report["config"]["params"]["I0"] == pytest.approx(10.0, rel=1e-4)

    # …and the same answer the step-by-step session path gives.
    sid = ctrl.open_dataset_from_data(**data)["session_id"]
    try:
        ctrl.select_simple_model(sid, "Guinier")
        ctrl.run_simple_fit(sid)
        session = ctrl.export_results(sid)
    finally:
        ctrl.close_session(sid)
    assert report["config"]["params"]["I0"] == pytest.approx(
        session["config"]["params"]["I0"], rel=1e-6)


def test_a_config_that_asks_for_smearing_keeps_its_own_slit_length():
    """A length the scientist set outranks the file's; neither is invented."""
    from pyirena.core.tool_config import apply_data_slit_settings

    model = SimpleFitModel()
    apply_data_slit_settings(model, {"use_slit_smearing": True,
                                     "slit_length": 0.02},
                             data_is_slit_smeared=True, data_slit_length=0.05)
    assert (model.use_slit_smearing, model.slit_length) == (True, 0.02)

    # A stale length beside an unticked box does not outrank the curve's.
    model = SimpleFitModel()
    apply_data_slit_settings(model, {"use_slit_smearing": False,
                                     "slit_length": 0.02},
                             data_is_slit_smeared=True, data_slit_length=0.05)
    assert (model.use_slit_smearing, model.slit_length) == (True, 0.05)

    # Pinhole data and a pinhole config leave the model alone.
    model = SimpleFitModel()
    apply_data_slit_settings(model, {"use_slit_smearing": False,
                                     "slit_length": 0.0})
    assert model.use_slit_smearing is False


# ── The Invariant integrates on top of its background ────────────────────

def test_the_invariant_replays_the_saved_background_prefit():
    """Nothing downstream refines a stale background out of an integration.

    An ordinary fit absorbs a bad starting background; the Invariant is a
    calculation, so whatever background it starts with is the one it
    integrates against. The batch path replayed the pre-fit and this one did
    not, and the two disagreed by more than two orders of magnitude.
    """
    q = np.linspace(0.01, 1.0, 150)
    intensity = 2.0 + 10.0 * np.exp(-q * q * 20.0**2 / 3.0)

    model = SimpleFitModel()
    model.set_model("Invariant")
    model.use_complex_bg = True
    model.params.update(BG_B=0.0, BG_P=4.0, BG_flat=0.0)
    model.bg_prefit = {"enabled": True,
                       "flat": {"use": True, "q_min": 0.8, "q_max": 1.0}}
    config = model.to_dict()

    reference = SimpleFitModel.from_dict(config)
    reference.prefit_background(q, intensity)
    expected = reference.fit(q, intensity)
    assert reference.params["BG_flat"] == pytest.approx(2.0, rel=0.05), \
        "the reference pre-fit did not recover the flat background"

    report = ctrl.analyze(_payload(q, intensity), {"simple_fits": config})
    assert "error" not in report, report
    assert report["config"]["params"]["BG_flat"] == pytest.approx(
        reference.params["BG_flat"], rel=1e-6)
    assert report["results"]["derived"]["Invariant"] == pytest.approx(
        expected["derived"]["Invariant"], rel=1e-6)
    assert any("Background pre-fit" in n for n in report["analyze"]["notes"])


def test_a_held_background_is_not_refit_by_the_invariant_prefit():
    """The pre-fit obeys the same held-parameter choices the fit does."""
    q = np.linspace(0.01, 1.0, 150)
    intensity = 2.0 + 10.0 * np.exp(-q * q * 20.0**2 / 3.0)

    model = SimpleFitModel()
    model.set_model("Invariant")
    model.use_complex_bg = True
    model.params.update(BG_B=0.0, BG_P=4.0, BG_flat=0.5)
    model.bg_prefit = {"enabled": True,
                       "flat": {"use": True, "q_min": 0.8, "q_max": 1.0}}

    report = ctrl.analyze(_payload(q, intensity),
                          {"simple_fits": {**model.to_dict(),
                                           "param_fixed": {"BG_flat": True}}})
    assert "error" not in report, report
    assert report["config"]["params"]["BG_flat"] == pytest.approx(0.5)


# ── A one-sided Q range is a range ───────────────────────────────────────

@pytest.mark.parametrize("bounds", [{"q_min": 1.0}, {"q_max": 1e-5}])
def test_a_one_sided_q_range_off_the_end_falls_back_instead_of_raising(bounds):
    """Half of the nonoverlap note's own bounds can legitimately be absent.

    Formatting the missing one threw ``TypeError`` out of the code path whose
    job was to explain the problem — so instead of the documented full-range
    fallback the caller got a traceback.
    """
    q, intensity, model = _guinier()
    report = ctrl.analyze(_payload(q, intensity),
                          {"simple_fits": {**model.to_dict(), **bounds}})
    assert "error" not in report, report
    assert report["analyze"]["fit_q_min"] is None
    assert report["analyze"]["fit_q_max"] is None
    assert any("does not overlap" in n for n in report["analyze"]["notes"])
