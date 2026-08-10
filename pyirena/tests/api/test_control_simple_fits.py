"""
Behaviour tests for :mod:`pyirena.api.control.simple_fits`.

The control surface is what an LLM agent drives, and an agent cannot read a
traceback — every failure has to come back as a ``{"error", "code",
"suggestion"}`` dict.  These tests therefore cover the error paths as carefully
as the happy path: a wrong model name, a parameter that does not exist in the
selected model, an empty Q range, everything held fixed.

Synthetic data is used so the numbers are known: a Guinier sphere-ish curve
with I0 and Rg we can fit back.
"""
from __future__ import annotations

import numpy as np
import pytest

import pyirena.api.control as ctrl
from pyirena.api.control.session import create_session, drop_session


@pytest.fixture
def session():
    """A session holding a clean Guinier curve (I0 = 120, Rg = 45 Å)."""
    q = np.logspace(np.log10(2e-3), np.log10(6e-2), 120)
    intensity = 120.0 * np.exp(-(q ** 2) * 45.0 ** 2 / 3.0)
    error = 0.02 * intensity
    s = create_session("/tmp/synthetic_guinier.h5", q, intensity, error,
                       label="synthetic")
    yield s
    drop_session(s.session_id)


@pytest.fixture
def guinier(session):
    """The same session with a Guinier model selected and a sensible Q range."""
    ctrl.select_simple_model(session.session_id, "Guinier")
    # Guinier applies for q*Rg < ~1.3 → q < 0.029 for Rg = 45
    ctrl.set_fit_q_range(session.session_id, 2e-3, 2.5e-2)
    return session


# ── Model listing and selection ──────────────────────────────────────────

def test_list_models_needs_no_session_and_describes_each_model():
    out = ctrl.list_simple_models()
    assert out["ok"]
    by_name = {m["name"]: m for m in out["models"]}
    assert "Guinier" in by_name and "Porod" in by_name and "Invariant" in by_name
    assert by_name["Guinier"]["params"] == ["I0", "Rg"]
    assert by_name["Guinier"]["linearizable"] is True
    assert by_name["Invariant"]["is_calculation"] is True


def test_select_model_reports_its_parameters(session):
    out = ctrl.select_simple_model(session.session_id, "Guinier")
    assert out["ok"] and out["model"] == "Guinier"
    names = [p["name"] for p in out["config"]["parameters"]]
    assert names == ["I0", "Rg"]
    assert all(p["fixed"] is False for p in out["config"]["parameters"])


def test_unknown_model_returns_an_actionable_error(session):
    out = ctrl.select_simple_model(session.session_id, "Not A Model")
    assert out["code"] == "BAD_MODEL"
    assert "list_simple_models" in out["suggestion"] or "Choose one" in out["suggestion"]


def test_calls_without_a_model_say_which_call_is_missing(session):
    out = ctrl.get_simple_parameters(session.session_id)
    assert out["code"] == "NO_SIMPLE_MODEL"
    assert "select_simple_model" in out["suggestion"]


def test_unknown_session_is_reported(session):
    assert ctrl.get_simple_config("nosuchid")["code"] == "NO_SESSION"


def test_switching_model_resets_parameters_and_fixed_state(guinier):
    sid = guinier.session_id
    ctrl.set_simple_parameter(sid, "Rg", 999.0)
    ctrl.fix_simple_parameter(sid, "I0")

    ctrl.select_simple_model(sid, "Porod")
    params = ctrl.get_simple_parameters(sid)
    assert [p["name"] for p in params["parameters"]] != ["I0", "Rg"]
    assert all(p["fixed"] is False for p in params["parameters"])


# ── Parameters ───────────────────────────────────────────────────────────

def test_set_parameter_value(guinier):
    out = ctrl.set_simple_parameter(guinier.session_id, "Rg", 50.0)
    assert out["ok"] and out["value"] == 50.0


def test_parameter_not_in_this_model_is_named_in_the_error(guinier):
    out = ctrl.set_simple_parameter(guinier.session_id, "Porod_constant", 1.0)
    assert out["code"] == "BAD_PARAM"
    assert "get_simple_parameters" in out["suggestion"]


def test_bounds_clamp_an_out_of_range_starting_value(guinier):
    sid = guinier.session_id
    ctrl.set_simple_parameter(sid, "Rg", 5000.0)
    out = ctrl.set_simple_parameter_bounds(sid, "Rg", lo=10.0, hi=100.0)
    assert out["ok"]
    # An out-of-bounds start point makes the fit fail immediately, so it is
    # pulled into range rather than left invalid.
    assert out["value"] == 100.0


def test_bounds_reject_an_inverted_range(guinier):
    out = ctrl.set_simple_parameter_bounds(guinier.session_id, "Rg",
                                           lo=100.0, hi=10.0)
    assert out["code"] == "BAD_BOUNDS"


def test_fix_and_free_round_trip(guinier):
    sid = guinier.session_id
    assert ctrl.fix_simple_parameter(sid, "I0")["fixed"] is True
    assert [p["fixed"] for p in ctrl.get_simple_parameters(sid)["parameters"]
            if p["name"] == "I0"] == [True]
    assert ctrl.free_simple_parameter(sid, "I0")["fixed"] is False
    assert [p["fixed"] for p in ctrl.get_simple_parameters(sid)["parameters"]
            if p["name"] == "I0"] == [False]


def test_reset_restores_defaults_and_frees_everything(guinier):
    sid = guinier.session_id
    ctrl.set_simple_parameter(sid, "Rg", 999.0)
    ctrl.fix_simple_parameter(sid, "Rg")
    out = ctrl.reset_simple_parameters(sid)
    params = {p["name"]: p for p in out["config"]["parameters"]}
    assert params["Rg"]["value"] != 999.0
    assert params["Rg"]["fixed"] is False


# ── Fitting ──────────────────────────────────────────────────────────────

def test_fit_recovers_the_known_parameters(guinier):
    out = ctrl.run_simple_fit(guinier.session_id)
    assert out["success"], out
    values = {p["name"]: p["value"] for p in out["parameters"]}
    assert values["Rg"] == pytest.approx(45.0, rel=0.02)
    assert values["I0"] == pytest.approx(120.0, rel=0.02)
    assert out["reduced_chi_squared"] is not None
    assert out["n_data"] > 10


def test_fit_reports_uncertainties_and_marks_fixed_parameters(guinier):
    sid = guinier.session_id
    ctrl.set_simple_parameter(sid, "I0", 120.0)
    ctrl.fix_simple_parameter(sid, "I0")
    out = ctrl.run_simple_fit(sid)
    assert out["success"]
    by_name = {p["name"]: p for p in out["parameters"]}
    assert by_name["I0"]["fixed"] is True
    assert by_name["I0"]["value"] == pytest.approx(120.0)
    assert by_name["Rg"]["fixed"] is False
    assert by_name["Rg"]["std"] is not None


def test_fitted_values_become_the_new_starting_point(guinier):
    sid = guinier.session_id
    ctrl.run_simple_fit(sid)
    params = {p["name"]: p["value"]
              for p in ctrl.get_simple_parameters(sid)["parameters"]}
    assert params["Rg"] == pytest.approx(45.0, rel=0.02)


def test_fit_with_everything_fixed_explains_itself(guinier):
    sid = guinier.session_id
    ctrl.fix_simple_parameter(sid, "I0")
    ctrl.fix_simple_parameter(sid, "Rg")
    out = ctrl.run_simple_fit(sid)
    assert out["code"] == "ALL_FIXED"
    assert "free_simple_parameter" in out["suggestion"]


def test_set_fit_q_range_outside_the_data_is_refused(guinier):
    """The shared Q-range tool rejects (and resets) an empty window itself."""
    out = ctrl.set_fit_q_range(guinier.session_id, 100.0, 200.0)
    assert out["code"] == "EMPTY_RANGE"


def test_fit_with_an_empty_q_range_explains_itself(guinier):
    """Belt and braces: a window that selects nothing must not reach the fitter."""
    guinier.fit_q_min, guinier.fit_q_max = 100.0, 200.0     # bypass the setter
    out = ctrl.run_simple_fit(guinier.session_id)
    assert out["code"] == "EMPTY_RANGE"
    assert "reset_fit_q_range" in out["suggestion"]


def test_results_are_unavailable_before_a_fit(guinier):
    out = ctrl.get_simple_results(guinier.session_id)
    assert out["code"] == "NO_FIT"


def test_results_match_the_fit(guinier):
    sid = guinier.session_id
    fit = ctrl.run_simple_fit(sid)
    res = ctrl.get_simple_results(sid)
    assert res["ok"]
    assert res["model"] == fit["model"]
    assert res["reduced_chi_squared"] == pytest.approx(fit["reduced_chi_squared"])


def test_non_finite_numbers_never_reach_json(guinier):
    """JSON has no NaN — every numeric field must be a number or None."""
    ctrl.run_simple_fit(guinier.session_id)
    res = ctrl.get_simple_results(guinier.session_id)
    for param in res["parameters"]:
        for key in ("value", "std"):
            assert param[key] is None or np.isfinite(param[key])
    for value in res["derived"].values():
        assert value is None or np.isfinite(value)


# ── Complex background ───────────────────────────────────────────────────

def test_enabling_the_complex_background_adds_its_parameters(guinier):
    sid = guinier.session_id
    out = ctrl.set_simple_background(sid, True)
    assert out["ok"] and out["use_complex_bg"] is True
    names = [p["name"] for p in out["parameters"]]
    assert {"BG_B", "BG_P", "BG_flat"} <= set(names)


def test_disabling_the_complex_background_removes_its_parameters(guinier):
    sid = guinier.session_id
    ctrl.set_simple_background(sid, True)
    out = ctrl.set_simple_background(sid, False)
    names = [p["name"] for p in out["parameters"]]
    assert "BG_B" not in names


def test_enabling_the_background_keeps_model_parameter_values(guinier):
    sid = guinier.session_id
    ctrl.set_simple_parameter(sid, "Rg", 77.0)
    out = ctrl.set_simple_background(sid, True)
    values = {p["name"]: p["value"] for p in out["parameters"]}
    assert values["Rg"] == 77.0


def test_models_with_their_own_background_reject_the_complex_one(session):
    sid = session.session_id
    ctrl.select_simple_model(sid, "Porod")
    out = ctrl.set_simple_background(sid, True)
    assert out["code"] == "NO_COMPLEX_BG"


# ── Images ───────────────────────────────────────────────────────────────

def test_fit_image_is_a_png(guinier):
    pytest.importorskip("matplotlib")
    import base64

    sid = guinier.session_id
    ctrl.run_simple_fit(sid)
    out = ctrl.get_simple_fit_image(sid, width=600, height=480, dpi=80)
    assert out["ok"]
    assert base64.b64decode(out["image_base64"])[:8] == b"\x89PNG\r\n\x1a\n"


def test_linearization_image_reports_the_straight_line_fit(guinier):
    """The linearized line comes from the *current* parameters.

    So it is only meaningful after a fit — before one it shows the default
    parameters against the data, which is exactly the visual the panel gives.
    """
    pytest.importorskip("matplotlib")
    sid = guinier.session_id
    ctrl.run_simple_fit(sid)
    out = ctrl.get_simple_linearization_image(sid, width=600, height=480, dpi=80)
    assert out["ok"]
    # Clean Guinier data linearizes almost perfectly once fitted.
    assert out["r_squared"] == pytest.approx(1.0, abs=0.01)
    assert out["image_base64"]


def test_linearization_says_so_when_the_model_has_none(session):
    sid = session.session_id
    ctrl.select_simple_model(sid, "Invariant")
    out = ctrl.get_simple_linearization_image(sid)
    assert out["code"] in ("NO_LINEARIZATION", "EMPTY_RANGE")


def test_fit_image_needs_a_fit_first(guinier):
    assert ctrl.get_simple_fit_image(guinier.session_id)["code"] == "NO_FIT"


# ── Persistence ──────────────────────────────────────────────────────────

def test_save_writes_the_results_group(guinier, tmp_path, monkeypatch):
    pytest.importorskip("h5py")
    import h5py

    from pyirena.io.nxcansas_unified import create_nxcansas_file

    target = tmp_path / "saved.h5"
    create_nxcansas_file(target, guinier.q, guinier.intensity, guinier.error,
                         sample_name="synthetic")
    guinier.file_path = str(target)

    ctrl.run_simple_fit(guinier.session_id)
    out = ctrl.save_simple_fit(guinier.session_id)
    assert out["ok"], out
    with h5py.File(target, "r") as f:
        assert "entry/simple_fit_results" in f


def test_save_needs_a_fit_first(guinier):
    assert ctrl.save_simple_fit(guinier.session_id)["code"] == "NO_FIT"
