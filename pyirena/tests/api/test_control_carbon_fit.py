"""
Behaviour tests for :mod:`pyirena.api.control.carbon_fit`.

This surface addresses parameters by the model's own dotted key rather than
through one setter per control, so what has to be honest is the *discovery*
loop an agent depends on: configure the shape, list the keys, set them, fit.
The tests below pin that loop, and in particular that the key list really does
track the shape — an agent that trusted a stale list would silently set
nothing.

The synthetic data comes from the model's own defaults, so a fit has a right
answer to recover.
"""
from __future__ import annotations

import numpy as np
import pytest

import pyirena.api.control as ctrl
from pyirena.api.control.session import create_session, drop_session


@pytest.fixture
def session():
    """A session holding a noise-free curve from the default Carbon model."""
    from pyirena.core.carbon_fit import CarbonFitModel

    q = np.logspace(-3, 0.65, 400)
    intensity = CarbonFitModel().evaluate(q)
    s = create_session("/tmp/synthetic_carbon.h5", q, intensity,
                       0.02 * intensity, label="synthetic carbon")
    yield s
    drop_session(s.session_id)


@pytest.fixture
def sid(session):
    ctrl.select_carbon_model(session.session_id)
    return session.session_id


# ── Discovery ───────────────────────────────────────────────────────────────

def test_options_need_no_session():
    options = ctrl.list_carbon_options()
    assert options["ok"]
    assert set(options["saxs_modes"]) == {"fractal", "teubner_strey"}
    assert set(options["waxs_envelopes"]) == {"none", "crumpled"}


def test_calls_before_select_explain_themselves(session):
    for result in (ctrl.list_carbon_parameters(session.session_id),
                   ctrl.run_carbon_fit(session.session_id),
                   ctrl.get_carbon_results(session.session_id)):
        assert result["code"] == "NO_CARBON_MODEL"
        assert "select_carbon_model" in result["suggestion"]


def test_unknown_session_is_reported(sid):
    assert ctrl.list_carbon_parameters("nope")["code"] == "NO_SESSION"


def test_select_reports_the_default_shape_and_peaks(sid):
    config = ctrl.get_carbon_config(sid)
    assert config["saxs_mode"] == "fractal"
    assert config["waxs_envelope"] == "none"
    assert [p["label"] for p in config["peaks"]] == ["002", "100", "004"]
    # The contrast chain is live from the moment the model exists.
    assert config["material"]["rho_struc"] == pytest.approx(2.26, rel=1e-6)
    assert config["material"]["contrast_micropore"] > 0


# ── The key list tracks the model's shape ───────────────────────────────────

def test_parameter_keys_follow_the_saxs_mode(sid):
    keys = {p["key"] for p in ctrl.list_carbon_parameters(sid)["parameters"]}
    assert "saxs.phi" in keys and "saxs.ts_C1" not in keys

    ctrl.configure_carbon_model(sid, saxs_mode="teubner_strey")
    keys = {p["key"] for p in ctrl.list_carbon_parameters(sid)["parameters"]}
    assert "saxs.ts_C1" in keys and "saxs.phi" not in keys


def test_parameter_keys_follow_the_switches(sid):
    keys = {p["key"] for p in ctrl.list_carbon_parameters(sid)["parameters"]}
    assert "background.S_rough" not in keys
    ctrl.configure_carbon_model(sid, use_roughness=True)
    keys = {p["key"] for p in ctrl.list_carbon_parameters(sid)["parameters"]}
    assert {"background.S_rough", "background.R_rough"} <= keys


def test_links_remove_parameters_from_the_list(sid):
    ctrl.configure_carbon_model(sid, waxs_envelope="crumpled")
    keys = {p["key"] for p in ctrl.list_carbon_parameters(sid)["parameters"]}
    assert {"waxs.R_layer", "waxs.fractal_D", "waxs.fractal_sigma"} <= keys

    ctrl.configure_carbon_model(sid, link_R_to_pore=True, link_D_to_saxs=True,
                                link_sigma_to_saxs=True)
    keys = {p["key"] for p in ctrl.list_carbon_parameters(sid)["parameters"]}
    assert not ({"waxs.R_layer", "waxs.fractal_D", "waxs.fractal_sigma"} & keys)


def test_inactive_parameters_are_visible_on_request(sid):
    active = {p["key"] for p in ctrl.list_carbon_parameters(sid)["parameters"]}
    everything = {p["key"] for p in
                  ctrl.list_carbon_parameters(sid, active_only=False)["parameters"]}
    assert "saxs.ts_C1" in everything - active


def test_bad_mode_names_the_alternatives(sid):
    result = ctrl.configure_carbon_model(sid, saxs_mode="nonsense")
    assert result["code"] == "BAD_SAXS_MODE"
    assert "teubner_strey" in result["suggestion"]


# ── Parameters ──────────────────────────────────────────────────────────────

def test_set_value_fit_and_bounds(sid):
    assert ctrl.set_carbon_parameter(sid, "saxs.pore_radius", 8.5)["value"] == 8.5
    assert ctrl.set_carbon_parameter_fit(sid, "saxs.globule_k", True)["fit"] is True
    bounds = ctrl.set_carbon_parameter_bounds(sid, "saxs.phi", lo=0.01, hi=0.6)
    assert (bounds["lo"], bounds["hi"]) == (0.01, 0.6)

    row = next(p for p in ctrl.list_carbon_parameters(sid)["parameters"]
               if p["key"] == "saxs.phi")
    assert (row["lo"], row["hi"]) == (0.01, 0.6)


def test_an_inactive_key_is_refused_with_the_live_list(sid):
    result = ctrl.set_carbon_parameter(sid, "saxs.ts_C1", 1.0)
    assert result["code"] == "BAD_PARAM"
    assert "list_carbon_parameters" in result["suggestion"]
    assert "saxs.phi" in result["suggestion"]


def test_reversed_bounds_are_refused(sid):
    assert ctrl.set_carbon_parameter_bounds(
        sid, "saxs.phi", lo=0.9, hi=0.1)["code"] == "BAD_BOUNDS"


def test_setting_a_value_reports_whether_it_is_inside_the_bounds(sid):
    ctrl.set_carbon_parameter_bounds(sid, "saxs.pore_radius", lo=1.0, hi=10.0)
    assert ctrl.set_carbon_parameter(sid, "saxs.pore_radius", 5.0)["within_bounds"]


# ── Peaks and the material chain ────────────────────────────────────────────

def test_add_and_remove_peaks(sid):
    before = len(ctrl.list_carbon_peaks(sid)["peaks"])
    result = ctrl.add_carbon_peak(sid, label="110", Q0=5.1)
    assert len(result["peaks"]) == before + 1
    assert result["peaks"][-1]["d_spacing"] == pytest.approx(2 * np.pi / 5.1)

    assert ctrl.add_carbon_peak(sid, label="110", Q0=5.2)["code"] == "DUPLICATE_PEAK"
    assert len(ctrl.remove_carbon_peak(sid, before)["peaks"]) == before
    assert ctrl.remove_carbon_peak(sid, 99)["code"] == "BAD_PEAK"


def test_moving_the_002_peak_moves_the_whole_contrast_chain(sid):
    before = ctrl.get_carbon_config(sid)["material"]
    ctrl.set_carbon_parameter(sid, "peak.002.Q0", 2 * np.pi / 3.80)
    after = ctrl.get_carbon_config(sid)["material"]
    assert after["d002"] == pytest.approx(3.80)
    assert after["rho_struc"] < before["rho_struc"]
    assert after["contrast_micropore"] < before["contrast_micropore"]


def test_material_overrides_take_effect_immediately(sid):
    result = ctrl.set_carbon_material(sid, contrast_mode="manual",
                                      contrast_porod=150.0,
                                      contrast_micropore=250.0)
    assert result["material"]["contrast_porod"] == 150.0
    assert result["material"]["contrast_micropore"] == 250.0

    bad = ctrl.set_carbon_material(sid, porosity_mode="sometimes")
    assert bad["code"] == "BAD_MODE"


def test_a_doped_formula_changes_the_sld(sid):
    plain = ctrl.get_carbon_config(sid)["material"]["sld_struc"]
    doped = ctrl.set_carbon_material(sid, formula="C0.9N0.1")["material"]
    assert doped["sld_struc"] != pytest.approx(plain)


# ── Fitting ─────────────────────────────────────────────────────────────────

def test_fit_recovers_the_generating_parameters(sid):
    truth = {p["key"]: p["value"]
             for p in ctrl.list_carbon_parameters(sid)["parameters"]}
    ctrl.set_carbon_parameter(sid, "saxs.phi", 0.2)
    ctrl.set_carbon_parameter(sid, "background.S_macro", 3.0e4)

    result = ctrl.run_carbon_fit(sid)
    assert result["ok"] and result["success"]
    assert result["parameters"]["saxs.phi"]["value"] == pytest.approx(
        truth["saxs.phi"], rel=1e-3)
    assert result["parameters"]["background.S_macro"]["value"] == pytest.approx(
        truth["background.S_macro"], rel=1e-3)


def test_results_carry_the_derived_layer(sid):
    ctrl.run_carbon_fit(sid)
    derived = ctrl.get_carbon_results(sid)["derived"]
    for key in ("S_part_m2_g", "S_mp_m2_g", "rho_struc", "rho_sample",
                "d002", "d100", "L_c", "N_layers", "L_a"):
        assert derived[key] is not None, key
    assert derived["S_mp_m2_g"] > 0


def test_results_are_json_safe(sid):
    """NaN is not JSON; a not-applicable derived value must come back as null."""
    import json

    ctrl.run_carbon_fit(sid)
    payload = json.dumps(ctrl.get_carbon_results(sid))
    assert "NaN" not in payload and "Infinity" not in payload
    # The Teubner-Strey quantities are undefined in fractal mode.
    assert ctrl.get_carbon_results(sid)["derived"]["ts_xi"] is None


def test_fit_with_nothing_ticked_says_so(sid):
    for row in ctrl.list_carbon_parameters(sid)["parameters"]:
        ctrl.set_carbon_parameter_fit(sid, row["key"], False)
    assert ctrl.run_carbon_fit(sid)["code"] == "NOTHING_TO_FIT"


def test_bad_weighting_is_refused(sid):
    assert ctrl.run_carbon_fit(sid, weighting="vibes")["code"] == "BAD_WEIGHTING"


def test_results_before_a_fit(sid):
    assert ctrl.get_carbon_results(sid)["code"] == "NO_FIT"
    assert ctrl.get_carbon_fit_image(sid)["code"] == "NO_FIT"


def test_fit_respects_the_shared_q_range(sid):
    ctrl.set_fit_q_range(sid, 0.01, 1.0)
    result = ctrl.run_carbon_fit(sid)
    assert result["n_points"] < 400


def test_monte_carlo_produces_uncertainties(sid):
    for row in ctrl.list_carbon_parameters(sid)["parameters"]:
        ctrl.set_carbon_parameter_fit(sid, row["key"], False)
    ctrl.set_carbon_parameter_fit(sid, "saxs.phi", True)
    result = ctrl.run_carbon_fit(sid, n_mc_runs=3)
    assert result["parameters"]["saxs.phi"]["std"] is not None


# ── Image and persistence ───────────────────────────────────────────────────

def test_fit_image_renders(sid, tmp_path, monkeypatch):
    monkeypatch.setenv("PYIRENA_PLOT_CACHE", str(tmp_path))
    ctrl.run_carbon_fit(sid)
    image = ctrl.get_carbon_fit_image(sid)
    assert image["ok"] and image["image_base64"]
    from pathlib import Path
    assert Path(image["image_path"]).exists()


def test_save_writes_the_group_and_the_setup(tmp_path, monkeypatch):
    from pyirena.io.nxcansas_carbon_fit import (
        load_carbon_fit_model,
        load_carbon_fit_results,
    )
    from pyirena.io.nxcansas_unified import create_nxcansas_file

    monkeypatch.setenv("PYIRENA_DATA_ROOT", str(tmp_path))
    from pyirena.core.carbon_fit import CarbonFitModel

    q = np.logspace(-3, 0.65, 300)
    intensity = CarbonFitModel().evaluate(q)
    path = tmp_path / "carbon.h5"
    create_nxcansas_file(path, q, intensity, error=0.02 * intensity,
                         sample_name="carbon")

    s = create_session(str(path), q, intensity, 0.02 * intensity)
    try:
        ctrl.select_carbon_model(s.session_id)
        ctrl.configure_carbon_model(s.session_id, use_roughness=True)
        ctrl.run_carbon_fit(s.session_id)
        saved = ctrl.save_carbon_fit(s.session_id)
        assert saved["ok"]
        assert saved["group"] == "entry/carbon_fit_results"

        back = load_carbon_fit_results(path)
        assert back["params"]
        # The embedded setup must reopen with the shape the agent chose.
        restored = load_carbon_fit_model(path)
        assert restored.background.use_roughness is True
    finally:
        drop_session(s.session_id)
