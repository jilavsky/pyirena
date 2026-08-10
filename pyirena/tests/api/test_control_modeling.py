"""
Behaviour tests for :mod:`pyirena.api.control.modeling`.

Modeling is the widest control surface: six population types, each with its own
parameters, some nested in ``dist_params`` / ``ff_params`` / ``sf_params``.  The
surface flattens all of that into one dotted namespace, so the tests here are
mostly about that translation being honest:

* the parameter list matches what the fitter would actually pack;
* changing the distribution or form factor re-derives it (a core-shell form
  factor must bring its SLDs, and dropping back must remove them);
* a name that is not active is refused with a message naming the alternative.

The synthetic data is a known two-parameter size distribution, so a fit has a
right answer to recover.
"""
from __future__ import annotations

import numpy as np
import pytest

import pyirena.api.control as ctrl
from pyirena.api.control.session import create_session, drop_session


@pytest.fixture
def session():
    """A session holding a lognormal sphere population plus flat background."""
    from pyirena.core.modeling import (
        ModelingConfig,
        ModelingEngine,
        SizeDistPopulation,
    )

    q = np.logspace(np.log10(0.005), np.log10(0.3), 80)
    pop = SizeDistPopulation()
    pop.dist_params = {"min_size": 10.0, "mean_size": 90.0, "sdeviation": 0.25}
    config = ModelingConfig(populations=[pop], background=0.002,
                            q_min=float(q.min()), q_max=float(q.max()))
    # total_intensity() already includes config.background — adding it again
    # would make the generating parameters *not* a solution, and every
    # "recovers the truth" assertion below would quietly be testing nothing.
    intensity = ModelingEngine().total_intensity(config, q, use_cache=False)[0]

    s = create_session("/tmp/synthetic_modeling.h5", q, intensity,
                       0.02 * intensity, label="synthetic")
    yield s
    drop_session(s.session_id)


@pytest.fixture
def modeling(session):
    """The same session with an empty Modeling configuration."""
    ctrl.select_modeling_model(session.session_id)
    return session


@pytest.fixture
def one_pop(modeling):
    """…plus a single size-distribution population at index 0."""
    ctrl.add_population(modeling.session_id, "size_dist", label="pores")
    return modeling


# ── Lifecycle ────────────────────────────────────────────────────────────

def test_select_starts_empty_with_the_full_q_range(session):
    out = ctrl.select_modeling_model(session.session_id)
    assert out["ok"]
    assert out["config"]["n_populations"] == 0
    assert out["config"]["q_min"] == pytest.approx(float(session.q.min()))
    assert out["config"]["q_max"] == pytest.approx(float(session.q.max()))


def test_calls_without_a_model_say_which_call_is_missing(session):
    out = ctrl.list_populations(session.session_id)
    assert out["code"] == "NO_MODELING_MODEL"
    assert "select_modeling_model" in out["suggestion"]


def test_population_types_describe_their_options():
    out = ctrl.list_population_types()
    by_type = {t["pop_type"]: t for t in out["population_types"]}
    assert set(by_type) == {
        "size_dist", "unified_level", "guinier_porod",
        "diffraction_peak", "mass_fractal", "surface_fractal",
    }
    size = by_type["size_dist"]
    assert "form_factor" in size["options"]
    # The form-factor catalogue must carry each shape's own parameters.
    assert size["form_factors"]["sphere"] == []
    assert "sld_core" in size["form_factors"]["cs_sphere_by_core"]
    assert "mean_size" in size["distributions"]["lognormal"]
    assert by_type["diffraction_peak"]["peak_types"] == [
        "gaussian", "lorentzian", "voigt",
    ]


# ── Populations ──────────────────────────────────────────────────────────

def test_add_population_returns_index_and_parameters(modeling):
    out = ctrl.add_population(modeling.session_id, "size_dist", label="pores")
    assert out["ok"] and out["index"] == 0
    assert out["population"]["label"] == "pores"
    names = [p["name"] for p in out["parameters"]]
    assert "dist.mean_size" in names and "scale" in names


def test_unknown_population_type_is_refused(modeling):
    out = ctrl.add_population(modeling.session_id, "not_a_population")
    assert out["code"] == "BAD_POPULATION_TYPE"


def test_each_population_type_can_be_added_and_lists_parameters(modeling):
    sid = modeling.session_id
    expected_first = {
        "unified_level": "G",
        "guinier_porod": "G",
        "diffraction_peak": "position",
        "mass_fractal": "Phi",
        "surface_fractal": "Surface",
    }
    for pop_type, first in expected_first.items():
        out = ctrl.add_population(sid, pop_type)
        assert out["ok"], pop_type
        names = [p["name"] for p in out["parameters"]]
        assert names and names[0] == first, (pop_type, names)


def test_remove_population_reindexes(modeling):
    sid = modeling.session_id
    ctrl.add_population(sid, "size_dist", label="first")
    ctrl.add_population(sid, "unified_level", label="second")
    out = ctrl.remove_population(sid, 0)
    assert out["ok"]
    assert [p["label"] for p in out["populations"]] == ["second"]
    assert out["populations"][0]["index"] == 0


def test_bad_population_index_is_reported(modeling):
    out = ctrl.get_population_parameters(modeling.session_id, 5)
    assert out["code"] == "BAD_POPULATION"
    assert "list_populations" in out["suggestion"]


def test_disabled_population_is_excluded_from_the_count(one_pop):
    sid = one_pop.session_id
    ctrl.set_population_enabled(sid, 0, False)
    config = ctrl.get_modeling_config(sid)["config"]
    assert config["n_populations"] == 1
    assert config["n_enabled"] == 0


# ── The flat parameter namespace ─────────────────────────────────────────

def test_parameters_use_dotted_names_for_nested_values(one_pop):
    out = ctrl.get_population_parameters(one_pop.session_id, 0)
    names = [p["name"] for p in out["parameters"]]
    assert "dist.mean_size" in names          # dist_params
    assert "scale" in names                   # plain attribute
    assert "contrast" in names                # plain sphere → contrast applies


def test_only_the_active_distribution_parameters_are_listed(one_pop):
    sid = one_pop.session_id
    ctrl.set_population_option(sid, 0, "dist_type", "gauss")
    names = [p["name"] for p in
             ctrl.get_population_parameters(sid, 0)["parameters"]]
    assert "dist.mean_size" in names and "dist.width" in names
    # lognormal-only parameter must be gone
    assert "dist.sdeviation" not in names


def test_switching_form_factor_adds_and_removes_its_parameters(one_pop):
    sid = one_pop.session_id
    out = ctrl.set_population_option(sid, 0, "form_factor", "cs_sphere_by_core")
    names = [p["name"] for p in out["parameters"]]
    assert {"ff.sld_core", "ff.sld_shell", "ff.sld_solvent", "ff.t_shell"} <= set(names)
    # SLDs encode the contrast, so the separate multiplier disappears
    assert "contrast" not in names

    back = ctrl.set_population_option(sid, 0, "form_factor", "sphere")
    names = [p["name"] for p in back["parameters"]]
    assert not any(n.startswith("ff.") for n in names)
    assert "contrast" in names


def test_structure_factor_parameters_appear_only_when_active(one_pop):
    sid = one_pop.session_id
    assert not any(p["name"].startswith("sf.") for p in
                   ctrl.get_population_parameters(sid, 0)["parameters"])
    out = ctrl.set_population_option(sid, 0, "structure_factor", "hard_sphere")
    names = [p["name"] for p in out["parameters"]]
    assert "sf.radius" in names and "sf.volume_fraction" in names
    assert "sf.eta" not in names          # that belongs to 'interferences'


def test_parameter_not_active_is_refused_with_guidance(one_pop):
    out = ctrl.set_population_parameter(one_pop.session_id, 0, "ff.sld_core", 1.0)
    assert out["code"] == "BAD_PARAM"
    assert "dist.mean_size" in out["suggestion"]


def test_set_value_fit_and_bounds_round_trip(one_pop):
    sid = one_pop.session_id
    assert ctrl.set_population_parameter(sid, 0, "dist.mean_size", 90.0)["value"] == 90.0
    assert ctrl.set_population_parameter_fit(sid, 0, "dist.mean_size", True)["fit"]
    out = ctrl.set_population_parameter_bounds(sid, 0, "dist.mean_size", 20.0, 300.0)
    assert (out["lo"], out["hi"]) == (20.0, 300.0)

    row = {p["name"]: p for p in
           ctrl.get_population_parameters(sid, 0)["parameters"]}["dist.mean_size"]
    assert row["value"] == 90.0 and row["fit"] is True
    assert (row["lo"], row["hi"]) == (20.0, 300.0)


def test_bounds_clamp_an_out_of_range_value(one_pop):
    sid = one_pop.session_id
    ctrl.set_population_parameter(sid, 0, "dist.mean_size", 5000.0)
    out = ctrl.set_population_parameter_bounds(sid, 0, "dist.mean_size", 10.0, 500.0)
    assert out["value"] == 500.0


def test_bounds_reject_an_inverted_range(one_pop):
    out = ctrl.set_population_parameter_bounds(one_pop.session_id, 0,
                                               "dist.mean_size", 500.0, 10.0)
    assert out["code"] == "BAD_BOUNDS"


def test_option_not_available_for_this_type_is_refused(modeling):
    sid = modeling.session_id
    ctrl.add_population(sid, "unified_level")
    out = ctrl.set_population_option(sid, 0, "form_factor", "sphere")
    assert out["code"] == "BAD_OPTION"
    assert "correlations" in out["suggestion"]


def test_unknown_option_value_is_refused(one_pop):
    out = ctrl.set_population_option(one_pop.session_id, 0,
                                     "form_factor", "banana")
    assert out["code"] == "BAD_OPTION_VALUE"


def test_correlations_option_adds_its_parameters(modeling):
    sid = modeling.session_id
    ctrl.add_population(sid, "unified_level")
    before = [p["name"] for p in
              ctrl.get_population_parameters(sid, 0)["parameters"]]
    assert "ETA" not in before
    after = ctrl.set_population_option(sid, 0, "correlations", "true")
    assert {"ETA", "PACK"} <= {p["name"] for p in after["parameters"]}


# ── Global settings ──────────────────────────────────────────────────────

def test_background_value_and_fit_flag(one_pop):
    out = ctrl.set_modeling_background(one_pop.session_id, 0.002, True)
    assert out["background"] == 0.002 and out["fit_background"] is True


def test_q_range_is_validated_against_the_data(one_pop):
    sid = one_pop.session_id
    out = ctrl.set_modeling_q_range(sid, 0.01, 0.2)
    assert out["ok"] and out["n_points"] > 0

    empty = ctrl.set_modeling_q_range(sid, 100.0, 200.0)
    assert empty["code"] == "EMPTY_RANGE"
    assert "covers" in empty["suggestion"]

    inverted = ctrl.set_modeling_q_range(sid, 0.2, 0.01)
    assert inverted["code"] == "BAD_RANGE"


# ── Fitting ──────────────────────────────────────────────────────────────

def test_fit_from_the_truth_reproduces_the_data(one_pop):
    """Started at the generating parameters, the model matches the data.

    Proves the values an agent sets through the flat namespace really reach
    the engine: a mis-wired parameter would show up as a large chi-squared.
    """
    sid = one_pop.session_id
    ctrl.set_population_parameter(sid, 0, "dist.mean_size", 90.0)
    ctrl.set_population_parameter(sid, 0, "dist.sdeviation", 0.25)
    ctrl.set_population_parameter_fit(sid, 0, "dist.mean_size", True)
    ctrl.set_modeling_background(sid, 0.002, False)

    out = ctrl.run_modeling_fit(sid)
    assert out["success"], out
    values = {p["name"]: p["value"] for p in out["populations"][0]["parameters"]}
    assert values["dist.mean_size"] == pytest.approx(90.0, rel=0.02)
    assert out["reduced_chi_squared"] < 1.0


def test_control_surface_fits_identically_to_a_hand_built_config(session):
    """The flat namespace must be a faithful view of the dataclasses.

    Same intent expressed two ways — through these tools, and by setting the
    core objects directly — has to reach the same optimum. This is the real
    guard on the dist./ff./sf. flattening: a fit flag written to the wrong
    dict would silently change the answer rather than raise.
    """
    from pyirena.core.modeling import (
        ModelingConfig,
        ModelingEngine,
        SizeDistPopulation,
    )

    # Hand-built
    pop = SizeDistPopulation()
    pop.dist_params = {"min_size": 10.0, "mean_size": 70.0, "sdeviation": 0.3}
    pop.dist_params_fit = {"min_size": False, "mean_size": True, "sdeviation": True}
    pop.fit_scale = False
    pop.scale = 0.001
    config = ModelingConfig(populations=[pop], background=0.002,
                            fit_background=False,
                            q_min=float(session.q.min()), q_max=float(session.q.max()))
    expected = ModelingEngine().fit(config, session.q, session.intensity,
                                    0.02 * session.intensity)

    # Through the control surface
    sid = session.session_id
    ctrl.select_modeling_model(sid)
    ctrl.add_population(sid, "size_dist")
    ctrl.set_population_parameter(sid, 0, "dist.mean_size", 70.0)
    ctrl.set_population_parameter(sid, 0, "dist.sdeviation", 0.3)
    ctrl.set_population_parameter_fit(sid, 0, "dist.mean_size", True)
    ctrl.set_population_parameter_fit(sid, 0, "dist.sdeviation", True)
    ctrl.set_population_parameter(sid, 0, "scale", 0.001)
    ctrl.set_population_parameter_fit(sid, 0, "scale", False)
    ctrl.set_modeling_background(sid, 0.002, False)
    out = ctrl.run_modeling_fit(sid)

    values = {p["name"]: p["value"] for p in out["populations"][0]["parameters"]}
    assert values["dist.mean_size"] == pytest.approx(
        expected.config.populations[0].dist_params["mean_size"], rel=1e-6)
    assert values["dist.sdeviation"] == pytest.approx(
        expected.config.populations[0].dist_params["sdeviation"], rel=1e-6)
    assert out["reduced_chi_squared"] == pytest.approx(
        expected.reduced_chi_squared, rel=1e-6)


def test_fit_improves_chi_squared_from_a_poor_start(one_pop):
    """A distant start converges to *a* minimum, not necessarily the truth.

    Modeling's chi-squared surface is multimodal — which is why run_modeling_fit
    offers fit_method='global'. The tool must report the minimum it reached
    honestly rather than pretend to have found the global one.
    """
    sid = one_pop.session_id
    ctrl.set_population_parameter(sid, 0, "dist.mean_size", 40.0)
    ctrl.set_population_parameter_fit(sid, 0, "dist.mean_size", True)
    ctrl.set_population_parameter_fit(sid, 0, "dist.sdeviation", True)
    ctrl.set_modeling_background(sid, 0.002, False)

    out = ctrl.run_modeling_fit(sid)
    assert out["success"]
    values = {p["name"]: p["value"] for p in out["populations"][0]["parameters"]}
    assert values["dist.mean_size"] > 40.0        # it moved towards the data
    assert out["reduced_chi_squared"] is not None


def test_fit_reports_derived_quantities(one_pop):
    sid = one_pop.session_id
    ctrl.set_population_parameter_fit(sid, 0, "dist.mean_size", True)
    out = ctrl.run_modeling_fit(sid)
    derived = out["populations"][0]["derived"]
    assert "volume_fraction" in derived and "Rg" in derived
    assert all(v is None or np.isfinite(v) for v in derived.values())


def test_fit_without_populations_explains_itself(modeling):
    out = ctrl.run_modeling_fit(modeling.session_id)
    assert out["code"] == "NO_POPULATIONS"
    assert "add_population" in out["suggestion"]


def test_fit_with_everything_held_explains_itself(one_pop):
    sid = one_pop.session_id
    for name in ("dist.mean_size", "dist.sdeviation", "scale"):
        ctrl.set_population_parameter_fit(sid, 0, name, False)
    ctrl.set_modeling_background(sid, 0.0, False)
    out = ctrl.run_modeling_fit(sid)
    assert out["code"] == "ALL_FIXED"
    assert "set_population_parameter_fit" in out["suggestion"]


def test_unknown_fit_method_is_refused(one_pop):
    out = ctrl.run_modeling_fit(one_pop.session_id, fit_method="magic")
    assert out["code"] == "BAD_METHOD"


def test_results_are_unavailable_before_a_fit(one_pop):
    assert ctrl.get_modeling_results(one_pop.session_id)["code"] == "NO_FIT"


def test_results_match_the_fit(one_pop):
    sid = one_pop.session_id
    ctrl.set_population_parameter_fit(sid, 0, "dist.mean_size", True)
    fit = ctrl.run_modeling_fit(sid)
    res = ctrl.get_modeling_results(sid)
    assert res["ok"]
    assert res["reduced_chi_squared"] == pytest.approx(fit["reduced_chi_squared"])
    assert len(res["populations"]) == len(fit["populations"])


def test_disabled_populations_stay_out_of_the_fit(one_pop):
    sid = one_pop.session_id
    ctrl.add_population(sid, "unified_level", label="unused")
    ctrl.set_population_enabled(sid, 1, False)
    ctrl.set_population_parameter_fit(sid, 0, "dist.mean_size", True)
    out = ctrl.run_modeling_fit(sid)
    assert out["success"]
    assert [p["index"] for p in out["populations"]] == [0]


# ── Image and persistence ────────────────────────────────────────────────

def test_fit_image_is_a_png(one_pop):
    pytest.importorskip("matplotlib")
    import base64

    sid = one_pop.session_id
    ctrl.set_population_parameter_fit(sid, 0, "dist.mean_size", True)
    ctrl.run_modeling_fit(sid)
    out = ctrl.get_modeling_fit_image(sid, width=600, height=480, dpi=80)
    assert out["ok"]
    assert base64.b64decode(out["image_base64"])[:8] == b"\x89PNG\r\n\x1a\n"


def test_fit_image_needs_a_fit_first(one_pop):
    assert ctrl.get_modeling_fit_image(one_pop.session_id)["code"] == "NO_FIT"


def test_save_writes_the_results_group(one_pop, tmp_path):
    pytest.importorskip("h5py")
    import h5py

    from pyirena.io.nxcansas_unified import create_nxcansas_file

    target = tmp_path / "saved.h5"
    create_nxcansas_file(target, one_pop.q, one_pop.intensity, one_pop.error,
                         sample_name="synthetic")
    one_pop.file_path = str(target)

    sid = one_pop.session_id
    ctrl.set_population_parameter_fit(sid, 0, "dist.mean_size", True)
    ctrl.run_modeling_fit(sid)
    out = ctrl.save_modeling_fit(sid)
    assert out["ok"], out
    with h5py.File(target, "r") as f:
        assert "entry/modeling_results" in f


def test_save_needs_a_fit_first(one_pop):
    assert ctrl.save_modeling_fit(one_pop.session_id)["code"] == "NO_FIT"
