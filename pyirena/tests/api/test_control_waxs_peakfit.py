"""
Behaviour tests for :mod:`pyirena.api.control.waxs_peakfit`.

A WAXS pattern is a background with peaks on it, so this surface is background
choice plus peak-list management.  The synthetic pattern below has two known
Gaussian peaks (Q0 = 2.10 and 3.55) on a broad halo, with noise — so the peak
finder has something real to find and the fit has a right answer.

As with the other control surfaces, the error paths matter as much as the happy
path: an agent cannot read a traceback, only a ``{"error", "code",
"suggestion"}`` dict.
"""
from __future__ import annotations

import numpy as np
import pytest

import pyirena.api.control as ctrl
from pyirena.api.control.session import create_session, drop_session

PEAK_Q0 = (2.10, 3.55)
PEAK_A = (800.0, 450.0)
PEAK_FWHM = (0.06, 0.08)


@pytest.fixture
def session():
    q = np.linspace(1.0, 6.0, 600)
    intensity = 200.0 + 20.0 * np.exp(-((q - 1.0) / 2.0) ** 2)      # halo
    for q0, amplitude, fwhm in zip(PEAK_Q0, PEAK_A, PEAK_FWHM):
        intensity = intensity + amplitude * np.exp(-((q - q0) / fwhm) ** 2)
    intensity = intensity + np.random.default_rng(0).normal(0, 3.0, size=q.size)

    s = create_session("/tmp/synthetic_waxs.h5", q, intensity,
                       np.sqrt(np.abs(intensity)), label="synthetic")
    yield s
    drop_session(s.session_id)


@pytest.fixture
def waxs(session):
    ctrl.select_waxs_model(session.session_id, "SNIP")
    return session


@pytest.fixture
def two_peaks(waxs):
    """…with the two real peaks placed by hand (no spurious detections)."""
    for q0, amplitude, fwhm in zip(PEAK_Q0, PEAK_A, PEAK_FWHM):
        ctrl.add_waxs_peak(waxs.session_id, q0, "Gauss",
                           amplitude=amplitude, fwhm=fwhm)
    return waxs


# ── Options and lifecycle ────────────────────────────────────────────────

def test_options_list_shapes_and_flag_adaptive_backgrounds():
    out = ctrl.list_waxs_options()
    assert "Gauss" in out["peak_shapes"] and "Pseudo-Voigt" in out["peak_shapes"]
    by_name = {b["name"]: b for b in out["background_shapes"]}
    assert by_name["SNIP"]["adaptive"] is True
    assert by_name["Linear"]["adaptive"] is False
    assert by_name["Linear"]["parameters"]            # polynomial coefficients
    assert out["weight_modes"] == ["standard", "equal", "relative"]


def test_select_model_starts_with_no_peaks(session):
    out = ctrl.select_waxs_model(session.session_id, "SNIP")
    assert out["ok"]
    assert out["config"]["n_peaks"] == 0
    assert out["config"]["background_is_adaptive"] is True


def test_unknown_background_is_refused(session):
    out = ctrl.select_waxs_model(session.session_id, "Chebyshev")
    assert out["code"] == "BAD_BG_SHAPE"


def test_calls_without_a_model_say_which_call_is_missing(session):
    out = ctrl.list_waxs_peaks(session.session_id)
    assert out["code"] == "NO_WAXS_MODEL"
    assert "select_waxs_model" in out["suggestion"]


# ── Background ───────────────────────────────────────────────────────────

def test_switching_background_keeps_the_peaks(two_peaks):
    sid = two_peaks.session_id
    out = ctrl.set_waxs_background(sid, "Linear")
    assert out["ok"]
    assert out["config"]["background_is_adaptive"] is False
    assert ctrl.get_waxs_config(sid)["config"]["n_peaks"] == 2


def test_polynomial_background_parameters_accept_fit_and_bounds(waxs):
    sid = waxs.session_id
    ctrl.set_waxs_background(sid, "Linear")
    name = ctrl.get_waxs_config(sid)["config"]["background_parameters"][0]["name"]
    out = ctrl.set_waxs_background_parameter(sid, name, value=100.0,
                                             fit=False, lo=0.0, hi=500.0)
    assert out["ok"]
    row = {r["name"]: r for r in out["background_parameters"]}[name]
    assert row["value"] == 100.0 and row["fit"] is False
    assert (row["lo"], row["hi"]) == (0.0, 500.0)


def test_adaptive_background_takes_a_value_but_not_a_fit_flag(waxs):
    sid = waxs.session_id
    name = ctrl.get_waxs_config(sid)["config"]["background_parameters"][0]["name"]
    out = ctrl.set_waxs_background_parameter(sid, name, value=0.2, fit=True)
    assert out["ok"]
    assert "note" in out and "not" in out["note"].lower()
    row = {r["name"]: r for r in out["background_parameters"]}[name]
    assert row["value"] == 0.2
    assert "fit" not in row          # adaptive rows carry a value only


def test_unknown_background_parameter_names_the_valid_ones(waxs):
    out = ctrl.set_waxs_background_parameter(waxs.session_id, "nonsense", value=1.0)
    assert out["code"] == "BAD_BG_PARAM"


# ── Peaks ────────────────────────────────────────────────────────────────

def test_find_peaks_locates_the_real_ones(waxs):
    out = ctrl.find_waxs_peaks(waxs.session_id)
    assert out["ok"] and out["n_found"] > 0
    found_q0 = [
        {p["name"]: p["value"] for p in peak["parameters"]}["Q0"]
        for peak in out["peaks"]
    ]
    for q0 in PEAK_Q0:
        assert any(abs(f - q0) < 0.05 for f in found_q0), (q0, found_q0)


def test_find_peaks_replaces_or_appends(waxs):
    sid = waxs.session_id
    first = ctrl.find_waxs_peaks(sid)["n_found"]
    same = ctrl.find_waxs_peaks(sid, replace=True)
    assert len(same["peaks"]) == first
    more = ctrl.find_waxs_peaks(sid, replace=False)
    assert len(more["peaks"]) == 2 * first


def test_raising_prominence_finds_fewer_peaks(waxs):
    sid = waxs.session_id
    loose = ctrl.find_waxs_peaks(sid, prominence_frac=0.01)["n_found"]
    strict = ctrl.find_waxs_peaks(sid, prominence_frac=0.4)["n_found"]
    assert strict <= loose


def test_add_peak_takes_its_amplitude_from_the_data(waxs):
    out = ctrl.add_waxs_peak(waxs.session_id, PEAK_Q0[0])
    values = {p["name"]: p["value"] for p in out["peak"]["parameters"]}
    # The peak sits ~1000 counts above a ~200 background; a generic default
    # of 1.0 would be useless as a starting point.
    assert values["A"] > 500.0


def test_peak_outside_the_data_range_is_refused(waxs):
    out = ctrl.add_waxs_peak(waxs.session_id, 99.0)
    assert out["code"] == "Q0_OUT_OF_RANGE"
    assert "covers" in out["suggestion"]


def test_remove_peak_reindexes(two_peaks):
    sid = two_peaks.session_id
    out = ctrl.remove_waxs_peak(sid, 0)
    assert out["ok"] and len(out["peaks"]) == 1
    assert out["peaks"][0]["index"] == 0


def test_bad_peak_index_is_reported(waxs):
    out = ctrl.get_waxs_peak_parameters(waxs.session_id, 3)
    assert out["code"] == "BAD_PEAK"
    assert "list_waxs_peaks" in out["suggestion"]


def test_peaks_report_a_derived_area(two_peaks):
    peaks = ctrl.list_waxs_peaks(two_peaks.session_id)["peaks"]
    assert all(p["area"] is not None and p["area"] > 0 for p in peaks)


# ── Peak parameters ──────────────────────────────────────────────────────

def test_set_value_fit_and_bounds_round_trip(two_peaks):
    sid = two_peaks.session_id
    assert ctrl.set_waxs_peak_parameter(sid, 0, "Q0", 2.11)["ok"]
    assert ctrl.set_waxs_peak_parameter_fit(sid, 0, "Q0", False)["fit"] is False
    out = ctrl.set_waxs_peak_parameter_bounds(sid, 0, "FWHM", 0.01, 0.2)

    row = {p["name"]: p for p in out["peak"]["parameters"]}
    assert row["Q0"]["value"] == 2.11 and row["Q0"]["fit"] is False
    assert (row["FWHM"]["lo"], row["FWHM"]["hi"]) == (0.01, 0.2)


def test_eta_exists_only_on_a_pseudo_voigt(two_peaks):
    sid = two_peaks.session_id
    out = ctrl.set_waxs_peak_parameter(sid, 0, "eta", 0.5)
    assert out["code"] == "BAD_PARAM"
    assert "Pseudo-Voigt" in out["suggestion"]

    changed = ctrl.set_waxs_peak_shape(sid, 0, "Pseudo-Voigt")
    assert "eta" in {p["name"] for p in changed["peak"]["parameters"]}
    assert ctrl.set_waxs_peak_parameter(sid, 0, "eta", 0.5)["ok"]


def test_changing_shape_keeps_position_and_width(two_peaks):
    sid = two_peaks.session_id
    before = {p["name"]: p["value"]
              for p in ctrl.get_waxs_peak_parameters(sid, 0)["peak"]["parameters"]}
    after = {p["name"]: p["value"]
             for p in ctrl.set_waxs_peak_shape(sid, 0, "Lorentz")["peak"]["parameters"]}
    assert after["Q0"] == before["Q0"]
    assert after["FWHM"] == before["FWHM"]


def test_unknown_peak_shape_is_refused(two_peaks):
    out = ctrl.set_waxs_peak_shape(two_peaks.session_id, 0, "Triangle")
    assert out["code"] == "BAD_PEAK_SHAPE"


def test_bounds_reject_an_inverted_range(two_peaks):
    out = ctrl.set_waxs_peak_parameter_bounds(two_peaks.session_id, 0,
                                              "FWHM", 1.0, 0.1)
    assert out["code"] == "BAD_BOUNDS"


# ── Fitting ──────────────────────────────────────────────────────────────

def test_fit_recovers_the_known_peaks(two_peaks):
    out = ctrl.run_waxs_fit(two_peaks.session_id)
    assert out["success"], out
    positions = sorted(
        {p["name"]: p["value"] for p in peak["parameters"]}["Q0"]
        for peak in out["peaks"]
    )
    for fitted, expected in zip(positions, PEAK_Q0):
        assert fitted == pytest.approx(expected, abs=0.01)


def test_fit_reports_uncertainties_and_areas(two_peaks):
    out = ctrl.run_waxs_fit(two_peaks.session_id)
    peak = out["peaks"][0]
    stds = [p["std"] for p in peak["parameters"] if p["fit"]]
    assert stds and all(s is None or np.isfinite(s) for s in stds)
    assert peak["area"] is not None and peak["area"] > 0


def test_held_position_stays_put(two_peaks):
    sid = two_peaks.session_id
    ctrl.set_waxs_peak_parameter(sid, 0, "Q0", 2.05)     # deliberately off
    ctrl.set_waxs_peak_parameter_fit(sid, 0, "Q0", False)
    out = ctrl.run_waxs_fit(sid)
    held = {p["name"]: p for p in out["peaks"][0]["parameters"]}["Q0"]
    assert held["value"] == pytest.approx(2.05)
    assert held["fit"] is False


def test_fit_without_peaks_explains_itself(waxs):
    out = ctrl.run_waxs_fit(waxs.session_id)
    assert out["code"] == "NO_PEAKS"
    assert "find_waxs_peaks" in out["suggestion"]


def test_unknown_weight_mode_is_refused(two_peaks):
    out = ctrl.run_waxs_fit(two_peaks.session_id, weight_mode="inverse_vibes")
    assert out["code"] == "BAD_WEIGHT_MODE"


def test_every_weight_mode_runs(two_peaks):
    for mode in ("standard", "equal", "relative"):
        out = ctrl.run_waxs_fit(two_peaks.session_id, weight_mode=mode)
        assert out["success"], (mode, out)


def test_fit_with_an_empty_q_range_explains_itself(two_peaks):
    two_peaks.fit_q_min, two_peaks.fit_q_max = 100.0, 200.0
    out = ctrl.run_waxs_fit(two_peaks.session_id)
    assert out["code"] == "EMPTY_RANGE"


def test_results_are_unavailable_before_a_fit(two_peaks):
    assert ctrl.get_waxs_results(two_peaks.session_id)["code"] == "NO_FIT"


def test_results_match_the_fit(two_peaks):
    sid = two_peaks.session_id
    fit = ctrl.run_waxs_fit(sid)
    res = ctrl.get_waxs_results(sid)
    assert res["ok"]
    assert res["reduced_chi_squared"] == pytest.approx(fit["reduced_chi_squared"])
    assert len(res["peaks"]) == len(fit["peaks"])


def test_non_finite_numbers_never_reach_json(two_peaks):
    ctrl.run_waxs_fit(two_peaks.session_id)
    res = ctrl.get_waxs_results(two_peaks.session_id)
    for peak in res["peaks"]:
        for param in peak["parameters"]:
            for key in ("value", "std", "lo", "hi"):
                assert param[key] is None or np.isfinite(param[key])


# ── Image and persistence ────────────────────────────────────────────────

def test_fit_image_is_a_png(two_peaks):
    pytest.importorskip("matplotlib")
    import base64

    sid = two_peaks.session_id
    ctrl.run_waxs_fit(sid)
    out = ctrl.get_waxs_fit_image(sid, width=600, height=480, dpi=80)
    assert out["ok"]
    assert base64.b64decode(out["image_base64"])[:8] == b"\x89PNG\r\n\x1a\n"


def test_fit_image_needs_a_fit_first(two_peaks):
    assert ctrl.get_waxs_fit_image(two_peaks.session_id)["code"] == "NO_FIT"


def test_save_writes_the_results_group(two_peaks, tmp_path):
    pytest.importorskip("h5py")
    import h5py

    from pyirena.io.nxcansas_unified import create_nxcansas_file

    target = tmp_path / "saved.h5"
    create_nxcansas_file(target, two_peaks.q, two_peaks.intensity,
                         two_peaks.error, sample_name="synthetic")
    two_peaks.file_path = str(target)

    ctrl.run_waxs_fit(two_peaks.session_id)
    out = ctrl.save_waxs_fit(two_peaks.session_id)
    assert out["ok"], out
    with h5py.File(target, "r") as f:
        assert "entry/waxs_peakfit_results" in f


def test_save_needs_a_fit_first(two_peaks):
    assert ctrl.save_waxs_fit(two_peaks.session_id)["code"] == "NO_FIT"
