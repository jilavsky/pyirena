"""Tests for ``export_results`` — the JSON counterpart to ``save_*``.

Three things have to hold for every one of the six fitting tools, because the
ZMQ service (``planning/zmq-service/``) sends this straight down a socket:

* **One envelope.** Same top-level keys whichever tool ran, so an agent learns
  the shape once.
* **Strictly serialisable.** ``json.dumps(..., allow_nan=False)`` must succeed.
  Python's own ``json`` emits bare ``NaN``/``Infinity`` and reads them back
  happily, so a pyIrena-to-pyIrena round trip hides the bug and a strict
  parser on the other end rejects the whole reply.
* **Replayable config.** ``config`` is the model's ``to_dict()``, which is what
  lets a result be re-applied to the next measurement.

The fits here are deliberately small (200 points) — this file tests the export
envelope, not fit quality.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

import pyirena.api.control as ctrl

N = 200


def _saxs():
    rng = np.random.default_rng(3)
    q = np.logspace(-3, 0, N)
    I = 1000.0 * np.exp(-(q**2) * 150.0**2 / 3) + 1e-6 * q**-4.0 + 0.01
    return q, I * (1 + 0.02 * rng.standard_normal(q.size))


def _waxs():
    rng = np.random.default_rng(4)
    q = np.linspace(0.5, 5.0, N)
    I = 10.0 + 2.0 / q
    for c, a, w in ((1.8, 50.0, 0.12), (3.1, 25.0, 0.15)):
        I = I + a * np.exp(-((q - c) ** 2) / (2 * w**2))
    return q, I * (1 + 0.01 * rng.standard_normal(q.size))


def _carbon():
    rng = np.random.default_rng(5)
    q = np.logspace(-3, 0.7, N)
    I = 1e-4 * q**-4 + 5.0 * np.exp(-(q**2) * 5.0**2 / 3) + 0.5
    return q, I * (1 + 0.02 * rng.standard_normal(q.size))


def _open(data):
    q, I = data
    r = ctrl.open_dataset_from_data(
        q=q.tolist(), intensity=I.tolist(), error=(0.02 * np.abs(I)).tolist(),
        label="export test",
    )
    assert "error" not in r, r
    return r["session_id"]


def _fit_unified(sid):
    ctrl.select_model(sid, "unified_fit")
    ctrl.add_unified_level(sid)
    ctrl.run_fit(sid)


def _fit_sizes(sid):
    ctrl.select_sizes_model(sid)
    ctrl.run_sizes_fit(sid)


def _fit_simple(sid):
    ctrl.select_simple_model(sid, "Guinier")
    ctrl.run_simple_fit(sid)


def _fit_modeling(sid):
    ctrl.select_modeling_model(sid)
    ctrl.add_population(sid, "size_dist", label="pores")
    ctrl.run_modeling_fit(sid)


def _fit_waxs(sid):
    ctrl.select_waxs_model(sid)
    ctrl.find_waxs_peaks(sid)
    ctrl.run_waxs_fit(sid)


def _fit_carbon(sid):
    ctrl.select_carbon_model(sid)
    ctrl.run_carbon_fit(sid)


#: tool name as reported by export_results → (data builder, fit driver)
TOOLS = {
    "unified_fit":  (_saxs,   _fit_unified),
    "sizes":        (_saxs,   _fit_sizes),
    "simple_fits":  (_saxs,   _fit_simple),
    "modeling":     (_saxs,   _fit_modeling),
    "waxs_peakfit": (_waxs,   _fit_waxs),
    "carbon_fit":   (_carbon, _fit_carbon),
}


@pytest.fixture(scope="module")
def fitted_sessions():
    """One fitted session per tool, shared across this module (fits are slow)."""
    sessions = {}
    for tool, (build, fit) in TOOLS.items():
        sid = _open(build())
        fit(sid)
        sessions[tool] = sid
    yield sessions
    for sid in sessions.values():
        ctrl.close_session(sid)


@pytest.mark.parametrize("tool", sorted(TOOLS))
def test_envelope_is_the_same_shape_for_every_tool(fitted_sessions, tool):
    report = ctrl.export_results(fitted_sessions[tool])
    assert "error" not in report, report

    assert report["ok"] is True
    assert report["tool"] == tool
    assert report["pyirena_version"]
    assert report["exported_at"].startswith("20")
    assert set(report) >= {
        "ok", "tool", "pyirena_version", "exported_at",
        "data", "fit_q_range", "quality", "results", "config",
    }
    assert set(report["quality"]) == {
        "chi_squared", "reduced_chi_squared", "dof", "n_points_fitted",
        "n_parameters", "success", "message", "metrics",
    }
    assert report["data"]["n_points"] == N
    assert report["data"]["file"] is None        # opened from arrays
    assert report["data"]["has_errors"] is True
    assert report["results"]


@pytest.mark.parametrize("tool", sorted(TOOLS))
def test_every_report_is_strict_json(fitted_sessions, tool):
    """allow_nan=False is the check that matters — see the module docstring."""
    for kwargs in ({}, {"include_arrays": True}):
        report = ctrl.export_results(fitted_sessions[tool], **kwargs)
        text = json.dumps(report, allow_nan=False)
        assert "NaN" not in text and "Infinity" not in text
        # Round-trips through a strict parser unchanged.
        assert json.loads(text)["tool"] == tool


@pytest.mark.parametrize("tool", sorted(TOOLS))
def test_config_is_the_models_to_dict(fitted_sessions, tool):
    from pyirena.api.control.session import get_session

    report = ctrl.export_results(fitted_sessions[tool])
    config = report["config"]
    assert isinstance(config, dict) and config, f"{tool}: empty config"

    model = get_session(fitted_sessions[tool]).model
    assert set(config) == set(model.to_dict())


@pytest.mark.parametrize("tool", sorted(TOOLS))
def test_arrays_are_off_by_default_and_aligned_when_on(fitted_sessions, tool):
    sid = fitted_sessions[tool]
    assert "arrays" not in ctrl.export_results(sid)

    arrays = ctrl.export_results(sid, include_arrays=True)["arrays"]
    n = arrays["n_points"]
    assert n > 0
    assert len(arrays["q"]) == n
    assert "intensity_model" in arrays
    for key, value in arrays.items():
        if isinstance(value, list) and value and isinstance(value[0], (int, float)):
            assert len(value) == n, f"{tool}: '{key}' is {len(value)} long, q is {n}"


def test_arrays_are_decimated_to_max_points(fitted_sessions):
    sid = fitted_sessions["unified_fit"]

    full = ctrl.export_results(sid, include_arrays=True, max_points=None)["arrays"]
    assert full["decimated"] is False

    small = ctrl.export_results(sid, include_arrays=True, max_points=50)["arrays"]
    assert small["decimated"] is True
    assert small["n_points"] <= 50
    assert small["decimation_stride"] > 1
    assert len(small["intensity_model"]) == small["n_points"]
    # Decimation must keep the Q axis and its companions in step.
    assert small["q"][0] == pytest.approx(full["q"][0])


def test_sizes_reports_the_distribution_without_include_arrays(fitted_sessions):
    """For a size distribution the histogram is the result, not an extra."""
    report = ctrl.export_results(fitted_sessions["sizes"])
    dist = report["results"]["distribution"]
    assert len(dist["r_grid"]) == len(dist["distribution"]) > 0


def test_modeling_reports_per_population_curves(fitted_sessions):
    arrays = ctrl.export_results(
        fitted_sessions["modeling"], include_arrays=True
    )["arrays"]
    pops = arrays["populations"]
    assert pops and len(pops[0]["intensity"]) == arrays["n_points"]


def test_carbon_reports_its_three_components(fitted_sessions):
    arrays = ctrl.export_results(
        fitted_sessions["carbon_fit"], include_arrays=True
    )["arrays"]
    assert {"intensity_porod", "intensity_micropore", "intensity_waxs"} <= set(arrays)


def test_uncertainties_travel_where_the_fit_computes_them(fitted_sessions):
    """Five tools compute a std per parameter; Unified Fit has none to report."""
    simple = ctrl.export_results(fitted_sessions["simple_fits"])
    assert all("std" in p for p in simple["results"]["parameters"])

    unified = ctrl.export_results(fitted_sessions["unified_fit"])
    assert unified["results"]["uncertainties_available"] is False


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

def test_unknown_session():
    assert ctrl.export_results("nope")["code"] == "NO_SESSION"


def test_session_with_no_model():
    sid = _open(_saxs())
    assert ctrl.export_results(sid)["code"] == "NO_MODEL"
    ctrl.close_session(sid)


def test_session_with_no_fit():
    sid = _open(_saxs())
    ctrl.select_model(sid, "unified_fit")
    assert ctrl.export_results(sid)["code"] == "NO_FIT"
    ctrl.close_session(sid)
