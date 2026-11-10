"""Monte-Carlo uncertainties on the agent surface (issues #28, #38).

``run_*_fit(n_mc_runs=…)`` for Unified Fit, Size Distribution, Simple Fits and
Modeling calls the same core implementation as the GUI and ``pyirena.batch``,
so these tests check the wiring and the contract an agent relies on — off by
default, refused without measured σ, reproducible with a seed, stored for
``get_parameter_uncertainties`` and ``export_results``, strictly JSON — not
the statistics, which ``pyirena/tests/test_mc_uncertainty.py`` covers.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

import pyirena.api.control as ctrl

N = 150


def _saxs(with_errors=True):
    rng = np.random.default_rng(3)
    q = np.logspace(-3, -0.5, N)
    I = 1000.0 * np.exp(-(q**2) * 150.0**2 / 3) + 1e-6 * q**-4.0 + 0.01
    I = I * (1 + 0.02 * rng.standard_normal(q.size))
    r = ctrl.open_dataset_from_data(
        q=q.tolist(), intensity=I.tolist(),
        error=(0.02 * np.abs(I)).tolist() if with_errors else None,
        label="mc test",
    )
    assert "error" not in r, r
    return r["session_id"]


def _unified(sid, **kw):
    ctrl.select_model(sid, "unified_fit")
    ctrl.add_unified_level(sid)
    return ctrl.run_fit(sid, **kw)


def _sizes(sid, **kw):
    ctrl.select_sizes_model(sid)
    return ctrl.run_sizes_fit(sid, **kw)


def _simple(sid, **kw):
    ctrl.select_simple_model(sid, "Guinier")
    ctrl.set_fit_q_range(sid, 1e-3, 0.01)
    return ctrl.run_simple_fit(sid, **kw)


def _modeling(sid, **kw):
    ctrl.select_modeling_model(sid)
    ctrl.add_population(sid, "size_dist", label="pores")
    return ctrl.run_modeling_fit(sid, **kw)


RUNNERS = {
    "unified_fit": _unified,
    "sizes": _sizes,
    "simple_fits": _simple,
    "modeling": _modeling,
}


@pytest.fixture
def session():
    sids = []

    def make(with_errors=True):
        sid = _saxs(with_errors)
        sids.append(sid)
        return sid

    yield make
    for sid in sids:
        ctrl.close_session(sid)


@pytest.mark.parametrize("tool", RUNNERS)
def test_off_by_default(session, tool):
    sid = session()
    r = RUNNERS[tool](sid)
    assert r.get("success"), r
    assert "uncertainty" not in r
    u = ctrl.get_parameter_uncertainties(sid)
    assert u["ok"] and u["tool"] == tool
    if tool != "simple_fits":                 # Simple Fits falls back to covariance
        assert u["available"] is False and "n_mc_runs" in u["reason"]


@pytest.mark.parametrize("tool", RUNNERS)
def test_mc_runs_and_is_stored(session, tool):
    sid = session()
    r = RUNNERS[tool](sid, n_mc_runs=4, random_seed=7)
    assert r.get("success"), r
    block = r["uncertainty"]
    assert block["method"] == "monte_carlo" and block["n_runs"] == 4
    assert block["available"], block
    u = ctrl.get_parameter_uncertainties(sid)
    assert u["available"] and u["method"] == "monte_carlo"
    if tool == "sizes":
        assert u["scalars"]["rg_std"] is not None
        assert "distribution_std" not in u
        arr = ctrl.get_parameter_uncertainties(sid, include_arrays=True)
        assert len(arr["distribution_std"]) > 0
    else:
        assert u["parameters"] and any(p["std"] for p in u["parameters"])
    report = ctrl.export_results(sid)
    assert report["uncertainty"]["method"] == "monte_carlo"
    json.dumps(report, allow_nan=False)
    json.dumps(u, allow_nan=False)


@pytest.mark.parametrize("tool", ["unified_fit", "sizes", "simple_fits"])
def test_seed_reproduces(session, tool):
    a = RUNNERS[tool](session(), n_mc_runs=3, random_seed=11)
    b = RUNNERS[tool](session(), n_mc_runs=3, random_seed=11)
    def key(r):
        if tool == "unified_fit":
            return [p.get("std") for p in r["parameters_updated"]]
        if tool == "sizes":
            return r["uncertainty"]["scalars"]
        return [p["std"] for p in r["parameters"]]

    assert key(a)
    assert key(a) == key(b)


@pytest.mark.parametrize("tool", RUNNERS)
def test_refused_without_measured_errors(session, tool):
    sid = session(with_errors=False)
    r = RUNNERS[tool](sid, n_mc_runs=3)
    assert r.get("success"), r
    assert r["uncertainty"]["available"] is False
    assert "uncertainties" in r["uncertainty"]["reason"]


@pytest.mark.parametrize("bad", [1, -2, "x"])
def test_bad_n_mc_runs(session, bad):
    sid = session()
    ctrl.select_model(sid, "unified_fit")
    ctrl.add_unified_level(sid)
    r = ctrl.run_fit(sid, n_mc_runs=bad)
    assert r["code"] == "BAD_ARGUMENTS"


def test_refit_clears_stale_uncertainty(session):
    sid = session()
    _unified(sid, n_mc_runs=3, random_seed=1)
    assert ctrl.get_parameter_uncertainties(sid)["available"]
    ctrl.run_fit(sid)
    assert ctrl.get_parameter_uncertainties(sid)["available"] is False


def test_simple_mc_replaces_covariance_std(session):
    sid = session()
    cov = _simple(sid)
    assert cov["uncertainty_method"] == "covariance"
    mc = ctrl.run_simple_fit(sid, n_mc_runs=5, random_seed=2)
    assert mc["uncertainty_method"] == "monte_carlo"
    assert ctrl.get_simple_results(sid)["parameters"]
    stds = {p["name"]: p["std"] for p in mc["parameters"]}
    assert stds == {p["name"]: p["std"]
                    for p in ctrl.get_parameter_uncertainties(sid)["parameters"]}


def test_unified_save_stores_mc_uncertainties(session, tmp_path):
    from pyirena.io.nxcansas_unified import load_unified_fit_results

    sid = session()
    r = _unified(sid, n_mc_runs=4, random_seed=3)
    rg_std = next(p["std"] for p in r["parameters_updated"] if p["name"] == "Rg_1")
    out = tmp_path / "unified_mc.h5"
    saved = ctrl.save_fit(sid, output_path=str(out))
    assert "error" not in saved, saved
    loaded = load_unified_fit_results(out)
    assert loaded["levels"][0]["Rg_err"] == pytest.approx(rg_std)


def test_analyze_passes_n_mc_runs(session):
    sid = session()
    _unified(sid)
    report = ctrl.export_results(sid)
    q = np.logspace(-3, -0.5, N)
    I = 1000.0 * np.exp(-(q**2) * 150.0**2 / 3) + 1e-6 * q**-4.0 + 0.01
    data = {"q": q.tolist(), "intensity": I.tolist(), "error": (0.02 * I).tolist()}
    out = ctrl.analyze(data, report, n_mc_runs=3)
    assert out["analyze"]["n_mc_runs"] == 3
    assert out["uncertainty"]["method"] == "monte_carlo"
    assert out["uncertainty"]["n_runs"] == 3
