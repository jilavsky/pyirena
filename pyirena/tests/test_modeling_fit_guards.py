"""Tests for the Modeling fit guards added after the bad1.h5 runaway-fit report.

Covers:
  * the relative restart threshold in ModelingEngine.fit() — a round that
    improves χ² by less than 1e-4·χ² must stop the restart loop instead of
    burning all five rounds (the ~10-minute fit),
  * `fit_warnings` diagnostics on ModelingResult — parameters pinned at fit
    limits, non-convergence, unphysically broad lognormal, huge reduced χ².
"""

import numpy as np

import pyirena.core.modeling as M
from pyirena.core.modeling import (
    ModelingConfig,
    SizeDistPopulation,
    _model_sanity_warnings,
    _pinned_parameter_warnings,
)


def _sphere_pop(**kw) -> SizeDistPopulation:
    pop = SizeDistPopulation()
    pop.enabled = True
    pop.form_factor = 'sphere'
    pop.ff_params = {}
    pop.dist_type = 'lognormal'
    pop.dist_params = {'min_size': 10.0, 'mean_size': 100.0, 'sdeviation': 0.2}
    pop.dist_params_fit = {'min_size': False, 'mean_size': False, 'sdeviation': False}
    pop.n_bins = 30
    for k, v in kw.items():
        setattr(pop, k, v)
    return pop


def _config(pop) -> ModelingConfig:
    cfg = ModelingConfig()
    cfg.populations = [pop]
    cfg.fit_background = False
    cfg.background = 0.0
    cfg.fit_method = 'local'
    return cfg


# ── _pinned_parameter_warnings ────────────────────────────────────────────────

def test_pinned_warning_detects_upper_and_lower_limits():
    x = np.array([0.125, 0.01, 50.0])
    lo = np.array([0.005, 0.01, 1.0])
    hi = np.array([0.125, 5.0, 1e6])
    keys = [('pop', 0, 'scale', 'scale'),
            ('pop', 0, 'dist', 'sdeviation'),
            ('pop', 0, 'dist', 'mean_size')]
    warns = _pinned_parameter_warnings(x, lo, hi, keys)
    assert len(warns) == 2
    assert any('scale' in w and 'upper' in w for w in warns)
    assert any('sdeviation' in w and 'lower' in w for w in warns)
    assert not any('mean_size' in w for w in warns)


def test_pinned_warning_no_false_positive_on_decade_spanning_limits():
    # A small value inside huge limits (0 … 1e10) is NOT at the lower limit —
    # a lo–hi span-relative test would wrongly flag it.
    warns = _pinned_parameter_warnings(
        np.array([0.26]), np.array([0.0]), np.array([1e10]), [('background',)])
    assert warns == []


# ── _model_sanity_warnings ────────────────────────────────────────────────────

def test_sanity_warning_for_extreme_lognormal_sigma():
    cfg = _config(_sphere_pop())
    cfg.populations[0].dist_params['sdeviation'] = 4.1
    warns = _model_sanity_warnings(cfg, rchi2=1.0)
    assert len(warns) == 1
    assert 'sdeviation' in warns[0]


def test_sanity_warning_for_huge_reduced_chi2():
    warns = _model_sanity_warnings(_config(_sphere_pop()), rchi2=1.2e4)
    assert any('reduced' in w for w in warns)


def test_sanity_no_warnings_for_reasonable_fit():
    assert _model_sanity_warnings(_config(_sphere_pop()), rchi2=1.2) == []


def test_sanity_ignores_disabled_populations():
    cfg = _config(_sphere_pop(enabled=False))
    cfg.populations[0].dist_params['sdeviation'] = 4.5
    assert _model_sanity_warnings(cfg, rchi2=1.0) == []


# ── Restart loop guard ────────────────────────────────────────────────────────

def test_restart_loop_stops_on_marginal_improvement(monkeypatch):
    """A round improving χ² by ≪ 1e-4·χ² must end the restart loop.

    This is the bad1.h5 failure mode: every round exhausted max_nfev in a
    flat χ² valley, improvements of a few counts on χ² ~ 10⁶ passed the old
    1e-8·χ² test, and all five restarts ran (~10 minutes). With the relative
    threshold the loop must stop after the second round.
    """
    pop = _sphere_pop()
    pop.dist_params_fit['mean_size'] = True
    cfg = _config(pop)
    q = np.linspace(0.001, 0.1, 20)
    cfg.q_min, cfg.q_max = float(q[0]), float(q[-1])

    engine = M.ModelingEngine()
    I, _, _, _ = engine.total_intensity_maybe_smeared(cfg, q, use_cache=False)
    dI = 0.05 * I

    calls = {'n': 0}
    # χ² per round: 1e6, then an improvement of only 10 (< 1e-4·1e6 = 100).
    chi2_seq = [1.0e6, 1.0e6 - 10.0, 1.0e6 - 20.0, 1.0e6 - 30.0, 1.0e6 - 40.0]

    def fake_least_squares(fun, x0, **kw):
        c = chi2_seq[min(calls['n'], len(chi2_seq) - 1)]
        calls['n'] += 1
        r = type('R', (), {})()
        r.x = np.asarray(x0, dtype=float)
        r.fun = np.full(len(q), np.sqrt(c / len(q)))
        r.status = 0          # max_nfev exhausted, never converged
        return r

    monkeypatch.setattr(M, 'least_squares', fake_least_squares)
    result = engine.fit(cfg, q, I, dI)

    assert calls['n'] == 2, (
        f"restart loop ran {calls['n']} rounds; marginal improvement should stop it at 2")
    # status 0 on the final round must be reported to the user
    assert any('without converging' in w for w in result.fit_warnings)


# ── End-to-end warning on a real (small) fit ──────────────────────────────────

def test_fit_warns_when_scale_pins_at_limit():
    """Data 10× above what the scale limit allows → scale pins, warning raised."""
    pop = _sphere_pop()
    pop.scale = 0.01
    pop.fit_scale = True
    pop.scale_limits = (0.001, 0.02)
    cfg = _config(pop)
    q = np.linspace(0.001, 0.1, 30)
    cfg.q_min, cfg.q_max = float(q[0]), float(q[-1])

    engine = M.ModelingEngine()
    I_model, _, _, _ = engine.total_intensity_maybe_smeared(cfg, q, use_cache=False)
    I_data = 10.0 * I_model                  # unreachable: needs scale = 0.1 > 0.02
    dI = 0.02 * I_data

    result = engine.fit(cfg, q, I_data, dI)
    assert any('scale' in w and 'upper fit limit' in w for w in result.fit_warnings)


def test_fit_no_warnings_on_well_posed_fit():
    pop = _sphere_pop()
    pop.scale = 0.012
    pop.fit_scale = True
    pop.scale_limits = (0.001, 1.0)
    cfg = _config(pop)
    q = np.linspace(0.001, 0.1, 30)
    cfg.q_min, cfg.q_max = float(q[0]), float(q[-1])

    engine = M.ModelingEngine()
    I_model, _, _, _ = engine.total_intensity_maybe_smeared(cfg, q, use_cache=False)
    # Data generated by the model itself at scale 0.01 — reachable, well posed.
    I_data = I_model * (0.01 / 0.012)
    dI = 0.02 * I_data

    result = engine.fit(cfg, q, I_data, dI)
    assert result.fit_warnings == []
