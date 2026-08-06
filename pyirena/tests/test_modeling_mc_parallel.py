"""
Unit tests for parallel Monte-Carlo uncertainty estimation in the Modeling tool
(``ModelingEngine.calculate_uncertainty_mc``), plus its config/state plumbing.

The central guarantee under test is that parallel and serial runs produce
*identical* uncertainties for a given seed. That holds because every
perturbation is drawn in the parent process before any pass starts, which also
rules out the failure mode where workers inherit one RNG state and return N
identical passes — silently collapsing every uncertainty towards zero.
"""

import json
import os
import pickle

import numpy as np
import pytest

import pyirena.core.modeling as M
from pyirena.core.modeling import (
    ModelingConfig,
    ModelingEngine,
    SizeDistPopulation,
    _MCRun,
    _resolve_mc_workers,
)


# ── Helpers ──────────────────────────────────────────────────────────────────

def _make_fit_config(start_r=60.0):
    """One lognormal sphere population with three free parameters + background."""
    p = SizeDistPopulation()
    p.enabled = True
    p.dist_type = 'lognormal'
    p.dist_params = {'min_size': 10.0, 'mean_size': start_r, 'sdeviation': 0.3}
    p.dist_params_fit = {'min_size': False, 'mean_size': True, 'sdeviation': True}
    p.dist_params_limits = {
        'min_size': (0.0, 1e6), 'mean_size': (1.0, 1e6), 'sdeviation': (0.01, 5.0),
    }
    p.form_factor = 'sphere'
    p.scale = 0.001
    p.fit_scale = True
    p.contrast = 1.0
    p.n_bins = 60                      # small grid — these tests are about plumbing
    return ModelingConfig(
        populations=[p], background=0.0, fit_background=True,
        q_min=1e-3, q_max=0.1,
    )


def _make_dataset(cfg, n_q=60, seed=1):
    eng = ModelingEngine()
    q = np.logspace(-3, -1, n_q)
    I_true, *_ = eng.total_intensity(cfg, q, use_cache=False)
    rng = np.random.default_rng(seed)
    dI = 0.03 * I_true + 1e-12
    return q, I_true + dI * rng.standard_normal(n_q), dI


@pytest.fixture
def mc_case():
    cfg = _make_fit_config()
    q, I_obs, dI = _make_dataset(cfg)
    return ModelingEngine(), cfg, q, I_obs, dI


# ── Config / worker-count defaults ───────────────────────────────────────────

class TestConfigDefaults:
    def test_mc_workers_default_is_auto(self):
        assert ModelingConfig().mc_workers == 0

    def test_resolve_auto_leaves_headroom(self):
        n_cpu = os.cpu_count() or 1
        assert _resolve_mc_workers(0) == max(1, n_cpu - 2)
        # None and negatives mean "auto" too.
        assert _resolve_mc_workers(None) == max(1, n_cpu - 2)
        assert _resolve_mc_workers(-4) == max(1, n_cpu - 2)

    def test_resolve_explicit_counts(self):
        n_cpu = os.cpu_count() or 1
        assert _resolve_mc_workers(1) == 1
        assert _resolve_mc_workers(2) == min(2, n_cpu)
        # Never hand out more workers than there are cores.
        assert _resolve_mc_workers(9999) == n_cpu


# ── The picklable pass ───────────────────────────────────────────────────────

class TestMCRun:
    def test_pickles_and_matches_an_inline_fit(self, mc_case):
        """_MCRun must survive the worker boundary and return exactly what an
        in-process fit of the same perturbed data would."""
        eng, cfg, q, I_obs, dI = mc_case
        _, _, _, keys = eng._pack_params(M.deepcopy(cfg))

        cfg_mc = M.deepcopy(cfg)
        cfg_mc.fit_method = 'local'
        cfg_mc.de_workers = 1
        run = _MCRun(cfg_mc, q, dI, len(keys))

        I_pert = I_obs + dI * np.random.default_rng(7).standard_normal(len(q))

        expected = run(I_pert)
        assert expected is not None and len(expected) == len(keys)

        revived = pickle.loads(pickle.dumps(run))
        assert revived(I_pert) == pytest.approx(expected)

    def test_failed_pass_returns_none_rather_than_raising(self, mc_case):
        """A pass that blows up is skipped, exactly as in the serial loop."""
        eng, cfg, q, I_obs, dI = mc_case
        cfg_mc = M.deepcopy(cfg)
        run = _MCRun(cfg_mc, q, dI, n_keys=4)
        # Wrong-length intensity → fit() raises inside the pass.
        assert run(I_obs[:5]) is None


# ── Serial vs parallel equivalence ───────────────────────────────────────────

class TestSerialParallelIdentity:
    def test_same_seed_gives_identical_sigmas(self, mc_case, monkeypatch):
        """The headline guarantee: with the pool forced on, parallel results are
        bit-identical to serial ones."""
        eng, cfg, q, I_obs, dI = mc_case
        monkeypatch.setattr(M, '_MC_PARALLEL_MIN_SECONDS', 0.0)   # force the pool

        serial = eng.calculate_uncertainty_mc(
            cfg, q, I_obs, dI, n_runs=4, workers=1, seed=42)
        parallel = eng.calculate_uncertainty_mc(
            cfg, q, I_obs, dI, n_runs=4, workers=2, seed=42)

        assert serial, "serial run produced no uncertainties"
        assert serial == parallel

    def test_progress_counts_every_pass_exactly_once(self, mc_case, monkeypatch):
        eng, cfg, q, I_obs, dI = mc_case
        monkeypatch.setattr(M, '_MC_PARALLEL_MIN_SECONDS', 0.0)

        seen = []
        eng.calculate_uncertainty_mc(
            cfg, q, I_obs, dI, n_runs=4, workers=2, seed=1,
            progress_cb=lambda done, total: seen.append((done, total)))

        assert seen == [(1, 4), (2, 4), (3, 4), (4, 4)]

    def test_different_seeds_give_different_sigmas(self, mc_case):
        """Guards the RNG plumbing: N passes must not share one noise draw."""
        eng, cfg, q, I_obs, dI = mc_case
        a = eng.calculate_uncertainty_mc(cfg, q, I_obs, dI, n_runs=4, workers=1, seed=1)
        b = eng.calculate_uncertainty_mc(cfg, q, I_obs, dI, n_runs=4, workers=1, seed=2)
        assert a and b and a != b
        # ...and every σ is a real spread, not a collapsed zero.
        assert all(v > 0 for v in a.values())


# ── Choosing between pool and serial ─────────────────────────────────────────

class TestPoolSelection:
    def test_short_run_stays_serial(self, mc_case, monkeypatch):
        """A job too small to repay pool startup must never build one, even when
        plenty of workers are requested."""
        eng, cfg, q, I_obs, dI = mc_case
        monkeypatch.setattr(M, '_MC_PARALLEL_MIN_SECONDS', 1e6)   # nothing qualifies

        called = []
        monkeypatch.setattr(M, '_run_mc_parallel',
                            lambda *a, **k: called.append(1) or [])

        stds = eng.calculate_uncertainty_mc(
            cfg, q, I_obs, dI, n_runs=3, workers=8, seed=3)

        assert called == [], "pool was built for a job that should stay serial"
        assert stds

    def test_single_worker_never_builds_a_pool(self, mc_case, monkeypatch):
        eng, cfg, q, I_obs, dI = mc_case
        monkeypatch.setattr(M, '_MC_PARALLEL_MIN_SECONDS', 0.0)

        called = []
        monkeypatch.setattr(M, '_run_mc_parallel',
                            lambda *a, **k: called.append(1) or [])

        eng.calculate_uncertainty_mc(cfg, q, I_obs, dI, n_runs=3, workers=1, seed=3)
        assert called == []


# ── Failure handling ─────────────────────────────────────────────────────────

class TestPoolFailureFallback:
    def test_unstartable_pool_falls_back_to_serial(self, mc_case, monkeypatch):
        """If the host cannot start workers the estimate must still complete —
        and must not double-count the passes it retries serially."""
        eng, cfg, q, I_obs, dI = mc_case
        monkeypatch.setattr(M, '_MC_PARALLEL_MIN_SECONDS', 0.0)

        def _boom(*a, **k):
            raise OSError("simulated: cannot start worker processes")

        monkeypatch.setattr('concurrent.futures.ProcessPoolExecutor', _boom)

        seen = []
        with pytest.warns(RuntimeWarning, match='running serially'):
            parallel = eng.calculate_uncertainty_mc(
                cfg, q, I_obs, dI, n_runs=4, workers=2, seed=42,
                progress_cb=lambda done, total: seen.append(done))

        # Every pass ran exactly once, and the answer matches a plain serial run.
        assert seen == [1, 2, 3, 4]
        serial = eng.calculate_uncertainty_mc(
            cfg, q, I_obs, dI, n_runs=4, workers=1, seed=42)
        assert parallel == serial


# ── Cancellation ─────────────────────────────────────────────────────────────

class TestCancel:
    def test_cancel_stops_early_and_keeps_completed_passes(self, mc_case):
        eng, cfg, q, I_obs, dI = mc_case

        seen = []
        stds = eng.calculate_uncertainty_mc(
            cfg, q, I_obs, dI, n_runs=20, workers=1, seed=5,
            progress_cb=lambda done, total: seen.append(done),
            cancel_cb=lambda: len(seen) >= 3)

        assert len(seen) < 20, "cancel did not stop the run"
        assert stds, "uncertainties from completed passes were thrown away"

    def test_cancel_before_two_passes_yields_no_estimate(self, mc_case):
        """Fewer than two passes cannot give a standard deviation."""
        eng, cfg, q, I_obs, dI = mc_case
        stds = eng.calculate_uncertainty_mc(
            cfg, q, I_obs, dI, n_runs=20, workers=1, seed=5,
            cancel_cb=lambda: True)
        assert stds == {}


# ── Degenerate inputs ────────────────────────────────────────────────────────

class TestDegenerate:
    def test_no_fittable_parameters_returns_empty(self, mc_case):
        eng, cfg, q, I_obs, dI = mc_case
        cfg.fit_background = False
        cfg.populations[0].fit_scale = False
        cfg.populations[0].dist_params_fit = {
            'min_size': False, 'mean_size': False, 'sdeviation': False,
        }
        assert eng.calculate_uncertainty_mc(cfg, q, I_obs, dI, n_runs=4) == {}

    def test_single_pass_cannot_give_a_spread(self, mc_case):
        eng, cfg, q, I_obs, dI = mc_case
        assert eng.calculate_uncertainty_mc(
            cfg, q, I_obs, dI, n_runs=1, workers=1, seed=1) == {}


# ── State migration (schema 4 → 5 adds mc_workers) ───────────────────────────

class TestStateMigration:
    def test_old_modeling_state_gains_mc_workers(self, tmp_path):
        from pyirena.state.state_manager import StateManager

        old = {
            'modeling': {
                'schema_version': 4,
                'background': 0.0,
                'fit_background': True,
                'no_limits': False,
                'fit_method': 'local',
                'de_workers': 1,
                'n_mc_runs': 10,
                'populations': [],
            }
        }
        f = tmp_path / 'state.json'
        f.write_text(json.dumps(old))

        sm = StateManager(state_file=f)
        assert sm.load() is True
        mod = sm.state['modeling']
        assert mod['mc_workers'] == 0          # auto
        assert mod['schema_version'] == 5
        assert mod['de_workers'] == 1          # untouched by the new migration
