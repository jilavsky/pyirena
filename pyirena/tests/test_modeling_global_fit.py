"""
Unit tests for the Modeling tool's Global (Differential Evolution → local
polish) fit option, plus its config/state plumbing.
"""

import json

import numpy as np
import pytest

from pyirena.core.modeling import (
    ModelingEngine,
    ModelingConfig,
    SizeDistPopulation,
)


# ── Helpers ──────────────────────────────────────────────────────────────────

def _make_sphere_dataset(true_r=80.0, width=4.0, scale=0.01, seed=1):
    """Synthetic monodisperse-ish sphere I(Q) with sharp form-factor lobes."""
    eng = ModelingEngine()
    q = np.logspace(np.log10(0.005), np.log10(0.35), 250)
    truth = SizeDistPopulation()
    truth.enabled = True
    truth.dist_type = 'gauss'
    truth.dist_params = {'mean_size': true_r, 'width': width}
    truth.form_factor = 'sphere'
    truth.scale = scale
    truth.contrast = 1.0
    cfg = ModelingConfig(populations=[truth], background=0.0,
                         fit_background=False, q_min=q.min(), q_max=q.max())
    I_true, *_ = eng.total_intensity(cfg, q, use_cache=False)
    rng = np.random.default_rng(seed)
    dI = 0.02 * I_true + 1e-7
    I_obs = I_true + dI * rng.standard_normal(len(q))
    return q, I_obs, dI


def _make_fit_config(method, start_r):
    p = SizeDistPopulation()
    p.enabled = True
    p.dist_type = 'gauss'
    p.dist_params = {'mean_size': start_r, 'width': 4.0}
    p.dist_params_fit = {'mean_size': True, 'width': True}
    p.dist_params_limits = {'mean_size': (20.0, 200.0), 'width': (0.5, 40.0)}
    p.form_factor = 'sphere'
    p.scale = 0.01
    p.fit_scale = True
    p.scale_limits = (1e-6, 1.0)
    p.contrast = 1.0
    p.fit_contrast = False
    return ModelingConfig(
        populations=[p], background=0.0, fit_background=False,
        q_min=0.005, q_max=0.35, no_limits=False, fit_method=method,
    )


# ── Config defaults ──────────────────────────────────────────────────────────

class TestConfigDefaults:
    def test_fit_method_default_is_local(self):
        cfg = ModelingConfig()
        assert cfg.fit_method == 'local'


# ── Global vs local behaviour ────────────────────────────────────────────────

class TestGlobalFit:
    def test_global_recovers_radius_when_local_fails(self):
        """From a bad starting radius the local fit sticks in the wrong Bessel
        lobe; the global search recovers the true radius."""
        q, I_obs, dI = _make_sphere_dataset(true_r=80.0)
        eng = ModelingEngine()

        local = eng.fit(_make_fit_config('local', start_r=150.0), q, I_obs, dI)
        glob = eng.fit(_make_fit_config('global', start_r=150.0), q, I_obs, dI)

        r_local = local.config.populations[0].dist_params['mean_size']
        r_glob = glob.config.populations[0].dist_params['mean_size']

        # Global lands on the truth; local does not.
        assert abs(r_glob - 80.0) < 6.0, f"global R={r_glob}"
        assert abs(r_local - 80.0) > 15.0, f"local unexpectedly succeeded R={r_local}"
        # And the global fit is dramatically better.
        assert glob.reduced_chi_squared < local.reduced_chi_squared

    def test_local_default_path_runs(self):
        """A good starting guess with the default local method fits fine."""
        q, I_obs, dI = _make_sphere_dataset(true_r=80.0)
        eng = ModelingEngine()
        res = eng.fit(_make_fit_config('local', start_r=78.0), q, I_obs, dI)
        r = res.config.populations[0].dist_params['mean_size']
        assert abs(r - 80.0) < 6.0


# ── No-limits forces local (global needs finite bounds) ──────────────────────

class TestNoLimitsForcesLocal:
    def test_global_not_invoked_in_no_limits_mode(self, monkeypatch):
        q, I_obs, dI = _make_sphere_dataset(true_r=80.0)
        eng = ModelingEngine()

        def _boom(*a, **k):
            raise AssertionError("global DE must not run in no_limits mode")

        monkeypatch.setattr(eng, '_run_global_de', _boom)

        cfg = _make_fit_config('global', start_r=78.0)
        cfg.no_limits = True   # must override fit_method='global'
        # Should complete via the Nelder-Mead path without calling DE.
        res = eng.fit(cfg, q, I_obs, dI)
        assert res is not None

    def test_global_invoked_when_bounded(self, monkeypatch):
        q, I_obs, dI = _make_sphere_dataset(true_r=80.0)
        eng = ModelingEngine()
        called = {'n': 0}
        orig = eng._run_global_de

        def _spy(*a, **k):
            called['n'] += 1
            return orig(*a, **k)

        monkeypatch.setattr(eng, '_run_global_de', _spy)
        eng.fit(_make_fit_config('global', start_r=78.0), q, I_obs, dI)
        assert called['n'] == 1


# ── Log-sampling transform inside the DE stage ───────────────────────────────

class TestLogSampling:
    def test_de_seed_within_bounds_with_wide_decade_ranges(self):
        """The DE helper must return a vector inside the (linear) bounds even
        when some parameters span many decades and are searched in log space."""
        q, I_obs, dI = _make_sphere_dataset(true_r=80.0)
        eng = ModelingEngine()
        cfg = _make_fit_config('global', start_r=78.0)
        x0, lo, hi, keys = eng._pack_params(cfg)
        lo = np.array(lo); hi = np.array(hi); x0 = np.array(x0)
        # scale spans 1e-6..1 → a >5-decade range that must be log-sampled.
        x_seed = eng._run_global_de(lo, hi, keys, cfg, q, I_obs, dI, x0, seed=0)
        assert x_seed.shape == x0.shape
        assert np.all(x_seed >= lo - 1e-9)
        assert np.all(x_seed <= hi + 1e-9)


# ── State migration (schema 2 → 3 adds fit_method) ───────────────────────────

class TestStateMigration:
    def test_old_modeling_state_gains_fit_method(self, tmp_path):
        from pyirena.state.state_manager import StateManager

        old = {
            'modeling': {
                'schema_version': 2,
                'background': 0.0,
                'fit_background': True,
                'no_limits': False,
                'n_mc_runs': 10,
                'populations': [],
            }
        }
        f = tmp_path / 'state.json'
        f.write_text(json.dumps(old))

        sm = StateManager(state_file=f)
        assert sm.load() is True
        mod = sm.state['modeling']
        assert mod['fit_method'] == 'local'
        assert mod['schema_version'] == 3


# ── Exported JSON config carries fit_method (batch/scripting) ────────────────

class TestExportJsonCarriesFitMethod:
    """The 'Save params to JSON' export must include fit_method so a headless
    batch run (fit_modeling) uses the method the user picked in the GUI."""

    def test_export_includes_selected_fit_method(self, tmp_path, monkeypatch):
        pytest.importorskip("pyirena.gui.modeling_panel")
        try:
            try:
                from PySide6.QtWidgets import QApplication
            except ImportError:
                from PyQt6.QtWidgets import QApplication
        except Exception:
            pytest.skip("Qt not available")

        import os
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from pyirena.gui import modeling_panel as mp

        app = QApplication.instance() or QApplication([])
        panel = mp.ModelingPanel()
        panel.no_limits_cb.setChecked(False)
        gi = panel.fit_method_combo.findData('global')
        panel.fit_method_combo.setCurrentIndex(gi)

        out = tmp_path / 'cfg.json'
        monkeypatch.setattr(
            mp.QFileDialog, 'getSaveFileName',
            staticmethod(lambda *a, **k: (str(out), '')),
        )
        panel.export_json()

        data = json.loads(out.read_text())
        assert data['modeling']['fit_method'] == 'global'

        # Under no-limits the export must record the effective (local) method.
        panel.no_limits_cb.setChecked(True)
        panel.export_json()
        data2 = json.loads(out.read_text())
        assert data2['modeling']['fit_method'] == 'local'
