"""Leftovers from the 2026-11 review (planning/unfinished, 11-11-2026).

* Modeling's global (DE) search takes a seed on every path.
* SAXS Morph and Data Manipulation use the "Generated uncertainty fraction"
  for data without uncertainties (and SAXS Morph repairs σ <= 0 with it).
* The Size Distribution GUI stores its Monte-Carlo per-bin σ.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from pyirena.tests.test_modeling_global_fit import _make_fit_config, _make_sphere_dataset


class _FakeState:
    """Just enough of StateManager for text_import_settings()."""

    def __init__(self, fraction):
        self.fraction = fraction

    def get(self, section, key, default=None):
        if (section, key) == ('data_selector', 'error_fraction'):
            return self.fraction
        return default


# ── Modeling: the DE seed ────────────────────────────────────────────────────

def _capture_de_seed(monkeypatch):
    import scipy.optimize

    import pyirena.core.modeling as modeling
    seen = []
    real = scipy.optimize.differential_evolution

    def spy(*args, **kwargs):
        seen.append(kwargs.get('seed'))
        kwargs['maxiter'] = 2               # keep the test fast
        return real(*args, **kwargs)

    monkeypatch.setattr(modeling, 'differential_evolution', spy)
    return seen


def test_engine_fit_passes_seed_to_global_search(monkeypatch):
    from pyirena.core.modeling import ModelingEngine
    seen = _capture_de_seed(monkeypatch)
    q, I, dI = _make_sphere_dataset()
    ModelingEngine().fit(_make_fit_config('global', 60.0), q, I, dI, seed=123)
    assert seen == [123]


def test_engine_fit_default_stays_unseeded(monkeypatch):
    from pyirena.core.modeling import ModelingEngine
    seen = _capture_de_seed(monkeypatch)
    q, I, dI = _make_sphere_dataset()
    ModelingEngine().fit(_make_fit_config('global', 60.0), q, I, dI)
    assert seen == [None]


def test_runner_passes_seed_to_global_search(monkeypatch):
    from pyirena.core.fit_runner import run_modeling
    seen = _capture_de_seed(monkeypatch)
    q, I, dI = _make_sphere_dataset()
    run_modeling(_make_fit_config('global', 60.0), q, I, dI, seed=7)
    assert seen == [7]


def test_seeded_global_fit_is_reproducible(monkeypatch):
    from pyirena.core.modeling import ModelingEngine
    _capture_de_seed(monkeypatch)
    q, I, dI = _make_sphere_dataset()
    r1 = ModelingEngine().fit(_make_fit_config('global', 60.0), q, I, dI, seed=5)
    r2 = ModelingEngine().fit(_make_fit_config('global', 60.0), q, I, dI, seed=5)
    assert r1.chi_squared == r2.chi_squared


# ── SAXS Morph: σ rule ───────────────────────────────────────────────────────

def test_saxs_morph_repair_sigma():
    from pyirena.core.saxs_morph import repair_sigma
    I = np.array([10.0, -4.0, 2.0, 1.0])
    np.testing.assert_allclose(repair_sigma(I, None, 0.1), [1.0, 0.4, 0.2, 0.1])
    dI = np.array([0.5, 0.0, -1.0, np.nan])
    np.testing.assert_allclose(repair_sigma(I, dI, 0.1), [0.5, 0.4, 0.2, 0.1])
    assert dI[1] == 0.0                     # the caller's array is not edited


def test_saxs_morph_compute_uses_error_fraction():
    """No dI: χ² scales with 1/fraction², i.e. the fraction is really used."""
    from pyirena.core.saxs_morph import SaxsMorphConfig, SaxsMorphEngine
    q = np.logspace(-2, -0.7, 60)
    I = 1e-3 * q ** -4 + 0.01
    cfg = SaxsMorphConfig(voxel_size_render=32, q_min=None, q_max=None, rng_seed=1)
    eng = SaxsMorphEngine()
    r05 = eng.compute_voxelgram(cfg, q, I, None, error_fraction=0.05)
    r10 = eng.compute_voxelgram(cfg, q, I, None, error_fraction=0.10)
    assert r05.chi_squared == pytest.approx(4.0 * r10.chi_squared, rel=1e-6)


# ── Data Manipulation: generated fraction ────────────────────────────────────

def test_data_manipulation_error_or_generated():
    pytest.importorskip("PySide6")
    pytest.importorskip("pyqtgraph")
    from pyirena.gui.data_manipulation_panel import DataManipulationPanel
    fake = SimpleNamespace(_sm=_FakeState(0.02))
    I = np.array([100.0, -50.0])
    out = DataManipulationPanel._error_or_generated(fake, {'Intensity': I})
    np.testing.assert_allclose(out, [2.0, 1.0])
    out = DataManipulationPanel._error_or_generated(fake, {'Intensity': I, 'Error': None})
    np.testing.assert_allclose(out, [2.0, 1.0])
    measured = np.array([3.0, 4.0])
    out = DataManipulationPanel._error_or_generated(fake, {'Intensity': I, 'Error': measured})
    assert out is measured


# ── Size Distribution GUI: MC σ is stored ────────────────────────────────────

def test_sizes_store_prefers_mc_distribution_std():
    pytest.importorskip("PySide6")
    pytest.importorskip("pyqtgraph")
    from pyirena.gui.sizes_panel import SizesFitPanel
    r = np.linspace(10, 100, 5)
    own = np.full(5, 9.0)
    mc = np.arange(5.0)
    pick = SizesFitPanel._distribution_std_to_store
    result = {'r_grid': r, 'distribution_std': own}
    assert pick(SimpleNamespace(_mc_dist_std=mc), result) is mc
    # No MC run, or an MC σ on another grid: the method's own spread (or None).
    assert pick(SimpleNamespace(_mc_dist_std=None), result) is own
    assert pick(SimpleNamespace(_mc_dist_std=np.arange(4.0)), result) is own
    assert pick(SimpleNamespace(_mc_dist_std=None), {'r_grid': r}) is None


# ── Size Distribution: MC σ of the scalars in the saved file ─────────────────

def _save_sizes(path, mc):
    from pyirena.io.nxcansas_sizes import save_sizes_results
    q = np.linspace(0.01, 0.2, 10)
    r = np.linspace(10, 100, 5)
    save_sizes_results(path, q, q * 0 + 1.0, q * 0 + 1.0, q * 0, r, r * 0 + 0.1,
                       params={'rg': 50.0, 'volume_fraction': 0.01},
                       distribution_std=r * 0 + 0.01, mc_uncertainty=mc)


def test_sizes_mc_scalars_round_trip(tmp_path):
    from pyirena.api.results import read_size_distribution
    from pyirena.io.nxcansas_sizes import load_sizes_results
    p = tmp_path / "s.h5"
    _save_sizes(p, {'available': True, 'n_runs': 20, 'n_ok': 19,
                    'sigma_source': 'measured', 'rg_std': 1.5,
                    'volume_fraction_std': 2e-4, 'peak_r_std': 3.0})
    res = load_sizes_results(p)
    assert res['rg_std'] == 1.5
    assert res['volume_fraction_std'] == pytest.approx(2e-4)
    assert res['peak_r_std'] == 3.0
    assert (res['mc_n_runs'], res['mc_n_ok']) == (20, 19)
    assert res['mc_sigma_source'] == 'measured'
    api = read_size_distribution(str(p))
    assert (api['rg_std'], api['peak_r_std'], api['mc_n_ok']) == (1.5, 3.0, 19)


def test_sizes_without_mc_reads_none(tmp_path):
    from pyirena.io.nxcansas_sizes import load_sizes_results
    p = tmp_path / "s.h5"
    _save_sizes(p, None)
    res = load_sizes_results(p)
    for k in ('rg_std', 'volume_fraction_std', 'peak_r_std', 'mc_n_runs', 'mc_n_ok'):
        assert res[k] is None
    # An MC run that failed stores nothing either.
    _save_sizes(p, {'available': False, 'n_runs': 20, 'n_ok': 1})
    assert load_sizes_results(p)['rg_std'] is None
