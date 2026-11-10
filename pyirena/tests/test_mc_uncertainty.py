"""Monte-Carlo uncertainty in core (issue #38).

One implementation per tool now serves the GUI, batch and api/control paths.
These tests pin the decisions made when the copies were merged:

* no measured σ → refused with a reason (no assumed 5 % of I), unless the
  caller passes an explicit ``noise_sigma`` (GUI/batch: I × the Data
  Selector's "Generated uncertainty fraction") or Sizes has fractional error
  on — in both the user chose σ;
* Sizes perturbs with the *effective* error (fractional / scaled);
* Unified refits the same model as the main fit: links and slit smearing kept;
* a seed reproduces the result exactly; ``n_ok`` of ``n_runs`` is reported.
"""

import copy

import numpy as np
import pytest

from pyirena.core import mc_uncertainty as mcu
from pyirena.core.simple_fits import SimpleFitModel
from pyirena.core.sizes import SizesDistribution
from pyirena.core.unified import UnifiedFitModel

Q = np.logspace(-3, -0.5, 120)


def _unified_data():
    m = UnifiedFitModel(num_levels=1)
    lv = m.levels[0]
    lv.G, lv.Rg, lv.P, lv.B = 1000.0, 50.0, 4.0, 1e-3
    m.background = 0.01
    I = m.calculate_intensity(Q)
    dI = 0.03 * I
    I_obs = I + dI * np.random.default_rng(1).standard_normal(len(Q))
    return m, I_obs, dI


# ── shared helpers ────────────────────────────────────────────────────────────

def test_check_sigma():
    assert mcu.check_sigma(None, 3)
    assert mcu.check_sigma(np.zeros(3), 3)
    assert mcu.check_sigma(np.ones(4), 3)
    assert mcu.check_sigma(np.ones(3), 3) is None


def test_perturbations_are_seeded_and_skip_bad_sigma():
    I = np.ones(5)
    s = np.array([0.1, 0.0, np.nan, 0.1, 0.1])
    a = mcu.perturbed_intensities(I, s, 3, seed=7)
    b = mcu.perturbed_intensities(I, s, 3, seed=7)
    assert all(np.array_equal(x, y) for x, y in zip(a, b, strict=True))
    for x in a:
        assert x[1] == 1.0 and x[2] == 1.0          # σ=0 / NaN left unperturbed
    c = mcu.perturbed_intensities(I, s, 3, seed=8)
    assert not np.array_equal(a[0], c[0])


def test_std_of_needs_two_finite_values():
    assert np.isnan(mcu.std_of([1.0]))
    assert np.isnan(mcu.std_of([1.0, np.nan]))
    assert mcu.std_of([1.0, 3.0]) == pytest.approx(np.std([1.0, 3.0], ddof=1))


# ── Unified ───────────────────────────────────────────────────────────────────

def test_unified_mc_reproducible_and_reports_counts():
    m, I_obs, dI = _unified_data()
    m.fit(Q, I_obs, dI)
    r1 = m.calculate_uncertainty_mc(Q, I_obs, dI, n_runs=6, seed=3)
    r2 = m.calculate_uncertainty_mc(Q, I_obs, dI, n_runs=6, seed=3)
    assert r1['available'] and r1['n_ok'] == 6 and r1['n_runs'] == 6
    assert r1['uncertainties'] == r2['uncertainties']
    rg_sigma = r1['uncertainties']['levels'][0]['Rg']
    assert 0 < rg_sigma < 5.0                       # 3 % noise on a clean Rg=50


def test_unified_mc_refuses_without_errors():
    m, I_obs, _ = _unified_data()
    r = m.calculate_uncertainty_mc(Q, I_obs, None, n_runs=5)
    assert r['available'] is False
    assert 'uncertainties' in r['reason']
    assert 'uncertainties' not in r or r.get('uncertainties') is None


def test_unified_mc_does_not_modify_model():
    m, I_obs, dI = _unified_data()
    m.fit(Q, I_obs, dI)
    before = m.to_dict()
    m.calculate_uncertainty_mc(Q, I_obs, dI, n_runs=4, seed=1)
    assert m.to_dict() == before


def test_unified_mc_refits_the_same_model(monkeypatch):
    """Links and slit smearing survive into every refit (the GUI dropped both)."""
    m, I_obs, dI = _unified_data()
    m.use_slit_smearing = True
    m.slit_length = 0.02
    m.levels[0].link_B = True
    seen = []
    real_fit = UnifiedFitModel.fit

    def spy(self, *a, **k):
        seen.append((self.use_slit_smearing, self.slit_length, self.levels[0].link_B,
                     k.get('no_limits')))
        return real_fit(self, *a, **k)

    monkeypatch.setattr(UnifiedFitModel, 'fit', spy)
    m.calculate_uncertainty_mc(Q, I_obs, dI, n_runs=3, seed=1, no_limits=True)
    assert seen == [(True, 0.02, True, True)] * 3


def test_unified_mc_counts_failures(monkeypatch):
    m, I_obs, dI = _unified_data()
    calls = {'n': 0}
    real_fit = UnifiedFitModel.fit

    def flaky(self, *a, **k):
        calls['n'] += 1
        if calls['n'] % 2:
            raise RuntimeError('boom')
        return real_fit(self, *a, **k)

    monkeypatch.setattr(UnifiedFitModel, 'fit', flaky)
    r = m.calculate_uncertainty_mc(Q, I_obs, dI, n_runs=6, seed=1)
    assert r['available'] and r['n_ok'] == 3 and r['n_runs'] == 6


def test_unified_mc_not_available_when_too_few_succeed(monkeypatch):
    m, I_obs, dI = _unified_data()
    monkeypatch.setattr(UnifiedFitModel, 'fit',
                        lambda self, *a, **k: {'success': False})
    r = m.calculate_uncertainty_mc(Q, I_obs, dI, n_runs=4, seed=1)
    assert r['available'] is False and r['n_ok'] == 0


def test_unified_mc_cancel_and_progress():
    m, I_obs, dI = _unified_data()
    progress = []
    r = m.calculate_uncertainty_mc(
        Q, I_obs, dI, n_runs=10, seed=1,
        progress_cb=lambda d, n: progress.append((d, n)),
        cancel_cb=lambda: len(progress) >= 3)
    assert progress == [(1, 10), (2, 10), (3, 10)]
    assert r['n_ok'] == 3 and r['available']


# ── Sizes ─────────────────────────────────────────────────────────────────────

def test_sizes_effective_error():
    s = SizesDistribution()
    I = np.array([1.0, -2.0])
    err = np.array([0.1, 0.2])
    assert np.array_equal(s.effective_error(I, err), err)
    s.error_scale = 2.0
    assert np.allclose(s.effective_error(I, err), [0.2, 0.4])
    assert s.effective_error(I, None) is None
    s.fractional_error, s.fractional_error_value = True, 0.05
    assert np.allclose(s.effective_error(I, err), [0.05, 0.10])
    assert np.allclose(s.effective_error(I, None), [0.05, 0.10])


def test_sizes_mc_perturbs_with_effective_error(monkeypatch):
    """With fractional error on, the noise is |I|·fraction, not the file σ."""
    _, I_obs, dI = _unified_data()
    s = SizesDistribution()
    s.fractional_error, s.fractional_error_value = True, 0.10
    drawn = []
    real = mcu.perturbed_intensities

    def spy(I, sigma, n, seed=None, floor=None):
        drawn.append(np.asarray(sigma))
        return real(I, sigma, n, seed, floor)

    monkeypatch.setattr(mcu, 'perturbed_intensities', spy)
    r = s.calculate_uncertainty_mc(Q, I_obs, dI, n_runs=3, seed=1)
    assert r['available']
    assert np.allclose(drawn[0], 0.10 * np.abs(I_obs))


def test_sizes_mc_refuses_without_errors_unless_fractional():
    _, I_obs, _ = _unified_data()
    s = SizesDistribution()
    assert s.calculate_uncertainty_mc(Q, I_obs, None, n_runs=3)['available'] is False
    s.fractional_error = True
    assert s.calculate_uncertainty_mc(Q, I_obs, None, n_runs=3, seed=1)['available']


def test_sizes_mc_reproducible_including_montecarlo_method():
    _, I_obs, dI = _unified_data()
    s = SizesDistribution()
    s.method = 'montecarlo'
    s.n_bins = 40
    s.montecarlo_n_repetitions = 4
    r1 = s.calculate_uncertainty_mc(Q, I_obs, dI, n_runs=3, seed=11)
    r2 = s.calculate_uncertainty_mc(Q, I_obs, dI, n_runs=3, seed=11)
    assert r1['available'] and r1['n_ok'] == 3
    assert np.array_equal(r1['distribution_std'], r2['distribution_std'])
    assert r1['rg_std'] == r2['rg_std']
    assert s.montecarlo_n_repetitions == 4          # the model itself is untouched


def test_sizes_fit_rng_seeds_montecarlo_inversion():
    """run_sizes_fit(random_seed=…) relied on np.random.seed, which the
    inversion's default_rng() ignored; the rng argument is what seeds it now."""
    _, I_obs, dI = _unified_data()
    s = SizesDistribution()
    s.method = 'montecarlo'
    s.n_bins = 40
    a = s.fit(Q, I_obs, dI, rng=5)['distribution']
    b = s.fit(Q, I_obs, dI, rng=5)['distribution']
    c = s.fit(Q, I_obs, dI, rng=6)['distribution']
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)


# ── Simple Fits ───────────────────────────────────────────────────────────────

def _guinier_data():
    sel = Q < 0.02
    m, I_obs, dI = _unified_data()
    return Q[sel], I_obs[sel], dI[sel]


def test_simple_mc_reproducible():
    q, I, dI = _guinier_data()
    sf = SimpleFitModel()
    sf.fit(q, I, dI)
    r1 = sf.calculate_uncertainty_mc(q, I, dI, n_runs=8, seed=2)
    r2 = sf.calculate_uncertainty_mc(q, I, dI, n_runs=8, seed=2)
    assert r1['available'] and r1['n_ok'] == 8
    assert r1['params_std'] == r2['params_std']
    assert set(r1['params_std']) == {'I0', 'Rg'}


def test_simple_mc_passes_fit_options(monkeypatch):
    q, I, dI = _guinier_data()
    sf = SimpleFitModel()
    seen = []
    real = SimpleFitModel.fit

    def spy(self, *a, **k):
        seen.append((k.get('fixed_params'), k.get('no_limits')))
        return real(self, *a, **k)

    monkeypatch.setattr(SimpleFitModel, 'fit', spy)
    sf.calculate_uncertainty_mc(q, I, dI, n_runs=2, fixed_params={'I0': 1.0},
                                no_limits=True, seed=1)
    assert seen == [({'I0': 1.0}, True)] * 2


def test_simple_mc_refusals():
    q, I, _ = _guinier_data()
    sf = SimpleFitModel()
    assert sf.calculate_uncertainty_mc(q, I, None, n_runs=4)['available'] is False
    calc = copy.deepcopy(sf)
    calc.model = 'Invariant'
    r = calc.calculate_uncertainty_mc(q, I, 0.03 * I, n_runs=4)
    assert r['available'] is False and 'calculation' in r['reason']


# ── explicit noise σ for data without uncertainties ──────────────────────────

def test_unified_noise_sigma_drives_noise_but_refits_keep_main_weighting(monkeypatch):
    m, I_obs, _ = _unified_data()
    seen = []
    real_fit = UnifiedFitModel.fit

    def spy(self, q, I, error=None, **k):
        seen.append(error)
        return real_fit(self, q, I, error, **k)

    monkeypatch.setattr(UnifiedFitModel, 'fit', spy)
    r = m.calculate_uncertainty_mc(Q, I_obs, None, n_runs=3, seed=1,
                                   noise_sigma=0.03 * I_obs)
    assert r['available'] and r['sigma_source'] == 'generated'
    assert seen == [None, None, None]          # unit weights, like the main fit


def test_measured_sigma_wins_over_noise_sigma():
    m, I_obs, dI = _unified_data()
    r = m.calculate_uncertainty_mc(Q, I_obs, dI, n_runs=3, seed=1,
                                   noise_sigma=10 * I_obs)
    assert r['sigma_source'] == 'measured'


def test_noise_sigma_scales_uncertainty():
    q, I, _ = _guinier_data()
    sf = SimpleFitModel()
    small = sf.calculate_uncertainty_mc(q, I, None, n_runs=8, seed=4, noise_sigma=0.01 * I)
    large = sf.calculate_uncertainty_mc(q, I, None, n_runs=8, seed=4, noise_sigma=0.04 * I)
    assert large['params_std']['Rg'] == pytest.approx(4 * small['params_std']['Rg'], rel=0.2)


def test_sizes_noise_sigma_and_source():
    _, I_obs, _ = _unified_data()
    s = SizesDistribution()
    r = s.calculate_uncertainty_mc(Q, I_obs, None, n_runs=3, seed=1, noise_sigma=0.03 * I_obs)
    assert r['available'] and r['sigma_source'] == 'generated'
    s.fractional_error = True
    r = s.calculate_uncertainty_mc(Q, I_obs, None, n_runs=3, seed=1)
    assert r['sigma_source'] == 'fractional_error'


def test_generated_noise_sigma_helper():
    from pyirena.state.state_manager import generated_noise_sigma

    class _State:
        def get(self, tool, key, default=None):
            return {'error_fraction': 0.02}.get(key, default)

    I = np.array([1.0, -2.0])
    assert np.allclose(generated_noise_sigma(I, None, _State()), [0.02, 0.04])
    assert generated_noise_sigma(I, np.ones(2), _State()) is None
