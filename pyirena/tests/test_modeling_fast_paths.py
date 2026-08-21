"""
Tests for the fast evaluation paths added to the Modeling hot loop.

These are pure-performance rewrites, so every test here asserts that the fast
path returns what the reference implementation returned, not that it is fast:

* ``schulz_zimm_cdf`` / ``gauss_cdf`` call ``scipy.special`` directly instead of
  going through ``scipy.stats``'s ``rv_continuous`` wrapper.
* ``_cs_g_from_pairs`` builds the core-shell G matrix as one vectorised
  expression instead of looping over radius bins.
* ``ModelingEngine`` memoises the radius grid and each population's I(Q) during
  a fit, so a step that moved none of a population's parameters reuses its
  curve.
"""

import numpy as np
import pytest
from scipy import stats

from pyirena.core import distributions as D
from pyirena.core.form_factors import (
    _GL_NODES,
    _GL_WEIGHTS,
    _build_g_css_sphere_by_core,
    _cs_g_from_pairs,
    _cs_spheroid_g_from_pairs,
    _sphere_amplitude,
)
from pyirena.core.modeling import (
    DiffractionPeakPopulation,
    ModelingConfig,
    ModelingEngine,
    PopulationEvaluationError,
    SizeDistPopulation,
)

# ── A. CDFs via scipy.special ────────────────────────────────────────────────

class TestCdfFastPaths:
    @pytest.mark.parametrize("mean_size, width", [
        (10.0, 1.0), (244.8, 41.36), (1e4, 10.0), (50.0, 49.0), (1e5, 1e4),
    ])
    def test_schulz_zimm_matches_scipy_stats(self, mean_size, width):
        """gammainc route is bit-identical to stats.gamma.cdf on the support."""
        r = mean_size * np.array([1e-6, 0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 1e4])
        a, scale = D._schulz_zimm_ab(mean_size, width)
        expected = stats.gamma.cdf(r, a=a, scale=scale)
        assert D.schulz_zimm_cdf(r, mean_size, width) == pytest.approx(
            expected, rel=0, abs=0)

    def test_schulz_zimm_non_positive_radius_is_zero(self):
        """gammainc returns NaN below zero; the Gamma CDF is 0 there.

        generate_radius_grid clamps at 1 A so it never asks, but ``cdf`` is
        public and Sizes calls it directly.
        """
        vals = D.schulz_zimm_cdf(np.array([-5.0, -1e-30, 0.0]), 244.8, 41.36)
        assert np.all(vals == 0.0)
        assert not np.any(np.isnan(vals))
        assert D.schulz_zimm_cdf(-2.0, 100.0, 10.0) == 0.0

    def test_schulz_zimm_scalar_in_scalar_out(self):
        """Callers index and compare the result; it must not become 0-d."""
        v = D.schulz_zimm_cdf(200.0, 244.8, 41.36)
        assert np.isscalar(v) or getattr(v, 'ndim', None) == 0
        assert float(v) == pytest.approx(
            float(stats.gamma.cdf(200.0, **dict(zip(
                ('a', 'scale'), D._schulz_zimm_ab(244.8, 41.36))))))

    @pytest.mark.parametrize("mean_size, width", [
        (0.0, 1.0), (200.0, 20.0), (-1e3, 1e-4), (1e4, 1e3),
    ])
    def test_gauss_matches_scipy_stats(self, mean_size, width):
        r = mean_size + width * np.array([-40, -10, -3, -1, 0, 1, 3, 10, 40.0])
        expected = stats.norm.cdf(r, loc=mean_size, scale=max(width, 1e-30))
        assert D.gauss_cdf(r, mean_size, width) == pytest.approx(
            expected, rel=0, abs=0)

    def test_radius_grid_unchanged(self):
        """The grid is what the CDF change was made for — it must not move."""
        params = {'mean_size': 244.8, 'width': 41.36}
        rg = D.generate_radius_grid('schulz_zimm', params, n_bins=199)
        assert rg.shape == (199,)
        assert np.all(np.diff(rg) > 0)
        # CDF at the grid ends brackets the 1 % tail precision it was built for.
        assert D.cdf('schulz_zimm', rg[0], params) < 0.02
        assert D.cdf('schulz_zimm', rg[-1], params) > 0.98


# ── C. vectorised core-shell G matrix ────────────────────────────────────────

def _cs_g_reference(q, r_c_arr, r_t_arr, sld_core, sld_shell, sld_solvent):
    """The per-bin loop this module replaced, kept as the reference."""
    d_rho_c = float(sld_core) - float(sld_shell)
    d_rho_s = float(sld_shell) - float(sld_solvent)
    G = np.empty((len(q), len(r_c_arr)), dtype=float)
    for j in range(len(r_c_arr)):
        r_t = float(r_t_arr[j])
        r_c = float(r_c_arr[j])
        V_t = (4.0 / 3.0) * np.pi * r_t ** 3
        V_c = (4.0 / 3.0) * np.pi * r_c ** 3
        F = (d_rho_c * V_c * _sphere_amplitude(q * r_c)
             + d_rho_s * V_t * _sphere_amplitude(q * r_t))
        G[:, j] = F ** 2 / V_t
    return G * 1e-4


class TestCoreShellGVectorised:
    @pytest.mark.parametrize("mode", ['by_total', 'by_core', 'by_shell'])
    def test_matches_per_bin_loop(self, mode):
        """All three polydispersity modes, including the varying-shell one."""
        rng = np.random.default_rng(20260821)
        for _ in range(25):
            q = np.sort(10.0 ** rng.uniform(-4, 0.5, 200))
            r = np.sort(10.0 ** rng.uniform(0, 4, 120))
            t = float(10.0 ** rng.uniform(-1, 3.5))
            if mode == 'by_total':
                r_c, r_t = np.maximum(r - t, 0.0), r
            elif mode == 'by_core':
                r_c, r_t = r, r + t
            else:                        # by_shell: shell thickness varies
                r_c = np.full_like(r, t)
                r_t = r_c + r
            slds = rng.uniform(-10, 15, 3)
            got = _cs_g_from_pairs(q, r_c, r_t, *slds)
            ref = _cs_g_reference(q, r_c, r_t, *slds)
            assert np.max(np.abs(got - ref)) <= 1e-11 * np.max(np.abs(ref))

    def test_zero_core_radius_is_pure_sphere(self):
        """R_core = 0 must drop the core term rather than divide by zero."""
        q = np.linspace(0.001, 0.5, 64)
        r_t = np.array([100.0, 200.0])
        G = _cs_g_from_pairs(q, np.zeros_like(r_t), r_t, 5.0, 3.0, 1.0)
        ref = _cs_g_reference(q, np.zeros_like(r_t), r_t, 5.0, 3.0, 1.0)
        assert np.all(np.isfinite(G))
        assert np.max(np.abs(G - ref)) <= 1e-11 * np.max(np.abs(ref))

    def test_single_bin_grid(self):
        """N = 1 skips the constant-shell branch; it must still be right."""
        q = np.linspace(0.002, 0.3, 32)
        G = _cs_g_from_pairs(q, np.array([80.0]), np.array([120.0]), 9.9, 9.4, 9.46)
        ref = _cs_g_reference(q, np.array([80.0]), np.array([120.0]), 9.9, 9.4, 9.46)
        assert G.shape == (32, 1)
        assert np.max(np.abs(G - ref)) <= 1e-11 * np.max(np.abs(ref))

    @pytest.mark.parametrize("qr", [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8])
    def test_guinier_regime_precision(self, qr):
        """sin(x) - x*cos(x) cancels as x -> 0 and needs the Taylor branch.

        Small particles at USAXS Q reach this regime; without the branch the
        result is only ~5 digits by qr = 1e-5 and meaningless by qr = 1e-7.
        Checked per element, not relative to the peak of the matrix.
        """
        r_t = 100.0
        q = np.array([qr / r_t])
        got = _cs_g_from_pairs(q, np.array([60.0]), np.array([r_t]), 9.9, 9.4, 9.46)
        ref = _cs_g_reference(q, np.array([60.0]), np.array([r_t]), 9.9, 9.4, 9.46)
        assert got[0, 0] == pytest.approx(ref[0, 0], rel=1e-12)


# ── C (cont.). the other two per-bin loops ───────────────────────────────────

def _css_g_reference(q, r_grid, sld_core, sld_shell1, sld_shell2, sld_solvent,
                     t_shell1, t_shell2):
    """The per-bin loop `_build_g_css_sphere_by_core` replaced."""
    r_c = np.asarray(r_grid, dtype=float)
    r_m = r_c + t_shell1
    r_t = r_m + t_shell2
    d1, d2, d3 = (sld_core - sld_shell1, sld_shell1 - sld_shell2,
                  sld_shell2 - sld_solvent)
    G = np.empty((len(q), len(r_c)), dtype=float)
    for j in range(len(r_c)):
        r_t_j = float(r_t[j])
        V_c = (4.0 / 3.0) * np.pi * float(r_c[j]) ** 3
        V_m = (4.0 / 3.0) * np.pi * float(r_m[j]) ** 3
        V_t = (4.0 / 3.0) * np.pi * r_t_j ** 3
        F = (d1 * V_c * _sphere_amplitude(q * float(r_c[j]))
             + d2 * V_m * _sphere_amplitude(q * float(r_m[j]))
             + d3 * V_t * _sphere_amplitude(q * r_t_j))
        G[:, j] = F ** 2 / V_t
    return G * 1e-4


def _cs_spheroid_g_reference(q, r_c_arr, r_t_arr, sld_core, sld_shell,
                             sld_solvent, aspect_ratio):
    """The per-bin loop `_cs_spheroid_g_from_pairs` replaced."""
    AR = float(aspect_ratio)
    cos_a = 0.5 * (_GL_NODES + 1.0)
    weights = 0.5 * _GL_WEIGHTS
    stretch = np.sqrt(1.0 + (AR ** 2 - 1.0) * cos_a ** 2)
    d_rho_c, d_rho_s = sld_core - sld_shell, sld_shell - sld_solvent
    G = np.empty((len(q), len(r_c_arr)), dtype=float)
    for j in range(len(r_c_arr)):
        R_c, R_t = float(r_c_arr[j]), float(r_t_arr[j])
        V_c = (4.0 / 3.0) * np.pi * R_c ** 3 * AR
        V_t = (4.0 / 3.0) * np.pi * R_t ** 3 * AR
        F_c = _sphere_amplitude(q[:, None] * (R_c * stretch))
        F_t = _sphere_amplitude(q[:, None] * (R_t * stretch))
        G[:, j] = ((d_rho_c * V_c * F_c + d_rho_s * V_t * F_t) ** 2 @ weights) / V_t
    return G * 1e-4


class TestOtherCoreShellBuildersVectorised:
    def test_css_sphere_matches_per_bin_loop(self):
        rng = np.random.default_rng(20260821)
        for _ in range(20):
            q = np.sort(10.0 ** rng.uniform(-4, 0.5, 150))
            r = np.sort(10.0 ** rng.uniform(0, 3.5, 90))
            t1 = float(10.0 ** rng.uniform(0, 2.5))
            t2 = float(10.0 ** rng.uniform(0, 2.5))
            slds = rng.uniform(-5, 15, 4)
            got = _build_g_css_sphere_by_core(q, r, 1.0, *slds, t1, t2)
            ref = _css_g_reference(q, r, *slds, t1, t2)
            assert np.max(np.abs(got - ref)) <= 1e-11 * np.max(np.abs(ref))

    @pytest.mark.parametrize("aspect_ratio", [0.2, 0.9, 1.0, 1.7, 8.0])
    def test_cs_spheroid_matches_per_bin_loop(self, aspect_ratio):
        """Oblate through prolate; AR = 1 must reduce to the sphere case."""
        rng = np.random.default_rng(31415)
        q = np.sort(10.0 ** rng.uniform(-3.5, 0.3, 120))
        r_t = np.sort(10.0 ** rng.uniform(1, 3.0, 70))
        r_c = np.maximum(r_t - 40.0, 0.0)
        slds = rng.uniform(-5, 15, 3)
        got = _cs_spheroid_g_from_pairs(q, r_c, r_t, *slds, aspect_ratio)
        ref = _cs_spheroid_g_reference(q, r_c, r_t, *slds, aspect_ratio)
        assert np.max(np.abs(got - ref)) <= 1e-11 * np.max(np.abs(ref))

    def test_cs_spheroid_unit_aspect_ratio_equals_cs_sphere(self):
        """AR = 1 is a sphere, so the two builders must agree."""
        q = np.logspace(-3, -0.5, 100)
        r_t = np.linspace(60.0, 300.0, 50)
        r_c = r_t - 25.0
        spheroid = _cs_spheroid_g_from_pairs(q, r_c, r_t, 9.9, 9.4, 9.46, 1.0)
        sphere = _cs_g_from_pairs(q, r_c, r_t, 9.9, 9.4, 9.46)
        assert np.max(np.abs(spheroid - sphere)) <= 1e-11 * np.max(np.abs(sphere))


# ── Strict evaluation during fitting ─────────────────────────────────────────

class TestStrictEvaluation:
    @staticmethod
    def _broken_config():
        cfg = _two_population_config()
        # An unknown form factor makes build_g_matrix raise ValueError, the same
        # shape of failure as a bad parameter combination or a coding error.
        cfg.populations[0].form_factor = 'no_such_form_factor'
        return cfg

    def test_display_path_still_warns_and_skips(self):
        """A half-typed parameter must not take the GUI panel down."""
        q = np.logspace(-2.4, -0.1, 64)
        with pytest.warns(RuntimeWarning, match="Population 1 failed"):
            I, pop_idx, _, _ = ModelingEngine().total_intensity(
                self._broken_config(), q, use_cache=True)
        assert pop_idx == [1]                 # only the peak survived
        assert np.all(np.isfinite(I))

    def test_strict_raises_instead_of_dropping(self):
        q = np.logspace(-2.4, -0.1, 64)
        with pytest.raises(PopulationEvaluationError, match="Population 1"):
            ModelingEngine().total_intensity(
                self._broken_config(), q, use_cache=False, strict=True)

    def test_fit_refuses_a_model_it_cannot_evaluate(self):
        """The regression this guards: the fit used to 'succeed' without it.

        With the population silently dropped, fit() converged on the remaining
        terms and returned a result whose reduced chi-squared bore no relation
        to the model the user asked for.
        """
        q = np.logspace(-2.4, -0.1, 200)
        I = np.full_like(q, 1e-3)
        dI = 0.05 * I
        with pytest.raises(PopulationEvaluationError):
            ModelingEngine().fit(self._broken_config(), q, I, dI)

    def test_healthy_fit_reports_no_rejected_trials(self):
        cfg = _two_population_config()
        q = np.logspace(-2.4, -0.1, 200)
        eng = ModelingEngine()
        I_true, _, _, _ = eng.total_intensity(cfg, q, use_cache=True)
        res = eng.fit(cfg, q, I_true, 0.01 * I_true)
        assert eng._eval_failures == 0
        assert not any('rejected' in w for w in (res.fit_warnings or []))


# ── B + D. radius-grid reuse and the per-population memo ─────────────────────

def _two_population_config():
    return ModelingConfig(
        populations=[
            SizeDistPopulation(
                dist_type='schulz_zimm',
                dist_params={'mean_size': 244.8, 'width': 41.36},
                form_factor='cs_sphere_by_total',
                ff_params={'sld_core': 9.911, 'sld_shell': 9.414,
                           'sld_solvent': 9.46, 't_shell': 75.68},
                scale=1.6e-3, n_bins=99,
            ),
            DiffractionPeakPopulation(
                peak_type='lorentzian', position=0.1027,
                amplitude=9.13e-4, width=0.0418,
            ),
        ],
        background=3.3e-4, q_min=0.004, q_max=0.8,
    )


class TestEngineMemo:
    def test_memo_returns_same_curve_as_uncached(self):
        cfg = _two_population_config()
        q = np.logspace(-2.4, -0.1, 300)
        cold = ModelingEngine()
        warm = ModelingEngine()
        expected, _, _, _ = cold.total_intensity(cfg, q, use_cache=False)
        warm.total_intensity(cfg, q, use_cache=False)          # fills the memo
        got, _, _, _ = warm.total_intensity(cfg, q, use_cache=False)
        assert got == pytest.approx(expected, rel=0, abs=0)

    def test_memo_invalidated_by_every_population_parameter(self):
        """Each field in the memo key must actually change the curve."""
        q = np.logspace(-2.4, -0.1, 200)
        eng = ModelingEngine()
        base_cfg = _two_population_config()
        base, _, _, _ = eng.total_intensity(base_cfg, q, use_cache=False)

        mutations = {
            'dist_params': lambda p: p.dist_params.update(mean_size=260.0),
            'n_bins':      lambda p: setattr(p, 'n_bins', 149),
            'ff_params':   lambda p: p.ff_params.update(t_shell=60.0),
            'scale':       lambda p: setattr(p, 'scale', 3.2e-3),
            'form_factor': lambda p: setattr(p, 'form_factor', 'cs_sphere_by_core'),
        }
        for name, mutate in mutations.items():
            cfg = _two_population_config()
            eng.total_intensity(cfg, q, use_cache=False)       # seed the memo
            mutate(cfg.populations[0])
            got, _, _, _ = eng.total_intensity(cfg, q, use_cache=False)
            assert not np.allclose(got, base), f"memo went stale on {name}"

    def test_memo_invalidated_by_contrast(self):
        """contrast is ignored by cs_* shapes, so check it on a plain sphere."""
        q = np.logspace(-2.4, -0.1, 200)
        eng = ModelingEngine()

        def sphere_cfg(contrast):
            cfg = _two_population_config()
            pop = cfg.populations[0]
            pop.form_factor, pop.ff_params, pop.contrast = 'sphere', {}, contrast
            return cfg

        base, _, _, _ = eng.total_intensity(sphere_cfg(1.0), q, use_cache=False)
        got, _, _, _ = eng.total_intensity(sphere_cfg(2.0), q, use_cache=False)
        assert not np.allclose(got, base)

    def test_memo_does_not_leak_across_q_grids(self):
        """The slit-smeared path evaluates on q_ext with the same populations."""
        cfg = _two_population_config()
        eng = ModelingEngine()
        q_a = np.logspace(-2.4, -0.1, 200)
        q_b = np.logspace(-2.4, -0.1, 200) * 1.05      # same length, different q
        ref = ModelingEngine()
        expected_a, _, _, _ = ref.total_intensity(cfg, q_a, use_cache=False)
        expected_b, _, _, _ = ModelingEngine().total_intensity(
            cfg, q_b, use_cache=False)
        got_a, _, _, _ = eng.total_intensity(cfg, q_a, use_cache=False)
        got_b, _, _, _ = eng.total_intensity(cfg, q_b, use_cache=False)
        assert got_a == pytest.approx(expected_a, rel=0, abs=0)
        assert got_b == pytest.approx(expected_b, rel=0, abs=0)

    def test_radius_grid_cache_tracks_dist_params(self):
        pop = _two_population_config().populations[0]
        eng = ModelingEngine()
        rg1 = eng.build_radius_grid(pop)
        assert eng.build_radius_grid(pop) is rg1          # reused, not rebuilt
        pop.dist_params['mean_size'] = 300.0
        rg2 = eng.build_radius_grid(pop)
        assert rg2 is not rg1
        assert rg2.mean() > rg1.mean()

    def test_clear_cache_drops_memo_and_grid(self):
        cfg = _two_population_config()
        eng = ModelingEngine()
        eng.total_intensity(cfg, np.logspace(-2.4, -0.1, 128), use_cache=False)
        assert eng._pop_memo and eng._rg_cache
        eng.clear_cache()
        assert not eng._pop_memo and not eng._rg_cache

    def test_calculate_pop_intensity_still_returns_three(self):
        """Public signature is unchanged by the internal 4-tuple split."""
        pop = _two_population_config().populations[0]
        out = ModelingEngine().calculate_pop_intensity(
            pop, np.logspace(-2.4, -0.1, 64))
        assert len(out) == 3
        I_pop, vol_dist, num_dist = out
        assert I_pop.shape == (64,)
        assert vol_dist.shape == num_dist.shape == (pop.n_bins,)

    def test_fit_is_unaffected_by_a_stale_engine(self):
        """fit() clears the caches, so reusing an engine cannot bias a fit."""
        cfg = _two_population_config()
        q = np.logspace(-2.4, -0.1, 250)
        eng = ModelingEngine()
        I_true, _, _, _ = eng.total_intensity(cfg, q, use_cache=True)
        dI = 0.01 * I_true

        other = _two_population_config()
        other.populations[0].dist_params['mean_size'] = 180.0
        eng.total_intensity(other, q, use_cache=False)     # poison the caches

        fit_cfg = _two_population_config()
        fit_cfg.populations[0].dist_params['mean_size'] = 235.0
        res = eng.fit(fit_cfg, q, I_true, dI)
        assert fit_cfg.populations[0].dist_params['mean_size'] == pytest.approx(
            244.8, rel=0.05)
        assert res.reduced_chi_squared < 1.0
