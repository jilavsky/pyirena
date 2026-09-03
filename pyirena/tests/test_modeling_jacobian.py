"""Validate the Modeling hybrid analytic/finite-difference Jacobian.

Modeling sums heterogeneous populations, and scipy needs one full Jacobian
matrix, so ``ModelingEngine._jacobian`` assembles it by hand: closed-form
columns for the cleanly-differentiable population types (unified_level,
diffraction_peak) and the constant background, finite-difference columns for the
rest (size-dist, guinier_porod, mass/surface fractal, scale, contrast).

These tests pin the *analytic* columns to a 5-point finite-difference reference
(the FD columns are excluded — comparing one finite difference against another
only measures truncation error), confirm the size-dist/fractal populations stay
on the finite-difference path, and check that analytic and finite-difference
fits reach the same optimum.
"""

import numpy as np
import pytest

from pyirena.core.modeling import (
    DiffractionPeakPopulation,
    GuinierPorodPopulation,
    ModelingConfig,
    ModelingEngine,
    UnifiedLevelPopulation,
)

# Groups whose Jacobian columns are analytic (plus the background column).
_ANALYTIC_GROUPS = {"uf", "peak"}


def _fd5_residual_column(eng, keys, cfg, x0, j, q, I, sigma):
    step = max(1e-4 * abs(x0[j]), 1e-5)

    def ev(mult):
        xv = x0.copy()
        xv[j] += mult * step
        return eng._residuals(xv, keys, cfg, q, I, sigma)

    col = (-ev(2) + 8 * ev(1) - 8 * ev(-1) + ev(-2)) / (12 * step)
    eng._residuals(x0, keys, cfg, q, I, sigma)   # restore config
    return col


def _assert_analytic_columns_match(cfg, tol=1e-6):
    """Every analytic (uf/peak/background) column equals the FD reference."""
    eng = ModelingEngine()
    q = np.logspace(-2.5, 0.0, 250)
    I = np.full_like(q, 5.0)          # nonzero so background is identifiable
    sigma = 0.05 * np.ones_like(q)

    x0, lo, hi, keys = eng._pack_params(cfg)
    x0 = np.array(x0, dtype=float)
    eng._fit_g_cache.clear()
    eng._rg_cache.clear()
    eng._pop_memo.clear()

    Ja = eng._jacobian(x0, keys, cfg, q, I, sigma)
    assert Ja.shape == (len(q), len(x0))
    assert np.all(np.isfinite(Ja))

    checked_analytic = 0
    for j, key in enumerate(keys):
        group = key[2] if key[0] == "pop" else "background"
        if group not in _ANALYTIC_GROUPS and group != "background":
            continue   # finite-difference column — skip
        checked_analytic += 1
        ref = _fd5_residual_column(eng, keys, cfg, x0, j, q, I, sigma)
        peak = np.abs(ref).max()
        big = np.abs(ref) > 1e-6 * max(peak, 1e-300)
        if not big.any():
            continue
        denom = np.maximum(np.abs(Ja[big, j]), np.abs(ref[big]))
        rel = np.median(np.abs(Ja[big, j] - ref[big]) / denom)
        assert rel < tol, f"{key}: median rel err {rel:.2e} >= {tol}"
    assert checked_analytic > 0


class TestModelingAnalyticJacobianColumns:

    def test_unified_level(self):
        c = ModelingConfig()
        p = UnifiedLevelPopulation(G=100, Rg=200, P=4.0, B=3e-7)
        p.fit_G = p.fit_Rg = p.fit_B = p.fit_P = True
        c.populations = [p]
        c.fit_background = True
        c.background = 0.5
        _assert_analytic_columns_match(c)

    def test_unified_level_low_power_and_rgco(self):
        c = ModelingConfig()
        p = UnifiedLevelPopulation(G=100, Rg=200, P=2.5, B=3e-7, RgCO=150.0)
        p.fit_G = p.fit_Rg = p.fit_B = p.fit_P = p.fit_RgCO = True
        c.populations = [p]
        c.fit_background = True
        c.background = 0.5
        _assert_analytic_columns_match(c)

    def test_unified_level_correlations(self):
        c = ModelingConfig()
        p = UnifiedLevelPopulation(G=100, Rg=200, P=4.0, B=3e-7,
                                   correlations=True, ETA=150.0, PACK=3.0)
        p.fit_G = p.fit_Rg = p.fit_B = True
        p.fit_ETA = p.fit_PACK = True
        c.populations = [p]
        c.fit_background = True
        c.background = 0.5
        _assert_analytic_columns_match(c)

    @pytest.mark.parametrize("peak_type", ["gaussian", "lorentzian", "voigt"])
    def test_diffraction_peak(self, peak_type):
        c = ModelingConfig()
        p = DiffractionPeakPopulation(peak_type=peak_type, position=0.3,
                                      amplitude=100.0, width=0.02, eta_voigt=0.4)
        p.fit_position = p.fit_amplitude = p.fit_width = True
        p.fit_eta_voigt = (peak_type == "voigt")
        c.populations = [p]
        c.fit_background = True
        c.background = 0.5
        _assert_analytic_columns_match(c)

    def test_mixed_analytic_and_fd_populations(self):
        """Unified (analytic) + Guinier-Porod (FD) in one model: the analytic
        columns must still be exact and the assembler must not confuse them."""
        c = ModelingConfig()
        p1 = UnifiedLevelPopulation(G=80, Rg=150, P=3.8, B=2e-7)
        p1.fit_G = p1.fit_Rg = p1.fit_B = True
        p2 = GuinierPorodPopulation(G=50, Rg1=80, P=3.5)
        p2.fit_G = p2.fit_Rg1 = True
        c.populations = [p1, p2]
        c.fit_background = True
        c.background = 0.5
        _assert_analytic_columns_match(c)


class TestModelingJacobianFit:

    def _fit_both(self, make_cfg):
        eng = ModelingEngine()
        q = np.logspace(-2.5, 0.0, 300)
        cfg_truth = make_cfg(truth=True)
        I_true, _, _, _ = eng.total_intensity(cfg_truth, q, use_cache=False,
                                              strict=True)
        rng = np.random.RandomState(5)
        I = I_true + 0.02 * I_true * rng.randn(len(q))
        dI = 0.05 * I_true

        def fit(analytic):
            eng.use_analytic_jacobian = analytic
            return eng.fit(make_cfg(truth=False), q, I, dI)

        return fit(True), fit(False)

    def test_unified_level_fit_agrees(self):
        def make(truth):
            c = ModelingConfig()
            c.fit_background = True
            if truth:
                p = UnifiedLevelPopulation(G=100, Rg=200, P=4.0, B=3e-7)
                c.background = 0.05
            else:
                p = UnifiedLevelPopulation(G=60, Rg=150, P=4.0, B=1e-6)
                c.background = 0.01
            p.fit_G = p.fit_Rg = p.fit_B = True
            c.populations = [p]
            return c

        ra, rf = self._fit_both(make)
        assert ra.reduced_chi_squared == pytest.approx(
            rf.reduced_chi_squared, rel=1e-3)
        pa, pf = ra.config.populations[0], rf.config.populations[0]
        for name in ("G", "Rg", "B"):
            assert getattr(pa, name) == pytest.approx(getattr(pf, name), rel=1e-2)

    def test_hybrid_model_fit_agrees(self):
        """A model mixing an analytic (unified) and FD (guinier-porod)
        population must converge identically with or without the analytic path."""
        def make(truth):
            c = ModelingConfig()
            c.fit_background = True
            if truth:
                p1 = UnifiedLevelPopulation(G=100, Rg=250, P=4.0, B=3e-7)
                p2 = GuinierPorodPopulation(G=50, Rg1=60, P=3.6)
                c.background = 0.05
            else:
                p1 = UnifiedLevelPopulation(G=70, Rg=200, P=4.0, B=1e-6)
                p2 = GuinierPorodPopulation(G=35, Rg1=50, P=3.6)
                c.background = 0.01
            p1.fit_G = p1.fit_Rg = p1.fit_B = True
            p2.fit_G = p2.fit_Rg1 = True
            c.populations = [p1, p2]
            return c

        ra, rf = self._fit_both(make)
        assert ra.reduced_chi_squared == pytest.approx(
            rf.reduced_chi_squared, rel=1e-3)

    def test_analytic_jacobian_fallback(self):
        """If _jacobian raises, the fit must still complete via finite
        differences and recover the parameters."""
        eng = ModelingEngine()
        q = np.logspace(-2.5, 0.0, 200)
        c_truth = ModelingConfig()
        c_truth.fit_background = True
        pt = UnifiedLevelPopulation(G=100, Rg=200, P=4.0, B=3e-7)
        c_truth.populations = [pt]
        c_truth.background = 0.05
        I_true, _, _, _ = eng.total_intensity(c_truth, q, use_cache=False,
                                              strict=True)
        dI = 0.05 * I_true

        def boom(*a, **k):
            raise RuntimeError("analytic Jacobian sabotaged")

        eng._jacobian = boom
        eng.use_analytic_jacobian = True

        c = ModelingConfig()
        c.fit_background = True
        p = UnifiedLevelPopulation(G=60, Rg=150, P=4.0, B=1e-6)
        p.fit_G = p.fit_Rg = p.fit_B = True
        c.populations = [p]
        c.background = 0.01
        res = eng.fit(c, q, I_true, dI)
        assert res.reduced_chi_squared < 1.0
        assert abs(res.config.populations[0].Rg - 200.0) / 200.0 < 0.1
