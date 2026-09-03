"""Validate the WAXS Peak Fit analytic Jacobian against finite differences.

The fit supplies ``curve_fit`` an analytic Jacobian (``_make_jac_func``) whenever
every peak has an analytic derivative — Gauss, Lorentz, Pseudo-Voigt — falling
back to finite differences when any peak is LogNormal.  These tests pin the
Jacobian to a 5-point finite-difference reference across peak shapes, polynomial
backgrounds, and free/fixed parameter mixes, and check the fallback + fit
behaviour end to end.

A 5-point stencil is the reference (not 2-point) so its O(h^4) truncation error
stays well below the tolerance — the same approach as the Unified Fit Jacobian
tests.
"""

import numpy as np
import pytest

from pyirena.core import waxs_peakfit as W


def _peak(shape, A, Q0, FWHM, eta=None, fit=True):
    d = {
        "shape": shape,
        "A": {"value": A, "fit": fit, "lo": 0.0, "hi": None},
        "Q0": {"value": Q0, "fit": fit, "lo": 0.0, "hi": None},
        "FWHM": {"value": FWHM, "fit": fit, "lo": 1e-6, "hi": None},
    }
    if eta is not None:
        d["eta"] = {"value": eta, "fit": fit, "lo": 0.0, "hi": 1.0}
    return d


def _bg(*coeffs, fit=True):
    return {f"bg{i}": {"value": c, "fit": fit, "lo": None, "hi": None}
            for i, c in enumerate(coeffs)}


def _assert_jac_matches(model, q, tol=1e-6):
    """``_make_jac_func`` equals a 5-point finite difference of the model."""
    p0, lb, ub, free_tags, fixed_vals = model._pack()
    f = model._make_model_func(free_tags, fixed_vals)
    jac = model._make_jac_func(free_tags, fixed_vals)

    Ja = jac(q, *p0)
    assert Ja.shape == (len(q), len(p0))
    assert np.all(np.isfinite(Ja))

    for j in range(len(p0)):
        step = 1e-5 * max(abs(p0[j]), 1e-9)

        def ev(mult, j=j, step=step):
            pv = list(p0)
            pv[j] += mult * step
            return f(q, *pv)

        ref = (-ev(2) + 8 * ev(1) - 8 * ev(-1) + ev(-2)) / (12 * step)
        peak = np.abs(ref).max()
        big = np.abs(ref) > 1e-6 * max(peak, 1e-300)
        if not big.any():
            continue
        denom = np.maximum(np.abs(Ja[big, j]), np.abs(ref[big]))
        rel = np.median(np.abs(Ja[big, j] - ref[big]) / denom)
        assert rel < tol, f"tag {free_tags[j]}: median rel err {rel:.2e} >= {tol}"


@pytest.fixture
def q():
    return np.linspace(0.5, 4.0, 600)


class TestWAXSAnalyticJacobian:
    """Analytic Jacobian must equal the finite-difference reference."""

    def test_gauss_linear_bg(self, q):
        m = W.WAXSPeakFitModel("Linear", [
            _peak("Gauss", 100, 1.85, 0.09),
            _peak("Gauss", 60, 2.55, 0.12),
            _peak("Gauss", 40, 3.05, 0.15)])
        m.bg_params = _bg(8.0, 1.5)
        _assert_jac_matches(m, q)

    def test_lorentz(self, q):
        m = W.WAXSPeakFitModel("Constant", [
            _peak("Lorentz", 50, 2.0, 0.08),
            _peak("Lorentz", 30, 2.6, 0.10)])
        m.bg_params = _bg(3.0)
        _assert_jac_matches(m, q)

    def test_pseudo_voigt_cubic_bg(self, q):
        m = W.WAXSPeakFitModel("Cubic", [
            _peak("Pseudo-Voigt", 200, 2.1, 0.07, eta=0.4),
            _peak("Pseudo-Voigt", 90, 2.8, 0.10, eta=0.6)])
        m.bg_params = _bg(10.0, -2.0, 0.8, -0.05)
        _assert_jac_matches(m, q)

    def test_mixed_shapes(self, q):
        m = W.WAXSPeakFitModel("Constant", [
            _peak("Gauss", 100, 1.5, 0.10),
            _peak("Lorentz", 80, 2.2, 0.09),
            _peak("Pseudo-Voigt", 60, 3.0, 0.12, eta=0.5)])
        m.bg_params = _bg(3.0)
        _assert_jac_matches(m, q)

    def test_some_params_fixed(self, q):
        """Fixed peak and background params must not get Jacobian columns."""
        m = W.WAXSPeakFitModel("Linear", [
            _peak("Lorentz", 50, 2.0, 0.08),
            _peak("Lorentz", 30, 2.6, 0.10)])
        m.peaks[0]["Q0"]["fit"] = False        # fix a peak param
        m.bg_params = _bg(5.0, 0.5)
        m.bg_params["bg1"]["fit"] = False       # fix a background coeff
        _assert_jac_matches(m, q)


class TestWAXSJacobianIntegration:
    """The Jacobian wiring must not change fit outcomes, and must fall back."""

    def _synth(self, q, peaks, bg_shape, bg_coeffs, seed=3):
        I = W.eval_background(q, bg_shape, bg_coeffs)
        for p in peaks:
            params = {k: p[k]["value"] for k in p if k != "shape"}
            I = I + W.eval_peak(q, p["shape"], params)
        rng = np.random.RandomState(seed)
        return I + 0.01 * I.max() * rng.randn(len(q))

    def test_lognormal_falls_back_to_fd(self, q):
        """A LogNormal peak disables the analytic path but the fit still runs."""
        peaks = [_peak("LogNormal", 100, 2.0, 0.15)]
        m = W.WAXSPeakFitModel("Constant", peaks)
        assert m._can_use_analytic_jac() is False
        I = self._synth(q, peaks, "Constant", [5.0])
        bg = _bg(3.0)
        res = m.fit(q, I, None, bg, [_peak("LogNormal", 70, 2.05, 0.20)])
        assert res["success"]

    def test_analytic_and_fd_agree(self, q):
        """Analytic and FD paths reach the same peak parameters."""
        truth = [_peak("Gauss", 100, 1.85, 0.09),
                 _peak("Gauss", 60, 2.55, 0.12)]
        I = self._synth(q, truth, "Linear", [8.0, 1.5])

        def fit(analytic):
            peaks = [_peak("Gauss", 70, 1.86, 0.11),
                     _peak("Gauss", 45, 2.56, 0.14)]
            bg = _bg(5.0, 1.0)
            m = W.WAXSPeakFitModel("Linear", peaks)
            m.use_analytic_jacobian = analytic
            return m.fit(q, I, None, bg, peaks)

        ra, rf = fit(True), fit(False)
        assert ra["success"] and rf["success"]
        assert ra["reduced_chi2"] == pytest.approx(rf["reduced_chi2"], rel=1e-3)
        for pa, pf in zip(ra["peaks"], rf["peaks"]):
            for name in ("A", "Q0", "FWHM"):
                assert pa[name]["value"] == pytest.approx(
                    pf[name]["value"], rel=1e-3)
