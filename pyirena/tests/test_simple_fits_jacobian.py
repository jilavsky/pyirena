"""Validate the Simple Fits analytic Jacobians against finite differences.

Eight registry models carry a closed-form ``'jacobian'`` (Guinier, Guinier Rod,
Guinier Sheet, Porod, Power Law, Debye Polymer Chain, Debye-Bueche,
Teubner-Strey); the rest fall back to scipy's finite-difference Jacobian.  These
tests pin each analytic Jacobian to a 5-point finite-difference reference —
including the complex-background columns and free/fixed parameter selection —
and check that analytic and finite-difference fits reach the same optimum.

The reference is a 5-point stencil evaluated at *realistic* parameter scales:
a plain difference of e.g. the Porod ``q**-4`` term at tiny q reaches ~1e12,
where adding a flat background and differencing loses it to float64 roundoff.
That is a limitation of the numerical reference, not the analytic derivative
(whose flat-background column is exactly 1); realistic scales avoid it.
"""

import numpy as np
import pytest

from pyirena.core.simple_fits import MODEL_REGISTRY, SimpleFitModel

# Models that carry an analytic Jacobian, with realistic parameter values where
# every parameter (including any flat background) is well determined.
_ANALYTIC_MODELS = {
    "Guinier": {"I0": 100.0, "Rg": 50.0},
    "Guinier Rod": {"I0": 100.0, "Rc": 20.0},
    "Guinier Sheet": {"I0": 100.0, "Rg": 20.0},
    "Porod": {"Kp": 0.001, "Background": 0.5},
    "Power Law": {"Prefactor": 0.01, "Exponent": 3.5, "Background": 0.5},
    "Debye Polymer Chain": {"Scale": 100.0, "Rg": 50.0},
    "Debye-Bueche": {"Prefactor": 1.0, "Eta": 0.05, "CorrLength": 80.0},
    "Teubner-Strey": {"Prefactor": 1.0, "A": 0.1, "C1": -30.0, "C2": 5000.0},
}


def _fd5(func, q, p0, j, step):
    """5-point central difference of func w.r.t. free parameter j."""
    def ev(mult):
        pv = list(p0)
        pv[j] += mult * step
        return func(q, *pv)
    return (-ev(2) + 8 * ev(1) - 8 * ev(-1) + ev(-2)) / (12 * step)


def _assemble(model, use_bg, params, fixed=None):
    """Return (func, jac, p0) with the same free/fixed wiring fit() uses."""
    m = SimpleFitModel()
    m.set_model(model)
    m.use_complex_bg = use_bg
    m.params.update(params)
    base_func = m._build_fit_func()
    base_jac = m._build_jac_func()
    assert base_jac is not None

    specs = m._active_param_specs()
    all_names = [s[0] for s in specs]
    fixed = fixed or {}
    free_specs = [(n, v, lo, hi) for (n, v, lo, hi) in specs if n not in fixed]
    free_names = [s[0] for s in free_specs]
    free_col_idx = [all_names.index(n) for n in free_names]
    p0 = [s[1] for s in free_specs]

    def func(q, *fv):
        fm = dict(zip(free_names, fv))
        av = [fixed[n] if n in fixed else fm[n] for n in all_names]
        return base_func(q, *av)

    def jac(q, *fv):
        fm = dict(zip(free_names, fv))
        av = [fixed[n] if n in fixed else fm[n] for n in all_names]
        return base_jac(q, *av)[:, free_col_idx]

    return func, jac, p0


def _assert_jac_matches(model, use_bg, params, fixed=None, tol=1e-6):
    q = np.logspace(-1.5, 0.0, 250)
    func, jac, p0 = _assemble(model, use_bg, params, fixed)
    Ja = jac(q, *p0)
    assert Ja.shape == (len(q), len(p0))
    assert np.all(np.isfinite(Ja))
    for j in range(len(p0)):
        step = max(1e-4 * abs(p0[j]), 1e-4)   # abs floor for zero-valued params
        ref = _fd5(func, q, p0, j, step)
        peak = np.abs(ref).max()
        big = np.abs(ref) > 1e-6 * max(peak, 1e-300)
        if not big.any():
            continue
        denom = np.maximum(np.abs(Ja[big, j]), np.abs(ref[big]))
        rel = np.median(np.abs(Ja[big, j] - ref[big]) / denom)
        assert rel < tol, f"{model} col {j}: median rel err {rel:.2e} >= {tol}"


class TestSimpleFitsAnalyticJacobian:
    """Each analytic Jacobian must equal the finite-difference reference."""

    @pytest.mark.parametrize("model", list(_ANALYTIC_MODELS))
    def test_model_jacobian(self, model):
        _assert_jac_matches(model, use_bg=False, params=_ANALYTIC_MODELS[model])

    @pytest.mark.parametrize("model", ["Guinier", "Debye-Bueche", "Teubner-Strey"])
    def test_with_complex_background(self, model):
        params = dict(_ANALYTIC_MODELS[model])
        params.update(BG_B=0.001, BG_P=3.0, BG_flat=0.5)
        _assert_jac_matches(model, use_bg=True, params=params)

    def test_with_fixed_parameter(self):
        _assert_jac_matches("Power Law", use_bg=False,
                            params=_ANALYTIC_MODELS["Power Law"],
                            fixed={"Exponent": 3.5})

    def test_hard_models_have_no_jacobian(self):
        """Models with quadratures/branches must not advertise an analytic jac."""
        for name in ("Sphere", "Spheroid", "Benedetti-Ciccariello",
                     "Hermans", "Hybrid Hermans", "Unified Born Green"):
            assert "jacobian" not in MODEL_REGISTRY[name]
            m = SimpleFitModel()
            m.set_model(name)
            assert m._build_jac_func() is None


class TestSimpleFitsJacobianFit:
    """Analytic and finite-difference fits must reach the same optimum."""

    def _fit_both(self, model, params, seed=0):
        q = np.logspace(-2.0, 0.0, 200)
        m = SimpleFitModel()
        m.set_model(model)
        m.params.update(params)
        I_true = m.compute(q)
        rng = np.random.RandomState(seed)
        I = I_true + 0.02 * I_true * rng.randn(len(q))
        dI = 0.05 * I_true

        def fit(analytic):
            mm = SimpleFitModel()
            mm.set_model(model)
            mm.params.update({k: v * 0.8 for k, v in params.items()})
            mm.use_analytic_jacobian = analytic
            return mm.fit(q, I, dI)

        return fit(True), fit(False)

    def test_guinier_agrees(self):
        ra, rf = self._fit_both("Guinier", {"I0": 100.0, "Rg": 50.0})
        assert ra["success"] and rf["success"]
        assert ra["reduced_chi2"] == pytest.approx(rf["reduced_chi2"], rel=1e-3)
        for k in ("I0", "Rg"):
            assert ra["params"][k] == pytest.approx(rf["params"][k], rel=1e-3)

    def test_power_law_agrees(self):
        ra, rf = self._fit_both("Power Law",
                                {"Prefactor": 0.01, "Exponent": 3.5,
                                 "Background": 0.5})
        assert ra["success"] and rf["success"]
        assert ra["reduced_chi2"] == pytest.approx(rf["reduced_chi2"], rel=1e-3)
        for k in ("Prefactor", "Exponent", "Background"):
            assert ra["params"][k] == pytest.approx(rf["params"][k], rel=2e-3)

    def test_degenerate_model_agrees_on_identifiable_combo(self):
        """Debye-Bueche: Prefactor and Eta only enter as Prefactor·Eta², so the
        two paths may split that product differently — the identifiable
        combination and chi-squared must still agree."""
        ra, rf = self._fit_both("Debye-Bueche",
                                {"Prefactor": 1.0, "Eta": 0.05,
                                 "CorrLength": 80.0})
        assert ra["success"] and rf["success"]
        assert ra["reduced_chi2"] == pytest.approx(rf["reduced_chi2"], rel=1e-3)
        assert ra["params"]["CorrLength"] == pytest.approx(
            rf["params"]["CorrLength"], rel=1e-3)
        combo_a = ra["params"]["Prefactor"] * ra["params"]["Eta"] ** 2
        combo_f = rf["params"]["Prefactor"] * rf["params"]["Eta"] ** 2
        assert combo_a == pytest.approx(combo_f, rel=1e-3)
