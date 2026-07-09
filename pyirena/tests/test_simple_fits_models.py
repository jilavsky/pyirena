"""
Fit-recovery tests for every model in the simple-fits registry.

Each model generates its own synthetic curve at the registry defaults;
the fit then starts from parameters perturbed by 10% and must reproduce
the curve. For the models without parameter degeneracies the fitted
parameters themselves must return to the defaults.

Treubner-Strey, Benedetti-Ciccariello and Hybrid Hermans have (near-)
degenerate parameter sets — different parameter values can produce the
same curve — so for those only curve reproduction is required.
"""

import numpy as np
import pytest

from pyirena.core.simple_fits import SimpleFitModel, MODEL_REGISTRY, MODEL_NAMES

Q = np.logspace(-2.3, -0.7, 200)

# Models whose parameters are uniquely identifiable from a clean curve
IDENTIFIABLE = [
    "Guinier", "Guinier Rod", "Guinier Sheet", "Porod", "Power Law",
    "Debye Polymer Chain", "Sphere", "Spheroid", "Debye-Bueche",
    "Hermans", "Unified Born Green",
]
DEGENERATE = ["Treubner-Strey", "Benedetti-Ciccariello", "Hybrid Hermans"]


def _self_consistency(name, perturbation=1.1):
    """Generate curve at defaults, fit from perturbed start, return
    (max relative curve error, params recovered?)."""
    m = SimpleFitModel()
    m.set_model(name)
    defaults = dict(m.params)
    I = m.compute(Q)
    assert np.isfinite(I).all()

    m2 = SimpleFitModel()
    m2.set_model(name)
    for k in m2.params:
        m2.params[k] = m2.params[k] * perturbation if m2.params[k] != 0 else 0.01
    result = m2.fit(Q, I, 0.01 * np.abs(I) + 1e-12)

    fitted = SimpleFitModel()
    fitted.set_model(name)
    fitted.params.update(
        {k: result["params"][k] for k in fitted.params if k in result["params"]}
    )
    I_fit = fitted.compute(Q)
    mask = I > 0
    curve_err = float(np.max(np.abs(I_fit[mask] - I[mask]) / I[mask]))

    recovered = all(
        abs(result["params"][k] - defaults[k]) <= 0.05 * max(abs(defaults[k]), 1e-12)
        for k in defaults if "background" not in k.lower()
    )
    return curve_err, recovered


class TestModelRegistry:
    def test_all_models_present(self):
        assert set(MODEL_NAMES) == set(MODEL_REGISTRY)
        assert set(IDENTIFIABLE) | set(DEGENERATE) == set(MODEL_NAMES)

    @pytest.mark.parametrize("name", MODEL_NAMES)
    def test_compute_finite_positive(self, name):
        m = SimpleFitModel()
        m.set_model(name)
        I = m.compute(Q)
        assert np.isfinite(I).all(), name
        assert (I > 0).any(), name


class TestFitSelfConsistency:
    @pytest.mark.parametrize("name", IDENTIFIABLE)
    def test_identifiable_models_recover_params(self, name):
        curve_err, recovered = _self_consistency(name)
        assert curve_err < 1e-2, f"{name}: curve error {curve_err:.2e}"
        assert recovered, f"{name}: parameters did not return to defaults"

    @pytest.mark.parametrize("name", DEGENERATE)
    def test_degenerate_models_reproduce_curve(self, name):
        curve_err, _ = _self_consistency(name)
        assert curve_err < 1e-2, f"{name}: curve error {curve_err:.2e}"


class TestSerialization:
    @pytest.mark.parametrize("name", MODEL_NAMES)
    def test_to_dict_from_dict_round_trip(self, name):
        m = SimpleFitModel()
        m.set_model(name)
        d = m.to_dict()
        m2 = SimpleFitModel.from_dict(d)
        assert m2.model == name
        assert m2.params == m.params
        np.testing.assert_allclose(m2.compute(Q), m.compute(Q))
