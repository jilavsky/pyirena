"""Validate the Unified Fit analytic Jacobian against finite differences.

The fit supplies an analytic Jacobian to ``least_squares`` instead of scipy's
finite-difference approximation (see ``UnifiedFitModel._jacobian``).  These tests
pin that Jacobian to a numerical reference across every model mode — base,
``link_B``, ``mass_fractal``, correlations, cross-level ``link_RGCO`` and slit
smearing — in the "fast path returns what the reference returns" spirit of
``test_modeling_fast_paths.py``.

The reference is a **5-point central difference** of ``_residuals``.  A plain
2-point difference is not accurate enough here: the Unified intensity spans many
decades, so differencing it suffers catastrophic cancellation and its own
truncation error (~1e-5) swamps any realistic Jacobian tolerance.  The 5-point
stencil's O(h^4) error lets the analytic derivative be checked to ~1e-7.
"""

import numpy as np
import pytest

from pyirena.core.unified import UnifiedFitModel


def _fd5_column(model, p0, j, rel_h=1e-4):
    """5-point central difference of ``_residuals`` w.r.t. parameter ``j``."""
    step = rel_h * max(abs(p0[j]), 1e-12)

    def f(mult):
        pv = p0.copy()
        pv[j] += mult * step
        return model._residuals(pv)

    col = (-f(2) + 8 * f(1) - 8 * f(-1) + f(-2)) / (12 * step)
    model._unpack_parameters(p0)  # restore model state
    return col


def _assert_jacobian_matches(model, q, tol=1e-6):
    """Analytic ``_jacobian`` equals the finite-difference reference.

    Compared only where the reference derivative is non-negligible relative to
    its own column peak — elsewhere both derivatives have decayed into
    floating-point underflow and a relative comparison is meaningless.
    """
    model.q_data = np.asarray(q, dtype=float)
    model.I_data = np.ones_like(model.q_data)          # residual weights/derivs
    model.error_data = 0.05 * np.ones_like(model.q_data)  # do not depend on I

    p0 = model._pack_parameters()
    J = model._jacobian(p0)
    assert J.shape == (len(q), len(p0))
    assert np.all(np.isfinite(J))

    for j in range(len(p0)):
        ref = _fd5_column(model, p0, j)
        peak = np.abs(ref).max()
        significant = np.abs(ref) > 1e-6 * max(peak, 1e-300)
        if not significant.any():
            continue
        denom = np.maximum(np.abs(J[significant, j]), np.abs(ref[significant]))
        rel = np.abs(J[significant, j] - ref[significant]) / denom
        assert np.median(rel) < tol, (
            f"column {j} median rel err {np.median(rel):.2e} exceeds {tol}"
        )


@pytest.fixture
def q():
    return np.logspace(np.log10(5e-4), np.log10(0.3), 300)


def _one_level(P=4.0, **flags):
    m = UnifiedFitModel(num_levels=1)
    lv = m.levels[0]
    lv.G, lv.Rg, lv.P, lv.B = 100.0, 200.0, P, 3e-7
    lv.fit_G = lv.fit_Rg = lv.fit_P = lv.fit_B = True
    m.background = 0.01
    m.fit_background = True
    for k, v in flags.items():
        setattr(lv, k, v)
    return m


class TestUnifiedAnalyticJacobian:
    """Analytic Jacobian must equal the finite-difference reference."""

    def test_base_porod(self, q):
        """Independent Rg, G, P, B, background with P > 3 (K = 1.0)."""
        _assert_jacobian_matches(_one_level(P=4.0), q)

    def test_base_low_power(self, q):
        """P < 3 exercises the K = 1.06 branch of the erf argument."""
        _assert_jacobian_matches(_one_level(P=2.5), q)

    def test_link_b(self, q):
        """link_B: B derived from G, Rg, P — chain rule into those columns."""
        _assert_jacobian_matches(_one_level(P=3.6, link_B=True), q)

    def test_mass_fractal(self, q):
        """mass_fractal: B from Gamma(P/2) — digamma term in dB/dP."""
        _assert_jacobian_matches(_one_level(P=2.5, mass_fractal=True), q)

    def test_correlations(self, q):
        """Born-Green correlation factor adds ETA and PACK columns."""
        m = _one_level(P=4.0, correlations=True, ETA=150.0, PACK=3.0)
        m.levels[0].fit_ETA = m.levels[0].fit_PACK = True
        _assert_jacobian_matches(m, q)

    def test_two_level_link_rgco(self, q):
        """link_RGCO couples an upper level's cutoff to the lower level's Rg."""
        m = UnifiedFitModel(num_levels=2)
        low, high = m.levels
        low.G, low.Rg, low.P, low.B = 8.0, 120.0, 4.0, 1.9e-7
        low.fit_G = low.fit_Rg = low.fit_P = low.fit_B = True
        high.G, high.Rg, high.P, high.B = 4000.0, 1200.0, 3.2, 1.4e-6
        high.link_RGCO = True
        high.fit_G = high.fit_Rg = high.fit_P = high.fit_B = True
        m.background = 0.02
        m.fit_background = True
        _assert_jacobian_matches(m, q)

    def test_slit_smeared(self, q):
        """Smearing is linear, so the smeared Jacobian must still match."""
        m = _one_level(P=4.0)
        m.use_slit_smearing = True
        m.slit_length = 0.03
        _assert_jacobian_matches(m, q)

    def test_column_order_matches_packing(self, q):
        """The Jacobian has exactly one column per packed free parameter."""
        m = _one_level(P=4.0)
        m.q_data = q
        m.I_data = np.ones_like(q)
        m.error_data = 0.05 * np.ones_like(q)
        assert m._jacobian(m._pack_parameters()).shape[1] == len(
            m._pack_parameters()
        )
