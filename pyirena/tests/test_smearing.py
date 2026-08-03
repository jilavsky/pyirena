"""
Tests for the core slit-smearing engine (:mod:`pyirena.core.smearing`).

Units (1-6) need no data files.  Integration test 7 uses the Matilda USAXS
file ``testData/BoehNaNO2_10m_34_84min_0923.h5``, which carries BOTH a
desmeared (DSM) entry and a slit-smeared (SMR) entry on the same Q grid — the
file IS the ground-truth pair: smearing the DSM curve must reproduce the SMR
curve.
"""

import os

import numpy as np
import pytest
from scipy.integrate import quad

from pyirena.core.smearing import (
    SlitSmearer,
    build_extended_q,
    build_smearing_matrix,
    smear_curve,
    smear_model,
)

_TESTDATA = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
    "testData",
    "BoehNaNO2_10m_34_84min_0923.h5",
)


def _log_q(qmin=1e-3, qmax=0.3, n=200):
    return np.geomspace(qmin, qmax, n)


# --------------------------------------------------------------------------- #
# 1. SL <= 0  =>  strict no-op (all APIs, incl. the matrix)
# --------------------------------------------------------------------------- #
def test_noop_identity_matrix():
    q = _log_q()
    W, q_ext = build_smearing_matrix(q, slit_length=0.0)
    assert q_ext.shape == q.shape
    np.testing.assert_array_equal(q_ext, q)
    I = q ** -3.0
    np.testing.assert_allclose(W @ I, I, rtol=0, atol=0)


def test_noop_smear_model_and_curve():
    q = _log_q()
    def fn(x):
        return x ** -4.0
    np.testing.assert_array_equal(smear_model(fn, q, 0.0), fn(q))
    np.testing.assert_array_equal(smear_curve(q, fn(q), 0.0), fn(q))


def test_noop_smearer():
    q = _log_q()
    s = SlitSmearer(q, -1.0)
    assert s.is_noop
    np.testing.assert_array_equal(s.q_ext, q)
    I = q ** -2.0
    np.testing.assert_array_equal(s.smear(I), I)


# --------------------------------------------------------------------------- #
# 2. Flat curve smears to itself (E4 — constant background)
# --------------------------------------------------------------------------- #
def test_flat_curve_invariant():
    q = _log_q()
    SL = 0.02
    bkg = 3.7
    # analytic model path
    out = smear_model(lambda x: np.full_like(x, bkg), q, SL)
    np.testing.assert_allclose(out, bkg, rtol=1e-6)
    # tabulated path
    out2 = smear_curve(q, np.full_like(q, bkg), SL, extrapolation="flat")
    np.testing.assert_allclose(out2, bkg, rtol=1e-6)


# --------------------------------------------------------------------------- #
# 3. Power law q^-P  ->  smeared slope ~ -(P-1) for q << SL, ~ -P for q >> SL
# --------------------------------------------------------------------------- #
def test_power_law_smeared_slopes():
    SL = 0.02
    P = 4.0
    # Wide grid so we have a regime well below and well above SL.
    q = np.geomspace(1e-4, 2.0, 400)
    I_sm = smear_model(lambda x: x ** -P, q, SL)

    def local_slope(qlo, qhi):
        sel = (q >= qlo) & (q <= qhi)
        s, _ = np.polyfit(np.log(q[sel]), np.log(I_sm[sel]), 1)
        return s

    # q << SL: infinite-slit result flattens the slope to -(P-1) = -3.
    assert local_slope(2e-4, 1e-3) == pytest.approx(-(P - 1), abs=0.15)
    # q >> SL: smearing negligible, recover -P = -4.
    assert local_slope(0.7, 1.8) == pytest.approx(-P, abs=0.15)


# --------------------------------------------------------------------------- #
# 4. Guinier vs direct numerical quadrature (per point) — the reference truth
# --------------------------------------------------------------------------- #
def test_guinier_vs_quadrature():
    SL = 0.015
    G, Rg = 100.0, 50.0
    def model(x):
        return G * np.exp(-(x ** 2) * Rg ** 2 / 3.0)
    q = np.geomspace(2e-3, 0.2, 60)
    I_sm = smear_model(model, q, SL)

    def quad_ref(qi):
        def integrand(t):
            return model(np.sqrt(qi ** 2 + t ** 2))
        val, _ = quad(integrand, 0.0, SL, limit=200)
        return val / SL

    ref = np.array([quad_ref(qi) for qi in q])
    rel = np.abs(I_sm - ref) / np.abs(ref)
    assert np.median(rel) < 1e-3
    # Max over the meaningful dynamic range (within ~6 decades of the peak).
    # The bare-Guinier tail decays as e^-(q²Rg²/3); below ~1e-6 of the peak it
    # is numerical dust (~1e-13 absolute) that no real fit reaches, and its
    # extreme curvature defeats any linear grid in *relative* terms.
    meaningful = ref > ref.max() * 1e-6
    assert np.max(rel[meaningful]) < 1e-2


# --------------------------------------------------------------------------- #
# 5. Matrix vs direct quadrature on a random smooth curve
# --------------------------------------------------------------------------- #
def test_matrix_vs_quadrature_smooth():
    SL = 0.01
    rng = np.random.default_rng(0)
    # A smooth positive curve: sum of a few power laws + a Guinier bump.
    a, b, c = rng.uniform(0.5, 2.0, 3)
    def model(x):
        return 10 * x ** -a + 0.1 * x ** -b + 5 * np.exp(-(x ** 2) * (30.0) ** 2 / 3.0) * c
    q = np.geomspace(3e-3, 0.15, 50)
    I_sm = smear_model(model, q, SL)

    def quad_ref(qi):
        def integrand(t):
            return model(np.sqrt(qi ** 2 + t ** 2))
        val, _ = quad(integrand, 0.0, SL, limit=200)
        return val / SL

    ref = np.array([quad_ref(qi) for qi in q])
    rel = np.abs(I_sm - ref) / np.abs(ref)
    assert np.median(rel) < 3e-3
    assert np.max(rel) < 2e-2


# --------------------------------------------------------------------------- #
# 6. Grid building (E1), subrange indexing (E2), NaN sanitising (E5)
# --------------------------------------------------------------------------- #
def test_extended_q_spans_slit_and_contains_original():
    q = np.geomspace(1e-3, 0.01, 50)   # qmax << 3*SL
    SL = 0.02
    q_ext = build_extended_q(q, SL)
    assert q_ext[-1] >= 3 * SL - 1e-12
    assert q_ext[-1] >= np.sqrt(q[-1] ** 2 + SL ** 2)
    # original nodes are a subset of the extended grid
    assert np.all(np.isin(q, q_ext))
    assert np.all(np.diff(q_ext) > 0)


def test_extended_q_minimal_when_long_enough():
    q = np.geomspace(1e-3, 2.0, 100)
    SL = 0.001
    q_ext = build_extended_q(q, SL, refine=1)
    # The l-integral at qmax reaches sqrt(qmax²+SL²) > qmax, so a minimal
    # extension is always needed for SL>0 — but only a point or two here.
    assert q_ext[-1] >= np.sqrt(q[-1] ** 2 + SL ** 2)
    assert len(q_ext) <= len(q) + 3


def test_subrange_indexing_matches_full_grid():
    # E2: smear on the FULL grid then index the sub-range == identical values
    # to what those points get in the full smear (never smear a sub-range alone).
    SL = 0.015
    q = np.geomspace(1e-3, 0.2, 120)
    def model(x):
        return 50 * x ** -3.5 + 0.02
    full = smear_model(model, q, SL)
    lo, hi = 20, 60
    # Re-deriving via the smearer and slicing must give the same numbers.
    s = SlitSmearer(q, SL)
    full2 = s.smear_model(model)
    np.testing.assert_allclose(full[lo:hi], full2[lo:hi], rtol=1e-12)


def test_smear_curve_rejects_nonpositive():
    q = _log_q()
    with pytest.raises(ValueError):
        smear_curve(q, np.zeros_like(q), 0.01)


def test_extrapolation_none_raises_when_needed():
    q = np.geomspace(1e-3, 0.1, 80)
    with pytest.raises(ValueError):
        smear_curve(q, q ** -4.0, 0.05, extrapolation="none")


# --------------------------------------------------------------------------- #
# 7. Integration: smear DSM data -> reproduce the file's SMR curve
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not os.path.exists(_TESTDATA), reason="USAXS test file absent")
def test_dsm_smear_reproduces_smr_file():
    import h5py

    with h5py.File(_TESTDATA, "r") as f:
        dsm = f["/entry/BoehNaNO2_10m_34_84min/sasdata"]
        smr = f["/entry/BoehNaNO2_10m_34_84min_SMR/sasdata"]
        q = dsm["Q"][()]
        I_dsm = dsm["I"][()]
        I_smr = smr["I"][()]
        SL = float(smr["dQl"][()])

    I_model = smear_curve(q, I_dsm, SL, extrapolation="power_law")

    # Compare over the trustworthy mid range (desmearing itself is
    # ill-conditioned at the extreme ends).
    qmin, qmax = q.min(), q.max()
    mid = (q > 5 * qmin) & (q < 0.5 * qmax) & (I_smr > 0)
    rel = np.abs(I_model[mid] - I_smr[mid]) / np.abs(I_smr[mid])
    assert np.median(rel) < 0.01, f"median rel dev {np.median(rel):.4f}"
    assert np.percentile(rel, 90) < 0.05, f"90th pct rel dev {np.percentile(rel, 90):.4f}"


# --------------------------------------------------------------------------- #
# 8. Unified Fit: fitting slit-smeared data recovers ideal-space parameters
# --------------------------------------------------------------------------- #
def test_unified_smeared_fit_recovers_ideal_params():
    from pyirena.core.unified import UnifiedFitModel, UnifiedLevel

    q = np.geomspace(1e-3, 0.3, 300)
    truth = UnifiedFitModel(num_levels=1)
    truth.levels[0] = UnifiedLevel(G=500.0, Rg=80.0, B=2e-3, P=3.8)
    truth.background = 0.05
    SL = 0.0177
    # Make slit-smeared "data" via the INDEPENDENT tabulated-curve path so the
    # test is not circular with the model-based smearing used in the fit.
    I_data = smear_curve(q, truth.calculate_intensity(q), SL)
    err = 0.02 * I_data

    m = UnifiedFitModel(num_levels=1)
    m.levels[0] = UnifiedLevel(
        G=300, Rg=60, B=1e-3, P=4.0,
        fit_G=True, fit_Rg=True, fit_B=True, fit_P=True,
        G_limits=(1, 1e6), Rg_limits=(1, 1e4), B_limits=(1e-8, 1), P_limits=(1, 5),
    )
    m.background = 0.01
    m.fit_background = True
    m.background_limits = (0, 1)
    m.use_slit_smearing = True
    m.slit_length = SL
    m.fit(q, I_data, err)

    assert m.levels[0].Rg == pytest.approx(80.0, rel=0.03)
    assert m.levels[0].G == pytest.approx(500.0, rel=0.03)
    assert m.levels[0].P == pytest.approx(3.8, abs=0.05)
    assert m.reduced_chi_squared < 1.0

    # Sanity: fitting the SAME smeared data WITHOUT smearing is biased.
    m2 = UnifiedFitModel(num_levels=1)
    m2.levels[0] = UnifiedLevel(
        G=300, Rg=60, B=1e-3, P=4.0,
        fit_G=True, fit_Rg=True, fit_B=True, fit_P=True,
        G_limits=(1, 1e6), Rg_limits=(1, 1e4), B_limits=(1e-8, 1), P_limits=(1, 5),
    )
    m2.background = 0.01
    m2.fit_background = True
    m2.background_limits = (0, 1)
    m2.fit(q, I_data, err)
    assert m2.reduced_chi_squared > m.reduced_chi_squared


def test_unified_smearer_is_noop_when_disabled():
    from pyirena.core.unified import UnifiedFitModel

    q = np.geomspace(1e-3, 0.3, 100)
    m = UnifiedFitModel(num_levels=1)
    # smearing off -> smeared intensity identical to ideal
    np.testing.assert_array_equal(
        m.calculate_intensity_smeared(q), m.calculate_intensity(q)
    )


def test_local_guinier_smearing_warns_below_slit():
    from pyirena.core.unified import fit_local_guinier

    # Fit window entirely below the slit length -> E2 warning present.
    q = np.geomspace(1e-3, 5e-3, 30)
    I = 100 * np.exp(-(q ** 2) * 200.0 ** 2 / 3.0)
    r = fit_local_guinier(q, I, slit_length=0.02)
    assert "warning" in r and "slit length" in r["warning"]


# --------------------------------------------------------------------------- #
# 9. Public input-contract validation (only enforced when SL > 0)
# --------------------------------------------------------------------------- #
def test_validation_noop_passes_bad_input_through_when_sl_zero():
    """SL <= 0 must stay a strict no-op even for NaN/unsorted q (contract)."""
    q = _log_q()
    bad = q.copy()
    bad[5] = np.nan
    I = q ** -3.0
    np.testing.assert_array_equal(smear_curve(bad, I, 0.0), I)
    W, q_ext = build_smearing_matrix(bad, 0.0)
    assert q_ext.shape == bad.shape


def test_validation_rejects_nan_q():
    q = _log_q()
    q[10] = np.nan
    I = _log_q() ** -3.0
    with pytest.raises(ValueError, match="finite"):
        SlitSmearer(q, 0.015)
    with pytest.raises(ValueError, match="finite"):
        smear_curve(q, I, 0.015)


def test_validation_rejects_unsorted_q():
    q = _log_q()
    with pytest.raises(ValueError, match="strictly increasing"):
        SlitSmearer(q[::-1], 0.015)
    swapped = q.copy()
    swapped[3], swapped[4] = swapped[4], swapped[3]
    with pytest.raises(ValueError, match="strictly increasing"):
        build_smearing_matrix(swapped, 0.015)


def test_validation_rejects_length_mismatch():
    q = _log_q(n=60)
    with pytest.raises(ValueError, match="length mismatch"):
        smear_curve(q, (q[:50]) ** -3.0, 0.015)


def test_validation_rejects_small_grid_and_n_l():
    q = _log_q()
    with pytest.raises(ValueError, match="n_l"):
        build_smearing_matrix(q, 0.015, n_l=1)
    with pytest.raises(ValueError, match="at least 2 points"):
        SlitSmearer(q[:1], 0.015)


def test_validation_rejects_bad_q_ext():
    q = _log_q()
    with pytest.raises(ValueError, match="q_ext"):
        build_smearing_matrix(q, 0.015, q_ext=q[::-1].copy())
