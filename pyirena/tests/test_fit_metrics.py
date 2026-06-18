"""
Unit tests for robust fit-quality diagnostics (pyirena.core.fit_metrics).

Covers the scenarios from the spec (planning/ai-agent/04-fit-quality-metrics.md
section 5): well-calibrated sigma, under-estimated sigma, localized misfit,
systematic misfit, missing sigma, the dof guard, and the I->0 fractional guard.
"""

import numpy as np
import pytest

from pyirena.core.fit_metrics import fit_quality_metrics


def _power_law(q, B=1.0e-3, P=4.0):
    """A smooth, strictly positive model intensity to fit against."""
    return B * q ** (-P) + 1.0


@pytest.fixture
def q():
    return np.logspace(-2, 0, 200)


@pytest.fixture
def rng():
    # Fixed seed -> deterministic test data (test code, not library code).
    return np.random.default_rng(12345)


# ---------------------------------------------------------------------------
#  Well-calibrated sigma
# ---------------------------------------------------------------------------

def test_well_calibrated_sigma(q, rng):
    """data = model + Gaussian noise at exactly sigma => s ~ 1, chi2 ~ 1, no outliers."""
    M = _power_law(q)
    sigma = 0.02 * M  # 2% noise
    I = M + rng.normal(0.0, 1.0, size=q.size) * sigma

    res = fit_quality_metrics(q, I, M, sigma, n_params=3)

    assert res["sigma_available"] is True
    assert res["n_valid"] == q.size
    assert res["robust_scale_s"] == pytest.approx(1.0, abs=0.25)
    assert res["reduced_chi2"] == pytest.approx(1.0, abs=0.25)
    assert res["realistic_reduced_chi2_floor"] == pytest.approx(1.0, abs=0.6)
    assert res["n_outliers_3s"] == 0
    # Random noise => short sign runs and near-zero sign autocorrelation.
    assert res["longest_same_sign_run"] < 15
    assert abs(res["sign_autocorr_lag1"]) < 0.3
    # median fractional uncertainty ~ 2%
    assert res["median_frac_uncertainty"] == pytest.approx(0.02, abs=0.01)


# ---------------------------------------------------------------------------
#  Under-estimated sigma (x3) -- the key "don't mistake scale for misfit" case
# ---------------------------------------------------------------------------

def test_underestimated_sigma_x3(q, rng):
    """Same data, sigma divided by 3 => s ~ 3, floor ~ 9, still NO 3-sigma outliers."""
    M = _power_law(q)
    true_sigma = 0.02 * M
    I = M + rng.normal(0.0, 1.0, size=q.size) * true_sigma
    reported_sigma = true_sigma / 3.0  # under-reported by 3x

    res = fit_quality_metrics(q, I, M, reported_sigma, n_params=3)

    assert res["robust_scale_s"] == pytest.approx(3.0, abs=0.7)
    assert res["sigma_misscale_factor"] == res["robust_scale_s"]
    assert res["realistic_reduced_chi2_floor"] == pytest.approx(9.0, rel=0.5)
    # Reduced chi2 should be ~ 9 (scale problem), NOT ~ 1.
    assert res["reduced_chi2"] > 4.0
    # Crucial: a pure scale problem produces no outliers above 3*s.
    assert res["n_outliers_3s"] == 0
    # And no spurious structure.
    assert res["longest_same_sign_run"] < 15


# ---------------------------------------------------------------------------
#  Localized misfit -- a bump the model can't follow in one Q band
# ---------------------------------------------------------------------------

def test_localized_misfit(q, rng):
    """Inject a localized bump => bulk s ~ 1, outliers > 0, hot band, large frac misfit."""
    M = _power_law(q)
    sigma = 0.02 * M
    I = M + rng.normal(0.0, 1.0, size=q.size) * sigma

    # Inject a strong localized excess around q ~ 0.1 (mid range).
    bump_center = 0.1
    bump = np.exp(-0.5 * ((np.log(q) - np.log(bump_center)) / 0.05) ** 2)
    I = I + 0.8 * M * bump  # up to 80% excess locally

    res = fit_quality_metrics(q, I, M, sigma, n_params=3, n_bands=4)

    # Bulk scatter still ~ 1 (MAD ignores the localized outliers).
    assert res["robust_scale_s"] == pytest.approx(1.0, abs=0.6)
    assert res["n_outliers_3s"] > 0
    # Gross fractional misfit, located inside the injected bump.
    assert res["max_abs_frac_misfit"] > 0.3
    assert res["q_at_max_frac_misfit"] == pytest.approx(bump_center, rel=0.5)
    # The band(s) containing the bump should be far hotter than the clean bands.
    band_chi2 = [b["reduced_chi2"] for b in res["bands"] if b["reduced_chi2"]]
    assert max(band_chi2) > 10.0 * min(band_chi2)


# ---------------------------------------------------------------------------
#  Systematic misfit -- wrong slope => long sign runs, high autocorrelation
# ---------------------------------------------------------------------------

def test_systematic_misfit_wrong_slope(q, rng):
    """A model with the wrong power-law slope => long same-sign run, high autocorr."""
    M_true = _power_law(q, P=4.0)
    sigma = 0.02 * M_true
    I = M_true + rng.normal(0.0, 1.0, size=q.size) * sigma
    M_wrong = _power_law(q, P=3.6)  # systematically wrong shape

    res = fit_quality_metrics(q, I, M_wrong, sigma, n_params=3)

    # Even where magnitudes are modest, structure betrays the misfit.
    assert res["longest_same_sign_run"] > q.size // 4
    assert res["sign_autocorr_lag1"] > 0.5


# ---------------------------------------------------------------------------
#  sigma missing / all <= 0
# ---------------------------------------------------------------------------

def test_sigma_missing(q, rng):
    """sigma=None => sigma_available False, sigma-quantities None, fractional still populated, no crash."""
    M = _power_law(q)
    I = M * (1.0 + 0.05 * rng.normal(0.0, 1.0, size=q.size))

    res = fit_quality_metrics(q, I, M, None, n_params=3)

    assert res["sigma_available"] is False
    assert res["norm_residual"] is None
    assert res["frac_uncertainty"] is None
    assert res["robust_scale_s"] is None
    assert res["reduced_chi2"] is None
    assert res["n_outliers_3s"] is None
    assert res["dof"] is None
    # Fractional residual still computed.
    assert np.any(np.isfinite(res["frac_residual"]))
    assert np.isfinite(res["max_abs_frac_misfit"])
    # Structure still computed from raw residuals.
    assert res["longest_same_sign_run"] >= 0


def test_sigma_all_nonpositive(q):
    """sigma all <= 0 => treated as unavailable, no crash."""
    M = _power_law(q)
    I = M.copy()
    sigma = np.zeros_like(q)

    res = fit_quality_metrics(q, I, M, sigma, n_params=3)
    assert res["sigma_available"] is False
    assert res["robust_scale_s"] is None


# ---------------------------------------------------------------------------
#  dof <= 0 guard
# ---------------------------------------------------------------------------

def test_dof_guard():
    """n_params >= n_valid => reduced_chi2 is None (no divide), other stats still fine."""
    q = np.logspace(-2, 0, 5)
    M = _power_law(q)
    sigma = 0.02 * M
    I = M + 0.01 * M

    res = fit_quality_metrics(q, I, M, sigma, n_params=10)
    assert res["dof"] == 5 - 10
    assert res["reduced_chi2"] is None
    # robust scale still computable
    assert res["robust_scale_s"] is not None


# ---------------------------------------------------------------------------
#  I -> 0 fractional guard
# ---------------------------------------------------------------------------

def test_fractional_guard_near_zero_intensity(rng):
    """Points where |I| ~ 0 must not blow up the fractional arrays."""
    q = np.logspace(-2, 0, 100)
    M = _power_law(q)
    sigma = 0.02 * M
    I = M + rng.normal(0.0, 1.0, size=q.size) * sigma
    # Force a few points to ~zero intensity (e.g. over-subtracted background).
    I[10] = 0.0
    I[20] = 1e-12

    res = fit_quality_metrics(q, I, M, sigma, n_params=3)
    # The near-zero points are excluded (NaN), not infinite.
    assert np.isnan(res["frac_residual"][10])
    assert np.isnan(res["frac_residual"][20])
    # Max fractional misfit remains finite.
    assert np.isfinite(res["max_abs_frac_misfit"])


# ---------------------------------------------------------------------------
#  Adaptive banding with few points
# ---------------------------------------------------------------------------

def test_adaptive_band_count_small_data():
    """With few points, n_bands is reduced so each band keeps enough points."""
    q = np.logspace(-2, 0, 8)
    M = _power_law(q)
    sigma = 0.02 * M
    I = M.copy()

    res = fit_quality_metrics(q, I, M, sigma, n_params=2, n_bands=4)
    assert res["n_bands_used"] <= 4
    assert res["n_bands_used"] >= 1
    # Sum of per-band counts equals number of valid points.
    assert sum(b["n"] for b in res["bands"]) == res["n_valid"]


def test_band_counts_sum_to_valid(q, rng):
    """Per-band point counts always partition the valid points."""
    M = _power_law(q)
    sigma = 0.02 * M
    I = M + rng.normal(0.0, 1.0, size=q.size) * sigma
    res = fit_quality_metrics(q, I, M, sigma, n_params=3, n_bands=4)
    assert sum(b["n"] for b in res["bands"]) == res["n_valid"]
