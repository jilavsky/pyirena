"""
Tests for pyirena.core.distributions — normalization, CDF consistency,
dispatcher routing, and the Igor-style radius-grid generator.
"""

import numpy as np
import pytest

from pyirena.core.distributions import (
    gauss_pdf, gauss_cdf,
    lognormal_pdf, lognormal_cdf,
    lsw_pdf, lsw_cdf,
    schulz_zimm_pdf, schulz_zimm_cdf,
    ardell_pdf, ardell_cdf,
    pdf, cdf,
    generate_radius_grid,
)


# (dist_type, params, integration range) for the parametrized tests.
# dist_type strings are the dispatcher keys used by core/modeling.py.
CASES = [
    ("gauss",       {"mean_size": 100.0, "width": 15.0},                        (1.0, 300.0)),
    ("lognormal",   {"min_size": 10.0, "mean_size": 80.0, "sdeviation": 0.3},   (10.0, 1000.0)),
    ("lsw",         {"location": 100.0},                                        (1.0, 200.0)),
    ("schulz_zimm", {"mean_size": 100.0, "width": 20.0},                        (1.0, 500.0)),
    ("ardell",      {"location": 100.0, "parameter": 2.5},                      (1.0, 400.0)),
]


class TestNormalization:
    @pytest.mark.parametrize("dist_type,params,rng", CASES)
    def test_pdf_integrates_to_one(self, dist_type, params, rng):
        r = np.linspace(*rng, 20000)
        area = np.trapezoid(pdf(dist_type, r, params), r)
        assert area == pytest.approx(1.0, abs=0.01), dist_type

    @pytest.mark.parametrize("dist_type,params,rng", CASES)
    def test_pdf_nonnegative(self, dist_type, params, rng):
        r = np.linspace(*rng, 2000)
        assert (pdf(dist_type, r, params) >= 0).all()

    @pytest.mark.parametrize("dist_type,params,rng", CASES)
    def test_cdf_monotonic_zero_to_one(self, dist_type, params, rng):
        r = np.linspace(*rng, 500)
        c = np.asarray(cdf(dist_type, r, params), dtype=float)
        assert (np.diff(c) >= -1e-9).all(), dist_type
        assert c[0] == pytest.approx(0.0, abs=0.01)
        assert c[-1] == pytest.approx(1.0, abs=0.01)


class TestKnownValues:
    def test_gauss_symmetric_about_mean(self):
        mean, width = 100.0, 15.0
        assert gauss_cdf(mean, mean, width) == pytest.approx(0.5, abs=1e-6)
        d = 10.0
        assert gauss_pdf(np.array([mean - d]), mean, width)[0] == pytest.approx(
            gauss_pdf(np.array([mean + d]), mean, width)[0])

    def test_gauss_matches_scipy(self):
        from scipy.stats import norm
        r = np.linspace(50, 150, 11)
        np.testing.assert_allclose(
            gauss_pdf(r, 100.0, 15.0), norm.pdf(r, 100.0, 15.0), rtol=1e-10)
        np.testing.assert_allclose(
            gauss_cdf(r, 100.0, 15.0), norm.cdf(r, 100.0, 15.0), rtol=1e-9)

    def test_lsw_zero_beyond_1p5_location(self):
        """LSW distribution has a hard cutoff at r = 1.5 * location."""
        loc = 100.0
        r = np.array([149.0, 150.5, 200.0])
        vals = lsw_pdf(r, loc)
        assert vals[0] > 0
        assert vals[1] == 0.0 and vals[2] == 0.0

    def test_lognormal_zero_below_min_size(self):
        vals = lognormal_pdf(np.array([5.0, 9.9]), min_size=10.0,
                             mean_size=80.0, sdeviation=0.3)
        assert (vals == 0.0).all()

    def test_lsw_cdf_continuous_at_cutoff(self):
        """Regression: lsw_cdf used to jump from ~0.226 to 1.0 at the
        1.5*location cutoff because the pdf was missing its (3/2)^(11/3)
        normalization factor."""
        loc = 100.0
        just_below = float(lsw_cdf(1.5 * loc - 0.5, loc))
        assert just_below == pytest.approx(1.0, abs=0.01)

    def test_schulz_zimm_mean(self):
        """First moment of the Schulz-Zimm pdf must equal mean_size."""
        r = np.linspace(0.1, 800.0, 40000)
        p = schulz_zimm_pdf(r, 100.0, 20.0)
        mean = np.trapezoid(r * p, r)
        assert mean == pytest.approx(100.0, rel=0.01)


class TestDispatchers:
    def test_pdf_routes_by_name(self):
        r = np.linspace(50, 150, 5)
        np.testing.assert_allclose(
            pdf("gauss", r, {"mean_size": 100.0, "width": 15.0}),
            gauss_pdf(r, 100.0, 15.0))

    def test_unknown_dist_raises(self):
        with pytest.raises((ValueError, KeyError)):
            pdf("NotADistribution", np.array([1.0]), {})


class TestGenerateRadiusGrid:
    def test_grid_properties(self):
        r_grid = generate_radius_grid("gauss",
                                      {"mean_size": 100.0, "width": 15.0},
                                      n_bins=50)
        r_grid = np.asarray(r_grid[0] if isinstance(r_grid, tuple) else r_grid)
        assert len(r_grid) == 50
        assert (np.diff(r_grid) > 0).all()
        # Grid must bracket the mode
        assert r_grid[0] < 100.0 < r_grid[-1]

    def test_grid_denser_near_mode(self):
        """Inverse-CDF binning: spacing near the mode must be tighter than
        at the tails."""
        r_grid = generate_radius_grid("gauss",
                                      {"mean_size": 100.0, "width": 15.0},
                                      n_bins=50)
        r_grid = np.asarray(r_grid[0] if isinstance(r_grid, tuple) else r_grid)
        spacing = np.diff(r_grid)
        mode_idx = np.argmin(np.abs(r_grid - 100.0))
        near_mode = spacing[max(mode_idx - 2, 0):mode_idx + 2].mean()
        tails = max(spacing[0], spacing[-1])
        assert near_mode < tails
