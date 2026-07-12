"""
Unit tests for the Sizes (size distribution) core module.

Tests cover:
  - Sphere form factor numerical correctness
  - Spheroid form factor degeneracy (AR=1 → sphere)
  - G-matrix shape and positivity
  - Bin-width helper
  - Each fitting method recovers a known synthetic distribution
  - TNNLS non-negativity
  - MaxEnt positivity
  - Regularization chi² convergence
  - HDF5 save/load round-trip
"""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from pyirena.core.form_factors import (
    sphere_ff,
    spheroid_ff,
    build_g_matrix,
    make_r_grid,
    bin_widths,
)
from pyirena.core.sizes import SizesDistribution


# ──────────────────────────────────────────────────────────────────────────────
# Helper: generate synthetic data from a known Gaussian distribution
# ──────────────────────────────────────────────────────────────────────────────

def _synthetic_data(r_true=200.0, r_sigma=30.0, vf_true=0.01,
                    q_min=1e-3, q_max=0.1, n_q=60,
                    r_min=50.0, r_max=600.0, n_bins=50,
                    contrast=1.0, noise_frac=0.02, seed=42):
    """
    Generate synthetic I(Q) from a Gaussian P(r) centred at r_true.

    Returns (q, I, err, r_grid, distribution_true).
    """
    rng = np.random.default_rng(seed)
    q = np.logspace(np.log10(q_min), np.log10(q_max), n_q)
    r_grid = make_r_grid(r_min, r_max, n_bins)

    # True P(r): Gaussian, normalised to vf_true
    dw = bin_widths(r_grid)
    dist_raw = np.exp(-0.5 * ((r_grid - r_true) / r_sigma) ** 2)
    dist_raw *= vf_true / float(np.trapezoid(dist_raw, r_grid))

    # G matrix
    G = build_g_matrix(q, r_grid, 'sphere', contrast)

    # True model = G × (dist × bin_widths)
    x_true = dist_raw * dw
    I_true = G @ x_true

    # Add noise
    err = I_true * noise_frac
    I = I_true + rng.normal(0.0, err)

    return q, I, err, r_grid, dist_raw


# ──────────────────────────────────────────────────────────────────────────────
# Form factor tests
# ──────────────────────────────────────────────────────────────────────────────

class TestSphereFF:
    def test_low_q_limit(self):
        """At q→0, sphere_ff → V (volume of sphere, not V²)."""
        r = 100.0
        V = (4.0 / 3.0) * np.pi * r ** 3
        q_tiny = np.array([1e-8])
        ff = sphere_ff(q_tiny, r)
        assert abs(ff[0] - V) / V < 1e-6

    def test_known_value(self):
        """Check against manually computed value at qr=1."""
        r = 100.0
        q = np.array([1.0 / r])   # qr = 1
        V = (4.0 / 3.0) * np.pi * r ** 3
        # F_norm(1) = 3[sin(1) - cos(1)] / 1
        F_expected = 3.0 * (np.sin(1.0) - np.cos(1.0)) / 1.0
        # sphere_ff returns V × F_norm²
        expected = V * F_expected ** 2
        ff = sphere_ff(q, r)
        assert abs(ff[0] - expected) / expected < 1e-9

    def test_positive(self):
        """Form factor must be non-negative everywhere."""
        q = np.logspace(-3, 0, 200)
        ff = sphere_ff(q, 100.0)
        assert np.all(ff >= 0)

    def test_shape(self):
        """Output shape matches q."""
        q = np.linspace(0.001, 0.1, 50)
        ff = sphere_ff(q, 50.0)
        assert ff.shape == (50,)


class TestSpheroidFF:
    def test_ar1_equals_sphere(self):
        """Spheroid with AR=1 must equal sphere form factor."""
        q = np.logspace(-3, -0.5, 40)
        r = 80.0
        ff_sph = sphere_ff(q, r)
        ff_sphd = spheroid_ff(q, r, aspect_ratio=1.0)
        np.testing.assert_allclose(ff_sphd, ff_sph, rtol=1e-4)

    def test_ar_not1_differs(self):
        """Spheroid with AR≠1 differs from sphere."""
        q = np.logspace(-2, -0.5, 30)
        r = 80.0
        ff_sph = sphere_ff(q, r)
        ff_sphd = spheroid_ff(q, r, aspect_ratio=2.0)
        # They should be different (not identical)
        assert not np.allclose(ff_sphd, ff_sph, rtol=1e-3)

    def test_positive(self):
        q = np.logspace(-3, -0.5, 80)
        ff = spheroid_ff(q, 100.0, aspect_ratio=3.0)
        assert np.all(ff >= 0)


# ──────────────────────────────────────────────────────────────────────────────
# G-matrix tests
# ──────────────────────────────────────────────────────────────────────────────

class TestGMatrix:
    def test_shape(self):
        q = np.logspace(-3, -1, 30)
        r_grid = make_r_grid(20, 500, 20)
        G = build_g_matrix(q, r_grid, 'sphere', contrast=1.0)
        assert G.shape == (30, 20)

    def test_positive(self):
        q = np.logspace(-3, -1, 30)
        r_grid = make_r_grid(20, 500, 20)
        G = build_g_matrix(q, r_grid, 'sphere', contrast=1.0)
        assert np.all(G >= 0)

    def test_contrast_scaling(self):
        """Doubling contrast should double G."""
        q = np.logspace(-3, -1, 20)
        r_grid = make_r_grid(50, 300, 15)
        G1 = build_g_matrix(q, r_grid, 'sphere', contrast=1.0)
        G2 = build_g_matrix(q, r_grid, 'sphere', contrast=2.0)
        np.testing.assert_allclose(G2, 2.0 * G1, rtol=1e-12)

    def test_unknown_shape_raises(self):
        q = np.array([0.01, 0.02])
        r = np.array([100.0])
        with pytest.raises(ValueError, match="Unknown shape"):
            build_g_matrix(q, r, 'banana', contrast=1.0)

    def test_spheroid_shape(self):
        q = np.logspace(-3, -1, 20)
        r_grid = make_r_grid(50, 300, 10)
        G = build_g_matrix(q, r_grid, 'spheroid', contrast=1.0, aspect_ratio=2.0)
        assert G.shape == (20, 10)
        assert np.all(G >= 0)


# ──────────────────────────────────────────────────────────────────────────────
# Bin-width tests
# ──────────────────────────────────────────────────────────────────────────────

class TestBinWidths:
    def test_uniform_grid(self):
        r = np.linspace(10, 100, 10)   # spacing = 10
        dw = bin_widths(r)
        # First and last: dw = spacing
        assert abs(dw[0] - (r[1] - r[0])) < 1e-12
        assert abs(dw[-1] - (r[-1] - r[-2])) < 1e-12
        # Interior: (r[i+1] - r[i-1]) / 2 = spacing
        expected_interior = (r[1] - r[0])   # uniform → all equal
        np.testing.assert_allclose(dw[1:-1], expected_interior, rtol=1e-12)

    def test_positive(self):
        r = make_r_grid(10, 1000, 50)
        dw = bin_widths(r)
        assert np.all(dw > 0)

    def test_single_point(self):
        r = np.array([100.0])
        dw = bin_widths(r)
        assert dw.shape == (1,)
        assert dw[0] == 1.0


# ──────────────────────────────────────────────────────────────────────────────
# SizesDistribution fitting tests
# ──────────────────────────────────────────────────────────────────────────────

def _make_model(**kwargs) -> SizesDistribution:
    s = SizesDistribution()
    s.r_min = 50.0
    s.r_max = 600.0
    s.n_bins = 50
    s.contrast = 1.0
    for k, v in kwargs.items():
        setattr(s, k, v)
    return s


class TestRegularization:
    def test_converges(self):
        q, I, err, _, _ = _synthetic_data()
        s = _make_model(method='regularization')
        result = s.fit(q, I, err)
        assert result['success'], result['message']

    def test_chi_squared_reasonable(self):
        """χ² should be in the ballpark of the number of data points."""
        q, I, err, _, _ = _synthetic_data(noise_frac=0.03)
        s = _make_model(method='regularization')
        result = s.fit(q, I, err)
        M = len(q)
        # Typically chi² / M ∈ [0.1, 10] for a reasonable fit
        assert result['chi_squared'] / M < 50.0

    def test_nonnegative(self):
        q, I, err, _, _ = _synthetic_data()
        s = _make_model(method='regularization')
        result = s.fit(q, I, err)
        assert np.all(result['distribution'] >= 0)

    def test_volume_fraction_positive(self):
        q, I, err, _, _ = _synthetic_data(vf_true=0.01)
        s = _make_model(method='regularization')
        result = s.fit(q, I, err)
        assert result['volume_fraction'] > 0

    def test_peak_near_true_r(self):
        """Distribution peak should be near the true radius (within ±50%)."""
        r_true = 200.0
        q, I, err, r_grid, _ = _synthetic_data(r_true=r_true)
        s = _make_model(method='regularization')
        result = s.fit(q, I, err)
        dist = result['distribution']
        r_peak = result['r_grid'][np.argmax(dist)]
        assert abs(r_peak - r_true) / r_true < 0.5, \
            f"Peak at {r_peak:.1f} Å, expected ~{r_true:.1f} Å"


class TestTNNLS:
    def test_converges(self):
        q, I, err, _, _ = _synthetic_data()
        s = _make_model(method='tnnls')
        result = s.fit(q, I, err)
        assert result['success'], result['message']

    def test_nonnegative(self):
        """All distribution values must be ≥ 0 (TNNLS constraint)."""
        q, I, err, _, _ = _synthetic_data()
        s = _make_model(method='tnnls')
        result = s.fit(q, I, err)
        assert np.all(result['distribution'] >= 0)

    def test_volume_fraction_positive(self):
        q, I, err, _, _ = _synthetic_data()
        s = _make_model(method='tnnls')
        result = s.fit(q, I, err)
        assert result['volume_fraction'] > 0

    def test_model_intensity_shape(self):
        q, I, err, _, _ = _synthetic_data(n_q=60)
        s = _make_model(method='tnnls')
        result = s.fit(q, I, err)
        assert result['model_intensity'].shape == (60,)


class TestMaxEnt:
    def test_converges(self):
        q, I, err, _, _ = _synthetic_data()
        s = _make_model(method='maxent', maxent_max_iter=500)
        result = s.fit(q, I, err)
        assert result['success'], result['message']

    def test_positive(self):
        """MaxEnt solution must be strictly positive."""
        q, I, err, _, _ = _synthetic_data()
        s = _make_model(method='maxent', maxent_max_iter=300)
        result = s.fit(q, I, err)
        assert np.all(result['distribution'] > 0), \
            "MaxEnt distribution must be strictly positive"

    def test_volume_fraction_positive(self):
        q, I, err, _, _ = _synthetic_data()
        s = _make_model(method='maxent')
        result = s.fit(q, I, err)
        assert result['volume_fraction'] > 0


# ──────────────────────────────────────────────────────────────────────────────
# SizesDistribution utility tests
# ──────────────────────────────────────────────────────────────────────────────

class TestUtilities:
    def test_to_from_dict(self):
        s = SizesDistribution()
        s.r_min = 20.0
        s.r_max = 800.0
        s.n_bins = 40
        s.shape = 'spheroid'
        s.shape_params = {'aspect_ratio': 2.5}
        s.method = 'tnnls'
        d = s.to_dict()
        s2 = SizesDistribution.from_dict(d)
        assert s2.r_min == 20.0
        assert s2.shape == 'spheroid'
        assert s2.shape_params == {'aspect_ratio': 2.5}
        assert s2.method == 'tnnls'

    def test_no_error_auto_generated(self):
        """Fit without providing error should succeed (auto 5% error)."""
        q, I, err, _, _ = _synthetic_data()
        s = _make_model(method='regularization')
        result = s.fit(q, I)   # no error provided
        assert result['success']

    def test_calculate_vf_and_rg(self):
        """After a fit, calculate_volume_fraction and calculate_rg return positive values."""
        q, I, err, _, _ = _synthetic_data()
        s = _make_model(method='regularization')
        s.fit(q, I, err)
        assert s.calculate_volume_fraction() > 0
        assert s.calculate_rg() > 0


# ──────────────────────────────────────────────────────────────────────────────
# HDF5 round-trip test
# ──────────────────────────────────────────────────────────────────────────────

class TestHDF5IO:
    def test_save_load_round_trip(self):
        """Save sizes results to HDF5, reload, check arrays match."""
        from pyirena.io.nxcansas_sizes import save_sizes_results, load_sizes_results
        from pyirena.io.nxcansas_unified import create_nxcansas_file

        q, I, err, _, _ = _synthetic_data(n_q=50)
        s = _make_model(method='regularization')
        result = s.fit(q, I, err)
        assert result['success']

        with tempfile.TemporaryDirectory() as td:
            fp = Path(td) / 'test_sizes.h5'
            # Create base NXcanSAS file
            create_nxcansas_file(fp, q, I, err, sample_name='test')

            params = {
                'shape':          s.shape,
                'method':         s.method,
                'contrast':       s.contrast,
                'chi_squared':    result['chi_squared'],
                'volume_fraction': result['volume_fraction'],
                'rg':             result['rg'],
                'n_iterations':   result['n_iterations'],
            }
            save_sizes_results(
                fp,
                q=q,
                intensity_data=I,
                intensity_model=result['model_intensity'],
                residuals=result['residuals'],
                r_grid=result['r_grid'],
                distribution=result['distribution'],
                params=params,
                intensity_error=err,
            )

            loaded = load_sizes_results(fp)

        np.testing.assert_allclose(loaded['Q'],            q,                        rtol=1e-9)
        np.testing.assert_allclose(loaded['distribution'], result['distribution'],   rtol=1e-9)
        np.testing.assert_allclose(loaded['r_grid'],       result['r_grid'],          rtol=1e-9)
        assert loaded['shape']  == s.shape
        assert loaded['method'] == s.method
        assert abs(loaded['chi_squared'] - result['chi_squared']) < 1e-6


# ──────────────────────────────────────────────────────────────────────────────
# Q-range → radius-grid bounds helper
# ──────────────────────────────────────────────────────────────────────────────

class TestRBoundsFromQRange:
    def test_exact_pi_over_q(self):
        from pyirena.core.form_factors import r_bounds_from_q_range
        r_min, r_max = r_bounds_from_q_range(1e-3, 0.1, pad=False)
        assert r_min == pytest.approx(np.pi / 0.1, rel=1e-9)
        assert r_max == pytest.approx(np.pi / 1e-3, rel=1e-9)

    def test_padded_brackets_exact(self):
        """Padded bounds must sit outside (or on) the exact π/Q range."""
        from pyirena.core.form_factors import r_bounds_from_q_range
        r_min_e, r_max_e = r_bounds_from_q_range(2.85e-4, 1.085e-2, pad=False)
        r_min_p, r_max_p = r_bounds_from_q_range(2.85e-4, 1.085e-2, pad=True)
        assert r_min_p <= r_min_e
        assert r_max_p >= r_max_e
        # nice 1-2-5 numbers
        assert r_min_p == pytest.approx(200.0)
        assert r_max_p == pytest.approx(20000.0)

    def test_swapped_inputs_ok(self):
        from pyirena.core.form_factors import r_bounds_from_q_range
        a = r_bounds_from_q_range(0.1, 1e-3, pad=True)
        b = r_bounds_from_q_range(1e-3, 0.1, pad=True)
        assert a == b
        assert a[0] < a[1]


# ──────────────────────────────────────────────────────────────────────────────
# recommend_sizes_setup — suitability for the common data classes
# ──────────────────────────────────────────────────────────────────────────────

def _broad_data(power_law_B=3e-9, power_law_P=4.0, flat=0.6,
                r0=500.0, log_sigma=0.5, vf=0.01, contrast=100.0,
                q_min=3e-4, q_max=0.3, n_q=400, noise=0.02, seed=0):
    """Broad log-normal size distribution on a power-law + flat background."""
    rng = np.random.default_rng(seed)
    q = np.logspace(np.log10(q_min), np.log10(q_max), n_q)
    rg = make_r_grid(20.0, 6000.0, 120, log_spacing=True)
    dw = bin_widths(rg)
    Dv = (1.0 / (rg * log_sigma * np.sqrt(2 * np.pi))
          * np.exp(-(np.log(rg / r0)) ** 2 / (2 * log_sigma ** 2)))
    Dv *= vf / float(np.trapezoid(Dv, rg))
    G = build_g_matrix(q, rg, 'sphere', contrast)
    I = G @ (Dv * dw) + power_law_B * q ** (-power_law_P) + flat
    err = noise * I
    I = I + err * rng.standard_normal(n_q)
    return q, I, err


class TestRecommendSizesSetup:
    def test_broad_distribution_on_background_is_suitable(self):
        """Broad distribution + power-law + flat bg — the everyday case."""
        from pyirena.core.sizes import recommend_sizes_setup
        q, I, err = _broad_data()
        rec = recommend_sizes_setup(q, I, sigma=err)
        assert rec["suitable"] is True
        r = rec["recommended"]
        # inversion window is a proper sub-range with the flat bg recovered
        assert r["inversion_q_min"] < r["inversion_q_max"]
        assert r["r_min"] is not None and r["r_max"] is not None
        assert r["r_min"] < r["r_max"]
        assert r["flat_background"] == pytest.approx(0.6, rel=0.25)
        assert rec["features"]["signal_to_bg_ratio"] is not None

    def test_solid_material_floating_slope_is_suitable(self):
        from pyirena.core.sizes import recommend_sizes_setup
        q, I, err = _broad_data(power_law_B=1e-7, power_law_P=3.5,
                                flat=0.3, r0=300.0, seed=2)
        rec = recommend_sizes_setup(q, I, sigma=err)
        assert rec["suitable"] is True

    def test_pure_background_is_not_suitable(self):
        """Flat noise with no particle signal must be flagged unsuitable."""
        from pyirena.core.sizes import recommend_sizes_setup
        rng = np.random.default_rng(1)
        q = np.logspace(-3, -0.5, 300)
        I = 0.5 + 0.5 * 0.02 * rng.standard_normal(300)
        err = 0.02 * np.abs(I)
        rec = recommend_sizes_setup(q, I, sigma=err)
        assert rec["suitable"] is False
        assert rec["warnings"]

    def test_multiple_knees_are_advisory_not_fatal(self):
        """A broad distribution that reads as several 'levels' is still suitable."""
        from pyirena.core.sizes import recommend_sizes_setup
        q, I, err = _broad_data(seed=5)
        rec = recommend_sizes_setup(q, I, sigma=err)
        # even if the detector reports >=3 levels, suitability holds
        assert rec["suitable"] is True


# ──────────────────────────────────────────────────────────────────────────────
# Monte Carlo — no r³ over-weighting; contributions stay in the resolvable band
# ──────────────────────────────────────────────────────────────────────────────

class TestMonteCarloWiring:
    def test_distribution_not_r3_shifted(self):
        """MC volume distribution must peak near the true size, not ~2× larger.

        The old r³ re-weighting bug pushed the volume-weighted mean radius to
        roughly twice the truth; the resolvable-band restriction also keeps the
        peak off the smallest bin.  Bounds are loose to stay robust to the RNG.
        """
        q, I, err, r_grid_true, dist_true = _synthetic_data(
            r_true=200.0, r_sigma=40.0, q_min=1e-3, q_max=0.1,
            r_min=50.0, r_max=600.0, n_bins=60,
        )
        s = SizesDistribution()
        s.r_min, s.r_max, s.n_bins, s.log_spacing = 50.0, 600.0, 60, True
        s.shape, s.contrast, s.method = 'sphere', 1.0, 'montecarlo'
        s.montecarlo_n_repetitions = 8
        s.montecarlo_max_iter = 40000
        res = s.fit(q, I, err)
        assert res["success"]
        rg = res["r_grid"]
        D = res["distribution"]
        vf = res["volume_fraction"]
        vmean = float(np.trapezoid(rg * D, rg) / vf)
        # true volume-mean radius ~200; the r³ bug gave ~2× (>350).
        assert 120.0 < vmean < 330.0
        # peak not pinned at the smallest bin
        assert rg[int(np.argmax(D))] > rg[0]
