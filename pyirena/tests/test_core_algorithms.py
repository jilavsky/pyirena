"""
Tranche-2 core tests: data_manipulation, data_merge, similarity,
morphology, waxs_peakfit peak math, fractal growth, and the simple-fit
model registry.

Optional-dependency modules (scattering_contrast -> xraydb/periodictable,
diffraction_lines -> Dans_Diffraction) are covered in
test_optional_dep_modules.py with import skips.
"""

import numpy as np
import pytest

from pyirena.core.data_manipulation import (
    DataManipulation, ScaleConfig, TrimConfig, RebinConfig, SubtractConfig,
)


@pytest.fixture
def synth():
    q = np.logspace(-3, -0.5, 300)
    I = 100.0 * np.exp(-q**2 * 50.0**2 / 3.0) + 1e-4 * q**-4 + 0.5
    dI = 0.02 * I
    return q, I, dI


# ---------------------------------------------------------------------------
# Data manipulation
# ---------------------------------------------------------------------------

class TestScale:
    def test_scale_and_background(self, synth):
        q, I, dI = synth
        res = DataManipulation.scale(q, I, dI, None, ScaleConfig(scale_I=2.0, background=0.1))
        np.testing.assert_allclose(res.I, 2.0 * I - 0.1)
        np.testing.assert_allclose(res.q, q)

    def test_uncertainty_scales_with_intensity(self, synth):
        q, I, dI = synth
        res = DataManipulation.scale(q, I, dI, None, ScaleConfig(scale_I=3.0))
        np.testing.assert_allclose(res.dI, 3.0 * dI)

    def test_identity(self, synth):
        q, I, dI = synth
        res = DataManipulation.scale(q, I, dI, None, ScaleConfig())
        np.testing.assert_allclose(res.I, I)
        np.testing.assert_allclose(res.dI, dI)


class TestTrim:
    def test_limits_respected(self, synth):
        q, I, dI = synth
        res = DataManipulation.trim(q, I, dI, None, TrimConfig(q_min=0.01, q_max=0.1))
        assert res.q.min() >= 0.01 and res.q.max() <= 0.1
        assert len(res.q) < len(q)

    def test_values_preserved_inside_window(self, synth):
        q, I, dI = synth
        res = DataManipulation.trim(q, I, dI, None, TrimConfig(q_min=0.01, q_max=0.1))
        mask = (q >= 0.01) & (q <= 0.1)
        np.testing.assert_allclose(res.I, I[mask])

    def test_dq_trimmed_alongside(self, synth):
        q, I, dI = synth
        dQ = 0.01 * q
        res = DataManipulation.trim(q, I, dI, dQ, TrimConfig(q_min=0.01, q_max=0.1))
        assert len(res.dQ) == len(res.q)


class TestRebin:
    def test_point_count_and_range(self, synth):
        q, I, dI = synth
        res = DataManipulation.rebin(q, I, dI, None, RebinConfig(n_points=50))
        assert len(res.q) <= 50
        assert res.q.min() >= q.min() * 0.999
        assert res.q.max() <= q.max() * 1.001

    def test_smooth_data_preserved(self, synth):
        """Rebinning smooth data must reproduce it (interpolation sanity)."""
        q, I, dI = synth
        res = DataManipulation.rebin(q, I, dI, None, RebinConfig(n_points=100))
        expected = np.interp(res.q, q, I)
        np.testing.assert_allclose(res.I, expected, rtol=0.05)

    def test_monotonic_output(self, synth):
        q, I, dI = synth
        res = DataManipulation.rebin(q, I, dI, None, RebinConfig(n_points=80))
        assert (np.diff(res.q) > 0).all()


class TestSubtract:
    def test_exact_subtraction_on_common_grid(self, synth):
        q, I, dI = synth
        buffer_I = np.full_like(q, 0.5)
        res = DataManipulation.subtract(
            q, I, dI, None, q, buffer_I, 0.02 * buffer_I,
            SubtractConfig(buffer_scale=1.0),
        )
        np.testing.assert_allclose(res.I, I - 0.5, rtol=1e-10)

    def test_buffer_scale_applied(self, synth):
        q, I, dI = synth
        buffer_I = np.full_like(q, 0.5)
        res = DataManipulation.subtract(
            q, I, dI, None, q, buffer_I, 0.02 * buffer_I,
            SubtractConfig(buffer_scale=2.0),
        )
        np.testing.assert_allclose(res.I, I - 1.0, rtol=1e-10)

    def test_uncertainties_add_in_quadrature(self, synth):
        q, I, dI = synth
        buffer_I = np.full_like(q, 0.5)
        dbuf = 0.1 * buffer_I
        res = DataManipulation.subtract(
            q, I, dI, None, q, buffer_I, dbuf,
            SubtractConfig(buffer_scale=1.0),
        )
        np.testing.assert_allclose(res.dI, np.sqrt(dI**2 + dbuf**2), rtol=1e-6)


class TestAverage:
    def test_average_of_identical_datasets(self, synth):
        q, I, dI = synth
        ds = [(q, I, dI, None) for _ in range(3)]
        res = DataManipulation.average(ds)
        np.testing.assert_allclose(res.I, I, rtol=1e-10)

    def test_average_of_two_levels(self, synth):
        q, I, dI = synth
        ds = [(q, I, dI, None), (q, 3.0 * I, dI, None)]
        res = DataManipulation.average(ds)
        np.testing.assert_allclose(res.I, 2.0 * I, rtol=1e-6)


# ---------------------------------------------------------------------------
# Data merge
# ---------------------------------------------------------------------------

class TestDataMerge:
    def _two_overlapping(self):
        """DS1 (low q) and DS2 (high q) from the same underlying curve,
        with DS2 mis-scaled by a known factor."""
        def curve(q):
            return 1000.0 * np.exp(-q**2 * 100.0**2 / 3.0) + 1e-3 * q**-4
        q1 = np.logspace(-3, -1.3, 150)          # 0.001 – 0.05
        q2 = np.logspace(-1.7, -0.3, 150)        # 0.02  – 0.5   (overlap 0.02–0.05)
        I1 = curve(q1)
        true_scale = 2.5                          # DS2 must be divided by this
        I2 = curve(q2) * true_scale
        return q1, I1, 0.01 * I1, q2, I2, 0.01 * I2, true_scale

    def test_optimize_recovers_scale(self):
        from pyirena.core.data_merge import DataMerge, MergeConfig
        q1, I1, dI1, q2, I2, dI2, true_scale = self._two_overlapping()
        cfg = MergeConfig(q_overlap_min=0.02, q_overlap_max=0.05,
                          fit_scale=True, scale_dataset="DS2", fit_qshift=False)
        result = DataMerge().optimize(q1, I1, dI1, q2, I2, dI2, cfg)
        assert result.success
        assert result.scale == pytest.approx(1.0 / true_scale, rel=0.02)

    def test_merged_output_is_monotonic_and_continuous(self):
        from pyirena.core.data_merge import DataMerge, MergeConfig
        q1, I1, dI1, q2, I2, dI2, true_scale = self._two_overlapping()
        cfg = MergeConfig(q_overlap_min=0.02, q_overlap_max=0.05,
                          fit_scale=True, scale_dataset="DS2", fit_qshift=False)
        dm = DataMerge()
        result = dm.optimize(q1, I1, dI1, q2, I2, dI2, cfg)
        qm, Im, dIm, dQm = dm.merge(q1, I1, dI1, None, q2, I2, dI2, None,
                                    result, cfg)
        # sorted; duplicate q values may occur where both datasets contribute
        assert (np.diff(qm) >= 0).all()
        assert qm.min() == pytest.approx(q1.min())
        assert qm.max() == pytest.approx(q2.max())
        # No scale jump: compare merged intensity against the true curve
        truth = 1000.0 * np.exp(-qm**2 * 100.0**2 / 3.0) + 1e-3 * qm**-4
        np.testing.assert_allclose(Im, truth, rtol=0.05)


# ---------------------------------------------------------------------------
# Similarity (CorMap)
# ---------------------------------------------------------------------------

class TestCormapStatistics:
    def test_p_value_exact_small_case(self):
        """n=4, c=2: 16 sequences, 2 alternating -> P(run>=2) = 14/16."""
        from pyirena.core.similarity import _cormap_p_value
        assert _cormap_p_value(4, 2) == pytest.approx(14.0 / 16.0)

    def test_no_recursion_error_for_long_datasets(self):
        """Regression: Schilling recurrence was recursive and blew the
        stack for n >~ 300 — an ordinary SAXS dataset length."""
        from pyirena.core.similarity import _cormap_p_value
        p = _cormap_p_value(2000, 10)
        assert 0.0 < p < 1.0

    def test_no_float_overflow_for_large_n(self):
        """Regression: float(delta) overflowed for n >~ 1024 and silently
        returned p=0, flagging long datasets as different regardless of
        the data."""
        from pyirena.core.similarity import _cormap_p_value
        # A short run in a long dataset is overwhelmingly probable
        assert _cormap_p_value(5000, 8) > 0.9

    def test_p_monotone_decreasing_in_run_length(self):
        from pyirena.core.similarity import _cormap_p_value
        ps = [_cormap_p_value(500, c) for c in (3, 6, 12, 25, 50)]
        assert all(a >= b for a, b in zip(ps, ps[1:]))


class TestSimilarity:
    def test_identical_datasets_accepted(self, synth):
        from pyirena.core.similarity import check_similarity
        q, I, dI = synth
        rng = np.random.default_rng(1)
        noisy = [I * (1.0 + 0.01 * rng.standard_normal(len(q))) for _ in range(3)]
        datasets = [(q, In, dI) for In in noisy]
        results = check_similarity(datasets, method="cormap")
        assert all(r.accepted for r in results)

    def test_outlier_rejected(self, synth):
        from pyirena.core.similarity import check_similarity
        q, I, dI = synth
        rng = np.random.default_rng(2)
        good = [I * (1.0 + 0.01 * rng.standard_normal(len(q))) for _ in range(3)]
        # Systematically different shape — steeper decay — must be flagged
        outlier = I * (q / q.min()) ** -0.5
        datasets = [(q, In, dI) for In in good + [outlier]]
        results = check_similarity(datasets, method="cormap", reference="first")
        assert all(r.accepted for r in results[:3])
        assert results[3].accepted is False


# ---------------------------------------------------------------------------
# Morphology metrics
# ---------------------------------------------------------------------------

class TestMorphology:
    def test_single_percolating_channel(self):
        from pyirena.core.morphology import compute_morphology_metrics
        # 10^3 solid block with a single open channel along x
        vox = np.ones((10, 10, 10), dtype=np.uint8)
        vox[:, 5, 5] = 0
        m = compute_morphology_metrics(vox, voxel_pitch_A=10.0)
        assert m.minority_phase_value == 0
        assert m.minority_volume_fraction == pytest.approx(10.0 / 1000.0)
        assert m.n_clusters == 1
        assert m.percolating_x
        assert not m.percolating_y and not m.percolating_z
        assert m.open_porosity_fraction == pytest.approx(1.0)

    def test_isolated_pore_is_closed(self):
        from pyirena.core.morphology import compute_morphology_metrics
        vox = np.ones((10, 10, 10), dtype=np.uint8)
        vox[4:6, 4:6, 4:6] = 0          # 8-voxel isolated pore
        m = compute_morphology_metrics(vox, voxel_pitch_A=10.0)
        assert m.n_clusters == 1
        assert not (m.percolating_x or m.percolating_y or m.percolating_z)
        assert m.closed_porosity_fraction == pytest.approx(1.0)
        assert m.open_porosity_fraction == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# WAXS peak math
# ---------------------------------------------------------------------------

class TestWaxsPeakMath:
    Q = np.linspace(1.0, 3.0, 4001)

    @pytest.mark.parametrize("shape", ["Gauss", "Lorentz", "Pseudo-Voigt", "LogNormal"])
    def test_peak_height_and_position(self, shape):
        from pyirena.core.waxs_peakfit import eval_peak
        params = {"A": 50.0, "Q0": 2.0, "FWHM": 0.05, "eta": 0.5}
        y = eval_peak(self.Q, shape, params)
        assert y.max() == pytest.approx(50.0, rel=1e-3), shape
        assert self.Q[np.argmax(y)] == pytest.approx(2.0, abs=0.01), shape

    @pytest.mark.parametrize("shape", ["Gauss", "Lorentz"])
    def test_fwhm_definition(self, shape):
        from pyirena.core.waxs_peakfit import eval_peak
        fwhm = 0.05
        y = eval_peak(self.Q, shape, {"A": 1.0, "Q0": 2.0, "FWHM": fwhm})
        above = self.Q[y >= 0.5]
        assert above[-1] - above[0] == pytest.approx(fwhm, rel=0.02), shape

    def test_gauss_area_analytic(self):
        """Gauss area = A * FWHM * sqrt(pi/(4 ln2)); numeric and helper agree."""
        from pyirena.core.waxs_peakfit import eval_peak, peak_area
        A, fwhm = 50.0, 0.05
        params = {"A": A, "Q0": 2.0, "FWHM": fwhm}
        analytic = A * fwhm * np.sqrt(np.pi / (4.0 * np.log(2.0)))
        numeric = np.trapezoid(eval_peak(self.Q, "Gauss", params), self.Q)
        assert numeric == pytest.approx(analytic, rel=1e-3)
        assert peak_area("Gauss", params) == pytest.approx(analytic, rel=1e-3)

    def test_constant_background(self):
        from pyirena.core.waxs_peakfit import eval_background
        y = eval_background(self.Q, "Constant", [1.5])
        np.testing.assert_allclose(y, 1.5)

    def test_snip_background_below_peaks(self):
        """SNIP background must stay at/below the data through a peak."""
        from pyirena.core.waxs_peakfit import eval_peak, snip_background
        base = 10.0 + 2.0 * self.Q
        data = base + eval_peak(self.Q, "Gauss", {"A": 100.0, "Q0": 2.0, "FWHM": 0.05})
        bg = snip_background(data, half_width_frac=0.05)
        peak_zone = np.abs(self.Q - 2.0) < 0.02
        assert (bg[peak_zone] < data[peak_zone]).all()
        # Away from the peak, background tracks the data closely
        flat_zone = np.abs(self.Q - 1.2) < 0.1
        np.testing.assert_allclose(bg[flat_zone], data[flat_zone], rtol=0.1)


# ---------------------------------------------------------------------------
# Fractal growth
# ---------------------------------------------------------------------------

class TestFractalGrowth:
    def test_growth_reproducible_with_seed(self):
        from pyirena.core.fractals import grow_aggregate, GrowthConfig
        a1 = grow_aggregate(GrowthConfig(z=40, seed=7))
        a2 = grow_aggregate(GrowthConfig(z=40, seed=7))
        np.testing.assert_array_equal(a1.positions, a2.positions)

    def test_aggregate_is_connected_and_sized(self):
        from pyirena.core.fractals import grow_aggregate, GrowthConfig
        z = 60
        agg = grow_aggregate(GrowthConfig(z=z, seed=3))
        assert agg.positions.shape == (z, 3)
        # every particle after the seed has at least one neighbor
        assert (agg.neighbor_count[1:] > 0).all()
        # positions are unique lattice sites
        assert len(np.unique(agg.positions, axis=0)) == z

    def test_params_physical_ranges(self):
        from pyirena.core.fractals import grow_aggregate, GrowthConfig
        agg = grow_aggregate(GrowthConfig(z=80, seed=11))
        p = agg.params
        assert p.z == 80
        assert 1.0 <= p.dmin <= 3.0
        assert p.rg_aggregate > p.rg_primary > 0

    def test_df_equals_dmin_times_c(self):
        """Beaucage relation d_f = d_min * c must hold by construction."""
        from pyirena.core.fractals import grow_aggregate, GrowthConfig
        p = grow_aggregate(GrowthConfig(z=80, seed=11)).params
        assert p.df == pytest.approx(p.dmin * p.c, rel=1e-10)

    def test_compute_params_consistent_on_same_aggregate(self):
        """Re-evaluating the same aggregate must give compatible dimensions
        (path sampling is stochastic, so allow a loose tolerance)."""
        from pyirena.core.fractals import (
            grow_aggregate, GrowthConfig, compute_fractal_params,
        )
        agg = grow_aggregate(GrowthConfig(z=80, seed=11))
        p2 = compute_fractal_params(agg.positions, agg.neighbor_list,
                                    agg.neighbor_count, rg_primary=10.0,
                                    num_test_paths=2500)
        assert p2.df == pytest.approx(agg.params.df, rel=0.15)
        assert p2.rg_aggregate == pytest.approx(agg.params.rg_aggregate, rel=1e-6)

    def test_tiny_aggregate_returns_nan_params(self):
        from pyirena.core.fractals import compute_fractal_params
        pos = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int32)
        nl = np.full((2, 26), -1, dtype=np.int32)
        nc = np.array([1, 1], dtype=np.uint8)
        p = compute_fractal_params(pos, nl, nc, rg_primary=10.0)
        assert np.isnan(p.df)
        assert p.z == 2


class TestFractalIntensity:
    def _params(self):
        from pyirena.core.fractals import grow_aggregate, GrowthConfig
        return grow_aggregate(GrowthConfig(z=80, seed=11)).params

    def test_intensity_finite_positive(self):
        from pyirena.core.fractals import intensity_unified
        p = self._params()
        q = np.logspace(-3, 0, 200)
        I = intensity_unified(p, q)
        assert np.isfinite(I).all()
        assert (I > 0).all()

    def test_low_q_plateau_equals_z(self):
        """Guinier plateau of the aggregate level: I(q->0) -> G = z (+1
        from the primary level)."""
        from pyirena.core.fractals import intensity_unified
        p = self._params()
        I0 = intensity_unified(p, np.array([1e-4]))[0]
        assert I0 == pytest.approx(p.z + 1, rel=0.05)

    def test_high_q_porod_slope(self):
        """Beyond the primary-particle knee the slope must approach -4."""
        from pyirena.core.fractals import intensity_unified
        p = self._params()
        q = np.logspace(0.0, 0.7, 50)   # q*Rg_primary >> 1  (Rg_prim = 10)
        I = intensity_unified(p, q)
        slope = np.polyfit(np.log(q), np.log(I), 1)[0]
        assert slope == pytest.approx(-4.0, abs=0.3)


# ---------------------------------------------------------------------------
# Simple-fit model registry (beyond the Guinier already covered)
# ---------------------------------------------------------------------------

class TestSimpleFitModels:
    def test_registry_names_consistent(self):
        from pyirena.core.simple_fits import MODEL_REGISTRY, MODEL_NAMES
        assert set(MODEL_NAMES) == set(MODEL_REGISTRY)

    def test_porod_fit_recovers_slope(self):
        from pyirena.core.simple_fits import SimpleFitModel
        q = np.logspace(-2, -0.5, 100)
        I = 5e-5 * q ** -4.0
        m = SimpleFitModel()
        m.model = "Power Law"
        m._reset_to_defaults()
        result = m.fit(q, I, 0.01 * I)
        # slope parameter should be ~4 whichever sign convention is used
        slopes = [v for k, v in result["params"].items()
                  if "slope" in k.lower() or "exponent" in k.lower()]
        assert slopes, f"no slope-like parameter in {list(result['params'])}"
        assert abs(slopes[0]) == pytest.approx(4.0, rel=0.02)

    def test_sphere_fit_recovers_radius(self):
        from pyirena.core.simple_fits import SimpleFitModel
        from pyirena.core.form_factors import sphere_ff
        r_true = 60.0
        q = np.logspace(-2.5, -1.0, 150)
        I = 1e-2 * sphere_ff(q, r_true)
        m = SimpleFitModel()
        m.model = "Sphere"
        m._reset_to_defaults()
        result = m.fit(q, I, 0.01 * I)
        radii = [v for k, v in result["params"].items() if k.lower().startswith("r")]
        assert radii and radii[0] == pytest.approx(r_true, rel=0.05)
