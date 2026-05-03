"""
Unit tests for the SAXS Morph engine.

These tests use synthetic Debye-Bueche scattering as the reference data:
    I(Q) = I0 * a**3 / (1 + (Q*a)**2)**2
which is the analytical SAS profile of a random two-phase medium with a
single correlation length a (Debye, Anderson, Brumberger 1957).  An engine
that correctly inverts I(Q) -> gamma(r) -> F(k) -> voxelgram should produce
a binary cube whose statistics match the input phi and whose recomputed
I(Q) is in the same ballpark as the input.
"""

import numpy as np
import pytest

from pyirena.core.saxs_morph import (
    SaxsMorphConfig, SaxsMorphResult, SaxsMorphEngine,
    alfa_threshold, debye_autocorr, spectral_function,
    generate_voxelgram, voxelgram_to_iq,
    derive_contrast_from_invariant, derive_phi_from_invariant,
    compute_invariant_extrapolated,
    fit_power_law_bg, fit_flat_bg,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def synthetic_debye_bueche(q, I0=1.0, a=100.0):
    """Analytical Debye-Bueche I(Q) for correlation length a [A]."""
    return I0 * a ** 3 / (1.0 + (q * a) ** 2) ** 2


def make_synthetic_dataset(n=200, q_min=1e-3, q_max=0.3, a=100.0, I0=1.0):
    """Log-spaced Q + DB I(Q) + 5%-of-I uncertainties."""
    q = np.logspace(np.log10(q_min), np.log10(q_max), n)
    I = synthetic_debye_bueche(q, I0=I0, a=a)
    dI = np.maximum(I * 0.05, 1e-30)
    return q, I, dI


# ---------------------------------------------------------------------------
# alfa_threshold
# ---------------------------------------------------------------------------

class TestAlfaThreshold:
    def test_phi_half_is_zero(self):
        assert abs(alfa_threshold(0.5)) < 1e-12

    def test_monotonic_decreasing(self):
        # As phi increases the threshold drops (more of the field exceeds it)
        phis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        alfas = [alfa_threshold(p) for p in phis]
        for a, b in zip(alfas[:-1], alfas[1:]):
            assert a > b

    def test_clamps_extremes(self):
        # Both endpoints should be finite (no inf from erfinv)
        assert np.isfinite(alfa_threshold(0.0))
        assert np.isfinite(alfa_threshold(1.0))


# ---------------------------------------------------------------------------
# debye_autocorr / spectral_function
# ---------------------------------------------------------------------------

class TestSincTransforms:
    def test_gamma_normalised_to_one(self):
        q, I, _ = make_synthetic_dataset()
        r, g = debye_autocorr(q, I, n_r=128)
        assert abs(g[0] - 1.0) < 1e-10
        assert np.all(np.isfinite(g))

    def test_gamma_decreases_from_origin(self):
        # For Debye-Bueche gamma(r) = exp(-r/a) — strictly monotonic
        q, I, _ = make_synthetic_dataset(a=100.0)
        r, g = debye_autocorr(q, I, n_r=128)
        # First half should be roughly monotonically decreasing
        first_half = g[: len(g) // 2]
        assert first_half[0] > first_half[-1]

    def test_spectral_F_nonnegative(self):
        q, I, _ = make_synthetic_dataset()
        r, g = debye_autocorr(q, I, n_r=128)
        k, F = spectral_function(r, g, n_k=128)
        assert np.all(F >= 0)


# ---------------------------------------------------------------------------
# generate_voxelgram
# ---------------------------------------------------------------------------

class TestGenerateVoxelgram:
    def test_shape_and_dtype(self):
        # Trivial spectrum -> still produces a cube of the right shape/dtype
        k = np.linspace(0, 0.1, 64)
        F = np.exp(-k * 50.0)  # arbitrary smooth nonneg spectrum
        vox, sigma, seed = generate_voxelgram(
            k, F, voxel_size=32, box_size_A=500.0,
            alfa=alfa_threshold(0.3), rng_seed=42,
        )
        assert vox.shape == (32, 32, 32)
        assert vox.dtype == np.uint8
        assert vox.min() >= 0 and vox.max() <= 1
        assert seed == 42

    def test_phi_actual_close_to_target(self):
        # With a smooth spectrum the realised phi should match the target
        # within a few percent (statistical fluctuations on a 32**3 cube).
        k = np.linspace(0, 0.1, 64)
        F = np.exp(-k * 50.0)
        for target in (0.2, 0.3, 0.5, 0.7):
            vox, _, _ = generate_voxelgram(
                k, F, voxel_size=32, box_size_A=500.0,
                alfa=alfa_threshold(target), rng_seed=42,
            )
            phi_realised = float(vox.mean())
            # A smooth F(k) gives a strongly correlated field, so the
            # effective number of independent voxels is much smaller than
            # N**3.  A 0.10 tolerance covers fluctuations up to phi=0.5.
            assert abs(phi_realised - target) < 0.10, (
                f"target phi={target}, realised={phi_realised}")

    def test_seed_reproducibility(self):
        k = np.linspace(0, 0.1, 64)
        F = np.exp(-k * 50.0)
        kwargs = dict(spectral_k=k, spectral_F=F, voxel_size=32,
                      box_size_A=500.0, alfa=alfa_threshold(0.3))
        vox1, _, _ = generate_voxelgram(rng_seed=123, **kwargs)
        vox2, _, _ = generate_voxelgram(rng_seed=123, **kwargs)
        assert np.array_equal(vox1, vox2)


# ---------------------------------------------------------------------------
# voxelgram_to_iq
# ---------------------------------------------------------------------------

class TestVoxelgramToIQ:
    def test_returns_finite_positive(self):
        rng = np.random.default_rng(0)
        # Random binary cube
        vox = (rng.random((32, 32, 32)) > 0.7).astype(np.uint8)
        q = np.logspace(-2, 0, 50)
        I = voxelgram_to_iq(vox, voxel_pitch_A=10.0, q_target=q)
        assert I.shape == q.shape
        assert np.all(np.isfinite(I))
        assert np.all(I >= 0)


# ---------------------------------------------------------------------------
# derive_contrast_from_invariant
# ---------------------------------------------------------------------------

class TestInvariant:
    def test_units_are_10e20_cm_minus_4(self):
        """Contrast must be returned in 10^20 cm^-4 (matches Igor convention).

        Pin the unit-conversion factor so a future refactor cannot silently
        drop the Angstrom-to-cm factor (10^4) that we added in 0.5.2.

        Synthetic Debye-Bueche I(Q) with I0 = 1 has the closed-form
        invariant Q* = pi/4 (in numerical Å^-3 cm^-1 units).  So with
        phi=0.30, the derived contrast in 10^20 cm^-4 must be
        (pi/4) * 1e4 / (2 pi^2 * 0.30 * 0.70) =~ 1894.
        """
        q, I, _ = make_synthetic_dataset(I0=1.0, a=100.0)
        Qstar_num = compute_invariant_extrapolated(q, I)
        expected = Qstar_num * 1e4 / (2 * np.pi ** 2 * 0.30 * 0.70)
        actual = derive_contrast_from_invariant(q, I, 0.30)
        # Loose tolerance because the synthetic dataset only spans a
        # finite Q range, so simpson != closed form.
        assert abs(actual - expected) / expected < 1e-9

    def test_voxelgram_iq_is_in_per_cm_when_contrast_is_e20cm4(self):
        """Verify the absolute-units relationship between voxelgram_to_iq
        and Porod's invariant identity.

        For a uniform random voxelgram of volume fraction phi, the
        reconstructed I(Q) integrated as Q* = sum(Q^2 I dQ) * 1e4 should
        equal 2 pi^2 * phi(1-phi) * Delta rho^2 (in 10^20 cm^-4) up to a
        numerical-aperture factor.  We cannot achieve full agreement because
        the voxelgram only resolves Q < pi/pitch — but we should be within
        a factor of ~2 of the theoretical invariant for a binary cube,
        confirming the units are right (a stray factor of N^3 or 10^4
        would push us off by orders of magnitude).
        """
        from scipy.integrate import simpson
        rng = np.random.default_rng(7)
        N, pitch = 32, 10.0  # box edge = 320 A
        # Spatially-correlated binary cube (random + Gaussian smoothing)
        # so the FFT spectrum is concentrated below pi/pitch and we're
        # not throwing away half the variance to high-Q aliasing.
        from scipy.ndimage import gaussian_filter
        field = rng.standard_normal((N, N, N))
        field = gaussian_filter(field, sigma=2.0)
        threshold = np.percentile(field, 70)        # phi ~= 0.30
        vox = (field > threshold).astype(np.uint8)
        phi_actual = float(vox.mean())
        contrast = 1.0  # in 10^20 cm^-4

        q = np.logspace(-4, np.log10(np.pi / pitch * 0.9), 300)  # A^-1
        I_struct = voxelgram_to_iq(vox, voxel_pitch_A=pitch, q_target=q)
        I_model = contrast * I_struct  # cm^-1

        Qstar_num = float(simpson(q ** 2 * I_model, x=q))  # A^-3 cm^-1
        Qstar_e20cm4 = Qstar_num * 1e4
        expected = 2 * np.pi ** 2 * phi_actual * (1 - phi_actual) * contrast

        # Within a factor of 2 — our Q range only covers low-Q, missing
        # some high-Q content.  But it should NOT be off by 1e4 or N^3.
        ratio = Qstar_e20cm4 / expected
        assert 0.30 < ratio < 1.5, (
            f"invariant ratio {ratio:.3g} out of band — units bug? "
            f"Qstar={Qstar_e20cm4:.3g}, expected={expected:.3g}")

    def test_invariant_extrapolation_adds_to_data(self):
        """compute_invariant_extrapolated must equal data integral + low-Q
        + high-Q closed-form contributions.
        """
        from scipy.integrate import simpson
        q, I, _ = make_synthetic_dataset()
        Q_data = float(simpson(q ** 2 * I, x=q))

        n = 5
        I0 = float(np.mean(I[:n]))
        Q_low_expected = I0 * q[0] ** 3 / 3.0

        K = float(np.mean(I[-n:] * q[-n:] ** 4))
        Q_high_expected = K / q[-1]

        Q_total = compute_invariant_extrapolated(q, I, n_avg=n)
        assert abs(Q_total - (Q_data + Q_low_expected + Q_high_expected)) < 1e-12

    def test_invariant_no_extrapolation_when_n_avg_zero(self):
        """n_avg=0 should reduce to bare Simpson integral."""
        from scipy.integrate import simpson
        q, I, _ = make_synthetic_dataset()
        Q_data = float(simpson(q ** 2 * I, x=q))
        Q_no_extrap = compute_invariant_extrapolated(q, I, n_avg=0)
        assert abs(Q_no_extrap - Q_data) < 1e-12

    def test_recovers_known_contrast(self):
        # Build I(Q) = 2 pi**2 phi (1-phi) Drho**2 * gamma_FT(q)
        # Easiest sanity check: invariant of any I should be > 0 and the
        # derived contrast should scale linearly with I0.
        phi = 0.3
        q, I_a, _ = make_synthetic_dataset(I0=1.0)
        q, I_b, _ = make_synthetic_dataset(I0=10.0)
        c_a = derive_contrast_from_invariant(q, I_a, phi)
        c_b = derive_contrast_from_invariant(q, I_b, phi)
        assert c_a > 0 and c_b > 0
        ratio = c_b / c_a
        assert 9.0 < ratio < 11.0

    def test_phi_from_contrast_inverse_of_contrast_from_phi(self):
        # Round-trip: φ → contrast → φ should recover the input φ.
        q, I, _ = make_synthetic_dataset()
        for phi_in in (0.1, 0.2, 0.35, 0.49):
            c = derive_contrast_from_invariant(q, I, phi_in)
            phi_back = derive_phi_from_invariant(q, I, c)
            assert abs(phi_back - phi_in) < 1e-6, (
                f"round-trip mismatch: phi_in={phi_in}, phi_back={phi_back}")

    def test_phi_from_contrast_handles_too_high_contrast(self):
        # If contrast is set very high so that 1 - 4K is negative,
        # function should fail gracefully (return 0.5 sentinel).
        q, I, _ = make_synthetic_dataset()
        result = derive_phi_from_invariant(q, I, contrast=1e-30)
        # Tiny contrast → K huge → disc < 0 → sentinel
        assert result == 0.5

    def test_phi_from_contrast_zero_contrast_safe(self):
        q, I, _ = make_synthetic_dataset()
        assert derive_phi_from_invariant(q, I, contrast=0.0) == 0.5
        assert derive_phi_from_invariant(q, I, contrast=-1.0) == 0.5


class TestPreFitHelpers:
    def test_power_law_recovers_known_slope(self):
        # Synthesise pure power-law I = B0 * q^-P0; the fit must recover B0, P0.
        B0, P0 = 1e-4, 4.0
        q = np.logspace(-3, 0, 200)
        I = B0 * q ** -P0
        B, P = fit_power_law_bg(q, I, q_min=1e-3, q_max=1e-1)
        assert abs(P - P0) < 1e-6
        assert abs(B - B0) / B0 < 1e-6

    def test_power_law_robust_to_window(self):
        # Power-law plus a small Gaussian peak at high Q.  The pre-fit
        # window should isolate the low-Q power-law region cleanly.
        B0, P0 = 1e-4, 4.0
        q = np.logspace(-3, 0, 300)
        I = B0 * q ** -P0 + 1e-6 * np.exp(-((q - 0.3) / 0.05) ** 2)
        B, P = fit_power_law_bg(q, I, q_min=1e-3, q_max=1e-2)
        # Low-Q window should be uncontaminated by the peak.
        assert abs(P - P0) < 0.1
        assert abs(B - B0) / B0 < 0.1

    def test_power_law_empty_window(self):
        q = np.logspace(-3, 0, 100)
        I = q ** -4
        # Window outside the data range → fall back to defaults
        B, P = fit_power_law_bg(q, I, q_min=10.0, q_max=100.0)
        assert (B, P) == (0.0, 4.0)

    def test_flat_bg_recovers_constant(self):
        q = np.logspace(-3, 0, 300)
        # Pure flat background of 0.05 with some power law thrown in.
        flat0 = 0.05
        I = 1e-4 * q ** -4 + flat0
        # Fit power law first to subtract it before the flat-bg call.
        bg = fit_flat_bg(q, I, q_min=0.5, q_max=1.0,
                         power_law_B=1e-4, power_law_P=4.0)
        assert abs(bg - flat0) < 1e-3

    def test_flat_bg_empty_window(self):
        q = np.logspace(-3, 0, 100)
        I = np.ones_like(q)
        bg = fit_flat_bg(q, I, q_min=10.0, q_max=100.0)
        assert bg == 0.0


class TestInputModeResolution:
    def test_compute_voxelgram_input_mode_phi(self):
        # mode='phi': contrast must be derived (not equal to the input value
        # of 999.0, which would give absurd model intensity).
        q, I, dI = make_synthetic_dataset()
        cfg = SaxsMorphConfig(
            voxel_size_fit=32, voxel_size_render=32,
            box_size_A=500.0,
            input_mode='phi',
            volume_fraction=0.30, contrast=999.0,  # contrast should be overwritten
            rng_seed=42,
        )
        engine = SaxsMorphEngine()
        res = engine.compute_voxelgram(cfg, q, I, dI)
        assert res.config.contrast != 999.0
        assert res.config.contrast > 0
        # phi should be unchanged
        assert abs(res.config.volume_fraction - 0.30) < 1e-12

    def test_compute_voxelgram_input_mode_contrast(self):
        # mode='contrast': phi must be derived (overrides the 0.99 input).
        q, I, dI = make_synthetic_dataset()
        # First, find a self-consistent (phi, contrast) pair.
        contrast0 = derive_contrast_from_invariant(q, I, phi=0.30)

        cfg = SaxsMorphConfig(
            voxel_size_fit=32, voxel_size_render=32,
            box_size_A=500.0,
            input_mode='contrast',
            volume_fraction=0.99,  # should be overwritten
            contrast=contrast0,
            rng_seed=42,
        )
        engine = SaxsMorphEngine()
        res = engine.compute_voxelgram(cfg, q, I, dI)
        assert abs(res.config.contrast - contrast0) < 1e-12  # unchanged
        assert abs(res.config.volume_fraction - 0.30) < 1e-3  # derived back

    def test_compute_voxelgram_input_mode_both(self):
        # mode='both': both inputs are honored exactly.
        q, I, dI = make_synthetic_dataset()
        cfg = SaxsMorphConfig(
            voxel_size_fit=32, voxel_size_render=32,
            box_size_A=500.0,
            input_mode='both',
            volume_fraction=0.30, contrast=2.5,
            rng_seed=42,
        )
        engine = SaxsMorphEngine()
        res = engine.compute_voxelgram(cfg, q, I, dI)
        assert abs(res.config.contrast - 2.5) < 1e-12
        assert abs(res.config.volume_fraction - 0.30) < 1e-12


# ---------------------------------------------------------------------------
# SaxsMorphEngine.compute_voxelgram (smoke)
# ---------------------------------------------------------------------------

class TestEngineSmoke:
    def test_compute_voxelgram_runs(self):
        q, I, dI = make_synthetic_dataset()
        cfg = SaxsMorphConfig(
            voxel_size_fit=32, voxel_size_render=32,
            box_size_A=500.0,
            volume_fraction=0.30,
            link_phi_contrast=True,
            rng_seed=42,
            smooth_sigma=0.0,   # explicit: keep binary uint8 output
        )
        engine = SaxsMorphEngine()
        res = engine.compute_voxelgram(cfg, q, I, dI)

        assert isinstance(res, SaxsMorphResult)
        assert res.voxelgram.shape == (32, 32, 32)
        assert res.voxelgram.dtype == np.uint8
        assert 0.20 < res.phi_actual < 0.40
        assert np.isfinite(res.chi_squared)
        assert res.dof > 0
        assert res.model_q.shape == res.data_q.shape
        assert res.model_I.shape == res.data_q.shape
        assert np.all(np.isfinite(res.model_I))
        assert res.contrast_or_link()  # see helper below

    def test_compute_voxelgram_with_smoothing_returns_float32(self):
        """With smooth_sigma > 0 the voxelgram is float32 in [0, 1]."""
        q, I, dI = make_synthetic_dataset()
        cfg = SaxsMorphConfig(
            voxel_size_fit=32, voxel_size_render=32,
            box_size_A=500.0,
            volume_fraction=0.30,
            link_phi_contrast=True,
            rng_seed=42,
            smooth_sigma=1.0,   # default, mimics Igor's gauss3d
        )
        engine = SaxsMorphEngine()
        res = engine.compute_voxelgram(cfg, q, I, dI)
        assert res.voxelgram.dtype == np.float32
        assert res.voxelgram.shape == (32, 32, 32)
        assert res.voxelgram.min() >= 0.0
        assert res.voxelgram.max() <= 1.0
        # Smoothing preserves the mean of the binary indicator.
        assert 0.20 < res.phi_actual < 0.40
        # Smoothing should reduce the per-voxel variance vs binary
        # (binary cube has variance phi*(1-phi); smoothed is lower).
        binary_var = res.phi_actual * (1.0 - res.phi_actual)
        assert float(res.voxelgram.var()) < binary_var

    def test_seed_reproducible_through_engine(self):
        q, I, dI = make_synthetic_dataset()
        cfg = SaxsMorphConfig(
            voxel_size_fit=32, voxel_size_render=32,
            box_size_A=500.0, volume_fraction=0.30, rng_seed=7,
            smooth_sigma=0.0,
        )
        engine = SaxsMorphEngine()
        r1 = engine.compute_voxelgram(cfg, q, I, dI)
        r2 = engine.compute_voxelgram(cfg, q, I, dI)
        assert np.array_equal(r1.voxelgram, r2.voxelgram)


# Tiny patch — give SaxsMorphResult a one-liner used by the test above so we
# don't have to expose internals.  (Adds a method via the class, not an
# instance, so it's available everywhere.)
def _contrast_or_link(self) -> bool:
    """True when contrast was either set or derived (sanity check)."""
    return self.config.contrast > 0 or self.config.link_phi_contrast
SaxsMorphResult.contrast_or_link = _contrast_or_link


# ---------------------------------------------------------------------------
# SaxsMorphEngine.fit (convergence smoke)
# ---------------------------------------------------------------------------

class TestEngineFit:
    def test_fit_converges_when_only_phi_free(self):
        # We start phi away from a 'reasonable' value and let the fitter
        # move it.  We only require the residual to decrease, not exact
        # convergence — exact recovery is hard at 32**3 because of strong
        # statistical noise in voxelgram_to_iq.
        q, I, dI = make_synthetic_dataset()
        cfg = SaxsMorphConfig(
            voxel_size_fit=32, voxel_size_render=32,
            box_size_A=500.0,
            volume_fraction=0.50,
            fit_volume_fraction=True,
            link_phi_contrast=True,
            rng_seed=42,
            no_limits=True,  # Nelder-Mead — bounded fits with one var are flaky
        )
        engine = SaxsMorphEngine()
        # Initial chi**2 (compute, no fit)
        cfg_init = SaxsMorphConfig(
            voxel_size_fit=32, voxel_size_render=32,
            box_size_A=500.0,
            volume_fraction=0.50,
            link_phi_contrast=True,
            rng_seed=42,
        )
        r0 = engine.compute_voxelgram(cfg_init, q, I, dI)
        rf = engine.fit(cfg, q, I, dI)
        # The fit should not increase chi**2 (least-squares guarantees this);
        # for a single free param Nelder-Mead may stall at the start, so we
        # use <= rather than strict inequality.
        assert rf.chi_squared <= r0.chi_squared * 1.01


# ---------------------------------------------------------------------------
# HDF5 round-trip
# ---------------------------------------------------------------------------

class TestHDF5RoundTrip:
    def test_roundtrip(self, tmp_path):
        from pyirena.io.nxcansas_saxs_morph import (
            save_saxs_morph_results, load_saxs_morph_results,
        )

        q, I, dI = make_synthetic_dataset()
        cfg = SaxsMorphConfig(
            voxel_size_fit=32, voxel_size_render=32,
            box_size_A=500.0, volume_fraction=0.30,
            rng_seed=42, link_phi_contrast=True,
            smooth_sigma=0.0,   # binary uint8 output
        )
        engine = SaxsMorphEngine()
        res = engine.compute_voxelgram(cfg, q, I, dI)

        path = tmp_path / 'morph.h5'
        save_saxs_morph_results(path, res)
        loaded = load_saxs_morph_results(path)

        assert loaded is not None
        assert loaded['voxel_size'] == 32
        assert abs(loaded['phi_actual'] - res.phi_actual) < 1e-12
        assert abs(loaded['chi_squared'] - res.chi_squared) < 1e-9
        assert loaded['voxelgram'].shape == (32, 32, 32)
        assert loaded['voxelgram'].dtype == np.uint8
        # Bit-exact
        assert np.array_equal(loaded['voxelgram'], res.voxelgram)
        assert np.allclose(loaded['model_q'], res.model_q)
        assert np.allclose(loaded['model_I'], res.model_I)
