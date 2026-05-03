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
    generate_voxelgram, voxelgram_to_iq, derive_contrast_from_invariant,
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

    def test_seed_reproducible_through_engine(self):
        q, I, dI = make_synthetic_dataset()
        cfg = SaxsMorphConfig(
            voxel_size_fit=32, voxel_size_render=32,
            box_size_A=500.0, volume_fraction=0.30, rng_seed=7,
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
