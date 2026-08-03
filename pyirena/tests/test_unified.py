"""
Unit tests for the Unified Fit model.
"""

import numpy as np
import pytest
from pyirena.core.unified import UnifiedFitModel, UnifiedLevel, compute_invariant_sv


class TestUnifiedLevel:
    """Test UnifiedLevel dataclass."""

    def test_default_initialization(self):
        """Test default parameter values."""
        level = UnifiedLevel()
        assert level.Rg == 10.0
        assert level.G == 1.0
        assert level.P == 4.0
        assert level.B == 1e-4
        assert level.K == 1.0
        assert level.correlations is False
        assert level.mass_fractal is False

    def test_auto_calculate_K(self):
        """Test automatic K calculation."""
        level = UnifiedLevel(P=4.5)
        level.auto_calculate_K()
        assert level.K == 1.0

        level.P = 2.5
        level.auto_calculate_K()
        assert level.K == 1.06

    def test_auto_calculate_B(self):
        """Test automatic B calculation from Hammouda relationship."""
        level = UnifiedLevel(Rg=50.0, G=1000.0, P=4.0, link_B=True)
        level.auto_calculate_B()
        assert level.B > 0
        assert np.isfinite(level.B)


class TestUnifiedFitModel:
    """Test UnifiedFitModel class."""

    def test_initialization(self):
        """Test model initialization."""
        model = UnifiedFitModel(num_levels=1)
        assert model.num_levels == 1
        assert len(model.levels) == 1
        assert model.background == 0.0

    def test_invalid_num_levels(self):
        """Test that invalid number of levels raises error."""
        with pytest.raises(ValueError):
            UnifiedFitModel(num_levels=0)
        with pytest.raises(ValueError):
            UnifiedFitModel(num_levels=6)

    def test_calculate_intensity_single_level(self):
        """Test intensity calculation for single level."""
        model = UnifiedFitModel(num_levels=1)
        model.levels[0].Rg = 50.0
        model.levels[0].G = 1000.0
        model.levels[0].P = 4.0
        model.levels[0].B = 1e-3
        model.background = 0.01

        q = np.logspace(-3, 0, 50)
        intensity = model.calculate_intensity(q)

        assert len(intensity) == len(q)
        assert np.all(intensity > 0)
        assert np.all(np.isfinite(intensity))

    def test_sphere_amplitude(self):
        """Test sphere amplitude function for correlations."""
        model = UnifiedFitModel(num_levels=1)
        q = np.array([0.01, 0.1, 1.0])
        eta = 50.0

        amp = model.sphere_amplitude(q, eta)
        assert len(amp) == len(q)
        assert np.all(np.isfinite(amp))
        assert np.all(np.abs(amp) <= 1.0)

    def test_fitting_synthetic_data(self):
        """Test fitting with synthetic data."""
        # Generate synthetic data
        q = np.logspace(-3, 0, 100)

        true_model = UnifiedFitModel(num_levels=1)
        true_model.levels[0].Rg = 50.0
        true_model.levels[0].G = 1000.0
        true_model.levels[0].P = 4.0
        true_model.levels[0].B = 1e-3
        true_model.background = 0.01

        I_true = true_model.calculate_intensity(q)

        # Add small noise
        np.random.seed(42)
        noise = 0.02 * I_true * np.random.randn(len(q))
        I_measured = I_true + noise
        I_error = 0.05 * I_true

        # Fit with initial guess
        fit_model = UnifiedFitModel(num_levels=1)
        fit_model.levels[0].Rg = 40.0
        fit_model.levels[0].G = 800.0
        fit_model.levels[0].P = 4.0
        fit_model.levels[0].B = 5e-4
        fit_model.levels[0].fit_P = False  # Fix P
        fit_model.background = 0.0

        results = fit_model.fit(q, I_measured, I_error, verbose=0)

        assert results['success']
        assert results['chi_squared'] > 0
        assert results['reduced_chi_squared'] > 0
        # Check that fitted parameters are close to true values
        assert abs(fit_model.levels[0].Rg - 50.0) / 50.0 < 0.2  # Within 20%
        assert abs(fit_model.levels[0].G - 1000.0) / 1000.0 < 0.2

    def test_calculate_invariant(self):
        """Test invariant calculation."""
        model = UnifiedFitModel(num_levels=1)
        model.levels[0].Rg = 50.0
        model.levels[0].G = 1000.0
        model.levels[0].P = 4.0
        model.levels[0].B = 1e-3

        inv_results = model.calculate_invariant(0, q_min=0.0, q_max=1.0)

        assert 'invariant' in inv_results
        assert 'invariant_cm4' in inv_results
        assert inv_results['invariant'] > 0
        assert np.isfinite(inv_results['invariant'])

    def test_parameter_packing_unpacking(self):
        """Test parameter packing and unpacking for fitting."""
        model = UnifiedFitModel(num_levels=2)
        model.levels[0].Rg = 20.0
        model.levels[0].G = 100.0
        model.levels[1].Rg = 200.0
        model.levels[1].G = 1000.0
        model.background = 0.05

        # Pack parameters
        packed = model._pack_parameters()
        assert len(packed) > 0

        # Modify and unpack
        packed_modified = packed * 1.1
        model._unpack_parameters(packed_modified)

        # Check that values changed
        assert model.levels[0].Rg != 20.0
        assert abs(model.levels[0].Rg - 22.0) < 0.1

    def test_get_parameter_summary(self):
        """Test parameter summary generation."""
        model = UnifiedFitModel(num_levels=1)
        summary = model.get_parameter_summary()

        assert isinstance(summary, str)
        assert "Level 1" in summary
        assert "Rg" in summary
        assert "=" in summary

    def test_eta_pack_not_fitted_without_correlations(self):
        """ETA/PACK must not be free parameters when correlations are off.

        With correlations disabled, ETA and PACK have no effect on the model,
        so fitting them makes the least-squares problem rank-deficient and the
        solver thrashes (badly amplified by slit smearing).  Regression guard:
        the packed parameter vector must exclude ETA/PACK unless correlations
        are on, regardless of the fit_ETA / fit_PACK flags.
        """
        model = UnifiedFitModel(num_levels=1)
        model.fit_background = False
        lvl = model.levels[0]
        lvl.fit_Rg = lvl.fit_G = lvl.fit_P = lvl.fit_B = False
        lvl.fit_ETA = lvl.fit_PACK = True

        # correlations off → ETA/PACK excluded
        lvl.correlations = False
        assert lvl.fit_ETA_effective is False
        assert lvl.fit_PACK_effective is False
        assert len(model._pack_parameters()) == 0
        lo, hi = model._get_bounds()
        assert len(lo) == 0 and len(hi) == 0

        # correlations on → ETA/PACK included
        lvl.correlations = True
        assert lvl.fit_ETA_effective is True
        assert lvl.fit_PACK_effective is True
        assert len(model._pack_parameters()) == 2


class TestComputeInvariantSv:
    """Regression tests for compute_invariant_sv().

    The expected values were pinned from the pre-refactor implementations
    (LevelParametersWidget.update_porod_surface_and_invariant in the GUI and
    batch._compute_invariant_sv), which were verified to agree with each
    other and with Igor Irena's IR1A_UpdatePorodSfcandInvariant.  If these
    numbers change, GUI-displayed and batch-saved Invariant/Sv change too.
    """

    def test_porod_level_pinned(self):
        inv, sv = compute_invariant_sv(100.0, 50.0, 1e-4, 4.0, 0.0, 0.0, 0.0, False)
        assert inv == pytest.approx(3.417756004764339e+21, rel=1e-10)
        assert sv == pytest.approx(919.1974644212241, rel=1e-10)

    def test_porod_level_with_cutoff_pinned(self):
        inv, sv = compute_invariant_sv(1000.0, 100.0, 2e-5, 4.0, 20.0, 0.0, 0.0, False)
        assert inv == pytest.approx(2.611245870080375e+21, rel=1e-10)
        assert sv == pytest.approx(240.62021042033047, rel=1e-10)

    def test_correlated_level_pinned(self):
        inv, sv = compute_invariant_sv(100.0, 50.0, 1e-4, 4.0, 0.0, 120.0, 2.0, True)
        assert inv == pytest.approx(3.2115542177814733e+21, rel=1e-10)
        assert sv == pytest.approx(978.215667727382, rel=1e-10)

    def test_sv_only_in_porod_regime(self):
        """Sv is None outside P = 3.95 - 4.05, invariant may still be valid."""
        inv, sv = compute_invariant_sv(100.0, 50.0, 1e-4, 3.8, 0.0, 0.0, 0.0, False)
        assert inv is not None and inv > 0
        assert sv is None

    def test_invalid_rg_returns_none(self):
        assert compute_invariant_sv(100.0, 0.0, 1e-4, 4.0) == (None, None)
        assert compute_invariant_sv(100.0, -5.0, 1e-4, 4.0) == (None, None)

    def test_negative_invariant_returns_none(self):
        """Bad extrapolation (large negative Porod tail) yields (None, None)."""
        assert compute_invariant_sv(500.0, 80.0, 0.05, 2.5) == (None, None)

    def test_guinier_only_level(self):
        """B=0 still yields a finite (Guinier-only) invariant, like the GUI did."""
        inv, sv = compute_invariant_sv(100.0, 50.0, 0.0, 4.0)
        assert inv is not None and inv > 0
        assert sv == 0.0

    def test_batch_wrapper_delegates(self):
        """batch._compute_invariant_sv must return identical values."""
        from pyirena.batch import _compute_invariant_sv
        args = (100.0, 50.0, 1e-4, 4.0, 0.0, 0.0, 0.0, False)
        assert _compute_invariant_sv(*args) == compute_invariant_sv(*args)

    def test_calculate_invariant_sv_units_match(self):
        """calculate_invariant()'s Sv now carries the 1e4 A^-1 -> m^2/cm^3
        conversion; over the same integration range it must equal
        1e4*pi*B/invariant."""
        model = UnifiedFitModel(num_levels=1)
        model.levels[0].Rg = 50.0
        model.levels[0].G = 100.0
        model.levels[0].P = 4.0
        model.levels[0].B = 1e-4
        res = model.calculate_invariant(0, q_min=0.0, q_max=1.0)
        assert res['surface_to_volume'] == pytest.approx(
            1e4 * np.pi * model.levels[0].B / res['invariant'], rel=1e-12
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
