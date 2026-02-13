"""
Unit tests for the Unified Fit model.
"""

import numpy as np
import pytest
from pyirena.core.unified import UnifiedFitModel, UnifiedLevel


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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
