"""
Unified Fit Model — minimal example usage.

Fits a one-level Unified model to synthetic Guinier + power-law data and
compares fitted with true parameters. Run with::

    python -m pyirena.examples.unified_fit_demo

Moved here from the ``__main__`` block of ``pyirena/core/unified.py``.
"""

import numpy as np

from pyirena.core.unified import UnifiedFitModel


def main():
    # Create a simple example with synthetic data
    print("Unified Fit Model - Example Usage")
    print("=" * 70)

    # Generate synthetic data (simple power law + Guinier)
    q = np.logspace(-3, 0, 100)  # 0.001 to 1.0 Å^-1

    # True parameters for level 1
    true_Rg = 50.0
    true_G = 1000.0
    true_P = 4.0
    true_B = 1e-3
    true_background = 0.01

    # Generate "measured" intensity
    model = UnifiedFitModel(num_levels=1)
    model.levels[0].Rg = true_Rg
    model.levels[0].G = true_G
    model.levels[0].P = true_P
    model.levels[0].B = true_B
    model.background = true_background

    I_true = model.calculate_intensity(q)

    # Add noise
    np.random.seed(42)
    noise = 0.05 * I_true * np.random.randn(len(q))
    I_measured = I_true + noise
    I_error = 0.05 * I_true

    # Fit with initial guess
    fit_model = UnifiedFitModel(num_levels=1)
    fit_model.levels[0].Rg = 40.0  # Initial guess
    fit_model.levels[0].G = 800.0
    fit_model.levels[0].P = 4.0
    fit_model.levels[0].B = 1e-4
    fit_model.levels[0].fit_P = False  # Fix P for this example
    fit_model.background = 0.0

    print("\nPerforming fit...")
    fit_model.fit(q, I_measured, I_error, verbose=1)

    print("\n" + fit_model.get_parameter_summary())

    print("\nTrue vs Fitted Parameters:")
    print(f"  Rg:  True = {true_Rg:.2f}, Fitted = {fit_model.levels[0].Rg:.2f}")
    print(f"  G:   True = {true_G:.2f}, Fitted = {fit_model.levels[0].G:.2f}")
    print(f"  B:   True = {true_B:.4e}, Fitted = {fit_model.levels[0].B:.4e}")
    print(f"  Bkg: True = {true_background:.4f}, Fitted = {fit_model.background:.4f}")

    # Calculate invariant
    inv_results = fit_model.calculate_invariant(0)
    print(f"\nInvariant Q = {inv_results['invariant']:.4e} cm⁻¹ Å⁻³")
    if inv_results['surface_to_volume'] is not None:
        print(f"Surface/Volume = {inv_results['surface_to_volume']:.4e} m²/cm³")


if __name__ == "__main__":
    main()
