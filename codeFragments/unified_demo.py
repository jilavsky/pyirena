"""
Demonstration script for Unified Fit Model.

Shows how to:
1. Load data from NXcanSAS file
2. Set up the model with multiple levels
3. Perform fitting
4. Visualize and export results
"""

import numpy as np
from unified import UnifiedFitModel, load_data_from_nxcansas
from unified_utils import (
    plot_fit_results, plot_guinier_analysis, plot_porod_analysis,
    export_fit_results, estimate_initial_parameters, apply_parameters_from_dict,
    calculate_size_distribution_moments
)


def demo_synthetic_data():
    """
    Demonstrate fitting with synthetic data.
    """
    print("\n" + "=" * 70)
    print("DEMO 1: Fitting Synthetic Data (Single Level)")
    print("=" * 70 + "\n")

    # Generate synthetic data
    q = np.logspace(-3, 0, 150)  # 0.001 to 1.0 Å^-1

    # True parameters
    true_Rg = 50.0
    true_G = 1000.0
    true_P = 4.0
    true_B = 2e-4
    true_background = 0.05

    # Create model with true parameters
    model_true = UnifiedFitModel(num_levels=1)
    model_true.levels[0].Rg = true_Rg
    model_true.levels[0].G = true_G
    model_true.levels[0].P = true_P
    model_true.levels[0].B = true_B
    model_true.background = true_background

    # Generate "measured" data
    I_true = model_true.calculate_intensity(q)

    # Add realistic noise (5% relative error)
    np.random.seed(42)
    I_error = 0.05 * I_true
    noise = I_error * np.random.randn(len(q))
    I_measured = I_true + noise

    # Estimate initial parameters
    print("Estimating initial parameters from data...")
    initial_params = estimate_initial_parameters(q, I_measured, num_levels=1)

    # Create fitting model
    fit_model = UnifiedFitModel(num_levels=1)
    apply_parameters_from_dict(fit_model, initial_params)

    # Configure fitting
    fit_model.levels[0].fit_Rg = True
    fit_model.levels[0].fit_G = True
    fit_model.levels[0].fit_P = True
    fit_model.levels[0].fit_B = True
    fit_model.fit_background = True

    # Set reasonable bounds
    fit_model.levels[0].Rg_limits = (1.0, 500.0)
    fit_model.levels[0].G_limits = (1.0, 1e5)
    fit_model.levels[0].P_limits = (2.0, 6.0)
    fit_model.levels[0].B_limits = (1e-10, 1.0)

    print("\nInitial parameters:")
    print(f"  Rg  = {fit_model.levels[0].Rg:.2f} Å")
    print(f"  G   = {fit_model.levels[0].G:.2e} cm⁻¹")
    print(f"  P   = {fit_model.levels[0].P:.2f}")
    print(f"  B   = {fit_model.levels[0].B:.2e}")
    print(f"  Bkg = {fit_model.background:.4f}")

    # Perform fit
    print("\nPerforming fit...")
    results = fit_model.fit(q, I_measured, I_error, verbose=0)

    # Display results
    print("\n" + fit_model.get_parameter_summary())

    print("\nComparison: True vs Fitted Parameters")
    print("-" * 50)
    print(f"{'Parameter':<10} {'True':<15} {'Fitted':<15} {'Error %':<10}")
    print("-" * 50)

    def print_comparison(name, true_val, fitted_val):
        error = 100 * abs(fitted_val - true_val) / true_val if true_val != 0 else 0
        print(f"{name:<10} {true_val:<15.4e} {fitted_val:<15.4e} {error:<10.2f}")

    print_comparison("Rg", true_Rg, fit_model.levels[0].Rg)
    print_comparison("G", true_G, fit_model.levels[0].G)
    print_comparison("P", true_P, fit_model.levels[0].P)
    print_comparison("B", true_B, fit_model.levels[0].B)
    print_comparison("Bkg", true_background, fit_model.background)
    print("-" * 50)

    # Calculate derived quantities
    print("\nDerived Quantities:")
    print("-" * 50)

    # Size distribution moments
    size_params = calculate_size_distribution_moments(fit_model, 0)
    print(f"Average sphere radius: {size_params['R_avg_sphere']:.2f} Å")
    print(f"Volume (sphere):       {size_params['volume_sphere']:.2e} Å³")
    print(f"Surface (sphere):      {size_params['surface_sphere']:.2e} Å²")

    # Invariant
    inv_results = fit_model.calculate_invariant(0)
    print(f"\nInvariant Q:           {inv_results['invariant']:.4e} cm⁻¹ Å⁻³")
    if inv_results['surface_to_volume'] is not None:
        print(f"Surface/Volume ratio:  {inv_results['surface_to_volume']:.4e} m²/cm³")

    # Plot results
    try:
        plot_fit_results(fit_model, show_residuals=True, log_scale=True)
        plot_guinier_analysis(fit_model, 0)
        plot_porod_analysis(fit_model, 0)
    except Exception as e:
        print(f"\nCould not create plots: {e}")

    # Export results
    export_fit_results(fit_model, "unified_fit_results.txt", include_levels=True)

    return fit_model


def demo_multi_level():
    """
    Demonstrate fitting with multiple structural levels.
    """
    print("\n" + "=" * 70)
    print("DEMO 2: Fitting Synthetic Data (Two Levels)")
    print("=" * 70 + "\n")

    # Generate synthetic data with two structural levels
    q = np.logspace(-3, 0, 200)

    # Level 1: Primary particles (small)
    model_true = UnifiedFitModel(num_levels=2)
    model_true.levels[0].Rg = 20.0
    model_true.levels[0].G = 100.0
    model_true.levels[0].P = 4.0
    model_true.levels[0].B = 1e-3

    # Level 2: Aggregates (large)
    model_true.levels[1].Rg = 200.0
    model_true.levels[1].G = 10000.0
    model_true.levels[1].P = 2.5  # Mass fractal
    model_true.levels[1].B = 1e-6
    model_true.levels[1].RgCO = 20.0  # Linked to level 1 Rg
    model_true.levels[1].link_RGCO = True

    model_true.background = 0.01

    # Generate data
    I_true = model_true.calculate_intensity(q)

    # Add noise
    np.random.seed(123)
    I_error = 0.05 * I_true
    noise = I_error * np.random.randn(len(q))
    I_measured = I_true + noise

    # Create fitting model
    fit_model = UnifiedFitModel(num_levels=2)

    # Initialize Level 1 (primary particles)
    fit_model.levels[0].Rg = 15.0
    fit_model.levels[0].G = 80.0
    fit_model.levels[0].P = 4.0
    fit_model.levels[0].B = 1e-4
    fit_model.levels[0].fit_Rg = True
    fit_model.levels[0].fit_G = True
    fit_model.levels[0].fit_P = False  # Fix Porod law
    fit_model.levels[0].fit_B = True

    # Initialize Level 2 (aggregates)
    fit_model.levels[1].Rg = 150.0
    fit_model.levels[1].G = 8000.0
    fit_model.levels[1].P = 2.5
    fit_model.levels[1].B = 1e-7
    fit_model.levels[1].link_RGCO = True  # Link to level 1
    fit_model.levels[1].fit_Rg = True
    fit_model.levels[1].fit_G = True
    fit_model.levels[1].fit_P = True
    fit_model.levels[1].fit_B = True

    fit_model.background = 0.0
    fit_model.fit_background = True

    print("Initial parameters:")
    print(fit_model.get_parameter_summary())

    # Perform fit
    print("\nPerforming fit...")
    results = fit_model.fit(q, I_measured, I_error, verbose=0)

    # Display results
    print("\n" + fit_model.get_parameter_summary())

    print("\nComparison: True vs Fitted Parameters")
    print("-" * 70)
    print("Level 1 (Primary Particles):")
    print(f"  Rg: True = {model_true.levels[0].Rg:.2f}, Fitted = {fit_model.levels[0].Rg:.2f}")
    print(f"  G:  True = {model_true.levels[0].G:.2f}, Fitted = {fit_model.levels[0].G:.2f}")
    print(f"  B:  True = {model_true.levels[0].B:.2e}, Fitted = {fit_model.levels[0].B:.2e}")

    print("\nLevel 2 (Aggregates):")
    print(f"  Rg: True = {model_true.levels[1].Rg:.2f}, Fitted = {fit_model.levels[1].Rg:.2f}")
    print(f"  G:  True = {model_true.levels[1].G:.2f}, Fitted = {fit_model.levels[1].G:.2f}")
    print(f"  P:  True = {model_true.levels[1].P:.2f}, Fitted = {fit_model.levels[1].P:.2f}")
    print(f"  B:  True = {model_true.levels[1].B:.2e}, Fitted = {fit_model.levels[1].B:.2e}")

    # Plot results
    try:
        plot_fit_results(fit_model, show_residuals=True, log_scale=True)
    except Exception as e:
        print(f"\nCould not create plots: {e}")

    return fit_model


def demo_with_correlations():
    """
    Demonstrate fitting with correlation function.
    """
    print("\n" + "=" * 70)
    print("DEMO 3: Fitting with Correlation Function")
    print("=" * 70 + "\n")

    q = np.logspace(-3, 0, 150)

    # Create model with correlations
    model_true = UnifiedFitModel(num_levels=1)
    model_true.levels[0].Rg = 30.0
    model_true.levels[0].G = 500.0
    model_true.levels[0].P = 4.0
    model_true.levels[0].B = 5e-4
    model_true.levels[0].correlations = True
    model_true.levels[0].ETA = 60.0  # Correlation distance
    model_true.levels[0].PACK = 0.3  # Packing factor
    model_true.background = 0.02

    # Generate data
    I_true = model_true.calculate_intensity(q)

    # Add noise
    np.random.seed(456)
    I_error = 0.05 * I_true
    I_measured = I_true + I_error * np.random.randn(len(q))

    # Fit without correlations first
    fit_model_no_corr = UnifiedFitModel(num_levels=1)
    fit_model_no_corr.levels[0].Rg = 25.0
    fit_model_no_corr.levels[0].G = 400.0
    fit_model_no_corr.levels[0].P = 4.0
    fit_model_no_corr.levels[0].B = 1e-4
    fit_model_no_corr.levels[0].fit_P = False

    print("Fitting WITHOUT correlations...")
    results_no_corr = fit_model_no_corr.fit(q, I_measured, I_error, verbose=0)

    # Fit with correlations
    fit_model_corr = UnifiedFitModel(num_levels=1)
    fit_model_corr.levels[0].Rg = 25.0
    fit_model_corr.levels[0].G = 400.0
    fit_model_corr.levels[0].P = 4.0
    fit_model_corr.levels[0].B = 1e-4
    fit_model_corr.levels[0].correlations = True
    fit_model_corr.levels[0].ETA = 50.0
    fit_model_corr.levels[0].PACK = 0.2
    fit_model_corr.levels[0].fit_P = False
    fit_model_corr.levels[0].fit_ETA = True
    fit_model_corr.levels[0].fit_PACK = True

    print("Fitting WITH correlations...")
    results_corr = fit_model_corr.fit(q, I_measured, I_error, verbose=0)

    # Compare results
    print("\n" + "=" * 70)
    print("COMPARISON OF FIT QUALITY")
    print("=" * 70)
    print(f"\nWithout correlations:")
    print(f"  χ² = {results_no_corr['chi_squared']:.4e}")
    print(f"  Reduced χ² = {results_no_corr['reduced_chi_squared']:.4f}")

    print(f"\nWith correlations:")
    print(f"  χ² = {results_corr['chi_squared']:.4e}")
    print(f"  Reduced χ² = {results_corr['reduced_chi_squared']:.4f}")

    print(f"\nTrue correlation parameters:")
    print(f"  ETA  = {model_true.levels[0].ETA:.2f} Å")
    print(f"  PACK = {model_true.levels[0].PACK:.3f}")

    print(f"\nFitted correlation parameters:")
    print(f"  ETA  = {fit_model_corr.levels[0].ETA:.2f} Å")
    print(f"  PACK = {fit_model_corr.levels[0].PACK:.3f}")

    return fit_model_corr


def demo_load_real_data(file_path: str):
    """
    Demonstrate loading and fitting real data from NXcanSAS file.

    Args:
        file_path: Path to NXcanSAS HDF5 file
    """
    print("\n" + "=" * 70)
    print("DEMO 4: Loading Real Data from NXcanSAS File")
    print("=" * 70 + "\n")

    try:
        # Load data
        print(f"Loading data from: {file_path}")
        data = load_data_from_nxcansas(file_path)

        print(f"Data loaded successfully:")
        print(f"  Number of points: {len(data['Q'])}")
        print(f"  Q range: {data['Q'][0]:.4e} to {data['Q'][-1]:.4e} Å⁻¹")
        print(f"  I range: {np.min(data['Intensity']):.4e} to {np.max(data['Intensity']):.4e} cm⁻¹")

        # Create model
        fit_model = UnifiedFitModel(num_levels=1)

        # Estimate initial parameters
        initial_params = estimate_initial_parameters(
            data['Q'], data['Intensity'], num_levels=1
        )
        apply_parameters_from_dict(fit_model, initial_params)

        print("\nEstimated initial parameters:")
        print(fit_model.get_parameter_summary())

        # Perform fit
        print("\nPerforming fit...")
        results = fit_model.fit(
            data['Q'], data['Intensity'], data['Error'],
            verbose=1
        )

        print("\n" + fit_model.get_parameter_summary())

        # Plot results
        try:
            plot_fit_results(fit_model, show_residuals=True, log_scale=True,
                           save_path="real_data_fit.png")
        except Exception as e:
            print(f"Could not create plots: {e}")

        # Export results
        export_fit_results(fit_model, "real_data_fit_results.txt", include_levels=True)

        return fit_model

    except Exception as e:
        print(f"Error loading or fitting data: {e}")
        return None


if __name__ == "__main__":
    print("\n")
    print("╔" + "═" * 68 + "╗")
    print("║" + " " * 15 + "UNIFIED FIT MODEL - DEMONSTRATION" + " " * 20 + "║")
    print("╚" + "═" * 68 + "╝")

    # Run demos
    model1 = demo_synthetic_data()
    model2 = demo_multi_level()
    model3 = demo_with_correlations()

    # Uncomment to test with real data
    # model4 = demo_load_real_data("/path/to/your/nxcansas/file.h5")

    print("\n" + "=" * 70)
    print("All demonstrations completed!")
    print("=" * 70 + "\n")
