"""
Utility functions for Unified Fit Model analysis.

Includes plotting, data export, and analysis helpers.
"""

import numpy as np
from typing import Optional, Dict, Tuple
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: matplotlib not available. Plotting functions disabled.")


def plot_fit_results(model, show_residuals: bool = True,
                     log_scale: bool = True,
                     save_path: Optional[str] = None):
    """
    Plot the fit results with data, model, and individual levels.

    Args:
        model: UnifiedFitModel instance with fit results
        show_residuals: Show residuals subplot
        log_scale: Use log-log scale
        save_path: Path to save the figure (if None, display only)
    """
    if not HAS_MATPLOTLIB:
        print("Matplotlib not available. Cannot plot.")
        return

    if model.q_data is None or model.fit_intensity is None:
        print("No fit data available. Run fit() first.")
        return

    # Create figure
    if show_residuals:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8),
                                        gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig, ax1 = plt.subplots(1, 1, figsize=(10, 6))

    # Plot data with error bars
    if model.error_data is not None:
        ax1.errorbar(model.q_data, model.I_data, yerr=model.error_data,
                    fmt='o', markersize=4, capsize=2, alpha=0.6,
                    label='Data', color='blue')
    else:
        ax1.plot(model.q_data, model.I_data, 'o', markersize=4,
                alpha=0.6, label='Data', color='blue')

    # Plot total fit
    ax1.plot(model.q_data, model.fit_intensity, '-', linewidth=2,
            label='Total Fit', color='red')

    # Plot individual levels
    colors = ['green', 'orange', 'purple', 'brown', 'pink']
    prev_Rg = 0.0
    for i in range(model.num_levels):
        level_intensity = model.calculate_level_intensity(model.q_data, i, prev_Rg)
        level_intensity += model.background / model.num_levels  # Distribute background
        ax1.plot(model.q_data, level_intensity, '--', linewidth=1.5,
                label=f'Level {i+1}', color=colors[i % len(colors)], alpha=0.7)
        prev_Rg = model.levels[i].Rg

    # Formatting
    if log_scale:
        ax1.set_xscale('log')
        ax1.set_yscale('log')

    ax1.set_xlabel('Q (Å⁻¹)', fontsize=12)
    ax1.set_ylabel('Intensity (cm⁻¹)', fontsize=12)
    ax1.set_title('Unified Fit Results', fontsize=14, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Plot residuals
    if show_residuals:
        residuals = model.I_data - model.fit_intensity
        if model.error_data is not None:
            normalized_residuals = residuals / model.error_data
            ylabel = 'Normalized Residuals'
        else:
            normalized_residuals = residuals
            ylabel = 'Residuals (cm⁻¹)'

        ax2.plot(model.q_data, normalized_residuals, 'o', markersize=3, color='black')
        ax2.axhline(0, color='red', linestyle='--', linewidth=1)
        ax2.set_xlabel('Q (Å⁻¹)', fontsize=12)
        ax2.set_ylabel(ylabel, fontsize=12)
        ax2.grid(True, alpha=0.3)

        if log_scale:
            ax2.set_xscale('log')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    plt.show()


def plot_guinier_analysis(model, level_idx: int = 0):
    """
    Plot Guinier analysis (ln(I) vs Q²) for a specific level.

    Args:
        model: UnifiedFitModel instance
        level_idx: Level index to analyze
    """
    if not HAS_MATPLOTLIB:
        print("Matplotlib not available. Cannot plot.")
        return

    if model.q_data is None:
        print("No data available.")
        return

    level = model.levels[level_idx]

    # Calculate Guinier region (typically Rg*Q < 1.3)
    q_max_guinier = 1.3 / level.Rg
    mask = model.q_data < q_max_guinier

    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot ln(I) vs Q²
    q_squared = model.q_data[mask] ** 2
    ln_I = np.log(model.I_data[mask])

    ax.plot(q_squared, ln_I, 'o', markersize=6, label='Data')

    # Plot model
    model_intensity = model.calculate_level_intensity(model.q_data[mask], level_idx)
    ln_I_model = np.log(model_intensity + model.background)
    ax.plot(q_squared, ln_I_model, '-', linewidth=2, color='red', label='Model')

    # Theoretical Guinier line
    guinier_theory = np.log(level.G) - q_squared * level.Rg ** 2 / 3.0
    ax.plot(q_squared, guinier_theory, '--', linewidth=1.5,
            color='green', label=f'Guinier: Rg={level.Rg:.2f} Å')

    ax.set_xlabel('Q² (Å⁻²)', fontsize=12)
    ax.set_ylabel('ln(I)', fontsize=12)
    ax.set_title(f'Guinier Analysis - Level {level_idx + 1}', fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


def plot_porod_analysis(model, level_idx: int = 0):
    """
    Plot Porod analysis (I*Q⁴ vs Q) for a specific level.

    Args:
        model: UnifiedFitModel instance
        level_idx: Level index to analyze
    """
    if not HAS_MATPLOTLIB:
        print("Matplotlib not available. Cannot plot.")
        return

    if model.q_data is None:
        print("No data available.")
        return

    level = model.levels[level_idx]

    fig, ax = plt.subplots(figsize=(8, 6))

    # Calculate I*Q⁴ (Porod representation)
    I_Q4_data = model.I_data * model.q_data ** 4
    I_Q4_model = model.fit_intensity * model.q_data ** 4

    ax.plot(model.q_data, I_Q4_data, 'o', markersize=4, alpha=0.6, label='Data')
    ax.plot(model.q_data, I_Q4_model, '-', linewidth=2, color='red', label='Model')

    # If P ≈ 4, show Porod plateau
    if abs(level.P - 4.0) < 0.1:
        porod_plateau = level.B * np.ones_like(model.q_data)
        ax.axhline(level.B, color='green', linestyle='--', linewidth=1.5,
                  label=f'Porod B={level.B:.4e}')

    ax.set_xlabel('Q (Å⁻¹)', fontsize=12)
    ax.set_ylabel('I·Q⁴ (cm⁻¹·Å⁴)', fontsize=12)
    ax.set_title(f'Porod Analysis - Level {level_idx + 1}', fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


def export_fit_results(model, filename: str, include_levels: bool = True):
    """
    Export fit results to a text file.

    Args:
        model: UnifiedFitModel instance with fit results
        filename: Output file path
        include_levels: Include individual level contributions
    """
    with open(filename, 'w') as f:
        # Header
        f.write("=" * 70 + "\n")
        f.write("UNIFIED FIT MODEL - FIT RESULTS\n")
        f.write("=" * 70 + "\n\n")

        # Parameters
        f.write(model.get_parameter_summary())
        f.write("\n\n")

        # Data table
        f.write("=" * 70 + "\n")
        f.write("DATA TABLE\n")
        f.write("=" * 70 + "\n")

        # Column headers
        headers = ["Q (1/Å)", "I_data (cm⁻¹)", "I_model (cm⁻¹)", "Residual"]
        if model.error_data is not None:
            headers.insert(2, "Error (cm⁻¹)")

        if include_levels:
            for i in range(model.num_levels):
                headers.append(f"Level_{i+1}")

        # Write headers
        f.write("\t".join(headers) + "\n")

        # Write data
        for i in range(len(model.q_data)):
            row = [
                f"{model.q_data[i]:.6e}",
                f"{model.I_data[i]:.6e}",
            ]

            if model.error_data is not None:
                row.append(f"{model.error_data[i]:.6e}")

            row.extend([
                f"{model.fit_intensity[i]:.6e}",
                f"{model.I_data[i] - model.fit_intensity[i]:.6e}"
            ])

            if include_levels:
                prev_Rg = 0.0
                for j in range(model.num_levels):
                    level_I = model.calculate_level_intensity(
                        np.array([model.q_data[i]]), j, prev_Rg
                    )[0]
                    row.append(f"{level_I:.6e}")
                    prev_Rg = model.levels[j].Rg

            f.write("\t".join(row) + "\n")

    print(f"Results exported to {filename}")


def calculate_size_distribution_moments(model, level_idx: int = 0) -> Dict[str, float]:
    """
    Calculate size distribution moments from Unified fit parameters.

    For spherical particles:
        R_avg = Rg * sqrt(5/3)
        Volume = 4/3 * π * R_avg³
        Surface = 4 * π * R_avg²

    Args:
        level_idx: Level index to analyze

    Returns:
        Dictionary with size parameters
    """
    level = model.levels[level_idx]

    # Average radius assuming spherical particles
    R_avg = level.Rg * np.sqrt(5.0 / 3.0)

    # Volume and surface
    volume = (4.0 / 3.0) * np.pi * R_avg ** 3
    surface = 4.0 * np.pi * R_avg ** 2

    # Number density (if G is available and assuming two-phase system)
    # G ≈ N * V² * Δρ² (for dilute systems)
    # This requires contrast information which we don't have, so it's qualitative

    return {
        'Rg': level.Rg,
        'R_avg_sphere': R_avg,
        'volume_sphere': volume,
        'surface_sphere': surface,
        'G': level.G,
        'B': level.B,
        'P': level.P
    }


def estimate_initial_parameters(q: np.ndarray, intensity: np.ndarray,
                               num_levels: int = 1) -> Dict:
    """
    Estimate initial parameters from data for fitting.

    Args:
        q: Scattering vector
        intensity: Measured intensity
        num_levels: Number of levels to estimate

    Returns:
        Dictionary with estimated parameters
    """
    params = {'levels': []}

    # Find low-Q and high-Q regions
    low_q_idx = np.argmin(q)
    high_q_idx = np.argmax(q)

    # Estimate background from high-Q plateau
    high_q_region = intensity[high_q_idx - 5:high_q_idx + 1]
    background = np.median(high_q_region)

    params['background'] = max(background, 0.0)

    # For single level
    if num_levels == 1:
        # Estimate G from low-Q intensity
        G = intensity[low_q_idx]

        # Estimate Rg from Guinier region
        # Try to find where ln(I) vs Q² is linear
        mask = (q < 0.1) & (intensity > background * 2)
        if np.sum(mask) > 3:
            q_guinier = q[mask]
            I_guinier = intensity[mask]

            # Fit ln(I) = ln(G) - q²Rg²/3
            ln_I = np.log(I_guinier)
            q_squared = q_guinier ** 2

            # Linear regression
            coeffs = np.polyfit(q_squared, ln_I, 1)
            Rg = np.sqrt(-3.0 * coeffs[0]) if coeffs[0] < 0 else 10.0
        else:
            Rg = 10.0

        # Estimate B from high-Q power law
        high_q_mask = q > 0.1
        if np.sum(high_q_mask) > 3:
            q_high = q[high_q_mask]
            I_high = intensity[high_q_mask] - background

            # Fit log(I) = log(B) - P*log(Q)
            log_I = np.log(I_high + 1e-10)
            log_q = np.log(q_high)

            coeffs = np.polyfit(log_q, log_I, 1)
            P = -coeffs[0]
            B = np.exp(coeffs[1])
        else:
            P = 4.0
            B = 1e-4

        params['levels'].append({
            'Rg': Rg,
            'G': G,
            'P': P,
            'B': B,
            'RgCO': 0.0,
            'ETA': Rg,
            'PACK': 0.0
        })

    return params


def apply_parameters_from_dict(model, params_dict: Dict):
    """
    Apply parameters from a dictionary to a model.

    Args:
        model: UnifiedFitModel instance
        params_dict: Dictionary with parameters (from estimate_initial_parameters)
    """
    if 'background' in params_dict:
        model.background = params_dict['background']

    if 'levels' in params_dict:
        for i, level_params in enumerate(params_dict['levels']):
            if i < len(model.levels):
                for key, value in level_params.items():
                    if hasattr(model.levels[i], key):
                        setattr(model.levels[i], key, value)
