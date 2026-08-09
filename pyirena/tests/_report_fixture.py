"""
Shared fixture inputs for the report formatter tests.

Kept in its own module (not a conftest fixture) so the same inputs can be fed
to the pre-move and post-move implementations when checking that relocating the
formatter changed nothing.
"""

import numpy as np


def data_info() -> dict:
    q = np.logspace(-3, -1, 50)
    intensity = 1e3 * q ** -3
    return {"Q": q, "I": intensity, "I_error": 0.05 * intensity}


def unified_results() -> dict:
    q = np.logspace(-3, -1, 50)
    intensity = 1e3 * q ** -3
    return {
        "Q": q,
        "intensity_data": intensity,
        "intensity_model": 1.01 * intensity,
        "residuals": np.linspace(-1.5, 1.5, 50),
        "intensity_error": 0.05 * intensity,
        "num_levels": 2,
        "background": 0.012,
        "chi_squared": 1.4321,
        "timestamp": "2026-08-09T10:00:00",
        "slit_length": 0.0,
        "fit_quality": {
            "robust_scale_s": 1.8,
            "realistic_reduced_chi2_floor": 3.24,
            "max_abs_frac_misfit": 0.087,
            "longest_same_sign_run": 7,
        },
        "levels": [
            {
                "G": 12345.0, "G_err": 210.0,
                "Rg": 185.4, "Rg_err": 3.1,
                "B": 1.2e-5, "B_err": 4e-7,
                "P": 4.01, "P_err": 0.05,
                "RgCutoff": 0.0,
                "correlated": False,
                "Sv": 1.7e5,
                "Invariant": 3.3e-4,
            },
            {
                "G": 55.0,
                "Rg": 21.7,
                "B": 3.4e-8,
                "P": 3.2,
                "RgCutoff": 190.0,
                "correlated": True,
                "ETA": 42.5,
                "PACK": 0.7412,
            },
        ],
    }


def sizes_results() -> dict:
    r_grid = np.linspace(10.0, 400.0, 40)
    dist = np.exp(-((r_grid - 120.0) / 30.0) ** 2)
    return {
        "r_grid": r_grid,
        "distribution": dist,
        "residuals": np.linspace(-0.8, 0.8, 50),
        "chi_squared": 2.13,
        "volume_fraction": 0.0123,
        "rg": 143.2,
        "n_iterations": 25,
        "method": "MaxEnt",
        "shape": "Spheroid",
        "contrast": 34.5,
        "aspect_ratio": 2.0,
        "n_bins": 40,
        "r_min": 10.0,
        "r_max": 400.0,
        "log_spacing": True,
        "background": 0.0021,
        "error_scale": 1.0,
        "power_law_B": 1.5e-6,
        "power_law_P": 3.85,
        "cursor_q_min": 0.0012,
        "cursor_q_max": 0.092,
        "maxent_sky_background": 1e-6,
        "maxent_stability": 0.01,
        "maxent_max_iter": 100,
        "timestamp": "2026-08-09T10:05:00",
        "slit_length": 0.0,
        "fit_quality": {"robust_scale_s": 1.2, "longest_same_sign_run": 4},
    }


def simple_fit_results() -> dict:
    return {
        "model": "Guinier",
        "chi_squared": 55.2,
        "reduced_chi_squared": 1.12,
        "dof": 47,
        "q_min": 0.0015,
        "q_max": 0.05,
        "use_complex_bg": False,
        "params": {"G": 1234.5, "Rg": 154.3},
        "params_std": {"G": 21.0, "Rg": 2.4},
        "derived": {"I0": 1234.5, "Rg_nm": 15.43},
        "timestamp": "2026-08-09T10:10:00",
        "slit_length": 0.0,
        "fit_quality": {"robust_scale_s": 0.9},
    }


def waxs_results() -> dict:
    return {
        "n_peaks": 2,
        "bg_shape": "SNIP",
        "chi_squared": 31.4,
        "reduced_chi_squared": 1.05,
        "dof": 30,
        "q_min": 1.0,
        "q_max": 5.0,
        "timestamp": "2026-08-09T10:15:00",
        "peaks": [
            # Peak rows carry {'value','fit','lo','hi'} dicts in the GUI…
            {"shape": "Gauss",
             "A": {"value": 500.0, "fit": True},
             "Q0": {"value": 2.5, "fit": True},
             "FWHM": {"value": 0.1, "fit": True}},
            # …and plain floats when loaded back from HDF5.
            {"shape": "Pseudo-Voigt", "A": 120.0, "Q0": 3.9,
             "FWHM": 0.22, "eta": 0.4},
        ],
        "peaks_std": [{"A": 12.0, "Q0": 0.002, "FWHM": 0.004}, {}],
        "fit_quality": {"robust_scale_s": 1.05},
    }


def modeling_results() -> dict:
    return {
        "chi_squared": 1.98,
        "reduced_chi_squared": 1.02,
        "dof": 55,
        "background": 0.0015,
        "q_min": 0.005,
        "q_max": 0.35,
        "slit_length": 0.0,
        "timestamp": "2026-08-09T10:20:00",
        "fit_quality": {"robust_scale_s": 1.3},
        "populations": [
            {
                "pop_type": "size_dist", "enabled": True, "population_index": 1,
                "label": "matrix pores", "dist_type": "lognormal",
                "dist_params": {"mean_size": 80.0, "sdeviation": 0.3},
                "form_factor": "sphere", "ff_params": {},
                "structure_factor": "none", "sf_params": {},
                "scale": 0.01, "contrast": 34.5,
                "derived": {"volume_fraction": 0.0101, "vol_mean_r": 82.4},
            },
            {
                "pop_type": "unified_level", "enabled": True,
                "population_index": 2, "label": "large scale",
                "G": 1200.0, "Rg": 350.0, "B": 2e-6, "P": 3.9,
                "RgCO": 0.0, "ETA": 0.0, "PACK": 0.0, "derived": {},
            },
            {
                "pop_type": "size_dist", "enabled": False,
                "population_index": 3, "label": "disabled", "derived": {},
            },
        ],
    }


def all_inputs() -> dict:
    """Every section at once — the widest path through the formatter."""
    return {
        "file_path": "/data/beamline/AeroGel_500C.h5",
        "data_info": data_info(),
        "fit_results": unified_results(),
        "sizes_results": sizes_results(),
        "simple_fit_results": simple_fit_results(),
    }


def strip_timestamp(report: str) -> str:
    """Drop the generation-time line so two runs can be compared."""
    return "\n".join(
        line for line in report.splitlines()
        if "Report generated" not in line
    )
