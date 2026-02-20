"""
NXcanSAS I/O for Sizes (size distribution) results.

Saves and loads size distribution fitting results in the same NXcanSAS HDF5
file that already contains experimental SAS data, mirroring the pattern used
by nxcansas_unified.py.  Results are stored under:

    entry/sizes_results   (NXprocess group)

Datasets stored
---------------
Q                  — Q vector [Å^-1]
intensity_data     — measured intensity [cm^-1]
intensity_error    — measurement uncertainty [cm^-1]  (if available)
intensity_model    — fitted model intensity [cm^-1]
residuals          — normalised residuals (I_data - I_model) / error
r_grid             — radius bin centres [Å]
distribution       — size distribution P(r) [volume fraction / Å]
distribution_std   — per-bin std of P(r) across McSAS repetitions (if available)

Group attributes (fit results)
-------------------------------
chi_squared, volume_fraction, rg, n_iterations, q_power

Group attributes (model / grid setup)
--------------------------------------
shape, contrast, aspect_ratio,
r_min, r_max, n_bins, log_spacing,
background, power_law_B, power_law_P,
method, error_scale,
maxent_sky_background, maxent_stability, maxent_max_iter,
regularization_evalue, regularization_min_ratio,
tnnls_approach_param, tnnls_max_iter,
mcsas_n_repetitions, mcsas_convergence, mcsas_max_iter,
power_law_q_min, power_law_q_max,
background_q_min, background_q_max,
cursor_q_min, cursor_q_max,
timestamp, program
"""

from __future__ import annotations

import h5py
import numpy as np
from datetime import datetime
from pathlib import Path
from typing import Optional


_GROUP = 'entry/sizes_results'
_PROGRAM = 'pyirena.core.sizes'


def save_sizes_results(
    filepath: Path,
    q: np.ndarray,
    intensity_data: np.ndarray,
    intensity_model: np.ndarray,
    residuals: np.ndarray,
    r_grid: np.ndarray,
    distribution: np.ndarray,
    params: dict,
    intensity_error: Optional[np.ndarray] = None,
    distribution_std: Optional[np.ndarray] = None,
) -> None:
    """
    Save size distribution fitting results to an NXcanSAS HDF5 file.

    The file must already exist (created by ``create_nxcansas_file`` or the
    data selector workflow).  If a ``sizes_results`` group already exists it
    is deleted and recreated.

    Args:
        filepath:        Path to the existing NXcanSAS HDF5 file.
        q:               Q vector [Å^-1].
        intensity_data:  Measured intensity [cm^-1].
        intensity_model: Fitted model intensity [cm^-1].
        residuals:       Normalised residuals (I_data - I_model) / error.
        r_grid:          Radius bin centres [Å].
        distribution:    Size distribution P(r) [volume fraction / Å].
        params:          Dict of scalar metadata to store as group attributes.
                         All keys are optional.  Fit results:
                           ``chi_squared``, ``volume_fraction``, ``rg``,
                           ``n_iterations``, ``q_power``.
                         Model / grid setup:
                           ``shape``, ``contrast``, ``aspect_ratio``,
                           ``r_min``, ``r_max``, ``n_bins``, ``log_spacing``,
                           ``background``, ``power_law_B``, ``power_law_P``,
                           ``method``, ``error_scale``,
                           ``maxent_sky_background``, ``maxent_stability``,
                           ``maxent_max_iter``, ``regularization_evalue``,
                           ``regularization_min_ratio``,
                           ``tnnls_approach_param``, ``tnnls_max_iter``,
                           ``mcsas_n_repetitions``,
                           ``mcsas_convergence``, ``mcsas_max_iter``,
                           ``power_law_q_min``, ``power_law_q_max``,
                           ``background_q_min``, ``background_q_max``,
                           ``cursor_q_min``, ``cursor_q_max``.
        intensity_error:  Measurement uncertainty [cm^-1]; stored if provided.
        distribution_std: Per-bin std of P(r) across McSAS repetitions;
                          stored if provided (McSAS only).
    """
    filepath = Path(filepath)
    timestamp = datetime.now().isoformat()

    mode = 'a' if filepath.exists() else 'w'
    with h5py.File(filepath, mode) as f:
        # Remove existing group if present
        if _GROUP in f:
            del f[_GROUP]

        grp = f.require_group(_GROUP)
        grp.attrs['NX_class'] = 'NXprocess'
        grp.attrs['program']   = _PROGRAM
        grp.attrs['timestamp'] = timestamp

        # Store scalar metadata as attributes
        _scalar_keys = (
            # Fit results
            'chi_squared', 'volume_fraction', 'rg', 'n_iterations', 'q_power',
            # Model / grid setup
            'shape', 'contrast', 'aspect_ratio',
            'r_min', 'r_max', 'n_bins', 'log_spacing',
            'background', 'power_law_B', 'power_law_P',
            'method', 'error_scale',
            # Method-specific parameters
            'maxent_sky_background', 'maxent_stability', 'maxent_max_iter',
            'regularization_evalue', 'regularization_min_ratio',
            'tnnls_approach_param', 'tnnls_max_iter',
            'mcsas_n_repetitions',
            'mcsas_convergence', 'mcsas_max_iter',
            # Q ranges used during fitting
            'power_law_q_min', 'power_law_q_max',
            'background_q_min', 'background_q_max',
            'cursor_q_min', 'cursor_q_max',
        )
        for k in _scalar_keys:
            if k in params and params[k] is not None:
                grp.attrs[k] = params[k]

        # Store arrays
        grp.create_dataset('Q',                 data=q.astype('f8'),              compression='gzip')
        grp.create_dataset('intensity_data',    data=intensity_data.astype('f8'), compression='gzip')
        grp.create_dataset('intensity_model',   data=intensity_model.astype('f8'), compression='gzip')
        grp.create_dataset('residuals',         data=residuals.astype('f8'),      compression='gzip')
        grp.create_dataset('r_grid',            data=r_grid.astype('f8'),         compression='gzip')
        grp.create_dataset('distribution',      data=distribution.astype('f8'),   compression='gzip')

        if intensity_error is not None:
            grp.create_dataset('intensity_error', data=intensity_error.astype('f8'), compression='gzip')

        if distribution_std is not None:
            grp.create_dataset('distribution_std', data=distribution_std.astype('f8'), compression='gzip')
            grp['distribution_std'].attrs['units'] = 'volume_fraction/angstrom'

        # Units annotations
        grp['Q'].attrs['units']              = '1/angstrom'
        grp['intensity_data'].attrs['units'] = '1/cm'
        grp['intensity_model'].attrs['units'] = '1/cm'
        grp['r_grid'].attrs['units']         = 'angstrom'
        grp['distribution'].attrs['units']   = 'volume_fraction/angstrom'


def load_sizes_results(filepath: Path) -> dict:
    """
    Load size distribution results from an NXcanSAS HDF5 file.

    Args:
        filepath: Path to the NXcanSAS HDF5 file.

    Returns:
        dict with keys:
            ``Q``, ``intensity_data``, ``intensity_model``, ``residuals``,
            ``r_grid``, ``distribution``,
            ``intensity_error`` (may be None),
            and all scalar metadata stored as group attributes
            (fit results: ``chi_squared``, ``volume_fraction``, ``rg``,
             ``n_iterations``, ``q_power``; model setup: ``shape``,
             ``contrast``, ``aspect_ratio``, ``r_min``, ``r_max``,
             ``n_bins``, ``log_spacing``, ``background``, ``power_law_B``,
             ``power_law_P``, ``method``, ``error_scale``,
             ``maxent_sky_background``, ``maxent_stability``,
             ``maxent_max_iter``, ``regularization_evalue``,
             ``regularization_min_ratio``, ``tnnls_approach_param``,
             ``tnnls_max_iter``, ``power_law_q_min``, ``power_law_q_max``,
             ``background_q_min``, ``background_q_max``,
             ``cursor_q_min``, ``cursor_q_max``,
             ``timestamp``, ``program``).

    Raises:
        KeyError:  if the file does not contain a ``sizes_results`` group.
        OSError:   if the file cannot be opened.
    """
    filepath = Path(filepath)
    with h5py.File(filepath, 'r') as f:
        if _GROUP not in f:
            raise KeyError(
                f"No sizes_results group found in {filepath.name}. "
                "Run the Sizes fitting tool first."
            )
        grp = f[_GROUP]

        result: dict = {}

        # Arrays
        result['Q']              = grp['Q'][:]
        result['intensity_data'] = grp['intensity_data'][:]
        result['intensity_model'] = grp['intensity_model'][:]
        result['residuals']      = grp['residuals'][:]
        result['r_grid']         = grp['r_grid'][:]
        result['distribution']   = grp['distribution'][:]

        if 'intensity_error' in grp:
            result['intensity_error'] = grp['intensity_error'][:]
        else:
            result['intensity_error'] = None

        if 'distribution_std' in grp:
            result['distribution_std'] = grp['distribution_std'][:]
        else:
            result['distribution_std'] = None

        # Scalar metadata from attributes
        for k in (
            # Fit results
            'chi_squared', 'volume_fraction', 'rg', 'n_iterations', 'q_power',
            # Model / grid setup
            'shape', 'contrast', 'aspect_ratio',
            'r_min', 'r_max', 'n_bins', 'log_spacing',
            'background', 'power_law_B', 'power_law_P',
            'method', 'error_scale',
            # Method-specific parameters
            'maxent_sky_background', 'maxent_stability', 'maxent_max_iter',
            'regularization_evalue', 'regularization_min_ratio',
            'tnnls_approach_param', 'tnnls_max_iter',
            'mcsas_n_repetitions',
            'mcsas_convergence', 'mcsas_max_iter',
            # Q ranges used during fitting
            'power_law_q_min', 'power_law_q_max',
            'background_q_min', 'background_q_max',
            'cursor_q_min', 'cursor_q_max',
            # Administrative
            'timestamp', 'program',
        ):
            result[k] = grp.attrs.get(k)

    return result


def print_sizes_results(results: dict) -> None:
    """Pretty-print size distribution results to the console."""
    print("=" * 60)
    print("  Size Distribution Fitting Results")
    print("=" * 60)
    print(f"  Shape:           {results.get('shape', 'N/A')}")
    print(f"  Method:          {results.get('method', 'N/A')}")
    print(f"  Contrast (Δρ)²:  {results.get('contrast', 'N/A'):.4g} ×10²⁰ cm⁻⁴")
    print(f"  Chi-squared:     {results.get('chi_squared', 'N/A'):.4f}")
    print(f"  Volume fraction: {results.get('volume_fraction', 'N/A'):.4g}")
    print(f"  Rg:              {results.get('rg', 'N/A'):.2f} Å")
    print(f"  Iterations:      {results.get('n_iterations', 'N/A')}")
    r = results.get('r_grid')
    if r is not None:
        print(f"  R range:         {r.min():.1f} – {r.max():.1f} Å")
    print(f"  Timestamp:       {results.get('timestamp', 'N/A')}")
    print("=" * 60)
