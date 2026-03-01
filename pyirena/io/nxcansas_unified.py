"""
NXcanSAS I/O for Unified Fit results.

This module handles saving and loading Unified Fit results to/from NXcanSAS HDF5 files.
It creates proper NeXus structure following NXcanSAS standard and stores Unified Fit
model parameters, calculated intensities, and residuals.
"""

import h5py
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple


def create_nxcansas_file(filepath: Path, q: np.ndarray, intensity: np.ndarray,
                         error: Optional[np.ndarray] = None,
                         sample_name: str = "data") -> None:
    """
    Create a new NXcanSAS HDF5 file with experimental data.

    Args:
        filepath: Path to output HDF5 file
        q: Q vector (1/Angstrom)
        intensity: Intensity data (1/cm)
        error: Error/uncertainty data (1/cm)
        sample_name: Name for the data entry
    """
    timestamp = datetime.now().isoformat()

    with h5py.File(filepath, "w") as f:
        # Root attributes
        f.attrs['default'] = 'entry'
        f.attrs['file_name'] = filepath.name
        f.attrs['file_time'] = timestamp
        f.attrs['creator'] = 'pyirena'
        f.attrs['NeXus_version'] = '4.3.0'
        f.attrs['HDF5_version'] = h5py.version.hdf5_version
        f.attrs['h5py_version'] = h5py.version.version

        # Create entry group
        nxentry = f.create_group('entry')
        nxentry.attrs['NX_class'] = 'NXentry'
        nxentry.attrs['canSAS_class'] = 'SASentry'
        nxentry.attrs['default'] = sample_name
        nxentry.create_dataset('definition', data='NXsas')

        # Create subentry for data
        data_path = f"entry/{sample_name}"
        nxdata_entry = f.create_group(data_path)
        nxdata_entry.attrs['NX_class'] = 'NXsubentry'
        nxdata_entry.attrs['canSAS_class'] = 'SASentry'
        nxdata_entry.attrs['default'] = 'sasdata'
        nxdata_entry.attrs['title'] = sample_name
        nxdata_entry.create_dataset('definition', data='NXcanSAS')
        nxdata_entry.create_dataset('title', data=sample_name)
        nxdata_entry.create_dataset('run', data='pyirena_run')

        # Create sasdata group
        sasdata = nxdata_entry.create_group('sasdata')
        sasdata.attrs['NX_class'] = 'NXdata'
        sasdata.attrs['canSAS_class'] = 'SASdata'
        sasdata.attrs['signal'] = 'I'
        sasdata.attrs['I_axes'] = 'Q'

        # Store intensity data
        ds_i = sasdata.create_dataset('I', data=intensity)
        ds_i.attrs['units'] = '1/cm'
        ds_i.attrs['long_name'] = 'Intensity'
        if error is not None:
            ds_i.attrs['uncertainties'] = 'Idev'

        # Store Q data
        ds_q = sasdata.create_dataset('Q', data=q)
        ds_q.attrs['units'] = '1/angstrom'
        ds_q.attrs['long_name'] = 'Q'

        # Store error data if provided
        if error is not None:
            ds_err = sasdata.create_dataset('Idev', data=error)
            ds_err.attrs['units'] = '1/cm'
            ds_err.attrs['long_name'] = 'Uncertainties'


def save_unified_fit_results(filepath: Path,
                             q: np.ndarray,
                             intensity_data: np.ndarray,
                             intensity_model: np.ndarray,
                             residuals: np.ndarray,
                             levels: List[Dict],
                             background: float,
                             chi_squared: float,
                             num_levels: int,
                             error: Optional[np.ndarray] = None,
                             uncertainties: Optional[Dict] = None) -> None:
    """
    Save Unified Fit results to NXcanSAS HDF5 file.

    Creates or appends to NXcanSAS file with Unified Fit results including
    model parameters, calculated intensity, and residuals.

    Args:
        filepath: Path to HDF5 file (will be created if doesn't exist)
        q: Q vector (1/Angstrom)
        intensity_data: Experimental intensity (1/cm)
        intensity_model: Calculated model intensity (1/cm)
        residuals: Fit residuals (normalized)
        levels: List of level dictionaries with parameters (G, Rg, B, P, etc.)
        background: Background value
        chi_squared: Chi-squared value from fit
        num_levels: Number of levels used
        error: Error/uncertainty data (1/cm), optional
        uncertainties: Monte Carlo parameter uncertainties (std devs), optional.
            Structure: {'levels': [{'G': float, 'Rg': float, ...}, ...],
                        'background': float}
            Non-zero values are stored as  <param>_err  attributes on each level
            group and as  background_err  on the unified_fit_results group.
    """
    timestamp = datetime.now().isoformat()

    # Open or create file
    with h5py.File(filepath, "a") as f:
        # Ensure entry exists
        if 'entry' not in f:
            # Create basic NXcanSAS structure
            f.attrs['default'] = 'entry'
            f.attrs['file_name'] = filepath.name
            f.attrs['file_time'] = timestamp
            f.attrs['creator'] = 'pyirena'
            f.attrs['NeXus_version'] = '4.3.0'
            f.attrs['HDF5_version'] = h5py.version.hdf5_version
            f.attrs['h5py_version'] = h5py.version.version

            nxentry = f.create_group('entry')
            nxentry.attrs['NX_class'] = 'NXentry'
            nxentry.attrs['canSAS_class'] = 'SASentry'
            nxentry.create_dataset('definition', data='NXsas')

        nxentry = f['entry']

        # Create or replace Unified Fit results group
        unified_path = 'entry/unified_fit_results'
        if unified_path in f:
            del f[unified_path]

        unified_group = f.create_group(unified_path)
        unified_group.attrs['NX_class'] = 'NXprocess'
        unified_group.attrs['analysis_type'] = 'Unified Fit'
        unified_group.attrs['program'] = 'pyirena'
        unified_group.attrs['timestamp'] = timestamp
        unified_group.attrs['num_levels'] = num_levels
        # Store fit results as scalar datasets (browseable/collectable in HDF5 viewer)
        unified_group.create_dataset('background', data=float(background))
        unified_group.create_dataset('chi_squared', data=float(chi_squared))
        if uncertainties is not None:
            bg_err = uncertainties.get('background', 0.0)
            if bg_err > 0.0:
                unified_group.create_dataset('background_err', data=float(bg_err))

        # Store Q vector
        ds_q = unified_group.create_dataset('Q', data=q)
        ds_q.attrs['units'] = '1/angstrom'
        ds_q.attrs['long_name'] = 'Q vector'

        # Store experimental data
        ds_data = unified_group.create_dataset('intensity_data', data=intensity_data)
        ds_data.attrs['units'] = '1/cm'
        ds_data.attrs['long_name'] = 'Experimental intensity'

        if error is not None:
            ds_err = unified_group.create_dataset('intensity_error', data=error)
            ds_err.attrs['units'] = '1/cm'
            ds_err.attrs['long_name'] = 'Experimental uncertainty'

        # Store model intensity
        ds_model = unified_group.create_dataset('intensity_model', data=intensity_model)
        ds_model.attrs['units'] = '1/cm'
        ds_model.attrs['long_name'] = 'Unified Fit model intensity'

        # Store residuals
        ds_resid = unified_group.create_dataset('residuals', data=residuals)
        ds_resid.attrs['long_name'] = 'Fit residuals (normalized)'

        # Store level parameters (and MC uncertainties if available)
        for i, level_params in enumerate(levels):
            level_num = i + 1
            level_group = unified_group.create_group(f'level_{level_num}')
            level_group.attrs['level_number'] = level_num

            for param_name, param_value in level_params.items():
                if isinstance(param_value, (int, float)) and not isinstance(param_value, bool):
                    # Numeric params as scalar datasets (browseable/collectable)
                    level_group.create_dataset(param_name, data=float(param_value))
                elif isinstance(param_value, str):
                    # String params stay as attributes
                    level_group.attrs[param_name] = param_value

            # Store MC uncertainties as <param>_err scalar datasets (non-zero only)
            if uncertainties is not None and i < len(uncertainties.get('levels', [])):
                ud = uncertainties['levels'][i]
                for param_key in ('G', 'Rg', 'B', 'P', 'ETA', 'PACK'):
                    err_val = ud.get(param_key, 0.0)
                    if err_val > 0.0:
                        level_group.create_dataset(f'{param_key}_err', data=float(err_val))

        print(f"Saved Unified Fit results to {filepath}")


def load_unified_fit_results(filepath: Path) -> Dict:
    """
    Load Unified Fit results from NXcanSAS HDF5 file.

    Args:
        filepath: Path to HDF5 file containing Unified Fit results

    Returns:
        Dictionary containing:
            - Q: Q vector
            - intensity_data: Experimental intensity
            - intensity_model: Model intensity
            - residuals: Fit residuals
            - intensity_error: Experimental uncertainty (if available)
            - num_levels: Number of levels
            - background: Background value
            - chi_squared: Chi-squared value
            - levels: List of dictionaries with level parameters
            - timestamp: When fit was saved
    """
    results = {}

    with h5py.File(filepath, "r") as f:
        if 'entry/unified_fit_results' not in f:
            raise ValueError(f"No Unified Fit results found in {filepath}")

        unified = f['entry/unified_fit_results']

        # Load global attributes
        results['num_levels'] = unified.attrs['num_levels']
        results['timestamp'] = unified.attrs['timestamp']
        results['program'] = unified.attrs.get('program', 'unknown')
        # background and chi_squared: try dataset (new format) then attr (old format)
        results['background'] = (
            float(unified['background'][()]) if 'background' in unified
            else unified.attrs.get('background', 0.0)
        )
        results['chi_squared'] = (
            float(unified['chi_squared'][()]) if 'chi_squared' in unified
            else unified.attrs.get('chi_squared', np.nan)
        )

        # Load arrays
        results['Q'] = unified['Q'][:]
        results['intensity_data'] = unified['intensity_data'][:]
        results['intensity_model'] = unified['intensity_model'][:]
        results['residuals'] = unified['residuals'][:]

        if 'intensity_error' in unified:
            results['intensity_error'] = unified['intensity_error'][:]
        else:
            results['intensity_error'] = None

        # Load level parameters
        results['levels'] = []
        for i in range(results['num_levels']):
            level_num = i + 1
            level_path = f'level_{level_num}'

            if level_path in unified:
                level_group = unified[level_path]
                level_params = {}

                # New format: numeric params stored as scalar datasets
                for key, item in level_group.items():
                    if isinstance(item, h5py.Dataset) and item.shape == ():
                        try:
                            level_params[key] = float(item[()])
                        except (TypeError, ValueError):
                            level_params[key] = item[()]

                # Backward compat: also read attrs not yet covered as datasets
                for key in level_group.attrs.keys():
                    if key not in level_params:
                        level_params[key] = level_group.attrs[key]

                results['levels'].append(level_params)

    print(f"Loaded Unified Fit results from {filepath}")
    print(f"  - {results['num_levels']} levels")
    print(f"  - χ² = {results['chi_squared']:.4f}")
    print(f"  - Timestamp: {results['timestamp']}")

    return results


def get_output_filepath(input_path: Path, is_nxcansas: bool = False) -> Path:
    """
    Determine output filepath for storing results.

    If input is already NXcanSAS, use it. Otherwise, create new file
    with "_NX.h5" suffix.

    Args:
        input_path: Original input file path
        is_nxcansas: Whether input is already NXcanSAS format

    Returns:
        Path to use for output file

    Examples:
        mydata.dat (not NXcanSAS) -> mydata_NX.h5
        mydata.txt (not NXcanSAS) -> mydata_NX.h5
        mydata.h5 (not NXcanSAS) -> mydata_NX.h5
        mydata.h5 (is NXcanSAS) -> mydata.h5 (same file)
    """
    if is_nxcansas:
        # Already NXcanSAS - use same file
        return input_path
    else:
        # Create new file with _NX.h5 extension
        # For .dat, .txt, or non-NXcanSAS .h5 files
        stem = input_path.stem
        new_name = f"{stem}_NX.h5"
        return input_path.parent / new_name


def print_unified_fit_results(results: Dict) -> None:
    """
    Pretty print Unified Fit results to console.

    Args:
        results: Dictionary from load_unified_fit_results()
    """
    print("\n" + "="*70)
    print("UNIFIED FIT RESULTS")
    print("="*70)
    print(f"Program: {results['program']}")
    print(f"Timestamp: {results['timestamp']}")
    print(f"Number of levels: {results['num_levels']}")
    print(f"Background: {results['background']:.3e} (1/cm)")
    print(f"Chi-squared: {results['chi_squared']:.4f}")
    print(f"Data points: {len(results['Q'])}")
    print()

    for i, level in enumerate(results['levels']):
        level_num = i + 1
        print(f"Level {level_num}:")
        print(f"  G  = {level.get('G', 0):.3e}")
        print(f"  Rg = {level.get('Rg', 0):.3e} Å")
        print(f"  B  = {level.get('B', 0):.3e}")
        print(f"  P  = {level.get('P', 0):.3f}")

        if level.get('RgCutoff', 0) > 0.01:
            print(f"  RgCutoff = {level['RgCutoff']:.3e} Å")

        if level.get('correlated', False):
            print(f"  ETA  = {level.get('ETA', 0):.1f} Å")
            print(f"  PACK = {level.get('PACK', 0):.2f}")

        if 'Sv' in level and level['Sv'] != 'N/A':
            print(f"  Sv = {level['Sv']:.3e} m²/cm³")

        if 'Invariant' in level and level['Invariant'] != 'N/A':
            print(f"  Invariant = {level['Invariant']:.3e} cm⁻⁴")

        print()

    print("="*70)


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python -m pyirena.io.nxcansas_unified <filepath>")
        print("       Loads and displays Unified Fit results from HDF5 file")
        sys.exit(1)

    filepath = Path(sys.argv[1])

    if not filepath.exists():
        print(f"Error: File not found: {filepath}")
        sys.exit(1)

    try:
        results = load_unified_fit_results(filepath)
        print_unified_fit_results(results)
    except Exception as e:
        print(f"Error loading file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
