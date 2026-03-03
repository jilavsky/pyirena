"""
nxcansas_data_merge.py — HDF5 / NXcanSAS I/O for the Data Merge tool.

Handles:
- Copying an input NXcanSAS file and replacing its Q/I/Idev/Qdev arrays with
  the merged data while stripping any existing pyirena result groups.
- Creating a fresh NXcanSAS file when the DS1 input is not already NXcanSAS.
- Appending a provenance NXprocess group recording the merge parameters.
"""
from __future__ import annotations

import shutil
from datetime import datetime
from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from pyirena.io.nxcansas_unified import create_nxcansas_file
from pyirena.io.hdf5 import find_matching_groups


# Known pyirena result group paths to strip when copying DS1
_PYIRENA_RESULT_GROUPS = [
    'entry/unified_fit_results',
    'entry/sizes_results',
    'entry/simple_fits_results',
    'entry/waxs_peakfit_results',
    'entry/data_merge_results',
]


def save_merged_data(
    output_folder: Path,
    ds1_path: Path,
    ds1_is_nxcansas: bool,
    q: np.ndarray,
    I: np.ndarray,
    dI: np.ndarray,
    dQ: Optional[np.ndarray],
    merge_result_dict: dict,
    ds2_path: Optional[Path] = None,
    output_stem_suffix: str = '_merged',
) -> Path:
    """Save merged SAS data to a NXcanSAS HDF5 file.

    Strategy
    --------
    - If *ds1_is_nxcansas*: copy DS1 → output dir (stripping pyirena results),
      then replace the Q/I/Idev/Qdev arrays in the existing sasdata group.
    - Otherwise: create a fresh NXcanSAS file via ``create_nxcansas_file()``.

    In both cases a ``data_merge_results`` NXprocess group is appended with all
    merge parameters for provenance.

    Parameters
    ----------
    output_folder : Path
        Directory where the output file will be written.  Must already exist.
    ds1_path : Path
        Full path to the DS1 input file.
    ds1_is_nxcansas : bool
        True if DS1 was read as a proper NXcanSAS HDF5 file.
    q, I, dI : ndarray
        Merged Q, intensity, and uncertainty arrays.
    dQ : ndarray or None
        Merged Q-resolution array, or None.
    merge_result_dict : dict
        Keys: scale, q_shift, background, chi_squared, q_overlap_min,
        q_overlap_max, scale_dataset, fit_scale, fit_qshift,
        split_at_left_cursor.  Values are scalars or booleans.
    ds2_path : Path or None
        Full path to DS2 input file (for provenance only).
    output_stem_suffix : str
        Appended to DS1 stem for the output filename.

    Returns
    -------
    Path
        Full path of the written output file.
    """
    output_folder = Path(output_folder)
    ds1_path = Path(ds1_path)

    # Determine output filename
    stem = ds1_path.stem
    ext = ds1_path.suffix if ds1_is_nxcansas else '.h5'
    out_name = f"{stem}{output_stem_suffix}{ext}"
    out_path = output_folder / out_name

    output_folder.mkdir(parents=True, exist_ok=True)

    if ds1_is_nxcansas:
        # Copy DS1 file, strip results, replace data arrays
        _copy_and_strip_results(ds1_path, out_path)
        _replace_nxcansas_data(out_path, q, I, dI, dQ)
    else:
        # Create fresh NXcanSAS file
        sample_name = stem
        create_nxcansas_file(out_path, q, I, error=dI, sample_name=sample_name)
        # Add dQ if available
        if dQ is not None:
            _append_dq(out_path, dQ, sample_name)

    # Append provenance group
    _append_merge_provenance(out_path, merge_result_dict, ds1_path, ds2_path)

    return out_path


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _copy_and_strip_results(src: Path, dst: Path) -> None:
    """Copy *src* to *dst* then delete known pyirena result groups."""
    shutil.copy2(src, dst)
    with h5py.File(dst, 'a') as f:
        for grp_path in _PYIRENA_RESULT_GROUPS:
            if grp_path in f:
                del f[grp_path]


def _replace_nxcansas_data(
    filepath: Path,
    q: np.ndarray,
    I: np.ndarray,
    dI: np.ndarray,
    dQ: Optional[np.ndarray],
) -> None:
    """Overwrite Q/I/Idev/Qdev arrays in the first NXcanSAS sasdata group."""
    with h5py.File(filepath, 'a') as f:
        # Locate the sasdata group dynamically (same approach as the reader)
        sasdata_paths = find_matching_groups(
            f,
            required_attributes={'canSAS_class': 'SASdata'},
            required_items={},
        )
        if not sasdata_paths:
            # Fallback: look for NXdata group
            sasdata_paths = find_matching_groups(
                f,
                required_attributes={'NX_class': 'NXdata'},
                required_items={},
            )
        if not sasdata_paths:
            raise RuntimeError(
                f"Could not locate a sasdata/NXdata group in {filepath}. "
                "Cannot replace merged data."
            )

        sasdata = f[sasdata_paths[0]]

        # Replace each array dataset; delete-then-recreate to allow size change
        for name, data, attrs in [
            ('Q',    q,  {'units': '1/angstrom', 'long_name': 'Q'}),
            ('I',    I,  {'units': '1/cm',       'long_name': 'Intensity'}),
            ('Idev', dI, {'units': '1/cm',       'long_name': 'Uncertainties'}),
        ]:
            if name in sasdata:
                del sasdata[name]
            ds = sasdata.create_dataset(name, data=data)
            for k, v in attrs.items():
                ds.attrs[k] = v

        # I.uncertainties attribute
        if 'I' in sasdata:
            sasdata['I'].attrs['uncertainties'] = 'Idev'

        # Qdev (Q resolution) — add only if provided
        if dQ is not None:
            if 'Qdev' in sasdata:
                del sasdata['Qdev']
            ds_qdev = sasdata.create_dataset('Qdev', data=dQ)
            ds_qdev.attrs['units'] = '1/angstrom'
            ds_qdev.attrs['long_name'] = 'Q resolution'
            sasdata['Q'].attrs['resolutions'] = 'Qdev'
        else:
            # Remove stale Qdev if the merged data has none
            if 'Qdev' in sasdata:
                del sasdata['Qdev']
            if 'Q' in sasdata and 'resolutions' in sasdata['Q'].attrs:
                del sasdata['Q'].attrs['resolutions']

        # Update file timestamp
        f.attrs['file_time'] = datetime.now().isoformat()


def _append_dq(filepath: Path, dQ: np.ndarray, sample_name: str) -> None:
    """Add a Qdev dataset to the sasdata group of a freshly-created NXcanSAS file."""
    with h5py.File(filepath, 'a') as f:
        sasdata_paths = find_matching_groups(
            f,
            required_attributes={'canSAS_class': 'SASdata'},
            required_items={},
        )
        if sasdata_paths:
            sasdata = f[sasdata_paths[0]]
            if 'Qdev' not in sasdata:
                ds = sasdata.create_dataset('Qdev', data=dQ)
                ds.attrs['units'] = '1/angstrom'
                ds.attrs['long_name'] = 'Q resolution'
                if 'Q' in sasdata:
                    sasdata['Q'].attrs['resolutions'] = 'Qdev'


def _append_merge_provenance(
    filepath: Path,
    merge_result_dict: dict,
    ds1_path: Path,
    ds2_path: Optional[Path],
) -> None:
    """Append ``entry/data_merge_results`` NXprocess group with merge parameters."""
    with h5py.File(filepath, 'a') as f:
        # Remove stale group if present (e.g. re-saving after re-merge)
        if 'entry/data_merge_results' in f:
            del f['entry/data_merge_results']

        grp = f.require_group('entry/data_merge_results')
        grp.attrs['NX_class'] = 'NXprocess'
        grp.attrs['program'] = 'pyirena'
        grp.attrs['version'] = '1.0'
        grp.attrs['timestamp'] = datetime.now().isoformat()

        # Store source filenames as string datasets (visible in HDF5 viewer)
        # and as group attributes (for quick programmatic access).
        ds1_str = str(ds1_path)
        ds2_str = str(ds2_path) if ds2_path else ''
        grp.create_dataset('ds1_file', data=ds1_str)
        grp.create_dataset('ds2_file', data=ds2_str)
        grp.attrs['ds1_path'] = ds1_str   # keep attr for backwards compatibility
        grp.attrs['ds2_path'] = ds2_str

        # Store merge parameters as scalar datasets (visible in HDF5 viewer)
        for key, value in merge_result_dict.items():
            if isinstance(value, bool):
                grp.create_dataset(key, data=int(value))
            elif value is None:
                grp.create_dataset(key, data=float('nan'))
            else:
                grp.create_dataset(key, data=value)
