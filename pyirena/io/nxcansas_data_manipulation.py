"""
nxcansas_data_manipulation.py — HDF5 / NXcanSAS I/O for the Data Manipulation tool.

Handles:
- Copying an input NXcanSAS file and replacing its Q/I/Idev/Qdev arrays with
  the manipulated data while stripping any existing pyirena result groups.
- Creating a fresh NXcanSAS file when the input is not NXcanSAS.
- Appending a provenance NXprocess group recording the operation parameters.
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


# Known pyirena result group paths to strip when copying the source file
_PYIRENA_RESULT_GROUPS = [
    'entry/unified_fit_results',
    'entry/sizes_results',
    'entry/simple_fits_results',
    'entry/waxs_peakfit_results',
    'entry/data_merge_results',
    'entry/modeling_results',
    'entry/data_manipulation_results',
]

# Default suffix per operation
_OPERATION_SUFFIXES = {
    'scaled': '_scaled',
    'trimmed': '_trimmed',
    'rebinned': '_rebinned',
    'avg': '_avg',
    'sub': '_sub',
    'div': '_div',
}


def save_manipulated_data(
    output_folder: Path,
    source_path: Path,
    source_is_nxcansas: bool,
    q: np.ndarray,
    I: np.ndarray,
    dI: np.ndarray,
    dQ: Optional[np.ndarray],
    operation: str,
    provenance: dict,
    output_stem_suffix: Optional[str] = None,
) -> Path:
    """Save manipulated SAS data to a NXcanSAS HDF5 file.

    Strategy
    --------
    - If *source_is_nxcansas*: copy source → output dir (stripping pyirena
      results), then replace Q/I/Idev/Qdev in the existing sasdata group.
    - Otherwise: create a fresh NXcanSAS file via ``create_nxcansas_file()``.

    In both cases a ``data_manipulation_results`` NXprocess group is appended
    with all operation parameters for provenance.

    Parameters
    ----------
    output_folder : Path
        Directory where the output file will be written.  Must already exist.
    source_path : Path
        Full path to the source input file.
    source_is_nxcansas : bool
        True if the source was read as a proper NXcanSAS HDF5 file.
    q, I, dI : ndarray
        Manipulated Q, intensity, and uncertainty arrays.
    dQ : ndarray or None
        Q-resolution array, or None.
    operation : str
        Operation key (e.g. 'scaled', 'trimmed', 'avg').
    provenance : dict
        Operation-specific parameters to store as datasets.
    output_stem_suffix : str or None
        Override the default suffix.  If None, looks up *operation* in
        ``_OPERATION_SUFFIXES``.

    Returns
    -------
    Path
        Full path of the written output file.
    """
    output_folder = Path(output_folder)
    source_path = Path(source_path)

    suffix = output_stem_suffix or _OPERATION_SUFFIXES.get(operation, f'_{operation}')
    stem = source_path.stem
    ext = source_path.suffix if source_is_nxcansas else '.h5'
    out_name = f"{stem}{suffix}{ext}"
    out_path = output_folder / out_name

    output_folder.mkdir(parents=True, exist_ok=True)

    if source_is_nxcansas:
        _copy_and_strip_results(source_path, out_path)
        _replace_nxcansas_data(out_path, q, I, dI, dQ)
    else:
        sample_name = stem
        create_nxcansas_file(out_path, q, I, error=dI, sample_name=sample_name)
        if dQ is not None:
            _append_dq(out_path, dQ, sample_name)

    _append_manipulation_provenance(out_path, operation, provenance, source_path)
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
        sasdata_paths = find_matching_groups(
            f,
            required_attributes={'canSAS_class': 'SASdata'},
            required_items={},
        )
        if not sasdata_paths:
            sasdata_paths = find_matching_groups(
                f,
                required_attributes={'NX_class': 'NXdata'},
                required_items={},
            )
        if not sasdata_paths:
            raise RuntimeError(
                f"Could not locate a sasdata/NXdata group in {filepath}. "
                "Cannot replace manipulated data."
            )

        sasdata = f[sasdata_paths[0]]

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

        if 'I' in sasdata:
            sasdata['I'].attrs['uncertainties'] = 'Idev'

        if dQ is not None:
            if 'Qdev' in sasdata:
                del sasdata['Qdev']
            ds_qdev = sasdata.create_dataset('Qdev', data=dQ)
            ds_qdev.attrs['units'] = '1/angstrom'
            ds_qdev.attrs['long_name'] = 'Q resolution'
            sasdata['Q'].attrs['resolutions'] = 'Qdev'
        else:
            if 'Qdev' in sasdata:
                del sasdata['Qdev']
            if 'Q' in sasdata and 'resolutions' in sasdata['Q'].attrs:
                del sasdata['Q'].attrs['resolutions']

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


def _append_manipulation_provenance(
    filepath: Path,
    operation: str,
    provenance: dict,
    source_path: Path,
) -> None:
    """Append ``entry/data_manipulation_results`` NXprocess group."""
    with h5py.File(filepath, 'a') as f:
        grp_path = 'entry/data_manipulation_results'
        if grp_path in f:
            del f[grp_path]

        grp = f.require_group(grp_path)
        grp.attrs['NX_class'] = 'NXprocess'
        grp.attrs['analysis_type'] = 'Data Manipulation'
        grp.attrs['program'] = 'pyirena'
        grp.attrs['version'] = '1.0'
        grp.attrs['timestamp'] = datetime.now().isoformat()

        grp.create_dataset('operation', data=operation)
        grp.create_dataset('source_file', data=str(source_path))

        for key, value in provenance.items():
            if isinstance(value, bool):
                grp.create_dataset(key, data=int(value))
            elif value is None:
                grp.create_dataset(key, data=float('nan'))
            elif isinstance(value, str):
                grp.create_dataset(key, data=value)
            elif isinstance(value, (int, float, np.integer, np.floating)):
                grp.create_dataset(key, data=float(value))
