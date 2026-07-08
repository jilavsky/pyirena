"""
_nxcansas_common.py — shared NXcanSAS output helpers.

Helpers used by both the Data Merge and Data Manipulation I/O modules
(``nxcansas_data_merge`` / ``nxcansas_data_manipulation``) when writing
processed 1-D data back to an NXcanSAS file:

- stripping non-positive / non-finite data points,
- copying an input file while removing existing pyirena result groups,
- replacing the Q/I/Idev/Qdev arrays in-place,
- appending a Qdev dataset to a freshly created file.
"""
from __future__ import annotations

import shutil
from datetime import datetime
from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from pyirena.io.hdf5 import find_matching_groups

# Known pyirena result group paths to strip when copying the source file
PYIRENA_RESULT_GROUPS = [
    'entry/unified_fit_results',
    'entry/sizes_results',
    'entry/simple_fit_results',
    'entry/waxs_peakfit_results',
    'entry/data_merge_results',
    'entry/modeling_results',
    'entry/data_manipulation_results',
]


def strip_nonpositive_intensities(
    q: np.ndarray,
    I: np.ndarray,
    dI: np.ndarray,
    dQ: Optional[np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray], int]:
    """Drop any points with non-positive or non-finite Q or I.

    Subtraction and merge operations can produce negative intensities (e.g.
    at sample-buffer boundaries, or at the edges of DS1 where a flat
    background is subtracted).  Such points are not physically meaningful
    and break downstream analysis (in particular log plots and the Size
    Distribution fit).  Igor Pro's equivalent routines delete these points
    before saving; we do the same.

    Returns the cleaned arrays and the number of points removed.
    """
    q  = np.asarray(q,  dtype=float)
    I  = np.asarray(I,  dtype=float)
    dI = np.asarray(dI, dtype=float)
    mask = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
    n_removed = int(np.sum(~mask))
    if n_removed == 0:
        return q, I, dI, dQ, 0
    q  = q[mask]
    I  = I[mask]
    dI = dI[mask]
    if dQ is not None:
        dQ = np.asarray(dQ, dtype=float)[mask]
    return q, I, dI, dQ, n_removed


def copy_and_strip_results(src: Path, dst: Path) -> None:
    """Copy *src* to *dst* then delete known pyirena result groups."""
    shutil.copy2(src, dst)
    with h5py.File(dst, 'a') as f:
        for grp_path in PYIRENA_RESULT_GROUPS:
            if grp_path in f:
                del f[grp_path]


def replace_nxcansas_data(
    filepath: Path,
    q: np.ndarray,
    I: np.ndarray,
    dI: np.ndarray,
    dQ: Optional[np.ndarray],
    context: str = "data",
) -> None:
    """Overwrite Q/I/Idev/Qdev arrays in the first NXcanSAS sasdata group.

    *context* is used only in the error message (e.g. "merged data",
    "manipulated data") so callers can identify which tool failed.
    """
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
                f"Cannot replace {context}."
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
            # Remove stale Qdev if the new data has none
            if 'Qdev' in sasdata:
                del sasdata['Qdev']
            if 'Q' in sasdata and 'resolutions' in sasdata['Q'].attrs:
                del sasdata['Q'].attrs['resolutions']

        # Update file timestamp
        f.attrs['file_time'] = datetime.now().isoformat()


def append_dq(filepath: Path, dQ: np.ndarray, sample_name: str) -> None:
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
