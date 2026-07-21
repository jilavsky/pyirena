"""
nxcansas_data_merge.py — HDF5 / NXcanSAS I/O for the Data Merge tool.

Handles:
- Copying an input NXcanSAS file and replacing its Q/I/Idev/Qdev arrays with
  the merged data while stripping any existing pyirena result groups.
- Creating a fresh NXcanSAS file when the DS1 input is not already NXcanSAS.
- Appending a provenance NXprocess group recording the merge parameters.
"""
from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from pyirena.io.nxcansas_unified import create_nxcansas_file
from pyirena.io._nxcansas_common import (
    strip_nonpositive_intensities as _strip_nonpositive_intensities,
    copy_and_strip_results as _copy_and_strip_results,
    replace_nxcansas_data,
    append_dq as _append_dq,
    append_dql as _append_dql,
    drop_smr_entries as _drop_smr_entries,
)

import logging

log = logging.getLogger(__name__)



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

    # Strip points with non-positive / non-finite Q or I before writing.
    # See _strip_nonpositive_intensities() for the rationale.
    q, I, dI, dQ, n_stripped = _strip_nonpositive_intensities(q, I, dI, dQ)
    if n_stripped:
        log.info(f"[data_merge] Stripped {n_stripped} non-positive/non-finite "
                 f"point(s) before saving merged data.")
    if len(q) < 2:
        raise ValueError(
            f"After stripping non-positive intensities only {len(q)} point(s) "
            f"remain — refusing to save an empty/degenerate merged dataset."
        )

    # Determine output filename
    stem = ds1_path.stem
    ext = ds1_path.suffix if ds1_is_nxcansas else '.h5'
    out_name = f"{stem}{output_stem_suffix}{ext}"
    out_path = output_folder / out_name

    output_folder.mkdir(parents=True, exist_ok=True)

    if ds1_is_nxcansas:
        # Copy DS1 file, strip results, replace data arrays
        _copy_and_strip_results(ds1_path, out_path)
        replace_nxcansas_data(out_path, q, I, dI, dQ, context="merged data")
        # The copied file may contain a stale slit-smeared (_SMR) twin of DS1
        # that no longer matches the merged default entry — drop it so a later
        # prefer_slit_smeared load can't silently return the wrong curve.
        n_smr = _drop_smr_entries(out_path)
        if n_smr:
            log.info(f"[data_merge] Removed {n_smr} stale _SMR entry(ies) "
                     f"from the merged output.")
    else:
        # Create fresh NXcanSAS file
        sample_name = stem
        create_nxcansas_file(out_path, q, I, error=dI, sample_name=sample_name)
        # Add dQ if available
        if dQ is not None:
            _append_dq(out_path, dQ, sample_name)

    # When the merged curve is slit smeared, mark the output so downstream
    # tools auto-detect it (writes scalar dQl + Q@resolutions).
    sl_merged = float(merge_result_dict.get('slit_length_merged', 0.0) or 0.0)
    if sl_merged > 0:
        _append_dql(out_path, sl_merged)

    # Append provenance group
    _append_merge_provenance(out_path, merge_result_dict, ds1_path, ds2_path)

    return out_path


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------




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
