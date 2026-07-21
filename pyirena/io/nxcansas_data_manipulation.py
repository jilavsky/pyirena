"""
nxcansas_data_manipulation.py — HDF5 / NXcanSAS I/O for the Data Manipulation tool.

Handles:
- Copying an input NXcanSAS file and replacing its Q/I/Idev/Qdev arrays with
  the manipulated data while stripping any existing pyirena result groups.
- Creating a fresh NXcanSAS file when the input is not NXcanSAS.
- Appending a provenance NXprocess group recording the operation parameters.
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
    slit_length: float = 0.0,
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

    # Strip points with non-positive / non-finite Q or I before writing.
    # See _strip_nonpositive_intensities() for the rationale.
    q, I, dI, dQ, n_stripped = _strip_nonpositive_intensities(q, I, dI, dQ)
    if n_stripped:
        log.info(f"[data_manipulation] Stripped {n_stripped} non-positive/non-finite "
                 f"point(s) before saving '{operation}' result.")
    if len(q) < 2:
        raise ValueError(
            f"After stripping non-positive intensities only {len(q)} point(s) "
            f"remain — refusing to save an empty/degenerate '{operation}' result."
        )

    suffix = output_stem_suffix or _OPERATION_SUFFIXES.get(operation, f'_{operation}')
    stem = source_path.stem
    ext = source_path.suffix if source_is_nxcansas else '.h5'
    out_name = f"{stem}{suffix}{ext}"
    out_path = output_folder / out_name

    output_folder.mkdir(parents=True, exist_ok=True)

    if source_is_nxcansas:
        _copy_and_strip_results(source_path, out_path)
        replace_nxcansas_data(out_path, q, I, dI, dQ, context="manipulated data")
        # A copied Matilda file may contain a stale slit-smeared (_SMR) twin
        # of the source that no longer matches the manipulated default entry —
        # drop it so a later prefer_slit_smeared load can't silently return the
        # wrong curve (F2).
        n_smr = _drop_smr_entries(out_path)
        if n_smr:
            log.info(f"[data_manipulation] Removed {n_smr} stale _SMR entry(ies) "
                     f"from the '{operation}' output.")
    else:
        sample_name = stem
        create_nxcansas_file(out_path, q, I, error=dI, sample_name=sample_name)
        if dQ is not None:
            _append_dq(out_path, dQ, sample_name)

    # When the manipulated curve is itself slit smeared (e.g. a single-entry
    # slit-smeared source, or subtract/divide of two matching smeared curves),
    # re-mark the output.  replace_nxcansas_data() dropped the source's stale
    # dQl, so this restores a *consistent* one for the new arrays.
    if slit_length and slit_length > 0:
        _append_dql(out_path, float(slit_length))

    _append_manipulation_provenance(out_path, operation, provenance, source_path)
    return out_path


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

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
