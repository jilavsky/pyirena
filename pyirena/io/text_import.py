"""
pyirena.io.text_import — clean and convert ASCII SAS files to NXcanSAS HDF5.

The canonical workflow for text data is:
1. Call ``ensure_nxcansas_sibling(text_path)`` once when a file is first used.
2. All subsequent operations (fitting, plotting, result saving) use the
   returned HDF5 path as if the data had always been in NXcanSAS format.

This confines all text-file awareness to the import step and removes the
ext-branch duplication that would otherwise appear in every consumer.

Cleaning rules (silent, logged, recorded in the HDF5)
------------------------------------------------------
- Points with Q ≤ 0 or non-finite Q are removed (occasional Q=0 rows).
- Points with I ≤ 0 or non-finite I are removed (beamstop zeros behind
  the beam-stop that the instrument should have stripped but often didn't;
  they are invisible on a log-intensity plot yet corrupt numerical fits).
- Surviving points with E ≤ 0 or non-finite E have their uncertainty
  replaced by  E = I × error_fraction  (default 5 %).
- dQ is truncated in parallel when present.

The cleaning report (counts per rule) is written into the HDF5 file under
``entry/notes/`` so provenance is never lost.

Naming convention
-----------------
The converted file is placed next to the original text file with the same
stem and a ``.h5`` extension: ``mydata.dat`` → ``mydata.h5`` (also applies
to ``.txt`` and ``.csv`` sources).

Collision guard: if ``mydata.h5`` already exists and was NOT created by
pyirena (no provenance marker), the fallback name ``mydata_NX.h5`` is used
and a one-line notice is printed.

Caching: if the sibling already exists, has the provenance marker, is newer
than the text file, and was converted with the same ``q_unit``, it is
returned as-is without re-reading the source. Pass ``force=True`` to always
regenerate.

Q units
-------
Text files carry no unit metadata, so pyIrena assumes Q is already in 1/Å
(its native unit) unless told otherwise via the ``q_unit`` parameter (one of
``pyirena.io.hdf5.Q_UNIT_TO_ANGSTROM``'s keys). The assumed unit is recorded
in the converted file's provenance (``source.q_unit_assumed`` metadata and
the ``pyirena_q_unit`` HDF5 attribute), and changing it for a given file
invalidates any previously cached sibling.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Marker attribute written to HDF5 root to identify files we created.
_PROVENANCE_ATTR = 'pyirena_converted_from'

# Q-unit assumed at conversion time, written alongside the provenance marker
# so a later change to the Q-unit setting invalidates a stale cached sibling
# even though its mtime still looks up to date.
_Q_UNIT_ATTR = 'pyirena_q_unit'


# ── Public API ────────────────────────────────────────────────────────────────

def clean_sas_arrays(
    Q: np.ndarray,
    I: np.ndarray,
    E: Optional[np.ndarray],
    dQ: Optional[np.ndarray] = None,
    error_fraction: float = 0.05,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray], dict]:
    """Remove bad points from SAS arrays and synthesize missing uncertainties.

    Parameters
    ----------
    Q : array-like   Q vector [Å⁻¹]
    I : array-like   Intensity [cm⁻¹]
    E : array-like or None   Uncertainty [cm⁻¹]; None → all synthesized
    dQ : array-like or None  Q resolution; truncated in parallel with Q/I/E
    error_fraction : float   Fraction of I used to synthesize missing E

    Returns
    -------
    Q_clean, I_clean, E_clean, dQ_clean : cleaned arrays (E always present)
    report : dict with counts removed/synthesized, for logging and provenance
    """
    Q = np.asarray(Q, dtype=float)
    I = np.asarray(I, dtype=float)
    E_in = np.asarray(E, dtype=float) if E is not None else None
    dQ_in = np.asarray(dQ, dtype=float) if dQ is not None else None

    n_original = len(Q)

    # ── Boolean mask: points to keep ──────────────────────────────────────
    keep = np.isfinite(Q) & (Q > 0) & np.isfinite(I) & (I > 0)

    n_q_removed = int(np.sum(~np.isfinite(Q) | (Q <= 0)))
    n_i_removed = int(np.sum(np.isfinite(Q) & (Q > 0) & (~np.isfinite(I) | (I <= 0))))

    Q_c = Q[keep]
    I_c = I[keep]
    dQ_c = dQ_in[keep] if dQ_in is not None and len(dQ_in) == len(Q) else None

    # ── Uncertainties ─────────────────────────────────────────────────────
    if E_in is not None and len(E_in) == len(Q):
        E_c = E_in[keep]
    else:
        E_c = np.zeros(len(Q_c), dtype=float)

    # Replace bad (≤0 or non-finite) uncertainties with fractional estimate
    bad_e = ~np.isfinite(E_c) | (E_c <= 0)
    n_e_synthesized = int(np.sum(bad_e))
    E_c[bad_e] = I_c[bad_e] * error_fraction

    report = {
        'n_original':      n_original,
        'n_kept':          len(Q_c),
        'n_q_removed':     n_q_removed,
        'n_i_removed':     n_i_removed,
        'n_e_synthesized': n_e_synthesized,
        'error_fraction':  error_fraction,
    }

    if n_q_removed or n_i_removed or n_e_synthesized:
        logger.info(
            "clean_sas_arrays: %d/%d points kept; "
            "removed %d Q≤0, %d I≤0; synthesized %d uncertainties (×%.4f)",
            len(Q_c), n_original, n_q_removed, n_i_removed,
            n_e_synthesized, error_fraction,
        )

    return Q_c, I_c, E_c, dQ_c, report


def converted_sibling_path(text_path: Path) -> Path:
    """Return the preferred HDF5 sibling path for a text SAS file.

    ``mydata.dat`` → ``mydata.h5`` (same directory).

    Does NOT check for existence or collisions — use
    ``ensure_nxcansas_sibling`` for the full conversion workflow.
    """
    return text_path.with_suffix('.h5')


def ensure_nxcansas_sibling(
    text_path: Path,
    error_fraction: float = 0.05,
    q_unit: str = '1/A',
    force: bool = False,
) -> Path:
    """Return an up-to-date NXcanSAS HDF5 sibling for a text SAS file.

    Converts, cleans, and caches.  All downstream tools should call this and
    then treat the returned path as a normal NXcanSAS file.

    Parameters
    ----------
    text_path : Path   Path to the source ``.dat`` / ``.txt`` / ``.csv`` file.
    error_fraction : float
        Fraction of I used to synthesize missing uncertainties (default 0.05).
    q_unit : str
        Unit of the Q (and dQ) column in the source file — one of
        ``pyirena.io.hdf5.Q_UNIT_TO_ANGSTROM``'s keys (default ``'1/A'``,
        i.e. already in pyIrena's native unit, no conversion).  A non-default
        unit forces reconversion of a previously cached sibling that was
        made with a different unit — see ``_is_valid_sibling``.
    force : bool
        Re-convert even if an up-to-date sibling already exists.

    Returns
    -------
    Path   Path to the cleaned NXcanSAS HDF5 file.

    Raises
    ------
    RuntimeError   If the text file cannot be read or converted.
    """
    from pyirena.io.hdf5 import readTextFile
    from pyirena.io.nxcansas_unified import create_nxcansas_file

    text_path = Path(text_path)
    if not text_path.exists():
        raise RuntimeError(f"Text file not found: '{text_path}'")

    # ── Determine output path (with collision guard) ───────────────────────
    preferred = converted_sibling_path(text_path)  # <stem>.h5
    out_path = _resolve_output_path(text_path, preferred)

    # ── Cache check ───────────────────────────────────────────────────────
    if not force and _is_valid_sibling(out_path, text_path, q_unit):
        logger.debug("ensure_nxcansas_sibling: reusing cached '%s'", out_path.name)
        return out_path

    # ── Read raw text ─────────────────────────────────────────────────────
    raw = readTextFile(str(text_path.parent), text_path.name,
                       error_fraction=error_fraction, q_unit=q_unit)
    if raw is None:
        raise RuntimeError(f"Could not read text file: '{text_path}'")

    Q_raw = np.asarray(raw['Q'],         dtype=float)
    I_raw = np.asarray(raw['Intensity'], dtype=float)
    E_raw = raw.get('Error')
    dQ_raw = raw.get('dQ')

    # ── Clean ─────────────────────────────────────────────────────────────
    Q, I, E, dQ, report = clean_sas_arrays(
        Q_raw, I_raw, E_raw, dQ_raw, error_fraction=error_fraction
    )

    if len(Q) == 0:
        raise RuntimeError(
            f"No valid data points remain after cleaning '{text_path.name}'. "
            f"Check that the file contains positive Q and positive I values."
        )

    # ── Build provenance metadata ──────────────────────────────────────────
    cleaning_lines = [
        f"Converted from: {text_path.name}",
        f"Q unit assumed: {q_unit}" + (
            "" if q_unit == '1/A' else " (converted to 1/Å)"
        ),
        f"Original points: {report['n_original']}",
        f"Kept: {report['n_kept']}",
        f"Removed (Q<=0 or non-finite Q): {report['n_q_removed']}",
        f"Removed (I<=0 or non-finite I): {report['n_i_removed']}",
        f"Uncertainties synthesized (E=I*{error_fraction:.4f}): {report['n_e_synthesized']}",
    ]
    metadata = {
        'source.converted_from':    text_path.name,
        'source.conversion_notes':  '; '.join(cleaning_lines),
        'source.error_fraction':    str(error_fraction),
        'source.q_unit_assumed':    q_unit,
    }

    # ── Write HDF5 ────────────────────────────────────────────────────────
    create_nxcansas_file(
        filepath=out_path,
        q=Q,
        intensity=I,
        error=E,
        dq=dQ,
        sample_name=text_path.stem,
        metadata=metadata,
    )

    # Write the provenance markers so the collision guard can identify our
    # files, and so a later q_unit change can invalidate this cached sibling.
    import h5py
    with h5py.File(out_path, 'a') as f:
        f.attrs[_PROVENANCE_ATTR] = text_path.name
        f.attrs[_Q_UNIT_ATTR] = q_unit

    logger.info("ensure_nxcansas_sibling: wrote '%s' (%d pts)", out_path.name, len(Q))
    return out_path


# ── Internal helpers ──────────────────────────────────────────────────────────

def _is_valid_sibling(sibling: Path, source: Path, q_unit: str = '1/A') -> bool:
    """True if sibling exists, has our marker, is newer than source, and was
    converted with the same ``q_unit`` requested now (older siblings without
    a recorded unit are treated as ``'1/A'``, matching pre-feature behavior).
    """
    if not sibling.exists():
        return False
    try:
        import h5py
        with h5py.File(sibling, 'r') as f:
            if _PROVENANCE_ATTR not in f.attrs:
                return False
            if f.attrs.get(_Q_UNIT_ATTR, '1/A') != q_unit:
                return False
        src_mtime = source.stat().st_mtime
        sib_mtime = sibling.stat().st_mtime
        return sib_mtime >= src_mtime
    except Exception:
        return False


def _resolve_output_path(text_path: Path, preferred: Path) -> Path:
    """Return preferred path unless it's an unrelated HDF5 (collision guard).

    If ``preferred`` exists without our provenance marker, fall back to
    ``<stem>_NX.h5`` and print a one-line notice.
    """
    if not preferred.exists():
        return preferred

    # Check whether it's ours
    try:
        import h5py
        with h5py.File(preferred, 'r') as f:
            if _PROVENANCE_ATTR in f.attrs:
                return preferred   # Our file — safe to overwrite/reuse
    except Exception:
        # Can't read it → treat as foreign
        logger.debug("Could not inspect %s; treating as foreign file", preferred,
                     exc_info=True)

    # Foreign file: use fallback name
    fallback = text_path.with_name(text_path.stem + '_NX.h5')
    logger.info(
        f"[pyirena] Note: '{preferred.name}' already exists and was not created by "
        f"pyirena. Using '{fallback.name}' for the converted copy."
    )
    return fallback
