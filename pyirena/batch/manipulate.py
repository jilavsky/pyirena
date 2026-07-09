"""
pyirena.batch.manipulate — Headless data manipulation (manipulate_data, average_data).

Split from the original monolithic pyirena/batch.py (no behavior change).
"""

from __future__ import annotations

import logging
import traceback
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional


from pyirena.logging_setup import ensure_console_output as _ensure_console

if TYPE_CHECKING:
    pass

log = logging.getLogger(__name__)

from pyirena.batch._common import _load_data


# ===========================================================================
# Data Manipulation (headless)
# ===========================================================================

def manipulate_data(
    data_file,
    operation: str,
    config: Optional[Dict] = None,
    buffer_file=None,
    output_folder=None,
    verbose: bool = True,
) -> Optional[Dict]:
    """Apply a data manipulation operation to a SAS dataset.

    Parameters
    ----------
    data_file : str or Path
        Path to the input data file.
    operation : str
        One of 'scale', 'trim', 'rebin', 'subtract', 'divide'.
    config : dict, optional
        Operation-specific parameters.  If None, uses defaults.
    buffer_file : str or Path, optional
        Path to the buffer/denominator file (for subtract/divide).
    output_folder : str or Path, optional
        Where to save the result.  If None, uses ``{data_folder}_manip``.
    verbose : bool
        Print progress messages.

    Returns
    -------
    dict or None
        Result dict with keys: success, operation, output_file, message.
    """
    _ensure_console()
    from pyirena.core.data_manipulation import (
        DataManipulation, ScaleConfig, TrimConfig, RebinConfig,
        SubtractConfig, DivideConfig,
    )
    from pyirena.io.nxcansas_data_manipulation import save_manipulated_data

    data_file = Path(data_file)
    data = _load_data(data_file)
    if data is None:
        return None

    q, I, dI = data['Q'], data['Intensity'], data.get('Error', data['Intensity'] * 0.05)
    dQ = data.get('dQ')
    cfg = config or {}
    engine = DataManipulation

    try:
        if operation == 'scale':
            result = engine.scale(q, I, dI, dQ, ScaleConfig(**cfg))
        elif operation == 'trim':
            result = engine.trim(q, I, dI, dQ, TrimConfig(**cfg))
        elif operation == 'rebin':
            result = engine.rebin(q, I, dI, dQ, RebinConfig(**cfg))
        elif operation == 'subtract':
            if buffer_file is None:
                log.info("[pyirena.batch] buffer_file required for subtract")
                return None
            buf = _load_data(Path(buffer_file))
            if buf is None:
                return None
            result = engine.subtract(
                q, I, dI, dQ,
                buf['Q'], buf['Intensity'],
                buf.get('Error', buf['Intensity'] * 0.05),
                SubtractConfig(**cfg),
            )
        elif operation == 'divide':
            if buffer_file is None:
                log.info("[pyirena.batch] buffer_file (denominator) required for divide")
                return None
            den = _load_data(Path(buffer_file))
            if den is None:
                return None
            result = engine.divide(
                q, I, dI, dQ,
                den['Q'], den['Intensity'],
                den.get('Error', den['Intensity'] * 0.05),
                DivideConfig(**cfg),
            )
        else:
            log.info(f"[pyirena.batch] Unknown operation: {operation}")
            return None
    except Exception:
        log.error(f"[pyirena.batch.manipulate_data] Error:\n{traceback.format_exc()}")
        return None

    # Save
    if output_folder is None:
        output_folder = str(data_file.parent) + '_manip'
    output_folder = Path(output_folder)

    try:
        out_path = save_manipulated_data(
            output_folder=output_folder,
            source_path=data_file,
            source_is_nxcansas=data.get('is_nxcansas', False),
            q=result.q, I=result.I, dI=result.dI, dQ=result.dQ,
            operation=result.operation,
            provenance=result.metadata,
        )
    except Exception:
        log.error(f"[pyirena.batch.manipulate_data] Save error:\n{traceback.format_exc()}")
        return None

    msg = f"{operation} complete: {len(result.q)} points → {out_path.name}"
    if verbose:
        log.info(f"[pyirena.batch] {msg}")

    return {
        'success': True,
        'operation': operation,
        'output_file': out_path,
        'message': msg,
    }


def average_data(
    data_files,
    output_folder=None,
    verbose: bool = True,
    similarity_check: bool = False,
    similarity_p_min: float = 0.01,
    similarity_method: str = 'cormap',
    similarity_reference: str = 'first',
    similarity_normalize_scale: bool = True,
) -> Optional[Dict]:
    """Average multiple SAS datasets.

    Parameters
    ----------
    data_files : list of str or Path
        Paths to data files to average.
    output_folder : str or Path, optional
        Where to save.  If None, uses ``{first_file_folder}_manip``.
    verbose : bool
        Print progress.
    similarity_check : bool
        When True, run a similarity analysis before averaging and discard
        outliers (likely radiation-damaged frames).
    similarity_p_min : float
        P-value threshold for rejection (0–1).  Frames with p < threshold are
        discarded.  Typical range: 0.001–0.05.
    similarity_method : str
        Similarity algorithm.  Currently ``'cormap'`` (Franke 2015).
    similarity_reference : str
        ``'first'`` — compare each frame vs. frame 0 (frame 0 always kept);
        ``'majority'`` — compare each frame vs. the median of all frames.
    similarity_normalize_scale : bool
        When True (default), rescale each frame to match the reference before
        comparing.  Removes flux-drift / absorption differences so that only
        shape differences are detected.

    Returns
    -------
    dict or None
        On success includes ``'rejected'``: list of ``(filename, p_value)``
        pairs for frames that were discarded by the similarity filter.
    """
    _ensure_console()
    from pyirena.core.data_manipulation import DataManipulation
    from pyirena.io.nxcansas_data_manipulation import save_manipulated_data

    files = [Path(f) for f in data_files]
    if not files:
        return None

    datasets = []
    loaded_files = []   # Path objects parallel to datasets
    loaded_data = []    # raw dicts parallel to datasets
    for fp in files:
        d = _load_data(fp)
        if d is None:
            continue
        q, I = d['Q'], d['Intensity']
        dI = d.get('Error', I * 0.05)
        dQ = d.get('dQ')
        datasets.append((q, I, dI, dQ))
        loaded_files.append(fp)
        loaded_data.append(d)

    if len(datasets) < 2:
        log.info("[pyirena.batch] Need at least 2 datasets to average")
        return None

    # --- Optional similarity filter ---
    rejected_pairs: list[tuple[str, float]] = []
    if similarity_check:
        from pyirena.core.similarity import check_similarity
        filenames = [fp.name for fp in loaded_files]
        sim_results = check_similarity(
            datasets,
            filenames=filenames,
            method=similarity_method,
            reference=similarity_reference,
            p_min=similarity_p_min,
            normalize_scale=similarity_normalize_scale,
        )
        if verbose:
            for r in sim_results:
                tag = "OK    " if r.accepted else "REJECT"
                log.info(f"  [{tag}] {r.filename}  p={r.p_value:.4f}"
                      f"  C={r.longest_run}  N={r.n_points}")
        accepted_idx = [r.idx for r in sim_results if r.accepted]
        rejected_pairs = [(r.filename, r.p_value) for r in sim_results if not r.accepted]
        datasets = [datasets[i] for i in accepted_idx]
        loaded_files = [loaded_files[i] for i in accepted_idx]
        loaded_data = [loaded_data[i] for i in accepted_idx]
        if len(datasets) < 2:
            log.info("[pyirena.batch] Too few datasets remain after similarity filtering.")
            return None

    result = DataManipulation.average(datasets, reference_index=0)

    if output_folder is None:
        output_folder = str(loaded_files[0].parent) + '_manip'
    output_folder = Path(output_folder)

    try:
        out_path = save_manipulated_data(
            output_folder=output_folder,
            source_path=loaded_files[0],
            source_is_nxcansas=loaded_data[0].get('is_nxcansas', False),
            q=result.q, I=result.I, dI=result.dI, dQ=result.dQ,
            operation=result.operation,
            provenance=result.metadata,
        )
    except Exception:
        log.error(f"[pyirena.batch.average_data] Save error:\n{traceback.format_exc()}")
        return None

    msg = f"Average of {len(datasets)} datasets → {out_path.name}"
    if verbose:
        log.info(f"[pyirena.batch] {msg}")

    return {
        'success': True,
        'operation': 'avg',
        'output_file': out_path,
        'message': msg,
        'n_datasets': len(datasets),
        'rejected': rejected_pairs,
    }
