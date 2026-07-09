"""
pyirena.batch.merge — Headless data merging (merge_data).

Split from the original monolithic pyirena/batch.py (no behavior change).
"""

from __future__ import annotations

import json
import logging
import traceback
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional, Union

import numpy as np

from pyirena.logging_setup import ensure_console_output as _ensure_console

if TYPE_CHECKING:
    pass

log = logging.getLogger(__name__)

from pyirena.batch._common import _load_data


# ---------------------------------------------------------------------------
# Data Merge
# ---------------------------------------------------------------------------

def merge_data(
    file1: Union[str, Path],
    file2: Union[str, Path],
    config_file: Optional[Union[str, Path]] = None,
    save_to_nexus: bool = True,
    output_folder: Optional[Union[str, Path]] = None,
    verbose: bool = True,
) -> Optional[Dict]:
    """Headless merge: load two SAS files, optimise merge parameters, save result.

    Parameters
    ----------
    file1 : str or Path
        DS1 file (lower Q, assumed absolute intensity scale).  Accepts .h5,
        .hdf5, .hdf (NXcanSAS or generic HDF5) or .dat/.txt (text).
    file2 : str or Path
        DS2 file (higher Q).  Same supported formats as file1.
    config_file : str or Path or None
        JSON file with merge parameters.  Expected keys (all optional):
        ``q_overlap_min``, ``q_overlap_max``, ``scale_dataset`` (1 or 2),
        ``fit_scale`` (bool), ``qshift_dataset`` (0/1/2), ``fit_qshift``
        (bool), ``split_at_left_cursor`` (bool).  Missing keys fall back to
        defaults (scale DS2, no Q-shift, include overlap in output).
    save_to_nexus : bool
        If True, write merged data to a NXcanSAS file.
    output_folder : str or Path or None
        Directory for the output file.  If None, a folder named
        ``{folder1_name}_merged`` is created next to DS1's parent folder.
    verbose : bool
        Print progress to stdout.

    Returns
    -------
    dict or None
        ``{'q', 'I', 'dI', 'dQ', 'scale', 'q_shift', 'background',
        'chi_squared', 'success', 'output_path'}``  or None on failure.
    """
    _ensure_console()
    from pyirena.core.data_merge import DataMerge, MergeConfig

    file1 = Path(file1)
    file2 = Path(file2)

    # ── Load data ────────────────────────────────────────────────────────────
    data1 = _load_data(file1)
    data2 = _load_data(file2)
    if data1 is None or data2 is None:
        return None

    q1 = data1['Q'];  I1 = data1['Intensity']
    _e1 = data1.get('Error');  dI1 = _e1 if _e1 is not None else I1 * 0.05
    q2 = data2['Q'];  I2 = data2['Intensity']
    _e2 = data2.get('Error');  dI2 = _e2 if _e2 is not None else I2 * 0.05
    dQ1 = data1.get('dQ')
    dQ2 = data2.get('dQ')

    # ── Validate data quality ─────────────────────────────────────────────────
    for label, q, I in (('DS1', q1, I1), ('DS2', q2, I2)):
        valid = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
        if not np.any(valid):
            log.error(f"[merge_data] {label} has no valid (finite, positive) data points — skipping.")
            return None
        if int(np.sum(valid)) < 3:
            log.error(f"[merge_data] {label} has fewer than 3 valid data points — skipping.")
            return None
    if float(q1.max()) <= float(q2.min()):
        log.error(f"[merge_data] No Q overlap: DS1 ends at {q1.max():.4g}, "
              f"DS2 starts at {q2.min():.4g} — skipping.")
        return None

    # ── Build MergeConfig ─────────────────────────────────────────────────────
    cfg_dict: Dict = {}
    if config_file is not None:
        config_file = Path(config_file)
        if not config_file.exists():
            if verbose:
                log.error(f"[merge_data] Config file not found: {config_file}")
        else:
            try:
                with open(config_file) as fh:
                    cfg_dict = json.load(fh)
            except Exception as exc:
                if verbose:
                    log.error(f"[merge_data] Cannot read config file: {exc}")

    # Determine Q overlap range
    q_overlap_min = cfg_dict.get('q_overlap_min')
    q_overlap_max = cfg_dict.get('q_overlap_max')
    if q_overlap_min is None or q_overlap_max is None:
        # Auto-detect: 10th–90th percentile of the Q intersection
        q_max1 = float(q1.max());  q_min2 = float(q2.min())
        if q_max1 <= q_min2:
            if verbose:
                log.info(f"[merge_data] No Q overlap detected between DS1 "
                      f"(max {q_max1:.4g}) and DS2 (min {q_min2:.4g}).")
            return None
        q_ov_lo = max(float(q1.min()), float(q2.min()))
        q_ov_hi = min(float(q1.max()), float(q2.max()))
        q_overlap_min = q_ov_lo + 0.1 * (q_ov_hi - q_ov_lo)
        q_overlap_max = q_ov_hi - 0.1 * (q_ov_hi - q_ov_lo)
        if verbose:
            log.info(f"[merge_data] Auto overlap range: "
                  f"[{q_overlap_min:.4g}, {q_overlap_max:.4g}] Å⁻¹")

    config = MergeConfig(
        q_overlap_min=float(q_overlap_min),
        q_overlap_max=float(q_overlap_max),
        fit_scale=bool(cfg_dict.get('fit_scale', True)),
        scale_dataset=int(cfg_dict.get('scale_dataset', 2)),
        fixed_scale_value=float(cfg_dict.get('fixed_scale_value', 1.0)),
        fit_qshift=bool(cfg_dict.get('fit_qshift', False)),
        fixed_qshift_value=float(cfg_dict.get('fixed_qshift_value', 0.0)),
        qshift_dataset=int(cfg_dict.get('qshift_dataset', 0)),
        split_at_left_cursor=bool(cfg_dict.get('split_at_left_cursor', False)),
    )

    # ── Optimise ──────────────────────────────────────────────────────────────
    engine = DataMerge()
    if verbose:
        log.info("[merge_data] Optimising merge …")
    try:
        result = engine.optimize(q1, I1, dI1, q2, I2, dI2, config)
    except Exception:
        log.error(f"[merge_data] Merge optimisation failed for "
              f"'{file1.name}' + '{file2.name}':\n{traceback.format_exc()}")
        return None
    if verbose:
        status = "OK" if result.success else "FAILED"
        log.info(f"[merge_data] {status}  scale={result.scale:.4g}  "
              f"bg={result.background:.4g}  χ²={result.chi_squared:.4g}  "
              f"n_pts={result.n_overlap_points}  msg={result.message}")

    # ── Merge arrays ──────────────────────────────────────────────────────────
    try:
        q_m, I_m, dI_m, dQ_m = engine.merge(
            q1, I1, dI1, dQ1, q2, I2, dI2, dQ2, result, config
        )
    except Exception:
        log.error(f"[merge_data] Merge array assembly failed for "
              f"'{file1.name}' + '{file2.name}':\n{traceback.format_exc()}")
        return None

    # ── Save ──────────────────────────────────────────────────────────────────
    out_path: Optional[Path] = None
    if save_to_nexus:
        from pyirena.io.nxcansas_data_merge import save_merged_data

        if output_folder is None:
            parent = file1.parent.parent if file1.parent.parent != file1.parent else file1.parent
            output_folder = parent / f"{file1.parent.name}_merged"

        merge_result_dict = {
            'scale': result.scale,
            'q_shift': result.q_shift,
            'background': result.background,
            'chi_squared': result.chi_squared,
            'n_overlap_points': result.n_overlap_points,
            'q_overlap_min': config.q_overlap_min,
            'q_overlap_max': config.q_overlap_max,
            'scale_dataset': config.scale_dataset,
            'fit_scale': config.fit_scale,
            'qshift_dataset': config.qshift_dataset,
            'fit_qshift': config.fit_qshift,
            'split_at_left_cursor': config.split_at_left_cursor,
        }
        try:
            out_path = save_merged_data(
                output_folder=Path(output_folder),
                ds1_path=file1,
                ds1_is_nxcansas=data1.get('is_nxcansas', False),
                q=q_m, I=I_m, dI=dI_m, dQ=dQ_m,
                merge_result_dict=merge_result_dict,
                ds2_path=file2,
            )
            if verbose:
                log.info(f"[merge_data] Saved → {out_path}")
        except Exception as exc:
            if verbose:
                log.error(f"[merge_data] Save error: {exc}")

    return {
        'q': q_m, 'I': I_m, 'dI': dI_m, 'dQ': dQ_m,
        'scale': result.scale,
        'q_shift': result.q_shift,
        'background': result.background,
        'chi_squared': result.chi_squared,
        'success': result.success,
        'output_path': str(out_path) if out_path else None,
    }
