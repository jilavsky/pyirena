"""
pyirena.batch.saxs_morph — Headless SAXS Morphology fitting (fit_saxs_morph).

Split from the original monolithic pyirena/batch.py (no behavior change).
"""

from __future__ import annotations

import logging
import traceback
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional, Union

import numpy as np

from pyirena.logging_setup import ensure_console_output as _ensure_console

if TYPE_CHECKING:
    pass

log = logging.getLogger(__name__)

from pyirena.batch._common import _load_config, _load_data


# ===========================================================================
# SAXS Morph (headless)
# ===========================================================================

def fit_saxs_morph(
    data_file: Union[str, Path],
    config_file: Union[str, Path],
    save_to_nexus: bool = True,
    with_uncertainty: bool = False,
    n_mc_runs: int = 10,
) -> Optional[Dict]:
    """Run a SAXS Morph (3D voxelgram) calculation driven by a pyIrena config file.

    The standard workflow has three sequential steps (matching what the GUI
    does when the user clicks "Fit Power-law Bckg" → "Fit Flat Bckg" →
    "Calculate 3D"):

      1. Power-law pre-fit over ``power_law_q_min/q_max``  → sets B, P.
      2. Flat-background pre-fit over ``background_q_min/q_max``  → sets bg.
      3. compute_voxelgram() over ``q_min/q_max`` with the resolved
         (φ, contrast) according to ``input_mode``.

    Each pre-fit step is skipped if its Q range is missing in the config
    (in which case the corresponding parameter value from the config is
    used as-is — useful when the user has already set good values and
    just wants to recompute the voxelgram).

    The ``with_uncertainty`` / ``n_mc_runs`` arguments are kept for API
    compatibility with the other fit_* functions but are no-ops here:
    SAXS Morph has no fittable parameters, so MC perturbation produces
    only the variation already captured by the RNG seed.

    Parameters
    ----------
    data_file       : NXcanSAS HDF5 (.h5/.hdf5/.nxs) or .dat/.txt SAS data.
    config_file     : pyIrena JSON config containing a ``'saxs_morph'`` block
                      (e.g. produced by Export Parameters in the GUI).
    save_to_nexus   : Save voxelgram + params into the HDF5 file (HDF5 only).
    with_uncertainty: Ignored (kept for API symmetry).
    n_mc_runs       : Ignored (kept for API symmetry).

    Returns
    -------
    dict or None    Same shape as :func:`fit_modeling`.
    """
    _ensure_console()
    from pyirena.core.saxs_morph import (
        SaxsMorphEngine, SaxsMorphConfig,
        fit_power_law_bg, fit_flat_bg,
    )
    from pyirena.io.nxcansas_saxs_morph import save_saxs_morph_results

    data_file = Path(data_file)
    config_file = Path(config_file)

    # --- Load config ---
    try:
        config = _load_config(config_file)
        if config is None:
            return None
        sm_cfg = config.get('saxs_morph')
        if sm_cfg is None:
            log.info(f"[pyirena.batch.fit_saxs_morph] No 'saxs_morph' section in '{config_file.name}'")
            return {'success': False, 'message': "No 'saxs_morph' section in config file"}
    except Exception:
        log.error(f"[pyirena.batch.fit_saxs_morph] Error reading config:\n{traceback.format_exc()}")
        return None

    # --- Load data ---
    try:
        data = _load_data(data_file)
        if data is None:
            return None
    except Exception:
        log.error(f"[pyirena.batch.fit_saxs_morph] Error loading data:\n{traceback.format_exc()}")
        return None

    q = np.asarray(data['Q'], dtype=float)
    I = np.asarray(data['Intensity'], dtype=float)
    _err = data.get('Error')
    dI = np.asarray(_err if _err is not None else 0.05 * np.abs(I), dtype=float)
    dI = np.where(dI <= 0, 0.05 * np.abs(I), dI)

    # --- Step 1: Power-law pre-fit (if Q range provided) ---
    pl_qmin = sm_cfg.get('power_law_q_min')
    pl_qmax = sm_cfg.get('power_law_q_max')
    pl_B = float(sm_cfg.get('power_law_B', 0.0))
    pl_P = float(sm_cfg.get('power_law_P', 4.0))
    if pl_qmin is not None and pl_qmax is not None and pl_qmax > pl_qmin:
        try:
            pl_B, pl_P = fit_power_law_bg(q, I, float(pl_qmin), float(pl_qmax))
            log.info(f"[pyirena.batch.fit_saxs_morph] Power-law pre-fit: "
                  f"B={pl_B:.4g}, P={pl_P:.4g}")
        except Exception:
            log.error(f"[pyirena.batch.fit_saxs_morph] Power-law pre-fit failed:\n"
                  f"{traceback.format_exc()}")

    # --- Step 2: Flat background pre-fit (if Q range provided) ---
    bg_qmin = sm_cfg.get('background_q_min')
    bg_qmax = sm_cfg.get('background_q_max')
    flat_bg = float(sm_cfg.get('background', 0.0))
    if bg_qmin is not None and bg_qmax is not None and bg_qmax > bg_qmin:
        try:
            flat_bg = fit_flat_bg(q, I, float(bg_qmin), float(bg_qmax),
                                  power_law_B=pl_B, power_law_P=pl_P)
            log.info(f"[pyirena.batch.fit_saxs_morph] Flat-bg pre-fit: "
                  f"background={flat_bg:.4g}")
        except Exception:
            log.error(f"[pyirena.batch.fit_saxs_morph] Flat-bg pre-fit failed:\n"
                  f"{traceback.format_exc()}")

    # --- Build SaxsMorphConfig with the (possibly updated) bg values ---
    try:
        cfg = SaxsMorphConfig(
            q_min=sm_cfg.get('q_min'),
            q_max=sm_cfg.get('q_max'),
            power_law_q_min=pl_qmin, power_law_q_max=pl_qmax,
            background_q_min=bg_qmin, background_q_max=bg_qmax,
            voxel_size_fit=int(sm_cfg.get('voxel_size_fit', 128)),
            voxel_size_render=int(sm_cfg.get('voxel_size_render', 256)),
            box_size_A=float(sm_cfg.get('box_size_A', 1000.0)),
            input_mode=str(sm_cfg.get('input_mode', 'phi')),
            volume_fraction=float(sm_cfg.get('volume_fraction', 0.30)),
            contrast=float(sm_cfg.get('contrast', 1.0)),
            power_law_B=pl_B, power_law_P=pl_P, background=flat_bg,
            rng_seed=sm_cfg.get('rng_seed'),
        )
    except Exception:
        log.error(f"[pyirena.batch.fit_saxs_morph] Error building SaxsMorphConfig:\n"
              f"{traceback.format_exc()}")
        return None

    # --- Step 3: Calculate 3D voxelgram at render resolution ---
    try:
        engine = SaxsMorphEngine()
        result = engine.compute_voxelgram(
            cfg, q, I, dI,
            voxel_size_override=cfg.voxel_size_render,
        )
    except Exception:
        log.error(f"[pyirena.batch.fit_saxs_morph] Calculate 3D failed:\n"
              f"{traceback.format_exc()}")
        return {'success': False, 'message': 'SAXS Morph calculation failed'}

    # --- Save ---
    out_path = None
    if save_to_nexus and data_file.suffix.lower() in ('.h5', '.hdf5', '.nxs'):
        try:
            save_saxs_morph_results(data_file, result)
            out_path = data_file
        except Exception:
            log.error(f"[pyirena.batch.fit_saxs_morph] Save error:\n"
                  f"{traceback.format_exc()}")

    chi2_str = f"χ²/dof={result.reduced_chi_squared:.4g}"
    msg = (f"SAXS Morph complete — {chi2_str}, "
           f"φ_actual={result.phi_actual:.3g}, "
           f"voxel={result.voxel_size}³")
    log.info(f"[pyirena.batch] {msg}")

    return {
        'success': True,
        'message': msg,
        'result': result,
        'output_file': out_path,
    }
