"""
pyirena.batch.unified — Headless Unified Fit (fit_unified) and its helpers.

Split from the original monolithic pyirena/batch.py (no behavior change).
"""

from __future__ import annotations

import logging
import traceback
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional, Union

import numpy as np

from pyirena.core.unified import UnifiedFitModel, UnifiedLevel
from pyirena.logging_setup import ensure_console_output as _ensure_console

if TYPE_CHECKING:
    pass

log = logging.getLogger(__name__)

from pyirena.batch._common import _load_config, _load_data


def _state_to_model(unified_state: Dict) -> UnifiedFitModel:
    """Convert the unified_fit config section into a configured UnifiedFitModel."""
    num_levels = unified_state.get('num_levels', 1)
    no_limits = unified_state.get('no_limits', False)
    bg = unified_state.get('background', {})

    model = UnifiedFitModel(num_levels=num_levels)
    model.background = float(bg.get('value', 0.0))
    model.fit_background = bool(bg.get('fit', False))

    for i, ls in enumerate(unified_state.get('levels', [])[:num_levels]):
        def _v(param, key='value', default=0.0):
            entry = ls.get(param, {})
            if isinstance(entry, dict):
                return float(entry.get(key, default))
            return float(entry)

        def _b(param, key, default=False):
            entry = ls.get(param, {})
            if isinstance(entry, dict):
                return bool(entry.get(key, default))
            return default

        if no_limits:
            level = UnifiedLevel(
                G=_v('G'), Rg=_v('Rg'), B=_v('B'), P=_v('P'),
                RgCO=float(ls.get('RgCutoff', 0.0)),
                ETA=_v('ETA'), PACK=_v('PACK'),
                correlations=bool(ls.get('correlated', False)),
                # Propagate parameter links so B (and RgCO) are recomputed at
                # every fit iteration from the live G/Rg/P, mirroring the GUI's
                # run_fit().  estimate_B == "Calculate B" maps to link_B.
                link_B=bool(ls.get('estimate_B', False)),
                link_RGCO=bool(ls.get('link_rgco', False)),
                fit_G=_b('G', 'fit'), fit_Rg=_b('Rg', 'fit'),
                fit_B=_b('B', 'fit'), fit_P=_b('P', 'fit'),
                fit_ETA=_b('ETA', 'fit'), fit_PACK=_b('PACK', 'fit'),
                # Wide default bounds
                G_limits=(1e-10, 1e10), Rg_limits=(0.1, 1e6),
                B_limits=(1e-20, 1e10), P_limits=(0.0, 6.0),
                ETA_limits=(0.1, 1e6), PACK_limits=(0.0, 16.0),
            )
        else:
            level = UnifiedLevel(
                G=_v('G'), Rg=_v('Rg'), B=_v('B'), P=_v('P'),
                RgCO=float(ls.get('RgCutoff', 0.0)),
                ETA=_v('ETA'), PACK=_v('PACK'),
                correlations=bool(ls.get('correlated', False)),
                # Propagate parameter links so B (and RgCO) are recomputed at
                # every fit iteration from the live G/Rg/P, mirroring the GUI's
                # run_fit().  estimate_B == "Calculate B" maps to link_B.
                link_B=bool(ls.get('estimate_B', False)),
                link_RGCO=bool(ls.get('link_rgco', False)),
                fit_G=_b('G', 'fit'), fit_Rg=_b('Rg', 'fit'),
                fit_B=_b('B', 'fit'), fit_P=_b('P', 'fit'),
                fit_ETA=_b('ETA', 'fit'), fit_PACK=_b('PACK', 'fit'),
                G_limits=(_v('G', 'low_limit'), _v('G', 'high_limit', 1e10)),
                Rg_limits=(_v('Rg', 'low_limit', 0.1), _v('Rg', 'high_limit', 1e6)),
                B_limits=(_v('B', 'low_limit'), _v('B', 'high_limit', 1e10)),
                P_limits=(_v('P', 'low_limit'), _v('P', 'high_limit', 6.0)),
                ETA_limits=(_v('ETA', 'low_limit', 0.1), _v('ETA', 'high_limit', 1e6)),
                PACK_limits=(_v('PACK', 'low_limit'), _v('PACK', 'high_limit', 16.0)),
            )

        model.levels[i] = level

    return model


def _compute_invariant_sv(G, Rg, B, P, RgCO, ETA, PACK, correlated):
    """
    Compute the scattering invariant and Sv for one Unified Fit level.

    Thin wrapper around :func:`pyirena.core.unified.compute_invariant_sv`,
    the single (Igor-faithful) implementation shared with the GUI, so the
    batch path produces the same values as the GUI.

    Returns (invariant_cm4, sv) where either value is None when not computable.
    Sv is only meaningful when P is in the Porod regime (3.95 – 4.05).
    """
    from pyirena.core.unified import compute_invariant_sv
    return compute_invariant_sv(G, Rg, B, P, RgCO, ETA, PACK, correlated)


def _mc_uncertainty_unified(
    model: UnifiedFitModel,
    q_fit: np.ndarray,
    intensity_fit: np.ndarray,
    error_fit: Optional[np.ndarray],
    n_runs: int = 10,
) -> Optional[Dict]:
    """Run MC noise-perturbation fits and return per-parameter std-devs.

    Each of *n_runs* iterations perturbs intensity by Gaussian noise scaled by
    the measurement errors, refits the model, and collects the fitted values.
    Returns an ``uncertainties`` dict (format expected by
    ``save_unified_fit_results``) or ``None`` when fewer than 2 runs succeed.
    """
    import copy

    num_levels = model.num_levels
    err = error_fit if error_fit is not None else intensity_fit * 0.05
    err = np.maximum(err, 1e-30)

    mc: Dict[str, list] = {
        'background': [],
        'G':    [[] for _ in range(num_levels)],
        'Rg':   [[] for _ in range(num_levels)],
        'B':    [[] for _ in range(num_levels)],
        'P':    [[] for _ in range(num_levels)],
        'ETA':  [[] for _ in range(num_levels)],
        'PACK': [[] for _ in range(num_levels)],
    }

    for _ in range(n_runs):
        I_noisy = intensity_fit + np.random.normal(0.0, 1.0, len(q_fit)) * err
        I_noisy = np.maximum(I_noisy, 1e-30)
        mc_model = copy.deepcopy(model)
        try:
            res = mc_model.fit(q_fit, I_noisy, error_fit)
        except Exception:
            continue
        if not res.get('success', False):
            continue
        mc['background'].append(res.get('background', 0.0))
        for i, lv in enumerate(res.get('levels', mc_model.levels)):
            mc['G'][i].append(lv.G)
            mc['Rg'][i].append(lv.Rg)
            mc['B'][i].append(lv.B)
            mc['P'][i].append(lv.P)
            mc['ETA'][i].append(lv.ETA)
            mc['PACK'][i].append(lv.PACK)

    n_ok = len(mc['background'])
    if n_ok < 2:
        return None

    ddof = 1
    uncertainties: Dict = {
        'background': float(np.std(mc['background'], ddof=ddof)),
        'levels': [],
    }
    for i in range(num_levels):
        ud: Dict[str, float] = {}
        for key in ('G', 'Rg', 'B', 'P', 'ETA', 'PACK'):
            vals = mc[key][i]
            ud[key] = float(np.std(vals, ddof=ddof)) if len(vals) > 1 else 0.0
        uncertainties['levels'].append(ud)

    return uncertainties


def _build_setup_state(unified_state: Dict, fit_result: Dict) -> Dict:
    """Build a GUI-compatible setup_state from the config + fitted parameter values.

    Starts from the original config state (which carries fit flags, bounds, links,
    cursor positions, etc.) and overwrites the parameter values with the fitted ones
    so that "Load Setup from File" in the GUI restores the fitted — not initial — values.
    """
    import copy
    state = copy.deepcopy(unified_state)

    for i, lv in enumerate(fit_result.get('levels', [])):
        if i >= len(state.get('levels', [])):
            break
        ls = state['levels'][i]
        for param, fitted_val in [('G', lv.G), ('Rg', lv.Rg), ('B', lv.B),
                                   ('P', lv.P), ('ETA', lv.ETA), ('PACK', lv.PACK)]:
            entry = ls.get(param)
            if isinstance(entry, dict):
                entry['value'] = float(fitted_val)
            else:
                ls[param] = float(fitted_val)
        # RgCutoff is stored as a plain float in the state dict
        ls['RgCutoff'] = float(lv.RgCO)

    bg_val = fit_result.get('background', 0.0)
    bg = state.get('background', {})
    if isinstance(bg, dict):
        bg['value'] = float(bg_val)
    else:
        state['background'] = {'value': float(bg_val), 'fit': False}

    return state


def _save_to_nexus(data: Dict, model: UnifiedFitModel,
                   fit_result: Dict, num_levels: int,
                   uncertainties: Optional[Dict] = None,
                   setup_state: Optional[Dict] = None) -> Optional[Path]:
    """Save fit results to NXcanSAS HDF5 file.  Returns output path, or None on error."""
    from pyirena.io.nxcansas_unified import (
        save_unified_fit_results, get_output_filepath, create_nxcansas_file
    )

    source_path = Path(data['filepath'])
    is_nxcansas = data.get('is_nxcansas', False)
    output_path = get_output_filepath(source_path, is_nxcansas)

    # Build level dicts expected by save_unified_fit_results
    levels = []
    for lv in fit_result['levels']:
        level_dict = {
            'G': lv.G, 'Rg': lv.Rg, 'B': lv.B, 'P': lv.P,
            'RgCutoff': lv.RgCO, 'ETA': lv.ETA, 'PACK': lv.PACK,
            'correlated': lv.correlations,
        }
        invariant, sv = _compute_invariant_sv(
            lv.G, lv.Rg, lv.B, lv.P, lv.RgCO, lv.ETA, lv.PACK, lv.correlations
        )
        if invariant is not None:
            level_dict['Invariant'] = invariant
        if sv is not None:
            level_dict['Sv'] = sv
        levels.append(level_dict)

    intensity_model = model.calculate_intensity(data['Q'])
    if data.get('Error') is not None:
        residuals = (data['Intensity'] - intensity_model) / data['Error']
    else:
        residuals = (data['Intensity'] - intensity_model) / data['Intensity']

    # Create base NXcanSAS file first if needed
    if not is_nxcansas and not output_path.exists():
        label = data.get('label', 'data')
        create_nxcansas_file(
            output_path, data['Q'], data['Intensity'],
            data.get('Error'), sample_name=Path(label).stem
        )

    save_unified_fit_results(
        filepath=output_path,
        q=data['Q'],
        intensity_data=data['Intensity'],
        intensity_model=intensity_model,
        residuals=residuals,
        levels=levels,
        background=fit_result['background'],
        chi_squared=fit_result.get('chi_squared', 0.0),
        num_levels=num_levels,
        error=data.get('Error'),
        uncertainties=uncertainties,
        setup_state=setup_state,
    )

    return output_path


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fit_unified(
    data_file: Union[str, Path],
    config_file: Union[str, Path],
    save_to_nexus: bool = True,
    with_uncertainty: bool = False,
    n_mc_runs: int = 10,
) -> Optional[Dict]:
    """Fit Unified Fit model to a data file using parameters from a pyIrena config file.

    Parameters
    ----------
    data_file : str or Path
        Path to SAS data file (.h5/.hdf5 NXcanSAS, or .dat/.txt text).
    config_file : str or Path
        Path to a pyIrena JSON configuration file (created by Export Parameters
        in the GUI, or written manually).
    save_to_nexus : bool, optional
        If True (default), save fit results to an NXcanSAS HDF5 file alongside
        the input data (same file for NXcanSAS input, or ``<name>_NX.h5`` for
        text input).
    with_uncertainty : bool, optional
        If True, run *n_mc_runs* Gaussian-noise-perturbed refits after the main
        fit and store per-parameter 1σ uncertainties (``G_err``, ``Rg_err``,
        ``B_err``, ``P_err``, ``ETA_err``, ``PACK_err``, ``background_err``)
        in the NXcanSAS HDF5 file.  Default False.
    n_mc_runs : int, optional
        Number of Monte Carlo perturbation runs when *with_uncertainty* is True.
        Default 10.

    Returns
    -------
    dict or None
        On success, a dictionary with:

        ``'success'``        bool — True if fit converged.
        ``'model'``          UnifiedFitModel — the fitted model object.
        ``'parameters'``     dict — structured fit summary (levels, chi_squared, …).
        ``'uncertainties'``  dict or None — MC parameter uncertainties (None when
                             ``with_uncertainty=False`` or fewer than 2 runs succeeded).
        ``'data'``           dict — Q, Intensity, Error, intensity_model, residuals.
        ``'output_file'``    Path or None — NXcanSAS output path (if save_to_nexus).
        ``'input_file'``     Path
        ``'config_file'``    Path
        ``'message'``        str — human-readable summary.

        Returns None if loading, configuration, or fitting fails fatally.
        A failed fit (optimizer did not converge) still returns a dict with
        ``'success': False`` so the caller can inspect the partial result.

    Examples
    --------
    >>> result = fit_unified("sample.h5", "pyirena_config.json")
    >>> if result and result['success']:
    ...     print(result['parameters']['chi_squared'])

    >>> # With MC uncertainty analysis:
    >>> result = fit_unified("sample.h5", "pyirena_config.json",
    ...                      with_uncertainty=True, n_mc_runs=20)
    >>> if result and result['uncertainties']:
    ...     print(result['uncertainties']['levels'][0]['Rg'])
    """
    _ensure_console()
    data_file = Path(data_file)
    config_file = Path(config_file)

    # --- Load config ---
    try:
        config = _load_config(config_file)
        if config is None:
            return None

        if 'unified_fit' not in config:
            log.info(f"[pyirena.batch] Config file '{config_file}' has no 'unified_fit' group.")
            return None

        unified_state = config['unified_fit']
    except Exception:
        log.error(f"[pyirena.batch] Unexpected error reading config:\n{traceback.format_exc()}")
        return None

    # --- Load data ---
    try:
        data = _load_data(data_file)
        if data is None:
            return None
    except Exception:
        log.error(f"[pyirena.batch] Unexpected error loading data:\n{traceback.format_exc()}")
        return None

    # --- Build model ---
    try:
        model = _state_to_model(unified_state)
        num_levels = model.num_levels
    except Exception:
        log.error(f"[pyirena.batch] Error building model from config:\n{traceback.format_exc()}")
        return None

    # --- Apply Q range from cursor positions ---
    try:
        q_min = unified_state.get('cursor_left')
        q_max = unified_state.get('cursor_right')
        Q = data['Q']
        Intensity = data['Intensity']
        Error = data.get('Error')

        if q_min is not None and q_max is not None:
            mask = (Q >= float(q_min)) & (Q <= float(q_max))
            q_fit = Q[mask]
            intensity_fit = Intensity[mask]
            error_fit = Error[mask] if Error is not None else None
        else:
            q_fit, intensity_fit, error_fit = Q, Intensity, Error
    except Exception:
        log.error(f"[pyirena.batch] Error applying Q range:\n{traceback.format_exc()}")
        return None

    # --- Run fit ---
    try:
        fit_result = model.fit(q_fit, intensity_fit, error_fit)
    except Exception:
        log.error(f"[pyirena.batch] Fitting failed for '{data_file}':\n{traceback.format_exc()}")
        return None

    # --- MC uncertainty analysis (optional) ---
    uncertainties = None
    if with_uncertainty and fit_result.get('success', False):
        log.info(f"[pyirena.batch] Running {n_mc_runs} MC uncertainty runs for '{data_file.name}' ...")
        try:
            uncertainties = _mc_uncertainty_unified(
                model, q_fit, intensity_fit, error_fit, n_runs=n_mc_runs
            )
            if uncertainties is None:
                log.error("[pyirena.batch] MC uncertainty: fewer than 2 runs succeeded; skipping.")
        except Exception:
            log.error(f"[pyirena.batch] MC uncertainty failed:\n{traceback.format_exc()}")

    # --- Build return structure ---
    try:
        intensity_model = model.calculate_intensity(Q)
        if Error is not None:
            residuals = (Intensity - intensity_model) / Error
        else:
            residuals = (Intensity - intensity_model) / Intensity

        level_summaries = []
        for i, lv in enumerate(fit_result.get('levels', model.levels)):
            level_summaries.append({
                'level': i + 1,
                'G': lv.G, 'Rg': lv.Rg, 'B': lv.B, 'P': lv.P,
                'RgCutoff': lv.RgCO, 'ETA': lv.ETA, 'PACK': lv.PACK,
                'correlated': lv.correlations,
            })

        parameters = {
            'num_levels': num_levels,
            'background': fit_result.get('background', model.background),
            'chi_squared': fit_result.get('chi_squared', 0.0),
            'reduced_chi_squared': fit_result.get('reduced_chi_squared', 0.0),
            'n_iterations': fit_result.get('n_iterations', 0),
            'levels': level_summaries,
        }

        success = bool(fit_result.get('success', False))
        result = {
            'success': success,
            'tool': 'unified_fit',
            'input_file': data_file,
            'config_file': config_file,
            'output_file': None,
            'model': model,
            'fit_result': fit_result,
            'parameters': parameters,
            'uncertainties': uncertainties,
            'data': {
                'Q': Q,
                'Intensity': Intensity,
                'Error': Error,
                'Q_fit': q_fit,
                'intensity_model': intensity_model,
                'residuals': residuals,
            },
            'message': (
                f"Unified Fit: {num_levels} level(s), "
                f"chi²={parameters['chi_squared']:.4g}, "
                f"success={success}"
            ),
        }
    except Exception:
        log.error(f"[pyirena.batch] Error building result structure:\n{traceback.format_exc()}")
        return None

    # --- Save to NXcanSAS ---
    if save_to_nexus:
        try:
            setup_state = _build_setup_state(unified_state, fit_result)
            output_path = _save_to_nexus(
                data, model, fit_result, num_levels,
                uncertainties=uncertainties, setup_state=setup_state,
            )
            result['output_file'] = output_path
        except Exception:
            log.warning(f"[pyirena.batch] Warning: could not save NXcanSAS file:\n"
                  f"{traceback.format_exc()}")
            # Non-fatal: return result without output_file

    log.info(f"[pyirena.batch] {result['message']}")
    return result
