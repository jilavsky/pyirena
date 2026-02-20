"""
pyirena.batch — headless (no-GUI) fitting API for scripting and automation.

Typical usage
-------------
Single file, one tool:
    from pyirena.batch import fit_unified
    result = fit_unified("data.h5", "pyirena_config.json")
    if result:
        print(result['parameters']['chi_squared'])

Single file, all configured tools:
    from pyirena.batch import fit_pyirena
    results = fit_pyirena("data.h5", "pyirena_config.json")

Batch over many files:
    results = [fit_pyirena(f, cfg) for f in data_files]
    results = [r for r in results if r is not None]  # filter failures
"""

from __future__ import annotations

import json
import os
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np

from pyirena import __version__
from pyirena.core.unified import UnifiedFitModel, UnifiedLevel

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_config(config_file: Union[str, Path]) -> Optional[Dict]:
    """Load and validate a pyIrena JSON config file.  Returns None on failure."""
    config_file = Path(config_file)
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
    except Exception as e:
        print(f"[pyirena.batch] Cannot read config file '{config_file}': {e}")
        return None

    if '_pyirena_config' not in config:
        print(f"[pyirena.batch] '{config_file}' is not a pyIrena configuration file "
              f"(missing '_pyirena_config' header).")
        return None

    return config


def _load_data(data_file: Union[str, Path]) -> Optional[Dict]:
    """Load SAS data from a text (.dat/.txt) or HDF5 (.h5/.hdf5) file.

    Returns a dict with keys: Q, Intensity, Error (may be None).
    """
    from pyirena.io.hdf5 import readGenericNXcanSAS, readTextFile

    data_file = Path(data_file)
    if not data_file.exists():
        print(f"[pyirena.batch] Data file not found: '{data_file}'")
        return None

    path = str(data_file.parent)
    filename = data_file.name
    ext = data_file.suffix.lower()

    try:
        if ext in ('.txt', '.dat'):
            data = readTextFile(path, filename)
            is_nxcansas = False
        else:
            data = readGenericNXcanSAS(path, filename)
            is_nxcansas = True
    except Exception as e:
        print(f"[pyirena.batch] Error reading '{data_file}': {e}")
        return None

    if data is None:
        print(f"[pyirena.batch] Could not read data from '{data_file}'")
        return None

    data['filepath'] = str(data_file)
    data['is_nxcansas'] = is_nxcansas
    data.setdefault('label', data_file.stem)
    return data


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


def _save_to_nexus(data: Dict, model: UnifiedFitModel,
                   fit_result: Dict, num_levels: int) -> Optional[Path]:
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
        levels.append({
            'G': lv.G, 'Rg': lv.Rg, 'B': lv.B, 'P': lv.P,
            'RgCutoff': lv.RgCO, 'ETA': lv.ETA, 'PACK': lv.PACK,
            'correlated': lv.correlations,
        })

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
    )

    return output_path


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fit_unified(
    data_file: Union[str, Path],
    config_file: Union[str, Path],
    save_to_nexus: bool = True,
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

    Returns
    -------
    dict or None
        On success, a dictionary with:

        ``'success'``      bool — True if fit converged.
        ``'model'``        UnifiedFitModel — the fitted model object.
        ``'parameters'``   dict — structured fit summary (levels, chi_squared, …).
        ``'data'``         dict — Q, Intensity, Error, intensity_model, residuals.
        ``'output_file'``  Path or None — NXcanSAS output path (if save_to_nexus).
        ``'input_file'``   Path
        ``'config_file'``  Path
        ``'message'``      str — human-readable summary.

        Returns None if loading, configuration, or fitting fails fatally.
        A failed fit (optimizer did not converge) still returns a dict with
        ``'success': False`` so the caller can inspect the partial result.

    Examples
    --------
    >>> result = fit_unified("sample.h5", "pyirena_config.json")
    >>> if result and result['success']:
    ...     print(result['parameters']['chi_squared'])
    """
    data_file = Path(data_file)
    config_file = Path(config_file)

    # --- Load config ---
    try:
        config = _load_config(config_file)
        if config is None:
            return None

        if 'unified_fit' not in config:
            print(f"[pyirena.batch] Config file '{config_file}' has no 'unified_fit' group.")
            return None

        unified_state = config['unified_fit']
    except Exception:
        print(f"[pyirena.batch] Unexpected error reading config:\n{traceback.format_exc()}")
        return None

    # --- Load data ---
    try:
        data = _load_data(data_file)
        if data is None:
            return None
    except Exception:
        print(f"[pyirena.batch] Unexpected error loading data:\n{traceback.format_exc()}")
        return None

    # --- Build model ---
    try:
        model = _state_to_model(unified_state)
        num_levels = model.num_levels
    except Exception:
        print(f"[pyirena.batch] Error building model from config:\n{traceback.format_exc()}")
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
        print(f"[pyirena.batch] Error applying Q range:\n{traceback.format_exc()}")
        return None

    # --- Run fit ---
    try:
        fit_result = model.fit(q_fit, intensity_fit, error_fit)
    except Exception:
        print(f"[pyirena.batch] Fitting failed for '{data_file}':\n{traceback.format_exc()}")
        return None

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

        result = {
            'success': bool(fit_result.get('success', False)),
            'tool': 'unified_fit',
            'input_file': data_file,
            'config_file': config_file,
            'output_file': None,
            'model': model,
            'fit_result': fit_result,
            'parameters': parameters,
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
                f"success={result['success'] if False else fit_result.get('success', False)}"
            ),
        }
        # Fix message after dict is built
        result['message'] = (
            f"Unified Fit: {num_levels} level(s), "
            f"chi²={parameters['chi_squared']:.4g}, "
            f"success={result['success']}"
        )
    except Exception:
        print(f"[pyirena.batch] Error building result structure:\n{traceback.format_exc()}")
        return None

    # --- Save to NXcanSAS ---
    if save_to_nexus:
        try:
            output_path = _save_to_nexus(data, model, fit_result, num_levels)
            result['output_file'] = output_path
        except Exception:
            print(f"[pyirena.batch] Warning: could not save NXcanSAS file:\n"
                  f"{traceback.format_exc()}")
            # Non-fatal: return result without output_file

    print(f"[pyirena.batch] {result['message']}")
    return result


def fit_sizes(
    data_file: Union[str, Path],
    config_file: Union[str, Path],
    save_to_nexus: bool = True,
) -> Optional[Dict]:
    """Fit a Size Distribution model to a data file using a pyIrena config file.

    Parameters
    ----------
    data_file : str or Path
        Path to SAS data file (.h5/.hdf5 NXcanSAS, or .dat/.txt text).
    config_file : str or Path
        Path to a pyIrena JSON configuration file (created by "Export Parameters"
        in the Sizes GUI, or written manually).  Must contain a ``'sizes'`` group.
    save_to_nexus : bool, optional
        If True (default), save fit results to an NXcanSAS HDF5 file alongside
        the input data (same file for NXcanSAS input, or ``<name>_NX.h5`` for
        text input).

    Returns
    -------
    dict or None
        On success, a dictionary with:

        ``'success'``       bool — True if fit converged.
        ``'parameters'``    dict — chi_squared, volume_fraction, rg, peak_r, etc.
        ``'data'``          dict — Q, Intensity, Error, Q_fit, r_grid,
                            distribution, intensity_model, residuals.
        ``'output_file'``   Path or None — NXcanSAS output path (if save_to_nexus).
        ``'input_file'``    Path
        ``'config_file'``   Path
        ``'message'``       str — human-readable summary.

        Returns None if loading, configuration, or fitting fails fatally.

    Examples
    --------
    >>> result = fit_sizes("sample.h5", "pyirena_config.json")
    >>> if result and result['success']:
    ...     print(result['parameters']['rg'])
    """
    from pyirena.core.sizes import SizesDistribution

    data_file = Path(data_file)
    config_file = Path(config_file)

    # --- Load config ---
    try:
        config = _load_config(config_file)
        if config is None:
            return None

        if 'sizes' not in config:
            print(f"[pyirena.batch] Config file '{config_file}' has no 'sizes' group.")
            return None

        sizes_state = config['sizes']
    except Exception:
        print(f"[pyirena.batch] Unexpected error reading config:\n{traceback.format_exc()}")
        return None

    # --- Load data ---
    try:
        data = _load_data(data_file)
        if data is None:
            return None
    except Exception:
        print(f"[pyirena.batch] Unexpected error loading data:\n{traceback.format_exc()}")
        return None

    # --- Build SizesDistribution from config ---
    try:
        s = SizesDistribution()
        s.r_min              = float(sizes_state.get('r_min', 10.0))
        s.r_max              = float(sizes_state.get('r_max', 1000.0))
        s.n_bins             = int(sizes_state.get('n_bins', 200))
        s.log_spacing        = bool(sizes_state.get('log_spacing', True))
        s.shape              = str(sizes_state.get('shape', 'sphere'))
        s.contrast           = float(sizes_state.get('contrast', 1.0))
        ar = sizes_state.get('aspect_ratio', 1.0)
        if s.shape == 'spheroid':
            s.shape_params = {'aspect_ratio': float(ar)}
        s.background         = float(sizes_state.get('background', 0.0))
        s.error_scale        = float(sizes_state.get('error_scale', 1.0))
        s.power_law_B        = float(sizes_state.get('power_law_B', 0.0))
        s.power_law_P        = float(sizes_state.get('power_law_P', 4.0))
        s.method             = str(sizes_state.get('method', 'regularization'))
        s.maxent_sky_background  = float(sizes_state.get('maxent_sky_background', 1e-6))
        s.maxent_stability       = float(sizes_state.get('maxent_stability', 0.01))
        s.maxent_max_iter        = int(sizes_state.get('maxent_max_iter', 1000))
        s.regularization_evalue  = float(sizes_state.get('regularization_evalue', 1.0))
        s.regularization_min_ratio = float(sizes_state.get('regularization_min_ratio', 1e-4))
        s.tnnls_approach_param   = float(sizes_state.get('tnnls_approach_param', 0.95))
        s.tnnls_max_iter         = int(sizes_state.get('tnnls_max_iter', 1000))
    except Exception:
        print(f"[pyirena.batch] Error building Sizes model from config:\n{traceback.format_exc()}")
        return None

    # --- Apply Q range from saved cursor positions ---
    try:
        q_min = sizes_state.get('cursor_q_min')
        q_max = sizes_state.get('cursor_q_max')
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
        print(f"[pyirena.batch] Error applying Q range:\n{traceback.format_exc()}")
        return None

    # --- Run fit ---
    try:
        fit_result = s.fit(q_fit, intensity_fit, error_fit)
    except Exception:
        print(f"[pyirena.batch] Fitting failed for '{data_file}':\n{traceback.format_exc()}")
        return None

    # --- Build return structure ---
    try:
        r_grid       = fit_result.get('r_grid')
        distribution = fit_result.get('distribution')
        chi2         = fit_result.get('chi_squared', float('nan'))
        vf           = fit_result.get('volume_fraction', float('nan'))
        rg           = fit_result.get('rg', float('nan'))
        intensity_model = fit_result.get('model_intensity', np.zeros_like(q_fit))
        residuals    = fit_result.get('residuals', np.zeros_like(q_fit))

        peak_r = float('nan')
        if distribution is not None and r_grid is not None and len(distribution) > 0:
            peak_r = float(r_grid[int(np.argmax(distribution))])

        # Build a complete params dict matching what the GUI saves via _get_current_state().
        # This is passed both to save_sizes_results() and returned as 'parameters'.
        ar = s.shape_params.get('aspect_ratio', 1.0) if s.shape == 'spheroid' else None
        save_params = {
            # Fit results
            'chi_squared':     chi2,
            'volume_fraction': vf,
            'rg':              rg,
            'n_iterations':    fit_result.get('n_iterations', 0),
            # Model / grid setup
            'method':          s.method,
            'shape':           s.shape,
            'contrast':        s.contrast,
            'aspect_ratio':    ar,
            'r_min':           s.r_min,
            'r_max':           s.r_max,
            'n_bins':          s.n_bins,
            'log_spacing':     s.log_spacing,
            'background':      s.background,
            'error_scale':     s.error_scale,
            'power_law_B':     s.power_law_B,
            'power_law_P':     s.power_law_P,
            # Method-specific parameters
            'maxent_sky_background':    s.maxent_sky_background,
            'maxent_stability':         s.maxent_stability,
            'maxent_max_iter':          s.maxent_max_iter,
            'regularization_evalue':    s.regularization_evalue,
            'regularization_min_ratio': s.regularization_min_ratio,
            'tnnls_approach_param':     s.tnnls_approach_param,
            'tnnls_max_iter':           s.tnnls_max_iter,
            'mcsas_n_repetitions':      getattr(s, 'mcsas_n_repetitions', 1),
            'mcsas_convergence':        s.mcsas_convergence,
            'mcsas_max_iter':           s.mcsas_max_iter,
            # Q ranges from config (cursor / background / power-law)
            'cursor_q_min':      sizes_state.get('cursor_q_min'),
            'cursor_q_max':      sizes_state.get('cursor_q_max'),
            'power_law_q_min':   sizes_state.get('power_law_q_min'),
            'power_law_q_max':   sizes_state.get('power_law_q_max'),
            'background_q_min':  sizes_state.get('background_q_min'),
            'background_q_max':  sizes_state.get('background_q_max'),
        }

        # Convenience summary for the caller (add computed peak_r)
        parameters = dict(save_params)
        parameters['peak_r'] = peak_r

        success = bool(fit_result.get('success', False))
        result = {
            'success':     success,
            'tool':        'sizes',
            'input_file':  data_file,
            'config_file': config_file,
            'output_file': None,
            'fit_result':  fit_result,
            'parameters':  parameters,
            'data': {
                'Q':               Q,
                'Intensity':       Intensity,
                'Error':           Error,
                'Q_fit':           q_fit,
                'r_grid':          r_grid,
                'distribution':    distribution,
                'intensity_model': intensity_model,
                'residuals':       residuals,
            },
            'message': (
                f"Sizes ({s.method}): chi²={chi2:.4g}, "
                f"Vf={vf:.4g}, Rg={rg:.4g} Å, peak_r={peak_r:.4g} Å, "
                f"success={success}"
            ),
        }
    except Exception:
        print(f"[pyirena.batch] Error building result structure:\n{traceback.format_exc()}")
        return None

    # --- Save to NXcanSAS ---
    if save_to_nexus and success:
        try:
            from pyirena.io.nxcansas_sizes import save_sizes_results
            from pyirena.io.nxcansas_unified import (
                get_output_filepath, create_nxcansas_file
            )

            source_path = Path(data['filepath'])
            is_nxcansas = data.get('is_nxcansas', False)
            output_path = get_output_filepath(source_path, is_nxcansas)

            if not is_nxcansas and not output_path.exists():
                label = data.get('label', 'data')
                create_nxcansas_file(
                    output_path, Q, Intensity, Error,
                    sample_name=Path(label).stem
                )

            save_sizes_results(
                filepath=output_path,
                q=q_fit,
                intensity_data=intensity_fit,
                intensity_model=intensity_model,
                residuals=residuals,
                r_grid=r_grid,
                distribution=distribution,
                params=save_params,           # complete parameter set
                distribution_std=fit_result.get('distribution_std'),
            )
            result['output_file'] = output_path
        except Exception:
            print(f"[pyirena.batch] Warning: could not save NXcanSAS file:\n"
                  f"{traceback.format_exc()}")

    print(f"[pyirena.batch] {result['message']}")
    return result


def fit_pyirena(
    data_file: Union[str, Path],
    config_file: Union[str, Path],
    save_to_nexus: bool = True,
) -> Optional[Dict]:
    """Run all analysis tools that have a configuration group in the config file.

    This is the top-level batch entry point.  It reads the config file, detects
    which tool sections are present, and dispatches to the appropriate fitting
    function for each tool.  Currently supported tools:

    ``unified_fit``
        Runs :func:`fit_unified`.
    ``sizes``
        Runs :func:`fit_sizes`.

    Additional tool sections will be dispatched automatically as they are added
    to pyIrena.  Unknown sections in the config file are silently skipped.

    Parameters
    ----------
    data_file : str or Path
        Path to SAS data file.
    config_file : str or Path
        Path to a pyIrena JSON configuration file.
    save_to_nexus : bool, optional
        Passed to each individual tool's fitting function (default True).

    Returns
    -------
    dict or None
        On success, a dictionary with:

        ``'input_file'``   Path
        ``'config_file'``  Path
        ``'tools_run'``    list of str — tool names that were executed.
        ``'results'``      dict mapping tool name → result dict (or None on failure).

        Returns None if the config file cannot be loaded or contains no known tools.

    Examples
    --------
    >>> results = fit_pyirena("sample.h5", "pyirena_config.json")
    >>> if results:
    ...     uf = results['results'].get('unified_fit')
    ...     if uf and uf['success']:
    ...         print(uf['parameters']['chi_squared'])
    """
    data_file = Path(data_file)
    config_file = Path(config_file)

    # Registry: config key → callable
    _TOOL_REGISTRY: Dict[str, callable] = {
        'unified_fit': lambda: fit_unified(data_file, config_file, save_to_nexus),
        'sizes':       lambda: fit_sizes(data_file, config_file, save_to_nexus),
    }

    # --- Load config to discover which tools are present ---
    try:
        config = _load_config(config_file)
        if config is None:
            return None
    except Exception:
        print(f"[pyirena.batch] Unexpected error reading config:\n{traceback.format_exc()}")
        return None

    tools_to_run = [key for key in _TOOL_REGISTRY if key in config]

    if not tools_to_run:
        print(f"[pyirena.batch] Config file '{config_file}' contains no recognised "
              f"tool sections. Known tools: {list(_TOOL_REGISTRY)}")
        return None

    # --- Run each tool ---
    all_results: Dict[str, Optional[Dict]] = {}
    for tool in tools_to_run:
        print(f"[pyirena.batch] Running '{tool}' on '{data_file.name}' ...")
        try:
            all_results[tool] = _TOOL_REGISTRY[tool]()
        except Exception:
            print(f"[pyirena.batch] Unhandled error in '{tool}':\n{traceback.format_exc()}")
            all_results[tool] = None

    return {
        'input_file': data_file,
        'config_file': config_file,
        'tools_run': tools_to_run,
        'results': all_results,
    }
