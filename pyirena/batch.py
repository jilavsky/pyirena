"""
pyirena.batch — headless (no-GUI) fitting API for scripting and automation.

Typical usage
-------------
Single file, one tool:
    from pyirena.batch import fit_unified
    result = fit_unified("data.h5", "pyirena_config.json")
    if result:
        print(result['parameters']['chi_squared'])

WAXS peak fitting (linear/linear, diffraction peaks):
    from pyirena.batch import fit_waxs
    result = fit_waxs("waxs_data.h5", "pyirena_config.json")
    if result and result['success']:
        print(result['n_peaks'], "peaks fitted")

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


def _compute_invariant_sv(G, Rg, B, P, RgCO, ETA, PACK, correlated):
    """
    Compute the scattering invariant and Sv for one Unified Fit level.

    Mirrors LevelParametersWidget.update_porod_surface_and_invariant() in
    unified_fit.py so that the batch path produces the same values as the GUI.

    Returns (invariant_cm4, sv) where either value is None when not computable.
    Sv is only meaningful when P is in the Porod regime (3.95 – 4.05).
    """
    from scipy.special import erf
    from scipy.integrate import simpson as _simpson

    if Rg <= 0 or B <= 0:
        return None, None
    try:
        maxQ = 2 * np.pi / (Rg / 10)
        surf_q = np.linspace(0, maxQ, 2000)

        # Unified intensity (matches _calculate_unified_intensity in unified_fit.py)
        K = 1.0 if P > 3 else 1.06
        q_safe = np.where(np.abs(surf_q) < 1e-10, 1e-10, surf_q)
        erf_cubed = erf(K * q_safe * Rg / np.sqrt(6)) ** 3
        erf_cubed = np.where(np.abs(erf_cubed) < 1e-10, 1e-10, erf_cubed)
        qstar = q_safe / erf_cubed
        qstar = np.where(np.isfinite(qstar), qstar, 1e-10)
        intensity = (G * np.exp(-surf_q**2 * Rg**2 / 3)
                     + (B / qstar**P) * np.exp(-RgCO**2 * surf_q**2 / 3))
        intensity = np.where(np.isfinite(intensity), intensity, 0.0)
        if correlated and PACK > 0 and ETA > 0:
            qr = np.where(surf_q * ETA == 0, 1e-10, surf_q * ETA)
            sphere_amp = 3 * (np.sin(qr) - qr * np.cos(qr)) / qr**3
            intensity = intensity / (1 + PACK * sphere_amp)
        if intensity[0] == 0 or np.isnan(intensity[0]):
            intensity[0] = intensity[1]

        # Invariant = ∫ I(Q)·Q² dQ
        invariant = _simpson(intensity * surf_q**2, x=surf_q)
        if RgCO < 0.1:
            invariant += -B * maxQ**(3 - abs(P)) / (3 - abs(P))
        if invariant <= 0:
            return None, None

        sv = 1e4 * np.pi * B / invariant if 3.95 <= P <= 4.05 else None
        return invariant * 1e24, sv   # convert to cm⁻⁴
    except Exception:
        return None, None


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


def _mc_uncertainty_sizes(
    s,
    q_fit: np.ndarray,
    intensity_fit: np.ndarray,
    error_fit: Optional[np.ndarray],
    n_runs: int = 10,
) -> Optional[np.ndarray]:
    """Run MC noise-perturbation fits and return per-bin std-dev of P(r).

    Returns a 1-D array with the same length as the r-grid, or ``None``
    when fewer than 2 runs succeed.
    """
    import copy

    err = error_fit if error_fit is not None else intensity_fit * 0.05
    err = np.maximum(err, 1e-30)

    # For Monte Carlo method: one internal MC run per perturbed fit (same as GUI)
    s_mc = copy.deepcopy(s)
    if s_mc.method == 'montecarlo':
        s_mc.montecarlo_n_repetitions = 1

    distributions = []
    for _ in range(n_runs):
        I_perturbed = intensity_fit + err * np.random.randn(len(intensity_fit))
        try:
            res = s_mc.fit(q_fit, I_perturbed, error_fit)
        except Exception:
            continue
        if not res.get('success', False):
            continue
        dist = res.get('distribution')
        if dist is not None:
            distributions.append(dist)

    if len(distributions) < 2:
        return None

    arr = np.array(distributions)          # (n_ok, n_bins)
    return np.std(arr, axis=0, ddof=1)


def _save_to_nexus(data: Dict, model: UnifiedFitModel,
                   fit_result: Dict, num_levels: int,
                   uncertainties: Optional[Dict] = None) -> Optional[Path]:
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

    # --- MC uncertainty analysis (optional) ---
    uncertainties = None
    if with_uncertainty and fit_result.get('success', False):
        print(f"[pyirena.batch] Running {n_mc_runs} MC uncertainty runs for '{data_file.name}' ...")
        try:
            uncertainties = _mc_uncertainty_unified(
                model, q_fit, intensity_fit, error_fit, n_runs=n_mc_runs
            )
            if uncertainties is None:
                print("[pyirena.batch] MC uncertainty: fewer than 2 runs succeeded; skipping.")
        except Exception:
            print(f"[pyirena.batch] MC uncertainty failed:\n{traceback.format_exc()}")

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
        print(f"[pyirena.batch] Error building result structure:\n{traceback.format_exc()}")
        return None

    # --- Save to NXcanSAS ---
    if save_to_nexus:
        try:
            output_path = _save_to_nexus(
                data, model, fit_result, num_levels, uncertainties=uncertainties
            )
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
    with_uncertainty: bool = False,
    n_mc_runs: int = 10,
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
    with_uncertainty : bool, optional
        If True, run *n_mc_runs* Gaussian-noise-perturbed refits after the main
        fit and store the per-bin 1σ standard deviation of the size distribution
        as ``distribution_std`` in the NXcanSAS HDF5 file.  For the Monte Carlo
        method each perturbed fit uses a single internal MC run.  Default False.
    n_mc_runs : int, optional
        Number of Monte Carlo perturbation runs when *with_uncertainty* is True.
        Default 10.

    Returns
    -------
    dict or None
        On success, a dictionary with:

        ``'success'``           bool — True if fit converged.
        ``'parameters'``        dict — chi_squared, volume_fraction, rg, peak_r, etc.
        ``'distribution_std'``  ndarray or None — per-bin 1σ uncertainty (None when
                                ``with_uncertainty=False`` or too few runs succeeded).
        ``'data'``              dict — Q, Intensity, Error, Q_fit, r_grid,
                                distribution, intensity_model, residuals.
        ``'output_file'``       Path or None — NXcanSAS output path (if save_to_nexus).
        ``'input_file'``        Path
        ``'config_file'``       Path
        ``'message'``           str — human-readable summary.

        Returns None if loading, configuration, or fitting fails fatally.

    Examples
    --------
    >>> result = fit_sizes("sample.h5", "pyirena_config.json")
    >>> if result and result['success']:
    ...     print(result['parameters']['rg'])

    >>> # With MC uncertainty analysis:
    >>> result = fit_sizes("sample.h5", "pyirena_config.json",
    ...                    with_uncertainty=True, n_mc_runs=20)
    >>> if result and result['distribution_std'] is not None:
    ...     import numpy as np
    ...     peak_r = result['data']['r_grid'][np.argmax(result['data']['distribution'])]
    ...     print(f"Peak radius = {peak_r:.1f} Å")
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

    # --- MC uncertainty analysis (optional) ---
    distribution_std = None
    if with_uncertainty and fit_result.get('success', False):
        print(f"[pyirena.batch] Running {n_mc_runs} MC uncertainty runs for '{data_file.name}' ...")
        try:
            distribution_std = _mc_uncertainty_sizes(
                s, q_fit, intensity_fit, error_fit, n_runs=n_mc_runs
            )
            if distribution_std is None:
                print("[pyirena.batch] MC uncertainty: fewer than 2 runs succeeded; skipping.")
        except Exception:
            print(f"[pyirena.batch] MC uncertainty failed:\n{traceback.format_exc()}")

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
            'montecarlo_n_repetitions': getattr(s, 'montecarlo_n_repetitions', 1),
            'montecarlo_convergence':   s.montecarlo_convergence,
            'montecarlo_max_iter':      s.montecarlo_max_iter,
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
            'success':          success,
            'tool':             'sizes',
            'input_file':       data_file,
            'config_file':      config_file,
            'output_file':      None,
            'fit_result':       fit_result,
            'parameters':       parameters,
            'distribution_std': distribution_std,
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

            # Use MC distribution_std when available; fall back to any std
            # produced by the method itself (e.g. Monte Carlo with n_repetitions > 1).
            std_to_save = distribution_std if distribution_std is not None \
                else fit_result.get('distribution_std')
            save_sizes_results(
                filepath=output_path,
                q=q_fit,
                intensity_data=intensity_fit,
                intensity_model=intensity_model,
                residuals=residuals,
                r_grid=r_grid,
                distribution=distribution,
                params=save_params,           # complete parameter set
                distribution_std=std_to_save,
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
    with_uncertainty: bool = False,
    n_mc_runs: int = 10,
) -> Optional[Dict]:
    """Run all analysis tools that have a configuration group in the config file.

    This is the top-level batch entry point.  It reads the config file, detects
    which tool sections are present, and dispatches to the appropriate fitting
    function for each tool.  Currently supported tools:

    ``unified_fit``
        Runs :func:`fit_unified`.
    ``sizes``
        Runs :func:`fit_sizes`.
    ``simple_fits``
        Runs :func:`fit_simple_from_config`.
    ``waxs_peakfit``
        Runs :func:`fit_waxs` (alias for :func:`fit_waxs_peaks_from_config`).

    Unknown sections in the config file are silently skipped.

    Parameters
    ----------
    data_file : str or Path
        Path to SAS data file.
    config_file : str or Path
        Path to a pyIrena JSON configuration file.
    save_to_nexus : bool, optional
        Passed to each individual tool's fitting function (default True).
    with_uncertainty : bool, optional
        Passed to each individual tool's fitting function (default False).
    n_mc_runs : int, optional
        Passed to each individual tool's fitting function (default 10).

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
        'unified_fit': lambda: fit_unified(
            data_file, config_file, save_to_nexus, with_uncertainty, n_mc_runs
        ),
        'sizes': lambda: fit_sizes(
            data_file, config_file, save_to_nexus, with_uncertainty, n_mc_runs
        ),
        'simple_fits': lambda: fit_simple_from_config(
            data_file, config_file, save_to_nexus, with_uncertainty, n_mc_runs
        ),
        'waxs_peakfit': lambda: fit_waxs(
            data_file, config_file, save_to_nexus, with_uncertainty, n_mc_runs
        ),
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


# ---------------------------------------------------------------------------
# fit_simple — headless Simple Fits API
# ---------------------------------------------------------------------------

def fit_simple(
    data_file: Union[str, Path],
    config: Union[Dict, 'SimpleFitModel'],
    with_uncertainty: bool = False,
    n_mc_runs: int = 50,
    q_min: Optional[float] = None,
    q_max: Optional[float] = None,
    verbose: bool = True,
) -> Optional[Dict]:
    """
    Fit a single analytical model to one SAS file and save the result.

    Parameters
    ----------
    data_file : str or Path
        Path to a text (.dat/.txt) or NXcanSAS HDF5 (.h5/.hdf5) file.
    config : dict or SimpleFitModel
        Model configuration.  If a dict, must have key ``'model'`` (model
        name string) and optionally ``'params'``, ``'limits'``,
        ``'use_complex_bg'``.  Example::

            {'model': 'Guinier', 'params': {'I0': 5.0, 'Rg': 80.0}}

        Alternatively pass a ``SimpleFitModel`` instance directly.
    with_uncertainty : bool
        If True, run a Monte Carlo uncertainty analysis after the initial
        fit.  Perturbs the data ``n_mc_runs`` times by one-sigma Gaussian
        noise and fits each perturbed dataset; the standard deviation of
        the resulting parameter distributions updates ``params_std`` in
        the returned result.
    n_mc_runs : int
        Number of Monte Carlo perturbation runs (used only when
        ``with_uncertainty=True``).
    q_min, q_max : float or None
        Q range [Å⁻¹] to use for fitting.  Data outside this range is
        excluded.  None means use all data.
    verbose : bool
        If True, print progress messages to stdout.

    Returns
    -------
    dict or None
        Fit result dict (same structure as ``SimpleFitModel.fit()``), with
        the result also saved to the HDF5 file.  Returns None on failure.
    """
    from pyirena.core.simple_fits import SimpleFitModel, MODEL_REGISTRY
    from pyirena.io.nxcansas_simple_fits import save_simple_fit_results

    data_file = Path(data_file)

    # ── Build model object ───────────────────────────────────────────────────
    if isinstance(config, SimpleFitModel):
        model = config
        if n_mc_runs != 50:
            model.n_mc_runs = n_mc_runs
    elif isinstance(config, dict):
        model = SimpleFitModel.from_dict(config)
        model.n_mc_runs = n_mc_runs
    else:
        print(f"[pyirena.batch.fit_simple] config must be a dict or SimpleFitModel, "
              f"got {type(config).__name__}")
        return None

    # ── Load data ────────────────────────────────────────────────────────────
    data = _load_data(data_file)
    if data is None:
        return None

    q = np.asarray(data.get('Q', data.get('q', [])), dtype=float)
    I = np.asarray(data.get('Intensity', data.get('intensity', [])), dtype=float)
    dI = data.get('Error', data.get('error', None))
    dI = np.asarray(dI, dtype=float) if dI is not None else None

    # ── Apply Q range mask ───────────────────────────────────────────────────
    mask = np.ones(len(q), dtype=bool)
    if q_min is not None:
        mask &= (q >= q_min)
    if q_max is not None:
        mask &= (q <= q_max)
    mask &= np.isfinite(q) & np.isfinite(I) & (q > 0)
    if dI is not None:
        mask &= np.isfinite(dI) & (dI > 0)

    if mask.sum() < 2:
        print(f"[pyirena.batch.fit_simple] Too few data points after Q masking "
              f"({mask.sum()} points, need ≥ 2).")
        return None

    qf = q[mask]
    If = I[mask]
    dIf = dI[mask] if dI is not None else None

    # ── Initial fit ──────────────────────────────────────────────────────────
    if verbose:
        print(f"[pyirena.batch] Fitting {model.model} to '{data_file.name}' ...")

    result = model.fit(qf, If, dIf)

    if not result.get('success'):
        msg = result.get('error', 'Unknown error')
        print(f"[pyirena.batch.fit_simple] Fit failed: {msg}")
        return result

    if verbose:
        rchi2 = result.get('reduced_chi2')
        print(f"[pyirena.batch] Fit succeeded.  Reduced χ² = {rchi2:.4g}")

    # ── Monte Carlo uncertainty ──────────────────────────────────────────────
    if with_uncertainty:
        if verbose:
            print(f"[pyirena.batch] Running {n_mc_runs} MC uncertainty runs ...")

        dIf_safe = (dIf if dIf is not None
                    else np.maximum(If * 0.05, 1e-30))
        mc_params: Dict[str, list] = {k: [] for k in result['params']}

        for _ in range(n_mc_runs):
            I_perturbed = If + dIf_safe * np.random.randn(len(If))
            mc_res = model.fit(qf, I_perturbed, dIf)
            if mc_res.get('success'):
                for k in mc_params:
                    mc_params[k].append(mc_res['params'].get(k, float('nan')))

        if any(len(v) > 1 for v in mc_params.values()):
            result['params_std'] = {
                k: float(np.std(v, ddof=1)) if len(v) > 1 else float('nan')
                for k, v in mc_params.items()
            }
            if verbose:
                print(f"[pyirena.batch] MC uncertainty done "
                      f"({sum(len(v) for v in mc_params.values()) // max(len(mc_params), 1)} "
                      f"successful runs).")

    # ── Save to HDF5 (only for NXcanSAS files) ───────────────────────────────
    ext = data_file.suffix.lower()
    if ext in ('.h5', '.hdf5', '.hdf', '.nx'):
        try:
            save_simple_fit_results(
                filepath=data_file,
                result=result,
                model_obj=model,
                intensity_data=I[mask] if mask is not None else I,
                intensity_error=dI[mask] if (dI is not None and mask is not None) else dI,
            )
            if verbose:
                print(f"[pyirena.batch] Results saved to '{data_file.name}'")
        except Exception as exc:
            print(f"[pyirena.batch.fit_simple] Could not save results: {exc}")

    return result


def fit_simple_from_config(
    data_file: Union[str, Path],
    config_file: Union[str, Path],
    save_to_nexus: bool = True,
    with_uncertainty: bool = False,
    n_mc_runs: int = 10,
) -> Optional[Dict]:
    """Fit a Simple Fits model using parameters from a pyIrena config file.

    Thin wrapper around :func:`fit_simple` that reads the ``simple_fits``
    section from a JSON config file and delegates to :func:`fit_simple`.
    This allows the same ``(data_file, config_file, save_to_nexus)`` calling
    convention used by :func:`fit_unified` and :func:`fit_sizes`.

    Parameters
    ----------
    data_file : str or Path
        Path to SAS data file (.h5/.hdf5 NXcanSAS, or .dat/.txt text).
    config_file : str or Path
        Path to a pyIrena JSON configuration file containing a
        ``'simple_fits'`` section (created by "Export Parameters" in the
        Simple Fits GUI panel).
    save_to_nexus : bool, optional
        If True (default), save fit results into the HDF5 file.
        Only HDF5 files are written; text-file results are not saved.
    with_uncertainty : bool, optional
        If True, run Monte Carlo uncertainty estimation.  Default False.
    n_mc_runs : int, optional
        Number of MC runs (used only when *with_uncertainty* is True).

    Returns
    -------
    dict or None
        Fit result dict (same structure as :meth:`SimpleFitModel.fit`) with
        ``'success'`` and ``'message'`` keys guaranteed.
        Returns None if the config cannot be read.
    """
    config_file = Path(config_file)
    config = _load_config(config_file)
    if config is None:
        print(f"[pyirena.batch.fit_simple_from_config] Cannot load config: {config_file}")
        return None

    sf_cfg = config.get('simple_fits')
    if sf_cfg is None:
        print(f"[pyirena.batch.fit_simple_from_config] No 'simple_fits' section in "
              f"'{config_file.name}'")
        return {'success': False, 'message': "No 'simple_fits' section in config file"}

    # Work on a shallow copy so we don't mutate the caller's config dict
    sf_cfg = dict(sf_cfg)

    # Extract optional Q range from config (stored by the GUI panel)
    q_min = sf_cfg.pop('q_min', None)
    q_max = sf_cfg.pop('q_max', None)

    # The GUI state uses 'param_limits'; SimpleFitModel.from_dict() expects 'limits'
    if 'param_limits' in sf_cfg and 'limits' not in sf_cfg:
        sf_cfg['limits'] = sf_cfg.pop('param_limits')
    else:
        sf_cfg.pop('param_limits', None)

    # Remove state-only keys that have no meaning for from_dict()
    for _k in ('schema_version', 'no_limits', 'param_fixed'):
        sf_cfg.pop(_k, None)

    result = fit_simple(
        data_file=data_file,
        config=sf_cfg,
        with_uncertainty=with_uncertainty,
        n_mc_runs=n_mc_runs,
        q_min=q_min,
        q_max=q_max,
        verbose=True,
    )

    if result is None:
        return {'success': False, 'message': 'fit_simple returned None (data load failure?)'}

    # Normalise error key so BatchWorker can always call result.get('message')
    if not result.get('success'):
        result.setdefault('message', result.get('error', 'fit failed'))

    return result


def fit_waxs_peaks_from_config(
    data_file: Union[str, Path],
    config_file: Union[str, Path],
    save_to_nexus: bool = True,
    with_uncertainty: bool = False,
    n_mc_runs: int = 10,
) -> Optional[Dict]:
    """Fit WAXS diffraction peaks using parameters from a pyIrena config file.

    Reads the ``'waxs_peakfit'`` section from a JSON config file and calls
    :func:`fit_waxs_peaks`.  The ``with_uncertainty`` and ``n_mc_runs``
    parameters are accepted for API consistency with :func:`fit_unified` and
    :func:`fit_sizes` but are not yet used (WAXS uncertainty comes from the
    ``curve_fit`` covariance matrix, not Monte Carlo).

    This function is also exported as :func:`fit_waxs` for convenience.

    Parameters
    ----------
    data_file : str or Path
        Path to a WAXS data file (HDF5/NXcanSAS or text).
    config_file : str or Path
        Path to a pyIrena JSON configuration file containing a
        ``'waxs_peakfit'`` section (written by the GUI's Export Parameters
        button or by hand).
    save_to_nexus : bool, optional
        Write fitted results to ``entry/waxs_peakfit_results`` in the HDF5
        file (default ``True``).
    with_uncertainty : bool, optional
        Accepted for API compatibility; not currently used.
    n_mc_runs : int, optional
        Accepted for API compatibility; not currently used.

    Returns
    -------
    dict
        Always returns a dict with ``'success'`` (bool) and ``'message'``
        (str) keys.  On success also contains ``'n_peaks'``, ``'bg_shape'``,
        ``'chi2'``, ``'reduced_chi2'``, ``'dof'``, ``'q_min'``, ``'q_max'``,
        ``'bg_params'``, ``'peaks'``, ``'I_fit'``, ``'I_bg'``, ``'q'``.
        Returns ``None`` only on a fatal pre-fit error (data unreadable).

    Examples
    --------
    >>> result = fit_waxs("waxs_data.h5", "pyirena_config.json")
    >>> if result and result['success']:
    ...     for pk in result['peaks']:
    ...         print(f"Q0={pk['Q0']['value']:.4f}  FWHM={pk['FWHM']['value']:.4f}")
    """
    config_file = Path(config_file)
    config = _load_config(config_file)
    if config is None:
        msg = f"Cannot load config: {config_file}"
        print(f"[pyirena.batch.fit_waxs_peaks_from_config] {msg}")
        return {'success': False, 'message': msg}

    wp_cfg = config.get('waxs_peakfit')
    if wp_cfg is None:
        msg = f"No 'waxs_peakfit' section in '{config_file.name}'"
        print(f"[pyirena.batch.fit_waxs_peaks_from_config] {msg}")
        return {'success': False, 'message': msg}

    q_min = wp_cfg.get('q_min')
    q_max = wp_cfg.get('q_max')

    result = fit_waxs_peaks(
        data_file=data_file,
        config=wp_cfg,
        q_min=q_min,
        q_max=q_max,
        save_to_nexus=save_to_nexus,
        verbose=True,
    )

    if result is None:
        return {'success': False, 'message': 'Data load or fit failure (fit_waxs_peaks returned None)'}

    if not result.get('success'):
        result.setdefault('message', result.get('error', 'fit failed'))

    return result


#: Short alias for :func:`fit_waxs_peaks_from_config`.
#: Follows the same ``(data_file, config_file, save_to_nexus)`` convention
#: as ``fit_unified``, ``fit_sizes``, and ``fit_simple_from_config``.
fit_waxs = fit_waxs_peaks_from_config


# ---------------------------------------------------------------------------
# fit_waxs_peaks — headless WAXS peak fitting
# ---------------------------------------------------------------------------

def fit_waxs_peaks(
    data_file: Union[str, Path],
    config: Dict,
    q_min: Optional[float] = None,
    q_max: Optional[float] = None,
    save_to_nexus: bool = True,
    verbose: bool = True,
) -> Optional[Dict]:
    """Fit WAXS diffraction peaks to one data file and (optionally) save.

    Parameters
    ----------
    data_file : str or Path
        Path to a text (.dat/.txt) or NXcanSAS HDF5 (.h5/.hdf5) file.
    config : dict
        waxs_peakfit state dict.  Must contain at least ``'bg_shape'`` and
        ``'peaks'`` (list of peak param dicts from
        ``pyirena.core.waxs_peakfit.default_peak_params()``).
        Optionally: ``'no_limits'``, ``'peak_find'`` (used only when
        ``peaks`` is empty and auto-finding is requested).
    q_min, q_max : float or None
        Restrict fitting to this Q range (Å⁻¹).  None = use all data.
    save_to_nexus : bool
        When True and *data_file* is an HDF5 file, save results back into
        ``entry/waxs_peakfit_results``.
    verbose : bool
        Print progress messages.

    Returns
    -------
    dict or None
        Result dict from ``WAXSPeakFitModel.fit()``, or None on failure.
    """
    from pyirena.core.waxs_peakfit import (
        WAXSPeakFitModel, find_peaks_in_data, default_bg_params,
    )
    from pyirena.io.hdf5 import readGenericNXcanSAS, readTextFile

    data_file = Path(data_file)
    if verbose:
        print(f"[pyirena.batch.fit_waxs_peaks] {data_file.name}")

    # ── Load data ────────────────────────────────────────────────────────────
    try:
        ext = data_file.suffix.lower()
        if ext in ('.txt', '.dat'):
            data = readTextFile(str(data_file.parent), data_file.name)
        else:
            data = readGenericNXcanSAS(str(data_file.parent), data_file.name)
        if data is None:
            if verbose:
                print(f"  [fit_waxs_peaks] Could not load data from {data_file.name}")
            return None
    except Exception as exc:
        if verbose:
            print(f"  [fit_waxs_peaks] Load error: {exc}")
        return None

    q   = np.asarray(data['Q'],        float)
    I   = np.asarray(data['Intensity'], float)
    dI  = np.asarray(data.get('Error', np.ones_like(I) * np.nan), float)

    # ── Q range crop ─────────────────────────────────────────────────────────
    mask = np.isfinite(q) & np.isfinite(I) & (q > 0)
    if q_min is not None:
        mask &= q >= q_min
    if q_max is not None:
        mask &= q <= q_max
    q_, I_ = q[mask], I[mask]
    dI_    = dI[mask] if dI is not None else None

    if len(q_) < 5:
        if verbose:
            print(f"  [fit_waxs_peaks] Too few data points in Q range.")
        return None

    # ── Build model config ────────────────────────────────────────────────────
    bg_shape  = config.get('bg_shape', 'Constant')
    no_limits = config.get('no_limits', False)
    bg_params = config.get('bg_params', default_bg_params(bg_shape))
    peaks     = config.get('peaks', [])

    # Auto-find peaks if none provided and peak_find config given
    if not peaks and 'peak_find' in config:
        pf = config['peak_find']
        peaks = find_peaks_in_data(
            q_, I_,
            prominence_frac=pf.get('prominence',    0.05),
            min_fwhm=       pf.get('min_fwhm',      0.001),
            max_fwhm=       pf.get('max_fwhm',      0.5),
            min_distance=   pf.get('min_distance',  0.005),
            sg_window_frac= pf.get('sg_window_frac', 0.15),
        )
        if verbose:
            print(f"  [fit_waxs_peaks] Auto-found {len(peaks)} peak(s).")

    if not peaks:
        if verbose:
            print(f"  [fit_waxs_peaks] No peaks defined; fitting background only.")
        peaks = []

    # ── Fit ───────────────────────────────────────────────────────────────────
    engine = WAXSPeakFitModel(bg_shape=bg_shape, peaks=peaks, no_limits=no_limits)
    try:
        result = engine.fit(q_, I_, dI_ if (dI_ is not None and np.any(np.isfinite(dI_))) else None,
                            bg_params=bg_params, peaks=peaks)
    except Exception as exc:
        if verbose:
            print(f"  [fit_waxs_peaks] Fit error: {exc}")
        return None

    if verbose:
        status = "OK" if result.get('success') else "FAILED"
        print(f"  [fit_waxs_peaks] {status}  "
              f"reduced-χ² = {result.get('reduced_chi2', float('nan')):.4g}")

    # ── Save to HDF5 ──────────────────────────────────────────────────────────
    if save_to_nexus and data_file.suffix.lower() in ('.h5', '.hdf5', '.hdf'):
        try:
            from pyirena.io.nxcansas_waxs_peakfit import save_waxs_peakfit_results
            save_waxs_peakfit_results(
                data_file, result, q_,
                intensity_data=I_,
                intensity_error=dI_ if dI_ is not None else None,
                q_min=float(q_.min()),
                q_max=float(q_.max()),
            )
            if verbose:
                print(f"  [fit_waxs_peaks] Saved to {data_file.name}.")
        except Exception as exc:
            if verbose:
                print(f"  [fit_waxs_peaks] Save error: {exc}")

    return result
