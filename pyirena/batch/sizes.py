"""
pyirena.batch.sizes — Headless Size Distribution fitting (fit_sizes).

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
    _ensure_console()
    from pyirena.core.sizes import SizesDistribution

    data_file = Path(data_file)
    config_file = Path(config_file)

    # --- Load config ---
    try:
        config = _load_config(config_file)
        if config is None:
            return None

        if 'sizes' not in config:
            log.info(f"[pyirena.batch] Config file '{config_file}' has no 'sizes' group.")
            return None

        sizes_state = config['sizes']
    except Exception:
        log.error(f"[pyirena.batch] Unexpected error reading config:\n{traceback.format_exc()}")
        return None

    # --- Load data ---
    # JSON keys (canonical): "load_slit_smeared" picks the file's _SMR dataset;
    # "slit_length" optionally overrides the file-derived value.
    load_slit_smeared = bool(sizes_state.get('load_slit_smeared', False))
    try:
        data = _load_data(data_file, load_slit_smeared=load_slit_smeared)
        if data is None:
            return None
    except Exception:
        log.error(f"[pyirena.batch] Unexpected error loading data:\n{traceback.format_exc()}")
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
        s.fractional_error   = bool(sizes_state.get('fractional_error', False))
        s.fractional_error_value = float(sizes_state.get('fractional_error_value', 0.03))
        s.power_law_B        = float(sizes_state.get('power_law_B', 0.0))
        s.power_law_P        = float(sizes_state.get('power_law_P', 4.0))
        s.method             = str(sizes_state.get('method', 'regularization'))
        s.maxent_sky_background  = float(sizes_state.get('maxent_sky_background', 1e-6))
        s.maxent_stability       = float(sizes_state.get('maxent_stability', 0.01))
        s.maxent_max_iter        = int(sizes_state.get('maxent_max_iter', 300))
        s.regularization_evalue  = float(sizes_state.get('regularization_evalue', 1.0))
        s.regularization_min_ratio = float(sizes_state.get('regularization_min_ratio', 1e-4))
        s.tnnls_approach_param   = float(sizes_state.get('tnnls_approach_param', 0.95))
        s.tnnls_max_iter         = int(sizes_state.get('tnnls_max_iter', 300))
        s.montecarlo_n_repetitions = 1  # main fit always uses a single MC run, matching GUI
        # Slit smearing: enable when the loaded data are slit smeared or the
        # config asks for it; slit length is file-derived unless overridden.
        cfg_sl = sizes_state.get('slit_length')
        sl = float(cfg_sl) if cfg_sl else float(data.get('slit_length', 0.0) or 0.0)
        if (bool(data.get('is_slit_smeared')) or bool(sizes_state.get('use_slit_smearing'))) and sl > 0:
            s.use_slit_smearing = True
            s.slit_length = sl
            log.info(f"[pyirena.batch] Sizes slit smearing enabled (SL={sl:.4g} 1/A).")
    except Exception:
        log.error(f"[pyirena.batch] Error building Sizes model from config:\n{traceback.format_exc()}")
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
        log.error(f"[pyirena.batch] Error applying Q range:\n{traceback.format_exc()}")
        return None

    # --- Optional background pre-fits (mirrors the GUI "Fit All" sequence) ---
    # 1. Power-law B·q^(-P): fit B and/or P over the power-law Q range when the
    #    config's fit flags are set, updating s.power_law_B / s.power_law_P.
    # 2. Flat background: averaged over the background Q range when one is set,
    #    updating s.background.
    # Both pre-fits use the full (un-masked) data Q so the chosen Q windows are
    # honoured independently of the size-fit cursor range, exactly as the GUI does.
    try:
        fit_B = bool(sizes_state.get('fit_power_law_B', False))
        fit_P = bool(sizes_state.get('fit_power_law_P', False))
        if fit_B or fit_P:
            pl_qmin = sizes_state.get('power_law_q_min')
            pl_qmax = sizes_state.get('power_law_q_max')
            if pl_qmin is None or pl_qmax is None:
                pl_qmin, pl_qmax = float(Q.min()), float(Q.max())
            pl_res = s.fit_power_law(
                Q, Intensity, float(pl_qmin), float(pl_qmax),
                fit_B=fit_B, fit_P=fit_P,
            )
            log.info(f"[pyirena.batch] Power-law pre-fit: {pl_res.get('message', '')}")

        bg_qmin = sizes_state.get('background_q_min')
        bg_qmax = sizes_state.get('background_q_max')
        if bg_qmin is not None and bg_qmax is not None:
            bg_res = s.fit_background_term(
                Q, Intensity, float(bg_qmin), float(bg_qmax),
            )
            log.info(f"[pyirena.batch] Background pre-fit: {bg_res.get('message', '')}")
    except Exception:
        log.error(f"[pyirena.batch] Background pre-fit failed (continuing with "
              f"configured values):\n{traceback.format_exc()}")

    # --- Run fit ---
    try:
        fit_result = s.fit(q_fit, intensity_fit, error_fit)
    except Exception:
        log.error(f"[pyirena.batch] Fitting failed for '{data_file}':\n{traceback.format_exc()}")
        return None

    # --- MC uncertainty analysis (optional) ---
    distribution_std = None
    if with_uncertainty and fit_result.get('success', False):
        log.info(f"[pyirena.batch] Running {n_mc_runs} MC uncertainty runs for '{data_file.name}' ...")
        try:
            distribution_std = _mc_uncertainty_sizes(
                s, q_fit, intensity_fit, error_fit, n_runs=n_mc_runs
            )
            if distribution_std is None:
                log.error("[pyirena.batch] MC uncertainty: fewer than 2 runs succeeded; skipping.")
        except Exception:
            log.error(f"[pyirena.batch] MC uncertainty failed:\n{traceback.format_exc()}")

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
            # Slit-smearing provenance
            'slit_length':          float(s.slit_length) if s.use_slit_smearing else 0.0,
            'data_is_slit_smeared': bool(s.use_slit_smearing),
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
            'fractional_error':       s.fractional_error,
            'fractional_error_value': s.fractional_error_value,
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
        log.error(f"[pyirena.batch] Error building result structure:\n{traceback.format_exc()}")
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

            # Robust fit-quality metrics on the consistent raw-basis triple
            # (observed I, model + complex background, error) — same as the GUI.
            fq_metrics = None
            try:
                from pyirena.core.fit_metrics import fit_quality_metrics
                _q_q = fit_result.get('q', q_fit)
                _I_obs = fit_result.get('I_data')
                _err = fit_result.get('err')
                if _I_obs is not None:
                    _bg = s.compute_complex_background(_q_q)
                    fq_metrics = fit_quality_metrics(
                        _q_q, _I_obs, intensity_model + _bg, _err, n_params=1)
            except Exception:
                fq_metrics = None

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
                fit_quality=fq_metrics,
                setup_state=sizes_state,
                intensity_model_ideal=fit_result.get('model_intensity_ideal'),
            )
            result['output_file'] = output_path
        except Exception:
            log.warning(f"[pyirena.batch] Warning: could not save NXcanSAS file:\n"
                  f"{traceback.format_exc()}")

    log.info(f"[pyirena.batch] {result['message']}")
    return result
