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
    seed=None,
) -> Dict:
    """MC uncertainties for a Sizes fit — a thin wrapper over
    :meth:`SizesDistribution.calculate_uncertainty_mc`, kept so existing
    callers keep their import.

    Returns the core result dict; ``distribution_std`` is the per-bin 1σ.
    Perturbs with the model's *effective* error (fractional error / error
    scale honoured), as the GUI always did; until 1.2.1 this path used the
    raw error. Without measured uncertainties, σ = I × the Data Selector's
    "Generated uncertainty fraction" (issue #38).
    """
    from pyirena.state.state_manager import generated_noise_sigma
    return s.calculate_uncertainty_mc(
        q_fit, intensity_fit, error_fit, n_runs=n_runs, seed=seed,
        noise_sigma=generated_noise_sigma(intensity_fit, error_fit))


def fit_sizes(
    data_file: Union[str, Path],
    config_file: Union[str, Path],
    save_to_nexus: bool = True,
    with_uncertainty: bool = False,
    n_mc_runs: Optional[int] = None,
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
        Default None, which takes the config's own ``unc_n_runs`` (the GUI's
        uncertainty-run spin box), falling back to 10 when it states none.

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

    data_file = Path(data_file)
    config_file = Path(config_file)

    # --- Load config (the one envelope reader, AGENTS.md invariant 6) ---
    from pyirena.core.tool_config import (
        ConfigError,
        apply_prefits,
        build_setup,
        unwrap_config,
    )
    try:
        config = _load_config(config_file)
        if config is None:
            return None
        _, sizes_state, _ = unwrap_config(config, 'sizes')
    except ConfigError as exc:
        log.info(f"[pyirena.batch] Config file '{config_file}': {exc}")
        return None
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

    # --- Build SizesDistribution + run settings (shared with analyze) ---
    try:
        setup = build_setup(
            config, 'sizes',
            data_is_slit_smeared=bool(data.get('is_slit_smeared')),
            data_slit_length=float(data.get('slit_length', 0.0) or 0.0))
        s = setup.model
        if s.use_slit_smearing:
            log.info(f"[pyirena.batch] Sizes slit smearing enabled "
                     f"(SL={s.slit_length:.4g} 1/A).")
    except Exception:
        log.error(f"[pyirena.batch] Error building Sizes model from config:\n"
                  f"{traceback.format_exc()}")
        return None

    Q = data['Q']
    Intensity = data['Intensity']
    Error = data.get('Error')

    # --- Background pre-fits (the GUI "Fit Cmplx. Bckg. & Sizes" sequence) ---
    # Shared with analyze: full-curve windows, and an empty power-law window
    # falls back to the cursor range like the GUI (decision D4).
    try:
        for note in apply_prefits(setup, sizes_state, Q, Intensity):
            log.info(f"[pyirena.batch] {note}")
    except Exception:
        log.error(f"[pyirena.batch] Background pre-fit failed (continuing with "
                  f"configured values):\n{traceback.format_exc()}")

    # The config's own 'unc_n_runs' is the number the user set in the GUI.
    # An explicit argument still wins, so existing callers are unaffected.
    if n_mc_runs is None:
        n_mc_runs = setup.mc_n_runs or 10

    # --- Fit (+ MC): the shared Sizes runner, as the API ---
    from pyirena.core.fit_runner import run_sizes
    from pyirena.state.state_manager import text_import_settings
    fraction = text_import_settings()["error_fraction"]
    try:
        run = run_sizes(s, Q, Intensity, Error,
                        q_min=setup.fit_q_min, q_max=setup.fit_q_max,
                        error_fraction=fraction,
                        n_mc_runs=(n_mc_runs if with_uncertainty else 0),
                        mc_noise_fraction=fraction)
        fit_result = run.result
        q_fit = run.data.q
        for note in run.notes:
            log.info(f"[pyirena.batch] {note}")
    except Exception:
        log.error(f"[pyirena.batch] Fitting failed for '{data_file}':\n{traceback.format_exc()}")
        return None

    # --- MC uncertainty (computed by the runner) ---
    distribution_std = None
    if run.mc is not None:
        if run.mc.get('available'):
            distribution_std = run.mc['distribution_std']
            log.info(f"[pyirena.batch] MC uncertainty: {run.mc['n_ok']}/{run.mc['n_runs']} "
                     "refits succeeded.")
        else:
            log.error(f"[pyirena.batch] MC uncertainty not computed: {run.mc['reason']}.")

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

        # Build a complete params dict matching what the GUI saves via _collect_state().
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
            from pyirena.io.nxcansas_unified import create_nxcansas_file, get_output_filepath

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
                    # Smeared when the model is, as the GUI's
                    # _display_background (it was always the pinhole one).
                    if s.use_slit_smearing and s.slit_length > 0:
                        from pyirena.core.smearing import SlitSmearer
                        _bg = SlitSmearer(np.asarray(_q_q, dtype=float),
                                          s.slit_length).smear_model(
                                              s.compute_complex_background)
                    else:
                        _bg = s.compute_complex_background(_q_q)
                    fq_metrics = fit_quality_metrics(
                        _q_q, _I_obs, intensity_model + _bg, _err, n_params=1)
            except Exception:
                log.warning("[pyirena.batch] Fit-quality metrics not computed; saving "
                            "without them.", exc_info=True)
                fq_metrics = None

            save_sizes_results(
                filepath=output_path,
                q=q_fit,
                intensity_data=run.data.intensity,
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
