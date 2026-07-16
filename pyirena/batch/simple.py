"""
pyirena.batch.simple — Headless simple-model fitting (fit_simple, fit_simple_from_config).

Split from the original monolithic pyirena/batch.py (no behavior change).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional, Union

import numpy as np

from pyirena.logging_setup import ensure_console_output as _ensure_console

if TYPE_CHECKING:
    from pyirena.core.simple_fits import SimpleFitModel

log = logging.getLogger(__name__)

from pyirena.batch._common import _load_config, _load_data


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
    fixed_params: Optional[Dict] = None,
    verbose: bool = True,
    setup_state: Optional[Dict] = None,
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
    _ensure_console()
    from pyirena.core.simple_fits import SimpleFitModel
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
        log.info(f"[pyirena.batch.fit_simple] config must be a dict or SimpleFitModel, "
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

    # ── Background prefit replay (Invariant) ────────────────────────────────
    # Re-determine the complex background from the Q ranges the GUI user
    # recorded (bg_prefit), using the FULL data — the background windows
    # usually lie outside the integration range applied below.
    if model.is_calculation and (model.bg_prefit or {}).get('enabled'):
        applied = model.prefit_background(q, I)
        if verbose and applied:
            vals = '  '.join(f'{k}={v:.4g}' for k, v in applied.items()
                             if k != 'warning')
            log.info(f"[pyirena.batch] Background prefit replayed: {vals}")
        if applied.get('warning'):
            log.warning(f"[pyirena.batch] Background prefit: {applied['warning']}")

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
        log.info(f"[pyirena.batch.fit_simple] Too few data points after Q masking "
              f"({mask.sum()} points, need ≥ 2).")
        return None

    qf = q[mask]
    If = I[mask]
    dIf = dI[mask] if dI is not None else None

    # ── Initial fit ──────────────────────────────────────────────────────────
    if verbose:
        log.info(f"[pyirena.batch] Fitting {model.model} to '{data_file.name}' ...")

    # Honor the GUI's "Fit?" checkboxes: parameters marked fixed are held
    # constant during fitting, matching SimpleFitsPanel._run_fit().
    fixed_params = fixed_params if fixed_params else None
    result = model.fit(qf, If, dIf, fixed_params=fixed_params)

    if not result.get('success'):
        msg = result.get('error', 'Unknown error')
        log.error(f"[pyirena.batch.fit_simple] Fit failed: {msg}")
        return result

    if verbose:
        rchi2 = result.get('reduced_chi2')
        if rchi2 is not None:
            log.info(f"[pyirena.batch] Fit succeeded.  Reduced χ² = {rchi2:.4g}")
        else:
            derived = result.get('derived', {})
            log.info(f"[pyirena.batch] Calculation done.  "
                     + '  '.join(f'{k}={v:.4g}' for k, v in derived.items()))

    # ── Monte Carlo uncertainty ──────────────────────────────────────────────
    if with_uncertainty and model.is_calculation:
        log.info("[pyirena.batch] MC uncertainty is not applicable to "
                 f"calculation model '{model.model}' — skipped.")
        with_uncertainty = False
    if with_uncertainty:
        if verbose:
            log.info(f"[pyirena.batch] Running {n_mc_runs} MC uncertainty runs ...")

        dIf_safe = (dIf if dIf is not None
                    else np.maximum(If * 0.05, 1e-30))
        mc_params: Dict[str, list] = {k: [] for k in result['params']}

        for _ in range(n_mc_runs):
            I_perturbed = If + dIf_safe * np.random.randn(len(If))
            mc_res = model.fit(qf, I_perturbed, dIf, fixed_params=fixed_params)
            if mc_res.get('success'):
                for k in mc_params:
                    mc_params[k].append(mc_res['params'].get(k, float('nan')))

        if any(len(v) > 1 for v in mc_params.values()):
            result['params_std'] = {
                k: float(np.std(v, ddof=1)) if len(v) > 1 else float('nan')
                for k, v in mc_params.items()
            }
            if verbose:
                log.info(f"[pyirena.batch] MC uncertainty done "
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
                setup_state=setup_state,
            )
            if verbose:
                log.info(f"[pyirena.batch] Results saved to '{data_file.name}'")
        except Exception as exc:
            log.error(f"[pyirena.batch.fit_simple] Could not save results: {exc}")

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
    _ensure_console()
    config_file = Path(config_file)
    config = _load_config(config_file)
    if config is None:
        log.error(f"[pyirena.batch.fit_simple_from_config] Cannot load config: {config_file}")
        return None

    sf_cfg = config.get('simple_fits')
    if sf_cfg is None:
        log.info(f"[pyirena.batch.fit_simple_from_config] No 'simple_fits' section in "
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

    # The GUI stores per-parameter "Fit?" state as param_fixed = {name: True if
    # held fixed}.  SimpleFitModel.fit() expects fixed_params = {name: value},
    # so build that from the params dict.  Without this the batch path would
    # refit every parameter, ignoring the user's fixed-parameter choices.
    param_fixed = sf_cfg.pop('param_fixed', {}) or {}
    params = sf_cfg.get('params', {}) or {}
    fixed_params = {
        name: params[name]
        for name, is_fixed in param_fixed.items()
        if is_fixed and name in params
    }

    # Remove state-only keys that have no meaning for from_dict()
    for _k in ('schema_version', 'no_limits'):
        sf_cfg.pop(_k, None)

    result = fit_simple(
        data_file=data_file,
        config=sf_cfg,
        with_uncertainty=with_uncertainty,
        n_mc_runs=n_mc_runs,
        q_min=q_min,
        q_max=q_max,
        fixed_params=fixed_params if fixed_params else None,
        verbose=True,
        setup_state=config.get('simple_fits'),
    )

    if result is None:
        return {'success': False, 'message': 'fit_simple returned None (data load failure?)'}

    # Normalise error key so BatchWorker can always call result.get('message')
    if not result.get('success'):
        result.setdefault('message', result.get('error', 'fit failed'))

    return result
