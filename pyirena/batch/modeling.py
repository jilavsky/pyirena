"""
pyirena.batch.modeling — Headless Modeling tool fitting (fit_modeling).

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

# ---------------------------------------------------------------------------
# fit_modeling — headless Modeling (parametric size distribution) API
# ---------------------------------------------------------------------------

def fit_modeling(
    data_file: Union[str, Path],
    config_file: Union[str, Path],
    save_to_nexus: bool = True,
    with_uncertainty: bool = False,
    n_mc_runs: Optional[int] = None,
    mc_workers: Optional[int] = None,
    random_seed: Optional[int] = None,
) -> Optional[Dict]:
    """Fit a Modeling (parametric size distribution) model using a pyIrena config file.

    Reads the ``'modeling'`` section from *config_file*, builds a
    :class:`~pyirena.core.modeling.ModelingConfig`, runs
    :meth:`~pyirena.core.modeling.ModelingEngine.fit`, and optionally saves
    results to HDF5 via :func:`~pyirena.io.nxcansas_modeling.save_modeling_results`.

    Parameters
    ----------
    data_file : str or Path
        Path to SAS data file (.h5/.hdf5 NXcanSAS, or .dat/.txt text).
    config_file : str or Path
        Path to a pyIrena JSON configuration file containing a ``'modeling'``
        section (created by "Export Parameters" in the Modeling GUI panel).
    save_to_nexus : bool, optional
        If True (default), save fit results into the HDF5 file.
        Text-file inputs are not saved.
    with_uncertainty : bool, optional
        If True, run Monte Carlo uncertainty estimation after the main fit.
        Default False.
    n_mc_runs : int, optional
        Number of MC runs (used only when *with_uncertainty* is True).
        Default None takes the config's own ``n_mc_runs`` (falling back to
        10); until 1.2.1 the argument's 10 always won over the config.
    mc_workers : int, optional
        Worker processes for the MC passes, overriding the config file's
        ``mc_workers``. 0 = auto (cpu_count - 2), 1 = serial, N > 1 = explicit.
        Each pass is an independent refit, so this scales nearly linearly for
        slow models. Short runs stay serial automatically, and the passes fall
        back to serial if the host cannot start worker processes.

        Note that parallelism belongs at one level only: when driving many files
        through :func:`pyirena.batch.pipeline`, leave this at 1 rather than
        nesting pools inside a per-file loop.
    random_seed : int, optional
        Seed for the global search (config ``fit_method: "global"``) and the
        Monte-Carlo perturbations. None (default) leaves them unseeded.

    Returns
    -------
    dict or None
        On success::

            {
                'success':      bool,
                'message':      str,
                'result':       ModelingResult,
                'output_file':  Path or None,
            }

        Returns None if data loading or config reading fails fatally.
    """
    _ensure_console()
    from pyirena.io.nxcansas_modeling import save_modeling_results

    data_file = Path(data_file)
    config_file = Path(config_file)

    # --- Load config (the one envelope reader, AGENTS.md invariant 6) ---
    from pyirena.core.tool_config import ConfigError, build_setup, unwrap_config
    try:
        config = _load_config(config_file)
        if config is None:
            return None
        _, mod_cfg, _ = unwrap_config(config, 'modeling')
    except ConfigError as exc:
        log.info(f"[pyirena.batch.fit_modeling] {exc}")
        return {'success': False, 'message': str(exc)}
    except Exception:
        log.error(f"[pyirena.batch.fit_modeling] Error reading config:\n{traceback.format_exc()}")
        return None

    # --- Load data ---
    # "load_slit_smeared" picks the file's _SMR dataset; smearing is then
    # auto-enabled at the file-derived (or config-overridden) slit length.
    load_slit_smeared = bool(mod_cfg.get('load_slit_smeared', False))
    try:
        data = _load_data(data_file, load_slit_smeared=load_slit_smeared)
        if data is None:
            return None
    except Exception:
        log.error(f"[pyirena.batch.fit_modeling] Error loading data:\n{traceback.format_exc()}")
        return None

    q  = np.asarray(data['Q'],         dtype=float)
    I  = np.asarray(data['Intensity'],  dtype=float)
    _err = data.get('Error')

    # --- Build ModelingConfig: the one translator (shared with analyze) ---
    try:
        setup = build_setup(
            config, 'modeling',
            data_is_slit_smeared=bool(data.get('is_slit_smeared')),
            data_slit_length=float(data.get('slit_length', 0.0) or 0.0))
        modeling_config = setup.model
        if mc_workers is not None:
            modeling_config.mc_workers = int(mc_workers)
    except Exception:
        log.error(f"[pyirena.batch.fit_modeling] Error building ModelingConfig:\n"
                  f"{traceback.format_exc()}")
        return None
    if modeling_config.use_slit_smearing:
        log.info(f"[pyirena.batch.fit_modeling] Slit smearing enabled "
                 f"(SL={modeling_config.slit_length:.4g} 1/A).")

    # The config's own n_mc_runs is what the user set; an explicit argument wins.
    if n_mc_runs is None:
        n_mc_runs = int(modeling_config.n_mc_runs or 10)

    # --- Fit (+ MC): the shared Modeling runner, as the API ---
    # No measured uncertainties: σ = I × the Data Selector's "Generated
    # uncertainty fraction" (the user's choice; issue #38).
    from pyirena.core.fit_runner import run_modeling
    from pyirena.state.state_manager import text_import_settings
    frac = text_import_settings()["error_fraction"]
    try:
        run = run_modeling(modeling_config, q, I, _err,
                           q_min=setup.fit_q_min, q_max=setup.fit_q_max,
                           error_fraction=frac,
                           n_mc_runs=(n_mc_runs if with_uncertainty else 0),
                           mc_noise_fraction=frac,
                           mc_workers=modeling_config.mc_workers,
                           seed=random_seed)
        fit_result = run.result
        for w in (getattr(fit_result, 'fit_warnings', None) or []) + run.notes:
            log.warning(f"[pyirena.batch.fit_modeling] {w}")
    except Exception:
        log.error(f"[pyirena.batch.fit_modeling] Fit error:\n{traceback.format_exc()}")
        return {'success': False, 'message': 'Modeling fit failed (see console for details)'}
    dI = run.data.error if run.data.error is not None else np.abs(run.data.intensity) * frac
    q, I = run.data.q, run.data.intensity

    if run.mc is not None:
        if run.mc.get('available'):
            fit_result.params_std.update(run.mc['stds'])
        else:
            log.error(f"[pyirena.batch.fit_modeling] MC uncertainty not computed: "
                      f"{run.mc['reason']}.")

    # --- Save to HDF5 ---
    out_path = None
    if save_to_nexus and data_file.suffix.lower() in ('.h5', '.hdf5', '.nxs'):
        try:
            # Robust fit-quality metrics over the fitted Q range.
            fq_metrics = None
            try:
                from pyirena.core.fit_metrics import fit_quality_metrics
                _mask = (q >= fit_result.config.q_min) & (q <= fit_result.config.q_max)
                if np.any(_mask) and fit_result.model_I is not None:
                    _n_free = max(1, int(np.count_nonzero(_mask)) - int(fit_result.dof))
                    fq_metrics = fit_quality_metrics(
                        q[_mask], I[_mask], fit_result.model_I,
                        dI[_mask] if dI is not None else None, n_params=_n_free)
            except Exception:
                log.warning("[pyirena.batch.fit_modeling] Fit-quality metrics not "
                            "computed; saving without them.", exc_info=True)
                fq_metrics = None
            save_modeling_results(data_file, fit_result, fit_quality=fq_metrics,
                                  setup_state=mod_cfg)
            out_path = data_file
        except Exception:
            log.error(f"[pyirena.batch.fit_modeling] Save error:\n{traceback.format_exc()}")

    chi2_str = f"χ²/dof={fit_result.reduced_chi_squared:.4g}"
    msg = f"Modeling fit complete — {chi2_str}, {len(fit_result.pop_indices)} active population(s)"
    log.info(f"[pyirena.batch] {msg}")

    return {
        'success': True,
        'message': msg,
        'result': fit_result,
        'output_file': out_path,
    }
