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
    n_mc_runs: int = 10,
    mc_workers: Optional[int] = None,
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
        Default 10.
    mc_workers : int, optional
        Worker processes for the MC passes, overriding the config file's
        ``mc_workers``. 0 = auto (cpu_count - 2), 1 = serial, N > 1 = explicit.
        Each pass is an independent refit, so this scales nearly linearly for
        slow models. Short runs stay serial automatically, and the passes fall
        back to serial if the host cannot start worker processes.

        Note that parallelism belongs at one level only: when driving many files
        through :func:`pyirena.batch.pipeline`, leave this at 1 rather than
        nesting pools inside a per-file loop.

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
    from pyirena.core.modeling import (
        ModelingConfig,
        ModelingEngine,
        population_from_dict,
    )
    from pyirena.io.nxcansas_modeling import save_modeling_results

    def _build_pop(pd):
        """Deserialize one population dict → the appropriate population dataclass.

        The field-by-field rebuild that used to live here (six branches, ninety
        lines, one per population type) is now
        :func:`pyirena.core.modeling.population_from_dict`, which walks the
        dataclass fields — so a field added to a population is understood by the
        batch path the moment it is declared.
        """
        pop = population_from_dict(pd)
        # Historical batch quirk, kept deliberately: a size-distribution
        # population whose config omits "enabled" has always been treated as
        # disabled here (every other type defaults to enabled).  GUI-written
        # configs always carry the key, so this only affects hand-written ones —
        # changing it would silently switch on a population in someone's
        # pipeline.
        if pd.get('pop_type', 'size_dist') == 'size_dist' and 'enabled' not in pd:
            pop.enabled = False
        return pop

    data_file = Path(data_file)
    config_file = Path(config_file)

    # --- Load config ---
    try:
        config = _load_config(config_file)
        if config is None:
            return None
        mod_cfg = config.get('modeling')
        if mod_cfg is None:
            log.info(f"[pyirena.batch.fit_modeling] No 'modeling' section in '{config_file.name}'")
            return {'success': False, 'message': "No 'modeling' section in config file"}
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
    dI = np.asarray(_err if _err is not None else np.ones_like(q) * 0.05 * I, dtype=float)
    dI = np.where(dI <= 0, 0.05 * np.abs(I), dI)

    # --- Build ModelingConfig from the dict stored by the GUI ---
    try:
        pops = [_build_pop(pd) for pd in mod_cfg.get('populations', [])]

        modeling_config = ModelingConfig(
            populations=pops,
            background=float(mod_cfg.get('background', 0.0)),
            fit_background=bool(mod_cfg.get('fit_background', True)),
            background_limits=tuple(mod_cfg.get('background_limits', [0.0, 1e10])),
            q_min=float(mod_cfg.get('q_min', float(q.min()))),
            q_max=float(mod_cfg.get('q_max', float(q.max()))),
            no_limits=bool(mod_cfg.get('no_limits', False)),
            n_mc_runs=int(mod_cfg.get('n_mc_runs', n_mc_runs)),
            fit_method=str(mod_cfg.get('fit_method', 'local')),
            de_workers=int(mod_cfg.get('de_workers', 1)),
            mc_workers=int(mc_workers if mc_workers is not None
                           else mod_cfg.get('mc_workers', 0)),
            use_slit_smearing=bool(data.get('is_slit_smeared'))
                              or bool(mod_cfg.get('use_slit_smearing', False)),
            slit_length=(float(mod_cfg['slit_length']) if mod_cfg.get('slit_length')
                         else float(data.get('slit_length', 0.0) or 0.0)),
        )
    except Exception:
        log.error(f"[pyirena.batch.fit_modeling] Error building ModelingConfig:\n{traceback.format_exc()}")
        return None
    # Guard: smearing requested but no slit length available.
    if modeling_config.use_slit_smearing and modeling_config.slit_length <= 0:
        modeling_config.use_slit_smearing = False
    if modeling_config.use_slit_smearing:
        log.info(f"[pyirena.batch.fit_modeling] Slit smearing enabled "
                 f"(SL={modeling_config.slit_length:.4g} 1/A).")

    # --- Run fit ---
    try:
        engine = ModelingEngine()
        fit_result = engine.fit(modeling_config, q, I, dI)
    except Exception:
        log.error(f"[pyirena.batch.fit_modeling] Fit error:\n{traceback.format_exc()}")
        return {'success': False, 'message': 'Modeling fit failed (see console for details)'}

    # --- Optional MC uncertainty ---
    if with_uncertainty:
        try:
            mc_runs = n_mc_runs or modeling_config.n_mc_runs
            stds = engine.calculate_uncertainty_mc(modeling_config, q, I, dI, mc_runs)
            fit_result.params_std.update(stds)
        except Exception:
            log.error(f"[pyirena.batch.fit_modeling] MC uncertainty error:\n{traceback.format_exc()}")

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
