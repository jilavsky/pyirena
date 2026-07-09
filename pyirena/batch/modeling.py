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
        ModelingEngine, ModelingConfig,
        SizeDistPopulation, UnifiedLevelPopulation, DiffractionPeakPopulation,
        GuinierPorodPopulation, MassFractalPopulation, SurfaceFractalPopulation,
    )
    from pyirena.io.nxcansas_modeling import save_modeling_results

    def _build_pop(pd):
        """Deserialize one population dict → the appropriate population dataclass."""
        pt = pd.get('pop_type', 'size_dist')
        if pt == 'unified_level':
            pop = UnifiedLevelPopulation()
            pop.enabled = bool(pd.get('enabled', True))
            pop.label = pd.get('label', '')
            for key in ['G', 'Rg', 'B', 'P', 'RgCO']:
                setattr(pop, key, float(pd.get(key, getattr(pop, key))))
                setattr(pop, f'fit_{key}', bool(pd.get(f'fit_{key}', getattr(pop, f'fit_{key}'))))
                lim = pd.get(f'{key}_limits', list(getattr(pop, f'{key}_limits')))
                setattr(pop, f'{key}_limits', tuple(lim))
            pop.correlations = bool(pd.get('correlations', False))
            for key in ['ETA', 'PACK']:
                setattr(pop, key, float(pd.get(key, getattr(pop, key))))
                setattr(pop, f'fit_{key}', bool(pd.get(f'fit_{key}', getattr(pop, f'fit_{key}'))))
                lim = pd.get(f'{key}_limits', list(getattr(pop, f'{key}_limits')))
                setattr(pop, f'{key}_limits', tuple(lim))
            return pop
        if pt == 'diffraction_peak':
            pop = DiffractionPeakPopulation()
            pop.enabled = bool(pd.get('enabled', True))
            pop.label = pd.get('label', '')
            pop.peak_type = pd.get('peak_type', 'gaussian')
            for key in ['position', 'amplitude', 'width', 'eta_voigt']:
                setattr(pop, key, float(pd.get(key, getattr(pop, key))))
                setattr(pop, f'fit_{key}', bool(pd.get(f'fit_{key}', getattr(pop, f'fit_{key}'))))
                lim = pd.get(f'{key}_limits', list(getattr(pop, f'{key}_limits')))
                setattr(pop, f'{key}_limits', tuple(lim))
            return pop
        if pt == 'guinier_porod':
            pop = GuinierPorodPopulation()
            pop.enabled = bool(pd.get('enabled', True))
            pop.label = pd.get('label', '')
            for key in ['G', 'Rg1', 's1', 'P', 'Rg2', 's2', 'RgCO']:
                setattr(pop, key, float(pd.get(key, getattr(pop, key))))
                setattr(pop, f'fit_{key}', bool(pd.get(f'fit_{key}', getattr(pop, f'fit_{key}'))))
                lim = pd.get(f'{key}_limits', list(getattr(pop, f'{key}_limits')))
                setattr(pop, f'{key}_limits', tuple(lim))
            pop.correlations = bool(pd.get('correlations', False))
            for key in ['ETA', 'PACK']:
                setattr(pop, key, float(pd.get(key, getattr(pop, key))))
                setattr(pop, f'fit_{key}', bool(pd.get(f'fit_{key}', getattr(pop, f'fit_{key}'))))
                lim = pd.get(f'{key}_limits', list(getattr(pop, f'{key}_limits')))
                setattr(pop, f'{key}_limits', tuple(lim))
            return pop
        if pt == 'mass_fractal':
            pop = MassFractalPopulation()
            pop.enabled = bool(pd.get('enabled', True))
            pop.label = pd.get('label', '')
            for key in ['Phi', 'Radius', 'Beta', 'Dv', 'Ksi', 'Eta', 'Contrast']:
                setattr(pop, key, float(pd.get(key, getattr(pop, key))))
                setattr(pop, f'fit_{key}', bool(pd.get(f'fit_{key}', getattr(pop, f'fit_{key}'))))
                lim = pd.get(f'{key}_limits', list(getattr(pop, f'{key}_limits')))
                setattr(pop, f'{key}_limits', tuple(lim))
            return pop
        if pt == 'surface_fractal':
            pop = SurfaceFractalPopulation()
            pop.enabled = bool(pd.get('enabled', True))
            pop.label = pd.get('label', '')
            for key in ['Surface', 'Ds', 'Ksi', 'Contrast', 'Qc', 'QcWidth']:
                setattr(pop, key, float(pd.get(key, getattr(pop, key))))
                setattr(pop, f'fit_{key}', bool(pd.get(f'fit_{key}', getattr(pop, f'fit_{key}'))))
                lim = pd.get(f'{key}_limits', list(getattr(pop, f'{key}_limits')))
                setattr(pop, f'{key}_limits', tuple(lim))
            pop.use_porod_transition = bool(pd.get('use_porod_transition', False))
            return pop
        # default: size_dist
        return SizeDistPopulation(
            enabled=pd.get('enabled', False),
            dist_type=pd.get('dist_type', 'lognormal'),
            dist_params=dict(pd.get('dist_params', {})),
            dist_params_fit=dict(pd.get('dist_params_fit', {})),
            dist_params_limits={k: tuple(v) for k, v in pd.get('dist_params_limits', {}).items()},
            form_factor=pd.get('form_factor', 'sphere'),
            ff_params=dict(pd.get('ff_params', {})),
            ff_params_fit=dict(pd.get('ff_params_fit', {})),
            ff_params_limits={k: tuple(v) for k, v in pd.get('ff_params_limits', {}).items()},
            structure_factor=pd.get('structure_factor', 'none'),
            sf_params=dict(pd.get('sf_params', {})),
            sf_params_fit=dict(pd.get('sf_params_fit', {})),
            sf_params_limits={k: tuple(v) for k, v in pd.get('sf_params_limits', {}).items()},
            contrast=float(pd.get('contrast', 1.0)),
            fit_contrast=bool(pd.get('fit_contrast', False)),
            contrast_limits=tuple(pd.get('contrast_limits', [0.0, 1e10])),
            scale=float(pd.get('scale', 0.001)),
            fit_scale=bool(pd.get('fit_scale', True)),
            scale_limits=tuple(pd.get('scale_limits', [0.0, 1.0])),
            use_number_dist=bool(pd.get('use_number_dist', False)),
            n_bins=int(pd.get('n_bins', 200)),
        )

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
    try:
        data = _load_data(data_file)
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
        )
    except Exception:
        log.error(f"[pyirena.batch.fit_modeling] Error building ModelingConfig:\n{traceback.format_exc()}")
        return None

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
