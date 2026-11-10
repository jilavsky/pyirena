"""
pyirena.batch.waxs — Headless WAXS peak fitting (fit_waxs_peaks, fit_waxs_peaks_from_config).

Split from the original monolithic pyirena/batch.py (no behavior change).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional, Union

import numpy as np

from pyirena.logging_setup import ensure_console_output as _ensure_console

if TYPE_CHECKING:
    pass

log = logging.getLogger(__name__)

from pyirena.batch._common import _load_config


def fit_waxs_peaks_from_config(
    data_file: Union[str, Path],
    config_file: Union[str, Path],
    save_to_nexus: bool = True,
    with_uncertainty: bool = False,
    n_mc_runs: int = 10,
) -> Optional[Dict]:
    """Fit WAXS diffraction peaks using parameters from a pyIrena config file.

    Reads the ``'waxs_peakfit'`` section from a JSON config file and calls
    :func:`fit_waxs_peaks`.

    **Uncertainties.** Every WAXS fit returns the ``curve_fit`` covariance
    1σ of each parameter (``peaks_std``, ``bg_params_std``, and the peak
    areas' σ); there is no Monte-Carlo estimate for WAXS. ``with_uncertainty``
    and ``n_mc_runs`` are in the signature because
    :func:`pyirena.batch.pipeline.fit_pyirena` passes them to every tool; asking for uncertainty here is
    logged and the result says which estimate it carries
    (``uncertainty_method: 'covariance'``).

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
        No effect on the numbers: the covariance σ are always returned. When
        True, that is logged so a caller expecting Monte-Carlo is told.
    n_mc_runs : int, optional
        Not used — WAXS has no Monte-Carlo uncertainty.

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
    _ensure_console()
    config_file = Path(config_file)
    config = _load_config(config_file)
    if config is None:
        msg = f"Cannot load config: {config_file}"
        log.info(f"[pyirena.batch.fit_waxs_peaks_from_config] {msg}")
        return {'success': False, 'message': msg}

    # The one envelope reader (AGENTS.md invariant 6): every wrapper pyIrena
    # writes, and the envelope's fit_q_range, as analyze() reads them.
    from pyirena.core.tool_config import ConfigError, build_setup, unwrap_config
    try:
        _, wp_cfg, _ = unwrap_config(config, 'waxs_peakfit')
        setup = build_setup(config, 'waxs_peakfit')
    except ConfigError as exc:
        msg = f"'{config_file.name}': {exc}"
        log.info(f"[pyirena.batch.fit_waxs_peaks_from_config] {msg}")
        return {'success': False, 'message': msg}

    q_min = setup.fit_q_min
    q_max = setup.fit_q_max

    result = fit_waxs_peaks(
        data_file=data_file,
        config=wp_cfg,
        q_min=q_min,
        q_max=q_max,
        save_to_nexus=save_to_nexus,
        verbose=True,
        setup_state=wp_cfg,
    )

    if result is None:
        return {'success': False, 'message': 'Data load or fit failure (fit_waxs_peaks returned None)'}

    if not result.get('success'):
        result.setdefault('message', result.get('error', 'fit failed'))
    else:
        result['uncertainty_method'] = 'covariance'
        if with_uncertainty:
            log.info("[pyirena.batch.fit_waxs] WAXS uncertainties are the curve_fit "
                     "covariance values already in the result; n_mc_runs does not apply.")

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
    setup_state: Optional[Dict] = None,
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
    _ensure_console()
    from pyirena.core.waxs_peakfit import (
        WAXSPeakFitModel,
    )
    from pyirena.io.hdf5 import readGenericNXcanSAS

    data_file = Path(data_file)
    if verbose:
        log.info(f"[pyirena.batch.fit_waxs_peaks] {data_file.name}")

    # ── Load data ────────────────────────────────────────────────────────────
    try:
        ext = data_file.suffix.lower()
        if ext in ('.txt', '.dat'):
            from pyirena.io.text_import import ensure_nxcansas_sibling
            from pyirena.state.state_manager import text_import_settings
            h5_file = ensure_nxcansas_sibling(data_file, **text_import_settings())
            data = readGenericNXcanSAS(str(h5_file.parent), h5_file.name)
        else:
            data = readGenericNXcanSAS(str(data_file.parent), data_file.name)
        if data is None:
            if verbose:
                log.error(f"  [fit_waxs_peaks] Could not load data from {data_file.name}")
            return None
    except Exception as exc:
        if verbose:
            log.error(f"  [fit_waxs_peaks] Load error: {exc}")
        return None

    q   = np.asarray(data['Q'],        float)
    I   = np.asarray(data['Intensity'], float)
    dI  = data.get('Error')
    dI  = np.asarray(dI, float) if dI is not None else None     # no Idev: fine

    # ── Fit: the shared WAXS runner, as the API ──────────────────────────────
    # Prunes peaks outside the range, auto-finds peaks when none are given
    # and the config has a peak_find block, runs the Q0 pre-search, applies
    # the shared point rule (it crashed on files without uncertainties).
    from pyirena.core.fit_runner import run_waxs
    from pyirena.state.state_manager import text_import_settings
    model = WAXSPeakFitModel.from_dict(config)
    try:
        run = run_waxs(model, q, I, dI, q_min=q_min, q_max=q_max,
                       error_fraction=text_import_settings()["error_fraction"],
                       weight_mode=config.get('weight_mode', 'standard'),
                       presearch=config.get('presearch') or None,
                       peak_find=config.get('peak_find') or None)
    except Exception as exc:
        if verbose:
            log.error(f"  [fit_waxs_peaks] Fit error: {exc}")
        return None
    if len(run.data.q) < 5:
        if verbose:
            log.info("  [fit_waxs_peaks] Too few data points in Q range.")
        return None
    if verbose:
        for note in run.notes:
            log.info(f"  [fit_waxs_peaks] {note}")
    result = run.result
    q_, I_, dI_ = run.data.q, run.data.intensity, run.data.error

    if verbose:
        status = "OK" if result.get('success') else "FAILED"
        log.info(f"  [fit_waxs_peaks] {status}  "
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
                setup_state=setup_state,
            )
            if verbose:
                log.info(f"  [fit_waxs_peaks] Saved to {data_file.name}.")
        except Exception as exc:
            if verbose:
                log.error(f"  [fit_waxs_peaks] Save error: {exc}")

    return result
