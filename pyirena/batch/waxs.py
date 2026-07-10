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
    _ensure_console()
    config_file = Path(config_file)
    config = _load_config(config_file)
    if config is None:
        msg = f"Cannot load config: {config_file}"
        log.info(f"[pyirena.batch.fit_waxs_peaks_from_config] {msg}")
        return {'success': False, 'message': msg}

    wp_cfg = config.get('waxs_peakfit')
    if wp_cfg is None:
        msg = f"No 'waxs_peakfit' section in '{config_file.name}'"
        log.info(f"[pyirena.batch.fit_waxs_peaks_from_config] {msg}")
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
        setup_state=wp_cfg,
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
        WAXSPeakFitModel, find_peaks_in_data, default_bg_params,
        cross_corr_q_shift, presearch_q0_per_peak, eval_model,
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
            h5_file = ensure_nxcansas_sibling(data_file)
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
            log.info("  [fit_waxs_peaks] Too few data points in Q range.")
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
            log.info(f"  [fit_waxs_peaks] Auto-found {len(peaks)} peak(s).")

    if not peaks:
        if verbose:
            log.info("  [fit_waxs_peaks] No peaks defined; fitting background only.")
        peaks = []

    # ── Q0 pre-search ─────────────────────────────────────────────────────────
    ps = config.get('presearch', {})
    run_cc   = bool(ps.get('cross_corr',    False))
    run_scan = bool(ps.get('per_peak_scan', False))
    if peaks and (run_cc or run_scan):
        ps_window = float(ps.get('search_window', 0.05))
        ps_steps  = int(  ps.get('n_steps',       50))
        if run_cc:
            I_model = eval_model(q_, bg_shape, bg_params, peaks, I=I_)
            shift   = cross_corr_q_shift(q_, I_, I_model, max_shift=ps_window)
            if abs(shift) > 1e-6:
                import copy as _cp
                peaks = _cp.deepcopy(peaks)
                # Only shift Q0 of peaks where Q0['fit'] is True; locked
                # Q0s stay where the user set them.
                for pk in peaks:
                    if bool(pk.get("Q0", {}).get("fit", True)):
                        pk["Q0"]["value"] = float(pk["Q0"]["value"]) + shift
                if verbose:
                    log.info(f"  [fit_waxs_peaks] Cross-corr. shift applied: {shift:+.4f} Å⁻¹")
        if run_scan:
            peaks = presearch_q0_per_peak(
                q_, I_, bg_shape, bg_params, peaks,
                search_window=ps_window, n_steps=ps_steps,
            )
            if verbose:
                log.info(f"  [fit_waxs_peaks] Per-peak Q0 scan done "
                      f"(window=±{ps_window:.3f} Å⁻¹, steps={ps_steps}).")

    # ── Fit ───────────────────────────────────────────────────────────────────
    weight_mode = config.get('weight_mode', 'standard')
    engine = WAXSPeakFitModel(bg_shape=bg_shape, peaks=peaks, no_limits=no_limits)
    try:
        result = engine.fit(q_, I_, dI_ if (dI_ is not None and np.any(np.isfinite(dI_))) else None,
                            bg_params=bg_params, peaks=peaks, weight_mode=weight_mode)
    except Exception as exc:
        if verbose:
            log.error(f"  [fit_waxs_peaks] Fit error: {exc}")
        return None

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
