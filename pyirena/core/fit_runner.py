"""
One fit runner per tool, shared by ``pyirena.batch`` and ``pyirena.api.control``.

Review of 09-11-2026 §3.2 / issue #38. Fit orchestration — which points, which
σ, Q range, how the fit is called, Monte-Carlo — used to be written separately
in the batch module, the control module (and so ``analyze``) and the GUI
panel of each tool, and an audit found ~35 places where the copies disagreed.
Each ``run_<tool>`` here is now the one place that decides; the decisions are
recorded in ``planning/fit-runner/README.md``. The GUI still has its own
orchestration (a later step) but follows the same data rule
(:mod:`pyirena.core.fit_data`).

A runner takes a ready model, the full data arrays and the run settings, and
returns a :class:`FitRun`. Loading files, building models from configs
(``core/tool_config.py``), saving and formatting replies stay with the
callers.

numpy/scipy only (through the core models).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np

from pyirena.core.fit_data import FitArrays, prepare_fit_arrays, q_range_mask


@dataclass
class FitRun:
    """What a runner did and got."""

    tool: str
    #: The fitted model (updated in place by the fit).
    model: Any
    #: The core fit's own result object (dict or dataclass, per tool).
    result: Any
    #: The points and σ the fit used (Q range applied).
    data: FitArrays
    #: Monte-Carlo result dict from the core, or None when not requested.
    mc: Optional[dict] = None
    #: Plain-language notes for the caller (repaired σ, dropped points, …).
    notes: list = field(default_factory=list)

    @property
    def success(self) -> bool:
        r = self.result
        if isinstance(r, dict):
            return bool(r.get("success", False))
        return bool(getattr(r, "success", False))


def _mc_noise(data: FitArrays, mc_noise_fraction: Optional[float]):
    """σ for Monte-Carlo when the data have no uncertainties (else None)."""
    if data.error is not None or mc_noise_fraction is None:
        return None
    return np.abs(data.intensity) * float(mc_noise_fraction)


# ---------------------------------------------------------------------------
# Unified Fit
# ---------------------------------------------------------------------------

def run_unified(model, q, intensity, error=None, *, q_min=None, q_max=None,
                error_fraction: float = 0.05, no_limits: bool = False,
                walk_limits: bool = True, max_iterations: Optional[int] = None,
                n_mc_runs: int = 0, seed=None,
                mc_noise_fraction: Optional[float] = None,
                progress_cb=None, cancel_cb=None) -> FitRun:
    """Fit a :class:`~pyirena.core.unified.UnifiedFitModel`.

    Args:
        model: The model to fit (levels, flags, limits, links, slit settings
            already set); updated in place.
        q, intensity, error: Full data arrays; *error* may be None.
        q_min, q_max: Inclusive fit range; None = open.
        error_fraction: σ for points whose σ ≤ 0 (D2).
        no_limits: Fit with the model's wide default bounds
            (``fit(no_limits=True)``); implies no limit walking.
        walk_limits: Re-centre limits a fitted parameter ends pinned at and
            refit (D1, the GUI's behaviour). Ignored with *no_limits*.
        max_iterations: Optional override of the optimiser budget.
        n_mc_runs: Monte-Carlo refits after a successful fit (0 = none).
        seed: Seed for the Monte-Carlo perturbations.
        mc_noise_fraction: When the data have no uncertainties, perturb with
            ``|I| ×`` this (the user's generated fraction); None = refuse.
        progress_cb, cancel_cb: Passed to the Monte-Carlo.
    """
    data = prepare_fit_arrays(q, intensity, error, error_fraction,
                              base_mask=q_range_mask(q, q_min, q_max))
    kwargs = {}
    if max_iterations is not None:
        kwargs["max_iterations"] = int(max_iterations)
    if no_limits:
        result = model.fit(data.q, data.intensity, data.error, no_limits=True, **kwargs)
    elif walk_limits:
        result = model.fit_with_limit_walking(data.q, data.intensity, data.error, **kwargs)
    else:
        result = model.fit(data.q, data.intensity, data.error, **kwargs)

    run = FitRun("unified_fit", model, result, data, notes=data.notes())
    if n_mc_runs and run.success:
        mc_kwargs = dict(kwargs)
        if no_limits:
            mc_kwargs["no_limits"] = True
        run.mc = model.calculate_uncertainty_mc(
            data.q, data.intensity, data.error, n_runs=n_mc_runs, seed=seed,
            noise_sigma=_mc_noise(data, mc_noise_fraction),
            progress_cb=progress_cb, cancel_cb=cancel_cb, **mc_kwargs)
    return run


# ---------------------------------------------------------------------------
# Size Distribution
# ---------------------------------------------------------------------------

def run_sizes(model, q, intensity, error=None, *, q_min=None, q_max=None,
              error_fraction: float = 0.05, seed=None, n_mc_runs: int = 0,
              mc_noise_fraction: Optional[float] = None,
              progress_cb=None, cancel_cb=None) -> FitRun:
    """Invert with a :class:`~pyirena.core.sizes.SizesDistribution`.

    Background pre-fits are the caller's (``core.tool_config.apply_prefits``
    on a replay); the model arrives with its background terms set. *seed*
    seeds the stochastic Monte Carlo inversion and the uncertainty refits.
    Points with σ ≤ 0 are repaired rather than dropped by the core (D2); the
    core's own I > 0 mask is unchanged.
    """
    data = prepare_fit_arrays(q, intensity, error, error_fraction,
                              base_mask=q_range_mask(q, q_min, q_max))
    result = model.fit(data.q, data.intensity, data.error, rng=seed)
    run = FitRun("sizes", model, result, data, notes=data.notes())
    if n_mc_runs and run.success:
        run.mc = model.calculate_uncertainty_mc(
            data.q, data.intensity, data.error, n_runs=n_mc_runs, seed=seed,
            noise_sigma=_mc_noise(data, mc_noise_fraction),
            progress_cb=progress_cb, cancel_cb=cancel_cb)
    return run


# ---------------------------------------------------------------------------
# Simple Fits
# ---------------------------------------------------------------------------

def run_simple(model, q, intensity, error=None, *, q_min=None, q_max=None,
               fixed_params: Optional[dict] = None, no_limits: bool = False,
               error_fraction: float = 0.05, n_mc_runs: int = 0, seed=None,
               mc_noise_fraction: Optional[float] = None,
               progress_cb=None, cancel_cb=None) -> FitRun:
    """Fit (or, for the Invariant, calculate) a
    :class:`~pyirena.core.simple_fits.SimpleFitModel`.

    The Invariant's complex-background pre-fit is replayed here on the FULL
    curve when the model records one (``bg_prefit``) — every path now does
    this; the API session path used to skip it. Points follow the shared rule
    (I ≤ 0 kept, σ ≤ 0 repaired); batch used to drop σ ≤ 0 points and the
    GUI I ≤ 0 points. Monte-Carlo does not apply to the Invariant.
    """
    from pyirena.core.tool_config import apply_simple_bg_prefit

    notes = list(apply_simple_bg_prefit(model, fixed_params, q, intensity))
    data = prepare_fit_arrays(q, intensity, error, error_fraction,
                              base_mask=q_range_mask(q, q_min, q_max))
    result = model.fit(data.q, data.intensity, data.error,
                       fixed_params=fixed_params or None, no_limits=bool(no_limits))
    run = FitRun("simple_fits", model, result, data, notes=notes + data.notes())
    if n_mc_runs and run.success and not model.is_calculation:
        run.mc = model.calculate_uncertainty_mc(
            data.q, data.intensity, data.error, n_runs=n_mc_runs,
            fixed_params=fixed_params or None, no_limits=bool(no_limits), seed=seed,
            noise_sigma=_mc_noise(data, mc_noise_fraction),
            progress_cb=progress_cb, cancel_cb=cancel_cb)
    return run


# ---------------------------------------------------------------------------
# Modeling
# ---------------------------------------------------------------------------

def run_modeling(config, q, intensity, error=None, *, q_min=None, q_max=None,
                 error_fraction: float = 0.05, no_limits: Optional[bool] = None,
                 n_mc_runs: int = 0, seed=None,
                 mc_noise_fraction: Optional[float] = None, mc_workers=None,
                 progress_cb=None, cancel_cb=None) -> FitRun:
    """Fit a :class:`~pyirena.core.modeling.ModelingConfig` (updated in place).

    Q range: *q_min*/*q_max*, a missing bound being the data's own (written
    onto the config, which the engine crops by). σ: the shared rule for
    σ ≤ 0; with no uncertainties at all, ``|I| × error_fraction`` — the fit
    needs a σ, and this is the user's generated fraction (GUI/batch) or the
    session's (API). *no_limits*, when given, applies to this fit only.
    Monte-Carlo starts from the fitted config; data without uncertainties
    are refused unless *mc_noise_fraction* is given. The core reports σ but
    not how many refits succeeded, so ``mc['n_ok']`` is None.
    """
    from pyirena.core import mc_uncertainty as mcu
    from pyirena.core.modeling import ModelingEngine

    q_arr = np.asarray(q, dtype=float)
    finite = q_arr[np.isfinite(q_arr) & (q_arr > 0)]
    lo = float(q_min) if q_min is not None else float(finite.min())
    hi = float(q_max) if q_max is not None else float(finite.max())
    config.q_min, config.q_max = lo, hi
    data = prepare_fit_arrays(q, intensity, error, error_fraction,
                              base_mask=q_range_mask(q, lo, hi))
    dI = data.error if data.error is not None else \
        np.abs(data.intensity) * float(error_fraction)

    engine = ModelingEngine()
    previous = config.no_limits
    if no_limits is not None:
        config.no_limits = bool(no_limits)
    try:
        result = engine.fit(config, data.q, data.intensity, dI)
    finally:
        config.no_limits = previous          # a property of this fit only
    result.config.no_limits = previous
    run = FitRun("modeling", result.config, result, data, notes=data.notes())

    if n_mc_runs:
        if data.error is None and mc_noise_fraction is None:
            run.mc = mcu.not_available(mcu.check_sigma(None, len(data.q)), n_mc_runs, 0, seed)
            run.mc["n_ok"] = None
        else:
            dI_mc = dI if data.error is not None else \
                np.abs(data.intensity) * float(mc_noise_fraction)
            cfg = result.config
            prev = cfg.no_limits
            if no_limits is not None:
                cfg.no_limits = bool(no_limits)
            try:
                stds = engine.calculate_uncertainty_mc(
                    cfg, data.q, data.intensity, dI_mc, n_runs=n_mc_runs, seed=seed,
                    workers=mc_workers, progress_cb=progress_cb, cancel_cb=cancel_cb)
            finally:
                cfg.no_limits = prev
            if stds:
                run.mc = mcu.base_result(n_mc_runs, 0, seed)
                run.mc["stds"] = stds
            else:
                run.mc = mcu.not_available("fewer than 2 refits succeeded",
                                           n_mc_runs, 0, seed)
            run.mc["n_ok"] = None
            run.mc["sigma_source"] = "measured" if data.error is not None else "generated"
    return run


# ---------------------------------------------------------------------------
# WAXS Peak Fit
# ---------------------------------------------------------------------------

def run_waxs(model, q, intensity, error=None, *, q_min=None, q_max=None,
             error_fraction: float = 0.05, weight_mode: str = "standard",
             no_limits: Optional[bool] = None, presearch: Optional[dict] = None,
             peak_find: Optional[dict] = None) -> FitRun:
    """Fit a :class:`~pyirena.core.waxs_peakfit.WAXSPeakFitModel` (updated in place).

    In order, as the GUI does:

    1. points: the shared rule within the Q range (files without
       uncertainties are fine — batch crashed on them);
    2. peaks whose Q0 lies outside the fitted range are removed (only the GUI
       did this; a peak outside the window is unconstrained but still bends
       the model) — listed in the notes;
    3. with no peaks and a ``peak_find`` block, peaks are found in the range;
    4. the Q0 pre-search (``presearch`` block: cross-correlation shift and/or
       per-peak scan) on the fitted range — the API session path never ran it;
    5. the fit, with *weight_mode*; *no_limits* None means the model's own.
    """
    from pyirena.core.tool_config import apply_waxs_presearch

    data = prepare_fit_arrays(q, intensity, error, error_fraction,
                              base_mask=q_range_mask(q, q_min, q_max))
    notes = data.notes()
    if len(data.q):
        lo, hi = float(data.q.min()), float(data.q.max())
        kept, dropped = [], []
        for pk in model.peaks or []:
            v = pk.get("Q0", {})
            q0 = float(v.get("value", 0.0) if isinstance(v, dict) else v)
            (kept if lo <= q0 <= hi else dropped).append(pk)
        if dropped:
            model.peaks = kept
            notes.append(f"Removed {len(dropped)} peak(s) with Q0 outside the fitted "
                         f"range [{lo:.4g}, {hi:.4g}] 1/A.")
    if not model.peaks and peak_find:
        from pyirena.core.waxs_peakfit import find_peaks_in_data
        model.peaks = find_peaks_in_data(
            data.q, data.intensity,
            prominence_frac=peak_find.get('prominence', 0.05),
            min_fwhm=peak_find.get('min_fwhm', 0.001),
            max_fwhm=peak_find.get('max_fwhm', 0.5),
            min_distance=peak_find.get('min_distance', 0.005),
            sg_window_frac=peak_find.get('sg_window_frac', 0.15))
        notes.append(f"Auto-found {len(model.peaks)} peak(s).")
    if presearch:
        notes.extend(apply_waxs_presearch(model, {"presearch": presearch},
                                          data.q, data.intensity))

    previous = model.no_limits
    if no_limits is not None:
        model.no_limits = bool(no_limits)
    try:
        result = model.fit(data.q, data.intensity, data.error,
                           model.bg_params, model.peaks, weight_mode=weight_mode)
    finally:
        model.no_limits = previous              # a property of this fit only
    if result.get("success", False):
        model.bg_params = result.get("bg_params", model.bg_params)
        model.peaks = result.get("peaks", model.peaks)
    return FitRun("waxs_peakfit", model, result, data, notes=notes)


# ---------------------------------------------------------------------------
# Carbon model
# ---------------------------------------------------------------------------

def run_carbon(model, q, intensity, error=None, *, q_min=None, q_max=None,
               error_fraction: float = 0.05, weighting: Optional[str] = None,
               n_mc_runs: int = 0, no_limits: bool = False,
               progress=None) -> FitRun:
    """Fit a :class:`~pyirena.core.carbon_fit.CarbonFitModel` (updated in place).

    *q_min*/*q_max* None = open (the model's ``0`` convention). *weighting*
    None keeps the model's own (API sessions used to overwrite it with
    ``'auto'``). *n_mc_runs* applies to this fit only and is restored
    afterwards, as the GUI does — so a config's stored count no longer turns
    MC on by itself in batch. Uncertainties (covariance, or MC) come back in
    ``result.errors``; ``FitRun.mc`` stays None. *progress* is the core's
    ``progress(iteration, chi2)`` callback (raise ``CarbonFitAborted`` from
    it to stop).
    """
    model.q_min = float(q_min) if q_min else 0.0
    model.q_max = float(q_max) if q_max else 0.0
    data = prepare_fit_arrays(q, intensity, error, error_fraction,
                              base_mask=q_range_mask(q, q_min or None, q_max or None))
    if weighting is not None:
        model.weighting = weighting
    previous_mc = model.n_mc_runs
    model.n_mc_runs = int(n_mc_runs or 0)
    try:
        result = model.fit(data.q, data.intensity, data.error,
                           progress=progress, no_limits=bool(no_limits))
    finally:
        model.n_mc_runs = previous_mc
    return FitRun("carbon_fit", model, result, data, notes=data.notes())
