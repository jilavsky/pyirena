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
