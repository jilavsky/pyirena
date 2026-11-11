"""
Shared Monte-Carlo uncertainty plumbing for the fitting tools.

Every tool estimates parameter uncertainties the same way, following Irena:
perturb the measured intensity with Gaussian noise of the measured 1σ,
``I_k = I + σ · N(0, 1)``, refit, repeat, and quote the standard deviation of
each fitted quantity over the successful refits. The refit itself is
tool-specific and lives on the tool's core model
(``UnifiedFitModel.calculate_uncertainty_mc``,
``SizesDistribution.calculate_uncertainty_mc``,
``SimpleFitModel.calculate_uncertainty_mc``; Modeling and Carbon have their
own). What is shared lives here, so the GUI, batch and api/control paths can
no longer drift apart (review of 09-11-2026, issue #38):

* **No invented σ.** MC needs a measured uncertainty. When there is none the
  estimate is refused with a reason rather than run on an assumed 5 % of I,
  which would make the quoted error bar a product of that assumption. (Sizes
  with "fractional error" is the exception: there the user chose σ.)
* **Reproducible.** All noise is drawn up front from
  ``numpy.random.default_rng(seed)``; ``seed=None`` keeps passes
  nondeterministic. Never the global ``np.random`` state.
* **Honest about failures.** The result always reports ``n_ok`` of
  ``n_runs``, and fewer than :data:`MIN_SUCCESSFUL_RUNS` successes means "not
  available", never a standard deviation of one value.

numpy only — no Qt, no matplotlib.
"""

from __future__ import annotations

from typing import Callable, Optional

import numpy as np

#: Fewer successful refits than this and no uncertainty is reported.
MIN_SUCCESSFUL_RUNS = 2

#: Floor applied to perturbed intensities by tools whose models need I > 0.
POSITIVE_FLOOR = 1e-30


def check_sigma(sigma, n_points: int) -> Optional[str]:
    """Return why *sigma* cannot drive a Monte-Carlo estimate, or None if it can.

    Args:
        sigma: Per-point 1σ uncertainty, or None.
        n_points: Number of intensity points it must match.
    """
    if sigma is None:
        return ("the data have no measured uncertainties; Monte-Carlo "
                "uncertainty needs them (no σ is assumed)")
    s = np.asarray(sigma, dtype=float)
    if s.shape != (n_points,):
        return f"uncertainty array has shape {s.shape}, expected ({n_points},)"
    finite = np.isfinite(s) & (s > 0)
    if not np.any(finite):
        return "the measured uncertainties are all zero or non-finite"
    return None


def resolve_noise_sigma(measured, noise_sigma, n_points: int):
    """Pick the σ that drives the perturbations, and say where it came from.

    *measured* is the data's own uncertainty (what the fit weights by);
    *noise_sigma* is an explicit choice from the caller — in pyIrena, I × the
    Data Selector's "Generated uncertainty fraction" for data that have no
    uncertainties (the user's stated σ; issue #38). Measured σ wins when both
    are given.

    Returns ``(sigma, source, reason)``: *source* is ``"measured"`` or
    ``"generated"``; *reason* is None, or why no σ is usable.
    """
    if measured is not None and check_sigma(measured, n_points) is None:
        return np.asarray(measured, dtype=float), "measured", None
    if noise_sigma is not None:
        why = check_sigma(noise_sigma, n_points)
        if why is None:
            return np.asarray(noise_sigma, dtype=float), "generated", None
        return None, None, why
    return None, None, check_sigma(measured, n_points)


def perturbed_intensities(I, sigma, n_runs: int, seed=None,
                          floor: Optional[float] = None) -> list[np.ndarray]:
    """Draw *n_runs* noisy copies ``I + sigma · N(0, 1)`` up front.

    Drawing everything before any refit keeps results identical however the
    refits are later scheduled, and a given *seed* reproduces them exactly.
    Points whose σ is zero or non-finite are left unperturbed.

    Args:
        I: Measured intensity.
        sigma: Per-point 1σ (already validated with :func:`check_sigma`).
        n_runs: Number of copies.
        seed: Anything ``numpy.random.default_rng`` accepts; None = fresh entropy.
        floor: Optional lower clip for the perturbed intensities.
    """
    I = np.asarray(I, dtype=float)
    s = np.asarray(sigma, dtype=float)
    s = np.where(np.isfinite(s) & (s > 0), s, 0.0)
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(int(n_runs)):
        Ik = I + s * rng.standard_normal(I.shape)
        if floor is not None:
            Ik = np.maximum(Ik, floor)
        out.append(Ik)
    return out


def child_seeds(seed, n: int) -> list:
    """Independent per-refit seeds derived from *seed* (for refits that are
    themselves stochastic, e.g. the Sizes Monte-Carlo inversion).

    With ``seed=None`` every child is None, i.e. fresh entropy.
    """
    if seed is None:
        return [None] * n
    ss = seed if isinstance(seed, np.random.SeedSequence) else np.random.SeedSequence(seed)
    return list(ss.spawn(n))


def std_of(values) -> float:
    """Sample standard deviation (ddof=1) over finite values; NaN if < 2."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < MIN_SUCCESSFUL_RUNS:
        return float("nan")
    return float(np.std(v, ddof=1))


def not_available(reason: str, n_runs: int = 0, n_ok: int = 0, seed=None) -> dict:
    """The common "no uncertainty" result shape."""
    return {
        "available": False,
        "reason": reason,
        "n_runs": int(n_runs),
        "n_ok": int(n_ok),
        "seed": _seed_repr(seed),
    }


def base_result(n_runs: int, n_ok: int, seed=None) -> dict:
    """The common header of a successful result; tools add their own keys."""
    return {
        "available": True,
        "reason": None,
        "n_runs": int(n_runs),
        "n_ok": int(n_ok),
        "seed": _seed_repr(seed),
    }


def cancelled(cancel_cb: Optional[Callable[[], bool]]) -> bool:
    """True when the caller asked to stop. A failing callback never stops a run."""
    if cancel_cb is None:
        return False
    try:
        return bool(cancel_cb())
    except Exception:
        return False


def report_progress(progress_cb, n_done: int, n_runs: int) -> None:
    """Call ``progress_cb(n_done, n_runs)``; a failing callback is ignored."""
    if progress_cb is None:
        return
    try:
        progress_cb(n_done, n_runs)
    except Exception:
        pass


def _seed_repr(seed):
    """JSON-safe echo of the seed the caller passed (int or None)."""
    if seed is None:
        return None
    try:
        return int(seed)
    except (TypeError, ValueError):
        return str(seed)
