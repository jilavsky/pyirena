"""
Robust fit-quality diagnostics for small-angle scattering fits.

This module provides a single, model-agnostic function, :func:`fit_quality_metrics`,
that operates only on numpy arrays (q, intensity, model, sigma). It therefore
serves every fit tool in pyirena (unified, sizes, simple_fits, modeling,
waxs_peakfit, ...) without any model-specific knowledge.

Motivation
----------
Fit quality is traditionally judged by reduced chi-squared plus the raw normalized
residual ``(I - M) / sigma``. That breaks down when the reported uncertainties
``sigma`` are mis-scaled (common in SAXS): a target of reduced chi-squared ~ 1 may
be physically unreachable, leading either to endless over-fitting or to genuine
misfits being dismissed as "sigma are unreliable".

The diagnostics here separate the two independent pieces of information carried by
sigma, which fail differently:

1. **Relative shape** (which points are noisier than which) -- usually trustworthy.
2. **Absolute scale** (overall noise level) -- often unreliable.

So we keep sigma for relative weighting, but let the residuals themselves recover
the true scale via a robust (MAD-based) estimate, and we add a sigma-independent
fractional-misfit backstop ``(I - M) / I`` that catches gross misfits no matter how
unreliable sigma is.

Design rule (deliberate)
------------------------
This function returns **facts only** -- no "good/bad" verdict, no thresholds, no
"stop fitting" flag. Decision logic lives in the consumer (the pyirena-ai agent
strategy, the GUI advisor, the user). Different workflows interpret the same
numbers differently.

See ``planning/ai-agent/04-fit-quality-metrics.md`` for the full rationale and
``planning/ai-agent/05-fit-quality-metrics-implementation-plan.md`` for the
implementation plan.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from pyirena.core.similarity import _longest_run

__all__ = ["fit_quality_metrics", "robust_residual_scale", "rescale_residuals"]

# Scale factor that turns the median absolute deviation into a consistent
# estimator of the Gaussian standard deviation.
_MAD_TO_SIGMA = 1.4826

# A point's fractional residual (I - M) / I is only meaningful where |I| is well
# above zero. Below this fraction of the median |I| the point is excluded from the
# fractional arrays (NaN) to avoid blow-ups after background subtraction.
_FRAC_I_FLOOR_FRACTION = 1e-3

# Minimum number of valid points required in a Q-band for its per-band statistics
# to be meaningful. n_bands is reduced adaptively so each band holds at least this.
_MIN_POINTS_PER_BAND = 5


def robust_residual_scale(norm_residual) -> float:
    """1.4826 * MAD(norm_residual) -- a Gaussian-sigma estimate robust to outliers.

    Centered on the median so a systematic offset does not inflate the scale; the
    offset itself shows up in the sign-structure diagnostics instead. Non-finite
    entries are ignored. Returns NaN if there are no finite points.

    This is the single source of truth for the robust residual scale ``s`` used
    both inside :func:`fit_quality_metrics` and by the GUI/plotting residual views.
    """
    a = np.asarray(norm_residual, dtype=float)
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return float("nan")
    center = np.median(finite)
    mad = np.median(np.abs(finite - center))
    return float(_MAD_TO_SIGMA * mad)


def rescale_residuals(norm_residual) -> tuple[np.ndarray, float]:
    """Rescale normalized residuals by the robust scale ``s``.

    Returns ``(r_prime, s)`` where ``r_prime = norm_residual / s`` re-references
    the scatter to the data's own robust noise floor instead of the (possibly
    mis-scaled) reported sigma. If ``s`` is not finite or not positive, the input
    is returned unchanged alongside that ``s`` (callers can fall back gracefully).
    """
    a = np.asarray(norm_residual, dtype=float)
    s = robust_residual_scale(a)
    if np.isfinite(s) and s > 0:
        return a / s, s
    return a, s


def _robust_scale(norm_residual: np.ndarray) -> float:
    """Internal alias kept for call sites that pass already-masked arrays."""
    return robust_residual_scale(norm_residual)


def _sign_autocorr_lag1(residual: np.ndarray) -> float:
    """Lag-1 autocorrelation of the *sign* of the residuals.

    Near +1 => long stretches of same-sign residuals (systematic misfit / wrong
    functional form). Near 0 => random scatter. Returns 0.0 when undefined.
    """
    if residual.size < 2:
        return 0.0
    signs = np.sign(residual)
    # Ignore exact zeros (carry no sign information).
    nz = signs != 0
    if np.count_nonzero(nz) < 2:
        return 0.0
    s = signs[nz].astype(float)
    s0, s1 = s[:-1], s[1:]
    denom = np.sum(s0 * s0)
    if denom == 0:
        return 0.0
    return float(np.sum(s0 * s1) / denom)


def _band_metrics(
    q: np.ndarray,
    resid_raw: np.ndarray,
    norm_residual: Optional[np.ndarray],
    frac_residual: np.ndarray,
    n_params: int,
    n_bands: int,
) -> tuple[int, list[dict]]:
    """Split the valid points into up to ``n_bands`` log-spaced Q windows.

    Returns ``(n_bands_used, bands)`` where each band is a dict with
    ``q_lo, q_hi, n, reduced_chi2, robust_scale_s, max_abs_frac_misfit``.
    The band count is reduced adaptively so each band holds at least
    ``_MIN_POINTS_PER_BAND`` points.
    """
    n = q.size
    if n == 0:
        return 0, []

    # Adaptive band count: at least _MIN_POINTS_PER_BAND points per band.
    n_bands_used = max(1, min(int(n_bands), n // _MIN_POINTS_PER_BAND))

    q_pos = q[q > 0]
    if q_pos.size == 0 or n_bands_used == 1:
        edges = np.array([np.min(q), np.max(q)], dtype=float)
        n_bands_used = 1
    else:
        edges = np.logspace(
            np.log10(np.min(q_pos)), np.log10(np.max(q_pos)), n_bands_used + 1
        )
        edges[0] = min(edges[0], np.min(q))
        edges[-1] = max(edges[-1], np.max(q))

    bands: list[dict] = []
    for b in range(n_bands_used):
        lo, hi = edges[b], edges[b + 1]
        if b == n_bands_used - 1:
            mask = (q >= lo) & (q <= hi)
        else:
            mask = (q >= lo) & (q < hi)
        nb = int(np.count_nonzero(mask))

        band: dict = {
            "q_lo": float(lo),
            "q_hi": float(hi),
            "n": nb,
            "reduced_chi2": None,
            "robust_scale_s": None,
            "max_abs_frac_misfit": None,
        }
        if nb > 0:
            if norm_residual is not None:
                nr = norm_residual[mask]
                dof_b = nb - n_params
                if dof_b > 0:
                    band["reduced_chi2"] = float(np.sum(nr ** 2) / dof_b)
                band["robust_scale_s"] = _robust_scale(nr)
            fr = frac_residual[mask]
            fr = fr[np.isfinite(fr)]
            if fr.size:
                band["max_abs_frac_misfit"] = float(np.max(np.abs(fr)))
        bands.append(band)

    return n_bands_used, bands


def fit_quality_metrics(
    q: np.ndarray,
    intensity: np.ndarray,
    model: np.ndarray,
    sigma: Optional[np.ndarray],
    n_params: int,
    n_bands: int = 4,
) -> dict:
    """Compute robust, sigma-scale-independent fit-quality diagnostics.

    Operates purely on arrays so it can be used by any fit tool. Returns *facts
    only* -- no thresholds or verdicts (those belong to the consumer).

    Parameters
    ----------
    q : np.ndarray
        Scattering vector, aligned with ``intensity``/``model``/``sigma``.
    intensity : np.ndarray
        Measured intensity ``I``.
    model : np.ndarray
        Fitted/model intensity ``M``.
    sigma : np.ndarray or None
        Reported uncertainty ``sigma`` on ``I``. May be ``None`` or contain
        non-positive / non-finite entries; such points are masked out. If no
        usable sigma remains, ``sigma_available`` is ``False`` and the
        sigma-dependent quantities are ``None`` (fractional quantities still
        computed). Sigma is never invented.
    n_params : int
        Number of free fit parameters, used only for degrees of freedom in
        reduced chi-squared. For regularized fits where this is ill-defined, pass
        the caller's best estimate; the robust metrics do not depend on it.
    n_bands : int, default 4
        Requested number of log-spaced Q bands. Reduced adaptively so each band
        keeps at least a few points; see ``n_bands_used`` in the result.

    Returns
    -------
    dict
        See module docs / the implementation plan for the full, stable schema.
        Key fields: ``sigma_available``, ``n_valid``, ``dof``, the per-point
        arrays (``q_valid``, ``norm_residual``, ``frac_residual``,
        ``frac_uncertainty``), the global scalars (``reduced_chi2``,
        ``robust_scale_s``, ``sigma_misscale_factor``,
        ``realistic_reduced_chi2_floor``, ``median_frac_uncertainty``,
        ``max_abs_frac_misfit``, ``q_at_max_frac_misfit``, ``n_outliers_3s``,
        ``frac_outliers_3s``), the structure scalars (``longest_same_sign_run``,
        ``sign_autocorr_lag1``), and the per-band list (``n_bands_used``,
        ``bands``).
    """
    q = np.asarray(q, dtype=float).ravel()
    intensity = np.asarray(intensity, dtype=float).ravel()
    model = np.asarray(model, dtype=float).ravel()
    n_params = int(n_params)

    # --- Determine sigma availability and align everything to a common length ---
    sigma_arr: Optional[np.ndarray]
    if sigma is None:
        sigma_arr = None
    else:
        sigma_arr = np.asarray(sigma, dtype=float).ravel()

    # Base finiteness mask on q, I, M.
    finite = np.isfinite(q) & np.isfinite(intensity) & np.isfinite(model)

    if sigma_arr is not None and sigma_arr.shape == q.shape:
        sigma_ok = np.isfinite(sigma_arr) & (sigma_arr > 0)
        sigma_available = bool(np.any(finite & sigma_ok))
    else:
        sigma_arr = None
        sigma_ok = np.zeros_like(finite)
        sigma_available = False

    if sigma_available:
        valid = finite & sigma_ok
    else:
        valid = finite

    qv = q[valid]
    Iv = intensity[valid]
    Mv = model[valid]
    sv = sigma_arr[valid] if sigma_available else None
    n_valid = int(qv.size)

    resid_raw = Iv - Mv  # (I - M), raw

    # --- Fractional arrays (sigma-independent backstop) ---
    # Guard against |I| ~ 0 (e.g. after background subtraction) where (I-M)/I and
    # sigma/I blow up: such points become NaN and are excluded from extrema.
    frac_residual = np.full(n_valid, np.nan)
    frac_uncertainty: Optional[np.ndarray] = None
    if n_valid:
        i_floor = _FRAC_I_FLOOR_FRACTION * np.median(np.abs(Iv)) if np.any(Iv) else 0.0
        good_I = np.isfinite(Iv) & (np.abs(Iv) > i_floor)
        frac_residual[good_I] = resid_raw[good_I] / Iv[good_I]
        if sigma_available:
            frac_uncertainty = np.full(n_valid, np.nan)
            frac_uncertainty[good_I] = sv[good_I] / Iv[good_I]

    # --- Normalized residual and sigma-dependent scalars ---
    norm_residual: Optional[np.ndarray] = None
    reduced_chi2: Optional[float] = None
    robust_scale_s: Optional[float] = None
    realistic_floor: Optional[float] = None
    n_outliers_3s: Optional[int] = None
    frac_outliers_3s: Optional[float] = None
    dof: Optional[int] = None

    if sigma_available and n_valid:
        norm_residual = resid_raw / sv
        dof = n_valid - n_params
        if dof > 0:
            reduced_chi2 = float(np.sum(norm_residual ** 2) / dof)
        robust_scale_s = _robust_scale(norm_residual)
        if robust_scale_s is not None and np.isfinite(robust_scale_s):
            realistic_floor = float(robust_scale_s ** 2)
            if robust_scale_s > 0:
                outlier_mask = np.abs(norm_residual) > 3.0 * robust_scale_s
                n_outliers_3s = int(np.count_nonzero(outlier_mask))
                frac_outliers_3s = float(n_outliers_3s / n_valid)

    # --- Fractional-misfit extrema (sigma-independent) ---
    median_frac_uncertainty = float("nan")
    if frac_uncertainty is not None:
        fu = frac_uncertainty[np.isfinite(frac_uncertainty)]
        if fu.size:
            median_frac_uncertainty = float(np.median(np.abs(fu)))

    max_abs_frac_misfit = float("nan")
    q_at_max_frac_misfit = float("nan")
    fr_finite_mask = np.isfinite(frac_residual)
    if np.any(fr_finite_mask):
        abs_fr = np.abs(frac_residual)
        abs_fr[~fr_finite_mask] = -np.inf
        idx = int(np.argmax(abs_fr))
        max_abs_frac_misfit = float(np.abs(frac_residual[idx]))
        q_at_max_frac_misfit = float(qv[idx])

    # --- Structure scalars ---
    longest_same_sign_run = int(_longest_run(resid_raw)) if n_valid else 0
    sign_autocorr_lag1 = _sign_autocorr_lag1(resid_raw) if n_valid else 0.0

    # --- Per-band ---
    n_bands_used, bands = _band_metrics(
        qv, resid_raw, norm_residual, frac_residual, n_params, n_bands
    )

    return {
        "sigma_available": sigma_available,
        "n_valid": n_valid,
        "n_params": n_params,
        "dof": dof,
        # per-point arrays (aligned to q_valid)
        "q_valid": qv,
        "norm_residual": norm_residual,
        "frac_residual": frac_residual,
        "frac_uncertainty": frac_uncertainty,
        # global scalars
        "reduced_chi2": reduced_chi2,
        "robust_scale_s": robust_scale_s,
        "sigma_misscale_factor": robust_scale_s,  # alias
        "realistic_reduced_chi2_floor": realistic_floor,
        "median_frac_uncertainty": median_frac_uncertainty,
        "max_abs_frac_misfit": max_abs_frac_misfit,
        "q_at_max_frac_misfit": q_at_max_frac_misfit,
        "n_outliers_3s": n_outliers_3s,
        "frac_outliers_3s": frac_outliers_3s,
        # structure scalars
        "longest_same_sign_run": longest_same_sign_run,
        "sign_autocorr_lag1": sign_autocorr_lag1,
        # per-band
        "n_bands_used": n_bands_used,
        "bands": bands,
    }
