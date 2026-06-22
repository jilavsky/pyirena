"""
Uniform fit-quality display helper for all pyirena fit GUIs.

Every fitting panel (unified, sizes, simple_fits, waxs_peakfit, modeling) shows a
residuals plot and a one-line quality summary. To keep that behaviour identical
everywhere, route it through this single helper rather than re-deriving the
robust-rescale and summary text in each panel.

The residual shown is the *rescaled* residual r' = r / s, where s is the robust
(MAD-based) noise scale — it compares the scatter to the data's own noise floor
instead of the (often mis-scaled) reported sigma. See pyirena.core.fit_metrics
and docs/fit_quality_metrics.md.
"""

from __future__ import annotations

import numpy as np

from pyirena.core.fit_metrics import fit_quality_metrics, rescale_residuals


def quality_suffix(metrics: dict) -> str:
    """Build the uniform ' | σ-scale: … | max|(I−M)/I|: …%' status suffix.

    Returns an empty string if no metric is available (e.g. no sigma and no
    finite fractional misfit).
    """
    parts: list[str] = []
    s = metrics.get("robust_scale_s")
    floor = metrics.get("realistic_reduced_chi2_floor")
    if s is not None and np.isfinite(s):
        if floor is not None and np.isfinite(floor):
            parts.append(f"σ-scale: {s:.2f}× (χ²ᵣ floor: {floor:.1f})")
        else:
            parts.append(f"σ-scale: {s:.2f}×")
    mfm = metrics.get("max_abs_frac_misfit")
    if mfm is not None and np.isfinite(mfm):
        parts.append(f"max|(I−M)/I|: {mfm * 100:.1f}%")
    return (" | " + " | ".join(parts)) if parts else ""


def compute_quality_display(q, intensity, model, sigma, n_params: int = 1):
    """Compute the rescaled residual view and quality summary for a fit.

    Parameters
    ----------
    q, intensity, model, sigma : array-like
        The fit's scattering vector, measured intensity, model intensity, and
        uncertainty (sigma may be None). All on a consistent basis (same length).
    n_params : int
        Free-parameter count for reduced-χ² dof. Robust metrics do not depend on
        it; pass 1 for regularized fits where it is ill-defined.

    Returns
    -------
    (q_plot, r_prime, suffix, metrics)
        q_plot : np.ndarray   q aligned with r_prime (valid points only)
        r_prime : np.ndarray  rescaled residual r' = r/s (or rescaled fractional
                              residual when sigma is unavailable)
        suffix : str          uniform status-line suffix (see quality_suffix)
        metrics : dict        full fit_quality_metrics result
    """
    metrics = fit_quality_metrics(q, intensity, model, sigma, n_params=n_params)
    q_plot = metrics["q_valid"]

    s = metrics.get("robust_scale_s")
    if metrics["norm_residual"] is not None and s is not None and np.isfinite(s) and s > 0:
        r_prime = metrics["norm_residual"] / s
    else:
        # No usable sigma: rescale the fractional residual by its own robust scale
        # so the plot still shows scatter relative to its noise floor.
        r_prime, _ = rescale_residuals(metrics["frac_residual"])

    return q_plot, r_prime, quality_suffix(metrics), metrics
