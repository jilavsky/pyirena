"""pyirena.api.control.sizes — agent-drivable Size Distribution fitting.

This module is the Sizes counterpart of
:mod:`pyirena.api.control.unified_fit`.  It lets an LLM agent (or any
automation script) configure a particle size-distribution inversion, fit the
complex background, run the inversion, inspect the result, and save it to
NXcanSAS — all through plain-dict tool calls.

The functions here operate on the *same* :class:`~pyirena.api.control.session.Session`
objects used by the Unified Fit control surface, so the shared session-lifecycle
and Q-range tools (``open_dataset``, ``list_open_sessions``, ``close_session``,
``get_data_q_range`` / ``get_fit_q_range`` / ``set_fit_q_range`` /
``reset_fit_q_range``) are reused as-is.  ``set_fit_q_range`` defines the
**inversion Q-range** (the batch's ``cursor_q_min`` / ``cursor_q_max``).

The Sizes model state lives in ``session.model`` (a
:class:`~pyirena.core.sizes.SizesDistribution`) with ``session.model_name ==
"sizes"``.

Typical workflow
----------------
>>> import pyirena.api.control as ctrl
>>> sid = ctrl.open_dataset("/data/sample.h5")["session_id"]
>>> ctrl.suggest_sizes_setup(sid)                 # data-driven recommendations
>>> ctrl.select_sizes_model(sid, method="maxent")
>>> ctrl.set_shape(sid, "sphere", contrast=...)
>>> ctrl.set_size_grid(sid, r_min=10, r_max=500, n_bins=200)
>>> ctrl.set_error_handling(sid, error_scale=1.0)
>>> ctrl.fit_power_law_background(sid, q_lo, q_hi, fit_B=True, fit_P=True)
>>> ctrl.fit_flat_background(sid, q_lo2, q_hi2)
>>> ctrl.set_fit_q_range(sid, q_min, q_max)       # inversion window (shared tool)
>>> ctrl.run_sizes_fit(sid)
>>> ctrl.get_sizes_fit_image(sid)
>>> ctrl.save_sizes_fit(sid)

All functions return plain dicts.  Errors are
``{"error": ..., "code": ..., "suggestion": ...}`` dicts rather than exceptions.
"""
from __future__ import annotations

import base64
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np

from pyirena.api.control.errors import make_error, no_session, no_model, no_fit
from pyirena.api.control.session import get_session, fit_mask, Session

# Inversion methods exposed to the agent.  MaxEnt is the recommended default.
SIZES_METHODS = ["maxent", "regularization", "tnnls", "montecarlo"]
SIZES_SHAPES = ["sphere", "spheroid"]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _require_sizes(session_id: str):
    """Return (session, None) when a Sizes model is ready, else (None, error)."""
    s = get_session(session_id)
    if s is None:
        return None, no_session(session_id)
    if s.model is None or s.model_name != "sizes":
        return None, make_error(
            f"Session '{session_id}' has no Sizes model selected.",
            suggestion="Call select_sizes_model() first.",
            code="NO_SIZES_MODEL",
        )
    return s, None


def _decimated(arr, max_points: int = 500) -> list:
    """Thin an array to at most max_points for JSON transport (NaN/inf -> None)."""
    a = np.asarray(arr, dtype=float)
    n = len(a)
    if n > max_points:
        step = max(1, n // max_points)
        a = a[::step]
    return [v if np.isfinite(v) else None for v in a.tolist()]


def _config_dict(model) -> dict:
    """Full current configuration of a SizesDistribution as a flat dict."""
    ar = model.shape_params.get("aspect_ratio") if model.shape == "spheroid" else None
    return {
        # grid
        "r_min": float(model.r_min),
        "r_max": float(model.r_max),
        "n_bins": int(model.n_bins),
        "log_spacing": bool(model.log_spacing),
        # shape
        "shape": str(model.shape),
        "contrast": float(model.contrast),
        "aspect_ratio": (float(ar) if ar is not None else None),
        # method
        "method": str(model.method),
        "maxent_sky_background": float(model.maxent_sky_background),
        "maxent_stability": float(model.maxent_stability),
        "maxent_max_iter": int(model.maxent_max_iter),
        "regularization_evalue": float(model.regularization_evalue),
        "regularization_min_ratio": float(model.regularization_min_ratio),
        "tnnls_approach_param": float(model.tnnls_approach_param),
        "tnnls_max_iter": int(model.tnnls_max_iter),
        "montecarlo_n_repetitions": int(model.montecarlo_n_repetitions),
        "montecarlo_convergence": float(model.montecarlo_convergence),
        "montecarlo_max_iter": int(model.montecarlo_max_iter),
        # error handling
        "error_scale": float(model.error_scale),
        "fractional_error": bool(model.fractional_error),
        "fractional_error_value": float(model.fractional_error_value),
        # complex background
        "power_law_B": float(model.power_law_B),
        "power_law_P": float(model.power_law_P),
        "background": float(model.background),
        # background-fit windows (tracked on the model for save/setup_state)
        "power_law_q_min": getattr(model, "_pl_q_min", None),
        "power_law_q_max": getattr(model, "_pl_q_max", None),
        "fit_power_law_B": bool(getattr(model, "_fit_pl_B", False)),
        "fit_power_law_P": bool(getattr(model, "_fit_pl_P", False)),
        "background_q_min": getattr(model, "_bg_q_min", None),
        "background_q_max": getattr(model, "_bg_q_max", None),
    }


def _render_image(fig, session_id: str, tag: str, dpi: int) -> str:
    import matplotlib.pyplot as plt  # noqa: PLC0415

    tmp = Path(tempfile.gettempdir()) / "pyirena-ctrl"
    tmp.mkdir(parents=True, exist_ok=True)
    out = tmp / f"sizes_{tag}_{session_id}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(out.read_bytes()).decode("ascii")


# ---------------------------------------------------------------------------
# Model lifecycle
# ---------------------------------------------------------------------------

def select_sizes_model(session_id: str, method: str = "maxent") -> dict:
    """Create a Size Distribution model for the session.

    Parameters
    ----------
    method : str
        Inversion method.  One of ``"maxent"`` (recommended default),
        ``"regularization"``, ``"tnnls"``, ``"montecarlo"``.

    Notes
    -----
    MaxEnt produces the smoothest distribution consistent with the data and is
    the safest default.  Regularization is a good alternative.  TNNLS enforces
    strict non-negativity without smoothing.  Monte Carlo (McSAS) is slower and
    stochastic.  Calling this again replaces the model and clears any prior fit.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    method = (method or "maxent").lower()
    if method not in SIZES_METHODS:
        return make_error(
            f"Unknown method '{method}'.",
            suggestion=f"Choose one of {SIZES_METHODS}.",
            code="BAD_METHOD",
        )

    from pyirena.core.sizes import SizesDistribution  # noqa: PLC0415

    model = SizesDistribution()
    model.method = method
    # Background-fit window bookkeeping (consumed by save / get_sizes_config)
    model._pl_q_min = None
    model._pl_q_max = None
    model._fit_pl_B = False
    model._fit_pl_P = False
    model._bg_q_min = None
    model._bg_q_max = None

    s.model = model
    s.model_name = "sizes"
    s.last_fit_result = None

    return {"ok": True, "model": "sizes", "config": _config_dict(model)}


def get_sizes_config(session_id: str) -> dict:
    """Return the current Sizes model configuration (grid, shape, method, errors,
    complex background)."""
    s, err = _require_sizes(session_id)
    if err:
        return err
    return {"ok": True, "config": _config_dict(s.model)}


# ---------------------------------------------------------------------------
# Size grid & shape
# ---------------------------------------------------------------------------

def set_size_grid(
    session_id: str,
    r_min: Optional[float] = None,
    r_max: Optional[float] = None,
    n_bins: Optional[int] = None,
    log_spacing: Optional[bool] = None,
) -> dict:
    """Set the radius grid for the inversion.

    Parameters
    ----------
    r_min, r_max : float
        Radius range [Å].  Should bracket the real particle sizes; a useful
        starting heuristic is r ≈ π/Q over the inversion Q-range.
    n_bins : int
        Number of radius bins (typical 100–300).
    log_spacing : bool
        Log-spaced bins (recommended for wide size ranges).
    """
    s, err = _require_sizes(session_id)
    if err:
        return err
    m = s.model
    if r_min is not None:
        m.r_min = float(r_min)
    if r_max is not None:
        m.r_max = float(r_max)
    if n_bins is not None:
        if int(n_bins) < 5:
            return make_error("n_bins must be at least 5.", code="BAD_NBINS")
        m.n_bins = int(n_bins)
    if log_spacing is not None:
        m.log_spacing = bool(log_spacing)
    if m.r_min >= m.r_max:
        return make_error(
            f"r_min ({m.r_min}) must be < r_max ({m.r_max}).",
            code="BAD_RANGE",
        )
    return {"ok": True, "config": _config_dict(m)}


def set_shape(
    session_id: str,
    shape: Optional[str] = None,
    contrast: Optional[float] = None,
    aspect_ratio: Optional[float] = None,
) -> dict:
    """Set the particle form factor.

    Parameters
    ----------
    shape : str
        ``"sphere"`` or ``"spheroid"``.
    contrast : float
        Scattering contrast (Δρ)² in units of 10²⁰ cm⁻⁴.  Sets the absolute
        scale of the recovered volume fraction; use 1.0 if unknown.
    aspect_ratio : float
        Spheroid axis ratio (only used when shape == "spheroid").
    """
    s, err = _require_sizes(session_id)
    if err:
        return err
    m = s.model
    if shape is not None:
        sh = str(shape).lower()
        if sh not in SIZES_SHAPES:
            return make_error(
                f"Unknown shape '{shape}'.",
                suggestion=f"Choose one of {SIZES_SHAPES}.",
                code="BAD_SHAPE",
            )
        m.shape = sh
    if contrast is not None:
        m.contrast = float(contrast)
    if aspect_ratio is not None:
        m.shape_params["aspect_ratio"] = float(aspect_ratio)
    # Ensure spheroid always carries an aspect ratio
    if m.shape == "spheroid" and "aspect_ratio" not in m.shape_params:
        m.shape_params["aspect_ratio"] = 1.5
    return {"ok": True, "config": _config_dict(m)}


# ---------------------------------------------------------------------------
# Inversion method & parameters
# ---------------------------------------------------------------------------

def set_method(
    session_id: str,
    method: str,
    maxent_sky_background: Optional[float] = None,
    maxent_max_iter: Optional[int] = None,
    regularization_evalue: Optional[float] = None,
    regularization_min_ratio: Optional[float] = None,
    tnnls_approach_param: Optional[float] = None,
    tnnls_max_iter: Optional[int] = None,
    montecarlo_n_repetitions: Optional[int] = None,
    montecarlo_convergence: Optional[float] = None,
    montecarlo_max_iter: Optional[int] = None,
) -> dict:
    """Choose the inversion method and (optionally) its tuning parameters.

    Only the parameters relevant to the chosen *method* are applied; the rest
    are ignored.  MaxEnt is the recommended default for most data.

    Method parameters
    ------------------
    maxent  : maxent_sky_background, maxent_max_iter
    regularization : regularization_evalue, regularization_min_ratio
    tnnls   : tnnls_approach_param, tnnls_max_iter
    montecarlo : montecarlo_n_repetitions, montecarlo_convergence, montecarlo_max_iter
    """
    s, err = _require_sizes(session_id)
    if err:
        return err
    m = s.model
    method = (method or "").lower()
    if method not in SIZES_METHODS:
        return make_error(
            f"Unknown method '{method}'.",
            suggestion=f"Choose one of {SIZES_METHODS}.",
            code="BAD_METHOD",
        )
    m.method = method

    if method == "maxent":
        if maxent_sky_background is not None:
            m.maxent_sky_background = float(maxent_sky_background)
        if maxent_max_iter is not None:
            m.maxent_max_iter = int(maxent_max_iter)
    elif method == "regularization":
        if regularization_evalue is not None:
            m.regularization_evalue = float(regularization_evalue)
        if regularization_min_ratio is not None:
            m.regularization_min_ratio = float(regularization_min_ratio)
    elif method == "tnnls":
        if tnnls_approach_param is not None:
            m.tnnls_approach_param = float(tnnls_approach_param)
        if tnnls_max_iter is not None:
            m.tnnls_max_iter = int(tnnls_max_iter)
    elif method == "montecarlo":
        if montecarlo_n_repetitions is not None:
            m.montecarlo_n_repetitions = int(montecarlo_n_repetitions)
        if montecarlo_convergence is not None:
            m.montecarlo_convergence = float(montecarlo_convergence)
        if montecarlo_max_iter is not None:
            m.montecarlo_max_iter = int(montecarlo_max_iter)

    return {"ok": True, "config": _config_dict(m)}


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

def set_error_handling(
    session_id: str,
    error_scale: Optional[float] = None,
    fractional_error: Optional[bool] = None,
    fractional_error_value: Optional[float] = None,
) -> dict:
    """Configure how measurement uncertainties are used during the inversion.

    Two mutually exclusive modes:

    - **Error scaling** (default): use the file/measured σ multiplied by
      ``error_scale`` (1.0 = unchanged).  Increase it when the reported error
      bars are too small and the fit chases noise.
    - **Fractional error**: set ``fractional_error=True`` to ignore the file σ
      entirely and use σ = |I| × ``fractional_error_value`` (e.g. 0.03 = 3%).
      Useful when the file uncertainties are unreliable.
    """
    s, err = _require_sizes(session_id)
    if err:
        return err
    m = s.model
    if error_scale is not None:
        m.error_scale = float(error_scale)
    if fractional_error is not None:
        m.fractional_error = bool(fractional_error)
    if fractional_error_value is not None:
        m.fractional_error_value = float(fractional_error_value)
    return {"ok": True, "config": _config_dict(m)}


# ---------------------------------------------------------------------------
# Complex background (power-law B·q^-P + flat term)
# ---------------------------------------------------------------------------

def set_background(
    session_id: str,
    power_law_B: Optional[float] = None,
    power_law_P: Optional[float] = None,
    background: Optional[float] = None,
) -> dict:
    """Set the complex-background terms directly (without fitting).

    The background subtracted before inversion is ``power_law_B·q^(-power_law_P)
    + background``.  Set ``power_law_B = 0`` to use a flat background only.
    """
    s, err = _require_sizes(session_id)
    if err:
        return err
    m = s.model
    if power_law_B is not None:
        m.power_law_B = float(power_law_B)
    if power_law_P is not None:
        m.power_law_P = float(power_law_P)
    if background is not None:
        m.background = float(background)
    return {"ok": True, "config": _config_dict(m)}


def fit_power_law_background(
    session_id: str,
    q_min: float,
    q_max: float,
    fit_B: bool = True,
    fit_P: bool = True,
) -> dict:
    """Fit the power-law background term B·q^(-P) over a Q-window.

    Typically the low-Q region where a steep slope dominates.  Updates
    ``power_law_B`` and/or ``power_law_P`` on the model.  Use the full data Q
    (not the inversion Q-range) — the window here is independent.

    Parameters
    ----------
    q_min, q_max : float
        Q-window [Å⁻¹] over which to fit the power law.
    fit_B, fit_P : bool
        Which parameters to vary (at least one must be True).
    """
    s, err = _require_sizes(session_id)
    if err:
        return err
    m = s.model
    res = m.fit_power_law(s.q, s.intensity, float(q_min), float(q_max),
                          fit_B=bool(fit_B), fit_P=bool(fit_P))
    if not res.get("success"):
        return make_error(res.get("message", "Power-law fit failed."),
                          code="POWERLAW_FIT_FAILED")
    # Remember the window/flags so save_sizes_fit() can reproduce the pre-fit.
    m._pl_q_min = float(q_min)
    m._pl_q_max = float(q_max)
    m._fit_pl_B = bool(fit_B)
    m._fit_pl_P = bool(fit_P)
    return {
        "ok": True,
        "power_law_B": res["B"],
        "power_law_P": res["P"],
        "message": res.get("message", ""),
    }


def fit_flat_background(session_id: str, q_min: float, q_max: float) -> dict:
    """Fit the flat background by averaging I − B·q^(-P) over a Q-window.

    Typically the high-Q region where the curve flattens to a constant.
    Updates ``background`` on the model.  Run this *after*
    ``fit_power_law_background`` if both terms are present.

    Parameters
    ----------
    q_min, q_max : float
        Q-window [Å⁻¹] over which to average the flat background.
    """
    s, err = _require_sizes(session_id)
    if err:
        return err
    m = s.model
    res = m.fit_background_term(s.q, s.intensity, float(q_min), float(q_max))
    if not res.get("success"):
        return make_error(res.get("message", "Background fit failed."),
                          code="BACKGROUND_FIT_FAILED")
    m._bg_q_min = float(q_min)
    m._bg_q_max = float(q_max)
    return {"ok": True, "background": res["background"], "message": res.get("message", "")}


def get_background_preview_image(
    session_id: str, width: int = 1024, height: int = 768, dpi: int = 120,
) -> dict:
    """Render the data with the current complex background overlaid (log-log PNG).

    Use this to visually confirm the background before running the inversion.
    Returns ``{"image_base64": ...}``.
    """
    s, err = _require_sizes(session_id)
    if err:
        return err
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt  # noqa: PLC0415

    m = s.model
    q, I = s.q, s.intensity
    valid = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
    qv, Iv = q[valid], I[valid]
    bg = m.compute_complex_background(qv)

    fig, ax = plt.subplots(1, 1, figsize=(width / dpi, height / dpi))
    if s.error is not None:
        ax.errorbar(qv, Iv, yerr=s.error[valid], fmt="o", markersize=3, capsize=2,
                    alpha=0.6, color="steelblue", label="Data", linewidth=0.5)
    else:
        ax.plot(qv, Iv, "o", markersize=3, alpha=0.6, color="steelblue", label="Data")
    ax.plot(qv, bg, "-", linewidth=2, color="darkorange", label="Complex background")
    ax.plot(qv, np.maximum(Iv - bg, 1e-300), "--", linewidth=1.2, color="green",
            alpha=0.7, label="Data − background")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Q  (Å⁻¹)", fontsize=10)
    ax.set_ylabel("I  (cm⁻¹)", fontsize=10)
    ax.set_title(
        f"Background preview  (B={m.power_law_B:.3g}, P={m.power_law_P:.3g}, "
        f"flat={m.background:.3g})",
        fontsize=11, fontweight="bold",
    )
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    return {"ok": True, "image_base64": _render_image(fig, session_id, "bg", dpi)}


# ---------------------------------------------------------------------------
# Fit execution & results
# ---------------------------------------------------------------------------

def run_sizes_fit(session_id: str, random_seed: Optional[int] = None) -> dict:
    """Run the size-distribution inversion on the current model and data.

    The inversion Q-range (set_fit_q_range) is applied first.  The complex
    background is subtracted internally by the model before inversion.

    Returns
    -------
    dict with keys: success, chi_squared, volume_fraction, rg, peak_r,
    n_iterations, n_data, message.
    """
    s, err = _require_sizes(session_id)
    if err:
        return err

    mask = fit_mask(s)
    if not np.any(mask):
        return make_error(
            "No data points in the current fit Q range.",
            suggestion="Call reset_fit_q_range() or widen set_fit_q_range().",
            code="EMPTY_RANGE",
        )

    q_fit = s.q[mask]
    I_fit = s.intensity[mask]
    err_fit = s.error[mask] if s.error is not None else None

    if random_seed is not None:
        np.random.seed(int(random_seed))

    try:
        result = s.model.fit(q_fit, I_fit, err_fit)
    except Exception as exc:  # pragma: no cover - defensive
        return make_error(
            f"Fit failed with exception: {exc}",
            suggestion="Check the size grid, contrast, and background settings.",
            code="FIT_EXCEPTION",
        )

    if not result.get("success", False):
        return make_error(
            result.get("message", "Inversion did not converge."),
            code="FIT_FAILED",
        )

    result["random_seed"] = random_seed
    s.last_fit_result = result

    dist = result.get("distribution")
    r_grid = result.get("r_grid")
    peak_r = None
    if dist is not None and r_grid is not None and len(dist) > 0:
        peak_r = float(r_grid[int(np.argmax(dist))])

    return {
        "success": True,
        "chi_squared": float(result.get("chi_squared", float("nan"))),
        "volume_fraction": float(result.get("volume_fraction", float("nan"))),
        "rg": float(result.get("rg", float("nan"))),
        "peak_r": peak_r,
        "n_iterations": int(result.get("n_iterations", 0)),
        "n_data": int(result.get("n_data", len(q_fit))),
        "random_seed": random_seed,
        "message": str(result.get("message", "")),
    }


def get_sizes_distribution(session_id: str, max_points: int = 500) -> dict:
    """Return the fitted size distribution arrays (decimated for transport).

    Returns ``r_grid`` [Å] and ``distribution`` P(r) [volume fraction / Å],
    plus ``distribution_std`` when available.
    """
    s, err = _require_sizes(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)
    res = s.last_fit_result
    out = {
        "ok": True,
        "r_grid": _decimated(res["r_grid"], max_points),
        "distribution": _decimated(res["distribution"], max_points),
    }
    std = res.get("distribution_std")
    out["distribution_std"] = _decimated(std, max_points) if std is not None else None
    return out


def get_sizes_results(session_id: str) -> dict:
    """Return the full scalar result + configuration for the last fit."""
    s, err = _require_sizes(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)
    res = s.last_fit_result
    dist = res.get("distribution")
    r_grid = res.get("r_grid")
    peak_r = None
    if dist is not None and r_grid is not None and len(dist) > 0:
        peak_r = float(r_grid[int(np.argmax(dist))])

    params = _config_dict(s.model)
    params.update({
        "chi_squared": float(res.get("chi_squared", float("nan"))),
        "volume_fraction": float(res.get("volume_fraction", float("nan"))),
        "rg": float(res.get("rg", float("nan"))),
        "peak_r": peak_r,
        "n_iterations": int(res.get("n_iterations", 0)),
        "n_data": int(res.get("n_data", 0)),
    })
    return {"ok": True, "parameters": params}


def get_sizes_fit_image(
    session_id: str, width: int = 1024, height: int = 900, dpi: int = 120,
) -> dict:
    """Render the fit as a two-panel PNG: (top) log-log data + model (+ background),
    (bottom) the size distribution P(r) vs r.  Returns ``{"image_base64": ...}``."""
    s, err = _require_sizes(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt  # noqa: PLC0415

    m = s.model
    res = s.last_fit_result
    q = np.asarray(res["q"], dtype=float)
    I_obs = np.asarray(res["I_data"], dtype=float)
    model_I = np.asarray(res["model_intensity"], dtype=float)
    bg = m.compute_complex_background(q)
    model_total = model_I + bg
    r_grid = np.asarray(res["r_grid"], dtype=float)
    dist = np.asarray(res["distribution"], dtype=float)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(width / dpi, height / dpi),
        gridspec_kw={"height_ratios": [3, 2]},
    )

    # Top: I(Q)
    err_arr = res.get("err")
    if err_arr is not None:
        ax1.errorbar(q, I_obs, yerr=np.asarray(err_arr, dtype=float), fmt="o",
                     markersize=3, capsize=2, alpha=0.6, color="steelblue",
                     label="Data", linewidth=0.5)
    else:
        ax1.plot(q, I_obs, "o", markersize=3, alpha=0.6, color="steelblue", label="Data")
    ax1.plot(q, model_total, "-", linewidth=2, color="red", label="Model + bg")
    if np.any(bg != 0.0):
        ax1.plot(q, bg, "--", linewidth=1.2, color="darkorange", alpha=0.7,
                 label="Background")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Q  (Å⁻¹)", fontsize=10)
    ax1.set_ylabel("I  (cm⁻¹)", fontsize=10)
    chi2 = res.get("chi_squared")
    n_data = res.get("n_data")
    title = f"Size Distribution ({m.method})"
    if chi2 is not None:
        title += f"   χ² = {chi2:.3g}"
        if n_data:
            title += f"  (n={n_data})"
    ax1.set_title(title, fontsize=11, fontweight="bold")
    ax1.legend(fontsize=8, loc="best")
    ax1.grid(True, which="both", alpha=0.25)

    # Bottom: P(r)
    ax2.plot(r_grid, dist, "-", linewidth=1.8, color="purple")
    std = res.get("distribution_std")
    if std is not None:
        std = np.asarray(std, dtype=float)
        ax2.fill_between(r_grid, np.maximum(dist - std, 0.0), dist + std,
                         color="purple", alpha=0.2)
    if m.log_spacing:
        ax2.set_xscale("log")
    ax2.set_xlabel("Radius r  (Å)", fontsize=10)
    ax2.set_ylabel("P(r)  (vol-frac / Å)", fontsize=10)
    ax2.grid(True, which="both", alpha=0.25)

    fig.tight_layout()
    return {"ok": True, "image_base64": _render_image(fig, session_id, "fit", dpi)}


# ---------------------------------------------------------------------------
# Guidance: suitability + auto-setup
# ---------------------------------------------------------------------------

def suggest_sizes_setup(session_id: str) -> dict:
    """Inspect the loaded data and recommend a Sizes setup, with suitability warnings.

    Returns
    -------
    dict with:
      - ``suitable``     : bool — whether the data looks like a viable single
                           size-distribution candidate.
      - ``recommended``  : {r_min, r_max, inversion_q_min, inversion_q_max,
                            power_law_q_min, power_law_q_max,
                            background_q_min, background_q_max}
      - ``warnings``     : list of human-readable caveats.
      - ``features``     : raw feature-detection summary (segments, knees).

    This is advisory only — it does not modify the model.  Apply the values via
    set_size_grid() / set_fit_q_range() / fit_power_law_background() /
    fit_flat_background() as appropriate.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    from pyirena.core.sizes import recommend_sizes_setup  # noqa: PLC0415

    rec = recommend_sizes_setup(s.q, s.intensity, sigma=s.error)
    return {"ok": True, **rec}


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def _model_to_gui_state(s: Session) -> dict:
    """Build a Sizes GUI-state dict (matches SizesFitPanel._get_current_state)
    so the saved file can be restored via 'Load Setup from File'."""
    cfg = _config_dict(s.model)
    cfg["aspect_ratio"] = cfg.get("aspect_ratio") or 1.5
    cfg["cursor_q_min"] = s.fit_q_min
    cfg["cursor_q_max"] = s.fit_q_max
    return cfg


def save_sizes_fit(session_id: str, output_path: Optional[str] = None) -> dict:
    """Save the fitted size distribution to a NXcanSAS HDF5 file.

    Parameters
    ----------
    output_path : str, optional
        Where to save.  Defaults to the original file (in-place).  Saving to a
        different path preserves the original.
    """
    s, err = _require_sizes(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    from pyirena.io.nxcansas_sizes import save_sizes_results  # noqa: PLC0415

    m = s.model
    res = s.last_fit_result
    target = Path(output_path) if output_path else Path(s.file_path)

    q_fit = np.asarray(res["q"], dtype=float)
    I_data = np.asarray(res["I_data"], dtype=float)
    model_I = np.asarray(res["model_intensity"], dtype=float)
    residuals = np.asarray(res["residuals"], dtype=float)
    r_grid = np.asarray(res["r_grid"], dtype=float)
    dist = np.asarray(res["distribution"], dtype=float)
    err_arr = res.get("err")
    err_arr = np.asarray(err_arr, dtype=float) if err_arr is not None else None

    peak_r = float(r_grid[int(np.argmax(dist))]) if len(dist) else float("nan")
    ar = m.shape_params.get("aspect_ratio") if m.shape == "spheroid" else None

    params = {
        "chi_squared": float(res.get("chi_squared", float("nan"))),
        "volume_fraction": float(res.get("volume_fraction", float("nan"))),
        "rg": float(res.get("rg", float("nan"))),
        "peak_r": peak_r,
        "n_iterations": int(res.get("n_iterations", 0)),
        "method": m.method,
        "shape": m.shape,
        "contrast": m.contrast,
        "aspect_ratio": ar,
        "r_min": m.r_min,
        "r_max": m.r_max,
        "n_bins": m.n_bins,
        "log_spacing": m.log_spacing,
        "background": m.background,
        "error_scale": m.error_scale,
        "fractional_error": m.fractional_error,
        "fractional_error_value": m.fractional_error_value,
        "power_law_B": m.power_law_B,
        "power_law_P": m.power_law_P,
        "maxent_sky_background": m.maxent_sky_background,
        "maxent_stability": m.maxent_stability,
        "maxent_max_iter": m.maxent_max_iter,
        "regularization_evalue": m.regularization_evalue,
        "regularization_min_ratio": m.regularization_min_ratio,
        "tnnls_approach_param": m.tnnls_approach_param,
        "tnnls_max_iter": m.tnnls_max_iter,
        "montecarlo_n_repetitions": m.montecarlo_n_repetitions,
        "montecarlo_convergence": m.montecarlo_convergence,
        "montecarlo_max_iter": m.montecarlo_max_iter,
        "cursor_q_min": s.fit_q_min,
        "cursor_q_max": s.fit_q_max,
        "power_law_q_min": getattr(m, "_pl_q_min", None),
        "power_law_q_max": getattr(m, "_pl_q_max", None),
        "background_q_min": getattr(m, "_bg_q_min", None),
        "background_q_max": getattr(m, "_bg_q_max", None),
    }

    # Robust fit-quality metrics on the raw-basis triple (observed I, model+bg, err)
    fq_metrics = None
    try:
        from pyirena.core.fit_metrics import fit_quality_metrics  # noqa: PLC0415
        bg = m.compute_complex_background(q_fit)
        fq_metrics = fit_quality_metrics(q_fit, I_data, model_I + bg, err_arr, n_params=1)
    except Exception:
        fq_metrics = None

    try:
        save_sizes_results(
            filepath=target,
            q=q_fit,
            intensity_data=I_data,
            intensity_model=model_I,
            residuals=residuals,
            r_grid=r_grid,
            distribution=dist,
            params=params,
            intensity_error=err_arr,
            distribution_std=res.get("distribution_std"),
            fit_quality=fq_metrics,
            setup_state=_model_to_gui_state(s),
        )
    except Exception as exc:
        return make_error(
            f"Could not save Sizes fit to '{target}': {exc}",
            suggestion="Check write permissions and that the file is a valid NXcanSAS HDF5.",
            code="SAVE_ERROR",
        )

    return {"ok": True, "file_path": str(target)}
