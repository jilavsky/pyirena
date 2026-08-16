"""pyirena.api.control.simple_fits — agent-drivable Simple Fits.

The Simple Fits counterpart of :mod:`pyirena.api.control.unified_fit` and
:mod:`pyirena.api.control.sizes`.  It lets an LLM agent (or any automation
script) pick one of the analytical models, set and constrain its parameters,
fit a Q sub-range, inspect the result, and save it to NXcanSAS — all through
plain-dict tool calls.

These functions operate on the *same*
:class:`~pyirena.api.control.session.Session` objects as the other control
surfaces, so the shared session-lifecycle and Q-range tools (``open_dataset``,
``list_open_sessions``, ``close_session``, ``get_data_q_range`` /
``get_fit_q_range`` / ``set_fit_q_range`` / ``reset_fit_q_range``) are reused
as-is; ``set_fit_q_range`` defines the fitted Q window.

The model state lives in ``session.model`` (a
:class:`~pyirena.core.simple_fits.SimpleFitModel`) with
``session.model_name == "simple_fits"``.

Which parameters are fitted is tracked per session rather than on the model:
:class:`~pyirena.core.simple_fits.SimpleFitModel` takes the held-fixed set as a
``fixed_params`` argument to :meth:`~pyirena.core.simple_fits.SimpleFitModel.fit`
(the GUI holds it in its "Fit?" checkboxes), so this module keeps the agent's
choices in ``session._simple_fixed``.

Typical workflow
----------------
>>> import pyirena.api.control as ctrl
>>> sid = ctrl.open_dataset("/data/sample.h5")["session_id"]
>>> ctrl.list_simple_models()                       # 13 models + Invariant
>>> ctrl.select_simple_model(sid, "Guinier")
>>> ctrl.set_simple_parameter(sid, "Rg", 150.0)
>>> ctrl.fix_simple_parameter(sid, "I0")            # hold I0, fit the rest
>>> ctrl.set_fit_q_range(sid, 0.002, 0.03)          # shared tool
>>> ctrl.run_simple_fit(sid)
>>> ctrl.get_simple_results(sid)
>>> ctrl.get_simple_fit_image(sid)
>>> ctrl.save_simple_fit(sid)

All functions return plain dicts.  Errors are
``{"error": ..., "code": ..., "suggestion": ...}`` dicts rather than exceptions.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from pyirena.api._paths import PathSecurityError, resolve_safe
from pyirena.api.control._images import render_png
from pyirena.api.control.errors import make_error, no_fit, no_session
from pyirena.api.control.session import fit_mask, get_session

__all__ = [
    "list_simple_models",
    "select_simple_model",
    "get_simple_config",
    "get_simple_parameters",
    "set_simple_parameter",
    "set_simple_parameter_bounds",
    "fix_simple_parameter",
    "free_simple_parameter",
    "reset_simple_parameters",
    "set_simple_background",
    "run_simple_fit",
    "get_simple_results",
    "get_simple_fit_image",
    "get_simple_linearization_image",
    "save_simple_fit",
]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _require_simple(session_id: str):
    """Return (session, None) when a Simple Fits model is ready, else (None, error)."""
    s = get_session(session_id)
    if s is None:
        return None, no_session(session_id)
    if s.model is None or s.model_name != "simple_fits":
        return None, make_error(
            f"Session '{session_id}' has no Simple Fits model selected.",
            suggestion="Call select_simple_model() first.",
            code="NO_SIMPLE_MODEL",
        )
    return s, None


def _fixed_set(session) -> set:
    """The agent's held-fixed parameter names for this session."""
    if not hasattr(session, "_simple_fixed") or session._simple_fixed is None:
        session._simple_fixed = set()
    return session._simple_fixed


def _bad_param(name: str, model) -> dict:
    return make_error(
        f"Parameter '{name}' not found in model '{model.model}'.",
        suggestion=(
            "Call get_simple_parameters() for the valid names of the selected "
            "model; each model has its own parameter set."
        ),
        code="BAD_PARAM",
    )


def _param_table(model, fixed: set) -> list:
    """Parameters as a list of dicts, in the model's own order."""
    rows = []
    for name, value in model.params.items():
        lo, hi = model.limits.get(name, (None, None))
        rows.append({
            "name": name,
            "value": float(value),
            "fixed": name in fixed,
            "lo": (None if lo is None else float(lo)),
            "hi": (None if hi is None else float(hi)),
        })
    return rows


def _config_dict(model, fixed: set) -> dict:
    """Full current configuration of a SimpleFitModel as a flat dict."""
    return {
        "model": str(model.model),
        "is_calculation": bool(model.is_calculation),
        "use_complex_bg": bool(model.use_complex_bg),
        "n_mc_runs": int(model.n_mc_runs),
        "invariant_porod_tail": bool(model.invariant_porod_tail),
        "use_slit_smearing": bool(model.use_slit_smearing),
        "slit_length": float(model.slit_length),
        "parameters": _param_table(model, fixed),
    }


def _render_image(fig, session_id: str, tag: str, dpi: int) -> tuple[str, str]:
    """Save the figure; return ``(base64 PNG, absolute path on disk)``."""
    return render_png(fig, f"simple_{tag}_{session_id}", dpi)


# ---------------------------------------------------------------------------
# Model lifecycle
# ---------------------------------------------------------------------------

def list_simple_models() -> dict:
    """List the available Simple Fits models and what each is for.

    No session needed — call this before ``select_simple_model`` to choose.

    Returns
    -------
    dict with ``models``: a list of ``{name, params, linearizable,
    supports_complex_background, is_calculation}``.  ``is_calculation`` marks
    models (Invariant) that are computed directly rather than least-squares
    fitted — for those, ``run_simple_fit`` evaluates instead of optimising.
    """
    from pyirena.core.simple_fits import MODEL_NAMES, MODEL_REGISTRY  # noqa: PLC0415

    models = []
    for name in MODEL_NAMES:
        entry = MODEL_REGISTRY[name]
        models.append({
            "name": name,
            "params": [p[0] for p in entry["params"]],
            "linearizable": entry.get("linearization") is not None,
            "supports_complex_background": bool(entry.get("complex_bg", False)),
            "is_calculation": bool(entry.get("calculation", False)),
        })
    return {"ok": True, "models": models}


def select_simple_model(session_id: str, model_name: str = "Guinier") -> dict:
    """Create a Simple Fits model for the session.

    Parameters
    ----------
    model_name : str
        One of the names from :func:`list_simple_models`.  Legacy spellings
        from before 1.1.0 are accepted and mapped transparently.

    Notes
    -----
    Calling this again — or switching model — resets every parameter to the
    new model's registry defaults, clears the held-fixed set, and discards any
    previous fit, because parameter names differ between models.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    from pyirena.core.simple_fits import MODEL_NAMES, SimpleFitModel  # noqa: PLC0415

    model = SimpleFitModel()
    try:
        model.set_model(model_name)
    except ValueError:
        return make_error(
            f"Unknown model '{model_name}'.",
            suggestion=f"Choose one of {list(MODEL_NAMES)}.",
            code="BAD_MODEL",
        )

    # Slit-smeared data: smear the analytic model before comparison so fitted
    # parameters stay ideal-space (parity with the GUI and batch layers).
    if s.is_slit_smeared and s.slit_length > 0:
        model.use_slit_smearing = True
        model.slit_length = float(s.slit_length)

    s.model = model
    s.model_name = "simple_fits"
    s.last_fit_result = None
    s._simple_fixed = set()

    return {"ok": True, "model": model.model,
            "config": _config_dict(model, _fixed_set(s))}


def get_simple_config(session_id: str) -> dict:
    """Return the current Simple Fits configuration (model, parameters, background)."""
    s, err = _require_simple(session_id)
    if err:
        return err
    return {"ok": True, "config": _config_dict(s.model, _fixed_set(s))}


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

def get_simple_parameters(session_id: str) -> dict:
    """Return every parameter of the selected model with value, bounds and fit state."""
    s, err = _require_simple(session_id)
    if err:
        return err
    return {
        "ok": True,
        "model": s.model.model,
        "parameters": _param_table(s.model, _fixed_set(s)),
    }


def set_simple_parameter(session_id: str, name: str, value: float) -> dict:
    """Set one parameter's value (its starting point for the fit)."""
    s, err = _require_simple(session_id)
    if err:
        return err
    if name not in s.model.params:
        return _bad_param(name, s.model)
    try:
        s.model.params[name] = float(value)
    except (TypeError, ValueError):
        return make_error(
            f"Value '{value}' for '{name}' is not a number.",
            suggestion="Pass a numeric value.",
            code="BAD_VALUE",
        )
    return {"ok": True, "name": name, "value": float(s.model.params[name])}


def set_simple_parameter_bounds(
    session_id: str,
    name: str,
    lo: Optional[float] = None,
    hi: Optional[float] = None,
) -> dict:
    """Set the fitting bounds of one parameter.

    Pass ``None`` for either side to leave it unbounded.  A value outside the
    new bounds is clamped into range, because an out-of-bounds start point
    makes the fit fail immediately.
    """
    s, err = _require_simple(session_id)
    if err:
        return err
    if name not in s.model.params:
        return _bad_param(name, s.model)

    lo_f = None if lo is None else float(lo)
    hi_f = None if hi is None else float(hi)
    if lo_f is not None and hi_f is not None and lo_f >= hi_f:
        return make_error(
            f"Lower bound {lo_f} is not below upper bound {hi_f}.",
            suggestion="Pass lo < hi.",
            code="BAD_BOUNDS",
        )

    s.model.limits[name] = (lo_f, hi_f)

    value = float(s.model.params[name])
    if lo_f is not None and value < lo_f:
        s.model.params[name] = lo_f
    elif hi_f is not None and value > hi_f:
        s.model.params[name] = hi_f

    return {
        "ok": True, "name": name, "lo": lo_f, "hi": hi_f,
        "value": float(s.model.params[name]),
    }


def fix_simple_parameter(session_id: str, name: str) -> dict:
    """Hold one parameter fixed at its current value during the fit."""
    s, err = _require_simple(session_id)
    if err:
        return err
    if name not in s.model.params:
        return _bad_param(name, s.model)
    _fixed_set(s).add(name)
    return {"ok": True, "name": name, "fixed": True,
            "value": float(s.model.params[name])}


def free_simple_parameter(session_id: str, name: str) -> dict:
    """Let one parameter vary during the fit (the default for every parameter)."""
    s, err = _require_simple(session_id)
    if err:
        return err
    if name not in s.model.params:
        return _bad_param(name, s.model)
    _fixed_set(s).discard(name)
    return {"ok": True, "name": name, "fixed": False}


def reset_simple_parameters(session_id: str) -> dict:
    """Reset every parameter to the model's registry defaults and free them all."""
    s, err = _require_simple(session_id)
    if err:
        return err
    s.model.set_model(s.model.model)     # re-applies registry defaults
    s._simple_fixed = set()
    s.last_fit_result = None
    return {"ok": True, "config": _config_dict(s.model, _fixed_set(s))}


# ---------------------------------------------------------------------------
# Complex background
# ---------------------------------------------------------------------------

def set_simple_background(session_id: str, enabled: bool = True) -> dict:
    """Enable or disable the complex background (power law + flat) for this model.

    Enabling adds ``BG_B``, ``BG_P`` and ``BG_flat`` to the parameter list;
    disabling removes them.  Models with an explicit Background parameter
    (Porod, Power Law) do not support it and return an error.
    """
    s, err = _require_simple(session_id)
    if err:
        return err

    from pyirena.core.simple_fits import MODEL_REGISTRY  # noqa: PLC0415

    entry = MODEL_REGISTRY[s.model.model]
    if enabled and not entry.get("complex_bg", False):
        return make_error(
            f"Model '{s.model.model}' has its own Background parameter and does "
            "not use the complex background.",
            suggestion="Set the model's Background parameter instead.",
            code="NO_COMPLEX_BG",
        )

    s.model.use_complex_bg = bool(enabled)
    # Re-derive params so the BG_* entries appear or disappear, keeping the
    # values the agent has already set for the model's own parameters.
    kept = dict(s.model.params)
    s.model.set_model(s.model.model)
    for name, value in kept.items():
        if name in s.model.params:
            s.model.params[name] = value
    _fixed_set(s).intersection_update(set(s.model.params))

    return {"ok": True, "use_complex_bg": s.model.use_complex_bg,
            "parameters": _param_table(s.model, _fixed_set(s))}


# ---------------------------------------------------------------------------
# Fit execution
# ---------------------------------------------------------------------------

def run_simple_fit(session_id: str, no_limits: bool = False) -> dict:
    """Fit the selected model to the data inside the current fit Q range.

    Parameters
    ----------
    no_limits : bool
        Fit unconstrained, ignoring every parameter bound.  Useful when a fit
        stalls against a bound that was only a guess.

    Returns
    -------
    dict with success, model, chi_squared, reduced_chi_squared, dof,
    q_min/q_max, n_data, parameters (value + 1σ) and derived quantities.
    Calculation models (Invariant) are evaluated rather than optimised.
    """
    s, err = _require_simple(session_id)
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

    fixed = _fixed_set(s)
    free = [n for n in s.model.params if n not in fixed]
    if not free and not s.model.is_calculation:
        return make_error(
            "Every parameter is held fixed — nothing to fit.",
            suggestion="Call free_simple_parameter() for at least one parameter.",
            code="ALL_FIXED",
        )

    fixed_params = {n: float(s.model.params[n]) for n in fixed}

    try:
        result = s.model.fit(q_fit, I_fit, err_fit,
                             fixed_params=fixed_params or None,
                             no_limits=bool(no_limits))
    except Exception as exc:  # pragma: no cover - defensive
        return make_error(
            f"Fit failed with exception: {exc}",
            suggestion="Check the parameter values, bounds and fit Q range.",
            code="FIT_EXCEPTION",
        )

    if not result.get("success", False):
        return make_error(
            result.get("error") or "Fit did not converge.",
            suggestion="Try different starting values, or widen the bounds.",
            code="FIT_FAILED",
        )

    # Fitted values become the new starting point, as in the GUI.
    s.model.params.update(result.get("params", {}))
    s.last_fit_result = result

    return {
        "success": True,
        "model": result.get("model"),
        "chi_squared": _f(result.get("chi2")),
        "reduced_chi_squared": _f(result.get("reduced_chi2")),
        "dof": (int(result["dof"]) if result.get("dof") is not None else None),
        "q_min": float(np.min(q_fit)),
        "q_max": float(np.max(q_fit)),
        "n_data": int(len(q_fit)),
        "parameters": [
            {
                "name": name,
                "value": _f(value),
                "std": _f((result.get("params_std") or {}).get(name)),
                "fixed": name in fixed,
            }
            for name, value in (result.get("params") or {}).items()
        ],
        "derived": {k: _f(v) for k, v in (result.get("derived") or {}).items()},
        "warning": result.get("warning"),
    }


def _f(value):
    """Float, or None for missing/non-finite — JSON has no NaN."""
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

def get_simple_results(session_id: str) -> dict:
    """Return the last fit's parameters, uncertainties, quality and derived values."""
    s, err = _require_simple(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    res = s.last_fit_result
    fixed = _fixed_set(s)
    return {
        "ok": True,
        "model": res.get("model"),
        "chi_squared": _f(res.get("chi2")),
        "reduced_chi_squared": _f(res.get("reduced_chi2")),
        "dof": (int(res["dof"]) if res.get("dof") is not None else None),
        "use_complex_bg": bool(s.model.use_complex_bg),
        "parameters": [
            {
                "name": name,
                "value": _f(value),
                "std": _f((res.get("params_std") or {}).get(name)),
                "fixed": name in fixed,
            }
            for name, value in (res.get("params") or {}).items()
        ],
        "derived": {k: _f(v) for k, v in (res.get("derived") or {}).items()},
    }


def get_simple_fit_image(
    session_id: str, width: int = 1024, height: int = 800, dpi: int = 120,
) -> dict:
    """Render the fit as a two-panel PNG: log-log data + model, and residuals.

    Returns ``{"image_base64": ...}``.
    """
    s, err = _require_simple(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt  # noqa: PLC0415

    res = s.last_fit_result
    q = np.asarray(res["q"], dtype=float)
    model_I = np.asarray(res["I_model"], dtype=float)
    residuals = np.asarray(res.get("residuals"), dtype=float)

    mask = fit_mask(s)
    I_obs = s.intensity[mask]
    err_obs = s.error[mask] if s.error is not None else None

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(width / dpi, height / dpi), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    if err_obs is not None:
        ax1.errorbar(q, I_obs, yerr=err_obs, fmt="o", markersize=3, capsize=2,
                     alpha=0.6, color="steelblue", label="Data", linewidth=0.5)
    else:
        ax1.plot(q, I_obs, "o", markersize=3, alpha=0.6, color="steelblue",
                 label="Data")
    ax1.plot(q, model_I, "-", linewidth=2, color="red", label=res.get("model", "Model"))
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_ylabel("I  (cm⁻¹)", fontsize=10)
    chi2r = _f(res.get("reduced_chi2"))
    title = f"Simple Fit — {res.get('model', '')}"
    if chi2r is not None:
        title += f"   reduced χ² = {chi2r:.3g}"
    ax1.set_title(title, fontsize=11, fontweight="bold")
    ax1.legend(fontsize=8, loc="best")
    ax1.grid(True, which="both", alpha=0.25)

    ax2.axhline(0.0, color="black", linewidth=1, linestyle="--")
    ax2.plot(q, residuals, "o", markersize=3, color="steelblue")
    ax2.set_xscale("log")
    ax2.set_xlabel("Q  (Å⁻¹)", fontsize=10)
    ax2.set_ylabel("Residuals", fontsize=10)
    ax2.grid(True, which="both", alpha=0.25)

    fig.tight_layout()
    b64, image_path = _render_image(fig, session_id, "fit", dpi)
    return {"ok": True, "image_base64": b64, "image_path": image_path}


def get_simple_linearization_image(
    session_id: str, width: int = 900, height: int = 700, dpi: int = 120,
) -> dict:
    """Render the model's linearized plot (Guinier plot, Porod plot, …) as a PNG.

    A straight line in the linearized view is the classic visual test that the
    model applies over the chosen Q range — worth checking before trusting a
    fitted Rg.  Models with no linearization return an error saying so.
    """
    s, err = _require_simple(session_id)
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

    lin = s.model.linearize(q_fit, I_fit, err_fit)
    if lin is None:
        return make_error(
            f"Model '{s.model.model}' has no linearized form.",
            suggestion="Use get_simple_fit_image() instead.",
            code="NO_LINEARIZATION",
        )

    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt  # noqa: PLC0415

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi))
    ax.plot(lin["X"], lin["Y"], "o", markersize=4, alpha=0.7, color="steelblue",
            label="Linearized data")
    ax.plot(lin["X_fit"], lin["Y_fit"], "-", linewidth=2, color="red",
            label="Linear fit")
    ax.set_xlabel(lin.get("x_label", "X"), fontsize=10)
    ax.set_ylabel(lin.get("y_label", "Y"), fontsize=10)
    ax.set_title(
        f"{s.model.model} linearization — "
        f"slope {lin['slope']:.4g}, intercept {lin['intercept']:.4g}, "
        f"R² {lin['r_squared']:.4f}",
        fontsize=10, fontweight="bold",
    )
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    _lin_b64, _lin_path = _render_image(fig, session_id, "lin", dpi)
    return {
        "ok": True,
        "slope": _f(lin.get("slope")),
        "intercept": _f(lin.get("intercept")),
        "r_squared": _f(lin.get("r_squared")),
        "image_base64": _lin_b64,
        "image_path": _lin_path,
    }


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def save_simple_fit(session_id: str, output_path: Optional[str] = None) -> dict:
    """Save the fit to a NXcanSAS HDF5 file under ``entry/simple_fit_results``.

    Parameters
    ----------
    output_path : str, optional
        Where to save.  Defaults to the original file (in-place).  Saving to a
        different path preserves the original.
    """
    s, err = _require_simple(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    from pyirena.io.nxcansas_simple_fits import save_simple_fit_results  # noqa: PLC0415

    try:
        src = resolve_safe(s.file_path, must_exist=False)
        target = resolve_safe(output_path, must_exist=False) if output_path else src
    except PathSecurityError as exc:
        return make_error(
            str(exc),
            suggestion="Save to a path inside PYIRENA_DATA_ROOT.",
            code="PATH_NOT_ALLOWED",
        )

    # Saving to a *new* location must yield a complete, re-openable NXcanSAS
    # file — not a results-only stub.
    if target != src and not target.exists():
        from pyirena.io._nxcansas_common import copy_and_strip_results  # noqa: PLC0415
        try:
            copy_and_strip_results(src, target)
        except Exception as exc:
            return make_error(
                f"Could not create output file '{target}' from source: {exc}",
                suggestion="Check the source file exists and the target is writable.",
                code="SAVE_ERROR",
            )

    mask = fit_mask(s)
    I_data = s.intensity[mask]
    I_err = s.error[mask] if s.error is not None else None

    # The same setup dict the GUI embeds, so a saved agent run reopens in the
    # panel with every control restored.
    setup_state = s.model.to_dict()
    setup_state["fixed_params"] = sorted(_fixed_set(s))
    setup_state["q_min"] = (None if s.fit_q_min is None else float(s.fit_q_min))
    setup_state["q_max"] = (None if s.fit_q_max is None else float(s.fit_q_max))

    try:
        save_simple_fit_results(
            filepath=target,
            result=s.last_fit_result,
            model_obj=s.model,
            intensity_data=I_data,
            intensity_error=I_err,
            setup_state=setup_state,
        )
    except Exception as exc:
        return make_error(
            f"Could not save results: {exc}",
            suggestion="Check the file is writable and not open elsewhere.",
            code="SAVE_ERROR",
        )

    return {"ok": True, "saved_to": str(target), "group": "entry/simple_fit_results"}
