"""Control-surface tools for Unified Fit model.

These functions give an LLM agent (or any automation script) the ability to:
  - Open a dataset and create a session
  - Select and configure the Unified Fit model
  - Read and modify per-parameter values, bounds, and fix/free status
  - Restrict the Q range used for fitting
  - Run fits and read results
  - Capture the current fit as a PNG image (base64-encoded)
  - Save results back to the HDF5 file

All public functions:
  - Accept plain JSON types (str, float, int, bool, list)
  - Return a JSON-serialisable dict
  - Never raise across the API boundary — errors are {"error": ..., ...} dicts
  - Are stateless themselves; state lives in the Session registry

Session model
-------------
Call open_dataset() to get a session_id.  Pass that id to every other call.
Multiple sessions can be open simultaneously.  Sessions are in-memory only
and are lost when the process exits.
"""
from __future__ import annotations

import base64
import os
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np

from pyirena.api.control.errors import (
    bad_param, make_error, no_fit, no_model, no_session,
)
from pyirena.api.control.session import (
    Session, all_sessions, create_session, drop_session, fit_mask, get_session,
)
from pyirena.core.fit_metrics import fit_quality_metrics

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Parameters that belong to each UnifiedLevel and their metadata.
# Order matters — it determines the parameter table order.
_LEVEL_PARAMS = [
    ("Rg",   "Å",             "Radius of gyration"),
    ("G",    "cm⁻¹",          "Guinier prefactor (low-Q amplitude)"),
    ("P",    "",              "Power law slope"),
    ("B",    "cm⁻¹·Å⁻ᴾ",     "Power law prefactor (Porod constant)"),
    ("ETA",  "Å",             "Correlation distance (Born-Green)"),
    ("PACK", "",              "Packing factor for correlations"),
    ("RgCO", "Å",             "Cutoff radius (high-Q exponential decay)"),
]


def _param_table(model) -> list[dict]:
    """Convert a UnifiedFitModel to a JSON-serialisable parameter table."""
    rows = []
    for i, level in enumerate(model.levels):
        level_num = i + 1
        for attr, units, desc in _LEVEL_PARAMS:
            value = getattr(level, attr)
            fit_flag = getattr(level, f"fit_{attr}", None)
            # fit_flag is None only for attrs that can't be fitted (K); skip those
            if fit_flag is None:
                continue
            limits = getattr(level, f"{attr}_limits", (None, None))
            rows.append({
                "name":        f"{attr}_{level_num}",
                "value":       float(value),
                "fixed":       not bool(fit_flag),
                "lo":          float(limits[0]) if limits[0] is not None else None,
                "hi":          float(limits[1]) if limits[1] is not None else None,
                "units":       units,
                "description": f"{desc} (level {level_num})",
            })
    # Model-wide background
    rows.append({
        "name":        "background",
        "value":       float(model.background),
        "fixed":       not bool(model.fit_background),
        "lo":          float(model.background_limits[0]),
        "hi":          float(model.background_limits[1]),
        "units":       "cm⁻¹",
        "description": "Flat incoherent background",
    })
    return rows


def _find_param(model, param_name: str):
    """Return (level_index_or_None, attr_name) for a parameter name, or None."""
    if param_name == "background":
        return ("background", None)
    for attr, _, _ in _LEVEL_PARAMS:
        for i in range(model.num_levels):
            if param_name == f"{attr}_{i + 1}":
                return (i, attr)
    return None


# Per-level boolean flags that materially change the level's intensity formula.
# These are *not* numeric parameters and cannot be set with set_parameter_value.
_LEVEL_OPTIONS = ("correlations", "mass_fractal", "link_B", "link_RGCO")


def _level_options(model) -> list[dict]:
    """Per-level boolean flag state — separate from numeric parameter rows."""
    return [
        {
            "level":        i + 1,
            "correlations": bool(lv.correlations),
            "mass_fractal": bool(lv.mass_fractal),
            "link_B":       bool(lv.link_B),
            "link_RGCO":    bool(lv.link_RGCO),
        }
        for i, lv in enumerate(model.levels)
    ]


def _decimated(arr: np.ndarray, max_points: int = 500) -> list:
    """Thin a numpy array to at most max_points for JSON transport."""
    n = len(arr)
    if n <= max_points:
        return [v if np.isfinite(v) else None for v in arr.tolist()]
    step = max(1, n // max_points)
    sub = arr[::step]
    return [v if np.isfinite(v) else None for v in sub.tolist()]


def _render_fit_image(session: Session, width: int, height: int, dpi: int = 120) -> str:
    """Render data + model (+ residuals if fit was run) as base64 PNG.

    Matches the visual layout of pyirena's UnifiedFitResultsWindow:
      - Top panel (log-log): data points + model line (+ per-level dashed)
      - Bottom panel (log-x linear-y): normalised residuals (only after fit)
    """
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    model = session.model
    has_fit = session.last_fit_result is not None and model.fit_intensity is not None

    # Choose layout
    if has_fit:
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(width / dpi, height / dpi),
            gridspec_kw={"height_ratios": [3, 1]},
        )
    else:
        fig, ax1 = plt.subplots(1, 1, figsize=(width / dpi, height / dpi))
        ax2 = None

    q = session.q
    I = session.intensity
    err = session.error

    # --- top panel ---
    if err is not None:
        ax1.errorbar(q, I, yerr=err, fmt="o", markersize=3, capsize=2,
                     alpha=0.6, label="Data", color="steelblue", linewidth=0.5)
    else:
        ax1.plot(q, I, "o", markersize=3, alpha=0.6, label="Data", color="steelblue")

    if model is not None:
        # Current model at full Q range (calculate_intensity works before fit too)
        model_I = model.calculate_intensity(q)
        ax1.plot(q, model_I, "-", linewidth=2, color="red", label="Model")

        # Individual levels (dashed)
        level_colors = ["#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#17becf"]
        prev_Rg = 0.0
        for i in range(model.num_levels):
            level_I = model.calculate_level_intensity(q, i, prev_Rg)
            level_I += model.background / model.num_levels
            ax1.plot(q, level_I, "--", linewidth=1.2, alpha=0.7,
                     color=level_colors[i % len(level_colors)],
                     label=f"Level {i + 1}")
            prev_Rg = model.levels[i].Rg

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Q  (Å⁻¹)", fontsize=10)
    ax1.set_ylabel("I  (cm⁻¹)", fontsize=10)
    chi_sq = (session.last_fit_result or {}).get("reduced_chi_squared")
    title = "Unified Fit"
    if chi_sq is not None:
        title += f"   χ²ᵣ = {chi_sq:.3f}"
    ax1.set_title(title, fontsize=11, fontweight="bold")
    ax1.legend(fontsize=8, loc="best")
    ax1.grid(True, which="both", alpha=0.25)

    # --- residuals panel (only after fit) ---
    if has_fit and ax2 is not None:
        resid = model.I_data - model.fit_intensity
        if model.error_data is not None and np.any(model.error_data > 0):
            y = resid / model.error_data
            ylabel = "Normalised residuals"
        else:
            y = resid
            ylabel = "Residuals  (cm⁻¹)"
        ax2.plot(model.q_data, y, "o", markersize=3, color="black", alpha=0.7)
        ax2.axhline(0, color="red", linestyle="--", linewidth=1)
        ax2.set_xscale("log")
        ax2.set_xlabel("Q  (Å⁻¹)", fontsize=10)
        ax2.set_ylabel(ylabel, fontsize=10)
        ax2.grid(True, which="both", alpha=0.25)

    fig.tight_layout()

    # Save to temp file and return base64
    tmp = Path(tempfile.gettempdir()) / "pyirena-ctrl"
    tmp.mkdir(parents=True, exist_ok=True)
    out = tmp / f"fit_{session.session_id}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(out.read_bytes()).decode("ascii")


# ---------------------------------------------------------------------------
# Category A — Session lifecycle
# ---------------------------------------------------------------------------

def open_dataset(file_path: str) -> dict:
    """Load a NXcanSAS HDF5 file and create a new fitting session.

    Returns
    -------
    dict with keys: session_id, summary (file, n_points, q_min, q_max,
    intensity_min, intensity_max, label).
    On failure returns an error dict.
    """
    from pyirena.io.hdf5 import readGenericNXcanSAS

    fp = Path(file_path)
    if not fp.exists():
        return make_error(
            f"File not found: '{file_path}'",
            suggestion="Check the path and try again.",
            code="FILE_NOT_FOUND",
        )

    try:
        data = readGenericNXcanSAS(str(fp.parent), fp.name)
    except Exception as exc:
        return make_error(
            f"Could not read '{file_path}': {exc}",
            code="READ_ERROR",
        )

    if data is None:
        return make_error(
            f"No readable SAS data found in '{file_path}'.",
            suggestion="Verify the file is a valid NXcanSAS HDF5 file.",
            code="NO_DATA",
        )

    q_raw = data.get("Q")
    I_raw = data.get("Intensity")
    q = np.asarray(q_raw if q_raw is not None else [], dtype=float)
    I = np.asarray(I_raw if I_raw is not None else [], dtype=float)
    err_raw = data.get("Error")
    err = np.asarray(err_raw, dtype=float) if err_raw is not None else None

    if len(q) == 0 or len(I) == 0:
        return make_error(
            f"Empty Q or Intensity array in '{file_path}'.",
            code="EMPTY_DATA",
        )

    label = data.get("label") or ""
    if isinstance(label, bytes):
        label = label.decode("utf-8", errors="replace")

    session = create_session(file_path=str(fp), q=q, intensity=I, error=err,
                             label=str(label))

    valid = np.isfinite(q) & np.isfinite(I)
    return {
        "session_id": session.session_id,
        "summary": {
            "file":          str(fp),
            "label":         session.label,
            "n_points":      int(len(q)),
            "q_min":         float(np.nanmin(q[valid])),
            "q_max":         float(np.nanmax(q[valid])),
            "intensity_min": float(np.nanmin(I[valid & (I > 0)])) if np.any(valid & (I > 0)) else None,
            "intensity_max": float(np.nanmax(I[valid])),
        },
    }


def list_open_sessions() -> dict:
    """Return a summary of all open sessions."""
    sessions = all_sessions()
    return {
        "sessions": [
            {
                "session_id": s.session_id,
                "file":       s.file_path,
                "label":      s.label,
                "model":      s.model_name,
                "has_fit":    s.last_fit_result is not None,
            }
            for s in sessions
        ],
        "count": len(sessions),
    }


def close_session(session_id: str) -> dict:
    """Close and discard a session, freeing its memory."""
    if not drop_session(session_id):
        return no_session(session_id)
    return {"ok": True, "session_id": session_id}


def get_session_summary(session_id: str) -> dict:
    """Return a summary of the current session state."""
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    valid = np.isfinite(s.q) & np.isfinite(s.intensity)
    chi = None
    if s.last_fit_result:
        chi = s.last_fit_result.get("reduced_chi_squared")

    return {
        "session_id":    session_id,
        "file":          s.file_path,
        "label":         s.label,
        "n_points":      int(len(s.q)),
        "q_min":         float(np.nanmin(s.q[valid])),
        "q_max":         float(np.nanmax(s.q[valid])),
        "model":         s.model_name,
        "nlevels":       s.model.num_levels if s.model is not None else None,
        "fit_q_min":     s.fit_q_min,
        "fit_q_max":     s.fit_q_max,
        "has_fit":       s.last_fit_result is not None,
        "reduced_chi_squared": chi,
    }


# ---------------------------------------------------------------------------
# Category B — Model selection and inspection
# ---------------------------------------------------------------------------

AVAILABLE_MODELS = ["unified_fit"]


def list_available_models() -> dict:
    """List the fitting models available in this version of the control API."""
    return {
        "models": AVAILABLE_MODELS,
        "note": "Additional models will be added in future phases.",
    }


def select_model(
    session_id: str,
    model_name: str = "unified_fit",
    nlevels: int = 1,
) -> dict:
    """Select and initialise a fitting model for the session.

    For unified_fit, nlevels sets the number of structural levels (1-5).
    Calling select_model() again replaces the model and clears any prior fit.

    Parameters
    ----------
    nlevels : int
        Number of Unified Fit levels (1-5).  Start with 1 and add levels
        as needed via add_unified_level().
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    if model_name not in AVAILABLE_MODELS:
        return make_error(
            f"Unknown model '{model_name}'.",
            suggestion=f"Available models: {AVAILABLE_MODELS}",
            code="UNKNOWN_MODEL",
        )

    if model_name == "unified_fit":
        if not (1 <= nlevels <= 5):
            return make_error(
                f"nlevels must be between 1 and 5, got {nlevels}.",
                code="BAD_NLEVELS",
            )
        from pyirena.core.unified import UnifiedFitModel
        s.model = UnifiedFitModel(num_levels=nlevels)
        s.model_name = "unified_fit"
        s.last_fit_result = None

    return {
        "ok":            True,
        "model":         s.model_name,
        "nlevels":       s.model.num_levels,
        "parameters":    _param_table(s.model),
        "level_options": _level_options(s.model),
    }


def get_model_parameters(session_id: str) -> dict:
    """Return the current parameter table for the session's model.

    Returns
    -------
    dict with keys:
      - parameters    : list of {name, value, fixed, lo, hi, units, description}
      - level_options : list of per-level boolean flags
                        {level, correlations, mass_fractal, link_B, link_RGCO}

    The numeric parameters live in `parameters`; the per-level boolean flags
    (which switch entire intensity-formula features on/off) live in
    `level_options`.  Use set_level_option() to toggle the booleans —
    set_parameter_value() only works on numeric parameters.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)

    return {
        "model":         s.model_name,
        "nlevels":       s.model.num_levels,
        "parameters":    _param_table(s.model),
        "level_options": _level_options(s.model),
    }


def get_model_description(session_id: str) -> dict:
    """Return a text description of the Unified Fit model and its parameters.

    Useful as context for a system prompt or to understand what each parameter
    physically means before setting values.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)

    return {
        "model": "unified_fit",
        "summary": (
            "The Unified Fit model (Beaucage, 1995/1996) describes hierarchical "
            "structures in SAXS/USAXS data.  Each structural level contributes a "
            "Guinier term (low-Q) and a power-law term (high-Q).  Multiple levels "
            "are summed to represent multi-scale structures."
        ),
        "parameter_guide": {
            "Rg":         "Radius of gyration [Å].  Controls where the Guinier knee appears.  Start with a value roughly matching the feature size expected from the data.",
            "G":          "Guinier prefactor [cm⁻¹].  Controls the low-Q intensity of this level.  Often estimated from the data value at the Guinier knee.",
            "P":          "Power-law slope.  P=4 → smooth surfaces (Porod).  P=3 → rough surfaces.  P<3 → mass fractal.  Usually fixed at 4 initially.",
            "B":          "Power-law prefactor [cm⁻¹·Å⁻ᴾ].  Rarely needs manual setting; often linked to G via the link_B level option.",
            "ETA":        "Correlation distance [Å].  HAS NO EFFECT unless the level's `correlations` option is enabled — use set_level_option(level, 'correlations', True).",
            "PACK":       "Packing factor.  HAS NO EFFECT unless the level's `correlations` option is enabled — use set_level_option(level, 'correlations', True).",
            "RgCO":       "Cutoff radius [Å].  High-Q exponential roll-off.  Useful when this level is bounded above by the next level's Rg.  Often linked to previous level's Rg via the link_RGCO level option.",
            "background": "Flat incoherent background [cm⁻¹].  Fit this first; fix once converged.",
        },
        "level_options_guide": {
            "_note": (
                "Per-level boolean flags that switch entire features of the "
                "intensity formula on or off.  These are NOT numeric parameters "
                "and cannot be changed with set_parameter_value() — use "
                "set_level_option(session_id, level, option, enabled) instead.  "
                "Read the current state from get_model_parameters()['level_options']."
            ),
            "correlations": (
                "Enable the Born-Green liquid-like-ordering correction for this level. "
                "Must be True for ETA and PACK to have any effect — without this flag, "
                "setting ETA/PACK does NOTHING.  Use only for concentrated/interacting "
                "systems where a clear interference peak is visible."
            ),
            "mass_fractal": (
                "Switch this level to mass-fractal mode.  B is auto-computed from G, "
                "Rg, and P; the manual B value is ignored."
            ),
            "link_B": (
                "Estimate B from G, Rg, P via the Porod invariant — reduces a free "
                "parameter.  Recommended for the smallest level when P ≈ 4."
            ),
            "link_RGCO": (
                "Link this level's RgCO to the previous level's Rg.  Use when the "
                "level represents an aggregate composed of the previous level's "
                "primary particles."
            ),
        },
        "fitting_tips": [
            "Start with one level.  Add more only if residuals show systematic structure.",
            "Fix P at 4 initially unless you have reason to believe otherwise.",
            "Fit Rg and G first (fix everything else) to get the gross shape right.",
            "Free background early — a wrong background biases all other parameters.",
            "The Q range matters: restrict if the data has a substrate or beam-stop artifact at the extremes.",
            "If chi_squared does not change after setting ETA/PACK values, the "
            "level's `correlations` option is OFF.  Toggle it on with "
            "set_level_option(session_id, level, 'correlations', True).",
        ],
    }


# ---------------------------------------------------------------------------
# Category C — Parameter control
# ---------------------------------------------------------------------------

def set_parameter_value(session_id: str, param_name: str, value: float) -> dict:
    """Set the starting/current value of a parameter.

    The change takes effect at the next run_fit() call.  If the new value lies
    outside the current bounds the bounds are silently expanded to contain it,
    so that a fixed parameter whose value is set to an extreme sentinel (e.g.
    G=0 to "remove" a level) does not cause scipy to reject the initial guess.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)

    found = _find_param(s.model, param_name)
    if found is None:
        return bad_param(param_name, s.model_name)

    v = float(value)
    level_idx, attr = found
    if level_idx == "background":
        s.model.background = v
        lo, hi = s.model.background_limits
        if v < lo or v > hi:
            s.model.background_limits = (min(lo, v), max(hi, v))
    else:
        setattr(s.model.levels[level_idx], attr, v)
        limits_attr = f"{attr}_limits"
        lo, hi = getattr(s.model.levels[level_idx], limits_attr)
        if v < lo or v > hi:
            setattr(s.model.levels[level_idx], limits_attr, (min(lo, v), max(hi, v)))

    return {"ok": True, "param": param_name, "value": v}


def set_parameter_bounds(
    session_id: str, param_name: str, lo: float, hi: float
) -> dict:
    """Set the lower and upper bounds for a parameter during fitting."""
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)

    found = _find_param(s.model, param_name)
    if found is None:
        return bad_param(param_name, s.model_name)

    lo_f, hi_f = float(lo), float(hi)
    if lo_f > hi_f:
        return make_error(
            f"lo ({lo_f}) must be ≤ hi ({hi_f}).",
            code="BAD_BOUNDS",
        )

    kind, idx = found
    if kind == "background":
        s.model.background_limits = (lo_f, hi_f)
    else:
        attr = param_name.rsplit("_", 1)[0]
        setattr(s.model.levels[kind], f"{attr}_limits", (lo_f, hi_f))

    return {"ok": True, "param": param_name, "lo": lo_f, "hi": hi_f}


def fix_parameter(session_id: str, param_name: str) -> dict:
    """Hold a parameter fixed at its current value during fitting."""
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)

    found = _find_param(s.model, param_name)
    if found is None:
        return bad_param(param_name, s.model_name)

    kind, idx = found
    if kind == "background":
        s.model.fit_background = False
    else:
        attr = param_name.rsplit("_", 1)[0]
        setattr(s.model.levels[kind], f"fit_{attr}", False)

    return {"ok": True, "param": param_name, "fixed": True}


def free_parameter(session_id: str, param_name: str) -> dict:
    """Allow a parameter to vary during fitting."""
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)

    found = _find_param(s.model, param_name)
    if found is None:
        return bad_param(param_name, s.model_name)

    kind, idx = found
    if kind == "background":
        s.model.fit_background = True
    else:
        attr = param_name.rsplit("_", 1)[0]
        setattr(s.model.levels[kind], f"fit_{attr}", True)

    return {"ok": True, "param": param_name, "fixed": False}


def fix_all_except(session_id: str, free_list: list[str]) -> dict:
    """Fix every parameter except those in free_list.

    Useful for staged fitting: fix everything, then release parameters one
    group at a time.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)

    # Validate free_list first
    for name in free_list:
        if _find_param(s.model, name) is None:
            return bad_param(name, s.model_name)

    all_params = [row["name"] for row in _param_table(s.model)]
    free_set = set(free_list)
    fixed_count = free_count = 0

    for pname in all_params:
        found = _find_param(s.model, pname)
        kind, idx = found
        if pname in free_set:
            if kind == "background":
                s.model.fit_background = True
            else:
                attr = pname.rsplit("_", 1)[0]
                setattr(s.model.levels[kind], f"fit_{attr}", True)
            free_count += 1
        else:
            if kind == "background":
                s.model.fit_background = False
            else:
                attr = pname.rsplit("_", 1)[0]
                setattr(s.model.levels[kind], f"fit_{attr}", False)
            fixed_count += 1

    return {
        "ok": True,
        "fixed_count": fixed_count,
        "free_count":  free_count,
        "free_params": list(free_set),
    }


def reset_parameters_to_defaults(session_id: str) -> dict:
    """Reset all parameters to their dataclass defaults and re-select the model.

    Equivalent to calling select_model() again with the same nlevels.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)

    nlevels = s.model.num_levels
    return select_model(session_id, model_name=s.model_name, nlevels=nlevels)


# ---------------------------------------------------------------------------
# Category C'' — Unified-specific level operations
# ---------------------------------------------------------------------------

def add_unified_level(session_id: str, position: int = -1) -> dict:
    """Add a structural level to the Unified Fit model.

    Parameters
    ----------
    position : int
        Where to insert the new level.  -1 (default) appends at the end.
        1 inserts before the current first level.  Uses 1-based indexing.

    Returns the updated parameter table so the agent can see the new names.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None or s.model_name != "unified_fit":
        return no_model(session_id)

    from pyirena.core.unified import UnifiedFitModel, UnifiedLevel

    current = s.model.num_levels
    if current >= 5:
        return make_error(
            "Cannot add a level: Unified Fit supports at most 5 levels.",
            code="MAX_LEVELS",
        )

    new_level = UnifiedLevel()
    levels = list(s.model.levels)

    if position == -1 or position >= current + 1:
        levels.append(new_level)
    else:
        insert_idx = max(0, position - 1)
        levels.insert(insert_idx, new_level)

    # Re-create model preserving settings
    new_model = UnifiedFitModel(num_levels=current + 1)
    new_model.levels = levels
    new_model.background = s.model.background
    new_model.fit_background = s.model.fit_background
    new_model.background_limits = s.model.background_limits
    s.model = new_model
    s.last_fit_result = None

    return {
        "ok":         True,
        "nlevels":    s.model.num_levels,
        "parameters": _param_table(s.model),
    }


def remove_unified_level(session_id: str, level: int) -> dict:
    """Remove a structural level (1-based index) from the Unified Fit model.

    After removal, levels are renumbered from 1.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None or s.model_name != "unified_fit":
        return no_model(session_id)

    from pyirena.core.unified import UnifiedFitModel

    current = s.model.num_levels
    if current <= 1:
        return make_error(
            "Cannot remove a level: model must have at least 1 level.",
            code="MIN_LEVELS",
        )
    if not (1 <= level <= current):
        return make_error(
            f"Level {level} does not exist (model has levels 1–{current}).",
            code="BAD_LEVEL",
        )

    levels = [lv for i, lv in enumerate(s.model.levels) if i != level - 1]
    new_model = UnifiedFitModel(num_levels=current - 1)
    new_model.levels = levels
    new_model.background = s.model.background
    new_model.fit_background = s.model.fit_background
    new_model.background_limits = s.model.background_limits
    s.model = new_model
    s.last_fit_result = None

    return {
        "ok":         True,
        "nlevels":    s.model.num_levels,
        "parameters": _param_table(s.model),
    }


# ---------------------------------------------------------------------------
# Category C''' — Per-level boolean options (correlations, mass_fractal, etc.)
# ---------------------------------------------------------------------------

def get_level_options(session_id: str, level: Optional[int] = None) -> dict:
    """Return per-level boolean flag state.

    Parameters
    ----------
    level : int, optional
        1-based level number.  If omitted, returns the state of all levels.

    Returns
    -------
    dict with key 'level_options' = list of
        {level, correlations, mass_fractal, link_B, link_RGCO}
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None or s.model_name != "unified_fit":
        return no_model(session_id)

    opts = _level_options(s.model)
    if level is not None:
        if not (1 <= level <= s.model.num_levels):
            return make_error(
                f"Level {level} does not exist (model has levels 1–{s.model.num_levels}).",
                code="BAD_LEVEL",
            )
        return {"level_options": [opts[level - 1]]}
    return {"level_options": opts}


def set_level_option(
    session_id: str,
    level: int,
    option: str,
    enabled: bool,
) -> dict:
    """Toggle a per-level boolean flag.

    These flags switch entire features of the intensity formula on or off.
    They are NOT numeric parameters and cannot be set with set_parameter_value().

    Parameters
    ----------
    level : int
        1-based level number to modify.
    option : str
        One of:
          - "correlations" — enable Born-Green liquid-like-ordering correction.
            **Required** for ETA and PACK to have any effect; without this
            flag set to True, ETA/PACK values are ignored.
          - "mass_fractal" — switch this level to mass-fractal mode (B is
            auto-computed from G, Rg, P; manual B is ignored).
          - "link_B" — estimate B from G, Rg, P via the Porod invariant.
          - "link_RGCO" — link this level's RgCO to the previous level's Rg.
    enabled : bool
        True to turn the option on, False to turn it off.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None or s.model_name != "unified_fit":
        return no_model(session_id)

    if not (1 <= level <= s.model.num_levels):
        return make_error(
            f"Level {level} does not exist (model has levels 1–{s.model.num_levels}).",
            code="BAD_LEVEL",
        )
    if option not in _LEVEL_OPTIONS:
        return make_error(
            f"Unknown level option '{option}'.",
            suggestion=f"Valid options: {list(_LEVEL_OPTIONS)}",
            code="BAD_OPTION",
        )

    setattr(s.model.levels[level - 1], option, bool(enabled))
    return {
        "ok":      True,
        "level":   level,
        "option":  option,
        "enabled": bool(enabled),
    }


def check_level_feasibility(
    session_id: str, level: Optional[int] = None
) -> dict:
    """Check whether each level's parameters are physically meaningful.

    A level is "feasible" when its Guinier and power-law regions connect
    smoothly at the Hammouda rollover Q point.  A large discontinuity at
    that point indicates an unphysical (G, Rg, B, P) combination — for
    example, a power-law tail that lies far above or below the Guinier
    plateau extrapolated to the same Q.

    Use this after run_fit to catch combinations that converged
    mathematically but are not physically interpretable.  If a level is
    not feasible, common culprits are B too small/large for the chosen P,
    or G inconsistent with B and Rg.

    Parameters
    ----------
    level : int, optional
        1-based level number.  If omitted, checks every level.

    Returns
    -------
    dict with keys:
      - feasibility : list of {level, feasible (bool)} for each checked level
      - all_feasible : True iff every checked level is feasible
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None or s.model_name != "unified_fit":
        return no_model(session_id)

    if level is not None and not (1 <= level <= s.model.num_levels):
        return make_error(
            f"Level {level} does not exist (model has levels 1–{s.model.num_levels}).",
            code="BAD_LEVEL",
        )

    levels_to_check = (
        [level] if level is not None else list(range(1, s.model.num_levels + 1))
    )
    results = [
        {"level": lv, "feasible": bool(s.model.levels[lv - 1].check_physical_feasibility())}
        for lv in levels_to_check
    ]
    return {
        "feasibility":   results,
        "all_feasible":  all(r["feasible"] for r in results),
    }


# ---------------------------------------------------------------------------
# Category C'' — Q range
# ---------------------------------------------------------------------------

def get_data_q_range(session_id: str) -> dict:
    """Return the Q range of the loaded dataset."""
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    valid = np.isfinite(s.q)
    return {
        "q_min":    float(np.nanmin(s.q[valid])),
        "q_max":    float(np.nanmax(s.q[valid])),
        "n_points": int(len(s.q)),
    }


def get_fit_q_range(session_id: str) -> dict:
    """Return the Q range currently used for fitting."""
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    mask = fit_mask(s)
    n_fit = int(np.sum(mask))
    valid = np.isfinite(s.q)
    return {
        "q_min":         s.fit_q_min if s.fit_q_min is not None else float(np.nanmin(s.q[valid])),
        "q_max":         s.fit_q_max if s.fit_q_max is not None else float(np.nanmax(s.q[valid])),
        "n_fit_points":  n_fit,
        "is_full_range": s.fit_q_min is None and s.fit_q_max is None,
    }


def set_fit_q_range(
    session_id: str,
    q_min: Optional[float] = None,
    q_max: Optional[float] = None,
) -> dict:
    """Restrict the Q range used for fitting.

    Pass q_min and/or q_max to set limits.  Values are clamped to the data
    range.  Pass None to leave that end unchanged.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    valid = np.isfinite(s.q)
    data_qmin = float(np.nanmin(s.q[valid]))
    data_qmax = float(np.nanmax(s.q[valid]))

    if q_min is not None:
        s.fit_q_min = float(max(q_min, data_qmin))
    if q_max is not None:
        s.fit_q_max = float(min(q_max, data_qmax))

    # Sanity: ensure min < max
    lo = s.fit_q_min if s.fit_q_min is not None else data_qmin
    hi = s.fit_q_max if s.fit_q_max is not None else data_qmax
    if lo >= hi:
        s.fit_q_min = None
        s.fit_q_max = None
        return make_error(
            f"Resulting Q range [{lo}, {hi}] is empty — range reset to full data range.",
            code="EMPTY_RANGE",
        )

    return get_fit_q_range(session_id)


def reset_fit_q_range(session_id: str) -> dict:
    """Restore the fit Q range to the full data range."""
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    s.fit_q_min = None
    s.fit_q_max = None
    return get_fit_q_range(session_id)


# ---------------------------------------------------------------------------
# Category C'''' — Local estimators (Guinier / power law on a Q sub-range)
# ---------------------------------------------------------------------------
#
# Standalone Q-range fits of single-term models, equivalent to the GUI's
# "Fit Rg/G btwn cursors" and "Fit P/B btwn cursors" buttons.  Useful for
# generating good starting values before the full multi-level fit runs.
#
# These tools do NOT mutate the model — they only return the fitted
# numbers.  The agent decides whether to apply them via set_parameter_value().

def _slice_q_range(session: Session, q_min: float, q_max: float):
    """Return q, I, err arrays restricted to [q_min, q_max].  May be empty."""
    mask = (session.q >= q_min) & (session.q <= q_max)
    if session.error is not None:
        return session.q[mask], session.intensity[mask], session.error[mask]
    return session.q[mask], session.intensity[mask], None


def fit_local_guinier(
    session_id: str,
    q_min: float,
    q_max: float,
) -> dict:
    """Fit a single Guinier term I(q) = G · exp(-q²·Rg²/3) on a Q sub-range.

    Equivalent to the GUI's "Fit Rg/G btwn cursors" button.  Useful for
    estimating starting values for Rg and G on one level before running
    the full multi-level fit.

    Does **not** modify the session's model — call set_parameter_value()
    afterwards if you want to apply the result.

    Parameters
    ----------
    q_min, q_max : float
        Q range to fit in [Å⁻¹].  Choose a window covering the Guinier
        knee of the level you are characterising.

    Returns
    -------
    dict with keys:
      - G, Rg                : fitted values
      - q_min, q_max         : actual Q range used (may be clipped to data)
      - n_points             : number of data points used
      - chi_squared          : sum of squared (weighted) residuals
      - reduced_chi_squared  : χ² per degree of freedom
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    if q_max <= q_min:
        return make_error(
            f"q_max ({q_max}) must be greater than q_min ({q_min}).",
            code="BAD_RANGE",
        )

    q_fit, I_fit, err_fit = _slice_q_range(s, float(q_min), float(q_max))
    if len(q_fit) < 3:
        return make_error(
            f"Not enough data points in [{q_min}, {q_max}] ({len(q_fit)} points). "
            "Need at least 3 points; widen the Q range.",
            code="TOO_FEW_POINTS",
        )

    # Starting estimates (same heuristic as the GUI button)
    q_avg = (float(q_fit[0]) + float(q_fit[-1])) / 2
    g0    = float((I_fit[0] + I_fit[-1]) / 2)
    rg0   = 2 * np.pi / q_avg if q_avg > 0 else 10.0

    from scipy.optimize import curve_fit  # noqa: PLC0415

    def guinier(q, G, Rg):
        return G * np.exp(-(q ** 2) * (Rg ** 2) / 3.0)

    sigma = err_fit if err_fit is not None and np.all(err_fit > 0) else None
    try:
        popt, _ = curve_fit(
            guinier, q_fit, I_fit,
            p0=[g0, rg0],
            sigma=sigma,
            absolute_sigma=False,
            maxfev=5000,
        )
    except Exception as exc:
        return make_error(
            f"Local Guinier fit failed: {exc}",
            suggestion="Try a different Q range that covers a visible Guinier knee.",
            code="FIT_FAILED",
        )

    g_fit  = float(abs(popt[0]))
    rg_fit = float(abs(popt[1]))

    model_I = guinier(q_fit, g_fit, rg_fit)
    resid   = (I_fit - model_I) / (err_fit if sigma is not None else 1.0)
    chi_sq  = float(np.sum(resid ** 2))
    dof     = max(len(q_fit) - 2, 1)

    return {
        "ok":                  True,
        "G":                   g_fit,
        "Rg":                  rg_fit,
        "q_min":               float(q_fit[0]),
        "q_max":               float(q_fit[-1]),
        "n_points":            int(len(q_fit)),
        "chi_squared":         chi_sq,
        "reduced_chi_squared": chi_sq / dof,
    }


def fit_local_power_law(
    session_id: str,
    q_min: float,
    q_max: float,
) -> dict:
    """Fit a single power-law term I(q) = B · q⁻ᴾ on a Q sub-range.

    Equivalent to the GUI's "Fit P/B btwn cursors" button.  Useful for
    estimating starting values for P and B on one level (typically the
    high-Q power-law region after the Guinier knee).

    Does **not** modify the session's model — call set_parameter_value()
    afterwards if you want to apply the result.

    Parameters
    ----------
    q_min, q_max : float
        Q range to fit in [Å⁻¹].  Pick a window over the linear power-law
        portion of the log-log plot, beyond the Guinier knee.

    Returns
    -------
    dict with keys:
      - P, B                 : fitted values
      - q_min, q_max         : actual Q range used (may be clipped to data)
      - n_points             : number of data points used
      - chi_squared          : sum of squared (weighted) residuals
      - reduced_chi_squared  : χ² per degree of freedom
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    if q_max <= q_min:
        return make_error(
            f"q_max ({q_max}) must be greater than q_min ({q_min}).",
            code="BAD_RANGE",
        )

    q_fit, I_fit, err_fit = _slice_q_range(s, float(q_min), float(q_max))
    if len(q_fit) < 3:
        return make_error(
            f"Not enough data points in [{q_min}, {q_max}] ({len(q_fit)} points). "
            "Need at least 3 points; widen the Q range.",
            code="TOO_FEW_POINTS",
        )

    # Reject non-positive intensities (log-slope estimate would crash)
    pos = (I_fit > 0) & (q_fit > 0)
    if int(np.sum(pos)) < 3:
        return make_error(
            "Not enough positive intensity values in range for a power-law fit.",
            code="TOO_FEW_POSITIVE",
        )
    q_pos = q_fit[pos]
    I_pos = I_fit[pos]

    # Starting estimates from endpoints (same as GUI button)
    log_q_first = float(np.log(q_pos[0]))
    log_q_last  = float(np.log(q_pos[-1]))
    log_I_first = float(np.log(I_pos[0]))
    log_I_last  = float(np.log(I_pos[-1]))
    denom = log_q_last - log_q_first
    p0 = abs((log_I_first - log_I_last) / denom) if denom != 0 else 4.0
    b0 = float(I_pos[0] * (q_pos[0] ** p0))

    from scipy.optimize import curve_fit  # noqa: PLC0415

    def power_law(q, B, P):
        return B * np.power(q, -P)

    err_pos = err_fit[pos] if err_fit is not None else None
    sigma   = err_pos if err_pos is not None and np.all(err_pos > 0) else None
    try:
        popt, _ = curve_fit(
            power_law, q_pos, I_pos,
            p0=[b0, p0],
            sigma=sigma,
            absolute_sigma=False,
            maxfev=5000,
        )
    except Exception as exc:
        return make_error(
            f"Local power-law fit failed: {exc}",
            suggestion="Try a different Q range over the linear power-law region.",
            code="FIT_FAILED",
        )

    b_fit = float(abs(popt[0]))
    p_fit = float(abs(popt[1]))

    model_I = power_law(q_pos, b_fit, p_fit)
    resid   = (I_pos - model_I) / (err_pos if sigma is not None else 1.0)
    chi_sq  = float(np.sum(resid ** 2))
    dof     = max(len(q_pos) - 2, 1)

    return {
        "ok":                  True,
        "P":                   p_fit,
        "B":                   b_fit,
        "q_min":               float(q_pos[0]),
        "q_max":               float(q_pos[-1]),
        "n_points":            int(len(q_pos)),
        "chi_squared":         chi_sq,
        "reduced_chi_squared": chi_sq / dof,
    }


# ---------------------------------------------------------------------------
# Category C''''' — Feature detection (slope-profile analysis)
# ---------------------------------------------------------------------------
#
# Analyses the loaded data in log-log space and identifies plateaus
# (Guinier knees), peaks (Guinier + structure factor), and power-law
# regions.  Intended to help the agent decide how many Unified Fit levels
# a curve needs and where to place Q-windows for local Guinier/Porod fits.
# Does not mutate the model.

def detect_features(
    session_id: str,
    q_min: Optional[float] = None,
    q_max: Optional[float] = None,
    q_max_clip: Optional[float] = 0.6,
    config_overrides: Optional[dict] = None,
) -> dict:
    """Segment the loaded I(Q) curve into power-law slope segments.

    Returns a piecewise classification of the curve in log-log space:
    each segment has a locally-constant slope; adjacent segments with
    different slopes imply Guinier knees between them (also returned).
    Intended to help the user / agent decide how many Unified Fit levels
    are needed and where to place Q-windows for ``fit_local_guinier`` /
    ``fit_local_power_law``.

    Does **not** modify the session's model or fit state.

    Parameters
    ----------
    session_id : str
        Session whose data should be analysed.
    q_min, q_max : float, optional
        Restrict analysis to this Q sub-range.  Defaults to full data range.
    q_max_clip : float, optional
        Drop data with Q greater than this value before segmentation.
        Default 0.6 Å⁻¹ — the practical upper limit of the small-angle
        approximation; above this, amorphous diffraction features should
        not be treated as SAS structure.  Pass ``None`` to disable.
    config_overrides : dict, optional
        Mapping of ``FeatureDetectConfig`` field names to override values
        (e.g. ``{"change_threshold_1": 0.60, "min_segment_decades": 0.15}``).

    Returns
    -------
    dict with keys:
      - ok                            : True on success
      - segments                      : list of {q_min, q_max, P, P_std, kind,
                                                 intensity_mid, width_decades}
        Ordered HIGH-Q → LOW-Q to match Unified Fit level numbering
        (segments[0] = Level 1, highest Q, smallest structure).
        kind is "background" | "guinier_plateau" | "power_law"
      - guinier_knees                 : list of {q_min, q_max, q_center,
                                                 P_low_q, P_high_q, delta_P}
        P_low_q and P_high_q are positive Porod exponents (I ∝ Q^-P).
        Only transitions where P_low_q < P_high_q are listed — exponent
        is smaller (shallower) at low Q, the physical Guinier-knee condition.
      - segments                      : each segment also carries "P" and "P_std"
        (positive Porod exponent, I ∝ Q^-P) in addition to "slope" / "slope_std".
      - background_q_min              : float | null
      - recommended_guinier_windows   : list of {feature_type, q_min_guinier,
                                                 q_max_guinier, q_min_powerlaw}
      - recommended_nlevels           : suggested number of Unified Fit levels
      - n_segments_found              : int
      - log_decades                   : total log10(Q) span analysed
      - q_min_analysed, q_max_analysed : actual Q range analysed
      - n_points                      : number of data points used
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    from pyirena.core.feature_detect import (  # noqa: PLC0415
        FeatureDetectConfig,
        detect_features as _detect,
    )

    if q_min is None:
        q_min = float(s.q.min())
    if q_max is None:
        q_max = float(s.q.max())
    if q_max <= q_min:
        return make_error(
            f"q_max ({q_max}) must be greater than q_min ({q_min}).",
            code="BAD_RANGE",
        )

    q_use, I_use, err_use = _slice_q_range(s, float(q_min), float(q_max))
    if q_use.size < 8:
        return make_error(
            f"Not enough data points in [{q_min}, {q_max}] ({q_use.size} points). "
            "Need at least 8 points; widen the Q range.",
            code="TOO_FEW_POINTS",
        )

    cfg = FeatureDetectConfig(q_max_clip=q_max_clip)
    if config_overrides:
        for k, v in config_overrides.items():
            if not hasattr(cfg, k):
                return make_error(
                    f"Unknown FeatureDetectConfig field: {k!r}",
                    suggestion="Valid fields include change_window_1, "
                               "change_threshold_1, change_window_2, "
                               "change_threshold_2, wide_region_decades, "
                               "min_segment_decades, edge_min_segment_decades, "
                               "merge_slope_tol, guinier_knee_min_delta_slope.",
                    code="BAD_OVERRIDE",
                )
            try:
                if v is None:
                    setattr(cfg, k, None)
                else:
                    setattr(cfg, k, float(v))
            except (TypeError, ValueError):
                return make_error(
                    f"Invalid value for config field {k!r}: {v!r}",
                    code="BAD_OVERRIDE",
                )

    result = _detect(q_use, I_use, sigma_I=err_use, config=cfg)
    out = result.to_dict()
    out["ok"] = True
    return out


# ---------------------------------------------------------------------------
# Category D — Fit execution
# ---------------------------------------------------------------------------

def run_fit(
    session_id: str,
    max_iter: Optional[int] = None,
    tolerance: Optional[float] = None,
    random_seed: Optional[int] = None,
) -> dict:
    """Run the fitting algorithm on the current session's model and data.

    The Q-range restriction (set_fit_q_range) is applied before fitting.
    Re-running fit() uses the current parameter values as the starting point
    (they are updated in place by the previous fit).

    Parameters
    ----------
    random_seed : int, optional
        Seed for numpy's random number generator before the fit.  The
        underlying scipy optimiser is deterministic, so this only matters for
        fitting methods that have stochastic components (e.g. future MCSaS
        integration).  Pass the same seed to reproduce a fit exactly from an
        audit trail.  Stored in the returned dict and in last_fit_result so
        pyirena-ai can stamp it into the audit JSON.

    Returns
    -------
    dict with keys: success, chi_squared, reduced_chi_squared, iterations,
    message, random_seed, parameters_updated (list of {name, value}).
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)

    mask = fit_mask(s)
    if not np.any(mask):
        return make_error(
            "No data points in the current fit Q range.",
            suggestion="Call reset_fit_q_range() or set_fit_q_range() with wider limits.",
            code="EMPTY_RANGE",
        )

    q_fit = s.q[mask]
    I_fit = s.intensity[mask]
    err_fit = s.error[mask] if s.error is not None else None

    if random_seed is not None:
        np.random.seed(int(random_seed))

    kwargs: dict = {}
    if max_iter is not None:
        kwargs["max_iterations"] = int(max_iter)

    try:
        result = s.model.fit(q_fit, I_fit, err_fit, **kwargs)
    except Exception as exc:
        return make_error(
            f"Fit failed with exception: {exc}",
            suggestion="Check that parameter bounds are reasonable and data is valid.",
            code="FIT_EXCEPTION",
        )

    result["random_seed"] = random_seed  # None if not supplied; stored for audit trail
    s.last_fit_result = result

    # Build per-parameter update list (current fitted values)
    updated = [
        {"name": row["name"], "value": row["value"]}
        for row in _param_table(s.model)
    ]

    # Robust fit-quality scalars (additive; complement reduced χ² when σ are
    # mis-scaled).  Full arrays/bands available via get_fit_quality().
    quality = _quality_scalars(_compute_quality(s.model))

    return {
        "success":             bool(result.get("success", False)),
        "chi_squared":         float(result.get("chi_squared", float("nan"))),
        "reduced_chi_squared": float(result.get("reduced_chi_squared", float("nan"))),
        "iterations":          int(result.get("n_iterations", 0)),
        "message":             str(result.get("message", "")),
        "random_seed":         random_seed,
        "parameters_updated":  updated,
        "quality":             quality,
    }


# ---------------------------------------------------------------------------
# Category E — Quality assessment
# ---------------------------------------------------------------------------

def get_chi_squared(session_id: str) -> dict:
    """Return χ² and reduced χ² from the last fit."""
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.last_fit_result is None:
        return no_fit(session_id)
    return {
        "chi_squared":         float(s.last_fit_result.get("chi_squared", float("nan"))),
        "reduced_chi_squared": float(s.last_fit_result.get("reduced_chi_squared", float("nan"))),
    }


def get_residuals(session_id: str) -> dict:
    """Return residuals array and summary statistics from the last fit."""
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.last_fit_result is None:
        return no_fit(session_id)

    model = s.model
    if model.fit_intensity is None:
        return no_fit(session_id)

    resid = model.I_data - model.fit_intensity
    if model.error_data is not None and np.any(model.error_data > 0):
        norm_resid = resid / model.error_data
    else:
        norm_resid = resid

    # Robust scale and alternative residual views (additive — the legacy
    # "residuals" key is unchanged).  See get_fit_quality() for the full set.
    metrics = _compute_quality(model)
    s = metrics["robust_scale_s"] if metrics else None
    has_s = s is not None and np.isfinite(s) and s > 0

    # Rescaled residual r' = r / s : compares scatter to the data's own robust
    # noise floor rather than to the (possibly mis-scaled) reported σ.
    rescaled = norm_resid / s if has_s else None
    # Fractional misfit (I−M)/I in %, σ-independent; guard |I| ~ 0.
    with np.errstate(divide="ignore", invalid="ignore"):
        frac_pct = np.where(np.abs(model.I_data) > 0, resid / model.I_data * 100.0, np.nan)

    max_pts = int(os.environ.get("PYIRENA_MAX_ARRAY_POINTS", 500))
    return {
        "q":          _decimated(model.q_data, max_pts),
        "residuals":  _decimated(norm_resid, max_pts),
        "rescaled_residual":   None if rescaled is None else _decimated(rescaled, max_pts),
        "frac_misfit_percent": _decimated(frac_pct, max_pts),
        "summary": {
            "rms":     float(np.sqrt(np.mean(norm_resid ** 2))),
            "max_abs": float(np.max(np.abs(norm_resid))),
            "mean":    float(np.mean(norm_resid)),
            "robust_scale_s": float(s) if has_s else None,
        },
        "normalised": model.error_data is not None,
    }


def _n_free_params(model) -> int:
    """Number of free (non-fixed) fit parameters — used for degrees of freedom."""
    return sum(1 for row in _param_table(model) if not row.get("fixed", False))


def _compute_quality(model, n_bands: int = 4) -> Optional[dict]:
    """Run fit_quality_metrics on a session's last fit. None if no usable fit."""
    if model is None or model.fit_intensity is None:
        return None
    return fit_quality_metrics(
        q=model.q_data,
        intensity=model.I_data,
        model=model.fit_intensity,
        sigma=model.error_data,
        n_params=_n_free_params(model),
        n_bands=n_bands,
    )


def _quality_scalars(metrics: Optional[dict]) -> dict:
    """The decision-relevant scalar subset, JSON-safe, for embedding in run_fit.

    Returns an empty dict if metrics could not be computed. None values are kept
    (they carry meaning: e.g. robust_scale_s is None when sigma is unavailable).
    """
    if not metrics:
        return {}

    def _f(x):
        return float(x) if x is not None and np.isfinite(x) else None

    return {
        "sigma_available":              bool(metrics["sigma_available"]),
        "robust_scale_s":               _f(metrics["robust_scale_s"]),
        "realistic_reduced_chi2_floor": _f(metrics["realistic_reduced_chi2_floor"]),
        "max_abs_frac_misfit":          _f(metrics["max_abs_frac_misfit"]),
        "q_at_max_frac_misfit":         _f(metrics["q_at_max_frac_misfit"]),
        "median_frac_uncertainty":      _f(metrics["median_frac_uncertainty"]),
        "n_outliers_3s":                metrics["n_outliers_3s"],
        "longest_same_sign_run":        int(metrics["longest_same_sign_run"]),
        "sign_autocorr_lag1":           _f(metrics["sign_autocorr_lag1"]),
    }


def _serialize_quality(metrics: Optional[dict], max_pts: int) -> Optional[dict]:
    """Convert a fit_quality_metrics dict to a fully JSON-serialisable form.

    Arrays are decimated to max_pts (NaN/inf -> None). Scalar None values are
    preserved. Returns None if metrics is None.
    """
    if metrics is None:
        return None

    def _f(x):
        return float(x) if x is not None and np.isfinite(x) else None

    def _arr(a):
        return None if a is None else _decimated(np.asarray(a, dtype=float), max_pts)

    bands = [
        {
            "q_lo":                float(b["q_lo"]),
            "q_hi":                float(b["q_hi"]),
            "n":                   int(b["n"]),
            "reduced_chi2":        _f(b["reduced_chi2"]),
            "robust_scale_s":      _f(b["robust_scale_s"]),
            "max_abs_frac_misfit": _f(b["max_abs_frac_misfit"]),
        }
        for b in metrics["bands"]
    ]

    return {
        "sigma_available":              bool(metrics["sigma_available"]),
        "n_valid":                      int(metrics["n_valid"]),
        "n_params":                     int(metrics["n_params"]),
        "dof":                          metrics["dof"],
        # per-point arrays (aligned to q_valid)
        "q":                            _arr(metrics["q_valid"]),
        "norm_residual":                _arr(metrics["norm_residual"]),
        "frac_residual":                _arr(metrics["frac_residual"]),
        "frac_uncertainty":             _arr(metrics["frac_uncertainty"]),
        # global scalars
        "reduced_chi2":                 _f(metrics["reduced_chi2"]),
        "robust_scale_s":               _f(metrics["robust_scale_s"]),
        "sigma_misscale_factor":        _f(metrics["sigma_misscale_factor"]),
        "realistic_reduced_chi2_floor": _f(metrics["realistic_reduced_chi2_floor"]),
        "median_frac_uncertainty":      _f(metrics["median_frac_uncertainty"]),
        "max_abs_frac_misfit":          _f(metrics["max_abs_frac_misfit"]),
        "q_at_max_frac_misfit":         _f(metrics["q_at_max_frac_misfit"]),
        "n_outliers_3s":                metrics["n_outliers_3s"],
        "frac_outliers_3s":             _f(metrics["frac_outliers_3s"]),
        # structure scalars
        "longest_same_sign_run":        int(metrics["longest_same_sign_run"]),
        "sign_autocorr_lag1":           _f(metrics["sign_autocorr_lag1"]),
        # per-band
        "n_bands_used":                 int(metrics["n_bands_used"]),
        "bands":                        bands,
    }


def get_fit_quality(session_id: str, n_bands: int = 4) -> dict:
    """Return robust, sigma-scale-independent fit-quality diagnostics.

    Complements get_chi_squared / get_residuals with metrics that stay
    informative when the reported uncertainties are mis-scaled — a common SAXS
    situation where chasing reduced χ² ≈ 1 is misleading.

    Key fields to read:
      - ``robust_scale_s``: MAD-based estimate of how many times the *actual*
        scatter exceeds the reported σ. s ≈ 1 → σ honest; s ≈ 3 → σ ~3× too
        small, so the realistic reduced-χ² floor is ~9 (``realistic_reduced_chi2_floor``).
      - ``max_abs_frac_misfit`` (+ ``q_at_max_frac_misfit``): the largest
        |(I−M)/I|, a σ-independent backstop. ≳ 0.3 means a gross local misfit no
        matter how unreliable σ is.
      - ``n_outliers_3s``: points beyond 3·robust_scale_s — genuine outliers
        even after accounting for a mis-scaled σ.
      - ``longest_same_sign_run`` / ``sign_autocorr_lag1``: structure that
        signals a wrong functional form (systematic), distinct from a pure
        σ-scale problem.
      - ``bands``: the same metrics per Q-decade; an uneven per-band χ² is itself
        a misfit signal.

    This function returns *facts only* — no good/bad verdict. Interpretation/
    thresholds belong to the caller.

    Returns
    -------
    dict
        The serialised metrics (see fields above), or an error dict if there is
        no session or no fit yet.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.last_fit_result is None or s.model is None or s.model.fit_intensity is None:
        return no_fit(session_id)

    max_pts = int(os.environ.get("PYIRENA_MAX_ARRAY_POINTS", 500))
    metrics = _compute_quality(s.model, n_bands=int(n_bands))
    return _serialize_quality(metrics, max_pts)


def get_fit_image(
    session_id: str,
    width: int = 1024,
    height: int = 768,
    format: str = "png",
) -> dict:
    """Capture the current fit as a PNG image (base64-encoded).

    Works before a fit is run (shows data + initial model) and after a fit
    (adds residuals subplot).  The layout mirrors pyirena's Unified Fit
    results window.

    Parameters
    ----------
    width, height : int
        Output image dimensions in pixels.
    format : str
        Only "png" is supported in Phase 1.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)

    try:
        b64 = _render_fit_image(s, width=width, height=height)
    except Exception as exc:
        return make_error(f"Image generation failed: {exc}", code="PLOT_ERROR")

    return {
        "image_base64": b64,
        "format":       "png",
        "width":        width,
        "height":       height,
        "has_residuals": s.last_fit_result is not None,
    }


def get_residuals_image(
    session_id: str,
    width: int = 1024,
    height: int = 768,
) -> dict:
    """Capture the fit+residuals panel as a PNG.

    Returns an error if no fit has been run yet (there are no residuals to show).
    The image is identical to get_fit_image() but requires a completed fit.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.last_fit_result is None:
        return no_fit(session_id)
    return get_fit_image(session_id, width=width, height=height)


def get_parameter_uncertainties(session_id: str) -> dict:
    """Return parameter uncertainties from the last fit.

    Note: scipy least_squares does not return a covariance matrix by default.
    Uncertainty estimates are not available in Phase 1.  This function returns
    a placeholder indicating they are unavailable.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.last_fit_result is None:
        return no_fit(session_id)

    return {
        "available": False,
        "note": (
            "Parameter uncertainties are not computed in Phase 1. "
            "A future phase will add covariance-based uncertainty estimation."
        ),
    }


# ---------------------------------------------------------------------------
# Category F — Persistence
# ---------------------------------------------------------------------------

def _session_to_gui_state(session: Session) -> dict:
    """Snapshot the session's UnifiedFitModel as a GUI-compatible state dict.

    The shape matches ``UnifiedFitPanel.get_current_state()`` exactly so the
    panel's ``apply_state`` can restore every control verbatim — the whole
    point of embedding this in the NXcanSAS file.
    """
    model = session.model
    state: dict = {
        "num_levels":     int(model.num_levels),
        "background": {
            "value": float(model.background),
            "fit":   bool(model.fit_background),
        },
        "cursor_left":    session.fit_q_min,
        "cursor_right":   session.fit_q_max,
        "update_auto":    False,
        "display_local":  False,
        "no_limits":      False,
        "skip_fit_check": False,
        "store_local":    False,
        "levels":         [],
    }

    def _lim(level, attr):
        lo, hi = getattr(level, f"{attr}_limits", (None, None))
        return (None if lo is None else float(lo),
                None if hi is None else float(hi))

    for i, lv in enumerate(model.levels):
        g_lo, g_hi   = _lim(lv, "G")
        rg_lo, rg_hi = _lim(lv, "Rg")
        b_lo, b_hi   = _lim(lv, "B")
        p_lo, p_hi   = _lim(lv, "P")
        e_lo, e_hi   = _lim(lv, "ETA")
        k_lo, k_hi   = _lim(lv, "PACK")

        state["levels"].append({
            "level":    i + 1,
            "G":   {"value": float(lv.G),   "fit": bool(lv.fit_G),   "low_limit": g_lo,  "high_limit": g_hi},
            "Rg":  {"value": float(lv.Rg),  "fit": bool(lv.fit_Rg),  "low_limit": rg_lo, "high_limit": rg_hi},
            "B":   {"value": float(lv.B),   "fit": bool(lv.fit_B),   "low_limit": b_lo,  "high_limit": b_hi},
            "P":   {"value": float(lv.P),   "fit": bool(lv.fit_P),   "low_limit": p_lo,  "high_limit": p_hi},
            "ETA": {"value": float(lv.ETA), "fit": bool(lv.fit_ETA), "low_limit": e_lo,  "high_limit": e_hi},
            "PACK":{"value": float(lv.PACK),"fit": bool(lv.fit_PACK),"low_limit": k_lo,  "high_limit": k_hi},
            # UnifiedLevel uses RgCO; the GUI key is RgCutoff
            "RgCutoff":   float(lv.RgCO),
            # GUI naming differs from model naming for the boolean flags
            "correlated": bool(lv.correlations),
            "estimate_B": bool(lv.link_B),
            "link_rgco":  bool(lv.link_RGCO),
        })

    return state


def save_fit(session_id: str, output_path: Optional[str] = None) -> dict:
    """Save the fitted result back to NXcanSAS HDF5.

    Parameters
    ----------
    output_path : str, optional
        Where to save.  Defaults to overwriting the original file.
        Saving to a different path preserves the original.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)
    if s.last_fit_result is None:
        return no_fit(session_id)

    from pyirena.io.nxcansas_unified import save_unified_fit_results

    model = s.model
    target = Path(output_path) if output_path else Path(s.file_path)

    # Convert UnifiedLevel objects to the dict format expected by nxcansas_unified
    levels_dicts = [
        {
            "Rg": lv.Rg, "G": lv.G, "P": lv.P, "B": lv.B,
            "ETA": lv.ETA, "PACK": lv.PACK, "RgCO": lv.RgCO,
        }
        for lv in model.levels
    ]
    resid = model.I_data - model.fit_intensity
    norm_resid = (resid / model.error_data
                  if model.error_data is not None else resid)

    try:
        save_unified_fit_results(
            filepath=target,
            q=model.q_data,
            intensity_data=model.I_data,
            intensity_model=model.fit_intensity,
            residuals=norm_resid,
            levels=levels_dicts,
            background=model.background,
            chi_squared=float(s.last_fit_result.get("chi_squared", float("nan"))),
            num_levels=model.num_levels,
            error=model.error_data,
            # Embed the full setup so the GUI can "Load Setup from File…"
            # and resume exactly where pyirena-ai left off.
            setup_state=_session_to_gui_state(s),
        )
    except Exception as exc:
        return make_error(
            f"Could not save fit to '{target}': {exc}",
            suggestion="Check write permissions and that the path is valid.",
            code="SAVE_ERROR",
        )

    return {"ok": True, "file_path": str(target)}


def export_fit_report(session_id: str, format: str = "json") -> dict:
    """Export a human-readable or machine-readable fit report.

    Parameters
    ----------
    format : "json" | "markdown"
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)
    if s.last_fit_result is None:
        return no_fit(session_id)

    params = _param_table(s.model)
    chi_sq = s.last_fit_result.get("reduced_chi_squared")
    chi_sq_str = f"{chi_sq:.4f}" if chi_sq is not None else "n/a"
    quality = _quality_scalars(_compute_quality(s.model))

    if format == "json":
        import json
        report = {
            "file":            s.file_path,
            "model":           s.model_name,
            "nlevels":         s.model.num_levels,
            "fit_q_min":       s.fit_q_min,
            "fit_q_max":       s.fit_q_max,
            "reduced_chi_squared": chi_sq,
            "quality":         quality,
            "parameters":      params,
        }
        return {"format": "json", "content": json.dumps(report, indent=2)}

    elif format == "markdown":
        lines = [
            f"# Unified Fit Report",
            f"",
            f"**File:** {s.file_path}",
            f"**Model:** Unified Fit — {s.model.num_levels} level(s)",
            f"**Reduced χ²:** {chi_sq_str}",
        ]
        if quality:
            s_val = quality.get("robust_scale_s")
            floor = quality.get("realistic_reduced_chi2_floor")
            mfm = quality.get("max_abs_frac_misfit")
            if s_val is not None:
                floor_s = f" (realistic reduced-χ² floor ≈ {floor:.1f})" if floor is not None else ""
                lines.append(f"**σ-scale (robust):** {s_val:.2f}×{floor_s}")
            if mfm is not None:
                lines.append(f"**Max |(I−M)/I|:** {mfm * 100:.1f}%")
            csr = quality.get("longest_same_sign_run")
            if csr is not None:
                lines.append(f"**Longest same-sign run:** {csr}")
        lines += [
            f"",
            f"## Parameters",
            f"",
            f"| Parameter | Value | Fixed | Lo | Hi | Units |",
            f"|-----------|-------|-------|----|----|-------|",
        ]
        for p in params:
            lo_s = f"{p['lo']:.3g}" if p["lo"] is not None else "—"
            hi_s = f"{p['hi']:.3g}" if p["hi"] is not None else "—"
            lines.append(
                f"| {p['name']} | {p['value']:.4g} | {'Yes' if p['fixed'] else 'No'} "
                f"| {lo_s} | {hi_s} | {p['units']} |"
            )
        return {"format": "markdown", "content": "\n".join(lines)}

    else:
        return make_error(
            f"Unknown format '{format}'.",
            suggestion="Use 'json' or 'markdown'.",
            code="BAD_FORMAT",
        )
