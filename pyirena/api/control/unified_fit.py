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

import os
from typing import Optional

import numpy as np

from pyirena.api._json import decimate_for_json as _decimated  # shared helper (#38)
from pyirena.api.control._images import render_png
from pyirena.api.control._save import resolve_save_target

# Session lifecycle, fit Q range, local estimators and feature detection serve
# every tool and live in data_session.py; re-exported here, where they used to
# be defined, so existing imports work.
from pyirena.api.control.data_session import (  # noqa: F401
    close_session,
    detect_features,
    fit_local_guinier,
    fit_local_power_law,
    get_data_q_range,
    get_fit_q_range,
    get_session_summary,
    list_open_sessions,
    open_dataset,
    open_dataset_from_data,
    reset_fit_q_range,
    set_fit_q_range,
)
from pyirena.api.control.errors import (
    bad_param,
    make_error,
    no_fit,
    no_model,
    no_session,
)
from pyirena.api.control.session import (
    Session,
    fit_mask,
    get_session,
)

# get_parameter_uncertainties serves every tool and lives in uncertainty.py;
# re-exported here, where it used to be defined, so existing imports work.
from pyirena.api.control.uncertainty import get_parameter_uncertainties  # noqa: F401
from pyirena.core.fit_metrics import fit_quality_metrics
from pyirena.core.unified import UNIFIED_STATE_SCHEMA_VERSION

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


def _render_fit_image(
    session: Session, width: int, height: int, dpi: int = 120,
) -> tuple[str, str]:
    """Render data + model (+ residuals if fit was run) as a PNG.

    Returns ``(base64 PNG, absolute path on disk)``.

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

    # Save to the plot cache and return both the base64 payload and the path
    return render_png(fig, f"fit_{session.session_id}", dpi)


# ---------------------------------------------------------------------------
# Category B — Model selection and inspection
# ---------------------------------------------------------------------------

AVAILABLE_MODELS = ["unified_fit", "sizes"]


def list_available_models() -> dict:
    """List the fitting models available in this version of the control API."""
    return {
        "models": AVAILABLE_MODELS,
        "details": {
            "unified_fit": (
                "Beaucage multi-level Unified Fit. Use select_model() and the "
                "parameter/level tools (set_parameter_value, fix_all_except, "
                "add_unified_level, …)."
            ),
            "sizes": (
                "Particle size-distribution inversion (MaxEnt / Regularization / "
                "TNNLS / Monte Carlo). Use select_sizes_model() and the sizes "
                "tools (set_size_grid, set_shape, set_method, set_error_handling, "
                "fit_power_law_background, fit_flat_background, run_sizes_fit, …). "
                "Best for dilute samples with a single identifiable particle "
                "population; call suggest_sizes_setup() first."
            ),
        },
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
        # Auto-enable slit smearing when the session's data are slit smeared
        # (SasView-style: slit-length presence drives model smearing).
        if s.is_slit_smeared and s.slit_length > 0:
            s.model.use_slit_smearing = True
            s.model.slit_length = s.slit_length

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
# Category D — Fit execution
# ---------------------------------------------------------------------------

def run_fit(
    session_id: str,
    max_iter: Optional[int] = None,
    tolerance: Optional[float] = None,
    random_seed: Optional[int] = None,
    walk_limits: bool = True,
    no_limits: bool = False,
    n_mc_runs: int = 0,
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
    walk_limits : bool, optional
        Default True: use
        :meth:`~pyirena.core.unified.UnifiedFitModel.fit_with_limit_walking`,
        matching the GUI Fit button — when a fitted parameter ends pinned at
        a limit, the limits are recentred on the fitted value
        (``panel_auto_limits``) and the fit rerun, so an optimum far outside
        the current limits window is still reached in one call.  Pass False
        to treat the current bounds (e.g. from set_parameter_bounds) as hard
        constraints; the fit then stops at the bound and reports the
        parameter in ``pinned_parameters``.
    no_limits : bool, optional
        Fit without the bounds set on the model — the GUI's "No limits?"
        check box. For Unified Fit that means the model's own wide defaults
        (Rg 0.1–1e6 Å, G 1e-10–1e10, …) rather than the limits you or the
        config supplied; the bounds themselves are untouched and come back
        after the fit. Implies ``walk_limits=False``, since there is nothing
        left to walk. Use it when a fit is stuck against a limit and you
        want to see where it actually wants to go.
    n_mc_runs : int, optional
        Monte-Carlo uncertainty refits after the fit; 0 (default) = none.
        The same estimate as the GUI's "Analyze Results" and
        ``batch.fit_unified(with_uncertainty=True)``: the fitted model —
        links, limits, slit smearing — is refitted to data perturbed by the
        measured σ. Needs measured uncertainties. *random_seed* makes it
        reproducible. Slow: n refits; mind a remote client's timeout.

    Returns
    -------
    dict with keys: success, chi_squared, reduced_chi_squared, iterations,
    message, random_seed, parameters_updated (list of {name, value}, plus
    ``std`` when Monte-Carlo ran), pinned_parameters, warnings, and
    ``uncertainty`` when n_mc_runs > 0.
    """
    from pyirena.api.control import uncertainty as unc

    n_mc_runs, bad = unc.check_n_mc_runs(n_mc_runs)
    if bad:
        return bad
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

    if random_seed is not None:
        np.random.seed(int(random_seed))

    # The shared Unified runner (core/fit_runner.py) — the same code batch
    # uses: data rule, walking (no walking with no_limits, since the wide
    # defaults leave nothing to walk), Monte-Carlo.
    from pyirena.core.fit_runner import run_unified
    try:
        run = run_unified(
            s.model, s.q, s.intensity, s.error,
            q_min=s.fit_q_min, q_max=s.fit_q_max,
            error_fraction=getattr(s, "error_fraction", 0.05),
            no_limits=bool(no_limits), walk_limits=bool(walk_limits),
            max_iterations=max_iter, n_mc_runs=n_mc_runs, seed=random_seed)
        result = run.result
    except Exception as exc:
        return make_error(
            f"Fit failed with exception: {exc}",
            suggestion="Check that parameter bounds are reasonable and data is valid.",
            code="FIT_EXCEPTION",
        )

    result["random_seed"] = random_seed  # None if not supplied; stored for audit trail
    s.last_fit_result = result
    s.last_uncertainty = None

    mc_block = None
    if run.mc is not None:
        mc = run.mc
        rows = []
        if mc.get("available"):
            u = mc["uncertainties"]
            for i, lv in enumerate(u["levels"]):
                for key, std in lv.items():
                    rows.append({"name": f"{key}_{i + 1}", "std": unc._f(std)})
            rows.append({"name": "background", "std": unc._f(u["background"])})
        mc_block = unc.summarize_mc(mc, rows)
        mc_block["_levels"] = mc.get("uncertainties")
        s.last_uncertainty = mc_block

    # Build per-parameter update list (current fitted values)
    std_by_name = ({row["name"]: row["std"] for row in mc_block["parameters"]}
                   if mc_block else {})
    updated = []
    for row in _param_table(s.model):
        entry = {"name": row["name"], "value": row["value"]}
        if row["name"] in std_by_name:
            entry["std"] = std_by_name[row["name"]]
        updated.append(entry)

    # Robust fit-quality scalars (additive; complement reduced χ² when σ are
    # mis-scaled).  Full arrays/bands available via get_fit_quality().
    quality = _quality_scalars(_compute_quality(s.model))

    # Parameters still sitting at a fit limit after the fit (after walking,
    # this list is normally empty; with walk_limits=False it shows which
    # bound stopped the fit — the difficult-to-debug case).
    pinned = [f"level {i + 1} {name}"
              for i, name in s.model.pinned_fitted_parameters()]
    warnings_list = list(run.notes)
    if pinned and no_limits:
        # The fit did not use these bounds, so nothing was pinned *at* one —
        # but the values landing outside them is the answer to the question
        # "where does this parameter actually want to go", which is why
        # no_limits was asked for. Say that, rather than "pinned at a limit".
        warnings_list.append(
            "with no_limits the fit went outside the bounds still set on the "
            "model, for: " + ", ".join(pinned)
            + " — widen or clear those bounds before fitting normally again, "
            "or the next bounded fit will stop at them")
        pinned = []
    elif pinned:
        warnings_list.append(
            "fitted parameter(s) pinned at a limit: " + ", ".join(pinned)
            + " — the reported value is dictated by the bound, not the data"
            + ("" if walk_limits else
               "; rerun with walk_limits=True or widen the bounds"))

    return {
        "no_limits":           bool(no_limits),
        "success":             bool(result.get("success", False)),
        "chi_squared":         float(result.get("chi_squared", float("nan"))),
        "reduced_chi_squared": float(result.get("reduced_chi_squared", float("nan"))),
        "iterations":          int(result.get("n_iterations", 0)),
        "message":             str(result.get("message", "")),
        "random_seed":         random_seed,
        "parameters_updated":  updated,
        "quality":             quality,
        "pinned_parameters":   pinned,
        "warnings":            warnings_list,
        # Report smearing state so the AI advisor can explain that fitted
        # parameters are ideal-space and the model was compared smeared.
        "slit_smearing_applied": bool(getattr(s.model, "use_slit_smearing", False)
                                       and getattr(s.model, "slit_length", 0.0) > 0),
        "slit_length":           float(getattr(s.model, "slit_length", 0.0)),
        **({"uncertainty": unc.mc_reply(mc_block)} if mc_block else {}),
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
    """Capture the current fit as a PNG image.

    Works before a fit is run (shows data + initial model) and after a fit
    (adds residuals subplot).  The layout mirrors pyirena's Unified Fit
    results window.

    Parameters
    ----------
    width, height : int
        Output image dimensions in pixels.
    format : str
        Only "png" is supported.

    Returns
    -------
    dict with ``image_base64`` (for inline display) and ``image_path`` (the
    same PNG on disk, for clients that cannot render images inline).
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)

    try:
        b64, image_path = _render_fit_image(s, width=width, height=height)
    except Exception as exc:
        return make_error(f"Image generation failed: {exc}", code="PLOT_ERROR")

    return {
        "image_base64": b64,
        "image_path":   image_path,
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


# ---------------------------------------------------------------------------
# Category F — Persistence
# ---------------------------------------------------------------------------

def _mc_levels(session: Session) -> Optional[dict]:
    """The MC uncertainties of the last fit in save_unified_fit_results' shape."""
    block = getattr(session, "last_uncertainty", None)
    if not block or not block.get("available"):
        return None
    return block.get("_levels")


def _session_to_gui_state(session: Session) -> dict:
    """Snapshot the session's UnifiedFitModel as the panel's setup state.

    The panel reads both dialects, and what is written is the core one —
    ``to_dict()``, the same shape the other five tools embed, which also
    carries ``K``, ``mass_fractal`` and ``RgCO``'s fit flag and bounds that
    the panel vocabulary has no room for. See ``planning/config-dialects/``
    §4: the panel dialect is read forever and written no longer.
    """
    model = session.model
    return {
        "schema_version":    UNIFIED_STATE_SCHEMA_VERSION,
        "num_levels":        int(model.num_levels),
        "levels":            [lv.to_dict() for lv in model.levels],
        "background":        float(model.background),
        "fit_background":    bool(model.fit_background),
        "background_limits": [float(v) for v in model.background_limits],
        "use_slit_smearing": bool(getattr(model, "use_slit_smearing", False)),
        "slit_length":       float(getattr(model, "slit_length", 0.0) or 0.0),
        "cursor_left":       session.fit_q_min,
        "cursor_right":      session.fit_q_max,
        "update_auto":       False,
        "display_local":     False,
        "no_limits":         False,
        "skip_fit_check":    False,
        "store_local":       False,
    }


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

    target, save_error = resolve_save_target(s, output_path)
    if save_error is not None:
        return save_error

    # Convert UnifiedLevel objects to the dict format expected by nxcansas_unified
    levels_dicts = [
        {
            "Rg": lv.Rg, "G": lv.G, "P": lv.P, "B": lv.B,
            "ETA": lv.ETA, "PACK": lv.PACK, "RgCO": lv.RgCO,
        }
        for lv in model.levels
    ]

    # Under slit smearing the curve that matches the data is the SMEARED model;
    # the ideal (pinhole) model is saved alongside it (parity with batch/GUI).
    slit_length = float(model.slit_length) if model.use_slit_smearing else 0.0
    intensity_model_ideal = None
    if slit_length > 0:
        intensity_model = model.calculate_intensity_smeared(model.q_data)
        intensity_model_ideal = model.calculate_intensity(model.q_data)
    else:
        intensity_model = model.calculate_intensity(model.q_data)

    resid = model.I_data - intensity_model
    norm_resid = (resid / model.error_data
                  if model.error_data is not None else resid)

    try:
        save_unified_fit_results(
            filepath=target,
            q=model.q_data,
            intensity_data=model.I_data,
            intensity_model=intensity_model,
            residuals=norm_resid,
            levels=levels_dicts,
            background=model.background,
            chi_squared=float(s.last_fit_result.get("chi_squared", float("nan"))),
            num_levels=model.num_levels,
            error=model.error_data,
            slit_length=slit_length,
            intensity_model_ideal=intensity_model_ideal,
            # Monte-Carlo uncertainties from run_fit(n_mc_runs=…), stored the
            # way the GUI's Store in File stores "Analyze Results".
            uncertainties=_mc_levels(s),
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


def _report_data_info(s) -> dict:
    """Session data arrays in the shape ``core.reporting`` expects."""
    return {"Q": s.q, "I": s.intensity, "I_error": s.error}


def _report_results_dict(s, chi_sq, quality) -> dict:
    """Session fit state as a ``load_unified_fit_results()``-shaped dict.

    Lets the control layer render its report with the same builder as the GUI
    instead of formatting a second, slightly different table of its own.
    """
    from datetime import datetime as _dt

    import numpy as _np

    model = s.model
    intensity_model = None
    residuals = None
    try:
        intensity_model = model.calculate_intensity(s.q)
        denom = s.error if s.error is not None else s.intensity
        residuals = (s.intensity - intensity_model) / denom
    except Exception:
        residuals = _np.zeros_like(s.q)

    levels = []
    for lv in model.levels[: model.num_levels]:
        level = {
            "G": lv.G, "Rg": lv.Rg, "B": lv.B, "P": lv.P,
            "RgCutoff": getattr(lv, "RgCO", 0.0),
            "correlated": bool(getattr(lv, "correlations", False)),
            "ETA": getattr(lv, "ETA", 0.0),
            "PACK": getattr(lv, "PACK", 0.0),
        }
        levels.append(level)

    return {
        "Q": s.q,
        "intensity_data": s.intensity,
        "intensity_model": intensity_model,
        "intensity_error": s.error,
        "residuals": residuals,
        "num_levels": model.num_levels,
        "background": float(getattr(model, "background", 0.0) or 0.0),
        "chi_squared": chi_sq if chi_sq is not None else 0.0,
        "levels": levels,
        "slit_length": (
            float(model.slit_length) if getattr(model, "use_slit_smearing", False) else 0.0
        ),
        "fit_quality": quality or None,
        "timestamp": _dt.now().strftime("%Y-%m-%d %H:%M:%S"),
    }


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
        # Same builder the GUI's "Copy results" button and the Data Selector's
        # Create Report use, so an agent and a user reading over their shoulder
        # see identical numbers and wording.
        from pyirena.core.reporting import build_report

        content = build_report(
            s.file_path or "results",
            data_info=_report_data_info(s),
            fit_results=_report_results_dict(s, chi_sq, quality),
        )
        # The per-parameter fit flags and bounds are an agent-only concern
        # (the GUI shows them as checkboxes), so they are appended rather than
        # pushed into the shared report.
        lines = [content, "## Fit Control", "",
                 "| Parameter | Value | Fixed | Lo | Hi | Units |",
                 "|-----------|-------|-------|----|----|-------|"]
        for p in params:
            lo_s = f"{p['lo']:.3g}" if p["lo"] is not None else "—"
            hi_s = f"{p['hi']:.3g}" if p["hi"] is not None else "—"
            lines.append(
                f"| {p['name']} | {p['value']:.4g} | {'Yes' if p['fixed'] else 'No'} "
                f"| {lo_s} | {hi_s} | {p['units']} |"
            )
        lines.append("")
        return {"format": "markdown", "content": "\n".join(lines)}

    else:
        return make_error(
            f"Unknown format '{format}'.",
            suggestion="Use 'json' or 'markdown'.",
            code="BAD_FORMAT",
        )
