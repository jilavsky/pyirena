"""Build a configured model from a pyIrena config section — one implementation.

A pyIrena config file (*Export Parameters* in the GUI, or hand-written) looks
like::

    {"_pyirena_config": {...metadata...},
     "modeling": {...the tool's settings...}}

Turning that section into a configured core model is something three callers
need: ``pyirena.batch`` (fit a file from a config), ``pyirena.api.control``
(``analyze()``, and the ZMQ service behind it), and anything else that replays
a saved setup. It used to live inline in each ``batch/<tool>.py``, which is
how the batch path and the GUI drifted apart once before — the comment in
``batch/unified.py`` saying so is why this module exists rather than a fourth
copy.

**A config is more than model state.** The fitted Q range, the slit settings
and (for Simple Fits) which parameters are held fixed are not fields of the
core model, and ``model.to_dict()`` rightly excludes them — but a fit replayed
without them is a different fit, silently. :class:`ToolSetup` carries the
model and those settings together, so a caller cannot use one and forget the
other.

Both config vocabularies are accepted. Five tools write the same names their
model's ``to_dict()`` uses; Unified Fit's config predates ``to_dict`` and
writes the panel's older names (``RgCutoff`` for ``RgCO``, and each parameter
as ``{"value": …, "fit": …}`` rather than a bare number). See
``planning/config-dialects/`` for why, and for the plan to converge them.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple

#: Config section name per tool, as written in a pyIrena config file. These are
#: also the ``model_name`` values used by ``pyirena.api.control`` sessions.
TOOL_SECTIONS: Tuple[str, ...] = (
    "unified_fit", "sizes", "simple_fits", "modeling", "waxs_peakfit", "carbon_fit",
)


class ConfigError(ValueError):
    """A config that cannot be turned into a model, with a readable reason."""


@dataclass
class ToolSetup:
    """A model plus the session settings a config carries alongside it."""

    tool: str
    model: Any

    #: Q range the fit should be restricted to (None = the full curve).
    fit_q_min: Optional[float] = None
    fit_q_max: Optional[float] = None

    #: Simple Fits only: {name: value} for parameters the user held fixed.
    fixed_params: Dict[str, float] = field(default_factory=dict)

    #: True when the config asks for the file's slit-smeared dataset.
    load_slit_smeared: bool = False

    #: Notes worth showing the caller (a clamped value, an ignored key).
    warnings: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# Locating the tool section
# ---------------------------------------------------------------------------

def detect_tool(config: Dict) -> Optional[str]:
    """Return the tool a config file is for, or None if it is unclear.

    A config written by the GUI holds exactly one tool section, so the tool
    need not be stated separately — which matters for ``analyze()``, where
    making the caller repeat it invites the two disagreeing.
    """
    if not isinstance(config, dict):
        return None
    present = [name for name in TOOL_SECTIONS if isinstance(config.get(name), dict)]
    if len(present) == 1:
        return present[0]
    stated = (config.get("_pyirena_config") or {}).get("tool")
    if stated in TOOL_SECTIONS:
        return stated
    if config.get("tool") in TOOL_SECTIONS:
        return config["tool"]
    return None


def section_for(config: Dict, tool: Optional[str] = None) -> Tuple[str, Dict]:
    """Return ``(tool, section)`` from a config file or a bare section.

    Accepts a whole config file, a bare tool section with *tool* given, or the
    ``{"tool": …, "model": …}`` envelope.
    """
    if not isinstance(config, dict):
        raise ConfigError(f"config must be an object, got {type(config).__name__}.")

    detected = detect_tool(config)
    if tool is None:
        tool = detected
    if tool is None:
        raise ConfigError(
            "Could not tell which tool this config is for. Expected one "
            f"section named one of: {', '.join(TOOL_SECTIONS)}."
        )
    if tool not in TOOL_SECTIONS:
        raise ConfigError(
            f"Unknown tool '{tool}'. One of: {', '.join(TOOL_SECTIONS)}."
        )

    section = config.get(tool)
    if isinstance(section, dict):
        return tool, section
    # A bare section was passed (no wrapper), or the envelope form.
    if isinstance(config.get("model"), dict) and detected is None:
        return tool, config
    if detected is None:
        return tool, config
    raise ConfigError(f"Config has no '{tool}' section.")


# ---------------------------------------------------------------------------
# Unified Fit — the one tool whose config speaks the older panel vocabulary
# ---------------------------------------------------------------------------

#: Defaults for a bound the config file leaves out, matching the wide bounds
#: the GUI ships with.
_LIMIT_DEFAULTS = {
    "G": (0.0, 1e10), "Rg": (0.1, 1e6), "B": (0.0, 1e10),
    "P": (0.0, 6.0), "ETA": (0.1, 1e6), "PACK": (0.0, 16.0),
}


def flatten_level_config(ls: Dict) -> Dict:
    """Turn one level of the JSON config into the panel's flat key names.

    Accepts both spellings a config file may use for a parameter: the nested
    ``{'value': …, 'fit': …, 'low_limit': …, 'high_limit': …}`` form written by
    the GUI, and a bare number for hand-written configs.
    """
    flat: Dict = {
        "RgCutoff":   float(ls.get("RgCutoff", ls.get("RgCO", 0.0)) or 0.0),
        "correlated": bool(ls.get("correlated", ls.get("correlations", False))),
        "estimate_B": bool(ls.get("estimate_B", ls.get("link_B", False))),
        "link_rgco":  bool(ls.get("link_rgco", ls.get("link_RGCO", False))),
    }
    for name, (lo_default, hi_default) in _LIMIT_DEFAULTS.items():
        entry = ls.get(name, {})
        if not isinstance(entry, dict):
            flat[name] = float(entry)
            continue
        flat[name] = float(entry.get("value", 0.0))
        flat[f"fit_{name}"] = bool(entry.get("fit", False))
        flat[f"{name}_low"] = float(entry.get("low_limit", lo_default))
        flat[f"{name}_high"] = float(entry.get("high_limit", hi_default))
    return flat


def unified_model_from_config(state: Dict):
    """Convert the ``unified_fit`` config section into a configured model."""
    from pyirena.core.unified import UnifiedFitModel, UnifiedLevel

    num_levels = state.get("num_levels", 1)
    no_limits = state.get("no_limits", False)

    model = UnifiedFitModel(num_levels=num_levels)

    # The panel writes {'value': …, 'fit': …}; the core dict writes a number.
    bg = state.get("background", {})
    if isinstance(bg, dict):
        model.background = float(bg.get("value", 0.0))
        model.fit_background = bool(bg.get("fit", False))
    else:
        model.background = float(bg or 0.0)
        model.fit_background = bool(state.get("fit_background", False))

    for i, ls in enumerate(state.get("levels", [])[:num_levels]):
        # The config file stores each parameter as a small dict; the core model
        # speaks the panel's flat key names. Flatten here, then let
        # UnifiedLevel.from_panel_params own the field mapping — it is the same
        # translation the GUI does, and keeping two copies is how the batch path
        # and the GUI drifted apart before.
        model.levels[i] = UnifiedLevel.from_panel_params(
            flatten_level_config(ls), with_limits=not no_limits)

    return model


# ---------------------------------------------------------------------------
# Size Distribution
# ---------------------------------------------------------------------------

def sizes_model_from_config(state: Dict, *, data_is_slit_smeared: bool = False,
                            data_slit_length: float = 0.0):
    """Convert the ``sizes`` config section into a configured SizesDistribution."""
    from pyirena.core.sizes import SizesDistribution

    s = SizesDistribution()
    s.r_min              = float(state.get("r_min", 10.0))
    s.r_max              = float(state.get("r_max", 1000.0))
    s.n_bins             = int(state.get("n_bins", 200))
    s.log_spacing        = bool(state.get("log_spacing", True))
    s.shape              = str(state.get("shape", "sphere"))
    s.contrast           = float(state.get("contrast", 1.0))
    ar = state.get("aspect_ratio", 1.0)
    if s.shape == "spheroid":
        s.shape_params = {"aspect_ratio": float(ar)}
    elif isinstance(state.get("shape_params"), dict):
        s.shape_params = dict(state["shape_params"])
    s.background         = float(state.get("background", 0.0))
    s.error_scale        = float(state.get("error_scale", 1.0))
    s.fractional_error   = bool(state.get("fractional_error", False))
    s.fractional_error_value = float(state.get("fractional_error_value", 0.03))
    s.power_law_B        = float(state.get("power_law_B", 0.0))
    s.power_law_P        = float(state.get("power_law_P", 4.0))
    s.method             = str(state.get("method", "regularization"))
    s.maxent_sky_background  = float(state.get("maxent_sky_background", 1e-6))
    s.maxent_stability       = float(state.get("maxent_stability", 0.01))
    s.maxent_max_iter        = int(state.get("maxent_max_iter", 300))
    s.regularization_evalue  = float(state.get("regularization_evalue", 1.0))
    s.regularization_min_ratio = float(state.get("regularization_min_ratio", 1e-4))
    s.tnnls_approach_param   = float(state.get("tnnls_approach_param", 0.95))
    s.tnnls_max_iter         = int(state.get("tnnls_max_iter", 300))
    # The main fit always uses a single MC run, matching the GUI.
    s.montecarlo_n_repetitions = 1

    # Slit smearing: enable when the loaded data are slit smeared or the config
    # asks for it; slit length is file-derived unless overridden.
    cfg_sl = state.get("slit_length")
    sl = float(cfg_sl) if cfg_sl else float(data_slit_length or 0.0)
    if (bool(data_is_slit_smeared) or bool(state.get("use_slit_smearing"))) and sl > 0:
        s.use_slit_smearing = True
        s.slit_length = sl
    return s


# ---------------------------------------------------------------------------
# Simple Fits
# ---------------------------------------------------------------------------

def simple_model_from_config(state: Dict) -> Tuple[Any, Dict[str, float]]:
    """Convert the ``simple_fits`` config section into ``(model, fixed_params)``.

    The GUI stores the per-parameter "Fit?" checkboxes as
    ``param_fixed = {name: True when held}``, while ``SimpleFitModel.fit()``
    wants ``fixed_params = {name: value}``. Without that translation a replay
    refits every parameter and quietly ignores the user's choices.
    """
    from pyirena.core.simple_fits import SimpleFitModel

    cfg = dict(state)
    cfg.pop("q_min", None)
    cfg.pop("q_max", None)

    # The GUI state uses 'param_limits'; SimpleFitModel.from_dict() wants 'limits'.
    if "param_limits" in cfg and "limits" not in cfg:
        cfg["limits"] = cfg.pop("param_limits")
    else:
        cfg.pop("param_limits", None)

    param_fixed = cfg.pop("param_fixed", {}) or {}
    params = cfg.get("params", {}) or {}
    fixed_params = {
        name: params[name]
        for name, is_fixed in param_fixed.items()
        if is_fixed and name in params
    }

    for key in ("schema_version", "no_limits"):
        cfg.pop(key, None)

    return SimpleFitModel.from_dict(cfg), fixed_params


# ---------------------------------------------------------------------------
# The three tools whose config already is their model's to_dict()
# ---------------------------------------------------------------------------

def modeling_model_from_config(state: Dict):
    from pyirena.core.modeling import ModelingConfig
    return ModelingConfig.from_dict(state)


def waxs_model_from_config(state: Dict):
    from pyirena.core.waxs_peakfit import WAXSPeakFitModel
    return WAXSPeakFitModel.from_dict(state)


def carbon_model_from_config(state: Dict):
    from pyirena.core.carbon_fit import CarbonFitModel
    # The GUI nests the model under 'model' alongside its own view state; a
    # bare model dict is accepted too.
    return CarbonFitModel.from_dict(state.get("model", state))


# ---------------------------------------------------------------------------
# One entry point
# ---------------------------------------------------------------------------

def build_setup(config: Dict, tool: Optional[str] = None, *,
                data_is_slit_smeared: bool = False,
                data_slit_length: float = 0.0) -> ToolSetup:
    """Build a :class:`ToolSetup` from a pyIrena config file or tool section.

    Parameters
    ----------
    config : dict
        A whole config file, or one tool's section with *tool* given.
    tool : str, optional
        Which tool; inferred from the config when it holds one tool section.
    data_is_slit_smeared, data_slit_length :
        Properties of the loaded curve. A config may enable slit smearing
        without stating a length, in which case the file's value is used.

    Raises
    ------
    ConfigError
        With a message naming what was wrong and what was expected.
    """
    tool, section = section_for(config, tool)
    setup = ToolSetup(tool=tool, model=None)

    try:
        if tool == "unified_fit":
            setup.model = unified_model_from_config(section)
            setup.fit_q_min = section.get("cursor_left")
            setup.fit_q_max = section.get("cursor_right")
        elif tool == "sizes":
            setup.model = sizes_model_from_config(
                section, data_is_slit_smeared=data_is_slit_smeared,
                data_slit_length=data_slit_length)
            setup.fit_q_min = section.get("cursor_q_min")
            setup.fit_q_max = section.get("cursor_q_max")
        elif tool == "simple_fits":
            setup.model, setup.fixed_params = simple_model_from_config(section)
            setup.fit_q_min = section.get("q_min")
            setup.fit_q_max = section.get("q_max")
        elif tool == "modeling":
            setup.model = modeling_model_from_config(section)
            setup.fit_q_min = section.get("q_min")
            setup.fit_q_max = section.get("q_max")
        elif tool == "waxs_peakfit":
            setup.model = waxs_model_from_config(section)
            setup.fit_q_min = section.get("q_min")
            setup.fit_q_max = section.get("q_max")
        elif tool == "carbon_fit":
            setup.model = carbon_model_from_config(section)
            setup.fit_q_min = section.get("q_min")
            setup.fit_q_max = section.get("q_max")
    except ConfigError:
        raise
    except Exception as exc:
        raise ConfigError(f"Could not build a {tool} model from this config: {exc}") from exc

    setup.load_slit_smeared = bool(section.get("load_slit_smeared", False))
    setup.fit_q_min = float(setup.fit_q_min) if setup.fit_q_min is not None else None
    setup.fit_q_max = float(setup.fit_q_max) if setup.fit_q_max is not None else None
    return setup


# ---------------------------------------------------------------------------
# Pre-fit steps a config asks for
# ---------------------------------------------------------------------------
#
# Two tools do preparatory work before the fit proper, driven by config keys
# rather than by model state, and skipping it changes the answer without
# changing anything visible. Size Distribution pre-fits the power-law and flat
# background over Q windows the user picked; WAXS Peak Fit re-centres peaks by
# scanning or cross-correlating Q0 before refining. A replay that omits them
# returns a plausible number, not an error — which is the worst kind.
#
# Both pre-fits deliberately use the FULL curve, not the fitted Q range: the
# windows are chosen independently of the fit cursors, exactly as the GUI does.

def apply_sizes_prefits(model, section: Dict, q, intensity) -> list:
    """Run the Size Distribution pre-fits a config requests. Returns notes."""
    import numpy as np

    notes: list = []
    fit_B = bool(section.get("fit_power_law_B", False))
    fit_P = bool(section.get("fit_power_law_P", False))
    if fit_B or fit_P:
        lo = section.get("power_law_q_min")
        hi = section.get("power_law_q_max")
        if lo is None or hi is None:
            lo, hi = float(np.min(q)), float(np.max(q))
        try:
            result = model.fit_power_law(q, intensity, float(lo), float(hi),
                                         fit_B=fit_B, fit_P=fit_P)
            notes.append(f"Power-law pre-fit: {result.get('message', '')}".strip())
        except Exception as exc:
            notes.append(f"Power-law pre-fit failed, using the configured values: {exc}")

    lo = section.get("background_q_min")
    hi = section.get("background_q_max")
    if lo is not None and hi is not None:
        try:
            result = model.fit_background_term(q, intensity, float(lo), float(hi))
            notes.append(f"Background pre-fit: {result.get('message', '')}".strip())
        except Exception as exc:
            notes.append(f"Background pre-fit failed, using the configured value: {exc}")
    return notes


def apply_waxs_presearch(model, section: Dict, q, intensity) -> list:
    """Re-centre WAXS peaks as the config's ``presearch`` block asks."""
    from pyirena.core.waxs_peakfit import (
        cross_corr_q_shift,
        eval_model,
        presearch_q0_per_peak,
    )

    ps = section.get("presearch") or {}
    run_cc = bool(ps.get("cross_corr", False))
    run_scan = bool(ps.get("per_peak_scan", False))
    if not model.peaks or not (run_cc or run_scan):
        return []

    window = float(ps.get("search_window", 0.05))
    steps = int(ps.get("n_steps", 50))
    notes: list = []

    if run_cc:
        I_model = eval_model(q, model.bg_shape, model.bg_params, model.peaks, I=intensity)
        shift = cross_corr_q_shift(q, intensity, I_model, max_shift=window)
        if abs(shift) > 1e-6:
            # Only shift peaks whose Q0 is free; a locked Q0 stays where the
            # user put it.
            for peak in model.peaks:
                if bool(peak.get("Q0", {}).get("fit", True)):
                    peak["Q0"]["value"] = float(peak["Q0"]["value"]) + shift
            notes.append(f"Cross-correlation shift applied: {shift:+.4f} 1/A.")

    if run_scan:
        model.peaks = presearch_q0_per_peak(
            q, intensity, model.bg_shape, model.bg_params, model.peaks,
            search_window=window, n_steps=steps,
        )
        notes.append(f"Per-peak Q0 scan done (window=±{window:.3f} 1/A, {steps} steps).")
    return notes


def apply_prefits(setup: "ToolSetup", section: Dict,
                  q_full, intensity_full, q_fit=None, intensity_fit=None) -> list:
    """Run whatever pre-fit steps *setup*'s tool defines. Returns notes.

    The two tools disagree about which data a pre-fit should see, and both are
    right. The Size Distribution background windows are chosen independently
    of the size-fit cursors, so they read the **full** curve. The WAXS peak
    presearch is re-centring the peaks that are about to be fitted, so it
    reads only the **fitted** range — scanning Q0 across data the fit will
    never see would move a peak towards a feature outside the window.
    """
    if q_fit is None:
        q_fit, intensity_fit = q_full, intensity_full
    if setup.tool == "sizes":
        return apply_sizes_prefits(setup.model, section, q_full, intensity_full)
    if setup.tool == "waxs_peakfit":
        return apply_waxs_presearch(setup.model, section, q_fit, intensity_fit)
    return []
