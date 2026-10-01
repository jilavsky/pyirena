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

**Both config vocabularies are read; only one is written.** Every tool now
writes the core dialect — the names its ``to_dict()`` uses. Unified Fit's
config predated ``to_dict`` and wrote the panel's older names (``RgCutoff``
for ``RgCO``, each parameter a ``{"value": …, "fit": …}`` block); those files
are read forever and never written again. ``is_panel_level`` tells the two
apart by shape, in one place, and ``unified_level_from_config`` is the only
thing that chooses a reader.

**Four envelopes, one reader.** ``unwrap_config`` accepts the *Export
Parameters* sidecar, the ``_pyirena_config`` attribute inside a result file,
an ``export_results()`` reply and the bare ``{tool, model}`` form. For most
of pyIrena's life only the first was accepted, so the two a remote caller
actually held were the two that could not be fed back in.

See ``planning/config-dialects/`` for the full inventory and
``docs/batch_api.md`` for the user-facing description.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple

#: Config section name per tool, as written in a pyIrena config file. These are
#: also the ``model_name`` values used by ``pyirena.api.control`` sessions.
TOOL_SECTIONS: Tuple[str, ...] = (
    "unified_fit", "sizes", "simple_fits", "modeling", "waxs_peakfit", "carbon_fit",
)


#: The header key every pyIrena config wrapper carries. Named here so the
#: reader and ``io/setup_config.py``'s writer cannot drift apart.
PYIRENA_CONFIG_HEADER = "_pyirena_config"


class ConfigError(ValueError):
    """A config that cannot be turned into a model, with a readable reason."""


@dataclass
class ToolSetup:
    """A model plus the session settings a config carries alongside it."""

    tool: str
    model: Any

    #: The tool's settings section, as resolved out of whichever wrapper the
    #: caller passed. Execution needs it as much as model construction does —
    #: ``no_limits``, the fit method and the weighting live beside the model,
    #: not inside it, and a caller that re-reads the wrapper with its own
    #: shallow lookup sees the wrapper instead of the section. Carrying the
    #: one ``unwrap_config`` already resolved is what keeps the two in step.
    section: Dict = field(default_factory=dict)

    #: Q range the fit should be restricted to (None = the full curve).
    fit_q_min: Optional[float] = None
    fit_q_max: Optional[float] = None

    #: Simple Fits only: {name: value} for parameters the user held fixed.
    fixed_params: Dict[str, float] = field(default_factory=dict)

    #: True when the config asks for the file's slit-smeared dataset.
    load_slit_smeared: bool = False

    #: Size Distribution only: the config's uncertainty-run count
    #: (``unc_n_runs``), or None when it does not state one. It is a GUI
    #: control the user sets and the replay path used to ignore outright.
    mc_n_runs: Optional[int] = None

    #: Notes worth showing the caller (a clamped value, an ignored key).
    warnings: list = field(default_factory=list)


#: Keys that describe the *view*, not the fit, and so do not belong in a
#: config file meant to be shared or replayed.
#:
#: ``last_folder`` is the one that is actively wrong: it is an absolute path
#: from whichever machine exported the file, and it travelled into every
#: shared config and every result file. The rest are harmless but noise — a
#: tab index and a couple of checkbox states cannot change a number, and a
#: reader has to satisfy itself of that one key at a time.
#:
#: These stay in the StateManager, which is where a user's view preferences
#: belong; only *Export Parameters* strips them. See
#: ``planning/config-dialects/`` §2.5 and §5 Step 3.4.
VIEW_ONLY_CONFIG_KEYS = frozenset({
    "last_folder", "active_tab", "waxs_zoom_visible", "show_components",
    "auto_update", "update_auto", "display_local", "store_local",
})


def without_view_state(section: Dict) -> Dict:
    """A tool section with the view-only keys of :data:`VIEW_ONLY_CONFIG_KEYS`
    removed. Non-destructive; the caller's dict is untouched."""
    if not isinstance(section, dict):
        return section
    return {k: v for k, v in section.items() if k not in VIEW_ONLY_CONFIG_KEYS}


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


@dataclass
class _Envelope:
    """What a wrapper carries besides the tool's own settings section.

    Only the fields a wrapper states; ``None`` means "this envelope says
    nothing about it", which is different from "it says zero".
    """

    fit_q_min: Optional[float] = None
    fit_q_max: Optional[float] = None
    is_slit_smeared: Optional[bool] = None
    slit_length: Optional[float] = None

    #: Simple Fits' held-parameter choices, in whichever spelling the wrapper
    #: used. They are not model state, so they only reach a replay if the
    #: wrapper states them (see :func:`resolve_fixed_params`).
    fixed_params: Any = None


def _q_range_of(block: Any) -> Tuple[Optional[float], Optional[float]]:
    """Read a ``{"q_min": …, "q_max": …}`` block, tolerating None and junk."""
    if not isinstance(block, dict):
        return None, None

    def _num(key):
        value = block.get(key)
        try:
            return None if value is None else float(value)
        except (TypeError, ValueError):
            return None

    return _num("q_min"), _num("q_max")


def unwrap_config(config: Dict, tool: Optional[str] = None
                  ) -> Tuple[str, Dict, _Envelope]:
    """Return ``(tool, section, envelope)`` for any wrapper pyIrena writes.

    pyIrena writes a tool's settings inside four different wrappers, and for
    a long time read only the first of them — so the two a remote caller
    actually holds, the one inside a result file and the one the service just
    handed back, were the two that could not be fed back in:

    ``{"_pyirena_config": …, "<tool>": {…}}``
        The *Export Parameters* JSON sidecar.
    ``{"_pyirena_config": …, "state": {…}}``
        The ``_pyirena_config`` attribute on every HDF5 results group.
    ``{"ok": …, "tool": …, "config": {…}, "fit_q_range": …, "data": …}``
        Every ``export_results()`` reply, over MCP and the ZMQ service.
    ``{"tool": …, "model": {…}}``
        The minimal envelope, for a caller assembling one by hand.

    A bare tool section with *tool* given is accepted too. The settings are
    the same object in all five cases; only the wrapper differs.
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

    env = _Envelope()

    # 1. The sidecar: the section is named after the tool.
    section = config.get(tool)
    if isinstance(section, dict):
        return tool, section, env

    # 2. The HDF5 attribute: {_pyirena_config, state}. Only a wrapper has a
    #    '_pyirena_config' header, so this cannot collide with a bare section.
    state = config.get("state")
    if isinstance(state, dict) and isinstance(config.get(PYIRENA_CONFIG_HEADER), dict):
        return tool, state, env

    # 3. An export_results() reply. Its 'config' is the model's to_dict(), and
    #    it states the fitted Q range and the curve's slit settings separately
    #    — which is the whole reason the reply is worth feeding back in.
    reply = config.get("config")
    if isinstance(reply, dict) and config.get("tool") in TOOL_SECTIONS:
        env.fit_q_min, env.fit_q_max = _q_range_of(config.get("fit_q_range"))
        if config.get("fixed_params") is not None:
            env.fixed_params = config["fixed_params"]
        data = config.get("data")
        if isinstance(data, dict):
            if data.get("is_slit_smeared") is not None:
                env.is_slit_smeared = bool(data["is_slit_smeared"])
            try:
                if data.get("slit_length") is not None:
                    env.slit_length = float(data["slit_length"])
            except (TypeError, ValueError):
                pass
        return tool, reply, env

    # 4. The minimal {tool, model} envelope. Guarded on an explicit 'tool'
    #    key: the Carbon section nests its own 'model', and reading that as
    #    this envelope would drop the Q range beside it.
    model = config.get("model")
    if isinstance(model, dict) and config.get("tool") in TOOL_SECTIONS:
        env.fit_q_min, env.fit_q_max = _q_range_of(config.get("fit_q_range"))
        return tool, model, env

    # 5. A bare section, passed with the tool named.
    if detected is None or detected == tool:
        return tool, config, env

    raise ConfigError(f"Config has no '{tool}' section.")


def section_for(config: Dict, tool: Optional[str] = None) -> Tuple[str, Dict]:
    """``(tool, section)`` from any wrapper — :func:`unwrap_config` without
    the envelope's own fields. Kept for callers that only want the section."""
    tool, section, _ = unwrap_config(config, tool)
    return tool, section


# ---------------------------------------------------------------------------
# Unified Fit — the one tool whose config speaks the older panel vocabulary
# ---------------------------------------------------------------------------

#: Defaults for a bound the config file leaves out, matching the wide bounds
#: the GUI ships with.
_LIMIT_DEFAULTS = {
    "G": (0.0, 1e10), "Rg": (0.1, 1e6), "B": (0.0, 1e10),
    "P": (0.0, 6.0), "ETA": (0.1, 1e6), "PACK": (0.0, 16.0),
}


#: Keys that only the panel dialect ever writes.  A level carrying one of
#: these is panel-shaped even when its parameters happen to be bare numbers.
_PANEL_ONLY_LEVEL_KEYS = (
    "RgCutoff", "correlated", "estimate_B", "link_rgco", "level",
    "G_low", "Rg_low", "B_low", "P_low", "ETA_low", "PACK_low",
    "G_high", "Rg_high", "B_high", "P_high", "ETA_high", "PACK_high",
)


def is_panel_level(ls: Dict) -> bool:
    """True when one level dict speaks the *panel* dialect, not the core one.

    The two are told apart by shape, which is reliable because they disagree
    about a parameter's type: the panel writes ``{"value": …, "fit": …}`` and
    the core writes a bare number beside ``fit_Rg`` and ``Rg_limits``.  A
    hand-written level with bare numbers and none of the panel's own key
    names reads as core, which is the dialect it is closest to and the one
    that carries flags and bounds.
    """
    if not isinstance(ls, dict):
        return False
    if any(isinstance(ls.get(name), dict) for name in _LIMIT_DEFAULTS):
        return True
    return any(key in ls for key in _PANEL_ONLY_LEVEL_KEYS)


def flatten_level_config(ls: Dict) -> Dict:
    """Turn one level of the JSON config into the panel's flat key names.

    Accepts every spelling a config file may use for a parameter: the nested
    ``{'value': …, 'fit': …, 'low_limit': …, 'high_limit': …}`` form written by
    the GUI, and a bare number — which, in the core dialect, is accompanied by
    its own ``fit_<name>`` flag and ``<name>_limits`` pair.  Reading a bare
    number and *skipping* those two is what silently freed parameters the
    scientist had pinned; see ``planning/config-dialects/`` §2.1.

    Levels that are wholly core-shaped go through :meth:`UnifiedLevel.from_dict`
    instead (see :func:`unified_level_from_config`), which also carries ``K``,
    ``mass_fractal`` and the ``RgCO`` flag and bounds that the panel vocabulary
    has no room for.  This function still has to handle the mixed case.
    """
    flat: Dict = {
        "RgCutoff":   float(ls.get("RgCutoff", ls.get("RgCO", 0.0)) or 0.0),
        "correlated": bool(ls.get("correlated", ls.get("correlations", False))),
        "estimate_B": bool(ls.get("estimate_B", ls.get("link_B", False))),
        "link_rgco":  bool(ls.get("link_rgco", ls.get("link_RGCO", False))),
    }
    for name, (lo_default, hi_default) in _LIMIT_DEFAULTS.items():
        entry = ls.get(name, {})
        if isinstance(entry, dict):
            flat[name] = float(entry.get("value", 0.0) or 0.0)
            fit = entry.get("fit", False)
            # A stated-but-null bound is how the shipped state file spells
            # "no bound here"; treat it as absent rather than crashing on
            # float(None).
            lo = entry.get("low_limit")
            hi = entry.get("high_limit")
            lo = lo_default if lo is None else lo
            hi = hi_default if hi is None else hi
        else:
            # Core dialect: the value is bare and its flag and bounds sit
            # beside it under their own names.
            flat[name] = float(entry)
            fit = ls.get(f"fit_{name}")
            pair = ls.get(f"{name}_limits")
            if isinstance(pair, (list, tuple)) and len(pair) == 2:
                lo, hi = pair
            else:
                lo, hi = ls.get(f"{name}_low"), ls.get(f"{name}_high")
            lo = lo_default if lo is None else lo
            hi = hi_default if hi is None else hi
            if fit is None:
                # Neither dialect stated a flag: leave the model's default
                # alone rather than inventing False.
                flat[f"{name}_low"] = float(lo)
                flat[f"{name}_high"] = float(hi)
                continue
        flat[f"fit_{name}"] = bool(fit)
        flat[f"{name}_low"] = float(lo)
        flat[f"{name}_high"] = float(hi)
    return flat


def unified_level_from_config(ls: Dict, *, with_limits: bool = True):
    """Build one :class:`UnifiedLevel` from a config level in either dialect.

    One translator, chosen by shape — the alternative is every reader
    hand-rolling its own detection, which is how the core dialect ended up
    half-read in the first place.
    """
    from pyirena.core.unified import UnifiedLevel

    if is_panel_level(ls):
        return UnifiedLevel.from_panel_params(
            flatten_level_config(ls), with_limits=with_limits)

    level = UnifiedLevel.from_dict(ls)
    if not with_limits:
        # "No limits" mode wants the wide defaults, whatever the file says.
        default = UnifiedLevel()
        for name in UnifiedLevel._LIMIT_FIELDS:
            setattr(level, name, getattr(default, name))
    return level


def unified_model_from_config(state: Dict):
    """Convert the ``unified_fit`` config section into a configured model."""
    from pyirena.core.unified import UnifiedFitModel

    num_levels = state.get("num_levels", 1)
    no_limits = state.get("no_limits", False)

    model = UnifiedFitModel(num_levels=num_levels)

    # The panel writes {'value': …, 'fit': …}; the core dict writes a number.
    bg = state.get("background", {})
    if isinstance(bg, dict):
        model.background = float(bg.get("value", 0.0))
        model.fit_background = bool(bg.get("fit", False))
        lo, hi = bg.get("low_limit"), bg.get("high_limit")
        if lo is not None and hi is not None and not no_limits:
            model.background_limits = (float(lo), float(hi))
    else:
        model.background = float(bg or 0.0)
        if state.get("fit_background") is not None:
            model.fit_background = bool(state["fit_background"])
        limits = state.get("background_limits")
        if isinstance(limits, (list, tuple)) and len(limits) == 2 and not no_limits:
            model.background_limits = (float(limits[0]), float(limits[1]))

    # Slit smearing is model state in the core dialect; the panel keeps it
    # beside the model, so the caller supplies it (see build_setup).
    if state.get("use_slit_smearing") is not None:
        model.use_slit_smearing = bool(state["use_slit_smearing"])
    if state.get("slit_length") is not None:
        model.slit_length = float(state["slit_length"])

    for i, ls in enumerate(state.get("levels", [])[:num_levels]):
        # Either dialect; unified_level_from_config picks the reader by shape
        # and both of them are the same translators the GUI uses. Keeping a
        # second copy here is how the batch path and the GUI drifted apart
        # before.
        model.levels[i] = unified_level_from_config(ls, with_limits=not no_limits)

    return model


# ---------------------------------------------------------------------------
# Slit smearing — a property of the measurement, not only of the config
# ---------------------------------------------------------------------------

def apply_data_slit_settings(model, section: Dict, *,
                             data_is_slit_smeared: bool = False,
                             data_slit_length: float = 0.0) -> None:
    """Enable model slit smearing from the config *and* the loaded curve.

    The session API does this in every ``select_*_model()``: a curve that is
    slit smeared must be compared against a smeared model, or the fitted
    parameters are wrong in a way nothing downstream flags. A replayed config
    is not where that decision lives — it is a property of the measurement in
    front of us, and the config may well have been exported from a pinhole
    run. So either source can switch smearing on, and neither switches it off:

    * the config's ``use_slit_smearing``, as the GUI exported it;
    * the data's own ``is_slit_smeared``, which is what the normal session
      path reads.

    Whichever source asked for the smearing names the length, because that is
    the one that knows it. A config that asks for smearing was written by a
    scientist who set the length, so it wins; a config that does not ask is
    only carrying a stale number beside an unticked box, and the curve the
    file declares is authoritative. Either falls back to the other when it
    states no length. Does nothing for a model that cannot smear (Carbon,
    WAXS).
    """
    if not hasattr(model, "use_slit_smearing"):
        return
    section = section if isinstance(section, dict) else {}

    def _length(value):
        try:
            return float(value) if value else 0.0
        except (TypeError, ValueError):
            return 0.0

    cfg_length = _length(section.get("slit_length"))
    data_length = _length(data_slit_length)

    if section.get("use_slit_smearing"):
        length = cfg_length or data_length
    elif data_is_slit_smeared:
        length = data_length or cfg_length
    else:
        return

    if length > 0:
        model.use_slit_smearing = True
        model.slit_length = length


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
    # Two dialects spell the same thing: the panel writes a flat
    # 'aspect_ratio', the core writes it inside 'shape_params'. Reading only
    # the flat key replayed a core-dialect spheroid config as a sphere — the
    # one shape where the aspect ratio matters. shape_params is splatted as
    # keyword arguments into the shape's G-matrix builder, so 'aspect_ratio'
    # is carried for the spheroid alone; the panel writes it whatever the
    # shape, and a sphere builder has no such argument.
    shape_params = dict(state["shape_params"]) if isinstance(
        state.get("shape_params"), dict) else {}
    if s.shape == "spheroid":
        ar = state.get("aspect_ratio", shape_params.get("aspect_ratio", 1.0))
        shape_params["aspect_ratio"] = float(ar)
    else:
        shape_params.pop("aspect_ratio", None)
    if shape_params:
        s.shape_params = shape_params
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
    apply_data_slit_settings(s, state,
                             data_is_slit_smeared=data_is_slit_smeared,
                             data_slit_length=data_slit_length)
    return s


# ---------------------------------------------------------------------------
# Simple Fits
# ---------------------------------------------------------------------------

def resolve_fixed_params(spec: Any, params: Dict) -> Dict[str, float]:
    """Normalise *which parameters are held fixed* into ``{name: value}``.

    Three parts of pyIrena spell the same choice three ways, and a reader that
    knows only one of them silently frees what the scientist pinned:

    ``{"Rg": true}``
        ``param_fixed`` — the GUI's per-parameter "Fit?" checkboxes, written
        into the *Export Parameters* sidecar.
    ``["Rg"]``
        ``fixed_params`` — the list ``save_simple_fit()`` embeds in the setup
        attribute of a result ``.h5``, and the one ``export_results()``
        reports. Already-written files use it, so it is read forever.
    ``{"Rg": 12.0}``
        The resolved form :meth:`SimpleFitModel.fit` itself takes, accepted so
        a caller can hand back what it was given.

    The held *value* always comes from the model's own parameters: a name with
    no parameter behind it is dropped rather than invented.
    """
    params = params or {}
    if not spec:
        return {}
    if isinstance(spec, dict):
        held = []
        for name, flag in spec.items():
            # {name: True} says "held"; {name: 12.0} is the resolved form, in
            # which every listed name is held and 0.0 is a real value.
            if isinstance(flag, bool):
                if flag:
                    held.append(name)
            elif flag is not None:
                held.append(name)
    elif isinstance(spec, (list, tuple, set, frozenset)):
        held = list(spec)
    else:
        return {}
    return {name: params[name] for name in held if name in params}


def simple_model_from_config(state: Dict) -> Tuple[Any, Dict[str, float]]:
    """Convert the ``simple_fits`` config section into ``(model, fixed_params)``.

    The held-parameter choices are not model state, so they travel beside the
    model in whichever spelling the writer used; :func:`resolve_fixed_params`
    reads all of them. Without that translation a replay refits every
    parameter and quietly ignores the user's choices.
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

    param_fixed = cfg.pop("param_fixed", None)
    saved_fixed = cfg.pop("fixed_params", None)
    params = cfg.get("params", {}) or {}
    fixed_params = resolve_fixed_params(
        param_fixed if param_fixed else saved_fixed, params)

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
    tool, section, envelope = unwrap_config(config, tool)
    setup = ToolSetup(tool=tool, model=None, section=section)

    # An envelope that states the curve's slit settings speaks for the data
    # the caller did not hand over separately (an export_results reply is
    # the case that matters: the data are not in the request).
    if envelope.is_slit_smeared is not None and not data_is_slit_smeared:
        data_is_slit_smeared = envelope.is_slit_smeared
    if envelope.slit_length is not None and not data_slit_length:
        data_slit_length = envelope.slit_length

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

    # Slit smearing is settled here for every tool that can smear, not in the
    # Sizes branch alone. The three analytic tools used to read the config
    # side only, so a slit-smeared curve replayed through a pinhole config
    # was fitted unsmeared — the one case where the normal session path and
    # a replay of the same setup disagreed about the physics.
    apply_data_slit_settings(setup.model, section,
                             data_is_slit_smeared=data_is_slit_smeared,
                             data_slit_length=data_slit_length)

    # The wrapper may carry the held-parameter choices the section does not:
    # an export_results() reply states them beside the model, because
    # model.to_dict() has no room for them.
    if tool == "simple_fits" and not setup.fixed_params and envelope.fixed_params:
        setup.fixed_params = resolve_fixed_params(
            envelope.fixed_params, getattr(setup.model, "params", {}))

    # The envelope's own fit_q_range wins over whatever the section spells
    # it: four tools keep the Q range beside the model and two keep it
    # inside to_dict(), so the envelope is the single place a reader looks.
    # See planning/config-dialects/ §4 rule 4.
    if envelope.fit_q_min is not None:
        setup.fit_q_min = envelope.fit_q_min
    if envelope.fit_q_max is not None:
        setup.fit_q_max = envelope.fit_q_max

    setup.load_slit_smeared = bool(section.get("load_slit_smeared", False))
    if section.get("unc_n_runs") is not None:
        try:
            setup.mc_n_runs = max(1, int(section["unc_n_runs"]))
        except (TypeError, ValueError):
            pass
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


def apply_simple_bg_prefit(model, fixed_params: Optional[Dict], q, intensity) -> list:
    """Replay the Simple Fits complex-background pre-fit. Returns notes.

    Only for a *calculation* model — the Invariant. An ordinary fit refines
    its background as part of the optimisation, so a stale starting value
    costs nothing; the Invariant is an integration, and nothing after this
    point ever touches the background again. Replaying the config without it
    integrates the curve on top of whatever background the file happened to
    be exported with, and the result is wrong by orders of magnitude while
    still reporting success.

    Reads the FULL curve: the saved background windows are normally outside
    the integration range, which is the whole reason they were recorded.
    """
    if not getattr(model, "is_calculation", False):
        return []
    if not (getattr(model, "bg_prefit", None) or {}).get("enabled"):
        return []
    try:
        applied = model.prefit_background(q, intensity,
                                          fixed_params=fixed_params or None)
    except Exception as exc:
        return [f"Background pre-fit failed, using the configured values: {exc}"]

    notes = []
    values = "  ".join(f"{k}={v:.4g}" for k, v in (applied or {}).items()
                       if k != "warning" and isinstance(v, (int, float)))
    if values:
        notes.append(f"Background pre-fit replayed: {values}")
    if (applied or {}).get("warning"):
        notes.append(f"Background pre-fit: {applied['warning']}")
    return notes


def apply_prefits(setup: "ToolSetup", section: Dict,
                  q_full, intensity_full, q_fit=None, intensity_fit=None) -> list:
    """Run whatever pre-fit steps *setup*'s tool defines. Returns notes.

    The tools disagree about which data a pre-fit should see, and each is
    right. The Size Distribution background windows are chosen independently
    of the size-fit cursors, so they read the **full** curve, and so does the
    Simple Fits background pre-fit that the Invariant depends on. The WAXS
    peak presearch is re-centring the peaks that are about to be fitted, so
    it reads only the **fitted** range — scanning Q0 across data the fit will
    never see would move a peak towards a feature outside the window.
    """
    if q_fit is None:
        q_fit, intensity_fit = q_full, intensity_full
    if setup.tool == "sizes":
        return apply_sizes_prefits(setup.model, section, q_full, intensity_full)
    if setup.tool == "waxs_peakfit":
        return apply_waxs_presearch(setup.model, section, q_fit, intensity_fit)
    if setup.tool == "simple_fits":
        return apply_simple_bg_prefit(setup.model, setup.fixed_params,
                                      q_full, intensity_full)
    return []
