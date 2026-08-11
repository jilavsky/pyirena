"""pyirena.api.control.modeling — agent-drivable Modeling.

The Modeling counterpart of :mod:`pyirena.api.control.unified_fit`,
:mod:`pyirena.api.control.sizes` and :mod:`pyirena.api.control.simple_fits`.
Modeling builds a curve from **several populations** — size distributions with
a form factor and optional structure factor, Beaucage unified levels,
Guinier-Porod levels, diffraction peaks, mass fractals, surface fractals — so
this surface is mostly population management: add one, configure it, fit, look.

These functions operate on the *same*
:class:`~pyirena.api.control.session.Session` objects as the other control
surfaces, so the shared session-lifecycle tools are reused as-is.  The
**fitted Q range comes from the model** here (``config.q_min`` / ``q_max``,
set by :func:`set_modeling_q_range`), not from the session's Q-range tools,
because the Modeling engine crops the data itself.

One flat parameter namespace
----------------------------
A population's parameters live in several places on the dataclass — plain
attributes (``scale``, ``G``, ``position``) and nested dicts (``dist_params``,
``ff_params``, ``sf_params``).  Rather than make an agent learn that layout,
every tool here takes a **flat, dotted name**:

======================  ================================================
``scale``, ``contrast`` plain attribute of the population
``dist.mean_size``      entry in ``dist_params``
``ff.sld_core``         entry in ``ff_params``
``sf.eta``              entry in ``sf_params``
======================  ================================================

:func:`get_population_parameters` lists exactly the names that are active for
the population's current distribution, form factor and structure factor — the
same set the fitter would pack — so an agent can enumerate rather than guess.

Typical workflow
----------------
>>> import pyirena.api.control as ctrl
>>> sid = ctrl.open_dataset("/data/sample.h5")["session_id"]
>>> ctrl.select_modeling_model(sid)
>>> ctrl.list_population_types()
>>> ctrl.add_population(sid, "size_dist", label="matrix pores")
>>> ctrl.set_population_option(sid, 0, "form_factor", "sphere")
>>> ctrl.set_population_parameter(sid, 0, "dist.mean_size", 80.0)
>>> ctrl.set_population_parameter_fit(sid, 0, "dist.mean_size", True)
>>> ctrl.set_modeling_q_range(sid, 0.005, 0.35)
>>> ctrl.run_modeling_fit(sid)
>>> ctrl.get_modeling_fit_image(sid)
>>> ctrl.save_modeling_fit(sid)

All functions return plain dicts.  Errors are
``{"error": ..., "code": ..., "suggestion": ...}`` dicts rather than exceptions.
"""
from __future__ import annotations

import base64
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np

from pyirena.api._paths import PathSecurityError, resolve_safe
from pyirena.api.control.errors import make_error, no_fit, no_session

__all__ = [
    "select_modeling_model",
    "get_modeling_config",
    "list_population_types",
    "add_population",
    "remove_population",
    "list_populations",
    "set_population_enabled",
    "get_population_parameters",
    "set_population_parameter",
    "set_population_parameter_fit",
    "set_population_parameter_bounds",
    "set_population_option",
    "set_modeling_background",
    "set_modeling_q_range",
    "run_modeling_fit",
    "get_modeling_results",
    "get_modeling_fit_image",
    "save_modeling_fit",
]

#: Population types an agent can add, with a one-line description.
POPULATION_TYPES = {
    "size_dist": (
        "Polydisperse particles: a size distribution with a form factor and an "
        "optional structure factor. The general-purpose population."
    ),
    "unified_level": (
        "One Beaucage unified level (G, Rg, B, P) — a structural level without "
        "assuming a particle shape."
    ),
    "guinier_porod": (
        "Guinier-Porod piecewise model (Hammouda 2010): a Guinier region "
        "joined to a power law, for non-spherical or rod/sheet-like objects."
    ),
    "diffraction_peak": (
        "A Gaussian, Lorentzian or pseudo-Voigt peak — lamellar spacings, "
        "crystalline reflections."
    ),
    "mass_fractal": (
        "Mass fractal aggregate (Teixeira 1988): primary particles of radius R "
        "aggregated with fractal dimension Dv up to correlation length Ksi."
    ),
    "surface_fractal": (
        "Surface fractal: rough interfaces with fractal dimension Ds, with an "
        "optional Porod transition."
    ),
}

#: Non-numeric switches, per population type, that ``set_population_option``
#: accepts.  Changing one re-derives which parameters exist.
POPULATION_OPTIONS = {
    "size_dist": ["dist_type", "form_factor", "structure_factor",
                  "use_number_dist", "n_bins", "label"],
    "unified_level": ["correlations", "label"],
    "guinier_porod": ["correlations", "label"],
    "diffraction_peak": ["peak_type", "label"],
    "mass_fractal": ["label"],
    "surface_fractal": ["use_porod_transition", "label"],
}

_PEAK_TYPES = ["gaussian", "lorentzian", "voigt"]
_STRUCTURE_FACTORS = ["none", "interferences", "hard_sphere"]
_SF_ACTIVE = {"interferences": ["eta", "pack"],
              "hard_sphere": ["radius", "volume_fraction"]}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _require_modeling(session_id: str):
    """Return (session, None) when a Modeling config is ready, else (None, error)."""
    from pyirena.api.control.session import get_session  # noqa: PLC0415

    s = get_session(session_id)
    if s is None:
        return None, no_session(session_id)
    if s.model is None or s.model_name != "modeling":
        return None, make_error(
            f"Session '{session_id}' has no Modeling configuration.",
            suggestion="Call select_modeling_model() first.",
            code="NO_MODELING_MODEL",
        )
    return s, None


def _require_population(session_id: str, index: int):
    """Return (session, population, None) or (None, None, error)."""
    s, err = _require_modeling(session_id)
    if err:
        return None, None, err
    pops = s.model.populations
    if not isinstance(index, int) or index < 0 or index >= len(pops):
        return None, None, make_error(
            f"No population at index {index} (there are {len(pops)}).",
            suggestion="Call list_populations() to see the current indices.",
            code="BAD_POPULATION",
        )
    return s, pops[index], None


def _dist_param_names(pop) -> list:
    from pyirena.core.distributions import DIST_PARAM_NAMES  # noqa: PLC0415

    return list(DIST_PARAM_NAMES.get(pop.dist_type, []))


def _active_parameters(pop) -> list:
    """Flat parameter names active for this population, in fitter order.

    Mirrors ``ModelingEngine._pack_params``: only the parameters that the
    current distribution / form factor / structure factor actually use, so an
    agent never sets something the fit will ignore.
    """
    pop_type = getattr(pop, "pop_type", "size_dist")

    if pop_type == "unified_level":
        names = ["G", "Rg", "B", "P", "RgCO"]
        if getattr(pop, "correlations", False):
            names += ["ETA", "PACK"]
        return names

    if pop_type == "guinier_porod":
        names = ["G", "Rg1", "s1", "P", "Rg2", "s2", "RgCO"]
        if getattr(pop, "correlations", False):
            names += ["ETA", "PACK"]
        return names

    if pop_type == "diffraction_peak":
        names = ["position", "amplitude", "width"]
        if getattr(pop, "peak_type", "gaussian") == "voigt":
            names.append("eta_voigt")
        return names

    if pop_type == "mass_fractal":
        return ["Phi", "Radius", "Beta", "Dv", "Ksi", "Eta", "Contrast"]

    if pop_type == "surface_fractal":
        names = ["Surface", "Ds", "Ksi", "Contrast"]
        if getattr(pop, "use_porod_transition", False):
            names += ["Qc", "QcWidth"]
        return names

    # size_dist
    names = [f"dist.{n}" for n in _dist_param_names(pop)]
    names += [f"ff.{n}" for n in pop.ff_params]
    names += [f"sf.{n}" for n in _SF_ACTIVE.get(
        str(pop.structure_factor).lower(), [])]
    names += ["scale"]
    from pyirena.core.form_factors import CONTRAST_FREE_FORM_FACTORS  # noqa: PLC0415
    if pop.form_factor not in CONTRAST_FREE_FORM_FACTORS:
        names.append("contrast")
    return names


def _resolve(pop, name: str):
    """Map a flat name to (kind, key); kind is 'attr', 'dist', 'ff' or 'sf'."""
    if "." in name:
        prefix, key = name.split(".", 1)
        if prefix in ("dist", "ff", "sf"):
            return prefix, key
        return None, None
    return "attr", name


def _get_param(pop, name: str):
    """Return (value, fit_flag, (lo, hi)) for a flat name, or None if absent."""
    kind, key = _resolve(pop, name)
    if kind == "attr":
        if not hasattr(pop, key):
            return None
        return (
            float(getattr(pop, key)),
            bool(getattr(pop, f"fit_{key}", False)),
            tuple(getattr(pop, f"{key}_limits", (None, None))),
        )
    if kind is None:
        return None
    values = getattr(pop, f"{kind}_params", {})
    if key not in values:
        return None
    fits = getattr(pop, f"{kind}_params_fit", {})
    limits = getattr(pop, f"{kind}_params_limits", {})
    return (
        float(values[key]),
        bool(fits.get(key, False)),
        tuple(limits.get(key, (None, None))),
    )


def _set_param(pop, name: str, value=None, fit=None, limits=None) -> bool:
    """Set value / fit flag / bounds for a flat name.  False if it does not exist."""
    kind, key = _resolve(pop, name)
    if kind == "attr":
        if not hasattr(pop, key):
            return False
        if value is not None:
            setattr(pop, key, float(value))
        if fit is not None:
            setattr(pop, f"fit_{key}", bool(fit))
        if limits is not None:
            setattr(pop, f"{key}_limits", tuple(limits))
        return True
    if kind is None:
        return False
    values = getattr(pop, f"{kind}_params", None)
    if values is None or key not in values:
        return False
    if value is not None:
        values[key] = float(value)
    if fit is not None:
        getattr(pop, f"{kind}_params_fit")[key] = bool(fit)
    if limits is not None:
        getattr(pop, f"{kind}_params_limits")[key] = tuple(limits)
    return True


def _param_table(pop) -> list:
    rows = []
    for name in _active_parameters(pop):
        got = _get_param(pop, name)
        if got is None:
            continue
        value, fit, (lo, hi) = got
        rows.append({
            "name": name,
            "value": value,
            "fit": fit,
            "lo": (None if lo is None else float(lo)),
            "hi": (None if hi is None else float(hi)),
        })
    return rows


def _population_summary(index: int, pop) -> dict:
    out = {
        "index": index,
        "pop_type": getattr(pop, "pop_type", "size_dist"),
        "label": getattr(pop, "label", "") or "",
        "enabled": bool(getattr(pop, "enabled", True)),
        "n_free_parameters": sum(1 for p in _param_table(pop) if p["fit"]),
    }
    if out["pop_type"] == "size_dist":
        out["dist_type"] = pop.dist_type
        out["form_factor"] = pop.form_factor
        out["structure_factor"] = pop.structure_factor
    elif out["pop_type"] == "diffraction_peak":
        out["peak_type"] = pop.peak_type
    return out


def _bad_param(name: str, pop) -> dict:
    return make_error(
        f"Parameter '{name}' is not active for this "
        f"{getattr(pop, 'pop_type', 'population')} population.",
        suggestion=(
            "Call get_population_parameters() for the active names. Nested "
            "parameters use a prefix: 'dist.mean_size', 'ff.sld_core', 'sf.eta'."
        ),
        code="BAD_PARAM",
    )


def _f(value):
    """Float, or None for missing/non-finite — JSON has no NaN."""
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def _render_image(fig, session_id: str, tag: str, dpi: int) -> str:
    import matplotlib.pyplot as plt  # noqa: PLC0415

    tmp = Path(tempfile.gettempdir()) / "pyirena-ctrl"
    tmp.mkdir(parents=True, exist_ok=True)
    out = tmp / f"modeling_{tag}_{session_id}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(out.read_bytes()).decode("ascii")


# ---------------------------------------------------------------------------
# Model lifecycle
# ---------------------------------------------------------------------------

def select_modeling_model(session_id: str) -> dict:
    """Start a Modeling configuration for the session (no populations yet).

    The Q range defaults to the data's full range; narrow it with
    :func:`set_modeling_q_range`.  Add populations with :func:`add_population`.
    Calling this again discards the current configuration and any fit.
    """
    from pyirena.api.control.session import get_session  # noqa: PLC0415
    from pyirena.core.modeling import ModelingConfig  # noqa: PLC0415

    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    config = ModelingConfig(
        populations=[],
        q_min=float(np.min(s.q)),
        q_max=float(np.max(s.q)),
    )
    if s.is_slit_smeared and s.slit_length > 0:
        config.use_slit_smearing = True
        config.slit_length = float(s.slit_length)

    s.model = config
    s.model_name = "modeling"
    s.last_fit_result = None

    return {"ok": True, "model": "modeling", "config": _config_summary(config)}


def _config_summary(config) -> dict:
    return {
        "q_min": float(config.q_min),
        "q_max": float(config.q_max),
        "background": float(config.background),
        "fit_background": bool(config.fit_background),
        "fit_method": str(config.fit_method),
        "use_slit_smearing": bool(config.use_slit_smearing),
        "slit_length": float(config.slit_length),
        "n_populations": len(config.populations),
        "n_enabled": sum(1 for p in config.populations if getattr(p, "enabled", True)),
    }


def get_modeling_config(session_id: str) -> dict:
    """Return the global Modeling settings and a summary of every population."""
    s, err = _require_modeling(session_id)
    if err:
        return err
    return {
        "ok": True,
        "config": _config_summary(s.model),
        "populations": [_population_summary(i, p)
                        for i, p in enumerate(s.model.populations)],
    }


def list_population_types() -> dict:
    """Describe the population types, their options and their parameters.

    No session needed.  Call before :func:`add_population` to choose a type,
    and to learn which non-numeric options (distribution, form factor,
    structure factor, peak type) that type accepts.
    """
    from pyirena.core.distributions import DIST_PARAM_NAMES  # noqa: PLC0415
    from pyirena.core.form_factors import (  # noqa: PLC0415
        FORM_FACTOR_PARAMS,
        list_form_factors,
    )

    types = []
    for name, description in POPULATION_TYPES.items():
        entry = {
            "pop_type": name,
            "description": description,
            "options": POPULATION_OPTIONS.get(name, []),
        }
        if name == "size_dist":
            entry["distributions"] = {
                d: list(params) for d, params in DIST_PARAM_NAMES.items()
            }
            entry["form_factors"] = {
                ff: list(FORM_FACTOR_PARAMS[ff]) for ff in list_form_factors()
            }
            entry["structure_factors"] = {
                sf: list(_SF_ACTIVE.get(sf, [])) for sf in _STRUCTURE_FACTORS
            }
        elif name == "diffraction_peak":
            entry["peak_types"] = list(_PEAK_TYPES)
        types.append(entry)
    return {"ok": True, "population_types": types}


# ---------------------------------------------------------------------------
# Populations
# ---------------------------------------------------------------------------

_POP_CLASSES = {
    "size_dist": "SizeDistPopulation",
    "unified_level": "UnifiedLevelPopulation",
    "guinier_porod": "GuinierPorodPopulation",
    "diffraction_peak": "DiffractionPeakPopulation",
    "mass_fractal": "MassFractalPopulation",
    "surface_fractal": "SurfaceFractalPopulation",
}


def add_population(session_id: str, pop_type: str = "size_dist",
                   label: str = "") -> dict:
    """Add a population to the model and return its index and parameters.

    Parameters
    ----------
    pop_type : str
        One of the types from :func:`list_population_types`.
    label : str
        Optional name shown in reports and plots ("matrix pores").

    Notes
    -----
    The population starts from the core defaults, which are deliberately
    generic — set the parameters that matter for your sample before fitting.
    """
    s, err = _require_modeling(session_id)
    if err:
        return err

    if pop_type not in _POP_CLASSES:
        return make_error(
            f"Unknown population type '{pop_type}'.",
            suggestion=f"Choose one of {list(_POP_CLASSES)}.",
            code="BAD_POPULATION_TYPE",
        )

    import pyirena.core.modeling as _modeling  # noqa: PLC0415

    pop = getattr(_modeling, _POP_CLASSES[pop_type])()
    if label:
        pop.label = str(label)

    # Form-factor parameters are empty by default; fill in the ones the
    # default form factor needs so the population is immediately evaluable.
    if pop_type == "size_dist":
        _sync_ff_params(pop)

    s.model.populations.append(pop)
    index = len(s.model.populations) - 1
    s.last_fit_result = None
    return {
        "ok": True,
        "index": index,
        "population": _population_summary(index, pop),
        "parameters": _param_table(pop),
    }


def _sync_ff_params(pop) -> None:
    """Make ``ff_params`` match the population's current form factor.

    Parameters the new form factor does not use are dropped; new ones arrive
    with their registry defaults and bounds.  Values the agent already set for
    a parameter shared by both form factors are kept.
    """
    from pyirena.core.form_factors import (  # noqa: PLC0415
        form_factor_param_defaults,
        form_factor_params,
    )

    wanted = form_factor_params(pop.form_factor)
    for key in list(pop.ff_params):
        if key not in wanted:
            pop.ff_params.pop(key, None)
            pop.ff_params_fit.pop(key, None)
            pop.ff_params_limits.pop(key, None)
    for key in wanted:
        if key not in pop.ff_params:
            value, lo, hi = form_factor_param_defaults(key)
            pop.ff_params[key] = value
            pop.ff_params_fit.setdefault(key, False)
            pop.ff_params_limits[key] = (lo, hi)


def remove_population(session_id: str, index: int) -> dict:
    """Remove the population at *index*; later populations shift down."""
    s, pop, err = _require_population(session_id, index)
    if err:
        return err
    s.model.populations.pop(index)
    s.last_fit_result = None
    return {"ok": True, "removed": index,
            "populations": [_population_summary(i, p)
                            for i, p in enumerate(s.model.populations)]}


def list_populations(session_id: str) -> dict:
    """List every population with its type, label, enabled state and free-parameter count."""
    s, err = _require_modeling(session_id)
    if err:
        return err
    return {"ok": True,
            "populations": [_population_summary(i, p)
                            for i, p in enumerate(s.model.populations)]}


def set_population_enabled(session_id: str, index: int, enabled: bool = True) -> dict:
    """Include or exclude a population from the model without deleting it.

    Disabling is the way to test what a population contributes: fit with it
    off, compare chi-squared, turn it back on.
    """
    s, pop, err = _require_population(session_id, index)
    if err:
        return err
    pop.enabled = bool(enabled)
    return {"ok": True, "population": _population_summary(index, pop)}


# ---------------------------------------------------------------------------
# Population parameters
# ---------------------------------------------------------------------------

def get_population_parameters(session_id: str, index: int) -> dict:
    """List the parameters active for this population, with value, fit flag and bounds.

    Only the parameters the current distribution / form factor / structure
    factor actually use are returned — the same set the fitter packs.
    """
    s, pop, err = _require_population(session_id, index)
    if err:
        return err
    return {
        "ok": True,
        "index": index,
        "population": _population_summary(index, pop),
        "parameters": _param_table(pop),
    }


def set_population_parameter(session_id: str, index: int,
                             name: str, value: float) -> dict:
    """Set one parameter's value.  Nested names are dotted: ``dist.mean_size``."""
    s, pop, err = _require_population(session_id, index)
    if err:
        return err
    if name not in _active_parameters(pop):
        return _bad_param(name, pop)
    try:
        ok = _set_param(pop, name, value=float(value))
    except (TypeError, ValueError):
        return make_error(
            f"Value '{value}' for '{name}' is not a number.",
            suggestion="Pass a numeric value.",
            code="BAD_VALUE",
        )
    if not ok:
        return _bad_param(name, pop)
    got = _get_param(pop, name)
    return {"ok": True, "index": index, "name": name, "value": got[0]}


def set_population_parameter_fit(session_id: str, index: int,
                                 name: str, fit: bool = True) -> dict:
    """Choose whether one parameter is fitted or held at its current value.

    Modeling has many parameters and few constraints, so fit a handful at a
    time — every population's defaults start mostly held.
    """
    s, pop, err = _require_population(session_id, index)
    if err:
        return err
    if name not in _active_parameters(pop):
        return _bad_param(name, pop)
    if not _set_param(pop, name, fit=bool(fit)):
        return _bad_param(name, pop)
    return {"ok": True, "index": index, "name": name, "fit": bool(fit)}


def set_population_parameter_bounds(
    session_id: str, index: int, name: str,
    lo: Optional[float] = None, hi: Optional[float] = None,
) -> dict:
    """Set the fitting bounds of one parameter.

    Both sides are required by the global fit method (differential evolution
    needs finite bounds), so passing ``None`` keeps the existing side rather
    than removing the bound.  A value outside the new range is clamped in.
    """
    s, pop, err = _require_population(session_id, index)
    if err:
        return err
    if name not in _active_parameters(pop):
        return _bad_param(name, pop)

    current = _get_param(pop, name)
    cur_lo, cur_hi = current[2]
    new_lo = cur_lo if lo is None else float(lo)
    new_hi = cur_hi if hi is None else float(hi)
    if new_lo is not None and new_hi is not None and new_lo >= new_hi:
        return make_error(
            f"Lower bound {new_lo} is not below upper bound {new_hi}.",
            suggestion="Pass lo < hi.",
            code="BAD_BOUNDS",
        )

    _set_param(pop, name, limits=(new_lo, new_hi))

    value = current[0]
    if new_lo is not None and value < new_lo:
        _set_param(pop, name, value=new_lo)
    elif new_hi is not None and value > new_hi:
        _set_param(pop, name, value=new_hi)

    got = _get_param(pop, name)
    return {"ok": True, "index": index, "name": name,
            "lo": _f(new_lo), "hi": _f(new_hi), "value": got[0]}


def set_population_option(session_id: str, index: int,
                          option: str, value: str) -> dict:
    """Set a non-numeric switch: distribution, form factor, structure factor, …

    Which options a population accepts is reported by
    :func:`list_population_types` and by ``get_population_parameters``.
    Changing one re-derives the active parameter list — switching form factor
    to ``cs_sphere_by_core``, for instance, adds the SLD and shell-thickness
    parameters with sensible defaults and drops any that no longer apply.

    Booleans accept true/false (or the strings "true"/"false"); ``n_bins``
    takes an integer; ``label`` any text.
    """
    s, pop, err = _require_population(session_id, index)
    if err:
        return err

    pop_type = getattr(pop, "pop_type", "size_dist")
    allowed = POPULATION_OPTIONS.get(pop_type, [])
    if option not in allowed:
        return make_error(
            f"'{option}' is not an option of a {pop_type} population.",
            suggestion=f"This type accepts {allowed}.",
            code="BAD_OPTION",
        )

    def _as_bool(v):
        if isinstance(v, bool):
            return v
        return str(v).strip().lower() in ("1", "true", "yes", "on")

    if option == "label":
        pop.label = str(value)
    elif option in ("correlations", "use_number_dist", "use_porod_transition"):
        setattr(pop, option, _as_bool(value))
    elif option == "n_bins":
        try:
            pop.n_bins = max(2, int(value))
        except (TypeError, ValueError):
            return make_error(
                f"n_bins '{value}' is not an integer.",
                suggestion="Pass an integer, e.g. 200.",
                code="BAD_VALUE",
            )
    elif option == "dist_type":
        from pyirena.core.distributions import DIST_PARAM_NAMES  # noqa: PLC0415
        if value not in DIST_PARAM_NAMES:
            return make_error(
                f"Unknown distribution '{value}'.",
                suggestion=f"Choose one of {list(DIST_PARAM_NAMES)}.",
                code="BAD_OPTION_VALUE",
            )
        pop.dist_type = str(value)
        _sync_dist_params(pop)
    elif option == "form_factor":
        from pyirena.core.form_factors import list_form_factors  # noqa: PLC0415
        if value not in list_form_factors():
            return make_error(
                f"Unknown form factor '{value}'.",
                suggestion=f"Choose one of {list_form_factors()}.",
                code="BAD_OPTION_VALUE",
            )
        pop.form_factor = str(value)
        _sync_ff_params(pop)
    elif option == "structure_factor":
        if value not in _STRUCTURE_FACTORS:
            return make_error(
                f"Unknown structure factor '{value}'.",
                suggestion=f"Choose one of {_STRUCTURE_FACTORS}.",
                code="BAD_OPTION_VALUE",
            )
        pop.structure_factor = str(value)
    elif option == "peak_type":
        if value not in _PEAK_TYPES:
            return make_error(
                f"Unknown peak type '{value}'.",
                suggestion=f"Choose one of {_PEAK_TYPES}.",
                code="BAD_OPTION_VALUE",
            )
        pop.peak_type = str(value)

    s.last_fit_result = None
    return {
        "ok": True,
        "index": index,
        "population": _population_summary(index, pop),
        "parameters": _param_table(pop),
    }


def _sync_dist_params(pop) -> None:
    """Give the new distribution its parameters, keeping shared ones."""
    from pyirena.core.distributions import DIST_PARAM_NAMES  # noqa: PLC0415

    defaults = {
        "min_size": (10.0, 0.0, 1e6),
        "mean_size": (100.0, 1.0, 1e6),
        "sdeviation": (0.3, 0.01, 5.0),
        "width": (20.0, 0.01, 1e6),
        "location": (100.0, 1.0, 1e6),
        "parameter": (2.5, 2.01, 2.99),
    }
    for key in DIST_PARAM_NAMES.get(pop.dist_type, []):
        if key not in pop.dist_params:
            value, lo, hi = defaults.get(key, (1.0, 0.0, 1e10))
            pop.dist_params[key] = value
            pop.dist_params_fit.setdefault(key, False)
            pop.dist_params_limits[key] = (lo, hi)


# ---------------------------------------------------------------------------
# Global settings
# ---------------------------------------------------------------------------

def set_modeling_background(session_id: str, value: Optional[float] = None,
                            fit: Optional[bool] = None) -> dict:
    """Set the flat background level and/or whether it is fitted."""
    s, err = _require_modeling(session_id)
    if err:
        return err
    if value is not None:
        try:
            s.model.background = float(value)
        except (TypeError, ValueError):
            return make_error(
                f"Background '{value}' is not a number.",
                suggestion="Pass a numeric value.",
                code="BAD_VALUE",
            )
    if fit is not None:
        s.model.fit_background = bool(fit)
    return {"ok": True, "background": float(s.model.background),
            "fit_background": bool(s.model.fit_background)}


def set_modeling_q_range(session_id: str, q_min: Optional[float] = None,
                         q_max: Optional[float] = None) -> dict:
    """Set the Q range the model is fitted over.

    Modeling crops the data itself, so this — not the shared
    ``set_fit_q_range`` — is what limits a Modeling fit.  Pass ``None`` to
    leave one end unchanged.
    """
    s, err = _require_modeling(session_id)
    if err:
        return err

    new_min = s.model.q_min if q_min is None else float(q_min)
    new_max = s.model.q_max if q_max is None else float(q_max)
    if new_min >= new_max:
        return make_error(
            f"q_min {new_min} is not below q_max {new_max}.",
            suggestion="Pass q_min < q_max.",
            code="BAD_RANGE",
        )
    n_points = int(np.sum((s.q >= new_min) & (s.q <= new_max)))
    if n_points == 0:
        return make_error(
            f"No data points between Q = {new_min} and {new_max}.",
            suggestion=(
                f"The data covers {float(np.min(s.q)):.4g} to "
                f"{float(np.max(s.q)):.4g} Å⁻¹."
            ),
            code="EMPTY_RANGE",
        )

    s.model.q_min, s.model.q_max = new_min, new_max
    return {"ok": True, "q_min": new_min, "q_max": new_max, "n_points": n_points}


# ---------------------------------------------------------------------------
# Fit execution
# ---------------------------------------------------------------------------

def run_modeling_fit(session_id: str, fit_method: str = "local") -> dict:
    """Fit the enabled populations to the data over the model's Q range.

    Parameters
    ----------
    fit_method : str
        ``"local"`` (default) is a least-squares refinement from the current
        values — fast, and what you want once the model is roughly right.
        ``"global"`` runs differential evolution first to find the basin, then
        polishes locally: slower, but the right choice for core-shell models
        whose chi-squared surface has many minima.  Global needs finite bounds
        on every fitted parameter.

    Returns
    -------
    dict with success, chi_squared, reduced_chi_squared, dof, n_data, and the
    fitted value of every free parameter grouped by population.
    """
    s, err = _require_modeling(session_id)
    if err:
        return err

    config = s.model
    enabled = [p for p in config.populations if getattr(p, "enabled", True)]
    if not enabled:
        return make_error(
            "No enabled populations to fit.",
            suggestion="Call add_population(), or set_population_enabled(index, True).",
            code="NO_POPULATIONS",
        )

    if fit_method not in ("local", "global"):
        return make_error(
            f"Unknown fit method '{fit_method}'.",
            suggestion="Use 'local' or 'global'.",
            code="BAD_METHOD",
        )
    config.fit_method = fit_method

    n_free = sum(1 for p in enabled for row in _param_table(p) if row["fit"])
    if config.fit_background:
        n_free += 1
    if n_free == 0:
        return make_error(
            "No free parameters — every parameter is held.",
            suggestion=(
                "Call set_population_parameter_fit(index, name, True) for the "
                "parameters you want fitted."
            ),
            code="ALL_FIXED",
        )

    error = s.error if s.error is not None else 0.05 * np.abs(s.intensity)

    from pyirena.core.modeling import ModelingEngine  # noqa: PLC0415

    try:
        result = ModelingEngine().fit(config, s.q, s.intensity, error)
    except ValueError as exc:
        return make_error(
            f"Fit failed: {exc}",
            suggestion="Check the Q range and the population parameters.",
            code="FIT_FAILED",
        )
    except Exception as exc:  # pragma: no cover - defensive
        return make_error(
            f"Fit failed with exception: {exc}",
            suggestion="Check the population setup; try fit_method='local'.",
            code="FIT_EXCEPTION",
        )

    s.last_fit_result = result
    # The engine returns the fitted config; adopt it so follow-up calls see the
    # optimised values (parity with the GUI, where the panel updates in place).
    s.model = result.config

    return {
        "success": True,
        "fit_method": fit_method,
        "chi_squared": _f(result.chi_squared),
        "reduced_chi_squared": _f(result.reduced_chi_squared),
        "dof": int(result.dof),
        "n_data": int(len(result.model_q)),
        "n_free_parameters": n_free,
        "background": _f(result.config.background),
        "populations": _fitted_population_rows(result),
    }


def _fitted_population_rows(result) -> list:
    """Per-population parameter values after a fit, plus derived quantities."""
    rows = []
    config = result.config
    for slot, pop_index in enumerate(result.pop_indices):
        pop = config.populations[pop_index]
        derived = result.derived[slot] if slot < len(result.derived) else {}
        rows.append({
            "index": pop_index,
            "pop_type": getattr(pop, "pop_type", "size_dist"),
            "label": getattr(pop, "label", "") or "",
            "parameters": [
                {"name": row["name"], "value": row["value"], "fit": row["fit"]}
                for row in _param_table(pop)
            ],
            "derived": {k: _f(v) for k, v in (derived or {}).items()},
        })
    return rows


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

def get_modeling_results(session_id: str) -> dict:
    """Return the last Modeling fit: quality, background and every population."""
    s, err = _require_modeling(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    result = s.last_fit_result
    return {
        "ok": True,
        "chi_squared": _f(result.chi_squared),
        "reduced_chi_squared": _f(result.reduced_chi_squared),
        "dof": int(result.dof),
        "background": _f(result.config.background),
        "q_min": _f(result.config.q_min),
        "q_max": _f(result.config.q_max),
        "populations": _fitted_population_rows(result),
    }


def get_modeling_fit_image(
    session_id: str, width: int = 1024, height: int = 800, dpi: int = 120,
) -> dict:
    """Render the fit as a PNG: log-log data, total model and each population.

    Seeing the per-population curves is the point — it shows which population
    carries which part of the curve, and whether one has collapsed to nothing.
    """
    s, err = _require_modeling(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt  # noqa: PLC0415

    result = s.last_fit_result
    q = np.asarray(result.model_q, dtype=float)
    total = np.asarray(result.model_I, dtype=float)
    mask = (s.q >= result.config.q_min) & (s.q <= result.config.q_max)
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
    ax1.plot(q, total, "-", linewidth=2, color="red", label="Total model")

    for slot, pop_index in enumerate(result.pop_indices):
        if slot >= len(result.pop_model_I):
            break
        pop = result.config.populations[pop_index]
        label = getattr(pop, "label", "") or f"P{pop_index + 1}"
        ax1.plot(q, np.asarray(result.pop_model_I[slot], dtype=float),
                 "--", linewidth=1.2, alpha=0.8,
                 label=f"{label} ({getattr(pop, 'pop_type', '')})")

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_ylabel("I  (cm⁻¹)", fontsize=10)
    chi2r = _f(result.reduced_chi_squared)
    title = "Modeling"
    if chi2r is not None:
        title += f" — reduced χ² = {chi2r:.3g}"
    ax1.set_title(title, fontsize=11, fontweight="bold")
    ax1.legend(fontsize=8, loc="best")
    ax1.grid(True, which="both", alpha=0.25)

    sigma = err_obs if err_obs is not None else np.abs(I_obs) * 0.05
    residuals = (I_obs - total) / np.maximum(sigma, 1e-30)
    ax2.axhline(0.0, color="black", linewidth=1, linestyle="--")
    ax2.plot(q, residuals, "o", markersize=3, color="steelblue")
    ax2.set_xscale("log")
    ax2.set_xlabel("Q  (Å⁻¹)", fontsize=10)
    ax2.set_ylabel("Residuals", fontsize=10)
    ax2.grid(True, which="both", alpha=0.25)

    fig.tight_layout()
    return {"ok": True, "image_base64": _render_image(fig, session_id, "fit", dpi)}


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def save_modeling_fit(session_id: str, output_path: Optional[str] = None) -> dict:
    """Save the Modeling fit to NXcanSAS HDF5 under ``entry/modeling_results``.

    Parameters
    ----------
    output_path : str, optional
        Where to save.  Defaults to the original file (in-place).  Saving to a
        different path preserves the original.
    """
    s, err = _require_modeling(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    from pyirena.io.nxcansas_modeling import save_modeling_results  # noqa: PLC0415

    try:
        src = resolve_safe(s.file_path, must_exist=False)
        target = resolve_safe(output_path, must_exist=False) if output_path else src
    except PathSecurityError as exc:
        return make_error(
            str(exc),
            suggestion="Save to a path inside PYIRENA_DATA_ROOT.",
            code="PATH_NOT_ALLOWED",
        )

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

    try:
        save_modeling_results(target, s.last_fit_result)
    except Exception as exc:
        return make_error(
            f"Could not save results: {exc}",
            suggestion="Check the file is writable and not open elsewhere.",
            code="SAVE_ERROR",
        )

    return {"ok": True, "saved_to": str(target), "group": "entry/modeling_results"}
