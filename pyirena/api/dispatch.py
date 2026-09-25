"""Fixed dispatcher over ``pyirena.api`` — one tool surface, many transports.

``pyirena-mcp`` used to register one MCP tool per ``pyirena.api.control``
function (~90 ``pyirena_ctrl_*`` tools). Combined with whatever else is
active in a client's session, that can exceed a hard provider-side cap on
the number of tools in a single request (observed via the ANL Argo gateway
proxy: 128 tools, ``array_above_max_length``). It is also a standing
per-turn token cost, since every enabled tool's schema is sent on every
turn regardless of use.

This module collapses that surface into four operations — list categories,
list tools in a category, describe one tool's schema, and call it by name —
built directly on top of ``pyirena.api.control.schemas.TOOL_SCHEMA_BY_NAME``,
which already has an Anthropic-style schema for every control function.
Session-lifecycle functions (``open_dataset``, ``open_dataset_from_data``,
``list_open_sessions``, ``close_session``, ``get_session_summary``) are
excluded here and stay top-level tools instead, since nearly every workflow
starts there.

The dispatcher is not control-only: ``_SOURCES`` below lists every schema
registry it serves. ``pyirena.api.calculators`` joins as the stateless
"calculators" category and ``pyirena.api.data_ops`` as the "data" category,
so support calculators and data operations (scattering contrast, averaging,
subtraction, merging) reach agents without adding a single top-level tool.

Profiles
--------
Not every transport can serve every tool. A stdio MCP server shares a
filesystem with its client, so paths are meaningful; a ZMQ service on
another machine has no shared filesystem, and a tool that takes
``file_path`` or returns ``image_path`` is worse than useless there — it
looks like it worked and points at a file the caller cannot read. Each tool
therefore carries two flags, ``touches_files`` and ``returns_image``, and
every public function takes a ``profile`` selecting which of those are
visible. ``PROFILE_ALL`` (the default, what MCP uses) hides nothing;
``PROFILE_JSON_ONLY`` hides both.

Design docs: AIDA's ``planning/mcp_tool_scaling.md`` (Tier 2) and
``planning/zmq-service/`` (profiles).

This module lives under ``pyirena/api/`` rather than ``pyirena/mcp/``
because two transports now use it; ``pyirena/mcp/dispatch.py`` remains as a
re-export shim. It has no dependency on the ``mcp`` package on purpose, so
it can be imported and tested without it installed. ``pyirena/mcp/server.py``
is the only place that turns a result into MCP content (e.g. image blocks).
"""
from __future__ import annotations

import difflib
import inspect
from dataclasses import dataclass
from typing import Any, Union

from pyirena.api import calculators as _calc
from pyirena.api import control as _ctrl
from pyirena.api import data_ops as _data_ops
from pyirena.api.calculator_schemas import CALCULATOR_SCHEMA_BY_NAME
from pyirena.api.control.schemas import TOOL_SCHEMA_BY_NAME
from pyirena.api.data_op_schemas import DATA_OP_SCHEMA_BY_NAME

SESSION_LIFECYCLE_NAMES = frozenset(
    {
        "open_dataset",
        "open_dataset_from_data",
        "list_open_sessions",
        "close_session",
        "get_session_summary",
    }
)

_CATEGORY_BY_MODULE = {
    "unified_fit": "unified",
    "sizes": "sizes",
    "simple_fits": "simple",
    "modeling": "modeling",
    "waxs_peakfit": "waxs",
    "carbon_fit": "carbon",
    "export": "results",
}

_CATEGORY_BLURBS = {
    "unified": "Unified Fit — Beaucage multi-level model for a whole multi-level SAXS/USAXS curve.",
    "sizes": "Size Distribution — invert a dilute single population to a size histogram P(r).",
    "simple": "Simple Fits — one analytical model (Guinier, Porod, Sphere, ...) over a Q sub-range.",
    "modeling": "Modeling — multi-population forward model with specific form/structure factors.",
    "waxs": "WAXS Peak Fit — peak position, width and area for wide-angle patterns.",
    "carbon": (
        "Carbon model — disordered carbons fitted across the whole SAXS+WAXS "
        "range at once: grain surface, micropores and turbostratic stacking, "
        "with density, contrast and pore/wall widths derived from the fit."
    ),
    "results": (
        "Results — export a finished fit as one JSON document, whichever "
        "of the six fitting tools produced it: parameters, quality, "
        "derived quantities and the model config needed to replay it."
    ),
    "calculators": (
        "Calculators — stateless support calculations that need no dataset "
        "and open no session (scattering contrast, SLDs, element lookup)."
    ),
    "data": (
        "Data operations — average, subtract, divide, scale, trim, rebin and "
        "merge datasets. The only tools that WRITE new data files."
    ),
}


# ---------------------------------------------------------------------------
# Profiles — which tools a given transport can honestly serve
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Profile:
    """What a transport is willing to expose.

    Attributes
    ----------
    name : str
        Identifier reported back to callers (``server_info``, error messages).
    allow_files : bool
        Expose tools that take or produce a path on the *server's* filesystem.
    allow_images : bool
        Expose tools that render a PNG. These return ``image_base64`` as well
        as ``image_path``, so a remote transport may enable them once it is
        prepared to drop the path.
    """

    name: str
    allow_files: bool = True
    allow_images: bool = True


PROFILE_ALL = Profile("all", allow_files=True, allow_images=True)
PROFILE_JSON_ONLY = Profile("json_only", allow_files=False, allow_images=False)
PROFILE_JSON_IMAGES = Profile("json_images", allow_files=False, allow_images=True)

_PROFILES_BY_NAME = {p.name: p for p in (PROFILE_ALL, PROFILE_JSON_ONLY, PROFILE_JSON_IMAGES)}

ProfileArg = Union[str, Profile, None]


def resolve_profile(profile: ProfileArg = None) -> Profile:
    """Coerce a profile name (or None, meaning ``all``) to a Profile."""
    if profile is None:
        return PROFILE_ALL
    if isinstance(profile, Profile):
        return profile
    try:
        return _PROFILES_BY_NAME[profile]
    except KeyError:
        raise ValueError(
            f"Unknown dispatcher profile '{profile}'. "
            f"Valid profiles: {', '.join(sorted(_PROFILES_BY_NAME))}."
        ) from None


# Tools whose arguments or results name a path on the machine running the
# service. Classified per source module rather than per schema so that a new
# tool is covered the day it is added: every data operation reads and writes
# data files, and in the control API the convention is that `save_*` writes
# one and `open_dataset` reads one. `pyirena/tests/api/test_dispatch_profiles.py`
# fails if a schema grows a path-like argument that these rules miss.
_FILE_TOOL_CATEGORIES = frozenset({"data"})


def _touches_files(name: str, category: str) -> bool:
    if category in _FILE_TOOL_CATEGORIES:
        return True
    return name.startswith("save_") or name in ("open_dataset", "match_merge_files")


def _returns_image(name: str, category: str) -> bool:
    return name.endswith("_image")


def _control_category_for(name: str) -> str:
    """Category of a control tool, from the api.control submodule owning it."""
    module_name = inspect.getmodule(getattr(_ctrl, name)).__name__.rsplit(".", 1)[-1]
    return _CATEGORY_BY_MODULE[module_name]


# Every schema registry the dispatcher serves, as
# (schema_by_name, owning_module, category_resolver). The owning module is
# stored per entry rather than the callable itself: functions are looked up
# by name at call time (see call_tool), same as every other pyirena_ctrl_*
# wrapper in server.py, so monkeypatching <module>.<name> (as the test suite
# does) is honoured.
_SOURCES: tuple[tuple[dict[str, dict], Any, Any], ...] = (
    (TOOL_SCHEMA_BY_NAME, _ctrl, _control_category_for),
    (CALCULATOR_SCHEMA_BY_NAME, _calc, lambda _name: "calculators"),
    (DATA_OP_SCHEMA_BY_NAME, _data_ops, lambda _name: "data"),
)


def _registry() -> dict[str, dict[str, Any]]:
    registry: dict[str, dict[str, Any]] = {}
    for schemas, module, category_for in _SOURCES:
        for name, schema in schemas.items():
            if name in SESSION_LIFECYCLE_NAMES:
                continue
            if name in registry:
                raise RuntimeError(
                    f"Duplicate dispatcher tool name '{name}' across schema "
                    f"registries; names must be unique."
                )
            category = category_for(name)
            registry[name] = {
                "category": category,
                "schema": schema,
                "module": module,
                "touches_files": _touches_files(name, category),
                "returns_image": _returns_image(name, category),
            }
    return registry


_REGISTRY = _registry()


def _visible(entry: dict[str, Any], profile: Profile) -> bool:
    if entry["touches_files"] and not profile.allow_files:
        return False
    if entry["returns_image"] and not profile.allow_images:
        return False
    return True


def _hidden_error(name: str, entry: dict[str, Any], profile: Profile) -> dict:
    reason = (
        "reads or writes a file on the server's filesystem"
        if entry["touches_files"] and not profile.allow_files
        else "returns a rendered image"
    )
    return {
        "error": f"Tool '{name}' {reason} and is not available in profile '{profile.name}'.",
        "code": "NOT_AVAILABLE_IN_PROFILE",
        "suggestion": (
            "Pass the data in the request instead (open_dataset_from_data) and "
            "read results with export_results."
            if entry["touches_files"]
            else "Ask the service operator to enable images, or use the numeric results instead."
        ),
    }


def _first_sentence(description: str) -> str:
    text = description.strip().split("\n", 1)[0]
    for sep in (". ", ".\n"):
        if sep in text:
            return text.split(sep, 1)[0] + "."
    return text


def _unknown_tool_error(name: str, profile: Profile) -> dict:
    if name in SESSION_LIFECYCLE_NAMES:
        return {
            "error": f"'{name}' is a session-lifecycle tool, not a dispatched one.",
            "code": "UNKNOWN_TOOL",
            "suggestion": f"Call the top-level 'pyirena_ctrl_{name}' tool directly.",
        }
    visible = [n for n, e in _REGISTRY.items() if _visible(e, profile)]
    candidates = difflib.get_close_matches(name, visible, n=3)
    suggestion = (
        f"Did you mean: {', '.join(candidates)}?"
        if candidates
        else "Call pyirena_list_categories() then pyirena_list_tools(category) to see valid names."
    )
    return {
        "error": f"Unknown tool '{name}'.",
        "code": "UNKNOWN_TOOL",
        "suggestion": suggestion,
    }


def tool_flags(name: str) -> dict:
    """Return ``{category, touches_files, returns_image}`` for *name*, or an error dict."""
    entry = _REGISTRY.get(name)
    if entry is None:
        return _unknown_tool_error(name, PROFILE_ALL)
    return {
        "name": name,
        "category": entry["category"],
        "touches_files": entry["touches_files"],
        "returns_image": entry["returns_image"],
    }


def list_categories(profile: ProfileArg = None) -> dict:
    """Return every dispatch category with its tool count in *profile*."""
    prof = resolve_profile(profile)
    counts: dict[str, int] = {}
    for entry in _REGISTRY.values():
        if not _visible(entry, prof):
            continue
        counts[entry["category"]] = counts.get(entry["category"], 0) + 1
    categories = [
        {"name": name, "description": _CATEGORY_BLURBS[name], "tool_count": counts.get(name, 0)}
        for name in _CATEGORY_BLURBS
        if counts.get(name, 0) > 0
    ]
    return {"categories": categories, "profile": prof.name}


def list_tools(category: str, profile: ProfileArg = None) -> dict:
    """Return the name and one-line summary of every tool in *category*."""
    prof = resolve_profile(profile)
    if category not in _CATEGORY_BLURBS:
        return {
            "error": f"Unknown category '{category}'.",
            "code": "UNKNOWN_CATEGORY",
            "suggestion": f"Valid categories: {', '.join(_CATEGORY_BLURBS)}.",
        }
    tools = sorted(
        (
            {"name": name, "summary": _first_sentence(entry["schema"]["description"])}
            for name, entry in _REGISTRY.items()
            if entry["category"] == category and _visible(entry, prof)
        ),
        key=lambda t: t["name"],
    )
    return {"category": category, "tools": tools, "profile": prof.name}


def describe_tool(name: str, profile: ProfileArg = None) -> dict:
    """Return the full Anthropic-style tool schema for *name*, plus its category."""
    prof = resolve_profile(profile)
    entry = _REGISTRY.get(name)
    if entry is None:
        return _unknown_tool_error(name, prof)
    if not _visible(entry, prof):
        return _hidden_error(name, entry, prof)
    return {**entry["schema"], "category": entry["category"]}


def call_tool(
    name: str, arguments: dict[str, Any] | None = None, profile: ProfileArg = None
) -> Any:
    """Dispatch to the real api function named *name*.

    Resolves against the module that owns the tool's schema — currently
    ``pyirena.api.control``, ``pyirena.api.calculators`` or
    ``pyirena.api.data_ops``.

    Returns whatever that function returns unchanged (a plain dict, or a
    dict shaped like an image result — ``pyirena/mcp/server.py`` is
    responsible for turning the latter into MCP image content).
    """
    prof = resolve_profile(profile)
    entry = _REGISTRY.get(name)
    if entry is None:
        return _unknown_tool_error(name, prof)
    if not _visible(entry, prof):
        return _hidden_error(name, entry, prof)
    try:
        return getattr(entry["module"], name)(**(arguments or {}))
    except TypeError as exc:
        return {
            "error": f"Bad arguments for '{name}': {exc}",
            "code": "BAD_ARGUMENTS",
            "suggestion": f"Call pyirena_describe_tool('{name}') to see the expected arguments.",
        }
