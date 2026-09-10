"""Fixed dispatcher over ``pyirena.api.control`` — the MCP tool-count fix.

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
Session-lifecycle functions (``open_dataset``, ``list_open_sessions``,
``close_session``, ``get_session_summary``) are excluded here and stay
top-level MCP tools instead, since nearly every workflow starts there.

The dispatcher is not control-only: ``_SOURCES`` below lists every schema
registry it serves. ``pyirena.api.calculators`` joins as the stateless
"calculators" category, so support calculators (scattering contrast today,
more later) reach agents without adding a single top-level MCP tool.

Design doc: AIDA's ``planning/mcp_tool_scaling.md`` (Tier 2) and pyIrena's
``planning/ai-agent/01-api-and-mcp-extensions.md`` ("Related work" section).

This module has no dependency on the ``mcp`` package on purpose, so it can
be imported and tested without it installed. ``pyirena/mcp/server.py`` is
the only place that turns a result into MCP content (e.g. image blocks).
"""
from __future__ import annotations

import difflib
import inspect
from typing import Any

from pyirena.api import calculators as _calc
from pyirena.api import control as _ctrl
from pyirena.api.calculator_schemas import CALCULATOR_SCHEMA_BY_NAME
from pyirena.api.control.schemas import TOOL_SCHEMA_BY_NAME

SESSION_LIFECYCLE_NAMES = frozenset(
    {"open_dataset", "list_open_sessions", "close_session", "get_session_summary"}
)

_CATEGORY_BY_MODULE = {
    "unified_fit": "unified",
    "sizes": "sizes",
    "simple_fits": "simple",
    "modeling": "modeling",
    "waxs_peakfit": "waxs",
}

_CATEGORY_BLURBS = {
    "unified": "Unified Fit — Beaucage multi-level model for a whole multi-level SAXS/USAXS curve.",
    "sizes": "Size Distribution — invert a dilute single population to a size histogram P(r).",
    "simple": "Simple Fits — one analytical model (Guinier, Porod, Sphere, ...) over a Q sub-range.",
    "modeling": "Modeling — multi-population forward model with specific form/structure factors.",
    "waxs": "WAXS Peak Fit — peak position, width and area for wide-angle patterns.",
    "calculators": (
        "Calculators — stateless support calculations that need no dataset "
        "and open no session (scattering contrast, SLDs, element lookup)."
    ),
}


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
            registry[name] = {
                "category": category_for(name),
                "schema": schema,
                "module": module,
            }
    return registry


_REGISTRY = _registry()


def _first_sentence(description: str) -> str:
    text = description.strip().split("\n", 1)[0]
    for sep in (". ", ".\n"):
        if sep in text:
            return text.split(sep, 1)[0] + "."
    return text


def _unknown_tool_error(name: str) -> dict:
    if name in SESSION_LIFECYCLE_NAMES:
        return {
            "error": f"'{name}' is a session-lifecycle tool, not a dispatched one.",
            "code": "UNKNOWN_TOOL",
            "suggestion": f"Call the top-level 'pyirena_ctrl_{name}' tool directly.",
        }
    candidates = difflib.get_close_matches(name, _REGISTRY.keys(), n=3)
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


def list_categories() -> dict:
    """Return every dispatch category with its tool count."""
    counts: dict[str, int] = {}
    for entry in _REGISTRY.values():
        counts[entry["category"]] = counts.get(entry["category"], 0) + 1
    categories = [
        {"name": name, "description": _CATEGORY_BLURBS[name], "tool_count": counts.get(name, 0)}
        for name in _CATEGORY_BLURBS
    ]
    return {"categories": categories}


def list_tools(category: str) -> dict:
    """Return the name and one-line summary of every tool in *category*."""
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
            if entry["category"] == category
        ),
        key=lambda t: t["name"],
    )
    return {"category": category, "tools": tools}


def describe_tool(name: str) -> dict:
    """Return the full Anthropic-style tool schema for *name*, plus its category."""
    entry = _REGISTRY.get(name)
    if entry is None:
        return _unknown_tool_error(name)
    return {**entry["schema"], "category": entry["category"]}


def call_tool(name: str, arguments: dict[str, Any] | None = None) -> Any:
    """Dispatch to the real api function named *name*.

    Resolves against the module that owns the tool's schema — currently
    ``pyirena.api.control`` or ``pyirena.api.calculators``.

    Returns whatever that function returns unchanged (a plain dict, or a
    dict shaped like an image result — ``pyirena/mcp/server.py`` is
    responsible for turning the latter into MCP image content).
    """
    entry = _REGISTRY.get(name)
    if entry is None:
        return _unknown_tool_error(name)
    try:
        return getattr(entry["module"], name)(**(arguments or {}))
    except TypeError as exc:
        return {
            "error": f"Bad arguments for '{name}': {exc}",
            "code": "BAD_ARGUMENTS",
            "suggestion": f"Call pyirena_describe_tool('{name}') to see the expected arguments.",
        }
