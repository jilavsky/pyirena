"""Tests for dispatcher profiles — what a transport is allowed to expose.

A stdio MCP server shares a filesystem with its client, so a tool taking
``file_path`` or returning ``image_path`` is meaningful. The ZMQ service
(``planning/zmq-service/``) runs on a different machine with no shared
storage, where the same tool would silently point at a file the caller can
never read. The ``json_only`` profile is what makes that difference explicit.

The load-bearing test here is
``test_every_path_taking_schema_is_flagged``: the flags are derived from
naming rules rather than hand-maintained, so the rules must be checked
against the actual schemas or a future tool will quietly leak a server path.
"""
from __future__ import annotations

import pytest

from pyirena.api import dispatch

PATHISH = ("path", "file", "folder", "dir", "output")


def _flagged(predicate) -> set[str]:
    return {n for n, e in dispatch._REGISTRY.items() if predicate(e)}


# ---------------------------------------------------------------------------
# Flag derivation
# ---------------------------------------------------------------------------

def test_every_path_taking_schema_is_flagged():
    """Any tool with a path-like argument must be touches_files.

    Guards the naming rules in ``_touches_files`` against a new tool that
    takes a path but is not called ``save_*`` and is not a data op.
    """
    leaks = []
    for name, entry in dispatch._REGISTRY.items():
        if entry["touches_files"]:
            continue
        props = (entry["schema"].get("input_schema") or {}).get("properties") or {}
        hits = [p for p in props if any(k in p.lower() for k in PATHISH)]
        if hits:
            leaks.append((name, hits))
    assert not leaks, (
        "these tools take a path-like argument but are not flagged "
        f"touches_files: {leaks}"
    )


def test_the_flagged_sets_are_non_empty_and_plausible():
    files = _flagged(lambda e: e["touches_files"])
    images = _flagged(lambda e: e["returns_image"])

    assert files and images
    assert {"save_fit", "save_sizes_fit", "save_carbon_fit"} <= files
    # Every data operation reads and writes files.
    assert {n for n, e in dispatch._REGISTRY.items() if e["category"] == "data"} <= files
    assert {"get_fit_image", "get_residuals_image"} <= images
    # A tool is not both by accident.
    assert not (files & images)


def test_pure_computation_is_not_flagged():
    for name in ("run_fit", "set_parameter_value", "get_chi_squared", "calc_contrast"):
        flags = dispatch.tool_flags(name)
        assert flags["touches_files"] is False, name
        assert flags["returns_image"] is False, name


# ---------------------------------------------------------------------------
# Profile behaviour
# ---------------------------------------------------------------------------

def test_json_only_hides_every_file_and_image_tool():
    hidden = _flagged(lambda e: e["touches_files"] or e["returns_image"])
    visible = set()
    for cat in (c["name"] for c in dispatch.list_categories("json_only")["categories"]):
        visible |= {t["name"] for t in dispatch.list_tools(cat, "json_only")["tools"]}

    assert not (visible & hidden)
    assert visible == set(dispatch._REGISTRY) - hidden


def test_json_images_hides_files_but_keeps_images():
    visible = set()
    for cat in (c["name"] for c in dispatch.list_categories("json_images")["categories"]):
        visible |= {t["name"] for t in dispatch.list_tools(cat, "json_images")["tools"]}

    assert "get_fit_image" in visible
    assert "save_fit" not in visible


def test_profile_all_is_the_default_and_hides_nothing():
    assert set(dispatch.list_tools("unified")["tools"][0]) == {"name", "summary"}
    for cat in ("unified", "data"):
        assert (
            dispatch.list_tools(cat)["tools"] == dispatch.list_tools(cat, "all")["tools"]
        )
    assert "save_fit" in {t["name"] for t in dispatch.list_tools("unified")["tools"]}


def test_the_data_category_disappears_entirely_under_json_only():
    names = {c["name"] for c in dispatch.list_categories("json_only")["categories"]}
    assert "data" not in names
    assert "unified" in names
    # ...but the six fitting categories and the calculators all survive.
    assert names == {"unified", "sizes", "simple", "modeling", "waxs", "carbon",
                     "results", "calculators"}


def test_describe_and_call_refuse_a_hidden_tool():
    described = dispatch.describe_tool("save_fit", "json_only")
    assert described["code"] == "NOT_AVAILABLE_IN_PROFILE"
    assert "json_only" in described["error"]

    called = dispatch.call_tool("save_fit", {"session_id": "nope"}, "json_only")
    assert called["code"] == "NOT_AVAILABLE_IN_PROFILE"
    # Refused before dispatch: no session lookup happened.
    assert "NO_SESSION" not in str(called)


def test_a_hidden_tool_is_not_offered_as_a_near_miss_suggestion():
    r = dispatch.describe_tool("save_ft", "json_only")
    assert r["code"] == "UNKNOWN_TOOL"
    assert "save_fit" not in r["suggestion"]
    # ...but it is under the permissive profile.
    assert "save_fit" in dispatch.describe_tool("save_ft", "all")["suggestion"]


def test_unknown_profile_name_is_a_programming_error():
    with pytest.raises(ValueError, match="Unknown dispatcher profile"):
        dispatch.resolve_profile("open_season")


def test_open_dataset_from_data_is_a_lifecycle_tool_not_a_dispatched_one():
    """It is the JSON transport's way in, so it stays top-level like open_dataset."""
    assert "open_dataset_from_data" in dispatch.SESSION_LIFECYCLE_NAMES
    assert "open_dataset_from_data" not in dispatch._REGISTRY


def test_the_shim_under_mcp_still_works():
    from pyirena.mcp import dispatch as shim

    assert shim.call_tool is dispatch.call_tool
    assert shim._REGISTRY is dispatch._REGISTRY
