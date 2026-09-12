"""Tests for pyirena.mcp.dispatch — the fixed dispatcher over
pyirena.api.control, pyirena.api.calculators and pyirena.api.data_ops
that keeps pyirena-mcp's
tool count small (see
pyirena/mcp/server.py's "Control API — dispatcher" section and
planning/ai-agent/01-api-and-mcp-extensions.md's "Related work" section).

Deliberately imports only pyirena.mcp.dispatch, not pyirena.mcp.server, so
these run without the optional ``mcp`` package installed.
"""
from __future__ import annotations

import pytest

from pyirena.api.calculator_schemas import CALCULATOR_SCHEMA_BY_NAME
from pyirena.api.control.schemas import TOOL_SCHEMA_BY_NAME
from pyirena.api.data_op_schemas import DATA_OP_SCHEMA_BY_NAME
from pyirena.mcp import dispatch

# Every category the dispatcher serves: the five fitting tools, the
# stateless calculators (pyirena.api.calculators) and the data operations
# (pyirena.api.data_ops).
ALL_CATEGORIES = [
    "unified", "sizes", "simple", "modeling", "waxs", "calculators", "data",
]


def test_session_lifecycle_excluded_from_dispatcher():
    assert dispatch.SESSION_LIFECYCLE_NAMES.isdisjoint(dispatch._REGISTRY)
    assert dispatch.SESSION_LIFECYCLE_NAMES <= set(TOOL_SCHEMA_BY_NAME)


def test_registry_covers_every_non_session_schema():
    expected = (
        set(TOOL_SCHEMA_BY_NAME)
        | set(CALCULATOR_SCHEMA_BY_NAME)
        | set(DATA_OP_SCHEMA_BY_NAME)
    ) - dispatch.SESSION_LIFECYCLE_NAMES
    assert set(dispatch._REGISTRY) == expected


def test_schema_registries_do_not_collide():
    """Dispatcher tool names are a flat namespace across all schema sources."""
    assert set(TOOL_SCHEMA_BY_NAME).isdisjoint(CALCULATOR_SCHEMA_BY_NAME)
    assert set(TOOL_SCHEMA_BY_NAME).isdisjoint(DATA_OP_SCHEMA_BY_NAME)
    assert set(CALCULATOR_SCHEMA_BY_NAME).isdisjoint(DATA_OP_SCHEMA_BY_NAME)


def test_list_categories_counts_match_registry():
    result = dispatch.list_categories()
    names = {c["name"] for c in result["categories"]}
    assert names == set(ALL_CATEGORIES)
    assert sum(c["tool_count"] for c in result["categories"]) == len(dispatch._REGISTRY)
    for c in result["categories"]:
        assert c["tool_count"] > 0, f"category {c['name']} has no tools"


@pytest.mark.parametrize("category", ALL_CATEGORIES)
def test_list_tools_names_round_trip_through_describe_tool(category):
    listed = dispatch.list_tools(category)
    assert listed["category"] == category
    assert listed["tools"], f"category {category} listed no tools"
    for tool in listed["tools"]:
        assert tool["summary"]
        schema = dispatch.describe_tool(tool["name"])
        assert schema["name"] == tool["name"]
        assert schema["category"] == category
        assert "input_schema" in schema


def test_list_tools_unknown_category_lists_valid_ones():
    result = dispatch.list_tools("bogus")
    assert result["code"] == "UNKNOWN_CATEGORY"
    for name in ALL_CATEGORIES:
        assert name in result["suggestion"]


def test_describe_tool_unknown_name_suggests_close_match():
    result = dispatch.describe_tool("run_ft")  # close to run_fit
    assert result["code"] == "UNKNOWN_TOOL"
    assert "run_fit" in result["suggestion"]


def test_describe_tool_session_lifecycle_name_points_at_top_level_tool():
    result = dispatch.describe_tool("open_dataset")
    assert result["code"] == "UNKNOWN_TOOL"
    assert "pyirena_ctrl_open_dataset" in result["suggestion"]


def test_call_tool_dispatches_to_the_real_control_function(monkeypatch):
    import pyirena.api.control as ctrl

    monkeypatch.setattr(ctrl, "run_fit", lambda session_id, **kw: {"ok": True, "session_id": session_id})
    result = dispatch.call_tool("run_fit", {"session_id": "s1"})
    assert result == {"ok": True, "session_id": "s1"}


def test_call_tool_bad_arguments_becomes_a_readable_error():
    result = dispatch.call_tool("run_fit", {"session_id": "s1", "not_a_real_kwarg": 1})
    assert result["code"] == "BAD_ARGUMENTS"
    assert "run_fit" in result["error"]


def test_call_tool_unknown_name():
    result = dispatch.call_tool("bogus_tool_name", {})
    assert result["code"] == "UNKNOWN_TOOL"


def test_call_tool_session_lifecycle_name_is_unreachable():
    result = dispatch.call_tool("open_dataset", {"file_path": "x.h5"})
    assert result["code"] == "UNKNOWN_TOOL"
