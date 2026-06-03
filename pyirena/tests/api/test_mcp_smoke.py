"""Smoke test for the MCP server module — only that it imports and
registers tools. Spawning a subprocess is overkill for CI."""
from __future__ import annotations

import pytest


def test_mcp_module_imports_and_registers_tools():
    pytest.importorskip("mcp")
    from pyirena.mcp.server import mcp
    # FastMCP exposes its registered tools via list_tools() coroutine.
    # We don't need to enter the event loop — just check the internal
    # tool registry has entries.
    # FastMCP stores tools in a _tools dict or similar (depends on version);
    # safest is to introspect the public interface if it exists.
    expected = {
        "pyirena_list_files", "pyirena_summarize_folder", "pyirena_inspect_file",
        "pyirena_read_reduced_data", "pyirena_read_metadata",
        "pyirena_read_simple_fit", "pyirena_read_unified_fit",
        "pyirena_read_size_distribution",
        "pyirena_read_modeling", "pyirena_read_saxs_morph",
        "pyirena_read_waxs_peakfit",
        "pyirena_read_fractals", "pyirena_read_merge_provenance",
        "pyirena_read_manipulation_provenance",
        "pyirena_tabulate_parameter", "pyirena_summarize_sample",
        "pyirena_plot_iq", "pyirena_plot_parameter_trend",
    }
    # Try the documented public method first
    if hasattr(mcp, "list_tools"):
        import asyncio
        tools = asyncio.run(mcp.list_tools())
        names = {t.name for t in tools}
    else:
        # Fallback: introspect a known private attribute
        registry = getattr(mcp, "_tool_manager", None) or getattr(mcp, "_tools", {})
        if hasattr(registry, "_tools"):
            names = set(registry._tools.keys())
        elif isinstance(registry, dict):
            names = set(registry.keys())
        else:
            pytest.skip("Could not introspect FastMCP tool registry")
    missing = expected - names
    assert not missing, f"MCP server is missing tools: {missing}"
