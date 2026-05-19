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
        "list_files", "summarize_folder", "inspect_file",
        "read_reduced_data", "read_metadata",
        "read_simple_fit", "read_unified_fit", "read_size_distribution",
        "read_modeling", "read_saxs_morph", "read_waxs_peakfit",
        "read_fractals", "read_merge_provenance",
        "read_manipulation_provenance",
        "tabulate_parameter", "summarize_sample",
        "plot_iq", "plot_parameter_trend",
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
