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


class _FastmcpBlocker:
    """Meta-path finder that makes ``import mcp.server.fastmcp`` fail.

    Lets us exercise the import guard in ``pyirena/mcp/server.py`` without
    actually installing a broken mcp.
    """

    def __init__(self, exc):
        self.exc = exc

    def find_spec(self, name, path=None, target=None):
        if name == "mcp.server.fastmcp" or name.startswith("mcp.server.fastmcp."):
            raise self.exc
        return None


def _import_server_with(monkeypatch, exc, reported_version):
    """Import pyirena.mcp.server with fastmcp broken; return the raised message."""
    import importlib
    import importlib.metadata
    import sys
    import types

    real_version = importlib.metadata.version

    def fake_version(name):
        if name == "mcp":
            if reported_version is None:
                raise importlib.metadata.PackageNotFoundError("mcp")
            return reported_version
        return real_version(name)

    monkeypatch.setattr(importlib.metadata, "version", fake_version)

    saved = {k: v for k, v in sys.modules.items() if k.startswith("pyirena.mcp")}
    for k in saved:
        monkeypatch.delitem(sys.modules, k, raising=False)
    monkeypatch.delitem(sys.modules, "mcp.server.fastmcp", raising=False)

    # The blocker below only fires for the leaf module, so the parent packages
    # have to resolve first -- otherwise, on a machine where mcp is not
    # installed at all (CI), the import dies at ``mcp`` with a plain
    # "No module named 'mcp'" and the injected exception is never seen.
    for parent in ("mcp", "mcp.server"):
        if parent not in sys.modules:
            stub = types.ModuleType(parent)
            stub.__path__ = []  # mark as a package so submodule import proceeds
            monkeypatch.setitem(sys.modules, parent, stub)

    blocker = _FastmcpBlocker(exc)
    sys.meta_path.insert(0, blocker)
    try:
        with pytest.raises(ImportError) as excinfo:
            importlib.import_module("pyirena.mcp.server")
        return str(excinfo.value)
    finally:
        sys.meta_path.remove(blocker)
        sys.modules.pop("pyirena.mcp.server", None)


def test_import_guard_when_mcp_missing(monkeypatch):
    msg = _import_server_with(
        monkeypatch, ModuleNotFoundError("No module named 'mcp'", name="mcp"), None
    )
    assert "pip install pyirena[mcp]" in msg


def test_import_guard_names_mcp_2_explicitly(monkeypatch):
    """mcp 2.x must not be misreported as 'mcp is not installed'.

    mcp 2.0 renamed FastMCP to MCPServer and left ``mcp.server.fastmcp`` as a
    stub raising ModuleNotFoundError (an ImportError subclass), so the guard
    has to distinguish it from a genuinely absent package.
    """
    stub_error = ModuleNotFoundError(
        "No module named 'mcp.server.fastmcp'. This is mcp 2.x, ...",
        name="mcp.server.fastmcp",
    )
    msg = _import_server_with(monkeypatch, stub_error, "2.0.1")
    assert "2.0.1" in msg
    assert "mcp>=1.0.0,<2.0" in msg
    assert "pip install pyirena[mcp]" not in msg  # the misleading advice


def test_import_guard_surfaces_other_failures(monkeypatch):
    msg = _import_server_with(
        monkeypatch, ImportError("cannot import name 'Image'"), "1.27.1"
    )
    assert "1.27.1" in msg
    assert "cannot import name 'Image'" in msg
