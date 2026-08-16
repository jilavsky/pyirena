"""Regression tests for how the MCP server returns images.

Two things used to be broken and both are guarded here.

1. FastMCP >= 1.10 derives a *structured output* JSON schema from a tool's
   return annotation and serialises the return value against it.  A
   ``mcp.server.fastmcp.Image`` is not JSON-serialisable, so every image tool
   died with "Unable to serialize unknown type: ... Image" *after* writing the
   PNG — Claude Desktop reported a serialization bug and no picture.  Image
   tools must therefore be registered with structured output disabled.

2. The control-API image tools returned base64 only.  Clients that do not
   render inline images (AnythingLLM in some modes, plain agent loops) had no
   way to reach the plot at all and told users no file existed.  Every image
   tool now emits a text block carrying the PNG's on-disk path *and* the
   inline image block.
"""
from __future__ import annotations

import asyncio
import base64

import pytest

pytest.importorskip("mcp")

# One-pixel PNG — enough to prove the bytes survive the round trip.
_PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg=="
)
_PNG_B64 = base64.b64encode(_PNG).decode("ascii")

# Every tool that returns a picture, with the pyirena.api.control function it
# delegates to. The control functions are stubbed so the test needs no data
# file and no matplotlib run.
IMAGE_TOOLS = {
    "pyirena_ctrl_get_fit_image": "get_fit_image",
    "pyirena_ctrl_get_residuals_image": "get_residuals_image",
    "pyirena_ctrl_sizes_get_background_image": "get_background_preview_image",
    "pyirena_ctrl_sizes_get_fit_image": "get_sizes_fit_image",
    "pyirena_ctrl_simple_get_fit_image": "get_simple_fit_image",
    "pyirena_ctrl_simple_get_linearization_image": "get_simple_linearization_image",
    "pyirena_ctrl_modeling_get_fit_image": "get_modeling_fit_image",
    "pyirena_ctrl_waxs_get_fit_image": "get_waxs_fit_image",
}
PLOT_TOOLS = ("pyirena_plot_iq", "pyirena_plot_parameter_trend")


def _content(result):
    """FastMCP returns either content or (content, structured_content)."""
    return result[0] if isinstance(result, tuple) else result


def _split(content):
    from mcp.types import ImageContent, TextContent

    texts = [c for c in content if isinstance(c, TextContent)]
    images = [c for c in content if isinstance(c, ImageContent)]
    return texts, images


def test_image_tools_have_no_structured_output_schema():
    """The registration bug: a structured schema means the Image cannot serialize."""
    from pyirena.mcp.server import mcp

    tools = asyncio.run(mcp.list_tools())
    expected = set(IMAGE_TOOLS) | set(PLOT_TOOLS)
    by_name = {t.name: t for t in tools}

    missing = expected - set(by_name)
    assert not missing, f"image tools missing from the server: {missing}"

    structured = sorted(n for n in expected if by_name[n].outputSchema)
    assert not structured, (
        "these image tools declare a structured output schema; FastMCP will try "
        f"to JSON-serialize the Image and fail at call time: {structured}"
    )


@pytest.mark.parametrize("tool_name,ctrl_func", sorted(IMAGE_TOOLS.items()))
def test_control_image_tool_returns_text_path_and_inline_image(
    tool_name, ctrl_func, tmp_path, monkeypatch
):
    from pyirena.mcp import server

    png_path = tmp_path / "stub.png"
    png_path.write_bytes(_PNG)

    monkeypatch.setattr(
        server._ctrl,
        ctrl_func,
        lambda *a, **k: {
            "ok": True,
            "image_base64": _PNG_B64,
            "image_path": str(png_path),
            "has_residuals": True,
        },
        raising=True,
    )

    content = _content(asyncio.run(server.mcp.call_tool(tool_name, {"session_id": "s1"})))
    texts, images = _split(content)

    assert len(images) == 1, f"{tool_name} returned no inline image block"
    assert images[0].mimeType == "image/png"
    assert base64.b64decode(images[0].data) == _PNG

    assert texts, f"{tool_name} returned no text block"
    assert str(png_path) in texts[0].text, (
        "the text block must carry the PNG path so clients that cannot render "
        "images inline can still open the file"
    )


@pytest.mark.parametrize("tool_name,ctrl_func", sorted(IMAGE_TOOLS.items()))
def test_control_image_tool_passes_errors_through_as_text(tool_name, ctrl_func, monkeypatch):
    """An error dict must stay a readable error, not become a broken image."""
    from mcp.types import ImageContent

    from pyirena.mcp import server

    monkeypatch.setattr(
        server._ctrl,
        ctrl_func,
        lambda *a, **k: {"error": "No session found.", "code": "NO_SESSION"},
        raising=True,
    )

    content = _content(asyncio.run(server.mcp.call_tool(tool_name, {"session_id": "nope"})))
    assert not any(isinstance(c, ImageContent) for c in content)
    assert "NO_SESSION" in "".join(getattr(c, "text", "") for c in content)


@pytest.mark.parametrize("tool_name", PLOT_TOOLS)
def test_plot_tools_return_text_path_and_inline_image(tool_name, tmp_path, monkeypatch):
    from pyirena.mcp import server

    png_path = tmp_path / "plot.png"
    png_path.write_bytes(_PNG)

    api_func = "plot_iq" if tool_name == "pyirena_plot_iq" else "plot_parameter_trend"
    monkeypatch.setattr(
        server.papi, api_func, lambda *a, **k: {"path": str(png_path)}, raising=True
    )

    args = (
        {"paths": ["x.h5"]}
        if tool_name == "pyirena_plot_iq"
        else {"folder": ".", "tool": "unified_fit", "parameter": "Rg"}
    )
    content = _content(asyncio.run(server.mcp.call_tool(tool_name, args)))
    texts, images = _split(content)

    assert len(images) == 1
    assert base64.b64decode(images[0].data) == _PNG
    assert str(png_path) in texts[0].text


def test_control_api_image_results_carry_a_path(monkeypatch, tmp_path):
    """render_png is the single place PNGs are written; it must report the path."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    from pyirena.api.control import _images

    monkeypatch.setenv("PYIRENA_PLOT_CACHE", str(tmp_path))

    fig, ax = plt.subplots()
    ax.plot([1, 2], [1, 2])
    b64, path = _images.render_png(fig, "unit_test_fig", dpi=50)

    assert path == str(tmp_path / "unit_test_fig.png")
    assert (tmp_path / "unit_test_fig.png").is_file()
    assert base64.b64decode(b64) == (tmp_path / "unit_test_fig.png").read_bytes()
