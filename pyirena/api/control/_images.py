"""Shared PNG rendering for the control API.

Every control-API image function returns *both* the base64 payload and the
absolute path of the same PNG on disk:

    {"ok": True, "image_base64": "...", "image_path": "/tmp/pyirena-ctrl/fit_s1.png"}

The base64 payload is what ``pyirena/mcp/server.py`` turns into an MCP image
content block, which is how Claude Desktop, ChatGPT and most other MCP clients
display the picture inline.  ``image_path`` is the fallback: clients that
cannot render inline images (or agents that want to hand the file to something
else) still get a real file they can open.  Reporting only the base64 payload
left such clients with no way to reach the plot at all.

The destination folder is ``$PYIRENA_PLOT_CACHE`` when set — the same override
``pyirena.api.plotting`` uses — otherwise a ``pyirena-ctrl`` folder inside the
system temp directory.
"""
from __future__ import annotations

import base64
import os
import tempfile
from pathlib import Path
from typing import Any


def image_cache_dir() -> Path:
    """Directory control-API PNGs are written to (created if missing)."""
    cache = os.environ.get("PYIRENA_PLOT_CACHE")
    if cache:
        path = Path(cache).expanduser()
    else:
        path = Path(tempfile.gettempdir()) / "pyirena-ctrl"
    path.mkdir(parents=True, exist_ok=True)
    return path


def render_png(fig, name: str, dpi: int = 120) -> tuple[str, str]:
    """Save *fig* as ``<cache>/<name>.png``.

    Args:
        fig: A matplotlib figure. It is closed before returning.
        name: File stem, e.g. ``"sizes_fit_s1"``. Must be unique per session.
        dpi: Output resolution.

    Returns:
        ``(base64_png, absolute_path)`` — the same image in both forms.
    """
    import matplotlib.pyplot as plt  # noqa: PLC0415

    out = image_cache_dir() / f"{name}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(out.read_bytes()).decode("ascii"), str(out)


def image_result(fig, name: str, dpi: int = 120, **extra: Any) -> dict:
    """Render *fig* and wrap it in the standard control-API image dict."""
    b64, path = render_png(fig, name, dpi)
    return {"ok": True, "image_base64": b64, "image_path": path, **extra}
