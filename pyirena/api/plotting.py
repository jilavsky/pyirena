"""Headless plotting for the api surface — matplotlib (Agg) only.

Both helpers save a PNG and optionally return base64 so an MCP client can
inline-display the image.
"""
from __future__ import annotations

import base64
import os
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np

from pyirena.api._paths import resolve_safe, resolve_safe_file
from pyirena.api.aggregate import tabulate_parameter
from pyirena.api.schemas import PlotResult


def _plot_cache_dir() -> Path:
    """Where plot PNGs are written by default. Overridable via PYIRENA_PLOT_CACHE."""
    cache = os.environ.get("PYIRENA_PLOT_CACHE")
    if cache:
        p = Path(cache).expanduser()
    else:
        p = Path(tempfile.gettempdir()) / "pyirena-mcp"
    p.mkdir(parents=True, exist_ok=True)
    return p


def _resolve_output(output_path: Optional[str], default_name: str) -> Path:
    if output_path:
        # User-supplied path: still subject to PYIRENA_DATA_ROOT if set
        p = resolve_safe(output_path, must_exist=False)
        p.parent.mkdir(parents=True, exist_ok=True)
        return p
    return _plot_cache_dir() / default_name


def _encode_base64(path: Path) -> str:
    return base64.b64encode(path.read_bytes()).decode("ascii")


def plot_iq(
    paths: list[str],
    overlay: bool = True,
    log_x: bool = True,
    log_y: bool = True,
    output_path: Optional[str] = None,
    return_base64: bool = True,
    dpi: int = 120,
) -> dict:
    """Plot I(Q) for one or more files. Saves a PNG, optionally inlines base64.

    Parameters
    ----------
    paths : list of str
        SAS files (NXcanSAS HDF5).
    overlay : bool
        True (default): one plot, all curves overlaid. False: one subplot per file.
    log_x, log_y : bool
        Log axes (default True). For WAXS data set both False.
    output_path : str, optional
        Where to save the PNG. Defaults to a temp file under PYIRENA_PLOT_CACHE.
    return_base64 : bool
        When True, include the PNG content as base64 in the result. Useful
        for MCP image content blocks.
    dpi : int
        Output resolution.
    """
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    if not paths:
        raise ValueError("plot_iq requires at least one path")
    files = [resolve_safe_file(p) for p in paths]

    # Load reduced data via existing api function (decimation off — full curves)
    from pyirena.api.data import read_reduced_data
    curves = []
    for fp in files:
        rd = read_reduced_data(str(fp), include_full=True)
        if not rd.get("found"):
            continue
        curves.append((fp.name, rd))
    if not curves:
        raise RuntimeError("None of the supplied files contain readable reduced data")

    if overlay:
        fig, ax = plt.subplots(figsize=(8, 6))
        for name, rd in curves:
            ax.plot(rd["Q"], rd["I"], marker="o", linestyle="", markersize=2.5,
                    label=name)
        if log_x:
            ax.set_xscale("log")
        if log_y:
            ax.set_yscale("log")
        ax.set_xlabel(r"Q  (Å$^{-1}$)")
        ax.set_ylabel(r"I  (cm$^{-1}$)")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8, loc="best")
        fig.tight_layout()
    else:
        n = len(curves)
        cols = min(2, n)
        rows = (n + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(8 * cols, 5 * rows),
                                  squeeze=False)
        for ax, (name, rd) in zip(axes.flat, curves):
            ax.plot(rd["Q"], rd["I"], marker="o", linestyle="", markersize=2.5)
            if log_x:
                ax.set_xscale("log")
            if log_y:
                ax.set_yscale("log")
            ax.set_xlabel(r"Q  (Å$^{-1}$)")
            ax.set_ylabel(r"I  (cm$^{-1}$)")
            ax.set_title(name, fontsize=10)
            ax.grid(True, which="both", alpha=0.3)
        # Hide unused axes
        for ax in axes.flat[len(curves):]:
            ax.axis("off")
        fig.tight_layout()

    out = _resolve_output(output_path, default_name="pyirena_iq.png")
    fig.savefig(out, dpi=dpi)
    w_px, h_px = fig.canvas.get_width_height()
    plt.close(fig)

    result = PlotResult(
        path=str(out),
        format="png",
        width=int(w_px),
        height=int(h_px),
        n_files=len(curves),
    )
    if return_base64:
        result.base64_png = _encode_base64(out)
    return result.to_dict()


def plot_parameter_trend(
    folder: str,
    tool: str,
    parameter: str,
    x_axis: str = "scan_number",
    subgroup_index: Optional[int] = None,
    sample_filter: Optional[str] = None,
    output_path: Optional[str] = None,
    return_base64: bool = True,
    dpi: int = 120,
) -> dict:
    """Plot a per-file parameter across many files in a folder.

    Internally calls :func:`tabulate_parameter` and renders the resulting
    rows as a line+marker plot with stddev error bars when available.
    """
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    tab = tabulate_parameter(
        folder=folder, tool=tool, parameter=parameter, x_axis=x_axis,
        subgroup_index=subgroup_index, sample_filter=sample_filter,
    )
    rows = tab.get("rows") or []
    if not rows:
        raise RuntimeError(
            f"No rows found for tool={tool} parameter={parameter} in {folder}"
        )

    xs: list[float] = []
    ys: list[float] = []
    es: list[float] = []
    labels: list[str] = []
    for i, r in enumerate(rows):
        v = r.get("value")
        if v is None:
            continue
        if x_axis == "scan_number" and r.get("scan_number") is not None:
            xs.append(float(r["scan_number"]))
        else:
            xs.append(float(i))
        ys.append(float(v))
        es.append(float(r["stddev"]) if r.get("stddev") is not None else float("nan"))
        labels.append(r.get("name") or "")
    if not xs:
        raise RuntimeError(
            f"All {len(rows)} rows had value=None for {tool}.{parameter}"
        )

    fig, ax = plt.subplots(figsize=(9, 5))
    has_err = any(np.isfinite(e) for e in es)
    if has_err:
        ax.errorbar(xs, ys, yerr=[None if np.isnan(e) else e for e in es],
                    marker="o", linestyle="-", capsize=3)
    else:
        ax.plot(xs, ys, marker="o", linestyle="-")
    ax.set_xlabel(x_axis if x_axis != "scan_number" else "Scan number")
    units = tab.get("units")
    ylab = tab.get("label") or parameter
    if units:
        ylab = f"{ylab}  ({units})"
    ax.set_ylabel(ylab)
    ax.set_title(f"{tool} — {parameter}")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    out = _resolve_output(output_path,
                          default_name=f"pyirena_trend_{tool}_{parameter}.png")
    fig.savefig(out, dpi=dpi)
    w_px, h_px = fig.canvas.get_width_height()
    plt.close(fig)

    result = PlotResult(
        path=str(out),
        format="png",
        width=int(w_px),
        height=int(h_px),
        n_files=len(rows),
    )
    if return_base64:
        result.base64_png = _encode_base64(out)
    return result.to_dict()
