"""MCP stdio server exposing pyirena.api tools.

Run via:    pyirena-mcp
Or in code: python -m pyirena.mcp.server

Environment overrides (see also pyirena.api):
    PYIRENA_DATA_ROOT       restrict file access to this subtree
    PYIRENA_MAX_ARRAY_POINTS default decimation cap (default 500)
    PYIRENA_PLOT_CACHE      where plot PNGs are saved (default tempdir)
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

try:
    from mcp.server.fastmcp import FastMCP, Image
except ImportError as exc:  # pragma: no cover - import guard
    raise ImportError(
        "The 'mcp' package is required to run pyirena-mcp. "
        "Install with: pip install pyirena[mcp]"
    ) from exc

from pyirena import api as papi


mcp = FastMCP(
    "pyirena",
    instructions=(
        "Tools for reading pyirena SAXS/USAXS analysis results from NXcanSAS "
        "HDF5 files. Start with summarize_folder() or list_files() to "
        "discover what is available, then drill in with inspect_file() or "
        "one of the read_<tool>() functions. Use plot_iq() / "
        "plot_parameter_trend() to visualize."
    ),
)


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

@mcp.tool()
def list_files(
    folder: str,
    pattern: str = "*.h5,*.hdf5,*.hdf,*.nx,*.nxs",
    sort: str = "mtime_desc",
    limit: int = 100,
    deep: bool = True,
) -> list[dict]:
    """List HDF5 files in *folder* with metadata.

    Each entry includes path, name, sample, scan_number, mtime, size, and
    the list of pyirena analyses present (e.g. ['unified_fit', 'modeling']).
    Use sort='mtime_desc' (default) to get the latest files first.
    """
    return papi.list_files(folder=folder, pattern=pattern, sort=sort,
                           limit=limit, deep=deep)


@mcp.tool()
def summarize_folder(folder: str, sample_filter: Optional[str] = None) -> dict:
    """Get an aggregate snapshot of a folder of SAS data.

    Returns file count, unique samples, per-analysis file counts, and mtime
    range. Cheap orientation call — use it BEFORE drilling into individual
    files. Optionally filter to one sample (case-insensitive substring).
    """
    return papi.summarize_folder(folder=folder, sample_filter=sample_filter)


@mcp.tool()
def inspect_file(path: str) -> dict:
    """Inspect a single file: sample name, analyses present, Q range, n_points."""
    return papi.inspect_file(path)


# ---------------------------------------------------------------------------
# Reading reduced data + metadata
# ---------------------------------------------------------------------------

@mcp.tool()
def read_reduced_data(path: str, decimate: int = 500,
                      include_full: bool = False) -> dict:
    """Read the raw reduced I(Q) curve from a SAS file.

    Arrays are decimated to ~*decimate* points by default to keep the
    response compact. Set include_full=True for full fidelity (avoid in
    LLM workflows — long arrays bloat context).
    """
    return papi.read_reduced_data(path=path, decimate=decimate,
                                   include_full=include_full)


@mcp.tool()
def read_metadata(path: str) -> dict:
    """Read sample / experiment metadata from a SAS file."""
    return papi.read_metadata(path)


# ---------------------------------------------------------------------------
# Per-tool results
# ---------------------------------------------------------------------------

@mcp.tool()
def read_simple_fit(path: str, include_arrays: bool = False,
                    max_points: int = 500) -> dict:
    """Read Simple Fits results (Guinier, Porod, etc.).

    Arrays (Q, I_model, residuals) are omitted by default. Set
    include_arrays=True to include them (decimated to max_points).
    """
    return papi.read_simple_fit(path=path, include_arrays=include_arrays,
                                 max_points=max_points)


@mcp.tool()
def read_unified_fit(path: str, include_arrays: bool = False,
                     max_points: int = 500) -> dict:
    """Read Unified Fit (Beaucage) results — multi-level Rg/G/B/P + correlations."""
    return papi.read_unified_fit(path=path, include_arrays=include_arrays,
                                  max_points=max_points)


@mcp.tool()
def read_size_distribution(path: str, include_arrays: bool = False,
                           max_points: int = 500) -> dict:
    """Read Size Distribution fit results — Vf, Rg, r_grid, distribution."""
    return papi.read_size_distribution(path=path, include_arrays=include_arrays,
                                        max_points=max_points)


@mcp.tool()
def read_modeling(path: str, include_arrays: bool = False,
                  max_points: int = 500) -> dict:
    """Read parametric Modeling results (size_dist / unified_level / diff_peak / fractal pops)."""
    return papi.read_modeling(path=path, include_arrays=include_arrays,
                               max_points=max_points)


@mcp.tool()
def read_saxs_morph(path: str, include_arrays: bool = False,
                    max_points: int = 500) -> dict:
    """Read SAXS Morph results — voxelgram-based forward modeling output.

    Note: the 3-D voxelgram itself is intentionally not returned; only
    derived scalar parameters and 1-D curves.
    """
    return papi.read_saxs_morph(path=path, include_arrays=include_arrays,
                                 max_points=max_points)


@mcp.tool()
def read_waxs_peakfit(path: str, include_arrays: bool = False,
                      max_points: int = 500) -> dict:
    """Read WAXS Peak Fit results — per-peak Q0, FWHM, A, eta, area."""
    return papi.read_waxs_peakfit(path=path, include_arrays=include_arrays,
                                   max_points=max_points)


@mcp.tool()
def read_fractals(path: str) -> dict:
    """List fractal aggregates stored in a file (Z, df, dmin, c, Rg per aggregate)."""
    return papi.read_fractals(path)


@mcp.tool()
def read_merge_provenance(path: str) -> dict:
    """Read Data Merge provenance: scale, q_shift, background, source files."""
    return papi.read_merge_provenance(path)


@mcp.tool()
def read_manipulation_provenance(path: str) -> dict:
    """Read Data Manipulation provenance: operation, parameters, source file."""
    return papi.read_manipulation_provenance(path)


# ---------------------------------------------------------------------------
# Cross-file aggregation
# ---------------------------------------------------------------------------

@mcp.tool()
def tabulate_parameter(
    folder: str,
    tool: str,
    parameter: str,
    x_axis: str = "scan_number",
    subgroup_index: Optional[int] = None,
    sample_filter: Optional[str] = None,
) -> dict:
    """Extract one scalar parameter across every file in *folder*.

    *tool* is a pyirena tool key (e.g. 'unified_fit', 'modeling',
    'size_distribution'). *parameter* is a scalar key declared by that
    tool's schema (e.g. 'Rg', 'background', 'volume_fraction'). For
    per-subgroup parameters (Unified Fit levels, modeling populations,
    WAXS peaks) supply subgroup_index (1-based, default 1).
    """
    return papi.tabulate_parameter(
        folder=folder, tool=tool, parameter=parameter, x_axis=x_axis,
        subgroup_index=subgroup_index, sample_filter=sample_filter,
    )


@mcp.tool()
def summarize_sample(folder: str, sample: str) -> dict:
    """Condense everything known about one sample across a folder.

    File list, per-analysis file count, and min/max/n of every top-level
    scalar parameter the file's analysis tools declare.
    """
    return papi.summarize_sample(folder=folder, sample=sample)


# ---------------------------------------------------------------------------
# Plotting — returns inline images
# ---------------------------------------------------------------------------

@mcp.tool()
def plot_iq(
    paths: list[str],
    overlay: bool = True,
    log_x: bool = True,
    log_y: bool = True,
    output_path: Optional[str] = None,
) -> Image:
    """Plot I(Q) for one or more files; returns the PNG inline.

    overlay=True puts all curves on one axes; False uses a grid.
    Set log_x=log_y=False for WAXS (linear-linear).
    The PNG is also saved to disk (under PYIRENA_PLOT_CACHE or
    *output_path* if given).
    """
    result = papi.plot_iq(
        paths=paths, overlay=overlay, log_x=log_x, log_y=log_y,
        output_path=output_path, return_base64=False,
    )
    return Image(data=Path(result["path"]).read_bytes(), format="png")


@mcp.tool()
def plot_parameter_trend(
    folder: str,
    tool: str,
    parameter: str,
    x_axis: str = "scan_number",
    subgroup_index: Optional[int] = None,
    sample_filter: Optional[str] = None,
    output_path: Optional[str] = None,
) -> Image:
    """Plot a parameter trend across many files; returns the PNG inline.

    See tabulate_parameter() for the *tool* / *parameter* / *subgroup_index*
    semantics. Useful for time-series questions like "how is Rg evolving
    across the latest 30 scans?"
    """
    result = papi.plot_parameter_trend(
        folder=folder, tool=tool, parameter=parameter, x_axis=x_axis,
        subgroup_index=subgroup_index, sample_filter=sample_filter,
        output_path=output_path, return_base64=False,
    )
    return Image(data=Path(result["path"]).read_bytes(), format="png")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Run the MCP server over stdio (default transport)."""
    mcp.run()


if __name__ == "__main__":
    main()
