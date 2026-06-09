"""MCP stdio server exposing pyirena.api tools.

Run via:    pyirena-mcp
Or in code: python -m pyirena.mcp.server

All MCP tools are prefixed ``pyirena_`` so they are globally unambiguous
when the client connects to multiple MCP servers and so small LLMs are
not confused by clients that render ``server-tool`` with a dash.

Environment overrides (see also pyirena.api):
    PYIRENA_DATA_ROOT       restrict file access to this subtree
    PYIRENA_MAX_ARRAY_POINTS default decimation cap (default 500)
    PYIRENA_PLOT_CACHE      where plot PNGs are saved (default tempdir)
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Any


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
        "Tools for reading and fitting SAXS/USAXS data via pyirena. "
        "All tool names are prefixed 'pyirena_'. "
        "\n\n"
        "READ-ONLY tools (pyirena_ prefix): read existing fit results from NXcanSAS "
        "HDF5 files. Start with pyirena_summarize_folder() or pyirena_list_files() "
        "to discover available data, then pyirena_inspect_file() or "
        "pyirena_read_<tool>() to retrieve results. Use pyirena_plot_iq() / "
        "pyirena_plot_parameter_trend() to visualize."
        "\n\n"
        "CONTROL tools (pyirena_ctrl_ prefix): drive fitting interactively. "
        "Typical workflow: pyirena_ctrl_open_dataset() → session_id → "
        "pyirena_ctrl_select_model() → pyirena_ctrl_fix_all_except() → "
        "pyirena_ctrl_run_fit() → pyirena_ctrl_get_fit_image() → "
        "pyirena_ctrl_save_fit(). Sessions are in-memory for this server process."
    ),
)


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

@mcp.tool()
def pyirena_list_files(
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
def pyirena_summarize_folder(folder: str, sample_filter: Optional[str] = None) -> dict:
    """Get an aggregate snapshot of a folder of SAS data.

    Returns file count, unique samples, per-analysis file counts, and mtime
    range. Cheap orientation call — use it BEFORE drilling into individual
    files. Optionally filter to one sample (case-insensitive substring).
    """
    return papi.summarize_folder(folder=folder, sample_filter=sample_filter)


@mcp.tool()
def pyirena_inspect_file(path: str) -> dict:
    """Inspect a single file: sample name, analyses present, Q range, n_points."""
    return papi.inspect_file(path)


# ---------------------------------------------------------------------------
# Reading reduced data + metadata
# ---------------------------------------------------------------------------

@mcp.tool()
def pyirena_read_reduced_data(path: str, decimate: int = 500,
                              include_full: bool = False) -> dict:
    """Read the raw reduced I(Q) curve from a SAS file.

    Arrays are decimated to ~*decimate* points by default to keep the
    response compact. Set include_full=True for full fidelity (avoid in
    LLM workflows — long arrays bloat context).
    """
    return papi.read_reduced_data(path=path, decimate=decimate,
                                   include_full=include_full)


@mcp.tool()
def pyirena_read_metadata(path: str) -> dict:
    """Read sample / experiment metadata from a SAS file."""
    return papi.read_metadata(path)


# ---------------------------------------------------------------------------
# Per-tool results
# ---------------------------------------------------------------------------

@mcp.tool()
def pyirena_read_simple_fit(path: str, include_arrays: bool = False,
                            max_points: int = 500) -> dict:
    """Read Simple Fits results (Guinier, Porod, etc.).

    Arrays (Q, I_model, residuals) are omitted by default. Set
    include_arrays=True to include them (decimated to max_points).
    """
    return papi.read_simple_fit(path=path, include_arrays=include_arrays,
                                 max_points=max_points)


@mcp.tool()
def pyirena_read_unified_fit(path: str, include_arrays: bool = False,
                             max_points: int = 500) -> dict:
    """Read Unified Fit (Beaucage) results — multi-level Rg/G/B/P + correlations."""
    return papi.read_unified_fit(path=path, include_arrays=include_arrays,
                                  max_points=max_points)


@mcp.tool()
def pyirena_read_size_distribution(path: str, include_arrays: bool = False,
                                   max_points: int = 500) -> dict:
    """Read Size Distribution fit results — Vf, Rg, r_grid, distribution."""
    return papi.read_size_distribution(path=path, include_arrays=include_arrays,
                                        max_points=max_points)


@mcp.tool()
def pyirena_read_modeling(path: str, include_arrays: bool = False,
                          max_points: int = 500) -> dict:
    """Read parametric Modeling results (size_dist / unified_level / diff_peak / fractal pops)."""
    return papi.read_modeling(path=path, include_arrays=include_arrays,
                               max_points=max_points)


@mcp.tool()
def pyirena_read_saxs_morph(path: str, include_arrays: bool = False,
                            max_points: int = 500) -> dict:
    """Read SAXS Morph results — voxelgram-based forward modeling output.

    Note: the 3-D voxelgram itself is intentionally not returned; only
    derived scalar parameters and 1-D curves.
    """
    return papi.read_saxs_morph(path=path, include_arrays=include_arrays,
                                 max_points=max_points)


@mcp.tool()
def pyirena_read_waxs_peakfit(path: str, include_arrays: bool = False,
                              max_points: int = 500) -> dict:
    """Read WAXS Peak Fit results — per-peak Q0, FWHM, A, eta, area."""
    return papi.read_waxs_peakfit(path=path, include_arrays=include_arrays,
                                   max_points=max_points)


@mcp.tool()
def pyirena_read_fractals(path: str) -> dict:
    """List fractal aggregates stored in a file (Z, df, dmin, c, Rg per aggregate)."""
    return papi.read_fractals(path)


@mcp.tool()
def pyirena_read_merge_provenance(path: str) -> dict:
    """Read Data Merge provenance: scale, q_shift, background, source files."""
    return papi.read_merge_provenance(path)


@mcp.tool()
def pyirena_read_manipulation_provenance(path: str) -> dict:
    """Read Data Manipulation provenance: operation, parameters, source file."""
    return papi.read_manipulation_provenance(path)


# ---------------------------------------------------------------------------
# Cross-file aggregation
# ---------------------------------------------------------------------------

@mcp.tool()
def pyirena_tabulate_parameter(
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
    tool's schema. Parameter names are case-insensitive. Common examples:
      - size_distribution: 'rg', 'volume_fraction', 'chi_squared', 'q_power'
      - unified_fit:       'Rg', 'G', 'B', 'P' (per-level; use subgroup_index)
      - modeling:          'pop_Rg', 'background', 'chi_squared'
      - waxs_peakfit:      'peak_Q0', 'peak_FWHM', 'peak_A' (per-peak)
    For per-subgroup parameters supply subgroup_index (1-based, default 1).
    """
    return papi.tabulate_parameter(
        folder=folder, tool=tool, parameter=parameter, x_axis=x_axis,
        subgroup_index=subgroup_index, sample_filter=sample_filter,
    )


@mcp.tool()
def pyirena_summarize_sample(folder: str, sample: str) -> dict:
    """Condense everything known about one sample across a folder.

    File list, per-analysis file count, and min/max/n of every top-level
    scalar parameter the file's analysis tools declare.
    """
    return papi.summarize_sample(folder=folder, sample=sample)


# ---------------------------------------------------------------------------
# Plotting — returns inline images
# ---------------------------------------------------------------------------

def _plot_result_as_mcp_content(result: dict) -> list[Any]:
    """Return plot output as mixed MCP content for clients such as AnythingLLM.

    AnythingLLM's MCP image renderer looks for image-bearing items in the
    MCP tool result content array. Returning a short text item followed by a
    FastMCP Image keeps the saved file path visible while still allowing
    FastMCP to serialize the PNG as:
        {"type": "image", "data": "...base64...", "mimeType": "image/png"}
    """
    path = Path(result["path"])
    return [
        f"Plot saved to: {path}",
        Image(data=path.read_bytes(), format="png"),
    ]


@mcp.tool()
def pyirena_plot_iq(
    paths: list[str],
    overlay: bool = True,
    log_x: bool = True,
    log_y: bool = True,
    output_path: Optional[str] = None,
) -> list[Any]:
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
    return _plot_result_as_mcp_content(result)


@mcp.tool()
def pyirena_plot_parameter_trend(
    folder: str,
    tool: str,
    parameter: str,
    x_axis: str = "scan_number",
    subgroup_index: Optional[int] = None,
    sample_filter: Optional[str] = None,
    output_path: Optional[str] = None,
) -> list[Any]:
    """Plot a parameter trend across many files; returns the PNG inline.

    See pyirena_tabulate_parameter() for the *tool* / *parameter* /
    *subgroup_index* semantics. Useful for time-series questions like
    "how is Rg evolving across the latest 30 scans?"
    """
    result = papi.plot_parameter_trend(
        folder=folder, tool=tool, parameter=parameter, x_axis=x_axis,
        subgroup_index=subgroup_index, sample_filter=sample_filter,
        output_path=output_path, return_base64=False,
    )
    return _plot_result_as_mcp_content(result)


# ---------------------------------------------------------------------------
# Control API — AI-driven fitting (stateful, session-based)
# ---------------------------------------------------------------------------
#
# These tools expose pyirena.api.control, which lets an AI agent drive the
# Unified Fit model end-to-end: open data → select model → configure params
# → run fit → evaluate → save.
#
# Workflow:
#   1. pyirena_ctrl_open_dataset()       → returns session_id
#   2. pyirena_ctrl_select_model()       → choose model + nlevels
#   3. pyirena_ctrl_fix_all_except()     → staged fitting setup
#   4. pyirena_ctrl_run_fit()            → run
#   5. pyirena_ctrl_get_fit_image()      → inspect visually
#   6. pyirena_ctrl_save_fit()           → persist to HDF5
#
# Sessions live in-memory for the lifetime of this server process.
# All tool names are prefixed pyirena_ctrl_ for global uniqueness.

import base64 as _base64

from pyirena.api import control as _ctrl


def _ctrl_image_result(result: dict, session_id: str) -> list[Any]:
    """Convert a get_fit_image / get_residuals_image result to MCP content."""
    if "error" in result:
        return [result]
    b64 = result.get("image_base64", "")
    label = (
        f"Fit image (session {session_id})"
        + (" — includes residuals subplot" if result.get("has_residuals") else " — pre-fit preview")
    )
    return [label, Image(data=_base64.b64decode(b64), format="png")]


# --- Session lifecycle ---

@mcp.tool()
def pyirena_ctrl_open_dataset(file_path: str) -> dict:
    """Load a NXcanSAS HDF5 file and open a fitting session.

    Returns a session_id you must pass to every other pyirena_ctrl_ tool.
    Also returns a data summary (Q range, n_points, intensity range).
    """
    return _ctrl.open_dataset(file_path)


@mcp.tool()
def pyirena_ctrl_list_open_sessions() -> dict:
    """List all currently open fitting sessions (session_id, file, model, fit status)."""
    return _ctrl.list_open_sessions()


@mcp.tool()
def pyirena_ctrl_close_session(session_id: str) -> dict:
    """Close and discard a fitting session, freeing its memory."""
    return _ctrl.close_session(session_id)


@mcp.tool()
def pyirena_ctrl_get_session_summary(session_id: str) -> dict:
    """Return a summary of the session: file, model, Q range, fit status, χ²."""
    return _ctrl.get_session_summary(session_id)


# --- Model selection ---

@mcp.tool()
def pyirena_ctrl_list_available_models() -> dict:
    """List the fitting models available in this version of the control API."""
    return _ctrl.list_available_models()


@mcp.tool()
def pyirena_ctrl_select_model(
    session_id: str,
    model_name: str = "unified_fit",
    nlevels: int = 1,
) -> dict:
    """Select and initialise a fitting model.

    For unified_fit, nlevels sets the number of structural levels (1–5).
    Start with 1; add levels with pyirena_ctrl_add_unified_level() if residuals
    show systematic structure after the initial fit.
    Returns the full initial parameter table.
    """
    return _ctrl.select_model(session_id, model_name=model_name, nlevels=nlevels)


@mcp.tool()
def pyirena_ctrl_get_model_parameters(session_id: str) -> dict:
    """Return the current parameter table.

    Each entry has: name, value, fixed, lo, hi, units, description.
    Level parameters are named <param>_<level>, e.g. Rg_1, G_2.
    """
    return _ctrl.get_model_parameters(session_id)


@mcp.tool()
def pyirena_ctrl_get_model_description(session_id: str) -> dict:
    """Return a text description of the model — physical meaning and fitting tips.

    Read this before starting a fit to understand what each parameter controls.
    """
    return _ctrl.get_model_description(session_id)


# --- Parameter control ---

@mcp.tool()
def pyirena_ctrl_set_parameter_value(
    session_id: str, param_name: str, value: float
) -> dict:
    """Set the starting/current value of a named parameter (e.g. 'Rg_1', 'background').

    Takes effect at the next pyirena_ctrl_run_fit() call.
    """
    return _ctrl.set_parameter_value(session_id, param_name, value)


@mcp.tool()
def pyirena_ctrl_set_parameter_bounds(
    session_id: str, param_name: str, lo: float, hi: float
) -> dict:
    """Set lower and upper bounds for a parameter during fitting."""
    return _ctrl.set_parameter_bounds(session_id, param_name, lo, hi)


@mcp.tool()
def pyirena_ctrl_fix_parameter(session_id: str, param_name: str) -> dict:
    """Hold a parameter fixed at its current value during fitting."""
    return _ctrl.fix_parameter(session_id, param_name)


@mcp.tool()
def pyirena_ctrl_free_parameter(session_id: str, param_name: str) -> dict:
    """Allow a parameter to vary during fitting."""
    return _ctrl.free_parameter(session_id, param_name)


@mcp.tool()
def pyirena_ctrl_fix_all_except(
    session_id: str, free_list: list[str]
) -> dict:
    """Fix every parameter except those in free_list.

    Core tool for staged fitting: fix everything, then free parameters
    in groups.  Example: fix_all_except(['Rg_1', 'G_1', 'background']).
    """
    return _ctrl.fix_all_except(session_id, free_list)


@mcp.tool()
def pyirena_ctrl_reset_parameters_to_defaults(session_id: str) -> dict:
    """Reset all parameters to factory defaults."""
    return _ctrl.reset_parameters_to_defaults(session_id)


# --- Unified Fit level management ---

@mcp.tool()
def pyirena_ctrl_add_unified_level(
    session_id: str, position: int = -1
) -> dict:
    """Add a structural level to the Unified Fit model.

    position=-1 appends; position=1 inserts before the first level.
    Returns the updated parameter table with new level names.
    """
    return _ctrl.add_unified_level(session_id, position=position)


@mcp.tool()
def pyirena_ctrl_remove_unified_level(session_id: str, level: int) -> dict:
    """Remove a structural level (1-based) from the Unified Fit model.

    Levels are renumbered from 1 after removal.
    """
    return _ctrl.remove_unified_level(session_id, level)


# --- Q range ---

@mcp.tool()
def pyirena_ctrl_get_data_q_range(session_id: str) -> dict:
    """Return the full Q range of the loaded dataset."""
    return _ctrl.get_data_q_range(session_id)


@mcp.tool()
def pyirena_ctrl_get_fit_q_range(session_id: str) -> dict:
    """Return the Q range currently used for fitting."""
    return _ctrl.get_fit_q_range(session_id)


@mcp.tool()
def pyirena_ctrl_set_fit_q_range(
    session_id: str,
    q_min: Optional[float] = None,
    q_max: Optional[float] = None,
) -> dict:
    """Restrict the Q range used for fitting.

    Useful to exclude low-Q beam stop artefacts or high-Q noise.
    Pass q_min and/or q_max; omit either to leave that end unchanged.
    """
    return _ctrl.set_fit_q_range(session_id, q_min=q_min, q_max=q_max)


@mcp.tool()
def pyirena_ctrl_reset_fit_q_range(session_id: str) -> dict:
    """Restore the fit Q range to the full data Q range."""
    return _ctrl.reset_fit_q_range(session_id)


# --- Fit execution ---

@mcp.tool()
def pyirena_ctrl_run_fit(
    session_id: str,
    max_iter: Optional[int] = None,
) -> dict:
    """Run the fitting algorithm.

    Uses current parameter values as starting point and the current Q range.
    Re-running after a partial fit continues from where it left off.
    Returns chi_squared, reduced_chi_squared, iterations, and updated parameter values.
    Good fit: reduced_chi_squared close to 1.0.
    """
    return _ctrl.run_fit(session_id, max_iter=max_iter)


# --- Quality assessment ---

@mcp.tool()
def pyirena_ctrl_get_chi_squared(session_id: str) -> dict:
    """Return χ² and reduced χ² from the last fit."""
    return _ctrl.get_chi_squared(session_id)


@mcp.tool()
def pyirena_ctrl_get_residuals(session_id: str) -> dict:
    """Return normalised residuals and summary statistics from the last fit.

    rms close to 1.0 means the fit matches the data within error bars.
    Systematic patterns in residuals suggest the model needs adjustment.
    """
    return _ctrl.get_residuals(session_id)


@mcp.tool()
def pyirena_ctrl_get_fit_image(
    session_id: str,
    width: int = 1024,
    height: int = 768,
) -> list[Any]:
    """Capture the current fit as an inline PNG image.

    Works before a fit (shows data + model at current parameter values)
    and after a fit (adds a residuals subplot).
    Examine the image to assess fit quality before deciding next steps.
    """
    result = _ctrl.get_fit_image(session_id, width=width, height=height)
    return _ctrl_image_result(result, session_id)


@mcp.tool()
def pyirena_ctrl_get_residuals_image(
    session_id: str,
    width: int = 1024,
    height: int = 768,
) -> list[Any]:
    """Capture fit + residuals panel as an inline PNG. Requires a completed fit."""
    result = _ctrl.get_residuals_image(session_id, width=width, height=height)
    return _ctrl_image_result(result, session_id)


# --- Persistence ---

@mcp.tool()
def pyirena_ctrl_save_fit(
    session_id: str, output_path: Optional[str] = None
) -> dict:
    """Save the fitted result to NXcanSAS HDF5.

    Defaults to overwriting the original file. Pass output_path to save
    to a different location and preserve the original.
    """
    return _ctrl.save_fit(session_id, output_path=output_path)


@mcp.tool()
def pyirena_ctrl_export_fit_report(
    session_id: str, format: str = "markdown"
) -> dict:
    """Export a human-readable fit report.

    format='markdown' (default) for display; format='json' for machine use.
    Returns the report text in the 'content' key.
    """
    return _ctrl.export_fit_report(session_id, format=format)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Run the MCP server over stdio (default transport)."""
    mcp.run()


if __name__ == "__main__":
    main()
