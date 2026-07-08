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
        "Two models are available — Unified Fit and Size Distribution (Sizes). "
        "Unified Fit workflow: pyirena_ctrl_open_dataset() → session_id → "
        "pyirena_ctrl_select_model() → pyirena_ctrl_fix_all_except() → "
        "pyirena_ctrl_run_fit() → pyirena_ctrl_get_fit_image() → "
        "pyirena_ctrl_save_fit(). "
        "Sizes workflow (pyirena_ctrl_sizes_ prefix): pyirena_ctrl_open_dataset() → "
        "pyirena_ctrl_sizes_suggest_setup() → pyirena_ctrl_sizes_select_model() → "
        "set_shape / set_size_grid / set_error_handling → "
        "fit_power_law_background + fit_flat_background → "
        "pyirena_ctrl_set_fit_q_range() (inversion window) → "
        "pyirena_ctrl_sizes_run_fit() → pyirena_ctrl_sizes_get_fit_image() → "
        "pyirena_ctrl_sizes_save_fit(). The session and Q-range tools are shared "
        "between both models. Sessions are in-memory for this server process."
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


# --- Per-level boolean options (correlations, mass_fractal, link_B, link_RGCO) ---

@mcp.tool()
def pyirena_ctrl_get_level_options(
    session_id: str, level: Optional[int] = None
) -> dict:
    """Return per-level boolean flag state.

    These flags (correlations, mass_fractal, link_B, link_RGCO) switch entire
    features of the intensity formula on or off and are SEPARATE from numeric
    parameters.  Omit `level` to get the state of all levels.
    """
    return _ctrl.get_level_options(session_id, level)


@mcp.tool()
def pyirena_ctrl_set_level_option(
    session_id: str, level: int, option: str, enabled: bool
) -> dict:
    """Toggle a per-level boolean flag.

    CRITICAL: setting numeric parameters ETA and PACK has NO EFFECT unless
    the level's `correlations` option is True — use this tool to enable it.

    option must be one of:
      - "correlations" — Born-Green liquid-like-ordering (uses ETA + PACK)
      - "mass_fractal" — auto-compute B from G, Rg, P; manual B is ignored
      - "link_B" — estimate B from G, Rg, P via the Porod invariant
      - "link_RGCO" — link RgCO to previous level's Rg
    """
    return _ctrl.set_level_option(session_id, level, option, enabled)


@mcp.tool()
def pyirena_ctrl_check_level_feasibility(
    session_id: str, level: Optional[int] = None
) -> dict:
    """Check whether each level's parameters are physically meaningful.

    A level is "feasible" when its Guinier and power-law regions connect
    smoothly at the Hammouda rollover Q point.  Use AFTER run_fit to catch
    combinations that converged mathematically but are not physically
    interpretable (e.g. B too small/large for the chosen P, or G inconsistent
    with B and Rg).  Omit `level` to check every level.
    """
    return _ctrl.check_level_feasibility(session_id, level)


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


# --- Local one-term estimators (good for starting values) ---

@mcp.tool()
def pyirena_ctrl_fit_local_guinier(
    session_id: str, q_min: float, q_max: float
) -> dict:
    """Fit Guinier I(q) = G·exp(-q²·Rg²/3) on a Q sub-range; return G and Rg.

    Equivalent to the GUI's 'Fit Rg/G btwn cursors' button.  Useful for
    estimating a level's starting Rg and G *before* running the full
    multi-level fit.  Does NOT modify the model — call set_parameter_value()
    afterwards if you want to apply the result.

    Pick q_min/q_max to cover the Guinier knee of the level you are
    characterising (where the log-log slope flattens).
    """
    return _ctrl.fit_local_guinier(session_id, q_min, q_max)


@mcp.tool()
def pyirena_ctrl_fit_local_power_law(
    session_id: str, q_min: float, q_max: float
) -> dict:
    """Fit power law I(q) = B·q⁻ᴾ on a Q sub-range; return P and B.

    Equivalent to the GUI's 'Fit P/B btwn cursors' button.  Useful for
    estimating a level's starting P and B from the linear portion of the
    log-log plot (typically just past the Guinier knee).  Does NOT modify
    the model — call set_parameter_value() afterwards if you want to apply
    the result.
    """
    return _ctrl.fit_local_power_law(session_id, q_min, q_max)


# --- Feature detection (slope-profile analysis) ---

@mcp.tool()
def pyirena_ctrl_detect_features(
    session_id: str,
    q_min: Optional[float] = None,
    q_max: Optional[float] = None,
    q_max_clip: Optional[float] = 0.6,
    config_overrides: Optional[dict] = None,
) -> dict:
    """Segment the loaded I(Q) curve into power-law slope segments.

    Returns a full piecewise classification of the curve in log-log space.
    Each segment has a locally-constant slope; adjacent segments with
    substantially different slopes imply Guinier knees between them
    (also returned).  Intended to help decide how many Unified Fit levels
    are needed and where to place Q-windows for local Guinier / Porod fits.

    Does NOT modify the model — purely diagnostic.

    Q clipping: by default data above q_max_clip=0.6 Å⁻¹ is dropped
    (the practical small-angle-approximation limit; amorphous diffraction
    above this should not be classified as SAS structure).  Pass None
    to disable.

    Returns dict with:
      segments — list of {q_min, q_max, P, P_std, kind, intensity_mid,
        width_decades} where P is the positive Porod exponent (I ∝ Q^-P) and
        kind is 'background' / 'guinier_plateau' / 'power_law'.
      guinier_knees — list of {q_min, q_max, q_center, P_low_q, P_high_q,
        delta_P} where P values are positive Porod exponents; P_low_q <
        P_high_q (shallower at low Q = physical Guinier knee condition).
      recommended_guinier_windows, recommended_nlevels, background_q_min,
      n_segments_found, log_decades, n_points, q_min_analysed, q_max_analysed.
    """
    return _ctrl.detect_features(
        session_id,
        q_min=q_min,
        q_max=q_max,
        q_max_clip=q_max_clip,
        config_overrides=config_overrides,
    )


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
    """Return residuals from the last fit (normalised, rescaled, and fractional).

    'residuals' = normalised (I-M)/sigma; rms close to 1.0 means the fit matches
    the data within error bars. Also returns 'rescaled_residual' (r/robust_scale_s)
    and 'frac_misfit_percent' ((I-M)/I in %, sigma-independent), plus
    summary.robust_scale_s. Systematic patterns in residuals suggest the model
    needs adjustment. For the full diagnostic set use pyirena_ctrl_get_fit_quality.
    """
    return _ctrl.get_residuals(session_id)


@mcp.tool()
def pyirena_ctrl_get_fit_quality(session_id: str, n_bands: int = 4) -> dict:
    """Robust, sigma-scale-independent fit-quality diagnostics for the last fit.

    Preferred over reduced chi-squared alone when reported uncertainties sigma may
    be mis-scaled (common in SAXS), where chasing reduced chi-squared ~ 1 is
    misleading. Key fields:
      - robust_scale_s: how many times the actual scatter exceeds reported sigma
        (~1 sigma honest; ~3 sigma ~3x too small, so realistic_reduced_chi2_floor ~9).
      - max_abs_frac_misfit (+ q_at_max_frac_misfit): largest |(I-M)/I|, a
        sigma-independent gross-misfit backstop (>~0.3 is a real local misfit).
      - n_outliers_3s: points beyond 3*robust_scale_s.
      - longest_same_sign_run / sign_autocorr_lag1: structure signalling a wrong
        functional form, distinct from a pure sigma-scale problem.
      - bands: the same metrics per Q-decade (uneven per-band chi2 is a misfit signal).

    Returns facts only — interpret thresholds yourself.
    """
    return _ctrl.get_fit_quality(session_id, n_bands=n_bands)


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
# Control API — Size Distribution (Sizes) fitting
# ---------------------------------------------------------------------------
#
# Sizes counterpart of the Unified Fit control tools above.  Reuses the shared
# session lifecycle (pyirena_ctrl_open_dataset / list_open_sessions /
# close_session / get_session_summary) and the Q-range tools
# (pyirena_ctrl_get_data_q_range / get_fit_q_range / set_fit_q_range /
# reset_fit_q_range) — set_fit_q_range defines the INVERSION Q-range.
#
# Workflow:
#   1. pyirena_ctrl_open_dataset()                 → session_id
#   2. pyirena_ctrl_sizes_suggest_setup()          → data-driven recommendations
#   3. pyirena_ctrl_sizes_select_model()           → choose inversion method
#   4. pyirena_ctrl_sizes_set_shape() / set_size_grid() / set_error_handling()
#   5. pyirena_ctrl_sizes_fit_power_law_background() + fit_flat_background()
#   6. pyirena_ctrl_set_fit_q_range()              → inversion window (shared tool)
#   7. pyirena_ctrl_sizes_run_fit()                → run inversion
#   8. pyirena_ctrl_sizes_get_fit_image()          → inspect visually
#   9. pyirena_ctrl_sizes_save_fit()               → persist to HDF5


@mcp.tool()
def pyirena_ctrl_sizes_select_model(session_id: str, method: str = "maxent") -> dict:
    """Create a Size Distribution model for the session.

    method: 'maxent' (recommended default), 'regularization', 'tnnls', or
    'montecarlo'. Best for dilute samples with a single particle population.
    Replaces any existing model and clears prior fit results.
    """
    return _ctrl.select_sizes_model(session_id, method=method)


@mcp.tool()
def pyirena_ctrl_sizes_get_config(session_id: str) -> dict:
    """Return the current Sizes configuration (grid, shape, method, error
    handling, complex background)."""
    return _ctrl.get_sizes_config(session_id)


@mcp.tool()
def pyirena_ctrl_sizes_suggest_setup(session_id: str) -> dict:
    """Inspect the data and recommend a Sizes setup, with suitability warnings.

    Returns a 'suitable' flag, 'recommended' r-range / inversion Q-range /
    background windows, and 'warnings' (e.g. no size scale, multiple
    populations). Advisory only — apply values with the set_* / fit_* tools.
    Call this before configuring the fit.
    """
    return _ctrl.suggest_sizes_setup(session_id)


@mcp.tool()
def pyirena_ctrl_sizes_set_size_grid(
    session_id: str,
    r_min: Optional[float] = None,
    r_max: Optional[float] = None,
    n_bins: Optional[int] = None,
    log_spacing: Optional[bool] = None,
) -> dict:
    """Set the radius grid [Å] for the inversion (r_min, r_max, n_bins,
    log_spacing). Heuristic: r ≈ π/Q over the inversion Q-range."""
    return _ctrl.set_size_grid(
        session_id, r_min=r_min, r_max=r_max, n_bins=n_bins, log_spacing=log_spacing
    )


@mcp.tool()
def pyirena_ctrl_sizes_set_shape(
    session_id: str,
    shape: Optional[str] = None,
    contrast: Optional[float] = None,
    aspect_ratio: Optional[float] = None,
) -> dict:
    """Set the form factor: shape ('sphere' or 'spheroid'), contrast (Δρ)² in
    10²⁰ cm⁻⁴ (use 1.0 if unknown), and aspect_ratio (spheroid only)."""
    return _ctrl.set_shape(
        session_id, shape=shape, contrast=contrast, aspect_ratio=aspect_ratio
    )


@mcp.tool()
def pyirena_ctrl_sizes_set_method(
    session_id: str,
    method: str,
    maxent_sky_background: Optional[float] = None,
    maxent_max_iter: Optional[int] = None,
    regularization_evalue: Optional[float] = None,
    regularization_min_ratio: Optional[float] = None,
    tnnls_approach_param: Optional[float] = None,
    tnnls_max_iter: Optional[int] = None,
    montecarlo_n_repetitions: Optional[int] = None,
    montecarlo_convergence: Optional[float] = None,
    montecarlo_max_iter: Optional[int] = None,
) -> dict:
    """Choose the inversion method and (optionally) its tuning parameters.

    Only parameters relevant to the chosen method are applied. MaxEnt is the
    recommended default. Method params: maxent_* / regularization_* / tnnls_* /
    montecarlo_*.
    """
    return _ctrl.set_method(
        session_id, method,
        maxent_sky_background=maxent_sky_background,
        maxent_max_iter=maxent_max_iter,
        regularization_evalue=regularization_evalue,
        regularization_min_ratio=regularization_min_ratio,
        tnnls_approach_param=tnnls_approach_param,
        tnnls_max_iter=tnnls_max_iter,
        montecarlo_n_repetitions=montecarlo_n_repetitions,
        montecarlo_convergence=montecarlo_convergence,
        montecarlo_max_iter=montecarlo_max_iter,
    )


@mcp.tool()
def pyirena_ctrl_sizes_set_error_handling(
    session_id: str,
    error_scale: Optional[float] = None,
    fractional_error: Optional[bool] = None,
    fractional_error_value: Optional[float] = None,
) -> dict:
    """Configure uncertainty handling. Either scale file errors (error_scale,
    1.0 = unchanged) or switch to fractional errors (fractional_error=True with
    fractional_error_value, e.g. 0.03 = 3%, which ignores file σ)."""
    return _ctrl.set_error_handling(
        session_id,
        error_scale=error_scale,
        fractional_error=fractional_error,
        fractional_error_value=fractional_error_value,
    )


@mcp.tool()
def pyirena_ctrl_sizes_set_background(
    session_id: str,
    power_law_B: Optional[float] = None,
    power_law_P: Optional[float] = None,
    background: Optional[float] = None,
) -> dict:
    """Set complex-background terms directly (no fitting). Background subtracted
    before inversion is power_law_B·q^(-power_law_P) + background. Set
    power_law_B=0 for a flat background only."""
    return _ctrl.set_background(
        session_id, power_law_B=power_law_B, power_law_P=power_law_P,
        background=background,
    )


@mcp.tool()
def pyirena_ctrl_sizes_fit_power_law_background(
    session_id: str,
    q_min: float,
    q_max: float,
    fit_B: bool = True,
    fit_P: bool = True,
) -> dict:
    """Fit the power-law background B·q^(-P) over [q_min, q_max] (typically the
    low-Q steep-slope region). Updates power_law_B/P. fit_B/fit_P select which
    vary (at least one True)."""
    return _ctrl.fit_power_law_background(
        session_id, q_min, q_max, fit_B=fit_B, fit_P=fit_P
    )


@mcp.tool()
def pyirena_ctrl_sizes_fit_flat_background(
    session_id: str, q_min: float, q_max: float
) -> dict:
    """Fit the flat background by averaging I − B·q^(-P) over [q_min, q_max]
    (typically the high-Q flat region). Updates background. Run after
    fit_power_law_background if both terms are present."""
    return _ctrl.fit_flat_background(session_id, q_min, q_max)


@mcp.tool()
def pyirena_ctrl_sizes_get_background_image(
    session_id: str, width: int = 1024, height: int = 768
) -> list[Any]:
    """Render the data with the current complex background overlaid (log-log).
    Use to visually confirm the background before inverting."""
    result = _ctrl.get_background_preview_image(session_id, width=width, height=height)
    if "error" in result:
        return [result]
    b64 = result.get("image_base64", "")
    return [f"Background preview (session {session_id})",
            Image(data=_base64.b64decode(b64), format="png")]


@mcp.tool()
def pyirena_ctrl_sizes_run_fit(
    session_id: str, random_seed: Optional[int] = None
) -> dict:
    """Run the size-distribution inversion. The inversion Q-range
    (pyirena_ctrl_set_fit_q_range) and complex background are applied first.

    Returns success, chi_squared, volume_fraction, rg, peak_r, n_iterations,
    n_data.
    """
    return _ctrl.run_sizes_fit(session_id, random_seed=random_seed)


@mcp.tool()
def pyirena_ctrl_sizes_get_distribution(
    session_id: str, max_points: int = 500
) -> dict:
    """Return the fitted distribution arrays: r_grid [Å] and distribution P(r)
    [vol-frac/Å] (decimated), plus distribution_std when available."""
    return _ctrl.get_sizes_distribution(session_id, max_points=max_points)


@mcp.tool()
def pyirena_ctrl_sizes_get_results(session_id: str) -> dict:
    """Return the full scalar results + configuration for the last Sizes fit
    (chi_squared, volume_fraction, rg, peak_r, plus all setup parameters)."""
    return _ctrl.get_sizes_results(session_id)


@mcp.tool()
def pyirena_ctrl_sizes_get_fit_image(
    session_id: str, width: int = 1024, height: int = 900
) -> list[Any]:
    """Render the Sizes fit as a two-panel PNG: (top) log-log data + model
    (+ background), (bottom) the size distribution P(r) vs r."""
    result = _ctrl.get_sizes_fit_image(session_id, width=width, height=height)
    if "error" in result:
        return [result]
    b64 = result.get("image_base64", "")
    return [f"Sizes fit image (session {session_id})",
            Image(data=_base64.b64decode(b64), format="png")]


@mcp.tool()
def pyirena_ctrl_sizes_save_fit(
    session_id: str, output_path: Optional[str] = None
) -> dict:
    """Save the fitted size distribution to NXcanSAS HDF5.

    Defaults to overwriting the original file. Pass output_path to save
    elsewhere and preserve the original.
    """
    return _ctrl.save_sizes_fit(session_id, output_path=output_path)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Run the MCP server over stdio (default transport)."""
    from pyirena.logging_setup import setup_logging, install_excepthook
    setup_logging("mcp")   # console handler writes to stderr; stdout stays clean for MCP
    install_excepthook()
    mcp.run()


if __name__ == "__main__":
    main()
