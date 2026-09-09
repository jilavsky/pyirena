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
from typing import Any, Optional

try:
    from mcp.server.fastmcp import FastMCP, Image
except ImportError as exc:  # pragma: no cover - import guard
    # Three different failures land here and they need three different fixes.
    # mcp 2.0 renamed ``FastMCP`` to ``MCPServer`` and replaced
    # ``mcp.server.fastmcp`` with a stub that raises ModuleNotFoundError --
    # a subclass of ImportError -- so a single "install pyirena[mcp]" message
    # would tell users to install a package they already have.
    try:
        from importlib.metadata import version as _pkg_version

        _installed_mcp = _pkg_version("mcp")
    except Exception:  # pragma: no cover - distribution metadata unavailable
        _installed_mcp = None

    _mcp_major = -1
    if _installed_mcp:
        try:
            _mcp_major = int(_installed_mcp.split(".", 1)[0])
        except ValueError:  # pragma: no cover - unparseable version string
            _mcp_major = -1

    if _installed_mcp is None:
        _hint = (
            "The 'mcp' package is required to run pyirena-mcp. "
            "Install with: pip install pyirena[mcp]"
        )
    elif _mcp_major >= 2:
        _hint = (
            f"pyirena-mcp requires the mcp 1.x SDK, but mcp {_installed_mcp} is "
            "installed. mcp 2.0 renamed FastMCP to MCPServer and removed "
            "'mcp.server.fastmcp', which pyirena/mcp/server.py is built on. "
            "Fix with:  pip install 'mcp>=1.0.0,<2.0'  -- or upgrade pyirena "
            "itself (pip install -U 'pyirena[mcp]'), which pins mcp<2 as of "
            "1.1.0b7. pyirena has not been migrated to the mcp 2.x API yet."
        )
    else:
        _hint = (
            f"Could not import 'mcp.server.fastmcp' from mcp {_installed_mcp}: {exc}"
        )
    raise ImportError(_hint) from exc

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
        "CONTROL tools: drive fitting interactively. Session lifecycle is "
        "top-level: pyirena_ctrl_open_dataset() → session_id, "
        "pyirena_ctrl_list_open_sessions(), pyirena_ctrl_close_session(), "
        "pyirena_ctrl_get_session_summary(). Everything else — model "
        "selection, parameters, fit execution, quality, persistence — goes "
        "through a small dispatcher instead of one MCP tool per function: "
        "pyirena_list_categories() to see what's available, "
        "pyirena_list_tools(category) to list names in one, "
        "pyirena_describe_tool(name) for a tool's full argument schema, and "
        "pyirena_call(name, arguments) to run it. This keeps the tool count "
        "small regardless of how many underlying functions exist. "
        "\n\n"
        "Five models are available — Unified Fit, Size Distribution (Sizes), "
        "Simple Fits, Modeling and WAXS Peak Fit — matching the dispatcher "
        "categories 'unified', 'sizes', 'simple', 'modeling', 'waxs'. Prefer "
        "Simple Fits when the question is about one feature over a "
        "restricted Q range (an Rg, a Porod slope, the invariant); Unified "
        "Fit for a whole multi-level curve; Sizes to invert a dilute single "
        "population to a size histogram; Modeling when the curve needs "
        "several components at once or a specific form factor (core-shell, "
        "cylinder); WAXS Peak Fit for wide-angle patterns where the "
        "questions are peak position, width and integrated area. "
        "Unified Fit workflow: pyirena_ctrl_open_dataset() → session_id → "
        "pyirena_call('select_model', {...}) → "
        "pyirena_call('fix_all_except', {...}) → "
        "pyirena_call('run_fit', {...}) → pyirena_call('get_fit_image', {...}) → "
        "pyirena_call('save_fit', {...}). "
        "Sizes workflow (category 'sizes'): pyirena_ctrl_open_dataset() → "
        "pyirena_call('suggest_sizes_setup', {...}) → "
        "pyirena_call('select_sizes_model', {...}) → "
        "set_shape / set_size_grid / set_error_handling → "
        "fit_power_law_background + fit_flat_background → "
        "pyirena_call('set_fit_q_range', {...}) (inversion window) → "
        "pyirena_call('run_sizes_fit', {...}) → "
        "pyirena_call('get_sizes_fit_image', {...}) → "
        "pyirena_call('save_sizes_fit', {...}). "
        "Simple Fits workflow (category 'simple'): pyirena_ctrl_open_dataset() → "
        "pyirena_call('list_simple_models', {}) → "
        "pyirena_call('select_simple_model', {...}) → "
        "pyirena_call('set_fit_q_range', {...}) → "
        "pyirena_call('run_simple_fit', {...}) → "
        "pyirena_call('get_simple_linearization_image', {...}) (validity check) → "
        "pyirena_call('save_simple_fit', {...}). "
        "Modeling workflow (category 'modeling'): pyirena_ctrl_open_dataset() → "
        "pyirena_call('select_modeling_model', {...}) → "
        "pyirena_call('list_population_types', {}) → "
        "pyirena_call('add_population', {...}) → set_population_option / "
        "set_population_parameter / set_population_parameter_fit → "
        "pyirena_call('set_modeling_q_range', {...}) → "
        "pyirena_call('run_modeling_fit', {...}) → "
        "pyirena_call('get_modeling_fit_image', {...}) → "
        "pyirena_call('save_modeling_fit', {...}). Modeling parameters use "
        "dotted names (dist.mean_size, ff.sld_core, sf.eta) — always list "
        "them with get_population_parameters rather than guessing. "
        "WAXS workflow (category 'waxs'): pyirena_ctrl_open_dataset() → "
        "pyirena_call('select_waxs_model', {...}) → "
        "pyirena_call('find_waxs_peaks', {...}) (data-driven starting positions) → "
        "pyirena_call('run_waxs_fit', {...}) → "
        "pyirena_call('get_waxs_results', {...}) (positions, widths, areas) → "
        "pyirena_call('save_waxs_fit', {...}). "
        "The session tools and Q-range tools (set_fit_q_range etc, category "
        "'unified') are shared between all five tools (Modeling has its own "
        "set_modeling_q_range). "
        "Sessions are in-memory for this server process."
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
# Images — how pyirena returns pictures over MCP
# ---------------------------------------------------------------------------
#
# Every image tool returns a two-item content list:
#
#   1. a text block naming the PNG's absolute path on disk
#   2. an MCP image content block  {"type": "image", "data": "<base64>",
#      "mimeType": "image/png"}
#
# That pair is the portable answer.  The image block is what Claude Desktop /
# claude.ai, ChatGPT and other MCP hosts render inline; the text block is the
# fallback for clients that drop image blocks (AnythingLLM in some versions,
# terminal clients, plain LLM agents), which can then open or hand off the
# file instead of being told the picture does not exist.
#
# Registration matters as much as the payload.  FastMCP >= 1.10 derives a
# *structured output* JSON schema from a tool's return annotation and then
# serialises the return value against it.  ``Image`` is not JSON-serialisable,
# so an image tool annotated ``-> list[Any]`` fails with
#     Unable to serialize unknown type: <class '...fastmcp.utilities.types.Image'>
# *after* the PNG has been written — the client sees only an error.  Images
# belong in the unstructured content array, so image tools opt out of
# structured output via ``_image_tool``.

def _image_tool(**kwargs: Any):
    """``@mcp.tool()`` for tools that return image content blocks.

    Disables FastMCP's structured-output schema (mcp >= 1.10). On older mcp
    releases there is no structured output and the plain decorator is already
    correct, so the ``TypeError`` fallback keeps those working.
    """
    try:
        return mcp.tool(structured_output=False, **kwargs)
    except TypeError:  # pragma: no cover - mcp < 1.10
        return mcp.tool(**kwargs)


def _image_content(
    label: str,
    png_bytes: bytes,
    path: Optional[str] = None,
) -> list[Any]:
    """Build the standard [text, image] content pair."""
    text = f"{label}\nPNG saved to: {path}" if path else label
    return [text, Image(data=png_bytes, format="png")]


def _plot_result_as_mcp_content(result: dict, label: str) -> list[Any]:
    """Turn a pyirena.api.plotting result into MCP content."""
    path = Path(result["path"])
    return _image_content(label, path.read_bytes(), str(path))


@_image_tool()
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
    return _plot_result_as_mcp_content(result, "I(Q) plot")


@_image_tool()
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
    return _plot_result_as_mcp_content(
        result, f"{parameter} trend across {tool} results"
    )


# ---------------------------------------------------------------------------
# Control API — AI-driven fitting (stateful, session-based)
# ---------------------------------------------------------------------------
#
# These tools expose pyirena.api.control, which lets an AI agent drive any
# of the five fitting tools end-to-end: open data → select model →
# configure params → run fit → evaluate → save. Session lifecycle
# (open_dataset / list_open_sessions / close_session / get_session_summary)
# is always top-level below; everything else is reached through the
# pyirena_call() dispatcher further down (see the comment above it).
#
# Unified Fit workflow, for example:
#   1. pyirena_ctrl_open_dataset()                              → session_id
#   2. pyirena_call("select_model", {...})                      → choose model + nlevels
#   3. pyirena_call("fix_all_except", {...})                    → staged fitting setup
#   4. pyirena_call("run_fit", {...})                           → run
#   5. pyirena_call("get_fit_image", {...})                     → inspect visually
#   6. pyirena_call("save_fit", {...})                          → persist to HDF5
#
# Sessions live in-memory for the lifetime of this server process.

import base64 as _base64

from pyirena.api import control as _ctrl


def _ctrl_image_result(result: dict, label: str) -> list[Any]:
    """Convert a control-API image result to MCP content.

    Control-API image functions return ``{"image_base64": ..., "image_path":
    ...}`` (or an error dict).  Both halves are passed through: the base64
    payload becomes the inline image block, the path becomes the text block.
    """
    if not isinstance(result, dict) or "error" in result:
        return [result]
    b64 = result.get("image_base64") or ""
    if not b64:
        return [{
            "error": "No image was produced.",
            "code": "NO_IMAGE",
            "result": result,
        }]
    return _image_content(label, _base64.b64decode(b64), result.get("image_path"))


def _dispatch_image_label(name: str, arguments: Optional[dict[str, Any]], result: dict) -> str:
    """Label for a pyirena_call() image result, noting the tool, session and residuals."""
    session_id = (arguments or {}).get("session_id", "?")
    suffix = (
        " — includes residuals subplot"
        if isinstance(result, dict) and result.get("has_residuals")
        else ""
    )
    return f"{name} (session {session_id}){suffix}"


# --- Session lifecycle ---

@mcp.tool()
def pyirena_ctrl_open_dataset(file_path: str, use_slit_smeared: bool = False) -> dict:
    """Load a NXcanSAS HDF5 file and open a fitting session.

    Returns a session_id you must pass to every other pyirena_ctrl_ tool.
    Also returns a data summary (Q range, n_points, intensity range) plus
    slit-smearing status: ``is_slit_smeared``, ``slit_length`` (1/Å), and
    ``has_slit_smeared_entry`` (True when the file also carries a slit-smeared
    ``_SMR`` dataset alongside the desmeared one).

    Set ``use_slit_smeared=True`` to load the file's slit-smeared dataset; the
    Unified Fit model then smears automatically at the file's slit length, and
    the fitted parameters are ideal-space (pinhole-equivalent).
    """
    return _ctrl.open_dataset(file_path, use_slit_smeared=use_slit_smeared)


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


# ---------------------------------------------------------------------------
# Control API — dispatcher (everything except session lifecycle)
# ---------------------------------------------------------------------------
#
# ~90 pyirena.api.control functions across Unified Fit, Sizes, Simple Fits,
# Modeling and WAXS Peak Fit share one JSON-schema registry
# (pyirena.api.control.schemas.TOOL_SCHEMA_BY_NAME) but used to each get
# their own MCP tool. Combined with whatever else is active in a client's
# session, that could exceed a provider-side cap on the number of tools in
# a single request (observed via the ANL Argo gateway proxy: 128 tools).
# pyirena.mcp.dispatch collapses them into the four tools below.
# Session-lifecycle tools stay top-level above since nearly every workflow
# starts there.
#
# Discovery: pyirena_list_categories() -> pyirena_list_tools(category) ->
# pyirena_describe_tool(name) -> pyirena_call(name, arguments).
#
# Categories, one per pyirena.api.control submodule: "unified" (model
# selection, parameter control, level management, Q range, local
# estimators, fit execution, quality, persistence for Unified Fit -- the Q
# range and level tools here are shared by every model, not Unified-Fit-only),
# "sizes", "simple", "modeling", "waxs". See pyirena/mcp/dispatch.py.

from pyirena.mcp import dispatch as _dispatch


@mcp.tool()
def pyirena_list_categories() -> dict:
    """List the pyirena_call() tool categories, with a one-line description
    and tool count each. Start here to discover what pyirena_call can do."""
    return _dispatch.list_categories()


@mcp.tool()
def pyirena_list_tools(category: str) -> dict:
    """List the tool names and one-line summaries in a pyirena_call() category.

    category is one of the names returned by pyirena_list_categories().
    """
    return _dispatch.list_tools(category)


@mcp.tool()
def pyirena_describe_tool(name: str) -> dict:
    """Return the full JSON schema (parameters, types, description) for one
    pyirena_call()-dispatched tool name."""
    return _dispatch.describe_tool(name)


@_image_tool()
def pyirena_call(name: str, arguments: Optional[dict[str, Any]] = None) -> Any:
    """Call one of the Unified Fit / Sizes / Simple Fits / Modeling / WAXS
    Peak Fit control tools by name (everything except session lifecycle).

    name must be a tool name from pyirena_list_tools(category); see
    pyirena_describe_tool(name) for its expected arguments. Session-lifecycle
    tools (open_dataset, list_open_sessions, close_session,
    get_session_summary) are NOT dispatched here -- call
    pyirena_ctrl_open_dataset() etc directly. Some tools return an inline
    PNG image instead of plain data; both cases are handled transparently.
    """
    result = _dispatch.call_tool(name, arguments)
    if isinstance(result, dict) and "image_base64" in result:
        return _ctrl_image_result(result, _dispatch_image_label(name, arguments, result))
    return result


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Run the MCP server over stdio (default transport)."""
    from pyirena.logging_setup import install_excepthook, setup_logging
    setup_logging("mcp")   # console handler writes to stderr; stdout stays clean for MCP
    install_excepthook()
    mcp.run()


if __name__ == "__main__":
    main()
