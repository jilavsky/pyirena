"""Result schemas for the pyirena.api surface.

Plain stdlib dataclasses. Each dataclass has a `.to_dict()` method that
returns a JSON-serializable dict — numpy arrays are decimated to a
bounded length and converted to Python lists; NaN/inf are replaced with
None so the result is strict-JSON-safe (the MCP protocol uses strict JSON).
"""
from __future__ import annotations

import dataclasses
import math
import os
from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np


# ---------------------------------------------------------------------------
# Array-handling helpers
# ---------------------------------------------------------------------------

def _default_max_points() -> int:
    """Default decimation cap; overridable via PYIRENA_MAX_ARRAY_POINTS."""
    raw = os.environ.get("PYIRENA_MAX_ARRAY_POINTS", "500")
    try:
        n = int(raw)
        return n if n > 0 else 500
    except ValueError:
        return 500


def _sanitize_scalar(v: Any) -> Any:
    """Convert numpy scalars and replace NaN/inf with None for strict JSON."""
    if v is None:
        return None
    if isinstance(v, (bool, str)):
        return v
    if isinstance(v, (np.bool_,)):
        return bool(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        f = float(v)
        return None if not math.isfinite(f) else f
    if isinstance(v, (np.ndarray,)) and v.ndim == 0:
        return _sanitize_scalar(v.item())
    if isinstance(v, bytes):
        try:
            return v.decode("utf-8")
        except UnicodeDecodeError:
            return v.decode("utf-8", errors="replace")
    if isinstance(v, dict):
        return {k: _sanitize_scalar(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_sanitize_scalar(x) for x in v]
    return v


def array_to_list(
    arr: Optional[np.ndarray | list],
    max_points: Optional[int] = None,
    include_full: bool = False,
) -> Optional[list[float | None]]:
    """Convert a numpy array to a JSON-safe list, decimating to max_points.

    Parameters
    ----------
    arr : np.ndarray or list or None
        Input array.
    max_points : int, optional
        Maximum number of points to return. Default: PYIRENA_MAX_ARRAY_POINTS
        env var (500 if unset). Ignored when include_full is True.
    include_full : bool
        If True, return the full array without decimation. Use sparingly —
        very long arrays bloat MCP responses and LLM context.
    """
    if arr is None:
        return None
    a = np.asarray(arr)
    if a.size == 0:
        return []
    if include_full:
        out = a
    else:
        cap = max_points if max_points is not None else _default_max_points()
        if a.size > cap:
            idx = np.linspace(0, a.size - 1, cap).astype(int)
            out = a[idx]
        else:
            out = a
    # Convert to list and sanitize NaN/inf -> None
    return [
        None if (isinstance(x, float) and not math.isfinite(x)) else
        (None if (isinstance(x, np.floating) and not np.isfinite(x)) else
         (float(x) if isinstance(x, (np.floating, np.integer, float, int)) else x))
        for x in out.tolist()
    ]


def asdict_clean(obj: Any) -> dict:
    """dataclasses.asdict() + scalar sanitisation."""
    if dataclasses.is_dataclass(obj):
        d = dataclasses.asdict(obj)
    elif isinstance(obj, dict):
        d = obj
    else:
        return _sanitize_scalar(obj)
    return _sanitize_scalar(d)


# ---------------------------------------------------------------------------
# Discovery schemas
# ---------------------------------------------------------------------------

@dataclass
class FileEntry:
    """One file row returned by list_files()."""
    path: str
    name: str
    sample: Optional[str] = None
    scan_number: Optional[int] = None
    mtime: Optional[str] = None       # ISO 8601
    size_bytes: int = 0
    analyses: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict_clean(self)


@dataclass
class FolderSummary:
    """Aggregate view returned by summarize_folder()."""
    folder: str
    n_files: int
    samples: list[str] = field(default_factory=list)
    analyses_count: dict[str, int] = field(default_factory=dict)
    mtime_min: Optional[str] = None
    mtime_max: Optional[str] = None
    sample_filter_applied: Optional[str] = None

    def to_dict(self) -> dict:
        return asdict_clean(self)


@dataclass
class FileInspection:
    """Single-file deep inspection returned by inspect_file()."""
    path: str
    name: str
    sample: Optional[str] = None
    scan_number: Optional[int] = None
    mtime: Optional[str] = None
    size_bytes: int = 0
    analyses_present: list[str] = field(default_factory=list)
    reduced_data_present: bool = False
    q_min: Optional[float] = None
    q_max: Optional[float] = None
    n_points: Optional[int] = None
    extra_metadata: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict_clean(self)


# ---------------------------------------------------------------------------
# Data schemas
# ---------------------------------------------------------------------------

@dataclass
class ReducedData:
    """Raw reduced I(Q) curve returned by read_reduced_data()."""
    path: str
    found: bool = False
    n_points: int = 0
    q_min: Optional[float] = None
    q_max: Optional[float] = None
    I_units: Optional[str] = None
    label: Optional[str] = None
    Q: Optional[list] = None
    I: Optional[list] = None
    dI: Optional[list] = None
    dQ: Optional[list] = None
    decimated: bool = False
    decimated_from: Optional[int] = None

    def to_dict(self) -> dict:
        return asdict_clean(self)


@dataclass
class SampleMetadata:
    """Sample / experiment metadata returned by read_metadata()."""
    path: str
    found: bool = False
    sample_name: Optional[str] = None
    label: Optional[str] = None
    thickness: Optional[float] = None
    blank: Optional[str] = None
    instrument: Optional[str] = None
    timestamp: Optional[str] = None
    extra: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict_clean(self)


# ---------------------------------------------------------------------------
# Result schemas — one per analysis tool
# ---------------------------------------------------------------------------

@dataclass
class _BaseToolResult:
    """Common fields for every per-tool result."""
    path: str
    found: bool = False
    tool: str = ""

    def to_dict(self) -> dict:
        return asdict_clean(self)


@dataclass
class SimpleFitResult(_BaseToolResult):
    tool: str = "simple_fits"
    fit_quality: Optional[dict] = None  # robust fit-quality metrics (see core.fit_metrics)
    model: Optional[str] = None
    success: Optional[bool] = None
    chi_squared: Optional[float] = None
    reduced_chi_squared: Optional[float] = None
    dof: Optional[int] = None
    q_min: Optional[float] = None
    q_max: Optional[float] = None
    n_mc_runs: Optional[int] = None
    timestamp: Optional[str] = None
    params: dict = field(default_factory=dict)
    params_std: dict = field(default_factory=dict)
    derived: dict = field(default_factory=dict)
    Q: Optional[list] = None
    I_model: Optional[list] = None
    residuals: Optional[list] = None
    intensity_data: Optional[list] = None
    intensity_error: Optional[list] = None


@dataclass
class UnifiedFitLevel:
    level_number: int
    G: Optional[float] = None
    Rg: Optional[float] = None
    B: Optional[float] = None
    P: Optional[float] = None
    RgCutoff: Optional[float] = None
    ETA: Optional[float] = None
    PACK: Optional[float] = None
    correlations: Optional[bool] = None
    G_err: Optional[float] = None
    Rg_err: Optional[float] = None
    B_err: Optional[float] = None
    P_err: Optional[float] = None
    ETA_err: Optional[float] = None
    PACK_err: Optional[float] = None


@dataclass
class UnifiedFitResult(_BaseToolResult):
    tool: str = "unified_fit"
    fit_quality: Optional[dict] = None  # robust fit-quality metrics (see core.fit_metrics)
    num_levels: Optional[int] = None
    background: Optional[float] = None
    background_err: Optional[float] = None
    chi_squared: Optional[float] = None
    reduced_chi_squared: Optional[float] = None
    timestamp: Optional[str] = None
    levels: list = field(default_factory=list)   # list of UnifiedFitLevel
    Q: Optional[list] = None
    intensity_data: Optional[list] = None
    intensity_model: Optional[list] = None
    intensity_error: Optional[list] = None
    residuals: Optional[list] = None


@dataclass
class SizeDistResult(_BaseToolResult):
    tool: str = "size_distribution"
    fit_quality: Optional[dict] = None  # robust fit-quality metrics (see core.fit_metrics)
    method: Optional[str] = None
    shape: Optional[str] = None
    chi_squared: Optional[float] = None
    reduced_chi_squared: Optional[float] = None
    volume_fraction: Optional[float] = None
    rg: Optional[float] = None
    n_iterations: Optional[int] = None
    q_power: Optional[float] = None
    contrast: Optional[float] = None
    aspect_ratio: Optional[float] = None
    r_min: Optional[float] = None
    r_max: Optional[float] = None
    n_bins: Optional[int] = None
    background: Optional[float] = None
    power_law_B: Optional[float] = None
    power_law_P: Optional[float] = None
    timestamp: Optional[str] = None
    Q: Optional[list] = None
    intensity_data: Optional[list] = None
    intensity_model: Optional[list] = None
    intensity_error: Optional[list] = None
    residuals: Optional[list] = None
    r_grid: Optional[list] = None
    distribution: Optional[list] = None
    distribution_std: Optional[list] = None
    number_dist: Optional[list] = None
    cumul_vol_dist: Optional[list] = None


@dataclass
class ModelingPopulation:
    population_index: int
    pop_type: str
    enabled: bool = True
    label: str = ""
    # Per-type parameters live in 'parameters' (size_dist uses dist_params/ff_params/sf_params)
    parameters: dict = field(default_factory=dict)
    derived: dict = field(default_factory=dict)
    # Size-dist only:
    dist_type: Optional[str] = None
    form_factor: Optional[str] = None
    structure_factor: Optional[str] = None
    radius_grid: Optional[list] = None
    volume_dist: Optional[list] = None
    number_dist: Optional[list] = None
    model_I: Optional[list] = None


@dataclass
class ModelingResult(_BaseToolResult):
    tool: str = "modeling"
    fit_quality: Optional[dict] = None  # robust fit-quality metrics (see core.fit_metrics)
    chi_squared: Optional[float] = None
    reduced_chi_squared: Optional[float] = None
    dof: Optional[int] = None
    background: Optional[float] = None
    q_min: Optional[float] = None
    q_max: Optional[float] = None
    timestamp: Optional[str] = None
    model_q: Optional[list] = None
    model_I: Optional[list] = None
    populations: list = field(default_factory=list)   # list of ModelingPopulation


@dataclass
class SAXSMorphResult(_BaseToolResult):
    tool: str = "saxs_morph"
    chi_squared: Optional[float] = None
    reduced_chi_squared: Optional[float] = None
    dof: Optional[int] = None
    volume_fraction: Optional[float] = None
    contrast: Optional[float] = None
    background: Optional[float] = None
    power_law_B: Optional[float] = None
    power_law_P: Optional[float] = None
    rg_A: Optional[float] = None
    phi_actual: Optional[float] = None
    voxel_size: Optional[int] = None
    box_size_A: Optional[float] = None
    timestamp: Optional[str] = None
    params_std: dict = field(default_factory=dict)
    morphology_metrics: Optional[dict] = None
    data_q: Optional[list] = None
    data_I: Optional[list] = None
    data_I_corr: Optional[list] = None
    model_q: Optional[list] = None
    model_I: Optional[list] = None
    r_grid: Optional[list] = None
    gamma_r: Optional[list] = None
    spectral_k: Optional[list] = None
    spectral_F: Optional[list] = None


@dataclass
class WAXSPeak:
    index: int
    shape: str
    params: dict = field(default_factory=dict)        # name -> value
    params_std: dict = field(default_factory=dict)    # name -> stddev
    area: Optional[float] = None
    area_std: Optional[float] = None


@dataclass
class WAXSPeakFitResult(_BaseToolResult):
    tool: str = "waxs_peakfit"
    fit_quality: Optional[dict] = None  # robust fit-quality metrics (see core.fit_metrics)
    n_peaks: Optional[int] = None
    bg_shape: Optional[str] = None
    chi_squared: Optional[float] = None
    reduced_chi_squared: Optional[float] = None
    dof: Optional[int] = None
    q_min: Optional[float] = None
    q_max: Optional[float] = None
    timestamp: Optional[str] = None
    bg_params: dict = field(default_factory=dict)
    peaks: list = field(default_factory=list)         # list of WAXSPeak
    Q: Optional[list] = None
    I_fit: Optional[list] = None
    I_bg: Optional[list] = None
    residuals: Optional[list] = None
    intensity_data: Optional[list] = None
    intensity_error: Optional[list] = None


@dataclass
class FractalAggregateEntry:
    group_path: str
    name: str
    timestamp: Optional[str] = None
    label: Optional[str] = None
    Z: Optional[int] = None
    df: Optional[float] = None
    dmin: Optional[float] = None
    c: Optional[float] = None
    RgPrimary: Optional[float] = None
    RgAggregate: Optional[float] = None


@dataclass
class FractalsResult(_BaseToolResult):
    tool: str = "fractals"
    n_aggregates: int = 0
    aggregates: list = field(default_factory=list)  # list of FractalAggregateEntry


@dataclass
class MergeProvenance(_BaseToolResult):
    tool: str = "data_merge"
    scale: Optional[float] = None
    q_shift: Optional[float] = None
    background: Optional[float] = None
    chi_squared: Optional[float] = None
    q_overlap_min: Optional[float] = None
    q_overlap_max: Optional[float] = None
    ds1_file: Optional[str] = None
    ds2_file: Optional[str] = None
    timestamp: Optional[str] = None
    extra: dict = field(default_factory=dict)


@dataclass
class ManipulationProvenance(_BaseToolResult):
    tool: str = "data_manipulation"
    operation: Optional[str] = None
    source_file: Optional[str] = None
    timestamp: Optional[str] = None
    parameters: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Aggregate / cross-file schemas
# ---------------------------------------------------------------------------

@dataclass
class TabulationRow:
    path: str
    name: str
    sample: Optional[str] = None
    scan_number: Optional[int] = None
    mtime: Optional[str] = None
    value: Optional[float] = None
    stddev: Optional[float] = None


@dataclass
class Tabulation:
    folder: str
    tool: str
    parameter: str
    x_axis: str
    units: Optional[str] = None
    label: Optional[str] = None
    n_rows: int = 0
    rows: list = field(default_factory=list)   # list of TabulationRow

    def to_dict(self) -> dict:
        return asdict_clean(self)


@dataclass
class SampleSummary:
    folder: str
    sample: str
    n_files: int = 0
    files: list = field(default_factory=list)            # list of FileEntry
    analyses_count: dict[str, int] = field(default_factory=dict)
    parameter_ranges: dict = field(default_factory=dict)  # tool -> {param: {min, max, n}}

    def to_dict(self) -> dict:
        return asdict_clean(self)


# ---------------------------------------------------------------------------
# Plot schemas
# ---------------------------------------------------------------------------

@dataclass
class PlotResult:
    path: str                           # absolute saved path
    format: str = "png"
    width: int = 0
    height: int = 0
    n_files: int = 0
    base64_png: Optional[str] = None    # filled in only when requested

    def to_dict(self) -> dict:
        return asdict_clean(self)
