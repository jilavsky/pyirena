"""Per-tool result readers — thin facades over pyirena.io.nxcansas_* loaders.

Each `read_<tool>()` returns a JSON-serializable dict (decimated arrays,
NaN -> None). Missing data returns a `found=False` dict; the caller does
not need a try/except.

All functions accept an `include_arrays: bool = False` kwarg controlling
whether the heavy I/Q/residual arrays are returned at full fidelity. By
default arrays are decimated to PYIRENA_MAX_ARRAY_POINTS (500 by default).
"""
from __future__ import annotations

import contextlib
import io
from typing import Optional

import h5py
import numpy as np

from pyirena.api._paths import resolve_safe_file
from pyirena.api.schemas import (
    FractalAggregateEntry,
    FractalsResult,
    ManipulationProvenance,
    MergeProvenance,
    ModelingPopulation,
    ModelingResult,
    SAXSMorphResult,
    SimpleFitResult,
    SizeDistResult,
    UnifiedFitLevel,
    UnifiedFitResult,
    WAXSPeak,
    WAXSPeakFitResult,
    array_to_list,
)


def _silent_call(fn, *args, **kwargs):
    """Call *fn* with stdout/stderr suppressed (some loaders print)."""
    with contextlib.redirect_stdout(io.StringIO()), \
         contextlib.redirect_stderr(io.StringIO()):
        return fn(*args, **kwargs)


def _arr(raw: dict, key: str, max_points: Optional[int], include: bool):
    """Return array_to_list(raw[key]) when *include* is True, else None.

    *max_points* tunes the decimation when arrays are included. Pass
    include=True to opt in (matches the include_arrays=True flag on each
    public read_<tool> function).
    """
    if not include:
        return None
    return array_to_list(raw.get(key), max_points=max_points, include_full=False)


# ---------------------------------------------------------------------------
# Simple Fits
# ---------------------------------------------------------------------------

def read_simple_fit(path: str, include_arrays: bool = False,
                    max_points: Optional[int] = None) -> dict:
    """Read Simple Fits results. Returns ``found=False`` if absent."""
    file_p = resolve_safe_file(path)
    result = SimpleFitResult(path=str(file_p))
    from pyirena.io.nxcansas_simple_fits import load_simple_fit_results
    try:
        raw = _silent_call(load_simple_fit_results, file_p)
    except (KeyError, OSError):
        return result.to_dict()
    if raw is None:
        return result.to_dict()

    result.found = True
    result.model = raw.get("model")
    result.success = raw.get("success")
    result.chi_squared = raw.get("chi_squared")
    result.reduced_chi_squared = raw.get("reduced_chi_squared")
    result.dof = raw.get("dof")
    result.q_min = raw.get("q_min")
    result.q_max = raw.get("q_max")
    result.n_mc_runs = raw.get("n_mc_runs")
    result.timestamp = raw.get("timestamp")
    result.params = {k: float(v) for k, v in (raw.get("params") or {}).items()}
    result.params_std = {k: float(v) for k, v in (raw.get("params_std") or {}).items()}
    result.derived = {k: float(v) for k, v in (raw.get("derived") or {}).items()}
    result.Q = _arr(raw, "Q", max_points, include_arrays)
    result.I_model = _arr(raw, "I_model", max_points, include_arrays)
    result.residuals = _arr(raw, "residuals", max_points, include_arrays)
    result.intensity_data = _arr(raw, "intensity_data", max_points, include_arrays)
    result.intensity_error = _arr(raw, "intensity_error", max_points, include_arrays)
    result.fit_quality = raw.get("fit_quality")
    return result.to_dict()


# ---------------------------------------------------------------------------
# Unified Fit
# ---------------------------------------------------------------------------

def read_unified_fit(path: str, include_arrays: bool = False,
                     max_points: Optional[int] = None) -> dict:
    """Read Unified Fit (Beaucage) results."""
    file_p = resolve_safe_file(path)
    result = UnifiedFitResult(path=str(file_p))
    from pyirena.io.nxcansas_unified import load_unified_fit_results
    try:
        raw = _silent_call(load_unified_fit_results, file_p)
    except (ValueError, KeyError, OSError):
        return result.to_dict()
    if raw is None:
        return result.to_dict()

    result.found = True
    result.num_levels = raw.get("num_levels")
    result.background = raw.get("background")
    result.background_err = raw.get("background_err")
    result.chi_squared = raw.get("chi_squared")
    result.reduced_chi_squared = raw.get("reduced_chi_squared")
    result.timestamp = raw.get("timestamp")

    levels: list = []
    for i, lv in enumerate(raw.get("levels", []) or [], 1):
        if not isinstance(lv, dict):
            continue
        level = UnifiedFitLevel(level_number=int(lv.get("level_number", i)))
        for key in ("G", "Rg", "B", "P", "RgCutoff", "ETA", "PACK"):
            v = lv.get(key)
            if v is not None:
                try:
                    setattr(level, key, float(v))
                except (TypeError, ValueError):
                    pass
        corr = lv.get("correlations")
        if corr is None:
            corr = lv.get("correlated")
        if corr is not None:
            try:
                level.correlations = bool(corr)
            except (TypeError, ValueError):
                pass
        for key in ("G_err", "Rg_err", "B_err", "P_err", "ETA_err", "PACK_err"):
            v = lv.get(key)
            if v is not None:
                try:
                    setattr(level, key, float(v))
                except (TypeError, ValueError):
                    pass
        levels.append(level)
    result.levels = levels

    result.Q = _arr(raw, "Q", max_points, include_arrays)
    result.intensity_data = _arr(raw, "intensity_data", max_points, include_arrays)
    result.intensity_model = _arr(raw, "intensity_model", max_points, include_arrays)
    result.intensity_model_ideal = _arr(raw, "intensity_model_ideal", max_points, include_arrays)
    result.intensity_error = _arr(raw, "intensity_error", max_points, include_arrays)
    result.residuals = _arr(raw, "residuals", max_points, include_arrays)
    result.fit_quality = raw.get("fit_quality")
    # Slit-smearing provenance (absent => legacy pinhole file).
    if raw.get("slit_length") is not None:
        result.slit_length = float(raw["slit_length"])
    if raw.get("data_is_slit_smeared") is not None:
        result.data_is_slit_smeared = bool(raw["data_is_slit_smeared"])
    return result.to_dict()


# ---------------------------------------------------------------------------
# Size Distribution
# ---------------------------------------------------------------------------

def read_size_distribution(path: str, include_arrays: bool = False,
                           max_points: Optional[int] = None) -> dict:
    """Read Size Distribution results."""
    file_p = resolve_safe_file(path)
    result = SizeDistResult(path=str(file_p))
    from pyirena.io.nxcansas_sizes import load_sizes_results
    try:
        raw = _silent_call(load_sizes_results, file_p)
    except (KeyError, OSError):
        return result.to_dict()
    if raw is None:
        return result.to_dict()

    result.found = True
    for key in ("method", "shape", "chi_squared", "reduced_chi_squared",
                "volume_fraction", "rg", "n_iterations", "q_power", "contrast",
                "aspect_ratio", "r_min", "r_max", "n_bins", "background",
                "power_law_B", "power_law_P", "timestamp"):
        if key in raw and raw[key] is not None:
            v = raw[key]
            if isinstance(v, bytes):
                v = v.decode("utf-8", errors="replace")
            try:
                if key in ("n_iterations", "n_bins"):
                    setattr(result, key, int(v))
                elif key in ("shape", "method", "timestamp"):
                    setattr(result, key, str(v))
                else:
                    setattr(result, key, float(v))
            except (TypeError, ValueError):
                pass

    for k in ("Q", "intensity_data", "intensity_model", "intensity_error",
              "residuals", "r_grid", "distribution", "distribution_std",
              "number_dist", "cumul_vol_dist"):
        setattr(result, k, _arr(raw, k, max_points, include_arrays))
    result.fit_quality = raw.get("fit_quality")
    return result.to_dict()


# ---------------------------------------------------------------------------
# Modeling
# ---------------------------------------------------------------------------

def read_modeling(path: str, include_arrays: bool = False,
                  max_points: Optional[int] = None) -> dict:
    """Read parametric Modeling results (size_dist, unified_level, etc.)."""
    file_p = resolve_safe_file(path)
    result = ModelingResult(path=str(file_p))
    from pyirena.io.nxcansas_modeling import load_modeling_results
    try:
        raw = _silent_call(load_modeling_results, file_p)
    except (OSError,):
        return result.to_dict()
    if raw is None:
        return result.to_dict()

    result.found = True
    for key in ("chi_squared", "reduced_chi_squared", "background",
                "q_min", "q_max", "timestamp"):
        result.__dict__[key] = raw.get(key)
    if raw.get("dof") is not None:
        try:
            result.dof = int(raw["dof"])
        except (TypeError, ValueError):
            pass
    result.model_q = _arr(raw, "model_q", max_points, include_arrays)
    result.model_I = _arr(raw, "model_I", max_points, include_arrays)

    pops_out: list = []
    for pd in raw.get("populations", []) or []:
        if not isinstance(pd, dict):
            continue
        pop = ModelingPopulation(
            population_index=int(pd.get("population_index", 0)),
            pop_type=str(pd.get("pop_type", "size_dist")),
            enabled=bool(pd.get("enabled", True)),
            label=str(pd.get("label", "")),
        )
        # Per-type parameters: collect scalar fields into `parameters`
        type_specific_keys = {
            "unified_level":   ["G", "Rg", "B", "P", "RgCO", "ETA", "PACK", "correlations"],
            "diffraction_peak": ["peak_type", "position", "amplitude", "width", "eta_voigt"],
            "guinier_porod":   ["G", "Rg1", "s1", "P", "Rg2", "s2", "RgCO", "ETA", "PACK", "correlations"],
            "mass_fractal":    ["Phi", "Radius", "Beta", "Dv", "Ksi", "Eta", "Contrast"],
            "surface_fractal": ["Surface", "Ds", "Ksi", "Contrast", "Qc", "QcWidth", "use_porod_transition"],
        }
        for k in type_specific_keys.get(pop.pop_type, []):
            if k in pd and pd[k] is not None:
                v = pd[k]
                if isinstance(v, bool):
                    pop.parameters[k] = bool(v)
                elif isinstance(v, str):
                    pop.parameters[k] = v
                else:
                    try:
                        pop.parameters[k] = float(v)
                    except (TypeError, ValueError):
                        pop.parameters[k] = v

        # Size-distribution specific
        if pop.pop_type == "size_dist":
            pop.dist_type = pd.get("dist_type")
            pop.form_factor = pd.get("form_factor")
            pop.structure_factor = pd.get("structure_factor")
            pop.parameters["dist_params"] = dict(pd.get("dist_params") or {})
            pop.parameters["ff_params"] = dict(pd.get("ff_params") or {})
            pop.parameters["sf_params"] = dict(pd.get("sf_params") or {})
            for k in ("scale", "contrast"):
                if k in pd and pd[k] is not None:
                    try:
                        pop.parameters[k] = float(pd[k])
                    except (TypeError, ValueError):
                        pass
            if include_arrays:
                pop.radius_grid = array_to_list(pd.get("radius_grid"),
                                                max_points=max_points)
                pop.volume_dist = array_to_list(pd.get("volume_dist"),
                                                max_points=max_points)
                pop.number_dist = array_to_list(pd.get("number_dist"),
                                                max_points=max_points)

        if include_arrays:
            pop.model_I = array_to_list(pd.get("model_I"),
                                        max_points=max_points)
        derived = pd.get("derived") or {}
        pop.derived = {k: float(v) for k, v in derived.items()
                       if v is not None}
        pops_out.append(pop)
    result.populations = pops_out
    result.fit_quality = raw.get("fit_quality")
    return result.to_dict()


# ---------------------------------------------------------------------------
# SAXS Morph
# ---------------------------------------------------------------------------

def read_saxs_morph(path: str, include_arrays: bool = False,
                    max_points: Optional[int] = None) -> dict:
    """Read SAXS Morph results.

    Note: the voxelgram (3-D uint8 array) is intentionally not returned —
    it is too large for an LLM context. Use h5py directly to access it.
    """
    file_p = resolve_safe_file(path)
    result = SAXSMorphResult(path=str(file_p))
    from pyirena.io.nxcansas_saxs_morph import load_saxs_morph_results
    try:
        raw = _silent_call(load_saxs_morph_results, file_p)
    except OSError:
        return result.to_dict()
    if raw is None:
        return result.to_dict()

    result.found = True
    for key in ("chi_squared", "reduced_chi_squared", "volume_fraction",
                "contrast", "background", "power_law_B", "power_law_P",
                "rg_A", "phi_actual", "box_size_A", "timestamp"):
        v = raw.get(key)
        if v is not None:
            try:
                if key == "timestamp":
                    result.__dict__[key] = str(v)
                else:
                    result.__dict__[key] = float(v)
            except (TypeError, ValueError):
                pass
    if raw.get("dof") is not None:
        try:
            result.dof = int(raw["dof"])
        except (TypeError, ValueError):
            pass
    if raw.get("voxel_size") is not None:
        try:
            result.voxel_size = int(raw["voxel_size"])
        except (TypeError, ValueError):
            pass
    result.params_std = {k: float(v) for k, v in (raw.get("params_std") or {}).items()
                         if v is not None}

    mm = raw.get("morphology_metrics")
    if mm is not None:
        if hasattr(mm, "to_dict"):
            result.morphology_metrics = mm.to_dict()
        elif isinstance(mm, dict):
            result.morphology_metrics = mm

    for k in ("data_q", "data_I", "data_I_corr", "model_q", "model_I",
              "r_grid", "gamma_r", "spectral_k", "spectral_F"):
        setattr(result, k, _arr(raw, k, max_points, include_arrays))
    return result.to_dict()


# ---------------------------------------------------------------------------
# WAXS Peak Fit
# ---------------------------------------------------------------------------

def read_waxs_peakfit(path: str, include_arrays: bool = False,
                      max_points: Optional[int] = None) -> dict:
    """Read WAXS Peak Fit results."""
    file_p = resolve_safe_file(path)
    result = WAXSPeakFitResult(path=str(file_p))
    from pyirena.io.nxcansas_waxs_peakfit import load_waxs_peakfit_results
    try:
        raw = _silent_call(load_waxs_peakfit_results, file_p)
    except (KeyError, OSError):
        return result.to_dict()
    if raw is None:
        return result.to_dict()

    result.found = True
    for key in ("n_peaks", "bg_shape", "chi_squared", "reduced_chi_squared",
                "dof", "q_min", "q_max", "timestamp"):
        v = raw.get(key)
        if v is not None:
            try:
                if key in ("n_peaks", "dof"):
                    setattr(result, key, int(v))
                elif key in ("bg_shape", "timestamp"):
                    setattr(result, key, str(v))
                else:
                    setattr(result, key, float(v))
            except (TypeError, ValueError):
                pass

    # Background params: collapse {name: {value, fit, lo, hi}} to {name: value}
    bg_in = raw.get("bg_params") or {}
    bg_std_in = raw.get("bg_params_std") or {}
    result.bg_params = {}
    for name, spec in bg_in.items():
        if isinstance(spec, dict) and "value" in spec:
            entry = {"value": float(spec["value"])}
            if name in bg_std_in:
                entry["std"] = float(bg_std_in[name])
            result.bg_params[name] = entry

    # Peaks list
    peaks_out: list = []
    raw_peaks = raw.get("peaks") or []
    raw_peaks_std = raw.get("peaks_std") or []
    for i, pk in enumerate(raw_peaks, 1):
        if not isinstance(pk, dict):
            continue
        shape = pk.get("shape", "Unknown")
        peak = WAXSPeak(index=i, shape=str(shape))
        for k, v in pk.items():
            if k == "shape":
                continue
            if isinstance(v, dict) and "value" in v:
                try:
                    peak.params[k] = float(v["value"])
                except (TypeError, ValueError):
                    pass
            elif k in ("area", "area_std") and v is not None:
                try:
                    setattr(peak, k, float(v))
                except (TypeError, ValueError):
                    pass
        if i - 1 < len(raw_peaks_std) and isinstance(raw_peaks_std[i - 1], dict):
            for k, v in raw_peaks_std[i - 1].items():
                try:
                    peak.params_std[k] = float(v)
                except (TypeError, ValueError):
                    pass
        peaks_out.append(peak)
    result.peaks = peaks_out

    for k in ("Q", "I_fit", "I_bg", "residuals", "intensity_data", "intensity_error"):
        setattr(result, k, _arr(raw, k, max_points, include_arrays))
    result.fit_quality = raw.get("fit_quality")
    return result.to_dict()


# ---------------------------------------------------------------------------
# Fractals (list of aggregates)
# ---------------------------------------------------------------------------

def read_fractals(path: str) -> dict:
    """Read summary of grown fractal aggregates stored in *path*.

    Returns metadata only (one row per aggregate). Use h5py directly to load
    full per-aggregate intensity arrays.
    """
    file_p = resolve_safe_file(path)
    result = FractalsResult(path=str(file_p))
    from pyirena.io.nxcansas_fractals import list_fractal_aggregates
    try:
        rows = _silent_call(list_fractal_aggregates, file_p)
    except OSError:
        return result.to_dict()
    if not rows:
        return result.to_dict()

    result.found = True
    result.n_aggregates = len(rows)
    result.aggregates = [
        FractalAggregateEntry(
            group_path=r.get("group_path", ""),
            name=r.get("name", ""),
            timestamp=r.get("timestamp") or None,
            label=r.get("label") or None,
            Z=int(r["Z"]) if r.get("Z") is not None else None,
            df=float(r["df"]) if r.get("df") is not None else None,
            dmin=float(r["dmin"]) if r.get("dmin") is not None else None,
            c=float(r["c"]) if r.get("c") is not None else None,
        )
        for r in rows
    ]
    return result.to_dict()


# ---------------------------------------------------------------------------
# Data Merge provenance
# ---------------------------------------------------------------------------

def read_merge_provenance(path: str) -> dict:
    """Read Data Merge provenance (scale, q_shift, background, source files)."""
    file_p = resolve_safe_file(path)
    result = MergeProvenance(path=str(file_p))
    try:
        with h5py.File(file_p, "r") as f:
            if "entry/data_merge_results" not in f:
                return result.to_dict()
            grp = f["entry/data_merge_results"]
            result.found = True
            result.timestamp = grp.attrs.get("timestamp")
            if isinstance(result.timestamp, bytes):
                result.timestamp = result.timestamp.decode("utf-8", errors="replace")
            for ds_name in ("scale", "q_shift", "background", "chi_squared",
                            "q_overlap_min", "q_overlap_max"):
                if ds_name in grp:
                    try:
                        v = float(grp[ds_name][()])
                        setattr(result, ds_name, v)
                    except (TypeError, ValueError):
                        pass
            for ds_name, target in (("ds1_file", "ds1_file"), ("ds2_file", "ds2_file")):
                if ds_name in grp:
                    v = grp[ds_name][()]
                    if isinstance(v, bytes):
                        v = v.decode("utf-8", errors="replace")
                    setattr(result, target, str(v))
                elif f"{ds_name.replace('_file', '')}_path" in grp.attrs:
                    v = grp.attrs[f"{ds_name.replace('_file', '')}_path"]
                    if isinstance(v, bytes):
                        v = v.decode("utf-8", errors="replace")
                    setattr(result, target, str(v))
            # Any other scalar datasets land in extra
            for key in grp:
                if key in ("scale", "q_shift", "background", "chi_squared",
                           "q_overlap_min", "q_overlap_max", "ds1_file", "ds2_file"):
                    continue
                try:
                    ds = grp[key]
                    if isinstance(ds, h5py.Dataset) and ds.shape == ():
                        val = ds[()]
                        if isinstance(val, bytes):
                            val = val.decode("utf-8", errors="replace")
                        result.extra[key] = val
                except Exception:
                    continue
    except OSError:
        pass
    return result.to_dict()


# ---------------------------------------------------------------------------
# Data Manipulation provenance
# ---------------------------------------------------------------------------

def read_manipulation_provenance(path: str) -> dict:
    """Read Data Manipulation provenance (operation, parameters, source)."""
    file_p = resolve_safe_file(path)
    result = ManipulationProvenance(path=str(file_p))
    try:
        with h5py.File(file_p, "r") as f:
            if "entry/data_manipulation_results" not in f:
                return result.to_dict()
            grp = f["entry/data_manipulation_results"]
            result.found = True
            result.timestamp = grp.attrs.get("timestamp")
            if isinstance(result.timestamp, bytes):
                result.timestamp = result.timestamp.decode("utf-8", errors="replace")
            if "operation" in grp:
                v = grp["operation"][()]
                if isinstance(v, bytes):
                    v = v.decode("utf-8", errors="replace")
                result.operation = str(v)
            if "source_file" in grp:
                v = grp["source_file"][()]
                if isinstance(v, bytes):
                    v = v.decode("utf-8", errors="replace")
                result.source_file = str(v)
            for key in grp:
                if key in ("operation", "source_file"):
                    continue
                try:
                    ds = grp[key]
                    if isinstance(ds, h5py.Dataset) and ds.shape == ():
                        val = ds[()]
                        if isinstance(val, bytes):
                            val = val.decode("utf-8", errors="replace")
                        elif isinstance(val, (np.floating, np.integer)):
                            val = float(val)
                        result.parameters[key] = val
                except Exception:
                    continue
    except OSError:
        pass
    return result.to_dict()
