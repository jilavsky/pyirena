"""Cross-file aggregation: parameter trends, sample summaries.

Uses the existing :mod:`pyirena.io.schema` TOOL_REGISTRY so that supported
tool/parameter pairs are discoverable rather than hard-coded.
"""
from __future__ import annotations

from datetime import datetime, timezone
from typing import Optional

import h5py
import numpy as np

from pyirena.api._paths import resolve_safe_folder
from pyirena.api.discovery import _glob_many, _sample_name, _scan_number_from_name, _DEFAULT_GLOBS
from pyirena.api.schemas import FileEntry, SampleSummary, Tabulation, TabulationRow
from pyirena.io.schema import TOOL_REGISTRY


def _read_scalar_dataset(grp: h5py.Group, rel_path: str) -> Optional[float]:
    """Read a scalar at *rel_path* inside *grp*; return None if missing/bad."""
    if rel_path not in grp:
        return None
    try:
        ds = grp[rel_path]
        if not isinstance(ds, h5py.Dataset):
            return None
        v = ds[()]
        return float(v)
    except Exception:
        return None


def _read_scalar_attr(grp: h5py.Group, name: str) -> Optional[float]:
    if name in grp.attrs:
        try:
            return float(grp.attrs[name])
        except (TypeError, ValueError):
            return None
    return None


def _fallback_lookup(
    grp: h5py.Group, parameter: str,
) -> tuple[Optional[float], Optional[float], Optional[str]]:
    """Resolve *parameter* against a tool group when no schema spec matched.

    Tries, in order:
      1. ``params/<leaf>`` (with sibling ``params_std/<leaf>`` for stddev)
      2. ``derived/<leaf>`` (model-specific derived quantities, e.g.
         Invariant, VolumeFraction, QmaxUsed, Thickness, CorrLength)
      3. ``<leaf>`` as a scalar dataset on the group root
      4. ``<leaf>`` as a group attribute

    where *leaf* = parameter with optional ``param_`` prefix stripped, so both
    ``Kp`` and ``param_Kp`` resolve to ``params/Kp``. Returns (None, None,
    None) when nothing is found. Units are always None for fallback hits
    because they aren't enumerated in the schema.
    """
    leaf = parameter[6:] if parameter.startswith("param_") else parameter

    # 1) params/<leaf> + params_std/<leaf>
    val = _read_scalar_dataset(grp, f"params/{leaf}")
    if val is not None:
        std_val = _read_scalar_dataset(grp, f"params_std/{leaf}")
        return val, std_val, None

    # 1b) derived/<leaf> (Invariant, VolumeFraction, Thickness, …)
    val = _read_scalar_dataset(grp, f"derived/{leaf}")
    if val is not None:
        return val, None, None

    # 2) Top-level scalar dataset on the group
    val = _read_scalar_dataset(grp, leaf)
    if val is not None:
        std_val = _read_scalar_dataset(grp, f"{leaf}_std")
        return val, std_val, None

    # 3) Group attribute
    val = _read_scalar_attr(grp, leaf)
    if val is not None:
        return val, None, None

    return None, None, None


def _extract_scalar(f: h5py.File, tool: str, parameter: str,
                    subgroup_index: Optional[int] = None) -> tuple[Optional[float], Optional[float], Optional[str]]:
    """Extract (value, stddev, units) for *parameter* in *tool* from open file.

    Looks the parameter up in TOOL_REGISTRY[tool]['scalars']. For per-subgroup
    parameters, *subgroup_index* selects the 1-based occurrence (e.g. level 1,
    pop_01, peak_01). Stddev is searched at sibling path with '_std' suffix.
    Returns (None, None, None) if not found.

    Runtime fallback: tools like simple_fits store model-dependent parameter
    names under ``params/<name>``. The schema only enumerates the most common
    ones, so if the spec lookup fails we also probe ``params/<leaf>`` and
    ``<leaf>`` directly on the group, where *leaf* is the parameter name with
    an optional ``param_`` prefix stripped. This lets the AI surface read
    parameters specific to a model (Porod Kp, Debye-Bueche Lc, etc.) without
    requiring a schema entry for every possibility.
    """
    schema = TOOL_REGISTRY.get(tool)
    if schema is None:
        return None, None, None
    group_path = schema["group"]
    if group_path not in f:
        return None, None, None
    grp = f[group_path]

    # Find matching scalar spec (case-insensitive so AI can pass 'Rg' or 'rg')
    spec = None
    parameter_lc = parameter.lower()
    for sc in schema["scalars"]:
        if sc["key"].lower() == parameter_lc:
            spec = sc
            break
    if spec is None:
        return _fallback_lookup(grp, parameter)
    units = spec.get("units") or None

    rel_path = spec["path"]
    if spec.get("per_subgroup"):
        n = subgroup_index if subgroup_index is not None else 1
        sub = schema.get("sub_groups") or {}
        fmt = sub.get("subgroup_fmt", "d")
        try:
            idx_str = format(n, fmt)
        except (ValueError, TypeError):
            idx_str = str(n)
        rel_path = rel_path.replace("{n}", idx_str)

    # Try dataset first, then group attribute
    val = _read_scalar_dataset(grp, rel_path)
    if val is None:
        # Some scalars are attrs of the top group
        attr_name = rel_path.rsplit("/", 1)[-1]
        val = _read_scalar_attr(grp, attr_name)

    # Look for stddev at sibling _std path
    std_rel = rel_path + "_std"
    std_val = _read_scalar_dataset(grp, std_rel)
    if std_val is None:
        # alternative location: params_std/<name>
        if "/" in rel_path:
            parent, leaf = rel_path.rsplit("/", 1)
            alt = f"{parent}_std/{leaf}"
            std_val = _read_scalar_dataset(grp, alt)
            if std_val is None:
                alt2 = f"{parent}/params_std/{leaf}"
                std_val = _read_scalar_dataset(grp, alt2)

    return val, std_val, units


def tabulate_parameter(
    folder: str,
    tool: str,
    parameter: str,
    x_axis: str = "scan_number",
    subgroup_index: Optional[int] = None,
    pattern: str = ",".join(_DEFAULT_GLOBS),
    sample_filter: Optional[str] = None,
) -> dict:
    """Extract a single parameter from every file in *folder* that has *tool*.

    Returns a Tabulation with one TabulationRow per file. Files without the
    requested tool result group are skipped.

    Parameters
    ----------
    folder, pattern, sample_filter : see summarize_folder().
    tool : str
        Tool key from TOOL_REGISTRY (e.g. 'unified_fit', 'modeling').
    parameter : str
        Scalar key from TOOL_REGISTRY[tool]['scalars'] (e.g. 'Rg', 'background').
    x_axis : str
        Sort key for rows: 'scan_number', 'mtime', or 'name'.
    subgroup_index : int, optional
        1-based index into per-subgroup parameters (level_1, pop_01, peak_01).
        Required for per_subgroup scalars; defaults to 1 if omitted.
    """
    if x_axis not in ("scan_number", "mtime", "name"):
        raise ValueError("x_axis must be 'scan_number', 'mtime', or 'name'")
    schema = TOOL_REGISTRY.get(tool)
    if schema is None:
        raise ValueError(f"Unknown tool '{tool}'. "
                         f"Choose from: {sorted(TOOL_REGISTRY.keys())}")
    # Determine label/units up front
    label = None
    units = None
    for sc in schema["scalars"]:
        if sc["key"] == parameter:
            label = sc.get("label")
            units = sc.get("units") or None
            break

    folder_p = resolve_safe_folder(folder)
    files = _glob_many(folder_p, pattern)
    sample_filter_lc = sample_filter.lower() if sample_filter else None
    rows: list[TabulationRow] = []

    for path in files:
        try:
            with h5py.File(path, "r") as f:
                sample = _sample_name(f)
                if sample_filter_lc:
                    if not sample or sample_filter_lc not in sample.lower():
                        continue
                if schema["group"] not in f:
                    continue
                val, std_val, _ = _extract_scalar(f, tool, parameter, subgroup_index)
                stat = path.stat()
                rows.append(TabulationRow(
                    path=str(path),
                    name=path.name,
                    sample=sample,
                    scan_number=_scan_number_from_name(path.stem),
                    mtime=datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
                    value=val,
                    stddev=std_val,
                ))
        except (OSError, KeyError):
            continue

    # Sort rows
    if x_axis == "scan_number":
        rows.sort(key=lambda r: (r.scan_number is None, r.scan_number if r.scan_number is not None else 0))
    elif x_axis == "mtime":
        rows.sort(key=lambda r: r.mtime or "")
    else:
        rows.sort(key=lambda r: r.name.lower())

    return Tabulation(
        folder=str(folder_p),
        tool=tool,
        parameter=parameter,
        x_axis=x_axis,
        units=units,
        label=label,
        n_rows=len(rows),
        rows=rows,
    ).to_dict()


def summarize_sample(
    folder: str,
    sample: str,
    pattern: str = ",".join(_DEFAULT_GLOBS),
) -> dict:
    """One-sample condensation: file list + analyses count + parameter ranges.

    For each tool present in any file, computes min/max/n of every scalar
    parameter the tool's schema declares (top-level scalars only; per-
    subgroup parameters are skipped because the index is not implied).
    """
    from pyirena.api.discovery import _file_entry
    folder_p = resolve_safe_folder(folder)
    sample_lc = sample.lower()
    files = _glob_many(folder_p, pattern)

    matched: list[FileEntry] = []
    analyses_count: dict[str, int] = {}
    # tool -> parameter -> list of values
    values: dict[str, dict[str, list[float]]] = {}

    for path in files:
        try:
            with h5py.File(path, "r") as f:
                sname = _sample_name(f)
                if not sname or sample_lc not in sname.lower():
                    continue
                fe = _file_entry(path, deep=True)
                matched.append(fe)
                for key, schema in TOOL_REGISTRY.items():
                    if schema["group"] not in f:
                        continue
                    analyses_count[key] = analyses_count.get(key, 0) + 1
                    for sc in schema["scalars"]:
                        if sc.get("per_subgroup"):
                            continue
                        val, _, _ = _extract_scalar(f, key, sc["key"])
                        if val is None:
                            continue
                        values.setdefault(key, {}).setdefault(sc["key"], []).append(val)
        except (OSError, KeyError):
            continue

    parameter_ranges: dict = {}
    for tool, pmap in values.items():
        parameter_ranges[tool] = {}
        for pname, lst in pmap.items():
            if not lst:
                continue
            arr = np.array(lst, dtype=float)
            parameter_ranges[tool][pname] = {
                "min": float(np.nanmin(arr)),
                "max": float(np.nanmax(arr)),
                "n": int(np.sum(np.isfinite(arr))),
            }

    return SampleSummary(
        folder=str(folder_p),
        sample=sample,
        n_files=len(matched),
        files=matched,
        analyses_count=dict(sorted(analyses_count.items())),
        parameter_ranges=parameter_ranges,
    ).to_dict()
