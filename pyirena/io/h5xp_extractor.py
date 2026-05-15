"""
pyirena/io/h5xp_extractor.py — Extract pyirena results from NXcanSAS HDF5
files and write them into an Igor Pro h5xp packed experiment.

This is the bridge between the pyirena HDF5 world and Igor Pro.  It reads
I(Q) data and every recognised result group present in a source file, then
calls :mod:`pyirena.io.h5xp_writer` to populate the h5xp.

Typical usage
-------------
Single file::

    from pyirena.io.h5xp_writer import create_h5xp, open_h5xp
    from pyirena.io.h5xp_extractor import extract_file_to_h5xp

    with create_h5xp("results.h5xp") as f:
        info = extract_file_to_h5xp("sample1.h5", f)
        print(info)

Batch (multiple files → one h5xp, also builds Results table)::

    from pyirena.io.h5xp_extractor import batch_extract_to_h5xp

    batch_extract_to_h5xp(["s1.h5", "s2.h5"], "batch.h5xp",
                           tools=["unified_fit", "size_distribution"])

Folder name in Igor
-------------------
Each source file becomes one subfolder inside ``Packed Data/SAXS/`` (or
``WAXS/``).  The folder name defaults to the file stem (no extension), but
can be overridden.  Igor folder names must not contain colons; any ``:`` in
the stem is replaced with ``_``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np

from pyirena.io.h5xp_writer import (
    make_wave_note,
    open_h5xp,
    create_h5xp,
    write_iq_data,
    write_result_wave,
    write_results_table,
)
from pyirena.io.igor_names import SIMPLE_FIT_MODEL_WAVE, TOOL_CROSS_REF
from pyirena.io.schema import TOOL_REGISTRY


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Suffix convention for Results-table waves that come from numbered sub-groups.
# These must produce unique, Igor-legal wave names.
_SUBGROUP_SUFFIX: dict[str, str] = {
    "unified_fit":  "L",   # Rg_L1, G_L1, B_L1, P_L1, …
    "waxs_peakfit": "P",   # peak_Q0_P1, peak_A_P1, …
    "modeling":     "",    # pop_Rg_1, pop_G_1, …
    "fractals":     "",    # agg_Z_1, agg_df_1, …
}


def _subgroup_result_key(tool_key: str, sc_key: str, n: int) -> str:
    """Build the Results-table wave name for a per-subgroup scalar."""
    suffix = _SUBGROUP_SUFFIX.get(tool_key, "")
    if suffix:
        return f"{sc_key}_{suffix}{n}"   # e.g. Rg_L1, peak_Q0_P1
    return f"{sc_key}_{n}"              # e.g. pop_Rg_1, agg_Z_1


def _safe_str(val: Any) -> str:
    """Decode bytes or numpy string scalar to plain str."""
    if isinstance(val, (bytes, np.bytes_)):
        return val.decode("utf-8", errors="replace")
    if isinstance(val, np.ndarray):
        val = val.flat[0]
        if isinstance(val, (bytes, np.bytes_)):
            return val.decode("utf-8", errors="replace")
    return str(val)


def _scalar(grp: h5py.Group, path: str, default: float = float("nan")) -> float:
    """Read a scalar from *grp*: tries dataset first, then group attribute.

    *path* is treated as a simple name (no slash) for the attribute lookup;
    nested paths (e.g. ``"params/Rg"``) only match datasets.
    """
    # Try as a dataset (new format)
    try:
        return float(grp[path][()])
    except (KeyError, TypeError, ValueError):
        pass
    # Fall back to a group attribute (old format, e.g. background/chi_squared on group)
    if "/" not in path:
        try:
            return float(grp.attrs[path])
        except (KeyError, TypeError, ValueError):
            pass
    return default


def _array(grp: h5py.Group, path: str) -> np.ndarray | None:
    """Read an array dataset; return None if absent."""
    if path in grp:
        return np.asarray(grp[path][:], dtype=np.float64)
    return None


def _igor_folder_name(stem: str) -> str:
    """Make a stem safe for use as an Igor data-folder name."""
    return stem.replace(":", "_").replace("/", "_").replace("\\", "_")


# ---------------------------------------------------------------------------
# NXcanSAS I(Q) reader
# ---------------------------------------------------------------------------

def _read_iq(source: h5py.File) -> dict[str, np.ndarray | None]:
    """Find the SASdata group and return Q, I, error, dQ arrays.

    Returns a dict with keys: Q, I, error (may be None), dq (may be None).
    Raises ValueError if no SASdata group is found.
    """
    from pyirena.io.hdf5 import find_matching_groups

    paths = find_matching_groups(source, {"canSAS_class": "SASdata"}, {})
    if not paths:
        # Fall back: look for any group named 'sasdata'
        paths = []
        def _find_sasdata(name, obj):
            if isinstance(obj, h5py.Group) and name.endswith("sasdata"):
                paths.append(name)
        source.visititems(_find_sasdata)

    if not paths:
        raise ValueError("No SASdata group found in file")

    sas = source[paths[0]]
    q     = np.asarray(sas["Q"][:],    dtype=np.float64)
    I     = np.asarray(sas["I"][:],    dtype=np.float64)
    error = np.asarray(sas["Idev"][:], dtype=np.float64) if "Idev" in sas else None
    dq    = np.asarray(sas["Qdev"][:], dtype=np.float64) if "Qdev" in sas else None
    return {"Q": q, "I": I, "error": error, "dq": dq}


# ---------------------------------------------------------------------------
# Metadata reader
# ---------------------------------------------------------------------------

_HPLANCK_EV_ANG = 12398.4   # h·c in eV·Å  (so E[keV] = 12.3984 / λ[Å])


def _read_group_scalars(grp: h5py.Group, prefix: str, out: dict[str, Any],
                        max_depth: int = 2) -> None:
    """Recursively collect all scalar datasets from *grp* into *out*.

    Keys are built as ``prefix + dataset_name`` (groups deepen the prefix with
    ``prefix + group_name + "_"``).  Only scalars (shape == () or (1,)) and
    short strings are collected; longer arrays are skipped.
    """
    for name, item in grp.items():
        key = prefix + name
        if isinstance(item, h5py.Dataset):
            if item.shape == () or item.shape == (1,):
                v = item[()]
                if isinstance(v, np.ndarray):
                    v = v.flat[0]
                if isinstance(v, (bytes, np.bytes_)):
                    s = v.decode("utf-8", errors="replace").strip()
                    if s:
                        out[key] = s
                elif isinstance(v, (int, float, np.integer, np.floating)):
                    out[key] = float(v)
                elif isinstance(v, str) and v.strip():
                    out[key] = v.strip()
        elif isinstance(item, h5py.Group) and max_depth > 0:
            _read_group_scalars(item, key + "_", out, max_depth - 1)


def _read_metadata(source: h5py.File, folder_name: str) -> dict[str, Any]:
    """Collect comprehensive wave-note metadata from a NXcanSAS entry group.

    Reads all scalar datasets from entry-level, sample/, instrument/ (depth 2)
    and computes X-ray energy from wavelength when available.
    """
    meta: dict[str, Any] = {"IgorFolder": folder_name}
    entry = source.get("entry")
    if entry is None:
        return meta

    # ── Top-level entry scalars ──────────────────────────────────────────────
    for key, path in [
        ("Title",       "title"),
        ("StartTime",   "start_time"),
        ("EndTime",     "end_time"),
        ("Definition",  "definition"),
        ("RunCycle",    "run_cycle"),
    ]:
        ds = entry.get(path)
        if ds is not None and ds.shape in ((), (1,)):
            v = ds[()]
            s = _safe_str(v).strip()
            if s:
                meta[key] = s

    # ── Sample group ─────────────────────────────────────────────────────────
    sample = entry.get("sample")
    if sample is not None:
        _read_group_scalars(sample, "sample_", meta, max_depth=1)
        # Friendly aliases so the most-used fields have clean names too
        for alias, raw_key in [
            ("SampleName",        "sample_name"),
            ("SampleDescription", "sample_description"),
        ]:
            if raw_key in meta and alias not in meta:
                meta[alias] = meta[raw_key]

    # ── Instrument group (depth 2: sub-groups are source/beam/detector/…) ───
    instrument = entry.get("instrument")
    if instrument is not None:
        _read_group_scalars(instrument, "instrument_", meta, max_depth=2)

    # ── Wavelength / energy ──────────────────────────────────────────────────
    # Check common NXcanSAS wavelength paths; compute energy in keV.
    wl: float | None = None
    for wl_path in (
        "instrument/beam/incident_wavelength",
        "instrument/monochromator/wavelength",
        "instrument/source/wavelength",
    ):
        ds = entry.get(wl_path)
        if ds is not None and ds.shape in ((), (1,)):
            v = ds[()]
            if isinstance(v, np.ndarray):
                v = float(v.flat[0])
            if isinstance(v, (int, float, np.integer, np.floating)) and float(v) > 0:
                wl = float(v)
                break
    if wl is not None:
        meta["Wavelength_A"]  = round(wl, 6)
        meta["Energy_keV"]    = round(_HPLANCK_EV_ANG / wl / 1000.0, 4)

    return meta


# ---------------------------------------------------------------------------
# Per-tool extraction functions
# ---------------------------------------------------------------------------

def _extract_unified_fit(grp: h5py.Group, h5xp: h5py.File,
                          folder: str, category: str) -> bool:
    """Write UnifiedFitIntensity wave from unified_fit_results group."""
    q_model = _array(grp, "Q")
    I_model = _array(grp, "intensity_model")
    if q_model is None or I_model is None:
        return False

    # Collect scalars for wave note
    params: dict[str, Any] = {
        "chi_squared": _scalar(grp, "chi_squared"),
        "background":  _scalar(grp, "background"),
    }
    # Per-level parameters (new format: scalar datasets; old format: group attrs)
    n = 1
    while f"level_{n}" in grp:
        lg = grp[f"level_{n}"]
        for pname in ("Rg", "G", "B", "P", "ETA", "PACK", "RgCutoff"):
            v = _scalar(lg, pname)   # _scalar already tries attrs as fallback
            if not np.isnan(v):
                params[f"{pname}_L{n}"] = v
        n += 1

    write_result_wave(h5xp, folder, "UnifiedFitIntensity",
                      q_model, I_model, params, category)
    return True


def _extract_size_distribution(grp: h5py.Group, h5xp: h5py.File,
                                 folder: str, category: str) -> bool:
    """Write SizesFitIntensity, SizesVolumeDistribution, SizesNumberDistribution."""
    wrote_any = False

    params: dict[str, Any] = {
        "chi_squared":     _scalar(grp, "chi_squared"),
        "volume_fraction": _scalar(grp, "volume_fraction"),
        "Rg":              _scalar(grp, "rg"),
        "q_power":         _scalar(grp, "q_power"),
    }

    q_model = _array(grp, "Q")
    I_model = _array(grp, "intensity_model")
    if q_model is not None and I_model is not None:
        write_result_wave(h5xp, folder, "SizesFitIntensity",
                          q_model, I_model, params, category)
        wrote_any = True

    r_grid = _array(grp, "r_grid")
    if r_grid is not None:
        for ds_name, igor_name in [
            ("distribution",   "SizesVolumeDistribution"),
            ("number_dist",    "SizesNumberDistribution"),
            ("cumul_vol_dist", "CumulativeSizeDist"),
        ]:
            y = _array(grp, ds_name)
            if y is not None:
                write_result_wave(h5xp, folder, igor_name,
                                  r_grid, y, params, category)
                wrote_any = True

    return wrote_any


def _extract_waxs_peakfit(grp: h5py.Group, h5xp: h5py.File,
                           folder: str, category: str) -> bool:
    """Write SADModelIntensity and per-peak SADModelIntPeak{n} waves."""
    wrote_any = False

    params_global: dict[str, Any] = {
        "chi_squared":         _scalar(grp, "chi_squared"),
        "reduced_chi_squared": _scalar(grp, "reduced_chi_squared"),
    }

    q_fit = _array(grp, "Q")
    I_fit = _array(grp, "I_fit")
    if q_fit is not None and I_fit is not None:
        write_result_wave(h5xp, folder, "SADModelIntensity",
                          q_fit, I_fit, params_global, category)
        wrote_any = True

    # Per-peak waves
    n = 1
    while f"peak_{n:02d}" in grp:
        pk = grp[f"peak_{n:02d}"]
        q_pk = _array(pk, "Q_peak")
        I_pk = _array(pk, "I_peak")
        if q_pk is not None and I_pk is not None:
            pk_params: dict[str, Any] = dict(params_global)
            pk_params["peak_n"] = n
            # Peak fit parameters
            if "params" in pk:
                for pname in ("Q0", "A", "FWHM", "eta"):
                    v = _scalar(pk["params"], pname)
                    if not np.isnan(v):
                        pk_params[pname] = v
            write_result_wave(h5xp, folder, f"SADModelIntPeak{n}",
                              q_pk, I_pk, pk_params, category)
            wrote_any = True
        n += 1

    return wrote_any


def _extract_simple_fits(grp: h5py.Group, h5xp: h5py.File,
                          folder: str, category: str) -> bool:
    """Write model wave with model-dependent Igor name."""
    q_model = _array(grp, "Q")
    I_model = _array(grp, "I_model")
    if q_model is None or I_model is None:
        return False

    model_name = _safe_str(grp.attrs.get("model", b""))
    igor_name  = SIMPLE_FIT_MODEL_WAVE.get(model_name, "SimFitUnknownI")

    params: dict[str, Any] = {
        "model":               model_name,
        "chi_squared":         _scalar(grp, "chi_squared"),
        "reduced_chi_squared": _scalar(grp, "reduced_chi_squared"),
    }
    if "params" in grp:
        for name, ds in grp["params"].items():
            if isinstance(ds, h5py.Dataset) and ds.shape == ():
                params[name] = float(ds[()])

    write_result_wave(h5xp, folder, igor_name, q_model, I_model, params, category)
    return True


def _extract_modeling(grp: h5py.Group, h5xp: h5py.File,
                       folder: str, category: str) -> bool:
    """Write ModelingIntensity and per-population distribution waves."""
    wrote_any = False

    params_global: dict[str, Any] = {
        "chi_squared":         _scalar(grp, "chi_squared"),
        "reduced_chi_squared": _scalar(grp, "reduced_chi_squared"),
        "background":          _scalar(grp, "background"),
    }

    q_model = _array(grp, "model_q")
    I_model = _array(grp, "model_I")
    if q_model is not None and I_model is not None:
        write_result_wave(h5xp, folder, "ModelingIntensity",
                          q_model, I_model, params_global, category)
        wrote_any = True

    # Per-population distributions (size_dist populations only)
    n = 1
    while f"pop_{n:02d}" in grp:
        pg = grp[f"pop_{n:02d}"]
        pop_type = _safe_str(pg.attrs.get("pop_type", b""))
        if pop_type == "size_dist":
            r_grid = _array(pg, "radius_grid")
            if r_grid is not None:
                for ds_name, igor_name in [
                    ("volume_dist", f"ModelingVolDist_Pop{n}"),
                    ("number_dist", f"ModelingNumDist_Pop{n}"),
                ]:
                    y = _array(pg, ds_name)
                    if y is not None:
                        pp: dict[str, Any] = dict(params_global)
                        pp["pop_n"] = n
                        pp["pop_type"] = pop_type
                        write_result_wave(h5xp, folder, igor_name,
                                          r_grid, y, pp, category)
                        wrote_any = True
        n += 1

    return wrote_any


def _extract_saxs_morph(grp: h5py.Group, h5xp: h5py.File,
                         folder: str, category: str) -> bool:
    """Write SAXS Morph model and correlation waves (pyirena-only names)."""
    wrote_any = False

    params: dict[str, Any] = {
        "chi_squared":     _scalar(grp, "chi_squared"),
        "volume_fraction": _scalar(grp, "volume_fraction"),
        "contrast":        _scalar(grp, "contrast"),
        "rg_A":            _scalar(grp, "rg_A"),
        "background":      _scalar(grp, "background"),
        "power_law_B":     _scalar(grp, "power_law_B"),
        "power_law_P":     _scalar(grp, "power_law_P"),
    }

    for x_ds, y_ds, igor_name in [
        ("model_q",    "model_I",    "SAXSMorphModelI"),
        ("data_q",     "data_I",     "SAXSMorphDataI"),
        ("r_grid",     "gamma_r",    "SAXSMorphGammaR"),
        ("spectral_k", "spectral_F", "SAXSMorphSpectralF"),
    ]:
        x = _array(grp, x_ds)
        y = _array(grp, y_ds)
        if x is not None and y is not None:
            write_result_wave(h5xp, folder, igor_name, x, y, params, category)
            wrote_any = True

    return wrote_any


def _extract_fractals(grp: h5py.Group, h5xp: h5py.File,
                       folder: str, category: str) -> bool:
    """Write FractFitIntensity and per-aggregate waves."""
    wrote_any = False

    # Total intensity
    q_int = _array(grp, "intensity/Q")
    I_uni = _array(grp, "intensity/I_unified")
    if q_int is not None and I_uni is not None:
        write_result_wave(h5xp, folder, "FractFitIntensity",
                          q_int, I_uni, {}, category)
        wrote_any = True

    # Per-aggregate components
    n = 1
    while f"aggregate_{n}" in grp:
        ag = grp[f"aggregate_{n}"]
        q_ag = _array(ag, "intensity/Q")
        params: dict[str, Any] = {"aggregate_n": n}
        if "parameters" in ag:
            for pname in ("Z", "df", "dmin", "c", "RgAggregate", "RgPrimary"):
                v = _scalar(ag["parameters"], pname)
                if not np.isnan(v):
                    params[pname] = v
        for ds_name, igor_name in [
            ("intensity/mass_fractal",    f"Mass{n}FractFitInt"),
            ("intensity/surface_fractal", f"Surf{n}FractFitInt"),
        ]:
            y = _array(ag, ds_name)
            if q_ag is not None and y is not None:
                write_result_wave(h5xp, folder, igor_name,
                                  q_ag, y, params, category)
                wrote_any = True
        n += 1

    return wrote_any


# ---------------------------------------------------------------------------
# Tool dispatch table
# ---------------------------------------------------------------------------

_TOOL_EXTRACTORS = {
    "unified_fit":      _extract_unified_fit,
    "size_distribution": _extract_size_distribution,
    "waxs_peakfit":     _extract_waxs_peakfit,
    "simple_fits":      _extract_simple_fits,
    "modeling":         _extract_modeling,
    "saxs_morph":       _extract_saxs_morph,
    "fractals":         _extract_fractals,
    # data_merge / data_manipulation have no plottable arrays → skip
}


# ---------------------------------------------------------------------------
# Public API — single file
# ---------------------------------------------------------------------------

def extract_file_to_h5xp(
    source_path: str | Path,
    h5xp_file: h5py.File,
    folder_name: str | None = None,
    tools: list[str] | None = None,
    category: str | None = None,
) -> dict[str, Any]:
    """Extract I(Q) data and results from one pyirena HDF5 file into *h5xp_file*.

    Parameters
    ----------
    source_path:
        Path to the NXcanSAS / pyirena HDF5 file to read.
    h5xp_file:
        Open h5py.File for the destination h5xp (from
        :func:`~pyirena.io.h5xp_writer.create_h5xp` or
        :func:`~pyirena.io.h5xp_writer.open_h5xp`).
    folder_name:
        Igor data-folder name for this file.  Defaults to the file stem
        (no extension).  Colons are replaced with underscores.
    tools:
        List of tool keys to export (e.g. ``["unified_fit"]``).  ``None``
        exports every result group that is present in the file.  Use
        ``pyirena.io.schema.TOOL_REGISTRY`` for the full list of keys.
    category:
        ``"SAXS"`` or ``"WAXS"``.  If ``None`` (default), auto-detected:
        ``"WAXS"`` when a ``waxs_peakfit_results`` group is present,
        otherwise ``"SAXS"``.

    Returns
    -------
    dict
        Summary with keys:

        * ``folder`` — Igor folder path written (e.g. ``"Packed Data/SAXS/sample1"``)
        * ``iq_written`` — True if I(Q) waves were written
        * ``tools_written`` — list of tool keys that were successfully exported
        * ``tools_skipped`` — list of tool keys that were present but had no data
        * ``tools_absent``  — list of tool keys requested but not in the file
    """
    src = Path(source_path)
    if folder_name is None:
        folder_name = _igor_folder_name(src.stem)

    result: dict[str, Any] = {
        "folder":        "",
        "iq_written":    False,
        "tools_written": [],
        "tools_skipped": [],
        "tools_absent":  [],
    }

    with h5py.File(src, "r") as f:
        # Auto-detect SAXS vs WAXS
        if category is None:
            category = "WAXS" if "entry/waxs_peakfit_results" in f else "SAXS"

        # Read metadata for wave note
        meta = _read_metadata(f, folder_name)

        # Read and write I(Q) data
        try:
            iq = _read_iq(f)
            write_iq_data(
                h5xp_file, folder_name,
                iq["Q"], iq["I"],
                error=iq["error"], dq=iq["dq"],
                wave_note=meta, category=category,
            )
            result["iq_written"] = True
        except (ValueError, KeyError) as exc:
            result["iq_error"] = str(exc)

        result["folder"] = f"Packed Data/{category}/{folder_name}"

        # Determine which tools to try
        tool_keys = tools if tools is not None else list(_TOOL_EXTRACTORS.keys())

        for key in tool_keys:
            if key not in _TOOL_EXTRACTORS:
                continue
            schema = TOOL_REGISTRY.get(key)
            if schema is None:
                continue
            group_path = schema["group"]
            if group_path not in f:
                result["tools_absent"].append(key)
                continue
            grp = f[group_path]
            ok = _TOOL_EXTRACTORS[key](grp, h5xp_file, folder_name, category)
            if ok:
                result["tools_written"].append(key)
            else:
                result["tools_skipped"].append(key)

    return result


# ---------------------------------------------------------------------------
# Public API — batch
# ---------------------------------------------------------------------------

def batch_extract_to_h5xp(
    source_paths: list[str | Path],
    h5xp_path: str | Path,
    tools: list[str] | None = None,
    category: str | None = None,
    overwrite: bool = False,
    build_results_table: bool = True,
) -> list[dict[str, Any]]:
    """Extract multiple pyirena files into one h5xp, optionally building a Results table.

    Parameters
    ----------
    source_paths:
        Ordered list of pyirena HDF5 file paths.
    h5xp_path:
        Destination h5xp path.  Created if it does not exist; appended to if
        it does (unless *overwrite* is True).
    tools:
        Tool keys to export from each file.  None = all present tools.
    category:
        Override SAXS/WAXS for all files.  None = auto-detect per file.
    overwrite:
        If True, overwrite an existing h5xp.  Default False (append).
    build_results_table:
        If True (default), after all files are processed, write a Results
        table for each scalar that was successfully collected.

    Returns
    -------
    list[dict]
        One summary dict per source file (from :func:`extract_file_to_h5xp`).
    """
    h5xp_p = Path(h5xp_path)
    open_fn = (lambda p: create_h5xp(p, overwrite=True)) if (overwrite or not h5xp_p.exists()) \
              else open_h5xp

    summaries: list[dict[str, Any]] = []
    # Accumulate scalars for the Results table: {tool_key: {param_key: [values...]}}
    results_accum: dict[str, dict[str, list[float]]] = {}
    sample_names: list[str] = []

    ctx = open_fn(h5xp_p) if not overwrite and h5xp_p.exists() else create_h5xp(h5xp_p, overwrite=overwrite)
    with ctx as h5xp:
        for src in source_paths:
            src = Path(src)
            info = extract_file_to_h5xp(src, h5xp, tools=tools, category=category)
            summaries.append(info)

            folder_name = _igor_folder_name(src.stem)
            sample_names.append(folder_name)

            if not build_results_table:
                continue

            # Collect scalars from each tool group for the Results table
            with h5py.File(src, "r") as f:
                for key in (tools or list(TOOL_REGISTRY.keys())):
                    schema = TOOL_REGISTRY.get(key)
                    if schema is None:
                        continue
                    grp_path = schema["group"]
                    if grp_path not in f:
                        continue
                    grp = f[grp_path]
                    tool_scalars = results_accum.setdefault(key, {})

                    for sc in schema["scalars"]:
                        if not sc["per_subgroup"]:
                            # Flat scalar: one value per file
                            val = _scalar(grp, sc["path"])
                            tool_scalars.setdefault(sc["key"], []).append(val)
                        else:
                            # Per-subgroup: enumerate actual sub-groups
                            sub_info = schema.get("sub_groups")
                            if sub_info is None:
                                continue
                            prefix = sub_info["prefix"]
                            fmt = sub_info.get("subgroup_fmt", "d")
                            n = 1
                            while True:
                                sg_name = prefix + format(n, fmt)
                                if sg_name not in grp:
                                    break
                                sg = grp[sg_name]
                                # Build the scalar path for this sub-group
                                # sc["path"] template uses "{n}" (bare int)
                                sc_path = sc["path"].replace("{n}", str(n))
                                # Strip the "prefix_{n}/" prefix that's already
                                # inside the subgroup object
                                sc_path_rel = sc_path.split("/", 1)[1] if "/" in sc_path else sc_path
                                val = _scalar(sg, sc_path_rel)
                                # Key name: append suffix for clarity
                                result_key = _subgroup_result_key(key, sc["key"], n)
                                tool_scalars.setdefault(result_key, []).append(val)
                                n += 1

        if build_results_table:
            tool_label_map = {k: v["label"] for k, v in TOOL_REGISTRY.items()}
            for tool_key, param_dict in results_accum.items():
                technique = tool_label_map.get(tool_key, tool_key)
                for param_key, values in param_dict.items():
                    # Pad missing entries with NaN so all waves have the same length
                    while len(values) < len(sample_names):
                        values.append(float("nan"))
                    # Find units from schema (flat scalars first, then per-subgroup)
                    schema_scalars = TOOL_REGISTRY[tool_key]["scalars"]
                    sc_spec = next(
                        (s for s in schema_scalars if s["key"] == param_key),
                        None,
                    )
                    units = sc_spec["units"] if sc_spec else ""
                    write_results_table(h5xp, technique, param_key,
                                        values, sample_names, units)

    return summaries
