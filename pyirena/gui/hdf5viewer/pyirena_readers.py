"""
PyirenaReaders — knows where pyirena tools store their results in HDF5.

All functions accept a file path and return a normalised dict suitable for
GraphWindow.add_curve() calls.  Every function is safe to call on a file
that doesn't have the expected data — it returns None in that case.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import h5py
import numpy as np


# ── Detection ──────────────────────────────────────────────────────────────

def detect_available_data(filepath: str | Path) -> list[str]:
    """
    Return a list of data-type keys available in *filepath*.

    Keys ∈ {"nxcansas", "unified_fit", "sizes", "waxs", "simple_fit"}.
    """
    available = []
    try:
        with h5py.File(str(filepath), "r") as f:
            if _has_nxcansas(f):
                available.append("nxcansas")
            if "entry/unified_fit_results" in f:
                available.append("unified_fit")
            if "entry/sizes_results" in f:
                available.append("sizes")
            if "entry/waxs_peakfit_results" in f:
                available.append("waxs")
            if "entry/simple_fit_results" in f:
                available.append("simple_fit")
    except Exception:
        pass
    return available


def _has_nxcansas(f: h5py.File) -> bool:
    """Return True if the file has at least one NXcanSAS entry."""
    found = []
    def _visitor(name, obj):
        if isinstance(obj, h5py.Group):
            css = obj.attrs.get("canSAS_class", b"")
            if isinstance(css, bytes):
                css = css.decode("utf-8", errors="ignore")
            if css == "SASentry":
                found.append(name)
                return True   # stop visiting
    try:
        f.visititems(_visitor)
    except Exception:
        pass
    return bool(found)


# ── NXcanSAS ───────────────────────────────────────────────────────────────

def read_nxcansas(filepath: str | Path) -> dict | None:
    """
    Read the first NXcanSAS entry.

    Returns {"Q", "I", "dI" (or None), "label"} or None on failure.
    """
    try:
        from pyirena.io.hdf5 import readGenericNXcanSAS
        result = readGenericNXcanSAS(
            str(Path(filepath).parent),
            Path(filepath).name,
        )
        if result is None:
            return None
        return {
            "Q":     np.asarray(result["Q"],         float),
            "I":     np.asarray(result["Intensity"],  float),
            "dI":    np.asarray(result["Error"],      float) if result.get("Error") is not None else None,
            "label": result.get("label", Path(filepath).stem),
        }
    except Exception:
        return None


# ── Unified Fit ────────────────────────────────────────────────────────────

def read_unified_fit(filepath: str | Path) -> dict | None:
    """
    Read Unified Fit results.

    Returns {
        "Q", "I_model", "I_data" (may be None),
        "levels": [{"level", "Rg", "G", "B", "P", "ETA", "PACK", …}],
        "chi2", "background", "label"
    } or None.
    """
    try:
        from pyirena.io.nxcansas_unified import load_unified_fit_results
        res = load_unified_fit_results(Path(filepath))
        if not res:
            return None
        _id = res.get("intensity_data")
        return {
            "Q":       np.asarray(res["Q"],              float),
            "I_model": np.asarray(res["intensity_model"], float),
            "I_data":  np.asarray(_id, float) if _id is not None and len(_id) > 0 else None,
            "levels":  res.get("levels", []),
            "chi2":    res.get("chi_squared"),
            "background": res.get("background"),
            "label":   Path(filepath).stem,
        }
    except Exception:
        return None


# ── Size Distribution ──────────────────────────────────────────────────────

def read_sizes(filepath: str | Path) -> dict | None:
    """
    Read Size Distribution results.

    Returns {
        "Q", "I_model", "r", "distribution", "distribution_std" (may be None),
        "chi2", "label"
    } or None.
    """
    try:
        from pyirena.io.nxcansas_sizes import load_sizes_results
        res = load_sizes_results(Path(filepath))
        if not res:
            return None
        dist_std = res.get("distribution_std")
        return {
            "Q":              np.asarray(res["Q"],              float),
            "I_model":        np.asarray(res["intensity_model"], float),
            "r":              np.asarray(res["r_grid"],          float),
            "distribution":   np.asarray(res["distribution"],   float),
            "distribution_std": np.asarray(dist_std, float) if dist_std is not None else None,
            "chi2":           res.get("chi_squared"),
            "label":          Path(filepath).stem,
        }
    except Exception:
        return None


# ── WAXS Peak Fit ──────────────────────────────────────────────────────────

def read_waxs(filepath: str | Path) -> dict | None:
    """
    Read WAXS Peak Fit results.

    Returns {
        "Q", "I_fit", "I_bg",
        "peaks": [{"shape", "Q0", "A", "FWHM", "Q_peak", "I_peak"}],
        "chi2", "label"
    } or None.
    """
    try:
        from pyirena.io.nxcansas_waxs_peakfit import load_waxs_peakfit_results
        res = load_waxs_peakfit_results(Path(filepath))
        if not res:
            return None
        return {
            "Q":     np.asarray(res.get("Q",      []), float),
            "I_fit": np.asarray(res.get("I_fit",  []), float),
            "I_bg":  np.asarray(res.get("I_bg",   []), float),
            "peaks": res.get("peaks", []),
            "chi2":  res.get("chi_squared"),
            "label": Path(filepath).stem,
        }
    except Exception:
        return None


# ── Simple Fits ────────────────────────────────────────────────────────────

def read_simple_fit(filepath: str | Path) -> dict | None:
    """
    Read Simple Fits results.

    Returns {
        "Q", "I_model", "params", "chi2", "model_name", "label"
    } or None.
    """
    try:
        from pyirena.io.nxcansas_simple_fits import load_simple_fit_results
        res = load_simple_fit_results(Path(filepath))
        if not res:
            return None
        return {
            "Q":         np.asarray(res.get("Q",       []), float),
            "I_model":   np.asarray(res.get("I_model", []), float),
            "params":    res.get("params",     {}),
            "chi2":      res.get("chi_squared"),
            "model_name": res.get("model", ""),
            "label":     Path(filepath).stem,
        }
    except Exception:
        return None


# ── Arbitrary dataset ──────────────────────────────────────────────────────

def read_dataset(filepath: str | Path, hdf5_path: str) -> np.ndarray | None:
    """
    Read any 1D or scalar dataset from an HDF5 file.

    Returns a numpy array, or None on failure.
    """
    try:
        with h5py.File(str(filepath), "r") as f:
            ds = f[hdf5_path]
            data = ds[()]
            if isinstance(data, (bytes, np.bytes_)):
                return None   # string dataset
            arr = np.asarray(data, float)
            return arr.flatten() if arr.ndim > 1 else arr
    except Exception:
        return None


# ── Value collection ───────────────────────────────────────────────────────

def collect_value(filepath: str | Path, spec: dict) -> float | None:
    """
    Extract a single scalar value from *filepath* according to *spec*.

    spec examples::

        {"type": "unified_fit", "item": "Rg", "level": 1}
        {"type": "unified_fit", "item": "G",  "level": 2}
        {"type": "unified_fit", "item": "chi2"}
        {"type": "sizes",      "item": "chi2"}
        {"type": "sizes",      "item": "volume_fraction"}
        {"type": "waxs",       "item": "Q0",  "peak": 0}
        {"type": "waxs",       "item": "A",   "peak": 0}
        {"type": "waxs",       "item": "chi2"}
        {"type": "simple_fit", "item": "param", "param_name": "Rg"}
        {"type": "simple_fit", "item": "chi2"}
        {"type": "custom",     "path": "/entry/unified_fit_results/level_1"}
                                                  # reads @Rg attr or dataset
    """
    data_type = spec.get("type", "")
    item      = spec.get("item", "")

    try:
        if data_type == "unified_fit":
            return _collect_unified(filepath, item, spec.get("level", 1))

        if data_type == "sizes":
            return _collect_sizes(filepath, item)

        if data_type == "waxs":
            return _collect_waxs(filepath, item, spec.get("peak", 0))

        if data_type == "simple_fit":
            return _collect_simple_fit(filepath, item, spec.get("param_name", ""))

        if data_type == "custom":
            return read_metadata_value(filepath, spec.get("path", ""))

    except Exception:
        pass
    return None


def _read_scalar_value(node, key: str) -> float | None:
    """Read a scalar from an HDF5 node: dataset first (new format), attr fallback (old)."""
    try:
        if key in node and isinstance(node[key], h5py.Dataset):
            return float(node[key][()])
        val = node.attrs.get(key)
        return float(val) if val is not None else None
    except Exception:
        return None


def _collect_unified(filepath, item: str, level: int) -> float | None:
    try:
        with h5py.File(str(filepath), "r") as f:
            grp = f["entry/unified_fit_results"]
            if item == "chi2":
                return _read_scalar_value(grp, "chi_squared")
            if item == "background":
                return _read_scalar_value(grp, "background")
            level_grp = grp[f"level_{level}"]
            return _read_scalar_value(level_grp, item)
    except Exception:
        return None


def _collect_sizes(filepath, item: str) -> float | None:
    try:
        with h5py.File(str(filepath), "r") as f:
            grp = f["entry/sizes_results"]
            return _read_scalar_value(grp, item)
    except Exception:
        return None


def _collect_waxs(filepath, item: str, peak: int) -> float | None:
    try:
        with h5py.File(str(filepath), "r") as f:
            grp = f["entry/waxs_peakfit_results"]
            if item == "chi2":
                return _read_scalar_value(grp, "chi_squared")
            peak_name = f"peak_{peak+1:02d}"
            params_grp = grp[peak_name]["params"]
            return float(params_grp[item][()])
    except Exception:
        return None


def _collect_simple_fit(filepath, item: str, param_name: str) -> float | None:
    try:
        with h5py.File(str(filepath), "r") as f:
            grp = f["entry/simple_fit_results"]
            if item == "chi2":
                return _read_scalar_value(grp, "chi_squared")
            if item == "param" and param_name:
                return float(grp["params"][param_name][()])
    except Exception:
        return None


def read_metadata_value(filepath: str | Path, hdf5_path: str) -> float | None:
    """
    Read a scalar (attribute or dataset) from *hdf5_path* in *filepath*.

    Supports both dataset paths (``/entry/foo/bar``) and attribute paths
    (``/entry/foo/bar@attr_name``).
    """
    if not hdf5_path:
        return None
    try:
        with h5py.File(str(filepath), "r") as f:
            if "@" in hdf5_path:
                path, attr = hdf5_path.rsplit("@", 1)
                path = path.rstrip("/") or "/"
                val = f[path].attrs[attr]
            else:
                val = f[hdf5_path][()]
            # Convert bytes → string → float where possible
            if isinstance(val, (bytes, np.bytes_)):
                val = val.decode("utf-8", errors="replace")
            return float(val)
    except Exception:
        return None


# ── X-axis helpers for Collect Values ─────────────────────────────────────

def extract_sort_key_value(filepath: str | Path, sort_index: int) -> float | None:
    """
    Extract a numeric X value from the filename using the sort key at *sort_index*.

    Uses the same _SORT_KEY functions as FileTreeWidget.
    """
    from .file_tree import _SORT_KEYS
    name = Path(filepath).name
    try:
        val = _SORT_KEYS[sort_index](name)
        return float(val) if val != float("inf") else None
    except Exception:
        return None
