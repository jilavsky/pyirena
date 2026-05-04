"""
NXcanSAS / HDF5 I/O for Fractal Aggregate results.

A grown aggregate is saved into an `entry/fractals_results/aggregate_{N}`
NXprocess group.  Multiple aggregates can coexist in the same file with
auto-incremented suffixes — the loader and lister enumerate all of them.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from pyirena.core.fractals import (
    FractalAggregate, FractalParams, GrowthConfig,
)


_ROOT_GROUP = "entry/fractals_results"


def _next_aggregate_index(parent: h5py.Group) -> int:
    """Return the next available aggregate_{N} index in `parent`."""
    used = []
    for name in parent.keys():
        if name.startswith("aggregate_"):
            try:
                used.append(int(name.split("_", 1)[1]))
            except ValueError:
                continue
    return max(used) + 1 if used else 1


def save_fractal_aggregate(
    filepath: Path,
    agg: FractalAggregate,
    group_index: Optional[int] = None,
) -> str:
    """Append `agg` to `entry/fractals_results/aggregate_{N}` in *filepath*.

    Parameters
    ----------
    filepath : Path
        Target HDF5 file.  Created if it does not exist.
    agg : FractalAggregate
    group_index : int or None
        Specific suffix to use; auto-increment when None.

    Returns
    -------
    str — the full HDF5 path of the written group, e.g. "entry/fractals_results/aggregate_3".
    """
    filepath = Path(filepath)
    timestamp = datetime.now().isoformat()

    with h5py.File(filepath, "a") as f:
        # Ensure /entry exists with NXcanSAS-friendly attributes
        if "entry" not in f:
            f.attrs["default"] = "entry"
            f.attrs["file_name"] = filepath.name
            f.attrs["file_time"] = timestamp
            f.attrs["creator"] = "pyirena"
            f.attrs["NeXus_version"] = "4.3.0"
            f.attrs["HDF5_version"] = h5py.version.hdf5_version
            f.attrs["h5py_version"] = h5py.version.version
            entry = f.create_group("entry")
            entry.attrs["NX_class"] = "NXentry"
            entry.attrs["canSAS_class"] = "SASentry"
            entry.create_dataset("definition", data="NXsas")

        # Ensure /entry/fractals_results parent group exists
        if _ROOT_GROUP not in f:
            parent = f.create_group(_ROOT_GROUP)
            parent.attrs["NX_class"] = "NXcollection"
        else:
            parent = f[_ROOT_GROUP]

        idx = group_index if group_index is not None else _next_aggregate_index(parent)
        gpath = f"{_ROOT_GROUP}/aggregate_{idx}"
        if gpath in f:
            del f[gpath]

        g = f.create_group(gpath)
        g.attrs["NX_class"] = "NXprocess"
        g.attrs["analysis_type"] = "Mass Fractal Aggregate"
        g.attrs["program"] = "pyirena"
        g.attrs["timestamp"] = timestamp
        g.attrs["uuid"] = agg.uuid
        g.attrs["label"] = agg.label or ""

        # Aggregate geometry
        g.create_dataset("positions", data=agg.positions.astype(np.int32),
                         compression="gzip", compression_opts=4)
        g.create_dataset("neighbor_list", data=agg.neighbor_list.astype(np.int32),
                         compression="gzip", compression_opts=4)
        g.create_dataset("neighbor_count", data=agg.neighbor_count.astype(np.uint8),
                         compression="gzip", compression_opts=4)
        g.create_dataset("attempt_value", data=int(agg.attempt_value))

        # Computed parameters (scalar datasets, browsable in HDF5 viewers)
        params_g = g.create_group("parameters")
        for name, value in (
            ("Z", agg.params.z),
            ("dmin", agg.params.dmin),
            ("c", agg.params.c),
            ("df", agg.params.df),
            ("R_dimensionless", agg.params.R_dimensionless),
            ("p", agg.params.p),
            ("s", agg.params.s),
            ("RgPrimary", agg.params.rg_primary),
            ("RgAggregate", agg.params.rg_aggregate),
            ("PrimaryDiameter", agg.params.primary_diameter),
            ("TrueStickingProbability", agg.params.true_sticking_prob),
            ("NumEndpoints", agg.params.num_endpoints),
            ("NumPathsUsed", agg.params.num_paths_used),
        ):
            params_g.create_dataset(name, data=float(value) if not isinstance(value, int) else int(value))

        # Input parameters
        ip_g = g.create_group("input_params")
        ip_g.create_dataset("StickingProbability", data=float(agg.config.sticking_prob))
        ip_g.create_dataset("NumberOfTestPaths", data=int(agg.config.num_test_paths))
        ip_g.create_dataset("AllowedNearDistance", data=int(agg.config.allowed_near_dist))
        ip_g.attrs["MultiParticleAttraction"] = str(agg.config.attraction)
        ip_g.create_dataset("Seed", data=int(agg.config.seed))

        # Optional intensity arrays
        if agg.q is not None and (agg.i_unified is not None or agg.i_montecarlo is not None):
            i_g = g.create_group("intensity")
            i_g.create_dataset("Q", data=np.asarray(agg.q, dtype=np.float64))
            if agg.i_unified is not None:
                i_g.create_dataset("I_unified",
                                    data=np.asarray(agg.i_unified, dtype=np.float64))
            if agg.i_montecarlo is not None:
                i_g.create_dataset("I_montecarlo",
                                    data=np.asarray(agg.i_montecarlo, dtype=np.float64))

    return gpath


def list_fractal_aggregates(filepath: Path) -> list[dict]:
    """Return one summary dict per aggregate_{N} group in *filepath*."""
    filepath = Path(filepath)
    if not filepath.exists():
        return []
    out = []
    with h5py.File(filepath, "r") as f:
        if _ROOT_GROUP not in f:
            return []
        parent = f[_ROOT_GROUP]
        for name in sorted(parent.keys()):
            if not name.startswith("aggregate_"):
                continue
            g = parent[name]
            params = g.get("parameters")
            entry = {
                "group_path": f"{_ROOT_GROUP}/{name}",
                "name": name,
                "timestamp": g.attrs.get("timestamp", ""),
                "label": g.attrs.get("label", ""),
                "Z": int(params["Z"][()]) if (params and "Z" in params) else 0,
                "dmin": float(params["dmin"][()]) if (params and "dmin" in params) else float("nan"),
                "c": float(params["c"][()]) if (params and "c" in params) else float("nan"),
                "df": float(params["df"][()]) if (params and "df" in params) else float("nan"),
            }
            out.append(entry)
    return out


def load_fractal_aggregate(filepath: Path, group_path: str) -> FractalAggregate:
    """Load a single aggregate from *filepath* at the given HDF5 path."""
    filepath = Path(filepath)
    with h5py.File(filepath, "r") as f:
        if group_path not in f:
            raise KeyError(f"{group_path} not found in {filepath}")
        g = f[group_path]

        positions = np.asarray(g["positions"][...], dtype=np.int32)
        neighbor_list = np.asarray(g["neighbor_list"][...], dtype=np.int32)
        neighbor_count = np.asarray(g["neighbor_count"][...], dtype=np.uint8)
        attempt_value = int(g["attempt_value"][()]) if "attempt_value" in g else 0

        params_g = g["parameters"]
        params = FractalParams(
            z=int(params_g["Z"][()]),
            dmin=float(params_g["dmin"][()]),
            c=float(params_g["c"][()]),
            df=float(params_g["df"][()]),
            R_dimensionless=float(params_g["R_dimensionless"][()]),
            p=float(params_g["p"][()]),
            s=float(params_g["s"][()]),
            rg_primary=float(params_g["RgPrimary"][()]),
            rg_aggregate=float(params_g["RgAggregate"][()]),
            primary_diameter=float(params_g["PrimaryDiameter"][()]),
            true_sticking_prob=float(params_g["TrueStickingProbability"][()]),
            num_endpoints=int(params_g["NumEndpoints"][()]),
            num_paths_used=int(params_g["NumPathsUsed"][()]),
        )
        ip_g = g["input_params"]
        config = GrowthConfig(
            z=int(params.z),
            sticking_prob=float(ip_g["StickingProbability"][()]),
            num_test_paths=int(ip_g["NumberOfTestPaths"][()]),
            rg_primary=float(params.rg_primary),
            allowed_near_dist=int(ip_g["AllowedNearDistance"][()]),
            attraction=str(ip_g.attrs.get("MultiParticleAttraction", "Neutral")),
            seed=int(ip_g["Seed"][()]),
        )

        q = i_u = i_mc = None
        if "intensity" in g:
            i_g = g["intensity"]
            if "Q" in i_g:
                q = np.asarray(i_g["Q"][...], dtype=np.float64)
            if "I_unified" in i_g:
                i_u = np.asarray(i_g["I_unified"][...], dtype=np.float64)
            if "I_montecarlo" in i_g:
                i_mc = np.asarray(i_g["I_montecarlo"][...], dtype=np.float64)

        agg = FractalAggregate(
            positions=positions,
            neighbor_list=neighbor_list,
            neighbor_count=neighbor_count,
            params=params,
            config=config,
            attempt_value=attempt_value,
            q=q, i_unified=i_u, i_montecarlo=i_mc,
            uuid=str(g.attrs.get("uuid", "")) or None or "",
            created_at=str(g.attrs.get("timestamp", datetime.now().isoformat())),
            label=str(g.attrs.get("label", "")),
        )
        return agg
