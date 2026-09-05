"""Stable, read-only scattering-data discovery and loading API.

This module is intentionally independent of Qt. Applications such as
Bernardyn can discover every NXcanSAS entry or read a text curve without
depending on pyIrena GUI internals or creating converted sibling files.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import h5py
import numpy as np

from pyirena.io.hdf5 import (
    Q_UNIT_TO_ANGSTROM,
    list_nxcansas_datasets,
    readGenericNXcanSAS,
    readTextFile,
)


def _decode(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        return _decode(value.item())
    if isinstance(value, np.ndarray):
        return [_decode(item) for item in value.tolist()]
    if isinstance(value, Mapping):
        return {str(key): _decode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_decode(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _fingerprint(path: Path) -> str:
    stat = path.stat()
    value = f"{path.resolve()}:{stat.st_size}:{stat.st_mtime_ns}".encode()
    return hashlib.sha256(value).hexdigest()


@dataclass(frozen=True)
class ScatteringLocation:
    """One loadable scattering dataset inside a file."""

    path: Path
    internal_path: str | None = None
    display_name: str = ""
    variant: str = "default"
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        path = Path(self.path).expanduser().resolve()
        object.__setattr__(self, "path", path)
        object.__setattr__(self, "display_name", self.display_name or path.name)
        object.__setattr__(self, "metadata", dict(self.metadata))


@dataclass(frozen=True)
class ScatteringRecord:
    """A loaded 1-D curve in linear physical units."""

    q: np.ndarray
    intensity: np.ndarray
    uncertainty: np.ndarray | None = None
    dq: np.ndarray | None = None
    q_unit: str = "1/angstrom"
    intensity_unit: str = "1/cm"
    label: str = "Dataset"
    metadata: Mapping[str, Any] = field(default_factory=dict)
    provenance: Mapping[str, Any] = field(default_factory=dict)
    source_fingerprint: str | None = None


def discover_scattering(path: str | Path) -> list[ScatteringLocation]:
    """Discover loadable curves without reading their numeric arrays."""
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if source.suffix.lower() in (".txt", ".dat", ".csv"):
        return [ScatteringLocation(source, display_name=source.name, variant="text")]
    if source.suffix.lower() not in (".h5", ".hdf5", ".hdf", ".nxs"):
        raise ValueError(f"unsupported scattering-data suffix: {source.suffix}")
    try:
        datasets = list_nxcansas_datasets(
            str(source.parent), source.name, include_smr=True
        )
    except TypeError:
        datasets = list_nxcansas_datasets(str(source.parent), source.name)
    if not datasets:
        datasets = _discover_simple_hdf5(source)
    if not datasets:
        return [ScatteringLocation(source, display_name=source.name, variant="simple-hdf5")]
    return [
        ScatteringLocation(
            source,
            internal_path=str(item["path"]),
            display_name=f"{source.name}: {item.get('name', item['path'])}",
            variant="slit-smeared" if "smr" in str(item["path"]).lower() else "default",
            metadata={"entry": item.get("entry"), "name": item.get("name")},
        )
        for item in datasets
    ]


def _discover_simple_hdf5(source: Path) -> list[dict[str, str | None]]:
    """Find conventional Q/I pairs in otherwise unlabelled HDF5 groups."""
    discovered: list[dict[str, str | None]] = []

    def inspect_group(name: str, group: h5py.Group) -> None:
        names = {key.lower(): key for key in group.keys()}
        q_name = next((names.get(key) for key in ("q", "qvec", "q_vector") if key in names), None)
        i_name = next((names.get(key) for key in ("i", "intensity", "r", "data") if key in names), None)
        if q_name and i_name:
            path = f"/{name}" if name else "/"
            discovered.append({"path": path, "entry": None, "name": path})

    with h5py.File(source, "r") as handle:
        inspect_group("", handle)

        def visitor(name: str, obj: h5py.Group | h5py.Dataset) -> None:
            if isinstance(obj, h5py.Group):
                inspect_group(name, obj)

        handle.visititems(visitor)
    return discovered


def load_scattering(
    location: ScatteringLocation | str | Path,
    *,
    q_unit: str = "1/A",
    error_fraction: float = 0.05,
) -> ScatteringRecord:
    """Load one curve without modifying its source file."""
    if not isinstance(location, ScatteringLocation):
        discovered = discover_scattering(location)
        if len(discovered) != 1:
            raise ValueError(
                f"{location} contains {len(discovered)} datasets; pass a ScatteringLocation"
            )
        location = discovered[0]
    source = location.path
    if source.suffix.lower() in (".txt", ".dat", ".csv"):
        if q_unit not in Q_UNIT_TO_ANGSTROM:
            raise ValueError(f"unknown Q unit {q_unit!r}")
        result = readTextFile(
            str(source.parent), source.name, error_fraction=error_fraction, q_unit=q_unit
        )
        if result is None:
            raise ValueError(f"could not load scattering data from {source.name}")
        metadata = {
            "error_fraction": error_fraction,
            "q_unit_assumed": q_unit,
            "read_only_import": True,
        }
        adapter = "pyirena.io.hdf5.readTextFile"
        label = source.stem
    else:
        result = readGenericNXcanSAS(
            str(source.parent), source.name, data_path=location.internal_path
        )
        if result is None:
            raise ValueError(f"could not load scattering data from {location.display_name}")
        metadata = {
            key: _decode(value)
            for key, value in result.items()
            if key.lower().endswith(("attrs", "attributes"))
        }
        adapter = "pyirena.io.hdf5.readGenericNXcanSAS"
        label = location.display_name
    q_attrs = result.get(
        "QAttrs", result.get("Q_attrs", result.get("Q_attributes", {}))
    ) or {}
    i_attrs = result.get(
        "IntensityAttrs", result.get("I_attrs", result.get("Int_attributes", {}))
    ) or {}
    return ScatteringRecord(
        q=np.asarray(result["Q"], dtype=float),
        intensity=np.asarray(result["Intensity"], dtype=float),
        uncertainty=(
            None if result.get("Error") is None else np.asarray(result["Error"], dtype=float)
        ),
        dq=None if result.get("dQ") is None else np.asarray(result["dQ"], dtype=float),
        q_unit=str(_decode(q_attrs.get("units", "1/angstrom"))),
        intensity_unit=str(_decode(i_attrs.get("units", "1/cm"))),
        label=label,
        metadata=metadata,
        provenance={
            "source_name": source.name,
            "source_path": str(source),
            "internal_path": location.internal_path,
            "variant": location.variant,
            "adapter": adapter,
        },
        source_fingerprint=_fingerprint(source),
    )
