"""File-level discovery: list files, summarize a folder, inspect one file.

These are the first calls an AI agent makes to orient itself.
"""
from __future__ import annotations

import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import h5py

from pyirena.api._paths import resolve_safe_file, resolve_safe_folder
from pyirena.api.schemas import FileEntry, FileInspection, FolderSummary
from pyirena.io.schema import TOOL_REGISTRY, available_tools


_DEFAULT_GLOBS = ("*.h5", "*.hdf5", "*.hdf", "*.nx", "*.nxs")
_VALID_SORT = frozenset({"name_asc", "name_desc", "mtime_asc", "mtime_desc",
                         "size_asc", "size_desc"})


def _glob_many(folder: Path, pattern: str) -> list[Path]:
    """Glob *pattern* in *folder*. Supports comma-separated patterns."""
    out: list[Path] = []
    for pat in (p.strip() for p in pattern.split(",")):
        if pat:
            out.extend(folder.glob(pat))
    # Dedup, preserving order
    seen: set[Path] = set()
    uniq: list[Path] = []
    for p in out:
        if p not in seen and p.is_file():
            seen.add(p)
            uniq.append(p)
    return uniq


_SCAN_RE = re.compile(r"(\d+)(?!.*\d)")  # last integer in a string


def _scan_number_from_name(stem: str) -> Optional[int]:
    m = _SCAN_RE.search(stem)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            return None
    return None


def _sample_name(f: h5py.File) -> Optional[str]:
    """Pull sample name from common NXcanSAS locations."""
    # Try standard paths in order
    for path in ("entry/sample/name", "entry/title", "entry/sample_name"):
        if path in f:
            try:
                v = f[path][()]
                if isinstance(v, bytes):
                    v = v.decode("utf-8", errors="replace")
                s = str(v).strip()
                if s:
                    return s
            except Exception:
                pass
    # Try a root attribute
    for attr in ("sample", "sample_name", "title"):
        if attr in f.attrs:
            try:
                v = f.attrs[attr]
                if isinstance(v, bytes):
                    v = v.decode("utf-8", errors="replace")
                s = str(v).strip()
                if s:
                    return s
            except Exception:
                pass
    # The data sub-entry name itself is often the sample name
    if "entry" in f:
        entry = f["entry"]
        for child in entry:
            obj = entry[child]
            if isinstance(obj, h5py.Group) and \
               obj.attrs.get("canSAS_class") == b"SASentry" or \
               obj.attrs.get("canSAS_class") == "SASentry":
                return str(child)
    return None


def _analyses_in_file(f: h5py.File) -> list[str]:
    """List the tool keys whose result group is present in *f*."""
    return [key for key, _ in available_tools(f)]


def _file_entry(path: Path, deep: bool = True) -> FileEntry:
    """Build a FileEntry. When *deep* is True, opens the file for analysis list."""
    stat = path.stat()
    mtime_iso = datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat()
    sample: Optional[str] = None
    analyses: list[str] = []
    if deep:
        try:
            with h5py.File(path, "r") as f:
                sample = _sample_name(f)
                analyses = _analyses_in_file(f)
        except (OSError, KeyError):
            # Not a valid HDF5 file — leave fields empty
            pass
    return FileEntry(
        path=str(path),
        name=path.name,
        sample=sample,
        scan_number=_scan_number_from_name(path.stem),
        mtime=mtime_iso,
        size_bytes=stat.st_size,
        analyses=analyses,
    )


def list_files(
    folder: str,
    pattern: str = ",".join(_DEFAULT_GLOBS),
    sort: str = "mtime_desc",
    limit: int = 100,
    deep: bool = True,
) -> list[dict]:
    """List HDF5 files in *folder* with metadata.

    Parameters
    ----------
    folder : str
        Folder to scan (non-recursive).
    pattern : str
        Glob pattern, or comma-separated patterns. Default: all HDF5 extensions.
    sort : str
        One of name_asc, name_desc, mtime_asc, mtime_desc, size_asc, size_desc.
    limit : int
        Max number of entries to return. 0 = no limit.
    deep : bool
        When True (default), opens each file to read sample name and detected
        analyses. Set False for faster shallow listing (path/size/mtime only).

    Returns
    -------
    list of dict — each is a FileEntry.to_dict().
    """
    if sort not in _VALID_SORT:
        raise ValueError(f"sort must be one of {sorted(_VALID_SORT)}")
    folder_p = resolve_safe_folder(folder)
    matches = _glob_many(folder_p, pattern)

    key_funcs = {
        "name_asc":   (lambda p: p.name.lower(), False),
        "name_desc":  (lambda p: p.name.lower(), True),
        "mtime_asc":  (lambda p: p.stat().st_mtime, False),
        "mtime_desc": (lambda p: p.stat().st_mtime, True),
        "size_asc":   (lambda p: p.stat().st_size, False),
        "size_desc":  (lambda p: p.stat().st_size, True),
    }
    keyfn, reverse = key_funcs[sort]
    matches.sort(key=keyfn, reverse=reverse)
    if limit and limit > 0:
        matches = matches[:limit]

    return [_file_entry(p, deep=deep).to_dict() for p in matches]


def summarize_folder(
    folder: str,
    pattern: str = ",".join(_DEFAULT_GLOBS),
    sample_filter: Optional[str] = None,
) -> dict:
    """Aggregate snapshot of a folder.

    Returns counts of files, distinct samples, and per-analysis file counts.
    Cheap orientation call: an AI agent should call this before drilling in.

    Parameters
    ----------
    folder : str
        Folder to scan (non-recursive).
    pattern : str
        Glob, comma-separated permitted.
    sample_filter : str, optional
        If set, only count files whose sample name matches (case-insensitive
        substring).
    """
    folder_p = resolve_safe_folder(folder)
    matches = _glob_many(folder_p, pattern)
    samples: set[str] = set()
    analyses_count: dict[str, int] = {}
    mtimes: list[float] = []
    n_files = 0
    sample_filter_lc = sample_filter.lower() if sample_filter else None

    for path in matches:
        try:
            with h5py.File(path, "r") as f:
                sample = _sample_name(f)
                if sample_filter_lc:
                    if not sample or sample_filter_lc not in sample.lower():
                        continue
                if sample:
                    samples.add(sample)
                for key in _analyses_in_file(f):
                    analyses_count[key] = analyses_count.get(key, 0) + 1
                mtimes.append(path.stat().st_mtime)
                n_files += 1
        except (OSError, KeyError):
            continue

    mtime_min = mtime_max = None
    if mtimes:
        mtime_min = datetime.fromtimestamp(min(mtimes), tz=timezone.utc).isoformat()
        mtime_max = datetime.fromtimestamp(max(mtimes), tz=timezone.utc).isoformat()

    return FolderSummary(
        folder=str(folder_p),
        n_files=n_files,
        samples=sorted(samples),
        analyses_count=dict(sorted(analyses_count.items())),
        mtime_min=mtime_min,
        mtime_max=mtime_max,
        sample_filter_applied=sample_filter,
    ).to_dict()


def inspect_file(path: str) -> dict:
    """Single-file inspection: metadata + analyses present + Q range.

    Always returns a FileInspection — `analyses_present` is empty for files
    that are not valid HDF5 or contain no recognised analyses.
    """
    file_p = resolve_safe_file(path)
    stat = file_p.stat()
    mtime_iso = datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat()
    inspection = FileInspection(
        path=str(file_p),
        name=file_p.name,
        mtime=mtime_iso,
        size_bytes=stat.st_size,
        scan_number=_scan_number_from_name(file_p.stem),
    )

    try:
        with h5py.File(file_p, "r") as f:
            inspection.sample = _sample_name(f)
            inspection.analyses_present = _analyses_in_file(f)
            # Reduced data Q range: cheaply read by sniffing the default NXcanSAS data
            try:
                from pyirena.io.hdf5 import find_matching_groups
                entries = find_matching_groups(
                    f,
                    {"canSAS_class": "SASentry", "NX_class": "NXsubentry"},
                    {"definition": "NXcanSAS"},
                )
                if entries:
                    cur = entries[0]
                    default = f[cur].attrs.get("default")
                    if default is not None:
                        cur = f"{cur}/{default}".strip("/")
                    if cur in f:
                        grp = f[cur]
                        signal = grp.attrs.get("signal")
                        i_axes = grp.attrs.get("I_axes")
                        if signal is not None and i_axes is not None:
                            i_ds = f"{cur}/{signal}"
                            q_ds = f"{cur}/{i_axes}"
                            if i_ds in f and q_ds in f:
                                q_arr = f[q_ds][:]
                                inspection.reduced_data_present = True
                                inspection.n_points = int(q_arr.size)
                                if q_arr.size:
                                    inspection.q_min = float(q_arr.min())
                                    inspection.q_max = float(q_arr.max())
            except Exception:
                pass

            # A few interesting root attrs (creator, instrument)
            for key in ("creator", "instrument", "file_time"):
                if key in f.attrs:
                    v = f.attrs[key]
                    if isinstance(v, bytes):
                        v = v.decode("utf-8", errors="replace")
                    inspection.extra_metadata[key] = str(v)
    except OSError:
        pass

    return inspection.to_dict()
