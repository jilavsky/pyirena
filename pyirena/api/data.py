"""Raw I(Q) reading and sample metadata for the api surface."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from pyirena.api._paths import resolve_safe_file
from pyirena.api.schemas import ReducedData, SampleMetadata, array_to_list
from pyirena.io.hdf5 import readGenericNXcanSAS


def read_reduced_data(
    path: str,
    decimate: Optional[int] = None,
    include_full: bool = False,
) -> dict:
    """Read the raw reduced I(Q) curve from a SAS file.

    Wraps :func:`pyirena.io.hdf5.readGenericNXcanSAS`. Arrays are decimated
    to *decimate* points (or PYIRENA_MAX_ARRAY_POINTS env default) so the
    response stays small. Pass ``include_full=True`` to keep full fidelity —
    use sparingly when feeding the result back to an LLM, since long arrays
    bloat context.

    Returns
    -------
    dict with keys: path, found, n_points, q_min, q_max, I_units, label,
    Q, I, dI, dQ, decimated, decimated_from.
    """
    file_p = resolve_safe_file(path)
    result = ReducedData(path=str(file_p))
    try:
        data = readGenericNXcanSAS(str(file_p.parent), file_p.name)
    except Exception:
        return result.to_dict()
    if data is None:
        return result.to_dict()

    q = data.get("Q")
    intensity = data.get("Intensity")
    error = data.get("Error")
    dq = data.get("dQ")
    if q is None or intensity is None:
        return result.to_dict()

    q_arr = np.asarray(q)
    n_orig = int(q_arr.size)
    if n_orig == 0:
        return result.to_dict()

    result.found = True
    result.n_points = n_orig
    result.q_min = float(np.nanmin(q_arr))
    result.q_max = float(np.nanmax(q_arr))
    units = data.get("units")
    if isinstance(units, bytes):
        units = units.decode("utf-8", errors="replace")
    result.I_units = str(units) if units else None
    label = data.get("label")
    if isinstance(label, bytes):
        label = label.decode("utf-8", errors="replace")
    result.label = str(label) if label else None

    cap = decimate if decimate is not None else None
    result.Q = array_to_list(q_arr, max_points=cap, include_full=include_full)
    result.I = array_to_list(intensity, max_points=cap, include_full=include_full)
    result.dI = array_to_list(error, max_points=cap, include_full=include_full)
    result.dQ = array_to_list(dq, max_points=cap, include_full=include_full)

    decimated_len = len(result.Q) if result.Q is not None else 0
    result.decimated = (not include_full) and (decimated_len < n_orig)
    if result.decimated:
        result.decimated_from = n_orig
    return result.to_dict()


def read_metadata(path: str) -> dict:
    """Read sample / experiment metadata from an NXcanSAS HDF5 file.

    Always returns a SampleMetadata; missing fields are None.
    """
    file_p = resolve_safe_file(path)
    meta = SampleMetadata(path=str(file_p))

    try:
        with h5py.File(file_p, "r") as f:
            meta.found = True
            # Sample name from common locations
            from pyirena.api.discovery import _sample_name
            meta.sample_name = _sample_name(f)

            # Sub-entry attributes (thickness, blank, label live on the I dataset)
            try:
                data = readGenericNXcanSAS(str(file_p.parent), file_p.name)
                if data:
                    thick = data.get("thickness")
                    if thick is not None and not isinstance(thick, bytes):
                        try:
                            meta.thickness = float(thick)
                        except (TypeError, ValueError):
                            pass
                    blank = data.get("blankname")
                    if isinstance(blank, bytes):
                        blank = blank.decode("utf-8", errors="replace")
                    if blank:
                        meta.blank = str(blank)
                    label = data.get("label")
                    if isinstance(label, bytes):
                        label = label.decode("utf-8", errors="replace")
                    if label:
                        meta.label = str(label)
            except Exception:
                pass

            # Top-level file attributes
            for attr_key, meta_field in (
                ("instrument", "instrument"),
                ("file_time", "timestamp"),
            ):
                if attr_key in f.attrs:
                    v = f.attrs[attr_key]
                    if isinstance(v, bytes):
                        v = v.decode("utf-8", errors="replace")
                    setattr(meta, meta_field, str(v))

            # Extra attrs (creator, NeXus_version, h5py_version, etc.) — keep as strings
            for key in f.attrs:
                if key in {"default", "instrument", "file_time", "sample"}:
                    continue
                try:
                    v = f.attrs[key]
                    if isinstance(v, bytes):
                        v = v.decode("utf-8", errors="replace")
                    meta.extra[key] = str(v)
                except Exception:
                    continue
    except OSError:
        pass

    return meta.to_dict()
