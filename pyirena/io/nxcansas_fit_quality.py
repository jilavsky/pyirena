"""
Shared NeXus/HDF5 persistence for robust fit-quality metrics.

Every fit tool stores its results in its own ``entry/<tool>_results`` NXprocess
group (see the ``nxcansas_*.py`` modules). To keep the new fit-quality
diagnostics uniform across all of them, this module provides a single
writer/reader pair that serializes the dict returned by
:func:`pyirena.core.fit_metrics.fit_quality_metrics` into a ``fit_quality``
sub-group, and reads it back.

Design notes
------------
- The metrics are computed at fit time (where the consistent (I, M, sigma) basis
  is known) and passed in; this module only serializes. That avoids the
  per-tool basis ambiguities you would hit trying to recompute from stored
  arrays (e.g. sizes stores raw data + background-subtracted model; modeling
  stores no data arrays at all).
- Only scalars + the per-band table are stored. The per-point arrays
  (norm_residual, frac_residual, …) are intentionally *not* duplicated here —
  the residuals array already lives in the parent group and the rest is cheap to
  re-derive if ever needed.
- ``None`` scalars are simply omitted; the reader returns ``None`` for anything
  absent. Backward compatible: older files without a ``fit_quality`` group make
  :func:`read_fit_quality` return ``None``.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

__all__ = ["write_fit_quality", "read_fit_quality", "QUALITY_GROUP"]

QUALITY_GROUP = "fit_quality"

# Scalar fields written when present and finite.
_FLOAT_FIELDS = (
    "reduced_chi2",
    "robust_scale_s",
    "sigma_misscale_factor",
    "realistic_reduced_chi2_floor",
    "median_frac_uncertainty",
    "max_abs_frac_misfit",
    "q_at_max_frac_misfit",
    "frac_outliers_3s",
    "sign_autocorr_lag1",
)
_INT_FIELDS = (
    "n_valid",
    "n_params",
    "dof",
    "n_outliers_3s",
    "longest_same_sign_run",
    "n_bands_used",
)
_BAND_FIELDS = ("q_lo", "q_hi", "n", "reduced_chi2", "robust_scale_s", "max_abs_frac_misfit")


def _is_num(x) -> bool:
    return x is not None and isinstance(x, (int, float, np.integer, np.floating)) and np.isfinite(x)


def write_fit_quality(parent_grp, metrics: Optional[dict], group_name: str = QUALITY_GROUP) -> bool:
    """Serialize a fit_quality_metrics dict under ``parent_grp/<group_name>``.

    Parameters
    ----------
    parent_grp : h5py.Group
        The tool's ``entry/<tool>_results`` group.
    metrics : dict or None
        The dict returned by ``fit_quality_metrics``. If None, nothing is
        written and the function returns False.
    group_name : str
        Sub-group name (default ``"fit_quality"``).

    Returns
    -------
    bool
        True if a group was written.
    """
    if not metrics:
        return False

    if group_name in parent_grp:
        del parent_grp[group_name]
    grp = parent_grp.create_group(group_name)
    grp.attrs["NX_class"] = "NXcollection"
    grp.attrs["description"] = (
        "Robust, sigma-scale-independent fit-quality diagnostics "
        "(see pyirena.core.fit_metrics)."
    )

    # sigma_available as an explicit 0/1 flag
    grp.create_dataset("sigma_available", data=np.int8(1 if metrics.get("sigma_available") else 0))

    for key in _INT_FIELDS:
        val = metrics.get(key)
        if val is not None and np.isfinite(val):
            grp.create_dataset(key, data=int(val))
    for key in _FLOAT_FIELDS:
        val = metrics.get(key)
        if _is_num(val):
            grp.create_dataset(key, data=float(val))

    # Per-band table as parallel arrays (NaN where a band value is None).
    bands = metrics.get("bands") or []
    if bands:
        bgrp = grp.create_group("bands")
        bgrp.attrs["NX_class"] = "NXcollection"
        for field in _BAND_FIELDS:
            col = []
            for b in bands:
                v = b.get(field)
                col.append(float(v) if _is_num(v) else np.nan)
            bgrp.create_dataset(field, data=np.asarray(col, dtype=float))
    return True


def read_fit_quality(parent_grp, group_name: str = QUALITY_GROUP) -> Optional[dict]:
    """Read a ``fit_quality`` sub-group back into a dict.

    Returns ``None`` if the group is absent (e.g. files written before this
    feature existed). Scalars absent from the file come back as ``None``.
    """
    if group_name not in parent_grp:
        return None
    grp = parent_grp[group_name]

    def _get(key):
        if key in grp:
            return grp[key][()]
        return None

    out: dict = {}
    sa = _get("sigma_available")
    out["sigma_available"] = bool(sa) if sa is not None else None
    for key in _INT_FIELDS:
        v = _get(key)
        out[key] = int(v) if v is not None else None
    for key in _FLOAT_FIELDS:
        v = _get(key)
        out[key] = float(v) if v is not None else None

    bands = []
    if "bands" in grp:
        bgrp = grp["bands"]
        cols = {f: bgrp[f][()] for f in _BAND_FIELDS if f in bgrp}
        n_rows = len(next(iter(cols.values()))) if cols else 0
        for i in range(n_rows):
            row = {}
            for f in _BAND_FIELDS:
                if f in cols:
                    val = cols[f][i]
                    if f == "n":
                        row[f] = int(val) if np.isfinite(val) else 0
                    else:
                        row[f] = float(val) if np.isfinite(val) else None
            bands.append(row)
    out["bands"] = bands
    return out
