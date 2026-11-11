"""Dataset sessions and model-independent data tools for the control API.

Split out of ``unified_fit.py`` (review of 09-11-2026, issue #38), where these
lived although every fitting tool uses them:

* **session lifecycle** — ``open_dataset``, ``open_dataset_from_data``,
  ``list_open_sessions``, ``close_session``, ``get_session_summary``;
* **fit Q range** — ``get_data_q_range``, ``get_fit_q_range``,
  ``set_fit_q_range``, ``reset_fit_q_range``;
* **local estimators** — ``fit_local_guinier``, ``fit_local_power_law``;
* **feature detection** — ``detect_features``.

No behaviour changed. ``pyirena.api.control`` and ``unified_fit`` re-export
every name, and the dispatcher keeps listing the Q-range / estimator /
feature tools under the ``unified`` category where agents already find them.
"""
from __future__ import annotations

import os
from typing import Optional

import numpy as np

from pyirena.api._paths import (
    PathSecurityError,
    resolve_safe_file,
)
from pyirena.api.control.errors import (
    make_error,
    no_session,
)
from pyirena.api.control.session import (
    Session,
    all_sessions,
    create_session,
    drop_session,
    fit_mask,
    get_session,
)

# get_parameter_uncertainties serves every tool and lives in uncertainty.py;
# re-exported here, where it used to be defined, so existing imports work.
from pyirena.api.control.uncertainty import get_parameter_uncertainties  # noqa: F401

# ---------------------------------------------------------------------------
# Category A — Session lifecycle
# ---------------------------------------------------------------------------

def open_dataset(file_path: str, use_slit_smeared: bool = False) -> dict:
    """Load a NXcanSAS HDF5 file and create a new fitting session.

    Parameters
    ----------
    use_slit_smeared : bool
        When True, load the file's slit-smeared (``_SMR``) dataset if present.
        Model smearing is then enabled automatically at the file-derived slit
        length when a Unified Fit model is created for the session.

    Returns
    -------
    dict with keys: session_id, summary (file, n_points, q_min, q_max,
    intensity_min, intensity_max, label, is_slit_smeared, slit_length,
    has_slit_smeared_entry).
    On failure returns an error dict.
    """
    from pyirena.io.hdf5 import file_has_smr_entry, readGenericNXcanSAS

    try:
        fp = resolve_safe_file(file_path)
    except PathSecurityError as exc:
        return make_error(
            str(exc),
            suggestion="Open a file inside PYIRENA_DATA_ROOT.",
            code="PATH_NOT_ALLOWED",
        )
    except (FileNotFoundError, IsADirectoryError):
        return make_error(
            f"File not found: '{file_path}'",
            suggestion="Check the path and try again.",
            code="FILE_NOT_FOUND",
        )

    try:
        data = readGenericNXcanSAS(
            str(fp.parent), fp.name, prefer_slit_smeared=use_slit_smeared)
    except Exception as exc:
        return make_error(
            f"Could not read '{file_path}': {exc}",
            code="READ_ERROR",
        )

    if data is None:
        return make_error(
            f"No readable SAS data found in '{file_path}'.",
            suggestion="Verify the file is a valid NXcanSAS HDF5 file.",
            code="NO_DATA",
        )

    q_raw = data.get("Q")
    I_raw = data.get("Intensity")
    q = np.asarray(q_raw if q_raw is not None else [], dtype=float)
    I = np.asarray(I_raw if I_raw is not None else [], dtype=float)
    err_raw = data.get("Error")
    err = np.asarray(err_raw, dtype=float) if err_raw is not None else None

    if len(q) == 0 or len(I) == 0:
        return make_error(
            f"Empty Q or Intensity array in '{file_path}'.",
            code="EMPTY_DATA",
        )

    label = data.get("label") or ""
    if isinstance(label, bytes):
        label = label.decode("utf-8", errors="replace")

    is_slit_smeared = bool(data.get("is_slit_smeared", False))
    slit_length = float(data.get("slit_length", 0.0) or 0.0)
    try:
        has_smr_entry = file_has_smr_entry(str(fp.parent), fp.name)
    except Exception:
        has_smr_entry = False

    session = create_session(file_path=str(fp), q=q, intensity=I, error=err,
                             label=str(label),
                             is_slit_smeared=is_slit_smeared,
                             slit_length=slit_length)

    valid = np.isfinite(q) & np.isfinite(I)
    return {
        "session_id": session.session_id,
        "summary": {
            "file":          str(fp),
            "label":         session.label,
            "n_points":      int(len(q)),
            "q_min":         float(np.nanmin(q[valid])),
            "q_max":         float(np.nanmax(q[valid])),
            "intensity_min": float(np.nanmin(I[valid & (I > 0)])) if np.any(valid & (I > 0)) else None,
            "intensity_max": float(np.nanmax(I[valid])),
            # Slit-smearing status so the AI advisor can explain/enable it.
            "is_slit_smeared":        is_slit_smeared,
            "slit_length":            slit_length,
            "has_slit_smeared_entry": has_smr_entry,
        },
    }


# ---------------------------------------------------------------------------
# Session lifecycle — creating a session from arrays (no file)
# ---------------------------------------------------------------------------
#
# A remote transport (pyirena/zmq/) has no shared filesystem with its caller,
# so the data arrives inside the request. Sessions created this way have
# file_path None; every reader of Session.file_path must cope with that, and
# the save_* tools refuse with NO_SOURCE_FILE unless given an output_path.

# A JSON array costs roughly 20 bytes per point, so 100k points is a ~6 MB
# request — already far beyond any real SAS curve (a USAXS+SAXS+WAXS merge is
# a few thousand points). The guard is against a malformed or runaway caller,
# not against legitimate data; raise PYIRENA_MAX_INPUT_POINTS if you have more.
_DEFAULT_MAX_INPUT_POINTS = 100_000

# Below this, every model in the package is under-determined and the fit
# metrics are meaningless.
_MIN_INPUT_POINTS = 5


def _max_input_points() -> int:
    raw = os.environ.get("PYIRENA_MAX_INPUT_POINTS", "").strip()
    if not raw:
        return _DEFAULT_MAX_INPUT_POINTS
    try:
        value = int(raw)
    except ValueError:
        return _DEFAULT_MAX_INPUT_POINTS
    return value if value > 0 else _DEFAULT_MAX_INPUT_POINTS


def open_dataset_from_data(
    q: list,
    intensity: list,
    error: Optional[list] = None,
    dq: Optional[list] = None,
    label: str = "",
    is_slit_smeared: bool = False,
    slit_length: float = 0.0,
    error_fraction: float = 0.05,
) -> dict:
    """Create a fitting session from in-memory arrays, with no file involved.

    The array-fed twin of :func:`open_dataset`, for callers that have the data
    but not a path the server can read — a remote service, a notebook holding
    arrays, a reduction pipeline that has not written a file yet. The returned
    ``session_id`` works with every other control tool exactly as one from
    ``open_dataset`` does, except that ``save_*`` needs an explicit
    ``output_path`` (there is no source file to write back into).

    The data get the fitting tools' shared point rule
    (:mod:`pyirena.core.fit_data`, issue #38): points with Q ≤ 0 or a
    non-finite Q or I are removed; **points with I ≤ 0 are kept** (real
    negative values of background-subtracted data — dropping only the
    negative half of the noise biases a fitted background); non-positive
    uncertainties are replaced by ``error_fraction × |I|``. Until 1.2.1
    I ≤ 0 points were removed here. The counts are reported back in
    ``summary["cleaning"]`` rather than applied silently.

    Parameters
    ----------
    q : list of float
        Q values in Å⁻¹. Need not be sorted; the session is sorted ascending.
    intensity : list of float
        Intensity, cm⁻¹ where the data is on an absolute scale.
    error : list of float, optional
        Intensity uncertainty. Omit it to fit unweighted, exactly as for a
        file with no uncertainty column — uncertainties are *not* synthesized
        when the argument is absent, only repaired when it is present.
    dq : list of float, optional
        Q resolution. Stored on the session for provenance; no fitting tool
        reads it yet.
    label : str
        Human-readable name for the curve, used in plots and reports.
    is_slit_smeared : bool
        Mark the data as slit-smeared; ``slit_length`` (Å⁻¹) is then applied
        to the model when one is created for this session.
    error_fraction : float
        Fraction of I used to repair non-positive uncertainties (default 0.05).

    Returns
    -------
    dict
        ``{session_id, summary}`` — the same shape as :func:`open_dataset`,
        with ``summary["file"] = None`` and an extra ``summary["cleaning"]``
        report. On failure, an error dict with one of the codes
        ``SHAPE_MISMATCH``, ``EMPTY_DATA``, ``TOO_FEW_POINTS``,
        ``TOO_MANY_POINTS``, ``BAD_VALUES`` or ``NO_VALID_POINTS``.
    """

    def _as_array(values, name):
        try:
            arr = np.asarray(values, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"'{name}' is not an array of numbers: {exc}") from None
        if arr.ndim != 1:
            raise ValueError(f"'{name}' must be a flat list, got {arr.ndim} dimensions.")
        return arr

    try:
        q_arr = _as_array(q, "q")
        I_arr = _as_array(intensity, "intensity")
        err_arr = _as_array(error, "error") if error is not None else None
        dq_arr = _as_array(dq, "dq") if dq is not None else None
    except ValueError as exc:
        return make_error(
            str(exc),
            suggestion="Pass plain lists of numbers, one value per data point.",
            code="BAD_VALUES",
        )

    if len(q_arr) == 0 or len(I_arr) == 0:
        return make_error(
            "Empty q or intensity array.",
            suggestion="Send at least "
                       f"{_MIN_INPUT_POINTS} points of real data.",
            code="EMPTY_DATA",
        )

    for name, arr in (("intensity", I_arr), ("error", err_arr), ("dq", dq_arr)):
        if arr is not None and len(arr) != len(q_arr):
            return make_error(
                f"'{name}' has {len(arr)} values but 'q' has {len(q_arr)}.",
                suggestion="All arrays must have one value per data point.",
                code="SHAPE_MISMATCH",
            )

    if len(q_arr) < _MIN_INPUT_POINTS:
        return make_error(
            f"Only {len(q_arr)} data points; at least {_MIN_INPUT_POINTS} are needed.",
            suggestion="Send the full measured curve.",
            code="TOO_FEW_POINTS",
        )

    max_points = _max_input_points()
    if len(q_arr) > max_points:
        return make_error(
            f"{len(q_arr)} data points exceeds the limit of {max_points}.",
            suggestion="Rebin the curve, or raise PYIRENA_MAX_INPUT_POINTS on the server.",
            code="TOO_MANY_POINTS",
        )

    # Sort ascending in Q. Unsorted input is common (concatenated USAXS + SAXS
    # + WAXS ranges) and silently breaks anything that assumes monotonic Q.
    was_sorted = bool(np.all(np.diff(q_arr) > 0))
    if not was_sorted:
        order = np.argsort(q_arr, kind="stable")
        q_arr = q_arr[order]
        I_arr = I_arr[order]
        err_arr = err_arr[order] if err_arr is not None else None
        dq_arr = dq_arr[order] if dq_arr is not None else None

    from pyirena.core.fit_data import prepare_fit_arrays
    prep = prepare_fit_arrays(q_arr, I_arr, err_arr, error_fraction)
    q_c, I_c, e_c = prep.q, prep.intensity, prep.error
    dq_c = dq_arr[prep.mask] if dq_arr is not None and len(dq_arr) == len(q_arr) else None
    n_q_removed = int(np.sum(~np.isfinite(q_arr) | (q_arr <= 0)))
    report = {
        "n_original": int(len(q_arr)),
        "n_kept": int(len(q_c)),
        "n_q_removed": n_q_removed,
        "n_i_removed": int(len(q_arr) - len(q_c) - n_q_removed),
        "n_e_synthesized": int(prep.n_sigma_repaired),
    }

    if len(q_c) < _MIN_INPUT_POINTS:
        return make_error(
            f"Only {len(q_c)} of {len(q_arr)} points survived cleaning "
            f"(removed {report['n_q_removed']} with Q<=0 or non-finite Q and "
            f"{report['n_i_removed']} with non-finite I or unusable σ).",
            suggestion="Check the Q and intensity units and sign conventions.",
            code="NO_VALID_POINTS",
        )

    # Uncertainties are repaired, never invented: a caller who sent no error
    # array gets an unweighted session, the same as a file with no error
    # column, rather than a silent 5%-of-I weighting they did not ask for.
    session_error = e_c if err_arr is not None else None

    session = create_session(
        file_path=None,
        q=q_c,
        intensity=I_c,
        error=session_error,
        label=str(label or ""),
        is_slit_smeared=bool(is_slit_smeared),
        slit_length=float(slit_length or 0.0),
        dq=dq_c,
    )
    session.error_fraction = float(error_fraction)

    return {
        "session_id": session.session_id,
        "summary": {
            "file":          None,
            "label":         session.label,
            "n_points":      int(len(q_c)),
            "q_min":         float(np.min(q_c)),
            "q_max":         float(np.max(q_c)),
            "intensity_min": float(np.min(I_c)),
            "intensity_max": float(np.max(I_c)),
            "has_errors":    session_error is not None,
            "is_slit_smeared":        bool(is_slit_smeared),
            "slit_length":            float(slit_length or 0.0),
            "has_slit_smeared_entry": False,
            "sorted_by_q":   not was_sorted,
            "cleaning": {
                "n_input":         int(report["n_original"]),
                "n_kept":          int(report["n_kept"]),
                "n_removed_q":     int(report["n_q_removed"]),
                "n_removed_i":     int(report["n_i_removed"]),
                "n_repaired_error": (
                    int(report["n_e_synthesized"]) if err_arr is not None else 0
                ),
            },
        },
    }


def list_open_sessions() -> dict:
    """Return a summary of all open sessions."""
    sessions = all_sessions()
    return {
        "sessions": [
            {
                "session_id": s.session_id,
                "file":       s.file_path,
                "label":      s.label,
                "model":      s.model_name,
                "has_fit":    s.last_fit_result is not None,
            }
            for s in sessions
        ],
        "count": len(sessions),
    }


def close_session(session_id: str) -> dict:
    """Close and discard a session, freeing its memory."""
    if not drop_session(session_id):
        return no_session(session_id)
    return {"ok": True, "session_id": session_id}


def _reduced_chi_squared(last_fit_result) -> Optional[float]:
    """Reduced chi-squared from a fit result, dict or dataclass, or None."""
    if last_fit_result is None:
        return None
    if isinstance(last_fit_result, dict):
        value = last_fit_result.get("reduced_chi_squared")
    else:
        value = getattr(last_fit_result, "reduced_chi_squared", None)
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def get_session_summary(session_id: str) -> dict:
    """Return a summary of the current session state."""
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    valid = np.isfinite(s.q) & np.isfinite(s.intensity)
    # last_fit_result is a plain dict for most tools but a result dataclass
    # for Modeling and the Carbon model, so read it both ways rather than
    # assuming .get() exists.
    chi = _reduced_chi_squared(s.last_fit_result)

    return {
        "session_id":    session_id,
        "file":          s.file_path,
        "label":         s.label,
        "n_points":      int(len(s.q)),
        "q_min":         float(np.nanmin(s.q[valid])),
        "q_max":         float(np.nanmax(s.q[valid])),
        "model":         s.model_name,
        # nlevels only applies to unified_fit; other models (e.g. sizes) lack it
        "nlevels":       getattr(s.model, "num_levels", None) if s.model is not None else None,
        "fit_q_min":     s.fit_q_min,
        "fit_q_max":     s.fit_q_max,
        "has_fit":       s.last_fit_result is not None,
        "reduced_chi_squared": chi,
    }


# ---------------------------------------------------------------------------
# Category C'' — Q range
# ---------------------------------------------------------------------------

def get_data_q_range(session_id: str) -> dict:
    """Return the Q range of the loaded dataset."""
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    valid = np.isfinite(s.q)
    return {
        "q_min":    float(np.nanmin(s.q[valid])),
        "q_max":    float(np.nanmax(s.q[valid])),
        "n_points": int(len(s.q)),
    }


def get_fit_q_range(session_id: str) -> dict:
    """Return the Q range currently used for fitting."""
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    mask = fit_mask(s)
    n_fit = int(np.sum(mask))
    valid = np.isfinite(s.q)
    return {
        "q_min":         s.fit_q_min if s.fit_q_min is not None else float(np.nanmin(s.q[valid])),
        "q_max":         s.fit_q_max if s.fit_q_max is not None else float(np.nanmax(s.q[valid])),
        "n_fit_points":  n_fit,
        "is_full_range": s.fit_q_min is None and s.fit_q_max is None,
    }


def set_fit_q_range(
    session_id: str,
    q_min: Optional[float] = None,
    q_max: Optional[float] = None,
) -> dict:
    """Restrict the Q range used for fitting.

    Pass q_min and/or q_max to set limits.  Values are clamped to the data
    range.  Pass None to leave that end unchanged.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    valid = np.isfinite(s.q)
    data_qmin = float(np.nanmin(s.q[valid]))
    data_qmax = float(np.nanmax(s.q[valid]))

    if q_min is not None:
        s.fit_q_min = float(max(q_min, data_qmin))
    if q_max is not None:
        s.fit_q_max = float(min(q_max, data_qmax))

    # Sanity: ensure min < max
    lo = s.fit_q_min if s.fit_q_min is not None else data_qmin
    hi = s.fit_q_max if s.fit_q_max is not None else data_qmax
    if lo >= hi:
        s.fit_q_min = None
        s.fit_q_max = None
        return make_error(
            f"Resulting Q range [{lo}, {hi}] is empty — range reset to full data range.",
            code="EMPTY_RANGE",
        )

    return get_fit_q_range(session_id)


def reset_fit_q_range(session_id: str) -> dict:
    """Restore the fit Q range to the full data range."""
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    s.fit_q_min = None
    s.fit_q_max = None
    return get_fit_q_range(session_id)


# ---------------------------------------------------------------------------
# Category C'''' — Local estimators (Guinier / power law on a Q sub-range)
# ---------------------------------------------------------------------------
#
# Standalone Q-range fits of single-term models, equivalent to the GUI's
# "Fit Rg/G btwn cursors" and "Fit P/B btwn cursors" buttons.  Useful for
# generating good starting values before the full multi-level fit runs.
#
# These tools do NOT mutate the model — they only return the fitted
# numbers.  The agent decides whether to apply them via set_parameter_value().

def _slice_q_range(session: Session, q_min: float, q_max: float):
    """Return q, I, err arrays restricted to [q_min, q_max].  May be empty."""
    mask = (session.q >= q_min) & (session.q <= q_max)
    if session.error is not None:
        return session.q[mask], session.intensity[mask], session.error[mask]
    return session.q[mask], session.intensity[mask], None


def fit_local_guinier(
    session_id: str,
    q_min: float,
    q_max: float,
) -> dict:
    """Fit a single Guinier term I(q) = G · exp(-q²·Rg²/3) on a Q sub-range.

    Equivalent to the GUI's "Fit Rg/G btwn cursors" button.  Useful for
    estimating starting values for Rg and G on one level before running
    the full multi-level fit.

    Does **not** modify the session's model — call set_parameter_value()
    afterwards if you want to apply the result.

    Parameters
    ----------
    q_min, q_max : float
        Q range to fit in [Å⁻¹].  Choose a window covering the Guinier
        knee of the level you are characterising.

    Returns
    -------
    dict with keys:
      - G, Rg                : fitted values
      - q_min, q_max         : actual Q range used (may be clipped to data)
      - n_points             : number of data points used
      - chi_squared          : sum of squared (weighted) residuals
      - reduced_chi_squared  : χ² per degree of freedom
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    if q_max <= q_min:
        return make_error(
            f"q_max ({q_max}) must be greater than q_min ({q_min}).",
            code="BAD_RANGE",
        )

    q_fit, I_fit, err_fit = _slice_q_range(s, float(q_min), float(q_max))
    if len(q_fit) < 3:
        return make_error(
            f"Not enough data points in [{q_min}, {q_max}] ({len(q_fit)} points). "
            "Need at least 3 points; widen the Q range.",
            code="TOO_FEW_POINTS",
        )

    from pyirena.core.unified import fit_local_guinier as _core_guinier  # noqa: PLC0415

    try:
        r = _core_guinier(q_fit, I_fit, error=err_fit,
                          slit_length=s.slit_length if s.is_slit_smeared else 0.0)
    except Exception as exc:
        return make_error(
            f"Local Guinier fit failed: {exc}",
            suggestion="Try a different Q range that covers a visible Guinier knee.",
            code="FIT_FAILED",
        )

    out = {
        "ok":                  True,
        "G":                   r["G"],
        "Rg":                  r["Rg"],
        "q_min":               float(q_fit[0]),
        "q_max":               float(q_fit[-1]),
        "n_points":            r["n_points"],
        "chi_squared":         r["chi_squared"],
        "reduced_chi_squared": r["reduced_chi_squared"],
    }
    if r.get("warning"):
        out["warning"] = r["warning"]
    return out


def fit_local_power_law(
    session_id: str,
    q_min: float,
    q_max: float,
) -> dict:
    """Fit a single power-law term I(q) = B · q⁻ᴾ on a Q sub-range.

    Equivalent to the GUI's "Fit P/B btwn cursors" button.  Useful for
    estimating starting values for P and B on one level (typically the
    high-Q power-law region after the Guinier knee).

    Does **not** modify the session's model — call set_parameter_value()
    afterwards if you want to apply the result.

    Parameters
    ----------
    q_min, q_max : float
        Q range to fit in [Å⁻¹].  Pick a window over the linear power-law
        portion of the log-log plot, beyond the Guinier knee.

    Returns
    -------
    dict with keys:
      - P, B                 : fitted values
      - q_min, q_max         : actual Q range used (may be clipped to data)
      - n_points             : number of data points used
      - chi_squared          : sum of squared (weighted) residuals
      - reduced_chi_squared  : χ² per degree of freedom
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    if q_max <= q_min:
        return make_error(
            f"q_max ({q_max}) must be greater than q_min ({q_min}).",
            code="BAD_RANGE",
        )

    q_fit, I_fit, err_fit = _slice_q_range(s, float(q_min), float(q_max))
    if len(q_fit) < 3:
        return make_error(
            f"Not enough data points in [{q_min}, {q_max}] ({len(q_fit)} points). "
            "Need at least 3 points; widen the Q range.",
            code="TOO_FEW_POINTS",
        )

    # Reject non-positive intensities (log-slope estimate would crash)
    if int(np.sum((I_fit > 0) & (q_fit > 0))) < 3:
        return make_error(
            "Not enough positive intensity values in range for a power-law fit.",
            code="TOO_FEW_POSITIVE",
        )

    from pyirena.core.unified import fit_local_power_law as _core_power_law  # noqa: PLC0415

    try:
        r = _core_power_law(q_fit, I_fit, error=err_fit,
                           slit_length=s.slit_length if s.is_slit_smeared else 0.0)
    except Exception as exc:
        return make_error(
            f"Local power-law fit failed: {exc}",
            suggestion="Try a different Q range over the linear power-law region.",
            code="FIT_FAILED",
        )

    q_pos = r["q"]
    out = {
        "ok":                  True,
        "P":                   r["P"],
        "B":                   r["B"],
        "q_min":               float(q_pos[0]),
        "q_max":               float(q_pos[-1]),
        "n_points":            r["n_points"],
        "chi_squared":         r["chi_squared"],
        "reduced_chi_squared": r["reduced_chi_squared"],
    }
    if r.get("warning"):
        out["warning"] = r["warning"]
    return out


# ---------------------------------------------------------------------------
# Category C''''' — Feature detection (slope-profile analysis)
# ---------------------------------------------------------------------------
#
# Analyses the loaded data in log-log space and identifies plateaus
# (Guinier knees), peaks (Guinier + structure factor), and power-law
# regions.  Intended to help the agent decide how many Unified Fit levels
# a curve needs and where to place Q-windows for local Guinier/Porod fits.
# Does not mutate the model.

def detect_features(
    session_id: str,
    q_min: Optional[float] = None,
    q_max: Optional[float] = None,
    q_max_clip: Optional[float] = 0.6,
    config_overrides: Optional[dict] = None,
) -> dict:
    """Segment the loaded I(Q) curve into power-law slope segments.

    Returns a piecewise classification of the curve in log-log space:
    each segment has a locally-constant slope; adjacent segments with
    different slopes imply Guinier knees between them (also returned).
    Intended to help the user / agent decide how many Unified Fit levels
    are needed and where to place Q-windows for ``fit_local_guinier`` /
    ``fit_local_power_law``.

    Does **not** modify the session's model or fit state.

    Parameters
    ----------
    session_id : str
        Session whose data should be analysed.
    q_min, q_max : float, optional
        Restrict analysis to this Q sub-range.  Defaults to full data range.
    q_max_clip : float, optional
        Drop data with Q greater than this value before segmentation.
        Default 0.6 Å⁻¹ — the practical upper limit of the small-angle
        approximation; above this, amorphous diffraction features should
        not be treated as SAS structure.  Pass ``None`` to disable.
    config_overrides : dict, optional
        Mapping of ``FeatureDetectConfig`` field names to override values
        (e.g. ``{"change_threshold_1": 0.60, "min_segment_decades": 0.15}``).

    Returns
    -------
    dict with keys:
      - ok                            : True on success
      - segments                      : list of {q_min, q_max, P, P_std, kind,
                                                 intensity_mid, width_decades}
        Ordered HIGH-Q → LOW-Q to match Unified Fit level numbering
        (segments[0] = Level 1, highest Q, smallest structure).
        kind is "background" | "guinier_plateau" | "power_law"
      - guinier_knees                 : list of {q_min, q_max, q_center,
                                                 P_low_q, P_high_q, delta_P}
        P_low_q and P_high_q are positive Porod exponents (I ∝ Q^-P).
        Only transitions where P_low_q < P_high_q are listed — exponent
        is smaller (shallower) at low Q, the physical Guinier-knee condition.
      - segments                      : each segment also carries "P" and "P_std"
        (positive Porod exponent, I ∝ Q^-P) in addition to "slope" / "slope_std".
      - background_q_min              : float | null
      - recommended_guinier_windows   : list of {feature_type, q_min_guinier,
                                                 q_max_guinier, q_min_powerlaw}
      - recommended_nlevels           : suggested number of Unified Fit levels
      - n_segments_found              : int
      - log_decades                   : total log10(Q) span analysed
      - q_min_analysed, q_max_analysed : actual Q range analysed
      - n_points                      : number of data points used
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    from pyirena.core.feature_detect import (  # noqa: PLC0415
        FeatureDetectConfig,
    )
    from pyirena.core.feature_detect import (
        detect_features as _detect,
    )

    if q_min is None:
        q_min = float(s.q.min())
    if q_max is None:
        q_max = float(s.q.max())
    if q_max <= q_min:
        return make_error(
            f"q_max ({q_max}) must be greater than q_min ({q_min}).",
            code="BAD_RANGE",
        )

    q_use, I_use, err_use = _slice_q_range(s, float(q_min), float(q_max))
    if q_use.size < 8:
        return make_error(
            f"Not enough data points in [{q_min}, {q_max}] ({q_use.size} points). "
            "Need at least 8 points; widen the Q range.",
            code="TOO_FEW_POINTS",
        )

    cfg = FeatureDetectConfig(q_max_clip=q_max_clip)
    if config_overrides:
        for k, v in config_overrides.items():
            if not hasattr(cfg, k):
                return make_error(
                    f"Unknown FeatureDetectConfig field: {k!r}",
                    suggestion="Valid fields include change_window_1, "
                               "change_threshold_1, change_window_2, "
                               "change_threshold_2, wide_region_decades, "
                               "min_segment_decades, edge_min_segment_decades, "
                               "merge_slope_tol, guinier_knee_min_delta_slope.",
                    code="BAD_OVERRIDE",
                )
            try:
                if v is None:
                    setattr(cfg, k, None)
                else:
                    setattr(cfg, k, float(v))
            except (TypeError, ValueError):
                return make_error(
                    f"Invalid value for config field {k!r}: {v!r}",
                    code="BAD_OVERRIDE",
                )

    result = _detect(q_use, I_use, sigma_I=err_use, config=cfg)
    out = result.to_dict()
    out["ok"] = True
    return out
