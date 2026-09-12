"""pyirena.api.data_ops — Data Manipulation and Data Merge for AI/scripting.

This is the **only** api group that creates data files. Everything else in
``pyirena.api`` reads; these functions average, subtract, divide, scale,
trim, rebin and merge SAS datasets, and write the result as a new
NXcanSAS HDF5 file.

Output location follows the convention the GUI and the batch layer already
use, so an agent's output lands where a user's own scripting would put it:

    <source folder>_manip/<stem><operation suffix><ext>     manipulation
    <source folder>_merged/<stem>_merged<ext>               merge

i.e. a **sibling** of the source folder, and a per-operation filename
suffix (``_avg``, ``_sub``, ``_div``, ``_scaled``, ``_trimmed``,
``_rebinned``, ``_merged``). Repeating the same operation on the same
input overwrites the previous result; different operations do not collide.
Pass ``output_folder`` to override.

Why this wraps ``core`` + ``io`` rather than ``pyirena.batch``, which does
the same orchestration:

1. ``pyirena.batch`` calls ``ensure_console_output()``, which attaches a
   ``StreamHandler(sys.stdout)`` when no handler is configured. MCP speaks
   JSON-RPC over stdout, so that would corrupt the transport.
2. The batch functions return bare ``None`` on every failure — bad
   operation, load failure, missing buffer, engine error, save error —
   which an agent cannot act on. They also return ``Path`` objects and raw
   numpy arrays, neither of which is JSON-safe.
3. CLAUDE.md's layer table allows the api layer to import ``core`` and
   ``io``; ``batch`` is a sibling layer, not a dependency.

Errors are **returned** as ``{"error", "suggestion", "code"}`` dicts rather
than raised, matching ``pyirena.api.control`` and ``pyirena.api.calculators``.

Q is in 1/Å and intensity in 1/cm (absolute) where the source file is
calibrated; these operations preserve whatever the input carries.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import numpy as np

from pyirena.api._paths import PathSecurityError, resolve_safe, resolve_safe_file

__all__ = [
    "average_data",
    "subtract_data",
    "divide_data",
    "scale_data",
    "trim_data",
    "rebin_data",
    "merge_datasets",
    "match_merge_files",
]

REBIN_MODES: tuple[str, ...] = ("log", "linear", "reference")
SIMILARITY_REFERENCES: tuple[str, ...] = ("first", "majority")


# ---------------------------------------------------------------------------
# Error helpers — same shape as pyirena.api.control.errors.make_error
# ---------------------------------------------------------------------------

def _error(message: str, suggestion: str = "", code: str = "ERROR") -> dict:
    return {"error": message, "suggestion": suggestion, "code": code}


def _path_error(exc: PathSecurityError) -> dict:
    return _error(
        str(exc),
        suggestion=(
            "Pass an output_folder inside PYIRENA_DATA_ROOT. The default "
            "output folder is a sibling of the source folder, which can "
            "fall outside the root when the root is the data folder itself."
        ),
        code="PATH_NOT_ALLOWED",
    )


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _load(path: str) -> tuple[Optional[dict], Optional[dict]]:
    """Load one dataset. Returns (data, None) or (None, error_dict).

    Mirrors pyirena.batch._common._load_data (text files are converted to a
    cached NXcanSAS sibling first) but reports failures instead of
    returning None, and never touches the logging configuration.
    """
    from pyirena.io.hdf5 import readGenericNXcanSAS

    try:
        file_p = resolve_safe_file(path)
    except PathSecurityError as exc:
        return None, _path_error(exc)
    except (FileNotFoundError, IsADirectoryError) as exc:
        return None, _error(
            str(exc),
            suggestion="Check the path with list_files() or inspect_file().",
            code="FILE_NOT_FOUND",
        )

    actual = file_p
    try:
        if file_p.suffix.lower() in (".txt", ".dat"):
            # Same conversion the GUI and batch use, including the
            # configured Q-unit assumption for text files.
            from pyirena.io.text_import import ensure_nxcansas_sibling
            from pyirena.state.state_manager import StateManager

            q_unit = StateManager().get("data_selector", "q_unit", "1/A")
            actual = ensure_nxcansas_sibling(file_p, q_unit=q_unit)
        data = readGenericNXcanSAS(str(actual.parent), actual.name)
    except Exception as exc:
        return None, _error(
            f"Could not read '{file_p.name}': {exc}",
            suggestion="Confirm the file is a readable NXcanSAS HDF5 or text dataset.",
            code="READ_FAILED",
        )

    if data is None or data.get("Q") is None or data.get("Intensity") is None:
        return None, _error(
            f"No I(Q) data found in '{file_p.name}'.",
            suggestion="Use inspect_file() to see what the file contains.",
            code="NO_DATA",
        )

    q = np.asarray(data["Q"], dtype=float)
    intensity = np.asarray(data["Intensity"], dtype=float)
    if q.size < 2:
        return None, _error(
            f"'{file_p.name}' has fewer than 2 data points.",
            suggestion="Check the file — it may be truncated.",
            code="NO_DATA",
        )
    # The core interpolation helpers use np.interp, which silently returns
    # garbage on an unsorted grid. Nothing downstream checks this.
    if np.any(np.diff(q) < 0):
        return None, _error(
            f"Q is not in ascending order in '{file_p.name}'.",
            suggestion="Re-export the dataset with Q sorted ascending.",
            code="UNSORTED_Q",
        )

    error = data.get("Error")
    data["Q"] = q
    data["Intensity"] = intensity
    data["Error"] = (
        np.asarray(error, dtype=float) if error is not None else intensity * 0.05
    )

    # For a slit-smeared file, readGenericNXcanSAS returns dQ as the scalar
    # slit length rather than a per-point array, and the core operations
    # index it (dQ[mask]) — which raises on a 0-d value. The slit length is
    # carried separately in slit_length / is_slit_smeared and re-written by
    # the saver, so drop a non-conforming dQ rather than propagating it.
    dq = data.get("dQ")
    if dq is not None:
        dq_arr = np.asarray(dq, dtype=float)
        data["dQ"] = dq_arr if dq_arr.ndim == 1 and dq_arr.size == q.size else None
    data["_source_path"] = file_p
    data.setdefault("is_nxcansas", True)
    data.setdefault("slit_length", 0.0)
    data.setdefault("is_slit_smeared", False)
    return data, None


def _smear(data: dict) -> tuple[bool, float]:
    """Slit-smearing flags for one dataset.

    These MUST be plumbed into subtract/divide: the core kwargs default to
    False/0.0, so forgetting them makes the compatibility check pass and
    produces a physically meaningless result with no warning.
    """
    return (
        bool(data.get("is_slit_smeared", False)),
        float(data.get("slit_length", 0.0) or 0.0),
    )


# ---------------------------------------------------------------------------
# Output paths
# ---------------------------------------------------------------------------

def _resolve_output_folder(
    output_folder: Optional[str], source: Path, suffix: str
) -> tuple[Optional[Path], Optional[dict]]:
    """Resolve the output folder, defaulting to the sibling batch/GUI use.

    The derived sibling is validated too: a sibling of an in-root input can
    itself be outside PYIRENA_DATA_ROOT when the root is the data folder.
    """
    raw = output_folder if output_folder else str(source.parent) + suffix
    try:
        return resolve_safe(raw, must_exist=False), None
    except PathSecurityError as exc:
        return None, _path_error(exc)


def _base_result(operation: str, out_path: Path, sources: list[Path]) -> dict:
    return {
        "success": True,
        "operation": operation,
        "output_path": str(out_path),
        "output_folder": str(out_path.parent),
        "output_name": out_path.name,
        "source_files": [str(p) for p in sources],
    }


def _diagnostics(q_in_total: int, result_q, result_I, out_path: Path) -> dict:
    """Point accounting, so a caller can see silent truncation and stripping.

    Two things quietly remove points: log-log interpolation does not
    extrapolate (out-of-range points become NaN and are masked), and the
    saver strips non-positive intensities before writing.
    """
    q_arr = np.asarray(result_q, dtype=float)
    i_arr = np.asarray(result_I, dtype=float)
    n_result = int(q_arr.size)

    finite = np.isfinite(q_arr) & np.isfinite(i_arr)
    n_nonfinite = int(np.count_nonzero(~finite))
    n_nonpositive = int(np.count_nonzero(finite & (i_arr <= 0)))

    n_written = 0
    q_min = q_max = None
    try:
        from pyirena.io.hdf5 import readGenericNXcanSAS

        written = readGenericNXcanSAS(str(out_path.parent), out_path.name)
        if written is not None and written.get("Q") is not None:
            wq = np.asarray(written["Q"], dtype=float)
            n_written = int(wq.size)
            if n_written:
                q_min, q_max = float(np.nanmin(wq)), float(np.nanmax(wq))
    except Exception:
        pass

    out = {
        "n_points_in": int(q_in_total),
        "n_points_result": n_result,
        "n_points_written": n_written,
        "n_dropped_nonpositive": n_nonpositive,
        "n_nonfinite": n_nonfinite,
        "q_min": q_min,
        "q_max": q_max,
    }
    notes = []
    if n_nonpositive:
        notes.append(
            f"{n_nonpositive} point(s) had intensity <= 0 and were stripped on "
            "save. For a subtraction this usually means over-subtraction — "
            "reduce buffer_scale."
        )
    if n_nonfinite:
        notes.append(
            f"{n_nonfinite} point(s) were non-finite (a zero denominator, or "
            "Q outside the other dataset's range — log-log interpolation does "
            "not extrapolate)."
        )
    if notes:
        out["warnings"] = notes
    return out


def _save_manip(
    out_folder: Path,
    source: dict,
    result: Any,
    slit_length: Optional[float] = None,
) -> tuple[Optional[Path], Optional[dict]]:
    from pyirena.io.nxcansas_data_manipulation import save_manipulated_data

    if slit_length is None:
        slit_length = float(
            result.metadata.get("slit_length", source.get("slit_length", 0.0)) or 0.0
        )
    try:
        out_path = save_manipulated_data(
            output_folder=out_folder,
            source_path=source["_source_path"],
            source_is_nxcansas=bool(source.get("is_nxcansas", True)),
            q=result.q, I=result.I, dI=result.dI, dQ=result.dQ,
            operation=result.operation,
            provenance=result.metadata,
            slit_length=slit_length,
        )
    except ValueError as exc:
        # save strips non-positive intensities and raises when <2 survive
        return None, _error(
            f"Nothing left to save: {exc}",
            suggestion=(
                "The result has fewer than 2 usable points. For a subtraction "
                "this usually means over-subtraction (lower buffer_scale); for "
                "a trim it means the Q window is too narrow."
            ),
            code="EMPTY_RESULT",
        )
    except Exception as exc:
        return None, _error(
            f"Could not write the output file: {exc}",
            suggestion="Check the output folder is writable.",
            code="SAVE_FAILED",
        )
    return out_path, None


# ---------------------------------------------------------------------------
# Single-dataset operations
# ---------------------------------------------------------------------------

def scale_data(
    file: str,
    scale_I: float = 1.0,
    background: float = 0.0,
    scale_uncertainty: Optional[float] = None,
    output_folder: Optional[str] = None,
) -> dict:
    """Scale a dataset's intensity and/or subtract a flat background.

    Computes ``I_out = scale_I * I - background`` — the background is
    subtracted *after* scaling. Writes ``<stem>_scaled.h5``.

    Parameters
    ----------
    file : str
        Input dataset.
    scale_I : float
        Multiplicative intensity factor.
    background : float
        Flat background in the intensity units of the file (1/cm when
        absolute), subtracted after scaling.
    scale_uncertainty : float, optional
        Factor applied to dI. Defaults to ``scale_I``.
    output_folder : str, optional
        Defaults to ``<source folder>_manip``.
    """
    from pyirena.core.data_manipulation import DataManipulation, ScaleConfig

    data, err = _load(file)
    if err:
        return err
    out_folder, err = _resolve_output_folder(output_folder, data["_source_path"], "_manip")
    if err:
        return err

    result = DataManipulation.scale(
        data["Q"], data["Intensity"], data["Error"], data.get("dQ"),
        ScaleConfig(
            scale_I=float(scale_I),
            background=float(background),
            scale_uncertainty=(
                None if scale_uncertainty is None else float(scale_uncertainty)
            ),
        ),
    )
    out_path, err = _save_manip(out_folder, data, result)
    if err:
        return err

    out = _base_result("scale", out_path, [data["_source_path"]])
    out["parameters"] = {
        "scale_I": float(scale_I),
        "background": float(background),
        "scale_uncertainty": (
            float(scale_uncertainty) if scale_uncertainty is not None else float(scale_I)
        ),
    }
    out.update(_diagnostics(data["Q"].size, result.q, result.I, out_path))
    return out


def trim_data(
    file: str,
    q_min: float = 0.0,
    q_max: Optional[float] = None,
    output_folder: Optional[str] = None,
) -> dict:
    """Keep only the points with q_min <= Q <= q_max. Writes ``<stem>_trimmed.h5``.

    Parameters
    ----------
    file : str
        Input dataset.
    q_min, q_max : float
        Q window in 1/Å, inclusive at both ends. ``q_max=None`` means no
        upper limit.
    output_folder : str, optional
        Defaults to ``<source folder>_manip``.
    """
    from pyirena.core.data_manipulation import DataManipulation, TrimConfig

    data, err = _load(file)
    if err:
        return err
    hi = float("inf") if q_max is None else float(q_max)
    if hi <= float(q_min):
        return _error(
            f"Empty Q window: q_min={q_min}, q_max={q_max}.",
            suggestion="Pass q_max greater than q_min, both in 1/Å.",
            code="BAD_Q_RANGE",
        )
    out_folder, err = _resolve_output_folder(output_folder, data["_source_path"], "_manip")
    if err:
        return err

    result = DataManipulation.trim(
        data["Q"], data["Intensity"], data["Error"], data.get("dQ"),
        TrimConfig(q_min=float(q_min), q_max=hi),
    )
    # core does not check the window caught anything
    if np.asarray(result.q).size < 2:
        q_arr = data["Q"]
        return _error(
            f"The Q window [{q_min}, {q_max}] contains fewer than 2 points.",
            suggestion=(
                f"The file spans Q = {float(q_arr.min()):.4g} to "
                f"{float(q_arr.max()):.4g} 1/Å — choose a window inside that."
            ),
            code="EMPTY_RESULT",
        )

    out_path, err = _save_manip(out_folder, data, result)
    if err:
        return err
    out = _base_result("trim", out_path, [data["_source_path"]])
    out["parameters"] = {"q_min": float(q_min), "q_max": hi if np.isfinite(hi) else None}
    out.update(_diagnostics(data["Q"].size, result.q, result.I, out_path))
    return out


def rebin_data(
    file: str,
    mode: str = "log",
    n_points: int = 200,
    q_min: Optional[float] = None,
    q_max: Optional[float] = None,
    reference_file: Optional[str] = None,
    output_folder: Optional[str] = None,
) -> dict:
    """Resample a dataset onto a new Q grid. Writes ``<stem>_rebinned.h5``.

    Interpolation is linear in log(I) vs log(Q) and does not extrapolate,
    so the result can be shorter than ``n_points``.

    Parameters
    ----------
    file : str
        Input dataset.
    mode : str
        'log' (geometric spacing), 'linear', or 'reference' (reuse the Q
        grid of ``reference_file``).
    n_points : int
        Target number of points for 'log' and 'linear'.
    q_min, q_max : float, optional
        Grid limits in 1/Å. Default to the data's own range.
    reference_file : str, optional
        Required for ``mode='reference'``; its Q grid is used verbatim.
    output_folder : str, optional
        Defaults to ``<source folder>_manip``.
    """
    from pyirena.core.data_manipulation import DataManipulation, RebinConfig

    if mode not in REBIN_MODES:
        return _error(
            f"Unknown rebin mode '{mode}'.",
            suggestion=f"Use one of: {', '.join(REBIN_MODES)}.",
            code="BAD_MODE",
        )
    data, err = _load(file)
    if err:
        return err

    reference_q = None
    if mode == "reference":
        if not reference_file:
            return _error(
                "mode='reference' needs a reference_file to take the Q grid from.",
                suggestion="Pass reference_file, or use mode='log'.",
                code="MISSING_REFERENCE",
            )
        ref, err = _load(reference_file)
        if err:
            return err
        reference_q = ref["Q"]
    elif int(n_points) < 2:
        return _error(
            f"n_points must be at least 2, got {n_points}.",
            suggestion="200 is a good default for a log grid.",
            code="BAD_N_POINTS",
        )

    lo = None if q_min is None else float(q_min)
    if mode == "log":
        # np.geomspace cannot start at or below zero
        if lo is not None and lo <= 0:
            return _error(
                f"mode='log' needs q_min > 0, got {q_min}.",
                suggestion="Use a positive q_min, or mode='linear'.",
                code="BAD_Q_RANGE",
            )
        if lo is None and float(data["Q"].min()) <= 0:
            return _error(
                "mode='log' needs positive Q, but the data starts at "
                f"{float(data['Q'].min()):.4g} 1/Å.",
                suggestion="Pass a positive q_min, or use mode='linear'.",
                code="BAD_Q_RANGE",
            )

    out_folder, err = _resolve_output_folder(output_folder, data["_source_path"], "_manip")
    if err:
        return err

    try:
        result = DataManipulation.rebin(
            data["Q"], data["Intensity"], data["Error"], data.get("dQ"),
            RebinConfig(
                mode=mode,
                n_points=int(n_points),
                q_min=lo,
                q_max=None if q_max is None else float(q_max),
                reference_q=reference_q,
            ),
        )
    except ValueError as exc:
        return _error(str(exc), suggestion="Check the rebin parameters.", code="REBIN_FAILED")

    out_path, err = _save_manip(out_folder, data, result)
    if err:
        return err
    out = _base_result("rebin", out_path, [data["_source_path"]])
    out["parameters"] = {
        "mode": mode,
        "n_points": int(n_points),
        "q_min": lo,
        "q_max": None if q_max is None else float(q_max),
        "reference_file": reference_file,
    }
    out.update(_diagnostics(data["Q"].size, result.q, result.I, out_path))
    return out


# ---------------------------------------------------------------------------
# Two-dataset operations
# ---------------------------------------------------------------------------

def subtract_data(
    sample_file: str,
    buffer_file: str,
    buffer_scale: float = 1.0,
    auto_scale: bool = False,
    auto_q_min: Optional[float] = None,
    auto_q_max: Optional[float] = None,
    output_folder: Optional[str] = None,
) -> dict:
    """Subtract a buffer/background dataset from a sample. Writes ``<stem>_sub.h5``.

    ``I_out = I_sample - buffer_scale * I_buffer``, with the buffer
    interpolated onto the sample's Q grid. The output keeps the sample's
    Q values, restricted to where the buffer could be interpolated — a
    buffer with a narrower Q range therefore truncates the result.

    Both datasets must have matching slit-smearing status; mixing a
    slit-smeared curve with a pinhole one is rejected.

    Parameters
    ----------
    sample_file, buffer_file : str
        The sample and the buffer/solvent/background to subtract.
    buffer_scale : float
        Multiplier applied to the buffer before subtraction.
    auto_scale : bool
        Fit ``buffer_scale`` from the integral ratio over
        [auto_q_min, auto_q_max]. Both bounds are required when this is on.
    auto_q_min, auto_q_max : float, optional
        Q window in 1/Å for the auto-scale fit.
    output_folder : str, optional
        Defaults to ``<source folder>_manip``.
    """
    from pyirena.core.data_manipulation import DataManipulation, SubtractConfig

    # core silently ignores auto_scale unless BOTH bounds are given
    if auto_scale and (auto_q_min is None or auto_q_max is None):
        return _error(
            "auto_scale needs both auto_q_min and auto_q_max.",
            suggestion=(
                "Pass both bounds (1/Å) for the region the sample and buffer "
                "should match in, or set auto_scale=false and give an "
                "explicit buffer_scale."
            ),
            code="MISSING_AUTO_RANGE",
        )

    sample, err = _load(sample_file)
    if err:
        return err
    buf, err = _load(buffer_file)
    if err:
        return err
    out_folder, err = _resolve_output_folder(output_folder, sample["_source_path"], "_manip")
    if err:
        return err

    s_smeared, s_slit = _smear(sample)
    b_smeared, b_slit = _smear(buf)
    try:
        result = DataManipulation.subtract(
            sample["Q"], sample["Intensity"], sample["Error"], sample.get("dQ"),
            buf["Q"], buf["Intensity"], buf["Error"],
            SubtractConfig(
                buffer_scale=float(buffer_scale),
                auto_scale=bool(auto_scale),
                auto_q_min=None if auto_q_min is None else float(auto_q_min),
                auto_q_max=None if auto_q_max is None else float(auto_q_max),
            ),
            is_slit_smeared_sample=s_smeared, slit_length_sample=s_slit,
            is_slit_smeared_buffer=b_smeared, slit_length_buffer=b_slit,
        )
    except ValueError as exc:
        return _error(
            str(exc),
            suggestion=(
                "Subtracting datasets with different slit smearing mixes two "
                "resolution functions. Use two pinhole or two matching "
                "slit-smeared datasets."
            ),
            code="SLIT_MISMATCH",
        )

    out_path, err = _save_manip(out_folder, sample, result)
    if err:
        return err
    out = _base_result("subtract", out_path, [sample["_source_path"], buf["_source_path"]])
    out["parameters"] = {
        "buffer_scale": float(result.metadata.get("buffer_scale", buffer_scale)),
        "auto_scale": bool(auto_scale),
        "auto_q_min": auto_q_min,
        "auto_q_max": auto_q_max,
        "is_slit_smeared": s_smeared,
        "slit_length": s_slit,
    }
    out.update(_diagnostics(sample["Q"].size, result.q, result.I, out_path))
    return out


def divide_data(
    numerator_file: str,
    denominator_file: str,
    denominator_scale: float = 1.0,
    denominator_background: float = 0.0,
    output_folder: Optional[str] = None,
) -> dict:
    """Divide one dataset by another. Writes ``<stem>_div.h5``.

    ``I_out = I_num / (denominator_scale * I_den - denominator_background)``,
    with the denominator interpolated onto the numerator's Q grid. Points
    where the denominator is exactly zero come back non-finite and are
    reported, then stripped on save.

    Slit-smearing status must match between the two datasets.

    Parameters
    ----------
    numerator_file, denominator_file : str
        The two datasets.
    denominator_scale : float
        Multiplier applied to the denominator before dividing.
    denominator_background : float
        Flat background subtracted from the scaled denominator.
    output_folder : str, optional
        Defaults to ``<source folder>_manip``.
    """
    from pyirena.core.data_manipulation import DataManipulation, DivideConfig

    num, err = _load(numerator_file)
    if err:
        return err
    den, err = _load(denominator_file)
    if err:
        return err
    out_folder, err = _resolve_output_folder(output_folder, num["_source_path"], "_manip")
    if err:
        return err

    n_smeared, n_slit = _smear(num)
    d_smeared, d_slit = _smear(den)
    try:
        result = DataManipulation.divide(
            num["Q"], num["Intensity"], num["Error"], num.get("dQ"),
            den["Q"], den["Intensity"], den["Error"],
            DivideConfig(
                denominator_scale=float(denominator_scale),
                denominator_background=float(denominator_background),
            ),
            is_slit_smeared_num=n_smeared, slit_length_num=n_slit,
            is_slit_smeared_den=d_smeared, slit_length_den=d_slit,
        )
    except ValueError as exc:
        return _error(
            str(exc),
            suggestion=(
                "Dividing datasets with different slit smearing mixes two "
                "resolution functions. Use two pinhole or two matching "
                "slit-smeared datasets."
            ),
            code="SLIT_MISMATCH",
        )

    out_path, err = _save_manip(out_folder, num, result)
    if err:
        return err
    out = _base_result("divide", out_path, [num["_source_path"], den["_source_path"]])
    out["parameters"] = {
        "denominator_scale": float(denominator_scale),
        "denominator_background": float(denominator_background),
    }
    out.update(_diagnostics(num["Q"].size, result.q, result.I, out_path))
    return out


# ---------------------------------------------------------------------------
# Many-dataset operation
# ---------------------------------------------------------------------------

def average_data(
    files: list,
    output_folder: Optional[str] = None,
    similarity_check: bool = False,
    similarity_p_min: float = 0.01,
    similarity_method: str = "cormap",
    similarity_reference: str = "first",
    similarity_normalize_scale: bool = True,
) -> dict:
    """Average two or more datasets onto the first one's Q grid.

    Writes ``<stem of first file>_avg.h5``. Uncertainties are the larger of
    the propagated error and the point-to-point standard deviation.

    Turn on ``similarity_check`` to screen for radiation damage: a cormap
    test compares the frames and discards outliers before averaging. The
    discarded frames come back in ``rejected``.

    Parameters
    ----------
    files : list of str
        Two or more datasets to average.
    output_folder : str, optional
        Defaults to ``<source folder>_manip``.
    similarity_check : bool
        Screen frames and drop outliers before averaging.
    similarity_p_min : float
        P-value threshold; frames below it are discarded. Typical 0.001-0.05.
    similarity_method : str
        Currently 'cormap'.
    similarity_reference : str
        'first' (compare each frame with frame 0, which is always kept) or
        'majority' (compare with the median of all frames).
    similarity_normalize_scale : bool
        Rescale each frame to the reference before comparing, so flux drift
        is not mistaken for a shape change.
    """
    from pyirena.core.data_manipulation import DataManipulation

    if not isinstance(files, (list, tuple)) or len(files) < 2:
        return _error(
            f"Averaging needs at least 2 files, got {len(files) if files else 0}.",
            suggestion="Pass a list of two or more dataset paths.",
            code="TOO_FEW_FILES",
        )
    if similarity_reference not in SIMILARITY_REFERENCES:
        return _error(
            f"Unknown similarity_reference '{similarity_reference}'.",
            suggestion=f"Use one of: {', '.join(SIMILARITY_REFERENCES)}.",
            code="BAD_REFERENCE",
        )

    loaded: list[dict] = []
    for f in files:
        data, err = _load(f)
        if err:
            return err
        loaded.append(data)

    # core's average() performs NO slit check — averaging a smeared curve
    # with a pinhole one would silently produce a meaningless mean.
    ref_smeared, ref_slit = _smear(loaded[0])
    for data in loaded[1:]:
        smeared, slit = _smear(data)
        ok, message = DataManipulation.check_slit_compatible(
            ref_smeared, ref_slit, smeared, slit, op="average"
        )
        if not ok:
            return _error(
                f"{message} ('{loaded[0]['_source_path'].name}' vs "
                f"'{data['_source_path'].name}')",
                suggestion=(
                    "Average only datasets with the same slit smearing — "
                    "mixing resolutions makes the mean meaningless."
                ),
                code="SLIT_MISMATCH",
            )

    datasets = [(d["Q"], d["Intensity"], d["Error"], d.get("dQ")) for d in loaded]
    n_in_total = sum(int(d["Q"].size) for d in loaded)

    rejected: list[dict] = []
    if similarity_check:
        from pyirena.core.similarity import check_similarity

        sim = check_similarity(
            datasets,
            filenames=[d["_source_path"].name for d in loaded],
            method=similarity_method,
            reference=similarity_reference,
            p_min=float(similarity_p_min),
            normalize_scale=bool(similarity_normalize_scale),
        )
        accepted = [r.idx for r in sim if r.accepted]
        rejected = [
            {"filename": r.filename, "p_value": float(r.p_value)}
            for r in sim
            if not r.accepted
        ]
        if len(accepted) < 2:
            return _error(
                f"Only {len(accepted)} of {len(datasets)} frames passed the "
                f"similarity filter (p_min={similarity_p_min}).",
                suggestion=(
                    "Lower similarity_p_min, or inspect the frames — this "
                    "much variation usually means radiation damage or a "
                    "sample change during the series."
                ),
                code="TOO_FEW_FILES",
            )
        datasets = [datasets[i] for i in accepted]
        loaded = [loaded[i] for i in accepted]

    result = DataManipulation.average(datasets, reference_index=0)

    out_folder, err = _resolve_output_folder(output_folder, loaded[0]["_source_path"], "_manip")
    if err:
        return err
    # batch.average_data omits slit_length here, which loses dQl on a
    # smeared average; pass it explicitly.
    out_path, err = _save_manip(out_folder, loaded[0], result, slit_length=ref_slit)
    if err:
        return err

    out = _base_result("average", out_path, [d["_source_path"] for d in loaded])
    out["n_datasets"] = len(datasets)
    out["rejected"] = rejected
    out["parameters"] = {
        "similarity_check": bool(similarity_check),
        "similarity_p_min": float(similarity_p_min),
        "similarity_method": similarity_method,
        "similarity_reference": similarity_reference,
        "similarity_normalize_scale": bool(similarity_normalize_scale),
        "is_slit_smeared": ref_smeared,
        "slit_length": ref_slit,
    }
    out.update(_diagnostics(n_in_total, result.q, result.I, out_path))
    return out


# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------

def merge_datasets(
    file1: str,
    file2: str,
    q_overlap_min: Optional[float] = None,
    q_overlap_max: Optional[float] = None,
    fit_scale: bool = True,
    scale_dataset: int = 2,
    fixed_scale_value: float = 1.0,
    fit_qshift: bool = False,
    fixed_qshift_value: float = 0.0,
    qshift_dataset: int = 0,
    split_at_left_cursor: bool = False,
    output_folder: Optional[str] = None,
) -> dict:
    """Merge two datasets that overlap in Q. Writes ``<stem of file1>_merged.h5``.

    ``file1`` is the lower-Q dataset and the absolute-intensity reference
    (typically USAXS); ``file2`` is the higher-Q dataset brought onto it
    (typically SAXS). A scale, and optionally a Q shift, are fitted in the
    overlap region.

    A background is **always** fitted and subtracted from file1 — there is
    no way to force it to zero. The fitted value is returned as
    ``background``; check it is small compared with your intensities.

    Mixing slit-smeared USAXS with pinhole SAXS is normal and allowed; any
    concern is returned in ``slit_warning`` rather than refused.

    Parameters
    ----------
    file1, file2 : str
        Lower-Q reference and higher-Q dataset.
    q_overlap_min, q_overlap_max : float, optional
        Overlap window in 1/Å. Omit both to auto-detect (intersect the two
        Q ranges and trim 10% off each side).
    fit_scale : bool
        Fit the scale factor. When False, ``fixed_scale_value`` is used.
    scale_dataset : int
        1 or 2 — which dataset is rescaled. 2 (the default) keeps file1 on
        its absolute scale and is much faster: it has a closed-form
        solution, whereas 1 falls back to iterative optimisation.
    fixed_scale_value : float
        Scale used when ``fit_scale`` is False.
    fit_qshift : bool
        Also fit a Q offset. Off by default.
    fixed_qshift_value : float
        Q shift in 1/Å used when ``fit_qshift`` is False.
    qshift_dataset : int
        0 (none), 1, or 2 — which dataset the Q shift applies to.
    split_at_left_cursor : bool
        True: hard split at q_overlap_min, no duplicate Q values. False
        (default): keep both datasets' points across the overlap.
    output_folder : str, optional
        Defaults to ``<parent of source folder>/<source folder>_merged``.
    """
    from pyirena.core.data_merge import DataMerge, MergeConfig

    data1, err = _load(file1)
    if err:
        return err
    data2, err = _load(file2)
    if err:
        return err

    q1, q2 = data1["Q"], data2["Q"]
    name1, name2 = data1["_source_path"].name, data2["_source_path"].name
    if float(q1.max()) <= float(q2.min()):
        return _error(
            f"The two datasets do not overlap in Q: '{name1}' ends at "
            f"{float(q1.max()):.4g} and '{name2}' starts at "
            f"{float(q2.min()):.4g} 1/Å.",
            suggestion=(
                "Pass file1 as the lower-Q dataset (usually USAXS) and file2 "
                "as the higher-Q one (usually SAXS)."
            ),
            code="NO_OVERLAP",
        )
    # file1 must be the lower-Q dataset: it is the absolute-intensity
    # reference and the background is subtracted from it. Swapped inputs
    # still "work" numerically but rescale the wrong curve.
    if float(q1.min()) > float(q2.min()):
        return _error(
            f"file1 and file2 look swapped: '{name1}' starts at "
            f"{float(q1.min()):.4g} 1/Å but '{name2}' starts lower, at "
            f"{float(q2.min()):.4g} 1/Å.",
            suggestion=(
                f"Pass the lower-Q dataset first: file1='{name2}', "
                f"file2='{name1}'. file1 is the absolute-intensity reference "
                "and keeps its scale."
            ),
            code="SWAPPED_INPUTS",
        )

    if q_overlap_min is None or q_overlap_max is None:
        # same auto-detect the batch layer uses
        lo = max(float(q1.min()), float(q2.min()))
        hi = min(float(q1.max()), float(q2.max()))
        span = hi - lo
        q_overlap_min = lo + 0.1 * span if q_overlap_min is None else float(q_overlap_min)
        q_overlap_max = hi - 0.1 * span if q_overlap_max is None else float(q_overlap_max)
    q_overlap_min, q_overlap_max = float(q_overlap_min), float(q_overlap_max)
    if q_overlap_max <= q_overlap_min:
        return _error(
            f"Empty overlap window: [{q_overlap_min}, {q_overlap_max}].",
            suggestion="Pass q_overlap_max greater than q_overlap_min, in 1/Å.",
            code="BAD_Q_RANGE",
        )

    try:
        # compared with == 2 in core, so the string "DS2" would silently
        # select the slower, less reliable iterative branch
        scale_dataset_i = int(scale_dataset)
        qshift_dataset_i = int(qshift_dataset)
    except (TypeError, ValueError):
        return _error(
            f"scale_dataset and qshift_dataset must be integers, got "
            f"{scale_dataset!r} and {qshift_dataset!r}.",
            suggestion="Use scale_dataset=2 (default) or 1; qshift_dataset 0, 1 or 2.",
            code="BAD_ARGUMENTS",
        )
    if scale_dataset_i not in (1, 2):
        return _error(
            f"scale_dataset must be 1 or 2, got {scale_dataset_i}.",
            suggestion="2 (rescale the higher-Q dataset) is the usual choice.",
            code="BAD_ARGUMENTS",
        )

    sl1 = float(data1.get("slit_length", 0.0) or 0.0)
    sl2 = float(data2.get("slit_length", 0.0) or 0.0)
    config = MergeConfig(
        q_overlap_min=q_overlap_min,
        q_overlap_max=q_overlap_max,
        fit_scale=bool(fit_scale),
        scale_dataset=scale_dataset_i,
        fixed_scale_value=float(fixed_scale_value),
        fit_qshift=bool(fit_qshift),
        fixed_qshift_value=float(fixed_qshift_value),
        qshift_dataset=qshift_dataset_i,
        split_at_left_cursor=bool(split_at_left_cursor),
        slit_length_ds1=sl1,
        slit_length_ds2=sl2,
    )

    engine = DataMerge()
    try:
        result = engine.optimize(
            q1, data1["Intensity"], data1["Error"],
            q2, data2["Intensity"], data2["Error"],
            config,
        )
    except Exception as exc:
        return _error(
            f"Merge optimisation failed: {exc}",
            suggestion="Check the overlap window contains data from both files.",
            code="MERGE_FAILED",
        )
    if not result.success:
        return _error(
            f"Merge optimisation did not converge: {result.message}",
            suggestion=(
                f"The overlap window [{q_overlap_min:.4g}, {q_overlap_max:.4g}] "
                "1/Å may hold too few points. Widen it, or omit both bounds to "
                "auto-detect."
            ),
            code="MERGE_FAILED",
        )

    try:
        # NOTE: merge() mutates `result`, stamping the slit fields — they
        # are only valid after this call.
        q_m, i_m, di_m, dq_m = engine.merge(
            q1, data1["Intensity"], data1["Error"], data1.get("dQ"),
            q2, data2["Intensity"], data2["Error"], data2.get("dQ"),
            result, config,
        )
    except Exception as exc:
        return _error(
            f"Merge assembly failed: {exc}",
            suggestion="Check both datasets have valid positive intensities.",
            code="MERGE_FAILED",
        )

    source1 = data1["_source_path"]
    if output_folder:
        try:
            out_folder = resolve_safe(output_folder, must_exist=False)
        except PathSecurityError as exc:
            return _path_error(exc)
    else:
        parent = source1.parent.parent if source1.parent.parent != source1.parent \
            else source1.parent
        try:
            out_folder = resolve_safe(
                str(parent / f"{source1.parent.name}_merged"), must_exist=False
            )
        except PathSecurityError as exc:
            return _path_error(exc)

    from pyirena.io.nxcansas_data_merge import save_merged_data

    merge_result_dict = {
        "scale": result.scale,
        "q_shift": result.q_shift,
        "background": result.background,
        "chi_squared": result.chi_squared,
        "n_overlap_points": result.n_overlap_points,
        "q_overlap_min": config.q_overlap_min,
        "q_overlap_max": config.q_overlap_max,
        "scale_dataset": config.scale_dataset,
        "fit_scale": config.fit_scale,
        "qshift_dataset": config.qshift_dataset,
        "fit_qshift": config.fit_qshift,
        "split_at_left_cursor": config.split_at_left_cursor,
        "slit_length_ds1": config.slit_length_ds1,
        "slit_length_ds2": config.slit_length_ds2,
        "slit_length_merged": result.slit_length_merged,
        "is_slit_smeared_merged": result.is_slit_smeared_merged,
    }
    try:
        out_path = save_merged_data(
            output_folder=out_folder,
            ds1_path=source1,
            ds1_is_nxcansas=bool(data1.get("is_nxcansas", True)),
            q=q_m, I=i_m, dI=di_m, dQ=dq_m,
            merge_result_dict=merge_result_dict,
            ds2_path=data2["_source_path"],
        )
    except ValueError as exc:
        return _error(
            f"Nothing left to save: {exc}",
            suggestion="The merged curve has fewer than 2 usable points.",
            code="EMPTY_RESULT",
        )
    except Exception as exc:
        return _error(
            f"Could not write the merged file: {exc}",
            suggestion="Check the output folder is writable.",
            code="SAVE_FAILED",
        )

    out = _base_result("merge", out_path, [source1, data2["_source_path"]])
    out["scale"] = float(result.scale)
    out["q_shift"] = float(result.q_shift)
    out["background"] = float(result.background)
    out["chi_squared"] = (
        float(result.chi_squared) if np.isfinite(result.chi_squared) else None
    )
    out["n_overlap_points"] = int(result.n_overlap_points)
    out["slit_length_merged"] = float(result.slit_length_merged)
    out["is_slit_smeared_merged"] = bool(result.is_slit_smeared_merged)
    if result.slit_warning:
        out["slit_warning"] = result.slit_warning
    out["parameters"] = {
        "q_overlap_min": q_overlap_min,
        "q_overlap_max": q_overlap_max,
        "fit_scale": bool(fit_scale),
        "scale_dataset": scale_dataset_i,
        "fit_qshift": bool(fit_qshift),
        "qshift_dataset": qshift_dataset_i,
        "split_at_left_cursor": bool(split_at_left_cursor),
    }
    out["note"] = (
        "A background is always fitted and subtracted from file1; it cannot "
        "be forced to zero. Check 'background' is small relative to your "
        "intensities."
    )
    out.update(_diagnostics(int(q1.size + q2.size), q_m, i_m, out_path))
    return out


def match_merge_files(folder1: str, folder2: str) -> dict:
    """Pair up files from two folders for batch merging. Reads only.

    Matches on a key of (text before the first underscore, last integer in
    the filename stem), so ``sampleA_usaxs_007.h5`` pairs with
    ``sampleA_saxs_007.h5``. Files with no partner are reported separately
    rather than silently dropped.

    Parameters
    ----------
    folder1, folder2 : str
        Folders holding the lower-Q and higher-Q datasets.

    Returns
    -------
    dict
        "pairs" (each with file1/file2 absolute paths ready for
        merge_datasets), plus "unmatched_1" and "unmatched_2".
    """
    from pyirena.api._paths import resolve_safe_folder
    from pyirena.core.data_merge import DataMerge
    from pyirena.core.file_sorting import sort_names

    try:
        dir1 = resolve_safe_folder(folder1)
        dir2 = resolve_safe_folder(folder2)
    except PathSecurityError as exc:
        return _path_error(exc)
    except (FileNotFoundError, NotADirectoryError) as exc:
        return _error(
            str(exc), suggestion="Check both folder paths.", code="FILE_NOT_FOUND"
        )

    exts = {".h5", ".hdf5", ".hdf", ".nx", ".nxs", ".dat", ".txt"}

    def _listing(d: Path) -> list[str]:
        names = [p.name for p in d.iterdir() if p.is_file() and p.suffix.lower() in exts]
        return sort_names(names)

    names1, names2 = _listing(dir1), _listing(dir2)
    pairs = DataMerge().match_files(names1, names2)

    matched1 = {a for a, _ in pairs}
    matched2 = {b for _, b in pairs}
    return {
        "folder1": str(dir1),
        "folder2": str(dir2),
        "pairs": [
            {"file1": str(dir1 / a), "file2": str(dir2 / b)} for a, b in pairs
        ],
        "n_pairs": len(pairs),
        "unmatched_1": [n for n in names1 if n not in matched1],
        "unmatched_2": [n for n in names2 if n not in matched2],
    }
