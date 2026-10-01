"""``export_results`` — one JSON report for a finished fit, whichever tool ran it.

Results have always left pyIrena as NXcanSAS HDF5 (``save_*``), which is the
right container for a file on disk and no use at all to a caller on another
machine with no shared filesystem. Unified Fit grew a JSON
``export_fit_report`` for that case; this is the same idea generalised to all
six fitting tools, as one tool rather than six, so an agent learns a single
name and gets the same envelope back every time.

What is in the envelope:

* **identification** — tool, pyIrena version, timestamp, the data's label and
  source file (null for a session opened from arrays)
* **quality** — chi-squared, reduced chi-squared, degrees of freedom and the
  fit-quality scalars, normalised across tools that spell them differently
* **results** — the tool's own answer, exactly as its ``get_*_results`` tool
  returns it, so there is one definition of "the result" and not two
* **config** — the model's ``to_dict()``, which is what ``pyirena.batch`` and
  the GUI's *Load Setup* consume: a result exported here can be reopened,
  re-applied to the next measurement, or diffed against another fit
* **arrays** — optional, off by default: the model curve, residuals and any
  per-component curves on the fitted Q grid

Uncertainties are included wherever the fit computes them (five of the six
tools report a ``std`` per parameter as part of fitting). Unified Fit has no
uncertainty estimation at all yet — see ``get_parameter_uncertainties`` — so
its parameters come back without one rather than with a fabricated zero.

Design doc: ``planning/zmq-service/`` (Phase 1.4).
"""
from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Optional

import numpy as np

from pyirena.api._json import to_strict_json
from pyirena.api.control.errors import make_error, no_fit, no_model, no_session
from pyirena.api.control.session import Session, fit_mask, get_session

# Default cap on any array returned by include_arrays. 2000 points is a full
# measured curve at the size real requests carry, so the default decimates
# nothing in practice; it is a ceiling for the pathological case.
_DEFAULT_MAX_POINTS = 2000


# ---------------------------------------------------------------------------
# Array packing
# ---------------------------------------------------------------------------

def _decimation_stride(n: int, max_points: Optional[int]) -> int:
    if not max_points or max_points <= 0 or n <= max_points:
        return 1
    return int(np.ceil(n / max_points))


def _pack_arrays(q, max_points: Optional[int], **named) -> dict:
    """Return q plus every companion array of the same length, decimated together.

    Arrays that do not line up with q are dropped rather than returned
    misaligned — a residual vector one point shorter than its Q axis is worse
    than no residual vector.
    """
    q = np.asarray(q, dtype=float)
    stride = _decimation_stride(len(q), max_points)

    out: dict[str, Any] = {"q": q[::stride]}
    for name, arr in named.items():
        if arr is None:
            continue
        arr = np.asarray(arr, dtype=float)
        if arr.shape != q.shape:
            continue
        out[name] = arr[::stride]

    out["n_points"] = int(len(out["q"]))
    out["decimated"] = stride > 1
    if stride > 1:
        out["decimation_stride"] = stride
    return out


def _masked(s: Session):
    """Session data restricted to the fitted Q range."""
    mask = fit_mask(s)
    return (
        s.q[mask],
        s.intensity[mask],
        s.error[mask] if s.error is not None else None,
    )


def _quality(chi2=None, reduced=None, dof=None, n_points=None, n_params=None,
             metrics=None, success=None, message=None) -> dict:
    """One quality block whatever the tool called its fields."""
    return {
        "chi_squared": chi2,
        "reduced_chi_squared": reduced,
        "dof": dof,
        "n_points_fitted": n_points,
        "n_parameters": n_params,
        "success": success,
        "message": message,
        "metrics": metrics,
    }


# ---------------------------------------------------------------------------
# Per-tool exporters — each returns {"results", "quality", "arrays"}
# ---------------------------------------------------------------------------

def _export_unified(s: Session, include_arrays: bool, max_points) -> dict:
    from pyirena.api.control.unified_fit import (
        _compute_quality,
        _param_table,
        _quality_scalars,
    )

    r = s.last_fit_result
    q, I, err = _masked(s)

    results = {
        "model": "unified_fit",
        "nlevels": int(s.model.num_levels),
        "background": float(getattr(s.model, "background", 0.0) or 0.0),
        "parameters": _param_table(s.model),
        "uncertainties_available": False,
    }
    quality = _quality(
        chi2=r.get("chi_squared"),
        reduced=r.get("reduced_chi_squared"),
        n_points=int(len(q)),
        metrics=_quality_scalars(_compute_quality(s.model)),
        success=r.get("success"),
        message=r.get("message"),
    )

    arrays = None
    if include_arrays:
        arrays = _pack_arrays(
            q, max_points,
            intensity=I, error=err,
            intensity_model=r.get("fit_intensity"),
            residuals=r.get("residuals"),
        )
    return {"results": results, "quality": quality, "arrays": arrays}


def _export_sizes(s: Session, include_arrays: bool, max_points) -> dict:
    from pyirena.api.control.sizes import get_sizes_distribution, get_sizes_results

    r = s.last_fit_result
    results = dict(get_sizes_results(s.session_id))
    results.pop("ok", None)

    # For a size distribution the histogram IS the result, not an optional
    # array, so it travels whether or not include_arrays is set. It honours
    # the caller's max_points rather than that tool's own 500-point default.
    dist = get_sizes_distribution(
        s.session_id, max_points=max_points if max_points else 10**9
    )
    if "error" not in dist:
        results["distribution"] = {
            k: v for k, v in dist.items() if k not in ("ok", "session_id")
        }

    quality = _quality(
        chi2=r.get("chi_squared"),
        n_points=r.get("n_data"),
        success=r.get("success"),
        message=r.get("message"),
    )

    arrays = None
    if include_arrays:
        arrays = _pack_arrays(
            r.get("q"), max_points,
            intensity=r.get("I_data"),
            error=r.get("err"),
            intensity_model=r.get("model_intensity"),
            intensity_model_ideal=r.get("model_intensity_ideal"),
            residuals=r.get("residuals"),
        )
    return {"results": results, "quality": quality, "arrays": arrays}


def _export_simple(s: Session, include_arrays: bool, max_points) -> dict:
    from pyirena.api.control.simple_fits import get_simple_results

    r = s.last_fit_result
    results = dict(get_simple_results(s.session_id))
    results.pop("ok", None)

    quality = _quality(
        chi2=r.get("chi2"),
        reduced=r.get("reduced_chi2"),
        dof=r.get("dof"),
        n_points=len(np.asarray(r.get("q", []))),
        success=r.get("success"),
        message=r.get("warning"),
    )

    arrays = None
    if include_arrays:
        _, I, err = _masked(s)
        arrays = _pack_arrays(
            r.get("q"), max_points,
            intensity=I, error=err,
            intensity_model=r.get("I_model"),
            intensity_model_ideal=r.get("I_model_ideal"),
            residuals=r.get("residuals"),
        )
    return {"results": results, "quality": quality, "arrays": arrays}


def _export_modeling(s: Session, include_arrays: bool, max_points) -> dict:
    from pyirena.api.control.modeling import get_modeling_results

    r = s.last_fit_result
    results = dict(get_modeling_results(s.session_id))
    results.pop("ok", None)

    quality = _quality(
        chi2=r.chi_squared,
        reduced=r.reduced_chi_squared,
        dof=r.dof,
        n_points=len(np.asarray(r.model_q)),
        message="; ".join(r.fit_warnings) if getattr(r, "fit_warnings", None) else None,
    )

    arrays = None
    if include_arrays:
        q = np.asarray(r.model_q, dtype=float)
        mask = (s.q >= r.config.q_min) & (s.q <= r.config.q_max)
        arrays = _pack_arrays(
            q, max_points,
            intensity=s.intensity[mask],
            error=s.error[mask] if s.error is not None else None,
            intensity_model=r.model_I,
            intensity_model_ideal=getattr(r, "model_I_ideal", None),
        )
        # Per-population curves are the point of a Modeling fit — they show
        # which population carries which decade.
        stride = _decimation_stride(len(q), max_points)
        populations = []
        for idx, curve in zip(getattr(r, "pop_indices", []) or [],
                              getattr(r, "pop_model_I", []) or []):
            curve = np.asarray(curve, dtype=float)
            if curve.shape == q.shape:
                populations.append({"index": int(idx), "intensity": curve[::stride]})
        if populations:
            arrays["populations"] = populations
    return {"results": results, "quality": quality, "arrays": arrays}


def _export_waxs(s: Session, include_arrays: bool, max_points) -> dict:
    from pyirena.api.control.waxs_peakfit import get_waxs_results

    r = s.last_fit_result
    results = dict(get_waxs_results(s.session_id))
    results.pop("ok", None)

    q, I, err = _masked(s)
    quality = _quality(
        chi2=r.get("chi2"),
        reduced=r.get("reduced_chi2"),
        dof=r.get("dof"),
        n_points=int(len(q)),
        success=r.get("success"),
        message=r.get("message"),
    )

    arrays = None
    if include_arrays:
        arrays = _pack_arrays(
            q, max_points,
            intensity=I, error=err,
            intensity_model=r.get("I_model"),
            intensity_background=r.get("I_bg"),
            residuals=r.get("residuals"),
        )
    return {"results": results, "quality": quality, "arrays": arrays}


def _export_carbon(s: Session, include_arrays: bool, max_points) -> dict:
    from pyirena.api.control.carbon_fit import get_carbon_results

    r = s.last_fit_result
    results = dict(get_carbon_results(s.session_id))
    results.pop("ok", None)
    results.pop("session_id", None)

    quality = _quality(
        chi2=r.chi_squared,
        reduced=r.reduced_chi_squared,
        n_points=r.n_points,
        n_params=r.n_params,
        metrics=results.get("fit_quality"),
        success=r.success,
        message=r.message,
    )

    arrays = None
    if include_arrays:
        arrays = _pack_arrays(
            r.q, max_points,
            intensity=r.I_data,
            error=r.I_error,
            intensity_model=r.I_model,
            residuals=r.residuals,
            # The three components: which decade each one owns is the whole
            # diagnostic value of a carbon fit.
            intensity_porod=getattr(r, "I_porod", None),
            intensity_micropore=getattr(r, "I_mp", None),
            intensity_waxs=getattr(r, "I_waxs", None),
        )
    return {"results": results, "quality": quality, "arrays": arrays}


_EXPORTERS = {
    "unified_fit": _export_unified,
    "sizes": _export_sizes,
    "simple_fits": _export_simple,
    "modeling": _export_modeling,
    "waxs_peakfit": _export_waxs,
    "carbon_fit": _export_carbon,
}


# ---------------------------------------------------------------------------
# Public tool
# ---------------------------------------------------------------------------

def export_results(
    session_id: str,
    include_arrays: bool = False,
    max_points: Optional[int] = _DEFAULT_MAX_POINTS,
) -> dict:
    """Export a finished fit as one JSON document, whichever tool produced it.

    The machine-readable counterpart to ``save_*``: everything a caller needs
    to record, compare or re-apply a fit, with no file written and no path
    involved. Use it when the caller is on another machine, when the result is
    going into a database or notebook, or when a fit needs to be replayed
    later against new data.

    Parameters
    ----------
    session_id : str
        A session with a completed fit, from any of the six fitting tools.
    include_arrays : bool
        Also return the curves: data, model, residuals and any per-component
        or per-population curves on the fitted Q grid. Off by default because
        the scalars are what most callers want and the arrays dominate the
        size of the reply.
    max_points : int or None
        Cap on the length of each returned array; longer ones are decimated by
        a constant stride and flagged with ``decimated``. None means no cap.
        Ignored unless ``include_arrays`` is set.

    Returns
    -------
    dict
        ``{ok, tool, pyirena_version, exported_at, data, quality, results,
        config, arrays?}``. ``config`` is the model's ``to_dict()``, the same
        shape ``pyirena.batch`` and the GUI's setup loader read, so a fit can
        be re-applied to the next measurement. Errors come back as the usual
        ``{error, code, suggestion}`` dict.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.model is None:
        return no_model(session_id)
    if s.last_fit_result is None:
        return no_fit(session_id)

    exporter = _EXPORTERS.get(s.model_name)
    if exporter is None:
        return make_error(
            f"No JSON exporter for model '{s.model_name}'.",
            suggestion=f"Supported: {', '.join(sorted(_EXPORTERS))}.",
            code="NO_EXPORTER",
        )

    from pyirena import __version__

    try:
        parts = exporter(s, include_arrays, max_points)
    except Exception as exc:                          # pragma: no cover - defensive
        return make_error(
            f"Could not export results for '{s.model_name}': {exc}",
            suggestion="Re-run the fit, or export without include_arrays.",
            code="EXPORT_ERROR",
        )

    try:
        config = s.model.to_dict()
    except Exception:
        config = None

    valid = np.isfinite(s.q) & np.isfinite(s.intensity)
    report = {
        "ok": True,
        "tool": s.model_name,
        "pyirena_version": __version__,
        "exported_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "data": {
            "label": s.label,
            "file": s.file_path,
            "n_points": int(len(s.q)),
            "q_min": float(np.nanmin(s.q[valid])),
            "q_max": float(np.nanmax(s.q[valid])),
            "has_errors": s.error is not None,
            "is_slit_smeared": bool(s.is_slit_smeared),
            "slit_length": float(s.slit_length or 0.0),
        },
        "fit_q_range": {"q_min": s.fit_q_min, "q_max": s.fit_q_max},
        "quality": parts["quality"],
        "results": parts["results"],
        "config": config,
    }

    # Simple Fits' held-parameter choices are session state, not model state,
    # so 'config' (a to_dict()) has no room for them. Without them here, a
    # report fed straight back into analyze() refits what the scientist
    # pinned and says nothing about it. Spelled the way save_simple_fit()
    # already writes it into the HDF5 setup attribute, so there is one name.
    if s.model_name == "simple_fits":
        from pyirena.api.control.simple_fits import _fixed_set
        report["fixed_params"] = sorted(_fixed_set(s))
    if parts["arrays"] is not None:
        report["arrays"] = parts["arrays"]

    # The whole point of this tool is that the result crosses a wire, so it
    # leaves strictly JSON-clean: no numpy scalars, no NaN, no Infinity.
    return to_strict_json(report)
