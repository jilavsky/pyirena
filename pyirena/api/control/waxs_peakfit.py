"""pyirena.api.control.waxs_peakfit — agent-drivable WAXS peak fitting.

The WAXS counterpart of the Unified Fit, Sizes, Simple Fits and Modeling
control surfaces.  A WAXS pattern is a smooth background with peaks on top, so
this surface is background choice plus peak-list management: find the peaks (or
place them by hand), set their shapes and starting values, fit, read back
position, width, amplitude and integrated area.

These functions operate on the *same*
:class:`~pyirena.api.control.session.Session` objects as the other control
surfaces, so the shared session-lifecycle and Q-range tools are reused as-is;
``set_fit_q_range`` limits the fitted range.

The model state lives in ``session.model`` (a
:class:`~pyirena.core.waxs_peakfit.WAXSPeakFitModel` carrying ``bg_params`` and
``peaks``) with ``session.model_name == "waxs_peakfit"``.

Peak parameters
---------------
Every peak is ``A`` (amplitude), ``Q0`` (position) and ``FWHM``, plus ``eta``
for a pseudo-Voigt.  Each carries a value, a fit flag and bounds, addressed by
peak index and parameter name.  The **area** under each peak is derived rather
than fitted, and is what most WAXS questions are really about — it is reported
with the results.

Typical workflow
----------------
>>> import pyirena.api.control as ctrl
>>> sid = ctrl.open_dataset("/data/waxs_scan.h5")["session_id"]
>>> ctrl.select_waxs_model(sid, bg_shape="SNIP")
>>> ctrl.find_waxs_peaks(sid)                  # data-driven starting point
>>> ctrl.list_waxs_peaks(sid)
>>> ctrl.set_waxs_peak_shape(sid, 0, "Pseudo-Voigt")
>>> ctrl.run_waxs_fit(sid)
>>> ctrl.get_waxs_results(sid)                 # positions, widths, areas
>>> ctrl.get_waxs_fit_image(sid)
>>> ctrl.save_waxs_fit(sid)

All functions return plain dicts.  Errors are
``{"error": ..., "code": ..., "suggestion": ...}`` dicts rather than exceptions.
"""
from __future__ import annotations

import base64
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np

from pyirena.api._paths import PathSecurityError, resolve_safe
from pyirena.api.control.errors import make_error, no_fit, no_session
from pyirena.api.control.session import fit_mask, get_session

__all__ = [
    "list_waxs_options",
    "select_waxs_model",
    "get_waxs_config",
    "set_waxs_background",
    "set_waxs_background_parameter",
    "find_waxs_peaks",
    "add_waxs_peak",
    "remove_waxs_peak",
    "list_waxs_peaks",
    "get_waxs_peak_parameters",
    "set_waxs_peak_shape",
    "set_waxs_peak_parameter",
    "set_waxs_peak_parameter_fit",
    "set_waxs_peak_parameter_bounds",
    "run_waxs_fit",
    "get_waxs_results",
    "get_waxs_fit_image",
    "save_waxs_fit",
]

WEIGHT_MODES = ["standard", "equal", "relative"]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _require_waxs(session_id: str):
    """Return (session, None) when a WAXS model is ready, else (None, error)."""
    s = get_session(session_id)
    if s is None:
        return None, no_session(session_id)
    if s.model is None or s.model_name != "waxs_peakfit":
        return None, make_error(
            f"Session '{session_id}' has no WAXS peak-fit model.",
            suggestion="Call select_waxs_model() first.",
            code="NO_WAXS_MODEL",
        )
    return s, None


def _require_peak(session_id: str, index: int):
    """Return (session, peak_dict, None) or (None, None, error)."""
    s, err = _require_waxs(session_id)
    if err:
        return None, None, err
    peaks = s.model.peaks
    if not isinstance(index, int) or index < 0 or index >= len(peaks):
        return None, None, make_error(
            f"No peak at index {index} (there are {len(peaks)}).",
            suggestion="Call list_waxs_peaks() to see the current indices.",
            code="BAD_PEAK",
        )
    return s, peaks[index], None


def _peak_param_names(peak) -> list:
    from pyirena.core.waxs_peakfit import _PEAK_PARAM_NAMES  # noqa: PLC0415

    return list(_PEAK_PARAM_NAMES.get(peak.get("shape", "Gauss"), []))


def _peak_row(index: int, peak, std: Optional[dict] = None) -> dict:
    from pyirena.core.waxs_peakfit import peak_area, peak_area_std  # noqa: PLC0415

    std = std or {}
    params = []
    for name in _peak_param_names(peak):
        entry = peak.get(name)
        if not isinstance(entry, dict):
            continue
        params.append({
            "name": name,
            "value": _f(entry.get("value")),
            "fit": bool(entry.get("fit", False)),
            "lo": _f(entry.get("lo")),
            "hi": _f(entry.get("hi")),
            "std": _f(std.get(name)),
        })

    # Area is derived, not fitted — and it is what most WAXS questions want.
    try:
        area = peak_area(peak.get("shape", "Gauss"), peak)
        area_err = peak_area_std(peak.get("shape", "Gauss"), peak, std)
    except (ValueError, KeyError, TypeError):
        area = area_err = None

    return {
        "index": index,
        "shape": peak.get("shape", "Gauss"),
        "parameters": params,
        "area": _f(area),
        "area_std": _f(area_err) or None,
    }


def _bg_rows(model) -> list:
    """Background parameters; adaptive shapes carry a value only."""
    from pyirena.core.waxs_peakfit import BG_ADAPTIVE, bg_param_names  # noqa: PLC0415

    adaptive = model.bg_shape in BG_ADAPTIVE
    rows = []
    for name in bg_param_names(model.bg_shape):
        entry = model.bg_params.get(name, {})
        row = {"name": name, "value": _f(entry.get("value"))}
        if not adaptive:
            row.update({
                "fit": bool(entry.get("fit", False)),
                "lo": _f(entry.get("lo")),
                "hi": _f(entry.get("hi")),
            })
        rows.append(row)
    return rows


def _config_dict(model) -> dict:
    from pyirena.core.waxs_peakfit import BG_ADAPTIVE  # noqa: PLC0415

    return {
        "bg_shape": model.bg_shape,
        "background_is_adaptive": model.bg_shape in BG_ADAPTIVE,
        "background_parameters": _bg_rows(model),
        "n_peaks": len(model.peaks),
    }


def _f(value):
    """Float, or None for missing/non-finite — JSON has no NaN."""
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def _render_image(fig, session_id: str, tag: str, dpi: int) -> str:
    import matplotlib.pyplot as plt  # noqa: PLC0415

    tmp = Path(tempfile.gettempdir()) / "pyirena-ctrl"
    tmp.mkdir(parents=True, exist_ok=True)
    out = tmp / f"waxs_{tag}_{session_id}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(out.read_bytes()).decode("ascii")


# ---------------------------------------------------------------------------
# Model lifecycle
# ---------------------------------------------------------------------------

def list_waxs_options() -> dict:
    """List the peak shapes, background shapes and weighting modes available.

    No session needed.  Adaptive backgrounds (SNIP, Rolling Quantile Spline,
    Rolling Ball) are estimated from the data rather than fitted — they have a
    tuning value but no fit flag, and are the right default for a WAXS pattern
    with a broad amorphous halo.
    """
    from pyirena.core.waxs_peakfit import (  # noqa: PLC0415
        BG_ADAPTIVE,
        BG_SHAPES,
        PEAK_SHAPES,
        bg_param_names,
    )

    return {
        "ok": True,
        "peak_shapes": list(PEAK_SHAPES),
        "background_shapes": [
            {
                "name": name,
                "adaptive": name in BG_ADAPTIVE,
                "parameters": bg_param_names(name),
            }
            for name in BG_SHAPES
        ],
        "weight_modes": list(WEIGHT_MODES),
    }


def select_waxs_model(session_id: str, bg_shape: str = "SNIP") -> dict:
    """Create a WAXS peak-fit model for the session, with no peaks yet.

    Parameters
    ----------
    bg_shape : str
        Background shape from :func:`list_waxs_options`.  ``"SNIP"`` (the
        default) estimates a smooth background from the data itself and suits
        most patterns; the polynomial shapes are fitted alongside the peaks.

    Notes
    -----
    Calling this again replaces the model and discards any peaks and fit.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    from pyirena.core.waxs_peakfit import (  # noqa: PLC0415
        BG_SHAPES,
        WAXSPeakFitModel,
        default_bg_params,
    )

    if bg_shape not in BG_SHAPES:
        return make_error(
            f"Unknown background shape '{bg_shape}'.",
            suggestion=f"Choose one of {list(BG_SHAPES)}.",
            code="BAD_BG_SHAPE",
        )

    model = WAXSPeakFitModel(bg_shape=bg_shape, peaks=[])
    model.bg_params = default_bg_params(bg_shape)

    s.model = model
    s.model_name = "waxs_peakfit"
    s.last_fit_result = None

    return {"ok": True, "model": "waxs_peakfit", "config": _config_dict(model)}


def get_waxs_config(session_id: str) -> dict:
    """Return the background setup and the current peak list."""
    s, err = _require_waxs(session_id)
    if err:
        return err
    return {
        "ok": True,
        "config": _config_dict(s.model),
        "peaks": [_peak_row(i, p) for i, p in enumerate(s.model.peaks)],
    }


# ---------------------------------------------------------------------------
# Background
# ---------------------------------------------------------------------------

def set_waxs_background(session_id: str, bg_shape: str) -> dict:
    """Switch the background shape, resetting its parameters to defaults.

    The peaks are kept — only the background changes, which is the usual way
    to test whether a stubborn residual is a background artefact.
    """
    s, err = _require_waxs(session_id)
    if err:
        return err

    from pyirena.core.waxs_peakfit import BG_SHAPES, default_bg_params  # noqa: PLC0415

    if bg_shape not in BG_SHAPES:
        return make_error(
            f"Unknown background shape '{bg_shape}'.",
            suggestion=f"Choose one of {list(BG_SHAPES)}.",
            code="BAD_BG_SHAPE",
        )

    s.model.bg_shape = bg_shape
    s.model.bg_params = default_bg_params(bg_shape)
    s.last_fit_result = None
    return {"ok": True, "config": _config_dict(s.model)}


def set_waxs_background_parameter(
    session_id: str,
    name: str,
    value: Optional[float] = None,
    fit: Optional[bool] = None,
    lo: Optional[float] = None,
    hi: Optional[float] = None,
) -> dict:
    """Set a background parameter's value, fit flag or bounds.

    Adaptive backgrounds (SNIP and friends) have a tuning *value* only — they
    are estimated from the data, not optimised, so ``fit`` and bounds do not
    apply and are ignored with a note.
    """
    s, err = _require_waxs(session_id)
    if err:
        return err

    from pyirena.core.waxs_peakfit import BG_ADAPTIVE, bg_param_names  # noqa: PLC0415

    if name not in bg_param_names(s.model.bg_shape):
        return make_error(
            f"'{name}' is not a parameter of the '{s.model.bg_shape}' background.",
            suggestion=(
                f"This background has {bg_param_names(s.model.bg_shape)}; "
                "call get_waxs_config() to see their values."
            ),
            code="BAD_BG_PARAM",
        )

    entry = s.model.bg_params.setdefault(name, {})
    note = None
    if value is not None:
        entry["value"] = float(value)
    if s.model.bg_shape in BG_ADAPTIVE:
        if fit is not None or lo is not None or hi is not None:
            note = ("Adaptive backgrounds are estimated from the data, not "
                    "fitted — fit flag and bounds ignored.")
    else:
        if fit is not None:
            entry["fit"] = bool(fit)
        if lo is not None:
            entry["lo"] = float(lo)
        if hi is not None:
            entry["hi"] = float(hi)

    out = {"ok": True, "name": name,
           "background_parameters": _bg_rows(s.model)}
    if note:
        out["note"] = note
    return out


# ---------------------------------------------------------------------------
# Peaks
# ---------------------------------------------------------------------------

def find_waxs_peaks(
    session_id: str,
    prominence_frac: float = 0.05,
    min_fwhm: float = 0.001,
    max_fwhm: float = 0.5,
    min_distance: float = 0.005,
    shape: str = "Gauss",
    replace: bool = True,
) -> dict:
    """Detect peaks in the data and add them, with starting values, to the model.

    The data-driven way to start: rather than guessing positions, this runs
    pyIrena's peak finder over the current fit Q range and creates one peak per
    detection with its position, amplitude and width already close.

    Parameters
    ----------
    prominence_frac : float
        Minimum prominence as a fraction of the intensity range.  Raise it to
        find fewer, stronger peaks; lower it to pick up shoulders.
    min_fwhm, max_fwhm : float
        Accepted width window [Å⁻¹].
    min_distance : float
        Minimum separation between detected peaks [Å⁻¹].
    shape : str
        Peak shape for the created peaks.
    replace : bool
        Replace the existing peak list (default) or append to it.
    """
    s, err = _require_waxs(session_id)
    if err:
        return err

    from pyirena.core.waxs_peakfit import PEAK_SHAPES, find_peaks_in_data  # noqa: PLC0415

    if shape not in PEAK_SHAPES:
        return make_error(
            f"Unknown peak shape '{shape}'.",
            suggestion=f"Choose one of {list(PEAK_SHAPES)}.",
            code="BAD_PEAK_SHAPE",
        )

    mask = fit_mask(s)
    if not np.any(mask):
        return make_error(
            "No data points in the current fit Q range.",
            suggestion="Call reset_fit_q_range() or widen set_fit_q_range().",
            code="EMPTY_RANGE",
        )

    try:
        # The finder returns ready-made peak dicts (value/fit/lo/hi per
        # parameter), so they go straight into the model.
        found = find_peaks_in_data(
            s.q[mask], s.intensity[mask],
            prominence_frac=prominence_frac,
            min_fwhm=min_fwhm, max_fwhm=max_fwhm,
            min_distance=min_distance,
            default_shape=shape,
        )
    except Exception as exc:  # pragma: no cover - defensive
        return make_error(
            f"Peak search failed: {exc}",
            suggestion="Try a different prominence_frac or Q range.",
            code="PEAK_SEARCH_FAILED",
        )

    if replace:
        s.model.peaks = list(found)
    else:
        s.model.peaks.extend(found)
    s.last_fit_result = None

    return {
        "ok": True,
        "n_found": len(found),
        "peaks": [_peak_row(i, p) for i, p in enumerate(s.model.peaks)],
        "message": (
            f"Found {len(found)} peak(s). Raise prominence_frac for fewer, "
            "lower it to pick up shoulders."
        ),
    }


def add_waxs_peak(session_id: str, q0: float, shape: str = "Gauss",
                  amplitude: Optional[float] = None,
                  fwhm: float = 0.01) -> dict:
    """Add one peak at position *q0*.

    When *amplitude* is omitted it is taken from the measured intensity at
    ``q0``, which is a far better starting point than a generic default.
    """
    s, err = _require_waxs(session_id)
    if err:
        return err

    from pyirena.core.waxs_peakfit import PEAK_SHAPES, default_peak_params  # noqa: PLC0415

    if shape not in PEAK_SHAPES:
        return make_error(
            f"Unknown peak shape '{shape}'.",
            suggestion=f"Choose one of {list(PEAK_SHAPES)}.",
            code="BAD_PEAK_SHAPE",
        )

    q0 = float(q0)
    if q0 < float(np.min(s.q)) or q0 > float(np.max(s.q)):
        return make_error(
            f"Q0 = {q0} is outside the data range.",
            suggestion=(
                f"The data covers {float(np.min(s.q)):.4g} to "
                f"{float(np.max(s.q)):.4g} Å⁻¹."
            ),
            code="Q0_OUT_OF_RANGE",
        )

    if amplitude is None:
        amplitude = float(s.intensity[int(np.argmin(np.abs(s.q - q0)))])

    s.model.peaks.append(
        default_peak_params(shape, Q0=q0, A=float(amplitude), FWHM=float(fwhm))
    )
    s.last_fit_result = None
    index = len(s.model.peaks) - 1
    return {"ok": True, "index": index, "peak": _peak_row(index, s.model.peaks[index])}


def remove_waxs_peak(session_id: str, index: int) -> dict:
    """Remove the peak at *index*; later peaks shift down."""
    s, peak, err = _require_peak(session_id, index)
    if err:
        return err
    s.model.peaks.pop(index)
    s.last_fit_result = None
    return {"ok": True, "removed": index,
            "peaks": [_peak_row(i, p) for i, p in enumerate(s.model.peaks)]}


def list_waxs_peaks(session_id: str) -> dict:
    """List every peak with its shape, parameters and derived area."""
    s, err = _require_waxs(session_id)
    if err:
        return err
    return {"ok": True,
            "peaks": [_peak_row(i, p) for i, p in enumerate(s.model.peaks)]}


def get_waxs_peak_parameters(session_id: str, index: int) -> dict:
    """Return one peak's parameters with values, fit flags, bounds and area."""
    s, peak, err = _require_peak(session_id, index)
    if err:
        return err
    return {"ok": True, "peak": _peak_row(index, peak)}


def set_waxs_peak_shape(session_id: str, index: int, shape: str) -> dict:
    """Change a peak's shape, keeping A, Q0 and FWHM.

    Switching to ``Pseudo-Voigt`` adds the ``eta`` mixing parameter (0 = pure
    Gaussian, 1 = pure Lorentzian); switching away drops it.
    """
    s, peak, err = _require_peak(session_id, index)
    if err:
        return err

    from pyirena.core.waxs_peakfit import PEAK_SHAPES, default_peak_params  # noqa: PLC0415

    if shape not in PEAK_SHAPES:
        return make_error(
            f"Unknown peak shape '{shape}'.",
            suggestion=f"Choose one of {list(PEAK_SHAPES)}.",
            code="BAD_PEAK_SHAPE",
        )

    fresh = default_peak_params(
        shape,
        Q0=float(peak["Q0"]["value"]),
        A=float(peak["A"]["value"]),
        FWHM=float(peak["FWHM"]["value"]),
    )
    # Keep the fit flags and bounds the agent already chose for shared params.
    for name in ("A", "Q0", "FWHM"):
        if isinstance(peak.get(name), dict):
            fresh[name].update({k: v for k, v in peak[name].items()
                                if k in ("fit", "lo", "hi")})
    s.model.peaks[index] = fresh
    s.last_fit_result = None
    return {"ok": True, "peak": _peak_row(index, fresh)}


def _peak_param(peak, name: str):
    entry = peak.get(name)
    return entry if isinstance(entry, dict) else None


def _bad_peak_param(name: str, peak) -> dict:
    return make_error(
        f"'{name}' is not a parameter of a {peak.get('shape', 'Gauss')} peak.",
        suggestion=(
            f"This shape has {_peak_param_names(peak)}; 'eta' exists only on "
            "a Pseudo-Voigt."
        ),
        code="BAD_PARAM",
    )


def set_waxs_peak_parameter(session_id: str, index: int,
                            name: str, value: float) -> dict:
    """Set one peak parameter's value (``A``, ``Q0``, ``FWHM`` or ``eta``)."""
    s, peak, err = _require_peak(session_id, index)
    if err:
        return err
    entry = _peak_param(peak, name)
    if entry is None or name not in _peak_param_names(peak):
        return _bad_peak_param(name, peak)
    try:
        entry["value"] = float(value)
    except (TypeError, ValueError):
        return make_error(
            f"Value '{value}' for '{name}' is not a number.",
            suggestion="Pass a numeric value.",
            code="BAD_VALUE",
        )
    return {"ok": True, "peak": _peak_row(index, peak)}


def set_waxs_peak_parameter_fit(session_id: str, index: int,
                                name: str, fit: bool = True) -> dict:
    """Choose whether one peak parameter is fitted or held at its value.

    Holding ``Q0`` at a known reflection position while fitting width and
    amplitude is the usual way to fit a phase you have already identified.
    """
    s, peak, err = _require_peak(session_id, index)
    if err:
        return err
    entry = _peak_param(peak, name)
    if entry is None or name not in _peak_param_names(peak):
        return _bad_peak_param(name, peak)
    entry["fit"] = bool(fit)
    return {"ok": True, "index": index, "name": name, "fit": bool(fit)}


def set_waxs_peak_parameter_bounds(
    session_id: str, index: int, name: str,
    lo: Optional[float] = None, hi: Optional[float] = None,
) -> dict:
    """Set bounds on one peak parameter.  ``None`` leaves that side unbounded."""
    s, peak, err = _require_peak(session_id, index)
    if err:
        return err
    entry = _peak_param(peak, name)
    if entry is None or name not in _peak_param_names(peak):
        return _bad_peak_param(name, peak)

    lo_f = None if lo is None else float(lo)
    hi_f = None if hi is None else float(hi)
    if lo_f is not None and hi_f is not None and lo_f >= hi_f:
        return make_error(
            f"Lower bound {lo_f} is not below upper bound {hi_f}.",
            suggestion="Pass lo < hi.",
            code="BAD_BOUNDS",
        )
    entry["lo"] = lo_f
    entry["hi"] = hi_f

    value = float(entry["value"])
    if lo_f is not None and value < lo_f:
        entry["value"] = lo_f
    elif hi_f is not None and value > hi_f:
        entry["value"] = hi_f

    return {"ok": True, "peak": _peak_row(index, peak)}


# ---------------------------------------------------------------------------
# Fit execution
# ---------------------------------------------------------------------------

def run_waxs_fit(session_id: str, weight_mode: str = "standard") -> dict:
    """Fit background and peaks to the data inside the current fit Q range.

    Parameters
    ----------
    weight_mode : str
        ``"standard"`` uses 1/σ² weighting; ``"equal"`` weights every point
        the same, which stops a low-noise background from dominating;
        ``"relative"`` (σ/I) emphasises peaks over the background.  Try
        ``"equal"`` when a strong background is pulling the fit away from
        weak peaks.

    Returns
    -------
    dict with success, chi_squared, reduced_chi_squared, dof, n_data, the
    background parameters, and every peak with its fitted values, 1σ
    uncertainties and integrated area.
    """
    s, err = _require_waxs(session_id)
    if err:
        return err

    if not s.model.peaks:
        return make_error(
            "No peaks to fit.",
            suggestion="Call find_waxs_peaks() or add_waxs_peak() first.",
            code="NO_PEAKS",
        )
    if weight_mode not in WEIGHT_MODES:
        return make_error(
            f"Unknown weight mode '{weight_mode}'.",
            suggestion=f"Choose one of {WEIGHT_MODES}.",
            code="BAD_WEIGHT_MODE",
        )

    mask = fit_mask(s)
    if not np.any(mask):
        return make_error(
            "No data points in the current fit Q range.",
            suggestion="Call reset_fit_q_range() or widen set_fit_q_range().",
            code="EMPTY_RANGE",
        )

    q_fit = s.q[mask]
    I_fit = s.intensity[mask]
    err_fit = s.error[mask] if s.error is not None else None

    try:
        result = s.model.fit(q_fit, I_fit, err_fit,
                             s.model.bg_params, s.model.peaks,
                             weight_mode=weight_mode)
    except Exception as exc:  # pragma: no cover - defensive
        return make_error(
            f"Fit failed with exception: {exc}",
            suggestion="Check the peak positions, widths and the Q range.",
            code="FIT_EXCEPTION",
        )

    if not result.get("success", False):
        return make_error(
            result.get("message") or "Fit did not converge.",
            suggestion="Try different starting positions, or weight_mode='equal'.",
            code="FIT_FAILED",
        )

    # Fitted values become the new state, as in the GUI.
    s.model.bg_params = result.get("bg_params", s.model.bg_params)
    s.model.peaks = result.get("peaks", s.model.peaks)
    s.last_fit_result = result

    peaks_std = result.get("peaks_std") or []
    return {
        "success": True,
        "weight_mode": weight_mode,
        "chi_squared": _f(result.get("chi2")),
        "reduced_chi_squared": _f(result.get("reduced_chi2")),
        "dof": (int(result["dof"]) if result.get("dof") is not None else None),
        "n_data": int(len(q_fit)),
        "bg_shape": s.model.bg_shape,
        "background_parameters": _bg_rows(s.model),
        "peaks": [
            _peak_row(i, p, peaks_std[i] if i < len(peaks_std) else {})
            for i, p in enumerate(s.model.peaks)
        ],
    }


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

def get_waxs_results(session_id: str) -> dict:
    """Return the last fit: quality, background, and every peak with its area."""
    s, err = _require_waxs(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    result = s.last_fit_result
    peaks_std = result.get("peaks_std") or []
    return {
        "ok": True,
        "chi_squared": _f(result.get("chi2")),
        "reduced_chi_squared": _f(result.get("reduced_chi2")),
        "dof": (int(result["dof"]) if result.get("dof") is not None else None),
        "bg_shape": s.model.bg_shape,
        "background_parameters": _bg_rows(s.model),
        "peaks": [
            _peak_row(i, p, peaks_std[i] if i < len(peaks_std) else {})
            for i, p in enumerate(s.model.peaks)
        ],
    }


def get_waxs_fit_image(
    session_id: str, width: int = 1024, height: int = 800, dpi: int = 120,
) -> dict:
    """Render the fit as a PNG: data, total model, background and each peak.

    Seeing the individual peaks over the background is how you spot a peak
    that has drifted onto its neighbour or collapsed to zero amplitude.
    """
    s, err = _require_waxs(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt  # noqa: PLC0415

    from pyirena.core.waxs_peakfit import eval_peak  # noqa: PLC0415

    result = s.last_fit_result
    mask = fit_mask(s)
    q = s.q[mask]
    I_obs = s.intensity[mask]
    I_model = np.asarray(result["I_model"], dtype=float)
    I_bg = np.asarray(result["I_bg"], dtype=float)
    residuals = np.asarray(result.get("residuals"), dtype=float)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(width / dpi, height / dpi), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    ax1.plot(q, I_obs, "o", markersize=3, alpha=0.6, color="steelblue",
             label="Data")
    ax1.plot(q, I_model, "-", linewidth=2, color="red", label="Total model")
    ax1.plot(q, I_bg, "--", linewidth=1.2, color="gray", label="Background")
    for i, peak in enumerate(s.model.peaks):
        try:
            values = {n: float(peak[n]["value"]) for n in _peak_param_names(peak)}
            ax1.plot(q, I_bg + eval_peak(q, peak["shape"], values),
                     ":", linewidth=1.0, alpha=0.8,
                     label=f"Peak {i + 1} ({peak['shape']})")
        except (KeyError, ValueError, TypeError):
            continue

    ax1.set_ylabel("Intensity", fontsize=10)
    chi2r = _f(result.get("reduced_chi2"))
    title = f"WAXS peak fit — {len(s.model.peaks)} peak(s), {s.model.bg_shape} background"
    if chi2r is not None:
        title += f"   reduced χ² = {chi2r:.3g}"
    ax1.set_title(title, fontsize=11, fontweight="bold")
    ax1.legend(fontsize=8, loc="best")
    ax1.grid(True, alpha=0.25)

    ax2.axhline(0.0, color="black", linewidth=1, linestyle="--")
    ax2.plot(q, residuals, "o", markersize=3, color="steelblue")
    ax2.set_xlabel("Q  (Å⁻¹)", fontsize=10)
    ax2.set_ylabel("Residuals", fontsize=10)
    ax2.grid(True, alpha=0.25)

    fig.tight_layout()
    return {"ok": True, "image_base64": _render_image(fig, session_id, "fit", dpi)}


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def save_waxs_fit(session_id: str, output_path: Optional[str] = None) -> dict:
    """Save the fit to NXcanSAS HDF5 under ``entry/waxs_peakfit_results``.

    Parameters
    ----------
    output_path : str, optional
        Where to save.  Defaults to the original file (in-place).  Saving to a
        different path preserves the original.
    """
    s, err = _require_waxs(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    from pyirena.io.nxcansas_waxs_peakfit import (  # noqa: PLC0415
        save_waxs_peakfit_results,
    )

    try:
        src = resolve_safe(s.file_path, must_exist=False)
        target = resolve_safe(output_path, must_exist=False) if output_path else src
    except PathSecurityError as exc:
        return make_error(
            str(exc),
            suggestion="Save to a path inside PYIRENA_DATA_ROOT.",
            code="PATH_NOT_ALLOWED",
        )

    if target != src and not target.exists():
        from pyirena.io._nxcansas_common import copy_and_strip_results  # noqa: PLC0415
        try:
            copy_and_strip_results(src, target)
        except Exception as exc:
            return make_error(
                f"Could not create output file '{target}' from source: {exc}",
                suggestion="Check the source file exists and the target is writable.",
                code="SAVE_ERROR",
            )

    mask = fit_mask(s)
    q_fit = s.q[mask]

    result = dict(s.last_fit_result)
    result.setdefault("bg_shape", s.model.bg_shape)
    result.setdefault("peaks", s.model.peaks)

    try:
        save_waxs_peakfit_results(
            filepath=target,
            result=result,
            q=q_fit,
            intensity_data=s.intensity[mask],
            intensity_error=(s.error[mask] if s.error is not None else None),
            q_min=float(np.min(q_fit)),
            q_max=float(np.max(q_fit)),
        )
    except Exception as exc:
        return make_error(
            f"Could not save results: {exc}",
            suggestion="Check the file is writable and not open elsewhere.",
            code="SAVE_ERROR",
        )

    return {"ok": True, "saved_to": str(target),
            "group": "entry/waxs_peakfit_results"}
