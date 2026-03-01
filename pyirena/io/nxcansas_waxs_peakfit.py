"""
HDF5 save / load for WAXS Peak-Fit results (NXcanSAS format extension).

Group path inside the NXcanSAS file:
    entry/waxs_peakfit_results   (NXprocess)

Structure
---------
waxs_peakfit_results/
    (attrs)  n_peaks, bg_shape, chi_squared, reduced_chi_squared, dof,
             q_min, q_max, timestamp, program, NX_class
    Q               — 1-D float64, Å⁻¹
    I_fit           — 1-D float64, total model
    I_bg            — 1-D float64, background-only model
    residuals       — 1-D float64, (I_data − I_fit) / I_error
    intensity_data  — 1-D float64, raw measured intensity
    intensity_error — 1-D float64, measurement uncertainty (optional)
    background/
        bg0 … bgN   — float64 scalar per coefficient
        (attrs on each) limit_low, limit_high
    background_std/
        bg0 … bgN   — float64 scalar, uncertainty
    peak_01/ … peak_NN/
        (attrs)     shape
        Q_peak      — 1-D float64, Q range ±5×FWHM around Q0
        I_peak      — 1-D float64, individual peak curve on Q_peak
        params/
            A, Q0, FWHM[, eta]  — float64 scalar each
            (attrs on each) limit_low, limit_high
        params_std/
            A, Q0, FWHM[, eta]  — float64 scalar, uncertainty
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

try:
    import h5py
except ImportError:
    h5py = None  # type: ignore

from pyirena.core.waxs_peakfit import (
    _PEAK_PARAM_NAMES,
    bg_param_names,
    eval_peak,
)

_GROUP = "entry/waxs_peakfit_results"
_PROGRAM = "pyIrena waxs_peakfit"


# ===========================================================================
# Save
# ===========================================================================

def save_waxs_peakfit_results(
    filepath: Path,
    result: Dict,
    q: np.ndarray,
    intensity_data: Optional[np.ndarray] = None,
    intensity_error: Optional[np.ndarray] = None,
    q_min: Optional[float] = None,
    q_max: Optional[float] = None,
) -> None:
    """Save WAXS peak-fit results into an NXcanSAS HDF5 file.

    Parameters
    ----------
    filepath : Path
        HDF5 file to write into (must already exist as a valid NXcanSAS file,
        or will be created if missing).
    result : dict
        Output dict from ``WAXSPeakFitModel.fit()``; must contain keys:
        ``success``, ``bg_shape``, ``bg_params``, ``bg_params_std``,
        ``peaks``, ``peaks_std``, ``chi2``, ``reduced_chi2``, ``dof``,
        ``I_model``, ``I_bg``, ``residuals``.
    q : array
        Full Q array in Å⁻¹ (same length as ``result['I_model']``).
    intensity_data : array or None
        Raw measured intensity.
    intensity_error : array or None
        Measurement uncertainty.
    q_min, q_max : float or None
        Q range used for fitting (written as metadata).
    """
    if h5py is None:
        raise ImportError("h5py is required to save results.")

    filepath = Path(filepath)
    peaks: List[Dict]  = result.get("peaks", [])
    bg_shape: str      = result.get("bg_shape", "Constant")
    bg_params: Dict    = result.get("bg_params", {})
    bg_params_std: Dict = result.get("bg_params_std", {})
    peaks_std: List    = result.get("peaks_std", [{} for _ in peaks])

    with h5py.File(filepath, "a") as f:
        # Remove old group if present
        if _GROUP in f:
            del f[_GROUP]

        grp = f.require_group(_GROUP)

        # ── Group attributes ──────────────────────────────────────────────
        grp.attrs["NX_class"]          = "NXprocess"
        grp.attrs["program"]           = _PROGRAM
        grp.attrs["timestamp"]         = datetime.now().isoformat()
        grp.attrs["n_peaks"]           = len(peaks)
        grp.attrs["bg_shape"]          = bg_shape
        grp.attrs["dof"]               = int(result.get("dof", 0))
        grp.attrs["q_min"]             = float(q_min) if q_min is not None else np.nan
        grp.attrs["q_max"]             = float(q_max) if q_max is not None else np.nan
        # Fit quality stored as scalar datasets (browseable/collectable in HDF5 viewer)
        grp.create_dataset("chi_squared",         data=float(result.get("chi2", np.nan)))
        grp.create_dataset("reduced_chi_squared",  data=float(result.get("reduced_chi2", np.nan)))

        # ── Main arrays ───────────────────────────────────────────────────
        q_arr = np.asarray(q, dtype=float)
        grp.create_dataset("Q", data=q_arr, dtype="float64")
        grp["Q"].attrs["units"] = "1/angstrom"

        I_model = result.get("I_model")
        if I_model is not None:
            grp.create_dataset("I_fit", data=np.asarray(I_model, float), dtype="float64")
            grp["I_fit"].attrs["units"] = "arb"

        I_bg = result.get("I_bg")
        if I_bg is not None:
            grp.create_dataset("I_bg", data=np.asarray(I_bg, float), dtype="float64")
            grp["I_bg"].attrs["units"] = "arb"

        resid = result.get("residuals")
        if resid is not None:
            grp.create_dataset("residuals", data=np.asarray(resid, float), dtype="float64")
            grp["residuals"].attrs["units"] = "dimensionless"

        if intensity_data is not None:
            grp.create_dataset("intensity_data",
                               data=np.asarray(intensity_data, float), dtype="float64")
            grp["intensity_data"].attrs["units"] = "arb"

        if intensity_error is not None:
            grp.create_dataset("intensity_error",
                               data=np.asarray(intensity_error, float), dtype="float64")
            grp["intensity_error"].attrs["units"] = "arb"

        # ── Background parameters ─────────────────────────────────────────
        bg_grp     = grp.create_group("background")
        bg_std_grp = grp.create_group("background_std")
        for name in bg_param_names(bg_shape):
            pd  = bg_params.get(name, {})
            val = float(pd.get("value", 0.0))
            ds  = bg_grp.create_dataset(name, data=val, dtype="float64")
            lo  = pd.get("lo")
            hi  = pd.get("hi")
            ds.attrs["limit_low"]  = float(lo) if lo is not None else np.nan
            ds.attrs["limit_high"] = float(hi) if hi is not None else np.nan
            std = float(bg_params_std.get(name, np.nan))
            bg_std_grp.create_dataset(name, data=std, dtype="float64")

        # ── Per-peak groups ───────────────────────────────────────────────
        for i, (peak, p_std) in enumerate(zip(peaks, peaks_std)):
            shape    = peak.get("shape", "Gauss")
            pnames   = _PEAK_PARAM_NAMES.get(shape, [])
            pk_grp   = grp.create_group(f"peak_{i+1:02d}")
            pk_grp.attrs["shape"] = shape

            # Individual peak curve (±5×FWHM around Q0)
            q0_val   = float(peak.get("Q0", {}).get("value", 0.0))
            fwhm_val = float(peak.get("FWHM", {}).get("value", 0.01))
            half_w   = 5.0 * max(fwhm_val, 1e-9)
            q_lo     = max(float(q_arr.min()), q0_val - half_w)
            q_hi     = min(float(q_arr.max()), q0_val + half_w)
            n_pts    = max(200, int((q_hi - q_lo) / max(fwhm_val / 20.0, 1e-10)))
            q_peak   = np.linspace(q_lo, q_hi, min(n_pts, 2000))

            p_vals   = {pn: float(peak[pn]["value"]) for pn in pnames if pn in peak}
            I_peak   = eval_peak(q_peak, shape, p_vals)

            pk_grp.create_dataset("Q_peak", data=q_peak, dtype="float64")
            pk_grp["Q_peak"].attrs["units"] = "1/angstrom"
            pk_grp.create_dataset("I_peak", data=I_peak, dtype="float64")
            pk_grp["I_peak"].attrs["units"] = "arb"

            # Fitted parameters
            p_grp   = pk_grp.create_group("params")
            ps_grp  = pk_grp.create_group("params_std")
            for pn in pnames:
                if pn not in peak:
                    continue
                pd  = peak[pn]
                val = float(pd.get("value", 0.0))
                ds  = p_grp.create_dataset(pn, data=val, dtype="float64")
                lo  = pd.get("lo")
                hi  = pd.get("hi")
                ds.attrs["limit_low"]  = float(lo) if lo is not None else np.nan
                ds.attrs["limit_high"] = float(hi) if hi is not None else np.nan
                std = float(p_std.get(pn, np.nan))
                ps_grp.create_dataset(pn, data=std, dtype="float64")


# ===========================================================================
# Load
# ===========================================================================

def load_waxs_peakfit_results(filepath: Path) -> Dict:
    """Load WAXS peak-fit results from an NXcanSAS HDF5 file.

    Returns
    -------
    dict with keys:
        n_peaks, bg_shape, chi_squared, reduced_chi_squared, dof,
        q_min, q_max, timestamp,
        Q, I_fit, I_bg, residuals, intensity_data, intensity_error,
        bg_params, bg_params_std,
        peaks (list of dicts), peaks_std (list of dicts),
        peak_curves (list of {Q_peak, I_peak} dicts)

    Raises
    ------
    KeyError
        If the ``entry/waxs_peakfit_results`` group is not found.
    """
    if h5py is None:
        raise ImportError("h5py is required to load results.")

    filepath = Path(filepath)

    with h5py.File(filepath, "r") as f:
        if _GROUP not in f:
            raise KeyError(f"No waxs_peakfit_results group in {filepath}")
        grp = f[_GROUP]

        attrs = dict(grp.attrs)
        n_peaks  = int(attrs.get("n_peaks", 0))
        bg_shape = str(attrs.get("bg_shape", "Constant"))

        def _arr(key):
            return np.array(grp[key], dtype=float) if key in grp else None

        def _scalar(key, default=np.nan):
            """Try dataset first (new format), fall back to attribute (old format)."""
            if key in grp and isinstance(grp[key], h5py.Dataset):
                try:
                    return float(grp[key][()])
                except Exception:
                    pass
            return float(attrs.get(key, default))

        result: Dict = {
            "n_peaks":           n_peaks,
            "bg_shape":          bg_shape,
            "chi_squared":       _scalar("chi_squared"),
            "reduced_chi_squared": _scalar("reduced_chi_squared"),
            "dof":               int(attrs.get("dof",   0)),
            "q_min":             float(attrs.get("q_min", np.nan)),
            "q_max":             float(attrs.get("q_max", np.nan)),
            "timestamp":         str(attrs.get("timestamp", "")),
            "Q":                 _arr("Q"),
            "I_fit":             _arr("I_fit"),
            "I_bg":              _arr("I_bg"),
            "residuals":         _arr("residuals"),
            "intensity_data":    _arr("intensity_data"),
            "intensity_error":   _arr("intensity_error"),
        }

        # Background
        bg_params     = {}
        bg_params_std = {}
        if "background" in grp:
            for name, ds in grp["background"].items():
                val = float(ds[()])
                lo  = ds.attrs.get("limit_low",  np.nan)
                hi  = ds.attrs.get("limit_high", np.nan)
                bg_params[name] = {
                    "value": val,
                    "fit":   True,
                    "lo":    None if np.isnan(lo) else float(lo),
                    "hi":    None if np.isnan(hi) else float(hi),
                }
        if "background_std" in grp:
            for name, ds in grp["background_std"].items():
                bg_params_std[name] = float(ds[()])
        result["bg_params"]     = bg_params
        result["bg_params_std"] = bg_params_std

        # Peaks
        peaks_list  = []
        peaks_std   = []
        peak_curves = []
        for i in range(1, n_peaks + 1):
            key = f"peak_{i:02d}"
            if key not in grp:
                continue
            pk_grp = grp[key]
            shape  = str(pk_grp.attrs.get("shape", "Gauss"))
            pnames = _PEAK_PARAM_NAMES.get(shape, [])

            peak_d  = {"shape": shape}
            p_std_d = {}
            if "params" in pk_grp:
                for pn in pnames:
                    if pn not in pk_grp["params"]:
                        continue
                    ds  = pk_grp["params"][pn]
                    val = float(ds[()])
                    lo  = ds.attrs.get("limit_low",  np.nan)
                    hi  = ds.attrs.get("limit_high", np.nan)
                    peak_d[pn] = {
                        "value": val,
                        "fit":   True,
                        "lo":    None if np.isnan(lo) else float(lo),
                        "hi":    None if np.isnan(hi) else float(hi),
                    }
            if "params_std" in pk_grp:
                for pn, ds in pk_grp["params_std"].items():
                    p_std_d[pn] = float(ds[()])

            q_pk  = np.array(pk_grp["Q_peak"], float) if "Q_peak" in pk_grp else None
            I_pk  = np.array(pk_grp["I_peak"], float) if "I_peak" in pk_grp else None

            peaks_list.append(peak_d)
            peaks_std.append(p_std_d)
            peak_curves.append({"Q_peak": q_pk, "I_peak": I_pk})

        result["peaks"]       = peaks_list
        result["peaks_std"]   = peaks_std
        result["peak_curves"] = peak_curves

    return result
