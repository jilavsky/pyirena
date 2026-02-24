"""
WAXS Peak-Fitting core: profile functions, background models, peak finder, fit engine.

All peak shapes are parameterized with **amplitude A** (= peak height at Q0) and
**FWHM** (full-width at half-maximum in Å⁻¹) so users see intuitive physical
quantities rather than the raw σ/γ fitting parameters.

Supported peak shapes
---------------------
* Gauss
* Lorentz  (Cauchy / Lorentzian)
* Pseudo-Voigt  (η·Lorentz + (1−η)·Gauss, η ∈ [0, 1])
* LogNormal  (asymmetric, mode at Q0, width set by FWHM)

Supported background shapes
---------------------------
* Constant   (1 parameter:  bg0)
* Linear     (2 parameters: bg0, bg1)
* Cubic      (4 parameters: bg0…bg3)
* 5th Poly   (6 parameters: bg0…bg5)

Peak finding
------------
Uses ``scipy.signal.find_peaks`` with a Savitzky-Golay baseline subtraction so
that peaks standing above a slowly-varying background are correctly identified.
"""

from __future__ import annotations

import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import optimize, signal


# ===========================================================================
# Constants
# ===========================================================================

PEAK_SHAPES  = ("Gauss", "Lorentz", "Pseudo-Voigt", "LogNormal")
BG_SHAPES    = ("Constant", "Linear", "Cubic", "5th Polynomial")

# Number of polynomial background coefficients per shape
_BG_NCOEFFS: Dict[str, int] = {
    "Constant":      1,
    "Linear":        2,
    "Cubic":         4,
    "5th Polynomial": 6,
}

# Parameter names for each peak shape (A, Q0, then shape-specific widths)
_PEAK_PARAM_NAMES: Dict[str, List[str]] = {
    "Gauss":        ["A", "Q0", "FWHM"],
    "Lorentz":      ["A", "Q0", "FWHM"],
    "Pseudo-Voigt": ["A", "Q0", "FWHM", "eta"],
    "LogNormal":    ["A", "Q0", "FWHM"],
}

# Default initial values and sensible fitting bounds per parameter name
_PARAM_DEFAULTS: Dict[str, Dict] = {
    "A":    {"value": 1.0,  "lo": 0.0,   "hi": None},
    "Q0":   {"value": 0.1,  "lo": 0.0,   "hi": None},
    "FWHM": {"value": 0.01, "lo": 1e-6,  "hi": None},
    "eta":  {"value": 0.5,  "lo": 0.0,   "hi": 1.0},
    "bg0":  {"value": 0.0,  "lo": None,  "hi": None},
    "bg1":  {"value": 0.0,  "lo": None,  "hi": None},
    "bg2":  {"value": 0.0,  "lo": None,  "hi": None},
    "bg3":  {"value": 0.0,  "lo": None,  "hi": None},
    "bg4":  {"value": 0.0,  "lo": None,  "hi": None},
    "bg5":  {"value": 0.0,  "lo": None,  "hi": None},
}


# ===========================================================================
# Peak profile functions
# ===========================================================================

_LN2 = np.log(2.0)


def gauss_peak(q: np.ndarray, A: float, Q0: float, FWHM: float) -> np.ndarray:
    """Gaussian peak: height A at Q0, full-width-half-max = FWHM."""
    sigma = FWHM / (2.0 * np.sqrt(2.0 * _LN2))
    return A * np.exp(-0.5 * ((q - Q0) / sigma) ** 2)


def lorentz_peak(q: np.ndarray, A: float, Q0: float, FWHM: float) -> np.ndarray:
    """Lorentzian (Cauchy) peak: height A at Q0, FWHM = FWHM."""
    gamma = FWHM / 2.0
    return A * gamma ** 2 / ((q - Q0) ** 2 + gamma ** 2)


def pseudo_voigt_peak(q: np.ndarray, A: float, Q0: float, FWHM: float,
                      eta: float) -> np.ndarray:
    """Pseudo-Voigt: linear mix η·Lorentz + (1−η)·Gauss, η ∈ [0, 1]."""
    eta = float(np.clip(eta, 0.0, 1.0))
    return eta * lorentz_peak(q, A, Q0, FWHM) + (1.0 - eta) * gauss_peak(q, A, Q0, FWHM)


def lognormal_peak(q: np.ndarray, A: float, Q0: float, FWHM: float) -> np.ndarray:
    """Log-normal peak with mode at Q0 and approximate FWHM.

    The log-normal is defined for q > 0.  For Q0 ≤ 0 or non-positive q
    values the function returns zero to avoid numerical issues.
    """
    q = np.asarray(q, dtype=float)
    out = np.zeros_like(q)
    if Q0 <= 0 or FWHM <= 0:
        return out
    # Relate sigma_ln to FWHM numerically: FWHM of log-normal in q-space.
    # For a log-normal with log-space sigma σ_ln, FWHM ≈ mode*(exp(σ_ln*√(2ln2)) − exp(−σ_ln*√(2ln2)))
    # We solve this implicitly via a Newton step approximation.
    # Good approximation valid for σ_ln < 1: FWHM ≈ Q0 * 2*sinh(σ_ln * √(2ln2))
    # => σ_ln ≈ arcsinh(FWHM / (2*Q0)) / √(2ln2)
    try:
        sigma_ln = float(np.arcsinh(FWHM / (2.0 * Q0)) / np.sqrt(2.0 * _LN2))
    except Exception:
        return out
    if sigma_ln <= 0:
        return out
    # Log-normal mode = exp(mu - sigma_ln^2), so mu = ln(Q0) + sigma_ln^2
    mu = np.log(Q0) + sigma_ln ** 2
    pos = q > 0
    out[pos] = (1.0 / (q[pos] * sigma_ln * np.sqrt(2.0 * np.pi))
                * np.exp(-0.5 * ((np.log(q[pos]) - mu) / sigma_ln) ** 2))
    # Normalise so peak height = A (value at Q0)
    peak_max = float(np.max(out)) if out.max() > 0 else 1.0
    out *= A / peak_max
    return out


def eval_peak(q: np.ndarray, shape: str, params: Dict) -> np.ndarray:
    """Evaluate a single peak on q from a parameter dict."""
    A, Q0, FWHM = float(params["A"]), float(params["Q0"]), float(params["FWHM"])
    if shape == "Gauss":
        return gauss_peak(q, A, Q0, FWHM)
    elif shape == "Lorentz":
        return lorentz_peak(q, A, Q0, FWHM)
    elif shape == "Pseudo-Voigt":
        return pseudo_voigt_peak(q, A, Q0, FWHM, float(params.get("eta", 0.5)))
    elif shape == "LogNormal":
        return lognormal_peak(q, A, Q0, FWHM)
    raise ValueError(f"Unknown peak shape: {shape!r}")


# ===========================================================================
# Background functions
# ===========================================================================

def eval_background(q: np.ndarray, bg_shape: str, coeffs: List[float]) -> np.ndarray:
    """Evaluate polynomial background on q."""
    q = np.asarray(q, dtype=float)
    n = _BG_NCOEFFS[bg_shape]
    coeffs = list(coeffs)[:n]
    result = np.zeros_like(q)
    for i, c in enumerate(coeffs):
        result += float(c) * q ** i
    return result


def bg_param_names(bg_shape: str) -> List[str]:
    """Return the coefficient names for the given background shape."""
    return [f"bg{i}" for i in range(_BG_NCOEFFS[bg_shape])]


# ===========================================================================
# Default state helpers
# ===========================================================================

def default_bg_params(bg_shape: str) -> Dict:
    """Return a fresh background params dict for bg_shape."""
    names = bg_param_names(bg_shape)
    return {
        name: {
            "value": _PARAM_DEFAULTS[name]["value"],
            "fit":   True,
            "lo":    _PARAM_DEFAULTS[name]["lo"],
            "hi":    _PARAM_DEFAULTS[name]["hi"],
        }
        for name in names
    }


def default_peak_params(shape: str, Q0: float = 0.1, A: float = 1.0,
                        FWHM: float = 0.01) -> Dict:
    """Return a fresh peak parameter dict with sensible initial guesses."""
    params: Dict = {
        "shape": shape,
        "A":     {"value": max(A, 1e-30),  "fit": True, "lo": 0.0,  "hi": None},
        "Q0":    {"value": Q0,             "fit": True, "lo": 0.0,  "hi": None},
        "FWHM":  {"value": max(FWHM, 1e-6),"fit": True, "lo": 1e-6, "hi": None},
    }
    if shape == "Pseudo-Voigt":
        params["eta"] = {"value": 0.5, "fit": True, "lo": 0.0, "hi": 1.0}
    return params


# ===========================================================================
# Peak finder
# ===========================================================================

def find_peaks_in_data(
    q: np.ndarray,
    I: np.ndarray,
    prominence_frac: float = 0.05,
    min_fwhm: float = 0.001,
    max_fwhm: float = 0.5,
    min_distance: float = 0.005,
    sg_window_frac: float = 0.15,
    default_shape: str = "Gauss",
) -> List[Dict]:
    """Identify peaks in I(q) data above a smooth background.

    Parameters
    ----------
    q, I : array
        Scattering vector (Å⁻¹) and intensity.  Must be 1-D and already
        sorted in ascending q order.
    prominence_frac : float
        Minimum peak prominence expressed as a fraction of (max − min) of I.
    min_fwhm, max_fwhm : float
        FWHM bounds in Å⁻¹.  Peaks narrower than *min_fwhm* or wider than
        *max_fwhm* are discarded.
    min_distance : float
        Minimum separation between adjacent peaks in Å⁻¹.
    sg_window_frac : float
        Savitzky-Golay baseline window as a fraction of len(q).  Larger
        values → smoother baseline.
    default_shape : str
        Peak shape assigned to every found peak.

    Returns
    -------
    list of dict
        Each dict has keys ``shape``, ``A``, ``Q0``, ``FWHM``, plus ``eta``
        for Pseudo-Voigt.  Suitable as input to ``default_peak_params``.
    """
    q = np.asarray(q, dtype=float)
    I = np.asarray(I, dtype=float)
    n = len(q)
    if n < 10:
        return []

    # ── 1. Savitzky-Golay baseline ────────────────────────────────────────
    sg_win = max(5, int(sg_window_frac * n))
    if sg_win % 2 == 0:
        sg_win += 1           # must be odd
    sg_win = min(sg_win, n if n % 2 == 1 else n - 1)
    try:
        baseline = signal.savgol_filter(I, window_length=sg_win, polyorder=3)
    except Exception:
        baseline = np.zeros_like(I)
    I_sub = I - baseline

    # ── 2. Convert physical params → sample-space params ─────────────────
    q_span = float(q[-1] - q[0])
    q_per_sample = q_span / max(n - 1, 1)
    min_width_samples = max(1, min_fwhm / q_per_sample)
    max_width_samples = max(min_width_samples + 1, max_fwhm / q_per_sample)
    min_dist_samples  = max(1, int(min_distance / q_per_sample))
    prominence_abs    = prominence_frac * (float(I.max()) - float(I.min()))

    # ── 3. Find peaks ─────────────────────────────────────────────────────
    try:
        peak_indices, props = signal.find_peaks(
            I_sub,
            prominence=max(prominence_abs, 0.0),
            width=(min_width_samples, max_width_samples),
            distance=min_dist_samples,
        )
    except Exception:
        return []

    # ── 4. Convert results back to physical units ─────────────────────────
    results: List[Dict] = []
    widths_samples = props.get("widths", np.ones(len(peak_indices)))
    for idx, w_s in zip(peak_indices, widths_samples):
        q0   = float(q[idx])
        amp  = max(float(I[idx] - baseline[idx]), 1e-30)
        fwhm = float(w_s * q_per_sample)
        fwhm = float(np.clip(fwhm, min_fwhm, max_fwhm))
        results.append(default_peak_params(default_shape, Q0=q0, A=amp, FWHM=fwhm))

    return results


# ===========================================================================
# Full model evaluation
# ===========================================================================

def eval_model(
    q: np.ndarray,
    bg_shape: str,
    bg_params: Dict,
    peaks: List[Dict],
) -> np.ndarray:
    """Evaluate the full model (background + all peaks) on q.

    Parameters are provided as nested dicts (GUI representation) with
    ``{'value': float, 'fit': bool, 'lo': ..., 'hi': ...}`` sub-dicts.
    """
    q = np.asarray(q, dtype=float)
    # Background
    coeffs = [float(bg_params[n]["value"]) for n in bg_param_names(bg_shape)]
    result = eval_background(q, bg_shape, coeffs)
    # Peaks
    for peak in peaks:
        shape = peak["shape"]
        p = {k: float(peak[k]["value"]) for k in _PEAK_PARAM_NAMES[shape]
             if k in peak}
        result += eval_peak(q, shape, p)
    return result


# ===========================================================================
# Fit engine
# ===========================================================================

class WAXSPeakFitModel:
    """Fit engine for WAXS peak fitting.

    Parameters
    ----------
    bg_shape : str
        One of ``BG_SHAPES``.
    peaks : list of dict
        Each dict from ``default_peak_params()`` / GUI state.
    no_limits : bool
        When True, fitting is unconstrained (bounds ignored).
    """

    def __init__(
        self,
        bg_shape: str,
        peaks: List[Dict],
        no_limits: bool = False,
    ):
        self.bg_shape   = bg_shape
        self.peaks      = peaks          # list of param dicts (shape + per-param sub-dicts)
        self.no_limits  = no_limits

    # ── Parameter vector helpers ──────────────────────────────────────────

    def _build_param_list(self) -> List[Tuple[str, float, bool, Optional[float], Optional[float]]]:
        """Return list of (name_tag, value, fit, lo, hi) for every free/fixed param."""
        items = []
        # Background
        for name in bg_param_names(self.bg_shape):
            pd = self.bg_params[name]
            items.append((f"bg:{name}", float(pd["value"]), bool(pd["fit"]),
                          pd.get("lo"), pd.get("hi")))
        # Peaks
        for i, peak in enumerate(self.peaks):
            shape = peak["shape"]
            for pname in _PEAK_PARAM_NAMES[shape]:
                if pname not in peak:
                    continue
                pd = peak[pname]
                items.append((f"p{i}:{pname}", float(pd["value"]), bool(pd["fit"]),
                              pd.get("lo"), pd.get("hi")))
        return items

    def _pack(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[Tuple]]:
        """Pack GUI state into p0, lb, ub, free_mask for curve_fit."""
        all_params = self._build_param_list()
        p0, lb, ub = [], [], []
        free_tags  = []
        fixed_vals = {}

        for tag, val, fit, lo, hi in all_params:
            if fit:
                p0.append(val)
                lb.append(-np.inf if (self.no_limits or lo is None) else lo)
                ub.append( np.inf if (self.no_limits or hi is None) else hi)
                free_tags.append(tag)
            else:
                fixed_vals[tag] = val

        return np.array(p0), np.array(lb), np.array(ub), free_tags, fixed_vals

    def _make_model_func(self, free_tags, fixed_vals):
        """Return a callable f(q, *free_params) → I_model."""
        bg_names   = bg_param_names(self.bg_shape)
        n_bg       = len(bg_names)
        n_peaks    = len(self.peaks)
        peak_shapes = [p["shape"] for p in self.peaks]
        peak_pnames = [_PEAK_PARAM_NAMES[s] for s in peak_shapes]

        def f(q, *free_vals):
            # Reconstruct full param dict from free + fixed
            full: Dict[str, float] = dict(fixed_vals)
            for tag, v in zip(free_tags, free_vals):
                full[tag] = float(v)

            # Background
            coeffs = [full[f"bg:{n}"] for n in bg_names]
            result = eval_background(q, self.bg_shape, coeffs)

            # Peaks
            for i in range(n_peaks):
                shape = peak_shapes[i]
                p = {pn: full[f"p{i}:{pn}"] for pn in peak_pnames[i]}
                result += eval_peak(q, shape, p)
            return result

        return f

    # ── Public API ────────────────────────────────────────────────────────

    def set_state(self, bg_params: Dict, peaks: List[Dict]):
        """Set the parameter state from GUI dicts."""
        self.bg_params = bg_params
        self.peaks     = peaks

    def predict(self, q: np.ndarray, bg_params: Dict, peaks: List[Dict]) -> np.ndarray:
        """Evaluate model at q without fitting."""
        self.set_state(bg_params, peaks)
        p0, lb, ub, free_tags, fixed_vals = self._pack()
        f = self._make_model_func(free_tags, fixed_vals)
        full_vals: Dict[str, float] = dict(fixed_vals)
        for tag, v in zip(free_tags, p0):
            full_vals[tag] = v
        bg_names   = bg_param_names(self.bg_shape)
        coeffs     = [full_vals[f"bg:{n}"] for n in bg_names]
        result     = eval_background(np.asarray(q, float), self.bg_shape, coeffs)
        for i, peak in enumerate(peaks):
            shape = peak["shape"]
            p     = {pn: full_vals[f"p{i}:{pn}"] for pn in _PEAK_PARAM_NAMES[shape]}
            result += eval_peak(np.asarray(q, float), shape, p)
        return result

    def fit(
        self,
        q: np.ndarray,
        I: np.ndarray,
        dI: Optional[np.ndarray],
        bg_params: Dict,
        peaks: List[Dict],
    ) -> Dict:
        """Fit the model to data.

        Returns
        -------
        dict with keys:
            success, message, chi2, reduced_chi2, dof,
            bg_params (updated), peaks (updated),
            bg_params_std, peaks_std,
            I_model, I_bg, residuals
        """
        self.set_state(bg_params, peaks)
        q  = np.asarray(q, float)
        I  = np.asarray(I, float)
        dI = np.asarray(dI, float) if dI is not None else None

        # Filter to finite, positive-q data
        mask = np.isfinite(q) & np.isfinite(I) & (q > 0)
        if dI is not None:
            mask &= np.isfinite(dI) & (dI > 0)
        q_, I_ = q[mask], I[mask]
        sigma_ = dI[mask] if dI is not None else np.ones_like(I_)

        if len(q_) < 3:
            return {"success": False, "message": "Too few valid data points."}

        p0, lb, ub, free_tags, fixed_vals = self._pack()
        if len(p0) == 0:
            # All parameters fixed — just evaluate
            f = self._make_model_func(free_tags, fixed_vals)
            I_model = f(q_, )
            resid   = (I_ - I_model) / sigma_
            chi2    = float(np.sum(resid ** 2))
            return self._package_result(
                q, I, dI, mask, free_tags, fixed_vals, p0,
                np.zeros((0, 0)), chi2, len(q_), True,
                "All parameters fixed — no fitting performed.",
            )

        f = self._make_model_func(free_tags, fixed_vals)

        # Clamp p0 inside bounds
        p0 = np.clip(p0, lb, ub)

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                popt, pcov = optimize.curve_fit(
                    f, q_, I_, p0=p0,
                    sigma=sigma_, absolute_sigma=True,
                    bounds=(lb, ub),
                    maxfev=50_000,
                )
            success = True
            message = "Fit converged."
        except optimize.OptimizeWarning as exc:
            popt    = p0
            pcov    = np.full((len(p0), len(p0)), np.inf)
            success = False
            message = f"OptimizeWarning: {exc}"
        except RuntimeError as exc:
            popt    = p0
            pcov    = np.full((len(p0), len(p0)), np.inf)
            success = False
            message = f"Fit did not converge: {exc}"

        I_model = f(q_, *popt)
        resid   = (I_ - I_model) / sigma_
        chi2    = float(np.sum(resid ** 2))
        dof     = max(1, len(q_) - len(popt))

        return self._package_result(
            q, I, dI, mask, free_tags, fixed_vals, popt, pcov,
            chi2, dof, success, message,
        )

    def _package_result(
        self, q, I, dI, mask, free_tags, fixed_vals, popt, pcov,
        chi2, dof, success, message,
    ) -> Dict:
        """Unpack fit result into updated param dicts."""
        import copy

        # Std deviations from covariance diagonal
        if pcov.size > 0 and np.all(np.isfinite(pcov.diagonal())):
            pstd = np.sqrt(np.abs(np.diag(pcov)))
        else:
            pstd = np.full(len(popt), np.nan)

        full_vals = dict(fixed_vals)
        full_stds: Dict[str, float] = {tag: 0.0 for tag in fixed_vals}
        for tag, v, s in zip(free_tags, popt, pstd):
            full_vals[tag] = float(v)
            full_stds[tag] = float(s)

        # Rebuild bg_params dict with fitted values
        bg_params_out = copy.deepcopy(self.bg_params)
        for n in bg_param_names(self.bg_shape):
            tag = f"bg:{n}"
            bg_params_out[n]["value"] = full_vals.get(tag, float(self.bg_params[n]["value"]))
        bg_params_std = {n: full_stds.get(f"bg:{n}", 0.0) for n in bg_param_names(self.bg_shape)}

        # Rebuild peaks list with fitted values
        peaks_out = copy.deepcopy(self.peaks)
        peaks_std = []
        for i, peak in enumerate(peaks_out):
            shape = peak["shape"]
            p_std = {}
            for pn in _PEAK_PARAM_NAMES[shape]:
                if pn not in peak:
                    continue
                tag = f"p{i}:{pn}"
                peak[pn]["value"] = full_vals.get(tag, float(self.peaks[i][pn]["value"]))
                p_std[pn] = full_stds.get(tag, 0.0)
            peaks_std.append(p_std)

        # Full-array model and residuals
        f_full  = self._make_model_func(free_tags, full_vals)
        q_arr   = np.asarray(q, float)
        I_model = np.full_like(q_arr, np.nan)
        I_bg    = np.full_like(q_arr, np.nan)
        resids  = np.full_like(q_arr, np.nan)
        bg_names = bg_param_names(self.bg_shape)
        coeffs   = [full_vals.get(f"bg:{n}", float(self.bg_params[n]["value"])) for n in bg_names]
        valid    = np.isfinite(q_arr) & (q_arr > 0)
        I_model[valid] = f_full(q_arr[valid], )
        I_bg[valid]    = eval_background(q_arr[valid], self.bg_shape, coeffs)
        if dI is not None:
            dI_arr = np.asarray(dI, float)
            s_arr  = np.where((dI_arr > 0) & np.isfinite(dI_arr), dI_arr, 1.0)
        else:
            s_arr = np.ones_like(q_arr)
        I_arr = np.asarray(I, float)
        resids[valid & mask] = ((I_arr - I_model) / s_arr)[valid & mask]

        return {
            "success":       success,
            "message":       message,
            "chi2":          chi2,
            "reduced_chi2":  chi2 / dof,
            "dof":           dof,
            "bg_params":     bg_params_out,
            "bg_params_std": bg_params_std,
            "peaks":         peaks_out,
            "peaks_std":     peaks_std,
            "I_model":       I_model,
            "I_bg":          I_bg,
            "residuals":     resids,
        }
