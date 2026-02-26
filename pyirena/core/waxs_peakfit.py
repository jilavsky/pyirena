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
Polynomial (optimised by curve_fit simultaneously with peaks):
* Constant   (1 parameter:  bg0)
* Linear     (2 parameters: bg0, bg1)
* Cubic      (4 parameters: bg0…bg3)
* 5th Poly   (6 parameters: bg0…bg5)

Adaptive / data-driven (pre-computed from data, not optimised):
* SNIP                  — iterative peak-clipping (standard for XRD/powder)
* Rolling Quantile Spline — percentile filter + cubic spline (Watier method)
* Rolling Ball          — morphological erosion + dilation

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
BG_SHAPES    = (
    "Constant", "Linear", "Cubic", "5th Polynomial",
    "SNIP", "Rolling Quantile Spline", "Rolling Ball",
)

# Adaptive background shapes: pre-computed from data, NOT fitted by curve_fit
BG_ADAPTIVE: frozenset = frozenset({"SNIP", "Rolling Quantile Spline", "Rolling Ball"})

# Number of polynomial background coefficients per shape
_BG_NCOEFFS: Dict[str, int] = {
    "Constant":      1,
    "Linear":        2,
    "Cubic":         4,
    "5th Polynomial": 6,
}

# Parameters for adaptive background shapes (value, bounds, display label)
_BG_ADAPTIVE_PARAMS: Dict[str, Dict[str, Dict]] = {
    "SNIP": {
        "half_width": {
            "value": 0.10, "lo": 0.02, "hi": 0.60,
            "label": "Half-width (fraction of data)",
        },
    },
    "Rolling Quantile Spline": {
        "window": {
            "value": 0.12, "lo": 0.02, "hi": 0.60,
            "label": "Window (fraction of data)",
        },
        "quantile": {
            "value": 0.10, "lo": 0.00, "hi": 0.49,
            "label": "Quantile  (0 = min)",
        },
    },
    "Rolling Ball": {
        "radius": {
            "value": 0.10, "lo": 0.02, "hi": 0.60,
            "label": "Radius (fraction of data)",
        },
    },
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
    """Return the parameter names for the given background shape."""
    if bg_shape in BG_ADAPTIVE:
        return list(_BG_ADAPTIVE_PARAMS[bg_shape].keys())
    return [f"bg{i}" for i in range(_BG_NCOEFFS[bg_shape])]


# ===========================================================================
# Adaptive (data-driven) background functions
# ===========================================================================

def snip_background(I: np.ndarray, half_width_frac: float) -> np.ndarray:
    """SNIP background estimation via iterative peak-clipping.

    The Statistics-sensitive Non-linear Iterative Peak-clipping (SNIP)
    algorithm iteratively replaces each point with the average of its
    neighbours at increasing separation, clipping from above so that
    peaks are progressively removed while the underlying background is
    retained.  It is the standard background estimator for XRD/powder
    diffraction and is particularly robust for data with many narrow
    peaks sitting on a slowly-varying continuum.

    Parameters
    ----------
    I : array_like
        Intensity values (1-D, linear scale).
    half_width_frac : float
        Maximum clip half-width expressed as a fraction of the data
        length.  For *n* data points the algorithm iterates p from 1 to
        ``max(2, int(half_width_frac * n))``.
    """
    I = np.asarray(I, float)
    n = len(I)
    hw = max(2, int(half_width_frac * n))
    bg = I.copy()
    for p in range(1, hw + 1):
        left           = np.empty_like(bg)
        right          = np.empty_like(bg)
        left[:p]       = bg[:p]
        left[p:]       = bg[:-p]
        right[:-p]     = bg[p:]
        right[-p:]     = bg[-p:]
        avg            = 0.5 * (left + right)
        bg             = np.minimum(bg, avg)
    return bg


def rolling_quantile_background(
    q: np.ndarray,
    I: np.ndarray,
    window_frac: float,
    quantile: float,
) -> np.ndarray:
    """Estimate background as a rolling low-quantile spline (Watier method).

    1. Apply ``scipy.ndimage.percentile_filter`` with a window of
       ``window_frac × n`` points and the given *quantile* (0 = minimum).
    2. Fit a ``CubicSpline`` through the resulting anchor points.

    Parameters
    ----------
    q, I : array_like
        Scattering vector and intensity (1-D).
    window_frac : float
        Rolling window width as a fraction of the data length.
    quantile : float
        Quantile fraction in [0, 0.5).  0 = rolling minimum.
    """
    from scipy.ndimage import percentile_filter
    from scipy.interpolate import CubicSpline

    q = np.asarray(q, float)
    I = np.asarray(I, float)
    n = len(I)
    window_pts = max(3, int(window_frac * n))
    if window_pts % 2 == 0:
        window_pts += 1          # percentile_filter works best with odd size
    pct = float(np.clip(quantile * 100.0, 0.0, 100.0))
    anchors = percentile_filter(I, pct, size=window_pts)
    # Prevent the spline from exceeding the actual data
    anchors = np.clip(anchors, I.min(), I.max())
    return CubicSpline(q, anchors)(q)


def rolling_ball_background(I: np.ndarray, radius_frac: float) -> np.ndarray:
    """Estimate background via morphological rolling-ball (erosion + dilation).

    Equivalent to rolling a ball of given radius along the underside of the
    spectrum.  Implemented as a grey erosion followed by grey dilation with a
    flat structuring element of size ``2*r+1`` where *r* = ``radius_frac × n``.

    Parameters
    ----------
    I : array_like
        Intensity values (1-D, linear scale).
    radius_frac : float
        Ball radius as a fraction of the data length.
    """
    from scipy.ndimage import grey_erosion, grey_dilation

    I = np.asarray(I, float)
    r = max(1, int(radius_frac * len(I)))
    size = 2 * r + 1
    eroded = grey_erosion(I,  size=size)
    return   grey_dilation(eroded, size=size)


def compute_adaptive_background(
    q: np.ndarray,
    I: np.ndarray,
    bg_shape: str,
    bg_params: Dict,
) -> np.ndarray:
    """Dispatch to the requested adaptive background algorithm.

    Parameters
    ----------
    q, I : array_like
        Scattering vector and intensity within the fit range.
    bg_shape : str
        One of the names in ``BG_ADAPTIVE``.
    bg_params : dict
        Parameter dict as returned by ``default_bg_params()`` or the GUI,
        e.g. ``{"half_width": {"value": 0.10}}``.

    Returns
    -------
    np.ndarray
        Background curve on the same grid as *q*.
    """
    q = np.asarray(q, float)
    I = np.asarray(I, float)
    if bg_shape == "SNIP":
        hw = float(bg_params.get("half_width", {}).get("value", 0.10))
        return snip_background(I, hw)
    if bg_shape == "Rolling Quantile Spline":
        win = float(bg_params.get("window",   {}).get("value", 0.12))
        qu  = float(bg_params.get("quantile", {}).get("value", 0.10))
        return rolling_quantile_background(q, I, win, qu)
    if bg_shape == "Rolling Ball":
        rad = float(bg_params.get("radius", {}).get("value", 0.10))
        return rolling_ball_background(I, rad)
    raise ValueError(f"Not an adaptive background: {bg_shape!r}")


# ===========================================================================
# Default state helpers
# ===========================================================================

def default_bg_params(bg_shape: str) -> Dict:
    """Return a fresh background params dict for bg_shape.

    For polynomial shapes the dict follows the standard GUI format::

        {name: {"value": float, "fit": bool, "lo": float|None, "hi": float|None}}

    For adaptive shapes only ``"value"`` is stored (no Fit?/lo/hi)::

        {name: {"value": float}}
    """
    if bg_shape in BG_ADAPTIVE:
        return {
            name: {"value": info["value"]}
            for name, info in _BG_ADAPTIVE_PARAMS[bg_shape].items()
        }
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
    I: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Evaluate the full model (background + all peaks) on q.

    Parameters are provided as nested dicts (GUI representation) with
    ``{'value': float, 'fit': bool, 'lo': ..., 'hi': ...}`` sub-dicts.

    For adaptive background shapes (*bg_shape* in ``BG_ADAPTIVE``) the raw
    intensity array *I* must be supplied so that the background can be
    computed from the data.
    """
    q = np.asarray(q, dtype=float)
    # Background
    if bg_shape in BG_ADAPTIVE:
        if I is None:
            raise ValueError(
                f"eval_model: I (intensity) must be provided for adaptive "
                f"background '{bg_shape}'."
            )
        result = compute_adaptive_background(q, np.asarray(I, float), bg_shape, bg_params)
    else:
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
        """Return list of (name_tag, value, fit, lo, hi) for every free/fixed param.

        For adaptive background shapes the background parameters are NOT
        included — the background is pre-computed and captured in the model
        function closure, so curve_fit only optimises the peak parameters.
        """
        items = []
        # Background — polynomial shapes only (adaptive bg is fixed in closure)
        if self.bg_shape not in BG_ADAPTIVE:
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
                lo_v = -np.inf if (self.no_limits or lo is None) else float(lo)
                hi_v =  np.inf if (self.no_limits or hi is None) else float(hi)
                # Guard: swap if limits are accidentally inverted (e.g. negative
                # background coefficients where lo*frac > hi*frac after sign flip).
                if lo_v > hi_v:
                    lo_v, hi_v = hi_v, lo_v
                lb.append(lo_v)
                ub.append(hi_v)
                free_tags.append(tag)
            else:
                fixed_vals[tag] = val

        return np.array(p0), np.array(lb), np.array(ub), free_tags, fixed_vals

    def _make_model_func(self, free_tags, fixed_vals, adaptive_bg=None):
        """Return a callable f(q, *free_params) → I_model.

        Parameters
        ----------
        adaptive_bg : np.ndarray or None
            Pre-computed background array (same length as the q passed to
            curve_fit).  When provided (*bg_shape* in BG_ADAPTIVE) the
            background is constant inside curve_fit — only peak parameters
            are optimised.
        """
        bg_shape    = self.bg_shape
        bg_names    = bg_param_names(bg_shape) if bg_shape not in BG_ADAPTIVE else []
        n_peaks     = len(self.peaks)
        peak_shapes = [p["shape"] for p in self.peaks]
        peak_pnames = [_PEAK_PARAM_NAMES[s] for s in peak_shapes]

        def f(q, *free_vals):
            # Reconstruct full param dict from free + fixed
            full: Dict[str, float] = dict(fixed_vals)
            for tag, v in zip(free_tags, free_vals):
                full[tag] = float(v)

            # Background
            if adaptive_bg is not None:
                result = np.array(adaptive_bg, dtype=float)
            else:
                coeffs = [full[f"bg:{n}"] for n in bg_names]
                result = eval_background(q, bg_shape, coeffs)

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

    def predict(
        self,
        q: np.ndarray,
        bg_params: Dict,
        peaks: List[Dict],
        I: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Evaluate model at q without fitting.

        Parameters
        ----------
        I : array_like or None
            Raw intensity data at *q*, required when *bg_shape* is an
            adaptive shape (in ``BG_ADAPTIVE``).
        """
        self.set_state(bg_params, peaks)
        p0, lb, ub, free_tags, fixed_vals = self._pack()

        # Precompute adaptive background if needed
        adaptive_bg = None
        if self.bg_shape in BG_ADAPTIVE:
            if I is None:
                raise ValueError(
                    f"predict(): I (intensity) required for adaptive background '{self.bg_shape}'."
                )
            adaptive_bg = compute_adaptive_background(
                np.asarray(q, float), np.asarray(I, float), self.bg_shape, bg_params
            )

        f = self._make_model_func(free_tags, fixed_vals, adaptive_bg=adaptive_bg)
        return f(np.asarray(q, float), *p0)

    def fit(
        self,
        q: np.ndarray,
        I: np.ndarray,
        dI: Optional[np.ndarray],
        bg_params: Dict,
        peaks: List[Dict],
        weight_mode: str = "standard",
    ) -> Dict:
        """Fit the model to data.

        Parameters
        ----------
        weight_mode : str
            How to weight the least-squares residuals:

            * ``"standard"`` — standard ``1/σ²`` weighting (default).
              Uses the measured uncertainties *dI* when available,
              otherwise uniform weights.
            * ``"equal"``    — all points weighted equally (σ_eff = 1).
              Prevents background points with small σ from dominating.
            * ``"relative"`` — relative weighting ``σ_eff = dI / I``.
              Emphasises peaks relative to the smooth background.

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
        dI_    = dI[mask] if dI is not None else None

        print(f"[WAXSModel.fit] bg={self.bg_shape!r}, n_in={len(q)}, "
              f"n_valid={mask.sum()}, has_dI={dI is not None}, weight={weight_mode!r}")

        if len(q_) < 3:
            print("[WAXSModel.fit] Too few valid points — returning failure")
            return {"success": False, "message": "Too few valid data points."}

        # ── Weight selection ─────────────────────────────────────────────
        if weight_mode == "equal":
            sigma_ = np.ones_like(I_)
        elif weight_mode == "relative":
            if dI_ is not None and np.all(I_ > 0):
                sigma_ = dI_ / np.clip(I_, 1e-30, None)
            else:
                sigma_ = 1.0 / np.clip(I_, 1e-30, None)
        else:  # "standard"
            sigma_ = dI_ if dI_ is not None else np.ones_like(I_)

        # ── Adaptive background pre-computation ──────────────────────────
        adaptive_bg_fit = None   # subset for curve_fit
        adaptive_bg_full = None  # full-array version for _package_result
        if self.bg_shape in BG_ADAPTIVE:
            print(f"[WAXSModel.fit] computing adaptive background ({self.bg_shape}) …")
            adaptive_bg_full_arr = compute_adaptive_background(
                q, I, self.bg_shape, bg_params
            )
            adaptive_bg_fit  = adaptive_bg_full_arr[mask]
            adaptive_bg_full = adaptive_bg_full_arr
            print("[WAXSModel.fit] adaptive background done")

        p0, lb, ub, free_tags, fixed_vals = self._pack()
        print(f"[WAXSModel.fit] n_free={len(p0)}, free_tags={free_tags}")
        if len(p0) == 0:
            # All parameters fixed — just evaluate
            print("[WAXSModel.fit] all params fixed — evaluating model only")
            f = self._make_model_func(free_tags, fixed_vals, adaptive_bg=adaptive_bg_fit)
            I_model = f(q_)
            resid   = (I_ - I_model) / sigma_
            chi2    = float(np.sum(resid ** 2))
            return self._package_result(
                q, I, dI, mask, free_tags, fixed_vals, p0,
                np.zeros((0, 0)), chi2, len(q_), True,
                "All parameters fixed — no fitting performed.",
                adaptive_bg_full=adaptive_bg_full,
            )

        f = self._make_model_func(free_tags, fixed_vals, adaptive_bg=adaptive_bg_fit)

        # Clamp p0 inside bounds
        p0 = np.clip(p0, lb, ub)

        sigma_min = float(sigma_.min())
        sigma_max = float(sigma_.max())
        sigma_ratio = sigma_max / max(sigma_min, 1e-30)
        print(f"[WAXSModel.fit] calling curve_fit: n_free={len(p0)}, "
              f"sigma=[{sigma_min:.3g}, {sigma_max:.3g}] ratio={sigma_ratio:.0f}x")
        if sigma_ratio > 200:
            print(f"[WAXSModel.fit] WARNING: sigma range ratio {sigma_ratio:.0f}x — "
                  f"extreme weights may slow convergence. "
                  f"Consider 'Equal' or 'Relative' weighting for faster fits.")

        # Wrap f to count calls
        _ncalls = [0]
        def _f_counted(q_, *args):
            _ncalls[0] += 1
            return f(q_, *args)

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                popt, pcov = optimize.curve_fit(
                    _f_counted, q_, I_, p0=p0,
                    sigma=sigma_, absolute_sigma=True,
                    bounds=(lb, ub),
                    maxfev=10_000,
                    # Relaxed convergence: scipy default ~1.5e-8 (machine precision)
                    # causes thousands of Jacobian iterations for large datasets.
                    # 1e-5 is more than adequate for WAXS peak positions/widths.
                    ftol=1e-5, xtol=1e-5, gtol=1e-5,
                )
            success = True
            message = "Fit converged."
            print(f"[WAXSModel.fit] curve_fit converged in {_ncalls[0]} calls")
        except optimize.OptimizeWarning as exc:
            popt    = p0
            pcov    = np.full((len(p0), len(p0)), np.inf)
            success = False
            message = f"OptimizeWarning: {exc}"
            print(f"[WAXSModel.fit] OptimizeWarning after {_ncalls[0]} calls: {exc}")
        except RuntimeError as exc:
            popt    = p0
            pcov    = np.full((len(p0), len(p0)), np.inf)
            success = False
            message = f"Fit did not converge: {exc}"
            print(f"[WAXSModel.fit] did not converge after {_ncalls[0]} calls: {exc}")

        print(f"[WAXSModel.fit] packaging result …")
        I_model = f(q_, *popt)
        resid   = (I_ - I_model) / sigma_
        chi2    = float(np.sum(resid ** 2))
        dof     = max(1, len(q_) - len(popt))
        print(f"[WAXSModel.fit] chi2={chi2:.4g}, dof={dof}, success={success}")

        return self._package_result(
            q, I, dI, mask, free_tags, fixed_vals, popt, pcov,
            chi2, dof, success, message,
            adaptive_bg_full=adaptive_bg_full,
        )

    def _package_result(
        self, q, I, dI, mask, free_tags, fixed_vals, popt, pcov,
        chi2, dof, success, message, adaptive_bg_full=None,
    ) -> Dict:
        """Unpack fit result into updated param dicts.

        Parameters
        ----------
        adaptive_bg_full : np.ndarray or None
            Pre-computed adaptive background on the full *q* grid.  When
            provided it is used directly for I_bg instead of calling
            ``eval_background()``.
        """
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

        # Rebuild bg_params dict with fitted values (adaptive params unchanged)
        bg_params_out = copy.deepcopy(self.bg_params)
        bg_params_std: Dict[str, float] = {}
        if self.bg_shape not in BG_ADAPTIVE:
            for n in bg_param_names(self.bg_shape):
                tag = f"bg:{n}"
                bg_params_out[n]["value"] = full_vals.get(tag, float(self.bg_params[n]["value"]))
            bg_params_std = {n: full_stds.get(f"bg:{n}", 0.0)
                             for n in bg_param_names(self.bg_shape)}
        # For adaptive shapes: bg_params_std is empty (no fitted bg params)

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
        q_arr   = np.asarray(q, float)
        I_model = np.full_like(q_arr, np.nan)
        I_bg    = np.full_like(q_arr, np.nan)
        resids  = np.full_like(q_arr, np.nan)
        valid   = np.isfinite(q_arr) & (q_arr > 0)

        # Reconstruct adaptive or polynomial background on valid points
        if self.bg_shape in BG_ADAPTIVE:
            if adaptive_bg_full is not None:
                I_bg[valid] = np.asarray(adaptive_bg_full, float)[valid]
            adaptive_bg_valid = I_bg[valid]  # may be NaN if not provided
        else:
            bg_names = bg_param_names(self.bg_shape)
            coeffs   = [full_vals.get(f"bg:{n}", float(self.bg_params[n]["value"]))
                        for n in bg_names]
            I_bg[valid] = eval_background(q_arr[valid], self.bg_shape, coeffs)
            adaptive_bg_valid = None

        # Full model = bg + peaks on valid points
        f_full = self._make_model_func(
            free_tags, full_vals,
            adaptive_bg=adaptive_bg_valid if self.bg_shape in BG_ADAPTIVE else None,
        )
        I_model[valid] = f_full(q_arr[valid])

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
