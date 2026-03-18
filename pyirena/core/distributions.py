"""
Parametric size-distribution functions for Modeling tool.

Five distribution families are provided (matching the Igor Pro IR2L_Modeling
package from Argonne National Laboratory):

  * Gauss       — symmetric normal distribution
  * Log-Normal  — 3-parameter shifted log-normal
  * LSW         — Lifshitz-Slyozov-Wagner (Ostwald-ripening)
  * Schulz-Zimm — gamma-based distribution (2 parameters)
  * Ardell      — extended LSW by Ardell (2 < m < 3)

Each distribution is available as:
  - ``{name}_pdf(r, ...)``  — probability density (vectorised over r)
  - ``{name}_cdf(r, ...)``  — cumulative distribution (scalar or vectorised)

Plus two dispatchers ``pdf()`` / ``cdf()`` that route by string ``dist_type``.

The key public function is ``generate_radius_grid()`` which implements the
Igor-style non-linear binning algorithm (``IR2L_GenerateRadiiDist``):
  - Find r_min / r_max where CDF crosses tail_precision / 1-tail_precision
  - Build a dense linear grid of 3*n_bins diameters over [r_min, r_max]
  - Compute CDF at every point
  - Interpolate inverse-CDF at n_bins equally-spaced cumulative targets
  - Result: bins are densely spaced near the mode, sparse at the tails
"""

from __future__ import annotations

import math
import warnings
from typing import Union

import numpy as np
from scipy import stats, integrate
from scipy.integrate import cumulative_trapezoid as _cumtrapz


# ──────────────────────────────────────────────────────────────────────────────
# Gauss  (normal distribution)
# Parameters: mean_size [Å], width [Å]
# ──────────────────────────────────────────────────────────────────────────────

def gauss_pdf(r: np.ndarray, mean_size: float, width: float) -> np.ndarray:
    """Normal probability density.

    f(r) = exp(-(r-mean)²/(2*width²)) / (width*sqrt(2π))

    Args:
        r: radius values [Å]
        mean_size: distribution mean [Å]
        width: standard deviation [Å]  (must be > 0)

    Returns:
        probability density [Å⁻¹], same shape as r
    """
    r = np.asarray(r, dtype=float)
    return stats.norm.pdf(r, loc=mean_size, scale=max(width, 1e-30))


def gauss_cdf(r: Union[float, np.ndarray], mean_size: float, width: float
              ) -> Union[float, np.ndarray]:
    """Normal cumulative distribution."""
    return stats.norm.cdf(r, loc=mean_size, scale=max(width, 1e-30))


# ──────────────────────────────────────────────────────────────────────────────
# Log-Normal  (3-parameter shifted)
# Parameters: min_size [Å], mean_size [Å], sdeviation (σ in log-space, >0)
# Distribution of (r - min_size) is log-normal with geometric mean mean_size
# and log-space standard deviation sdeviation.
# ──────────────────────────────────────────────────────────────────────────────

def lognormal_pdf(r: np.ndarray, min_size: float, mean_size: float,
                  sdeviation: float) -> np.ndarray:
    """Shifted 3-parameter log-normal probability density.

    Let x = r - min_size.  Then x follows a log-normal distribution with:
        μ_log = ln(mean_size)  and  σ_log = sdeviation.

    f(r) = 1 / ((r - min_size) * sdeviation * √2π)
           * exp(-(ln((r - min_size)/mean_size))² / (2*sdeviation²))

    Args:
        r: radius values [Å]
        min_size: minimum radius offset [Å]; f=0 for r ≤ min_size
        mean_size: geometric mean of (r - min_size) [Å]  (must be > 0)
        sdeviation: log-space standard deviation  (must be > 0)

    Returns:
        probability density [Å⁻¹], same shape as r
    """
    r = np.asarray(r, dtype=float)
    x = r - float(min_size)
    out = np.zeros_like(r)
    mask = x > 0
    if np.any(mask):
        out[mask] = stats.lognorm.pdf(
            x[mask],
            s=max(sdeviation, 1e-10),
            scale=max(mean_size, 1e-30),
        )
    return out


def lognormal_cdf(r: Union[float, np.ndarray], min_size: float,
                  mean_size: float, sdeviation: float
                  ) -> Union[float, np.ndarray]:
    """Shifted 3-parameter log-normal cumulative distribution."""
    scalar = np.ndim(r) == 0
    r = np.atleast_1d(np.asarray(r, dtype=float))
    x = r - float(min_size)
    out = np.zeros_like(r)
    mask = x > 0
    if np.any(mask):
        out[mask] = stats.lognorm.cdf(
            x[mask],
            s=max(sdeviation, 1e-10),
            scale=max(mean_size, 1e-30),
        )
    return float(out[0]) if scalar else out


# ──────────────────────────────────────────────────────────────────────────────
# LSW  (Lifshitz-Slyozov-Wagner)
# Parameters: location [Å]  — the critical radius (= mean radius in LSW theory)
# Distribution is defined only for  0 < r < 1.5 * location
# ──────────────────────────────────────────────────────────────────────────────

def _lsw_pdf_scalar(r: float, location: float) -> float:
    """LSW probability density for a single r value."""
    if location <= 0:
        return 0.0
    u = r / location
    if u <= 0 or u >= 1.5:
        return 0.0
    # Classical LSW distribution (matrix-diffusion-controlled coarsening):
    # f(u) = (4/9)*u²*(3/(3+u))^(7/3) * exp(-u/(1.5-u)) / (1.5-u)^(11/3)
    # This is the form matched to the Igor IR1_LSWProbability convention.
    try:
        a = (4.0 / 9.0) * u ** 2
        b = (3.0 / (3.0 + u)) ** (7.0 / 3.0)
        denom = 1.5 - u
        if denom < 1e-12:
            return 0.0
        c = math.exp(-u / denom) / denom ** (11.0 / 3.0)
        return a * b * c / location   # Jacobian 1/location converts u→r
    except (OverflowError, ZeroDivisionError, ValueError):
        return 0.0


def lsw_pdf(r: np.ndarray, location: float) -> np.ndarray:
    """Lifshitz-Slyozov-Wagner probability density.

    Defined for 0 < r < 1.5 * location; zero outside this range.

    Args:
        r: radius values [Å]
        location: critical / mean radius [Å]

    Returns:
        probability density [Å⁻¹], same shape as r
    """
    r = np.asarray(r, dtype=float)
    return np.vectorize(_lsw_pdf_scalar)(r, location)


def lsw_cdf(r: Union[float, np.ndarray], location: float
            ) -> Union[float, np.ndarray]:
    """LSW cumulative distribution — numerical integration up to r."""
    scalar = np.ndim(r) == 0
    r_arr = np.atleast_1d(np.asarray(r, dtype=float))
    out = np.empty_like(r_arr)
    r_max = 1.5 * location
    for i, ri in enumerate(r_arr):
        if ri <= 0:
            out[i] = 0.0
        elif ri >= r_max:
            out[i] = 1.0
        else:
            val, _ = integrate.quad(
                _lsw_pdf_scalar, 0.0, ri,
                args=(location,),
                limit=100,
            )
            out[i] = min(max(val, 0.0), 1.0)
    return float(out[0]) if scalar else out


# ──────────────────────────────────────────────────────────────────────────────
# Schulz-Zimm  (gamma distribution parametrized by mean and width)
# Parameters: mean_size [Å], width [Å]
# Equivalent to a Gamma distribution with Z = (mean/width)² - 1
# ──────────────────────────────────────────────────────────────────────────────

def _schulz_zimm_ab(mean_size: float, width: float):
    """Return (a, scale) for scipy.stats.gamma.

    Z + 1 = a = (mean/width)²  (Z ≥ 0 requires width ≤ mean)
    scale = mean / a
    """
    a = max((mean_size / max(width, 1e-10)) ** 2, 0.01)  # Z+1
    scale = mean_size / a
    return a, scale


def schulz_zimm_pdf(r: np.ndarray, mean_size: float, width: float
                    ) -> np.ndarray:
    """Schulz-Zimm (gamma) probability density.

    Uses Z = (mean/width)² - 1 with the Gamma distribution:
      f(r) = Gamma(Z+1, scale=mean/(Z+1)).pdf(r)

    Args:
        r: radius values [Å]
        mean_size: distribution mean [Å]
        width: standard deviation [Å]

    Returns:
        probability density [Å⁻¹], same shape as r
    """
    r = np.asarray(r, dtype=float)
    a, scale = _schulz_zimm_ab(mean_size, width)
    return stats.gamma.pdf(r, a=a, scale=scale)


def schulz_zimm_cdf(r: Union[float, np.ndarray], mean_size: float,
                    width: float) -> Union[float, np.ndarray]:
    """Schulz-Zimm cumulative distribution."""
    a, scale = _schulz_zimm_ab(mean_size, width)
    return stats.gamma.cdf(r, a=a, scale=scale)


# ──────────────────────────────────────────────────────────────────────────────
# Ardell  (extended LSW by Ardell 1972)
# Parameters: location [Å], parameter m  (must satisfy 2 < m < 3)
#
# Implements Igor's IR1_ArdellProbability / IR1_ArdellProbNormalized exactly:
#
#   F(zp, m) = (zp*(m-1))^(m-1) / [m^m*(zp-1) − zp^m*(m-1)^(m-1)]
#   P(zp, m) = ∫₀^zp F(t, m) dt          (F < 0 everywhere → P < 0)
#   f(zp)    = −3 · F(zp) · exp(P(zp))   (unnormalized; > 0)
#
# The distribution is defined for 0 < zp < zp_max = m/(m-1):
#   m=2 → zp_max=2.0,  m=2.5 → zp_max=1.667,  m=3 → zp_max=1.5 (LSW limit)
#
# The modal radius is ~1.3–1.4 × location; the hard cutoff is at
# zp_max × location with a sharp (exponential) drop just below it.
# ──────────────────────────────────────────────────────────────────────────────

def _ardell_umax(m: float) -> float:
    """Upper limit of normalised radius: zp_max = m / (m-1).

    Derived from the root of the denominator of Ardell_F.
    Checked against Igor: m=3 → 1.5 (LSW limit), m=2 → 2.0.
    """
    return m / (m - 1.0)


def _ardell_F_vec(zp: np.ndarray, m: float) -> np.ndarray:
    """Vectorised Igor Ardell_F(zp, m).

    F(zp, m) = (zp*(m-1))^(m-1) / [m^m*(zp-1) − zp^m*(m-1)^(m-1)]

    F is negative throughout 0 < zp < zp_max; F(0)=0; F→−∞ at zp→zp_max.
    """
    mm  = m ** m
    cm1 = (m - 1.0) ** (m - 1.0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        num   = (zp * (m - 1.0)) ** (m - 1.0)
        denom = mm * (zp - 1.0) - zp ** m * cm1
        F = np.where(np.abs(denom) > 1e-200, num / denom, 0.0)
    return np.where(np.isfinite(F), F, 0.0)


# Module-level cache keyed on rounded m value — avoids recomputing the grid
# on every PDF/CDF call during fitting and radius-grid generation.
_ARDELL_GRID_CACHE: dict[float, tuple] = {}


def _ardell_build_grid(m: float, n: int = 3000):
    """Pre-compute normalised Ardell PDF and CDF on a fine zp-grid.

    Returns (zp_grid, pdf_norm, cdf_norm).  Results are cached by m.
    """
    m_key = round(m, 5)
    if m_key in _ARDELL_GRID_CACHE:
        return _ARDELL_GRID_CACHE[m_key]

    zp_max = _ardell_umax(m)
    # Fine grid from 0 to just below the singularity at zp_max
    zp_grid = np.linspace(0.0, zp_max * (1.0 - 1e-6), n)

    F_grid = _ardell_F_vec(zp_grid, m)

    # P(zp) = ∫₀^zp F(t) dt  (cumulative trapezoid, initial value = 0)
    P_grid = _cumtrapz(F_grid, zp_grid, initial=0.0)

    # Unnormalized PDF: f(zp) = −3·F(zp)·exp(P(zp))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        f_unnorm = -3.0 * F_grid * np.exp(P_grid)
    f_unnorm = np.where(np.isfinite(f_unnorm) & (f_unnorm >= 0.0), f_unnorm, 0.0)

    # Normalise so ∫ f(zp) dzp = 1
    norm = float(np.trapezoid(f_unnorm, zp_grid))
    pdf_norm = f_unnorm / norm if norm > 0.0 else np.ones_like(zp_grid) / zp_grid[-1]

    # CDF on the same grid
    cdf_norm = _cumtrapz(pdf_norm, zp_grid, initial=0.0)
    cdf_norm = np.clip(cdf_norm, 0.0, 1.0)

    result = (zp_grid, pdf_norm, cdf_norm)
    _ARDELL_GRID_CACHE[m_key] = result
    return result


def ardell_pdf(r: np.ndarray, location: float, parameter: float
               ) -> np.ndarray:
    """Ardell extended-LSW probability density (port of Igor IR1_ArdellProbNormalized).

    f(r) = −3·F(r/r₀)·exp(P(r/r₀)) / r₀ / normalization

    where F and P are defined in _ardell_F_vec and _ardell_build_grid.

    Args:
        r: radius values [Å]
        location: scale radius r₀ [Å];  distribution cuts off at r₀·m/(m-1)
        parameter: Ardell m-parameter (strictly between 2 and 3)

    Returns:
        probability density [Å⁻¹], same shape as r
    """
    r   = np.asarray(r, dtype=float)
    loc = max(float(location), 1e-10)
    m   = max(2.001, min(2.999, float(parameter)))

    zp_grid, pdf_norm, _ = _ardell_build_grid(m)

    # Interpolate; apply Jacobian 1/loc to convert zp-density to r-density
    zp_r = r / loc
    return np.interp(zp_r, zp_grid, pdf_norm, left=0.0, right=0.0) / loc


def ardell_cdf(r: Union[float, np.ndarray], location: float,
               parameter: float) -> Union[float, np.ndarray]:
    """Ardell cumulative distribution (vectorised, grid-based)."""
    scalar = np.ndim(r) == 0
    r_arr  = np.atleast_1d(np.asarray(r, dtype=float))
    loc    = max(float(location), 1e-10)
    m      = max(2.001, min(2.999, float(parameter)))

    zp_grid, _, cdf_norm = _ardell_build_grid(m)

    zp_r = r_arr / loc
    out  = np.interp(zp_r, zp_grid, cdf_norm, left=0.0, right=1.0)
    return float(out[0]) if scalar else out


# ──────────────────────────────────────────────────────────────────────────────
# Dispatchers
# ──────────────────────────────────────────────────────────────────────────────

#: Supported distribution type strings and their required parameter keys.
DIST_PARAM_NAMES: dict[str, list[str]] = {
    'gauss':       ['mean_size', 'width'],
    'lognormal':   ['min_size', 'mean_size', 'sdeviation'],
    'lsw':         ['location'],
    'schulz_zimm': ['mean_size', 'width'],
    'ardell':      ['location', 'parameter'],
}

#: Human-readable labels for each distribution.
DIST_LABELS: dict[str, str] = {
    'gauss':       'Gauss',
    'lognormal':   'Log-Normal',
    'lsw':         'LSW',
    'schulz_zimm': 'Schulz-Zimm',
    'ardell':      'Ardell',
}

#: Default parameter values for each distribution.
DIST_DEFAULTS: dict[str, dict] = {
    'gauss':       {'mean_size': 100.0, 'width': 20.0},
    'lognormal':   {'min_size': 10.0, 'mean_size': 100.0, 'sdeviation': 0.3},
    'lsw':         {'location': 100.0},
    'schulz_zimm': {'mean_size': 100.0, 'width': 30.0},
    'ardell':      {'location': 100.0, 'parameter': 2.5},
}


def pdf(dist_type: str, r: np.ndarray, params: dict) -> np.ndarray:
    """Dispatch to the appropriate PDF by distribution type string.

    Args:
        dist_type: one of 'gauss', 'lognormal', 'lsw', 'schulz_zimm', 'ardell'
        r: radius array [Å]
        params: dict of distribution parameters (see DIST_PARAM_NAMES)

    Returns:
        probability density array [Å⁻¹]
    """
    r = np.asarray(r, dtype=float)
    t = dist_type.lower()
    if t == 'gauss':
        return gauss_pdf(r, params['mean_size'], params['width'])
    if t == 'lognormal':
        return lognormal_pdf(r, params['min_size'], params['mean_size'],
                             params['sdeviation'])
    if t == 'lsw':
        return lsw_pdf(r, params['location'])
    if t == 'schulz_zimm':
        return schulz_zimm_pdf(r, params['mean_size'], params['width'])
    if t == 'ardell':
        return ardell_pdf(r, params['location'], params['parameter'])
    raise ValueError(f"Unknown distribution type: {dist_type!r}. "
                     f"Supported: {sorted(DIST_PARAM_NAMES)}")


def cdf(dist_type: str, r, params: dict):
    """Dispatch to the appropriate CDF by distribution type string.

    Args:
        dist_type: one of 'gauss', 'lognormal', 'lsw', 'schulz_zimm', 'ardell'
        r: scalar or array of radius values [Å]
        params: dict of distribution parameters

    Returns:
        cumulative probability (scalar or array matching r)
    """
    t = dist_type.lower()
    if t == 'gauss':
        return gauss_cdf(r, params['mean_size'], params['width'])
    if t == 'lognormal':
        return lognormal_cdf(r, params['min_size'], params['mean_size'],
                             params['sdeviation'])
    if t == 'lsw':
        return lsw_cdf(r, params['location'])
    if t == 'schulz_zimm':
        return schulz_zimm_cdf(r, params['mean_size'], params['width'])
    if t == 'ardell':
        return ardell_cdf(r, params['location'], params['parameter'])
    raise ValueError(f"Unknown distribution type: {dist_type!r}. "
                     f"Supported: {sorted(DIST_PARAM_NAMES)}")


# ──────────────────────────────────────────────────────────────────────────────
# Igor-style non-linear radius grid
# ──────────────────────────────────────────────────────────────────────────────

def _mode_estimate(dist_type: str, params: dict) -> float:
    """Return a starting estimate of the distribution mode (peak location)."""
    t = dist_type.lower()
    if t == 'gauss':
        return params['mean_size']
    if t == 'lognormal':
        mn = params['min_size']
        ms = params['mean_size']
        s = params['sdeviation']
        # Mode of shifted log-normal: min_size + mean_size * exp(-s²)
        return mn + ms * math.exp(-s ** 2)
    if t == 'lsw':
        return params['location']          # peak is near location
    if t == 'schulz_zimm':
        # Mode of Gamma(a, scale): (a-1)*scale = mean*(1 - 1/a) ≈ mean for large a
        a, scale = _schulz_zimm_ab(params['mean_size'], params['width'])
        return max((a - 1.0) * scale, params['mean_size'] * 0.5)
    if t == 'ardell':
        return params['location']
    raise ValueError(f"Unknown distribution type: {dist_type!r}")


def _step_size(dist_type: str, params: dict) -> float:
    """Return an initial step size for the CDF-walking algorithm."""
    t = dist_type.lower()
    if t == 'gauss':
        return params['width'] * 0.02
    if t == 'lognormal':
        s = params['sdeviation']
        ms = params['mean_size']
        # Standard deviation of the log-normal: ms*sqrt(exp(s²)-1)*exp(s²/2)
        std = ms * math.sqrt(max(math.exp(s ** 2) - 1, 1e-10)) * math.exp(s ** 2 / 2)
        return max(std * 0.02, params['min_size'] * 0.01 + 0.1)
    if t == 'lsw':
        return params['location'] * 0.05
    if t == 'schulz_zimm':
        return params['width'] * 0.02
    if t == 'ardell':
        return params['location'] * 0.02
    raise ValueError(f"Unknown distribution type: {dist_type!r}")


def generate_radius_grid(
    dist_type: str,
    params: dict,
    n_bins: int = 200,
    tail_precision: float = 0.01,
) -> np.ndarray:
    """Generate a non-linear radius bin grid for size distribution calculations.

    Implements the Igor Pro ``IR2L_GenerateRadiiDist`` algorithm:

    1. Start at the distribution mode.
    2. Walk downward in steps until CDF < ``tail_precision``  → r_min
    3. Walk upward in steps until CDF > 1 - ``tail_precision`` → r_max
    4. Create a dense grid of 3*n_bins equally-spaced points in [r_min, r_max].
    5. Evaluate the CDF at every point.
    6. Generate n_bins equally-spaced cumulative probability targets
       from p_min to 1 - tail_precision.
    7. Invert (interpolate) the CDF to find the radius for each target.

    The result has bins that are densely spaced near the peak and sparse at
    the tails — efficient for accurate intensity integration.

    Args:
        dist_type:       Distribution type string (see DIST_PARAM_NAMES)
        params:          Distribution parameter dict
        n_bins:          Number of radius bins (default 200)
        tail_precision:  Fraction of probability to ignore at each tail (default 0.01 = 1 %)

    Returns:
        1-D array of radius bin centres [Å], shape (n_bins,)
    """
    if n_bins < 2:
        raise ValueError(f"n_bins must be ≥ 2, got {n_bins}")

    mode = _mode_estimate(dist_type, params)
    step = max(_step_size(dist_type, params), 1e-3)

    r_min_hard = 1.0          # minimum possible radius (1 Å)
    r_max_hard = 1e15         # maximum physical radius

    # ── Find r_min (walk downward from mode) ────────────────────────────────
    x = mode
    p_min = tail_precision
    while True:
        x -= step
        if x < r_min_hard:
            x = r_min_hard
            break
        c = cdf(dist_type, x, params)
        if c < p_min or x <= r_min_hard:
            break
    r_min = x

    # Recalculate actual p_min in case we hit the hard floor
    p_start = float(cdf(dist_type, r_min, params))
    if p_start >= tail_precision:
        p_start = tail_precision        # use nominal value

    # ── Find r_max (walk upward from mode) ──────────────────────────────────
    x = mode
    p_max_target = 1.0 - tail_precision
    while True:
        x += step
        if x > r_max_hard:
            x = r_max_hard
            break
        c = cdf(dist_type, x, params)
        if c >= p_max_target or x >= r_max_hard:
            break
        # Extra safety for LSW/Ardell that have hard upper bounds
        if dist_type.lower() == 'lsw' and x > 1.499 * params['location']:
            x = 1.499 * params['location']
            break
        if dist_type.lower() == 'ardell':
            u_max = _ardell_umax(float(params['parameter']))
            if x > 0.999 * u_max * params['location']:
                x = 0.999 * u_max * params['location']
                break
    r_max = x

    if r_max <= r_min:
        # Fallback: linear grid centred on mode
        r_min = max(mode * 0.01, r_min_hard)
        r_max = mode * 10.0

    # ── Dense grid + CDF evaluation ─────────────────────────────────────────
    n_dense = max(3 * n_bins, 300)
    r_dense = np.linspace(r_min, r_max, n_dense)
    cdf_dense = np.asarray(cdf(dist_type, r_dense, params), dtype=float)

    # Ensure CDF is monotonically non-decreasing (numerical noise fix)
    cdf_dense = np.maximum.accumulate(cdf_dense)
    # Clamp to [0, 1]
    cdf_dense = np.clip(cdf_dense, 0.0, 1.0)

    # ── Equally-spaced probability targets ──────────────────────────────────
    p_end = float(cdf_dense[-1])
    if p_end <= p_start:
        # Degenerate: return linear grid
        return np.linspace(max(r_min, r_min_hard), r_max, n_bins)

    targets = np.linspace(p_start, p_end, n_bins)

    # ── Invert CDF: find radius for each probability target ─────────────────
    r_grid = np.interp(targets, cdf_dense, r_dense)

    # Guard against NaN / non-finite
    if not np.all(np.isfinite(r_grid)):
        warnings.warn(
            f"generate_radius_grid({dist_type}): some bins are non-finite; "
            "falling back to linear grid",
            RuntimeWarning,
        )
        r_grid = np.linspace(r_min, r_max, n_bins)

    return r_grid
