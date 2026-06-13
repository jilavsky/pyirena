"""Feature detection for small-angle scattering I(Q) curves.

This module identifies characteristic features in a log-log I(Q) profile —
Guinier plateaus, structure-factor peaks, power-law regions, and high-Q
background — by analysing the smoothed slope d(log I)/d(log Q) in log(Q)
space.  All thresholds are expressed in *log decades* so the algorithm is
independent of point count and works equally for SAXS and USAXS data.

The detector is intended to:

* Help users decide how many Unified Fit levels a dataset needs.
* Suggest initial Q-windows for Guinier and power-law sub-fits.
* Give the pyirena-ai LLM agent a structured, numerical view of the curve
  before it issues `set_parameter_value` / `run_fit` calls.

The module has **no GUI dependencies** — it returns a plain dataclass result
that is JSON-serialisable for the MCP layer.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class FeatureDetectConfig:
    """Tunable parameters for `detect_features`.

    All width / span values are in *log10(Q) decades*; slopes are
    d(log10 I)/d(log10 Q) (dimensionless).
    """

    sigma_smooth: float = 0.15
    """Gaussian smoothing width for log(I), in log(Q) decades."""

    span_deriv: float = 0.30
    """Total span (±span/2) for central-difference slope estimation, decades."""

    plateau_slope_max: float = 0.3
    """A plateau is any region with |slope| < this value."""

    plateau_min_width: float = 0.2
    """Minimum plateau width in decades to be reported."""

    peak_min_slope_mag: float = 0.25
    """Required |slope| on both sides of a peak zero-crossing."""

    peak_min_side_span: float = 0.1
    """Currently unused (reserved for future side-span check)."""

    power_law_max_slope: float = -0.5
    """A power-law region is any region with slope < this value."""

    power_law_min_width: float = 0.3
    """Minimum power-law region width in decades to be reported."""

    power_law_slope_stability: float = 0.25
    """Maximum std-dev of slope within a power-law region.

    Slightly looser than the original spec (0.2) so brief perturbations
    (e.g. from filtered-out noise points creating a small gap) do not
    fragment an otherwise stable power-law region.
    """

    noisy_point_threshold: float = 0.5
    """Drop points with sigma_I/I > this before smoothing."""

    @classmethod
    def saxs_preset(cls) -> "FeatureDetectConfig":
        """Defaults — tuned for SAXS data covering ~1-2 log decades."""
        return cls()

    @classmethod
    def usaxs_preset(cls) -> "FeatureDetectConfig":
        """Relaxed thresholds for USAXS data (>2.5 log decades, steeper slopes)."""
        return cls(
            plateau_slope_max=0.5,
            plateau_min_width=0.3,
            power_law_max_slope=-1.0,
        )

    @classmethod
    def auto(cls, q: np.ndarray) -> "FeatureDetectConfig":
        """Pick saxs_preset or usaxs_preset based on log10(q.max()/q.min())."""
        q = np.asarray(q, dtype=float)
        q_pos = q[(q > 0) & np.isfinite(q)]
        if q_pos.size < 2:
            return cls.saxs_preset()
        decades = float(np.log10(q_pos.max() / q_pos.min()))
        return cls.usaxs_preset() if decades > 2.5 else cls.saxs_preset()


# ---------------------------------------------------------------------------
# Result
# ---------------------------------------------------------------------------

@dataclass
class FeatureDetectResult:
    plateaus: list[dict] = field(default_factory=list)
    peaks: list[dict] = field(default_factory=list)
    power_law_regions: list[dict] = field(default_factory=list)
    background_q_min: Optional[float] = None
    recommended_guinier_windows: list[dict] = field(default_factory=list)
    recommended_nlevels: int = 0
    log_decades: float = 0.0
    n_points: int = 0
    preset_used: str = "custom"

    def to_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# Algorithm
# ---------------------------------------------------------------------------

def _gaussian_smooth_log(logQ: np.ndarray, logI: np.ndarray,
                         sigma: float) -> np.ndarray:
    """Gaussian-weighted smoothing of logI along logQ axis.

    Uses reflection padding across both ends so points near the boundary are
    smoothed against a symmetric kernel — without this, the smoother biases
    edge values toward the interior and corrupts slope estimates there.

    Vectorised O(N^2) — fine for typical SAS data (N ~ 200-2000).
    """
    n = logQ.size
    if n == 0:
        return logI
    # Reflection padding: extend logQ by mirroring across logQ[0] and logQ[-1],
    # and mirror logI to match.  Use a pad width of 3σ measured in points,
    # capped at n-1, so the kernel sees enough support.
    if n >= 3:
        dlog_med = float(np.median(np.diff(logQ)))
        pad = min(max(int(np.ceil(3.0 * sigma / max(dlog_med, 1e-9))), 1), n - 1)
    else:
        pad = 0
    if pad > 0:
        # logQ: reflect across the boundary point
        left_logQ  = 2 * logQ[0]  - logQ[1:pad + 1][::-1]
        right_logQ = 2 * logQ[-1] - logQ[-pad - 1:-1][::-1]
        # logI: linear extrapolation in (logQ, logI) — fit a line to the
        # first/last `pad` interior points and extend it across the boundary.
        # This preserves the local slope at the edges, which simple reflection
        # destroys (reflection makes the smoother see a symmetric V/Λ shape
        # at the boundary, biasing the slope toward zero).
        nfit = min(pad + 1, n)
        slope_l, intercept_l = np.polyfit(logQ[:nfit], logI[:nfit], 1)
        slope_r, intercept_r = np.polyfit(logQ[-nfit:], logI[-nfit:], 1)
        left_logI  = slope_l * left_logQ  + intercept_l
        right_logI = slope_r * right_logQ + intercept_r
        ext_logQ = np.concatenate([left_logQ, logQ, right_logQ])
        ext_logI = np.concatenate([left_logI, logI, right_logI])
    else:
        ext_logQ = logQ
        ext_logI = logI
    # Smooth each original point against the extended array
    d = logQ[:, None] - ext_logQ[None, :]
    w = np.exp(-(d ** 2) / (2.0 * sigma ** 2))
    w_sum = w.sum(axis=1)
    return (w @ ext_logI) / w_sum


def _local_slope(logQ: np.ndarray, logI_s: np.ndarray,
                 span: float) -> np.ndarray:
    """Central-difference slope over ±span/2 in log(Q) space.

    Falls back to numpy.gradient at indices where the symmetric window
    cannot be formed (near the data ends).
    """
    n = logQ.size
    slope = np.full(n, np.nan)
    half = span / 2.0
    for i in range(n):
        lo_target = logQ[i] - half
        hi_target = logQ[i] + half
        lo = int(np.searchsorted(logQ, lo_target, side="left"))
        hi = int(np.searchsorted(logQ, hi_target, side="right")) - 1
        lo = max(lo, 0)
        hi = min(hi, n - 1)
        if hi > lo and logQ[hi] > logQ[lo]:
            slope[i] = (logI_s[hi] - logI_s[lo]) / (logQ[hi] - logQ[lo])
    # Fill any remaining NaN with numpy.gradient as fallback
    if np.any(np.isnan(slope)) and n >= 2:
        grad = np.gradient(logI_s, logQ)
        nan_mask = np.isnan(slope)
        slope[nan_mask] = grad[nan_mask]
    return slope


def _find_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Return [(start, end_exclusive)] for each contiguous True run in mask."""
    if mask.size == 0:
        return []
    runs = []
    in_run = False
    start = 0
    for i, m in enumerate(mask):
        if m and not in_run:
            start = i
            in_run = True
        elif not m and in_run:
            runs.append((start, i))
            in_run = False
    if in_run:
        runs.append((start, mask.size))
    return runs


def detect_features(
    q: np.ndarray,
    I: np.ndarray,
    sigma_I: Optional[np.ndarray] = None,
    config: Optional[FeatureDetectConfig] = None,
) -> FeatureDetectResult:
    """Detect plateaus, peaks, and power-law regions in I(Q).

    Parameters
    ----------
    q : array of Å⁻¹
    I : array of intensity in same units as the user wants reported
    sigma_I : optional array of intensity uncertainty
    config : FeatureDetectConfig; defaults to FeatureDetectConfig.auto(q)

    Returns
    -------
    FeatureDetectResult
    """
    if config is None:
        config = FeatureDetectConfig.auto(q)

    q = np.asarray(q, dtype=float)
    I = np.asarray(I, dtype=float)

    # Filter: positive Q, positive I, finite, optionally drop noisy points
    mask = (q > 0) & (I > 0) & np.isfinite(q) & np.isfinite(I)
    if sigma_I is not None:
        sigma_I = np.asarray(sigma_I, dtype=float)
        if sigma_I.shape == I.shape:
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = np.where(I > 0, sigma_I / I, np.inf)
            mask &= rel < config.noisy_point_threshold

    q = q[mask]
    I = I[mask]
    if q.size < 5:
        return FeatureDetectResult(n_points=int(q.size))

    # Ensure ascending Q
    order = np.argsort(q)
    q = q[order]
    I = I[order]

    logQ = np.log10(q)
    logI = np.log10(I)

    logI_s = _gaussian_smooth_log(logQ, logI, config.sigma_smooth)
    slope = _local_slope(logQ, logI_s, config.span_deriv)

    log_decades = float(logQ[-1] - logQ[0])

    # Detect plateaus
    plateau_mask = np.abs(slope) < config.plateau_slope_max
    plateaus: list[dict] = []
    for s, e in _find_runs(plateau_mask):
        width = logQ[e - 1] - logQ[s]
        if width >= config.plateau_min_width:
            plateaus.append({
                "q_min":      float(10 ** logQ[s]),
                "q_max":      float(10 ** logQ[e - 1]),
                "q_center":   float(10 ** (0.5 * (logQ[s] + logQ[e - 1]))),
                "mean_slope": float(np.mean(slope[s:e])),
                "width_decades": float(width),
            })

    # Detect power-law regions.  Strategy:
    #   1. Find candidate runs where slope < power_law_max_slope.
    #   2. Inside each candidate, find the longest contiguous sub-run where
    #      every point is within ±slope_stability of the run's local median.
    #      This separates a clean Porod tail from the curved knee transition
    #      that precedes it.
    pl_mask = slope < config.power_law_max_slope
    power_law_regions: list[dict] = []
    tol = config.power_law_slope_stability
    for s, e in _find_runs(pl_mask):
        seg = slope[s:e]
        if seg.size < 3:
            continue
        med = float(np.median(seg[np.isfinite(seg)]))
        # Stable sub-mask within this candidate
        stable = np.abs(seg - med) < tol
        # Find longest run of `stable` within [s, e)
        best_start, best_end = -1, -1
        run_s = -1
        for k in range(stable.size + 1):
            if k < stable.size and stable[k]:
                if run_s < 0:
                    run_s = k
            else:
                if run_s >= 0:
                    if (k - run_s) > (best_end - best_start):
                        best_start, best_end = run_s, k
                    run_s = -1
        if best_start < 0:
            continue
        sub_s = s + best_start
        sub_e = s + best_end
        width = logQ[sub_e - 1] - logQ[sub_s]
        if width < config.power_law_min_width:
            continue
        sub_seg = slope[sub_s:sub_e]
        power_law_regions.append({
            "q_min":         float(10 ** logQ[sub_s]),
            "q_max":         float(10 ** logQ[sub_e - 1]),
            "slope":         float(np.mean(sub_seg)),
            "slope_std":     float(np.std(sub_seg)),
            "width_decades": float(width),
        })

    # Detect peaks — zero-crossings of slope from positive to negative
    # with |slope| > peak_min_slope_mag on both sides
    peaks: list[dict] = []
    for i in range(1, slope.size - 1):
        if (np.isfinite(slope[i - 1]) and np.isfinite(slope[i + 1])
                and slope[i - 1] > config.peak_min_slope_mag
                and slope[i + 1] < -config.peak_min_slope_mag):
            # Walk left until slope drops below peak_min_slope_mag
            j_lo = i - 1
            while j_lo > 0 and np.isfinite(slope[j_lo - 1]) \
                    and slope[j_lo - 1] > config.peak_min_slope_mag:
                j_lo -= 1
            # Walk right until slope rises above -peak_min_slope_mag
            j_hi = i + 1
            while j_hi < slope.size - 1 and np.isfinite(slope[j_hi + 1]) \
                    and slope[j_hi + 1] < -config.peak_min_slope_mag:
                j_hi += 1
            peaks.append({
                "q_peak":         float(q[i]),
                "q_low":          float(q[j_lo]),
                "q_high":         float(q[j_hi]),
                "slope_rising":   float(slope[i - 1]),
                "slope_falling":  float(slope[i + 1]),
            })

    # Deduplicate adjacent peaks (multiple zero-crossings in same neighbourhood)
    if peaks:
        deduped = [peaks[0]]
        for p in peaks[1:]:
            prev = deduped[-1]
            if p["q_peak"] > prev["q_peak"] * 1.05:
                deduped.append(p)
        peaks = deduped

    # Background: highest-Q plateau (if it touches the data's max Q)
    background_q_min: Optional[float] = None
    if plateaus:
        last = plateaus[-1]
        if abs(np.log10(last["q_max"]) - logQ[-1]) < 0.05:
            background_q_min = last["q_min"]

    # Recommended Guinier windows
    recommended: list[dict] = []
    for p in plateaus:
        if background_q_min is not None and p["q_min"] >= background_q_min:
            continue  # skip the background plateau
        q_max_guinier = p["q_max"] * 1.5
        recommended.append({
            "feature_type":     "plateau",
            "q_min_guinier":    float(p["q_min"]),
            "q_max_guinier":    float(q_max_guinier),
            "q_min_powerlaw":   float(q_max_guinier),
        })
    for pk in peaks:
        q_max_guinier = pk["q_high"] * 0.8
        recommended.append({
            "feature_type":     "peak",
            "q_min_guinier":    float(pk["q_peak"]),
            "q_max_guinier":    float(q_max_guinier),
            "q_min_powerlaw":   float(q_max_guinier),
        })

    return FeatureDetectResult(
        plateaus=plateaus,
        peaks=peaks,
        power_law_regions=power_law_regions,
        background_q_min=background_q_min,
        recommended_guinier_windows=recommended,
        recommended_nlevels=len(recommended),
        log_decades=log_decades,
        n_points=int(q.size),
        preset_used=_classify_preset(config),
    )


def _classify_preset(cfg: FeatureDetectConfig) -> str:
    """Best-effort label for a config object."""
    saxs = FeatureDetectConfig.saxs_preset()
    usaxs = FeatureDetectConfig.usaxs_preset()
    if (cfg.plateau_slope_max == saxs.plateau_slope_max
            and cfg.power_law_max_slope == saxs.power_law_max_slope
            and cfg.plateau_min_width == saxs.plateau_min_width):
        return "saxs"
    if (cfg.plateau_slope_max == usaxs.plateau_slope_max
            and cfg.power_law_max_slope == usaxs.power_law_max_slope
            and cfg.plateau_min_width == usaxs.plateau_min_width):
        return "usaxs"
    return "custom"


__all__ = [
    "FeatureDetectConfig",
    "FeatureDetectResult",
    "detect_features",
]
