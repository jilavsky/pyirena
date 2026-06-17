"""Feature detection for small-angle scattering I(Q) curves.

Segments an I(Q) log-log curve into power-law-slope regions
(`d(log I)/d(log Q)` locally constant) and derives Guinier knees as the
transitions between adjacent stable-slope segments with substantially
different slopes.

Conceptual model
----------------
Every region of a SAS curve has a *local power-law slope*.  Steep slopes
(P ≈ 2-4) correspond to Porod tails; slope ≈ 0 corresponds to a Guinier
plateau or, at the highest-Q end, a constant background.  Therefore a
complete segmentation of the curve into piecewise-constant-slope segments
recovers all the information a Unified Fit user needs to choose levels and
Q-windows; Guinier knees fall naturally between adjacent slope segments.

Algorithm
---------
1. Filter positive-Q / positive-I points, optionally drop noisy points
   (`σ_I / I > noisy_point_threshold`).
2. Clip Q to ``config.q_max_clip`` (defaults to 0.6 Å⁻¹, the SAS-approximation
   limit; diffraction features beyond this should not be segmented).
3. Smooth ``log10(I)`` along ``log10(Q)`` with a Gaussian kernel of width
   ``sigma_smooth`` (in decades), using reflection + linear-extrapolation
   padding so edge slopes are preserved.
4. Estimate the local slope ``d(log I) / d(log Q)`` with a central difference
   over ``±span_deriv / 2`` decades.
5. Scan a sliding window of width ``stability_window`` decades along the
   slope profile; mark each point as **stable** if the slope std within the
   window is below ``stability_std_max``.
6. Group runs of stable points into raw segments; merge adjacent runs whose
   mean slopes agree within ``merge_slope_tol`` and whose gap is small.
7. Drop segments narrower than ``min_segment_decades``.
8. Classify each remaining segment by its mean slope: ``background``
   (|slope| small AND segment is rightmost AND touches the high-Q end),
   ``guinier_plateau`` (|slope| small elsewhere), or ``power_law``.
9. Derive Guinier knees as the unstable gaps between adjacent segments
   where the slope changes by at least ``guinier_knee_min_delta_slope``.

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

    # --- Smoothing & slope estimation ---
    sigma_smooth: float = 0.15
    """Gaussian smoothing width for log(I), in log(Q) decades."""

    span_deriv: float = 0.30
    """Total span (±span/2) for central-difference slope estimation, decades."""

    # --- Stability-window segmentation ---
    stability_window: float = 0.40
    """Width of the sliding window over which slope std is evaluated."""

    stability_std_max: float = 0.40
    """A window is "stable" if std(slope) within it is below this value.

    Tuned empirically against 31 hand-labelled SAXS curves: 0.40 hits
    ~95% of human-identified power-law-slope regions; tighter values
    (0.15-0.25) miss regions where the slope drifts smoothly.
    """

    merge_slope_tol: float = 0.40
    """Adjacent stable runs are merged if |mean_slope_left - mean_slope_right|
    is below this threshold (and the gap between them is small)."""

    merge_max_gap_decades: float = 0.20
    """Adjacent stable runs are only merged if the unstable gap between them
    is narrower than this (in decades)."""

    min_segment_decades: float = 0.20
    """Minimum width of a reported segment.  Narrower stable runs are dropped.

    Note: when the user sees Q-decades-based width here, it includes the
    smoothing window of ±sigma_smooth/2 on each side.  Practical minimum
    where slope can be confidently constant is ~2×sigma_smooth.
    """

    # --- Q clipping (SAS approximation limit) ---
    q_max_clip: Optional[float] = 0.6
    """Drop data with Q > this before segmentation (Å⁻¹).

    Defaults to 0.6, the practical upper limit of the small-angle
    approximation for X-rays.  Beyond this the curve typically contains
    amorphous diffraction features that should not be classified as SAS
    structure.  Pass ``None`` to disable.
    """

    # --- Classification thresholds (applied to mean slope per segment) ---
    background_slope_max: float = 0.25
    """|slope| < this AND segment touches data's high-Q end → background."""

    guinier_plateau_slope_max: float = 0.4
    """|slope| < this and not background → guinier_plateau."""

    guinier_knee_min_delta_slope: float = 0.5
    """Minimum slope jump between adjacent segments to emit a Guinier knee."""

    # --- Noise filtering ---
    noisy_point_threshold: float = 0.5
    """Drop points with sigma_I/I > this before smoothing."""


# ---------------------------------------------------------------------------
# Result
# ---------------------------------------------------------------------------

@dataclass
class FeatureDetectResult:
    segments: list[dict] = field(default_factory=list)
    """One entry per stable-slope segment, sorted by ascending Q.

    Each: ``q_min``, ``q_max``, ``slope``, ``slope_std``, ``kind``,
    ``intensity_mid``, ``width_decades``.

    ``kind`` is one of ``"background"``, ``"guinier_plateau"``, ``"power_law"``.

    Note: segments do not necessarily cover the full Q range — unstable
    transition regions (Guinier knees) are reported separately under
    ``guinier_knees``.
    """

    guinier_knees: list[dict] = field(default_factory=list)
    """Derived list of inferred Guinier knees between adjacent segments.

    Each: ``q_min``, ``q_max``, ``q_center``, ``slope_left``, ``slope_right``,
    ``delta_slope``.
    """

    background_q_min: Optional[float] = None
    recommended_guinier_windows: list[dict] = field(default_factory=list)
    recommended_nlevels: int = 0
    log_decades: float = 0.0
    n_points: int = 0
    q_min_analysed: Optional[float] = None
    q_max_analysed: Optional[float] = None
    n_segments_found: int = 0

    def to_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# Smoothing & slope (kept from v1 — unchanged)
# ---------------------------------------------------------------------------

def _gaussian_smooth_log(logQ: np.ndarray, logI: np.ndarray,
                         sigma: float) -> np.ndarray:
    """Gaussian-weighted smoothing of logI along logQ axis.

    Uses reflection + linear-extrapolation padding so points near the
    boundary are smoothed against a symmetric kernel — without this, the
    smoother biases edge values toward the interior and corrupts slope
    estimates there.

    Vectorised O(N²) — fine for typical SAS data (N ~ 200-2000).
    """
    n = logQ.size
    if n == 0:
        return logI
    if n >= 3:
        dlog_med = float(np.median(np.diff(logQ)))
        pad = min(max(int(np.ceil(3.0 * sigma / max(dlog_med, 1e-9))), 1), n - 1)
    else:
        pad = 0
    if pad > 0:
        left_logQ  = 2 * logQ[0]  - logQ[1:pad + 1][::-1]
        right_logQ = 2 * logQ[-1] - logQ[-pad - 1:-1][::-1]
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
    d = logQ[:, None] - ext_logQ[None, :]
    w = np.exp(-(d ** 2) / (2.0 * sigma ** 2))
    w_sum = w.sum(axis=1)
    return (w @ ext_logI) / w_sum


def _local_slope(logQ: np.ndarray, logI_s: np.ndarray,
                 span: float) -> np.ndarray:
    """Central-difference slope over ±span/2 in log(Q) space."""
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
    if np.any(np.isnan(slope)) and n >= 2:
        grad = np.gradient(logI_s, logQ)
        nan_mask = np.isnan(slope)
        slope[nan_mask] = grad[nan_mask]
    return slope


# ---------------------------------------------------------------------------
# Stability-window segmentation
# ---------------------------------------------------------------------------

def _stable_runs(logQ: np.ndarray, slope: np.ndarray,
                 window_decades: float, std_max: float) -> tuple[np.ndarray, np.ndarray]:
    """For each point, compute (is_stable, local_mean_slope).

    A point is stable if std(slope) within ``±window_decades/2`` < std_max.
    local_mean_slope is the mean within that window (NaN for unstable points).
    """
    n = logQ.size
    is_stable = np.zeros(n, dtype=bool)
    local_mean = np.full(n, np.nan)
    half = window_decades / 2.0
    for i in range(n):
        lo = int(np.searchsorted(logQ, logQ[i] - half, side="left"))
        hi = int(np.searchsorted(logQ, logQ[i] + half, side="right"))
        if hi - lo >= 3:
            seg = slope[lo:hi]
            std = float(np.std(seg))
            if std < std_max:
                is_stable[i] = True
                local_mean[i] = float(np.mean(seg))
    return is_stable, local_mean


def _runs_from_mask(mask: np.ndarray) -> list[tuple[int, int]]:
    """Return [(start, end_exclusive)] for each contiguous True run."""
    runs = []
    in_run = False
    start = 0
    for i in range(mask.size):
        if mask[i] and not in_run:
            start = i
            in_run = True
        elif not mask[i] and in_run:
            runs.append((start, i))
            in_run = False
    if in_run:
        runs.append((start, mask.size))
    return runs


def _merge_adjacent_runs(runs: list[tuple[int, int]], logQ: np.ndarray,
                         local_mean: np.ndarray,
                         slope_tol: float, max_gap_decades: float
                         ) -> list[tuple[int, int]]:
    """Merge adjacent runs whose mean slopes are similar and gap is small."""
    if not runs:
        return []
    merged = [runs[0]]
    for s, e in runs[1:]:
        prev_s, prev_e = merged[-1]
        prev_mean = float(np.nanmean(local_mean[prev_s:prev_e]))
        cur_mean = float(np.nanmean(local_mean[s:e]))
        gap = float(logQ[s] - logQ[prev_e - 1])
        if (abs(prev_mean - cur_mean) < slope_tol
                and gap < max_gap_decades):
            merged[-1] = (prev_s, e)
        else:
            merged.append((s, e))
    return merged


def _filter_by_width(runs: list[tuple[int, int]], logQ: np.ndarray,
                     min_decades: float) -> list[tuple[int, int]]:
    """Drop runs narrower than min_decades."""
    out = []
    for s, e in runs:
        width = float(logQ[e - 1] - logQ[s])
        if width >= min_decades:
            out.append((s, e))
    return out


# ---------------------------------------------------------------------------
# Classification & assembly
# ---------------------------------------------------------------------------

def _classify_segment(mean_slope: float, q_max_seg: float,
                      data_q_max: float, cfg: FeatureDetectConfig) -> str:
    """Classify a segment by its mean slope and position.

    A segment is *background* iff its |slope| is below background_slope_max
    AND it touches the high-Q end of the analysed data (last few percent).
    """
    abs_s = abs(mean_slope)
    touches_end = (np.log10(data_q_max) - np.log10(q_max_seg)) < 0.05
    if abs_s < cfg.background_slope_max and touches_end:
        return "background"
    if abs_s < cfg.guinier_plateau_slope_max:
        return "guinier_plateau"
    return "power_law"


def _build_segments(seg_ranges: list[tuple[int, int]],
                    q: np.ndarray, I: np.ndarray, slope: np.ndarray,
                    cfg: FeatureDetectConfig) -> list[dict]:
    if not seg_ranges:
        return []
    data_q_max = float(q[-1])
    out: list[dict] = []
    for s, e in seg_ranges:
        seg_slope = slope[s:e]
        seg_q = q[s:e]
        seg_I = I[s:e]
        mean_slope = float(np.mean(seg_slope))
        q_min_seg = float(seg_q[0])
        q_max_seg = float(seg_q[-1])
        kind = _classify_segment(mean_slope, q_max_seg, data_q_max, cfg)
        log_q_mid = 0.5 * (np.log10(q_min_seg) + np.log10(q_max_seg))
        mid_idx = int(np.argmin(np.abs(np.log10(seg_q) - log_q_mid)))
        out.append({
            "q_min":         q_min_seg,
            "q_max":         q_max_seg,
            "slope":         mean_slope,
            "slope_std":     float(np.std(seg_slope)),
            "kind":          kind,
            "intensity_mid": float(seg_I[mid_idx]),
            "width_decades": float(np.log10(q_max_seg) - np.log10(q_min_seg)),
        })
    return out


def _derive_knees(segments: list[dict], cfg: FeatureDetectConfig) -> list[dict]:
    knees: list[dict] = []
    for i in range(len(segments) - 1):
        left = segments[i]
        right = segments[i + 1]
        delta = abs(left["slope"] - right["slope"])
        if delta < cfg.guinier_knee_min_delta_slope:
            continue
        q_lo = float(left["q_max"])
        q_hi = float(right["q_min"])
        if q_hi < q_lo:
            q_hi = q_lo
        q_center = float(np.sqrt(max(q_lo * q_hi, 1e-30)))
        knees.append({
            "q_min":       q_lo,
            "q_max":       q_hi,
            "q_center":    q_center,
            "slope_left":  float(left["slope"]),
            "slope_right": float(right["slope"]),
            "delta_slope": float(delta),
        })
    return knees


def _recommended_guinier_windows(segments: list[dict],
                                 knees: list[dict]) -> list[dict]:
    """For each Guinier knee, propose a Q-window covering the transition.

    Window covers the lower-Q side of the knee (where the Guinier function
    dominates) up to slightly past the start of the higher-Q power-law segment.
    """
    out: list[dict] = []
    for knee in knees:
        right_seg = None
        for seg in segments:
            if seg["q_min"] >= knee["q_max"] and seg["kind"] == "power_law":
                right_seg = seg
                break
        if right_seg is None:
            continue
        left_seg = None
        for seg in segments:
            if seg["q_max"] <= knee["q_min"]:
                left_seg = seg
            else:
                break
        if left_seg is None:
            continue
        q_min_g = float(left_seg["q_min"])
        q_max_g = float(right_seg["q_min"] * 1.2)
        out.append({
            "feature_type":   "knee",
            "q_min_guinier":  q_min_g,
            "q_max_guinier":  q_max_g,
            "q_min_powerlaw": float(right_seg["q_min"]),
        })
    return out


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def detect_features(
    q: np.ndarray,
    I: np.ndarray,
    sigma_I: Optional[np.ndarray] = None,
    config: Optional[FeatureDetectConfig] = None,
) -> FeatureDetectResult:
    """Segment an I(Q) curve into power-law slope regions and derive knees.

    Parameters
    ----------
    q : array of Å⁻¹
    I : array of intensity
    sigma_I : optional array of intensity uncertainty
    config : FeatureDetectConfig; defaults to FeatureDetectConfig()

    Returns
    -------
    FeatureDetectResult
    """
    if config is None:
        config = FeatureDetectConfig()

    q = np.asarray(q, dtype=float)
    I = np.asarray(I, dtype=float)

    # 1. Mask & noise filter
    mask = (q > 0) & (I > 0) & np.isfinite(q) & np.isfinite(I)
    if sigma_I is not None:
        sigma_I = np.asarray(sigma_I, dtype=float)
        if sigma_I.shape == I.shape:
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = np.where(I > 0, sigma_I / I, np.inf)
            mask &= rel < config.noisy_point_threshold
            sigma_I = sigma_I[mask]
        else:
            sigma_I = None
    q = q[mask]
    I = I[mask]

    # 2. Clip to q_max_clip
    if config.q_max_clip is not None:
        clip_mask = q <= config.q_max_clip
        q = q[clip_mask]
        I = I[clip_mask]
        if sigma_I is not None:
            sigma_I = sigma_I[clip_mask]

    if q.size < 8:
        return FeatureDetectResult(n_points=int(q.size))

    # 3. Sort
    order = np.argsort(q)
    q = q[order]
    I = I[order]
    if sigma_I is not None:
        sigma_I = sigma_I[order]

    # 4. Smooth & slope
    logQ = np.log10(q)
    logI = np.log10(I)
    logI_s = _gaussian_smooth_log(logQ, logI, config.sigma_smooth)
    slope = _local_slope(logQ, logI_s, config.span_deriv)
    if np.any(~np.isfinite(slope)):
        idx_ok = np.where(np.isfinite(slope))[0]
        if idx_ok.size == 0:
            return FeatureDetectResult(n_points=int(q.size))
        slope = np.interp(np.arange(slope.size), idx_ok, slope[idx_ok])

    # 5. Stability windows → raw stable runs
    is_stable, local_mean = _stable_runs(
        logQ, slope,
        window_decades=config.stability_window,
        std_max=config.stability_std_max,
    )
    runs = _runs_from_mask(is_stable)

    # 6. Merge adjacent runs with similar mean slope
    runs = _merge_adjacent_runs(
        runs, logQ, local_mean,
        slope_tol=config.merge_slope_tol,
        max_gap_decades=config.merge_max_gap_decades,
    )

    # 7. Filter by minimum width
    runs = _filter_by_width(runs, logQ, config.min_segment_decades)

    # 8. Build segments + classify
    segments = _build_segments(runs, q, I, slope, config)

    # 9. Derive knees, windows, etc.
    knees = _derive_knees(segments, config)
    windows = _recommended_guinier_windows(segments, knees)
    background_q_min = None
    for seg in segments:
        if seg["kind"] == "background":
            background_q_min = seg["q_min"]
            break
    n_levels = sum(1 for seg in segments if seg["kind"] != "background")

    return FeatureDetectResult(
        segments=segments,
        guinier_knees=knees,
        background_q_min=background_q_min,
        recommended_guinier_windows=windows,
        recommended_nlevels=n_levels,
        log_decades=float(logQ[-1] - logQ[0]),
        n_points=int(q.size),
        q_min_analysed=float(q[0]),
        q_max_analysed=float(q[-1]),
        n_segments_found=len(segments),
    )


__all__ = [
    "FeatureDetectConfig",
    "FeatureDetectResult",
    "detect_features",
]
