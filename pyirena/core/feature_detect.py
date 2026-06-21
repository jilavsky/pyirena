"""Feature detection for small-angle scattering I(Q) curves.

Segments an I(Q) log-log curve into power-law-slope regions using a
two-pass change-point detector on the smoothed slope profile, and derives
Guinier knees from adjacent segments whose slopes change in the expected
direction (Porod tail giving way to Guinier plateau as Q decreases).

Conceptual model
----------------
Every region of a SAS curve has a *local power-law slope*.  Steep slopes
(P ≈ 2-4) correspond to Porod tails; slope ≈ 0 corresponds to a Guinier
plateau or, at the highest-Q end, a constant background.  Therefore a
complete segmentation of the curve into piecewise-constant-slope segments
recovers all the information a Unified Fit user needs to choose levels and
Q-windows; Guinier knees fall naturally between adjacent slope segments.

Algorithm (v0.8.5)
------------------
The detector is **change-point based**: at each candidate boundary point,
compute the difference between the mean slope in a window to its left and
the mean slope in a window to its right.  Local maxima of this difference
that exceed a threshold are boundaries.  Segments are the intervals between
boundaries.

Two passes:

  Pass 1 (loose): wider window, larger threshold — find the major
                  transitions cleanly.
  Pass 2 (tight): narrower window, smaller threshold — re-scan any segment
                  wider than ``wide_region_decades`` for sub-structure.
                  Catches the case where a single broad region actually
                  contains multiple distinct power-law slopes that drift
                  smoothly between each other.

Post-processing: merge adjacent segments whose mean slopes are within
``merge_slope_tol`` (suppresses spurious splits at noise spikes); width
filter that uses a looser ``edge_min_segment_decades`` for segments touching
the data extremes (lets narrow low-Q PLSes and high-Q backgrounds through).

This replaces the variance-based stability check from v0.8.4, which could
not distinguish "slope is constant" from "slope is slowly changing".  The
change-point statistic answers a different question — "does the mean differ
between the two sides?" — which catches both sharp transitions and smooth
drifts.

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

    # --- Change-point detection: Pass 1 (coarse) ---
    change_window_1: float = 0.30
    """Pass-1 half-window width for change-point detection (decades on each
    side of the candidate boundary point)."""

    change_threshold_1: float = 0.40
    """Pass-1 minimum |mean_left − mean_right| to declare a change-point."""

    # --- Change-point detection: Pass 2 (refine wide regions) ---
    change_window_2: float = 0.20
    """Pass-2 half-window width."""

    change_threshold_2: float = 0.20
    """Pass-2 minimum |Δmean| — tighter than pass 1 to catch sub-structure
    inside an apparently wide region."""

    wide_region_decades: float = 1.0
    """A segment wider than this gets a Pass-2 re-scan to look for hidden
    sub-segments (smooth slope drifts that pass-1 missed)."""

    # --- Segment filtering ---
    min_segment_decades: float = 0.10
    """Minimum width of an interior segment (not touching either data
    extreme)."""

    edge_min_segment_decades: float = 0.05
    """Minimum width for segments touching the lowest-Q or highest-Q data
    point.  Looser than ``min_segment_decades`` so that narrow low-Q
    Guinier plateaus and high-Q backgrounds are not silently discarded."""

    merge_slope_tol: float = 0.15
    """Adjacent segments are merged into one if their mean slopes differ by
    less than this — suppresses spurious sub-divisions at noise spikes."""

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

    guinier_plateau_slope_max: float = 0.40
    """|slope| < this and not background → guinier_plateau."""

    guinier_knee_min_delta_slope: float = 0.10
    """Minimum slope jump between adjacent segments to emit a Guinier knee.

    Per the physical direction check in `_derive_knees`, a knee is only
    reported when the lower-Q segment has a *shallower* slope than the
    higher-Q segment (Porod tail giving way to a Guinier plateau).
    """

    # --- Noise filtering ---
    noisy_point_threshold: float = 0.5
    """Drop points with sigma_I/I > this before smoothing."""


# ---------------------------------------------------------------------------
# Result
# ---------------------------------------------------------------------------

@dataclass
class FeatureDetectResult:
    segments: list[dict] = field(default_factory=list)
    """One entry per detected segment, sorted by ascending Q.

    Each: ``q_min``, ``q_max``, ``P``, ``P_std``, ``kind``,
    ``intensity_mid``, ``width_decades``.

    ``P`` is the positive Porod exponent (I ∝ Q^-P); ``kind`` is one of
    ``"background"``, ``"guinier_plateau"``, ``"power_law"``.
    """

    guinier_knees: list[dict] = field(default_factory=list)
    """Derived list of inferred Guinier knees between adjacent segments.

    Each: ``q_min``, ``q_max``, ``q_center``, ``P_low_q``, ``P_high_q``,
    ``delta_P``.

    ``P_low_q`` and ``P_high_q`` are positive Porod exponents (I ∝ Q^-P).
    Only transitions where P_low_q < P_high_q are listed — i.e. the
    exponent is smaller (shallower slope) at low Q, the physical Guinier-knee
    condition.
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
# Smoothing & slope (unchanged from v0.8.4)
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
# Change-point detection
# ---------------------------------------------------------------------------

def _changepoint_z(logQ: np.ndarray, slope: np.ndarray,
                   half_window_decades: float,
                   lo: int = 0, hi: Optional[int] = None) -> np.ndarray:
    """Compute |mean_L − mean_R| at each point.

    For each interior index `i`, the left window covers slope[max(lo, i-W):i]
    and the right window covers slope[i:min(hi, i+W)] where W is the index
    range covered by ``half_window_decades`` on each side of ``logQ[i]``.

    Windows are clipped to ``[lo, hi)`` so this can be called on sub-ranges
    (used by the two-pass refinement to scan inside a single segment without
    leaking statistics from the rest of the curve).

    Returns an array of length ``len(logQ)`` with z=0 outside ``[lo, hi)``
    and at points where either side window has fewer than 2 samples.
    """
    n = logQ.size
    if hi is None:
        hi = n
    z = np.zeros(n)
    for i in range(max(lo + 1, 1), min(hi - 1, n - 1)):
        lo_L = int(np.searchsorted(logQ, logQ[i] - half_window_decades,
                                   side="left"))
        lo_L = max(lo_L, lo)
        hi_L = i
        lo_R = i
        hi_R = int(np.searchsorted(logQ, logQ[i] + half_window_decades,
                                   side="right"))
        hi_R = min(hi_R, hi)
        if hi_L - lo_L < 2 or hi_R - lo_R < 2:
            continue
        mean_L = float(np.mean(slope[lo_L:hi_L]))
        mean_R = float(np.mean(slope[lo_R:hi_R]))
        z[i] = abs(mean_L - mean_R)
    return z


def _find_changepoints(z: np.ndarray, logQ: np.ndarray,
                       threshold: float, min_separation_decades: float,
                       lo: int = 0, hi: Optional[int] = None) -> list[int]:
    """Return indices of local maxima of `z` in ``[lo, hi)`` exceeding
    `threshold`, with no two selections closer than `min_separation_decades`.

    Resolves overlaps by keeping the higher z; output is sorted ascending.
    """
    n = z.size
    if hi is None:
        hi = n
    # Collect candidate local maxima
    candidates: list[tuple[int, float]] = []
    for i in range(max(lo + 1, 1), min(hi - 1, n - 1)):
        if z[i] < threshold:
            continue
        if z[i] > z[i - 1] and z[i] >= z[i + 1]:
            candidates.append((i, z[i]))
    if not candidates:
        return []
    # Sort by z descending, then greedily accept respecting separation
    candidates.sort(key=lambda x: -x[1])
    selected: list[int] = []
    for idx, _ in candidates:
        too_close = False
        for s in selected:
            if abs(logQ[idx] - logQ[s]) < min_separation_decades:
                too_close = True
                break
        if not too_close:
            selected.append(idx)
    return sorted(selected)


def _two_pass_segment(logQ: np.ndarray, slope: np.ndarray,
                      cfg: FeatureDetectConfig) -> list[tuple[int, int]]:
    """Run two-pass change-point segmentation, return [(start, end_exclusive)]."""
    n = logQ.size
    if n < 4:
        return [(0, n)] if n > 0 else []

    # Pass 1: coarse boundaries on the full curve
    z1 = _changepoint_z(logQ, slope, cfg.change_window_1 / 2.0)
    cps1 = _find_changepoints(
        z1, logQ,
        threshold=cfg.change_threshold_1,
        min_separation_decades=cfg.min_segment_decades,
    )

    coarse_boundaries = [0] + cps1 + [n]

    # Pass 2: re-scan any segment wider than wide_region_decades
    refined_boundaries: list[int] = [0]
    for k in range(len(coarse_boundaries) - 1):
        s, e = coarse_boundaries[k], coarse_boundaries[k + 1]
        if e <= s:
            continue
        width = float(logQ[e - 1] - logQ[s]) if e > s else 0.0
        if width > cfg.wide_region_decades and e - s >= 8:
            z2 = _changepoint_z(logQ, slope, cfg.change_window_2 / 2.0,
                                lo=s, hi=e)
            sub_cps = _find_changepoints(
                z2, logQ,
                threshold=cfg.change_threshold_2,
                min_separation_decades=cfg.min_segment_decades,
                lo=s, hi=e,
            )
            for c in sub_cps:
                if c != s and c != e:
                    refined_boundaries.append(c)
        if e < n:
            refined_boundaries.append(e)
    refined_boundaries.append(n)
    refined_boundaries = sorted(set(refined_boundaries))

    seg_ranges = []
    for k in range(len(refined_boundaries) - 1):
        s, e = refined_boundaries[k], refined_boundaries[k + 1]
        if e - s < 2:
            continue
        seg_ranges.append((s, e))
    return seg_ranges


# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------

def _merge_similar_segments(seg_ranges: list[tuple[int, int]],
                            slope: np.ndarray,
                            tol: float) -> list[tuple[int, int]]:
    """Merge adjacent segments whose mean slopes differ by less than ``tol``."""
    if not seg_ranges:
        return []
    merged = [seg_ranges[0]]
    for s, e in seg_ranges[1:]:
        ps, pe = merged[-1]
        prev_mean = float(np.mean(slope[ps:pe]))
        cur_mean = float(np.mean(slope[s:e]))
        if abs(prev_mean - cur_mean) < tol:
            merged[-1] = (ps, e)
        else:
            merged.append((s, e))
    return merged


def _filter_segments(seg_ranges: list[tuple[int, int]], logQ: np.ndarray,
                     n_total: int, min_decades: float,
                     edge_min_decades: float) -> list[tuple[int, int]]:
    """Drop segments narrower than the applicable threshold.

    Edge-touching segments (those whose s==0 or e==n_total) use the looser
    ``edge_min_decades``; all others must satisfy ``min_decades``.
    """
    out = []
    for s, e in seg_ranges:
        width = float(logQ[e - 1] - logQ[s])
        is_edge = (s == 0) or (e == n_total)
        threshold = edge_min_decades if is_edge else min_decades
        if width >= threshold:
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
            "slope":         mean_slope,   # kept for internal use (negative)
            "slope_std":     float(np.std(seg_slope)),
            "P":             -mean_slope,  # positive Porod exponent (I ∝ Q^-P)
            "P_std":         float(np.std(seg_slope)),
            "kind":          kind,
            "intensity_mid": float(seg_I[mid_idx]),
            "width_decades": float(np.log10(q_max_seg) - np.log10(q_min_seg)),
        })
    return out


def _derive_knees(segments: list[dict], cfg: FeatureDetectConfig) -> list[dict]:
    """Derive Guinier knees between adjacent segments.

    A Guinier knee is physically meaningful only when, going from high-Q to
    low-Q, the power-law slope becomes *shallower* (smaller absolute value).
    Equivalently, in the segment list (sorted low-Q → high-Q), a knee between
    segments[i] (lower-Q) and segments[i+1] (higher-Q) is valid when:

        |slope_low_Q| < |slope_high_Q|

    i.e. the lower-Q segment has a *shallower* slope.  If the low-Q segment
    is *steeper* than the high-Q segment, the transition is not a Guinier
    knee — it could be a mass-fractal transition, aggregate crossover, or
    structure factor; the caller should not label it as a knee.
    """
    knees: list[dict] = []
    for i in range(len(segments) - 1):
        low_q_seg  = segments[i]
        high_q_seg = segments[i + 1]
        p_low_q  = abs(low_q_seg["slope"])   # P = -slope (positive)
        p_high_q = abs(high_q_seg["slope"])
        delta_p  = abs(p_low_q - p_high_q)
        if delta_p < cfg.guinier_knee_min_delta_slope:
            continue
        # Guinier knee: P is smaller (shallower) at low-Q side
        if p_low_q >= p_high_q:
            continue
        q_lo = float(low_q_seg["q_max"])
        q_hi = float(high_q_seg["q_min"])
        if q_hi < q_lo:
            q_hi = q_lo
        q_center = float(np.sqrt(max(q_lo * q_hi, 1e-30)))
        knees.append({
            "q_min":    q_lo,
            "q_max":    q_hi,
            "q_center": q_center,
            "P_low_q":  float(p_low_q),   # Porod exponent on the low-Q side
            "P_high_q": float(p_high_q),  # Porod exponent on the high-Q side
            "delta_P":  float(delta_p),
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

    # 5. Two-pass change-point segmentation
    seg_ranges = _two_pass_segment(logQ, slope, config)

    # 6. Merge segments with very similar slopes (suppress noise splits)
    seg_ranges = _merge_similar_segments(seg_ranges, slope,
                                         tol=config.merge_slope_tol)

    # 7. Width filter (asymmetric: edge segments use looser threshold)
    n_total = logQ.size
    seg_ranges = _filter_segments(seg_ranges, logQ, n_total,
                                  min_decades=config.min_segment_decades,
                                  edge_min_decades=config.edge_min_segment_decades)

    # 8. Build & classify segments
    segments = _build_segments(seg_ranges, q, I, slope, config)

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
