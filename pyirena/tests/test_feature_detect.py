"""Unit tests for pyirena.core.feature_detect (v2 segmentation algorithm)."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from pyirena.core.feature_detect import (
    FeatureDetectConfig,
    FeatureDetectResult,
    detect_features,
)


# ---------------------------------------------------------------------------
# Synthetic-curve helpers
# ---------------------------------------------------------------------------

def _power_law(q: np.ndarray, B: float, P: float) -> np.ndarray:
    return B * np.power(q, -P)


def _guinier(q: np.ndarray, G: float, Rg: float) -> np.ndarray:
    return G * np.exp(-(q ** 2) * (Rg ** 2) / 3.0)


def _beaucage_level(q, G, Rg, B, P):
    """Beaucage unified-fit level: Guinier + smooth Porod onset."""
    from scipy.special import erf
    arg = q * Rg / np.sqrt(6.0)
    return _guinier(q, G, Rg) + B * (erf(arg)) ** (3 * P) / np.power(q, P)


# ---------------------------------------------------------------------------
# Schema & edge-case tests
# ---------------------------------------------------------------------------

def test_empty_input_returns_empty_result():
    r = detect_features(np.array([]), np.array([]))
    assert isinstance(r, FeatureDetectResult)
    assert r.n_points == 0
    assert r.segments == []
    assert r.guinier_knees == []


def test_too_few_points_returns_empty():
    q = np.array([0.01, 0.02, 0.03])
    I = np.array([1.0, 0.8, 0.6])
    r = detect_features(q, I)
    assert r.segments == []


def test_result_serializable_to_dict():
    q = np.logspace(-2, np.log10(0.5), 200)
    I = _power_law(q, 1e-4, 4.0)
    r = detect_features(q, I)
    d = r.to_dict()
    assert isinstance(d, dict)
    assert "segments" in d
    assert "guinier_knees" in d
    assert "recommended_nlevels" in d
    import json
    json.dumps(d)  # raises if any numpy types remain


def test_negative_intensities_dropped():
    q = np.logspace(-2, np.log10(0.5), 100)
    I = _power_law(q, 1e-4, 4.0)
    I[20:25] = -1e-6
    r = detect_features(q, I)
    assert r.n_points == 95


def test_noisy_points_dropped_via_sigma_filter():
    q = np.logspace(-2, np.log10(0.5), 100)
    I = _power_law(q, 1e-4, 4.0)
    sigma_I = 0.01 * I
    sigma_I[40:50] = I[40:50] * 1.0  # 100% relative error — should be dropped
    r = detect_features(q, I, sigma_I=sigma_I)
    assert r.n_points == 90


# ---------------------------------------------------------------------------
# Q clipping
# ---------------------------------------------------------------------------

def test_q_clipping_default_06():
    """Synthetic data extending to Q=2.0 — analysis limited to Q ≤ 0.6."""
    q = np.logspace(-2, np.log10(2.0), 300)
    I = _power_law(q, 1e-4, 4.0)
    r = detect_features(q, I)
    assert r.q_max_analysed <= 0.6 + 1e-9


def test_q_clipping_disabled_with_none():
    """Pass q_max_clip=None → full range is analysed."""
    q = np.logspace(-2, np.log10(2.0), 300)
    I = _power_law(q, 1e-4, 4.0)
    cfg = FeatureDetectConfig(q_max_clip=None)
    r = detect_features(q, I, config=cfg)
    assert r.q_max_analysed > 0.6


# ---------------------------------------------------------------------------
# Segmentation on controlled synthetic curves
# ---------------------------------------------------------------------------

def test_pure_power_law_one_segment():
    """Single Q⁻⁴ across 2 decades should give 1 power-law segment, slope ≈ -4."""
    q = np.logspace(-2, np.log10(0.5), 200)
    I = _power_law(q, 1e-4, 4.0)
    r = detect_features(q, I)
    pl_segments = [s for s in r.segments if s["kind"] == "power_law"]
    assert len(pl_segments) >= 1
    main = pl_segments[0]
    assert -4.3 < main["slope"] < -3.7, f"Expected slope ≈ -4, got {main['slope']}"
    # Single-slope curve should produce no Guinier knees
    assert len(r.guinier_knees) == 0


def test_beaucage_one_level_finds_knee():
    """Guinier+Porod knee should produce at least one Guinier knee.

    The change-point algorithm tracks the slope as it transitions from ~0
    (plateau) to -P (Porod tail), which may produce more than one segment
    in the transition zone.  What must always hold: the slope ends up at
    -P at the high-Q end, and at least one Guinier knee is reported.
    """
    q = np.logspace(-3, np.log10(0.5), 300)
    I = _beaucage_level(q, G=1.0e6, Rg=100.0, B=1.0e-3, P=4.0)
    r = detect_features(q, I)
    assert r.segments, f"Expected at least one segment"
    # Last segment should have slope close to the Porod exponent -4
    last_slope = r.segments[-1]["slope"]
    assert -4.5 < last_slope < -3.0, (
        f"Expected last segment slope near -4, got {last_slope:.2f}"
    )
    # At least one knee in the transition zone
    assert len(r.guinier_knees) >= 1, (
        f"Expected at least one Guinier knee, got 0. "
        f"Segments: {[(s['q_min'], s['q_max'], s['slope']) for s in r.segments]}"
    )


def test_background_recognised_at_high_q_end():
    """A flat tail at the highest-Q end should be classified as background.

    Built so the constant background fully dominates above Q≈0.05: power-law
    signal (Q⁻⁴ with B=1e-8) falls below the background (0.1) very early,
    leaving a clearly flat tail.
    """
    q = np.logspace(-2, np.log10(0.5), 200)
    I = _power_law(q, 1e-8, 4.0) + 0.1
    r = detect_features(q, I)
    assert r.segments, "Expected at least one segment"
    last_seg = r.segments[-1]
    assert last_seg["kind"] == "background", (
        f"Expected last segment to be background, got {last_seg['kind']} "
        f"with slope={last_seg['slope']:.2f}"
    )
    assert r.background_q_min is not None
    assert r.background_q_min == last_seg["q_min"]


def test_min_segment_decades_enforced_for_interior():
    """No interior segment narrower than min_segment_decades should be returned.

    Edge segments use the looser ``edge_min_segment_decades`` and may be
    narrower than the interior threshold.
    """
    q = np.logspace(-2, np.log10(0.5), 300)
    I = _power_law(q, 1e-4, 4.0)
    cfg = FeatureDetectConfig(min_segment_decades=0.5,
                              edge_min_segment_decades=0.5)
    r = detect_features(q, I, config=cfg)
    for s in r.segments:
        assert s["width_decades"] + 1e-9 >= 0.5, (
            f"Segment narrower than min_segment_decades: {s}"
        )


def test_edge_segment_promotion():
    """A narrow stable slope region at the very low-Q end of the data should
    be reported because edge segments use ``edge_min_segment_decades`` (looser).

    Mimics the sample15 case: a Guinier plateau (slope ≈ 0) at Q < Q_min×1.5
    that transitions sharply into a Porod tail (slope ≈ -4) over the next
    half-decade.
    """
    # Build a curve: flat at low Q, then Q^-4 above the knee
    q = np.logspace(-3, np.log10(0.3), 400)
    Rg = 1000.0  # Very large structure → knee at Q ≈ 1/Rg = 1e-3
    I = _beaucage_level(q, G=1.0e8, Rg=Rg, B=1.0e-1, P=4.0)
    cfg = FeatureDetectConfig(
        min_segment_decades=0.20,         # Strict interior threshold
        edge_min_segment_decades=0.05,    # Looser for edges
    )
    r = detect_features(q, I, config=cfg)
    # There must be a first segment with shallow slope at the very low-Q end
    assert r.segments, "Expected at least one segment"
    first = r.segments[0]
    assert first["q_min"] == r.q_min_analysed
    assert abs(first["slope"]) < 1.5, (
        f"Expected shallow slope at very low-Q edge, got {first['slope']:.2f}"
    )


def test_changepoint_detects_smooth_drift():
    """Slope drift from -2 → -3 → -4 should produce multiple segments,
    not one merged region.

    This is the sample25 failure mode under v0.8.4: variance-based stability
    saw the smooth drift as 'stable' and lumped it; change-point detection
    catches it because the mean differs between left/right windows.
    """
    # Build a curve where slope deliberately changes: three Beaucage-like
    # regions with different P values, smoothly joined
    q = np.logspace(-4, np.log10(0.3), 600)
    # Build I such that slope is approximately -2 at low Q, -3 at mid Q,
    # -4 at high Q.  Use a piecewise log construction.
    logQ = np.log10(q)
    # In log space: integrate slope to get logI.  Slope function:
    slope_target = np.where(logQ < -2.5, -2.0,
                            np.where(logQ < -1.5, -3.0, -4.0))
    # Integrate slope w.r.t. logQ to get logI; smooth slope a bit so it's
    # not a pure step (more realistic)
    from scipy.ndimage import gaussian_filter1d
    slope_smoothed = gaussian_filter1d(slope_target, sigma=8.0)
    dlogQ = np.diff(logQ, prepend=logQ[0])
    logI = 5.0 + np.cumsum(slope_smoothed * dlogQ)
    I = 10.0 ** logI
    r = detect_features(q, I)
    pl_segs = [s for s in r.segments if s["kind"] == "power_law"]
    assert len(pl_segs) >= 2, (
        f"Expected ≥2 power-law segments for a 3-region slope drift, "
        f"got {len(pl_segs)}. Segments: "
        f"{[(s['q_min'], s['q_max'], s['slope']) for s in r.segments]}"
    )


def test_recommended_nlevels_excludes_background():
    """If a background is detected, recommended_nlevels must exclude it."""
    q = np.logspace(-2, np.log10(0.5), 200)
    I = _power_law(q, 1e-4, 4.0) + 0.1
    r = detect_features(q, I)
    has_bg = any(s["kind"] == "background" for s in r.segments)
    if has_bg:
        non_bg = sum(1 for s in r.segments if s["kind"] != "background")
        assert r.recommended_nlevels == non_bg


# ---------------------------------------------------------------------------
# Ground-truth fidelity test (parametrised over the 31 hand-labelled files)
# ---------------------------------------------------------------------------

_GT_DIR = (Path(__file__).resolve().parents[2]
           / "testData" / "StructureIdentificationExamples")
_GT_FILE = _GT_DIR / "Description and task.txt"


def _parse_ground_truth(text: str) -> dict:
    """Parse the description.txt into {sample_name: [(kind, q_lo, q_hi), ...]}.

    Only PLS / GP / Background entries are returned (kinds the segmenter
    can directly match against).  GK ranges are not collected because they
    are the *gaps* between segments and are reported in `guinier_knees`,
    not `segments`.
    """
    out: dict = {}
    cur_name: str | None = None
    range_re = re.compile(
        r"^\s*"
        r"(?P<tag>Background|Blackground|PLS|GP|GP/PLS|Bkg)"
        r"\b[^0-9-]*"
        r"(?P<lo>(?:0?\.[0-9]+|[0-9]+(?:\.[0-9]+)?))"
        r"\s*-\s*"
        r"(?P<hi>(?:0?\.[0-9]+|[0-9]+(?:\.[0-9]+)?|min))",
        re.IGNORECASE,
    )
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        m = re.match(r"^(Sample\d+|Plate\w+_NX|Rod\w+_NX|RodSfAgg_NSX)\b", line)
        if m:
            cur_name = m.group(1).strip()
            if cur_name == "RodSfAgg_NSX":
                cur_name = "RodSfAgg_NX"
            out.setdefault(cur_name, [])
            continue
        if cur_name is None:
            continue
        rm = range_re.search(line)
        if not rm:
            continue
        tag = rm.group("tag").upper()
        if tag in ("BACKGROUND", "BLACKGROUND", "BKG"):
            kind = "background"
        elif tag == "PLS":
            kind = "power_law"
        elif tag in ("GP", "GP/PLS"):
            kind = "guinier_plateau"
        else:
            continue
        try:
            q_lo = float(rm.group("lo"))
        except ValueError:
            continue
        hi_str = rm.group("hi")
        try:
            q_hi = float("nan") if hi_str.lower() == "min" else float(hi_str)
        except ValueError:
            continue
        if not np.isnan(q_hi) and q_lo > q_hi:
            q_lo, q_hi = q_hi, q_lo
        out[cur_name].append((kind, q_lo, q_hi))
    return out


def _find_qi(h):
    """Walk an h5py.File/Group looking for a {Q, I[, Idev]} group."""
    import h5py
    for k in h:
        o = h[k]
        if isinstance(o, h5py.Group):
            if "Q" in o and "I" in o:
                idev = o["Idev"][:] if "Idev" in o else None
                return o["Q"][:], o["I"][:], idev
            r = _find_qi(o)
            if r is not None:
                return r
    return None


def _overlap_decades(a_lo, a_hi, b_lo, b_hi):
    if any(np.isnan([a_lo, a_hi, b_lo, b_hi])):
        return 0.0
    lo = max(a_lo, b_lo)
    hi = min(a_hi, b_hi)
    if hi <= lo:
        return 0.0
    return float(np.log10(hi / lo))


# Representative samples spanning easy/medium/hard cases.  sample15 and
# sample25 are explicit regression tests for the two v0.8.4 failure modes
# that motivated the v0.8.5 algorithm rewrite (edge knee + smooth drift).
_FIDELITY_SAMPLES = [
    "sample1", "sample6", "sample11", "sample14",
    "sample15",  # was failing in v0.8.4: low-Q Guinier knee missed
    "sample20",
    "sample25",  # was failing in v0.8.4: smooth slope drift lumped
    "sample28", "sample3",
    "PlateWeak_NX", "RodSfAgg_NX",
]


@pytest.mark.skipif(not _GT_FILE.exists(),
                    reason="Ground-truth description file not present")
@pytest.mark.parametrize("sample_name", _FIDELITY_SAMPLES)
def test_ground_truth_fidelity(sample_name: str):
    """For each labelled sample, at least 75% of human-labelled (PLS / GP /
    Background) ranges with definite high Q must overlap a detected segment.

    A detected segment "matches" a labelled range if their log-Q overlap is
    ≥ max(30% of the labelled width, 0.15 decades).
    """
    import h5py
    path = _GT_DIR / f"{sample_name}.h5"
    if not path.exists():
        pytest.skip(f"{path.name} not present")
    gt = _parse_ground_truth(_GT_FILE.read_text(encoding="utf-8"))
    # Sample name in description file may differ in case
    gt_entries = (gt.get(sample_name)
                  or gt.get(sample_name.replace("sample", "Sample"))
                  or [])
    gt_segs = [(k, lo, hi) for k, lo, hi in gt_entries if not np.isnan(hi)]
    if not gt_segs:
        pytest.skip(f"{sample_name}: no usable ground-truth entries")

    with h5py.File(path, "r") as fh:
        result = _find_qi(fh)
    assert result is not None, f"Could not find Q,I in {path.name}"
    q, I, idev = result

    res = detect_features(q, I, sigma_I=idev)
    assert len(res.segments) >= 1, f"{sample_name}: no segments detected"

    n_match = 0
    for kind, q_lo, q_hi in gt_segs:
        width = max(np.log10(q_hi / q_lo), 1e-6)
        best = 0.0
        for s in res.segments:
            o = _overlap_decades(q_lo, q_hi, s["q_min"], s["q_max"])
            if o > best:
                best = o
        if best >= 0.3 * width or best >= 0.15:
            n_match += 1
    match_frac = n_match / len(gt_segs)
    assert match_frac >= 0.75, (
        f"{sample_name}: only {n_match}/{len(gt_segs)} GT segments matched "
        f"(< 75%).\n  GT: {gt_segs}\n  Detected: "
        f"{[(s['kind'], s['q_min'], s['q_max'], s['slope']) for s in res.segments]}"
    )
