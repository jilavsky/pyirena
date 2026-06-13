"""Unit tests for pyirena.core.feature_detect."""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from pyirena.core.feature_detect import (
    FeatureDetectConfig,
    FeatureDetectResult,
    detect_features,
)


# ---------------------------------------------------------------------------
# Synthetic data tests — deterministic and independent of test data files
# ---------------------------------------------------------------------------

def _guinier(q: np.ndarray, G: float, Rg: float) -> np.ndarray:
    return G * np.exp(-(q ** 2) * (Rg ** 2) / 3.0)


def _power_law(q: np.ndarray, B: float, P: float) -> np.ndarray:
    return B * np.power(q, -P)


def test_empty_input_returns_empty_result():
    """No data should produce an empty result without raising."""
    r = detect_features(np.array([]), np.array([]))
    assert isinstance(r, FeatureDetectResult)
    assert r.n_points == 0
    assert r.plateaus == []
    assert r.peaks == []
    assert r.power_law_regions == []


def test_too_few_points_returns_empty():
    """Fewer than 5 points returns empty result."""
    q = np.array([0.01, 0.02, 0.03])
    I = np.array([1.0, 0.8, 0.6])
    r = detect_features(q, I)
    assert r.n_points == 3
    assert r.plateaus == []


def test_pure_power_law_detected():
    """Single Q⁻⁴ power-law should produce a power-law region and no plateau."""
    q = np.logspace(-2, 0, 200)  # 2 decades
    I = _power_law(q, B=1e-4, P=4.0)
    r = detect_features(q, I, config=FeatureDetectConfig.saxs_preset())
    assert len(r.power_law_regions) >= 1
    pl = r.power_law_regions[0]
    assert -4.5 < pl["slope"] < -3.5, f"Expected slope near -4, got {pl['slope']}"
    assert pl["width_decades"] > 1.0


def _beaucage_level(q, G, Rg, B, P):
    """Beaucage unified-fit level: smooth Guinier→Porod transition.

        I(q) = G·exp(-q²Rg²/3) + B·[erf(q·Rg/√6)]^(3P) / q^P

    The error-function factor switches the Porod term on only above the
    Guinier knee — without it the Porod tail diverges as q→0.
    """
    from scipy.special import erf
    guinier = G * np.exp(-(q ** 2) * (Rg ** 2) / 3.0)
    arg = q * Rg / np.sqrt(6.0)
    porod = B * (erf(arg)) ** (3 * P) / np.power(q, P)
    return guinier + porod


def test_single_guinier_with_porod_detected():
    """Single Beaucage level (Guinier knee + cutoff Porod tail) — should
    produce ≥1 plateau (the Guinier region) and ≥1 power-law (the Porod tail).
    """
    q = np.logspace(-3, 0, 300)  # 3 decades, USAXS-like
    I = _beaucage_level(q, G=1.0e6, Rg=100.0, B=1.0e-3, P=4.0)
    r = detect_features(q, I)  # auto preset → USAXS
    assert r.preset_used == "usaxs"
    assert len(r.plateaus) >= 1, f"Expected ≥1 plateau, got: {r.plateaus}"
    assert len(r.power_law_regions) >= 1, (
        f"Expected ≥1 power-law region, got: {r.power_law_regions}"
    )


def test_two_level_guinier():
    """Two well-separated Beaucage levels — large aggregate (mass-fractal
    P=2.5) + small primary particles (P=4) — should give ≥2 plateaus.

    A mass-fractal aggregate has P<3 so its Porod tail decays slowly enough
    that the smaller primary-particle plateau remains visible at intermediate Q.
    This is the canonical pyirena two-level test case.
    """
    q = np.logspace(-4, 0, 500)
    I = (
        _beaucage_level(q, G=1.0e8, Rg=500.0, B=2.0e-2, P=2.0)
        + _beaucage_level(q, G=1.0e5, Rg=10.0,  B=1.0e-3, P=4.0)
    )
    r = detect_features(q, I)
    assert len(r.plateaus) >= 2, (
        f"Expected ≥2 plateaus for two-level curve, got {len(r.plateaus)}: "
        f"{r.plateaus}"
    )
    assert r.recommended_nlevels >= 2


def test_auto_preset_selection_saxs():
    """≤2 log decades → SAXS preset."""
    q = np.logspace(-2, 0, 100)  # 2 decades
    cfg = FeatureDetectConfig.auto(q)
    # SAXS preset has plateau_slope_max=0.3, USAXS has 0.5
    assert cfg.plateau_slope_max == 0.3


def test_auto_preset_selection_usaxs():
    """>2.5 log decades → USAXS preset."""
    q = np.logspace(-4, 0, 200)  # 4 decades
    cfg = FeatureDetectConfig.auto(q)
    assert cfg.plateau_slope_max == 0.5
    assert cfg.power_law_max_slope == -1.0


def test_noisy_points_dropped_via_sigma_filter():
    """Points with σ_I/I > threshold should be filtered out."""
    q = np.logspace(-2, 0, 100)
    I = _power_law(q, 1e-4, 4.0)
    sigma_I = np.full_like(I, 0.01) * I  # 1% — well below threshold
    sigma_I[50:60] = I[50:60] * 1.0       # 100% — dropped
    r = detect_features(q, I, sigma_I=sigma_I)
    # 10 noisy points dropped → 90 points used
    assert r.n_points == 90


def test_negative_intensities_handled():
    """Negative I values (e.g. from subtraction) should not crash."""
    q = np.logspace(-2, 0, 100)
    I = _power_law(q, 1e-4, 4.0)
    I[20:25] = -1e-6  # spurious negatives
    r = detect_features(q, I)
    # 5 dropped → 95 points
    assert r.n_points == 95
    assert len(r.power_law_regions) >= 1


def test_result_serializable_to_dict():
    """Result must be JSON-friendly for MCP layer."""
    q = np.logspace(-2, 0, 100)
    I = _power_law(q, 1e-4, 4.0)
    r = detect_features(q, I)
    d = r.to_dict()
    assert isinstance(d, dict)
    assert "plateaus" in d
    assert "recommended_nlevels" in d
    # All values must be plain Python types (not numpy)
    import json
    json.dumps(d)  # raises if any np scalars remain


# ---------------------------------------------------------------------------
# Optional: real test-data spot check, skipped if files absent
# ---------------------------------------------------------------------------

_REAL_DATA = Path(__file__).resolve().parents[2] / "temp" / "test_data_out" / "SAXS" / "Sa1_sphNP_5mg_0004.h5"


def _find_qi_in_h5(h):
    """Walk h5py file/group looking for a child with Q and I datasets."""
    import h5py
    for k in h:
        o = h[k]
        if isinstance(o, h5py.Group):
            if "Q" in o and "I" in o:
                return o["Q"][:], o["I"][:]
            found = _find_qi_in_h5(o)
            if found is not None:
                return found
    return None


@pytest.mark.skipif(not _REAL_DATA.exists(),
                    reason="real test data not present in checkout")
def test_real_saxs_file_runs_without_error():
    """Real SAXS file should produce a valid result (any feature count).

    The exact feature counts depend on the sample's actual structure —
    not all SAXS files contain a clean Guinier knee.  This test just
    confirms the algorithm runs end-to-end on real data without raising.
    """
    import h5py
    with h5py.File(_REAL_DATA, "r") as f:
        result = _find_qi_in_h5(f)
    assert result is not None, "Could not find Q,I in test file"
    q, I = result
    r = detect_features(q, I)
    assert r.n_points > 0
    assert r.log_decades > 0
    # Result must be JSON-serialisable
    import json
    json.dumps(r.to_dict())
