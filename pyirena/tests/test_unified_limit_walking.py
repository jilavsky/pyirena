"""Tests for Unified Fit limit walking (fit_with_limit_walking).

The Unified panel keeps per-parameter fit limits centred on the current value
(B: value×0.2 … ×5, etc. — `panel_auto_limits`). When the optimum lies far
outside that window — e.g. changing P by ~1 moves the matching B by orders of
magnitude — a single bounded fit ends with B pinned at a limit and users had
to press Fit repeatedly. `fit_with_limit_walking` automates that loop.

Reproduces the user report on temp/Example2.h5 (P changed 3.8 → 3 and B set
to 1.41e-4: three Fit presses were needed to reach B ≈ 1.4e-6, P ≈ 3.86).
"""

import numpy as np

from pyirena.core.unified import UnifiedFitModel, panel_auto_limits


def _power_law_model(B_start, P_start, *, panel_limits=True):
    m = UnifiedFitModel(1)
    L = m.levels[0]
    L.G = 0.0
    L.fit_G = False
    L.Rg = 1e10
    L.fit_Rg = False
    L.B = B_start
    L.fit_B = True
    L.P = P_start
    L.fit_P = True
    L.fit_ETA = False
    L.fit_PACK = False
    if panel_limits:
        L.B_limits = panel_auto_limits('B', B_start)
        L.P_limits = panel_auto_limits('P', P_start)
    m.background = 0.2
    m.fit_background = True
    return m


def _synthetic_data(B_true=1.4e-6, P_true=3.86, background=0.24):
    q = np.logspace(-4, -0.7, 120)
    I = B_true * q ** (-P_true) + background
    dI = 0.01 * I
    return q, I, dI


# ── panel_auto_limits policy ─────────────────────────────────────────────────

def test_panel_auto_limits_policy():
    assert panel_auto_limits('B', 1.41e-4) == (1.41e-4 * 0.2, 1.41e-4 * 5)
    assert panel_auto_limits('P', 3.0) == (1.5, 4.5)
    assert panel_auto_limits('P', 4.0) == (2.0, 5.0)      # hi clamped to 5
    lo, hi = panel_auto_limits('Rg', 0.3)
    assert lo == 0.1                                       # lo clamped to 0.1
    lo, hi = panel_auto_limits('PACK', 5.0)
    assert hi == 12.0                                      # hi clamped to 12


# ── pinned detection ─────────────────────────────────────────────────────────

def test_pinned_fitted_parameters_detects_bounds():
    m = _power_law_model(1e-4, 3.0)
    L = m.levels[0]
    L.B = L.B_limits[0]                      # sitting exactly at the low limit
    assert (0, 'B') in m.pinned_fitted_parameters()
    L.B = 5e-5                               # interior again
    assert m.pinned_fitted_parameters() == []


def test_pinned_ignores_fixed_parameters():
    m = _power_law_model(1e-4, 3.0)
    L = m.levels[0]
    L.fit_B = False
    L.B = L.B_limits[0]
    assert m.pinned_fitted_parameters() == []


# ── The walking fit itself ───────────────────────────────────────────────────

def test_plain_fit_ends_pinned_when_optimum_is_outside_window():
    q, I, dI = _synthetic_data()
    m = _power_law_model(1.41e-4, 3.0)       # true B is 100× below the window
    m.fit(q, I, dI)
    assert (0, 'B') in m.pinned_fitted_parameters()


def test_walking_fit_reaches_distant_optimum_in_one_call():
    B_true, P_true, bkg_true = 1.4e-6, 3.86, 0.24
    q, I, dI = _synthetic_data(B_true, P_true, bkg_true)
    m = _power_law_model(1.41e-4, 3.0)
    m.fit_with_limit_walking(q, I, dI)
    L = m.levels[0]
    assert m.pinned_fitted_parameters() == []
    assert abs(L.P - P_true) < 0.05
    assert abs(L.B - B_true) / B_true < 0.25
    assert abs(m.background - bkg_true) < 0.05


def test_walking_fit_no_op_when_nothing_pinned():
    """A well-posed start must behave exactly like a plain fit."""
    B_true, P_true = 1.4e-6, 3.86
    q, I, dI = _synthetic_data(B_true, P_true)
    m = _power_law_model(1.5e-6, 3.8)        # optimum inside the window
    m.fit_with_limit_walking(q, I, dI)
    L = m.levels[0]
    assert abs(L.P - P_true) < 0.05
    assert abs(L.B - B_true) / B_true < 0.25
