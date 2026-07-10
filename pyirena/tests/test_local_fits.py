"""Regression tests for the unified local Guinier / power-law fitters.

These pin the numeric output of :func:`pyirena.core.unified.fit_local_guinier`
and :func:`pyirena.core.unified.fit_local_power_law`, which are the single
shared implementation behind both the GUI "btwn cursors" buttons and the
``api.control`` ``fit_local_guinier`` / ``fit_local_power_law`` tools.

The expected values were captured from the original inline implementations
before they were unified, so any drift here signals a behaviour change.
"""
from __future__ import annotations

import numpy as np
import pytest

from pyirena.core.unified import fit_local_guinier, fit_local_power_law


def _synthetic():
    """Guinier knee + power-law tail + flat background."""
    q = np.logspace(-3, 0, 300)
    intensity = 1e-2 * np.exp(-q**2 * 200**2 / 3) + 5e-6 * q**(-4) + 1e-4
    err = 0.05 * intensity
    return q, intensity, err


def test_guinier_both_unweighted():
    q, I, _ = _synthetic()
    m = (q >= 2e-3) & (q <= 1.5e-2)
    r = fit_local_guinier(q[m], I[m], error=None)
    # pinned from the original GUI inline curve_fit (unweighted)
    assert r["G"] == pytest.approx(949614.7485, rel=1e-6)
    assert r["Rg"] == pytest.approx(959.8475972, rel=1e-6)
    assert r["n_points"] == int(np.sum(m))
    assert r["model_I"].shape == q[m].shape


def test_guinier_both_weighted():
    q, I, err = _synthetic()
    m = (q >= 2e-3) & (q <= 1.5e-2)
    r = fit_local_guinier(q[m], I[m], error=err[m])
    # pinned from the original api.control weighted curve_fit
    assert r["G"] == pytest.approx(16992.33361, rel=1e-6)
    assert r["Rg"] == pytest.approx(328.2490252, rel=1e-6)


def test_guinier_fix_G():
    q, I, _ = _synthetic()
    m = (q >= 2e-3) & (q <= 1.5e-2)
    g0 = (I[m][0] + I[m][-1]) / 2.0
    r = fit_local_guinier(q[m], I[m], error=None, fit_G=False, G=g0)
    assert r["G"] == pytest.approx(g0)
    assert r["Rg"] == pytest.approx(452.8012963, rel=1e-6)


def test_guinier_fix_Rg():
    q, I, _ = _synthetic()
    m = (q >= 2e-3) & (q <= 1.5e-2)
    r = fit_local_guinier(q[m], I[m], error=None, fit_Rg=False, Rg=300.0)
    assert r["Rg"] == pytest.approx(300.0)
    assert r["G"] > 0


def test_guinier_errors():
    q, I, _ = _synthetic()
    with pytest.raises(ValueError):
        fit_local_guinier(q[:2], I[:2])
    with pytest.raises(ValueError):
        fit_local_guinier(q, I, fit_G=False, fit_Rg=False, G=1.0, Rg=1.0)


def test_power_law_both_unweighted():
    q, I, _ = _synthetic()
    m = (q >= 0.1) & (q <= 1.0)
    r = fit_local_power_law(q[m], I[m], error=None)
    # pinned from the original GUI inline power-law curve_fit
    assert r["B"] == pytest.approx(5.387858497e-06, rel=1e-6)
    assert r["P"] == pytest.approx(3.967488127, rel=1e-6)


def test_power_law_both_weighted():
    q, I, err = _synthetic()
    m = (q >= 0.1) & (q <= 1.0)
    r = fit_local_power_law(q[m], I[m], error=err[m])
    # pinned from the original api.control weighted power-law curve_fit
    assert r["B"] == pytest.approx(2.185025847e-05, rel=1e-6)
    assert r["P"] == pytest.approx(3.108311571, rel=1e-6)


def test_power_law_drops_nonpositive():
    q, I, _ = _synthetic()
    m = (q >= 0.1) & (q <= 1.0)
    q_win, I_win = q[m].copy(), I[m].copy()
    # inject a couple of non-positive points that must be filtered out
    I_win[1] = -1.0
    I_win[5] = 0.0
    r = fit_local_power_law(q_win, I_win, error=None)
    assert r["n_points"] == int(np.sum(m)) - 2
    assert np.all(r["q"] > 0)


def test_power_law_errors():
    q, I, _ = _synthetic()
    with pytest.raises(ValueError):
        fit_local_power_law(q[:2], I[:2])
    with pytest.raises(ValueError):
        fit_local_power_law(q, I, fit_B=False, fit_P=False, B=1.0, P=4.0)


# --- api.control delegation: the tool must return the same numbers as core ---

def test_api_local_guinier_delegates():
    from pyirena.api.control.session import create_session, drop_session
    from pyirena.api.control.unified_fit import fit_local_guinier as api_guinier

    q, I, err = _synthetic()
    s = create_session(file_path="synthetic", q=q, intensity=I, error=err)
    try:
        m = (q >= 2e-3) & (q <= 1.5e-2)
        core = fit_local_guinier(q[m], I[m], error=err[m])
        api = api_guinier(s.session_id, q_min=2e-3, q_max=1.5e-2)
        assert api["ok"] is True
        assert api["G"] == pytest.approx(core["G"], rel=1e-9)
        assert api["Rg"] == pytest.approx(core["Rg"], rel=1e-9)
        assert api["chi_squared"] == pytest.approx(core["chi_squared"], rel=1e-9)
    finally:
        drop_session(s.session_id)


def test_api_local_power_law_delegates():
    from pyirena.api.control.session import create_session, drop_session
    from pyirena.api.control.unified_fit import fit_local_power_law as api_power

    q, I, err = _synthetic()
    s = create_session(file_path="synthetic", q=q, intensity=I, error=err)
    try:
        m = (q >= 0.1) & (q <= 1.0)
        core = fit_local_power_law(q[m], I[m], error=err[m])
        api = api_power(s.session_id, q_min=0.1, q_max=1.0)
        assert api["ok"] is True
        assert api["B"] == pytest.approx(core["B"], rel=1e-9)
        assert api["P"] == pytest.approx(core["P"], rel=1e-9)
        assert api["chi_squared"] == pytest.approx(core["chi_squared"], rel=1e-9)
    finally:
        drop_session(s.session_id)
