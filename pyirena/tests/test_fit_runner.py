"""Shared fit runners and data rule (issue #38, planning/fit-runner/README.md)."""
import numpy as np
import pytest

import pyirena.api.control as ctrl
from pyirena.core.fit_data import prepare_fit_arrays, q_range_mask
from pyirena.core.tool_config import apply_sizes_prefits, sizes_model_from_config

Q = np.logspace(-3, -0.5, 120)


def _curve():
    from pyirena.core.unified import UnifiedFitModel
    m = UnifiedFitModel(1)
    lv = m.levels[0]
    lv.G, lv.Rg, lv.P, lv.B = 1000.0, 50.0, 4.0, 1e-3
    m.background = 0.01
    I = m.calculate_intensity(Q) * (1 + 0.02 * np.random.default_rng(2).standard_normal(Q.size))
    return I, 0.02 * I


# ── data rule ────────────────────────────────────────────────────────────────

def test_prepare_keeps_negative_I_and_repairs_sigma():
    q = np.array([0.0, 0.01, 0.02, 0.03, 0.04, np.nan])
    I = np.array([1.0, -2.0, 3.0, 0.0, 5.0, 1.0])
    e = np.array([0.1, 0.1, 0.0, 0.0, np.inf, 0.1])
    p = prepare_fit_arrays(q, I, e, error_fraction=0.1)
    assert list(p.q) == [0.01, 0.02, 0.04]           # Q<=0, NaN Q, I=0&σ=0 dropped
    assert list(p.intensity) == [-2.0, 3.0, 5.0]     # negative I kept (D3)
    assert np.allclose(p.error, [0.1, 0.3, 0.5])     # σ<=0 / inf repaired (D2)
    assert p.n_sigma_repaired == 2 and p.n_dropped == 3
    assert prepare_fit_arrays(q, I, None).error is None


def test_q_range_mask_inclusive_and_open():
    q = np.array([1.0, 2.0, 3.0])
    assert list(q_range_mask(q, 2.0, None)) == [False, True, True]
    assert list(q_range_mask(q, None, 2.0)) == [True, True, False]


# ── Sizes ────────────────────────────────────────────────────────────────────

def test_sizes_config_reads_mc_settings_and_mcsas_alias():
    s = sizes_model_from_config({"method": "mcsas", "mcsas_convergence": 2.5,
                                 "montecarlo_max_iter": 1234})
    assert s.method == "montecarlo"
    assert s.montecarlo_convergence == 2.5 and s.montecarlo_max_iter == 1234
    assert s.maxent_max_iter == 1000 and s.tnnls_max_iter == 1000
    assert s.montecarlo_n_repetitions == 1


def test_sizes_power_law_prefit_falls_back_to_cursor_range(monkeypatch):
    s = sizes_model_from_config({})
    seen = {}
    monkeypatch.setattr(type(s), "fit_power_law",
                        lambda self, q, I, lo, hi, **k: seen.update(lo=lo, hi=hi) or {})
    apply_sizes_prefits(s, {"fit_power_law_P": True}, Q, Q, fit_q_min=0.01, fit_q_max=0.1)
    assert (seen["lo"], seen["hi"]) == (0.01, 0.1)
    apply_sizes_prefits(s, {"fit_power_law_P": True}, Q, Q)
    assert (seen["lo"], seen["hi"]) == (Q.min(), Q.max())


def test_control_sizes_session_smears_and_uses_one_repetition():
    I, e = _curve()
    sid = ctrl.open_dataset_from_data(q=Q.tolist(), intensity=I.tolist(), error=e.tolist(),
                                      is_slit_smeared=True, slit_length=0.02)["session_id"]
    try:
        ctrl.select_sizes_model(sid)
        from pyirena.api.control.session import get_session
        m = get_session(sid).model
        assert m.use_slit_smearing and m.slit_length == pytest.approx(0.02)
        assert m.montecarlo_n_repetitions == 1
    finally:
        ctrl.close_session(sid)


def test_control_sizes_fit_reports_repaired_sigma():
    I, e = _curve()
    e = e.copy()
    e[40] = 0.0
    sid = ctrl.open_dataset_from_data(q=Q.tolist(), intensity=I.tolist(),
                                      error=e.tolist())["session_id"]
    try:
        ctrl.select_sizes_model(sid)
        r = ctrl.run_sizes_fit(sid)
        assert r["success"]
    finally:
        ctrl.close_session(sid)
