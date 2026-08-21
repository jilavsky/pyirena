"""API-level tests for Unified run_fit limit walking and pinned reporting.

run_fit(walk_limits=True) (the default, matching the GUI Fit button) must
reach an optimum far outside the current limits window; walk_limits=False
must respect the bounds and report the pinned parameter so an agent can see
which bound stopped the fit.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import pyirena.api.control as ctrl
from pyirena.io.nxcansas_unified import create_nxcansas_file

B_TRUE, P_TRUE, BKG_TRUE = 1.4e-6, 3.86, 0.24


def _make_powerlaw_file(folder: Path) -> Path:
    q = np.logspace(-4, -0.7, 150)
    intensity = B_TRUE * q ** (-P_TRUE) + BKG_TRUE
    fp = folder / "powerlaw.h5"
    create_nxcansas_file(fp, q, intensity, error=intensity * 0.01,
                         sample_name="powerlaw")
    return fp


def _setup_pinned_start(sid: str) -> None:
    """Pure power-law level starting 100× above the true B, panel-style limits.

    One level only.  select_model already supplies level 1, so calling
    add_unified_level here used to leave a second, entirely default level free
    in the fit — a duplicate power law the data cannot distinguish from the
    first.  Its G and B then drift to their 1e-10 / 1e-20 lower bounds, and
    whether they land *on* the bound or just above it is numerical luck that
    differs between platforms.  That is what made this test fail on Windows
    only, reporting level 2 as pinned.
    """
    sel = ctrl.select_model(sid, "unified_fit")
    assert "error" not in sel
    assert sel["nlevels"] == 1, "select_model no longer starts with one level"
    assert "error" not in ctrl.set_parameter_value(sid, "G_1", 0.0)
    assert "error" not in ctrl.fix_parameter(sid, "G_1")
    assert "error" not in ctrl.set_parameter_value(sid, "Rg_1", 1e10)
    assert "error" not in ctrl.fix_parameter(sid, "Rg_1")
    assert "error" not in ctrl.set_parameter_value(sid, "B_1", 1.41e-4)
    assert "error" not in ctrl.free_parameter(sid, "B_1")
    assert "error" not in ctrl.set_parameter_bounds(
        sid, "B_1", 1.41e-4 * 0.2, 1.41e-4 * 5)
    assert "error" not in ctrl.set_parameter_value(sid, "P_1", 3.0)
    assert "error" not in ctrl.free_parameter(sid, "P_1")
    assert "error" not in ctrl.set_parameter_bounds(sid, "P_1", 1.5, 4.5)


def _fitted(res, name):
    return next(row["value"] for row in res["parameters_updated"]
                if row["name"] == name)


def test_run_fit_walks_to_distant_optimum_by_default(tmp_path):
    fp = _make_powerlaw_file(tmp_path)
    sid = ctrl.open_dataset(str(fp))["session_id"]
    _setup_pinned_start(sid)

    res = ctrl.run_fit(sid)
    assert "error" not in res
    assert res["pinned_parameters"] == []
    assert res["warnings"] == []
    assert abs(_fitted(res, "P_1") - P_TRUE) < 0.05
    assert abs(_fitted(res, "B_1") - B_TRUE) / B_TRUE < 0.25
    ctrl.close_session(sid)


def test_run_fit_without_walking_reports_pinned_bound(tmp_path):
    fp = _make_powerlaw_file(tmp_path)
    sid = ctrl.open_dataset(str(fp))["session_id"]
    _setup_pinned_start(sid)

    res = ctrl.run_fit(sid, walk_limits=False)
    assert "error" not in res
    # B cannot leave its window: it must stop at the lower bound and say so.
    assert abs(_fitted(res, "B_1") - 1.41e-4 * 0.2) / (1.41e-4 * 0.2) < 1e-3
    assert "level 1 B" in res["pinned_parameters"]
    assert any("pinned" in w for w in res["warnings"])
    ctrl.close_session(sid)
