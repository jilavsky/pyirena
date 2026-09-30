"""Tests for ``analyze()`` — fit a curve with a saved config, in one call.

Driven by `testData/Scripting/`: one configuration exported from the GUI per
tool, with the data file it was set up on. That matters more than it sounds.
A config I wrote from the spec would test my reading of the spec; these
files test what the GUI actually writes, which is the contract `analyze`
exists to honour.

The load-bearing test is `test_analyze_agrees_with_the_batch_path`: the same
data and the same config, fitted through `analyze` and through
`pyirena.batch`, must produce the same numbers. They share their
config-to-model translation (`core/tool_config.py`) precisely so that they
cannot drift, and this is what would notice if they did.

Skips cleanly when the fixtures are absent (an installed package rather than
a checkout).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import pyirena.api.control as ctrl

SCRIPTING = Path(__file__).resolve().parents[3] / "testData" / "Scripting"

#: tool → (data file, config file). One exported config per tool.
FIXTURES = {
    "unified_fit":  ("UnifiedFit_PP15.h5",  "UnifiedFit"),
    "sizes":        ("SizeDis_PP15.h5",     "SizeDis.json"),
    "simple_fits":  ("SimpleFits_Porod.h5", "SimpleFits_Porod.json"),
    "modeling":     ("Modeling_PP15.h5",    "modeling.json"),
    "waxs_peakfit": ("WAXS_Al_7075.hdf",    "WAXS"),
    "carbon_fit":   ("Carbon_CE_1400.h5",   "carbon_CE_1400.json"),
}

pytestmark = pytest.mark.skipif(
    not SCRIPTING.exists(), reason="testData/Scripting fixtures not present"
)


def _load(tool: str):
    """Return (data payload, config dict) for one tool's fixture."""
    from pyirena.io.hdf5 import readGenericNXcanSAS

    data_file, config_file = FIXTURES[tool]
    raw = readGenericNXcanSAS(str(SCRIPTING), data_file)
    payload = {
        "q": np.asarray(raw["Q"], dtype=float).tolist(),
        "intensity": np.asarray(raw["Intensity"], dtype=float).tolist(),
        "label": data_file,
    }
    if raw.get("Error") is not None:
        payload["error"] = np.asarray(raw["Error"], dtype=float).tolist()
    return payload, json.loads((SCRIPTING / config_file).read_text(encoding="utf-8"))


def _chi(report: dict) -> float:
    q = report["quality"]
    return q["reduced_chi_squared"] if q["reduced_chi_squared"] is not None \
        else q["chi_squared"]


@pytest.fixture(autouse=True)
def _no_sessions_left_behind():
    yield
    assert ctrl.list_open_sessions()["count"] == 0, "analyze leaked a session"


# ---------------------------------------------------------------------------
# One call, every tool
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tool", sorted(FIXTURES))
def test_a_real_exported_config_fits_in_one_call(tool):
    data, config = _load(tool)
    report = ctrl.analyze(data, config)
    assert "error" not in report, report

    assert report["ok"] is True
    assert report["tool"] == tool
    assert np.isfinite(_chi(report))
    assert report["results"]
    assert report["config"]

    info = report["analyze"]
    assert info["tool"] == tool
    assert info["config_applied"] is True
    assert info["n_points_input"] > 0


@pytest.mark.parametrize("tool", sorted(FIXTURES))
def test_the_tool_is_inferred_from_the_config(tool):
    """The GUI writes one tool section, so the caller should not repeat it."""
    data, config = _load(tool)
    assert ctrl.analyze(data, config)["tool"] == tool
    # ...and naming it explicitly agrees.
    assert ctrl.analyze(data, config, tool=tool)["tool"] == tool


@pytest.mark.parametrize("tool", sorted(FIXTURES))
def test_the_reply_is_strict_json(tool):
    data, config = _load(tool)
    report = ctrl.analyze(data, config, include_arrays=True, max_points=200)
    json.dumps(report, allow_nan=False)      # must not raise


# ---------------------------------------------------------------------------
# The one that matters: analyze and batch must not drift
# ---------------------------------------------------------------------------

#: How to run each fixture through pyirena.batch, and which quality figure to
#: compare — batch reports absolute chi-squared for some tools and reduced for
#: others, so the analyze side is read with the matching key.
_BATCH = {
    "unified_fit":  ("fit_unified",            lambda r: r["fit_result"]["chi_squared"],
                     "chi_squared"),
    "sizes":        ("fit_sizes",              lambda r: r["fit_result"]["chi_squared"],
                     "chi_squared"),
    "simple_fits":  ("fit_simple_from_config", lambda r: r["reduced_chi2"],
                     "reduced_chi_squared"),
    "modeling":     ("fit_modeling",           lambda r: r["result"].reduced_chi_squared,
                     "reduced_chi_squared"),
    "waxs_peakfit": ("fit_waxs",               lambda r: r["reduced_chi2"],
                     "reduced_chi_squared"),
    "carbon_fit":   ("fit_carbon",             lambda r: r["reduced_chi_squared"],
                     "reduced_chi_squared"),
}


@pytest.mark.parametrize("tool", sorted(FIXTURES))
def test_analyze_agrees_with_the_batch_path(tool):
    """Same data, same config, same numbers — or the two have drifted."""
    import pyirena.batch as batch

    data_file, config_file = FIXTURES[tool]
    fn_name, chi_of, quality_key = _BATCH[tool]
    expected = chi_of(getattr(batch, fn_name)(
        str(SCRIPTING / data_file), str(SCRIPTING / config_file),
        save_to_nexus=False,
    ))

    data, config = _load(tool)
    got = ctrl.analyze(data, config)["quality"][quality_key]

    # They share core/tool_config.py, so this should be exact bar the data
    # round-tripping through JSON lists.
    assert got == pytest.approx(expected, rel=1e-9), (
        f"{tool}: analyze got {got}, pyirena.batch got {expected}"
    )


# ---------------------------------------------------------------------------
# A config is more than model state
# ---------------------------------------------------------------------------

def test_the_configs_q_range_is_applied_not_ignored():
    """Sizes fits a narrow window; using the whole curve is a different fit."""
    data, config = _load("sizes")
    section = config["sizes"]
    assert section["cursor_q_min"] and section["cursor_q_max"]

    report = ctrl.analyze(data, config)
    info = report["analyze"]
    assert info["fit_q_min"] == pytest.approx(section["cursor_q_min"])
    assert info["fit_q_max"] == pytest.approx(section["cursor_q_max"])
    # Materially narrower than the data, or this test proves nothing.
    assert info["fit_q_max"] < report["data"]["q_max"]


def test_a_q_range_that_misses_the_curve_is_reported_not_silently_applied():
    data, config = _load("sizes")
    config = json.loads(json.dumps(config))
    config["sizes"]["cursor_q_min"] = 500.0
    config["sizes"]["cursor_q_max"] = 900.0

    report = ctrl.analyze(data, config)
    assert "error" not in report
    assert report["analyze"]["fit_q_min"] is None
    assert any("does not overlap" in n for n in report["analyze"]["notes"])


def test_the_sizes_background_prefits_run():
    """The config asks for them; skipping them changes the answer silently."""
    data, config = _load("sizes")
    notes = ctrl.analyze(data, config)["analyze"]["notes"]
    assert any("Power-law pre-fit" in n for n in notes)
    assert any("Background pre-fit" in n for n in notes)


def test_the_waxs_peak_presearch_runs():
    data, config = _load("waxs_peakfit")
    notes = ctrl.analyze(data, config)["analyze"]["notes"]
    assert any("Q0 scan" in n or "Cross-correlation" in n for n in notes)


def test_held_parameters_stay_held():
    data, config = _load("simple_fits")
    section = config["simple_fits"]
    fixed = [n for n, held in (section.get("param_fixed") or {}).items() if held]
    if not fixed:
        pytest.skip("this exported config holds no parameter fixed")

    report = ctrl.analyze(data, config)
    assert set(report["analyze"]["fixed_parameters"]) == set(fixed)
    values = {p["name"]: p["value"] for p in report["results"]["parameters"]}
    for name in fixed:
        assert values[name] == pytest.approx(section["params"][name], rel=1e-9)


# ---------------------------------------------------------------------------
# Errors, and never leaking a session
# ---------------------------------------------------------------------------

def test_a_config_for_the_wrong_shape_of_object():
    data, _ = _load("sizes")
    assert ctrl.analyze(data, "not a config")["code"] == "BAD_CONFIG"
    assert ctrl.analyze("not data", {"sizes": {}})["code"] == "BAD_ARGUMENTS"


def test_a_config_with_no_recognisable_tool_section():
    data, _ = _load("sizes")
    result = ctrl.analyze(data, {"_pyirena_config": {}, "not_a_tool": {}})
    assert result["code"] == "BAD_CONFIG"
    assert "unified_fit" in result["suggestion"]


def test_an_unknown_tool_name():
    data, config = _load("sizes")
    assert ctrl.analyze(data, config, tool="telepathy")["code"] == "BAD_CONFIG"


def test_bad_data_is_refused_before_any_fitting():
    _, config = _load("sizes")
    result = ctrl.analyze({"q": [1, 2, 3], "intensity": [1, 2]}, config)
    assert result["code"] == "SHAPE_MISMATCH"


def test_a_failing_fit_still_closes_its_session():
    """The autouse fixture checks the leak; this makes the failure happen."""
    data, config = _load("modeling")
    config = json.loads(json.dumps(config))
    config["modeling"]["populations"] = []
    ctrl.analyze(data, config)          # error or not, no session may survive


def test_a_q_range_matching_the_data_is_not_reported_as_clipped():
    """A config usually carries the data's own limits, a few ulps out."""
    data, config = _load("modeling")
    notes = ctrl.analyze(data, config)["analyze"]["notes"]
    assert not any("clipped" in n for n in notes), notes


def test_a_genuinely_narrower_range_is_still_reported():
    data, config = _load("modeling")
    config = json.loads(json.dumps(config))
    config["modeling"]["q_min"] = 1e-9          # far below the data
    config["modeling"]["q_max"] = 1e3           # far above
    notes = ctrl.analyze(data, config)["analyze"]["notes"]
    assert any("clipped" in n for n in notes), notes


# ---------------------------------------------------------------------------
# The fixtures are read-only inputs
# ---------------------------------------------------------------------------

def test_save_to_nexus_false_leaves_the_data_file_alone():
    """It was advertised and ignored, so fits wrote into the input file.

    ``fit_simple_from_config`` took ``save_to_nexus`` and never passed it on,
    and ``fit_simple`` had no such parameter at all and always wrote. Running
    these fixtures during development silently rewrote the
    ``simple_fit_results`` group of a checked-in data file.
    """
    import hashlib

    import pyirena.batch as batch

    data_file = SCRIPTING / FIXTURES["simple_fits"][0]
    before = hashlib.sha256(data_file.read_bytes()).hexdigest()

    result = batch.fit_simple_from_config(
        str(data_file), str(SCRIPTING / FIXTURES["simple_fits"][1]),
        save_to_nexus=False,
    )
    assert result and result["success"]
    assert hashlib.sha256(data_file.read_bytes()).hexdigest() == before, (
        f"{data_file.name} was modified despite save_to_nexus=False"
    )


@pytest.mark.parametrize("tool", sorted(FIXTURES))
def test_analyze_never_touches_a_file(tool):
    """analyze works from arrays, so it has no file to write to at all."""
    import hashlib

    data_file = SCRIPTING / FIXTURES[tool][0]
    before = hashlib.sha256(data_file.read_bytes()).hexdigest()
    data, config = _load(tool)
    assert "error" not in ctrl.analyze(data, config)
    assert hashlib.sha256(data_file.read_bytes()).hexdigest() == before
