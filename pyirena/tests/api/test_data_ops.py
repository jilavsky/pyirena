"""Tests for pyirena.api.data_ops — the data manipulation / merge surface.

These are the only api functions that write data files, so every test uses
a function-scoped ``tmp_path``. Covers the output-path convention (the
sibling ``_manip`` / ``_merged`` folder plus the per-operation filename
suffix, previously untested anywhere), PYIRENA_DATA_ROOT confinement, and
the guards this layer adds on top of core.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pyirena.api import (
    average_data,
    divide_data,
    match_merge_files,
    merge_datasets,
    rebin_data,
    scale_data,
    subtract_data,
    trim_data,
)
from pyirena.api.data_op_schemas import DATA_OP_SCHEMA_BY_NAME, DATA_OP_TOOL_SCHEMAS
from pyirena.io.hdf5 import readGenericNXcanSAS
from pyirena.io.nxcansas_unified import create_nxcansas_file


def _curve(q: np.ndarray, mult: float = 1.0) -> np.ndarray:
    return mult * (1000 * np.exp(-((q * 80) ** 2) / 3) + 5 * q**-4 + 0.5)


def _make(folder: Path, name: str, mult: float = 1.0, flat_only: bool = False) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    q = np.logspace(-3, -0.5, 300)
    intensity = (5 * q**-4 + 0.5) if flat_only else _curve(q, mult)
    fp = folder / name
    create_nxcansas_file(fp, q, intensity, error=intensity * 0.02, sample_name=name)
    return fp


@pytest.fixture
def run(tmp_path: Path) -> Path:
    """A folder of three near-identical frames plus a buffer."""
    folder = tmp_path / "run"
    for i, mult in enumerate((1.0, 1.02, 0.98), start=1):
        _make(folder, f"s{i}.h5", mult)
    _make(folder, "buffer.h5", flat_only=True)
    return folder


def _assert_json_safe(result):
    json.dumps(result)


def _n_written(path: str) -> int:
    p = Path(path)
    data = readGenericNXcanSAS(str(p.parent), p.name)
    return int(np.asarray(data["Q"]).size)


# ---------------------------------------------------------------------------
# Output-path convention
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "call,suffix",
    [
        (lambda run: scale_data(str(run / "s1.h5"), scale_I=2.0), "_scaled"),
        (lambda run: trim_data(str(run / "s1.h5"), q_min=0.01, q_max=0.1), "_trimmed"),
        (lambda run: rebin_data(str(run / "s1.h5"), n_points=50), "_rebinned"),
        (
            lambda run: subtract_data(str(run / "s1.h5"), str(run / "buffer.h5")),
            "_sub",
        ),
        (
            lambda run: divide_data(str(run / "s1.h5"), str(run / "buffer.h5")),
            "_div",
        ),
        (
            lambda run: average_data([str(run / f"s{i}.h5") for i in (1, 2, 3)]),
            "_avg",
        ),
    ],
)
def test_output_lands_in_sibling_manip_folder_with_operation_suffix(run, call, suffix):
    result = call(run)
    assert "error" not in result, result
    out = Path(result["output_path"])
    assert out.name == f"s1{suffix}.h5"
    # sibling of the source folder, not a subfolder
    assert out.parent.name == f"{run.name}_manip"
    assert out.parent.parent == run.parent
    assert out.is_file()
    _assert_json_safe(result)


def test_different_operations_do_not_collide(run):
    a = scale_data(str(run / "s1.h5"), scale_I=2.0)
    b = trim_data(str(run / "s1.h5"), q_min=0.01, q_max=0.1)
    assert a["output_path"] != b["output_path"]
    assert Path(a["output_path"]).is_file() and Path(b["output_path"]).is_file()


def test_repeating_an_operation_overwrites(run):
    first = scale_data(str(run / "s1.h5"), scale_I=2.0)
    second = scale_data(str(run / "s1.h5"), scale_I=4.0)
    assert first["output_path"] == second["output_path"]
    out = Path(second["output_path"])
    data = readGenericNXcanSAS(str(out.parent), out.name)
    # the second call's factor won
    source = readGenericNXcanSAS(str(run), "s1.h5")
    ratio = np.asarray(data["Intensity"])[0] / np.asarray(source["Intensity"])[0]
    assert ratio == pytest.approx(4.0, rel=1e-6)


def test_explicit_output_folder_is_honoured(run, tmp_path):
    dest = tmp_path / "elsewhere"
    result = scale_data(str(run / "s1.h5"), output_folder=str(dest))
    assert Path(result["output_path"]).parent == dest


# ---------------------------------------------------------------------------
# Numerical behaviour
# ---------------------------------------------------------------------------

def test_scale_applies_factor_then_background(run):
    result = scale_data(str(run / "s1.h5"), scale_I=2.0, background=0.1)
    out = Path(result["output_path"])
    written = readGenericNXcanSAS(str(out.parent), out.name)
    source = readGenericNXcanSAS(str(run), "s1.h5")
    expected = 2.0 * np.asarray(source["Intensity"]) - 0.1
    # save strips non-positive points, so compare the surviving head
    n = int(np.asarray(written["Q"]).size)
    np.testing.assert_allclose(np.asarray(written["Intensity"])[:n], expected[:n], rtol=1e-9)


def test_trim_keeps_only_the_window(run):
    result = trim_data(str(run / "s1.h5"), q_min=0.01, q_max=0.1)
    out = Path(result["output_path"])
    q = np.asarray(readGenericNXcanSAS(str(out.parent), out.name)["Q"])
    assert q.min() >= 0.01 and q.max() <= 0.1
    assert result["n_points_written"] < result["n_points_in"]


def test_rebin_respects_n_points(run):
    result = rebin_data(str(run / "s1.h5"), n_points=50)
    assert result["n_points_written"] <= 50


def test_average_of_identical_frames_is_the_frame(tmp_path):
    folder = tmp_path / "same"
    for i in (1, 2, 3):
        _make(folder, f"s{i}.h5", mult=1.0)
    result = average_data([str(folder / f"s{i}.h5") for i in (1, 2, 3)])
    out = Path(result["output_path"])
    written = readGenericNXcanSAS(str(out.parent), out.name)
    source = readGenericNXcanSAS(str(folder), "s1.h5")
    n = int(np.asarray(written["Q"]).size)
    np.testing.assert_allclose(
        np.asarray(written["Intensity"])[:n], np.asarray(source["Intensity"])[:n], rtol=1e-9
    )
    assert result["n_datasets"] == 3


def test_subtract_reports_stripped_nonpositive_points(run):
    """Over-subtraction silently loses points on save — it must be reported."""
    result = subtract_data(str(run / "s1.h5"), str(run / "buffer.h5"))
    assert result["n_dropped_nonpositive"] > 0
    assert any("over-subtraction" in w for w in result["warnings"])
    assert result["n_points_written"] < result["n_points_in"]


def test_merge_recovers_a_known_scale(tmp_path):
    """DS2 generated 1.15x too high must be scaled back by 1/1.15."""
    lo, hi = tmp_path / "usaxs", tmp_path / "saxs"
    lo.mkdir(); hi.mkdir()

    def curve(q):
        return 1000 * np.exp(-((q * 300) ** 2) / 3) + 2e-3 * q**-3.2 + 0.05

    q1 = np.logspace(-4, -1.4, 200)
    q2 = np.logspace(-2.2, -0.5, 200)
    create_nxcansas_file(lo / "s_usaxs_001.h5", q1, curve(q1), curve(q1) * 0.02, sample_name="u")
    i2 = curve(q2) * 1.15
    create_nxcansas_file(hi / "s_saxs_001.h5", q2, i2, i2 * 0.02, sample_name="s")

    result = merge_datasets(str(lo / "s_usaxs_001.h5"), str(hi / "s_saxs_001.h5"))
    assert "error" not in result, result
    assert result["scale"] == pytest.approx(1 / 1.15, rel=0.02)
    out = Path(result["output_path"])
    assert out.name == "s_usaxs_001_merged.h5"
    assert out.parent.name == "usaxs_merged"
    assert "background" in result and "note" in result
    _assert_json_safe(result)


def test_match_merge_files_pairs_by_index(tmp_path):
    lo, hi = tmp_path / "usaxs", tmp_path / "saxs"
    for i in (1, 2):
        _make(lo, f"sample_usaxs_{i:03d}.h5")
        _make(hi, f"sample_saxs_{i:03d}.h5")
    _make(lo, "sample_usaxs_099.h5")  # no partner

    result = match_merge_files(str(lo), str(hi))
    assert result["n_pairs"] == 2
    assert result["unmatched_1"] == ["sample_usaxs_099.h5"]
    assert Path(result["pairs"][0]["file1"]).name.startswith("sample_usaxs")
    assert Path(result["pairs"][0]["file2"]).name.startswith("sample_saxs")
    _assert_json_safe(result)


# ---------------------------------------------------------------------------
# Guards this layer adds
# ---------------------------------------------------------------------------

def test_trim_to_an_empty_window_errors(run):
    result = trim_data(str(run / "s1.h5"), q_min=5.0, q_max=9.0)
    assert result["code"] == "EMPTY_RESULT"
    assert "1/Å" in result["suggestion"] or "1/A" in result["suggestion"]


def test_rebin_rejects_unknown_mode(run):
    assert rebin_data(str(run / "s1.h5"), mode="bogus")["code"] == "BAD_MODE"


def test_rebin_reference_mode_needs_a_reference_file(run):
    assert rebin_data(str(run / "s1.h5"), mode="reference")["code"] == "MISSING_REFERENCE"


def test_rebin_reference_mode_uses_the_reference_grid(run):
    ref = rebin_data(str(run / "s2.h5"), n_points=40)
    result = rebin_data(
        str(run / "s1.h5"), mode="reference", reference_file=ref["output_path"]
    )
    assert "error" not in result
    assert result["n_points_written"] <= 40


def test_log_rebin_rejects_nonpositive_q_min(run):
    assert rebin_data(str(run / "s1.h5"), mode="log", q_min=0.0)["code"] == "BAD_Q_RANGE"


def test_auto_scale_without_a_window_errors(run):
    """core silently ignores auto_scale unless BOTH bounds are given."""
    result = subtract_data(
        str(run / "s1.h5"), str(run / "buffer.h5"), auto_scale=True
    )
    assert result["code"] == "MISSING_AUTO_RANGE"


def test_auto_scale_with_a_window_is_accepted(run):
    result = subtract_data(
        str(run / "s1.h5"), str(run / "buffer.h5"),
        auto_scale=True, auto_q_min=0.1, auto_q_max=0.3,
    )
    assert "error" not in result
    assert result["parameters"]["auto_scale"] is True


def test_average_needs_two_files(run):
    assert average_data([str(run / "s1.h5")])["code"] == "TOO_FEW_FILES"


def test_missing_file_reports_not_found(run):
    assert scale_data(str(run / "nope.h5"))["code"] == "FILE_NOT_FOUND"


def test_merge_rejects_swapped_inputs(tmp_path):
    lo, hi = tmp_path / "usaxs", tmp_path / "saxs"
    lo.mkdir(); hi.mkdir()
    q1 = np.logspace(-4, -1.4, 100)
    q2 = np.logspace(-2.2, -0.5, 100)
    create_nxcansas_file(lo / "u.h5", q1, _curve(q1), _curve(q1) * 0.02, sample_name="u")
    create_nxcansas_file(hi / "s.h5", q2, _curve(q2), _curve(q2) * 0.02, sample_name="s")
    result = merge_datasets(str(hi / "s.h5"), str(lo / "u.h5"))
    assert result["code"] == "SWAPPED_INPUTS"
    assert "u.h5" in result["suggestion"]


def test_merge_rejects_non_integer_scale_dataset(tmp_path):
    lo, hi = tmp_path / "usaxs", tmp_path / "saxs"
    lo.mkdir(); hi.mkdir()
    q1 = np.logspace(-4, -1.4, 100)
    q2 = np.logspace(-2.2, -0.5, 100)
    create_nxcansas_file(lo / "u.h5", q1, _curve(q1), _curve(q1) * 0.02, sample_name="u")
    create_nxcansas_file(hi / "s.h5", q2, _curve(q2), _curve(q2) * 0.02, sample_name="s")
    result = merge_datasets(str(lo / "u.h5"), str(hi / "s.h5"), scale_dataset="DS2")
    assert result["code"] == "BAD_ARGUMENTS"


def test_errors_are_json_safe(run):
    _assert_json_safe(trim_data(str(run / "s1.h5"), q_min=5.0, q_max=9.0))


# ---------------------------------------------------------------------------
# Slit smearing
# ---------------------------------------------------------------------------

def _make_smeared(folder: Path, name: str, slit: float = 0.018, mult: float = 1.0) -> Path:
    from pyirena.io._nxcansas_common import append_dql

    fp = _make(folder, name, mult)
    append_dql(fp, slit)
    return fp


def test_average_of_smeared_data_keeps_the_slit_length(tmp_path):
    """A slit-smeared file returns dQ as a scalar, which core would index."""
    folder = tmp_path / "smeared"
    for i in (1, 2, 3):
        _make_smeared(folder, f"s{i}.h5", mult=1.0 + 0.01 * i)
    result = average_data([str(folder / f"s{i}.h5") for i in (1, 2, 3)])
    assert "error" not in result, result
    out = Path(result["output_path"])
    back = readGenericNXcanSAS(str(out.parent), out.name)
    assert back.get("is_slit_smeared") is True
    assert float(back.get("slit_length")) == pytest.approx(0.018)


def test_trim_of_smeared_data_succeeds(tmp_path):
    folder = tmp_path / "smeared"
    _make_smeared(folder, "s1.h5")
    result = trim_data(str(folder / "s1.h5"), q_min=0.01, q_max=0.1)
    assert "error" not in result, result


@pytest.mark.parametrize("op", ["subtract", "average"])
def test_mixed_slit_status_is_refused(tmp_path, op):
    """core guards subtract/divide but NOT average — this layer adds it."""
    folder = tmp_path / "mixed"
    smeared = _make_smeared(folder, "smeared.h5")
    pinhole = _make(folder, "pinhole.h5")
    if op == "subtract":
        result = subtract_data(str(smeared), str(pinhole))
    else:
        result = average_data([str(smeared), str(pinhole)])
    assert result["code"] == "SLIT_MISMATCH"


# ---------------------------------------------------------------------------
# Sandbox confinement
# ---------------------------------------------------------------------------

def test_default_sibling_folder_blocked_by_data_root(run, monkeypatch):
    """The derived sibling of an in-root input can itself be out of root."""
    monkeypatch.setenv("PYIRENA_DATA_ROOT", str(run))
    result = scale_data(str(run / "s1.h5"))
    assert result["code"] == "PATH_NOT_ALLOWED"
    assert "output_folder" in result["suggestion"]


def test_explicit_in_root_output_folder_is_allowed(run, monkeypatch):
    monkeypatch.setenv("PYIRENA_DATA_ROOT", str(run))
    result = scale_data(str(run / "s1.h5"), output_folder=str(run / "out"))
    assert "error" not in result, result
    assert Path(result["output_path"]).parent == run / "out"


def test_input_outside_data_root_is_rejected(run, tmp_path, monkeypatch):
    inside = tmp_path / "inside"
    inside.mkdir()
    monkeypatch.setenv("PYIRENA_DATA_ROOT", str(inside))
    assert scale_data(str(run / "s1.h5"))["code"] == "PATH_NOT_ALLOWED"


# ---------------------------------------------------------------------------
# Schema parity
# ---------------------------------------------------------------------------

def _sig_params(fn) -> dict[str, bool]:
    import inspect

    return {
        name: p.default is not inspect.Parameter.empty
        for name, p in inspect.signature(fn).parameters.items()
    }


@pytest.mark.parametrize("schema", DATA_OP_TOOL_SCHEMAS, ids=lambda s: s["name"])
def test_data_op_schema_matches_callable_signature(schema):
    import pyirena.api.data_ops as data_ops

    fn = getattr(data_ops, schema["name"], None)
    assert callable(fn), f"schema '{schema['name']}' has no data_ops callable"

    params = _sig_params(fn)
    mandatory = {n for n, has_default in params.items() if not has_default}
    ins = schema["input_schema"]
    props = set(ins.get("properties", {}))
    required = set(ins.get("required", []))

    assert ins["type"] == "object"
    assert schema["description"].strip()
    assert props <= set(params), \
        f"{schema['name']}: schema exposes non-parameters {props - set(params)}"
    assert required == mandatory, \
        f"{schema['name']}: required {required} != mandatory params {mandatory}"
    assert set(params) <= props, \
        f"{schema['name']}: parameters missing from schema {set(params) - props}"


def test_every_data_op_is_exported_from_the_api_facade():
    import pyirena.api as papi

    for name in DATA_OP_SCHEMA_BY_NAME:
        assert name in papi.__all__, f"{name} missing from pyirena.api.__all__"
        assert callable(getattr(papi, name))


def test_importing_data_ops_writes_nothing_to_stdout(capfd):
    """batch attaches a stdout log handler; the api layer must not.

    MCP speaks JSON-RPC over stdout, so a stray handler corrupts it.
    """
    import importlib

    import pyirena.api.data_ops as data_ops

    importlib.reload(data_ops)
    out, _ = capfd.readouterr()
    assert out == ""
