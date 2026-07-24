"""Behaviour tests for the control-surface mutating tools (save + open).

Covers three review findings on ``pyirena.api.control`` (also the MCP write
surface):

* **Path confinement** — ``open_dataset`` (read) and ``save_fit`` /
  ``save_sizes_fit`` (write) must honour ``PYIRENA_DATA_ROOT``.
* **Complete output files** — saving to a *new* ``output_path`` must produce a
  re-openable NXcanSAS file (reduced data + results), leaving the original
  untouched — not a results-only stub.
* **Slit-smearing provenance** — a slit-smeared fit must be saved with its
  ``slit_length`` and an ideal (pinhole) model curve, not silently as pinhole.
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

import pyirena.api.control as ctrl
from pyirena.api.control.session import get_session
from pyirena.io.hdf5 import readGenericNXcanSAS
from pyirena.io.nxcansas_unified import create_nxcansas_file


def _make_source(folder: Path, name: str = "samp.h5") -> Path:
    q = np.logspace(-3, 0, 200)
    intensity = 0.01 / (1 + (q * 50) ** 2) ** 2 + 1e-4
    fp = folder / name
    create_nxcansas_file(fp, q, intensity, error=intensity * 0.03, sample_name=name)
    return fp


def _fit_unified(sid):
    assert "error" not in ctrl.select_model(sid, "unified_fit")
    assert "error" not in ctrl.add_unified_level(sid)
    assert "error" not in ctrl.run_fit(sid)


def _fit_sizes(sid):
    assert "error" not in ctrl.select_sizes_model(sid, method="maxent")
    assert "error" not in ctrl.set_size_grid(sid, r_min=10, r_max=500, n_bins=40)
    assert "error" not in ctrl.set_shape(sid, "sphere", contrast=1.0)
    assert "error" not in ctrl.run_sizes_fit(sid)


# ── Path confinement ─────────────────────────────────────────────────────────

def test_open_dataset_rejects_path_outside_root(tmp_path, monkeypatch):
    root = tmp_path / "root"
    root.mkdir()
    outside = _make_source(tmp_path, "outside.h5")   # not under root
    monkeypatch.setenv("PYIRENA_DATA_ROOT", str(root))

    res = ctrl.open_dataset(str(outside))
    assert res.get("code") == "PATH_NOT_ALLOWED"
    assert "session_id" not in res


def test_open_dataset_rejects_traversal(tmp_path, monkeypatch):
    root = tmp_path / "root"
    root.mkdir()
    monkeypatch.setenv("PYIRENA_DATA_ROOT", str(root))
    assert ctrl.open_dataset("../escape.h5").get("code") == "PATH_NOT_ALLOWED"


def test_save_fit_rejects_output_outside_root(tmp_path, monkeypatch):
    root = tmp_path / "root"
    root.mkdir()
    src = _make_source(root)                          # source inside root
    sid = ctrl.open_dataset(str(src))["session_id"]
    _fit_unified(sid)
    monkeypatch.setenv("PYIRENA_DATA_ROOT", str(root))

    res = ctrl.save_fit(sid, output_path=str(tmp_path / "escape.h5"))
    assert res.get("code") == "PATH_NOT_ALLOWED"


def test_save_sizes_fit_rejects_output_outside_root(tmp_path, monkeypatch):
    root = tmp_path / "root"
    root.mkdir()
    src = _make_source(root)
    sid = ctrl.open_dataset(str(src))["session_id"]
    _fit_sizes(sid)
    monkeypatch.setenv("PYIRENA_DATA_ROOT", str(root))

    res = ctrl.save_sizes_fit(sid, output_path=str(tmp_path / "escape.h5"))
    assert res.get("code") == "PATH_NOT_ALLOWED"


# ── Complete output files (output_path preserves original & is re-openable) ───

def test_unified_output_path_is_complete_and_preserves_original(tmp_path):
    src = _make_source(tmp_path)
    sid = ctrl.open_dataset(str(src))["session_id"]
    _fit_unified(sid)

    out = tmp_path / "out_new.h5"
    assert "error" not in ctrl.save_fit(sid, output_path=str(out))

    # New file is a complete, re-openable NXcanSAS file.
    data = readGenericNXcanSAS(str(out.parent), out.name)
    assert data is not None and len(data.get("Q", [])) > 0
    with h5py.File(out) as f:
        assert "entry/unified_fit_results" in f

    # Original is untouched (no results written into it).
    with h5py.File(src) as f:
        assert "entry/unified_fit_results" not in f


def test_sizes_output_path_is_complete_and_preserves_original(tmp_path):
    src = _make_source(tmp_path)
    sid = ctrl.open_dataset(str(src))["session_id"]
    _fit_sizes(sid)

    out = tmp_path / "sizes_out.h5"
    assert "error" not in ctrl.save_sizes_fit(sid, output_path=str(out))

    data = readGenericNXcanSAS(str(out.parent), out.name)
    assert data is not None and len(data.get("Q", [])) > 0
    with h5py.File(out) as f:
        assert "entry/sizes_results" in f
    with h5py.File(src) as f:
        assert "entry/sizes_results" not in f


# ── Slit-smearing provenance survives a control-driven save ──────────────────

def test_unified_save_records_slit_provenance(tmp_path):
    src = _make_source(tmp_path)
    sid = ctrl.open_dataset(str(src))["session_id"]
    _fit_unified(sid)
    # Simulate a slit-smeared session (open_dataset(use_slit_smeared=True) sets
    # these from the file; here we set them directly to isolate the save path).
    s = get_session(sid)
    s.model.use_slit_smearing = True
    s.model.slit_length = 0.018
    assert "error" not in ctrl.run_fit(sid)

    out = tmp_path / "smeared_unified.h5"
    assert "error" not in ctrl.save_fit(sid, output_path=str(out))
    with h5py.File(out) as f:
        g = f["entry/unified_fit_results"]
        assert float(g.attrs["slit_length"]) == pytest.approx(0.018)
        assert "intensity_model_ideal" in g


def test_sizes_save_records_slit_provenance(tmp_path):
    src = _make_source(tmp_path)
    sid = ctrl.open_dataset(str(src))["session_id"]
    ctrl.select_sizes_model(sid, method="maxent")
    ctrl.set_size_grid(sid, r_min=10, r_max=500, n_bins=40)
    ctrl.set_shape(sid, "sphere", contrast=1.0)
    s = get_session(sid)
    s.model.use_slit_smearing = True
    s.model.slit_length = 0.02
    assert "error" not in ctrl.run_sizes_fit(sid)

    out = tmp_path / "smeared_sizes.h5"
    assert "error" not in ctrl.save_sizes_fit(sid, output_path=str(out))
    with h5py.File(out) as f:
        g = f["entry/sizes_results"]
        assert float(g.attrs["slit_length"]) == pytest.approx(0.02)
        assert bool(g.attrs["data_is_slit_smeared"]) is True
        assert "intensity_model_ideal" in g
