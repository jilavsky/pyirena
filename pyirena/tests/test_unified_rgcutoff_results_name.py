"""Saved Unified Fit results use one name for the cutoff radius: ``RgCutoff``.

Settings (``model.to_dict()``, config JSON) call it ``RgCO``; the saved
results group ``entry/unified_fit_results/level_<n>/`` is a separate, frozen
vocabulary that calls it ``RgCutoff`` (docs/developer/adding-features.md).
The GUI and ``pyirena.batch`` always wrote ``RgCutoff``, but
``pyirena.api.control.save_fit`` wrote ``RgCO`` until 2026-11-16, so every
results reader (the Fractals target detection, the Data Selector table, Igor
export, MCP tabulation) silently lost the cutoff for agent-written files.

These tests pin both halves of the fix: the writer now uses ``RgCutoff``, and
the readers still understand files that carry ``RgCO``.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import h5py
import pytest

from pyirena.api.aggregate import _extract_scalar
from pyirena.io.nxcansas_unified import load_unified_fit_results

DATA = Path(__file__).resolve().parents[2] / "validationData" / "unified" / "unified_two_level.h5"

pytestmark = pytest.mark.skipif(not DATA.is_file(), reason="validationData not present")


def _save_two_level_fit(tmp_path: Path) -> Path:
    """Fit the two-level validation curve with RgCO linked and save via the api."""
    from pyirena.api.control import (
        close_session,
        open_dataset,
        run_fit,
        save_fit,
        select_model,
        set_level_option,
        set_parameter_value,
    )

    src = tmp_path / DATA.name
    shutil.copy2(DATA, src)
    sid = open_dataset(str(src))["session_id"]
    try:
        select_model(sid, nlevels=2)
        set_level_option(sid, 2, "link_RGCO", True)
        for name, value in {"G_1": 8, "Rg_1": 120, "B_1": 1.9e-7, "P_1": 4,
                            "G_2": 4000, "Rg_2": 1200, "B_2": 1.4e-6, "P_2": 3.2}.items():
            set_parameter_value(sid, name, value)
        assert run_fit(sid)["success"]
        out = tmp_path / "fit.h5"
        save_fit(sid, output_path=str(out))
    finally:
        close_session(sid)
    return out


def _rename_to_old_name(path: Path) -> None:
    """Make *path* look like a file written by save_fit before the fix."""
    with h5py.File(path, "r+") as f:
        for n in (1, 2):
            lg = f[f"entry/unified_fit_results/level_{n}"]
            lg.move("RgCutoff", "RgCO")


def test_api_save_fit_writes_rgcutoff(tmp_path):
    out = _save_two_level_fit(tmp_path)
    with h5py.File(out, "r") as f:
        lg = f["entry/unified_fit_results/level_2"]
        assert "RgCutoff" in lg and "RgCO" not in lg
        rgco = float(lg["RgCutoff"][()])
        rg1 = float(f["entry/unified_fit_results/level_1/Rg"][()])
    assert rgco == pytest.approx(rg1, rel=1e-6)      # linked to Rg of level 1


def test_loader_reads_old_rgco_files_as_rgcutoff(tmp_path):
    out = _save_two_level_fit(tmp_path)
    expected = load_unified_fit_results(out)["levels"][1]["RgCutoff"]
    _rename_to_old_name(out)
    levels = load_unified_fit_results(out)["levels"]
    assert levels[1]["RgCutoff"] == pytest.approx(expected)
    assert "RgCO" not in levels[1]


def test_schema_readers_find_old_rgco_dataset(tmp_path):
    out = _save_two_level_fit(tmp_path)
    _rename_to_old_name(out)
    with h5py.File(out, "r") as f:
        val, _std, units = _extract_scalar(f, "unified_fit", "RgCutoff", subgroup_index=2)
        rg1 = float(f["entry/unified_fit_results/level_1/Rg"][()])
    assert val == pytest.approx(rg1, rel=1e-6)
    assert units == "Å"


def test_fractals_targets_from_api_written_file(tmp_path):
    """The Fractals panel finds the mass-fractal pair (RgCutoff_2 ≈ Rg_1)."""
    # Not pytest.importorskip: since pytest 9.1 it skips only on
    # ModuleNotFoundError, and pyirena.gui._qt raises a plain ImportError
    # (with an install hint) when no Qt binding is installed.
    try:
        from pyirena.gui.fractals_panel import FractalsPanel
    except ImportError as exc:
        pytest.skip(f"Qt GUI not available: {exc}")

    out = _save_two_level_fit(tmp_path)
    for rename in (False, True):          # new file, then an old RgCO file
        if rename:
            _rename_to_old_name(out)
        levels = load_unified_fit_results(out)["levels"]
        targets = FractalsPanel._derive_targets_from_levels(levels)
        assert targets, f"no targets (old-name file: {rename})"
        assert targets["level_low"] == 1 and targets["level_high"] == 2
