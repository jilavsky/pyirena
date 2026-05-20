"""Shared fixtures for pyirena.api tests.

Builds a small synthetic NXcanSAS HDF5 file containing:
- reduced data (sasdata: Q, I, Idev),
- simple_fit_results group,
- unified_fit_results group (1 level),
- sizes_results group,
- data_merge_results group (provenance only).

The fixture is session-scoped — one file used by every test in this dir.
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest


@pytest.fixture(scope="session")
def synth_nxcansas_file(tmp_path_factory) -> Path:
    """Create a synthetic NXcanSAS file with several pyirena result groups."""
    folder = tmp_path_factory.mktemp("synth_data")
    fp = folder / "sample_A_scan_42.h5"

    # Reduced data
    q = np.logspace(-3, 0, 300)
    intensity = 1e-2 / (1 + (q * 50) ** 2) ** 2 + 1e-4
    err = 0.05 * intensity

    with h5py.File(fp, "w") as f:
        # Root
        f.attrs["default"] = "entry"
        f.attrs["file_name"] = fp.name
        f.attrs["file_time"] = "2026-05-19T10:00:00"
        f.attrs["creator"] = "pyirena-test"
        f.attrs["NeXus_version"] = "4.3.0"

        # /entry NXsas wrapper
        entry = f.create_group("entry")
        entry.attrs["NX_class"] = "NXentry"
        entry.attrs["canSAS_class"] = "SASentry"
        entry.attrs["default"] = "sample_A"
        entry.create_dataset("definition", data="NXsas")

        # /entry/sample_A NXsubentry → NXcanSAS data
        sub = entry.create_group("sample_A")
        sub.attrs["NX_class"] = "NXsubentry"
        sub.attrs["canSAS_class"] = "SASentry"
        sub.attrs["default"] = "sasdata"
        sub.attrs["title"] = "sample_A"
        sub.create_dataset("definition", data="NXcanSAS")
        sub.create_dataset("title", data="sample_A")
        sub.create_dataset("run", data="test_run")

        sasdata = sub.create_group("sasdata")
        sasdata.attrs["NX_class"] = "NXdata"
        sasdata.attrs["canSAS_class"] = "SASdata"
        sasdata.attrs["signal"] = "I"
        sasdata.attrs["I_axes"] = "Q"

        ds_i = sasdata.create_dataset("I", data=intensity)
        ds_i.attrs["units"] = "1/cm"
        ds_i.attrs["uncertainties"] = "Idev"
        ds_i.attrs["thickness"] = 0.1
        ds_i.attrs["blankname"] = "blank.h5"
        ds_i.attrs["label"] = "sample_A_scan_42"
        ds_q = sasdata.create_dataset("Q", data=q)
        ds_q.attrs["units"] = "1/angstrom"
        sasdata.create_dataset("Idev", data=err).attrs["units"] = "1/cm"

        # Sample group (so _sample_name() finds it)
        sample_grp = entry.create_group("sample")
        sample_grp.attrs["NX_class"] = "NXsample"
        sample_grp.create_dataset("name", data="sample_A")

        # ── simple_fit_results ────────────────────────────────────────────
        sfg = f.create_group("entry/simple_fit_results")
        sfg.attrs["NX_class"] = "NXprocess"
        sfg.attrs["model"] = "Guinier"
        sfg.attrs["success"] = True
        sfg.attrs["dof"] = 297
        sfg.attrs["q_min"] = float(q.min())
        sfg.attrs["q_max"] = float(q.max())
        sfg.attrs["timestamp"] = "2026-05-19T10:05:00"
        sfg.create_dataset("chi_squared", data=1.23)
        sfg.create_dataset("reduced_chi_squared", data=0.0041)
        sfg.create_dataset("Q", data=q)
        sfg.create_dataset("I_model", data=intensity * 0.98)
        sfg.create_dataset("residuals", data=np.random.RandomState(0).randn(q.size))
        sfg.create_dataset("intensity_data", data=intensity)
        sfg.create_dataset("intensity_error", data=err)
        pgrp = sfg.create_group("params")
        pgrp.create_dataset("I0", data=0.01)
        pgrp.create_dataset("Rg", data=50.0)
        # Kp is a model-specific Porod param NOT enumerated in TOOL_REGISTRY —
        # tests use it to exercise the runtime-fallback parameter resolution.
        pgrp.create_dataset("Kp", data=3.7e-5)
        psg = sfg.create_group("params_std")
        psg.create_dataset("I0", data=0.0005)
        psg.create_dataset("Rg", data=1.2)
        psg.create_dataset("Kp", data=2.1e-6)
        dg = sfg.create_group("derived")
        dg.create_dataset("Q_min_used", data=float(q.min()))

        # ── unified_fit_results (1 level) ─────────────────────────────────
        ufg = f.create_group("entry/unified_fit_results")
        ufg.attrs["NX_class"] = "NXprocess"
        ufg.attrs["num_levels"] = 1
        ufg.attrs["timestamp"] = "2026-05-19T10:10:00"
        ufg.attrs["program"] = "pyirena-test"
        ufg.create_dataset("background", data=1e-4)
        ufg.create_dataset("chi_squared", data=2.1)
        ufg.create_dataset("Q", data=q)
        ufg.create_dataset("intensity_data", data=intensity)
        ufg.create_dataset("intensity_model", data=intensity * 0.97)
        ufg.create_dataset("residuals", data=np.random.RandomState(1).randn(q.size))
        ufg.create_dataset("intensity_error", data=err)
        lvg = ufg.create_group("level_1")
        for k, v in [("G", 0.012), ("Rg", 52.0), ("B", 5e-6),
                     ("P", 4.0), ("RgCutoff", 0.0),
                     ("ETA", 0.0), ("PACK", 1.0)]:
            lvg.create_dataset(k, data=float(v))

        # ── sizes_results ─────────────────────────────────────────────────
        szg = f.create_group("entry/sizes_results")
        szg.attrs["NX_class"] = "NXprocess"
        szg.attrs["timestamp"] = "2026-05-19T10:15:00"
        szg.attrs["shape"] = "sphere"
        szg.attrs["method"] = "maxent"
        szg.attrs["r_min"] = 10.0
        szg.attrs["r_max"] = 1000.0
        szg.attrs["n_bins"] = 50
        szg.attrs["log_spacing"] = True
        szg.attrs["aspect_ratio"] = 1.0
        szg.attrs["contrast"] = 1.0
        szg.attrs["background"] = 1e-4
        szg.attrs["power_law_B"] = 0.0
        szg.attrs["power_law_P"] = 0.0
        szg.create_dataset("chi_squared", data=1.5)
        szg.create_dataset("volume_fraction", data=0.05)
        szg.create_dataset("rg", data=55.0)
        szg.create_dataset("n_iterations", data=42)
        szg.create_dataset("q_power", data=4.0)
        szg.create_dataset("Q", data=q)
        szg.create_dataset("intensity_data", data=intensity)
        szg.create_dataset("intensity_model", data=intensity * 0.99)
        szg.create_dataset("residuals", data=np.random.RandomState(2).randn(q.size))
        r_grid = np.linspace(10, 1000, 50)
        szg.create_dataset("r_grid", data=r_grid)
        szg.create_dataset("distribution", data=np.exp(-((r_grid - 100) / 30) ** 2))
        szg.create_dataset("number_dist", data=np.exp(-((r_grid - 100) / 30) ** 2) / r_grid ** 3)
        szg.create_dataset("cumul_vol_dist", data=np.linspace(0, 1, 50))

        # ── data_merge_results (provenance only) ──────────────────────────
        dmg = f.create_group("entry/data_merge_results")
        dmg.attrs["NX_class"] = "NXprocess"
        dmg.attrs["timestamp"] = "2026-05-19T10:20:00"
        dmg.create_dataset("scale", data=1.02)
        dmg.create_dataset("q_shift", data=0.0)
        dmg.create_dataset("background", data=2e-5)
        dmg.create_dataset("chi_squared", data=0.8)
        dmg.create_dataset("q_overlap_min", data=0.05)
        dmg.create_dataset("q_overlap_max", data=0.2)
        dmg.create_dataset("ds1_file", data="usaxs_run.h5")
        dmg.create_dataset("ds2_file", data="saxs_run.h5")

    return fp


@pytest.fixture(scope="session")
def synth_folder(synth_nxcansas_file: Path) -> Path:
    """The folder containing the synthetic file (for list_files / summarize)."""
    return synth_nxcansas_file.parent


@pytest.fixture(scope="session")
def synth_folder_multi(tmp_path_factory, synth_nxcansas_file: Path) -> Path:
    """A folder with a few files of the same sample, varying scan numbers,
    used to exercise tabulate_parameter() across multiple files."""
    out_dir = tmp_path_factory.mktemp("synth_multi")
    import shutil
    # Copy the synth file three times under different scan-number names with
    # slightly different Rg values so the trend is non-trivial.
    for i, rg in enumerate([50.0, 52.0, 55.0], start=10):
        dest = out_dir / f"sample_A_scan_{i}.h5"
        shutil.copy(synth_nxcansas_file, dest)
        with h5py.File(dest, "a") as f:
            # Tweak Rg in level_1 so each file has a different value
            if "entry/unified_fit_results/level_1/Rg" in f:
                del f["entry/unified_fit_results/level_1/Rg"]
                f.create_dataset("entry/unified_fit_results/level_1/Rg",
                                  data=float(rg))
    return out_dir
