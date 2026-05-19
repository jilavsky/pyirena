"""Tests for pyirena.api.results — per-tool readers."""
from __future__ import annotations

from pyirena.api import (
    read_manipulation_provenance,
    read_merge_provenance,
    read_simple_fit,
    read_size_distribution,
    read_unified_fit,
)


def test_read_simple_fit(synth_nxcansas_file):
    r = read_simple_fit(str(synth_nxcansas_file))
    assert r["found"] is True
    assert r["tool"] == "simple_fits"
    assert r["model"] == "Guinier"
    assert r["chi_squared"] == 1.23
    assert r["params"]["Rg"] == 50.0
    assert r["params"]["I0"] == 0.01
    assert r["params_std"]["Rg"] == 1.2
    # Arrays omitted by default
    assert r["Q"] is None
    assert r["I_model"] is None


def test_read_simple_fit_with_arrays(synth_nxcansas_file):
    r = read_simple_fit(str(synth_nxcansas_file), include_arrays=True,
                        max_points=20)
    assert r["Q"] is not None
    assert len(r["Q"]) <= 20


def test_read_unified_fit(synth_nxcansas_file):
    r = read_unified_fit(str(synth_nxcansas_file))
    assert r["found"] is True
    assert r["num_levels"] == 1
    assert r["background"] == 1e-4
    assert r["chi_squared"] == 2.1
    assert len(r["levels"]) == 1
    lv = r["levels"][0]
    assert lv["Rg"] == 52.0
    assert lv["G"] == 0.012
    assert lv["P"] == 4.0


def test_read_size_distribution(synth_nxcansas_file):
    r = read_size_distribution(str(synth_nxcansas_file))
    assert r["found"] is True
    assert r["volume_fraction"] == 0.05
    assert r["rg"] == 55.0
    assert r["shape"] == "sphere"
    assert r["method"] == "maxent"
    assert r["n_bins"] == 50
    # Arrays omitted by default
    assert r["r_grid"] is None


def test_read_size_distribution_with_arrays(synth_nxcansas_file):
    r = read_size_distribution(str(synth_nxcansas_file), include_arrays=True,
                                max_points=10)
    assert r["r_grid"] is not None
    assert len(r["r_grid"]) <= 10
    assert r["distribution"] is not None


def test_read_merge_provenance(synth_nxcansas_file):
    r = read_merge_provenance(str(synth_nxcansas_file))
    assert r["found"] is True
    assert r["scale"] == 1.02
    assert r["q_shift"] == 0.0
    assert r["chi_squared"] == 0.8
    assert r["ds1_file"] == "usaxs_run.h5"
    assert r["ds2_file"] == "saxs_run.h5"


def test_read_manipulation_provenance_absent(synth_nxcansas_file):
    """File has no manipulation group — should return found=False."""
    r = read_manipulation_provenance(str(synth_nxcansas_file))
    assert r["found"] is False
    assert r["operation"] is None
