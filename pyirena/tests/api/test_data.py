"""Tests for pyirena.api.data."""
from __future__ import annotations

from pyirena.api import read_metadata, read_reduced_data


def test_read_reduced_data(synth_nxcansas_file):
    rd = read_reduced_data(str(synth_nxcansas_file), decimate=50)
    assert rd["found"] is True
    assert rd["n_points"] == 300
    assert rd["q_min"] > 0
    assert rd["I_units"] == "1/cm"
    assert len(rd["Q"]) <= 50
    assert len(rd["I"]) == len(rd["Q"])
    assert rd["decimated"] is True
    assert rd["decimated_from"] == 300


def test_read_reduced_data_full(synth_nxcansas_file):
    rd = read_reduced_data(str(synth_nxcansas_file), include_full=True)
    assert rd["found"] is True
    assert len(rd["Q"]) == 300
    assert rd["decimated"] is False


def test_read_metadata(synth_nxcansas_file):
    meta = read_metadata(str(synth_nxcansas_file))
    assert meta["found"] is True
    assert meta["sample_name"] == "sample_A"
    assert meta["thickness"] == 0.1
    assert meta["blank"] == "blank.h5"
    assert meta["label"] == "sample_A_scan_42"
    assert meta["timestamp"] is not None
