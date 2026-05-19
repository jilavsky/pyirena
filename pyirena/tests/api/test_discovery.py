"""Tests for pyirena.api.discovery."""
from __future__ import annotations

from pyirena.api import inspect_file, list_files, summarize_folder


def test_list_files_returns_entries(synth_folder):
    rows = list_files(str(synth_folder))
    assert isinstance(rows, list)
    assert len(rows) == 1
    r = rows[0]
    assert r["name"].endswith(".h5")
    assert r["sample"] == "sample_A"
    assert r["scan_number"] == 42
    assert set(r["analyses"]) >= {
        "simple_fits", "unified_fit", "size_distribution", "data_merge",
    }
    assert r["size_bytes"] > 0
    assert "mtime" in r


def test_list_files_shallow(synth_folder):
    rows = list_files(str(synth_folder), deep=False)
    assert rows
    # Shallow: analyses empty, sample None
    assert rows[0]["analyses"] == []
    assert rows[0]["sample"] is None


def test_summarize_folder(synth_folder):
    s = summarize_folder(str(synth_folder))
    assert s["n_files"] == 1
    assert "sample_A" in s["samples"]
    assert s["analyses_count"]["unified_fit"] == 1
    assert s["analyses_count"]["simple_fits"] == 1
    assert s["analyses_count"]["size_distribution"] == 1
    assert s["mtime_min"] is not None


def test_summarize_folder_sample_filter(synth_folder):
    s = summarize_folder(str(synth_folder), sample_filter="A")
    assert s["n_files"] == 1
    s = summarize_folder(str(synth_folder), sample_filter="nonexistent")
    assert s["n_files"] == 0


def test_inspect_file(synth_nxcansas_file):
    info = inspect_file(str(synth_nxcansas_file))
    assert info["sample"] == "sample_A"
    assert info["scan_number"] == 42
    assert info["reduced_data_present"] is True
    assert info["n_points"] == 300
    assert info["q_min"] > 0
    assert info["q_max"] > info["q_min"]
    assert "unified_fit" in info["analyses_present"]
