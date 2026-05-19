"""Tests for pyirena.api.aggregate."""
from __future__ import annotations

from pyirena.api import summarize_sample, tabulate_parameter


def test_tabulate_parameter_across_multi(synth_folder_multi):
    tab = tabulate_parameter(
        folder=str(synth_folder_multi),
        tool="unified_fit",
        parameter="Rg",
        subgroup_index=1,
        x_axis="scan_number",
    )
    assert tab["tool"] == "unified_fit"
    assert tab["parameter"] == "Rg"
    assert tab["n_rows"] == 3
    # Rg values should be {50, 52, 55}, ordered by scan number 10,11,12
    rg_values = [r["value"] for r in tab["rows"]]
    assert rg_values == [50.0, 52.0, 55.0]
    scan_numbers = [r["scan_number"] for r in tab["rows"]]
    assert scan_numbers == [10, 11, 12]
    assert tab["units"] == "Å"


def test_tabulate_parameter_top_level_scalar(synth_folder_multi):
    tab = tabulate_parameter(
        folder=str(synth_folder_multi),
        tool="unified_fit",
        parameter="chi_squared",
    )
    assert tab["n_rows"] == 3
    for row in tab["rows"]:
        assert row["value"] == 2.1


def test_tabulate_parameter_unknown_tool_raises(synth_folder_multi):
    import pytest
    with pytest.raises(ValueError):
        tabulate_parameter(folder=str(synth_folder_multi),
                            tool="bogus_tool", parameter="x")


def test_summarize_sample(synth_folder_multi):
    s = summarize_sample(folder=str(synth_folder_multi), sample="sample_A")
    assert s["n_files"] == 3
    assert s["analyses_count"]["unified_fit"] == 3
    assert "unified_fit" in s["parameter_ranges"]
    chi = s["parameter_ranges"]["unified_fit"]["chi_squared"]
    assert chi["n"] == 3
    assert chi["min"] == chi["max"] == 2.1
