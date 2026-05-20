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


def test_tabulate_parameter_runtime_simple_fits_param(synth_folder_multi):
    """Regression: simple_fits has model-specific param names (Kp, Lc, ...)
    not enumerated in TOOL_REGISTRY. _extract_scalar must fall back to
    probing params/<leaf> directly. Reported by the AI consumer."""
    # 'Kp' is in the fixture's params/ but NOT in the static scalar list
    tab = tabulate_parameter(
        folder=str(synth_folder_multi),
        tool="simple_fits",
        parameter="Kp",
    )
    assert tab["n_rows"] == 3
    for row in tab["rows"]:
        assert row["value"] == 3.7e-5
        assert row["stddev"] == 2.1e-6


def test_tabulate_parameter_runtime_param_with_prefix(synth_folder_multi):
    """The fallback strips an optional 'param_' prefix, so 'param_Kp'
    resolves to the same params/Kp as bare 'Kp'."""
    tab = tabulate_parameter(
        folder=str(synth_folder_multi),
        tool="simple_fits",
        parameter="param_Kp",
    )
    assert tab["n_rows"] == 3
    for row in tab["rows"]:
        assert row["value"] == 3.7e-5


def test_tabulate_parameter_truly_missing_returns_null(synth_folder_multi):
    """A name that exists in NO schema spec and NO params/ entry still
    returns rows with value=None (it does not raise)."""
    tab = tabulate_parameter(
        folder=str(synth_folder_multi),
        tool="simple_fits",
        parameter="DoesNotExistAnywhere",
    )
    assert tab["n_rows"] == 3
    for row in tab["rows"]:
        assert row["value"] is None
        assert row["stddev"] is None


def test_summarize_sample(synth_folder_multi):
    s = summarize_sample(folder=str(synth_folder_multi), sample="sample_A")
    assert s["n_files"] == 3
    assert s["analyses_count"]["unified_fit"] == 3
    assert "unified_fit" in s["parameter_ranges"]
    chi = s["parameter_ranges"]["unified_fit"]["chi_squared"]
    assert chi["n"] == 3
    assert chi["min"] == chi["max"] == 2.1
