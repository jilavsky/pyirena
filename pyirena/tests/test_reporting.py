"""
Tests for :mod:`pyirena.core.reporting` — the one Markdown report builder.

Needs no Qt: the builder is pure numpy + stdlib, which is the point of moving
it out of ``gui/`` (``api/`` and ``core/`` may not import from ``gui/``).

The report is user-facing text that ends up in logbooks and manuscripts, so
these tests pin what it must contain and, just as importantly, that a section
is *absent* when its tool was not run — a report claiming a Size Distribution
that was never fitted is worse than no report.
"""

import pytest

from pyirena.core.reporting import _build_report, build_report
from pyirena.tests import _report_fixture as fx


def test_public_and_private_names_are_the_same_builder():
    """`_build_report` is kept because the Data Selector imported it for years."""
    args = fx.all_inputs()
    assert fx.strip_timestamp(build_report(**args)) == fx.strip_timestamp(
        _build_report(**args)
    )


def test_report_is_reachable_from_the_old_gui_location():
    """The Data Selector's import path must keep working after the move."""
    pytest.importorskip("pyirena.gui._qt", reason="Qt binding not installed")
    from pyirena.gui.data_selector.report import _build_report as gui_build

    assert gui_build is _build_report


# ── Header and data summary ──────────────────────────────────────────────

def test_header_shows_the_file_name_only():
    text = build_report("/data/beamline/AeroGel_500C.h5")
    assert "| **File** | `AeroGel_500C.h5` |" in text
    assert "/data/beamline" not in text        # no leaking of local paths


def test_data_summary_reports_ranges_and_point_count():
    text = build_report("s.h5", data_info=fx.data_info())
    assert "## Data Summary" in text
    assert "| Data points | 50 |" in text
    assert "Q range" in text and "Å⁻¹" in text


# ── Per-tool sections ────────────────────────────────────────────────────

def test_unified_section_lists_levels_and_uncertainties():
    text = build_report("s.h5", fit_results=fx.unified_results())
    assert "## Fit Quality" in text
    assert "| Number of levels | 2 |" in text
    assert "## Level 1 Parameters" in text
    assert "## Level 2 Parameters" in text
    # Monte-Carlo errors present on level 1 → uncertainty column appears
    assert "Uncertainty (1σ)" in text
    # Correlated level 2 reports its packing parameters
    assert "| ETA |" in text and "| PACK |" in text


def test_unified_section_reports_robust_quality_metrics():
    text = build_report("s.h5", fit_results=fx.unified_results())
    assert "σ-scale (robust)" in text
    assert "Longest same-sign run" in text


def test_sizes_section_reports_method_and_grid():
    text = build_report("s.h5", sizes_results=fx.sizes_results())
    assert "## Size Distribution" in text
    assert "| Method | MaxEnt |" in text
    assert "| Particle shape | Spheroid |" in text
    assert "**Grid / model setup:**" in text
    assert "**MaxEnt parameters:**" in text        # method-specific block


def test_simple_fits_section_reports_parameters_with_uncertainties():
    text = build_report("s.h5", simple_fit_results=fx.simple_fit_results())
    assert "## Simple Fits" in text
    assert "| Model | Guinier |" in text
    assert "Uncertainty (1σ)" in text
    assert "| Rg |" in text


def test_waxs_section_reports_peaks_in_both_parameter_forms():
    """Peak parameters arrive as {'value': …} dicts live, floats from HDF5."""
    text = build_report("s.h5", waxs_peakfit_results=fx.waxs_results())
    assert "## WAXS Peak Fit" in text
    assert "| Number of peaks | 2 |" in text
    assert "### Peak 1 (Gauss)" in text
    assert "### Peak 2 (Pseudo-Voigt)" in text
    assert "| A | 500 |" in text.replace(" | ± 12 |", " |")   # dict form read
    assert "| Q0 Å⁻¹ | 3.9 |" in text or "3.9" in text        # float form read
    assert "Area (∫ peak dq)" in text


def test_waxs_section_survives_a_fit_that_has_not_been_run():
    """Peaks placed by hand but not yet fitted must still report."""
    partial = {"n_peaks": 1, "bg_shape": "SNIP", "peaks": [
        {"shape": "Gauss", "A": 10.0, "Q0": 2.0, "FWHM": 0.1}]}
    text = build_report("s.h5", waxs_peakfit_results=partial)
    assert "| Chi-squared (χ²) | N/A |" in text
    assert "| DOF | N/A |" in text          # None must not print as "None"


def test_unknown_peak_shape_costs_one_cell_not_the_whole_report():
    """A file written by a different pyIrena build must still report."""
    odd = {"n_peaks": 1, "bg_shape": "SNIP", "peaks": [
        {"shape": "SomeFutureShape", "A": 10.0, "Q0": 2.0, "FWHM": 0.1}]}
    text = build_report("s.h5", waxs_peakfit_results=odd)
    assert "### Peak 1 (SomeFutureShape)" in text
    assert "| Area (∫ peak dq) | N/A | — |" in text


def test_uncomputed_peak_area_uncertainty_is_a_dash_not_zero():
    partial = {"n_peaks": 1, "bg_shape": "SNIP", "peaks": [
        {"shape": "Gauss", "A": 10.0, "Q0": 2.0, "FWHM": 0.1}]}
    area_row = [ln for ln in build_report("s.h5", waxs_peakfit_results=partial)
                .splitlines() if "Area" in ln][0]
    assert area_row.endswith("| — |"), area_row


def test_modeling_section_lists_only_enabled_populations():
    text = build_report("s.h5", modeling_results=fx.modeling_results())
    assert "## Modeling Results" in text
    assert "| Populations (enabled) | 2 |" in text
    assert "### Population 1 — matrix pores (size_dist)" in text
    assert "### Population 2 — large scale (unified_level)" in text
    assert "disabled" not in text


def test_modeling_section_reports_type_specific_parameters():
    text = build_report("s.h5", modeling_results=fx.modeling_results())
    # size_dist population: distribution, form factor, derived quantities
    assert "| Dist: mean_size | 80 |" in text
    assert "| Form factor | sphere |" in text
    assert "| Vol. fraction |" in text
    # unified_level population: Beaucage parameters
    assert "| G [cm⁻¹] |" in text and "| Rg [Å] |" in text


def test_sections_are_omitted_when_their_tool_was_not_run():
    text = build_report("s.h5", fit_results=fx.unified_results())
    assert "## Size Distribution" not in text
    assert "## Simple Fits" not in text
    assert "## Data Summary" not in text


def test_full_report_contains_every_requested_section():
    text = build_report(**fx.all_inputs())
    for heading in ("## Data Summary", "## Fit Quality", "## Level 1 Parameters",
                    "## Size Distribution", "## Simple Fits"):
        assert heading in text, f"{heading} missing"
    assert text.rstrip().endswith("*Generated by pyIrena*")


def test_missing_values_render_as_na_not_a_traceback():
    """Half-filled results are normal — an aborted fit must still report."""
    sparse = {"model": "Porod", "params": {"B": 1e-5}}
    text = build_report("s.h5", simple_fit_results=sparse)
    assert "## Simple Fits" in text
    assert "N/A" in text


def test_report_is_stable_across_calls():
    args = fx.all_inputs()
    assert fx.strip_timestamp(build_report(**args)) == fx.strip_timestamp(
        build_report(**args)
    )
