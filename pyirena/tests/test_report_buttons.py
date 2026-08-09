"""
Tests for the fit panels' "Copy results" / "Save report…" buttons.

Skipped without Qt; the `test-gui` CI job runs them offscreen.

What matters here is the join between a live panel and the shared formatter:
each panel's ``results_for_report()`` must produce the dict shape
:mod:`pyirena.core.reporting` consumes, with the key names the *saved* results
use — the fit objects use different names internally (``chi2`` vs
``chi_squared``), and a mismatch shows up as "N/A" in the user's report rather
than as an error.
"""

import numpy as np
import pytest

pytest.importorskip("pyirena.gui._qt", reason="Qt binding not installed")
pytest.importorskip("pyqtgraph")

from pyirena.gui._qt import QApplication  # noqa: E402
from pyirena.gui.report_buttons import (  # noqa: E402
    build_panel_report,
    copy_report,
    make_report_buttons,
    save_report,
)

_APP = None


@pytest.fixture(scope="session")
def qapp():
    global _APP
    _APP = QApplication.instance() or QApplication([])
    yield _APP


@pytest.fixture
def unified_panel(qapp):
    from pyirena.gui.unified_fit import UnifiedFitPanel

    panel = UnifiedFitPanel()
    q = np.logspace(-3, -1, 40)
    intensity = 1e3 * q ** -3
    panel.data = {
        "Q": q, "Intensity": intensity, "Error": 0.05 * intensity,
        "filepath": "/data/AeroGel_500C.h5", "label": "AeroGel_500C.h5",
        "is_nxcansas": True,
    }
    yield panel
    panel.close()


# ── The panel → formatter join ───────────────────────────────────────────

def test_unified_panel_produces_a_report_dict(unified_panel):
    results = unified_panel.results_for_report()
    for key in ("Q", "intensity_data", "residuals", "num_levels",
                "background", "chi_squared", "levels"):
        assert key in results, f"{key} missing from the report dict"
    assert len(results["levels"]) == results["num_levels"]
    level = results["levels"][0]
    # Report-shaped level keys, not the model's internal names (RgCO/correlations)
    for key in ("G", "Rg", "B", "P", "RgCutoff", "correlated"):
        assert key in level, f"{key} missing from level dict"


def test_unified_panel_report_text_has_real_numbers(unified_panel):
    text = build_panel_report(
        unified_panel.results_for_report(), "fit_results",
        "/data/AeroGel_500C.h5", unified_panel._data_info_for_report(),
    )
    assert "# pyIrena Report" in text
    assert "| **File** | `AeroGel_500C.h5` |" in text
    assert "## Data Summary" in text
    assert "## Level 1 Parameters" in text
    assert "N/A" not in text.split("## Level 1")[1]


def test_report_dict_matches_what_the_panel_saves(unified_panel):
    """The report and the HDF5 save must describe the same levels."""
    report_levels = unified_panel.results_for_report()["levels"]
    saved_levels = unified_panel._collect_level_dicts()
    assert report_levels == saved_levels


# ── Clipboard and file ───────────────────────────────────────────────────

def test_copy_report_puts_markdown_on_the_clipboard(qapp):
    assert copy_report(None, "# report\n\nbody") is True
    assert QApplication.clipboard().text() == "# report\n\nbody"


def test_copy_report_refuses_empty_text(qapp, monkeypatch):
    from pyirena.gui import report_buttons as rb

    shown = []
    monkeypatch.setattr(rb.QMessageBox, "information",
                        staticmethod(lambda *a, **k: shown.append(a)))
    assert copy_report(None, "") is False
    assert shown, "the user must be told why nothing was copied"


def test_save_report_writes_markdown(qapp, tmp_path, monkeypatch):
    from pyirena.gui import report_buttons as rb

    target = tmp_path / "fit_report.md"
    monkeypatch.setattr(rb.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(target), "")))
    monkeypatch.setattr(rb, "save_report", rb.save_report)   # keep real function
    written = save_report(None, "# report\n", "unified_fit")
    assert written == str(target)
    assert target.read_text() == "# report\n"


def test_save_report_appends_the_extension(qapp, tmp_path, monkeypatch):
    from pyirena.gui import report_buttons as rb

    monkeypatch.setattr(rb.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(tmp_path / "r"), "")))
    written = save_report(None, "# report\n", "unified_fit")
    assert written.endswith(".md")


def test_build_panel_report_returns_empty_string_without_results():
    assert build_panel_report(None, "fit_results") == ""
    assert build_panel_report({}, "fit_results") == ""


# ── The button row ───────────────────────────────────────────────────────

def test_make_report_buttons_builds_two_wired_buttons(qapp):
    row = make_report_buttons(None, lambda: None, tool_key="fit_results")
    assert row.copy_button.text() == "Copy results"
    assert row.save_button.text() == "Save report…"
    assert row.copy_button.toolTip() and row.save_button.toolTip()


def test_button_report_reflects_the_current_results(qapp):
    """The provider is called on click, not when the buttons are built."""
    state = {"results": None}
    row = make_report_buttons(None, lambda: state["results"],
                              tool_key="simple_fit_results")
    assert row.report_text() == ""
    state["results"] = {"model": "Guinier", "params": {"Rg": 150.0}}
    text = row.report_text()
    assert "## Simple Fits" in text and "Guinier" in text


def test_panels_expose_the_report_hook(qapp):
    """Every fit panel must keep the method the buttons call."""
    from pyirena.gui.modeling_panel import ModelingPanel
    from pyirena.gui.simple_fits_panel import SimpleFitsPanel
    from pyirena.gui.sizes_panel import SizesFitPanel
    from pyirena.gui.unified_fit import UnifiedFitPanel
    from pyirena.gui.waxs_peakfit_panel import WAXSPeakFitPanel

    for cls in (UnifiedFitPanel, SizesFitPanel, SimpleFitsPanel,
                ModelingPanel, WAXSPeakFitPanel):
        assert callable(getattr(cls, "results_for_report", None)), cls.__name__
        assert callable(getattr(cls, "_data_info_for_report", None)), cls.__name__


def test_panels_report_nothing_before_a_fit(qapp):
    """Copy results on a fresh panel must offer no report, not crash.

    WAXS is excluded on purpose: its report describes the peak rows on screen,
    which exist before any fit is run.
    """
    from pyirena.gui.modeling_panel import ModelingPanel
    from pyirena.gui.simple_fits_panel import SimpleFitsPanel
    from pyirena.gui.sizes_panel import SizesFitPanel

    for cls in (SizesFitPanel, SimpleFitsPanel, ModelingPanel):
        panel = cls()
        try:
            assert panel.results_for_report() is None, cls.__name__
            assert panel._data_info_for_report() is None, cls.__name__
        finally:
            panel.close()


def test_waxs_panel_reports_the_peaks_on_screen(qapp):
    """The report follows the peak rows, not the last fit's snapshot.

    Peak counts are read from the panel rather than assumed: a panel restores
    the peaks of the previous session from the state file, so a fixed number
    would make this test depend on the developer's machine.
    """
    import numpy as np

    from pyirena.gui.waxs_peakfit_panel import WAXSPeakFitPanel

    panel = WAXSPeakFitPanel()
    try:
        q = np.linspace(1.0, 5.0, 200)
        intensity = 100 + 500 * np.exp(-((q - 2.5) / 0.05) ** 2)
        panel._q, panel._I, panel._dI = q, intensity, np.sqrt(intensity)
        before = len(panel._get_peaks())
        panel._add_peak_row({"shape": "Gauss", "A": 500.0, "Q0": 2.5, "FWHM": 0.1})
        expected = before + 1

        results = panel.results_for_report()
        assert results["n_peaks"] == expected
        # Padded with empty dicts for peaks the last fit never saw.
        assert len(results["peaks_std"]) == expected

        text = build_panel_report(results, "waxs_peakfit_results", "waxs.h5")
        assert f"| Number of peaks | {expected} |" in text
        assert f"### Peak {expected} (Gauss)" in text
    finally:
        panel.close()


def test_modeling_report_dict_comes_from_the_core_converter(qapp):
    """The panel must not hand-roll the population flattening."""
    import numpy as np

    from pyirena.core.modeling import (
        ModelingConfig,
        ModelingEngine,
        SizeDistPopulation,
        result_to_report_dict,
    )

    q = np.logspace(np.log10(0.01), np.log10(0.3), 40)
    pop = SizeDistPopulation()
    pop.dist_params_fit = {k: False for k in pop.dist_params}
    pop.ff_params_fit = {k: False for k in pop.ff_params}
    cfg = ModelingConfig(populations=[pop], background=0.001,
                         fit_background=False, q_min=q.min(), q_max=q.max())
    engine = ModelingEngine()
    intensity = engine.total_intensity(cfg, q, use_cache=False)[0]
    result = engine.fit(cfg, q, intensity, 0.02 * intensity + 1e-9)

    as_dict = result_to_report_dict(result)
    assert as_dict["background"] == pytest.approx(0.001)
    assert len(as_dict["populations"]) == 1
    # Flattened with asdict → every population field is present, so a new
    # parameter shows up in reports without touching this converter.
    assert "dist_params" in as_dict["populations"][0]
    assert as_dict["populations"][0]["population_index"] == 1

    text = build_panel_report(as_dict, "modeling_results", "s.h5")
    assert "## Modeling Results" in text
    assert "### Population 1" in text
