"""
Rules about which "Fit?" boxes may be checked, and which bounds are shown.

Two user-reported confusions, both about parameters that cannot meaningfully
be fitted:

* Modeling → size-distribution population: ``Contrast`` and ``Scale`` enter the
  model only as their product, so freeing both leaves the least-squares problem
  with a flat direction — the solver wanders and the covariance is singular.
  At most one may be fitted, and the box the user just ticked wins.
* Simple Fits → Invariant: the Invariant is a direct calculation with no
  least-squares step, so its ``Contrast`` is never fitted.  Showing lo/hi
  bound fields next to it only prompted users to ask what they were for.
"""

from __future__ import annotations

import numpy as np
import pytest


def _require_qt() -> None:
    try:
        import pyirena.gui._qt  # noqa: F401
    except ImportError:
        pytest.skip("Qt (PySide6/PyQt6) not available", allow_module_level=True)


_require_qt()
pytest.importorskip("pyqtgraph")

from pyirena.gui._qt import QApplication  # noqa: E402

_APP = None


@pytest.fixture(scope="session")
def qapp():
    global _APP
    _APP = QApplication.instance() or QApplication([])
    yield _APP


# ── Modeling: Contrast / Scale are mutually exclusive ───────────────────────

@pytest.fixture
def population(qapp):
    from pyirena.gui.modeling_panel import PopulationTab

    return PopulationTab(0)


def _flags(pop):
    return pop.contrast_fit_cb.isChecked(), pop.scale_fit_cb.isChecked()


def test_default_fits_scale_only(population):
    assert _flags(population) == (False, True)


def test_checking_contrast_unchecks_scale(population):
    population.contrast_fit_cb.setChecked(True)
    assert _flags(population) == (True, False)


def test_checking_scale_unchecks_contrast(population):
    population.contrast_fit_cb.setChecked(True)
    population.scale_fit_cb.setChecked(True)
    assert _flags(population) == (False, True)


def test_last_box_ticked_wins_repeatedly(population):
    for _ in range(3):
        population.contrast_fit_cb.setChecked(True)
        assert _flags(population) == (True, False)
        population.scale_fit_cb.setChecked(True)
        assert _flags(population) == (False, True)


def test_neither_fitted_is_allowed(population):
    """Holding both fixed is a legitimate choice — only *both fitted* is not."""
    population.scale_fit_cb.setChecked(False)
    assert _flags(population) == (False, False)
    population.contrast_fit_cb.setChecked(True)
    population.contrast_fit_cb.setChecked(False)
    assert _flags(population) == (False, False)


def test_loaded_state_with_both_checked_is_normalised(population):
    """Setups written before this rule can carry both flags true."""
    population._building = True          # emulate the loader's guard
    population.contrast_fit_cb.setChecked(True)
    population.scale_fit_cb.setChecked(True)
    population._building = False

    population._enforce_contrast_scale_exclusive()
    assert _flags(population) == (False, True)


# ── Simple Fits: the Invariant's Contrast has no bounds ─────────────────────

@pytest.fixture
def simple_panel(qapp):
    from pyirena.gui.simple_fits_panel import SimpleFitsPanel

    panel = SimpleFitsPanel()
    panel.show()          # visibility flags are only meaningful once shown
    q = np.logspace(-3, 0, 200)
    intensity = 1e3 * q ** -3 + 0.1
    panel.set_data(q, intensity, 0.02 * intensity, label="t")
    yield panel
    panel.hide()


def test_fit_model_shows_bounds_for_every_parameter(simple_panel):
    simple_panel.model_combo.setCurrentText("Sphere")
    assert simple_panel._param_lo_edits, "Sphere should expose parameters"
    for name, edit in simple_panel._param_lo_edits.items():
        assert edit.isVisible(), f"{name} lost its lower bound field"
        assert simple_panel._param_hi_edits[name].isVisible()


def test_invariant_contrast_has_no_fit_box_and_no_bounds(simple_panel):
    simple_panel.complex_bg_check.setChecked(False)
    simple_panel.model_combo.setCurrentText("Invariant")

    assert not simple_panel._param_fit_checks["Contrast"].isVisible()
    assert not simple_panel._param_lo_edits["Contrast"].isVisible()
    assert not simple_panel._param_hi_edits["Contrast"].isVisible()


def test_invariant_without_bg_hides_the_bound_column_headers(simple_panel):
    simple_panel.complex_bg_check.setChecked(False)
    simple_panel.model_combo.setCurrentText("Invariant")

    assert not simple_panel._lo_header_lbl.isVisible()
    assert not simple_panel._hi_header_lbl.isVisible()


def test_invariant_keeps_bounds_for_the_background_terms(simple_panel):
    """The BG_* terms *are* fitted by the background prefit — keep their bounds."""
    simple_panel.model_combo.setCurrentText("Invariant")
    simple_panel.complex_bg_check.setChecked(True)

    assert not simple_panel._param_lo_edits["Contrast"].isVisible()
    for name in ("BG_B", "BG_P", "BG_flat"):
        assert simple_panel._param_lo_edits[name].isVisible(), f"{name} lost its bounds"
        assert simple_panel._param_hi_edits[name].isVisible()
    assert simple_panel._lo_header_lbl.isVisible()
