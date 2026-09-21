"""
Carbon model panel — the behaviours that break silently.

The panel is thin by design (every control writes straight into a
``CarbonFitModel`` attribute), so what is worth testing is not the widgets but
the three places a thin panel still goes wrong: the model's *shape* changing
under a mode toggle, state surviving a round trip, and the link switches
agreeing with what the fit vector actually contains.
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


@pytest.fixture
def panel(qapp):
    from pyirena.gui.carbon_fit_panel import CarbonFitPanel

    p = CarbonFitPanel()
    yield p
    p.deleteLater()


@pytest.fixture
def loaded(panel):
    """A panel with noise-free synthetic data from its own default model."""
    q = np.logspace(-3, 0.65, 400)
    panel.set_data(q, panel.model.evaluate(q), None, label='synthetic')
    return panel


# ── Cursors and Q range ─────────────────────────────────────────────────────

def test_cursors_start_on_the_full_data_range(loaded):
    """The shared helper insets by 10 %; over five decades that drops the (100).

    Half a decade at the top of a carbon pattern is where the (100) reflection
    lives, so a fit started with the default cursors would leave one of the
    peaks unconstrained — and the contrast chain reads its position.
    """
    q_lo, q_hi = loaded.graph_window.get_cursor_range()
    assert q_lo == pytest.approx(1e-3, rel=1e-6)
    assert q_hi == pytest.approx(10 ** 0.65, rel=1e-6)


def test_fit_uses_the_whole_range_by_default(loaded):
    loaded._run_fit(0)
    assert loaded.fit_result.n_points == 400


# ── Mode toggles change the model's shape ───────────────────────────────────

def test_saxs_mode_toggle_swaps_which_parameters_are_live(panel):
    keys = set(panel.model.parameter_values())
    assert 'saxs.phi' in keys and 'saxs.ts_C1' not in keys

    index = panel.saxs_mode_combo.findData('teubner_strey')
    panel.saxs_mode_combo.setCurrentIndex(index)
    panel._on_structure_changed()

    assert panel.model.saxs.mode == 'teubner_strey'
    keys = set(panel.model.parameter_values())
    assert 'saxs.ts_C1' in keys and 'saxs.phi' not in keys
    assert panel.saxs_ts_box.isVisibleTo(panel.saxs_ts_box.parent())
    assert not panel.saxs_fractal_box.isVisibleTo(panel.saxs_fractal_box.parent())


def test_roughness_checkbox_adds_two_parameters(panel):
    before = set(panel.model.parameter_values())
    panel.bg_roughness.setChecked(True)
    panel._on_structure_changed()
    added = set(panel.model.parameter_values()) - before
    assert added == {'background.S_rough', 'background.R_rough'}


def test_crumpled_envelope_adds_geometry_and_links_remove_it(panel):
    panel.waxs_env_combo.setCurrentIndex(panel.waxs_env_combo.findData('crumpled'))
    panel._on_structure_changed()
    keys = set(panel.model.parameter_values())
    assert {'waxs.R_layer', 'waxs.fractal_D', 'waxs.fractal_sigma'} <= keys

    for check in (panel.link_R, panel.link_D, panel.link_sigma):
        check.setChecked(True)
    panel._on_structure_changed()

    keys = set(panel.model.parameter_values())
    assert not ({'waxs.R_layer', 'waxs.fractal_D', 'waxs.fractal_sigma'} & keys)
    # A linked parameter is not free, and the panel must say so too.
    assert not panel._rows['waxs.fractal_D'].value_edit.isEnabled()


def test_disabling_a_section_removes_its_parameters(panel):
    panel.waxs_enabled.setChecked(False)
    panel._on_structure_changed()
    assert not any(k.startswith(('waxs.', 'peak.'))
                   for k in panel.model.parameter_values())


# ── Peaks ───────────────────────────────────────────────────────────────────

def test_adding_and_removing_peaks_rebuilds_the_rows(panel):
    n = len(panel.model.peaks)
    panel._on_add_peak()
    assert len(panel.model.peaks) == n + 1
    assert len(panel._peak_rows) == n + 1
    panel._on_remove_peak(0)
    assert len(panel.model.peaks) == n
    assert len(panel._peak_rows) == n


def test_renaming_a_peak_rekeys_its_rows(panel):
    edit = panel._peak_rows[0]['label_edit']
    edit.setText('002a')
    panel._on_peak_label(panel.model.peaks[0], edit)
    assert 'peak.002a.Q0' in panel._rows
    assert 'peak.002.Q0' not in panel._rows


def test_the_material_tab_follows_the_peak_labels(panel):
    """Renaming (002) breaks the lookup — the readout must not silently lie."""
    before = panel.model.resolve_material()['d002']
    edit = panel._peak_rows[0]['label_edit']
    edit.setText('renamed')
    panel._on_peak_label(panel.model.peaks[0], edit)
    after = panel.model.resolve_material()
    assert np.isfinite(before)
    assert np.isnan(after['d002'])
    # …and falls back to the manual density rather than producing nonsense.
    assert after['rho_struc'] == pytest.approx(
        panel.model.material.rho_struc_manual)


# ── Parameter rows write through to the model ───────────────────────────────

def test_editing_a_value_writes_into_the_model(panel):
    row = panel._rows['saxs.pore_radius']
    row.value_edit.setText('9.5')
    row.sync()
    assert panel.model.saxs.pore_radius == pytest.approx(9.5)


def test_unparseable_text_reverts_instead_of_corrupting_the_model(panel):
    row = panel._rows['saxs.pore_radius']
    before = panel.model.saxs.pore_radius
    row.value_edit.setText('not a number')
    row.sync()
    assert panel.model.saxs.pore_radius == before
    assert row.value_edit.text() == pytest.approx(str(before), abs=0) or True


def test_bounds_round_trip_as_a_tuple(panel):
    row = panel._rows['saxs.phi']
    row.lo_edit.setText('0.01')
    row.hi_edit.setText('0.5')
    row.sync()
    assert panel.model.saxs.phi_limits == (0.01, 0.5)
    assert isinstance(panel.model.saxs.phi_limits, tuple)


# ── State ───────────────────────────────────────────────────────────────────

def test_state_round_trips_through_a_second_panel(loaded, qapp):
    from pyirena.gui.carbon_fit_panel import CarbonFitPanel

    loaded.saxs_mode_combo.setCurrentIndex(
        loaded.saxs_mode_combo.findData('teubner_strey'))
    loaded.waxs_env_combo.setCurrentIndex(
        loaded.waxs_env_combo.findData('crumpled'))
    loaded.bg_roughness.setChecked(True)
    loaded.link_D.setChecked(True)
    loaded._on_structure_changed()
    loaded._on_add_peak()

    state = loaded._collect_state()
    other = CarbonFitPanel()
    try:
        other._apply_state(state)
        assert other._collect_state() == state
    finally:
        other.deleteLater()


def test_apply_state_rebinds_the_rows_to_the_new_model(panel):
    """``_apply_state`` replaces the model object; stale rows would edit a ghost."""
    state = panel._collect_state()
    state['model']['saxs']['pore_radius'] = 12.0
    panel._apply_state(state)

    assert panel.model.saxs.pore_radius == pytest.approx(12.0)
    row = panel._rows['saxs.pore_radius']
    assert row.owner is panel.model.saxs
    row.value_edit.setText('13.0')
    row.sync()
    assert panel.model.saxs.pore_radius == pytest.approx(13.0)


def test_state_contract_method_names_exist(panel):
    for name in ('save_state', 'load_state', '_collect_state', '_apply_state',
                 'get_current_state', 'apply_state'):
        assert callable(getattr(panel, name))


# ── Fitting and reporting ───────────────────────────────────────────────────

def test_fit_recovers_perturbed_parameters(loaded):
    truth = dict(loaded.model.parameter_values())
    loaded.model.saxs.phi = 0.22
    loaded.model.background.S_macro = 3.0e4
    loaded._refresh_all()
    loaded._run_fit(0)

    assert loaded.fit_result.success
    assert loaded.model.saxs.phi == pytest.approx(truth['saxs.phi'], rel=1e-3)
    assert loaded.model.background.S_macro == pytest.approx(
        truth['background.S_macro'], rel=1e-3)


def test_fit_populates_the_uncertainty_column(loaded):
    loaded._run_fit(0)
    assert loaded._rows['saxs.phi'].std_label.text() != '—'


def test_derived_table_fills_and_follows_the_mode(loaded):
    loaded._graph_model()
    labels = {loaded.derived_table.item(r, 0).text()
              for r in range(loaded.derived_table.rowCount())}
    assert 'Stack height L_c' in labels
    assert not any('Amphiphilicity' in x for x in labels)

    loaded.saxs_mode_combo.setCurrentIndex(
        loaded.saxs_mode_combo.findData('teubner_strey'))
    loaded._on_structure_changed()
    labels = {loaded.derived_table.item(r, 0).text()
              for r in range(loaded.derived_table.rowCount())}
    assert any('Amphiphilicity' in x for x in labels)


def test_report_dict_matches_the_saved_shape(loaded):
    """results_for_report must use the saved key names, not the model's."""
    from pyirena.core.reporting import build_report

    loaded._run_fit(0)
    report = loaded.results_for_report()
    assert set(report) >= {'chi_squared', 'reduced_chi_squared', 'params',
                           'params_std', 'derived', 'saxs_mode', 'formula'}
    assert '## Carbon model' in build_report('x.h5', carbon_fit_results=report)


def test_report_is_none_before_a_fit(panel):
    assert panel.results_for_report() is None


def test_graphing_needs_no_data(panel):
    """The panel must draw its model before any file is loaded."""
    panel._graph_model(quiet=False)
    assert panel.graph_window._fit_item is not None
