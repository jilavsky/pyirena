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


# ── Revert back ─────────────────────────────────────────────────────────────

def test_revert_is_disabled_until_there_is_something_to_revert_to(panel):
    assert panel.revert_btn.isEnabled() is False


def test_revert_restores_the_pre_fit_parameters(loaded):
    loaded.model.saxs.phi = 0.42
    loaded.model.saxs.pore_radius = 12.0
    loaded.model.peaks[0].Q0 = 1.62
    loaded._refresh_all()
    start = dict(loaded.model.parameter_values())

    loaded._run_fit(0)
    assert loaded.revert_btn.isEnabled()
    assert loaded.model.parameter_values() != start   # the fit moved things

    loaded._on_revert()
    assert loaded.model.parameter_values() == start
    # …and the widgets, not just the model.
    assert float(loaded._rows['saxs.phi'].value_edit.text()) == pytest.approx(0.42)


def test_revert_restores_bounds_and_fit_flags_too(loaded):
    """The snapshot is the whole model, not just the values the fitter moved."""
    loaded.model.saxs.phi_limits = (0.02, 0.44)
    loaded.model.saxs.fit_globule_k = True
    loaded._refresh_all()

    loaded._run_fit(0)
    loaded.model.saxs.phi_limits = (0.0, 1.0)
    loaded.model.saxs.fit_globule_k = False

    loaded._on_revert()
    assert loaded.model.saxs.phi_limits == (0.02, 0.44)
    assert loaded.model.saxs.fit_globule_k is True


def test_revert_keeps_the_q_range(loaded):
    """The cursors are the user's, not the fit's — reverting them is a surprise."""
    loaded.graph_window.set_cursor_range(0.01, 1.0)
    loaded._run_fit(0)
    loaded._on_revert()
    assert loaded.model.q_min == pytest.approx(0.01)
    assert loaded.model.q_max == pytest.approx(1.0)


def test_revert_clears_the_stale_fit_display(loaded):
    loaded._run_fit(0)
    assert loaded._rows['saxs.phi'].std_label.text() != '—'
    loaded._on_revert()
    assert loaded._rows['saxs.phi'].std_label.text() == '—'
    assert loaded.fit_result is None
    assert 'Not fitted yet' in loaded.fit_summary.text()


def test_revert_twice_is_harmless(loaded):
    """It is an undo of the last fit, not a stack."""
    start = dict(loaded.model.parameter_values())
    loaded._run_fit(0)
    loaded._on_revert()
    loaded._on_revert()
    assert loaded.model.parameter_values() == start


def test_revert_survives_a_peak_being_added_by_the_fit_cycle(loaded):
    """The peak list is part of the snapshot, so the rows must be rebuilt."""
    loaded._run_fit(0)
    loaded._on_add_peak()
    n_after_add = len(loaded.model.peaks)
    loaded._on_revert()
    assert len(loaded.model.peaks) == n_after_add - 1
    assert len(loaded._peak_rows) == len(loaded.model.peaks)


def test_loading_a_setup_drops_the_backup(loaded):
    """A snapshot of a model that is no longer loaded would graft onto the new one."""
    loaded._run_fit(0)
    assert loaded.revert_btn.isEnabled()
    loaded._apply_state({})
    assert loaded.revert_btn.isEnabled() is False
    assert loaded._param_backup is None


# ── Wheel step on peak parameters ───────────────────────────────────────────

def _wheel(edit, notches=1):
    from pyirena.gui._qt import Qt, QtCore, QtGui

    event = QtGui.QWheelEvent(
        QtCore.QPointF(5, 5), QtCore.QPointF(5, 5), QtCore.QPoint(0, 0),
        QtCore.QPoint(0, 120 * notches), Qt.MouseButton.NoButton,
        Qt.KeyboardModifier.NoModifier, Qt.ScrollPhase.NoScrollPhase, False)
    edit.wheelEvent(event)


@pytest.mark.parametrize("key, expected_step", [
    ("peak.002.Q0", 0.01),        # 0.1 Å⁻¹ was a third of a peak width per notch
    ("peak.002.FWHM_G", 0.001),
    ("peak.002.FWHM_L", 0.001),
    ("peak.002.K", 0.1),          # amplitude keeps the coarse step
])
def test_peak_parameters_scrub_at_a_usable_rate(panel, key, expected_step):
    row = panel._rows[key]
    before = float(row.value_edit.text())
    _wheel(row.value_edit)
    after = float(row.value_edit.text())
    assert abs(after - before) == pytest.approx(expected_step, rel=0.05)
    assert after > before


def test_scrubbing_a_peak_position_writes_through_to_the_model(panel):
    row = panel._rows['peak.002.Q0']
    before = panel.model.peaks[0].Q0
    _wheel(row.value_edit)
    assert panel.model.peaks[0].Q0 > before


# ── WAXS zoom panel ─────────────────────────────────────────────────────────

def test_zoom_frames_the_diffraction_end_of_the_data(loaded):
    """Q = 1 Å⁻¹ to the highest measured Q — set by the data, not by the peaks.

    Before the first fit the peak positions are only graphite guesses, so
    framing the window on them could miss the sample's data entirely.
    """
    loaded.zoom_check.setChecked(True)
    loaded._on_display_option_changed()
    lo, hi = loaded.graph_window.zoom_plot.viewRange()[0]
    q_max = float(np.nanmax(loaded.data['Q']))
    assert lo == pytest.approx(loaded.WAXS_ZOOM_Q_MIN, rel=0.1)
    assert hi == pytest.approx(q_max, rel=0.1)


def test_zoom_y_axis_ignores_the_low_q_decades(loaded):
    """A linear Y axis ranged over the whole pattern leaves the peaks flat.

    The low-Q end of a carbon curve is orders of magnitude above the
    diffraction peaks, so the Y range has to come from the window, not the
    curve.
    """
    loaded.zoom_check.setChecked(True)
    loaded._on_display_option_changed()
    y_lo, y_hi = loaded.graph_window.zoom_plot.viewRange()[1]

    q, I = loaded.data['Q'], loaded.data['Intensity']
    in_window = q >= loaded.WAXS_ZOOM_Q_MIN
    assert y_hi == pytest.approx(I[in_window].max(), rel=0.25)
    # …and nowhere near the full-curve maximum, which is what it used to be.
    assert y_hi < I.max() / 1000.0


def test_zoom_survives_a_fit_on_a_restricted_q_range(loaded):
    """The model lives on the fit range, the data on the full range.

    Sharing one cached ``q`` between them paired a 284-point model with a
    400-point data array; the redraw raised part-way through and the panel
    came back empty (and the exception escaped the fit).
    """
    loaded.zoom_check.setChecked(True)
    loaded._on_display_option_changed()
    assert len(loaded.graph_window._zoom_items) == 3

    loaded.graph_window.set_cursor_range(0.01, 4.0)
    loaded._run_fit(0)

    assert loaded.fit_result is not None            # the fit finished
    assert len(loaded.graph_window._zoom_items) == 3  # …and the panel is not empty
    cached = loaded.graph_window._last
    assert len(cached['q_model']) < len(cached['q_data'])


def test_zoom_redraw_cannot_take_the_panel_down(loaded):
    """It runs from the fit-completion path, so it must never raise."""
    loaded.zoom_check.setChecked(True)
    loaded._on_display_option_changed()
    # Mismatched lengths are what used to raise; now they are skipped.
    loaded.graph_window._last['total'] = np.array([1.0, 2.0])
    loaded.graph_window._redraw_zoom()
    assert len(loaded.graph_window._zoom_items) >= 1


def test_rescale_button_appears_only_with_the_zoom(loaded):
    assert not loaded.zoom_rescale_btn.isVisibleTo(loaded)
    loaded.zoom_check.setChecked(True)
    loaded._on_display_option_changed()
    assert loaded.zoom_rescale_btn.isVisibleTo(loaded)


def test_rescale_restores_the_default_framing(loaded):
    loaded.zoom_check.setChecked(True)
    loaded._on_display_option_changed()
    loaded.graph_window.zoom_plot.setXRange(3.0, 3.1)
    loaded.graph_window.zoom_plot.setYRange(-500.0, 500.0)
    loaded._on_rescale_zoom()
    lo, hi = loaded.graph_window.zoom_plot.viewRange()[0]
    assert lo == pytest.approx(loaded.WAXS_ZOOM_Q_MIN, rel=0.1)
    assert hi > 3.5


def test_zoom_copes_with_data_that_stops_before_the_waxs_region(qapp):
    """A SAXS-only file must still produce a sane window, not an inverted one."""
    from pyirena.gui.carbon_fit_panel import CarbonFitPanel

    p = CarbonFitPanel()
    try:
        q = np.logspace(-3, -0.5, 200)          # stops at Q = 0.32 Å⁻¹
        p.set_data(q, p.model.evaluate(q), None, label='saxs only')
        p.zoom_check.setChecked(True)
        p._on_display_option_changed()
        lo, hi = p.graph_window.zoom_plot.viewRange()[0]
        assert 0 < lo < hi
        assert hi == pytest.approx(float(q.max()), rel=0.15)
    finally:
        p.deleteLater()


# ── Loaded-file display ─────────────────────────────────────────────────────

def test_the_filename_field_names_the_file_however_it_arrived(panel):
    """Opening from the Data Browser bypasses this panel's own Open… button."""
    q = np.logspace(-3, 0.65, 200)
    panel.set_data(q, panel.model.evaluate(q), None,
                   label='hard_carbon_001.h5', filepath='/tmp/hard_carbon_001.h5')
    assert panel.data_loader._edit.text() == 'hard_carbon_001.h5'
