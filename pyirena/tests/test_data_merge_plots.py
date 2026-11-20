"""
Data Merge display modes and hand-tuning (GitHub issue #35).

Two things here are easy to get wrong and expensive when wrong:

* **A display mode must not reach the maths.** Porod is log(I·Q⁴) vs log(Q) on screen
  and nothing else — Optimize and merge still see the original Q and I. A
  transform that leaked into the merge would silently change every saved
  curve.
* **The cursors carry physical Q.** They are what defines the overlap region,
  so switching view must not move them, in either direction.

The hand-tuning path is the other half: with poor data Optimize does not
converge to anything useful, so the user nudges the three numbers by hand.
That path must produce exactly what Save then writes.
"""

from __future__ import annotations

import numpy as np
import pytest


def _qt_or_skip():
    try:
        from pyirena.gui._qt import QApplication
    except ImportError:
        pytest.skip("Qt (PySide6/PyQt6) not available")
    return QApplication.instance() or QApplication([])


@pytest.fixture
def panel():
    _qt_or_skip()
    from pyirena.gui.data_merge_panel import DataMergePanel

    p = DataMergePanel()
    q1 = np.logspace(-3, -1.3, 120)
    I1 = 1e-4 * q1 ** -4 + 0.5
    q2 = np.logspace(-1.6, -0.3, 90)
    I2 = (1e-4 * q2 ** -4) * 0.8
    p._data1 = {'Q': q1, 'Intensity': I1, 'Error': 0.03 * I1, 'slit_length': 0.0}
    p._data2 = {'Q': q2, 'Intensity': I2, 'Error': 0.03 * I2, 'slit_length': 0.0}
    p._graph.plot_ds1(q1, I1, 0.03 * I1)
    p._graph.plot_ds2(q2, I2, 0.03 * I2)
    p._graph.init_cursors(0.02, 0.05)
    return p


# ── the Porod view ──────────────────────────────────────────────────────────

class TestPorodMode:
    def test_all_three_modes_are_offered(self, panel):
        keys = [panel._mode_combo.itemData(i)
                for i in range(panel._mode_combo.count())]
        assert keys == ['saxs', 'waxs', 'porod']

    def test_transform_is_i_q4_versus_q(self, panel):
        """y is scaled by Q⁴; x stays plain Q."""
        from pyirena.gui.data_merge_panel import MODE_POROD

        panel._graph.set_mode(MODE_POROD)
        g = panel._graph
        q = np.array([0.01, 0.05, 0.2])
        I = np.array([100.0, 2.0, 0.1])
        np.testing.assert_allclose(g._tx(q), q)
        np.testing.assert_allclose(g._ty(q, I), I * q ** 4)

    @pytest.mark.parametrize("mode", ["saxs", "waxs", "porod"])
    def test_axis_position_round_trips(self, panel, mode):
        """Cursor placement must be exactly invertible, or dragging drifts."""
        panel._graph.set_mode(mode)
        g = panel._graph
        for q in (1e-4, 0.001, 0.02, 0.3, 1.0):
            assert g._q_from_axis_pos(g._axis_pos(q)) == pytest.approx(q, rel=1e-12)

    def test_switching_mode_does_not_move_the_cursors(self, panel):
        """The overlap region is the user's choice; a view change is not."""
        before = panel._graph.get_overlap_range()
        for mode in ("porod", "waxs", "saxs", "porod", "saxs"):
            panel._graph.set_mode(mode)
            assert panel._graph.get_overlap_range() == pytest.approx(before)

    def test_default_axis_scales_per_mode(self, panel):
        """Porod is log(I·Q⁴) vs log(Q): the Q⁴ factor stretches the y axis in
        the overlap region, which is what makes a mismatch visible."""
        g = panel._graph
        for mode, expect in (("saxs", True), ("porod", True), ("waxs", False)):
            g.set_mode(mode)
            assert g._x_log is expect and g._y_log is expect

    @pytest.mark.parametrize("mode", ["saxs", "waxs", "porod"])
    def test_menu_transform_keeps_cursors_on_the_same_q(self, panel, mode):
        """Right-click → Plot options → Transforms flips the axis scale under
        the cursors. They must stay on the same physical Q, and the overlap
        handed to the merge must not change (it used to read 10**0.05)."""
        g = panel._graph
        g.set_mode(mode)
        g.init_cursors(0.03, 0.08)
        for _ in range(2):  # lin→log→lin, or log→lin→log
            box = g._plot.ctrl.logXCheck
            box.setChecked(not box.isChecked())
            assert g._x_log is box.isChecked()
            assert g.get_overlap_range() == pytest.approx((0.03, 0.08))
        # Dragged cursor position survives a toggle too.
        g._cursor_a.setPos(g._axis_pos(0.02))
        box.setChecked(not box.isChecked())
        assert g.get_overlap_range() == pytest.approx((0.02, 0.08))

    def test_axis_labels_follow_the_mode(self, panel):
        panel._graph.set_mode("porod")
        x, y = panel._graph._axis_labels()
        assert "Q⁴" not in x and "I·Q⁴" in y
        panel._graph.set_mode("saxs")
        x, y = panel._graph._axis_labels()
        assert x.startswith("Q") and "Q⁴" not in x


class TestModeNeverReachesTheMaths:
    """The whole safety argument for the feature, in one class."""

    def test_merged_curve_is_identical_in_every_mode(self, panel):
        results = {}
        for mode in ("saxs", "waxs", "porod"):
            panel._graph.set_mode(mode)
            panel._scale_result.setText("0.8")
            panel._bg_result.setText("0.5")
            panel._qshift_result.setText("0.0")
            panel._live_pending = True
            panel._do_live_update()
            results[mode] = (panel._last_q_merged.copy(),
                             panel._last_I_merged.copy())

        q_ref, I_ref = results["saxs"]
        for mode in ("waxs", "porod"):
            q, I = results[mode]
            np.testing.assert_array_equal(q, q_ref)
            np.testing.assert_array_equal(I, I_ref)

    def test_config_carries_physical_q(self, panel):
        """Overlap bounds handed to the engine are Q, never Q⁴."""
        panel._graph.set_mode("porod")
        panel._live_pending = True
        panel._do_live_update()
        cfg = panel._last_config
        assert cfg.q_overlap_min == pytest.approx(0.02)
        assert cfg.q_overlap_max == pytest.approx(0.05)


# ── hand tuning ─────────────────────────────────────────────────────────────

class TestManualTuning:
    def test_fields_stay_editable_while_fit_is_checked(self, panel):
        """The whole point of #35 item 2: tweak what Optimize produced."""
        for chk, field in ((panel._fit_scale_chk, panel._scale_result),
                           (panel._fit_qshift_chk, panel._qshift_result),
                           (panel._fit_bg_chk, panel._bg_result)):
            chk.setChecked(True)
            assert not field.isReadOnly()

    def test_live_update_does_not_optimise(self, panel):
        """Every fit flag off — apply these numbers and merge, nothing else."""
        panel._scale_result.setText("0.77")
        panel._live_pending = True
        panel._do_live_update()
        cfg = panel._last_config
        assert (cfg.fit_scale, cfg.fit_qshift, cfg.fit_background) == (False, False, False)
        assert panel._last_result.scale == pytest.approx(0.77)

    def test_the_typed_value_is_the_one_used(self, panel):
        for scale in (0.5, 1.0, 2.5):
            panel._scale_result.setText(str(scale))
            panel._live_pending = True
            panel._do_live_update()
            assert panel._last_result.scale == pytest.approx(scale)

    def test_mid_typing_garbage_is_ignored(self, panel):
        """A half-typed number must not blank the preview."""
        panel._scale_result.setText("1.0")
        panel._live_pending = True
        panel._do_live_update()
        good = panel._last_q_merged.copy()

        panel._scale_result.setText("1.0e")      # not yet a number
        panel._live_pending = True
        panel._do_live_update()
        np.testing.assert_array_equal(panel._last_q_merged, good)

    def test_preview_is_what_save_would_write(self, panel):
        """_last_* is the single source for both the preview and the file."""
        panel._scale_result.setText("0.9")
        panel._live_pending = True
        panel._do_live_update()
        assert panel._last_q_merged is not None
        assert len(panel._last_q_merged) == len(panel._last_I_merged)
        assert np.all(np.diff(panel._last_q_merged) >= 0), "merged Q must be sorted"

    def test_split_at_left_cursor_reaches_the_live_path(self, panel):
        """Item 3: trimming still follows the checkbox when tuning by hand."""
        panel._split_chk.setChecked(False)
        panel._live_pending = True
        panel._do_live_update()
        n_all = len(panel._last_q_merged)

        panel._split_chk.setChecked(True)
        panel._live_pending = True
        panel._do_live_update()
        assert panel._last_config.split_at_left_cursor is True
        # A hard split drops the overlapping points, so it cannot be longer.
        assert len(panel._last_q_merged) <= n_all


# ── the nudge field ─────────────────────────────────────────────────────────

class TestNudgeField:
    def test_step_is_a_percentage_of_the_value(self):
        _qt_or_skip()
        from pyirena.gui.nudge_field import NudgeField

        f = NudgeField("2.0")
        f.nudge(+1)
        assert float(f.text()) == pytest.approx(2.02)
        f.setText("1e-5")
        f.nudge(+1)
        assert float(f.text()) == pytest.approx(1.01e-5)

    def test_modifiers_choose_coarse_and_fine(self):
        _qt_or_skip()
        from pyirena.gui.nudge_field import STEP_COARSE, STEP_FINE, NudgeField

        f = NudgeField("100.0")
        f.nudge(+1, STEP_COARSE)
        assert float(f.text()) == pytest.approx(110.0)
        f.setText("100.0")
        f.nudge(+1, STEP_FINE)
        assert float(f.text()) == pytest.approx(100.1)

    def test_zero_uses_the_absolute_fallback(self):
        """A background starts at 0, and 1% of 0 is 0 — the field must still move."""
        _qt_or_skip()
        from pyirena.gui.nudge_field import NudgeField

        f = NudgeField("0.0")
        f.set_fallback_step(0.002)
        f.nudge(+1)
        assert float(f.text()) == pytest.approx(0.002)
        f.setText("0.0")
        f.nudge(-1)
        assert float(f.text()) == pytest.approx(-0.002)

    def test_non_numeric_text_is_left_alone(self):
        _qt_or_skip()
        from pyirena.gui.nudge_field import NudgeField

        f = NudgeField("not a number")
        f.nudge(+1)
        assert f.text() == "not a number"

    def test_nudging_emits_exactly_once(self):
        _qt_or_skip()
        from pyirena.gui.nudge_field import NudgeField

        f = NudgeField("1.0")
        seen = []
        f.nudged.connect(lambda: seen.append(f.text()))
        f.nudge(+1)
        f.nudge(-1)
        assert len(seen) == 2

    def test_buttons_do_nothing_on_a_read_only_field(self):
        _qt_or_skip()
        from pyirena.gui._qt import QPushButton
        from pyirena.gui.nudge_field import NudgeField, make_nudge_buttons

        f = NudgeField("1.0")
        f.setReadOnly(True)
        holder = make_nudge_buttons(f)
        for btn in holder.findChildren(QPushButton):
            btn.click()
        assert f.text() == "1.0"


# ── the right-hand axis ─────────────────────────────────────────────────────

class TestRightAxis:
    def test_mirror_is_on_by_default(self, panel):
        assert panel._graph._plot.getAxis('right').isVisible()
        assert panel._graph._ds2_vb is None

    def test_second_axis_can_be_turned_on_and_off(self, panel):
        panel._ds2_right_chk.setChecked(True)
        assert panel._graph._ds2_vb is not None
        assert "DS2" in panel._graph._plot.getAxis('right').labelText

        panel._ds2_right_chk.setChecked(False)
        assert panel._graph._ds2_vb is None

    def test_it_survives_a_mode_change(self, panel):
        """set_mode() rebuilds the PlotItem, which drops the second ViewBox."""
        panel._ds2_right_chk.setChecked(True)
        panel._on_mode_changed(0)
        panel._mode_combo.setCurrentIndex(2)       # Porod
        assert panel._graph._ds2_vb is not None

    def test_second_axis_does_not_change_the_merge(self, panel):
        panel._scale_result.setText("0.8")
        panel._live_pending = True
        panel._do_live_update()
        before = panel._last_I_merged.copy()

        panel._ds2_right_chk.setChecked(True)
        panel._live_pending = True
        panel._do_live_update()
        np.testing.assert_array_equal(panel._last_I_merged, before)


# ── saved state ─────────────────────────────────────────────────────────────

class TestStateRoundTrip:
    def test_porod_and_right_axis_survive_save_load(self, panel):
        panel._mode_combo.setCurrentIndex(2)        # Porod
        panel._ds2_right_chk.setChecked(True)
        panel.save_state()

        panel._mode_combo.setCurrentIndex(0)
        panel._ds2_right_chk.setChecked(False)
        panel.load_state()
        assert panel._mode_combo.currentData() == 'porod'
        assert panel._ds2_right_chk.isChecked()

    def test_a_pre_issue35_state_still_loads(self, panel):
        """Old states hold only 'saxs'/'waxs' and no ds2_right_axis key."""
        # Replace the section rather than update() it: update() merges,
        # which would leave a key an old state file could not have had.
        panel._sm.state['data_merge'] = {'plot_mode': 'waxs'}
        panel.load_state()
        assert panel._mode_combo.currentData() == 'waxs'
        assert not panel._ds2_right_chk.isChecked()

    def test_an_unknown_mode_falls_back_to_saxs(self, panel):
        panel._sm.state['data_merge'] = {'plot_mode': 'something_newer'}
        panel.load_state()
        assert panel._mode_combo.currentData() == 'saxs'


# ── recovering from unusable cursor positions ───────────────────────────────
#
# Before the menu-transform fix, switching an axis from the right-click menu
# made the panel read log10(Q) as Q and save it (q_overlap_min = -1.39). On the
# next start log10 of that gave NaN: invisible cursors, "nan" in both boxes,
# and no control that moved them.

def _damage(g, *which):
    """Set cursor positions to NaN without letting the panel react yet."""
    for c in which:
        c.blockSignals(True)
        c.setPos(float('nan'))
        c.blockSignals(False)


class TestCursorRecovery:
    CORRUPT = (-1.388269756027181, -0.8961919290668727)   # real values from a state file

    def test_graph_refuses_unusable_positions(self, panel):
        g = panel._graph
        before = g.get_overlap_range()
        for bad in (self.CORRUPT, (float('nan'), 0.05), (0.02, float('inf')), (0.0, 0.05)):
            assert g.init_cursors(*bad) is False
            assert g.get_overlap_range() == pytest.approx(before)

    def test_nan_cursor_is_reported_as_none_not_nan(self, panel):
        """A damaged cursor must never reach the merge or the state file."""
        g = panel._graph
        _damage(g, g._cursor_a)
        assert g.get_overlap_range() == (None, None)

    @pytest.mark.parametrize("saved", [CORRUPT, (float('nan'), float('nan')),
                                       (5.0, 9.0), (0.05, 0.02)])
    def test_bad_saved_range_falls_back_to_default(self, panel, saved):
        """Fresh start with damaged state: cursors appear at the default."""
        panel._sm.update('data_merge', {'q_overlap_min': saved[0],
                                        'q_overlap_max': saved[1]})
        panel._graph._cursor_a = panel._graph._cursor_b = None
        panel._init_cursors_from_data()
        default = panel._default_overlap_range()
        assert panel._graph.get_overlap_range() == pytest.approx(default)
        assert panel._q_min_edit.text() != "" and "nan" not in panel._q_min_edit.text()

    def test_good_saved_range_is_still_honoured(self, panel):
        panel._sm.update('data_merge', {'q_overlap_min': 0.03, 'q_overlap_max': 0.06})
        panel._graph._cursor_a = panel._graph._cursor_b = None
        panel._init_cursors_from_data()
        assert panel._graph.get_overlap_range() == pytest.approx((0.03, 0.06))

    def test_reset_button_recovers_damaged_cursors(self, panel):
        g = panel._graph
        _damage(g, g._cursor_a, g._cursor_b)
        assert g.get_overlap_range() == (None, None)
        panel._reset_cursors_btn.click()
        assert panel._graph.get_overlap_range() == pytest.approx(
            panel._default_overlap_range())

    def test_display_heals_itself(self, panel):
        """No button press needed: nan cursors are replaced when noticed."""
        g = panel._graph
        _damage(g, g._cursor_a)
        panel._update_cursor_display()
        assert panel._graph.get_overlap_range()[0] is not None
        assert "nan" not in panel._q_min_edit.text()

    def test_reset_without_data_is_harmless(self, panel):
        panel._data1 = None
        panel._reset_cursors()           # must not raise
        assert "Load DS1" in panel._status.text()

    @pytest.mark.parametrize("mode", ["porod", "waxs", "saxs"])
    def test_q_boxes_follow_a_dragged_cursor_after_a_mode_switch(self, panel, mode):
        """The cursors are rebuilt on a mode switch; the display must stay wired."""
        panel._graph.set_mode(mode)
        panel._graph._cursor_a.setPos(panel._graph._axis_pos(0.025))
        assert panel._q_min_edit.text() == "0.025"
