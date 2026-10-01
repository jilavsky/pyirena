"""
Data Merge display modes and hand-tuning (GitHub issue #35).

Two things here are easy to get wrong and expensive when wrong:

* **A display mode must not reach the maths.** Porod is I·Q⁴ vs Q⁴ on screen
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

    def test_porod_is_linear_not_log(self, panel):
        """A Porod plateau is only flat on linear axes — that is the point."""
        panel._graph.set_mode("porod")
        assert panel._graph._log_mode is False
        panel._graph.set_mode("saxs")
        assert panel._graph._log_mode is True

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
