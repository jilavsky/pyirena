"""
Tests for :mod:`pyirena.gui.plot_export` — the shared plot export used by every
pyqtgraph plot in the GUI.

Skipped when no Qt binding is installed; the `test-gui` CI job runs them with
``QT_QPA_PLATFORM=offscreen``.

The behaviours pinned here are the ones that break silently:
  * curves come back in **linear** units even when the plot is in log mode
    (``getData()`` returns log10 values and would double-log in Igor);
  * error-bar segments are not mistaken for data curves;
  * CSV pads curves of unequal length instead of truncating;
  * a save dialog with no typed extension still writes the chosen format;
  * the ITX file keeps its IGOR header, log-axis and data-folder commands.
"""

import numpy as np
import pytest


def _require_qt() -> None:
    """Skip this module when no Qt binding is installed.

    ``pytest.importorskip`` is not enough: ``pyirena.gui._qt`` raises a plain
    ImportError carrying an installation hint, and pytest >= 8.2 treats that as
    a *broken* module rather than a missing one — so the plain (no-GUI) CI job
    would report a collection error instead of a skip.
    """
    try:
        import pyirena.gui._qt  # noqa: F401
    except ImportError:
        pytest.skip("Qt (PySide6/PyQt6) not available", allow_module_level=True)


_require_qt()

pytest.importorskip("pyqtgraph")

import pyqtgraph as pg  # noqa: E402

from pyirena.gui import plot_export as pe  # noqa: E402
from pyirena.gui._qt import QApplication  # noqa: E402
from pyirena.gui.sas_plot import make_sas_plot, plot_iq_data, plot_iq_model  # noqa: E402

# Module-level reference: a QApplication collected while pyqtgraph scenes are
# still alive segfaults the interpreter on exit, which would fail CI even
# though every test passed.
_APP = None


@pytest.fixture(scope="session")
def qapp():
    global _APP
    _APP = QApplication.instance() or QApplication([])
    yield _APP


@pytest.fixture
def layout(qapp):
    """A GraphicsLayoutWidget torn down before the test ends.

    pyqtgraph scenes must be closed and deleted explicitly; leaving them to the
    cyclic garbage collector at interpreter shutdown crashes Qt.
    """
    widget = pg.GraphicsLayoutWidget()
    yield widget
    widget.clear()
    widget.close()
    widget.deleteLater()
    qapp.processEvents()


@pytest.fixture
def iq_plot(layout):
    """A log-log I(Q) plot with a data scatter (+dI) and a model curve."""
    plot = make_sas_plot(layout, 0, 0, parent_widget=layout,
                         jpeg_default_name="test_graph")
    q = np.logspace(-3, -1, 40)
    intensity = 1e3 * q ** -3
    plot_iq_data(plot, q, intensity, 0.05 * intensity, label="Data")
    plot_iq_model(plot, q, 1.05 * intensity, label="Model")
    return plot


@pytest.fixture(autouse=True)
def _no_folder_persistence(monkeypatch):
    """Never read or write the user's real state file from a test.

    Without stubbing the *read* as well, folder-resolution tests pass or fail
    depending on where the developer last saved a graph on that machine.
    """
    monkeypatch.setattr(pe, "_LAST_EXPORT_FOLDER", None, raising=False)
    monkeypatch.setattr(pe, "remember_export_folder", lambda path: None)
    monkeypatch.setattr(pe, "_saved_folders", lambda: (None, None))


def _accept_dialog(monkeypatch, path, selected_filter=""):
    """Make the next QFileDialog.getSaveFileName return *path*."""
    monkeypatch.setattr(
        pe.QFileDialog, "getSaveFileName",
        staticmethod(lambda *a, **k: (str(path), selected_filter)),
    )


# ── Curve collection ─────────────────────────────────────────────────────

def test_collect_returns_named_curves_only(iq_plot):
    curves = pe.collect_plot_curves(iq_plot)
    assert [c.name for c in curves] == ["Data", "Model"]


def test_collect_returns_linear_values_from_a_log_plot(iq_plot):
    """The plot is in log mode; exported values must still be linear."""
    data = pe.collect_plot_curves(iq_plot)[0]
    assert data.x.min() == pytest.approx(1e-3)
    assert data.x.max() == pytest.approx(1e-1)
    # I = 1e3 q^-3 → at q = 1e-3 that is 1e12, not log10 of it.
    assert data.y.max() == pytest.approx(1e12, rel=1e-6)


def test_collect_keeps_uncertainties_for_data_and_none_for_model(iq_plot):
    data, model = pe.collect_plot_curves(iq_plot)
    assert data.dy is not None and len(data.dy) == len(data.x)
    assert model.dy is None


def test_collect_skips_error_bar_segments(iq_plot):
    """Error bars are unnamed NaN-separated lines — never a curve."""
    n_items = len(iq_plot.listDataItems())
    assert n_items > len(pe.collect_plot_curves(iq_plot))


def test_unnamed_curves_are_exported_with_an_axis_derived_name(layout):
    """A residual panel labels nothing — exporting must still work.

    Regression: Unified Fit, Modeling and Simple Fits residuals, the Contrast
    (Δρ)² plot and the Simple Fits linearization all plot without ``name=``,
    and every export reported "No named data curves found to export" while the
    curve was plainly visible on screen.
    """
    plot = layout.addPlot()
    plot.setLabel("left", "Residuals r' (rescaled)")
    plot.plot([1.0, 2.0, 3.0], [0.1, -0.2, 0.3], pen=None, symbol="o")

    curves = pe.collect_plot_curves(plot)
    assert [c.name for c in curves] == ["Residuals r' (rescaled)"]
    assert list(curves[0].y) == [0.1, -0.2, 0.3]


def test_several_unnamed_curves_get_distinct_names(layout):
    plot = layout.addPlot()
    plot.setLabel("left", "Y")
    plot.plot([1.0, 2.0], [1.0, 2.0])
    plot.plot([1.0, 2.0], [3.0, 4.0])
    assert [c.name for c in pe.collect_plot_curves(plot)] == ["Y 1", "Y 2"]


def test_named_curves_win_over_unnamed_ones(layout):
    """The fallback must not smuggle decoration into a normal export."""
    plot = layout.addPlot()
    plot.plot([1.0, 2.0], [1.0, 2.0], name="Data")
    plot.plot([1.0, 2.0], [9.0, 9.0])          # unnamed guide line
    assert [c.name for c in pe.collect_plot_curves(plot)] == ["Data"]


def test_scatter_plot_items_are_collected(layout):
    """The Simple Fits linearization and collect windows use ScatterPlotItem.

    These are not ``PlotDataItem`` and were skipped entirely, so those plots
    could not be exported at all.
    """
    plot = layout.addPlot()
    plot.addItem(pg.ScatterPlotItem(x=[1.0, 2.0, 3.0], y=[4.0, 5.0, 6.0],
                                    name="Linearized data"))
    curves = pe.collect_plot_curves(plot)
    assert [c.name for c in curves] == ["Linearized data"]
    assert list(curves[0].x) == [1.0, 2.0, 3.0]


def test_error_bar_segments_are_never_exported_as_a_curve(layout):
    """Even in the unnamed fallback, NaN-separated bars must not become data."""
    plot = layout.addPlot()
    plot.setLabel("left", "I")
    x = [1.0, 1.0, np.nan, 2.0, 2.0, np.nan, 3.0, 3.0, np.nan]
    y = [0.9, 1.1, np.nan, 1.9, 2.1, np.nan, 2.9, 3.1, np.nan]
    plot.plot(x, y, connect="finite")
    assert pe.collect_plot_curves(plot) == []


def test_tag_curve_uncertainty_puts_dy_into_the_export(layout):
    """Panels that draw their own error bars must still export a dY column.

    Regression: only ``plot_iq_data`` recorded uncertainties, so Unified Fit,
    Sizes and WAXS — which draw their own scatter — exported no dY even though
    error bars were on screen.
    """
    plot = layout.addPlot()
    item = plot.plot([1.0, 2.0, 3.0], [10.0, 20.0, 30.0], name="Data")
    pe.tag_curve_uncertainty(item, [1.0, 2.0, 3.0])

    curve = pe.collect_plot_curves(plot)[0]
    assert curve.dy is not None
    assert list(curve.dy) == [1.0, 2.0, 3.0]
    assert "dY (Data)" in pe.curves_to_csv_text([curve]).splitlines()[0]


def test_collect_picks_up_explicit_itx_export_items(layout):
    """Step-mode bar charts stash their data as ``_itx_export``."""
    plot = layout.addPlot()
    item = pg.PlotCurveItem([0, 1, 2, 3], [1.0, 2.0, 3.0], stepMode=True)
    item._itx_export = {"name": "P(r)", "x": [0.5, 1.5, 2.5], "y": [1.0, 2.0, 3.0]}
    plot.addItem(item)
    curves = pe.collect_plot_curves(plot)
    assert [c.name for c in curves] == ["P(r)"]
    assert list(curves[0].y) == [1.0, 2.0, 3.0]


# ── CSV of curve data (A4) ───────────────────────────────────────────────

def test_curves_to_csv_has_one_column_group_per_curve(iq_plot):
    text = pe.curves_to_csv_text(pe.collect_plot_curves(iq_plot))
    assert text.splitlines()[0] == "X (Data),Y (Data),dY (Data),X (Model),Y (Model)"


def test_curves_to_csv_pads_unequal_curves(qapp):
    curves = [
        pe.PlotCurve("short", np.array([1.0, 2.0]), np.array([3.0, 4.0]), "#000", None),
        pe.PlotCurve("long", np.arange(3.0), np.arange(3.0), "#000", None),
    ]
    lines = pe.curves_to_csv_text(curves).splitlines()
    assert lines[0] == "X (short),Y (short),X (long),Y (long)"
    assert lines[3] == ",,2,2"          # short curve padded, long one continues


def test_save_curves_as_csv_writes_the_file(iq_plot, tmp_path, monkeypatch):
    target = tmp_path / "curves.csv"
    _accept_dialog(monkeypatch, target)
    written = pe.save_curves_as_csv(iq_plot, None, "test_graph")
    assert written == str(target)
    text = target.read_text(encoding="utf-8")
    assert text.startswith("X (Data),Y (Data),dY (Data)")
    assert len(text.splitlines()) == 41       # header + 40 points


def test_save_curves_as_csv_appends_the_extension(iq_plot, tmp_path, monkeypatch):
    _accept_dialog(monkeypatch, tmp_path / "curves")
    written = pe.save_curves_as_csv(iq_plot, None, "test_graph")
    assert written.endswith(".csv")


def test_save_curves_as_csv_warns_when_nothing_is_plotted(layout, monkeypatch):
    empty = layout.addPlot()
    warned = []
    monkeypatch.setattr(pe.QMessageBox, "warning",
                        staticmethod(lambda *a, **k: warned.append(a)))
    assert pe.save_curves_as_csv(empty, None, "empty") is None
    assert warned


# ── Images and clipboard (A3) ────────────────────────────────────────────

def test_copy_plot_to_clipboard_puts_an_image_on_the_clipboard(iq_plot):
    assert pe.copy_plot_to_clipboard(iq_plot) is True
    image = QApplication.clipboard().image()
    assert not image.isNull()
    assert image.width() == pe.EXPORT_WIDTH


def test_save_plot_image_writes_png(iq_plot, tmp_path, monkeypatch):
    target = tmp_path / "graph.png"
    _accept_dialog(monkeypatch, target)
    written = pe.save_plot_image(iq_plot, None, "test_graph")
    assert written == str(target)
    assert target.stat().st_size > 0
    assert target.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_save_plot_image_defaults_to_png_without_an_extension(
        iq_plot, tmp_path, monkeypatch):
    _accept_dialog(monkeypatch, tmp_path / "graph")
    written = pe.save_plot_image(iq_plot, None, "test_graph")
    assert written.endswith(".png")
    assert (tmp_path / "graph.png").exists()


def test_save_plot_image_honours_the_selected_jpeg_filter(
        iq_plot, tmp_path, monkeypatch):
    _accept_dialog(monkeypatch, tmp_path / "graph", "JPEG image (*.jpg *.jpeg)")
    written = pe.save_plot_image(iq_plot, None, "test_graph")
    assert written.endswith(".jpg")
    assert (tmp_path / "graph.jpg").exists()


def test_save_widget_image_captures_the_whole_window(
        iq_plot, layout, tmp_path, monkeypatch):
    layout.resize(400, 300)
    _accept_dialog(monkeypatch, tmp_path / "window.png")
    written = pe.save_widget_image(layout, None, "test_window")
    assert written and (tmp_path / "window.png").exists()


# ── Igor ITX ─────────────────────────────────────────────────────────────

def test_itx_export_writes_waves_and_log_commands(iq_plot, tmp_path, monkeypatch):
    target = tmp_path / "graph.itx"
    _accept_dialog(monkeypatch, target)
    pe.save_itx_from_plot(iq_plot, None, "test_graph")

    text = target.read_text(encoding="utf-8")
    assert text.startswith("IGOR\n")
    assert "WAVES/D  X_Data_01" in text
    assert "WAVES/D  Y_Data_01" in text
    assert "WAVES/D  Yerr_Data_01" in text     # dI came through
    assert "WAVES/D  Y_Model_02" in text
    assert "X ModifyGraph log(bottom)=1" in text
    assert "X ModifyGraph log(left)=1" in text
    assert "X ErrorBars Y_Data_01" in text
    # Linear values, not log10 — 1e12 would be 12 if double-logged.
    assert "1e+12" in text


def test_itx_marks_scatter_curves_as_markers_and_lines_as_lines(
        layout, tmp_path, monkeypatch):
    """Igor draws every imported wave as a line unless told otherwise.

    Regression: data points arrived in Igor as a zig-zag line, indistinguishable
    from a model curve — and an exported grey "outside fit range" scatter looked
    like a mystery trace across the graph.
    """
    plot = layout.addPlot()
    plot.plot([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], pen=None, symbol="o", name="Data")
    plot.plot([1.0, 2.0, 3.0], [1.1, 2.1, 3.1], name="Fit")

    _accept_dialog(monkeypatch, tmp_path / "styled.itx")
    pe.save_itx_from_plot(plot, None, "styled")
    text = (tmp_path / "styled.itx").read_text(encoding="utf-8")

    assert "X ModifyGraph mode(Y_Data_01)=3,marker(Y_Data_01)=19" in text
    assert "X ModifyGraph mode(Y_Fit_02)=0,lsize(Y_Fit_02)=1.5" in text


def test_itx_wave_names_stay_unique_when_labels_are_long(
        layout, tmp_path, monkeypatch):
    """Igor wave names cap at 31 characters — truncate the label, not the index.

    Two curves whose labels share a long prefix would otherwise get the same
    wave name and Igor would silently overwrite the first one.
    """
    plot = layout.addPlot()
    plot.plot([1.0, 2.0], [1.0, 2.0],
              name="Volume distribution population one")
    plot.plot([1.0, 2.0], [3.0, 4.0],
              name="Volume distribution population two")

    _accept_dialog(monkeypatch, tmp_path / "long.itx")
    pe.save_itx_from_plot(plot, None, "long")
    text = (tmp_path / "long.itx").read_text(encoding="utf-8")

    waves = [ln.split()[-1] for ln in text.splitlines() if ln.startswith("WAVES")]
    assert len(waves) == len(set(waves)), f"duplicate wave names: {waves}"
    assert all(len(w) <= 31 for w in waves), waves
    assert any(w.endswith("_01") for w in waves)
    assert any(w.endswith("_02") for w in waves)


def test_curve_style_is_detected_from_the_item(layout):
    plot = layout.addPlot()
    plot.plot([1.0, 2.0], [1.0, 2.0], pen=None, symbol="o", name="scatter")
    plot.plot([1.0, 2.0], [1.0, 2.0], name="line")
    plot.plot([1.0, 2.0], [1.0, 2.0], symbol="o", name="both")
    plot.addItem(pg.ScatterPlotItem(x=[1.0, 2.0], y=[1.0, 2.0], name="item"))
    styles = {c.name: c.style for c in pe.collect_plot_curves(plot)}
    assert styles == {
        "scatter": "markers", "line": "line", "both": "both", "item": "markers",
    }


def test_itx_export_routes_into_an_igor_data_folder(iq_plot, tmp_path, monkeypatch):
    _accept_dialog(monkeypatch, tmp_path / "graph.itx")
    pe.save_itx_from_plot(iq_plot, None, "test_graph",
                          technique="Sizes", sample="AeroGel_500C.h5")
    text = (tmp_path / "graph.itx").read_text(encoding="utf-8")
    assert "X NewDataFolder/O/S root:Sizes" in text
    assert "X NewDataFolder/O/S root:Sizes:'AeroGel_500C'" in text   # extension dropped
    assert "X SetDataFolder root:" in text


def test_itx_folder_cmds_sanitises_names():
    open_lines, close_lines = pe._itx_folder_cmds("3D saxsMorph", "a'b.h5")
    assert open_lines[0] == "X NewDataFolder/O/S root:f_3D_saxsMorph"
    assert open_lines[1] == "X NewDataFolder/O/S root:f_3D_saxsMorph:'a_b'"
    assert close_lines == ["X SetDataFolder root:"]
    assert pe._itx_folder_cmds(None, "sample") == ([], [])


# ── Remembered folder (A8) ───────────────────────────────────────────────

def test_export_folder_uses_the_hint_before_anything_was_exported(tmp_path):
    """First export of the session lands next to the data."""
    assert pe.export_folder(tmp_path) == str(tmp_path)


def test_last_export_folder_beats_the_panel_data_folder(tmp_path, monkeypatch):
    """Where the user last saved wins over the panel's data folder.

    Regression: Unified Fit and Modeling passed their data folder as a hint,
    which took priority, so those two panels always reopened in the data
    folder and looked like the folder memory was broken — while panels that
    passed no hint remembered correctly.
    """
    desktop = tmp_path / "Desktop"
    data = tmp_path / "data"
    desktop.mkdir()
    data.mkdir()
    monkeypatch.setattr(pe, "_LAST_EXPORT_FOLDER", str(desktop), raising=False)
    assert pe.export_folder(data) == str(desktop)


def test_export_folder_accepts_a_file_and_returns_its_directory(tmp_path):
    a_file = tmp_path / "data.h5"
    a_file.write_text("")
    assert pe.export_folder(a_file) == str(tmp_path)


def test_export_folder_falls_back_to_home_for_a_missing_folder():
    """A remembered folder that has since been deleted must not be offered."""
    result = pe.export_folder("/no/such/folder/anywhere")
    assert result and "/no/such/folder" not in result


def test_saved_export_folder_beats_the_hint_after_a_restart(tmp_path, monkeypatch):
    saved = tmp_path / "saved"
    data = tmp_path / "data"
    saved.mkdir()
    data.mkdir()
    monkeypatch.setattr(pe, "_saved_folders", lambda: (str(saved), None))
    assert pe.export_folder(data) == str(saved)


# ── The menu itself ──────────────────────────────────────────────────────

def _export_entries(plot):
    return [
        a.text() for a in plot.getViewBox().menu.actions()
        if a.text() and ("Save" in a.text() or "Copy" in a.text())
    ]


def test_attach_plot_export_adds_the_standard_block(layout):
    plot = layout.addPlot()
    pe.attach_plot_export(plot, None, "graph", window=layout)
    assert _export_entries(plot) == [
        "Copy graph to clipboard",
        "Save graph as image…",
        "Save whole window as image…",
        "Save curve data as CSV…",
        "Save as Igor Pro ITX…",
    ]


def test_whole_window_entry_is_omitted_without_a_window(layout):
    plot = layout.addPlot()
    pe.attach_plot_export(plot, None, "graph")
    assert "Save whole window as image…" not in _export_entries(plot)


def test_itx_entry_can_be_suppressed(layout):
    """Windows with a purpose-built ITX writer must not show two ITX entries."""
    plot = layout.addPlot()
    pe.attach_plot_export(plot, None, "graph", itx=False)
    assert "Save as Igor Pro ITX…" not in _export_entries(plot)
    assert "Save curve data as CSV…" in _export_entries(plot)


def test_make_sas_plot_attaches_the_shared_menu(iq_plot):
    assert len(_export_entries(iq_plot)) == 5


def test_make_sas_plot_without_parent_has_no_export_menu(layout):
    """parent_widget=None suppresses the menu — a documented opt-out."""
    plot = make_sas_plot(layout, 0, 0, parent_widget=None)
    assert _export_entries(plot) == []


def test_folder_callable_is_evaluated_at_click_time(layout, tmp_path, monkeypatch):
    """Panels pass a callable because the data folder changes as files load."""
    plot = layout.addPlot()
    plot.plot([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], name="curve")

    folder = {"value": None}
    pe.attach_plot_export(plot, None, "graph", folder=lambda: folder["value"])
    folder["value"] = str(tmp_path)

    seen = {}

    def _fake_dialog(parent, caption, directory, filt):
        seen["directory"] = directory
        return "", ""

    monkeypatch.setattr(pe.QFileDialog, "getSaveFileName", staticmethod(_fake_dialog))
    csv_action = next(a for a in plot.getViewBox().menu.actions()
                      if a.text() == "Save curve data as CSV…")
    csv_action.trigger()
    assert seen["directory"].startswith(str(tmp_path))
