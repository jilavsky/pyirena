"""
Source-level guard: every window a user can see answers ⌘W.

pyIrena's windows are independent top-level widgets with no common base class,
so a shortcut cannot be inherited — each window has to ask for one.  A window
that forgets is indistinguishable from a window that opted out, and the user
just finds that ⌘W does nothing *on that one window* (GitHub issue #16).

So the table below is the canonical list of **windows** versus **panes**, and
every Qt widget class in ``pyirena/gui`` must appear in exactly one of:

``WINDOWS``
    Shown as its own top-level window.  Must call
    ``install_standard_shortcuts`` (see ``gui/shortcuts.py``).
``EMBEDDED``
    Lives inside another widget's layout — a splitter pane, a tab, a row.  Its
    shortcuts come from the window it sits in.  Needs a written reason, because
    the ``*GraphWindow`` names make several of these look like windows: all of
    them except Contrast's are ``splitter.addWidget()`` panes.

A ``QDialog`` is exempt from both: Qt already closes one on Esc, and a dialog
that is up is the thing the user is answering.

Source text and import-only: no Qt, no display, so this runs in every CI job.
"""

from __future__ import annotations

import re
from pathlib import Path

GUI = Path(__file__).resolve().parents[1] / "gui"

#: Classes shown as their own top-level window → must install shortcuts.
WINDOWS: set[str] = {
    "carbon_fit_panel.py::CarbonFitPanel",
    "contrast_panel.py::ContrastGraphWindow",
    "contrast_panel.py::ContrastPanel",
    "data_manipulation_panel.py::DataManipulationPanel",
    "data_merge_panel.py::DataMergePanel",
    "data_selector/panel.py::DataSelectorPanel",
    "data_selector/results_windows.py::CarbonFitResultsWindow",
    "data_selector/results_windows.py::GraphWindow",
    "data_selector/results_windows.py::SimpleFitResultsWindow",
    "data_selector/results_windows.py::SizeDistResultsWindow",
    "data_selector/results_windows.py::TabulateResultsWindow",
    "data_selector/results_windows.py::UnifiedFitResultsWindow",
    "data_selector/results_windows.py::WAXSPeakFitResultsWindow",
    "feature_identifier.py::FeatureIdentifierDialog",
    "fractals_panel.py::FractalsGraphWindow",
    "hdf5viewer/collect_window.py::CollectWindow",
    "hdf5viewer/graph_window.py::GraphWindow",
    "hdf5viewer/main_window.py::HDF5ViewerWindow",
    "hdf5viewer/multi_collect_window.py::MultiCollectWindow",
    "modeling_panel.py::ModelingPanel",
    "saxs_morph_3d.py::VoxelViewerWindow",
    "saxs_morph_panel.py::SaxsMorphPanel",
    "simple_fits_panel.py::SimpleFitsPanel",
    "sizes_panel.py::SizesFitPanel",
    "unified_fit.py::UnifiedFitPanel",
    "waxs_peakfit_panel.py::WAXSPeakFitPanel",
}

#: Classes that live inside another widget → why they are not windows.
EMBEDDED: dict[str, str] = {
    "carbon_fit_panel.py::CarbonFitGraphWindow":
        "splitter pane of CarbonFitPanel despite the name",
    "data_loading.py::DataFileLoaderRow": "a row of controls",
    "data_manipulation_panel.py::DataManipulationGraphWindow":
        "splitter pane of DataManipulationPanel",
    "data_manipulation_panel.py::_ManipFileBrowser": "file list inside the panel",
    "data_merge_panel.py::DataMergeGraphWindow": "splitter pane of DataMergePanel",
    "data_merge_panel.py::_DatasetSelectorWidget": "DS1/DS2 picker inside the panel",
    "diffraction_lines_panel.py::DiffractionLinesPanel":
        "embedded in the WAXS Peak Fit panel",
    "fractals_panel.py::FractalsPanel":
        "central widget of FractalsGraphWindow, which is the window",
    "hdf5viewer/export_to_igor_tab.py::ExportToIgorTab": "a tab of the Data Explorer",
    "hdf5viewer/file_tree.py::FileTreeWidget": "tree pane of the Data Explorer",
    "hdf5viewer/hdf5_browser.py::HDF5BrowserWidget": "browser pane of the Data Explorer",
    "hdf5viewer/plot_controls.py::PlotControlsPanel": "controls pane of the graph window",
    "modeling_panel.py::ModelingGraphWindow": "splitter pane of ModelingPanel",
    "modeling_panel.py::PopulationTab": "one tab per population inside the panel",
    "q_range_ui.py::QRangeFields": "a pair of Q entry fields",
    "saxs_morph_3d.py::Slice2DViewer": "2-D slice pane; pops out via _PopoutDialog",
    "saxs_morph_3d.py::Voxel3DViewer": "3-D view pane; pops out via _PopoutDialog",
    "saxs_morph_panel.py::SaxsMorphGraphWindow": "splitter pane of SaxsMorphPanel",
    "simple_fits_panel.py::SimpleFitsGraphWindow": "splitter pane of SimpleFitsPanel",
    "sizes_panel.py::SizesFitGraphWindow": "splitter pane of SizesFitPanel",
    "unified_fit.py::LevelParametersWidget": "one level's parameter rows",
    "unified_fit.py::UnifiedFitGraphWindow": "splitter pane of UnifiedFitPanel",
    "waxs_peakfit_panel.py::PeakRowWidget": "one peak's row of controls",
    "waxs_peakfit_panel.py::WAXSPeakFitGraphWindow": "splitter pane of WAXSPeakFitPanel",
}

_CLASS = re.compile(r"^class\s+(\w+)\s*\(([^)]*)\)\s*:", re.M)
_QT_BASE = re.compile(r"\bQWidget\b|\bQMainWindow\b|\bQDialog\b")


def _qt_classes() -> dict[str, tuple[str, str]]:
    """``key -> (bases, class body)`` for every Qt widget class under gui/."""
    found: dict[str, tuple[str, str]] = {}
    for path in sorted(GUI.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        src = path.read_text(encoding="utf-8")
        matches = list(_CLASS.finditer(src))
        for i, m in enumerate(matches):
            bases = m.group(2)
            if not _QT_BASE.search(bases):
                continue
            end = matches[i + 1].start() if i + 1 < len(matches) else len(src)
            key = f"{path.relative_to(GUI).as_posix()}::{m.group(1)}"
            found[key] = (bases, src[m.start():end])
    return found


def test_every_qt_class_is_classified():
    """A new window must be declared a window or a pane — silence is not an option."""
    unclassified = sorted(
        key for key, (bases, _) in _qt_classes().items()
        if key not in WINDOWS and key not in EMBEDDED and "QDialog" not in bases
    )
    assert not unclassified, (
        "New Qt widget class(es) not in WINDOWS or EMBEDDED in this file:\n  "
        + "\n  ".join(unclassified)
        + "\n\nIf it is shown as its own window, add it to WINDOWS and call "
          "install_standard_shortcuts(self) in it so ⌘W works.\n"
          "If it lives inside another widget, add it to EMBEDDED with the reason."
    )


def test_declared_windows_install_shortcuts():
    """Every window in the table actually makes the call."""
    classes = _qt_classes()
    missing = sorted(
        key for key in WINDOWS
        if "install_standard_shortcuts(" not in classes.get(key, ("", ""))[1]
    )
    assert not missing, (
        "Window class(es) that never call install_standard_shortcuts, so ⌘W "
        "does nothing on them:\n  " + "\n  ".join(missing)
    )


def test_table_has_no_stale_entries():
    """A renamed or deleted class must not leave a row behind."""
    known = set(_qt_classes())
    stale = sorted((WINDOWS | set(EMBEDDED)) - known)
    assert not stale, (
        "Entries naming a class that no longer exists (renamed or removed?):\n  "
        + "\n  ".join(stale)
    )


def test_windows_and_embedded_do_not_overlap():
    overlap = sorted(WINDOWS & set(EMBEDDED))
    assert not overlap, f"Classified as both window and pane: {overlap}"


def test_embedded_reasons_are_written():
    blank = sorted(k for k, why in EMBEDDED.items() if not why.strip())
    assert not blank, f"EMBEDDED entries with no reason: {blank}"


# ── Live Qt (offscreen) ─────────────────────────────────────────────────────
#
# The tests above read source text, which proves the call is written, not that
# a key reaches anything. These build real widgets.

import pytest  # noqa: E402


def _qt_or_skip():
    try:
        from pyirena.gui._qt import QApplication
    except ImportError:
        pytest.skip("Qt (PySide6/PyQt6) not available")
    return QApplication.instance() or QApplication([])


def test_close_shortcut_closes_the_window():
    _qt_or_skip()
    from pyirena.gui._qt import QWidget
    from pyirena.gui.shortcuts import install_standard_shortcuts

    w = QWidget()
    w.show()
    shortcuts = install_standard_shortcuts(w)
    assert shortcuts, "no close shortcut was installed"
    shortcuts[0].activated.emit()
    assert not w.isVisible()


def test_close_binding_is_the_platform_one():
    """⌘W on macOS, Ctrl+W on Windows/Linux — from StandardKey, not a literal.

    Hard-coding "Ctrl+W" would give Mac users Control+W, which is a different
    key from the one the issue is about.
    """
    _qt_or_skip()
    from pyirena.gui._qt import QKeySequence, QWidget
    from pyirena.gui.shortcuts import install_standard_shortcuts

    expected = {k.toString() for k in
                QKeySequence.keyBindings(QKeySequence.StandardKey.Close)}
    assert expected, "Qt reports no Close binding on this platform"
    installed = {s.key().toString() for s in install_standard_shortcuts(QWidget())}
    assert installed == expected


def test_save_and_help_press_the_button_they_were_given():
    _qt_or_skip()
    from pyirena.gui._qt import QPushButton, QWidget
    from pyirena.gui.shortcuts import install_standard_shortcuts

    w = QWidget()
    pressed = []
    save_btn = QPushButton("Save", w)
    save_btn.clicked.connect(lambda: pressed.append("save"))
    shortcuts = install_standard_shortcuts(
        w, close=False, on_save=save_btn, on_help=lambda: pressed.append("help"),
    )
    for sc in shortcuts:
        sc.activated.emit()
    assert "save" in pressed and "help" in pressed


def test_a_disabled_button_does_nothing():
    """Going through click() means a disabled action stays disabled."""
    _qt_or_skip()
    from pyirena.gui._qt import QPushButton, QWidget
    from pyirena.gui.shortcuts import install_standard_shortcuts

    w = QWidget()
    pressed = []
    btn = QPushButton("Save", w)
    btn.setEnabled(False)
    btn.clicked.connect(lambda: pressed.append("save"))
    for sc in install_standard_shortcuts(w, close=False, on_save=btn):
        sc.activated.emit()
    assert pressed == []


def test_shortcuts_are_scoped_to_their_own_window():
    """Two open tool windows must not fight over ⌘W."""
    _qt_or_skip()
    from pyirena.gui._qt import Qt, QWidget
    from pyirena.gui.shortcuts import install_standard_shortcuts

    for sc in install_standard_shortcuts(QWidget()):
        assert sc.context() == Qt.ShortcutContext.WindowShortcut


def test_a_real_panel_gets_them():
    """End to end on a window that exists, not a bare QWidget."""
    _qt_or_skip()
    from pyirena.gui._qt import QKeySequence
    from pyirena.gui.contrast_panel import ContrastPanel

    panel = ContrastPanel()
    close_keys = {k.toString() for k in
                  QKeySequence.keyBindings(QKeySequence.StandardKey.Close)}
    help_keys = {k.toString() for k in
                 QKeySequence.keyBindings(QKeySequence.StandardKey.HelpContents)}
    installed = {s.key().toString() for s in panel.findChildren(_shortcut_type())}
    assert close_keys <= installed, f"no close shortcut: {installed}"
    assert help_keys <= installed, f"no help shortcut: {installed}"


def _shortcut_type():
    from pyirena.gui._qt import QShortcut
    return QShortcut
