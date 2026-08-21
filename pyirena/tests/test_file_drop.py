"""
Tests for :mod:`pyirena.gui.file_drop` — opening data by dropping it.

The part worth testing is not the Qt plumbing but what a drop *means*: a folder
should become the files inside it, anything pyIrena cannot open should be
dropped silently, and the result should arrive in the order the file browsers
would show it, because it is fed straight into a file list.

:func:`collect_dropped_paths` is pure, so all of that is checked with real files
in a tmp directory and no display.  The Qt wiring gets one smoke test.
"""

from __future__ import annotations

import os

import pytest

from pyirena.gui.file_drop import (
    DATA_EXTENSIONS,
    collect_dropped_paths,
    drop_hint,
    first_folder,
)


def _make(tmp_path, *names):
    """Create empty files and return their absolute paths, in the given order."""
    out = []
    for name in names:
        p = tmp_path / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("")
        out.append(str(p))
    return out


# ── What comes back ──────────────────────────────────────────────────────

def test_dropped_files_come_back_as_absolute_paths(tmp_path):
    made = _make(tmp_path, "a.h5", "b.dat")
    assert collect_dropped_paths(made) == sorted(made)


def test_files_pyirena_cannot_open_are_ignored(tmp_path):
    good, _bad = _make(tmp_path, "data.h5", "notes.pdf")
    got = collect_dropped_paths([good, str(tmp_path / "notes.pdf")])
    assert got == [good]


def test_dropping_only_unopenable_files_yields_nothing(tmp_path):
    _make(tmp_path, "notes.pdf", "image.png")
    got = collect_dropped_paths([str(tmp_path / "notes.pdf"), str(tmp_path / "image.png")])
    assert got == []


def test_a_path_that_does_not_exist_is_ignored(tmp_path):
    real, = _make(tmp_path, "real.h5")
    got = collect_dropped_paths([real, str(tmp_path / "ghost.h5")])
    assert got == [real]


def test_empty_and_blank_entries_are_ignored():
    assert collect_dropped_paths(["", None if False else ""]) == []


def test_an_empty_extension_list_accepts_anything(tmp_path):
    made = _make(tmp_path, "notes.pdf")
    assert collect_dropped_paths(made, extensions=()) == made


def test_a_panel_can_narrow_the_accepted_extensions(tmp_path):
    h5, dat = _make(tmp_path, "a.h5", "b.dat")
    assert collect_dropped_paths([h5, dat], extensions=(".h5",)) == [h5]


def test_extension_matching_ignores_case(tmp_path):
    made = _make(tmp_path, "LOUD.H5")
    assert collect_dropped_paths(made) == made


# ── Folders ──────────────────────────────────────────────────────────────

def test_dropping_a_folder_yields_the_data_files_inside(tmp_path):
    _make(tmp_path, "run/one.h5", "run/two.h5", "run/readme.pdf")
    got = collect_dropped_paths([str(tmp_path / "run")])
    assert [os.path.basename(p) for p in got] == ["one.h5", "two.h5"]


def test_a_dropped_folder_does_not_reach_into_subfolders_by_default(tmp_path):
    _make(tmp_path, "run/top.h5", "run/deeper/inner.h5")
    got = collect_dropped_paths([str(tmp_path / "run")])
    assert [os.path.basename(p) for p in got] == ["top.h5"]


def test_a_caller_can_ask_for_one_more_level(tmp_path):
    _make(tmp_path, "run/top.h5", "run/deeper/inner.h5")
    got = collect_dropped_paths([str(tmp_path / "run")], folder_depth=2)
    assert sorted(os.path.basename(p) for p in got) == ["inner.h5", "top.h5"]


def test_an_empty_folder_yields_nothing(tmp_path):
    (tmp_path / "empty").mkdir()
    assert collect_dropped_paths([str(tmp_path / "empty")]) == []


# ── Order and duplicates ─────────────────────────────────────────────────

def test_the_same_file_dropped_twice_appears_once(tmp_path):
    one, = _make(tmp_path, "a.h5")
    assert collect_dropped_paths([one, one]) == [one]


def test_a_file_and_its_folder_together_do_not_duplicate_it(tmp_path):
    _make(tmp_path, "run/a.h5")
    inside = str(tmp_path / "run" / "a.h5")
    got = collect_dropped_paths([inside, str(tmp_path / "run")])
    assert got == [inside]


def test_results_are_ordered_the_way_the_file_browsers_order_them(tmp_path):
    # The browsers' default is the case-insensitive filename order from
    # core.file_sorting, not the order the file manager happened to send.
    _make(tmp_path, "run/Beta.h5", "run/alpha.h5", "run/Gamma.h5")
    got = collect_dropped_paths([str(tmp_path / "run")])
    assert [os.path.basename(p) for p in got] == ["alpha.h5", "Beta.h5", "Gamma.h5"]


def test_files_of_the_same_name_from_two_folders_are_both_kept(tmp_path):
    _make(tmp_path, "a/scan.h5", "b/scan.h5")
    got = collect_dropped_paths([str(tmp_path / "a" / "scan.h5"), str(tmp_path / "b" / "scan.h5")])
    assert len(got) == 2


# ── Small helpers ────────────────────────────────────────────────────────

def test_first_folder_is_where_a_browser_should_go(tmp_path):
    made = _make(tmp_path, "run/a.h5", "run/b.h5")
    assert first_folder(made) == str(tmp_path / "run")


def test_first_folder_of_nothing_is_none():
    assert first_folder([]) is None


def test_the_tooltip_hint_mentions_dragging():
    assert "drag" in drop_hint("a data file").lower()


def test_the_default_extensions_cover_hdf5_and_text():
    assert ".h5" in DATA_EXTENSIONS and ".dat" in DATA_EXTENSIONS


# ── Qt wiring (smoke) ────────────────────────────────────────────────────

def _require_qt():
    try:
        from pyirena.gui._qt import QApplication  # noqa: F401
    except Exception:
        pytest.skip("Qt not available", allow_module_level=True)


def test_enabling_drop_marks_the_widget_and_survives_repeat_calls(tmp_path):
    _require_qt()
    from pyirena.gui._qt import QApplication, QWidget
    from pyirena.gui.file_drop import enable_file_drop

    app = QApplication.instance() or QApplication([])
    assert app is not None
    widget = QWidget()
    try:
        got = []
        assert enable_file_drop(widget, got.append) is not None
        assert widget.acceptDrops()
        # Re-wiring (a panel rebuilt after a settings change) must not raise.
        assert enable_file_drop(widget, got.append) is not None
    finally:
        widget.deleteLater()


def test_every_file_browser_accepts_drops_and_selects_what_was_dropped(tmp_path):
    """The four drop targets a user actually aims at must all be live.

    Wiring is one call per panel and easy to lose in a refactor, so this checks
    the behaviour end to end rather than the presence of the call: each browser
    must switch to the dropped file's folder and leave it selected.
    """
    _require_qt()
    pytest.importorskip("pyqtgraph")
    import h5py

    from pyirena.gui._qt import QApplication

    app = QApplication.instance() or QApplication([])
    assert app is not None

    folder = tmp_path / "run"
    folder.mkdir()
    for name in ("a_scan.h5", "b_scan.h5"):
        h5py.File(folder / name, "w").close()
    dropped = [str(folder / "a_scan.h5")]

    from pyirena.gui.data_manipulation_panel import _ManipFileBrowser
    from pyirena.gui.data_merge_panel import _DatasetSelectorWidget
    from pyirena.gui.hdf5viewer.file_tree import FileTreeWidget

    browsers = [
        (FileTreeWidget(), lambda w: [os.path.basename(p) for p in w.get_selected_paths()]),
        (_DatasetSelectorWidget(1), lambda w: w.get_selected_filenames()),
        (_ManipFileBrowser(), lambda w: w.get_selected_filenames()),
    ]
    try:
        for widget, selection in browsers:
            assert widget.acceptDrops(), f"{type(widget).__name__} does not accept drops"
            widget.open_dropped_paths(dropped)
            assert selection(widget) == ["a_scan.h5"], type(widget).__name__
    finally:
        for widget, _ in browsers:
            widget.deleteLater()


def test_paths_from_mime_reads_both_urls_and_plain_text(tmp_path):
    _require_qt()
    from pyirena.gui._qt import QApplication, QMimeData, QUrl
    from pyirena.gui.file_drop import paths_from_mime

    app = QApplication.instance() or QApplication([])
    assert app is not None
    one, = _make(tmp_path, "a.h5")

    mime = QMimeData()
    mime.setUrls([QUrl.fromLocalFile(one)])
    assert paths_from_mime(mime) == [one]

    # Applications that send text instead of URLs are handled too.
    text_mime = QMimeData()
    text_mime.setText(f"file://{one}\n")
    assert paths_from_mime(text_mime) == [one]


def test_text_drop_payloads_normalise_to_native_paths():
    """Every spelling of a dropped path must come back with native separators.

    Windows CI caught this: ``QUrl.toLocalFile()`` returns ``C:/Users/...``
    with forward slashes, and Explorer sends ``file:///C:/Users/...``.  Both
    used to reach the panels unnormalised, where they compare unequal to the
    same file found by a folder scan.  The Windows half is exercised here by
    swapping in ntpath, so a POSIX-only CI run still catches a regression.
    """
    import ntpath
    import nturl2path

    from pyirena.gui import file_drop as fd

    posix_cases = {
        "/data/a.h5": "/data/a.h5",
        "file:///data/a.h5": "/data/a.h5",
        "file:///data/my%20file.h5": "/data/my file.h5",
        "": "",
    }
    for payload, want in posix_cases.items():
        assert fd._path_from_text_line(payload) == want, payload

    real_ospath, real_u2p = fd.os.path, fd.url2pathname
    fd.os.path, fd.url2pathname = ntpath, nturl2path.url2pathname
    try:
        win_cases = {
            # Explorer / Qt: a proper file URL with an empty authority.
            "file:///C:/Users/me/a.h5": r"C:\Users\me\a.h5",
            # Applications that drop the third slash send a bare path.
            r"file://C:\Users\me\a.h5": r"C:\Users\me\a.h5",
            # Already a path, but with the separators Qt prefers.
            "C:/Users/me/a.h5": r"C:\Users\me\a.h5",
            "file:///C:/Users/me/my%20file.h5": r"C:\Users\me\my file.h5",
        }
        for payload, want in win_cases.items():
            assert fd._path_from_text_line(payload) == want, payload
    finally:
        fd.os.path, fd.url2pathname = real_ospath, real_u2p
