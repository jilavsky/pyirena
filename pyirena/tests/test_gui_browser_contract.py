"""
Source-level guard for the shared file-browser behaviour.

The filename filter and the sort keys have each been duplicated across the four
file browsers at some point, and each time a fix landed in one and was missed in
the others.  These checks read source text only — no Qt, no display — so they
run in every CI job and fail the build if a fifth copy appears.

If this test fails on a browser you added:

    from pyirena.core.file_sorting import SORT_LABELS, SORT_TOOLTIP, sort_names
    from pyirena.gui.file_filter import FILTER_PLACEHOLDER, FILTER_TOOLTIP

    self.sort_combo.addItems(SORT_LABELS)
    self.sort_combo.setToolTip(SORT_TOOLTIP)
    ...
    files = sort_names(files, self.sort_combo.currentIndex())
"""

from __future__ import annotations

import re
from pathlib import Path

GUI_DIR = Path(__file__).resolve().parents[1] / "gui"
REPO = GUI_DIR.parent.parent

#: Re-export shims are allowed to name the private helpers.
_SHIMS = {"sorting.py"}

#: The four file browsers, all of which must share filter and sort behaviour.
BROWSERS = [
    "pyirena/gui/data_selector/panel.py",
    "pyirena/gui/hdf5viewer/file_tree.py",
    "pyirena/gui/data_manipulation_panel.py",
    "pyirena/gui/data_merge_panel.py",
]


def _gui_sources() -> list[Path]:
    return sorted(GUI_DIR.rglob("*.py"))


def test_no_module_redefines_a_sort_key():
    """The regexes live in pyirena.core.file_sorting and nowhere else."""
    pattern = re.compile(r"^def _?sort_key_(name|temperature|time|order|pressure)",
                         re.MULTILINE)
    offenders = [
        str(p.relative_to(REPO))
        for p in _gui_sources()
        if p.name not in _SHIMS and pattern.search(p.read_text(encoding="utf-8"))
    ]
    assert not offenders, (
        "These modules define their own filename sort keys; import them from "
        "pyirena.core.file_sorting instead:\n  " + "\n  ".join(offenders)
    )


def test_no_module_hard_codes_the_sort_labels():
    """One label list, so the dropdowns cannot drift apart again."""
    offenders = [
        str(p.relative_to(REPO))
        for p in _gui_sources()
        if re.search(r'["\']Temperature\s+↑', p.read_text(encoding="utf-8"))
    ]
    assert not offenders, (
        "These modules hard-code sort labels; use SORT_LABELS from "
        "pyirena.core.file_sorting:\n  " + "\n  ".join(offenders)
    )


def test_every_browser_uses_the_shared_filter():
    """Each file browser must offer the shared regex filter box."""
    missing = [
        rel for rel in BROWSERS
        if "FILTER_TOOLTIP" not in (REPO / rel).read_text(encoding="utf-8")
    ]
    assert not missing, (
        "These browsers do not use pyirena.gui.file_filter:\n  "
        + "\n  ".join(missing)
    )


def test_no_module_redefines_the_file_type_table():
    """The extensions per data-file type live in pyirena.core.file_types."""
    offenders = [
        str(p.relative_to(REPO))
        for p in _gui_sources()
        if re.search(r"^_FILE_TYPE_EXTS\s*=", p.read_text(encoding="utf-8"),
                     re.MULTILINE)
    ]
    assert not offenders, (
        "These modules keep their own file-type table; import FILE_TYPES / "
        "files_in_folder from pyirena.core.file_types:\n  " + "\n  ".join(offenders)
    )


def test_every_browser_uses_the_shared_sort_modes():
    """Each file browser must offer the same sort dropdown."""
    missing = []
    for rel in BROWSERS:
        text = (REPO / rel).read_text(encoding="utf-8")
        if "SORT_LABELS" not in text or "SORT_TOOLTIP" not in text:
            missing.append(rel)
    assert not missing, (
        "These browsers lack the shared sort dropdown (SORT_LABELS + "
        "SORT_TOOLTIP from pyirena.core.file_sorting):\n  " + "\n  ".join(missing)
    )


def test_every_browser_accepts_dropped_files():
    """Drag-and-drop is part of "a file browser", not a per-panel extra."""
    missing = [
        rel for rel in BROWSERS
        if "enable_file_drop" not in (REPO / rel).read_text(encoding="utf-8")
    ]
    assert not missing, (
        "These browsers do not accept dropped files; call enable_file_drop() "
        "from pyirena.gui.file_drop:\n  " + "\n  ".join(missing)
    )


def test_no_browser_lists_a_folder_by_hand():
    """One listing implementation, in pyirena.core.file_types.

    Every hand-rolled ``os.listdir`` loop answers the awkward questions its own
    way — sub-directories, case, an unreadable folder — and the four browsers
    disagreed about all three before this was shared.
    """
    pattern = re.compile(r"for \w+ in os\.listdir\(")
    offenders = []
    for rel in BROWSERS:
        text = (REPO / rel).read_text(encoding="utf-8")
        if pattern.search(text):
            offenders.append(rel)
    assert not offenders, (
        "These browsers walk a folder themselves; use files_in_folder() or "
        "files_with_extensions() from pyirena.core.file_types:\n  "
        + "\n  ".join(offenders)
    )


def test_no_browser_keeps_its_own_extension_list():
    """The extensions per data type live in one table.

    A private ``['.hdf', '.h5', '.hdf5']`` in a panel is how a browser ends up
    not listing the file type someone added last month.
    """
    pattern = re.compile(r"=\s*\[\s*['\"]\.(?:hdf|h5|txt|dat|csv)['\"]")
    offenders = []
    for rel in BROWSERS:
        text = (REPO / rel).read_text(encoding="utf-8")
        if pattern.search(text):
            offenders.append(rel)
    assert not offenders, (
        "These browsers hard-code file extensions; take them from "
        "pyirena.core.file_types (FILE_TYPE_EXTS / extensions_for):\n  "
        + "\n  ".join(offenders)
    )


def test_the_list_browsers_share_the_drop_selection_step():
    """"Switch folder, select what was dropped" is written once.

    The tree browser is excluded: its items are tree nodes carrying full paths,
    so it selects them its own way.
    """
    list_browsers = [rel for rel in BROWSERS if "file_tree" not in rel]
    missing = [
        rel for rel in list_browsers
        if "select_dropped_in_list" not in (REPO / rel).read_text(encoding="utf-8")
    ]
    assert not missing, (
        "These browsers re-implement the post-drop selection; use "
        "select_dropped_in_list() from pyirena.gui.file_drop:\n  "
        + "\n  ".join(missing)
    )
