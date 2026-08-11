"""
Source-level guard for the "every table is copyable" contract.

This is the check the feature-parity review asked for
(``grep -rn "QTableWidget(" pyirena/gui`` — confirm each hit uses the shared
helper), automated so the next table added to the GUI cannot quietly ship
without clipboard support.  It reads source text only: no Qt, no display, so it
runs in every CI job.

If this test fails on a table you added, the fix is one line:

    from pyirena.gui.table_utils import attach_table_copy
    ...
    attach_table_copy(my_table)          # + enable_table_sorting(...) if the
                                         #   row order carries no meaning

If a table genuinely must not be copyable, add its module to
``_EXEMPT_MODULES`` below with a comment saying why.
"""

from __future__ import annotations

import re
from pathlib import Path

GUI_DIR = Path(__file__).resolve().parents[1] / "gui"

# Modules that construct a QTableWidget but deliberately do not attach the
# shared copy behaviour.  Keep this list empty if at all possible.
_EXEMPT_MODULES: dict[str, str] = {}

# ``_qt.py`` re-exports the Qt names; it does not build tables.
_SHIM_MODULES = {"_qt.py"}


def _modules_constructing_tables() -> list[Path]:
    out = []
    for path in sorted(GUI_DIR.rglob("*.py")):
        if path.name in _SHIM_MODULES:
            continue
        text = path.read_text(encoding="utf-8")
        if re.search(r"\bQTableWidget\s*\(", text):
            out.append(path)
    return out


def test_every_table_module_uses_the_shared_helper():
    offenders = []
    for path in _modules_constructing_tables():
        if path.name in _EXEMPT_MODULES:
            continue
        text = path.read_text(encoding="utf-8")
        if "attach_table_copy" not in text:
            offenders.append(str(path.relative_to(GUI_DIR.parent.parent)))
    assert not offenders, (
        "These modules create a QTableWidget without calling "
        "attach_table_copy() from pyirena.gui.table_utils:\n  "
        + "\n  ".join(offenders)
    )


def test_no_bespoke_clipboard_copy_left_in_table_modules():
    """Table modules must not hand-roll ``clipboard().setText`` for a table.

    Two divergent implementations is how this feature drifted the first time.
    """
    offenders = []
    for path in _modules_constructing_tables():
        text = path.read_text(encoding="utf-8")
        if "clipboard()" in text:
            offenders.append(str(path.relative_to(GUI_DIR.parent.parent)))
    assert not offenders, (
        "Bespoke clipboard code found in a table module; use "
        "pyirena.gui.table_utils.copy_table_selection instead:\n  "
        + "\n  ".join(offenders)
    )


def test_helper_module_exists_and_exports_the_public_api():
    src = (GUI_DIR / "table_utils.py").read_text(encoding="utf-8")
    for name in (
        "attach_table_copy",
        "enable_table_sorting",
        "rows_to_csv_text",
        "save_rows_as_csv",
        "make_numeric_item",
        "populating",
    ):
        assert f'"{name}"' in src, f"{name} missing from table_utils.__all__"
