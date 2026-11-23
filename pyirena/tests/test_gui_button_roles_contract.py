"""
Buttons are coloured by function, through one table — keep it that way.

pyIrena's push buttons get their look from the application stylesheet that
:func:`pyirena.gui.theme.apply_theme` builds from ``BUTTON_FUNCTIONS``:
panels call ``set_button_role(btn, "primary")`` and never name a colour. That
is what lets one line in ``theme.py`` recolour every button of a function, and
what keeps a per-widget stylesheet from bringing the old flat rectangle back
(see ``planning/raised-buttons/PLAN.md``).

Source text only; no display needed.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import pyirena

GUI_DIR = Path(pyirena.__file__).resolve().parent / "gui"

#: Files that define the button look itself.
_DEFINERS = {"theme.py", "raised_buttons.py"}

#: Colour-swatch buttons: the background *is* the value shown (a line or
#: marker colour), so it cannot come from a role.
SWATCH_FILES: set[str] = {
    "data_selector/config_dialogs.py",
    "diffraction_lines_panel.py",
    "hdf5viewer/graph_window.py",
}

#: A button selector followed, within the same rule, by a background.
_BUTTON_BG = re.compile(
    r"Q(?:Push|Tool)Button[^{}\"']*\{[^}]*?(?<![-\w])background(?:-color)?\s*:", re.S)
#: Per-widget helpers removed in favour of set_button_role; keep them gone.
_OLD_HELPER = re.compile(r"\b(accent_button_css|soft_button_css|chip_button_css|CHIP_BUTTON_CSS)\b")
#: set_button_role(<widget>, "<role>"
_ROLE_CALL = re.compile(r"set_button_role\(\s*[^,()]+,\s*['\"]([^'\"]+)['\"]")


def _gui_files():
    for path in sorted(GUI_DIR.rglob("*.py")):
        rel = path.relative_to(GUI_DIR).as_posix()
        if rel in _DEFINERS:
            continue
        yield rel, path.read_text(encoding="utf-8")


def test_no_inline_button_background():
    offenders = []
    for rel, text in _gui_files():
        if rel in SWATCH_FILES:
            continue
        for m in _BUTTON_BG.finditer(text):
            line = text.count("\n", 0, m.start()) + 1
            offenders.append(f"{rel}:{line}")
    assert not offenders, (
        "Inline button stylesheet with a background — use "
        "theme.set_button_role(btn, '<function>') instead:\n  " + "\n  ".join(offenders))


def test_old_button_css_helpers_unused():
    offenders = [rel for rel, text in _gui_files() if _OLD_HELPER.search(text)]
    assert not offenders, (
        "accent/soft/chip_button_css are superseded by set_button_role: "
        + ", ".join(offenders))


def test_every_role_literal_is_a_known_function():
    pytest.importorskip("PySide6", reason="theme imports the Qt shim")
    from pyirena.gui.theme import BUTTON_ROLES

    unknown = []
    for rel, text in _gui_files():
        for m in _ROLE_CALL.finditer(text):
            if m.group(1) not in BUTTON_ROLES:
                line = text.count("\n", 0, m.start()) + 1
                unknown.append(f"{rel}:{line} {m.group(1)!r}")
    assert not unknown, "Unknown button role(s):\n  " + "\n  ".join(unknown)
