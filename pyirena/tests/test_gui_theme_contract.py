"""
Contrast contract for the GUI's inline stylesheets.

The bug this guards against, reported from a beamline user's machine:
``pyirena/gui`` carries a few hundred inline Qt stylesheets written for a light
colour scheme, and Qt overrides only the properties a stylesheet actually
names.  Under a *dark* system scheme the unnamed half came from the platform
palette — so a rule that set only ``background-color: #ecf0f1`` painted the
desktop's near-white text on our near-white box.  The Modeling panel's
"Fit B/P btwn cursors" and "Fit Flat btwn cursors" buttons became blank
rectangles the user never realised were buttons.

Two defences, one test module:

1. :func:`test_theme_pins_both_colours_for_every_control` — the shipped
   baseline stylesheet must name a text colour wherever it names a background,
   because that stylesheet is what makes pyIrena look the same on every
   platform (``pyirena.gui.theme.apply_theme``).
2. :func:`test_no_inline_background_without_text_colour` — no *new* inline
   stylesheet may set a background without also setting ``color``.  Existing
   offenders are listed in ``KNOWN_BARE_BACKGROUNDS``; the list may shrink,
   never grow.

Source text only where possible: no display needed.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import pyirena

GUI_DIR = Path(pyirena.__file__).resolve().parent / "gui"

#: ``background`` / ``background-color`` declaration.
_BACKGROUND = re.compile(r"(?<![-\w])background(?:-color)?\s*:")
#: ``color:`` but not ``background-color:`` / ``border-color:`` / ``selection-color:``.
_TEXT_COLOUR = re.compile(r"(?<![-\w])color\s*:")
#: A QSS rule block: ``SELECTOR { declarations }``.
_RULE = re.compile(r"([^{}]*)\{([^{}]*)\}")
#: A ``/* ... */`` comment.  Stripped before rules are split, or a comment
#: sitting above a rule is swallowed into that rule's selector — which once
#: let a ``::indicator`` mentioned in prose exempt the rule below it.
_COMMENT = re.compile(r"/\*.*?\*/", re.DOTALL)

#: Selectors that paint a surface which never renders text — a scrollbar
#: groove, a splitter grip, a progress-bar chunk.  There is no text colour to
#: pair with, so requiring one would be noise.
TEXTLESS_SELECTORS = re.compile(
    r"(QScrollBar|QSplitter)\b"
    r"|::(handle|chunk|groove|add-line|sub-line|add-page|sub-page"
    r"|drop-down|up-button|down-button|indicator)\b"
)

#: Files whose inline styles set a background with no paired text colour, and
#: which the forced light palette already makes safe.  Shrink this list when
#: you convert one to :mod:`pyirena.gui.theme`; never add to it.
KNOWN_BARE_BACKGROUNDS: set[str] = {
    # colour-swatch buttons: the background *is* the value being shown
    "data_selector/config_dialogs.py",
    "diffraction_lines_panel.py",
    "hdf5viewer/graph_window.py",
    # tinted section headers, safe under the shipped palette
    "hdf5viewer/file_tree.py",
    "hdf5viewer/hdf5_browser.py",
}


def _qss_rules(css: str):
    """Yield ``(selector, declarations)`` for *css*, treating a bare
    declaration list as one rule with an empty selector."""
    css = _COMMENT.sub(" ", css)
    blocks = _RULE.findall(css)
    return blocks or [("", css)]


def _bare_background_selectors(css: str) -> list[str]:
    """Selectors in *css* that set a background with no text colour in force.

    A pseudo-state rule (``QPushButton:hover``) inherits ``color`` from the
    base rule for the same selector, so both are considered together.
    """
    rules = _qss_rules(css)
    with_colour = {
        sel.split(":")[0].strip()
        for sel, decl in rules
        if _TEXT_COLOUR.search(decl)
    }
    bare = []
    for sel, decl in rules:
        if not _BACKGROUND.search(decl):
            continue
        if _TEXT_COLOUR.search(decl):
            continue
        if sel.split(":")[0].strip() in with_colour:
            continue
        if TEXTLESS_SELECTORS.search(sel):
            continue
        bare.append(sel.strip() or "<bare declarations>")
    return bare


def test_theme_pins_both_colours_for_every_control():
    """The shipped baseline stylesheet never sets a background alone."""
    pytest.importorskip("PySide6", reason="Qt binding required to build the theme")
    from pyirena.gui import theme

    bare = _bare_background_selectors(theme._base_stylesheet())
    assert not bare, (
        "pyirena.gui.theme._base_stylesheet() sets a background with no text "
        "colour for: " + ", ".join(bare) + ".\nUnder a dark system scheme the "
        "text colour would come from the platform palette and the control "
        "would be unreadable — name `color:` in the same rule."
    )


def test_theme_helpers_pin_both_colours():
    """Every CSS helper in :mod:`pyirena.gui.theme` names both colours."""
    pytest.importorskip("PySide6", reason="Qt binding required to build the theme")
    from pyirena.gui import theme

    helpers = {
        "CHIP_BUTTON_CSS": theme.CHIP_BUTTON_CSS,
        "READONLY_FIELD_CSS": theme.READONLY_FIELD_CSS,
        "accent_button_css": theme.accent_button_css(theme.ACCENT_GREEN),
        "soft_button_css": theme.soft_button_css(theme.SOFT_GREEN),
        "readout_label_css": theme.readout_label_css(),
        "status_css": theme.status_css("warning"),
    }
    for name, css in helpers.items():
        assert not _bare_background_selectors(css), (
            f"theme.{name} sets a background without a text colour: {css}"
        )


def _inline_stylesheets(path: Path):
    """Yield ``(line_number, css)`` for each ``setStyleSheet(...)`` literal.

    Only string literals passed directly are inspected — that is where the
    hardcoded colours live.  Styles built into a variable first are caught by
    the theme-helper tests instead.
    """
    src = path.read_text(encoding="utf-8")
    for m in re.finditer(r"setStyleSheet\s*\(", src):
        i, depth, j = m.end(), 1, m.end()
        while j < len(src) and depth:
            ch = src[j]
            if ch in "([{":
                depth += 1
            elif ch in ")]}":
                depth -= 1
            elif ch in "\"'":
                quote = src[j:j + 3] if src[j:j + 3] in ('"""', "'''") else ch
                j += len(quote)
                while j < len(src) and src[j:j + len(quote)] != quote:
                    j += 2 if src[j] == "\\" else 1
                j += len(quote) - 1
            j += 1
        arg = src[i:j - 1]
        literals = re.findall(
            r'"""((?:.|\n)*?)"""|\'\'\'((?:.|\n)*?)\'\'\''
            r'|"((?:[^"\\\n]|\\.)*)"|\'((?:[^\'\\\n]|\\.)*)\'',
            arg,
        )
        css = "\n".join(part for group in literals for part in group if part)
        if css:
            yield src[: m.start()].count("\n") + 1, css


def test_no_inline_background_without_text_colour():
    """New inline styles must set ``color`` wherever they set a background."""
    offenders: dict[str, list[str]] = {}
    for path in sorted(GUI_DIR.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        rel = str(path.relative_to(GUI_DIR))
        # theme.py is the one place that *defines* these strings; its own
        # helpers are checked by test_theme_helpers_pin_both_colours, and its
        # docstrings quote the old broken styles on purpose.
        if rel == "theme.py" or rel in KNOWN_BARE_BACKGROUNDS:
            continue
        for line, css in _inline_stylesheets(path):
            for sel in _bare_background_selectors(css):
                offenders.setdefault(rel, []).append(f"line {line}: {sel}")

    assert not offenders, (
        "Inline stylesheets set a background with no text colour:\n"
        + "\n".join(f"  {name}: {'; '.join(hits)}" for name, hits in offenders.items())
        + "\nA background-only rule inherits the *system* text colour, which is "
        "light on a dark desktop — the control disappears.  Use a helper from "
        "pyirena.gui.theme, or name `color:` alongside the background."
    )


# ---------------------------------------------------------------------------
# The indicator actually reaches the screen
# ---------------------------------------------------------------------------
# A source-text check cannot see this one.  ``QCheckBox { background-color:
# transparent }`` looked harmless and passed every rule above, but Qt copies a
# QSS background into the palette it hands the base style, and Fusion derives
# the indicator's outline from ``palette.window()`` — so the box around every
# "Fit?" check was painted in a fully transparent pen and simply was not there.
# These tests render a check box and look at the pixels.

#: WCAG minimum contrast for a user-interface component against its surround.
_MIN_COMPONENT_CONTRAST = 3.0


def _relative_luminance(hex_colour: str) -> float:
    """WCAG relative luminance of ``#rrggbb``."""
    channels = []
    for i in (1, 3, 5):
        c = int(hex_colour[i:i + 2], 16) / 255.0
        channels.append(c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4)
    r, g, b = channels
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _contrast(a: str, b: str) -> float:
    """WCAG contrast ratio between two ``#rrggbb`` colours."""
    la, lb = _relative_luminance(a), _relative_luminance(b)
    lo, hi = sorted((la, lb))
    return (hi + 0.05) / (lo + 0.05)


def _render_check_box(checked: bool):
    """Return the set of ``#rrggbb`` colours a themed check box paints."""
    pytest.importorskip("PySide6", reason="Qt binding required to render")
    from pyirena.gui import theme
    from pyirena.gui._qt import QApplication, QCheckBox, QColor

    app = QApplication.instance() or QApplication([])
    theme.apply_theme(app)
    box = QCheckBox("")
    box.setChecked(checked)
    box.resize(box.sizeHint())
    image = box.grab().toImage()
    return {
        QColor(image.pixel(x, y)).name()
        for y in range(image.height())
        for x in range(image.width())
    }


def test_theme_never_makes_selection_controls_transparent():
    """``background-color`` on QCheckBox/QRadioButton erases the indicator.

    Qt writes a QSS background into ``QPalette::Window`` and ``::Base`` for the
    base style, and Fusion builds the indicator's outline and interior out of
    exactly those two roles.  Style ``::indicator`` instead.
    """
    pytest.importorskip("PySide6", reason="Qt binding required to build the theme")
    from pyirena.gui import theme

    for sel, decl in _qss_rules(theme._base_stylesheet()):
        if "::" in sel:  # ::indicator is styled on purpose
            continue
        if re.search(r"\bQ(CheckBox|RadioButton)\b", sel):
            assert not _BACKGROUND.search(decl), (
                f"`{sel.strip()}` sets a background; that reaches Fusion as a "
                "palette override and the check indicator disappears."
            )


def test_check_box_indicator_is_visible():
    """An unchecked box is drawn, with an outline the user can see."""
    from pyirena.gui import theme

    painted = _render_check_box(checked=False)
    assert theme.BASE_BG in painted, (
        "the check box interior was never painted — the indicator is missing, "
        "not merely faint"
    )
    best = max(_contrast(c, theme.WINDOW_BG) for c in painted)
    assert best >= _MIN_COMPONENT_CONTRAST, (
        f"the strongest colour in an unchecked check box reaches only "
        f"{best:.2f}:1 against the panel background; a hairline like that "
        f"rounds away on a fractionally scaled Windows or Linux display"
    )


def test_checked_box_is_distinguishable_from_unchecked():
    """Checked and unchecked must differ by more than a thin glyph."""
    from pyirena.gui import theme

    off = _render_check_box(checked=False)
    on = _render_check_box(checked=True)
    assert theme.ACCENT_BLUE in on and theme.ACCENT_BLUE not in off, (
        "the checked state is not filled with the accent colour"
    )
    assert _contrast(theme.ACCENT_BLUE, theme.BASE_BG) >= _MIN_COMPONENT_CONTRAST
