"""
Contract for :mod:`pyirena.gui.raised_buttons` — the shared, vendored
raised-button stylesheet (see ``planning/raised-buttons/PLAN.md``).

The module is copied verbatim into Matilda, Bernardyn and AIDA, so its own
promises are tested here: it imports without Qt, every role pins a text
colour, the pressed state does not change the button's size, and a themed
button really paints a raised (top-light, bottom-dark) face.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

from pyirena.gui import raised_buttons as rb

_HEX = re.compile(r"^#[0-9a-fA-F]{6}$")
_RULE = re.compile(r"([^{}]*)\{([^{}]*)\}")


def test_module_imports_no_gui_toolkit():
    """The vendored file must drop into any package regardless of its shim."""
    tree = ast.parse(Path(rb.__file__).read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            names = [node.module or ""]
        else:
            continue
        for name in names:
            assert not re.match(r"(PySide|PyQt|pyirena)", name), (
                f"raised_buttons imports {name!r}; it must stay pure Python")


@pytest.mark.parametrize("roles", [rb.DEFAULT_ROLES, rb.PALETTE], ids=["roles", "palette"])
def test_every_role_has_face_and_text_colour(roles):
    for name, (bg, fg, bold) in roles.items():
        assert _HEX.match(bg) and _HEX.match(fg), name
        assert isinstance(bold, bool), name


def test_function_table_resolves_through_palette():
    """Each function names a palette colour; a typo must fail loudly, not
    silently fall back to a neutral button."""
    assert set(rb.DEFAULT_ROLES) == set(rb.FUNCTIONS)
    for role, colour in rb.FUNCTIONS.items():
        assert rb.DEFAULT_ROLES[role] == rb.PALETTE[colour]
    with pytest.raises(KeyError):
        rb.resolve_roles({"primary": "gren"})
    assert rb.resolve_roles({"x": "#123456"})["x"][0] == "#123456"


def test_recolouring_a_function_is_one_table_entry():
    css = rb.stylesheet(rb.resolve_roles({**rb.FUNCTIONS, "help": "blue"}))
    help_rule = re.search(r'QPushButton\[btnRole="help"\] \{([^}]*)\}', css).group(1)
    assert rb.PALETTE["blue"][0] in help_rule


def test_theme_rejects_unknown_role():
    pytest.importorskip("PySide6", reason="Qt binding required")
    from pyirena.gui import theme
    from pyirena.gui._qt import QApplication, QPushButton

    QApplication.instance() or QApplication([])
    with pytest.raises(KeyError):
        theme.set_button_role(QPushButton("x"), "green")   # a colour, not a function


def test_theme_roles_cover_the_shared_role_names():
    """pyIrena may recolour roles, but every shared name must exist, so a
    ``set_role(btn, 'teal')`` written in any package works in every package."""
    pytest.importorskip("PySide6", reason="theme imports the Qt shim")
    from pyirena.gui import theme

    assert set(rb.DEFAULT_ROLES) <= set(theme.BUTTON_ROLES)


def test_every_background_rule_pins_text_colour():
    css = rb.stylesheet()
    for selector, decl in _RULE.findall(css):
        if re.search(r"(?<![-\w])background\s*:", decl):
            assert re.search(r"(?<![-\w])color\s*:", decl), selector.strip()


@pytest.mark.parametrize("size", list(rb.SIZES))
def test_pressed_state_keeps_outer_size(size):
    """The label sinks one pixel; total vertical padding must not change."""
    pad_v, pad_h, font = rb.SIZES[size]
    css = rb.size_qss("QPushButton", pad_v, pad_h, font)
    pressed = re.search(r":pressed[^{]*\{\s*padding:\s*(\d+)px \d+px (\d+)px", css)
    assert pressed, css
    top, bottom = int(pressed.group(1)), int(pressed.group(2))
    if pad_v >= 1:
        assert top + bottom == 2 * pad_v


def test_colour_arithmetic():
    assert rb.lighten("#000000", 1.0) == "#ffffff"
    assert rb.darken("#ffffff", 1.0) == "#000000"
    assert rb.mix("#000000", "#ffffff", 0.5) == "#808080"
    assert rb.mix("#abc", "#abc", 0.3) == "#aabbcc"
    with pytest.raises(ValueError):
        rb.mix("#12345", "#000000", 0.5)


def _render(role, *, pressed=False, enabled=True):
    pytest.importorskip("PySide6", reason="Qt binding required to render")
    from pyirena.gui import theme
    from pyirena.gui._qt import QApplication, QColor, QPushButton

    app = QApplication.instance() or QApplication([])
    theme.apply_theme(app)
    btn = QPushButton("Fit")
    theme.set_button_role(btn, role)
    btn.setEnabled(enabled)
    btn.setDown(pressed)
    btn.resize(btn.sizeHint().width() + 20, btn.sizeHint().height())
    img = btn.grab().toImage()
    x = img.width() // 2 - 15  # left of the label

    def lum(y):
        c = QColor(img.pixel(x, y))
        return 0.299 * c.red() + 0.587 * c.green() + 0.114 * c.blue()

    h = img.height()
    return lum(3), lum(h - 4)


@pytest.mark.parametrize("role", ["neutral", "primary", "setup"])
def test_button_face_is_raised_and_pressed_is_sunken(role):
    top, bottom = _render(role)
    assert top > bottom + 8, f"{role}: face is not lit from above ({top:.0f} vs {bottom:.0f})"
    ptop, pbottom = _render(role, pressed=True)
    assert ptop < pbottom, f"{role}: pressed face is not sunken ({ptop:.0f} vs {pbottom:.0f})"


def test_disabled_button_is_flat():
    top, bottom = _render("primary", enabled=False)
    assert abs(top - bottom) < 3, "disabled buttons should lose the bevel"
