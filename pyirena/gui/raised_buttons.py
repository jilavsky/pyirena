"""
raised_buttons — one raised, Igor-like push-button look for every Qt package.

VENDORED FILE — keep byte-identical across packages.
    Source of truth: pyirena/gui/raised_buttons.py
    Copies:          matilda/gui/, bernardyn/gui/, aida/ui/qt/, ...
    Edit the source, bump ``__version__``, then run the sync script.

Why this exists
---------------
Any Qt stylesheet that sets a button's ``background`` (and, under Fusion, any
``border``) makes Qt stop asking the platform style to draw the button. The
result is the flat coloured rectangle that is hard to recognise as a button.
This module puts the 3-D cues back *inside* the stylesheet:

* a vertical gradient (light top, darker bottom) — the raised face,
* a two-tone border (lighter top/left, darker bottom/right) — the bevel,
* a pressed state with the gradient and bevel inverted and the label shifted
  down one pixel — the sunken look,
* a flat, desaturated disabled state — "not available" reads at a glance.

Design rules
------------
* **No Qt import.** Everything here builds strings. The only function that
  touches a widget (:func:`set_role`) does so by duck typing, so the file
  drops into a package regardless of its ``_qt.py`` shim.
* **Performance.** Build the stylesheet once and install it once on the
  ``QApplication``. Buttons pick a colour with a dynamic property
  (``btnRole``), not their own ``setStyleSheet`` — one parsed stylesheet
  instead of hundreds. Gradients are painted by QPainter and cost nothing
  measurable. Never use ``QGraphicsDropShadowEffect`` for buttons: it renders
  each widget off-screen on every repaint.
* **Never write ``QWidget { background: ... }``** in a host stylesheet: it
  matches every button inside that widget and flattens them all again. Set a
  panel's background through its ``QPalette`` or a narrower selector.

Usage
-----
App level (once, right after ``QApplication``)::

    app.setStyle("Fusion")          # recommended, not required
    app.setStyleSheet(host_qss + raised_buttons.stylesheet())

Per button — name what the button *does*, never its colour::

    raised_buttons.set_role(fit_button, "primary", size="large")
    raised_buttons.set_role(cancel_button, "stop")
    raised_buttons.set_role(helper_button, None, size="small")   # neutral

Colours are looked up once, when the stylesheet is built: :data:`FUNCTIONS`
maps each function to a :data:`PALETTE` colour. Recolouring every "Revert"
button in a package is a one-line change to that table (or to the host's
copy of it, see :func:`resolve_roles`) — no button code changes.

Fallback for a widget that lives in an app which did not install the
stylesheet::

    btn.setStyleSheet(raised_buttons.inline_qss("revert"))
"""

from __future__ import annotations

__version__ = "0.2.0"

__all__ = [
    "DEFAULT_BEVEL", "DEFAULT_RADIUS", "ROLE_PROPERTY", "SIZE_PROPERTY",
    "PALETTE", "FUNCTIONS", "DEFAULT_ROLES", "SIZES",
    "resolve_roles", "button_qss", "size_qss", "stylesheet", "inline_qss",
    "set_role", "get_role",
    "lighten", "darken", "mix",
]

#: Depth of the bevel, 0 (flat) … ~0.45 (strong). Chosen by Jan 2026-11-23.
DEFAULT_BEVEL = 0.35
#: Corner radius in px. 7 px reads clearly as "button"; chosen 2026-11-23.
DEFAULT_RADIUS = 7

#: Dynamic property selecting the colour role, ``QPushButton[btnRole="green"]``.
ROLE_PROPERTY = "btnRole"
#: Dynamic property selecting the size, ``QPushButton[btnSize="small"]``.
SIZE_PROPERTY = "btnSize"

#: colour name -> (face colour, text colour, bold label). The shared button
#: palette: every package draws "green" with the same green.
PALETTE: dict[str, tuple[str, str, bool]] = {
    "neutral":    ("#e6eaed", "#20303c", False),
    "green":      ("#27ae60", "#ffffff", True),
    "mint":       ("#52c77a", "#ffffff", True),
    "blue":       ("#2980b9", "#ffffff", True),
    "teal":       ("#16a085", "#ffffff", True),
    "purple":     ("#8e44ad", "#ffffff", True),
    "orange":     ("#e67e22", "#ffffff", True),
    "red":        ("#c0392b", "#ffffff", True),
    "grey":       ("#7f8c8d", "#ffffff", True),
    "soft-green": ("#cdeccd", "#20303c", False),
    "soft-amber": ("#ffe082", "#20303c", False),
}

#: button function -> colour name. The shared vocabulary: a button is given
#: one of these roles by what it *does*, so the same action has the same
#: colour in every tool and every package. A host may add functions or
#: re-point one (see :func:`resolve_roles`); keep the names.
FUNCTIONS: dict[str, str] = {
    "neutral":  "neutral",     # utility, navigation, dialogs: Open…, Browse, All/None
    "primary":  "green",       # the main action: Fit, Calculate, Grow, Save merged data
    "graph":    "mint",        # put something on the graph: Graph model, Results to graph
    "tool":     "blue",        # secondary tools and output files: Identify features, Report
    "analysis": "teal",        # uncertainty, Monte Carlo, tabular results: Calc. MC, Copy
    "helper":   "purple",      # helper sub-fits, special analyses: Fit background, Ask AI
    "revert":   "orange",      # undo / reset: Revert back, Reset to defaults
    "stop":     "red",         # stop or destroy: Cancel fit, Remove, ×
    "settings": "soft-green",  # settings to/from file: Store in File, Save/Load JSON
    "setup":    "soft-amber",  # load a saved setup: Load setup from file…
    "help":     "red",         # "? Help"
}


def resolve_roles(functions: dict[str, str] | None = None,
                  palette: dict[str, tuple[str, str, bool]] | None = None,
                  ) -> dict[str, tuple[str, str, bool]]:
    """Turn a function table and a palette into ``role -> (face, text, bold)``.

    Args:
        functions: function -> colour name; defaults to :data:`FUNCTIONS`.
            A value may also be a literal ``#rrggbb`` (white bold text).
        palette: colour name -> (face, text, bold); defaults to
            :data:`PALETTE`.

    Raises:
        KeyError: a function names a colour missing from the palette — a
            typo would otherwise silently fall back to neutral.
    """
    functions = FUNCTIONS if functions is None else functions
    palette = PALETTE if palette is None else palette
    out: dict[str, tuple[str, str, bool]] = {}
    for role, colour in functions.items():
        if colour.startswith("#"):
            out[role] = (colour, "#ffffff", True)
        else:
            out[role] = palette[colour]
    return out


#: role -> (face, text, bold), resolved from the two tables above.
DEFAULT_ROLES: dict[str, tuple[str, str, bool]] = resolve_roles()

#: size -> (vertical padding px, horizontal padding px, font size or None).
SIZES: dict[str, tuple[int, int, str | None]] = {
    # Same padding the flat look used, so no existing layout reflows.
    "normal": (3, 8, None),
    "small":  (1, 6, "10px"),
    "large":  (3, 8, "12px"),
}

# Disabled look: flat and grey on purpose — no bevel means "cannot press".
_DISABLED_BG = "#e3e6e8"
_DISABLED_TEXT = "#9aa3ab"
_DISABLED_BORDER = "#cfd4d8"
# Checked toggle (e.g. "Log X"): sunken, in the selection colour.
_CHECKED_BG = "#2f6fb5"
_CHECKED_TEXT = "#ffffff"


# ---------------------------------------------------------------------------
# Colour arithmetic (pure Python; runs once when the stylesheet is built)
# ---------------------------------------------------------------------------

def _rgb(colour: str) -> tuple[int, int, int]:
    h = colour.strip().lstrip("#")
    if len(h) == 3:
        h = "".join(c * 2 for c in h)
    if len(h) != 6:
        raise ValueError(f"expected #rgb or #rrggbb, got {colour!r}")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def mix(a: str, b: str, t: float) -> str:
    """Blend colour *a* toward *b* by fraction *t* (0 = a, 1 = b)."""
    t = min(max(t, 0.0), 1.0)
    ra, ga, ba = _rgb(a)
    rb, gb, bb = _rgb(b)
    return "#{:02x}{:02x}{:02x}".format(
        round(ra + (rb - ra) * t), round(ga + (gb - ga) * t), round(ba + (bb - ba) * t))


def lighten(colour: str, t: float) -> str:
    """Blend *colour* toward white by fraction *t*."""
    return mix(colour, "#ffffff", t)


def darken(colour: str, t: float) -> str:
    """Blend *colour* toward black by fraction *t*."""
    return mix(colour, "#000000", t)


def _gradient(top: str, mid: str, bottom: str) -> str:
    return (f"qlineargradient(x1:0, y1:0, x2:0, y2:1, "
            f"stop:0 {top}, stop:0.45 {mid}, stop:1 {bottom})")


# ---------------------------------------------------------------------------
# Stylesheet builders
# ---------------------------------------------------------------------------

def button_qss(bg: str, text: str, *, selector: str = "QPushButton",
               bold: bool = False, bevel: float = DEFAULT_BEVEL, radius: int = DEFAULT_RADIUS) -> str:
    """Colour and bevel rules for one button selector, all states.

    Sets no padding or font size — those come from :func:`size_qss`, so a
    colour role and a size combine freely.

    Args:
        bg: Face colour.
        text: Label colour (always set explicitly — never inherit it).
        selector: QSS selector, e.g. ``'QPushButton[btnRole="green"]'``.
        bold: Bold label. ``False`` leaves the widget's own font alone.
        bevel: Depth of the effect, 0 (flat) … ~0.45 (strong).
        radius: Corner radius in px.

    Returns:
        QSS text for the normal, hover, pressed, checked and disabled states.
    """
    s = bevel
    face = _gradient(lighten(bg, s), bg, darken(bg, s * 0.45))
    hover = _gradient(lighten(bg, s + 0.10), lighten(bg, 0.06), darken(bg, s * 0.30))
    sunk = _gradient(darken(bg, s * 0.55), darken(bg, s * 0.25), bg)
    # border-color order: top right bottom left
    raised_edge = (f"{darken(bg, 0.12 + s * 0.2)} {darken(bg, 0.28 + s * 0.4)} "
                   f"{darken(bg, 0.36 + s * 0.5)} {darken(bg, 0.14 + s * 0.2)}")
    sunk_edge = (f"{darken(bg, 0.36 + s * 0.5)} {darken(bg, 0.14 + s * 0.2)} "
                 f"{darken(bg, 0.12 + s * 0.2)} {darken(bg, 0.30 + s * 0.4)}")
    chk = _CHECKED_BG
    checked = _gradient(darken(chk, s * 0.55), darken(chk, s * 0.25), chk)
    checked_edge = (f"{darken(chk, 0.45)} {darken(chk, 0.2)} "
                    f"{darken(chk, 0.15)} {darken(chk, 0.4)}")
    # Only ever *add* bold: writing `font-weight: normal` would override a
    # bold font a panel set in code (QSS beats QWidget.setFont).
    weight = " font-weight: bold;" if bold else ""
    return (
        f"{selector} {{ color: {text};{weight} background: {face}; "
        f"border: 1px solid; border-color: {raised_edge}; border-radius: {radius}px; }}\n"
        f"{selector}:hover {{ color: {text}; background: {hover}; }}\n"
        f"{selector}:pressed {{ color: {text}; background: {sunk}; "
        f"border-color: {sunk_edge}; }}\n"
        f"{selector}:checked {{ color: {_CHECKED_TEXT}; background: {checked}; "
        f"border-color: {checked_edge}; }}\n"
        f"{selector}:disabled {{ color: {_DISABLED_TEXT}; background: {_DISABLED_BG}; "
        f"border-color: {_DISABLED_BORDER}; }}\n"
    )


def size_qss(selector: str, pad_v: int, pad_h: int, font_size: str | None = None) -> str:
    """Padding / font rules for one selector, including the pressed text shift.

    The pressed state moves one pixel of padding from bottom to top, so the
    label sinks with the face while the button's outer size stays the same.
    """
    font = f" font-size: {font_size};" if font_size else ""
    down_v = max(pad_v - 1, 0)
    return (
        f"{selector} {{ padding: {pad_v}px {pad_h}px;{font} }}\n"
        f"{selector}:pressed, {selector}:checked "
        f"{{ padding: {pad_v + 1}px {pad_h}px {down_v}px {pad_h}px; }}\n"
    )


def stylesheet(roles: dict[str, tuple[str, str, bool]] | None = None, *,
               bevel: float = DEFAULT_BEVEL, radius: int = DEFAULT_RADIUS,
               selectors: tuple[str, ...] = ("QPushButton", "QToolButton")) -> str:
    """The complete application-level button stylesheet.

    Unroled buttons get the ``neutral`` role; ``btnRole`` / ``btnSize``
    properties select the others. Append the result to the host's own QSS
    and install it once with ``QApplication.setStyleSheet``.

    Args:
        roles: Role table; defaults to :data:`DEFAULT_ROLES`. A host may
            extend or recolour it (pass its theme tokens), keeping the names.
        bevel: Depth of the effect, see :func:`button_qss`.
        radius: Corner radius in px.
        selectors: Widget classes to style.
    """
    roles = roles or DEFAULT_ROLES
    bg, fg, bold = roles.get("neutral", DEFAULT_ROLES["neutral"])
    parts: list[str] = []
    for sel in selectors:
        parts.append(button_qss(bg, fg, selector=sel, bold=bold, bevel=bevel, radius=radius))
        pad_v, pad_h, font = SIZES["normal"]
        parts.append(size_qss(sel, pad_v, pad_h, font))
    for name, (rbg, rfg, rbold) in roles.items():
        sel = f'QPushButton[{ROLE_PROPERTY}="{name}"]'
        parts.append(button_qss(rbg, rfg, selector=sel, bold=rbold, bevel=bevel, radius=radius))
    for name, (pad_v, pad_h, font) in SIZES.items():
        if name != "normal":
            sel = f'QPushButton[{SIZE_PROPERTY}="{name}"]'
            parts.append(size_qss(sel, pad_v, pad_h, font))
    return "".join(parts)


def inline_qss(role_or_colour: str, text: str | None = None, *, size: str = "normal",
               bold: bool | None = None, bevel: float = DEFAULT_BEVEL, radius: int = DEFAULT_RADIUS) -> str:
    """Self-contained stylesheet for one widget's ``setStyleSheet``.

    Use only where the app-level stylesheet may be absent (a panel embedded
    in someone else's application). Prefer :func:`set_role` otherwise.

    Args:
        role_or_colour: A role (:data:`FUNCTIONS`), a :data:`PALETTE`
            colour name, or a ``#hex``.
        text: Label colour when a hex colour is given (default white).
        size: Key of :data:`SIZES`.
        bold: Override the role's bold flag.
    """
    if role_or_colour in DEFAULT_ROLES:
        bg, fg, rbold = DEFAULT_ROLES[role_or_colour]
    elif role_or_colour in PALETTE:
        bg, fg, rbold = PALETTE[role_or_colour]
    else:
        bg, fg, rbold = role_or_colour, text or "#ffffff", True
    fg = text or fg
    pad_v, pad_h, font = SIZES[size]
    return (button_qss(bg, fg, bold=rbold if bold is None else bold,
                       bevel=bevel, radius=radius)
            + size_qss("QPushButton", pad_v, pad_h, font))


def set_role(button, role: str | None, size: str | None = None) -> None:
    """Give *button* a colour role (and optionally a size) at any time.

    Re-polishes the one widget so a change at runtime (Fit → Cancel) takes
    effect immediately. Duck-typed: works with any Qt binding.

    Args:
        button: A ``QPushButton``.
        role: Key of the installed role table, or ``None`` for neutral.
        size: ``'small'``, ``'large'``, or ``None`` for normal.
    """
    button.setProperty(ROLE_PROPERTY, role or "")
    button.setProperty(SIZE_PROPERTY, size or "")
    style = button.style()
    style.unpolish(button)
    style.polish(button)
    button.update()


def get_role(button) -> tuple[str | None, str | None]:
    """Return ``(role, size)`` previously set by :func:`set_role`.

    Use it to restore a button after a temporary change (Fit → Cancel → Fit).
    """
    role = button.property(ROLE_PROPERTY) or None
    size = button.property(SIZE_PROPERTY) or None
    return role, size
