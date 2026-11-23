"""
pyirena.gui.theme — one deterministic look for every platform and OS setting.

Why this module exists
----------------------
pyIrena's panels carry several hundred inline Qt stylesheets written for a
light colour scheme (``background-color: #ecf0f1``, ``background-color:
lightgreen``, …).  Qt only overrides what a stylesheet actually names, so on a
machine whose *system* scheme is dark those widgets ended up with the dark
palette's near-white **text** painted on our hardcoded near-white
**background** — the "Fit B/P btwn cursors" / "Fit Flat btwn cursors" helper
buttons became invisible white boxes, and unstyled controls (the Modeling
panel's *Population type* and *Distribution* combo boxes) picked up whatever
the platform style happened to do with a dark palette.

Chasing every OS × theme combination is not winnable.  Instead pyIrena now
ships its own palette: :func:`apply_theme` installs the Fusion style plus an
explicit :class:`QPalette` and a baseline stylesheet on the ``QApplication``,
so the application looks the same whether the desktop is set to light or dark.
That also keeps the GUI consistent with the plots, which are drawn on a hard
white background throughout (see :mod:`pyirena.gui.sas_plot`).

Using it
--------
Entry points call :func:`apply_theme` right after creating the
``QApplication``::

    app = QApplication(sys.argv)
    apply_theme(app)

Panel code should take its colours from the tokens and helpers here rather
than hardcoding hex values, so a future dark variant is a matter of swapping
the token table::

    from pyirena.gui.theme import readonly_field_css, set_button_role

    set_button_role(fit_btn, 'primary')      # buttons: by function, never by colour
    field.setStyleSheet(readonly_field_css())

Every helper emits **both** ``background-color`` and ``color`` — that pairing
is the actual bug fix, and the rule to follow in new code.

One exception, learned the hard way: never put ``background-color`` on a
``QCheckBox`` or ``QRadioButton`` selector.  Qt copies a QSS background into
the :class:`QPalette` it hands the base style — ``Window`` *and* ``Base`` —
and Fusion builds the indicator's outline out of ``palette.window()`` and its
interior out of ``palette.base()``.  ``background-color: transparent`` on those
two selectors therefore made Fusion paint the box in a transparent pen, and
every "Fit?" check box in the Unified Fit panel became an empty gap.  Style
``::indicator`` instead; :func:`_base_stylesheet` does, and
``pyirena/tests/test_gui_theme_contract.py`` renders a check box to keep it
that way.
"""

from __future__ import annotations

import os
from pathlib import Path

from pyirena.gui import raised_buttons as _rb
from pyirena.gui._qt import QColor, QPalette

# ---------------------------------------------------------------------------
# Colour tokens
# ---------------------------------------------------------------------------
# Semantic names, not "blue"/"grey" — a dark variant can redefine the table
# without every call site changing.  Values are the light scheme pyIrena's
# panels were written against, so adopting the tokens is behaviour-preserving
# on a light desktop and *corrective* on a dark one.

# Base surfaces
WINDOW_BG      = '#f0f0f0'   # dialog / panel background
BASE_BG        = '#ffffff'   # text-entry and list background
ALT_BASE_BG    = '#f7f7f7'   # alternating rows
TEXT           = '#202020'   # primary text on WINDOW_BG / BASE_BG
TEXT_MUTED     = '#5a6570'   # secondary text, hints, units
TEXT_DISABLED  = '#9aa3ab'
BORDER         = '#b8c0c6'
BORDER_LIGHT   = '#d6dbdf'
# Outline for controls whose whole meaning is the outline — the box of a
# check box, the ring of a radio button.  BORDER is a 1.6:1 hairline there,
# which survives a 2x Retina pixel grid and vanishes on a 125%-scaled
# Windows or Linux display; this one clears 4.5:1 against every surface.
BORDER_STRONG  = '#6d767c'

# Raised "chip" surfaces — small helper buttons, read-only readouts
CHIP_BG        = '#e6eaed'
CHIP_BG_HOVER  = '#d5dbe0'
CHIP_TEXT      = '#20303c'

# Read-only field (a value the user cannot edit)
READONLY_BG    = '#ecf0f1'
READONLY_TEXT  = '#4d5a63'

# Selection
HIGHLIGHT_BG   = '#2f6fb5'
HIGHLIGHT_TEXT = '#ffffff'

# Status tints (background, text) — always specified as a readable pair
OK_BG,   OK_TEXT   = '#e7f6ec', '#1d6b38'
WARN_BG, WARN_TEXT = '#fdf3e0', '#8a5a08'
ERR_BG,  ERR_TEXT  = '#fdecea', '#a02c22'
INFO_BG, INFO_TEXT = '#eaf2fb', '#1c4f82'

# Saturated accents used for primary action buttons.  These are paired with
# white text, which is legible in either scheme, so they need no dark variant.
ACCENT_GREEN   = '#27ae60'
ACCENT_GREEN_D = '#1e8449'
ACCENT_MINT    = '#52c77a'
ACCENT_MINT_D  = '#3eb56a'
ACCENT_BLUE    = '#2980b9'
ACCENT_BLUE_D  = '#1f618d'
ACCENT_TEAL    = '#16a085'
ACCENT_TEAL_D  = '#138d75'
ACCENT_PURPLE  = '#8e44ad'
ACCENT_PURPLE_D = '#7d3c98'
ACCENT_ORANGE  = '#e67e22'
ACCENT_ORANGE_D = '#d35400'
ACCENT_RED     = '#c0392b'
ACCENT_RED_D   = '#a93226'
ACCENT_GREY    = '#95a5a6'

# "Soft" accents that are light enough to need dark text
SOFT_GREEN     = '#cdeccd'   # replaces the old bare `lightgreen`
SOFT_AMBER     = '#ffe082'
SOFT_TEXT      = '#20303c'

# ---------------------------------------------------------------------------
# Button colours — by function, from one table
# ---------------------------------------------------------------------------
# Panels never name a button colour.  They name what the button *does*:
#
#     set_button_role(self.fit_button, "primary", size="large")
#
# BUTTON_FUNCTIONS maps each function to a colour of BUTTON_PALETTE, and the
# stylesheet is generated from the result once, in apply_theme().  To recolour
# every "Revert back" / "Reset to defaults" button in pyIrena, change the
# 'revert' line below — nothing else.  The function names are shared with
# Jan's other Qt packages (raised_buttons.FUNCTIONS); keep them.

#: colour name -> (face, text, bold), filled from pyIrena's tokens.
BUTTON_PALETTE: dict[str, tuple[str, str, bool]] = {
    'neutral':    (CHIP_BG, CHIP_TEXT, False),
    'green':      (ACCENT_GREEN, '#ffffff', True),
    'mint':       (ACCENT_MINT, '#ffffff', True),
    'blue':       (ACCENT_BLUE, '#ffffff', True),
    'teal':       (ACCENT_TEAL, '#ffffff', True),
    'purple':     (ACCENT_PURPLE, '#ffffff', True),
    'orange':     (ACCENT_ORANGE, '#ffffff', True),
    'red':        (ACCENT_RED, '#ffffff', True),
    'grey':       (_rb.PALETTE['grey'][0], '#ffffff', True),
    'soft-green': (SOFT_GREEN, SOFT_TEXT, False),
    'soft-amber': (SOFT_AMBER, SOFT_TEXT, False),
    # Darker shades for the launcher's "(script)" twins, plus Carbon's slate.
    'green-dark':  (ACCENT_GREEN_D, '#ffffff', True),
    'blue-dark':   (ACCENT_BLUE_D, '#ffffff', True),
    'orange-dark': (ACCENT_ORANGE_D, '#ffffff', True),
    'purple-dark': ('#6c3483', '#ffffff', True),
    'slate':       ('#34495e', '#ffffff', True),
    'slate-dark':  ('#2c3e50', '#ffffff', True),
}

#: button function -> colour name.  Starts from the shared vocabulary; add a
#: pyIrena-only function here if one is ever needed.
BUTTON_FUNCTIONS: dict[str, str] = {
    **_rb.FUNCTIONS,
    # 'help': 'blue',   # e.g. — one line recolours every "? Help" button

    # Main-window launchers (pyIrena only).  Each tool has an identity colour
    # so users find it at a glance; its "(script)" button is a darker shade.
    'launch-unified':            'green',
    'launch-unified-script':     'green-dark',
    'launch-sizes':              'blue',
    'launch-sizes-script':       'blue-dark',
    'launch-modeling':           'orange',
    'launch-modeling-script':    'orange-dark',
    'launch-simple-fits':        'green',
    'launch-simple-fits-script': 'green-dark',
    'launch-waxs':               'blue',
    'launch-waxs-script':        'blue-dark',
    'launch-carbon':             'slate',
    'launch-carbon-script':      'slate-dark',
    'launch-saxs-morph':         'purple',
    'launch-saxs-morph-script':  'purple-dark',
    'launch-utility':            'teal',   # Data Merge, Manipulation, Contrast, Fractals
    'launch-explorer':           'red',    # Data Explorer: a different kind of tool
}

#: role -> (face, text, bold) — what the stylesheet is generated from.
BUTTON_ROLES: dict[str, tuple[str, str, bool]] = _rb.resolve_roles(
    BUTTON_FUNCTIONS, BUTTON_PALETTE)

__all__ = [
    'apply_theme', 'is_theme_applied', 'set_button_role', 'button_role',
    'BUTTON_ROLES', 'BUTTON_FUNCTIONS', 'BUTTON_PALETTE',
    'readonly_field_css', 'readout_label_css', 'status_css', 'muted_css',
    'READONLY_FIELD_CSS',
    'WINDOW_BG', 'BASE_BG', 'TEXT', 'TEXT_MUTED', 'TEXT_DISABLED',
    'BORDER', 'BORDER_LIGHT', 'BORDER_STRONG',
    'CHIP_BG', 'CHIP_BG_HOVER', 'CHIP_TEXT',
    'READONLY_BG', 'READONLY_TEXT',
    'OK_BG', 'OK_TEXT', 'WARN_BG', 'WARN_TEXT', 'ERR_BG', 'ERR_TEXT',
    'INFO_BG', 'INFO_TEXT',
    'ACCENT_GREEN', 'ACCENT_GREEN_D', 'ACCENT_MINT', 'ACCENT_MINT_D',
    'ACCENT_BLUE', 'ACCENT_BLUE_D', 'ACCENT_TEAL', 'ACCENT_TEAL_D',
    'ACCENT_PURPLE', 'ACCENT_PURPLE_D', 'ACCENT_ORANGE', 'ACCENT_ORANGE_D',
    'ACCENT_RED', 'ACCENT_RED_D', 'ACCENT_GREY',
    'SOFT_GREEN', 'SOFT_AMBER', 'SOFT_TEXT',
]


# ---------------------------------------------------------------------------
# CSS helpers — every one of these emits background *and* text colour
# ---------------------------------------------------------------------------

def set_button_role(button, role: str | None, size: str | None = None) -> None:
    """Colour *button* by role — the preferred way to style a push button.

    Thin wrapper over :func:`pyirena.gui.raised_buttons.set_role`; the role
    rules themselves are part of the application stylesheet installed by
    :func:`apply_theme`, so this costs one property set and one re-polish.

    Args:
        button: A ``QPushButton``.
        role: A button function (key of :data:`BUTTON_ROLES`: ``'primary'``,
            ``'graph'``, ``'stop'``, …), or ``None`` for neutral.
        size: ``'small'``, ``'large'`` or ``None``.
    """
    if role is not None and role not in BUTTON_ROLES:
        raise KeyError(f'unknown button role {role!r}; use one of {sorted(BUTTON_ROLES)}')
    _rb.set_role(button, role, size)


def button_role(button) -> tuple[str | None, str | None]:
    """``(role, size)`` of *button*, for restoring it after a temporary change."""
    return _rb.get_role(button)


def readonly_field_css() -> str:
    """Stylesheet for a read-only ``QLineEdit`` used as a value readout."""
    return (
        f'QLineEdit {{ background-color: {READONLY_BG}; color: {READONLY_TEXT}; '
        f'border: 1px solid {BORDER_LIGHT}; border-radius: 2px; }}'
    )


#: Ready-made read-only field stylesheet.
READONLY_FIELD_CSS = readonly_field_css()


def readout_label_css() -> str:
    """Stylesheet for a boxed read-only ``QLabel`` (derived results, Vf, …)."""
    return (
        f'color: {READONLY_TEXT}; background-color: {READONLY_BG}; '
        f'padding: 1px 4px; border-radius: 2px;'
    )


def status_css(kind: str = 'info') -> str:
    """Stylesheet for a status banner.

    Args:
        kind: One of ``'ok'``, ``'warning'``, ``'error'``, ``'info'``.
    """
    table = {
        'ok': (OK_BG, OK_TEXT),
        'success': (OK_BG, OK_TEXT),
        'warning': (WARN_BG, WARN_TEXT),
        'warn': (WARN_BG, WARN_TEXT),
        'error': (ERR_BG, ERR_TEXT),
        'info': (INFO_BG, INFO_TEXT),
    }
    bg, fg = table.get(kind, (INFO_BG, INFO_TEXT))
    return (f'background-color: {bg}; color: {fg}; '
            f'border: 1px solid {BORDER_LIGHT}; border-radius: 3px; '
            f'padding: 3px 6px;')


def muted_css(font_size: str = '10px', *, italic: bool = False) -> str:
    """Stylesheet for a secondary hint label."""
    style = ' font-style: italic;' if italic else ''
    return f'font-size: {font_size}; color: {TEXT_MUTED};{style}'


# ---------------------------------------------------------------------------
# Application-level theme
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Indicator glyphs
# ---------------------------------------------------------------------------
# A check box drawn by a stylesheet has to supply its own tick: once a QSS rule
# makes ``::indicator`` drawable, Qt renders the rule and never calls the base
# style that would otherwise paint the mark.  The two SVGs are shipped in
# ``pyirena/gui/assets`` (see ``package-data`` in pyproject.toml); Qt rasterises
# them at device resolution, so they stay crisp on a Retina or 150%-scaled
# display.

_ASSET_DIR = Path(__file__).resolve().parent / 'assets'


def _glyph(name: str) -> str:
    """``image: url(...)`` for a shipped indicator glyph.

    Returns an empty declaration when the file is missing — an install that
    dropped the assets then shows a plain accent-filled box, which still reads
    as "checked", rather than raising.

    Args:
        name: File name inside ``pyirena/gui/assets``.
    """
    path = _ASSET_DIR / name
    if not path.is_file():
        return ''
    # QSS wants a filesystem path with forward slashes, on Windows too.
    return f'image: url("{path.as_posix()}");'


_APPLIED = False


def is_theme_applied() -> bool:
    """True once :func:`apply_theme` has run in this process."""
    return _APPLIED


def _palette() -> QPalette:
    """Build the explicit light palette pyIrena renders against."""
    p = QPalette()
    role = QPalette.ColorRole
    group = QPalette.ColorGroup

    p.setColor(role.Window,          QColor(WINDOW_BG))
    p.setColor(role.WindowText,      QColor(TEXT))
    p.setColor(role.Base,            QColor(BASE_BG))
    p.setColor(role.AlternateBase,   QColor(ALT_BASE_BG))
    p.setColor(role.Text,            QColor(TEXT))
    p.setColor(role.PlaceholderText, QColor(TEXT_DISABLED))
    p.setColor(role.Button,          QColor(WINDOW_BG))
    p.setColor(role.ButtonText,      QColor(TEXT))
    p.setColor(role.BrightText,      QColor('#ffffff'))
    p.setColor(role.ToolTipBase,     QColor('#ffffdc'))
    p.setColor(role.ToolTipText,     QColor(TEXT))
    p.setColor(role.Highlight,       QColor(HIGHLIGHT_BG))
    p.setColor(role.HighlightedText, QColor(HIGHLIGHT_TEXT))
    p.setColor(role.Link,            QColor(ACCENT_BLUE))
    p.setColor(role.LinkVisited,     QColor(ACCENT_PURPLE))
    p.setColor(role.Light,           QColor('#ffffff'))
    p.setColor(role.Midlight,        QColor(BORDER_LIGHT))
    p.setColor(role.Mid,             QColor(BORDER))
    p.setColor(role.Dark,            QColor('#8d979e'))
    p.setColor(role.Shadow,          QColor(BORDER_STRONG))

    for r in (role.WindowText, role.Text, role.ButtonText,
              role.HighlightedText, role.PlaceholderText):
        p.setColor(group.Disabled, r, QColor(TEXT_DISABLED))
    p.setColor(group.Disabled, role.Base, QColor(WINDOW_BG))
    p.setColor(group.Disabled, role.Button, QColor(WINDOW_BG))
    # Inactive windows keep readable text (macOS otherwise greys it heavily).
    for r in (role.WindowText, role.Text, role.ButtonText):
        p.setColor(group.Inactive, r, QColor(TEXT))
    p.setColor(group.Inactive, role.Highlight, QColor('#c9d6e5'))
    p.setColor(group.Inactive, role.HighlightedText, QColor(TEXT))
    return p


def _base_stylesheet() -> str:
    """Baseline QSS applied to the whole application.

    Covers the controls that carried no inline stylesheet at all and were
    therefore at the mercy of the platform style under a dark palette — combo
    boxes (``Population type``, ``Distribution``), spin boxes, line edits,
    check boxes, tabs and headers.  Each rule names foreground *and*
    background so the result cannot depend on the desktop's scheme.
    """
    tick = _glyph('checkbox-tick.svg')
    dot = _glyph('radio-dot.svg')
    return f"""
    QWidget {{ color: {TEXT}; }}
    QMainWindow, QDialog, QScrollArea, QTabWidget::pane, QStackedWidget {{
        background-color: {WINDOW_BG}; color: {TEXT};
    }}
    QToolTip {{
        background-color: #ffffdc; color: {TEXT};
        border: 1px solid {BORDER}; padding: 3px;
    }}

    /* --- text entry ---------------------------------------------------- */
    QLineEdit, QPlainTextEdit, QTextEdit, QTextBrowser {{
        background-color: {BASE_BG}; color: {TEXT};
        border: 1px solid {BORDER}; border-radius: 3px;
        selection-background-color: {HIGHLIGHT_BG};
        selection-color: {HIGHLIGHT_TEXT};
    }}
    QLineEdit:disabled, QPlainTextEdit:disabled, QTextEdit:disabled {{
        background-color: {WINDOW_BG}; color: {TEXT_DISABLED};
        border-color: {BORDER_LIGHT};
    }}
    QLineEdit:read-only {{
        background-color: {READONLY_BG}; color: {READONLY_TEXT};
    }}

    /* --- combo boxes and spin boxes ------------------------------------ */
    /* These carried no inline style, so a dark system palette decided how
       they looked.  Pin both colours, and pin the popup list too — the popup
       is a separate top-level widget and does not inherit the editor's. */
    QComboBox {{
        background-color: {BASE_BG}; color: {TEXT};
        border: 1px solid {BORDER}; border-radius: 3px;
        padding: 2px 4px; min-height: 18px;
    }}
    QComboBox:hover {{ border-color: {HIGHLIGHT_BG}; }}
    QComboBox:disabled {{
        background-color: {WINDOW_BG}; color: {TEXT_DISABLED};
        border-color: {BORDER_LIGHT};
    }}
    /* NOTE: ::drop-down is deliberately left unstyled.  Styling it makes
       Fusion stop drawing its built-in arrow (a QSS-styled sub-control must
       supply its own image), which leaves the combo looking like a plain
       text field.  The colours above are enough to guarantee contrast. */
    QComboBox QAbstractItemView {{
        background-color: {BASE_BG}; color: {TEXT};
        border: 1px solid {BORDER};
        selection-background-color: {HIGHLIGHT_BG};
        selection-color: {HIGHLIGHT_TEXT};
        outline: none;
    }}
    QSpinBox, QDoubleSpinBox {{
        background-color: {BASE_BG}; color: {TEXT};
        border: 1px solid {BORDER}; border-radius: 3px; padding: 1px 2px;
    }}
    QSpinBox:disabled, QDoubleSpinBox:disabled {{
        background-color: {WINDOW_BG}; color: {TEXT_DISABLED};
    }}

    /* --- buttons: see pyirena.gui.raised_buttons (appended below) ---- */

    /* --- selection controls -------------------------------------------- */
    /* Do NOT add `background-color: transparent` to QCheckBox or
       QRadioButton.  Qt copies a QSS background into the *palette* it hands
       the base style (QRenderRule::configurePalette writes it to both Window
       and Base), and Fusion derives the indicator's outline from
       palette.window() and its interior from palette.base() — so a
       transparent background made Fusion paint the box in a fully
       transparent pen.  Every "Fit?" check box in the Unified Fit panel was
       an empty gap with, at most, a floating tick.  Neither widget paints a
       background of its own, so leaving it unset is also what we want. */
    /* `padding` restores the height Fusion's own metrics gave these: a
       styled ::indicator makes Qt size the widget from the rule, which
       would otherwise shrink every check box row by 5 px. */
    QCheckBox, QRadioButton {{
        color: {TEXT}; spacing: 6px; padding: 2px 0;
    }}
    QGroupBox, QLabel {{ background-color: transparent; color: {TEXT}; }}
    QCheckBox:disabled, QRadioButton:disabled, QLabel:disabled {{
        color: {TEXT_DISABLED};
    }}

    /* The indicator is pinned here rather than left to Fusion, which derives
       it from palette.window().darker(140) — #bdbdbd against WINDOW_BG, a
       1.6:1 hairline that a 125% or 150% display scale rounds away to
       nothing.  A drawable ::indicator rule means Qt stops calling the base
       style, so the tick and the dot have to come from our own glyphs. */
    QCheckBox::indicator, QRadioButton::indicator, QGroupBox::indicator {{
        width: 12px; height: 12px;
        background-color: {BASE_BG};
        border: 1px solid {BORDER_STRONG};
    }}
    QCheckBox::indicator, QGroupBox::indicator {{ border-radius: 3px; }}
    QRadioButton::indicator {{ border-radius: 7px; }}
    QCheckBox::indicator:hover, QRadioButton::indicator:hover,
    QGroupBox::indicator:hover {{ border-color: {HIGHLIGHT_BG}; }}
    QCheckBox::indicator:checked, QGroupBox::indicator:checked {{
        background-color: {ACCENT_BLUE}; border-color: {ACCENT_BLUE_D};
        {tick}
    }}
    QRadioButton::indicator:checked {{
        background-color: {ACCENT_BLUE}; border-color: {ACCENT_BLUE_D};
        {dot}
    }}
    QCheckBox::indicator:disabled, QRadioButton::indicator:disabled,
    QGroupBox::indicator:disabled {{
        background-color: {WINDOW_BG}; border-color: {BORDER};
    }}
    QCheckBox::indicator:checked:disabled,
    QRadioButton::indicator:checked:disabled,
    QGroupBox::indicator:checked:disabled {{
        background-color: {TEXT_DISABLED}; border-color: {TEXT_DISABLED};
    }}
    QGroupBox {{
        border: 1px solid {BORDER_LIGHT}; border-radius: 4px;
        margin-top: 8px; padding-top: 4px;
    }}
    QGroupBox::title {{
        subcontrol-origin: margin; subcontrol-position: top left;
        left: 8px; padding: 0 3px; color: {TEXT};
    }}

    /* --- containers ----------------------------------------------------- */
    QTabWidget::pane {{ border: 1px solid {BORDER_LIGHT}; }}
    QTabBar::tab {{
        background-color: {WINDOW_BG}; color: {TEXT};
        border: 1px solid {BORDER_LIGHT}; border-bottom: none;
        padding: 4px 10px; margin-right: 1px;
    }}
    QTabBar::tab:selected {{ background-color: {BASE_BG}; color: {TEXT}; }}
    QTabBar::tab:!selected:hover {{ background-color: {CHIP_BG}; }}
    QTabBar::tab:disabled {{ color: {TEXT_DISABLED}; }}

    QTableWidget, QTableView, QTreeWidget, QTreeView,
    QListWidget, QListView {{
        background-color: {BASE_BG}; color: {TEXT};
        alternate-background-color: {ALT_BASE_BG};
        border: 1px solid {BORDER_LIGHT};
        selection-background-color: {HIGHLIGHT_BG};
        selection-color: {HIGHLIGHT_TEXT};
    }}
    QHeaderView::section {{
        background-color: {CHIP_BG}; color: {CHIP_TEXT};
        border: none; border-right: 1px solid {BORDER_LIGHT};
        border-bottom: 1px solid {BORDER_LIGHT}; padding: 3px 4px;
    }}
    QMenuBar {{ background-color: {WINDOW_BG}; color: {TEXT}; }}
    QMenuBar::item:selected {{
        background-color: {HIGHLIGHT_BG}; color: {HIGHLIGHT_TEXT};
    }}
    QMenu {{
        background-color: {BASE_BG}; color: {TEXT};
        border: 1px solid {BORDER};
    }}
    QMenu::item:selected {{
        background-color: {HIGHLIGHT_BG}; color: {HIGHLIGHT_TEXT};
    }}
    QMenu::item:disabled {{ color: {TEXT_DISABLED}; }}
    QStatusBar {{ background-color: {WINDOW_BG}; color: {TEXT}; }}
    QProgressBar {{
        background-color: {BASE_BG}; color: {TEXT};
        border: 1px solid {BORDER}; border-radius: 3px; text-align: center;
    }}
    QProgressBar::chunk {{ background-color: {ACCENT_BLUE}; }}
    QSplitter::handle {{ background-color: {BORDER_LIGHT}; }}
    QScrollBar:vertical, QScrollBar:horizontal {{
        background-color: {WINDOW_BG}; border: none;
    }}
    QScrollBar::handle:vertical, QScrollBar::handle:horizontal {{
        background-color: {BORDER}; border-radius: 4px;
    }}
    QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; width: 0; }}
    """ + _button_stylesheet()


def _button_stylesheet() -> str:
    """Raised push/tool buttons, one rule set per :data:`BUTTON_ROLES` entry.

    Replaces the old flat ``QPushButton { background; border }`` rule, which
    made every button a coloured rectangle.  The checked state (graph-window
    X/Y log toggles) is sunken in the selection colour, so on/off differ.
    """
    return _rb.stylesheet(BUTTON_ROLES)


def apply_theme(app=None) -> None:
    """Install pyIrena's own light theme on *app*.

    Sets the Fusion style, an explicit :class:`QPalette` and a baseline
    stylesheet, so the GUI renders identically whether the desktop is set to
    light or dark and on any platform.  Idempotent — calling it again (for
    instance when a sub-tool is launched from the main window and also has a
    standalone ``main()``) is a no-op.

    Set ``PYIRENA_NATIVE_THEME=1`` to skip this and fall back to the platform
    style and system palette, for debugging or for users who prefer it.

    Args:
        app: The ``QApplication``; the running instance is used when omitted.
    """
    global _APPLIED
    if _APPLIED:
        return
    if os.environ.get('PYIRENA_NATIVE_THEME', '').strip().lower() not in (
            '', '0', 'false', 'no'):
        _APPLIED = True
        return

    from pyirena.gui._qt import QApplication
    app = app or QApplication.instance()
    if app is None:
        return

    app.setStyle('Fusion')
    app.setPalette(_palette())
    # Fusion caches per-widget palettes for widgets created before the call;
    # setting the style *first* and the palette after covers the normal
    # "apply right after QApplication()" order used by every entry point.
    app.setStyleSheet(_base_stylesheet())
    _APPLIED = True
