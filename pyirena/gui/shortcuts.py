"""
Standard keyboard shortcuts for pyIrena windows.

Tool windows are independent top-level widgets with no common base class, so
until now none of them answered ⌘W — the one shortcut every Mac user tries
first (GitHub issue #16). This module is the single place that decides what
the standard keys are, and :func:`install_standard_shortcuts` is the one call
a window makes to get them.

Platform correctness comes from ``QKeySequence.StandardKey`` rather than from
hard-coded strings: ``StandardKey.Close`` is ⌘W on macOS and Ctrl+W (plus
Ctrl+F4) on Windows, and Qt picks per platform. Writing ``"Ctrl+W"`` would
give Mac users Control+W, which is not the same key and is not what they
press. ``keyBindings()`` can return more than one sequence for a key, so each
one gets its own shortcut.

What is deliberately *not* here
-------------------------------
Nothing invents behaviour. ``on_save`` and ``on_help`` wire ⌘S and F1 to a
button the window already shows; a window without that button passes nothing
and simply has no such shortcut. ⌘C already works in tables (see
``table_utils``) and is left alone. ⌘Q is left to Qt and the OS — a quit key
that fires mid-fit is worse than no quit key.

Scope: every shortcut is a ``WindowShortcut``, so it fires only when its own
window is the active one, and two open tool windows do not fight over ⌘W.

Coverage is enforced by ``pyirena/tests/test_window_shortcuts.py``: a new
top-level window class must either call this or carry a written reason.
"""
from __future__ import annotations

import logging
from typing import Callable, Optional, Union

from pyirena.gui._qt import QKeySequence, QShortcut, Qt

log = logging.getLogger(__name__)

#: A shortcut target: a button to click, or something to call.
Target = Union[Callable[[], None], object]


def _activate(target: Target) -> None:
    """Press *target* if it is a button, otherwise call it."""
    click = getattr(target, "click", None)
    if callable(click):
        # A QPushButton: go through click() rather than the slot so the window
        # sees exactly what it would see from a mouse press — including any
        # enabled/disabled state, which is the point. A disabled button does
        # nothing, which is the correct answer for "Save" with nothing to save.
        click()
    elif callable(target):
        target()
    else:
        log.debug("shortcut target %r is neither a button nor callable", target)


def _bind(widget, standard_key, target: Target) -> list:
    """Create one shortcut per platform binding of *standard_key*."""
    made = []
    try:
        sequences = QKeySequence.keyBindings(standard_key)
    except Exception:
        log.debug("could not read key bindings for %r", standard_key, exc_info=True)
        return made

    for seq in sequences:
        try:
            sc = QShortcut(seq, widget)
            sc.setContext(Qt.ShortcutContext.WindowShortcut)
            sc.activated.connect(lambda t=target: _activate(t))
            made.append(sc)
        except Exception:
            log.debug("could not install %s on %r", seq.toString(), widget, exc_info=True)
    return made


def install_standard_shortcuts(
    widget,
    *,
    close: bool = True,
    on_save: Optional[Target] = None,
    on_help: Optional[Target] = None,
) -> list:
    """Give *widget* the standard window shortcuts. Call once, after the UI is built.

    ::

        install_standard_shortcuts(self, on_save=self.btn_save,
                                   on_help=self.help_btn)

    Parameters
    ----------
    widget :
        The top-level window. Shortcuts are children of it and die with it.
    close :
        Install Close (⌘W / Ctrl+W). True for a window; pass False for one
        that must not be closed from the keyboard.
    on_save, on_help :
        A ``QPushButton`` the window already shows, or a zero-argument
        callable, for Save (⌘S / Ctrl+S) and Help (F1). Omit either when the
        window has no such action — no shortcut is installed for it.

    Returns
    -------
    list
        The ``QShortcut`` objects created, for tests and for callers that want
        to disable one later. Never raises: a window that cannot take a
        shortcut is a degraded window, not a broken one.
    """
    shortcuts: list = []

    if close:
        shortcuts += _bind(widget, QKeySequence.StandardKey.Close, widget.close)
    if on_save is not None:
        shortcuts += _bind(widget, QKeySequence.StandardKey.Save, on_save)
    if on_help is not None:
        shortcuts += _bind(widget, QKeySequence.StandardKey.HelpContents, on_help)

    return shortcuts
