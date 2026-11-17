"""Pattern for a page module (not loaded: the name starts with ``_``).

Copy to ``scripts/doc_shots/<page>.py`` (``gui-overview`` -> ``gui_overview.py``)
for ``docs/**/<page>.md``. One function per shot; the docstring at the top
names the dataset and its known answer, if there is one.
"""

from __future__ import annotations

from ._common import Context, draw_callouts, shot, widget_rect

SIZE = (1200, 900)  # fixed window size: same image on every run


def _open(ctx: Context):
    from PySide6.QtWidgets import QLabel  # replace: the package's real window

    window = QLabel("main window")
    ctx.show(window, SIZE)
    # window.load(ctx.dataset("some/file.h5"))   # always a scratch copy
    ctx.settle()
    return window


@shot("example-01-window")
def window(ctx: Context):
    w = _open(ctx)  # keep this reference until the grab: a shot that grabs only a
    image = w.grab()  # child widget can otherwise lose the window to garbage collection
    # Numbered callouts from the widgets' own geometry; the page's legend uses the same numbers.
    return draw_callouts(image, [(1, widget_rect(w, w), "c")])
