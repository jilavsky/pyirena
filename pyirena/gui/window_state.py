"""
Remember where each pyIrena window was, and put it back safely.

Igor and Origin users expect a tool to reopen where they left it, at the size
they left it, with the control panel the width they dragged it to.  The part
that makes this harder than it looks is **display changes**: laptops get
undocked, a second monitor disappears, resolution changes, a window that was on
the right-hand screen is now at coordinates no screen covers.  Restoring such a
position blindly puts the window somewhere the user cannot see or drag it back
from — a bug far worse than not remembering at all.

So the decision of *where a saved rectangle may actually go* is isolated in one
pure function, :func:`resolve_geometry`, which takes plain ``(x, y, w, h)``
tuples and the list of current screen rectangles and returns either an adjusted
rectangle or ``None`` (meaning "give up, use the window's default").  It has no
Qt dependency, so its behaviour can be tested for every awkward case —
off-screen, half-off-screen, larger-than-screen, screen removed — without a
display.  Expect to tune the constants below; they are deliberately all in one
place.

Escape hatches, in order of convenience
---------------------------------------
1. Hold **Shift** while clicking a tool's button — *that* tool forgets its
   saved geometry and opens at its default size, centred.  This is the Irena
   gesture and it works every time, not only on the first launch of a session:
   :func:`install_window_state` records the window's coded default before it
   applies anything saved, so :func:`reset_window` can put an already-built
   window back to it.
2. Set ``PYIRENA_RESET_WINDOWS=1`` in the environment — *every* window opens at
   its default and the saved geometry is cleared.  Use this when a window has
   ended up somewhere you cannot reach it.
3. Call :func:`clear_window_state` from code or a console.

Turning the feature off
-----------------------
Set :data:`GEOMETRY_PERSISTENCE_ENABLED` to ``False``.  Every function then
becomes a no-op and windows use their coded defaults, with no other change
needed anywhere in the GUI.
"""

from __future__ import annotations

import logging
import os
from typing import Dict, List, Optional, Sequence, Tuple

log = logging.getLogger(__name__)

__all__ = [
    "GEOMETRY_PERSISTENCE_ENABLED",
    "clear_window_state",
    "install_window_state",
    "reset_requested",
    "reset_window",
    "reset_window_if_shift",
    "resolve_geometry",
    "restore_window_state",
    "save_window_state",
    "screen_rects",
]

# ──────────────────────────────────────────────────────────────────────────
#  Feature switch — flip to False to ship without geometry persistence
# ──────────────────────────────────────────────────────────────────────────

#: When False every public function is a no-op and windows open at their
#: coded defaults.  One switch, no other edits needed.
GEOMETRY_PERSISTENCE_ENABLED = True

# ──────────────────────────────────────────────────────────────────────────
#  Tunables — the placement policy lives here
# ──────────────────────────────────────────────────────────────────────────

#: A restored window must overlap some screen by at least this many pixels in
#: each direction, otherwise it counts as "lost" and is relocated.
MIN_VISIBLE_PX = 120

#: Height of the strip at the top of the window that carries the title bar.
#: It has to be on-screen or the user cannot drag the window anywhere.
TITLE_BAR_GRAB_PX = 32

#: Gap left between a relocated window and the screen edge.
EDGE_MARGIN_PX = 24

#: Smallest window we will ever restore; anything smaller is treated as junk.
MIN_WINDOW_PX = (240, 180)

#: File name (beside the main state file) holding one entry per window key.
#:
#: Deliberately *not* a section of the shared state file: every panel owns its
#: own ``StateManager`` instance and ``save()`` rewrites the whole file from
#: that instance's in-memory copy, so geometry written while closing one panel
#: is erased when the next panel saves its (older) copy.  Window geometry is
#: also written far more often than tool state — on every window close — and is
#: pure UI preference, so a small separate file keeps it out of harm's way and
#: makes "delete it to reset" a one-liner.
GEOMETRY_FILE_NAME = "window_geometry.json"

#: Environment variable that forces every window back to its default.
RESET_ENV_VAR = "PYIRENA_RESET_WINDOWS"

Rect = Tuple[int, int, int, int]     # (x, y, w, h)


# ──────────────────────────────────────────────────────────────────────────
#  The placement policy (pure — no Qt, fully testable)
# ──────────────────────────────────────────────────────────────────────────

def _intersection(a: Rect, b: Rect) -> Tuple[int, int]:
    """Overlap of two rectangles as (width, height); zero or negative if none."""
    ax, ay, aw, ah = a
    bx, by, bw, bh = b
    return (min(ax + aw, bx + bw) - max(ax, bx),
            min(ay + ah, by + bh) - max(ay, by))


def _is_comfortably_visible(rect: Rect, screen: Rect) -> bool:
    """True when *rect* overlaps *screen* enough to be seen and grabbed."""
    overlap_w, overlap_h = _intersection(rect, screen)
    if overlap_w < MIN_VISIBLE_PX or overlap_h < MIN_VISIBLE_PX:
        return False
    # The title bar must be on-screen, or the window cannot be moved by hand.
    title_bar = (rect[0], rect[1], rect[2], TITLE_BAR_GRAB_PX)
    bar_w, bar_h = _intersection(title_bar, screen)
    return bar_w >= MIN_VISIBLE_PX and bar_h > 0


def _centre(rect: Rect) -> Tuple[float, float]:
    return rect[0] + rect[2] / 2.0, rect[1] + rect[3] / 2.0


def _nearest_screen(rect: Rect, screens: Sequence[Rect]) -> Rect:
    """The screen whose centre is closest to *rect*'s centre."""
    rx, ry = _centre(rect)
    return min(screens,
               key=lambda s: (_centre(s)[0] - rx) ** 2 + (_centre(s)[1] - ry) ** 2)


def resolve_geometry(saved: Optional[Rect],
                     screens: Sequence[Rect]) -> Optional[Rect]:
    """Decide where a saved window rectangle may actually be placed.

    This is the whole placement policy, kept pure so every awkward case can be
    tested without a display.

    The rules, in order:

    1. No saved rectangle, no screens, or a nonsense size → ``None``
       ("use the window's own default").
    2. The rectangle is comfortably visible on some screen → use it unchanged.
       This is the ordinary case and must stay cheap.
    3. Otherwise the display setup has changed: shrink the window to fit the
       nearest screen if it is now too big, then slide it inside that screen's
       bounds.  The window keeps its size where possible and always lands
       somewhere the user can see and drag.
    4. If even that cannot produce a sane rectangle → ``None``.

    Args:
        saved: The remembered ``(x, y, w, h)``, or None.
        screens: Available screen rectangles, in the same coordinate space.

    Returns:
        A rectangle to apply, or None to fall back to the window's default.
    """
    if not saved or not screens:
        return None

    try:
        x, y, w, h = (int(v) for v in saved)
    except (TypeError, ValueError):
        return None

    if w < MIN_WINDOW_PX[0] or h < MIN_WINDOW_PX[1]:
        # A collapsed or corrupt entry — better to use the coded default.
        return None

    rect: Rect = (x, y, w, h)

    # 2 — still fits where it was.
    for screen in screens:
        if _is_comfortably_visible(rect, screen):
            return rect

    # 3 — relocate onto the nearest screen, shrinking only if we must.
    screen = _nearest_screen(rect, screens)
    sx, sy, sw, sh = screen

    new_w = min(w, max(MIN_WINDOW_PX[0], sw - 2 * EDGE_MARGIN_PX))
    new_h = min(h, max(MIN_WINDOW_PX[1], sh - 2 * EDGE_MARGIN_PX))

    # Keep the window's relative position on the screen where sensible, but
    # never outside it.
    new_x = min(max(x, sx + EDGE_MARGIN_PX), sx + sw - new_w - EDGE_MARGIN_PX)
    new_y = min(max(y, sy + EDGE_MARGIN_PX), sy + sh - new_h - EDGE_MARGIN_PX)
    # A screen narrower than the margins would invert the clamp above.
    new_x = max(new_x, sx)
    new_y = max(new_y, sy)

    relocated: Rect = (int(new_x), int(new_y), int(new_w), int(new_h))
    if not _is_comfortably_visible(relocated, screen):
        return None
    log.info("Window geometry %s was off-screen; relocated to %s", rect, relocated)
    return relocated


# ──────────────────────────────────────────────────────────────────────────
#  Qt glue
# ──────────────────────────────────────────────────────────────────────────

def screen_rects() -> List[Rect]:
    """Available geometry of every connected screen, as plain tuples.

    "Available" excludes the menu bar and taskbar, so a window placed inside
    one of these rectangles is genuinely usable.
    """
    try:
        from pyirena.gui._qt import QApplication  # noqa: PLC0415

        app = QApplication.instance()
        if app is None:
            return []
        rects = []
        for screen in QApplication.screens():
            geo = screen.availableGeometry()
            rects.append((geo.x(), geo.y(), geo.width(), geo.height()))
        return rects
    except Exception:
        log.debug("could not read screen geometry", exc_info=True)
        return []


def is_top_level(widget) -> bool:
    """True when *widget* is its own window rather than a pane inside one.

    Geometry only means something for a window.  Several pyIrena classes are
    called ``…GraphWindow`` but are actually added to a splitter as the right
    pane, and setting their geometry fights the layout: the pane jumps to the
    remembered rectangle and the control panel next to it changes width.  Every
    apply/save path checks this, so registering the wrong widget is harmless.
    """
    try:
        return bool(widget.isWindow())
    except Exception:
        # Not a QWidget (a test double, say) — assume the caller meant it.
        return True


def _screen_fingerprint(screens: Sequence[Rect]) -> str:
    """A short signature of the display setup, for logging a change."""
    return "|".join(f"{x},{y},{w},{h}" for x, y, w, h in sorted(screens))


def shift_held() -> bool:
    """True when the Shift key is down right now.

    ``queryKeyboardModifiers()`` asks the window system for the *current* state,
    where ``keyboardModifiers()`` reports the state at the last event Qt
    processed.  Both are tried, because the second is what a button-click slot
    sees, and the first is what a window built a moment later sees.
    """
    try:
        from pyirena.gui._qt import QApplication, Qt  # noqa: PLC0415

        if QApplication.instance() is None:
            return False
        shift = Qt.KeyboardModifier.ShiftModifier
        try:
            if bool(QApplication.queryKeyboardModifiers() & shift):
                return True
        except Exception:
            log.debug("queryKeyboardModifiers unavailable", exc_info=True)
        return bool(QApplication.keyboardModifiers() & shift)
    except Exception:
        log.debug("could not read keyboard modifiers", exc_info=True)
        return False


def reset_requested() -> bool:
    """True when the user is asking for default geometry.

    Either the Shift key is held right now (the gesture: hold Shift while
    opening the tool) or :data:`RESET_ENV_VAR` is set.
    """
    if os.environ.get(RESET_ENV_VAR, "").strip() not in ("", "0", "false", "False"):
        return True
    return shift_held()


def geometry_file():
    """Path of the geometry file — beside the main pyIrena state file."""
    from pyirena.state.state_manager import get_default_state_file  # noqa: PLC0415

    return get_default_state_file().parent / GEOMETRY_FILE_NAME


def _load_all() -> dict:
    """Every saved window entry; an unreadable or corrupt file reads as empty."""
    import json  # noqa: PLC0415

    try:
        path = geometry_file()
        if not path.exists():
            return {}
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        return data if isinstance(data, dict) else {}
    except Exception:
        log.debug("could not read window geometry", exc_info=True)
        return {}


def _write_all(entries: dict) -> None:
    import json  # noqa: PLC0415

    try:
        path = geometry_file()
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(entries, fh, indent=2)
    except Exception:
        log.debug("could not save window geometry", exc_info=True)


def _store(key: str, entry: Optional[dict]) -> None:
    """Write (or delete, when *entry* is None) one window's saved geometry."""
    entries = _load_all()
    if entry is None:
        entries.pop(key, None)
    else:
        entries[key] = entry
    _write_all(entries)


def clear_window_state(key: Optional[str] = None) -> None:
    """Forget one window's geometry, or all of them when *key* is None."""
    if key is None:
        _write_all({})
    else:
        _store(key, None)


def restore_window_state(widget, key: str,
                         splitters: Optional[Dict[str, object]] = None) -> bool:
    """Put *widget* back where it was, if that is still a sensible place.

    Args:
        widget: The window to position.
        key: Stable identifier for this window, e.g. ``"unified_fit_panel"``.
        splitters: Optional ``{name: QSplitter}`` whose pane sizes — the width
            of the control panel — are restored too.

    Returns:
        True when saved geometry was applied; False when the window should
        keep its coded default (nothing saved, reset requested, or the saved
        position is no longer usable).
    """
    if not GEOMETRY_PERSISTENCE_ENABLED:
        return False

    if reset_requested():
        log.info("Window geometry reset requested — using defaults for '%s'", key)
        clear_window_state(key)
        return False

    entry = _load_all().get(key)
    if not entry:
        return False

    screens = screen_rects()
    saved_fingerprint = entry.get("screens")
    now_fingerprint = _screen_fingerprint(screens)
    if saved_fingerprint and saved_fingerprint != now_fingerprint:
        log.info("Display setup changed since '%s' was saved (%s → %s)",
                 key, saved_fingerprint, now_fingerprint)

    rect = resolve_geometry(entry.get("geometry"), screens)
    applied = False
    if rect is not None and is_top_level(widget):
        try:
            widget.setGeometry(*rect)
            applied = True
        except Exception:
            log.debug("could not apply geometry for '%s'", key, exc_info=True)

    return applied


def _run_after_first_layout(widget, func) -> None:
    """Call *func* once *widget* has been shown and actually laid out.

    Splitter sizes are meaningless until then.  A panel sets its panes in the
    constructor while the splitter still has no width, so both *reading* the
    coded default and *writing* a remembered width have to wait for the first
    real layout — doing either early is what made a reset come back with the
    control panel the wrong width.
    """
    try:
        from pyirena.gui._qt import QEvent, QtCore, QTimer  # noqa: PLC0415

        def _later():
            # One more turn of the event loop: the Show event arrives before
            # the layout it triggers has run.
            QTimer.singleShot(0, lambda: func() if _still_alive(widget) else None)

        if widget.isVisible():
            _later()
            return

        class _OnFirstShow(QtCore.QObject):
            def eventFilter(self, obj, event):
                if event.type() == QEvent.Type.Show:
                    obj.removeEventFilter(self)
                    _later()
                return False

        watcher = _OnFirstShow(widget)
        widget.installEventFilter(watcher)
        widget._pyirena_first_layout_watcher = watcher
    except Exception:
        log.debug("could not defer to the first layout", exc_info=True)


def _restore_splitters(entry: dict, splitters, key: str) -> None:
    """Restore pane sizes, ignoring entries whose pane count has changed."""
    if not splitters:
        return
    saved = entry.get("splitters") or {}
    for name, splitter in splitters.items():
        sizes = saved.get(name)
        if not sizes or splitter is None:
            continue
        try:
            if len(sizes) != splitter.count() or sum(sizes) <= 0:
                # The panel gained or lost a pane since this was saved.
                continue
            splitter.setSizes([int(v) for v in sizes])
        except Exception:
            log.debug("could not restore splitter '%s' of '%s'", name, key,
                      exc_info=True)


def save_window_state(widget, key: str,
                      splitters: Optional[Dict[str, object]] = None) -> None:
    """Remember *widget*'s position, size and pane widths under *key*."""
    if not GEOMETRY_PERSISTENCE_ENABLED:
        return
    if not is_top_level(widget):
        # An embedded pane's geometry belongs to its layout, not to us.
        return
    try:
        if widget.isMinimized() or widget.isFullScreen():
            # Saving a minimised or full-screen rectangle would restore the
            # window in a state the user never chose.
            return
        geo = widget.geometry()
        entry = {
            "geometry": [geo.x(), geo.y(), geo.width(), geo.height()],
            "screens": _screen_fingerprint(screen_rects()),
        }
    except Exception:
        log.debug("could not read geometry of '%s'", key, exc_info=True)
        return

    if splitters:
        sizes = {}
        for name, splitter in splitters.items():
            try:
                if splitter is not None and sum(splitter.sizes()) > 0:
                    sizes[name] = list(splitter.sizes())
            except Exception:
                log.debug("could not read splitter '%s'", name, exc_info=True)
        if sizes:
            entry["splitters"] = sizes

    _store(key, entry)


#: Attribute holding everything needed to reset a window later: its key, its
#: splitters and the size/pane widths it had *before* any saved geometry was
#: applied — i.e. the tool's own coded default.
_INFO_ATTR = "_pyirena_window_state"


def centred_default(widget, size: Tuple[int, int]) -> Optional[Rect]:
    """Where a window of *size* should sit when it has no remembered position.

    Centred on the screen the window is currently on (the one the user is
    looking at), clamped so it always fits.  Returns None when there is no
    screen information — the caller then leaves the window alone.
    """
    screens = screen_rects()
    if not screens:
        return None

    current: Optional[Rect] = None
    try:
        screen = widget.screen()
        if screen is not None:
            geo = screen.availableGeometry()
            current = (geo.x(), geo.y(), geo.width(), geo.height())
    except Exception:
        log.debug("could not read the window's screen", exc_info=True)
    if current is None or current not in screens:
        current = screens[0]

    sx, sy, sw, sh = current
    w = min(max(int(size[0]), MIN_WINDOW_PX[0]), sw - 2 * EDGE_MARGIN_PX)
    h = min(max(int(size[1]), MIN_WINDOW_PX[1]), sh - 2 * EDGE_MARGIN_PX)
    return (int(sx + (sw - w) / 2), int(sy + (sh - h) / 2), w, h)


def _default_size(widget, info: dict) -> Tuple[int, int]:
    """The size a tool opens at when nothing is remembered.

    The size the panel asked for — ``self.resize(...)`` in its constructor,
    recorded at install time — is the answer whenever it is usable.  The size
    hint is only a fallback: for a panel full of nested group boxes it can be
    far larger than anything the author intended, which would make "reset"
    produce a window bigger than the one being reset.
    """
    recorded = tuple(int(v) for v in (info.get("default_size") or (0, 0)))
    if recorded[0] >= MIN_WINDOW_PX[0] and recorded[1] >= MIN_WINDOW_PX[1]:
        return (recorded[0], recorded[1])

    width = height = 0
    try:
        hint = widget.sizeHint()
        width = max(hint.width(), widget.minimumWidth())
        height = max(hint.height(), widget.minimumHeight())
    except Exception:
        log.debug("could not read size hint", exc_info=True)
    return (max(width, MIN_WINDOW_PX[0]), max(height, MIN_WINDOW_PX[1]))


def _splitter_fractions(splitter) -> Optional[List[float]]:
    """A splitter's pane sizes as fractions of its total — or None if unusable.

    Fractions rather than pixels on purpose: a reset re-centres the window at
    its default *size*, which is rarely the size the splitter has right now, and
    pixel sizes applied to a differently-sized splitter land wrong.  Qt scales
    what you pass to ``setSizes`` anyway, so the ratio is the real information.
    """
    try:
        sizes = [float(v) for v in splitter.sizes()]
    except Exception:
        log.debug("could not read splitter sizes", exc_info=True)
        return None
    total = sum(sizes)
    if total <= 0 or any(s < 0 for s in sizes):
        return None
    return [s / total for s in sizes]


def _apply_fractions(splitter, fractions: Sequence[float]) -> bool:
    """Set pane sizes from fractions of the splitter's current extent."""
    try:
        if splitter is None or len(fractions) != splitter.count():
            return False
        extent = splitter.width() if splitter.orientation() == _horizontal() \
            else splitter.height()
        if extent <= 0:
            extent = sum(splitter.sizes())
        if extent <= 0:
            return False
        splitter.setSizes([max(1, int(round(f * extent))) for f in fractions])
        return True
    except Exception:
        log.debug("could not apply splitter fractions", exc_info=True)
        return False


def _horizontal():
    from pyirena.gui._qt import Qt  # noqa: PLC0415

    return Qt.Orientation.Horizontal


def _linked_windows(widget) -> List[object]:
    """Other windows this one owns — a panel's separate graph window.

    A tool can be more than one window (Unified Fit and WAXS each have a
    detachable graph), and resetting the panel while its graph stays lost on a
    vanished monitor would only half-answer the user.  Attributes are scanned
    rather than the Qt parent tree because those windows are deliberately
    parentless.
    """
    found = []
    try:
        for value in list(vars(widget).values()):
            if value is widget or value is None:
                continue
            if getattr(value, _INFO_ATTR, None) is None:
                continue
            if not is_top_level(value):
                # A pane, not a window — see is_top_level().
                continue
            found.append(value)
    except Exception:
        log.debug("could not scan for linked windows", exc_info=True)
    return found


def reset_window(widget, include_linked: bool = True) -> bool:
    """Forget this window's saved geometry and put it back to its default.

    Works on a window that is already built and open — which is the point: tool
    panels are created once and reused, so a reset gesture that only worked at
    construction time would work once per session.

    Args:
        widget: A window that went through :func:`install_window_state`.
        include_linked: Also reset windows the panel owns (its graph window).

    Returns:
        True when something was reset.
    """
    if not GEOMETRY_PERSISTENCE_ENABLED:
        return False

    info = getattr(widget, _INFO_ATTR, None)
    if not info:
        return False

    key = info.get("key")
    if key:
        clear_window_state(key)

    _apply_defaults(widget, info, key)

    # Apply again once the event loop has run.  Showing a window triggers a
    # layout pass, and on some window managers a placement pass, either of
    # which can undo what we just set — the pane widths in particular are only
    # meaningful once the splitter has its real width.  The second pass is
    # cheap and makes the result the same whether the window was already open
    # or is about to be shown for the first time.
    try:
        from pyirena.gui._qt import QTimer  # noqa: PLC0415

        QTimer.singleShot(0, lambda: _apply_defaults(widget, info, key))
    except Exception:
        log.debug("could not schedule the second reset pass", exc_info=True)

    if include_linked:
        for other in _linked_windows(widget):
            reset_window(other, include_linked=False)

    log.info("Window geometry reset to default for '%s'", key)
    return True


def _apply_defaults(widget, info: dict, key: Optional[str]) -> None:
    """Put one window back to its recorded default size, position and panes."""
    try:
        if not _still_alive(widget):
            return
        if is_top_level(widget):
            rect = centred_default(widget, _default_size(widget, info))
            if rect is not None:
                widget.setGeometry(*rect)
        for name, fractions in (info.get("default_splitters") or {}).items():
            _apply_fractions((info.get("splitters") or {}).get(name), fractions)
    except Exception:
        log.debug("could not apply defaults for '%s'", key, exc_info=True)


def _still_alive(widget) -> bool:
    """False once Qt has destroyed the C++ object behind *widget*."""
    try:
        widget.objectName()
        return True
    except RuntimeError:
        return False
    except Exception:
        return True


def reset_window_if_shift(widget) -> bool:
    """Reset *widget* when the user is holding Shift as the tool is launched.

    One call at the launch site is the whole gesture::

        reset_window_if_shift(self.unified_fit_window)
        self.unified_fit_window.show()

    Safe to call with anything: a window that never installed geometry
    persistence, or None, simply returns False.
    """
    if widget is None or not GEOMETRY_PERSISTENCE_ENABLED:
        return False
    if not shift_held():
        return False
    return reset_window(widget)


def install_window_state(widget, key: str,
                         splitters: Optional[Dict[str, object]] = None) -> bool:
    """Restore now and arrange to save on close — the one call a panel needs.

    ::

        install_window_state(self, "unified_fit_panel",
                             splitters={"main": main_splitter})

    Saving happens on the window's close event through an event filter, so a
    panel does not need a ``closeEvent`` (and an existing one is untouched).

    Returns:
        Whether saved geometry was applied.
    """
    if not GEOMETRY_PERSISTENCE_ENABLED:
        return False

    # The size the panel asked for in its constructor is the default to reset
    # to; unlike pane widths it is readable straight away.
    try:
        size = widget.size()
        default_size = (size.width(), size.height())
    except Exception:
        default_size = (0, 0)
    info = {
        "key": key,
        "splitters": splitters or {},
        "default_size": default_size,
        "default_splitters": {},     # filled in at the first layout, below
    }
    setattr(widget, _INFO_ATTR, info)

    applied = restore_window_state(widget, key, splitters)

    def _first_layout() -> None:
        """Record the coded pane widths, then apply the remembered ones."""
        if not info["default_splitters"]:
            for name, splitter in (splitters or {}).items():
                fractions = _splitter_fractions(splitter) if splitter is not None else None
                if fractions:
                    info["default_splitters"][name] = fractions
        # Re-read rather than closing over the entry: a Shift launch may have
        # cleared it between construction and the first layout, and in that
        # case there is deliberately nothing to restore.
        entry = _load_all().get(key)
        if entry:
            _restore_splitters(entry, splitters, key)

    if splitters:
        _run_after_first_layout(widget, _first_layout)

    try:
        from pyirena.gui._qt import QEvent, QtCore  # noqa: PLC0415

        class _SaveOnClose(QtCore.QObject):
            def eventFilter(self, obj, event):
                if event.type() == QEvent.Type.Close:
                    save_window_state(obj, key, splitters)
                return False

        keeper = _SaveOnClose(widget)
        widget.installEventFilter(keeper)
        # Keep a reference so the filter lives as long as the window.
        widget._pyirena_window_state_filter = keeper
    except Exception:
        log.debug("could not install save-on-close for '%s'", key, exc_info=True)

    return applied
