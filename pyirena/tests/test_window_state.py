"""
Tests for :mod:`pyirena.gui.window_state` — remembering window geometry safely.

The dangerous case is not "did it remember", it is "the display setup changed
and the window came back somewhere unreachable".  Undocking a laptop, unplugging
a second monitor or changing resolution all leave saved coordinates that no
screen covers, and a window restored there cannot be seen *or* dragged back.

:func:`resolve_geometry` is the whole placement policy and is pure — plain
``(x, y, w, h)`` tuples in, a rectangle or None out — so every awkward case is
tested here with no display at all.  The Qt glue is smoke-tested separately.
"""

from __future__ import annotations

import pytest

from pyirena.gui.window_state import (
    EDGE_MARGIN_PX,
    MIN_VISIBLE_PX,
    resolve_geometry,
)

# A common two-monitor setup: laptop plus an external screen to its right.
LAPTOP = (0, 0, 1680, 1050)
EXTERNAL = (1680, 0, 2560, 1440)
BOTH = [LAPTOP, EXTERNAL]


# ── The ordinary case ────────────────────────────────────────────────────

def test_a_window_still_on_screen_is_left_alone():
    saved = (100, 100, 1200, 800)
    assert resolve_geometry(saved, [LAPTOP]) == saved


def test_a_window_on_the_second_screen_is_left_alone():
    saved = (2000, 200, 1200, 900)
    assert resolve_geometry(saved, BOTH) == saved


def test_a_window_spanning_two_screens_is_left_alone():
    """Users park a window across the seam; that is their choice to keep."""
    saved = (1400, 100, 900, 700)
    assert resolve_geometry(saved, BOTH) == saved


# ── Nothing to restore ───────────────────────────────────────────────────

@pytest.mark.parametrize("saved", [None, (), (1, 2, 3)])
def test_missing_or_malformed_geometry_falls_back_to_the_default(saved):
    assert resolve_geometry(saved, [LAPTOP]) is None


def test_no_screens_falls_back_to_the_default():
    """Headless or mid-reconfiguration — never guess."""
    assert resolve_geometry((100, 100, 800, 600), []) is None


def test_a_collapsed_window_is_not_restored():
    """A zero or tiny rectangle is corruption, not a user preference."""
    assert resolve_geometry((100, 100, 10, 10), [LAPTOP]) is None


def test_non_numeric_entries_are_ignored():
    assert resolve_geometry(("x", "y", "w", "h"), [LAPTOP]) is None


# ── The display changed ──────────────────────────────────────────────────

def test_window_on_a_removed_screen_comes_back_to_the_remaining_one():
    """The external monitor is gone; the window was entirely on it."""
    saved = (2400, 300, 1200, 900)
    out = resolve_geometry(saved, [LAPTOP])
    assert out is not None
    x, y, w, h = out
    # Fully inside the laptop screen…
    assert x >= LAPTOP[0] and y >= LAPTOP[1]
    assert x + w <= LAPTOP[0] + LAPTOP[2]
    assert y + h <= LAPTOP[1] + LAPTOP[3]
    # …and no bigger than it was.
    assert w <= saved[2] and h <= saved[3]


def test_window_larger_than_the_remaining_screen_is_shrunk_to_fit():
    saved = (2000, 0, 2400, 1300)          # sized for the external monitor
    x, y, w, h = resolve_geometry(saved, [LAPTOP])
    assert w <= LAPTOP[2] - 2 * EDGE_MARGIN_PX
    assert h <= LAPTOP[3] - 2 * EDGE_MARGIN_PX
    assert x + w <= LAPTOP[0] + LAPTOP[2]
    assert y + h <= LAPTOP[1] + LAPTOP[3]


def test_a_window_hanging_off_the_bottom_is_pulled_back():
    saved = (200, 1040, 900, 700)          # only a sliver visible
    out = resolve_geometry(saved, [LAPTOP])
    assert out is not None
    _, y, _, h = out
    assert y + h <= LAPTOP[1] + LAPTOP[3]


def test_a_window_whose_title_bar_is_above_the_screen_is_pulled_back():
    """The body is visible but the title bar is not — it cannot be dragged."""
    saved = (300, -200, 900, 700)
    out = resolve_geometry(saved, [LAPTOP])
    assert out is not None and out[1] >= LAPTOP[1]


def test_a_barely_overlapping_window_is_relocated():
    """A few pixels on screen is not enough to grab."""
    saved = (LAPTOP[2] - (MIN_VISIBLE_PX // 2), 100, 900, 700)
    out = resolve_geometry(saved, [LAPTOP])
    assert out is not None
    assert out[0] + out[2] <= LAPTOP[0] + LAPTOP[2]


def test_relocation_prefers_the_nearest_screen():
    """A window just off the right-hand screen returns to that one."""
    saved = (4300, 200, 800, 600)          # just beyond the external monitor
    x, _, w, _ = resolve_geometry(saved, BOTH)
    assert x >= EXTERNAL[0]
    assert x + w <= EXTERNAL[0] + EXTERNAL[2]


def test_negative_screen_origins_are_handled():
    """A monitor placed left of the primary gives negative coordinates."""
    left_screen = (-1920, 0, 1920, 1080)
    saved = (-1800, 100, 900, 700)
    assert resolve_geometry(saved, [left_screen, LAPTOP]) == saved

    lost = (-5000, 100, 900, 700)
    out = resolve_geometry(lost, [left_screen, LAPTOP])
    assert out is not None and out[0] >= left_screen[0]


def test_every_relocation_lands_somewhere_visible():
    """The property that actually matters, over a spread of bad positions."""
    from pyirena.gui.window_state import _is_comfortably_visible

    screens = BOTH
    for x in (-6000, -900, 0, 900, 4000, 9000):
        for y in (-3000, -100, 0, 500, 3000):
            out = resolve_geometry((x, y, 1000, 800), screens)
            if out is None:
                continue
            assert any(_is_comfortably_visible(out, s) for s in screens), (
                f"({x}, {y}) resolved to {out}, which is not usable"
            )


# ── Storage ──────────────────────────────────────────────────────────────

@pytest.fixture
def isolated_store(tmp_path, monkeypatch):
    """Point the geometry file at a temp dir — never touch the real one."""
    from pyirena.gui import window_state as ws

    target = tmp_path / "window_geometry.json"
    monkeypatch.setattr(ws, "geometry_file", lambda: target)
    monkeypatch.delenv(ws.RESET_ENV_VAR, raising=False)
    return target


def test_geometry_lives_outside_the_shared_state_file(isolated_store):
    """Its own file, deliberately.

    Every panel owns a StateManager instance whose save() rewrites the whole
    state file from an in-memory copy, so geometry written while closing one
    panel was erased by the next panel's save. A separate file also makes
    "delete this to reset" a one-liner.
    """
    from pyirena.gui import window_state as ws

    ws._store("some_panel", {"geometry": [10, 20, 900, 700]})
    assert isolated_store.exists()
    assert ws._load_all()["some_panel"]["geometry"] == [10, 20, 900, 700]

    from pyirena.state.state_manager import get_default_state_file

    assert isolated_store != get_default_state_file()


def test_a_corrupt_geometry_file_reads_as_empty(isolated_store):
    from pyirena.gui import window_state as ws

    isolated_store.write_text("{not json at all")
    assert ws._load_all() == {}


def test_clear_removes_one_entry_or_all(isolated_store):
    from pyirena.gui import window_state as ws

    ws._store("a", {"geometry": [0, 0, 800, 600]})
    ws._store("b", {"geometry": [0, 0, 800, 600]})
    ws.clear_window_state("a")
    assert set(ws._load_all()) == {"b"}
    ws.clear_window_state()
    assert ws._load_all() == {}


# ── Qt glue (skipped without a Qt binding) ───────────────────────────────

def _require_qt() -> None:
    try:
        import pyirena.gui._qt  # noqa: F401
    except ImportError:
        pytest.skip("Qt (PySide6/PyQt6) not available", allow_module_level=True)


def test_reset_via_environment_variable(monkeypatch):
    from pyirena.gui import window_state as ws

    monkeypatch.setenv(ws.RESET_ENV_VAR, "1")
    assert ws.reset_requested() is True
    monkeypatch.setenv(ws.RESET_ENV_VAR, "0")
    assert ws.reset_requested() is False
    monkeypatch.delenv(ws.RESET_ENV_VAR, raising=False)


def test_feature_switch_disables_every_entry_point(monkeypatch):
    """One flag must be enough to ship without the feature."""
    from pyirena.gui import window_state as ws

    monkeypatch.setattr(ws, "GEOMETRY_PERSISTENCE_ENABLED", False)

    calls = []
    monkeypatch.setattr(ws, "_store", lambda *a, **k: calls.append(a))

    class _FakeWindow:
        def setGeometry(self, *args):
            calls.append(("setGeometry", args))

    assert ws.restore_window_state(_FakeWindow(), "k") is False
    ws.save_window_state(_FakeWindow(), "k")
    assert ws.install_window_state(_FakeWindow(), "k") is False
    assert calls == []


# ── Live windows (offscreen Qt) ──────────────────────────────────────────

def _qt_or_skip():
    try:
        from pyirena.gui._qt import QApplication
    except ImportError:
        pytest.skip("Qt (PySide6/PyQt6) not available")
    return QApplication.instance() or QApplication([])


def test_a_window_reopens_where_it_was_left(isolated_store):
    """The whole point, end to end: move, close, reopen."""
    app = _qt_or_skip()
    from pyirena.gui._qt import QSplitter, QWidget
    from pyirena.gui.window_state import install_window_state

    first = QWidget()
    first.resize(640, 480)
    install_window_state(first, "round_trip")
    first.move(70, 55)
    first.close()
    app.processEvents()

    second = QWidget()
    second.resize(300, 200)          # a different default
    install_window_state(second, "round_trip")
    geo = second.geometry()
    assert (geo.x(), geo.y()) == (70, 55)
    assert (geo.width(), geo.height()) == (640, 480)
    first.deleteLater()
    second.deleteLater()
    del QSplitter                    # keep the import meaningful to linters


def _make_split_window(pane_sizes=(300, 600), panes=2):
    """A window laid out the way a real panel is: splitter inside a layout.

    Built and *shown* through the same path a panel takes, because splitter
    sizes only become real at the first layout — a test that skips the show
    tests nothing about panes.
    """
    from pyirena.gui._qt import QSplitter, Qt, QVBoxLayout, QWidget

    window = QWidget()
    window.resize(900, 600)
    layout = QVBoxLayout(window)
    layout.setContentsMargins(0, 0, 0, 0)
    splitter = QSplitter(Qt.Orientation.Horizontal)
    for _ in range(panes):
        splitter.addWidget(QWidget())
    layout.addWidget(splitter)
    splitter.setSizes(list(pane_sizes))
    return window, splitter


def _settle(app, times=6):
    """Let the deferred first-layout pass run."""
    for _ in range(times):
        app.processEvents()


def test_control_panel_width_is_remembered(isolated_store):
    """The splitter position is the third thing users expect to persist."""
    app = _qt_or_skip()
    from pyirena.gui.window_state import install_window_state

    window, splitter = _make_split_window((300, 600))
    install_window_state(window, "splitter_case", splitters={"main": splitter})
    window.show()
    _settle(app)
    splitter.setSizes([300, 600])            # the width the user chose
    _settle(app)
    window.close()
    _settle(app)

    window2, splitter2 = _make_split_window((450, 450))
    install_window_state(window2, "splitter_case", splitters={"main": splitter2})
    window2.show()
    _settle(app)
    # Qt reserves a few pixels for the handle, so the restored width is close
    # rather than exact — what matters is that it is the saved width, not the
    # 450 this window was built with.
    assert abs(splitter2.sizes()[0] - 300) <= 5
    window.deleteLater()
    window2.deleteLater()


def test_a_splitter_that_gained_a_pane_is_left_alone(isolated_store):
    """Saved sizes from an older layout must not be forced onto a new one."""
    app = _qt_or_skip()
    from pyirena.gui._qt import QSplitter, Qt, QWidget
    from pyirena.gui.window_state import install_window_state, save_window_state

    old = QWidget()
    old_splitter = QSplitter(Qt.Orientation.Horizontal, old)
    for _ in range(2):
        old_splitter.addWidget(QWidget())
    old_splitter.setSizes([200, 400])
    save_window_state(old, "pane_count", splitters={"main": old_splitter})
    app.processEvents()

    new, new_splitter = _make_split_window((300, 300, 300), panes=3)  # a pane was added
    install_window_state(new, "pane_count", splitters={"main": new_splitter})
    new.show()
    _settle(app)
    assert len(new_splitter.sizes()) == 3    # untouched, not crashed
    old.deleteLater()
    new.deleteLater()


def test_shift_at_open_restores_the_default_and_forgets(isolated_store, monkeypatch):
    app = _qt_or_skip()
    from pyirena.gui import window_state as ws
    from pyirena.gui._qt import QWidget

    saved = QWidget()
    saved.resize(700, 500)
    ws.install_window_state(saved, "reset_case")
    saved.move(200, 150)
    saved.close()
    app.processEvents()
    assert "reset_case" in ws._load_all()

    monkeypatch.setattr(ws, "reset_requested", lambda: True)
    fresh = QWidget()
    fresh.resize(400, 300)
    assert ws.install_window_state(fresh, "reset_case") is False
    assert "reset_case" not in ws._load_all()   # and forgotten, not just skipped
    saved.deleteLater()
    fresh.deleteLater()


# ── Shift-click on a tool button resets that one tool ─────────────────────

def test_shift_clicking_a_tool_button_resets_an_already_open_window(isolated_store,
                                                                    monkeypatch):
    """The Irena gesture, on the case that actually matters.

    Tool panels are built once and reused, so a reset that only ran at
    construction would work once per session.  This resets a window that is
    already open: the saved entry goes, and the window returns to the size it
    was coded with — not the size the user had dragged it to.
    """
    app = _qt_or_skip()
    from pyirena.gui import window_state as ws
    from pyirena.gui._qt import QWidget

    win = QWidget()
    win.resize(700, 500)            # the tool's coded default
    ws.install_window_state(win, "shift_reset")

    win.setGeometry(120, 90, 1100, 820)   # the user drags and resizes it
    win.close()
    app.processEvents()
    assert "shift_reset" in ws._load_all()

    monkeypatch.setattr(ws, "shift_held", lambda: True)
    assert ws.reset_window_if_shift(win) is True

    assert "shift_reset" not in ws._load_all()
    assert (win.width(), win.height()) == (700, 500)
    win.deleteLater()


def test_without_shift_a_launch_leaves_the_window_alone(isolated_store, monkeypatch):
    app = _qt_or_skip()
    from pyirena.gui import window_state as ws
    from pyirena.gui._qt import QWidget

    win = QWidget()
    win.resize(700, 500)
    ws.install_window_state(win, "no_shift")
    win.setGeometry(120, 90, 1100, 820)
    win.close()
    app.processEvents()

    monkeypatch.setattr(ws, "shift_held", lambda: False)
    assert ws.reset_window_if_shift(win) is False
    assert "no_shift" in ws._load_all()
    assert (win.width(), win.height()) == (1100, 820)
    win.deleteLater()


def test_resetting_is_safe_on_anything_a_launcher_might_hold(monkeypatch):
    """Launch sites call this on whatever they have; it must never raise."""
    _qt_or_skip()
    from pyirena.gui import window_state as ws
    from pyirena.gui._qt import QWidget

    monkeypatch.setattr(ws, "shift_held", lambda: True)

    assert ws.reset_window_if_shift(None) is False
    plain = QWidget()                      # never installed
    assert ws.reset_window_if_shift(plain) is False
    plain.deleteLater()


def test_reset_restores_the_control_panel_width_too(isolated_store, monkeypatch):
    app = _qt_or_skip()
    from pyirena.gui import window_state as ws

    win, splitter = _make_split_window((300, 600))   # the tool's coded panes
    ws.install_window_state(win, "reset_panes", splitters={"main": splitter})
    win.show()
    _settle(app)
    default_left = splitter.sizes()[0]

    splitter.setSizes([550, 350])          # user drags the divider
    _settle(app)
    win.close()
    _settle(app)
    assert splitter.sizes()[0] != default_left

    monkeypatch.setattr(ws, "shift_held", lambda: True)
    ws.reset_window_if_shift(win)
    _settle(app)
    assert abs(splitter.sizes()[0] - default_left) <= 5
    win.deleteLater()


def test_resetting_a_panel_also_resets_its_graph_window(isolated_store, monkeypatch):
    """Unified Fit and WAXS are two windows; a half reset is not a reset."""
    app = _qt_or_skip()
    from pyirena.gui import window_state as ws
    from pyirena.gui._qt import QWidget

    panel = QWidget()
    panel.resize(800, 600)
    graph = QWidget()                       # parentless, like the real ones
    graph.resize(500, 400)
    panel.graph = graph

    ws.install_window_state(graph, "linked_graph")
    ws.install_window_state(panel, "linked_panel")
    panel.setGeometry(50, 50, 1200, 900)
    graph.setGeometry(2000, 40, 900, 700)
    panel.close()
    graph.close()
    app.processEvents()
    assert {"linked_panel", "linked_graph"} <= set(ws._load_all())

    monkeypatch.setattr(ws, "shift_held", lambda: True)
    ws.reset_window_if_shift(panel)

    saved = set(ws._load_all())
    assert "linked_panel" not in saved and "linked_graph" not in saved
    assert (graph.width(), graph.height()) == (500, 400)
    panel.deleteLater()
    graph.deleteLater()


def test_a_reset_window_lands_fully_on_a_real_screen(isolated_store, monkeypatch):
    _qt_or_skip()
    from pyirena.gui import window_state as ws
    from pyirena.gui._qt import QWidget

    win = QWidget()
    win.resize(600, 450)
    ws.install_window_state(win, "reset_placement")
    win.setGeometry(9000, 7000, 600, 450)   # off on a monitor that is gone
    monkeypatch.setattr(ws, "shift_held", lambda: True)
    ws.reset_window_if_shift(win)

    rect = (win.x(), win.y(), win.width(), win.height())
    screens = ws.screen_rects()
    assert any(ws._is_comfortably_visible(rect, s) for s in screens)
    win.deleteLater()


def test_a_window_bigger_than_the_screen_is_shrunk_when_reset(isolated_store,
                                                              monkeypatch):
    """A default coded for a large monitor must still fit a laptop screen."""
    _qt_or_skip()
    from pyirena.gui import window_state as ws

    class _Fake:
        def __init__(self):
            self._geo = None

        def size(self):
            class _S:
                def width(self_inner): return 4000
                def height(self_inner): return 3000
            return _S()

        def sizeHint(self):
            return self.size()

        def minimumWidth(self): return 0
        def minimumHeight(self): return 0
        def screen(self): return None
        def setGeometry(self, *rect): self._geo = rect

    monkeypatch.setattr(ws, "screen_rects", lambda: [(0, 0, 1440, 900)])
    fake = _Fake()
    setattr(fake, ws._INFO_ATTR, {"key": "huge", "splitters": {},
                                  "default_size": (4000, 3000),
                                  "default_splitters": {}})
    assert ws.reset_window(fake) is True
    x, y, w, h = fake._geo
    assert w <= 1440 and h <= 900 and x >= 0 and y >= 0


def test_every_tool_launcher_honours_the_shift_gesture():
    """A new tool must not quietly miss the reset gesture.

    The wiring is one line per launch site, which is exactly the kind of line
    that gets forgotten when a tool is added, so this reads the source rather
    than trusting review.
    """
    import re
    from pathlib import Path

    import pyirena

    source = (Path(pyirena.__file__).parent / "gui" / "data_selector" / "panel.py").read_text(encoding="utf-8")
    shown = set(re.findall(r"self\.(\w+_window)\.show\(\)", source))
    reset = set(re.findall(r"reset_window_if_shift\(self\.(\w+_window)\)", source))

    # The plot window opened by double-clicking a file is deliberately not in
    # this list: Shift-click in a file list means "extend the selection".
    shown.discard("graph_window")

    assert shown, "no tool launchers found — has the panel been restructured?"
    assert shown <= reset, f"tool windows shown without a Shift reset: {sorted(shown - reset)}"


# ── Panes are not windows ────────────────────────────────────────────────

def test_an_embedded_pane_is_never_moved_or_saved(isolated_store):
    """Registering a pane by mistake must be harmless.

    Several pyIrena classes are named ``…GraphWindow`` but are added to a
    splitter as the right-hand pane.  Applying a remembered rectangle to one of
    those moves it inside the layout, which shows up as the *control panel*
    changing width — the bug this test exists for.
    """
    app = _qt_or_skip()
    from pyirena.gui import window_state as ws

    window, splitter = _make_split_window((300, 600))
    pane = splitter.widget(1)
    window.show()
    _settle(app)

    ws._store("stray_pane", {"geometry": [40, 40, 500, 400]})
    assert ws.restore_window_state(pane, "stray_pane") is False   # not moved

    before = splitter.sizes()
    ws.save_window_state(pane, "stray_pane_2")
    assert "stray_pane_2" not in ws._load_all()                   # not saved
    assert splitter.sizes() == before
    window.deleteLater()


def test_resetting_a_panel_leaves_its_panes_at_the_coded_default(isolated_store,
                                                                 monkeypatch):
    """End to end on the real panel, which is where this went wrong.

    Unified Fit builds its graph as the right pane of a splitter and sets the
    control panel to a third of the width.  After a Shift launch the panes must
    be back at that third — not at whatever the user last dragged, and not at
    some ratio measured before the window had ever been laid out.
    """
    app = _qt_or_skip()
    pytest.importorskip("pyqtgraph")
    from pyirena.gui import window_state as ws
    from pyirena.gui.unified_fit import UnifiedFitPanel

    panel = UnifiedFitPanel()
    panel.show()
    _settle(app)
    default_panes = list(panel.main_splitter.sizes())
    assert default_panes[0] < default_panes[1], "control panel should be the narrow one"

    panel.main_splitter.setSizes([900, 300])      # user widens the controls
    panel.setGeometry(60, 40, 1400, 1000)
    _settle(app)
    panel.close()
    _settle(app)
    assert "unified_fit_panel" in ws._load_all()

    monkeypatch.setattr(ws, "shift_held", lambda: True)
    ws.reset_window_if_shift(panel)
    panel.show()
    _settle(app)

    assert "unified_fit_panel" not in ws._load_all()
    assert abs(panel.main_splitter.sizes()[0] - default_panes[0]) <= 5
    panel.deleteLater()
    _settle(app)
