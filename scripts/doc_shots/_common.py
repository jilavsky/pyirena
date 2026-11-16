"""Shared helpers for the documentation screenshot script.

Used by ``scripts/make_doc_screenshots.py`` and the per-page modules in this
folder. Nothing here is imported by pyIrena itself. The few Qt classes that
``pyirena/gui/_qt.py`` does not re-export (QPoint, QRect, QPainter) come from
PySide6 directly: this folder is outside the package, and pyIrena's GUI is
PySide6 only.

The script must not touch the user's settings or the tracked data:

* ``PYIRENA_STATE_DIR`` points at a throw-away folder, so panels start from
  their coded defaults and nothing is written to ``~/.pyirena``;
* window-geometry persistence is switched off;
* every dataset is copied into a scratch folder before a panel opens it, so
  "Store in File" or an auto-created NXcanSAS sibling never lands in
  ``validationData/`` or ``testData/``.
"""

from __future__ import annotations

import os
import shutil
import tempfile
import time
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DOCS_IMG = REPO / "docs" / "tools" / "img"
VALIDATION = REPO / "validationData"
TEST_DATA = REPO / "testData"

#: Callout colour: pyIrena's orange accent, readable on white and on plots.
CALLOUT_FILL = "#c2410c"
CALLOUT_RADIUS = 12

#: Shot id -> registered function.
REGISTRY: dict[str, Shot] = {}


@dataclass
class Shot:
    """One registered screenshot."""

    shot_id: str
    func: Callable[[Context], object]
    page: str = field(init=False)

    def __post_init__(self) -> None:
        # "unified-fit-03-guinier-fit" -> page "unified-fit": the page is
        # everything before the two-digit sequence number.
        parts = self.shot_id.split("-")
        for i, part in enumerate(parts):
            if len(part) == 2 and part.isdigit():
                self.page = "-".join(parts[:i])
                break
        else:
            raise ValueError(f"shot id without <nn>: {self.shot_id}")

    @property
    def path(self) -> Path:
        return DOCS_IMG / self.page / f"{self.shot_id}.png"


def shot(shot_id: str):
    """Register a function that returns the image (QPixmap/QImage) for *shot_id*."""

    def deco(func):
        if shot_id in REGISTRY:
            raise ValueError(f"duplicate shot id {shot_id}")
        REGISTRY[shot_id] = Shot(shot_id, func)
        return func

    return deco


class Context:
    """Per-run state handed to every shot function."""

    def __init__(self, app, scratch: Path):
        self.app = app
        self.scratch = scratch

    def dataset(self, relpath: str, root: Path = VALIDATION) -> str:
        """Copy a tracked dataset into the scratch folder; return the copy's path.

        A fresh copy for every call, so a shot that stores results cannot
        change what the next shot sees.
        """
        src = root / relpath
        sub = Path(tempfile.mkdtemp(dir=self.scratch))
        dst = sub / src.name
        shutil.copy2(src, dst)
        return str(dst)

    def settle(self, rounds: int = 30) -> None:
        """Let Qt process pending events (layouts, timers, repaints)."""
        for _ in range(rounds):
            self.app.processEvents()
            time.sleep(0.005)

    def wait_until(self, cond: Callable[[], bool], timeout: float = 120.0) -> None:
        """Process events until *cond()* is true; raise on timeout."""
        t0 = time.monotonic()
        while not cond():
            self.app.processEvents()
            time.sleep(0.01)
            if time.monotonic() - t0 > timeout:
                raise TimeoutError("condition not reached")
        self.settle()

    def wait_jobs(self, panel, timeout: float = 300.0) -> None:
        """Wait for a panel's background fit/MC job (FitJobRunner) to finish."""
        jobs = getattr(panel, "_jobs", None)
        self.settle(5)
        if jobs is not None:
            self.wait_until(lambda: not jobs.busy, timeout)

    def show(self, widget, size: tuple[int, int]) -> None:
        """Show *widget* at a fixed size and let it lay out."""
        widget.resize(*size)
        widget.show()
        self.settle()


def widget_rect(widget, root):
    """Rectangle of *widget* in *root*'s coordinates (QRect)."""
    from PySide6.QtCore import QPoint, QRect

    top_left = widget.mapTo(root, QPoint(0, 0))
    return QRect(top_left, widget.size())


def scene_item_rect(view, item, root):
    """Rectangle of a pyqtgraph scene item (e.g. a PlotItem) in *root* coordinates."""
    from PySide6.QtCore import QPoint, QRect

    poly = view.mapFromScene(item.sceneBoundingRect())
    r = poly.boundingRect()
    top_left = view.mapTo(root, QPoint(r.x(), r.y()))
    return QRect(top_left, r.size())


def draw_callouts(image, marks: Sequence[tuple]):
    """Draw numbered circles on *image* (QPixmap) and return it.

    Args:
        image: The grabbed window.
        marks: ``(number, QRect, anchor)`` or ``(number, QRect, anchor, (dx, dy))``;
            *anchor* is where on the rectangle the circle sits: ``"tl"``,
            ``"tr"``, ``"bl"``, ``"br"``, ``"l"``, ``"r"``, ``"t"``, ``"c"``;
            the optional offset moves it, e.g. ``(16, 0)`` puts an ``"r"``
            circle just outside the widget. It is kept inside the image.
    """
    from PySide6.QtGui import QPainter

    from pyirena.gui._qt import QBrush, QColor, QFont, QPen, QPointF, Qt

    p = QPainter(image)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    font = QFont()
    font.setBold(True)
    font.setPixelSize(13)
    p.setFont(font)
    rad = CALLOUT_RADIUS
    w, h = image.width() / image.devicePixelRatio(), image.height() / image.devicePixelRatio()
    for number, rect, anchor, *offset in marks:
        dx, dy = offset[0] if offset else (0, 0)
        x = {"l": rect.left(), "tl": rect.left(), "bl": rect.left(),
             "r": rect.right(), "tr": rect.right(), "br": rect.right()}.get(
            anchor, rect.center().x())
        y = {"t": rect.top(), "tl": rect.top(), "tr": rect.top(),
             "bl": rect.bottom(), "br": rect.bottom()}.get(anchor, rect.center().y())
        x, y = x + dx, y + dy
        x = min(max(x, rad + 1), w - rad - 1)
        y = min(max(y, rad + 1), h - rad - 1)
        centre = QPointF(x, y)
        p.setPen(QPen(QColor("white"), 2))
        p.setBrush(QBrush(QColor(CALLOUT_FILL)))
        p.drawEllipse(centre, rad, rad)
        p.setPen(QColor("white"))
        p.drawText(int(x - rad), int(y - rad), 2 * rad, 2 * rad,
                   int(Qt.AlignmentFlag.AlignCenter), str(number))
    p.end()
    return image


def make_app():
    """Create the offscreen QApplication with pyIrena's light theme."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ.pop("PYIRENA_NATIVE_THEME", None)
    from pyirena.gui import window_state

    window_state.GEOMETRY_PERSISTENCE_ENABLED = False
    from pyirena.gui._qt import QApplication
    from pyirena.gui.theme import apply_theme

    app = QApplication.instance() or QApplication(["make_doc_screenshots"])
    apply_theme(app)
    return app


def save_png(image, path: Path) -> int:
    """Save *image* as PNG, recompress with oxipng if installed; return bytes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(image, "save"):
        if not image.save(str(path), "PNG"):
            raise OSError(f"could not write {path}")
    oxipng = shutil.which("oxipng")
    if oxipng:
        import subprocess

        subprocess.run([oxipng, "-q", "-o", "2", "--strip", "safe", str(path)], check=False)
    return path.stat().st_size


def selected(ids: Iterable[str] | None) -> list[Shot]:
    """Registered shots, in id order, optionally restricted to *ids* (or page names)."""
    shots = [REGISTRY[k] for k in sorted(REGISTRY)]
    if not ids:
        return shots
    wanted = set(ids)
    unknown = wanted - set(REGISTRY) - {s.page for s in shots}
    if unknown:
        raise SystemExit(f"unknown shot id or page: {', '.join(sorted(unknown))}")
    return [s for s in shots if s.shot_id in wanted or s.page in wanted]
