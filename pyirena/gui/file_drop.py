"""
Drop a file onto a pyIrena window to open it.

Igor, Origin and every modern analysis GUI accept a file dragged onto the
window; pyIrena accepted none, so opening data always meant Select Folder →
find it in the list.  This module is the one implementation, in the same spirit
as :mod:`pyirena.gui.file_filter` and :mod:`pyirena.gui.table_utils`: a panel
calls one function and gets the whole behaviour, including the parts that are
easy to get wrong.

What it handles so panels do not have to:

* **Accepting only what the panel can open.**  The drag is rejected (the cursor
  shows "no") while it is still over the window if none of the dragged files
  have a usable extension, so the user finds out before letting go.
* **Folders.**  Dropping a folder yields the matching files inside it, which is
  what a user dragging a run directory means.
* **Both drop payloads.**  Finder and Explorer send file *URLs*; some tools
  send plain text paths.  Both are read.
* **Not subclassing.**  The behaviour is installed through an event filter, so
  it works on an existing widget without changing its class — several pyIrena
  panels build their file lists inline.

A panel wires it in one call::

    from pyirena.gui.file_drop import enable_file_drop

    enable_file_drop(self.file_list, self._on_files_dropped)

and receives a list of absolute paths, already filtered and sorted.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Callable, Iterable, List, Optional, Sequence
from urllib.request import url2pathname

log = logging.getLogger(__name__)

__all__ = [
    "DATA_EXTENSIONS",
    "collect_dropped_paths",
    "enable_file_drop",
    "paths_from_mime",
    "select_dropped_in_list",
]

#: Extensions pyIrena can open, as a sensible default for any drop target.
#: Panels that read only one kind pass their own list.
DATA_EXTENSIONS: tuple = (
    ".h5", ".hdf5", ".hdf", ".nxs", ".nx",      # NXcanSAS / HDF5
    ".dat", ".txt", ".csv",                      # text data
)

#: How deep to look when a folder is dropped.  One level: dropping a run
#: directory should offer its files, not everything in a tree of sub-runs.
FOLDER_SCAN_DEPTH = 1

#: Sort mode used for the returned list — index 0 in
#: :data:`pyirena.core.file_sorting.SORT_LABELS` ("Filename A→Z").
NAME_SORT_INDEX = 0


# ──────────────────────────────────────────────────────────────────────────
#  Path collection (pure — no Qt, so it can be tested headlessly)
# ──────────────────────────────────────────────────────────────────────────

def collect_dropped_paths(
    raw_paths: Iterable[str],
    extensions: Sequence[str] = DATA_EXTENSIONS,
    folder_depth: int = FOLDER_SCAN_DEPTH,
) -> List[str]:
    """Turn what was dropped into a list of openable files.

    Args:
        raw_paths: Filesystem paths from the drop (files and/or folders).
        extensions: Acceptable extensions, lower-case with a leading dot.
            Pass an empty sequence to accept anything.
        folder_depth: How many directory levels to scan for a dropped folder.

    Returns:
        Absolute paths of existing, matching files — de-duplicated, and sorted
        so a dropped folder arrives in the order the browsers would show it.
    """
    accepted = {e.lower() for e in extensions}
    found: List[str] = []

    for raw in raw_paths:
        if not raw:
            continue
        path = Path(raw).expanduser()
        try:
            if path.is_dir():
                found.extend(_scan_folder(path, accepted, folder_depth))
            elif path.is_file() and _matches(path, accepted):
                found.append(str(path.resolve()))
        except OSError:
            log.debug("could not inspect dropped path %s", raw, exc_info=True)

    # De-duplicate while keeping first-seen order, then sort by file name the
    # way the file browsers do, so a dropped folder is not in random order.
    # Filename A→Z (index 0) on purpose rather than the browsers' default sort
    # mode: the panel receiving these re-sorts them with its own setting, and a
    # mode that finds nothing to sort on (no temperature, no order number in the
    # names) would leave them in whatever order the file manager happened to send.
    unique = list(dict.fromkeys(found))
    try:
        from pyirena.core.file_sorting import sort_names  # noqa: PLC0415

        by_name = {os.path.basename(p): p for p in unique}
        if len(by_name) == len(unique):
            return [by_name[n] for n in sort_names(list(by_name), NAME_SORT_INDEX)]
    except Exception:
        log.debug("could not sort dropped paths", exc_info=True)
    return unique


def _matches(path: Path, accepted: set) -> bool:
    return not accepted or path.suffix.lower() in accepted


def _scan_folder(folder: Path, accepted: set, depth: int) -> List[str]:
    """Files inside a dropped folder, at most *depth* levels down."""
    out: List[str] = []
    try:
        for entry in sorted(folder.iterdir()):
            if entry.is_file() and _matches(entry, accepted):
                out.append(str(entry.resolve()))
            elif entry.is_dir() and depth > 1:
                out.extend(_scan_folder(entry, accepted, depth - 1))
    except OSError:
        log.debug("could not list dropped folder %s", folder, exc_info=True)
    return out


def _path_from_text_line(line: str) -> str:
    """Turn one line of a text drop payload into a native filesystem path.

    Qt hands back ``C:/Users/...`` on Windows and applications send several
    spellings of a file URL, so both are normalised here rather than in every
    caller.  A path that reaches a panel with the wrong separator compares
    unequal to the same file discovered by a folder scan, which is how a
    dropped file ends up loaded twice.
    """
    if not line:
        return ""
    if line.startswith("file://"):
        rest = line[7:]
        if rest.startswith("/"):
            # A well-formed file URL carries an empty authority, so the
            # remainder is ``/tmp/x`` on POSIX and ``/C:/Users/x`` on Windows.
            # url2pathname knows the difference (and un-escapes ``%20``).
            line = url2pathname(rest)
        else:
            # Some applications drop the third slash and send ``file://C:\x``.
            # That is not a URL; take it literally rather than mangling it.
            line = rest
    if not line:
        return ""
    return os.path.normpath(line)


def paths_from_mime(mime) -> List[str]:
    """Extract filesystem paths from a Qt drop payload.

    File managers send URLs; some applications send newline-separated text.
    Both are read, because which one arrives depends on the source app rather
    than on anything pyIrena controls.
    """
    paths: List[str] = []
    try:
        if mime.hasUrls():
            for url in mime.urls():
                local = url.toLocalFile()
                if local:
                    paths.append(os.path.normpath(local))
        if not paths and mime.hasText():
            for line in mime.text().splitlines():
                local = _path_from_text_line(line.strip())
                if local:
                    paths.append(local)
    except Exception:
        log.debug("could not read drop payload", exc_info=True)
    return paths


# ──────────────────────────────────────────────────────────────────────────
#  Qt wiring
# ──────────────────────────────────────────────────────────────────────────

def enable_file_drop(
    widget,
    on_files: Callable[[List[str]], None],
    extensions: Sequence[str] = DATA_EXTENSIONS,
    folder_depth: int = FOLDER_SCAN_DEPTH,
    accept_folders: bool = True,
):
    """Let *widget* accept dropped data files, calling *on_files* with them.

    Args:
        widget: Any ``QWidget`` — a file list, a graph, or a whole panel.
        on_files: Called with a list of absolute paths once the user drops.
            Never called with an empty list.
        extensions: Extensions this target can open.
        folder_depth: Directory levels scanned when a folder is dropped.
        accept_folders: Whether dropping a folder is meaningful here.

    Returns:
        The installed event filter (kept alive on the widget), or None if the
        wiring failed — a panel should never break because drag-and-drop could
        not be set up.
    """
    try:
        from pyirena.gui._qt import QEvent, QtCore  # noqa: PLC0415
    except Exception:
        log.debug("Qt unavailable; drag-and-drop not installed", exc_info=True)
        return None

    class _FileDropFilter(QtCore.QObject):
        """Accepts the drag only when it carries something we can open."""

        def eventFilter(self, obj, event):
            kind = event.type()
            if kind in (QEvent.Type.DragEnter, QEvent.Type.DragMove):
                if self._wanted(event):
                    event.acceptProposedAction()
                    return True
                event.ignore()
                return True
            if kind == QEvent.Type.Drop:
                paths = self._paths(event)
                if not paths:
                    event.ignore()
                    return True
                event.acceptProposedAction()
                try:
                    on_files(paths)
                except Exception:
                    log.exception("handler for dropped files failed")
                return True
            return False

        def _paths(self, event) -> List[str]:
            raw = paths_from_mime(event.mimeData())
            if not accept_folders:
                raw = [p for p in raw if not Path(p).is_dir()]
            return collect_dropped_paths(raw, extensions, folder_depth)

        def _wanted(self, event) -> bool:
            """Cheap check during the drag — is anything here openable?

            Deliberately does not scan folders: that runs on every mouse move
            over the window.  A folder is assumed to be worth accepting and is
            expanded on drop, where an empty result is simply ignored.
            """
            accepted = {e.lower() for e in extensions}
            for raw in paths_from_mime(event.mimeData()):
                path = Path(raw)
                if accept_folders and path.is_dir():
                    return True
                if not accepted or path.suffix.lower() in accepted:
                    return True
            return False

    try:
        widget.setAcceptDrops(True)
        drop_filter = _FileDropFilter(widget)
        widget.installEventFilter(drop_filter)
        # Hold a reference so the filter lives as long as the widget.
        widget._pyirena_file_drop_filter = drop_filter
        return drop_filter
    except Exception:
        log.debug("could not enable drag-and-drop", exc_info=True)
        return None


def drop_hint(target: str = "files") -> str:
    """Tooltip line advertising the feature, so it is discoverable."""
    return f"You can also drag {target} here from Finder/Explorer."


def select_dropped_in_list(list_widget, paths: Sequence[str], folder: str) -> int:
    """Select the dropped files in a ``QListWidget`` and scroll to the first.

    The second half of every drop handler: the browser has switched folder, and
    now the files the user actually dragged should be the selection.  Only
    paths from *folder* are considered — several folders can be dropped at once
    and a file list shows one.

    Args:
        list_widget: The ``QListWidget`` showing file *names*.
        paths: Absolute paths from the drop.
        folder: The folder the browser is now showing.

    Returns:
        How many items were selected — zero means the drop landed on files the
        current type/filter settings hide, which is worth telling the user
        rather than leaving the drop looking like it did nothing.
    """
    wanted = {os.path.basename(p) for p in paths if os.path.dirname(p) == folder}
    if not wanted:
        return 0
    try:
        list_widget.clearSelection()
        matched = 0
        for row in range(list_widget.count()):
            item = list_widget.item(row)
            if item.text() in wanted:
                item.setSelected(True)
                matched += 1
                if matched == 1:
                    list_widget.scrollToItem(item)
        return matched
    except Exception:
        log.debug("could not select dropped files", exc_info=True)
        return 0


def first_folder(paths: Sequence[str]) -> Optional[str]:
    """Folder of the first dropped path — what a browser should switch to."""
    for path in paths:
        parent = Path(path).parent
        if parent.is_dir():
            return str(parent)
    return None
