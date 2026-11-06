"""
Open HDF5 files so that two different files can never be confused.

HDF5 recognises an already-open file by the identity the filesystem reports
for it -- device and inode on POSIX, volume serial number and file index on
Windows -- not by its name.  Opening a file whose identity matches a file that
is already open in the process silently hands back the *already open* file.

Local disks give every file its own identity, so this never shows.  Network
filesystems need not: SSHFS-Win reports file index 0 for every file.  There,
any two HDF5 files open at the same time are one file as far as HDF5 is
concerned, and pyIrena reads -- or writes -- the wrong one without any error.
Seen at the beamline as Collect Values returning the transmission of the file
shown in the Data Explorer's browser for every selected file, and as
"HDF5 group not found" when Bernardyn opened a series in parallel.

:func:`open_h5` is a drop-in replacement for ``h5py.File(path, mode)``.  It
opens natively, exactly like ``h5py.File``, unless HDF5 would confuse the file
with a *different* file that is already open; then it opens through a Python
file object, which HDF5 cannot confuse with anything else.  With no other HDF5
file open -- the usual case -- it costs nothing.

All of pyIrena's HDF5 opens go through here.  ``pyirena/tests/test_h5open.py``
fails the build if a direct ``h5py.File(`` call appears in the package.
"""

from __future__ import annotations

import logging
import os
import threading
from pathlib import Path

import h5py

log = logging.getLogger(__name__)

# The check for a confusable open file and the open itself must happen as one
# step, or two threads can both decide "nothing else is open" and open natively.
_OPEN_LOCK = threading.RLock()


class _FileObjectH5(h5py.File):
    """An ``h5py.File`` opened on a Python file object that it owns.

    Closing it closes the Python file too, so no OS handle is left behind for
    the garbage collector (on Windows an open handle blocks other programs
    from replacing the file).
    """

    def __init__(self, fileobj, mode: str, path: str) -> None:
        super().__init__(fileobj, mode)
        self._pyirena_fileobj = fileobj
        self._pyirena_path = path

    @property
    def filename(self) -> str:
        """The path on disk, as for a natively opened file."""
        return self._pyirena_path

    def close(self) -> None:
        try:
            super().close()
        finally:
            self._pyirena_fileobj.close()


def _identity(path: str) -> tuple[int, int]:
    """(device, file id) as the filesystem reports them -- what HDF5 compares.

    On Windows these are the volume serial number and the file index, the
    same fields HDF5's Windows file drivers use.  Raises OSError if *path*
    does not exist.
    """
    st = os.stat(path)
    return st.st_dev, st.st_ino


def _same_path(a: str, b: str) -> bool:
    return os.path.normcase(os.path.abspath(a)) == os.path.normcase(os.path.abspath(b))


def _open_native_files() -> list[str]:
    """Paths of the HDF5 files open natively in this process.

    Files opened through a Python file object are named after the object
    (``"<_io.BufferedReader name=...>"``); HDF5 never confuses those with a
    file on disk, so they are left out.
    """
    names = []
    for fid in h5py.h5f.get_obj_ids(h5py.h5f.OBJ_ALL, h5py.h5f.OBJ_FILE):
        try:
            name = os.fsdecode(fid.name)
        except Exception:
            continue
        if not name.startswith("<"):
            names.append(name)
    return names


def confusable_open_file(path: str | Path) -> str | None:
    """Return an open HDF5 file that HDF5 would confuse with *path*, or None.

    That is an open file with a different path but the same filesystem
    identity.  For a file that does not exist yet (about to be created), any
    open file on the same device whose identity is 0 counts: such a device
    gives the new file identity 0 as well.
    """
    path = os.fspath(path)
    others = [n for n in _open_native_files() if not _same_path(n, path)]
    if not others:
        return None
    try:
        device, file_id = _identity(path)
        creating = False
    except OSError:
        try:
            device, _ = _identity(os.path.dirname(os.path.abspath(path)))
        except OSError:
            return None
        file_id = 0
        creating = True
    for other in others:
        try:
            other_device, other_id = _identity(other)
        except OSError:
            continue
        if other_device != device:
            continue
        if other_id == file_id or (creating and other_id == 0):
            return other
    return None


def _open_through_file_object(path: str, mode: str) -> _FileObjectH5:
    # h5py mode -> Python file mode, and the h5py mode to use on that object.
    if mode == "r":
        fileobj, h5_mode = open(path, "rb"), "r"
    elif mode == "r+" or (mode == "a" and os.path.exists(path)):
        fileobj, h5_mode = open(path, "r+b"), "r+"
    elif mode in ("w", "a"):
        fileobj, h5_mode = open(path, "w+b"), "w"
    elif mode in ("w-", "x"):
        fileobj, h5_mode = open(path, "x+b"), "w"
    else:
        raise ValueError(f"Invalid HDF5 file mode {mode!r}")
    try:
        return _FileObjectH5(fileobj, h5_mode, path)
    except Exception:
        fileobj.close()
        raise


def open_h5(path: str | Path, mode: str = "r") -> h5py.File:
    """Open an HDF5 file; a drop-in replacement for ``h5py.File(path, mode)``.

    Opens natively unless HDF5 would confuse *path* with a different file that
    is already open (see the module docstring), in which case the file is
    opened through a Python file object instead.  Either way the result is an
    ``h5py.File`` with the same behaviour, usable as a context manager.

    Args:
        path: File on disk.
        mode: Any ``h5py.File`` mode: ``"r"``, ``"r+"``, ``"w"``, ``"w-"``,
            ``"x"``, ``"a"``.

    Returns:
        The open ``h5py.File``.  Close it, or use it in a ``with`` block.
    """
    path = os.fspath(path)
    with _OPEN_LOCK:
        other = confusable_open_file(path)
        if other is None:
            return h5py.File(path, mode)
        log.debug(
            "HDF5 cannot tell %s from open file %s (same filesystem identity); "
            "opening through a Python file object", path, other,
        )
        return _open_through_file_object(path, mode)
