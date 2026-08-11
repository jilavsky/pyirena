"""
Shared data-file type table and folder listing.

The "Type:" dropdown that filters a folder down to the files a tool can open
(HDF5 NXcanSAS, generic HDF5, text) existed as a verbatim copy in Data
Manipulation and Data Merge, alongside their own ``os.listdir`` loops.  Adding
a format — or fixing which extensions count as text — meant editing both and
missing one, the pattern that produced the filter and sort-key incidents.

Kept Qt-free (and in ``core`` next to :mod:`pyirena.core.file_sorting`) so
batch scripts can enumerate a folder exactly as the GUI shows it: the same
extensions, the same exclusions, and — combined with ``sort_names`` — the same
order.

Note that the Data Selector deliberately offers a *different* set of labels
("HDF5 Files", "Text Files", "All Supported Files") with its own combo
user-data keys, because it also lists mixed folders; it is not a consumer of
this table.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

__all__ = [
    "FILE_TYPES",
    "FILE_TYPE_EXTS",
    "DEFAULT_FILE_TYPE",
    "extensions_for",
    "files_in_folder",
    "files_with_extensions",
]

#: Dropdown entries, in display order.
FILE_TYPES: List[str] = [
    "HDF5 Nexus",
    "HDF5 Generic",
    "Text (.dat/.txt/.csv)",
]

#: Extensions accepted for each entry (lower-case, leading dot).
#:
#: "HDF5 Nexus" and "HDF5 Generic" share extensions on purpose — the choice
#: tells the *loader* how to read the file, not which files to list.
FILE_TYPE_EXTS: Dict[str, Tuple[str, ...]] = {
    "HDF5 Nexus":            (".h5", ".hdf5", ".hdf"),
    "HDF5 Generic":          (".h5", ".hdf5", ".hdf"),
    "Text (.dat/.txt/.csv)": (".dat", ".txt", ".csv"),
}

DEFAULT_FILE_TYPE = FILE_TYPES[0]


def extensions_for(file_type: str) -> Tuple[str, ...]:
    """Extensions for a dropdown entry; unknown entries list nothing.

    Returning an empty tuple rather than raising means a state file written by
    a build with an extra type leaves the list empty instead of crashing the
    panel on open.
    """
    return FILE_TYPE_EXTS.get(file_type, ())


def files_with_extensions(folder: str, extensions: Sequence[str]) -> List[str]:
    """List file *names* in *folder* whose extension is in *extensions*.

    The lower-level half of :func:`files_in_folder`, for a browser whose
    dropdown is not :data:`FILE_TYPES` — the Data Selector offers an "All
    supported files" entry that spans two of them.  Sharing this much still
    gives every browser the same answers to the awkward questions: skip
    sub-directories, match the extension case-insensitively, and treat an
    unreadable folder as empty rather than raising.

    Args:
        folder: Directory to list.  A missing or non-directory path yields [].
        extensions: Lower-case extensions with a leading dot.

    Returns:
        Unsorted file names — pass through
        :func:`pyirena.core.file_sorting.sort_names` for display order.
    """
    if not folder or not os.path.isdir(folder) or not extensions:
        return []
    wanted = {e.lower() for e in extensions}
    try:
        entries = os.listdir(folder)
    except OSError:
        return []
    return [
        name for name in entries
        if os.path.isfile(os.path.join(folder, name))
        and Path(name).suffix.lower() in wanted
    ]


def files_in_folder(folder: str, file_type: str = DEFAULT_FILE_TYPE) -> List[str]:
    """List the file *names* in *folder* matching *file_type*.

    Names, not paths, because the sort keys and the filename filter both match
    on the name.  Sub-directories are skipped.  An unreadable folder yields an
    empty list rather than raising: a browser pointed at a directory whose
    permissions changed should show nothing, not fall over.

    Args:
        folder: Directory to list.  A missing or non-directory path yields [].
        file_type: One of :data:`FILE_TYPES`.

    Returns:
        Unsorted file names — pass through
        :func:`pyirena.core.file_sorting.sort_names` for display order.
    """
    return files_with_extensions(folder, extensions_for(file_type))
