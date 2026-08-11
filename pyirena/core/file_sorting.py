"""
Shared filename sort keys for every file browser in pyIrena.

Beamline filenames carry their own metadata — ``AeroGel_500C_10min_03.h5`` says
temperature, time and run order without opening the file — so the browsers sort
on patterns pulled out of the name.  That logic existed **four times**: in the
Data Selector, the HDF5 Viewer (Data Explorer), Data Manipulation, and a lone
order-number copy in Data Merge.  A new sort mode, or a fix to one of the
regexes, had to be applied in all four and would silently be missed in some —
the same failure mode as the filename filter before
:mod:`pyirena.gui.file_filter`.

This module is the single implementation.  It is deliberately Qt-free so batch
and api code can order a list of input files exactly as the GUI displays them:
a batch run over a temperature series produces results in the same order the
user saw on screen.

Recognised patterns (case-insensitive, matched on the file name only)

======================  ==========================  ======================
Sort mode               Pattern                     Example
======================  ==========================  ======================
Temperature             ``_<number>C``              ``sample_500C.h5``
Time                    ``_<number>min``            ``sample_10min.h5``
Pressure                ``_<number>PSI``            ``sample_100PSI.h5``
Order number            trailing bare ``_<int>``    ``sample_10min_03.h5``
======================  ==========================  ======================

Files whose name lacks the pattern sort **last** (the key is ``+inf``), so a
stray file never lands in the middle of an ordered series.
"""

from __future__ import annotations

import re
from typing import Callable, Iterable, List, Union

__all__ = [
    "SORT_LABELS",
    "SORT_TOOLTIP",
    "SORT_KEYS",
    "DEFAULT_SORT_INDEX",
    "sort_key_name",
    "sort_key_temperature",
    "sort_key_time",
    "sort_key_order",
    "sort_key_pressure",
    "sort_key_for_index",
    "is_descending",
    "sort_names",
]

SortKey = Union[str, float]


# ── Key extractors ────────────────────────────────────────────────────────

def sort_key_name(name: str) -> str:
    """Case-insensitive filename key."""
    return name.lower()


def sort_key_temperature(name: str) -> float:
    """Temperature from ``_500C`` / ``_-20C``; +inf when absent."""
    m = re.search(r'_(-?\d+(?:\.\d+)?)C(?=_|\.|$)', name, re.IGNORECASE)
    return float(m.group(1)) if m else float('inf')


def sort_key_time(name: str) -> float:
    """Elapsed time from ``_10min``; +inf when absent."""
    m = re.search(r'_(\d+(?:\.\d+)?)min(?=_|\.|$)', name, re.IGNORECASE)
    return float(m.group(1)) if m else float('inf')


def sort_key_order(name: str) -> float:
    """Run/order number: the last bare integer segment of the stem.

    Strips the extension then scans ``_``-separated segments right to left for
    a digits-only part.  Any segment containing letters is skipped, so
    ``_merged``, ``_mrg``, ``_scaled`` and unit-bearing tokens such as
    ``_10min`` or ``_5C`` do not masquerade as the run number:
    ``sample_500C_10min_03_merged.h5`` sorts as 3.
    """
    stem = re.sub(r'\.[^.]+$', '', name)
    for part in reversed(stem.split('_')):
        if re.fullmatch(r'\d+', part):
            return float(part)
    return float('inf')


def sort_key_pressure(name: str) -> float:
    """Pressure from ``_100PSI``; +inf when absent."""
    m = re.search(r'_(\d+(?:\.\d+)?)PSI(?=_|\.|$)', name, re.IGNORECASE)
    return float(m.group(1)) if m else float('inf')


# ── The mode list shared by every browser's dropdown ──────────────────────
#
# Ascending/descending alternate, so index >> 1 selects the quantity and the
# low bit selects the direction.  Keep the two lists in step.

SORT_LABELS: List[str] = [
    "Filename A→Z",
    "Filename Z→A",
    "Temperature ↑",
    "Temperature ↓",
    "Time ↑",
    "Time ↓",
    "Order number ↑",
    "Order number ↓",
    "Pressure ↑",
    "Pressure ↓",
]

SORT_KEYS: List[Callable[[str], SortKey]] = [
    sort_key_name,         # 0 Filename A→Z
    sort_key_name,         # 1 Filename Z→A
    sort_key_temperature,  # 2 Temperature ↑
    sort_key_temperature,  # 3 Temperature ↓
    sort_key_time,         # 4 Time ↑
    sort_key_time,         # 5 Time ↓
    sort_key_order,        # 6 Order number ↑
    sort_key_order,        # 7 Order number ↓
    sort_key_pressure,     # 8 Pressure ↑
    sort_key_pressure,     # 9 Pressure ↓
]

#: Order number ↑ — what a measurement series is usually in already.
DEFAULT_SORT_INDEX = 6

SORT_TOOLTIP = (
    "Sort order for the file list.\n"
    "Patterns recognised in filenames:\n"
    "  Temperature  : _25C      (e.g. sample_500C.h5)\n"
    "  Time         : _50min    (e.g. sample_10min.h5)\n"
    "  Pressure     : _100PSI   (e.g. sample_100PSI.h5)\n"
    "  Order number : trailing _12  (e.g. sample_10min_03.h5)\n"
    "Matching is case-insensitive.  Files without the pattern sort last."
)


# ── Applying a mode ───────────────────────────────────────────────────────

def sort_key_for_index(index: int) -> Callable[[str], SortKey]:
    """Key function for a dropdown index, clamped to the valid range."""
    if index < 0:
        index = 0
    return SORT_KEYS[min(index, len(SORT_KEYS) - 1)]


def is_descending(index: int) -> bool:
    """True when the dropdown index selects the descending direction."""
    return bool(index % 2)


def sort_names(names: Iterable[str], index: int = DEFAULT_SORT_INDEX) -> List[str]:
    """Sort file names by the mode at dropdown *index*.

    Files whose name lacks the pattern sort last in **both** directions.  A
    plain ``reverse=True`` would flip them to the top, so a descending
    temperature list used to open with every unmatched file — logs, notes, a
    stray average — above the hottest measurement.  Among themselves they keep
    the order they arrived in.

    Args:
        names: File names (not paths — the patterns are matched on the name).
        index: Position in :data:`SORT_LABELS`.

    Returns:
        A new sorted list; the input is not modified.
    """
    key_fn = sort_key_for_index(index)
    descending = is_descending(index)

    if key_fn is sort_key_name:          # plain text sort, nothing can be missing
        return sorted(names, key=key_fn, reverse=descending)

    def ordering(name: str):
        value = key_fn(name)
        if value == float('inf'):        # pattern not present in this filename
            return (1, 0.0)
        # Sorting ascending on a negated value gives descending order while
        # leaving the "missing" flag above untouched.
        return (0, -value if descending else value)

    return sorted(names, key=ordering)
