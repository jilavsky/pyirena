"""
Shared filename-filter logic for every file browser in the GUI.

All of pyIrena's file-listing widgets (Data Selector, Data Manipulation, Data
Merge, HDF5 Viewer) present a "Filter" box.  Historically only the HDF5 Viewer
interpreted that text as a regular expression; the others did a plain
case-insensitive substring test, which surprised users who expected the
grep-like behaviour the documentation promised and that Igor Irena provided.

This module is the single implementation.  It is deliberately Qt-free and
matplotlib-free so it can be unit-tested headlessly and reused from any layer.

Semantics
---------
- The filter text is treated as a Python regular expression matched with
  ``re.search`` (i.e. it matches anywhere in the name, like grep) and
  ``re.IGNORECASE``.
- A plain fragment such as ``60C`` is a valid regex, so simple substring
  filtering keeps working exactly as before.
- If the text is *not* a valid regex — most commonly a half-typed pattern such
  as ``sample(`` while the user is still typing — the filter silently falls
  back to a case-insensitive substring test rather than hiding everything or
  raising.
- Empty/whitespace-only text matches everything.
"""

from __future__ import annotations

import re
from typing import Callable, Iterable, List

__all__ = ["make_file_matcher", "filter_names", "is_valid_filter", "FILTER_TOOLTIP",
           "FILTER_PLACEHOLDER"]


# Shared UI strings so every panel advertises the same behaviour.
FILTER_PLACEHOLDER = "Filter… (regex OK, e.g. 60C|0[12]min)"

FILTER_TOOLTIP = (
    "Filter filenames.  Regular expressions are supported:\n"
    "  60C          — names containing '60C'\n"
    "  0[12]min     — names with '01min' or '02min'\n"
    "  60C|100C     — names with '60C' or '100C'\n"
    "  ^sample      — names starting with 'sample'\n"
    "  \\.h5$        — names ending with '.h5'\n"
    "  ^(?!.*bkg)   — names NOT containing 'bkg'\n"
    "Matching is case-insensitive.  Plain text fragments also work; an\n"
    "incomplete or invalid pattern falls back to a substring match."
)


def make_file_matcher(filter_text: str) -> Callable[[str], bool]:
    """Build a predicate that decides whether a filename passes the filter.

    Args:
        filter_text: Raw contents of a Filter box.  May be empty, a plain
            substring, or any Python regular expression.

    Returns:
        A callable taking a filename and returning True if it should be shown.
        Compilation happens once here, not per filename, so the returned
        predicate is cheap to apply across a large listing.
    """
    if filter_text is None:
        return lambda name: True

    text = filter_text.strip()
    if not text:
        return lambda name: True

    try:
        pattern = re.compile(text, re.IGNORECASE)
    except re.error:
        # Half-typed or malformed pattern — degrade to substring matching so
        # the list stays useful while the user is still typing.
        lowered = text.lower()
        return lambda name: lowered in name.lower()

    return lambda name: bool(pattern.search(name))


def filter_names(names: Iterable[str], filter_text: str) -> List[str]:
    """Return the subset of ``names`` matching ``filter_text``.

    Convenience wrapper around :func:`make_file_matcher` for the panels that
    rebuild their list widget from scratch rather than hiding rows.
    """
    matches = make_file_matcher(filter_text)
    return [n for n in names if matches(n)]


def is_valid_filter(filter_text: str) -> bool:
    """Report whether ``filter_text`` compiles as a regular expression.

    Empty text counts as valid.  Useful if a panel wants to give the user
    visual feedback about a malformed pattern.
    """
    if not filter_text or not filter_text.strip():
        return True
    try:
        re.compile(filter_text.strip())
    except re.error:
        return False
    return True
