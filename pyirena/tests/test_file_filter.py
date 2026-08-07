"""Tests for the shared filename-filter helper (pyirena.gui.file_filter).

The module is deliberately Qt-free, so these run headless without importing
PySide6.
"""

import pytest

from pyirena.gui.file_filter import (
    filter_names,
    is_valid_filter,
    make_file_matcher,
)

NAMES = [
    "sample_60C_01min.h5",
    "sample_60C_02min.h5",
    "sample_100C_01min.h5",
    "bkg_60C_01min.h5",
    "Sample_25C_10min.dat",
    "spheres_Rg50.h5",
]


# ── Empty / trivial filters ────────────────────────────────────────────────

@pytest.mark.parametrize("text", ["", "   ", None])
def test_empty_filter_matches_everything(text):
    assert filter_names(NAMES, text) == NAMES


def test_plain_substring_still_works():
    assert filter_names(NAMES, "Rg50") == ["spheres_Rg50.h5"]


def test_matching_is_case_insensitive():
    assert filter_names(NAMES, "sample") == [
        "sample_60C_01min.h5",
        "sample_60C_02min.h5",
        "sample_100C_01min.h5",
        "Sample_25C_10min.dat",
    ]


# ── Real regex features — the point of the change ──────────────────────────

def test_alternation():
    got = filter_names(NAMES, "100C|25C")
    assert got == ["sample_100C_01min.h5", "Sample_25C_10min.dat"]


def test_character_class():
    got = filter_names(NAMES, "0[12]min")
    assert "sample_60C_01min.h5" in got
    assert "sample_60C_02min.h5" in got
    assert "Sample_25C_10min.dat" not in got


def test_anchor_start():
    got = filter_names(NAMES, "^sample")
    # Case-insensitive, so the capitalised variant is included too.
    assert "bkg_60C_01min.h5" not in got
    assert "Sample_25C_10min.dat" in got


def test_anchor_end():
    got = filter_names(NAMES, r"\.h5$")
    assert "Sample_25C_10min.dat" not in got
    assert len(got) == 5


def test_negative_lookahead_excludes():
    got = filter_names(NAMES, "^(?!.*bkg)")
    assert "bkg_60C_01min.h5" not in got
    assert len(got) == len(NAMES) - 1


def test_search_not_fullmatch():
    """Patterns match anywhere in the name, like grep — not just the whole name."""
    assert filter_names(NAMES, "60C") == [
        "sample_60C_01min.h5",
        "sample_60C_02min.h5",
        "bkg_60C_01min.h5",
    ]


def test_quantifier_and_group():
    got = filter_names(NAMES, r"_(\d+)C_0\dmin")
    assert "Sample_25C_10min.dat" not in got
    assert "sample_100C_01min.h5" in got


# ── Invalid patterns degrade gracefully ────────────────────────────────────

def test_invalid_regex_falls_back_to_substring():
    # A half-typed group: must not raise, and must not blank the list.
    got = filter_names(["a(b.h5", "plain.h5"], "a(b")
    assert got == ["a(b.h5"]


def test_invalid_regex_unbalanced_bracket():
    got = filter_names(["x[1].h5", "y.h5"], "x[1")
    assert got == ["x[1].h5"]


def test_make_file_matcher_returns_callable():
    m = make_file_matcher("60C")
    assert callable(m)
    assert m("sample_60C_01min.h5") is True
    assert m("sample_25C_01min.h5") is False


# ── Validity reporting ─────────────────────────────────────────────────────

@pytest.mark.parametrize("text", ["", "  ", "60C", "0[12]min", "^(?!.*bkg)", r"\.h5$"])
def test_is_valid_filter_true(text):
    assert is_valid_filter(text) is True


@pytest.mark.parametrize("text", ["a(b", "x[1", "*bad", "(?P<"])
def test_is_valid_filter_false(text):
    assert is_valid_filter(text) is False
