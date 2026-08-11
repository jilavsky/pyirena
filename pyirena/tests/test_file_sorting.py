"""
Tests for :mod:`pyirena.core.file_sorting` — the shared filename sort keys.

No Qt needed: the module is pure regex, which is why it can live in ``core``
and be reused by batch and api as well as by the four GUI file browsers.

Beamline users rely on these patterns daily (``AeroGel_500C_10min_03.h5``), so
the cases pinned here are the ones that would quietly reorder someone's series:
negative temperatures, decimals, suffixed filenames, and files that carry no
pattern at all.
"""

import pytest

from pyirena.core.file_sorting import (
    DEFAULT_SORT_INDEX,
    SORT_KEYS,
    SORT_LABELS,
    SORT_TOOLTIP,
    is_descending,
    sort_key_for_index,
    sort_key_name,
    sort_key_order,
    sort_key_pressure,
    sort_key_temperature,
    sort_key_time,
    sort_names,
)

INF = float("inf")


# ── Individual keys ──────────────────────────────────────────────────────

@pytest.mark.parametrize("name,expected", [
    ("sample_500C.h5", 500.0),
    ("sample_500c.h5", 500.0),            # case-insensitive
    ("sample_-20C.h5", -20.0),            # cryo runs
    ("sample_25.5C.h5", 25.5),
    ("sample_500C_10min.h5", 500.0),      # followed by another segment
    ("sample.h5", INF),                   # no pattern → last
    ("sample_500Celsius.h5", INF),        # C must end the segment
    ("500C.h5", INF),                     # needs the leading underscore
])
def test_temperature_key(name, expected):
    assert sort_key_temperature(name) == expected


@pytest.mark.parametrize("name,expected", [
    ("sample_10min.h5", 10.0),
    ("sample_0.5min.h5", 0.5),
    ("sample_10MIN.h5", 10.0),
    ("sample_10minutes.h5", INF),
    ("sample.h5", INF),
])
def test_time_key(name, expected):
    assert sort_key_time(name) == expected


@pytest.mark.parametrize("name,expected", [
    ("sample_100PSI.h5", 100.0),
    ("sample_12.5psi.h5", 12.5),
    ("sample.h5", INF),
])
def test_pressure_key(name, expected):
    assert sort_key_pressure(name) == expected


@pytest.mark.parametrize("name,expected", [
    ("sample_03.h5", 3.0),
    ("sample_500C_10min_03.h5", 3.0),          # last bare integer wins
    ("sample_500C_10min_03_merged.h5", 3.0),   # word suffixes are skipped…
    ("sample_03_mrg.h5", 3.0),
    ("sample_03_scaled.h5", 3.0),
    ("sample_10min.h5", INF),                  # …and unit tokens are not numbers
    ("sample_5C.h5", INF),
    ("sample.h5", INF),
])
def test_order_key(name, expected):
    assert sort_key_order(name) == expected


def test_name_key_is_case_insensitive():
    assert sort_key_name("Sample.H5") == "sample.h5"


# ── Mode table ───────────────────────────────────────────────────────────

def test_labels_and_keys_stay_in_step():
    assert len(SORT_LABELS) == len(SORT_KEYS)
    # Ascending/descending alternate, so the low bit is the direction.
    for i, label in enumerate(SORT_LABELS):
        if i % 2:
            assert label.endswith("↓") or "Z→A" in label, label
        else:
            assert label.endswith("↑") or "A→Z" in label, label


def test_direction_comes_from_the_index_parity():
    assert not is_descending(0)
    assert is_descending(1)
    assert not is_descending(DEFAULT_SORT_INDEX)


def test_index_is_clamped_rather_than_raising():
    """A state file from a build with more modes must not crash the browser."""
    assert sort_key_for_index(99) is SORT_KEYS[-1]
    assert sort_key_for_index(-1) is SORT_KEYS[0]


def test_default_is_order_number_ascending():
    assert SORT_LABELS[DEFAULT_SORT_INDEX] == "Order number ↑"


def test_tooltip_documents_every_pattern():
    for token in ("_25C", "_50min", "_100PSI", "_12"):
        assert token in SORT_TOOLTIP


# ── Sorting whole lists ──────────────────────────────────────────────────

SERIES = [
    "s_500C_10min_03.h5",
    "s_100C_05min_01.h5",
    "s_300C_20min_02.h5",
    "notes.h5",              # carries no pattern at all
]


def test_sort_by_temperature_ascending():
    assert sort_names(SERIES, 2) == [
        "s_100C_05min_01.h5", "s_300C_20min_02.h5", "s_500C_10min_03.h5",
        "notes.h5",
    ]


def test_files_without_the_pattern_sort_last_in_both_directions():
    """Regression: descending used to put unmatched files at the top.

    A plain reverse=True flips the +inf sentinel to the front, so a descending
    temperature list opened with every log and stray file above the hottest
    measurement.
    """
    ascending = sort_names(SERIES, 2)
    descending = sort_names(SERIES, 3)
    assert ascending[-1] == "notes.h5"
    assert descending[-1] == "notes.h5"
    assert descending[:3] == [
        "s_500C_10min_03.h5", "s_300C_20min_02.h5", "s_100C_05min_01.h5",
    ]


def test_descending_reverses_only_the_matched_files():
    assert sort_names(SERIES, 5) == [           # Time ↓
        "s_300C_20min_02.h5", "s_500C_10min_03.h5", "s_100C_05min_01.h5",
        "notes.h5",
    ]


def test_order_number_is_the_default_view():
    assert sort_names(SERIES) == [
        "s_100C_05min_01.h5", "s_300C_20min_02.h5", "s_500C_10min_03.h5",
        "notes.h5",
    ]


def test_unmatched_files_keep_their_relative_order():
    names = ["zeta.h5", "alpha.h5", "s_100C.h5"]
    assert sort_names(names, 2) == ["s_100C.h5", "zeta.h5", "alpha.h5"]


def test_sort_names_does_not_modify_the_input():
    names = list(SERIES)
    sort_names(names, 3)
    assert names == SERIES


def test_filename_sort_still_reverses_wholesale():
    assert sort_names(["b.h5", "a.h5", "c.h5"], 1) == ["c.h5", "b.h5", "a.h5"]
