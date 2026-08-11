"""
Tests for :mod:`pyirena.core.file_types` — the shared data-file type table.

No Qt needed, which is the point: batch scripts enumerate a folder with the
same extensions the GUI browsers show, so a batch run cannot silently process a
different set of files than the user saw.

The behaviours pinned here are the ones that would otherwise surface as "the
browser is empty" or "the panel crashed on open".
"""

import pytest

from pyirena.core.file_types import (
    DEFAULT_FILE_TYPE,
    FILE_TYPE_EXTS,
    FILE_TYPES,
    extensions_for,
    files_in_folder,
)


@pytest.fixture
def folder(tmp_path):
    for name in ("a_01.h5", "b_02.hdf5", "c_03.hdf",
                 "d_04.dat", "e_05.txt", "f_06.csv",
                 "g_07.pxp", "notes"):
        (tmp_path / name).touch()
    (tmp_path / "subfolder").mkdir()
    (tmp_path / "subfolder" / "buried.h5").touch()
    return tmp_path


# ── The table ────────────────────────────────────────────────────────────

def test_every_dropdown_entry_has_extensions():
    for file_type in FILE_TYPES:
        assert FILE_TYPE_EXTS[file_type], file_type


def test_default_is_the_first_entry():
    assert DEFAULT_FILE_TYPE == FILE_TYPES[0] == "HDF5 Nexus"


def test_extensions_are_lower_case_with_a_leading_dot():
    for exts in FILE_TYPE_EXTS.values():
        for ext in exts:
            assert ext.startswith(".") and ext == ext.lower(), ext


def test_unknown_type_lists_nothing_instead_of_raising():
    """A state file from a build with an extra type must not crash a panel."""
    assert extensions_for("Some Future Format") == ()


# ── Listing ──────────────────────────────────────────────────────────────

def test_hdf5_types_list_all_hdf5_extensions(folder):
    assert sorted(files_in_folder(str(folder), "HDF5 Nexus")) == [
        "a_01.h5", "b_02.hdf5", "c_03.hdf",
    ]
    # Generic reads the same files differently — same listing.
    assert (sorted(files_in_folder(str(folder), "HDF5 Generic"))
            == sorted(files_in_folder(str(folder), "HDF5 Nexus")))


def test_text_type_lists_text_extensions(folder):
    assert sorted(files_in_folder(str(folder), "Text (.dat/.txt/.csv)")) == [
        "d_04.dat", "e_05.txt", "f_06.csv",
    ]


def test_unrelated_and_extensionless_files_are_ignored(folder):
    listed = files_in_folder(str(folder), "HDF5 Nexus")
    assert "g_07.pxp" not in listed
    assert "notes" not in listed


def test_subfolders_are_not_listed(folder):
    listed = files_in_folder(str(folder), "HDF5 Nexus")
    assert "subfolder" not in listed
    assert "buried.h5" not in listed          # no recursion


def test_names_are_returned_not_paths(folder):
    """Sort keys and the filename filter both match on the name."""
    for name in files_in_folder(str(folder), "HDF5 Nexus"):
        assert "/" not in name and "\\" not in name


def test_uppercase_extensions_are_matched(tmp_path):
    (tmp_path / "LOUD.H5").touch()
    assert files_in_folder(str(tmp_path), "HDF5 Nexus") == ["LOUD.H5"]


def test_missing_folder_yields_empty_list():
    assert files_in_folder("/no/such/folder", "HDF5 Nexus") == []
    assert files_in_folder("", "HDF5 Nexus") == []


def test_unreadable_folder_yields_empty_list(tmp_path, monkeypatch):
    """A folder whose permissions changed shows nothing, it does not raise."""
    def _boom(_path):
        raise PermissionError("nope")

    monkeypatch.setattr("os.listdir", _boom)
    assert files_in_folder(str(tmp_path), "HDF5 Nexus") == []


def test_listing_is_unsorted_by_design(folder):
    """Order comes from file_sorting, so the browsers stay consistent."""
    from pyirena.core.file_sorting import sort_names

    listed = files_in_folder(str(folder), "HDF5 Nexus")
    assert sorted(sort_names(listed, 6)) == sorted(listed)
