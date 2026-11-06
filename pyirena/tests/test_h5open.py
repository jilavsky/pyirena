"""
HDF5 must never hand back one file when another was asked for.

HDF5 recognises an already-open file by its filesystem identity, not its
name.  On SSHFS-Win every file reports identity 0, so with one HDF5 file open
any other file opened natively *is* that open file: Collect Values at the
beamline returned the browsed sample's transmission for every selected file.
``pyirena/io/h5open.py`` explains the mechanism; these tests pin the fix.

Hard links reproduce the real HDF5 behaviour on any POSIX disk -- two names,
one identity.  A patched ``_identity`` reproduces SSHFS-Win's identity 0 for
files whose contents differ, which is the case that returned wrong numbers.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import h5py
import pytest

import pyirena
from pyirena.io import h5open
from pyirena.io.h5open import open_h5

PACKAGE = Path(pyirena.__file__).parent


def _readers():
    """The Data Explorer readers; their package imports Qt, so skip without it."""
    return pytest.importorskip("pyirena.gui.hdf5viewer.pyirena_readers")


def _write_sample(path: Path, sample: str, transmission: float) -> None:
    with h5py.File(path, "w") as f:
        f.create_group(f"entry/{sample}/sastransmission").attrs["T"] = transmission


@pytest.fixture
def same_identity_everywhere(monkeypatch):
    """Make every file report device 1, file id 0 -- as SSHFS-Win does."""
    monkeypatch.setattr(h5open, "_identity", lambda path: (1, 0))


def test_hdf5_itself_confuses_files_with_the_same_identity(tmp_path):
    """The premise: a native open of a same-identity file returns the open one."""
    a, b = tmp_path / "a.h5", tmp_path / "b.h5"
    _write_sample(a, "sampleA", 0.11)
    try:
        os.link(a, b)
    except (OSError, NotImplementedError):
        pytest.skip("hard links not supported here")
    with h5py.File(a, "r"):
        with pytest.raises(OSError, match="already open"):
            h5py.File(b, "r+")


def test_open_h5_opens_a_same_identity_file_independently(tmp_path):
    a, b = tmp_path / "a.h5", tmp_path / "b.h5"
    _write_sample(a, "sampleA", 0.11)
    try:
        os.link(a, b)
    except (OSError, NotImplementedError):
        pytest.skip("hard links not supported here")
    with open_h5(a, "r"):
        assert h5open.confusable_open_file(b) == str(a)
        with open_h5(b, "r+") as fb:       # natively this is refused, see above
            assert fb.filename == str(b)


def test_nothing_else_open_means_a_plain_native_open(tmp_path, same_identity_everywhere):
    a = tmp_path / "a.h5"
    _write_sample(a, "sampleA", 0.11)
    with open_h5(a, "r") as f:
        assert type(f) is h5py.File


def test_reopening_the_same_file_stays_native(tmp_path, same_identity_everywhere):
    a = tmp_path / "a.h5"
    _write_sample(a, "sampleA", 0.11)
    with open_h5(a, "r"), open_h5(a, "r") as again:
        assert type(again) is h5py.File


def test_collect_reads_each_file_while_another_is_browsed(tmp_path, same_identity_everywhere):
    """The beamline case: the browser holds one sample open, Collect reads others."""
    readers = _readers()
    values = {"S1": 0.11, "S2": 0.22, "S3": 0.33}
    paths = []
    for sample, t in values.items():
        p = tmp_path / f"{sample}.h5"
        _write_sample(p, sample, t)
        paths.append(p)

    with open_h5(paths[0], "r"):           # the file shown in the Data Explorer
        got = []
        for p in paths:
            with readers.opened_or_path(p) as h5:
                got.append(readers.read_metadata_value(h5, "/entry/*/sastransmission@T"))
    assert got == list(values.values())


def test_writes_land_in_the_file_asked_for(tmp_path, same_identity_everywhere):
    a, b, new = tmp_path / "a.h5", tmp_path / "b.h5", tmp_path / "new.h5"
    _write_sample(a, "sampleA", 0.11)
    _write_sample(b, "sampleB", 0.22)
    with open_h5(a, "r"):
        with open_h5(b, "a") as fb:
            fb.attrs["saved"] = "here"
        with open_h5(new, "w") as fn:      # does not exist yet: created
            fn.attrs["saved"] = "new"
    with h5py.File(a, "r") as f:
        assert "saved" not in f.attrs
    with h5py.File(b, "r") as f:
        assert f.attrs["saved"] == "here"
    with h5py.File(new, "r") as f:
        assert f.attrs["saved"] == "new"


def test_fallback_file_releases_its_os_handle(tmp_path, same_identity_everywhere):
    a, b = tmp_path / "a.h5", tmp_path / "b.h5"
    _write_sample(a, "sampleA", 0.11)
    _write_sample(b, "sampleB", 0.22)
    with open_h5(a, "r"):
        fb = open_h5(b, "r")
        assert isinstance(fb, h5open._FileObjectH5)
        fb.close()
        assert fb._pyirena_fileobj.closed


def test_readers_leave_a_passed_open_file_open(tmp_path):
    readers = _readers()
    a = tmp_path / "a.h5"
    _write_sample(a, "sampleA", 0.11)
    with open_h5(a, "r") as f:
        assert readers.read_metadata_value(f, "/entry/*/sastransmission@T") == 0.11
        assert readers.read_hdf5_item(f, "/entry/sampleA/sastransmission@T") == (0.11, "scalar")
        assert f.id.valid


def test_every_hdf5_open_goes_through_open_h5():
    """A direct ``h5py.File(`` call reintroduces the wrong-file reads."""
    call = re.compile(r"\bh5py\.File\(")
    offenders = sorted(
        str(p.relative_to(PACKAGE))
        for p in PACKAGE.rglob("*.py")
        if "__pycache__" not in p.parts
        and "tests" not in p.relative_to(PACKAGE).parts
        and p != PACKAGE / "io" / "h5open.py"
        and call.search(p.read_text(encoding="utf-8"))
    )
    assert not offenders, f"open HDF5 files with pyirena.io.h5open.open_h5: {offenders}"
