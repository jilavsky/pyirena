"""Tests for pyirena.api._paths."""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from pyirena.api._paths import (
    PathSecurityError, resolve_safe, resolve_safe_file, resolve_safe_folder,
)


def test_resolve_safe_absolute(tmp_path):
    f = tmp_path / "x.txt"
    f.write_text("hi")
    out = resolve_safe(str(f))
    assert out == f.resolve()


def test_resolve_safe_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        resolve_safe(str(tmp_path / "nope.txt"))


def test_resolve_safe_must_exist_false(tmp_path):
    p = resolve_safe(str(tmp_path / "future.png"), must_exist=False)
    assert p == (tmp_path / "future.png").resolve()


def test_resolve_safe_file_rejects_dir(tmp_path):
    with pytest.raises(IsADirectoryError):
        resolve_safe_file(str(tmp_path))


def test_resolve_safe_folder_rejects_file(tmp_path):
    f = tmp_path / "x.txt"
    f.write_text("hi")
    with pytest.raises(NotADirectoryError):
        resolve_safe_folder(str(f))


def test_data_root_enforced(tmp_path, monkeypatch):
    inside = tmp_path / "inside"
    inside.mkdir()
    target = inside / "ok.txt"
    target.write_text("hi")

    outside = tmp_path / "outside"
    outside.mkdir()
    forbidden = outside / "no.txt"
    forbidden.write_text("nope")

    monkeypatch.setenv("PYIRENA_DATA_ROOT", str(inside))

    # Inside: OK
    out = resolve_safe(str(target))
    assert out == target.resolve()

    # Outside: rejected
    with pytest.raises(PathSecurityError):
        resolve_safe(str(forbidden))


def test_data_root_relative_resolution(tmp_path, monkeypatch):
    inside = tmp_path / "inside"
    inside.mkdir()
    (inside / "f.txt").write_text("hi")
    monkeypatch.setenv("PYIRENA_DATA_ROOT", str(inside))

    # Relative path resolves against PYIRENA_DATA_ROOT
    out = resolve_safe("f.txt")
    assert out == (inside / "f.txt").resolve()


def test_data_root_traversal_rejected(tmp_path, monkeypatch):
    inside = tmp_path / "inside"
    inside.mkdir()
    monkeypatch.setenv("PYIRENA_DATA_ROOT", str(inside))
    with pytest.raises(PathSecurityError):
        resolve_safe("../escape.txt", must_exist=False)
