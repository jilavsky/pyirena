"""Text files are read and written as UTF-8, explicitly, everywhere.

Opening a text file without saying which encoding gets the *locale's* codec.
On Linux and macOS that is UTF-8 and nothing goes wrong; on Windows it is
cp1252, and a file holding any character outside it — an angstrom sign, a
chi, a superscript, an em dash, a warning emoji in a doc — raises
``UnicodeDecodeError`` on read, or writes mojibake.

This is not hypothetical and not only a CI problem. pyIrena's config files,
state file and exported reports are written on one machine and read on
another as a matter of routine: a beamline workstation writes a
``pyirena_config.json`` with an angstrom sign in a label, and a Windows
laptop cannot open it. The failure lands on the user, not on the developer
who left the argument out, which is why this is a test and not a convention.

The check parses the source rather than grepping it, so prose that mentions
a call — including this docstring — is not mistaken for one. Binary modes
are exempt: they have no encoding to declare.
"""

from __future__ import annotations

import ast
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]


def _source_files():
    for path in sorted(PACKAGE.rglob("*.py")):
        if "__pycache__" not in path.parts:
            yield path


def _is_binary(call: ast.Call) -> bool:
    """True when this ``open()`` was given a binary mode."""
    mode = None
    if len(call.args) >= 2 and isinstance(call.args[1], ast.Constant):
        mode = call.args[1].value
    for kw in call.keywords:
        if kw.arg == "mode" and isinstance(kw.value, ast.Constant):
            mode = kw.value.value
    return isinstance(mode, str) and "b" in mode


def _offenders(predicate) -> list[str]:
    found = []
    for path in _source_files():
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except SyntaxError:                                # pragma: no cover
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if not predicate(node):
                continue
            if any(kw.arg == "encoding" for kw in node.keywords):
                continue
            found.append(f"{path.relative_to(PACKAGE.parent)}:{node.lineno}")
    return found


def _is_builtin_open(node: ast.Call) -> bool:
    """A bare ``open(...)`` — not ``h5py.File`` and not ``something.open()``."""
    return (isinstance(node.func, ast.Name) and node.func.id == "open"
            and not _is_binary(node))


def _is_path_text_io(node: ast.Call) -> bool:
    return (isinstance(node.func, ast.Attribute)
            and node.func.attr in ("read_text", "write_text"))


def test_every_text_open_declares_utf8():
    offenders = _offenders(_is_builtin_open)
    assert not offenders, (
        "text-mode open() without encoding='utf-8'. These read as cp1252 on "
        "Windows and raise UnicodeDecodeError on the first non-Latin-1 "
        "character:\n  " + "\n  ".join(offenders)
    )


def test_every_read_text_and_write_text_declares_utf8():
    offenders = _offenders(_is_path_text_io)
    assert not offenders, (
        "Path.read_text()/write_text() without encoding='utf-8' — same "
        "problem, same platform:\n  " + "\n  ".join(offenders)
    )


def test_the_documentation_really_is_utf8():
    """The docs carry angstroms, chi-squareds and an emoji.

    ``docs/batch_api.md`` broke the Windows CI job exactly this way: a test
    read it with the locale codec and hit the 0x8f byte of a warning sign.
    """
    import pytest

    docs = PACKAGE.parent / "docs"
    if not docs.is_dir():
        pytest.skip("docs/ not present")
    for path in sorted(docs.glob("*.md")):
        path.read_text(encoding="utf-8")      # must not raise
