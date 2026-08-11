"""
One Qt import point, enforced.

``pyirena/gui/_qt.py`` exists so that the PySide6/PyQt6 fallback is written
once.  The rule erodes the same way every time: a function needs one more
class, the author writes another local try/except around the two bindings
because it is only three lines, and the project drifts back to a dozen
half-maintained shims — one of which still had a stale PyQt5 branch that would
have been an outright bug on a PyQt6-only install.

These tests read the source rather than importing it, so they run on a machine
with no Qt at all and fail the build the moment a new binding import appears.
"""

from __future__ import annotations

import re
from pathlib import Path

import pyirena

PACKAGE = Path(pyirena.__file__).parent
SHIM = PACKAGE / "gui" / "_qt.py"

#: ``from PySide6…`` / ``import PyQt6…`` at any indentation.  Matching the
#: statement (not the bare word) lets diagnostics and skip messages keep naming
#: the bindings in strings, which is exactly what they are for.
BINDING_IMPORT = re.compile(r"^\s*(?:from|import)\s+(PySide[26]|PyQt[456])\b", re.M)


def _python_files():
    return sorted(p for p in PACKAGE.rglob("*.py") if "__pycache__" not in p.parts)


def test_only_the_shim_imports_a_qt_binding():
    offenders = {}
    for path in _python_files():
        if path == SHIM:
            continue
        hits = BINDING_IMPORT.findall(path.read_text(encoding="utf-8"))
        if hits:
            offenders[str(path.relative_to(PACKAGE))] = sorted(set(hits))

    assert not offenders, (
        "Qt bindings must be imported through pyirena.gui._qt, not directly:\n"
        + "\n".join(f"  {name}: {found}" for name, found in offenders.items())
        + "\nAdd any missing class to pyirena/gui/_qt.py and import it from there."
    )


def test_there_is_exactly_one_shim():
    """A second ``_qt.py`` in a subpackage is the other way this drifts."""
    shims = [p for p in _python_files() if p.name == "_qt.py"]
    assert shims == [SHIM], f"expected one Qt shim, found {[str(p) for p in shims]}"


def test_no_module_defines_its_own_binding_fallback():
    """Catch the pattern even if someone reaches the binding another way."""
    pattern = re.compile(r"importlib\.import_module\(\s*['\"](PySide[26]|PyQt[456])",
                         re.M)
    offenders = [str(p.relative_to(PACKAGE)) for p in _python_files()
                 if p != SHIM and p.name != "diagnostics.py"
                 and pattern.search(p.read_text(encoding="utf-8"))]
    assert not offenders, f"dynamic Qt binding imports outside the shim: {offenders}"


def test_the_shim_exports_everything_the_gui_imports_from_it():
    """``__all__`` is what readers trust; keep it honest.

    A name added as a module attribute but left out of ``__all__`` still works
    for ``from … import X`` and so goes unnoticed until someone uses a star
    import or reads the list to see what is available.
    """
    source = SHIM.read_text(encoding="utf-8")
    assigned = set(re.findall(r"^(Q\w+|Signal|Qt) = Qt(?:Widgets|Core|Gui)\.", source, re.M))
    listed = set(re.findall(r'"(\w+)"', source.split("__all__")[1]))
    missing = sorted(assigned - listed)
    assert not missing, f"names in _qt.py but absent from __all__: {missing}"


def test_every_name_the_gui_asks_for_actually_exists():
    """The shim must cover what the package imports from it.

    Without this a typo or a name that was never added surfaces as an
    ImportError at GUI start-up on a user's machine, not here.
    """
    wanted: dict[str, set[str]] = {}
    block = re.compile(r"from pyirena\.gui\._qt import \(([^)]*)\)|"
                       r"from pyirena\.gui\._qt import ([^\n(]+)")
    for path in _python_files():
        if path == SHIM:
            continue
        for parens, plain in block.findall(path.read_text(encoding="utf-8")):
            names = (parens or plain).replace("\n", " ").split(",")
            for raw in names:
                name = raw.split("#")[0].split(" as ")[0].strip()
                if name:
                    wanted.setdefault(name, set()).add(str(path.relative_to(PACKAGE)))

    source = SHIM.read_text(encoding="utf-8")
    defined = set(re.findall(r"^(\w+) = Qt(?:Widgets|Core|Gui)\.", source, re.M))
    defined |= {"QtWidgets", "QtCore", "QtGui", "Signal", "QT_BINDING"}

    missing = {n: sorted(f) for n, f in wanted.items() if n not in defined}
    assert not missing, f"imported from _qt.py but not defined there: {missing}"
