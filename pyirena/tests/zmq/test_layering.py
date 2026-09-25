"""Layering guards for ``pyirena/zmq/``.

The service is an optional extra, so two things must stay true:

* ``import pyirena`` works with no pyzmq installed, and so does importing the
  protocol itself — that is what lets the envelope be tested without a socket
  and keeps ``pyirena/tests/zmq/test_protocol.py`` running everywhere.
* No Qt and no matplotlib at module level. The service exists precisely
  because the caller cannot see this machine's screen or filesystem.

Qt is already covered package-wide by ``test_gui_qt_contract.py``; these add
the pyzmq and matplotlib halves.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pyirena

ZMQ_DIR = Path(pyirena.__file__).parent / "zmq"

#: A module-level (unindented) import of a heavy optional dependency.
TOP_LEVEL_IMPORT = re.compile(
    r"^(?:from|import)\s+(matplotlib|PySide6|PyQt6|zmq)\b", re.M
)


def _modules():
    return sorted(p for p in ZMQ_DIR.rglob("*.py") if "__pycache__" not in p.parts)


def test_no_module_level_heavy_imports():
    offenders = {}
    for path in _modules():
        hits = TOP_LEVEL_IMPORT.findall(path.read_text(encoding="utf-8"))
        if hits:
            offenders[path.name] = sorted(set(hits))
    assert not offenders, (
        "pyirena/zmq must import zmq, Qt and matplotlib lazily inside "
        f"functions, not at module level: {offenders}"
    )


def test_the_protocol_imports_and_answers_without_pyzmq():
    """The logic lives where no optional dependency is needed to reach it."""
    script = """
import sys

class Blocked:
    def find_module(self, name, path=None):
        return self if name == "zmq" or name.startswith("zmq.") else None
    def load_module(self, name):
        raise ImportError("pyzmq is blocked for this test")

sys.meta_path.insert(0, Blocked())

import pyirena
import pyirena.zmq
from pyirena.zmq.protocol import handle_request
from pyirena.zmq.options import ServerOptions

assert "zmq" not in sys.modules, "something imported pyzmq"
reply = handle_request("ping", ServerOptions())
assert '"pong": true' in reply, reply
print("OK")
"""
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True,
        cwd=str(Path(pyirena.__file__).parent.parent),
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


def test_importing_the_package_does_not_pull_in_the_service():
    """`import pyirena` must not cost a socket library."""
    result = subprocess.run(
        [sys.executable, "-c",
         "import sys, pyirena; "
         "print('zmq' in sys.modules, 'pyirena.zmq.server' in sys.modules)"],
        capture_output=True, text=True,
        cwd=str(Path(pyirena.__file__).parent.parent),
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "False False"


def test_the_server_explains_itself_when_pyzmq_is_missing(monkeypatch, capsys):
    """A missing extra must name the fix, not raise ModuleNotFoundError."""
    import builtins

    from pyirena.zmq.options import ServerOptions
    from pyirena.zmq.server import run

    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "zmq":
            raise ModuleNotFoundError("No module named 'zmq'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    code = run(ServerOptions())
    monkeypatch.undo()

    assert code == 2
    assert "pip install 'pyirena[zmq]'" in capsys.readouterr().err
