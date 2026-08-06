"""
Tests for pyirena.diagnostics — the dependency/environment troubleshooter.

These cover the distinction that motivated the module: a package that is
*absent* and a package that is *present but unloadable* both raise
``ImportError``, and pyIrena must not report the second as the first.

Headless and Qt-free: everything is exercised through injected fake modules.
"""
import builtins
import importlib
import sys

import pytest

from pyirena import diagnostics


# --------------------------------------------------------------------------
# probe()
# --------------------------------------------------------------------------

def test_probe_reports_ok_for_stdlib_module():
    status = diagnostics.probe("json")
    assert status.status == "ok"
    assert status.ok
    assert status.path is not None


def test_probe_reports_missing_for_absent_module():
    status = diagnostics.probe("pyirena_definitely_not_a_real_module")
    assert status.status == "missing"
    assert not status.ok
    assert "NOT INSTALLED" in status.summary()


def test_probe_reports_broken_when_import_raises(monkeypatch):
    """A module whose spec exists but whose import fails is 'broken'."""
    boom = ImportError("dlopen(...): incompatible architecture")

    monkeypatch.setattr(diagnostics, "_find_spec", lambda name: object())

    def fake_import(name):
        raise boom

    monkeypatch.setattr(diagnostics.importlib, "import_module", fake_import)

    status = diagnostics.probe("PySide6")
    assert status.status == "broken"
    assert status.error is boom
    assert "FAILS TO LOAD" in status.summary()


def test_probe_survives_non_importerror(monkeypatch):
    """A module that raises something exotic must not crash the probe."""
    monkeypatch.setattr(diagnostics, "_find_spec", lambda name: object())

    def fake_import(name):
        raise RuntimeError("Qt platform plugin exploded")

    monkeypatch.setattr(diagnostics.importlib, "import_module", fake_import)

    status = diagnostics.probe("PySide6")
    assert status.status == "broken"
    assert isinstance(status.error, RuntimeError)


# --------------------------------------------------------------------------
# explain_import_error()
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "message, expected_fragment",
    [
        ("mach-o file, but is an incompatible architecture", "Architecture mismatch"),
        ("libGL.so.1: cannot open shared object file", "libgl1"),
        ("libxkbcommon-x11.so.0: cannot open shared object", "libxkbcommon-x11-0"),
        ("DLL load failed while importing QtCore", "Visual C++"),
        ("Symbol not found: __ZN9QMetaType", "mixed install"),
        ("shiboken6 version mismatch", "shiboken6"),
        ("/lib/x86_64-linux-gnu/libc.so.6: version `GLIBC_2.28' not found", "C library"),
    ],
)
def test_explain_import_error_recognises_common_causes(message, expected_fragment):
    hint = diagnostics.explain_import_error(ImportError(message))
    assert hint is not None
    assert expected_fragment.lower() in hint.lower()


def test_explain_import_error_returns_none_for_unknown_message():
    assert diagnostics.explain_import_error(ImportError("something novel")) is None


# --------------------------------------------------------------------------
# format_qt_import_failure()
# --------------------------------------------------------------------------

def test_qt_failure_message_when_truly_missing():
    statuses = [
        diagnostics.DependencyStatus("PySide6", "missing"),
        diagnostics.DependencyStatus("PyQt6", "missing"),
    ]
    msg = diagnostics.format_qt_import_failure(statuses)
    assert "No Qt bindings found" in msg
    assert "python -m pip install" in msg
    # The interpreter must be named: this is what resolves wrong-env reports.
    assert sys.executable in msg


def test_qt_failure_message_when_installed_but_broken():
    """The key regression: do not tell users to install what they have."""
    broken = diagnostics.DependencyStatus(
        "PySide6",
        "broken",
        path="/some/site-packages/PySide6/__init__.py",
        error=ImportError("incompatible architecture (have 'x86_64', need 'arm64')"),
    )
    msg = diagnostics.format_qt_import_failure(
        [broken, diagnostics.DependencyStatus("PyQt6", "missing")]
    )
    assert "IS installed" in msg
    assert "reinstalling it will probably not help" in msg
    assert "/some/site-packages/PySide6/__init__.py" in msg
    assert "Architecture mismatch" in msg
    # Must NOT give the misleading install instruction.
    assert "No Qt bindings found" not in msg


def test_qt_failure_message_mentions_doctor():
    msg = diagnostics.format_qt_import_failure(
        [
            diagnostics.DependencyStatus("PySide6", "missing"),
            diagnostics.DependencyStatus("PyQt6", "missing"),
        ]
    )
    assert "pyirena-doctor" in msg


# --------------------------------------------------------------------------
# format_gui_import_failure()
# --------------------------------------------------------------------------

def test_gui_failure_names_the_optional_dependency():
    exc = ModuleNotFoundError("No module named 'pyqtgraph'", name="pyqtgraph")
    msg = diagnostics.format_gui_import_failure(exc)
    assert "pyqtgraph" in msg
    assert "pyirena[gui]" in msg


def test_gui_failure_treats_unknown_module_as_a_bug():
    exc = ModuleNotFoundError(
        "No module named 'pyirena.gui.typo'", name="pyirena.gui.typo"
    )
    msg = diagnostics.format_gui_import_failure(exc)
    assert "bug" in msg.lower()
    assert "issues" in msg
    # Must not blame the user's dependencies for an internal error.
    assert "pip install" not in msg


def test_gui_failure_diagnoses_a_raw_qt_importerror():
    """A bare PySide6 ImportError gets the full Qt diagnosis attached."""
    exc = ImportError("No module named 'PySide6'", name="PySide6")
    msg = diagnostics.format_gui_import_failure(exc)
    assert "could not load a Qt binding" in msg
    assert sys.executable in msg


def test_gui_failure_passes_existing_qt_diagnosis_through_unchanged():
    """An already-diagnosed message is not wrapped in a second header."""
    diagnosed = (
        "pyIrena could not load a Qt binding.\n\nPySide6 IS installed but ..."
    )
    msg = diagnostics.format_gui_import_failure(ImportError(diagnosed))
    assert msg == diagnosed
    assert msg.count("could not load a Qt binding") == 1


# --------------------------------------------------------------------------
# environment + report
# --------------------------------------------------------------------------

def test_environment_lines_include_interpreter_and_platform():
    lines = "\n".join(diagnostics.environment_lines())
    assert sys.executable in lines
    assert "python" in lines
    assert "machine" in lines


def test_report_runs_and_mentions_required_dependencies():
    text = "\n".join(diagnostics._report())
    for name in diagnostics.REQUIRED_DEPENDENCIES:
        assert name in text
    assert "pyIrena environment report" in text


def test_main_returns_exit_code(capsys):
    code = diagnostics.main()
    out = capsys.readouterr().out
    assert "pyIrena environment report" in out
    # 0 when a Qt binding loads, 1 otherwise — both are valid here since CI
    # may or may not have the gui extra installed.
    assert code in (0, 1)
    assert code == (0 if any(s.ok for s in diagnostics.probe_qt()) else 1)


# --------------------------------------------------------------------------
# _qt.py integration: the shim must not claim "not installed" wrongly
# --------------------------------------------------------------------------

def test_qt_shim_raises_diagnosed_error_when_bindings_absent(monkeypatch):
    """Simulate an environment with neither binding and re-import the shim."""
    real_import = builtins.__import__

    def blocked_import(name, *args, **kwargs):
        if name.split(".")[0] in ("PySide6", "PyQt6"):
            raise ImportError(f"No module named {name!r}", name=name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked_import)
    for mod in list(sys.modules):
        if mod.startswith(("PySide6", "PyQt6", "pyirena.gui._qt")):
            monkeypatch.delitem(sys.modules, mod, raising=False)

    with pytest.raises(ImportError) as excinfo:
        importlib.import_module("pyirena.gui._qt")

    message = str(excinfo.value)
    assert "could not load a Qt binding" in message
    assert "pyirena-doctor" in message
    # The old message must be gone — it was the source of the confusion.
    assert "Neither PySide6 nor PyQt6 found" not in message
