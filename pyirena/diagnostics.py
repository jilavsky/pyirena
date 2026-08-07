"""
Environment and dependency diagnostics for pyIrena.

This module exists because the single most common support request is some
variant of *"pip says PySide6 is installed but pyIrena says it isn't"*.  That
report is usually accurate: ``ImportError`` covers two very different
failures, and the old error message conflated them.

    ModuleNotFoundError   the package really is absent from *this* interpreter
    ImportError (other)   the package is present but its binaries will not
                          load — wrong CPU architecture, missing system
                          libraries, a half-finished install, a shiboken6 /
                          PySide6 version mismatch

Telling a user to reinstall something that is already installed sends them in
a circle.  The helpers here identify which failure actually occurred, attach
the interpreter that hit it, and translate the cryptic dynamic-loader message
into an actionable hint.

Standard library only, deliberately: this code has to run *when the
dependencies are broken*.  Nothing here may import numpy, h5py, Qt, or
matplotlib at module scope.

Entry point: ``pyirena-doctor`` (see :func:`main`).
"""
from __future__ import annotations

import importlib
import importlib.util
import os
import platform
import shutil
import sys
import sysconfig
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

__all__ = [
    "DependencyStatus",
    "probe",
    "probe_qt",
    "explain_import_error",
    "environment_lines",
    "format_qt_import_failure",
    "format_gui_import_failure",
    "main",
]

# Python versions this project actually supports.  pyproject.toml allows 3.9+;
# the GUI stack (PySide6 wheels) is only reliable through 3.13.
PYTHON_MIN = (3, 9)
PYTHON_GUI_MAX_TESTED = (3, 13)

# module name -> (extra that provides it, what it is for)
OPTIONAL_DEPENDENCIES: Dict[str, Tuple[str, str]] = {
    "PySide6": ("gui", "Qt bindings (preferred)"),
    "PyQt6": ("gui", "Qt bindings (fallback)"),
    "pyqtgraph": ("gui", "interactive plots"),
    "matplotlib": ("plotting", "static/headless plots"),
    "xraydb": ("gui", "X-ray scattering contrast"),
    "periodictable": ("gui", "neutron scattering contrast"),
    "Dans_Diffraction": ("gui", "diffraction line markers"),
    "pyvista": ("gui3d", "3D voxelgram viewer"),
    "pyvistaqt": ("gui3d", "3D viewer Qt embedding"),
    "mcp": ("mcp", "MCP server"),
    "anthropic": ("gui", "AI advisor (optional)"),
    "openai": ("gui", "AI advisor (optional)"),
}

# Required (non-optional) runtime dependencies.
REQUIRED_DEPENDENCIES: Tuple[str, ...] = ("numpy", "scipy", "h5py", "igor2")

# Substring -> advice.  Matched case-insensitively against the exception text
# of a *broken* (not missing) import.  Ordered: first match wins.
_ERROR_HINTS: Sequence[Tuple[str, str]] = (
    (
        "incompatible architecture",
        "Architecture mismatch: the installed wheel was built for a different "
        "CPU than this Python. On Apple Silicon this usually means an x86_64 "
        "wheel under an arm64 interpreter (or a terminal running under "
        "Rosetta). Check `python -c \"import platform; print(platform.machine())\"` "
        "and reinstall with a matching interpreter.",
    ),
    (
        "mach-o",
        "The macOS dynamic loader rejected the installed binary — usually an "
        "architecture mismatch or a truncated download. Try: "
        "python -m pip install --force-reinstall --no-cache-dir PySide6",
    ),
    (
        "libgl.so.1",
        "Linux is missing OpenGL runtime libraries. On Debian/Ubuntu: "
        "sudo apt install libgl1 libglib2.0-0",
    ),
    (
        "libxkbcommon",
        "Linux is missing X11/keyboard libraries needed by the Qt platform "
        "plugin. On Debian/Ubuntu: sudo apt install libxkbcommon-x11-0 "
        "libxcb-cursor0 libxcb-icccm4 libxcb-keysyms1 libxcb-shape0",
    ),
    (
        "libegl",
        "Linux is missing EGL libraries. On Debian/Ubuntu: "
        "sudo apt install libegl1",
    ),
    (
        "glibc",
        "The system C library is older than the wheel requires. Install an "
        "older PySide6 (e.g. PySide6==6.5.*) or update the OS.",
    ),
    (
        "dll load failed",
        "Windows could not load the Qt DLLs. Install the Microsoft Visual C++ "
        "Redistributable (x64), then: "
        "python -m pip install --force-reinstall PySide6",
    ),
    (
        "shiboken",
        "shiboken6 does not match the installed PySide6. They must be "
        "installed together: "
        "python -m pip install --force-reinstall PySide6 shiboken6",
    ),
    (
        "symbol not found",
        "A partially upgraded or mixed install. Reinstall the binding cleanly: "
        "python -m pip uninstall -y PySide6 shiboken6 && "
        "python -m pip install PySide6",
    ),
)


@dataclass
class DependencyStatus:
    """Result of probing one importable module.

    Attributes:
        name: Module name as imported.
        status: ``"ok"``, ``"missing"`` or ``"broken"``.
        version: ``__version__`` if the module imported and exposes one.
        path: ``__file__`` of the imported module, or the spec origin when the
            module was found but failed to load.
        error: The exception raised, for ``"broken"``.
    """

    name: str
    status: str
    version: Optional[str] = None
    path: Optional[str] = None
    error: Optional[BaseException] = None

    @property
    def ok(self) -> bool:
        return self.status == "ok"

    def summary(self) -> str:
        """One-line human-readable status."""
        if self.status == "ok":
            return f"{self.name} {self.version or '(version unknown)'}"
        if self.status == "missing":
            return f"{self.name} NOT INSTALLED"
        return f"{self.name} INSTALLED BUT FAILS TO LOAD: {self.error}"


def _find_spec(name: str):
    """``importlib.util.find_spec`` that never raises."""
    try:
        return importlib.util.find_spec(name)
    except (ImportError, AttributeError, ValueError):
        return None


def probe(name: str) -> DependencyStatus:
    """Determine whether *name* is absent, importable, or present-but-broken.

    Args:
        name: Module name, e.g. ``"PySide6"``.

    Returns:
        A :class:`DependencyStatus`. ``"missing"`` means no distribution
        provides the module for this interpreter; ``"broken"`` means the files
        are on disk but importing them raised.
    """
    spec = _find_spec(name)
    if spec is None:
        return DependencyStatus(name=name, status="missing")

    origin = getattr(spec, "origin", None)
    try:
        module = importlib.import_module(name)
    except BaseException as exc:  # noqa: BLE001 - report anything, never crash
        return DependencyStatus(
            name=name, status="broken", path=origin, error=exc
        )

    return DependencyStatus(
        name=name,
        status="ok",
        version=getattr(module, "__version__", None),
        path=getattr(module, "__file__", origin),
    )


def probe_qt() -> List[DependencyStatus]:
    """Probe both supported Qt bindings, PySide6 first."""
    return [probe("PySide6"), probe("PyQt6")]


def explain_import_error(exc: BaseException) -> Optional[str]:
    """Translate a dynamic-loader error message into actionable advice.

    Args:
        exc: The exception raised by a failed import.

    Returns:
        A hint string, or ``None`` when the message is not recognised.
    """
    text = str(exc).lower()
    for needle, advice in _ERROR_HINTS:
        if needle in text:
            return advice
    return None


def _python_version_warning() -> Optional[str]:
    """Warn when the interpreter is outside the range the GUI is tested on."""
    v = sys.version_info[:2]
    if v < PYTHON_MIN:
        return (
            f"Python {v[0]}.{v[1]} is below the minimum supported "
            f"{PYTHON_MIN[0]}.{PYTHON_MIN[1]}."
        )
    if v > PYTHON_GUI_MAX_TESTED:
        return (
            f"Python {v[0]}.{v[1]} is newer than the newest version the GUI "
            f"stack is tested against ({PYTHON_GUI_MAX_TESTED[0]}."
            f"{PYTHON_GUI_MAX_TESTED[1]}). Qt wheels are often unavailable or "
            f"unstable on brand-new Python releases — a "
            f"{PYTHON_GUI_MAX_TESTED[0]}.{PYTHON_GUI_MAX_TESTED[1]} "
            f"environment is the safe choice."
        )
    return None


def environment_lines() -> List[str]:
    """Interpreter/platform facts worth pasting into a bug report."""
    try:
        from pyirena import __version__ as pyirena_version
    except Exception:  # pragma: no cover - only if the package is broken
        pyirena_version = "unknown"

    lines = [
        f"pyirena          : {pyirena_version}",
        f"python           : {platform.python_version()} ({sys.executable})",
        f"platform         : {platform.platform()}",
        f"machine          : {platform.machine()}",
    ]

    # A venv/conda env is the usual reason pip and the launcher disagree.
    in_venv = sys.prefix != getattr(sys, "base_prefix", sys.prefix)
    env_name = os.environ.get("CONDA_DEFAULT_ENV")
    if env_name:
        lines.append(f"conda env        : {env_name}")
    lines.append(
        f"environment      : {'virtualenv/venv' if in_venv else 'system/base'}"
        f" (prefix {sys.prefix})"
    )
    return lines


def _launcher_mismatch() -> Optional[str]:
    """Detect the classic 'pip and the launcher use different Pythons' bug.

    Reads the shebang of the installed ``pyirena-gui`` console script and
    compares it with the running interpreter.

    Returns:
        A warning string when they differ, else ``None``.
    """
    script = shutil.which("pyirena-gui")
    if not script:
        return None
    # Windows console scripts are .exe wrappers with no readable shebang.
    if os.name == "nt" or script.lower().endswith(".exe"):
        return None
    try:
        with open(script, "r", encoding="utf-8", errors="replace") as fh:
            first = fh.readline().strip()
    except OSError:
        return None
    if not first.startswith("#!"):
        return None
    interpreter = first[2:].strip().strip('"')
    # "#!/usr/bin/env python" tells us nothing definite.
    if not interpreter or interpreter.endswith("env python"):
        return None
    if os.path.realpath(interpreter) == os.path.realpath(sys.executable):
        return None
    return (
        "The 'pyirena-gui' launcher runs a DIFFERENT Python than the one you "
        "are using now:\n"
        f"    launcher uses : {interpreter}\n"
        f"    you are using : {sys.executable}\n"
        "Packages installed with one are invisible to the other. Use\n"
        "    python -m pip install 'pyirena[gui]'\n"
        "    python -m pyirena.gui.launch\n"
        "so that installation and launch always share one interpreter."
    )


def format_qt_import_failure(
    statuses: Optional[Sequence[DependencyStatus]] = None,
) -> str:
    """Build the message raised when no Qt binding can be imported.

    Distinguishes "not installed" from "installed but will not load" — the
    distinction the old message threw away.

    Args:
        statuses: Result of :func:`probe_qt`; probed here when omitted.

    Returns:
        A multi-line message suitable for an ``ImportError``.
    """
    statuses = list(statuses) if statuses is not None else probe_qt()
    broken = [s for s in statuses if s.status == "broken"]
    parts: List[str] = []

    if broken:
        first = broken[0]
        parts.append(
            f"{first.name} IS installed but could not be loaded — so "
            f"reinstalling it will probably not help."
        )
        parts.append(f"    location : {first.path}")
        parts.append(f"    error    : {type(first.error).__name__}: {first.error}")
        hint = explain_import_error(first.error) if first.error else None
        if hint:
            parts.append(f"    likely cause: {hint}")
    else:
        parts.append(
            "No Qt bindings found for this Python interpreter. Install with:"
        )
        parts.append("    python -m pip install 'pyirena[gui]'")
        parts.append(
            "(Use 'python -m pip', not bare 'pip' — that guarantees the "
            "packages land in the interpreter shown below. Quote the brackets; "
            "zsh expands them otherwise.)"
        )

    version_warning = _python_version_warning()
    if version_warning:
        parts.append(f"\nNote: {version_warning}")

    mismatch = _launcher_mismatch()
    if mismatch:
        parts.append("\n" + mismatch)

    parts.append("\nEnvironment:")
    parts.extend("    " + line for line in environment_lines())
    parts.append("\nFull dependency report: run 'pyirena-doctor'")
    return "\n".join(parts)


def format_gui_import_failure(exc: BaseException) -> str:
    """Build the message shown when launching the GUI fails on an import.

    The failing module may be a Qt binding, another optional GUI dependency,
    or an internal module (i.e. a genuine bug). These need different advice,
    so identify which happened rather than blaming "GUI dependencies".

    Args:
        exc: The ``ImportError`` raised while importing the GUI.

    Returns:
        A multi-line message for the terminal.
    """
    name = getattr(exc, "name", None) or ""
    root = name.split(".")[0]

    # Already diagnosed by pyirena.gui._qt — pass it through unchanged rather
    # than wrapping a second header around a complete message.
    if "Qt binding" in str(exc):
        return str(exc)

    if root in ("PySide6", "PyQt6"):
        return "pyIrena could not load a Qt binding.\n\n" + format_qt_import_failure()

    if root in OPTIONAL_DEPENDENCIES:
        extra, purpose = OPTIONAL_DEPENDENCIES[root]
        status = probe(root)
        lines = [
            f"pyIrena could not start the GUI: '{root}' ({purpose}) is "
            f"unavailable.",
            "",
        ]
        if status.status == "broken":
            lines.append(
                f"{root} is installed at {status.path} but fails to load:"
            )
            lines.append(f"    {type(status.error).__name__}: {status.error}")
            hint = explain_import_error(status.error) if status.error else None
            if hint:
                lines.append(f"    likely cause: {hint}")
        else:
            lines.append(f"Install it with:  python -m pip install 'pyirena[{extra}]'")
        lines.append("")
        lines.append("Environment:")
        lines.extend("    " + line for line in environment_lines())
        mismatch = _launcher_mismatch()
        if mismatch:
            lines.append("")
            lines.append(mismatch)
        lines.append("")
        lines.append("Full dependency report: run 'pyirena-doctor'")
        return "\n".join(lines)

    # Not a known dependency: most likely a bug inside pyirena itself.
    lines = [
        "pyIrena could not start the GUI because of an unexpected import "
        "error.",
        "",
        f"    {type(exc).__name__}: {exc}",
        "",
        "This looks like a bug rather than a missing package. The full "
        "traceback follows, and is also written to:",
        f"    {_log_path()}",
        "",
        "Please report it at https://github.com/jilavsky/pyirena/issues",
    ]
    return "\n".join(lines)


def _log_path() -> str:
    """Best-effort path to the GUI log file."""
    try:
        from pyirena.logging_setup import get_log_dir

        return str(get_log_dir() / "gui.log")
    except Exception:  # pragma: no cover - only if logging setup is broken
        return "~/.pyirena/logs/gui.log"


# --------------------------------------------------------------------------
# pyirena-doctor
# --------------------------------------------------------------------------

def _report() -> List[str]:
    """Assemble the full diagnostic report."""
    out: List[str] = []
    out.append("pyIrena environment report")
    out.append("=" * 60)
    out.extend(environment_lines())

    try:
        import pyirena

        out.append(f"package path     : {os.path.dirname(pyirena.__file__)}")
    except Exception as exc:  # pragma: no cover
        out.append(f"package path     : UNAVAILABLE ({exc})")

    out.append(f"scripts dir      : {sysconfig.get_path('scripts')}")
    for tool in ("pyirena-gui", "pyirena-viewer", "pyirena-mcp"):
        found = shutil.which(tool)
        out.append(f"  {tool:<16}: {found or 'not on PATH'}")

    out.append("")
    out.append("Required dependencies")
    out.append("-" * 60)
    problems: List[DependencyStatus] = []
    for name in REQUIRED_DEPENDENCIES:
        status = probe(name)
        out.append(f"  {status.summary()}")
        if not status.ok:
            problems.append(status)

    out.append("")
    out.append("Optional dependencies")
    out.append("-" * 60)
    for name, (extra, purpose) in OPTIONAL_DEPENDENCIES.items():
        status = probe(name)
        marker = "  " if status.ok else "! "
        out.append(f"{marker}{status.summary()}   [{extra}: {purpose}]")
        if status.status == "broken":
            problems.append(status)

    out.append("")
    out.append("Diagnosis")
    out.append("-" * 60)

    qt = probe_qt()
    if any(s.ok for s in qt):
        binding = next(s for s in qt if s.ok)
        out.append(f"  OK: Qt bindings available ({binding.name} {binding.version}).")
    else:
        out.append("  PROBLEM: no usable Qt binding — the GUI cannot start.")
        out.append("")
        out.extend("  " + line for line in format_qt_import_failure(qt).splitlines())

    for status in problems:
        if status.status != "broken":
            continue
        hint = explain_import_error(status.error) if status.error else None
        if hint:
            out.append(f"  {status.name}: {hint}")

    version_warning = _python_version_warning()
    if version_warning:
        out.append(f"  WARNING: {version_warning}")

    mismatch = _launcher_mismatch()
    if mismatch:
        out.append("")
        out.extend("  " + line for line in mismatch.splitlines())

    if not problems and any(s.ok for s in qt) and not mismatch:
        out.append("  No problems detected.")

    out.append("")
    out.append(f"Log files: {os.path.dirname(_log_path())}")
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point for ``pyirena-doctor``.

    Prints a dependency and environment report. Users can paste the output
    into a bug report; maintainers can ask for it instead of guessing.

    Returns:
        Process exit code: 0 when the GUI stack looks usable, 1 otherwise.
    """
    del argv  # no options yet; accepted for symmetry with other entry points
    lines = _report()
    print("\n".join(lines))
    return 0 if any(s.ok for s in probe_qt()) else 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
