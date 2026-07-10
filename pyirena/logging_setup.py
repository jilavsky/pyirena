"""
Central logging configuration for pyirena.

Log files live in ``~/.pyirena/logs/`` (next to ``state.json``), one file
per entry point (``gui.log``, ``mcp.log``, ``batch.log``), rotated at 2 MB
with 5 backups — a 10 MB hard cap per tool.  The file handler records
DEBUG and above; the console handler shows INFO and above (plain messages,
so console output looks the same as the old ``print()`` style).

Design rules
------------
- Importing pyirena NEVER configures logging (standard library-package
  behavior).  Only entry points call :func:`setup_logging`.
- Scripts that want file logs can call it themselves::

      import pyirena
      pyirena.setup_logging("batch")

- The batch API calls :func:`ensure_console_output` so progress messages
  remain visible to scripts that never configure logging — without
  writing any files behind the user's back.
- :func:`install_excepthook` routes uncaught exceptions into the log,
  so GUI crashes are captured even when the user never saw a traceback.
"""
from __future__ import annotations

import logging
import logging.handlers
import platform
import sys
from pathlib import Path

_FILE_FORMAT = "%(asctime)s %(levelname)-7s %(name)s %(funcName)s: %(message)s"
_CONSOLE_FORMAT = "%(message)s"

MAX_BYTES = 2_000_000   # rotate each log file at ~2 MB
BACKUP_COUNT = 5        # keep 5 rotated files -> 10 MB cap per tool

_configured_tools: set = set()
_console_attached = False


def get_log_dir() -> Path:
    """Return (and create) the pyirena log directory: ~/.pyirena/logs/."""
    log_dir = Path.home() / ".pyirena" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir


def get_logger(name: str = "") -> logging.Logger:
    """Return a logger under the ``pyirena`` hierarchy.

    ``get_logger("io.hdf5")`` -> logger ``pyirena.io.hdf5``.  Modules may
    equally use ``logging.getLogger(__name__)`` (identical result for
    modules inside the package).
    """
    return logging.getLogger(f"pyirena.{name}" if name else "pyirena")


def setup_logging(
    tool: str = "pyirena",
    file_level: int = logging.DEBUG,
    console_level: int = logging.INFO,
    console: bool = True,
) -> logging.Logger:
    """
    Configure pyirena logging for an entry point.  Idempotent per *tool*.

    Args:
        tool: Log file stem — ``gui``, ``mcp``, ``batch``, ...
              (-> ``~/.pyirena/logs/<tool>.log``).
        file_level: Level captured in the rotating file (default DEBUG —
              rotation caps the cost, and post-mortem detail is the point).
        console_level: Level echoed to the console (default INFO).
        console: Attach a console (stderr) handler.  Safe for the MCP
              server: stdout carries the protocol, stderr does not.

    Returns:
        The package root logger ``pyirena``.
    """
    global _console_attached
    logger = logging.getLogger("pyirena")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    if tool in _configured_tools:
        return logger
    _configured_tools.add(tool)

    try:
        file_handler = logging.handlers.RotatingFileHandler(
            get_log_dir() / f"{tool}.log",
            maxBytes=MAX_BYTES,
            backupCount=BACKUP_COUNT,
            encoding="utf-8",
            delay=True,
        )
        file_handler.setLevel(file_level)
        file_handler.setFormatter(logging.Formatter(_FILE_FORMAT))
        logger.addHandler(file_handler)
        # Route warnings.warn() output into the same file
        logging.captureWarnings(True)
        logging.getLogger("py.warnings").addHandler(file_handler)
    except OSError:
        # Home dir not writable (sandboxes, read-only accounts): carry on
        # console-only rather than refusing to start.
        file_handler = None

    if console and not _console_attached:
        console_handler = logging.StreamHandler()  # stderr
        console_handler.setLevel(console_level)
        console_handler.setFormatter(logging.Formatter(_CONSOLE_FORMAT))
        logger.addHandler(console_handler)
        _console_attached = True

    # Session header — file only (DEBUG), keeps the console clean
    try:
        from pyirena import __version__
    except Exception:
        __version__ = "unknown"
    logger.debug(
        "=== pyirena %s session start | tool=%s | python %s | %s %s ===",
        __version__, tool, platform.python_version(),
        platform.system(), platform.release(),
    )
    return logger


def ensure_console_output() -> None:
    """
    Attach a plain-message INFO console handler iff logging is completely
    unconfigured (no pyirena handlers, no root handlers).

    Called by the batch API so its progress messages stay visible to
    scripts exactly like the old print() output — while users who did
    configure logging keep full control, and no files are written.
    """
    global _console_attached
    logger = logging.getLogger("pyirena")
    if logger.handlers or logging.getLogger().handlers:
        return
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter(_CONSOLE_FORMAT))
    logger.addHandler(console_handler)
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    _console_attached = True


def install_excepthook() -> None:
    """
    Log uncaught exceptions (CRITICAL, full traceback) before the default
    handler runs.  In the GUI this is the difference between a user report
    of "it crashed" and an actionable traceback in ~/.pyirena/logs/gui.log.
    """
    previous_hook = sys.excepthook

    def _hook(exc_type, exc_value, exc_tb):
        if not issubclass(exc_type, KeyboardInterrupt):
            logging.getLogger("pyirena").critical(
                "Uncaught exception", exc_info=(exc_type, exc_value, exc_tb)
            )
        previous_hook(exc_type, exc_value, exc_tb)

    sys.excepthook = _hook
