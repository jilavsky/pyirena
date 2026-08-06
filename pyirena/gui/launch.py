#!/usr/bin/env python
"""
Launcher script for pyIrena GUI.

Usage:
    python -m pyirena.gui.launch
    or
    pyirena-gui (if installed)
"""

import logging
import sys


def main():
    """Launch the pyIrena data selector GUI."""
    from pyirena.logging_setup import setup_logging, install_excepthook
    setup_logging("gui")
    install_excepthook()
    try:
        from pyirena.gui.data_selector import main as gui_main
        gui_main()
    except ImportError as e:
        # An ImportError here is not automatically a missing dependency: it
        # may be an installed-but-unloadable binary, or a bug in pyirena.
        # pyirena.diagnostics identifies which, and the traceback always goes
        # to the log so a user report is actionable.
        from pyirena.diagnostics import format_gui_import_failure

        # DEBUG, not ERROR: the file handler records DEBUG and above, so the
        # traceback is preserved for support, while the console shows only the
        # formatted diagnosis instead of a wall of stack frames.
        logging.getLogger("pyirena").debug(
            "GUI failed to start", exc_info=True
        )
        print(format_gui_import_failure(e), file=sys.stderr)
        if _verbose():
            raise
        print(
            "\n(Re-run with PYIRENA_DEBUG=1 to see the full traceback.)",
            file=sys.stderr,
        )
        sys.exit(1)


def _verbose() -> bool:
    """True when the user asked for raw tracebacks via PYIRENA_DEBUG."""
    import os

    return os.environ.get("PYIRENA_DEBUG", "").strip().lower() not in (
        "", "0", "false", "no",
    )


if __name__ == "__main__":
    main()
