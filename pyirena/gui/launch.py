#!/usr/bin/env python
"""
Launcher script for pyIrena GUI.

Usage:
    python -m pyirena.gui.launch
    or
    pyirena-gui (if installed)
"""

import sys


def main():
    """Launch the pyIrena data selector GUI."""
    try:
        from pyirena.gui.data_selector import main as gui_main
        gui_main()
    except ImportError as e:
        print("Error: GUI dependencies not installed.")
        print("Install with: pip install pyirena[gui]")
        print(f"\nDetails: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
