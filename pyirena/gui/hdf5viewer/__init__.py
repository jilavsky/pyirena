"""
HDF5 Viewer — standalone tool for browsing, plotting, and extracting data
from HDF5/NXcanSAS files produced by pyirena and compatible software.

Can be launched:
  - Standalone:  pyirena-viewer  (entry point in pyproject.toml)
  - From Data Selector: HDF5ViewerWindow(initial_folder=...)
"""

from .main_window import HDF5ViewerWindow

__all__ = ["HDF5ViewerWindow", "main"]


def main(initial_folder: str | None = None) -> None:
    """Entry point for standalone launch via ``pyirena-viewer``."""
    import sys

    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        from PyQt6.QtWidgets import QApplication  # type: ignore[no-redef]

    app = QApplication.instance() or QApplication(sys.argv)
    app.setStyle("Fusion")

    window = HDF5ViewerWindow(initial_folder=initial_folder)
    window.show()

    # Only call exec() when we created the QApplication (standalone mode).
    if not hasattr(app, "_hdf5viewer_embedded"):
        sys.exit(app.exec())
