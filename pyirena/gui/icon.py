"""Application icon for pyIrena: Dock, taskbar, Cmd+Tab switcher and window title bars.

GUI layer. ``assets/app_icon.png`` is the 1024 px render of ``icons/pyIrena_icon.svg`` at
the repo root, one of the shared icon family for pyIrena, Matilda, Anisa, Bernardyn and AIDA
(the generator is ``pyirena/icons/build_icons.py``). The PNG ships in the wheel through
``[tool.setuptools.package-data]``. On macOS, Qt forwards ``QApplication.setWindowIcon`` to
the Dock tile, so no ``.app`` bundle is needed; ``icons/pyIrena.icns`` is there for one.
"""

from __future__ import annotations

from pathlib import Path

from pyirena.gui._qt import QIcon

_ICON_PATH = Path(__file__).resolve().parent / "assets" / "app_icon.png"


def app_icon() -> QIcon:
    """Return the pyIrena application icon (a null ``QIcon`` if the PNG is missing)."""
    return QIcon(str(_ICON_PATH))


def set_app_icon(app) -> None:
    """Make the pyIrena icon the default for every window of *app* (the ``QApplication``)."""
    app.setWindowIcon(app_icon())


__all__ = ["app_icon", "set_app_icon"]
