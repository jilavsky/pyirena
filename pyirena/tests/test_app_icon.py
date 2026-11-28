"""The pyIrena application icon ships inside the package and installs on the QApplication."""

from __future__ import annotations

import os

import pytest

# Skip cleanly when no Qt binding is installed (headless / core-only environments).
pytest.importorskip("pyirena.gui._qt")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from pyirena.gui._qt import QApplication, QIcon  # noqa: E402
from pyirena.gui.icon import _ICON_PATH, app_icon, set_app_icon  # noqa: E402


def test_icon_png_is_inside_the_package():
    assert _ICON_PATH.is_file()
    assert "pyirena/gui/assets" in _ICON_PATH.as_posix()


def test_icon_loads_and_is_installed_on_the_app():
    app = QApplication.instance() or QApplication([])
    icon = app_icon()
    assert isinstance(icon, QIcon)
    assert not icon.isNull()
    set_app_icon(app)
    assert not app.windowIcon().isNull()
