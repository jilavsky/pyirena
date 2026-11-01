"""
Tests for the user-zoom preservation in :mod:`pyirena.gui.sas_plot`.

A user who zooms an I(Q) plot to inspect the fit in a small Q region must keep
that view when the model is recalculated (Graph model, autoupdate, fit done).
New data, the "A" button and "View All" give the full view back.

Skipped when no Qt binding is installed; the `test-gui` CI job runs them with
``QT_QPA_PLATFORM=offscreen``.
"""

import numpy as np
import pytest


def _require_qt() -> None:
    """Skip this module when no Qt binding is installed (see test_plot_export)."""
    try:
        import pyirena.gui._qt  # noqa: F401
    except ImportError:
        pytest.skip("Qt (PySide6/PyQt6) not available", allow_module_level=True)


_require_qt()

pytest.importorskip("pyqtgraph")

import pyqtgraph as pg  # noqa: E402

from pyirena.gui._qt import QApplication  # noqa: E402
from pyirena.gui.sas_plot import (  # noqa: E402
    make_sas_plot,
    plot_iq_data,
    set_robust_y_range,
    user_view_active,
)

# Module-level reference: a QApplication collected while pyqtgraph scenes are
# still alive segfaults the interpreter on exit.
_APP = None

ZOOM = ([-2.0, -1.5], [-1.0, 1.0])     # log10 view coordinates


@pytest.fixture(scope="session")
def qapp():
    global _APP
    _APP = QApplication.instance() or QApplication([])
    yield _APP


@pytest.fixture
def layout(qapp):
    widget = pg.GraphicsLayoutWidget()
    yield widget
    widget.clear()
    widget.close()
    widget.deleteLater()
    qapp.processEvents()


@pytest.fixture
def data():
    q = np.logspace(-3, -0.5, 200)
    I = 1e3 * np.exp(-((q * 100) ** 2) / 3) + 1e-3 * q**-4 + 0.01
    return q, I, 0.05 * I


def _user_zoom(plot):
    """Do what a mouse wheel / drag does: set a range, emit the manual signal."""
    vb = plot.getViewBox()
    vb.setRange(xRange=ZOOM[0], yRange=ZOOM[1], padding=0)
    vb.sigRangeChangedManually.emit(vb.state["mouseEnabled"])


def _view(plot):
    return [list(np.round(r, 6)) for r in plot.getViewBox().viewRange()]


def _redraw(plot, q, I, dI):
    """A model recalculation as the fitting panels do it."""
    plot.clear()
    plot_iq_data(plot, q, I, dI)
    set_robust_y_range(plot, I)
    plot.plot(q, I * 1.01)


def test_redraw_of_same_data_keeps_user_zoom(layout, data):
    plot = make_sas_plot(layout, 0, 0)
    plot_iq_data(plot, *data)
    _user_zoom(plot)
    _redraw(plot, *data)
    assert user_view_active(plot)
    assert _view(plot) == [ZOOM[0], ZOOM[1]]


def test_without_user_zoom_redraw_still_autoranges(layout, data):
    plot = make_sas_plot(layout, 0, 0)
    plot_iq_data(plot, *data)
    full = _view(plot)
    _redraw(plot, *data)
    assert not user_view_active(plot)
    assert _view(plot) == full


def test_new_data_resets_user_zoom(layout, data):
    q, I, dI = data
    plot = make_sas_plot(layout, 0, 0)
    plot_iq_data(plot, q, I, dI)
    _user_zoom(plot)
    _redraw(plot, q, 3 * I, dI)
    assert not user_view_active(plot)
    assert _view(plot) != [ZOOM[0], ZOOM[1]]


def test_view_all_hands_axes_back(layout, data):
    plot = make_sas_plot(layout, 0, 0)
    plot_iq_data(plot, *data)
    _user_zoom(plot)
    plot.getViewBox().menu.viewAll.trigger()
    assert not user_view_active(plot)


def test_zoom_in_x_linked_plot_protects_shared_axis(layout, data):
    iq = make_sas_plot(layout, 0, 0)
    resid = make_sas_plot(layout, 1, 0, log_y=False, x_link=iq)
    plot_iq_data(iq, *data)
    _user_zoom(resid)
    assert user_view_active(iq)
