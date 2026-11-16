"""Unified Fit panel fields that depend on the values must follow a fit.

Two bugs of the same kind, found in the documentation screenshots
(2026-11-16):

1. The "Level feasible / not feasible" banner follows the values shown.

The banner was recalculated only when the user edited a field
(``editingFinished``) or a setup was loaded. A fit, a local Guinier/Porod fit,
Revert back or a level copy writes the new values through
``LevelParametersWidget.set_parameters()`` with ``setText()``, which emits no
signal, so the banner kept describing the previous values: after a good fit
from the coded defaults it still said "Level not feasible".

2. A linked RgCutoff (Link RgCutoff: RgCutoff of level N = Rg of level N-1)
was synced before a fit but not after it. The fit itself used the linked
value, but the field kept the pre-fit Rg, and "Results to graphs" and
"Store in File" (which read the fields) showed and saved the stale cutoff.
"""

from __future__ import annotations

import time
from pathlib import Path

import pytest


def _require_qt() -> None:
    try:
        import pyirena.gui._qt  # noqa: F401
    except ImportError:
        pytest.skip("Qt (PySide6/PyQt6) not available", allow_module_level=True)


_require_qt()
pytest.importorskip("pyqtgraph")

from pyirena.gui._qt import QApplication  # noqa: E402
from pyirena.gui.unified_fit import LevelParametersWidget, UnifiedFitPanel  # noqa: E402

DATA = Path(__file__).resolve().parents[2] / "validationData" / "unified" / "unified_two_level.h5"

_APP = None


@pytest.fixture(scope="session")
def qapp():
    global _APP
    _APP = QApplication.instance() or QApplication([])
    yield _APP


def test_set_parameters_refreshes_the_feasibility_banner(qapp):
    lw = LevelParametersWidget(1)
    # The coded defaults (G 100, Rg 100, B 0.01, P 4) are not a feasible level.
    lw.update_feasibility_status()
    assert lw.feasibility_label.text() == "Level not feasible"

    # Level 1 of validationData/unified/unified_two_level.h5: feasible.
    lw.set_parameters({"G": 8.0, "Rg": 120.0, "B": 1.88e-7, "P": 4.0})
    assert lw.feasibility_label.text().startswith("Level feasible")

    # And back: a B far above the Guinier-consistent value.
    lw.set_parameters({"B": 1e-3})
    assert lw.feasibility_label.text() == "Level not feasible"


@pytest.mark.skipif(not DATA.is_file(), reason="validationData not present")
def test_linked_rgcutoff_follows_the_fitted_rg(qapp, tmp_path, monkeypatch):
    from pyirena.io.hdf5 import readGenericNXcanSAS

    # The panel saves its state on close; keep that out of the session-wide
    # settings folder, which later panel tests load their defaults from.
    monkeypatch.setenv("PYIRENA_STATE_DIR", str(tmp_path))

    d = readGenericNXcanSAS(str(DATA.parent), DATA.name)
    panel = UnifiedFitPanel()
    try:
        panel.set_data(d["Q"], d["Intensity"], d.get("Error"), label=DATA.name)
        panel.num_levels_spin.setValue(2)
        start = {1: {"G": 10, "Rg": 100, "B": 1e-7, "P": 4},
                 2: {"G": 3000, "Rg": 1000, "B": 1e-6, "P": 3}}
        for level, values in start.items():
            lw = panel.level_widgets[level - 1]
            lw.set_parameters({**values, "fit_G": True, "fit_Rg": True, "fit_B": True,
                               "fit_P": level == 2})
            lw.fix_limits()
        panel.level_widgets[1].link_rgco_check.setChecked(True)
        panel.fit_background_check.setChecked(True)

        panel.run_fit()
        t0 = time.monotonic()
        while panel._jobs.busy and time.monotonic() - t0 < 120:
            qapp.processEvents()
            time.sleep(0.01)
        qapp.processEvents()

        rg1 = float(panel.level_widgets[0].rg_value.text())
        assert rg1 == pytest.approx(120, rel=0.02)          # the fit moved Rg_1 off 100
        assert float(panel.level_widgets[1].rg_cutoff.text()) == pytest.approx(rg1, rel=1e-3)
        stored = panel._collect_level_dicts()[1]["RgCutoff"]  # what Store in File writes
        assert stored == pytest.approx(rg1, rel=1e-3)
        assert panel.level_widgets[0].feasibility_label.text().startswith("Level feasible")
    finally:
        panel.close()
