"""Saved/exported panel state carries the slit keys a replay needs (issue #38).

Without ``load_slit_smeared`` a batch / analyze replay loaded the file's
desmeared copy and smeared the model against it; Simple Fits and Modeling
wrote no slit keys at all.
"""
import os

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app(tmp_path_factory):
    os.environ["HOME"] = str(tmp_path_factory.mktemp("home"))
    from pyirena.gui._qt import QApplication
    return QApplication.instance() or QApplication([])


@pytest.mark.parametrize("module, cls", [
    ("pyirena.gui.unified_fit", "UnifiedFitPanel"),
    ("pyirena.gui.sizes_panel", "SizesFitPanel"),
    ("pyirena.gui.simple_fits_panel", "SimpleFitsPanel"),
])
def test_state_has_slit_keys(app, module, cls):
    import importlib
    panel = getattr(importlib.import_module(module), cls)()
    q = np.logspace(-3, -1, 50)
    panel.set_data(q, 1.0 / q, 0.05 / q, "s", slit_length=0.02, is_slit_smeared=True)
    state = panel._collect_state()
    assert state["use_slit_smearing"] is True
    assert state["slit_length"] == pytest.approx(0.02)
    assert state["load_slit_smeared"] is False      # no SMR copy in a file here
    panel.close()


def test_replayed_simple_config_with_smearing_builds_a_smeared_model():
    from pyirena.core.tool_config import build_setup
    setup = build_setup({"simple_fits": {"model": "Guinier", "use_slit_smearing": True,
                                         "slit_length": 0.03}})
    assert setup.model.use_slit_smearing and setup.model.slit_length == pytest.approx(0.03)
