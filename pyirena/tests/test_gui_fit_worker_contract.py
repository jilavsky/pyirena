"""
Fits and Monte-Carlo run off the GUI thread, enforced (issue #38).

Until 1.2.1 only Modeling ran its fit in a QThread. The other five fitting
panels froze the window, and two kept it painting with
``QApplication.processEvents()`` — which also let a second click start a
second fit inside the first. They now share ``pyirena.gui.fit_worker``.

The source-scan half runs without Qt and fails the build if a panel goes back
to pumping the event loop or stops using the shared runner. The behaviour half
(skipped without Qt) checks the runner's thread contract: the job runs off the
GUI thread, every callback runs on it, buttons are locked and restored, and a
cancel reaches the job.
"""

from __future__ import annotations

import threading
import time
from pathlib import Path

import pytest

import pyirena

GUI = Path(pyirena.__file__).parent / "gui"

#: Panels whose fit / MC must go through FitJobRunner. Modeling keeps its own
#: two QThread workers (the pattern this was modelled on).
PANELS = [
    "unified_fit.py",
    "sizes_panel.py",
    "simple_fits_panel.py",
    "carbon_fit_panel.py",
    "waxs_peakfit_panel.py",
]


@pytest.mark.parametrize("name", PANELS)
def test_panel_uses_the_shared_runner(name):
    src = (GUI / name).read_text(encoding="utf-8")
    assert "FitJobRunner(" in src, f"{name} no longer creates a FitJobRunner"
    assert "self._jobs.shutdown()" in src, (
        f"{name}: closeEvent must join a running fit — a QThread destroyed "
        "while running aborts the process")


@pytest.mark.parametrize("name", PANELS)
def test_panel_does_not_pump_the_event_loop(name):
    src = (GUI / name).read_text(encoding="utf-8")
    assert "processEvents(" not in src, (
        f"{name} calls processEvents(): long work belongs in a FitJobRunner job, "
        "not on the GUI thread with the event loop pumped from inside it")


# ── behaviour (needs Qt) ─────────────────────────────────────────────────────

def _qt_app():
    pytest.importorskip("PySide6")
    import os

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from pyirena.gui._qt import QApplication

    return QApplication.instance() or QApplication([])


def test_runner_threads_locks_and_cancels():
    _qt_app()
    from pyirena.gui._qt import QPushButton, QWidget
    from pyirena.gui.fit_worker import FitJobRunner

    host = QWidget()
    go, other = QPushButton("Go", host), QPushButton("Other", host)
    runner = FitJobRunner(host, lock=[go, other])
    gui_thread = threading.get_ident()
    seen = {}

    def job(progress_cb, cancel_cb):
        seen["job"] = threading.get_ident()
        for i in range(200):
            if cancel_cb():
                return f"stopped at {i}"
            progress_cb(i + 1, 200)
            time.sleep(0.005)
        return "done"

    progress_threads = set()
    assert runner.start(job, lambda r: seen.update(result=r, done=threading.get_ident()),
                        on_progress=lambda d, n: progress_threads.add(threading.get_ident()),
                        button=go, cancellable=True)
    assert go.text() == "Cancel" and go.isEnabled() and not other.isEnabled()
    assert runner.busy
    assert runner.start(job, lambda r: None) is False      # one job at a time
    time.sleep(0.05)
    runner.cancel()
    assert runner.wait(10)
    assert seen["job"] != gui_thread
    assert seen["done"] == gui_thread
    assert progress_threads <= {gui_thread}
    assert seen["result"].startswith("stopped at")
    assert runner.cancelled
    assert go.text() == "Go" and go.isEnabled() and other.isEnabled()
    assert not runner.busy


def test_runner_reports_errors_on_gui_thread():
    _qt_app()
    from pyirena.gui._qt import QWidget
    from pyirena.gui.fit_worker import FitJobRunner

    host = QWidget()              # keep a reference: the runner is its child
    runner = FitJobRunner(host)
    got = {}

    def boom(progress_cb, cancel_cb):
        raise ValueError("no good")

    runner.start(boom, lambda r: None,
                 on_error=lambda m: got.update(msg=m, thread=threading.get_ident()))
    assert runner.wait(10)
    assert got["msg"] == "no good"
    assert got["thread"] == threading.get_ident()
    assert not runner.busy
