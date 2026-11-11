"""
pyirena.gui.fit_worker — run fits and Monte-Carlo off the GUI thread.

Until 1.2.1 only Modeling did this. Unified Fit, Size Distribution, Simple
Fits, the Carbon model and WAXS Peak Fit called the core fit on the GUI
thread, so the window froze for the length of a fit; Unified and WAXS kept it
painting with ``QApplication.processEvents()``, which also let a second click
start a second fit inside the first. Monte-Carlo uncertainty — the longest
operation in those tools — had no progress and no way to stop it (review of
09-11-2026, issue #38).

One shared runner now gives every tool the same behaviour:

* the computation runs in a :class:`FitWorker` (a ``QThread``);
* while it runs, the panel's action buttons are disabled, so nothing can be
  started twice; the button that started a *cancellable* job turns into a
  red "Cancel" button;
* progress and the result come back to the GUI thread as queued signals, so
  every widget update stays on the GUI thread.

**What can be cancelled.** Monte-Carlo runs: the core
``calculate_uncertainty_mc`` functions poll a cancel callback between refits
and return what finished. A single fit cannot be interrupted — scipy's
optimisers take no cancel hook — but fits take seconds, and the window stays
responsive while they run.

**The job contract.** A job is a plain callable ``job(progress_cb, cancel_cb)``
returning the result object. It runs on the worker thread, so it must not
touch a widget, and should work on copies of anything the GUI may redraw
from meanwhile (the panel's model, its data arrays). Collect inputs on the
GUI thread *before* starting; apply the result in ``on_done``, which runs on
the GUI thread.

Usage in a panel::

    self._jobs = FitJobRunner(self, lock=[self.fit_button, self.mc_button])
    ...
    self._jobs.start(job, on_done=self._finish_mc, button=self.mc_button,
                     cancellable=True, on_progress=self._show_mc_progress)
    # a second click on mc_button while it runs:  self._jobs.cancel()

``wait()`` blocks until the job is done and its callbacks have run — for
tests and for code that genuinely needs the result before going on.
"""

from __future__ import annotations

import logging
import time
from typing import Callable, Optional, Sequence

from pyirena.gui._qt import QApplication, QObject, QThread, Signal, Slot

log = logging.getLogger(__name__)

#: Red "Cancel" look for the button of a running, cancellable job.
_CANCEL_STYLE = (
    "QPushButton {background: #e74c3c; color: white; font-weight: bold;"
    " border-radius: 4px;}"
    "QPushButton:hover {background: #c0392b;}"
    "QPushButton:disabled {background: #bdc3c7;}"
)


class FitWorker(QThread):
    """Runs one ``job(progress_cb, cancel_cb)`` on a worker thread.

    Signals carry the outcome back; connect them to methods of a QObject that
    lives on the GUI thread (as :class:`FitJobRunner` does) so the slots run
    there, never here.
    """

    progress = Signal(object, object)   # (done, total) — or any pair a job reports
    succeeded = Signal(object)      # the job's return value
    failed = Signal(str)            # message

    def __init__(self, job: Callable, parent=None):
        super().__init__(parent)
        self._job = job
        self._cancelled = False

    def cancel(self) -> None:
        """Ask the job to stop at its next cancel check."""
        self._cancelled = True

    @property
    def cancelled(self) -> bool:
        return self._cancelled

    def run(self) -> None:  # noqa: D401 - QThread API
        try:
            result = self._job(lambda a, b: self.progress.emit(a, b),
                               lambda: self._cancelled)
        except Exception as exc:
            log.exception("Background fit job failed")
            self.failed.emit(str(exc) or exc.__class__.__name__)
            return
        self.succeeded.emit(result)


class FitJobRunner(QObject):
    """One background job at a time for a panel, with button handling.

    Args:
        parent: The panel (a QWidget on the GUI thread).
        lock: Buttons/widgets to disable while any job runs.
    """

    def __init__(self, parent, lock: Sequence = ()):
        super().__init__(parent)
        self._lock = list(lock)
        self._worker: Optional[FitWorker] = None
        self._button = None
        self._button_text = ""
        self._button_style = ""
        self._on_done = None
        self._on_error = None
        self._on_progress = None
        self._pending = False           # callbacks not yet delivered

    # ── state ────────────────────────────────────────────────────────────

    @property
    def busy(self) -> bool:
        """True from start() until the job's callback has run."""
        return self._pending

    @property
    def cancelled(self) -> bool:
        """Whether the current (or last) job was asked to cancel."""
        return bool(self._worker is not None and self._worker.cancelled)

    def add_lock(self, *widgets) -> None:
        """Register more widgets to disable while a job runs."""
        self._lock.extend(w for w in widgets if w is not None)

    # ── control ──────────────────────────────────────────────────────────

    def start(self, job: Callable, on_done: Callable, *, on_error: Optional[Callable] = None,
              on_progress: Optional[Callable] = None, button=None,
              cancellable: bool = False, busy_text: Optional[str] = None) -> bool:
        """Start *job* in the background. Returns False if a job is already running.

        Args:
            job: ``job(progress_cb, cancel_cb) -> result``; runs off the GUI thread.
            on_done: ``on_done(result)`` on the GUI thread.
            on_error: ``on_error(message)`` on the GUI thread (default: log it).
            on_progress: ``on_progress(a, b)`` on the GUI thread, with whatever
                pair the job passes to its ``progress_cb`` — usually
                ``(done, total)``.
            button: The button that started the job.
            cancellable: Turn *button* into a Cancel button while the job runs.
            busy_text: Label for *button* while running (default "Cancel" when
                cancellable, otherwise the label is kept and the button disabled).
        """
        if self._pending:
            return False
        self._on_done, self._on_error, self._on_progress = on_done, on_error, on_progress
        self._button = button
        self._set_locked(True)
        if button is not None:
            self._button_text = button.text()
            self._button_style = button.styleSheet()
            if cancellable:
                button.setText(busy_text or "Cancel")
                button.setStyleSheet(_CANCEL_STYLE)
                button.setEnabled(True)
            else:
                if busy_text:
                    button.setText(busy_text)
                button.setEnabled(False)

        worker = FitWorker(job, parent=self)
        worker.progress.connect(self._deliver_progress)
        worker.succeeded.connect(self._deliver_success)
        worker.failed.connect(self._deliver_failure)
        self._worker = worker
        self._pending = True
        worker.start()
        return True

    def cancel(self) -> None:
        """Ask the running job to stop (Monte-Carlo keeps the passes it finished)."""
        if self._worker is not None and self._pending:
            self._worker.cancel()
            if self._button is not None:
                self._button.setText("Cancelling…")
                self._button.setEnabled(False)

    def wait(self, timeout_s: float = 600.0) -> bool:
        """Block until the current job finished and its callback ran.

        Returns True when it did, False on timeout.
        """
        deadline = time.monotonic() + timeout_s
        while self._pending:
            if self._worker is not None:
                self._worker.wait(50)
            QApplication.processEvents()
            if time.monotonic() > deadline:
                return False
        return True

    def shutdown(self, timeout_s: float = 30.0) -> None:
        """Cancel and join the worker — call from the panel's ``closeEvent``.

        A QThread destroyed while running aborts the process, so a window
        closed mid-fit must wait for it. Callbacks are dropped: the panel is
        going away.
        """
        worker = self._worker
        if worker is None or not worker.isRunning():
            return
        worker.cancel()
        self._on_done = self._on_error = self._on_progress = None
        if not worker.wait(int(timeout_s * 1000)):
            log.warning("Background fit still running after %.0f s at window close", timeout_s)

    # ── delivery (GUI thread) ────────────────────────────────────────────

    def _set_locked(self, locked: bool) -> None:
        for w in self._lock:
            try:
                w.setEnabled(not locked)
            except RuntimeError:        # widget already deleted
                pass

    def _release(self) -> None:
        if self._button is not None:
            try:
                self._button.setText(self._button_text)
                self._button.setStyleSheet(self._button_style)
                self._button.setEnabled(True)
            except RuntimeError:
                pass
        self._set_locked(False)
        self._pending = False

    @Slot(object, object)
    def _deliver_progress(self, done, total) -> None:
        if self._on_progress is not None:
            self._on_progress(done, total)

    @Slot(object)
    def _deliver_success(self, result) -> None:
        self._release()
        if self._on_done is not None:
            try:
                self._on_done(result)
            except Exception:
                log.exception("Applying a background fit result failed")

    @Slot(str)
    def _deliver_failure(self, message: str) -> None:
        self._release()
        if self._on_error is not None:
            self._on_error(message)
        else:
            log.error("Background fit job failed: %s", message)
