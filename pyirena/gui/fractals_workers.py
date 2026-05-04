"""
QThread workers for the Fractals tool.

Two workers:
  * GrowthQueueWorker  — single background thread that consumes a job queue.
                         Jobs may be a single growth, a batch (Grow Many),
                         or an optimizer run.  The UI stays responsive and
                         multiple aggregates can be enqueued while one is
                         running.
  * MCIntensityWorker  — single-shot Monte-Carlo I(Q) computation for the
                         currently active aggregate.

Both workers expose Qt signals for completion / progress / failure that the
panel binds to via Qt.QueuedConnection (the default for cross-thread signals).
"""

from __future__ import annotations

import queue
import time
import traceback
import uuid
from dataclasses import dataclass
from typing import Optional

import numpy as np

try:
    from PySide6.QtCore import QThread, Signal
except ImportError:
    try:
        from PyQt6.QtCore import QThread, pyqtSignal as Signal
    except ImportError:
        from PyQt5.QtCore import QThread, pyqtSignal as Signal

from pyirena.core.fractals import (
    FractalAggregate, GrowthConfig, OptimizerConfig,
    grow_aggregate, optimize_growth, intensity_unified, intensity_montecarlo,
)


# ---------------------------------------------------------------------------
# Job descriptors
# ---------------------------------------------------------------------------

@dataclass
class _GrowJob:
    job_id: str
    label: str
    config: GrowthConfig


@dataclass
class _OptimizeJob:
    job_id: str
    label: str
    config: OptimizerConfig


class _ShutdownToken:
    """Sentinel pushed onto the queue to signal worker exit."""


# ---------------------------------------------------------------------------
# Growth queue worker
# ---------------------------------------------------------------------------

class GrowthQueueWorker(QThread):
    """Consumes growth/optimize jobs sequentially in a background thread."""

    # job_id, human-readable label
    job_started = Signal(str, str)
    # job_id, percent (0-100), short message
    job_progress = Signal(str, int, str)
    # job_id, FractalAggregate (one stick per Grow / per optimizer trial)
    aggregate_complete = Signal(str, object)
    # job_id, best_aggregate, list_of_all_attempted
    optimizer_complete = Signal(str, object, object)
    # job_id, error message
    job_failed = Signal(str, str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._queue: queue.Queue = queue.Queue()
        self._cancel_flags: dict[str, bool] = {}
        self._current_job_id: Optional[str] = None
        self._stop = False

    # ── Public API ───────────────────────────────────────────────────────

    def enqueue_grow(self, config: GrowthConfig, label: str = "") -> str:
        job_id = uuid.uuid4().hex[:8]
        self._cancel_flags[job_id] = False
        self._queue.put(_GrowJob(job_id=job_id, label=label or f"Grow Z={config.z}",
                                  config=config))
        return job_id

    def enqueue_grow_many(self, base_config: GrowthConfig, n: int) -> list[str]:
        job_ids = []
        for i in range(int(n)):
            cfg = GrowthConfig(
                z=base_config.z,
                sticking_prob=base_config.sticking_prob,
                num_test_paths=base_config.num_test_paths,
                rg_primary=base_config.rg_primary,
                allowed_near_dist=base_config.allowed_near_dist,
                attraction=base_config.attraction,
                # Vary seed deterministically when user pinned a non-zero seed,
                # else stay at 0 (fully random).
                seed=base_config.seed + i if base_config.seed != 0 else 0,
            )
            job_ids.append(self.enqueue_grow(cfg, label=f"Grow {i+1}/{n}"))
        return job_ids

    def enqueue_optimize(self, opt_config: OptimizerConfig) -> str:
        job_id = uuid.uuid4().hex[:8]
        self._cancel_flags[job_id] = False
        label = (f"Optimize → dmin={opt_config.target_dmin:.2f}, "
                 f"c={opt_config.target_c:.2f}")
        self._queue.put(_OptimizeJob(job_id=job_id, label=label, config=opt_config))
        return job_id

    def request_cancel(self, job_id: str) -> None:
        if job_id in self._cancel_flags:
            self._cancel_flags[job_id] = True

    def shutdown(self) -> None:
        """Signal the worker to drain the current job and exit."""
        self._stop = True
        self._queue.put(_ShutdownToken())
        # Cancel everything still pending
        for k in list(self._cancel_flags.keys()):
            self._cancel_flags[k] = True
        if self.isRunning():
            self.wait(5000)

    # ── Thread main loop ─────────────────────────────────────────────────

    def run(self) -> None:
        while not self._stop:
            try:
                job = self._queue.get(timeout=0.25)
            except queue.Empty:
                continue
            if isinstance(job, _ShutdownToken):
                break
            if isinstance(job, _GrowJob):
                self._run_grow(job)
            elif isinstance(job, _OptimizeJob):
                self._run_optimize(job)

    # ── Job handlers ─────────────────────────────────────────────────────

    def _check_cancel(self, job_id: str) -> bool:
        return self._cancel_flags.get(job_id, False) or self._stop

    def _run_grow(self, job: _GrowJob) -> None:
        self._current_job_id = job.job_id
        self.job_started.emit(job.job_id, job.label)
        try:
            t0 = time.perf_counter()

            def _progress(current, total):
                pct = int(100 * current / max(1, total))
                self.job_progress.emit(job.job_id, pct,
                                        f"{current}/{total}")
            agg = grow_aggregate(
                job.config,
                progress_cb=_progress,
                cancel_check=lambda: self._check_cancel(job.job_id),
            )
            agg.label = job.label
            elapsed = time.perf_counter() - t0
            self.job_progress.emit(job.job_id, 100,
                                    f"done in {elapsed:.1f}s")
            self.aggregate_complete.emit(job.job_id, agg)
        except Exception as exc:
            traceback.print_exc()
            self.job_failed.emit(job.job_id, str(exc))
        finally:
            self._cancel_flags.pop(job.job_id, None)
            self._current_job_id = None

    def _run_optimize(self, job: _OptimizeJob) -> None:
        self._current_job_id = job.job_id
        self.job_started.emit(job.job_id, job.label)
        try:
            t0 = time.perf_counter()
            iter_count = [0]

            def _on_iter(it: int, agg: FractalAggregate, obj: float) -> None:
                iter_count[0] += 1
                pct = int(100 * iter_count[0] / max(1, 3 * job.config.max_iter))
                self.job_progress.emit(job.job_id, pct,
                                        f"iter {it}, obj={obj:.4f}")
                self.aggregate_complete.emit(job.job_id, agg)

            best, all_aggs = optimize_growth(
                job.config,
                on_iter_complete=_on_iter,
                cancel_check=lambda: self._check_cancel(job.job_id),
            )
            elapsed = time.perf_counter() - t0
            self.job_progress.emit(job.job_id, 100,
                                    f"best in {elapsed:.1f}s")
            self.optimizer_complete.emit(job.job_id, best, all_aggs)
        except Exception as exc:
            traceback.print_exc()
            self.job_failed.emit(job.job_id, str(exc))
        finally:
            self._cancel_flags.pop(job.job_id, None)
            self._current_job_id = None


# ---------------------------------------------------------------------------
# Monte-Carlo intensity worker
# ---------------------------------------------------------------------------

class MCIntensityWorker(QThread):
    """Single-shot MC scattering computation for a given aggregate."""

    progress = Signal(int)             # percent 0-100
    finished_ok = Signal(str, object, object)   # agg_uuid, q, I_mc
    failed = Signal(str, str)          # agg_uuid, error message

    def __init__(self, parent=None):
        super().__init__(parent)
        self._aggregate: Optional[FractalAggregate] = None
        self._q: Optional[np.ndarray] = None
        self._cancel = False
        self._oversample = 10
        self._sphere_voxel_radius = 10
        self._max_pairs = 10_000_000
        self._time_budget_s = 20.0

    def configure(self, aggregate: FractalAggregate, q: np.ndarray,
                  *, oversample: int = 10, sphere_voxel_radius: int = 10,
                  max_pairs: int = 10_000_000, time_budget_s: float = 20.0) -> None:
        self._aggregate = aggregate
        self._q = np.asarray(q, dtype=np.float64)
        self._oversample = int(oversample)
        self._sphere_voxel_radius = int(sphere_voxel_radius)
        self._max_pairs = int(max_pairs)
        self._time_budget_s = float(time_budget_s)
        self._cancel = False

    def request_cancel(self) -> None:
        self._cancel = True

    def shutdown(self) -> None:
        self._cancel = True
        if self.isRunning():
            self.wait(5000)

    def run(self) -> None:
        if self._aggregate is None or self._q is None:
            self.failed.emit("", "MC worker not configured.")
            return
        agg_uuid = self._aggregate.uuid
        try:
            I_mc = intensity_montecarlo(
                self._aggregate, self._q,
                oversample=self._oversample,
                sphere_voxel_radius=self._sphere_voxel_radius,
                max_pairs=self._max_pairs,
                time_budget_s=self._time_budget_s,
                progress_cb=self.progress.emit,
                cancel_check=lambda: self._cancel,
            )
            self.finished_ok.emit(agg_uuid, self._q, I_mc)
        except Exception as exc:
            traceback.print_exc()
            self.failed.emit(agg_uuid, str(exc))
