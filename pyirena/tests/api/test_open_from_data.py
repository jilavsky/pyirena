"""Tests for ``open_dataset_from_data`` — sessions built from arrays, no file.

This is the entry point the ZMQ service (``planning/zmq-service/``) needs: the
orchestrator and the analysis workstation share no filesystem, so the data
arrives inside the request. The same function is useful from any script that
holds arrays rather than a path.

Two properties matter beyond "it runs a fit":

* **Cleaning parity with files.** Arrays go through the same
  ``clean_sas_arrays`` rules as text-file import, so a curve behaves
  identically whichever way it arrived — and the counts are reported, not
  applied silently.
* **``file_path is None`` is handled everywhere.** Every reader of
  ``Session.file_path`` (summaries, ``save_*``) must cope; ``save_*`` without
  an ``output_path`` must fail with a readable error rather than a traceback.
"""
from __future__ import annotations

import numpy as np
import pytest

import pyirena.api.control as ctrl
from pyirena.api.control.session import get_session


def _curve(n: int = 200, seed: int = 0):
    """A Guinier + Porod curve with realistic dynamic range."""
    rng = np.random.default_rng(seed)
    q = np.logspace(-3, 0, n)
    I = 1000.0 * np.exp(-(q**2) * 150.0**2 / 3) + 1e-6 * q**-4.0 + 0.01
    I = I * (1 + 0.02 * rng.standard_normal(q.size))
    return q, I, 0.02 * I


@pytest.fixture
def arrays():
    return _curve()


def _open(q, I, **kw):
    return ctrl.open_dataset_from_data(q=list(q), intensity=list(I), **kw)


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------

def test_opens_a_session_with_no_file_behind_it(arrays):
    q, I, e = arrays
    r = _open(q, I, error=list(e), label="curve A")
    assert "error" not in r

    summary = r["summary"]
    assert summary["file"] is None
    assert summary["label"] == "curve A"
    assert summary["n_points"] == len(q)
    assert summary["has_errors"] is True
    assert summary["q_min"] == pytest.approx(q.min())
    assert summary["q_max"] == pytest.approx(q.max())

    s = get_session(r["session_id"])
    assert s.file_path is None
    assert len(s.q) == len(q)
    ctrl.close_session(r["session_id"])


def test_session_summary_and_listing_tolerate_a_null_file(arrays):
    q, I, _ = arrays
    sid = _open(q, I)["session_id"]

    assert ctrl.get_session_summary(sid)["file"] is None
    listed = [s for s in ctrl.list_open_sessions()["sessions"] if s["session_id"] == sid]
    assert listed and listed[0]["file"] is None
    ctrl.close_session(sid)


def test_errors_are_repaired_when_given_but_never_invented(arrays):
    """Omitting `error` means unweighted, as for a file with no error column."""
    q, I, e = arrays

    sid = _open(q, I)["session_id"]
    assert get_session(sid).error is None
    ctrl.close_session(sid)

    bad = e.copy()
    bad[:5] = 0.0
    r = _open(q, I, error=list(bad))
    assert r["summary"]["cleaning"]["n_repaired_error"] == 5
    s = get_session(r["session_id"])
    assert s.error is not None and np.all(s.error > 0)
    ctrl.close_session(r["session_id"])


def test_unsorted_input_is_sorted_and_reported(arrays):
    q, I, e = arrays
    order = np.argsort(-q)                       # descending
    r = _open(q[order], I[order], error=list(e[order]))

    assert r["summary"]["sorted_by_q"] is True
    s = get_session(r["session_id"])
    assert np.all(np.diff(s.q) > 0)
    # The pairing must survive the sort, not just the Q axis.
    assert s.intensity[0] == pytest.approx(I[np.argmin(q)])
    ctrl.close_session(r["session_id"])


def test_dq_is_stored_for_provenance(arrays):
    q, I, e = arrays
    dq = 0.1 * q
    r = _open(q, I, error=list(e), dq=list(dq))
    s = get_session(r["session_id"])
    assert s.dq is not None and len(s.dq) == len(s.q)
    ctrl.close_session(r["session_id"])


# ---------------------------------------------------------------------------
# Cleaning — same rules as text-file import
# ---------------------------------------------------------------------------

def test_nonpositive_and_nonfinite_points_are_removed_and_counted(arrays):
    q, I, e = arrays
    q = q.copy()
    I = I.copy()
    q[0] = 0.0            # Q <= 0     → removed
    q[1] = np.nan         # non-finite → removed
    I[10] = 0.0           # I <= 0     → removed
    I[11] = -5.0          # I <= 0     → removed

    r = _open(q, I, error=list(e))
    cleaning = r["summary"]["cleaning"]
    assert cleaning["n_input"] == len(q)
    assert cleaning["n_removed_q"] == 2
    assert cleaning["n_removed_i"] == 2
    assert cleaning["n_kept"] == len(q) - 4
    assert r["summary"]["n_points"] == len(q) - 4
    ctrl.close_session(r["session_id"])


# ---------------------------------------------------------------------------
# Validation — every rejection has its own code an agent can branch on
# ---------------------------------------------------------------------------

def test_mismatched_lengths(arrays):
    q, I, _ = arrays
    r = ctrl.open_dataset_from_data(q=list(q), intensity=list(I[:-3]))
    assert r["code"] == "SHAPE_MISMATCH"
    assert "intensity" in r["error"]


def test_mismatched_error_length(arrays):
    q, I, e = arrays
    r = ctrl.open_dataset_from_data(q=list(q), intensity=list(I), error=list(e[:10]))
    assert r["code"] == "SHAPE_MISMATCH"


def test_empty_arrays():
    assert ctrl.open_dataset_from_data(q=[], intensity=[])["code"] == "EMPTY_DATA"


def test_too_few_points():
    r = ctrl.open_dataset_from_data(q=[0.1, 0.2, 0.3], intensity=[1.0, 2.0, 3.0])
    assert r["code"] == "TOO_FEW_POINTS"


def test_too_many_points_is_configurable(monkeypatch, arrays):
    q, I, _ = arrays
    monkeypatch.setenv("PYIRENA_MAX_INPUT_POINTS", "50")
    r = _open(q, I)
    assert r["code"] == "TOO_MANY_POINTS"
    assert "50" in r["error"]

    monkeypatch.setenv("PYIRENA_MAX_INPUT_POINTS", "100000")
    r = _open(q, I)
    assert "error" not in r
    ctrl.close_session(r["session_id"])


def test_non_numeric_values():
    r = ctrl.open_dataset_from_data(q=[0.1, "x", 0.3, 0.4, 0.5],
                                    intensity=[1.0, 2.0, 3.0, 4.0, 5.0])
    assert r["code"] == "BAD_VALUES"


def test_nested_arrays_are_rejected():
    r = ctrl.open_dataset_from_data(q=[[0.1, 0.2], [0.3, 0.4]],
                                    intensity=[[1.0, 2.0], [3.0, 4.0]])
    assert r["code"] == "BAD_VALUES"


def test_everything_cleaned_away():
    r = ctrl.open_dataset_from_data(q=[-1.0, -2.0, -3.0, -4.0, -5.0, -6.0],
                                    intensity=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    assert r["code"] == "NO_VALID_POINTS"


# ---------------------------------------------------------------------------
# A full fit, and persistence without a source file
# ---------------------------------------------------------------------------

def test_fit_runs_on_an_in_memory_session(arrays):
    q, I, e = arrays
    sid = _open(q, I, error=list(e))["session_id"]

    assert "error" not in ctrl.select_model(sid, "unified_fit")
    assert "error" not in ctrl.add_unified_level(sid)
    fit = ctrl.run_fit(sid)
    assert "error" not in fit
    assert np.isfinite(fit["reduced_chi_squared"])
    ctrl.close_session(sid)


def test_save_without_output_path_is_refused_clearly(arrays):
    q, I, e = arrays
    sid = _open(q, I, error=list(e))["session_id"]
    ctrl.select_model(sid, "unified_fit")
    ctrl.add_unified_level(sid)
    ctrl.run_fit(sid)

    r = ctrl.save_fit(sid)
    assert r["code"] == "NO_SOURCE_FILE"
    assert "output_path" in r["suggestion"]
    ctrl.close_session(sid)


def test_save_with_output_path_writes_a_reopenable_file(tmp_path, arrays):
    q, I, e = arrays
    sid = _open(q, I, error=list(e), label="from_arrays")["session_id"]
    ctrl.select_model(sid, "unified_fit")
    ctrl.add_unified_level(sid)
    ctrl.run_fit(sid)

    out = tmp_path / "from_arrays.h5"
    assert ctrl.save_fit(sid, output_path=str(out))["ok"] is True
    assert out.exists()
    ctrl.close_session(sid)

    # The written file must be complete NXcanSAS — reduced data included, not
    # a results-only stub — so it reopens as a normal session.
    reopened = ctrl.open_dataset(str(out))
    assert "error" not in reopened
    assert reopened["summary"]["n_points"] > 0
    ctrl.close_session(reopened["session_id"])
