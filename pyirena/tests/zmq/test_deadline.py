"""Tests for the request deadline — the thing that makes a synchronous
contract safe.

The orchestrator's client waits 60 s and then raises, and a ZMQ REQ socket
that has timed out is unusable afterwards. So the server must answer inside
its own, shorter budget no matter what the fit is doing. A slow fit must
produce a ``TIMEOUT`` reply, not silence, and the session it abandoned must be
marked so nothing downstream reads a half-fitted model as if it were a result.

These drive ``serve_one`` directly rather than through a socket: the deadline
is server logic, and this keeps the test fast and deterministic.
"""
from __future__ import annotations

import json
import time
from concurrent.futures import ThreadPoolExecutor

import pytest

from pyirena.zmq import protocol, server
from pyirena.zmq.options import ServerOptions


@pytest.fixture
def pool():
    executor = ThreadPoolExecutor(max_workers=1)
    yield executor
    # Do not wait: a worker abandoned by the deadline is still sleeping, and
    # the real server does not wait for it either.
    executor.shutdown(wait=False, cancel_futures=True)


@pytest.fixture(autouse=True)
def _clean_stale():
    protocol.reset_stale()
    yield
    protocol.reset_stale()


def _slow_op(seconds: float):
    def handler(envelope, options=None):
        time.sleep(seconds)
        return protocol.make_reply({"finished": True}, envelope.get("id"))
    return handler


def _serve(payload, pool, options):
    return json.loads(server.serve_one(json.dumps(payload), options, pool))


def test_a_fit_that_overruns_gets_a_timeout_reply(monkeypatch, pool):
    monkeypatch.setattr(protocol, "handle_parsed", _slow_op(1.0))
    options = ServerOptions(request_budget_s=0.3)

    started = time.perf_counter()
    reply = _serve({"op": "call", "tool": "run_fit",
                    "args": {"session_id": "abc"}, "id": "r1"}, pool, options)
    elapsed = time.perf_counter() - started

    assert reply["ok"] is False
    assert reply["error"]["code"] == "TIMEOUT"
    assert reply["id"] == "r1"
    # The whole point: the caller hears back promptly, not after the fit ends.
    assert elapsed < 2.0
    assert "run_fit" in reply["error"]["message"]
    assert "stale" in reply["error"]["suggestion"]


def test_the_abandoned_session_is_marked_stale(monkeypatch, pool):
    monkeypatch.setattr(protocol, "handle_parsed", _slow_op(1.0))
    options = ServerOptions(request_budget_s=0.3)

    _serve({"op": "call", "tool": "run_fit", "args": {"session_id": "abc"}},
           pool, options)
    assert protocol.is_session_stale("abc")


def test_a_fast_call_is_unaffected(pool):
    options = ServerOptions(request_budget_s=30)
    reply = _serve({"op": "ping", "id": "r2"}, pool, options)
    assert reply["ok"] and reply["result"]["pong"] is True
    assert protocol.is_session_stale("") is False


def test_a_call_may_shorten_its_own_deadline(monkeypatch, pool):
    monkeypatch.setattr(protocol, "handle_parsed", _slow_op(1.0))
    options = ServerOptions(request_budget_s=30)

    started = time.perf_counter()
    reply = _serve({"op": "call", "tool": "run_fit", "args": {"session_id": "s"},
                    "options": {"request_budget_s": 0.3}}, pool, options)
    assert reply["error"]["code"] == "TIMEOUT"
    assert time.perf_counter() - started < 2.0


def test_a_call_cannot_extend_the_deadline(pool):
    options = ServerOptions(request_budget_s=1)
    reply = _serve({"op": "ping", "options": {"request_budget_s": 600}}, pool, options)
    assert reply["error"]["code"] == "BAD_OPTION"


def test_a_malformed_request_never_reaches_the_worker(pool):
    options = ServerOptions(request_budget_s=30)
    reply = json.loads(server.serve_one("{{{", options, pool))
    assert reply["error"]["code"] == "BAD_JSON"


# ---------------------------------------------------------------------------
# Session hygiene
# ---------------------------------------------------------------------------

def _open_session():
    import numpy as np

    import pyirena.api.control as ctrl
    q = np.logspace(-3, 0, 50)
    I = 1.0 + q**-1
    return ctrl.open_dataset_from_data(q=q.tolist(), intensity=I.tolist())["session_id"]


def test_idle_sessions_are_evicted_past_their_ttl():
    import pyirena.api.control as ctrl

    sid = _open_session()
    server._touch_sessions()
    # Pretend it was last used an hour ago.
    server._LAST_TOUCHED[sid] = time.time() - 3600

    evicted = server.evict_sessions(ServerOptions(session_ttl_min=1))
    assert sid in evicted
    assert ctrl.get_session_summary(sid)["code"] == "NO_SESSION"


def test_the_oldest_idle_session_is_evicted_when_over_the_cap():
    import pyirena.api.control as ctrl
    from pyirena.api.control.session import all_sessions

    for s in all_sessions():                   # start from a clean registry
        ctrl.close_session(s.session_id)

    first, second, third = _open_session(), _open_session(), _open_session()
    now = time.time()
    server._LAST_TOUCHED.update({first: now - 300, second: now - 200, third: now})

    evicted = server.evict_sessions(ServerOptions(max_sessions=2, session_ttl_min=0))
    assert evicted == [first]
    assert ctrl.get_session_summary(third)["session_id"] == third

    for sid in (second, third):
        ctrl.close_session(sid)


def test_eviction_clears_a_stale_marking():
    sid = _open_session()
    protocol.mark_session_stale(sid)
    server._LAST_TOUCHED[sid] = time.time() - 3600

    server.evict_sessions(ServerOptions(session_ttl_min=1))
    assert protocol.is_session_stale(sid) is False
