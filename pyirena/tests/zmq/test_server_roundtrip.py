"""End-to-end over a real socket: the workflow the orchestrator will run.

Skipped when pyzmq is not installed, so the suite still passes on an install
without the ``[zmq]`` extra.

Half of these use the shipped :class:`PyIrenaClient`; the rest drive a bare
``send_string`` / ``recv_string`` REQ socket, because that is the client the
orchestrator actually has. Testing only our own client would prove the two
halves of our own code agree and nothing about the contract we were given.
"""
from __future__ import annotations

import json
import threading
import time

import numpy as np
import pytest

zmq = pytest.importorskip("zmq")

import pyirena.api.control as ctrl  # noqa: E402
from pyirena.zmq import protocol  # noqa: E402
from pyirena.zmq.client import PyIrenaClient, PyIrenaError, PyIrenaTimeout  # noqa: E402
from pyirena.zmq.options import ServerOptions  # noqa: E402
from pyirena.zmq.server import _Stopping, run  # noqa: E402


def _free_port() -> int:
    import socket as pysocket

    with pysocket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture(scope="module")
def service():
    """A real REP server on a random port, in a thread, for this module."""
    options = ServerOptions(bind=f"tcp://127.0.0.1:{_free_port()}", request_budget_s=30)
    stop = _Stopping()
    thread = threading.Thread(target=run, args=(options, stop), daemon=True)
    thread.start()

    # Wait for the bind rather than sleeping a guessed interval.
    deadline = time.time() + 10
    while time.time() < deadline:
        try:
            with PyIrenaClient(options.bind, timeout_ms=300) as probe:
                probe.ping()
            break
        except (PyIrenaTimeout, PyIrenaError):
            time.sleep(0.05)
    else:                                                    # pragma: no cover
        stop.requested = True
        pytest.fail("the service never came up")

    yield options
    stop.requested = True
    thread.join(timeout=10)


@pytest.fixture
def client(service):
    with PyIrenaClient(service.bind, timeout_ms=30_000) as c:
        yield c


@pytest.fixture
def raw_socket(service):
    """A bare REQ socket — the orchestrator's own client shape."""
    socket = zmq.Context.instance().socket(zmq.REQ)
    socket.setsockopt(zmq.LINGER, 0)
    socket.connect(service.bind)
    yield socket
    socket.close()


def _curve(n=300):
    rng = np.random.default_rng(13)
    q = np.logspace(-3, 0, n)
    I = 1000.0 * np.exp(-(q**2) * 150.0**2 / 3) + 1e-6 * q**-4.0 + 0.01
    I = I * (1 + 0.02 * rng.standard_normal(q.size))
    return q, I, 0.02 * np.abs(I)


# ---------------------------------------------------------------------------
# The orchestrator's own client shape: send_string / recv_string
# ---------------------------------------------------------------------------

def test_a_bare_req_socket_can_drive_a_whole_fit(raw_socket):
    def send(payload):
        raw_socket.send_string(
            payload if isinstance(payload, str) else json.dumps(payload)
        )
        reply = json.loads(raw_socket.recv_string())
        assert reply["ok"], reply.get("error")
        return reply["result"]

    assert send("ping")["pong"] is True

    q, I, e = _curve()
    sid = send({"op": "open_dataset",
                "args": {"q": q.tolist(), "intensity": I.tolist(),
                         "error": e.tolist(), "label": "orchestrator"}})["session_id"]

    send({"op": "call", "tool": "select_model",
          "args": {"session_id": sid, "model_name": "unified_fit"}})
    send({"op": "call", "tool": "add_unified_level", "args": {"session_id": sid}})
    send({"op": "call", "tool": "run_fit", "args": {"session_id": sid}})

    results = send({"op": "call", "tool": "export_results", "args": {"session_id": sid}})
    assert results["tool"] == "unified_fit"
    assert results["data"]["label"] == "orchestrator"
    assert np.isfinite(results["quality"]["reduced_chi_squared"])

    assert send({"op": "close_session", "args": {"session_id": sid}})["ok"] is True


def test_the_server_survives_everything_malformed(raw_socket):
    """Each of these must come back as a reply, not a dropped connection."""
    for payload in ("", "   ", "not json", "[1,2,3]", '{"op": "nope"}',
                    '{"op": "call"}', '{"protocol": "other/9", "op": "ping"}',
                    '{"op": "open_dataset", "args": {"q": [1], "intensity": []}}'):
        raw_socket.send_string(payload)
        reply = json.loads(raw_socket.recv_string())
        assert reply["ok"] is False
        assert reply["error"]["code"]

    # ...and it is still healthy afterwards.
    raw_socket.send_string("ping")
    assert json.loads(raw_socket.recv_string())["ok"] is True


# ---------------------------------------------------------------------------
# The shipped reference client
# ---------------------------------------------------------------------------

def test_discovery_then_a_fit(client):
    info = client.server_info()
    assert info["options"]["profile"] == "json_only"
    assert info["capabilities"]["async_jobs"] is False

    categories = [c["name"] for c in client.list_categories()["categories"]]
    assert "data" not in categories                     # hidden: it writes files
    assert "unified" in categories

    schema = client.describe_tool("run_fit")
    assert schema["input_schema"]["properties"]["session_id"]

    q, I, e = _curve()
    sid = client.open_dataset(q, I, e, label="client")["session_id"]
    try:
        client.call("select_model", session_id=sid, model_name="unified_fit")
        client.call("add_unified_level", session_id=sid)
        client.call("run_fit", session_id=sid)

        results = client.call("export_results", session_id=sid)
        assert "arrays" not in results
        with_arrays = client.call("export_results", session_id=sid,
                                  options={"include_arrays": True})
        assert with_arrays["arrays"]["n_points"] == len(q)
    finally:
        client.close_session(sid)


def test_errors_arrive_as_typed_exceptions(client):
    with pytest.raises(PyIrenaError) as excinfo:
        client.call("run_fit", session_id="no-such-session")
    assert excinfo.value.code == "NO_SESSION"
    assert excinfo.value.suggestion

    with pytest.raises(PyIrenaError) as excinfo:
        client.call("save_fit", session_id="whatever")
    assert excinfo.value.code == "NOT_AVAILABLE_OVER_ZMQ"


def test_a_2000_point_curve_round_trips(client):
    """The payload size the orchestrator actually sends."""
    q, I, e = _curve(2000)
    opened = client.open_dataset(q, I, e)
    sid = opened["session_id"]
    try:
        assert opened["summary"]["n_points"] == 2000
        summary = client.session_summary(sid)
        assert summary["q_min"] == pytest.approx(q.min())
        assert summary["q_max"] == pytest.approx(q.max())
    finally:
        client.close_session(sid)


def test_the_client_recovers_from_a_timeout(service):
    """A timed-out REQ socket is dead; the client must rebuild it silently."""
    with PyIrenaClient(service.bind, timeout_ms=30_000) as c:
        # Nothing is listening here, so this must time out and reconnect.
        dead = PyIrenaClient(f"tcp://127.0.0.1:{_free_port()}", timeout_ms=200)
        with pytest.raises(PyIrenaTimeout):
            dead.ping()
        with pytest.raises(PyIrenaTimeout):
            dead.ping()          # would raise ZMQError(EFSM) without the rebuild
        dead.close()

        assert c.ping()["pong"] is True


def test_sessions_are_visible_and_closable_across_connections(service):
    with PyIrenaClient(service.bind, timeout_ms=30_000) as first:
        q, I, e = _curve(100)
        sid = first.open_dataset(q, I, e, label="shared")["session_id"]

    with PyIrenaClient(service.bind, timeout_ms=30_000) as second:
        listed = second.list_sessions()["sessions"]
        assert sid in [s["session_id"] for s in listed]
        second.close_session(sid)
        assert sid not in [s["session_id"] for s in second.list_sessions()["sessions"]]


def test_a_cleaned_curve_reports_what_was_removed(client):
    q, I, e = _curve(100)
    I = I.copy()
    I[5] = 0.0                       # beamstop zero: removed
    q = q.copy()
    q[0] = 0.0                       # direct beam: removed

    opened = client.open_dataset(q, I, e)
    try:
        cleaning = opened["summary"]["cleaning"]
        assert cleaning["n_removed_q"] == 1
        assert cleaning["n_removed_i"] == 1
        assert opened["summary"]["n_points"] == 98
    finally:
        client.close_session(opened["session_id"])


def test_module_state_is_left_clean(service):
    protocol.reset_stale()
    for s in list(ctrl.list_open_sessions()["sessions"]):
        ctrl.close_session(s["session_id"])
