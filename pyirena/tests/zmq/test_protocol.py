"""Tests for the ZMQ envelope — no pyzmq needed, no socket involved.

``pyirena.zmq.protocol`` is deliberately the whole service minus the socket,
so this is where the protocol is actually pinned down: what a malformed
request gets back, how a tool-level error is reported, and — the one that
matters most for a machine on the other end — that every single reply is
strictly serialisable JSON.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

import pyirena.api.control as ctrl
from pyirena.zmq import protocol
from pyirena.zmq.options import ServerOptions

OPTIONS = ServerOptions()


def send(payload, options: ServerOptions = OPTIONS) -> dict:
    """Round-trip one request through the protocol, asserting strict JSON."""
    raw = payload if isinstance(payload, (str, bytes)) else json.dumps(payload)
    text = protocol.handle_request(raw, options)
    # The point of the service: a strict parser on the far end must accept it.
    assert "NaN" not in text and "Infinity" not in text
    reply = json.loads(text)
    assert reply["protocol"] == protocol.PROTOCOL
    assert isinstance(reply["ok"], bool)
    assert ("result" in reply) != ("error" in reply)
    return reply


@pytest.fixture(autouse=True)
def _clean_stale():
    protocol.reset_stale()
    yield
    protocol.reset_stale()


def _curve(n=200):
    rng = np.random.default_rng(11)
    q = np.logspace(-3, 0, n)
    I = 1000.0 * np.exp(-(q**2) * 150.0**2 / 3) + 1e-6 * q**-4.0 + 0.01
    I = I * (1 + 0.02 * rng.standard_normal(q.size))
    return q.tolist(), I.tolist(), (0.02 * np.abs(I)).tolist()


@pytest.fixture
def session():
    q, I, e = _curve()
    reply = send({"op": "open_dataset", "args": {"q": q, "intensity": I, "error": e}})
    sid = reply["result"]["session_id"]
    yield sid
    ctrl.close_session(sid)


# ---------------------------------------------------------------------------
# Malformed input — the server must survive anything and always answer
# ---------------------------------------------------------------------------

def test_bare_ping_string_is_answered():
    """The orchestrator's other services take bare string commands."""
    reply = send("ping")
    assert reply["ok"] and reply["result"]["pong"] is True


def test_garbage_is_a_readable_error_not_a_crash():
    reply = send("not json at all")
    assert reply["error"]["code"] == "BAD_JSON"
    assert "op" in reply["error"]["suggestion"]


def test_json_that_is_not_an_object():
    assert send("[1, 2, 3]")["error"]["code"] == "BAD_ENVELOPE"


def test_missing_op():
    assert send({"id": "x"})["error"]["code"] == "BAD_ENVELOPE"


def test_unknown_op_lists_the_real_ones():
    reply = send({"op": "teleport"})
    assert reply["error"]["code"] == "UNKNOWN_OP"
    for op in ("ping", "call", "open_dataset"):
        assert op in reply["error"]["suggestion"]


def test_wrong_protocol_version():
    reply = send({"protocol": "some-other-service/3", "op": "ping"})
    assert reply["error"]["code"] == "UNSUPPORTED_PROTOCOL"


def test_matching_protocol_family_is_accepted():
    assert send({"protocol": "pyirena-zmq/1", "op": "ping"})["ok"]


def test_args_must_be_an_object():
    assert send({"op": "call", "tool": "run_fit", "args": [1]})["error"]["code"] == "BAD_ENVELOPE"


def test_oversize_message_is_refused_before_parsing():
    tiny = ServerOptions(max_message_mb=1)
    reply = send({"op": "ping", "pad": "x" * (2 * 1024 * 1024)}, tiny)
    assert reply["error"]["code"] == "MESSAGE_TOO_LARGE"


def test_invalid_utf8_bytes():
    reply = json.loads(protocol.handle_request(b"\xff\xfe not utf8", OPTIONS))
    assert reply["error"]["code"] == "BAD_JSON"


def test_request_id_is_echoed_on_success_and_on_failure():
    assert send({"op": "ping", "id": "abc123"})["id"] == "abc123"
    assert send({"op": "nope", "id": "abc123"})["id"] == "abc123"


def test_an_unexpected_exception_becomes_internal_error(monkeypatch):
    def boom(*_a, **_kw):
        raise RuntimeError("something in the fit blew up")

    monkeypatch.setattr(protocol, "_dispatch_op", boom)
    reply = send({"op": "ping"})
    assert reply["error"]["code"] == "INTERNAL_ERROR"
    assert "something in the fit blew up" in reply["error"]["message"]


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

def test_ping_reports_liveness_and_capacity():
    result = send({"op": "ping"})["result"]
    assert result["pong"] is True
    assert result["uptime_s"] >= 0
    assert set(protocol.OPS) == set(result["ops"])


def test_server_info_describes_what_this_deployment_allows():
    result = send({"op": "server_info"})["result"]
    assert result["options"]["profile"] == "json_only"
    assert result["capabilities"]["async_jobs"] is False
    assert result["capabilities"]["images"] is False
    assert result["tool_count"] > 0
    assert "include_arrays" in result["per_call_options"]


def test_list_and_describe_round_trip():
    categories = send({"op": "list_categories"})["result"]["categories"]
    assert {"unified", "results"} <= {c["name"] for c in categories}

    tools = send({"op": "list_tools", "args": {"category": "results"}})["result"]["tools"]
    assert [t["name"] for t in tools] == ["analyze", "export_results"]

    schema = send({"op": "describe_tool", "args": {"name": "export_results"}})["result"]
    assert schema["name"] == "export_results"
    assert "session_id" in schema["input_schema"]["properties"]


def test_list_tools_without_a_category_says_so():
    assert send({"op": "list_tools"})["error"]["code"] == "BAD_ARGUMENTS"


# ---------------------------------------------------------------------------
# The json_only profile
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tool", ["save_fit", "get_fit_image", "average_data"])
def test_file_and_image_tools_are_refused(tool, session):
    reply = send({"op": "call", "tool": tool, "args": {"session_id": session}})
    assert reply["error"]["code"] == "NOT_AVAILABLE_OVER_ZMQ"


def test_images_can_be_enabled_by_the_operator(session):
    permissive = ServerOptions(allow_images=True)
    # Still refused for files...
    assert send({"op": "call", "tool": "save_fit", "args": {"session_id": session}},
                permissive)["error"]["code"] == "NOT_AVAILABLE_OVER_ZMQ"
    # ...but an image tool is now reachable (it errors on NO_MODEL, not on profile).
    reply = send({"op": "call", "tool": "get_fit_image", "args": {"session_id": session}},
                 permissive)
    assert reply["error"]["code"] != "NOT_AVAILABLE_OVER_ZMQ"


def test_open_dataset_ignores_a_file_path(session):
    """Over this transport 'open a dataset' can only mean 'here is the data'."""
    q, I, _ = _curve()
    reply = send({"op": "open_dataset",
                  "args": {"q": q, "intensity": I, "file_path": "/etc/passwd"}})
    assert reply["ok"]
    assert reply["result"]["summary"]["file"] is None
    ctrl.close_session(reply["result"]["session_id"])


# ---------------------------------------------------------------------------
# Tool errors and options
# ---------------------------------------------------------------------------

def test_a_tool_level_error_is_promoted_to_the_envelope():
    reply = send({"op": "call", "tool": "run_fit", "args": {"session_id": "nope"}})
    assert reply["ok"] is False
    assert reply["error"]["code"] == "NO_SESSION"
    assert "suggestion" in reply["error"]


def test_bad_tool_name_suggests_a_real_one():
    reply = send({"op": "call", "tool": "run_ft", "args": {}})
    assert reply["error"]["code"] == "UNKNOWN_TOOL"


def test_call_without_a_tool_name():
    assert send({"op": "call"})["error"]["code"] == "BAD_ARGUMENTS"


def test_an_unknown_per_call_option_is_refused_not_ignored():
    reply = send({"op": "ping", "options": {"include_arrys": True}})
    assert reply["error"]["code"] == "BAD_OPTION"
    assert "include_arrays" in reply["error"]["suggestion"]


def test_a_call_cannot_widen_the_servers_policy():
    reply = send({"op": "ping", "options": {"allow_images": True}})
    assert reply["error"]["code"] == "BAD_OPTION"


def test_server_defaults_fill_in_tool_arguments(session):
    """include_arrays is a server option AND an export_results argument."""
    send({"op": "call", "tool": "select_model",
          "args": {"session_id": session, "model_name": "unified_fit"}})
    send({"op": "call", "tool": "add_unified_level", "args": {"session_id": session}})
    send({"op": "call", "tool": "run_fit", "args": {"session_id": session}})

    default = send({"op": "call", "tool": "export_results",
                    "args": {"session_id": session}})["result"]
    assert "arrays" not in default

    per_call = send({"op": "call", "tool": "export_results",
                     "args": {"session_id": session},
                     "options": {"include_arrays": True}})["result"]
    assert "arrays" in per_call

    server_on = send({"op": "call", "tool": "export_results",
                      "args": {"session_id": session}},
                     ServerOptions(include_arrays=True))["result"]
    assert "arrays" in server_on

    # An explicit argument always beats the server's default.
    explicit_off = send({"op": "call", "tool": "export_results",
                         "args": {"session_id": session, "include_arrays": False}},
                        ServerOptions(include_arrays=True))["result"]
    assert "arrays" not in explicit_off


# ---------------------------------------------------------------------------
# Sessions
# ---------------------------------------------------------------------------

def test_session_lifecycle_over_the_wire(session):
    listed = send({"op": "list_sessions"})["result"]
    assert session in [s["session_id"] for s in listed["sessions"]]
    assert all("stale" in s for s in listed["sessions"])

    summary = send({"op": "get_session_summary",
                    "args": {"session_id": session}})["result"]
    assert summary["file"] is None and summary["stale"] is False


def test_session_summary_works_after_a_modeling_fit(session):
    """Modeling stores a result dataclass; the summary must not assume a dict."""
    send({"op": "call", "tool": "select_modeling_model", "args": {"session_id": session}})
    send({"op": "call", "tool": "add_population",
          "args": {"session_id": session, "population_type": "size_dist"}})
    send({"op": "call", "tool": "run_modeling_fit", "args": {"session_id": session}})

    reply = send({"op": "get_session_summary", "args": {"session_id": session}})
    assert reply["ok"], reply
    assert reply["result"]["model"] == "modeling"


def test_a_stale_session_is_refused_until_it_is_closed(session):
    protocol.mark_session_stale(session)

    reply = send({"op": "call", "tool": "run_fit", "args": {"session_id": session}})
    assert reply["error"]["code"] == "SESSION_STALE"
    assert send({"op": "get_session_summary",
                 "args": {"session_id": session}})["result"]["stale"] is True

    send({"op": "close_session", "args": {"session_id": session}})
    assert protocol.is_session_stale(session) is False


def test_close_session_without_an_id():
    assert send({"op": "close_session"})["error"]["code"] == "BAD_ARGUMENTS"


# ---------------------------------------------------------------------------
# A whole fit, through the envelope only
# ---------------------------------------------------------------------------

def test_a_full_workflow_never_produces_an_unserialisable_reply(session):
    steps = [
        {"op": "call", "tool": "select_model",
         "args": {"session_id": session, "model_name": "unified_fit"}},
        {"op": "call", "tool": "add_unified_level", "args": {"session_id": session}},
        {"op": "call", "tool": "get_model_parameters", "args": {"session_id": session}},
        {"op": "call", "tool": "run_fit", "args": {"session_id": session}},
        {"op": "call", "tool": "get_fit_quality", "args": {"session_id": session}},
        {"op": "call", "tool": "get_residuals", "args": {"session_id": session}},
        {"op": "call", "tool": "export_results",
         "args": {"session_id": session}, "options": {"include_arrays": True}},
    ]
    for step in steps:                       # send() asserts strict JSON on each
        reply = send(step)
        assert reply["ok"], (step["tool"], reply.get("error"))
    assert reply["result"]["tool"] == "unified_fit"


# ---------------------------------------------------------------------------
# analyze — the one-call path the orchestrator will mostly use
# ---------------------------------------------------------------------------

def test_analyze_fits_a_curve_from_a_config_in_one_request():
    """One round trip: data in, fitted results out, no session to manage."""
    q, I, e = _curve(300)
    config = {
        "_pyirena_config": {"tool": "unified_fit"},
        "unified_fit": {
            "num_levels": 1,
            "background": {"value": 0.01, "fit": True},
            "levels": [{
                "G":  {"value": 900.0, "fit": True,  "low_limit": 1.0, "high_limit": 1e6},
                "Rg": {"value": 140.0, "fit": True,  "low_limit": 10.0, "high_limit": 1e4},
                "B":  {"value": 1e-6,  "fit": True,  "low_limit": 0.0, "high_limit": 1.0},
                "P":  {"value": 4.0,   "fit": False, "low_limit": 0.0, "high_limit": 6.0},
            }],
        },
    }
    reply = send({"op": "call", "tool": "analyze",
                  "args": {"data": {"q": q, "intensity": I, "error": e},
                           "config": config}})
    assert reply["ok"], reply.get("error")
    result = reply["result"]
    assert result["tool"] == "unified_fit"
    assert result["analyze"]["config_applied"] is True
    assert result["quality"]["reduced_chi_squared"] is not None

    # The point of analyze: it owns its session and leaves nothing behind.
    assert send({"op": "list_sessions"})["result"]["count"] == 0


def test_analyze_is_reachable_in_the_json_only_profile():
    """It takes no path and returns no image, so a remote caller may use it."""
    tools = send({"op": "list_tools", "args": {"category": "results"}})["result"]["tools"]
    assert "analyze" in [t["name"] for t in tools]


def test_a_bad_config_is_refused_without_opening_a_session():
    q, I, _ = _curve(50)
    reply = send({"op": "call", "tool": "analyze",
                  "args": {"data": {"q": q, "intensity": I},
                           "config": {"nonsense": {}}}})
    assert reply["error"]["code"] == "BAD_CONFIG"
    assert send({"op": "list_sessions"})["result"]["count"] == 0


def test_a_profile_refusal_has_one_name_whatever_op_asked():
    """`call` and `describe_tool` must not disagree about the same refusal."""
    for request in ({"op": "call", "tool": "save_fit", "args": {"session_id": "x"}},
                    {"op": "describe_tool", "args": {"name": "save_fit"}}):
        assert send(request)["error"]["code"] == "NOT_AVAILABLE_OVER_ZMQ", request["op"]
