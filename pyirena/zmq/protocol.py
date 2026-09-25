"""Request/reply envelope for the ZMQ service. No ``zmq`` import lives here.

All of the service's logic is in this module — parse, validate, route, call,
wrap, serialise — so that it can be tested without a socket, and so that
``server.py`` is only a loop around it.

The wire format is one UTF-8 JSON document per message, in both directions,
because the orchestrator's existing client is ``send_string`` /
``recv_string``. There is no multipart and no binary framing.

Request::

    {"protocol": "pyirena-zmq/1", "id": "c7f3", "op": "call",
     "tool": "run_fit", "args": {"session_id": "a1b2c3d4"},
     "options": {"include_arrays": false}}

Reply::

    {"protocol": "pyirena-zmq/1", "id": "c7f3", "ok": true,
     "result": {...}, "elapsed_s": 1.84, "server": {"pyirena": "1.2.0"}}

Failure::

    {"protocol": "pyirena-zmq/1", "id": "c7f3", "ok": false,
     "error": {"code": "BAD_ARGUMENTS", "message": "...", "suggestion": "..."}}

Two design points worth keeping:

* **A tool-level error is ``ok: false``.** ``pyirena.api`` never raises; it
  returns ``{"error", "code", "suggestion"}``. Those are promoted to the
  envelope's error slot so a caller has exactly one place to look, rather than
  having to know that a successful transport can still carry a failed fit.
  (Still to confirm with the project owner — see planning/zmq-service/ open
  question 6 — and the one line that decides it is ``_reply_for_result``.)
* **Every reply is strictly serialisable.** ``to_strict_json`` runs over the
  whole envelope and ``json.dumps`` uses ``allow_nan=False``, so a NaN in one
  uncertainty cannot produce a reply the caller's parser rejects wholesale.
"""
from __future__ import annotations

import json
import logging
import time
from typing import Any, Optional, Tuple

from pyirena.api._json import to_strict_json
from pyirena.zmq.options import PER_CALL_OPTIONS, TOOL_ARGUMENT_OPTIONS, ServerOptions

logger = logging.getLogger(__name__)

PROTOCOL = "pyirena-zmq/1"
PROTOCOL_FAMILY = "pyirena-zmq/"

#: Operations the service answers. Kept as data so `ping` can report them and
#: an unknown op can list them.
OPS = (
    "ping",
    "server_info",
    "list_categories",
    "list_tools",
    "describe_tool",
    "call",
    "open_dataset",
    "close_session",
    "list_sessions",
    "get_session_summary",
)

_STARTED_AT = time.time()

#: Sessions abandoned mid-fit by a deadline. The orphaned computation may
#: still be mutating the model, so the session's state is untrustworthy until
#: the caller closes it. Kept here rather than in pyirena.api.control because
#: it is a property of this transport's deadline, not of the session itself.
_STALE_SESSIONS: set[str] = set()


def mark_session_stale(session_id: str) -> None:
    if session_id:
        _STALE_SESSIONS.add(session_id)


def is_session_stale(session_id: str) -> bool:
    return session_id in _STALE_SESSIONS


def clear_stale(session_id: str) -> None:
    _STALE_SESSIONS.discard(session_id)


def reset_stale() -> None:
    """Forget every stale marking. For tests."""
    _STALE_SESSIONS.clear()


# ---------------------------------------------------------------------------
# Envelope construction
# ---------------------------------------------------------------------------

def _server_block() -> dict:
    from pyirena import __version__
    return {"pyirena": __version__, "protocol": PROTOCOL}


def make_error(code: str, message: str, suggestion: str = "",
               request_id: Any = None, **extra: Any) -> dict:
    error = {"code": code, "message": message}
    if suggestion:
        error["suggestion"] = suggestion
    error.update(extra)
    return {
        "protocol": PROTOCOL,
        "id": request_id,
        "ok": False,
        "error": error,
        "server": _server_block(),
    }


def make_reply(result: Any, request_id: Any = None, elapsed_s: Optional[float] = None) -> dict:
    reply = {
        "protocol": PROTOCOL,
        "id": request_id,
        "ok": True,
        "result": result,
        "server": _server_block(),
    }
    if elapsed_s is not None:
        reply["elapsed_s"] = round(elapsed_s, 4)
    return reply


def encode(reply: dict) -> str:
    """Serialise a reply as JSON a strict parser will accept.

    Never raises: if the reply itself cannot be serialised that is a bug here,
    and answering with a minimal INTERNAL_ERROR beats leaving the caller's REQ
    socket hanging until it times out.
    """
    try:
        return json.dumps(to_strict_json(reply), allow_nan=False)
    except Exception as exc:                                  # pragma: no cover
        logger.exception("could not serialise reply")
        return json.dumps({
            "protocol": PROTOCOL,
            "id": reply.get("id") if isinstance(reply, dict) else None,
            "ok": False,
            "error": {
                "code": "INTERNAL_ERROR",
                "message": f"Reply could not be serialised: {exc}",
            },
        })


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_request(raw: Any, options: Optional[ServerOptions] = None
                  ) -> Tuple[Optional[dict], Optional[dict]]:
    """Parse and validate one request. Returns (envelope, error_reply).

    Exactly one of the two is not None. Split out from :func:`handle_request`
    so the server can learn the request id — and which session it names —
    before starting work that may hit the deadline.
    """
    options = options or ServerOptions()

    if isinstance(raw, (bytes, bytearray)):
        if len(raw) > options.max_message_bytes:
            return None, make_error(
                "MESSAGE_TOO_LARGE",
                f"Request is {len(raw)} bytes; this server accepts "
                f"{options.max_message_bytes}.",
                "Send fewer points, or ask the operator to raise max_message_mb.",
            )
        try:
            raw = raw.decode("utf-8")
        except UnicodeDecodeError as exc:
            return None, make_error(
                "BAD_JSON", f"Request is not valid UTF-8: {exc}",
                "Send one UTF-8 JSON document per message.",
            )

    if not isinstance(raw, str):
        return None, make_error(
            "BAD_JSON", f"Expected a JSON string, got {type(raw).__name__}.",
            "Send one UTF-8 JSON document per message.",
        )

    if len(raw.encode("utf-8")) > options.max_message_bytes:
        return None, make_error(
            "MESSAGE_TOO_LARGE",
            f"Request exceeds the {options.max_message_mb} MB limit.",
            "Send fewer points, or ask the operator to raise max_message_mb.",
        )

    # The orchestrator's other services take bare string commands, so a bare
    # "ping" is answered rather than rejected on a pedantry. Anything else
    # gets told what the envelope looks like.
    stripped = raw.strip()
    if stripped.lower() in ("ping", '"ping"'):
        return {"op": "ping", "id": None}, None

    try:
        envelope = json.loads(stripped)
    except (json.JSONDecodeError, ValueError) as exc:
        return None, make_error(
            "BAD_JSON", f"Could not parse request as JSON: {exc}",
            'Send an object such as {"op": "ping"} or '
            '{"op": "call", "tool": "run_fit", "args": {"session_id": "..."}}.',
        )

    if not isinstance(envelope, dict):
        return None, make_error(
            "BAD_ENVELOPE",
            f"Request must be a JSON object, got {type(envelope).__name__}.",
            'Example: {"op": "list_categories"}.',
        )

    request_id = envelope.get("id")

    protocol = envelope.get("protocol")
    if protocol is not None and not str(protocol).startswith(PROTOCOL_FAMILY):
        return None, make_error(
            "UNSUPPORTED_PROTOCOL",
            f"This server speaks '{PROTOCOL}', the request said '{protocol}'.",
            "Omit 'protocol' or set it to " + PROTOCOL + ".",
            request_id=request_id,
        )

    op = envelope.get("op")
    if not op:
        return None, make_error(
            "BAD_ENVELOPE", "Request has no 'op'.",
            f"One of: {', '.join(OPS)}.",
            request_id=request_id,
        )
    if op not in OPS:
        return None, make_error(
            "UNKNOWN_OP", f"Unknown op '{op}'.",
            f"One of: {', '.join(OPS)}.",
            request_id=request_id,
        )

    args = envelope.get("args")
    if args is not None and not isinstance(args, dict):
        return None, make_error(
            "BAD_ENVELOPE", "'args' must be an object.",
            'Example: {"op": "call", "tool": "run_fit", "args": {"session_id": "..."}}.',
            request_id=request_id,
        )

    return envelope, None


def session_id_of(envelope: Optional[dict]) -> Optional[str]:
    """The session a request is about, if it names one. Used by the deadline."""
    if not isinstance(envelope, dict):
        return None
    args = envelope.get("args")
    if isinstance(args, dict):
        sid = args.get("session_id")
        if isinstance(sid, str):
            return sid
    sid = envelope.get("session_id")
    return sid if isinstance(sid, str) else None


# ---------------------------------------------------------------------------
# Routing
# ---------------------------------------------------------------------------

def _reply_for_result(result: Any, request_id: Any, elapsed_s: float) -> dict:
    """Wrap an api result, promoting a tool-level error to the envelope."""
    if isinstance(result, dict) and "error" in result:
        return make_error(
            str(result.get("code") or "TOOL_ERROR"),
            str(result["error"]),
            str(result.get("suggestion") or ""),
            request_id=request_id,
        )
    return make_reply(result, request_id=request_id, elapsed_s=elapsed_s)


def _ping_result(options: ServerOptions) -> dict:
    from pyirena.api.control import list_open_sessions
    from pyirena.api.control.session import all_sessions

    return {
        "pong": True,
        "uptime_s": round(time.time() - _STARTED_AT, 1),
        "open_sessions": len(all_sessions()),
        "stale_sessions": len(_STALE_SESSIONS),
        "max_sessions": options.max_sessions,
        "ops": list(OPS),
        **{k: v for k, v in list_open_sessions().items() if k == "count"},
    }


def _server_info_result(options: ServerOptions) -> dict:
    from pyirena.api import dispatch

    categories = dispatch.list_categories(options.profile)["categories"]
    return {
        "server": _server_block(),
        "options": options.to_dict(),
        "per_call_options": sorted(PER_CALL_OPTIONS),
        "capabilities": {
            "async_jobs": False,
            "binary_arrays": False,
            "images": options.allow_images,
            "files": options.allow_files,
            "request_budget_s": options.request_budget_s,
            "max_input_points": options.max_input_points,
            "max_message_mb": options.max_message_mb,
        },
        "categories": categories,
        "tool_count": sum(c["tool_count"] for c in categories),
    }


def _apply_option_defaults(tool: str, args: dict, options: ServerOptions) -> dict:
    """Fill a tool's include_arrays / max_points from the effective options.

    Only when the caller did not name them, and only for tools whose schema
    actually has them, so the server's defaults are real defaults rather than
    an override of what the caller asked for.
    """
    from pyirena.api import dispatch

    entry = dispatch._REGISTRY.get(tool)
    if entry is None:
        return args
    properties = (entry["schema"].get("input_schema") or {}).get("properties") or {}
    filled = dict(args)
    for name in TOOL_ARGUMENT_OPTIONS:
        if name in properties and name not in filled:
            filled[name] = getattr(options, name)
    return filled


def _dispatch_op(envelope: dict, options: ServerOptions) -> Any:
    """Run one operation and return the raw api result."""
    from pyirena.api import control as ctrl
    from pyirena.api import dispatch

    op = envelope["op"]
    args = envelope.get("args") or {}

    if op == "ping":
        return _ping_result(options)

    if op == "server_info":
        return _server_info_result(options)

    if op == "list_categories":
        return dispatch.list_categories(options.profile)

    if op == "list_tools":
        category = args.get("category") or envelope.get("category")
        if not category:
            return {"error": "list_tools needs a category.",
                    "code": "BAD_ARGUMENTS",
                    "suggestion": "Call list_categories first."}
        return dispatch.list_tools(category, options.profile)

    if op == "describe_tool":
        name = args.get("name") or envelope.get("tool") or envelope.get("name")
        if not name:
            return {"error": "describe_tool needs a tool name.",
                    "code": "BAD_ARGUMENTS",
                    "suggestion": "Call list_tools(category) first."}
        return dispatch.describe_tool(name, options.profile)

    if op == "open_dataset":
        # Deliberately routed to the array form: this transport has no shared
        # filesystem, so "open a dataset" can only mean "here is the data".
        args = dict(args)
        args.pop("file_path", None)
        return ctrl.open_dataset_from_data(**args)

    if op == "close_session":
        sid = session_id_of(envelope)
        clear_stale(sid or "")
        return ctrl.close_session(sid) if sid else {
            "error": "close_session needs a session_id.", "code": "BAD_ARGUMENTS",
            "suggestion": "Pass args.session_id."}

    if op == "list_sessions":
        result = ctrl.list_open_sessions()
        for row in result.get("sessions", []):
            row["stale"] = is_session_stale(row.get("session_id", ""))
        return result

    if op == "get_session_summary":
        sid = session_id_of(envelope)
        if not sid:
            return {"error": "get_session_summary needs a session_id.",
                    "code": "BAD_ARGUMENTS", "suggestion": "Pass args.session_id."}
        result = ctrl.get_session_summary(sid)
        if isinstance(result, dict) and "error" not in result:
            result["stale"] = is_session_stale(sid)
        return result

    # op == "call"
    tool = envelope.get("tool") or args.get("tool")
    if not tool:
        return {"error": "call needs a tool name.", "code": "BAD_ARGUMENTS",
                "suggestion": "Set 'tool'; see list_tools(category)."}

    sid = session_id_of(envelope)
    if sid and is_session_stale(sid):
        return {
            "error": (
                f"Session '{sid}' was abandoned mid-fit when a call hit the "
                "server's time budget, so its model state cannot be trusted."
            ),
            "code": "SESSION_STALE",
            "suggestion": (
                "Close it and open a new session, or re-run the fit on a "
                "narrower Q range."
            ),
        }

    call_args = _apply_option_defaults(tool, args, options)
    result = dispatch.call_tool(tool, call_args, options.profile)

    if isinstance(result, dict) and result.get("code") == "NOT_AVAILABLE_IN_PROFILE":
        result = dict(result)
        result["code"] = "NOT_AVAILABLE_OVER_ZMQ"
    return result


def handle_parsed(envelope: dict, options: Optional[ServerOptions] = None) -> dict:
    """Run one already-parsed request and return the reply envelope."""
    options = options or ServerOptions()
    request_id = envelope.get("id")

    effective, option_error = options.with_overrides(envelope.get("options"))
    if option_error:
        return make_error(
            option_error["code"], option_error["message"],
            option_error.get("suggestion", ""), request_id=request_id,
        )

    started = time.perf_counter()
    try:
        result = _dispatch_op(envelope, effective)
    except Exception as exc:
        # The api layer returns errors rather than raising, so reaching here
        # is a bug. The traceback goes to the log; the caller gets the message
        # only, since it may name server-side paths.
        logger.exception("unhandled exception in op %r", envelope.get("op"))
        return make_error(
            "INTERNAL_ERROR", f"{type(exc).__name__}: {exc}",
            "This is a server-side bug; the traceback is in the service log.",
            request_id=request_id,
        )

    return _reply_for_result(result, request_id, time.perf_counter() - started)


def handle_request(raw: Any, options: Optional[ServerOptions] = None) -> str:
    """Parse, run and serialise one request. Never raises.

    The whole service in one function; ``server.py`` adds only a socket and a
    deadline around it.
    """
    envelope, error_reply = parse_request(raw, options)
    if error_reply is not None:
        return encode(error_reply)
    return encode(handle_parsed(envelope, options))
