"""A small reference client for the pyIrena ZMQ service.

Ships so the orchestrator team has something known-good to read and copy. It
depends on ``pyzmq`` and nothing else in pyIrena, so it can be lifted into
another codebase whole.

    from pyirena.zmq.client import PyIrenaClient

    with PyIrenaClient("tcp://workstation:9865") as pyirena:
        sid = pyirena.open_dataset(q, intensity, error)["session_id"]
        pyirena.call("select_model", session_id=sid, model_name="unified_fit")
        pyirena.call("add_unified_level", session_id=sid)
        pyirena.call("run_fit", session_id=sid)
        results = pyirena.call("export_results", session_id=sid)
        pyirena.close_session(sid)

**The one thing to copy if you copy nothing else** is the reconnect in
:meth:`PyIrenaClient.request`. A ZMQ REQ socket enforces strict
send/recv alternation: once a request times out, that socket can never be used
again — the next ``send`` raises ``EFSM``. The socket must be closed and
rebuilt. This is the "lazy pirate" pattern from the ZMQ guide, and a client
that only raises on timeout (as the orchestrator's current one does) will work
exactly once before every later call fails for a reason that has nothing to do
with pyIrena.

The service is designed never to make you need it: it answers ``TIMEOUT``
inside its own budget rather than going silent. The reconnect is for the cases
the server cannot cover — a network drop, a killed process, a machine reboot.
"""
from __future__ import annotations

import json
import uuid
from typing import Any, Optional, Sequence

DEFAULT_TIMEOUT_MS = 60_000
PROTOCOL = "pyirena-zmq/1"


class PyIrenaError(RuntimeError):
    """A reply with ``ok: false``. Carries the service's error code."""

    def __init__(self, code: str, message: str, suggestion: str = "") -> None:
        super().__init__(f"[{code}] {message}" + (f"  ({suggestion})" if suggestion else ""))
        self.code = code
        self.message = message
        self.suggestion = suggestion


class PyIrenaTimeout(TimeoutError):
    """No reply within the client's timeout. The socket has been rebuilt."""


class PyIrenaClient:
    """Synchronous REQ client for one pyIrena ZMQ service."""

    def __init__(self, address: str, timeout_ms: int = DEFAULT_TIMEOUT_MS) -> None:
        self.address = address
        self.timeout_ms = timeout_ms
        self._socket = None
        self._poller = None
        self._connect()

    # --- plumbing ----------------------------------------------------------

    def _connect(self) -> None:
        import zmq

        context = zmq.Context.instance()
        self._socket = context.socket(zmq.REQ)
        # Drop unsent messages on close instead of blocking the process exit.
        self._socket.setsockopt(zmq.LINGER, 0)
        self._socket.connect(self.address)
        self._poller = zmq.Poller()
        self._poller.register(self._socket, zmq.POLLIN)

    def _reconnect(self) -> None:
        if self._socket is not None:
            self._poller.unregister(self._socket)
            self._socket.close()
        self._connect()

    def close(self) -> None:
        if self._socket is not None:
            self._poller.unregister(self._socket)
            self._socket.close()
            self._socket = None

    def __enter__(self) -> "PyIrenaClient":
        return self

    def __exit__(self, *_exc) -> None:
        self.close()

    # --- the protocol ------------------------------------------------------

    def request(self, op: str, *, timeout_ms: Optional[int] = None, **fields: Any) -> dict:
        """Send one request and return its reply envelope.

        Raises :class:`PyIrenaTimeout` if nothing comes back in time, having
        first rebuilt the socket so the client stays usable.
        """
        envelope = {"protocol": PROTOCOL, "id": uuid.uuid4().hex[:8], "op": op}
        envelope.update({k: v for k, v in fields.items() if v is not None})

        self._socket.send_string(json.dumps(envelope))
        if not dict(self._poller.poll(timeout_ms or self.timeout_ms)):
            # A timed-out REQ socket is unusable; rebuild before raising so
            # the caller can retry without knowing any of this.
            self._reconnect()
            raise PyIrenaTimeout(
                f"no reply within {timeout_ms or self.timeout_ms} ms "
                f"from {self.address} (op={op})"
            )
        return json.loads(self._socket.recv_string())

    def _result(self, reply: dict) -> Any:
        if not reply.get("ok"):
            error = reply.get("error") or {}
            raise PyIrenaError(
                error.get("code", "ERROR"),
                error.get("message", "unknown error"),
                error.get("suggestion", ""),
            )
        return reply.get("result")

    # --- operations --------------------------------------------------------

    def ping(self) -> dict:
        return self._result(self.request("ping", timeout_ms=5_000))

    def server_info(self) -> dict:
        """What this deployment allows: options, limits, categories."""
        return self._result(self.request("server_info"))

    def list_categories(self) -> dict:
        return self._result(self.request("list_categories"))

    def list_tools(self, category: str) -> dict:
        return self._result(self.request("list_tools", args={"category": category}))

    def describe_tool(self, name: str) -> dict:
        """The tool's JSON schema — hand it straight to an LLM as a tool definition."""
        return self._result(self.request("describe_tool", args={"name": name}))

    def open_dataset(
        self,
        q: Sequence[float],
        intensity: Sequence[float],
        error: Optional[Sequence[float]] = None,
        **kwargs: Any,
    ) -> dict:
        """Create a session from arrays. Returns ``{session_id, summary}``."""
        args = {"q": list(q), "intensity": list(intensity)}
        if error is not None:
            args["error"] = list(error)
        args.update(kwargs)
        return self._result(self.request("open_dataset", args=args))

    def call(self, tool: str, *, options: Optional[dict] = None, **args: Any) -> Any:
        """Call any dispatched tool by name, with its arguments as keywords."""
        return self._result(self.request("call", tool=tool, args=args, options=options))

    def list_sessions(self) -> dict:
        return self._result(self.request("list_sessions"))

    def session_summary(self, session_id: str) -> dict:
        return self._result(
            self.request("get_session_summary", args={"session_id": session_id})
        )

    def close_session(self, session_id: str) -> dict:
        return self._result(self.request("close_session", args={"session_id": session_id}))
