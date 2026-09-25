"""ZMQ service exposing ``pyirena.api`` to a remote orchestrator as JSON.

One UTF-8 JSON document per message, plain REQ/REP, synchronous. See
``docs/zmq_service.md`` for the protocol and ``planning/zmq-service/`` for why
it is shaped this way.

Nothing here imports ``zmq`` at module level: ``pyirena.zmq.protocol`` and
``pyirena.zmq.options`` are pure Python and hold the logic and most of the
tests, while ``pyirena.zmq.server`` imports pyzmq lazily inside ``main()``.
That keeps ``import pyirena`` working with no extras installed and lets the
protocol be tested without a socket.
"""
from __future__ import annotations

__all__ = ["ServerOptions", "PROTOCOL", "handle_request"]


def __getattr__(name: str):
    # Lazy so that importing the package costs nothing and pulls in no pyzmq.
    if name == "ServerOptions":
        from pyirena.zmq.options import ServerOptions
        return ServerOptions
    if name in ("PROTOCOL", "handle_request"):
        from pyirena.zmq import protocol
        return getattr(protocol, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
