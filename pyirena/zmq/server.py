"""``pyirena-zmq`` — the REP loop, a deadline, and nothing else.

Everything interesting is in :mod:`pyirena.zmq.protocol`; this module binds a
socket, reads a message, hands it over, and writes the answer back.

Two things it does own:

**The deadline.** The contract is synchronous and the orchestrator's client
waits 60 s and then raises — leaving a REQ socket that must be thrown away.
So the server must never be the reason the client times out. Each request runs
on a worker thread with ``request_budget_s`` to finish; if it overruns, the
caller gets a ``TIMEOUT`` reply while the orphaned computation is left to
finish and be discarded. Its session is marked stale, because its model is
half-fitted and anything read from it afterwards would be a lie.

**Session hygiene.** Sessions live in a process-wide dict with no owner but
this loop, so the loop evicts them: idle longer than ``session_ttl_min``, or
oldest-first once ``max_sessions`` is reached. Checked between requests, never
during one.

Only one request runs at a time. The session registry is a plain dict and
scipy fits are CPU-bound, so serialising is both the simple choice and the
safe one — with the consequence that ``ping`` is not answered while a fit
runs. The deadline bounds how long that can last.
"""
from __future__ import annotations

import logging
import signal
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import TimeoutError as FutureTimeout
from typing import Optional

from pyirena.zmq import protocol
from pyirena.zmq.options import ServerOptions, options_from_args

logger = logging.getLogger(__name__)

#: How often the poll loop wakes to notice a signal. Short enough that Ctrl+C
#: feels immediate; long enough to be free.
_POLL_MS = 250


class _Stopping:
    """Set by SIGINT/SIGTERM; the loop finishes the current request and exits."""

    def __init__(self) -> None:
        self.requested = False

    def __call__(self, signum, _frame) -> None:
        logger.info("signal %s received; shutting down after the current request", signum)
        self.requested = True


# ---------------------------------------------------------------------------
# Session hygiene
# ---------------------------------------------------------------------------

_LAST_TOUCHED: dict[str, float] = {}


def _touch_sessions() -> None:
    """Record 'seen now' for every live session, and forget dead ones."""
    from pyirena.api.control.session import all_sessions

    now = time.time()
    live = {s.session_id for s in all_sessions()}
    for sid in live:
        _LAST_TOUCHED.setdefault(sid, now)
    for sid in list(_LAST_TOUCHED):
        if sid not in live:
            del _LAST_TOUCHED[sid]


def _touch(session_id: Optional[str]) -> None:
    if session_id:
        _LAST_TOUCHED[session_id] = time.time()


def evict_sessions(options: ServerOptions) -> list[str]:
    """Close idle and surplus sessions. Returns the ids closed."""
    from pyirena.api.control import close_session
    from pyirena.api.control.session import all_sessions

    _touch_sessions()
    now = time.time()
    evicted: list[str] = []

    ttl_s = options.session_ttl_min * 60.0
    if ttl_s > 0:
        for sid, seen in sorted(_LAST_TOUCHED.items(), key=lambda kv: kv[1]):
            if now - seen > ttl_s:
                close_session(sid)
                protocol.clear_stale(sid)
                evicted.append(sid)

    # Oldest-idle first, so a caller working through several curves keeps the
    # one it is using.
    surplus = len(all_sessions()) - options.max_sessions
    if surplus > 0:
        by_age = sorted(
            (s.session_id for s in all_sessions()),
            key=lambda sid: _LAST_TOUCHED.get(sid, 0.0),
        )
        for sid in by_age[:surplus]:
            close_session(sid)
            protocol.clear_stale(sid)
            evicted.append(sid)

    for sid in evicted:
        _LAST_TOUCHED.pop(sid, None)
    if evicted:
        logger.info("evicted %d session(s): %s", len(evicted), ", ".join(evicted))
    return evicted


# ---------------------------------------------------------------------------
# One request
# ---------------------------------------------------------------------------

def serve_one(raw, options: ServerOptions, pool: ThreadPoolExecutor) -> str:
    """Handle one raw request under the deadline. Always returns a reply."""
    started = time.perf_counter()
    envelope, error_reply = protocol.parse_request(raw, options)
    if error_reply is not None:
        logger.info("rejected: %s", error_reply["error"]["code"])
        return protocol.encode(error_reply)

    request_id = envelope.get("id")
    op = envelope.get("op")
    tool = envelope.get("tool")
    session_id = protocol.session_id_of(envelope)

    # A per-call budget may shorten the deadline but never extend it; the
    # options object has already refused anything longer.
    effective, _ = options.with_overrides(envelope.get("options"))
    budget = effective.request_budget_s

    future = pool.submit(protocol.handle_parsed, envelope, options)
    try:
        reply = future.result(timeout=budget)
    except FutureTimeout:
        # The worker is still running and will be abandoned. Do NOT cancel and
        # do NOT touch the session's data from here: the point is to answer
        # the caller before its socket dies, and to make clear that whatever
        # that session holds afterwards is not a result.
        protocol.mark_session_stale(session_id or "")
        logger.warning("op=%s tool=%s exceeded the %.0fs budget; session %s marked stale",
                       op, tool, budget, session_id)
        reply = protocol.make_error(
            "TIMEOUT",
            f"{tool or op} exceeded this server's {budget:g} s budget.",
            (
                "Narrow the Q range, fix parameters, or fit a simpler model. "
                "The session is marked stale; close it and start again."
            ),
            request_id=request_id,
        )
    except Exception as exc:                                    # pragma: no cover
        logger.exception("worker failed")
        reply = protocol.make_error(
            "INTERNAL_ERROR", f"{type(exc).__name__}: {exc}",
            "This is a server-side bug; the traceback is in the service log.",
            request_id=request_id,
        )

    _touch(session_id)
    # A reply that creates a session should start that session's clock.
    result = reply.get("result")
    if isinstance(result, dict):
        _touch(result.get("session_id"))

    elapsed = time.perf_counter() - started
    logger.info(
        "id=%s op=%s tool=%s session=%s %s %.3fs",
        request_id, op, tool, session_id,
        "ok" if reply.get("ok") else reply.get("error", {}).get("code", "error"),
        elapsed,
    )
    return protocol.encode(reply)


# ---------------------------------------------------------------------------
# The loop
# ---------------------------------------------------------------------------

def run(options: ServerOptions, stop: Optional[_Stopping] = None) -> int:
    """Bind and serve until interrupted. Returns a process exit code."""
    try:
        import zmq
    except ModuleNotFoundError:
        print(
            "pyirena-zmq needs pyzmq.\n"
            "    pip install 'pyirena[zmq]'",
            file=sys.stderr,
        )
        return 2

    stop = stop or _Stopping()
    context = zmq.Context.instance()
    socket = context.socket(zmq.REP)
    socket.setsockopt(zmq.MAXMSGSIZE, options.max_message_bytes)
    # Do not let a queued reply hold the process open at shutdown.
    socket.setsockopt(zmq.LINGER, 0)

    try:
        socket.bind(options.bind)
    except Exception as exc:
        logger.error("could not bind %s: %s", options.bind, exc)
        print(f"pyirena-zmq: could not bind {options.bind}: {exc}", file=sys.stderr)
        socket.close()
        return 1

    logger.warning(
        "pyirena-zmq listening on %s (profile=%s, budget=%.0fs, max_message=%dMB)",
        options.bind, options.profile, options.request_budget_s, options.max_message_mb,
    )
    if options.bind.startswith("tcp://0.0.0.0") or options.bind.startswith("tcp://*"):
        logger.warning(
            "bound to every interface with no authentication — restrict this "
            "port to the orchestrator host at the firewall"
        )

    poller = zmq.Poller()
    poller.register(socket, zmq.POLLIN)

    # One worker: requests are serialised by design, and the thread exists so
    # the deadline can be enforced, not for concurrency.
    pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix="pyirena-zmq")
    served = 0
    try:
        while not stop.requested:
            if not dict(poller.poll(_POLL_MS)):
                continue
            try:
                raw = socket.recv(copy=True)
            except zmq.ZMQError as exc:                         # pragma: no cover
                logger.error("recv failed: %s", exc)
                continue

            reply = serve_one(raw, options, pool)
            try:
                socket.send_string(reply)
            except zmq.ZMQError as exc:                         # pragma: no cover
                logger.error("send failed (client gone?): %s", exc)
            served += 1
            evict_sessions(options)
    finally:
        poller.unregister(socket)
        socket.close()
        # A worker abandoned by the deadline may still be running; do not wait
        # for it to finish a fit nobody is waiting for.
        pool.shutdown(wait=False, cancel_futures=True)
        logger.warning("pyirena-zmq stopped after %d request(s)", served)
    return 0


def configure_logging(options: ServerOptions) -> None:
    """Log the way the rest of pyIrena does, plus an operator-chosen file.

    ``~/.pyirena/logs/zmq.log`` always (rotating, like gui/mcp/batch), stderr
    only when asked, and never stdout. ``--log-file`` adds a second rotating
    handler for deployments that collect logs elsewhere.
    """
    import logging.handlers

    from pyirena.logging_setup import (
        BACKUP_COUNT,
        MAX_BYTES,
        install_excepthook,
        setup_logging,
    )

    level = getattr(logging, str(options.log_level).upper(), logging.INFO)
    root = setup_logging("zmq", console_level=level, console=options.log_stderr)
    install_excepthook()

    if options.log_file:
        handler = logging.handlers.RotatingFileHandler(
            options.log_file, maxBytes=MAX_BYTES, backupCount=BACKUP_COUNT,
            encoding="utf-8", delay=True,
        )
        handler.setLevel(level)
        handler.setFormatter(logging.Formatter(
            "%(asctime)s %(levelname)-7s %(name)s: %(message)s"))
        root.addHandler(handler)


def main(argv: Optional[list[str]] = None) -> int:
    import os

    options, error = options_from_args(argv)
    if error:
        print(f"pyirena-zmq: {error}", file=sys.stderr)
        return 2

    configure_logging(options)

    if options.data_root:
        # api/_paths reads this to confine every file operation.
        os.environ.setdefault("PYIRENA_DATA_ROOT", options.data_root)
    os.environ.setdefault("PYIRENA_MAX_INPUT_POINTS", str(options.max_input_points))

    stop = _Stopping()
    signal.signal(signal.SIGINT, stop)
    signal.signal(signal.SIGTERM, stop)
    return run(options, stop)


if __name__ == "__main__":                                      # pragma: no cover
    raise SystemExit(main())
