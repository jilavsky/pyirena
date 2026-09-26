#!/usr/bin/env python3
"""Smoke-test a running pyIrena ZMQ service, from anywhere.

Run this on the *client* machine — the one the orchestrator runs on — to
check that a deployed service is reachable and behaving. It needs **pyzmq and
nothing else**: no pyIrena, no numpy, no h5py, so it works on a bare
orchestrator box.

    python zmq_smoke_test.py tcp://workstation:9865
    python zmq_smoke_test.py tcp://workstation:9865 --config modeling.json \
                             --data curve.json

Every check prints PASS or FAIL and the script exits non-zero if any failed,
so it can go straight into a deployment script.

``--config`` additionally replays a real exported pyIrena configuration
through ``analyze``. Without ``--data`` it invents a curve, which exercises
the transport but tells you nothing about the science; pass a JSON file of
``{"q": [...], "intensity": [...], "error": [...]}`` to use real data.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
import uuid

try:
    import zmq
except ModuleNotFoundError:
    sys.exit("This script needs pyzmq:  pip install pyzmq")

PROTOCOL = "pyirena-zmq/1"

_passed = 0
_failed = 0


def check(label: str, ok: bool, detail: str = "") -> bool:
    global _passed, _failed
    if ok:
        _passed += 1
        print(f"  PASS  {label}" + (f"   ({detail})" if detail else ""))
    else:
        _failed += 1
        print(f"  FAIL  {label}" + (f"   ({detail})" if detail else ""))
    return ok


class Client:
    """A REQ client that rebuilds its socket after a timeout (lazy pirate)."""

    def __init__(self, address: str, timeout_ms: int = 60_000):
        self.address = address
        self.timeout_ms = timeout_ms
        self._connect()

    def _connect(self):
        self.ctx = zmq.Context.instance()
        self.sock = self.ctx.socket(zmq.REQ)
        self.sock.setsockopt(zmq.LINGER, 0)
        self.sock.connect(self.address)
        self.poller = zmq.Poller()
        self.poller.register(self.sock, zmq.POLLIN)

    def send(self, op, timeout_ms=None, **fields):
        envelope = {"protocol": PROTOCOL, "id": uuid.uuid4().hex[:8], "op": op}
        envelope.update({k: v for k, v in fields.items() if v is not None})
        return self.raw(json.dumps(envelope), timeout_ms)

    def raw(self, text, timeout_ms=None):
        self.sock.send_string(text)
        if not dict(self.poller.poll(timeout_ms or self.timeout_ms)):
            # A timed-out REQ socket can never be used again — rebuild it.
            self.sock.close()
            self.poller.unregister(self.sock)
            self._connect()
            raise TimeoutError(f"no reply within {timeout_ms or self.timeout_ms} ms")
        return json.loads(self.sock.recv_string())

    def close(self):
        self.sock.close()


def synthetic_curve(n=500):
    """A Guinier + Porod curve. Pure Python so this script needs no numpy."""
    q, intensity, error = [], [], []
    for i in range(n):
        qi = 10 ** (-3 + 3 * i / (n - 1))
        val = 1000.0 * math.exp(-(qi**2) * 150.0**2 / 3) + 1e-6 * qi**-4 + 0.01
        q.append(qi)
        intensity.append(val)
        error.append(0.02 * val)
    return {"q": q, "intensity": intensity, "error": error, "label": "smoke test"}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("address", help="e.g. tcp://workstation:9865")
    p.add_argument("--timeout-ms", type=int, default=60_000,
                   help="client timeout; must exceed the server's budget (default 60000)")
    p.add_argument("--config", help="a pyIrena config JSON to replay through analyze")
    p.add_argument("--data", help='JSON {"q": [...], "intensity": [...], "error": [...]}')
    p.add_argument("--points", type=int, default=500,
                   help="points in the synthetic curve (default 500)")
    args = p.parse_args(argv)

    print(f"\npyIrena ZMQ smoke test -> {args.address}\n")
    c = Client(args.address, args.timeout_ms)

    # --- 1. reachable ------------------------------------------------------
    print("1. Reachability")
    try:
        reply = c.send("ping", timeout_ms=5000)
    except TimeoutError:
        check("server answers ping", False,
              "no reply in 5 s — wrong host/port, service down, or firewall")
        print("\nCannot continue without a reachable server.")
        return 1
    check("server answers ping", reply.get("ok") is True)
    result = reply.get("result", {})
    print(f"        uptime {result.get('uptime_s')}s, "
          f"{result.get('open_sessions')} open session(s)")
    server = reply.get("server", {})
    print(f"        pyirena {server.get('pyirena')}, protocol {server.get('protocol')}")

    # A bare string, as the other services in the project send.
    try:
        check("bare 'ping' string is accepted", c.raw("ping", 5000).get("ok") is True)
    except TimeoutError:
        check("bare 'ping' string is accepted", False, "timed out")

    # --- 2. what this deployment allows ------------------------------------
    print("\n2. Configuration")
    info = c.send("server_info")["result"]
    opts = info["options"]
    print(f"        profile={opts['profile']}  budget={opts['request_budget_s']}s  "
          f"max_message={opts['max_message_mb']}MB  max_sessions={opts['max_sessions']}")
    check("client timeout exceeds the server budget",
          args.timeout_ms / 1000.0 > opts["request_budget_s"],
          f"client {args.timeout_ms/1000:.0f}s vs server {opts['request_budget_s']}s")
    check("tools are exposed", info["tool_count"] > 0, f"{info['tool_count']} tools")
    categories = [cat["name"] for cat in info["categories"]]
    print(f"        categories: {', '.join(categories)}")

    # --- 3. discovery ------------------------------------------------------
    print("\n3. Discovery")
    tools = c.send("list_tools", args={"category": "results"})["result"]["tools"]
    check("'analyze' is listed", "analyze" in [t["name"] for t in tools])
    schema = c.send("describe_tool", args={"name": "run_fit"})["result"]
    check("describe_tool returns a usable schema",
          isinstance(schema.get("input_schema"), dict))

    # --- 4. a real fit, step by step ---------------------------------------
    print("\n4. A fit, the step-by-step way")
    curve = json.load(open(args.data)) if args.data else synthetic_curve(args.points)
    n = len(curve["q"])
    payload = json.dumps({"op": "open_dataset", "args": curve})
    print(f"        {n} points, request {len(payload)/1024:.1f} KB")

    session = None
    try:
        opened = c.send("open_dataset", args=curve)
        if check("open_dataset from arrays", opened.get("ok") is True,
                 opened.get("error", {}).get("message", "")[:70]):
            session = opened["result"]["session_id"]
            summary = opened["result"]["summary"]
            print(f"        session {session}, {summary['n_points']} points kept")

        if session:
            for tool, kwargs in [("select_model", {"model_name": "unified_fit"}),
                                 ("add_unified_level", {})]:
                c.send("call", tool=tool, args={"session_id": session, **kwargs})
            t0 = time.time()
            fit = c.send("call", tool="run_fit", args={"session_id": session})
            elapsed = time.time() - t0
            check("run_fit", fit.get("ok") is True,
                  f"{elapsed:.2f}s" if fit.get("ok") else
                  fit.get("error", {}).get("message", "")[:70])

            exported = c.send("call", tool="export_results",
                              args={"session_id": session})
            check("export_results", exported.get("ok") is True)
            if exported.get("ok"):
                q = exported["result"]["quality"]
                print(f"        chi2={q.get('chi_squared')}  "
                      f"reduced={q.get('reduced_chi_squared')}")
    finally:
        if session:
            closed = c.send("close_session", args={"session_id": session})
            check("close_session", closed.get("ok") is True)

    check("no sessions left open",
          c.send("ping")["result"]["open_sessions"] == 0)

    # --- 5. one-call analyze -----------------------------------------------
    if args.config:
        print("\n5. One call: analyze with a real config")
        config = json.load(open(args.config))
        t0 = time.time()
        reply = c.send("call", tool="analyze",
                       args={"data": curve, "config": config})
        elapsed = time.time() - t0
        if check("analyze", reply.get("ok") is True,
                 f"{elapsed:.2f}s" if reply.get("ok")
                 else f"{reply.get('error',{}).get('code')}: "
                      f"{reply.get('error',{}).get('message','')[:60]}"):
            r = reply["result"]
            q = r["quality"]
            print(f"        tool={r['tool']}  chi2={q.get('chi_squared')}  "
                  f"reduced={q.get('reduced_chi_squared')}")
            print(f"        fitted Q range: {r['analyze']['fit_q_min']} .. "
                  f"{r['analyze']['fit_q_max']}")
            for note in r["analyze"]["notes"]:
                print(f"        note: {note}")
            check("analyze left no session behind",
                  c.send("ping")["result"]["open_sessions"] == 0)
            check("analyze stayed inside the budget",
                  elapsed < opts["request_budget_s"], f"{elapsed:.1f}s")
    else:
        print("\n5. One call: analyze        (skipped — pass --config to test it)")

    # --- 6. the service survives bad input ---------------------------------
    print("\n6. Robustness")
    for label, payload, expected in [
        ("malformed JSON",        "not json at all",              "BAD_JSON"),
        ("unknown op",            '{"op": "teleport"}',           "UNKNOWN_OP"),
        ("unknown tool",          '{"op": "call", "tool": "nope"}', "UNKNOWN_TOOL"),
        ("a file tool is refused",
         '{"op": "call", "tool": "save_fit", "args": {"session_id": "x"}}',
         "NOT_AVAILABLE_OVER_ZMQ"),
    ]:
        try:
            reply = c.raw(payload, 10_000)
            code = reply.get("error", {}).get("code")
            check(f"{label} -> {expected}", code == expected, f"got {code}")
        except TimeoutError:
            check(f"{label} -> {expected}", False, "timed out (server may have died)")

    try:
        check("still healthy after bad input", c.send("ping", 5000).get("ok") is True)
    except TimeoutError:
        check("still healthy after bad input", False, "server stopped answering")

    c.close()
    print(f"\n{'-' * 56}\n{_passed} passed, {_failed} failed\n")
    return 1 if _failed else 0


if __name__ == "__main__":
    sys.exit(main())
