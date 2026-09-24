# ZMQ Service — Planning

Internal planning artifact, not user-facing documentation. Treat it as a
starting point to argue with, not a spec. Delete it (git history keeps it)
once the phases it describes have shipped — a shipped plan is worse than no
plan.

Written: 24-09-2026, against pyIrena 1.1.1.
Phase 0 answered by the project owner: 24-09-2026.
**Status: not started; implementation begins 25-09-2026.**

---

## 1. Why this exists

An experiment project runs an AI **orchestrator** on a modest machine that
drives instrument and analysis **services** over ZMQ. pyIrena is to be one of
those services, running on a powerful analysis workstation. The machines are a
mix of operating systems, Windows included, and **share no filesystem**.

The project owner's contract for services is deliberately simple:

- **JSON in, JSON out.** Data arrives inside the request; results leave inside
  the reply. The orchestrator unpacks and stores them itself.
- The orchestrator's **agent** is the caller (it issues the ZMQ requests
  through the orchestrator), not hand-written orchestration code.

The client the orchestrator already uses for its other services is exactly
this, and pyIrena must fit it without modification:

```python
sock = zmq.Context.instance().socket(zmq.REQ)
sock.connect("tcp://%s:%d" % (host, port))
sock.send_string(command)
if poller.poll(timeout):
    return sock.recv_string()
raise ZMQError("no reply within %d ms from %s:%d" % ...)
```

That fixes more of the design than anything else in this document: a **single
UTF-8 text frame** each way (no multipart), plain **REQ/REP**, and the client
gives up after a timeout the owner sets to **60 s**. Everything below is
shaped to always answer inside that window.

Importing pyIrena into the orchestrator would have been simpler, but the
service boundary buys compute placement, crash isolation (a runaway fit
cannot kill the orchestrator), and no numpy/scipy/pyIrena dependency on the
orchestrator machine.

## 2. Decisions taken

| Decision | Rationale |
|---|---|
| The service lives **inside pyIrena** as `pyirena/zmq/`, an optional `[zmq]` extra with a `pyirena-zmq` entry point — the same shape as `[mcp]` | One package to maintain. Explicitly chosen over a companion package. |
| **JSON in / JSON out** for v1; no file paths cross the wire | Owner's requirement; mixed OS, no shared filesystem. |
| **One UTF-8 text frame per message**, JSON document as the payload | Matches the owner's existing client (`send_string` / `recv_string`). No multipart in v1 — see Phase 5 for what that costs binary arrays later. |
| **Plain REQ/REP**, synchronous | Confirmed by the owner. |
| **Default port 9865**, TCP | Chosen by the owner (another service already holds a different port on `purple`). |
| **No async job API in v1** | Confirmed: fits are expected to finish well inside the client's 60 s timeout; anything slower should time out rather than be queued. Phase 4 stays on the shelf. |
| **The server always replies inside a wall-clock budget** (default 55 s) — a fit that overruns gets a `TIMEOUT` reply, not silence | A REQ socket that times out is unusable and must be recreated; the owner's client raises instead. Replying keeps the client's socket alive and tells the agent *why* it failed. |
| **Server platform: Red Hat Linux, any Python we choose** (target 3.12/3.13) | Owner: "we can make any python version on the server." Removes the Windows-service question for the server; the *client* may still be Windows. |
| **Typical payload ≈ 2000 points per curve** | ~100 KB of JSON text. JSON is comfortably adequate; binary framing is not needed and is not worth the multipart break from the owner's client. |
| The service is a **thin transport over the existing dispatcher**, like MCP | `pyirena/mcp/dispatch.py` is already protocol-agnostic, returns JSON dicts, and serves Anthropic-style tool schemas the agent can use directly. |
| **Server options are a first-class, discoverable surface** (§6, Phase 2.6) rather than an afterthought | Owner: "in the long run I hope to support reasonable variation of options for the server… simple useful options should be implemented." Cheap ones ship in v1; the rest get a defined place to land. |

Note this partly overrides the standing ai-agent decision "separate package"
— that decision was about the *standalone AI app* (LLM SDKs, config). The ZMQ
service carries no LLM dependencies and adds only `pyzmq`, so it fits the
extras model instead.

## 3. What already exists and is reused unchanged

- `pyirena/mcp/dispatch.py` — `list_categories`, `list_tools`,
  `describe_tool`, `call_tool` over three schema registries (control,
  calculators, data ops). No `mcp` import.
- `pyirena/api/control/` — stateful fitting sessions for Unified Fit, Sizes,
  Simple Fits, Modeling, WAXS Peak Fit, Carbon model; errors returned as
  `{error, suggestion, code}` dicts, never raised.
- `pyirena/api/calculators` — stateless, no files; usable over JSON as is.
- Layering invariant 3 — api returns JSON-serialisable dicts only.
- Fit-quality metrics, residuals, parameter tables — already JSON.

## 4. Gaps v1 must close

1. **Loading data from arrays.** Every entry point today is
   `open_dataset(file_path)` → `readGenericNXcanSAS`. There is no way to
   create a session from q/I/dI. `Session.file_path` is a required `str`.
2. **Getting results out as JSON.** Results leave via `save_*_fit(session_id,
   output_path)` → NXcanSAS HDF5. Only Unified Fit has a JSON
   `export_fit_report`. The other five fitting tools need an equivalent.
3. **Tools that touch the filesystem** must be hidden from a JSON-only
   service: `open_dataset` (path), all `save_*`, all of `data_ops` (read and
   write files), and the file-reading `api` functions. Image tools return
   `image_path` on the *server's* disk, which is meaningless to the caller.
4. **The dispatcher lives under `mcp/`**, which is the wrong home once two
   transports use it.
5. **Strict JSON.** Invariant 3 is enforced in spirit; NaN/Inf handling is
   done per reader. The service needs `json.dumps(..., allow_nan=False)` to
   succeed on every reply, which no test checks today.
6. **Nothing bounds how long a call runs.** With a synchronous contract and a
   60 s client, a slow fit is now a protocol problem, not just a nuisance
   (§6, Phase 2.7).

## 5. Architecture after v1

```
                orchestrator (any OS)                 analysis workstation (RHEL)
  agent ──> orchestrator ──ZMQ REQ──JSON text──>  pyirena/zmq/server.py  (REP :9865)
                                                       │
                                               pyirena/zmq/protocol.py   (envelope, no zmq import)
                                                       │
                                               pyirena/api/dispatch.py   (moved from mcp/)
                                                       │
                            ┌──────────────────────────┼──────────────────┐
                     api/control (sessions)     api/calculators     (data_ops, file readers:
                                                                      hidden in JSON-only profile)
```

MCP keeps working unchanged: `pyirena/mcp/dispatch.py` becomes a re-export of
`pyirena/api/dispatch.py`.

New layer row for AGENTS.md §2:

| Layer | Path | Responsibility | May import |
|---|---|---|---|
| **zmq** | `pyirena/zmq/` | ZMQ transport exposing `pyirena.api` dispatcher as JSON request/reply | api, pyzmq (lazy) |

## 6. Phased plan

Sizes: **S** ≈ one focused session, **M** ≈ two to three, **L** ≈ more.

### Phase 0 — Confirm the contract with the project owner — **DONE (24-09-2026)**

| Question | Answer | Consequence |
|---|---|---|
| Common request/reply envelope in the project? | Not specified; the client just sends and receives a string | Use §6.2's envelope, carried as one UTF-8 JSON text frame. Revisit if the owner later standardises one. |
| Socket pattern | Plain REQ/REP | As planned; no ROUTER in v1. |
| Port / bind interface | **9865** | Default `tcp://0.0.0.0:9865`; firewall it to the orchestrator host. |
| Longest acceptable call / client timeout | **60 s**; "your fitting should be fast, then timeout" | No async jobs (Phase 4 deferred). Server-side budget of 55 s, always replies. |
| Typical data size | **~2000 points** per curve | JSON text is fine; no binary framing needed. |
| Tool error at envelope level or inside `result`? | Not answered | Keep v1 proposal: `ok: false` with an `error` object. Cheap to change before first deployment; confirm when the orchestrator side is written. |
| Python / OS / service management | **Red Hat; any Python** | Target 3.12–3.13, ship a `systemd` unit example. No Windows-service work for the server. |
| Longer-term intent | "Support reasonable variation of options for the server" | Phase 2.6: a real options surface, discoverable via `ping`/`server_info`, plus per-call overrides. |

Remaining to confirm with the owner *during* implementation, not blocking:
error placement (row 6 above) and whether the orchestrator wants
`analyze` (Phase 3) as its primary entry point.

### Phase 1 — Transport-neutral groundwork inside `pyirena.api` (M)

No `pyzmq` anywhere in this phase. Everything here is useful to scripting and
MCP too, and can ship on its own.

**1.1 Move the dispatcher.** `pyirena/mcp/dispatch.py` → `pyirena/api/dispatch.py`.
Leave a re-export shim at the old path so nothing breaks. Update
`test_mcp_dispatch.py` imports (or keep them, via the shim).

**1.2 Tool metadata: `touches_files`.** Add a per-tool flag (in the dispatcher
registry, derived from a set per source module, not by hand-editing every
schema) marking tools that take or produce server-side paths:
`open_dataset`, `save_*`, every `data_ops` tool, image-returning tools. Add
`list_tools(category, profile="all" | "json_only")` so a transport can ask
only for what it can serve. MCP keeps `profile="all"`.

**1.3 `open_dataset_from_data`** in `api/control` (lifecycle group, next to
`open_dataset`):

```python
open_dataset_from_data(
    q: list[float], intensity: list[float], error: list[float] | None = None,
    label: str = "", is_slit_smeared: bool = False, slit_length: float = 0.0,
) -> dict   # same return shape as open_dataset: {session_id, summary}
```

- Validation, each with its own error code: equal lengths (`SHAPE_MISMATCH`),
  ≥ some minimum points, finite values, q > 0, sorted by q (sort, don't
  reject, and say so in the summary), error > 0 where given, a max-points
  guard (`PYIRENA_MAX_INPUT_POINTS`, default 100 000 — 50× the expected 2000).
- `Session.file_path` becomes `Optional[str]` (`None` = in-memory). Every
  reader of `s.file_path` (`save_*`, summaries, `export_fit_report`) must
  handle `None`; `save_*` without `output_path` returns a clear
  `NO_SOURCE_FILE` error instead of failing on `resolve_safe(None)`.
- Schema entry in `control/schemas.py`; it is a lifecycle tool, so add it to
  `SESSION_LIFECYCLE_NAMES` handling in MCP too (useful there as well).
- dQ / resolution: not stored by `Session` today. Accept it as optional and
  store it, but do not use it until a tool needs it. Document that.

**1.4 JSON result export for all six fitting tools.** Generalise
`unified_fit.export_fit_report(session_id, format="json")` into one generic
`export_results(session_id)` that switches on `session.model_name` — one tool
for the agent to learn. Content:

- model name, pyIrena version, timestamp
- parameter table: name, value, uncertainty (if computed), fixed/free, bounds
- derived quantities (invariant, S/V, volume fraction, … as each tool has)
- χ², reduced χ², fit-quality scalars
- fit Q range and point count
- optionally the model curve and residuals on the data Q grid
  (`include_arrays`, `max_points`; allow `max_points=None` = full)
- **the model's `to_dict()` config**, so a result can later be re-opened in
  the GUI or re-applied (all six fitting tools already have core `to_dict`)
- **Monte-Carlo uncertainties are off by default** (`with_uncertainty=False`).
  They are the one thing that reliably blows a 60 s budget; the agent must ask
  for them, and the docs must say they may time out.

**1.5 Strict-JSON guard.** One helper `api/_json.py: to_strict_json(obj)`
(numpy scalars → python, NaN/Inf → `None`, arrays → lists) and a test that
drives each fitting tool through open → select model → fit → export on
`testData` and asserts `json.dumps(result, allow_nan=False)` succeeds.

**1.6 Timing reality check (new, small but do it first).** Before Phase 2,
measure each of the six tools on a 2000-point curve on the target hardware:
open → select model → fit → export. Record the numbers in this file. If any
default path exceeds ~30 s, either it gets a cheaper default (fewer MC
samples, looser convergence) or Phase 4 comes back onto the table. This is
the assumption the whole synchronous design rests on, so it should be
measured, not assumed.

**Phase 1 exit criteria:** existing test suite green; MCP unchanged for users;
a pure-Python script can fit any of the six tools from arrays and get a JSON
report without touching the disk; timings recorded.

### Phase 2 — `pyirena/zmq/` v1: synchronous REQ/REP (M)

**2.1 Packaging.**

```toml
zmq = [
    # ZMQ service exposing pyirena.api to remote orchestrators (JSON in/out).
    # No Qt, no matplotlib required.
    "pyzmq>=25",
    "pyirena[contrast]",   # so calculators work, as for [mcp]
]
```

Add to `all`. Entry point `pyirena-zmq = "pyirena.zmq.server:main"`.
Conda recipe: add `pyzmq` to the optional run requirements if the recipe
lists extras.

**2.2 Envelope.** One UTF-8 JSON document per frame, both directions.

Request:

```json
{"protocol": "pyirena-zmq/1", "id": "c7f3", "op": "call",
 "tool": "run_fit", "args": {"session_id": "a1b2c3d4"},
 "options": {"include_arrays": false}}
```

Reply:

```json
{"protocol": "pyirena-zmq/1", "id": "c7f3", "ok": true,
 "result": {"...": "..."}, "elapsed_s": 1.84, "server": {"pyirena": "1.2.0"}}
```

Failure (transport *or* tool-level — a result dict containing `"error"` is
promoted to `ok: false`):

```json
{"protocol": "pyirena-zmq/1", "id": "c7f3", "ok": false,
 "error": {"code": "BAD_ARGUMENTS", "message": "...", "suggestion": "..."}}
```

`id` is echoed verbatim and may be absent; with REQ/REP it is for the
caller's logs, not for correlation.

Operations (`op`):

| op | Purpose |
|---|---|
| `ping` | liveness; versions, uptime, open-session count, busy/idle |
| `server_info` | effective options and capabilities (§2.6) — how the agent discovers what this deployment allows |
| `list_categories` | dispatcher, current profile |
| `list_tools` | dispatcher, current profile |
| `describe_tool` | dispatcher — the Anthropic-style schema the agent can use directly |
| `call` | dispatcher `call_tool`, refusing `touches_files` tools with `NOT_AVAILABLE_OVER_ZMQ` |
| `open_dataset` | → `open_dataset_from_data` (lifecycle, top level as in MCP) |
| `close_session`, `list_sessions`, `get_session_summary` | lifecycle |

Transport-level error codes: `BAD_JSON`, `BAD_ENVELOPE`,
`UNSUPPORTED_PROTOCOL`, `UNKNOWN_OP`, `MESSAGE_TOO_LARGE`,
`NOT_AVAILABLE_OVER_ZMQ`, `TIMEOUT`, `SERVER_BUSY`, `INTERNAL_ERROR`
(unexpected exception — logged with traceback server-side, message only on
the wire).

A bare non-JSON string (the owner's other services may send plain commands)
gets `BAD_JSON` with a `suggestion` naming the expected envelope — except
`ping`, which is accepted as a bare word for convenience.

**2.3 Module layout.**

- `pyirena/zmq/__init__.py` — nothing that imports `zmq`.
- `pyirena/zmq/protocol.py` — `handle_request(raw: str | bytes) -> str`.
  Parse, validate envelope, route op, call dispatcher, wrap reply, enforce
  strict JSON, catch everything. **No `zmq` import** — this is where all
  logic and most tests live.
- `pyirena/zmq/options.py` — the options dataclass, defaults, config-file and
  CLI parsing, and the `server_info` payload. Shared by server and tests.
- `pyirena/zmq/server.py` — `main()`: argparse, logging, bind a REP socket,
  loop `recv → handle_request (with deadline) → send`. Lazy `import zmq`
  with a friendly "install pyirena[zmq]" message.
- `pyirena/zmq/client.py` — a ~40-line reference client
  (`PyIrenaClient(address, timeout_ms=60000).call(tool, **args)`) with the
  "lazy pirate" reconnect (a REQ socket that times out must be closed and
  recreated — the owner's snippet raises instead, so note this in the docs).
  Ships so the orchestrator team has a working example; it depends only on
  `pyzmq`.

**2.4 Server behaviour.**

- `--bind` default `tcp://0.0.0.0:9865` (`--port` as a shorthand). The
  service is useless on loopback, so binding wide is the default; there is no
  auth in v1, so the deployment note is "firewall 9865 to the orchestrator
  host." The startup log states the effective bind loudly.
- Poll loop with a short timeout (e.g. 250 ms) so Ctrl+C works.
- SIGTERM / Ctrl+C: finish the current request, close socket, exit 0.
- Logging to a rotating file (`--log-file`, default in a user log dir) and
  optionally stderr; **never stdout**. One line per request: id, op, tool,
  elapsed, ok/error code. Reuse `pyirena/logging_setup.py` if it fits.
- `zmq.MAXMSGSIZE` from `--max-message-mb` (default 16 — 2000 points is
  ~0.1 MB, so this is already 100× headroom).
- Session hygiene: idle-session TTL (`--session-ttl-min`, default 60) and a
  cap on open sessions (`--max-sessions`, default 32), evicting oldest idle.
  Checked between requests. Reported in `ping`.
- Images: image-returning tools are hidden by default; `--allow-images`
  switches them to base64-in-JSON (§2.6).
- Single-threaded request handling by design: one request at a time. The
  session registry is a plain dict and scipy fits are CPU-bound, so serial
  execution is also the safe choice. The work itself runs on a watchdog
  thread only so the deadline in §2.7 can be enforced; there is still exactly
  one reply per request and never two fits at once.
- A request arriving while another is still running is impossible with a
  single REQ client; with two clients the second simply waits. If a second
  client's request is picked up while the deadline thread is still cleaning
  up, reply `SERVER_BUSY` rather than interleaving.

**2.5 Agent-facing notes.** The orchestrator can fetch `describe_tool`
schemas and hand them to its LLM as tools, wrapping each call into a `call`
envelope. A short "how an agent should use this service" section (open from
arrays → select model → set/fix parameters → set Q range → fit → check
quality → export results → close) goes into the user doc; it is the same
workflow MCP agents follow today, minus files. It must also say: **one
session per curve, close it when done**, and **`analyze` (Phase 3) is the
one-call shortcut**.

**2.6 Server options — the "reasonable variation" surface (S, ships in v1).**

The owner asked for this explicitly, so it gets designed once rather than
accreting flags. Three layers, each overriding the one above:

1. **Config file** — `--config <file>` (TOML), or
   `$PYIRENA_ZMQ_CONFIG`. Same keys as the CLI flags, underscored.
2. **CLI flags** — override the file.
3. **Per-call `options` in the envelope** — override for that request only,
   and only within what the server allows (a request cannot turn on
   `allow_files`; it can turn off `include_arrays`).

| Option | Default | Per-call? | Notes |
|---|---|---|---|
| `bind` / `port` | `tcp://0.0.0.0:9865` | no | |
| `profile` | `json_only` | no | `all` requires `allow_files` |
| `max_message_mb` | 16 | no | |
| `request_budget_s` | 55 | yes (may only shorten) | §2.7 |
| `max_input_points` | 100 000 | no | |
| `session_ttl_min` | 60 | no | |
| `max_sessions` | 32 | no | |
| `include_arrays` | false | **yes** | model curve + residuals in `export_results` |
| `max_points` | 2000 | **yes** | decimation cap for returned arrays |
| `with_uncertainty` | false | **yes** | MC uncertainties; may hit the budget |
| `allow_images` | false | yes (off only) | base64 PNG instead of `image_path` |
| `allow_files` | false | no | opens `open_dataset`, `save_*`, `data_ops` under `data_root` |
| `data_root` | unset | no | required when `allow_files` |
| `log_file`, `log_level`, `log_stderr` | rotating file, INFO, off | no | |

`server_info` returns the effective values plus a `capabilities` list, so the
agent discovers what a given deployment permits instead of guessing — and so
the same orchestrator code works against a locked-down and a permissive
server. Unknown keys in `options` are an error, not silently ignored:
silently dropping an option the agent thinks it set is worse than failing.

**2.7 The request deadline (S, but it is what makes the sync contract safe).**

The owner's client waits 60 s and then raises; its REQ socket is then dead and
must be recreated. So the server must never be the reason the client times
out. The dispatcher call runs on a worker thread with
`future.result(timeout=request_budget_s)`; on overrun the server replies:

```json
{"ok": false, "error": {"code": "TIMEOUT",
 "message": "run_fit exceeded the 55 s server budget",
 "suggestion": "narrow the Q range, fix parameters, or set with_uncertainty=false"}}
```

The orphaned computation is left to finish and its result discarded; the
session it was fitting is marked `stale` and reported as such by
`get_session_summary`, because its model state is mid-fit and untrustworthy.
The next request on a stale session gets a clear error telling the agent to
re-run or close it. Simple, honest, and it never leaves the client guessing.

**Phase 2 exit criteria:** from the orchestrator machine, the reference client
opens a dataset from arrays, fits it with each of the six tools, exports JSON
results and closes the session; the server survives malformed input, bare
strings, oversize messages, unknown tools, a killed client, and a fit that
overruns the budget.

### Phase 3 — One-shot "recipe" call (S–M, now recommended, not optional)

With a synchronous 60 s contract, every extra round trip costs latency and
another chance to strand a session. A single coarse tool for the common case:

```python
analyze(data: {q, intensity, error, ...}, tool: str, config: dict,
        include_arrays: bool = False) -> dict   # the export_results() payload
```

`config` is the same JSON the GUI's *Export Parameters* writes and
`pyirena.batch` reads, so a scientist can set up a fit in the GUI once and the
agent can replay it on every new measurement. Implemented in `api` on top of
Phase 1 (open from data → build model from config → fit → export → close),
**not** by calling `pyirena.batch` (stdout logging, `None` on failure, file
I/O). Reuse batch's config-to-model helpers (`batch/unified.py:_state_to_model`
and siblings) by moving them to `core`/`api` where needed. Exposed through the
dispatcher, therefore automatically through MCP too. It opens and closes its
own session, so it cannot leak one, and it is the path the docs should lead
with; the fine-grained session tools stay for hard fits that need iteration.

### Phase 4 — Asynchronous jobs (M) — **deferred by the owner**

Not in scope. Revisit only if Phase 1.6's timings, or real use, show fits that
genuinely need more than the budget. Sketch kept for that day:

- Switch the server to ROUTER; a job table plus a worker (thread first;
  process pool if the GIL bites — note sessions then must live in the worker).
- New ops: `submit` (returns `job_id` immediately), `status`, `result`,
  `cancel`, `list_jobs`. `ping` answers during fits.
- Cancellation needs a cooperative check in the fit loops (a callback or
  flag checked per iteration); scope that per tool.
- Optional PUB socket for progress / completion events.
- Keep synchronous `call` working for short operations; clients choose.

The §2.7 deadline is deliberately the *same seam* a job API would use, so this
is an extension, not a rewrite.

### Phase 5 — Richer data paths (deferred; pick by need)

| Option | When | Sketch |
|---|---|---|
| **File paths on a shared filesystem** | If some deployment does share storage | Already reserved as `allow_files` + `data_root` in §2.6; flipping it switches the profile to `all`. Mostly already built. |
| **Zarr input** | When the project's Zarr layout stabilises | `pyirena/io/zarr_sas.py` behind a `[zarr]` extra, returning the same dict as `readGenericNXcanSAS`; an `open_dataset` that dispatches on file type. Layout mapping configurable, not hard-coded. |
| **Binary arrays** | If curves get far larger than 2000 points, or many per call | Multipart: JSON header + raw little-endian float64 frames. **This breaks the owner's `send_string`/`recv_string` client**, so it needs their buy-in, not just ours. A base64 array block inside the JSON is the compatible middle road if it ever comes up. |
| **Images over the wire** | If the agent is multimodal | Already reserved as `allow_images` in §2.6: `image_base64`, no `image_path`. |
| **Results written server-side** | If the orchestrator wants NXcanSAS sidecars | `save_*` with an `output_path` under `data_root`; return path + JSON. |
| **Authentication / encryption** | If the port leaves a trusted subnet | ZMQ CURVE with keys in a config file; `--curve-keys`. Until then, firewall. |
| **Multiple workers / clients** | Several orchestrators or parallel fits | ROUTER–DEALER broker, session affinity per worker. |

## 7. Effect on existing users

- Nothing changes without `pip install pyirena[zmq]`; `import pyirena`
  must still work with no extras (`test_optional_dep_modules.py`).
- MCP: the dispatcher move is invisible through the shim. New tools
  (`open_dataset_from_data`, `export_results`, `analyze`) appear to MCP
  clients as well — additive only.
- `Session.file_path` becoming optional touches every `save_*`; covered by
  `test_control_save.py` plus new in-memory-session cases.
- No GUI changes. No change to HDF5 formats.

## 8. Tests

- `pyirena/tests/api/test_open_from_data.py` — validation cases, sorting,
  in-memory session through a full fit for each tool, `save_*` without path.
- `pyirena/tests/api/test_export_results.py` — per tool; strict-JSON check;
  `to_dict` config round-trips back into a model.
- `pyirena/tests/api/test_dispatch_profiles.py` — `json_only` hides every
  `touches_files` tool; the set is non-empty and stable.
- `pyirena/tests/zmq/test_protocol.py` — `handle_request` without pyzmq:
  bad JSON, bare string, bad envelope, unknown op/tool, tool error promotion,
  internal exception, oversize, unknown `options` key, strict JSON of every
  reply.
- `pyirena/tests/zmq/test_options.py` — precedence file < CLI < per-call;
  per-call cannot widen a server restriction; `server_info` matches effective
  options.
- `pyirena/tests/zmq/test_deadline.py` — a deliberately slow fake tool returns
  `TIMEOUT` within the budget and marks the session stale.
- `pyirena/tests/zmq/test_server_roundtrip.py` —
  `pytest.importorskip("zmq")`; server in a thread on a random port, driven
  **by a raw `send_string`/`recv_string` REQ socket** (the owner's client
  shape, not only ours) through open → fit → export → close; plus client
  timeout and reconnect.
- Layering: `pyirena/zmq` must not import Qt or matplotlib at module level;
  extend the grep in AGENTS.md invariant 1 and the Qt contract test.
- Manual: one real cross-machine run against the orchestrator.

## 9. Documentation to touch when implementing

- `docs/zmq_service.md` (new, user-facing): install, run, flags and config
  file, envelope, op table, error codes, reference client, agent workflow,
  the 60 s budget and what `TIMEOUT` means, systemd unit, firewall note,
  limits.
- `AGENTS.md`: layer-table row (§5 above), `pyirena-zmq` in commands,
  "Where to look" row, invariant 1 grep includes `pyirena/zmq`.
- `docs/module_map.md`, `docs/developer_adding_features.md` (new fitting tool
  → needs `export_results` support and a `touches_files` decision),
  `pyirena/api/README.md`, `docs/ai_integration.md`, `CHANGELOG.md`.

## 10. Open questions

1. ~~Envelope and socket pattern~~ — settled (Phase 0).
2. ~~Is the longest realistic fit shorter than the orchestrator's timeout?~~ —
   assumed yes; **Phase 1.6 measures it** rather than leaving it assumed.
3. ~~`export_results` generic vs one per tool~~ — generic, confirmed here.
4. ~~Should `analyze` be the primary documented path?~~ — yes; session tools
   are the fallback for hard fits.
5. ~~Uncertainties in `export_results`~~ — opt-in only (`with_uncertainty`),
   because of the 60 s budget.
6. Tool-level errors: `ok: false` (current plan) or `ok: true` with an error
   inside `result`? Owner did not say; confirm before the orchestrator side
   is written.
7. Does anything else in the project (other services) already solve
   logging/service management on the RHEL box, so pyIrena should match it?
8. Is 9865 reachable from the orchestrator through the beamline firewall, and
   who owns that rule?
