# ZMQ Service — Planning

Internal planning artifact, not user-facing documentation. Treat it as a
starting point to argue with, not a spec. Delete it (git history keeps it)
once the phases it describes have shipped — a shipped plan is worse than no
plan.

Written: 24-09-2026, against pyIrena 1.1.1. **Status: not started.**

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

Importing pyIrena into the orchestrator would have been simpler, but the
service boundary buys compute placement, crash isolation (a runaway fit
cannot kill the orchestrator), and no numpy/scipy/pyIrena dependency on the
orchestrator machine.

## 2. Decisions taken

| Decision | Rationale |
|---|---|
| The service lives **inside pyIrena** as `pyirena/zmq/`, an optional `[zmq]` extra with a `pyirena-zmq` entry point — the same shape as `[mcp]` | One package to maintain. Explicitly chosen over a companion package. |
| **JSON in / JSON out** for v1; no file paths cross the wire | Owner's requirement; mixed OS, no shared filesystem. |
| **No Zarr support in v1** | The Zarr layout is still changing; JSON arrays sidestep it entirely. |
| The service is a **thin transport over the existing dispatcher**, like MCP | `pyirena/mcp/dispatch.py` is already protocol-agnostic, returns JSON dicts, and serves Anthropic-style tool schemas the agent can use directly. |
| Richer options (async jobs, file/Zarr input, binary arrays, auth) are **planned but deferred** | See Phases 4–5 in §6. v1 must not paint them into a corner. |

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

## 5. Architecture after v1

```
                orchestrator (any OS)                 analysis workstation
  agent ──> orchestrator ──ZMQ REQ──JSON──>  pyirena/zmq/server.py  (REP socket)
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

### Phase 0 — Confirm the contract with the project owner (no code)

Settle before writing code; each answer changes something below.

- [ ] Envelope: does the project have a common request/reply shape (field
      names for id, operation, status, error)? If yes, adopt it verbatim in
      `protocol.py`; if not, use §6.2's.
- [ ] Socket pattern: plain REQ/REP acceptable? (v1 assumes yes.)
- [ ] Port and bind interface on the workstation; firewall owner.
- [ ] Longest acceptable call duration / the orchestrator's timeout. Decides
      whether Phase 4 (async jobs) is needed before first real use.
- [ ] Typical data size (points per curve, curves per request).
- [ ] Should a tool-level error (bad parameter) be `ok: false` at the
      envelope level, or `ok: true` with an error dict inside `result`?
      (v1 proposes `ok: false`; see §6.2.)
- [ ] Python version and OS of the workstation; how the service is kept
      running (systemd, Windows service/NSSM, a terminal).

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
  guard (`PYIRENA_MAX_INPUT_POINTS`, generous default).
- `Session.file_path` becomes `Optional[str]` (`None` = in-memory). Every
  reader of `s.file_path` (`save_*`, summaries, `export_fit_report`) must
  handle `None`; `save_*` without `output_path` returns a clear
  `NO_SOURCE_FILE` error instead of failing on `resolve_safe(None)`.
- Schema entry in `control/schemas.py`; it is a lifecycle tool, so add it to
  `SESSION_LIFECYCLE_NAMES` handling in MCP too (useful there as well).
- dQ / resolution: not stored by `Session` today. Accept it as optional and
  store it, but do not use it until a tool needs it. Document that.

**1.4 JSON result export for all six fitting tools.** Generalise
`unified_fit.export_fit_report(session_id, format="json")` into one function
per tool (or one generic `export_results(session_id)` that switches on
`session.model_name` — preferred, one tool for the agent to learn). Content:

- model name, pyIrena version, timestamp
- parameter table: name, value, uncertainty (if computed), fixed/free, bounds
- derived quantities (invariant, S/V, volume fraction, … as each tool has)
- χ², reduced χ², fit-quality scalars
- fit Q range and point count
- optionally the model curve and residuals on the data Q grid
  (`include_arrays`, `max_points`; allow `max_points=None` = full)
- **the model's `to_dict()` config**, so a result can later be re-opened in
  the GUI or re-applied (all six fitting tools already have core `to_dict`)

**1.5 Strict-JSON guard.** One helper `api/_json.py: to_strict_json(obj)`
(numpy scalars → python, NaN/Inf → `None`, arrays → lists) and a test that
drives each fitting tool through open → select model → fit → export on
`testData` and asserts `json.dumps(result, allow_nan=False)` succeeds.

**Phase 1 exit criteria:** existing test suite green; MCP unchanged for users;
a pure-Python script can fit any of the six tools from arrays and get a JSON
report without touching the disk.

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

**2.2 Envelope (proposal — replace with the project's if Phase 0 gives one).**

Request:

```json
{"protocol": "pyirena-zmq/1", "id": "c7f3", "op": "call",
 "tool": "run_fit", "args": {"session_id": "a1b2c3d4"}}
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

Operations (`op`):

| op | Purpose |
|---|---|
| `ping` | liveness; returns versions, uptime, open-session count |
| `list_categories` | dispatcher, `profile="json_only"` |
| `list_tools` | dispatcher, `profile="json_only"` |
| `describe_tool` | dispatcher — the Anthropic-style schema the agent can use directly |
| `call` | dispatcher `call_tool`, refusing `touches_files` tools with `NOT_AVAILABLE_OVER_ZMQ` |
| `open_dataset` | → `open_dataset_from_data` (lifecycle, top level as in MCP) |
| `close_session`, `list_sessions`, `get_session_summary` | lifecycle |

Transport-level error codes: `BAD_JSON`, `BAD_ENVELOPE`,
`UNSUPPORTED_PROTOCOL`, `UNKNOWN_OP`, `MESSAGE_TOO_LARGE`,
`NOT_AVAILABLE_OVER_ZMQ`, `INTERNAL_ERROR` (unexpected exception — logged
with traceback server-side, message only on the wire).

**2.3 Module layout.**

- `pyirena/zmq/__init__.py` — nothing that imports `zmq`.
- `pyirena/zmq/protocol.py` — `handle_request(raw: bytes) -> bytes`. Parse,
  validate envelope, route op, call dispatcher, wrap reply, enforce strict
  JSON, catch everything. **No `zmq` import** — this is where all logic and
  most tests live.
- `pyirena/zmq/server.py` — `main()`: argparse, logging, bind a REP socket,
  loop `recv → handle_request → send`. Lazy `import zmq` with a friendly
  "install pyirena[zmq]" message.
- `pyirena/zmq/client.py` — a ~40-line reference client
  (`PyIrenaClient(address).call(tool, **args)`) with timeout and the "lazy
  pirate" reconnect (a REQ socket that times out must be closed and
  recreated). Ships so the orchestrator team has a working example; it
  depends only on `pyzmq`.

**2.4 Server behaviour.**

- `--bind` default `tcp://127.0.0.1:5577`; binding to the network
  (`tcp://0.0.0.0:5577` or an interface address) must be explicit.
- Poll loop with a short timeout (e.g. 250 ms) so Ctrl+C works on Windows
  (a blocking `recv` ignores it there).
- SIGTERM / Ctrl+C: finish the current request, close socket, exit 0.
- Logging to a rotating file (`--log-file`, default in a user log dir) and
  optionally stderr; **never stdout**. One line per request: id, op, tool,
  elapsed, ok/error code. Reuse `pyirena/logging_setup.py` if it fits.
- `zmq.MAXMSGSIZE` from `--max-message-mb` (default e.g. 50).
- Session hygiene: idle-session TTL (`--session-ttl-min`, default 60) and a
  cap on open sessions (`--max-sessions`), evicting oldest idle. Checked
  between requests. Report evictions in `ping`.
- Images: image-returning tools are hidden in v1 (Phase 5 lists the option to
  return base64 PNGs on request).
- Single-threaded by design: one request at a time. The session registry is
  a plain dict and scipy fits are CPU-bound, so serial execution is also the
  safe choice. Consequence: `ping` is not answered while a fit runs — the
  client timeout must exceed the longest fit. Phase 4 removes this.

**2.5 Agent-facing notes.** The orchestrator can fetch `describe_tool`
schemas and hand them to its LLM as tools, wrapping each call into a `call`
envelope. A short "how an agent should use this service" section (open from
arrays → select model → set/fix parameters → set Q range → fit → check
quality → export results → close) goes into the user doc; it is the same
workflow MCP agents follow today, minus files.

**Phase 2 exit criteria:** from a second machine (ideally Windows), the
reference client opens a dataset from arrays, fits it with each of the six
tools, exports JSON results and closes the session; the server survives
malformed input, oversize messages, unknown tools and a killed client.

### Phase 3 — One-shot "recipe" call (S–M, optional but recommended)

A single coarse tool for the common case, where the agent (or plain
orchestrator code) already knows what it wants:

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
dispatcher, therefore automatically through MCP too.

### Phase 4 — Asynchronous jobs (M) — when fits outgrow the client timeout

- Switch the server to ROUTER; a job table plus a worker (thread first;
  process pool if the GIL bites — note sessions then must live in the worker).
- New ops: `submit` (returns `job_id` immediately), `status`, `result`,
  `cancel`, `list_jobs`. `ping` answers during fits.
- Cancellation needs a cooperative check in the fit loops (a callback or
  flag checked per iteration); scope that per tool.
- Optional PUB socket for progress / completion events.
- Keep synchronous `call` working for short operations; clients choose.

### Phase 5 — Richer data paths (deferred; pick by need)

| Option | When | Sketch |
|---|---|---|
| **File paths on a shared filesystem** | If some deployment does share storage | `--allow-files` + `PYIRENA_DATA_ROOT`; switches the profile back to `all`, exposing `open_dataset`, `save_*`, `data_ops`. Mostly already built. |
| **Zarr input** | When the project's Zarr layout stabilises | `pyirena/io/zarr_sas.py` behind a `[zarr]` extra, returning the same dict as `readGenericNXcanSAS`; an `open_dataset` that dispatches on file type. Layout mapping configurable, not hard-coded. |
| **Binary arrays** | If JSON size/latency matters (≫10⁴ points or many curves per call) | Multipart messages: JSON header + raw little-endian float64 frames referenced by index. JSON stays the default. |
| **Images over the wire** | If the agent is multimodal and benefits from seeing fits | `include_images: true` → keep `image_base64`, drop `image_path`. |
| **Results written server-side** | If the orchestrator wants NXcanSAS sidecars | `save_*` with an `output_path` under a server data root; return path + JSON. |
| **Authentication / encryption** | If the port leaves a trusted subnet | ZMQ CURVE with keys in a config file; `--curve-keys`. |
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
  bad JSON, bad envelope, unknown op/tool, tool error promotion, internal
  exception, oversize, strict JSON of every reply.
- `pyirena/tests/zmq/test_server_roundtrip.py` —
  `pytest.importorskip("zmq")`; server in a thread on a random port, the
  reference client drives open → fit → export → close; client timeout and
  reconnect.
- Layering: `pyirena/zmq` must not import Qt or matplotlib at module level;
  extend the grep in AGENTS.md invariant 1 and the Qt contract test.
- Manual: one real cross-machine run, Windows client ↔ workstation server.

## 9. Documentation to touch when implementing

- `docs/zmq_service.md` (new, user-facing): install, run, flags, envelope,
  op table, error codes, reference client, agent workflow, Windows and
  firewall notes, limits.
- `AGENTS.md`: layer-table row (§5 above), `pyirena-zmq` in commands,
  "Where to look" row, invariant 1 grep includes `pyirena/zmq`.
- `docs/module_map.md`, `docs/developer_adding_features.md` (new fitting tool
  → needs `export_results` support and a `touches_files` decision),
  `pyirena/api/README.md`, `docs/ai_integration.md`, `CHANGELOG.md`.

## 10. Open questions

1. Envelope and socket pattern — adopt the project's if it has one (Phase 0).
2. Is the longest realistic fit shorter than the orchestrator's timeout? If
   not, Phase 4 moves ahead of first deployment.
3. `export_results` as one generic tool vs one per fitting tool — generic
   preferred; confirm when implementing.
4. Should `analyze` (Phase 3) be the *primary* documented path for the agent,
   with the fine-grained session tools as the fallback for hard fits?
5. Uncertainties: include Monte-Carlo parameter uncertainties in
   `export_results` (slow) or only on request (`with_uncertainty=True`)?
6. Does anything else in the project (other services) already solve
   logging/service management, so pyIrena should match it?
