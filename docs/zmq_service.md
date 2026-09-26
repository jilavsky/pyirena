# pyIrena ZMQ service

`pyirena-zmq` serves pyIrena's fitting tools to a remote caller over ZMQ:
**JSON in, JSON out, one request one reply**. It exists for the case where the
caller — an experiment orchestrator, an agent, a reduction pipeline — runs on a
different machine from the analysis workstation and **shares no filesystem with
it**. The data travels inside the request; the results travel inside the reply.

If your caller is on the same machine and can read the same files, you
probably want [the MCP server](ai_integration.md) instead. Same tools, same
schemas; MCP can also open and save files.

**Deploying or testing it for the first time?** Work through
[zmq_service_testing.md](zmq_service_testing.md) — the same material as a
tickable procedure, with a smoke-test script that needs only pyzmq.

---

## Install and run

```bash
pip install 'pyirena[zmq]'
pyirena-zmq                        # binds tcp://0.0.0.0:9865
pyirena-zmq --port 9999 --log-stderr
pyirena-zmq --config /etc/pyirena/zmq.json
```

The service listens on **port 9865** by default, on every interface, **with no
authentication**. Restrict that port to the orchestrator's host at the
firewall. On the analysis workstation:

```bash
sudo firewall-cmd --permanent --add-rich-rule \
  'rule family="ipv4" source address="10.0.0.5" port port="9865" protocol="tcp" accept'
sudo firewall-cmd --reload
```

Logs go to `~/.pyirena/logs/zmq.log` (rotating), never to stdout. One line per
request: id, op, tool, session, outcome, elapsed.

### Keeping it running (systemd)

```ini
# /etc/systemd/system/pyirena-zmq.service
[Unit]
Description=pyIrena ZMQ analysis service
After=network.target

[Service]
Type=simple
User=pyirena
ExecStart=/opt/pyirena/bin/pyirena-zmq --config /etc/pyirena/zmq.json
Restart=on-failure
RestartSec=5

[Install]
WantedBy=multi-user.target
```

```bash
sudo systemctl enable --now pyirena-zmq
journalctl -u pyirena-zmq -f
```

---

## The protocol

One **UTF-8 JSON document per message**, both directions, on a plain REQ/REP
socket. No multipart, no binary frames — `send_string` and `recv_string` are
all a client needs.

### Request

```json
{"protocol": "pyirena-zmq/1", "id": "c7f3", "op": "call",
 "tool": "run_fit", "args": {"session_id": "a1b2c3d4"},
 "options": {"include_arrays": true}}
```

| Field | |
|---|---|
| `op` | **required** — the operation (table below) |
| `tool` | the tool name, for `op: "call"` |
| `args` | object of arguments for the op or tool |
| `id` | echoed back verbatim; for your logs, optional |
| `protocol` | optional; rejected if it names a different protocol |
| `options` | per-request overrides (see [Options](#options)) |

A bare `ping` string (not JSON) is also answered, for quick liveness checks.

### Reply

```json
{"protocol": "pyirena-zmq/1", "id": "c7f3", "ok": true,
 "result": { ... }, "elapsed_s": 1.84,
 "server": {"pyirena": "1.1.1", "protocol": "pyirena-zmq/1"}}
```

```json
{"protocol": "pyirena-zmq/1", "id": "c7f3", "ok": false,
 "error": {"code": "NO_SESSION", "message": "...", "suggestion": "..."}}
```

**Check `ok` and nothing else.** A failed fit and a malformed request both
arrive as `ok: false` with a code — there is no second place to look. Every
reply is strict JSON: no `NaN`, no `Infinity`, so `JSON.parse` and Go's
`encoding/json` accept it.

### Operations

| `op` | What it does |
|---|---|
| `ping` | Liveness: uptime, open sessions, the op list |
| `server_info` | **Start here.** Effective options, limits, capabilities and categories — what *this* deployment allows |
| `list_categories` | Tool categories visible in the current profile |
| `list_tools` | `args: {category}` — names and one-line summaries |
| `describe_tool` | `args: {name}` — the full JSON schema; hand it to an LLM as a tool definition |
| `call` | `tool` + `args` — call any listed tool |
| `open_dataset` | `args: {q, intensity, error?, ...}` — create a session **from arrays** |
| `get_session_summary` | `args: {session_id}` |
| `list_sessions` | Every open session, with its `stale` flag |
| `close_session` | `args: {session_id}` — do this when you are done |

### Error codes

| Code | Meaning |
|---|---|
| `BAD_JSON` | Not parseable as a JSON document |
| `BAD_ENVELOPE` | Parsed, but not a valid request object |
| `UNSUPPORTED_PROTOCOL` | `protocol` names something this server does not speak |
| `UNKNOWN_OP` | No such `op` |
| `MESSAGE_TOO_LARGE` | Over `max_message_mb` |
| `BAD_OPTION` | An unknown per-call option, or one trying to widen server policy |
| `BAD_ARGUMENTS` | The op or tool was called without an argument it needs |
| `NOT_AVAILABLE_OVER_ZMQ` | The tool needs a shared filesystem or returns an image |
| `SESSION_STALE` | A previous call on this session hit the time budget |
| `TIMEOUT` | This call exceeded the server's budget (see below) |
| `INTERNAL_ERROR` | A bug; the traceback is in the server log |

Anything else is a pyIrena tool error passed through — `NO_SESSION`,
`NO_MODEL`, `NO_FIT`, `BAD_PARAM`, `SHAPE_MISMATCH`, `TOO_MANY_POINTS`,
`NO_VALID_POINTS`, … Each carries a `suggestion` saying what to do instead.

---

## Timeouts — the one thing to get right

Fits are fast (a 2000-point curve fits in well under a second for most tools;
Size Distribution and Modeling take a few seconds), but they are not bounded,
so the service is:

* The server gives each request **55 seconds** (`request_budget_s`). If the
  work overruns, you get a `TIMEOUT` **reply** — the server never just goes
  quiet.
* Set your client's timeout **above** that. 60 s is the intended pairing.
* On `TIMEOUT` the session is marked **stale**: the abandoned fit may still be
  mutating its model, so every later call on it is refused with
  `SESSION_STALE`. Close it and start again.

### A REQ socket that times out is dead

This is a property of ZMQ, not of pyIrena, and it catches everyone once. A REQ
socket enforces strict send/recv alternation: after a timeout the *next* send
raises `EFSM`, and it will look like a pyIrena failure. **Close the socket and
build a new one** — the "lazy pirate" pattern:

```python
if poller.poll(timeout_ms):
    return sock.recv_string()
sock.close()                 # <-- without this the client is broken from here on
sock = ctx.socket(zmq.REQ)
sock.connect(address)
raise TimeoutError(...)
```

`pyirena/zmq/client.py` does this for you and is ~150 lines you can copy.

---

## One call: `analyze`

Most of the time the fit is already known — a scientist set it up in the GUI,
exported the parameters, and every new measurement should be fitted the same
way. That is one request, not six:

```json
{"op": "call", "tool": "analyze",
 "args": {"data": {"q": [...], "intensity": [...], "error": [...]},
          "config": { ...the exported pyirena_config.json... }}}
```

```python
results = pyirena.call("analyze", data={"q": q, "intensity": I, "error": e},
                       config=json.load(open("modeling.json")))
```

It opens its own session, applies the configuration, fits, returns the same
payload `export_results` does (plus an `analyze` block), and closes the
session in a `finally` — so a client that dies mid-call cannot leak one.

**The tool is inferred** from the config, because a GUI export holds exactly
one tool section. Pass `tool` only to override it.

**The whole configuration is applied, not just the model.** A config carries
the fitted Q range, the slit settings, which parameters are held fixed, and
per-tool preparatory steps — the Size Distribution's power-law and flat
background pre-fits, the WAXS peak Q0 presearch. Those change the answer, so
skipping them would return a plausible wrong number rather than an error.
Anything adjusted is reported:

```json
"analyze": {"tool": "sizes", "config_applied": true,
            "fit_q_min": 0.00493, "fit_q_max": 0.0736,
            "fixed_parameters": [], "n_points_input": 516,
            "cleaning": {"n_input": 516, "n_kept": 516, "...": "..."},
            "notes": ["Power-law pre-fit: B = 0.0004554, P = 3.183",
                      "Background pre-fit: 0.4996 cm⁻¹"]}
```

Read `notes` — it is where a Q range that missed the curve, or a pre-fit that
failed, is reported rather than swallowed.

`analyze` and `pyirena.batch` share their config-to-model translation
(`pyirena/core/tool_config.py`), so a fit replayed over the wire and the same
fit run from the command line produce identical numbers. The test suite
asserts that on real exported configs for all six tools.

### When not to use it

Use the session tools when the fit is *not* already known — when an agent
needs to look at the data, try a model, free a parameter and refit. `analyze`
is one shot: it cannot iterate.

## A session, end to end

```python
from pyirena.zmq.client import PyIrenaClient

with PyIrenaClient("tcp://workstation:9865", timeout_ms=60_000) as pyirena:
    sid = pyirena.open_dataset(q, intensity, error, label="scan_007")["session_id"]
    try:
        pyirena.call("select_model", session_id=sid, model_name="unified_fit")
        pyirena.call("add_unified_level", session_id=sid)
        pyirena.call("set_fit_q_range", session_id=sid, q_min=1e-3, q_max=0.1)
        pyirena.call("run_fit", session_id=sid)

        quality = pyirena.call("get_fit_quality", session_id=sid)
        results = pyirena.call("export_results", session_id=sid)
    finally:
        pyirena.close_session(sid)
```

Without the helper client, each of those is one `send_string` of

```json
{"op": "call", "tool": "run_fit", "args": {"session_id": "a1b2c3d4"}}
```

### Getting data in

`open_dataset` takes arrays, never a path:

| Argument | |
|---|---|
| `q` | **required** — Q in Å⁻¹. Need not be sorted; the session sorts it |
| `intensity` | **required** — I, cm⁻¹ where the data is on an absolute scale |
| `error` | uncertainty on I. **Omit it for an unweighted fit** — uncertainties are repaired when you send them, never invented when you do not |
| `dq` | Q resolution; stored for provenance, not yet used in fitting |
| `label` | name for the curve, used in reports |
| `is_slit_smeared`, `slit_length` | slit-smeared data; `slit_length` in Å⁻¹ |

The data is cleaned exactly as a text file would be: points with Q ≤ 0 or
I ≤ 0 (beamstop zeros, direct-beam rows) and any non-finite value are removed,
and non-positive uncertainties are replaced by `error_fraction × I`. **The
counts come back** in `summary.cleaning`, so nothing is dropped silently:

```json
{"session_id": "a1b2c3d4",
 "summary": {"file": null, "n_points": 1998, "q_min": 0.00013, "q_max": 0.61,
             "has_errors": true, "sorted_by_q": false,
             "cleaning": {"n_input": 2000, "n_kept": 1998,
                          "n_removed_q": 1, "n_removed_i": 1,
                          "n_repaired_error": 0}}}
```

### Getting results out

`export_results` is the one tool for all six fitting tools:

```json
{"op": "call", "tool": "export_results",
 "args": {"session_id": "a1b2c3d4", "include_arrays": false}}
```

```json
{"ok": true, "tool": "unified_fit", "pyirena_version": "1.1.1",
 "exported_at": "2026-09-25T14:02:11+00:00",
 "data": {"label": "scan_007", "file": null, "n_points": 1998, "...": "..."},
 "fit_q_range": {"q_min": 0.001, "q_max": 0.1},
 "quality": {"chi_squared": 1841.2, "reduced_chi_squared": 1.02,
             "dof": null, "n_points_fitted": 1204, "n_parameters": null,
             "success": true, "message": "...", "metrics": {"...": "..."}},
 "results": {"...": "tool-specific: parameters, peaks, populations, derived"},
 "config": {"...": "the model's to_dict() — replay this fit on the next scan"}}
```

`config` is the model's own `to_dict()` — what you need to rebuild the model
and fit again. For Modeling, the Carbon model and Simple Fits it is already
the shape `pyirena.batch` and the GUI's *Export Parameters* use, so a fit can
be replayed as-is. **Unified Fit is the exception**: the GUI and batch write a
second, older vocabulary there (`RgCutoff` for `RgCO`, and each parameter as
`{"value": …, "fit": …}` rather than a bare number), so the two are not
interchangeable yet. Note also that `config` is *model* state only — the
fitted Q range and slit settings live on the session and are reported
separately under `fit_q_range` and `data`.

Set `include_arrays: true` for the curves — `q`, `intensity`, `error`,
`intensity_model`, `residuals`, plus per-population curves (Modeling) or the
three components (Carbon model). Arrays are capped at `max_points` and flagged
`decimated` if they were thinned. A report is 1–10 KB without arrays and
200–350 KB with them.

---

## Options

Three layers; each overrides the one above.

1. **Config file** — `--config FILE`, or `$PYIRENA_ZMQ_CONFIG`. JSON, or TOML
   on Python 3.11+. Keys are the option names below.
2. **CLI flags** — `--max-sessions 4`, `--allow-images`, …
3. **Per-call `options`** in the request envelope, for the four marked below.

```json
{
  "bind": "tcp://0.0.0.0:9865",
  "request_budget_s": 55,
  "max_sessions": 8,
  "allow_images": false
}
```

| Option | Default | Per-call | |
|---|---|:-:|---|
| `bind` / `--port` | `tcp://0.0.0.0:9865` | | ZMQ endpoint |
| `max_message_mb` | 16 | | Larger requests are refused |
| `request_budget_s` | 55 | ✓ | Reply `TIMEOUT` rather than exceed this |
| `max_input_points` | 100000 | | Largest curve accepted |
| `session_ttl_min` | 60 | | Evict sessions idle this long |
| `max_sessions` | 32 | | Cap; oldest idle is evicted first |
| `include_arrays` | false | ✓ | Default for `export_results` |
| `max_points` | 2000 | ✓ | Cap on each returned array; `null` = no cap |
| `allow_images` | false | ✓ | Expose image tools as base64 PNG (no paths) |
| `allow_files` | false | | Expose file tools; **requires** `data_root` |
| `data_root` | — | | Confines all file access to this directory |
| `log_file`, `log_level`, `log_stderr` | rotating, INFO, off | | |

A per-call option may only **narrow** what the operator allowed: a request can
turn images off but never on, and can shorten its own deadline but never
extend it. Anything else is `BAD_OPTION`. An unknown option name is an error
too — never silently ignored.

Call `server_info` to see the effective values, so the same client code works
against a locked-down and a permissive deployment.

### Profiles — why some tools are missing

The visible tool set follows from `allow_files` and `allow_images`:

| `allow_files` | `allow_images` | Profile | Hidden |
|:-:|:-:|---|---|
| off | off | `json_only` *(default)* | `save_*`, `open_dataset` (the path form), all data operations, all `*_image` tools |
| off | on | `json_images` | file tools only |
| on | any | `all` | nothing |

By default the whole `data` category is gone, because every data operation
reads and writes files on the server. A hidden tool called anyway returns
`NOT_AVAILABLE_OVER_ZMQ` rather than pretending to work and handing back a
path you cannot open.

---

## For an agent driving the service

1. `server_info` — what this deployment allows.
2. **If a saved configuration exists for this kind of measurement, call
   `analyze` and stop here.** It is one request, it applies the whole
   configuration, and it cannot leak a session.
3. Otherwise: `list_categories` → `list_tools(category)` →
   `describe_tool(name)`. The schemas are Anthropic tool-use shaped; pass
   them to the model as tools and wrap each call it makes in
   `{"op": "call", "tool": ..., "args": ...}`.
4. `open_dataset` with the arrays → `session_id`.
5. `select_model` → set parameters / fix what should not float → `set_fit_q_range`.
6. `run_fit` → `get_fit_quality`.
7. `export_results` → keep the JSON.
8. `close_session`. **Always** — sessions are evicted eventually, but a
   forgotten one holds memory until then.

One session per curve. Two callers can use the service at once, but requests
are served one at a time, so a long fit delays everyone — including `ping`.

---

## Limits in this version

* **Synchronous only.** No job queue; a call either finishes in the budget or
  returns `TIMEOUT`. In practice a fit takes well under a second to a few
  seconds; the one setting that reliably outruns the budget is Monte-Carlo
  uncertainties, which is why they are never run unless asked for.
* **One request at a time.** Fits are CPU-bound and the session registry is
  not thread-safe, so the service serialises deliberately.
* **JSON arrays only.** No binary framing — it would break the one-frame
  string contract. A 2000-point curve is ~100 KB of JSON, comfortably inside
  the 16 MB cap.
* **No authentication or encryption.** Firewall the port.
* **Sessions do not survive a restart.** They live in the process.

See [`planning/zmq-service/`](../planning/zmq-service/) for what is deferred
and why (async jobs, Zarr input, binary arrays, CURVE authentication).
