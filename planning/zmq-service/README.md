# ZMQ Service — decisions and what is left

Internal planning artifact, not user-facing documentation. Written
24-09-2026; Phases 1–3 shipped on `feature/zmq-service` (25–26-09-2026) and
were released with 1.2.0. The protocol, profiles, options and error codes are
documented in **`docs/zmq_service.md`**; the cross-machine test run
(29-09-2026, `usaxscontrol:9865`) in `docs/zmq_service_testing.md`.

Shrunk 11-11-2026 to the decisions and the open work; the 750-line build plan
for the shipped phases is in git history.

---

## 1. Why this exists

An experiment project runs an AI **orchestrator** on a modest machine that
drives instrument and analysis **services** over ZMQ. pyIrena is one of those
services, running on an analysis workstation. The machines are a mix of
operating systems, Windows included, and **share no filesystem**: JSON in,
JSON out, data inside the request, results inside the reply.

## 2. Decisions taken (24–26-09-2026, with the project owner)

| Decision | Rationale |
|---|---|
| The service lives **inside pyIrena** as `pyirena/zmq/`, an optional `[zmq]` extra with a `pyirena-zmq` entry point — the same shape as `[mcp]` | One package to maintain. Explicitly chosen over a companion package. |
| **JSON in / JSON out**; no file paths cross the wire | Owner's requirement; mixed OS, no shared filesystem. |
| **One UTF-8 text frame per message**, JSON document as the payload | Matches the owner's client (`send_string` / `recv_string`). No multipart — see §4 for what that costs binary arrays. |
| **Plain REQ/REP**, synchronous | Confirmed by the owner. |
| **Default port 9865**, TCP | Chosen by the owner. |
| **No async job API** | Fits finish well inside the client's 60 s timeout; anything slower should time out rather than be queued. See §3.2. |
| **The server always replies inside a wall-clock budget** (default 55 s) — an overrunning fit gets a `TIMEOUT` reply, not silence | A REQ socket that times out is unusable and must be recreated. Replying keeps the client's socket alive and says why it failed. |
| **A tool-level error is `ok: false`** (settled 25-09-2026) | One place for the caller to look; `pyirena.zmq.protocol._reply_for_result` is the one line that decides it. |
| **Server platform: Red Hat Linux**, Python 3.12/3.13 | Owner: any Python on the server. The *client* may still be Windows. |
| **Typical payload ≈ 2000 points per curve** | ~100 KB of JSON; binary framing is not worth the break from the owner's client. |
| The service is a **thin transport over the existing dispatcher** (`pyirena/api/dispatch.py`), like MCP | Protocol-agnostic, returns JSON dicts, serves tool schemas an agent can use directly. **Profiles** (`PROFILE_JSON_ONLY`) hide tools that take a `file_path` or return an `image_path`. |
| **Server options are a discoverable surface** | Owner: support reasonable variation of options over time. Cheap ones shipped; the rest have a defined place to land (§4). |
| **`analyze` (one-shot recipe) is the primary documented path**; session tools are the fallback for hard fits | One round trip, replayable from an exported config. |
| **Uncertainties travel in `export_results`** with an `uncertainty` block saying which estimate they are; Monte-Carlo is opt-in (`n_mc_runs`, default 0) because it can outrun the budget | Decided 25-09, extended to all tools by #38 (11-11-2026). |

This partly overrides the ai-agent decision "separate package" — that was
about the standalone AI app (LLM SDKs). The ZMQ service adds only `pyzmq`, so
it fits the extras model.

## 3. What is left before / after deployment

### 3.1 Deployment mechanics — issue #32

A systemd unit (drafted in `docs/zmq_service.md`), a service account, the
conda/venv the service runs from, and the firewall rule. Open questions:

1. Does anything else in the project (other services) already solve
   logging/service management on the RHEL box, so pyIrena should match it?
2. Is 9865 reachable from the orchestrator through the beamline firewall, and
   who owns that rule?

Carried forward from the 29-09 cross-machine run:

- The server is **~3× slower per fit than a laptop** (13.4 s vs 4.8 s for
  `ModelingSF_SaD`), so the 55 s budget has less headroom than local timings
  suggest.
- **An abandoned request still runs to completion**, and because the service
  is deliberately serial the next caller waits it out (measured 12.3 s). A
  client that times out and retries at once queues behind itself — stated in
  `docs/zmq_service.md`; the orchestrator must back off between retries.

### 3.2 Watch Modeling; asynchronous jobs only if triggered — issue #33

Modeling is the one tool whose fit can plausibly reach the 55 s budget with
several populations. If real use hits `TIMEOUT` regularly, that is the signal
for async jobs — not a reason to raise the budget past the client's timeout.
Sketch kept for that day:

- Switch the server to ROUTER; a job table plus a worker (thread first;
  process pool if the GIL bites — sessions then live in the worker).
- New ops: `submit` (returns `job_id`), `status`, `result`, `cancel`,
  `list_jobs`. `ping` answers during fits.
- Cancellation needs a cooperative check in the fit loops — the `cancel_cb`
  hooks added to the Monte-Carlo loops for #38 are a start; a single scipy
  optimiser call still cannot be interrupted.
- Optional PUB socket for progress / completion events.
- Keep synchronous `call` working for short operations.

The wall-clock deadline is deliberately the same seam a job API would use, so
this is an extension, not a rewrite.

## 4. Richer data paths (deferred; pick by need)

| Option | When | Sketch |
|---|---|---|
| **File paths on a shared filesystem** | If a deployment shares storage | Reserved as `allow_files` + `data_root`; flipping it switches the profile to `all`. Mostly built. |
| **Zarr input** | When the project's Zarr layout stabilises | `pyirena/io/zarr_sas.py` behind a `[zarr]` extra, returning the same dict as `readGenericNXcanSAS`; `open_dataset` dispatches on file type. Layout mapping configurable. |
| **Binary arrays** | If curves get far larger than 2000 points, or many per call | Multipart JSON header + raw float64 frames. **Breaks the owner's `send_string` client** — needs their buy-in. A base64 block inside the JSON is the compatible middle road. |
| **Images over the wire** | If the agent is multimodal | Reserved as `allow_images`: `image_base64`, no `image_path`. |
| **Results written server-side** | If the orchestrator wants NXcanSAS sidecars | `save_*` with an `output_path` under `data_root`; return path + JSON. |
| **Authentication / encryption** | If the port leaves a trusted subnet | ZMQ CURVE with keys in a config file. Until then, firewall. |
| **Multiple workers / clients** | Several orchestrators or parallel fits | ROUTER–DEALER broker, session affinity per worker. |
