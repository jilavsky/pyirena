# 01 — API & MCP Extensions (Foundation)

**Status:** Draft
**Last updated:** 2026-06-09
**Depends on:** nothing
**Consumed by:** [02 Standalone AI app](02-standalone-ai-app.md), [03 In-GUI AI advisor](03-ai-advisor-in-gui.md)

---

## Goal

Extend `pyirena/api/` from a read-only result-summarization layer into a
complete **control surface** that an AI agent (or any external automation)
can drive end-to-end: load data → configure model → set/fix/free parameters
→ run fit → read residuals and quality → iterate.

This is the **foundation**. Both subprojects 2 and 3 consume these tools.
MCP exposure of the new tools is a thin wrapper added at the end.

---

## Current state

Today `pyirena/api/` contains:
- `discovery.py` — find NXcanSAS files in a folder
- `data.py` — read raw data arrays
- `results.py` — read existing fit results from HDF5
- `aggregate.py` — combine parameters across files
- `plotting.py` — generate plot images
- `schemas.py` — JSON schemas for MCP tools

All of this is **read-only**: it tells the AI what *was* fit, not how to fit.

The MCP server (`pyirena/mcp/`) exposes 18 tools, all observation-side.

---

## Tool categories to add

### Category A — Session / state management
Functions an agent calls to bring data into a fitting session, before
touching any model.

- `open_dataset(path)` → session_id; loads NXcanSAS file
- `list_open_sessions()` → list of active sessions
- `close_session(session_id)`
- `get_session_summary(session_id)` → metadata, q range, intensity range

### Category B — Model selection and inspection
- `list_available_models()` → `['unified', 'sizes', 'fractals', 'reflectivity', ...]`
- `select_model(session_id, model_name, **model_options)` — e.g. for Unified, set number of levels
- `get_model_parameters(session_id)` → full parameter table: name, value, fixed/free, lo/hi bounds
- `get_model_description(session_id)` → text description of what this model does and what each parameter means (so the AI doesn't need to know SAXS conventions a priori)

### Category C — Parameter control
The granularity the AI needs to do multi-stage fitting properly.

- `set_parameter_value(session_id, param_name, value)`
- `set_parameter_bounds(session_id, param_name, lo, hi)`
- `fix_parameter(session_id, param_name)` — hold at current value
- `free_parameter(session_id, param_name)` — release for fitting
- `fix_all_except(session_id, free_list)` — convenience for staged fitting
- `reset_parameters_to_defaults(session_id)`

### Category D — Fit execution
- `run_fit(session_id, max_iter=None, tolerance=None)` → returns fit result summary (success, χ², iterations, final values)
- `cancel_fit(session_id)` — for long-running fits
- `get_fit_status(session_id)` — for async fits if we go that way

### Category E — Quality assessment
The most important readback for AI decision-making.

- `get_chi_squared(session_id)`
- `get_residuals(session_id)` → numeric residuals array + summary stats
- `get_residuals_image(session_id)` → PNG bytes of fit + residuals plot (so AI vision can evaluate)
- `get_fit_image(session_id)` → PNG bytes of data + model overlay
- `get_parameter_uncertainties(session_id)`

### Category F — Persistence
- `save_fit(session_id, output_path)` — write NXcanSAS HDF5 with fit result
- `export_fit_report(session_id, format='json'|'markdown')`

---

## Design principles

1. **Session-based state, not global.** The AI may have multiple datasets
   open. Every tool takes a `session_id`. Avoid hidden module-level state.

2. **Synchronous fits to start.** `run_fit` blocks until done. Async only if
   fits get long enough to matter (Subproject 2 phase 3 concern).

3. **JSON-only types in tool signatures.** Strings, numbers, booleans,
   lists, dicts, PNG bytes (base64). No numpy arrays, no DFREFs, no objects
   the LLM can't reason about.

4. **Self-describing tools.** `get_model_description` returns enough text
   that an LLM with no SAXS knowledge could still fit. The system prompt
   doesn't need to encode every model's parameter conventions.

5. **Idempotent where possible.** `set_parameter_value` twice with the same
   value is a no-op. `run_fit` always produces a result object even if it
   doesn't converge.

6. **Errors are data.** Tools return `{"error": "...", "suggestion": "..."}`
   rather than throwing exceptions across the MCP boundary. The agent reads
   the error and decides what to do.

---

## Wiring pattern

Following the existing convention you described:

1. Add the function to `pyirena/api/<area>.py`
2. Add a wrapper in `services/pyirena_tools.py` (in the consuming AI app)
3. Register it in `_TOOL_FUNCS` dict and add its JSON schema in
   `services/llm_service.py`

For MCP exposure (later):

4. Add a corresponding MCP tool in `pyirena/mcp/` that calls the same
   `pyirena.api` function

The contract is the **pyirena.api function**, not the MCP tool or the
JSON schema — those are two different wrappers around the same thing.

---

## Phasing

### Phase 1 — Minimum for advisor (small)
The in-GUI advisor (Subproject 3) needs only:
- Read current model parameters
- Read fit quality (χ², residuals image)
- Capture a plot image

This may already be partially in place via the existing read-only tools.
Audit and fill gaps.

### Phase 2 — Read-modify-write parameter control (medium → large)
Add categories B (model selection), C (parameter control), and basic D (synchronous `run_fit`).

This is the bulk of the work. Per-model functions are needed for each
supported tool. Each tool has different parameter conventions and is
expected to be its own implementation + debugging project.

**Supported tools (final list — see [00-overall-plan.md](00-overall-plan.md)):**
analysis — Unified Fit, Modeling, Size distribution, Simple fits, WAXS;
utilities — Scattering contrast calculator, Data merge, Data manipulation.

**Strategy:** ship **Unified Fit first as the test case** for all
infrastructure (session model, tool granularity, error reporting, image
capture). Once that pattern is proven, replicate to other tools in an order
to be decided based on observed per-tool cost. **Modeling is expected to
be the largest single effort** because of its parameter richness.

### Phase 3 — Quality and persistence (small)
Polish Category E (uncertainties, residual stats) and add Category F
(save fit, export report).

### Phase 4 — MCP wrappers (small)
Once API tools are stable, expose them via MCP. The existing
`pyirena/mcp/` infrastructure makes this a mechanical mapping.

---

## Resolved decisions

| Decision | Rationale |
|----------|-----------|
| **Unified Fit = one model with `nlevels` parameter** (not one model per level configuration) | Matches the way a human approaches the problem — start with one level, add more as needed. Same instruction patterns will work for AI. Makes "add a level mid-fit" a single tool call rather than a model switch. |
| **Sessions in-memory only for v1** | Per-fit working state (current values, fixed/free, bounds) lives in the AI app's process. No need to persist sessions across app restarts initially. |
| **Fit cancellation deferred** | Existing pyirena fits are already iteration-limited and complete in seconds at most. Add cooperative cancellation later if real users hit minutes-long fits. |

## Open questions

| Question | Notes |
|----------|-------|
| **Per-tool wrapping cost** | Unknown until Unified Fit is wrapped end-to-end as the test case. Modeling expected to be the largest. Each tool will likely be its own implementation + debugging project. Estimate after Unified Fit lands. |
| **Interaction with pyirena's existing per-tool state persistence** | pyirena already saves per-tool / per-user GUI parameter state in a JSON file, so re-opening a tool restores prior parameters. This is GUI-side machinery and may be bypassed when running headless via the API. **Action before implementation:** audit how that persistence interacts with the API layer — does an API session see those defaults, ignore them, or could it inadvertently mutate the user's GUI defaults? Decide whether the AI session should (a) start from per-tool defaults, (b) start from a clean slate, (c) be configurable. |
| **Numerical reproducibility** | If the agent re-runs the same fit twice, results must match. Audit any random-seed usage in MCSaS / Reg / etc. |

---

## Out of scope

- New fitting algorithms (we wrap existing ones, don't invent)
- GUI changes (this layer is pure API; GUI work lives in Subprojects 2/3)
- Schema migration of NXcanSAS HDF5 files (uses existing pyirena schema)
- Multi-user / network access (everything is local, single-process)

---

## Success criteria

A Python script using only `pyirena.api` can:
1. Open a NXcanSAS file
2. Select Unified Fit with 2 levels
3. Fix all parameters except scale and background
4. Run fit
5. Read χ² and residuals image
6. Decide (programmatically) to free Rg and run again
7. Save the final fit back to HDF5

If that script works end-to-end, Subproject 1 is done.
