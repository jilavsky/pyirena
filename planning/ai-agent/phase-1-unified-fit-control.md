# Phase 1 — Detailed Implementation Plan: Unified Fit Control Surface

**Status:** Draft — for review before implementation
**Last updated:** 2026-06-09
**Branch:** `feature/ai-api-mcp`
**Parent docs:** [00-overall-plan.md](00-overall-plan.md), [01-api-and-mcp-extensions.md](01-api-and-mcp-extensions.md)

---

## What "Phase 1 done" looks like

End state — a single Python script using only `pyirena` can:

1. Open a NXcanSAS file
2. Select Unified Fit with 1 level
3. Read the current parameter table (names, values, fixed/free, bounds)
4. Set starting values, fix all but scale and background, run a fit
5. Read χ² and capture a residuals image
6. Programmatically decide to free Rg, run again
7. Add a 2nd level mid-session, repeat
8. Save the fit to NXcanSAS and export a JSON report

Plus: every public function has a JSON schema registered, so an LLM agent
(or the in-GUI advisor in Subproject 3) can call them through tool-use.

If that script works end-to-end on a real USAXS dataset and the schemas
validate against Anthropic's tool-use format, **Phase 1 is done**.

---

## Codebase landscape (current state)

Confirmed by inspection:

| Component | Location | What it gives us |
|-----------|----------|------------------|
| `UnifiedFitModel` class | [pyirena/core/unified.py:98](../../pyirena/core/unified.py) (636 lines total) | The actual fitting engine. Already has per-parameter `fit_*` flags and `*_limits` tuples. |
| `UnifiedLevel` dataclass | [pyirena/core/unified.py:27](../../pyirena/core/unified.py) | Per-level parameter container with fit-control fields already in place |
| Headless config-driven fit | [pyirena/batch.py:383](../../pyirena/batch.py) (`fit_unified`) | Existing single-shot API: takes a config JSON + data file, returns result. **Not stateful — runs once, returns.** |
| Read-only API facade | [pyirena/api/__init__.py](../../pyirena/api/__init__.py) | `read_unified_fit` etc. Returns JSON-safe dicts. No Qt. |
| MCP server | [pyirena/mcp/server.py](../../pyirena/mcp/server.py) | Single-file MCP server wrapping the read API |
| Per-tool state persistence | [pyirena/state/state_manager.py](../../pyirena/state/state_manager.py) | `~/.pyirena/state.json`. **GUI-side only.** batch.py uses its own JSON config format. |
| GUI Unified Fit panel | [pyirena/gui/unified_fit.py](../../pyirena/gui/unified_fit.py) (4095 lines) | Where most parameter-control UX lives; we will mirror its semantics, not depend on it |

### Key implication
The fitting engine (`UnifiedFitModel`) already has the per-parameter control
hooks we need (`fit_Rg`, `Rg_limits`, etc.). **Phase 1 is mostly a stateful
session wrapper + JSON schemas + image capture, not new fitting code.**

The GUI state persistence and batch.py config format are *separate* from
each other — confirming the clean separation. We can add a third path
(session-based API) without disturbing either.

---

## Pre-implementation investigations (do FIRST)

Before writing any API code, do these three audits. Each is a small
focused task whose result changes the design. Expect 1-2 hours each.

### I-1 — State persistence interaction audit
**Open question from [01](01-api-and-mcp-extensions.md).**

Walk through: when a script calls `open_dataset()` and then `select_model('unified_fit')`,
does anything read from `~/.pyirena/state.json`? Does running fits via
`UnifiedFitModel` directly (bypassing the GUI) mutate that file?

**Output:** a short note in this doc or in `01` confirming one of:
- (a) Defaults come from `state.json` → AI session inherits user GUI prefs
- (b) Defaults come from hardcoded model defaults → AI session is clean
- (c) Behavior is currently inconsistent → we must pick one

**Decision needed before I-2.** Recommendation likely (b) — clean slate per
AI session — but verify before locking in.

### I-2 — `UnifiedFitModel` introspection
Build a complete picture of what the model exposes:
- All level parameters and their meanings
- All model-wide parameters (background, scale, q range)
- Which `fit_*` flags exist
- Which `*_limits` tuples exist
- How results come back from `.fit()` (success status, χ², residuals, uncertainties)
- Whether re-running `.fit()` reuses prior state or resets

**Output:** a one-page reference table of every parameter the AI will see,
with current default values, units, and physical meaning. Becomes the
basis for `get_model_description()`'s output text — and feeds the system
prompt for the advisor.

### I-3 — Plot/image generation path
The advisor (Subproject 3) and the agent (Subproject 2) both need
PNG bytes of the current fit. Today, GUI panels render via pyqtgraph;
headless rendering is via matplotlib in [pyirena/api/plotting.py](../../pyirena/api/plotting.py)
and in [pyirena/batch.py](../../pyirena/batch.py).

**Question:** can we use the existing matplotlib path to render a
"data + model + residuals" image from a `UnifiedFitModel` instance, or do
we need a new renderer?

**Output:** decision — reuse existing or write thin new renderer in
`pyirena/api/control/plotting.py`.

---

## Architecture

### New module: `pyirena/api/control/`

Sibling to the existing `pyirena/api/` (which stays read-only). Keeping
them separate avoids accidentally importing Qt or stateful machinery into
the read-only facade.

```
pyirena/api/control/
├── __init__.py           # Public API surface; re-exports tool functions
├── README.md             # AI-agent-facing reference (mirrors api/README.md style)
├── session.py            # Session class + in-memory registry
├── unified_fit.py        # Unified Fit-specific control tools
├── plotting.py           # Capture fit/residuals as PNG bytes (or import from api/plotting)
├── schemas.py            # JSON schemas for each tool (for LLM tool-use)
└── errors.py             # Structured error shape: {"error": "...", "suggestion": "..."}
```

### Session model

A `Session` is an in-memory object holding:
- The loaded data (Q, I, error arrays)
- The current `UnifiedFitModel` instance (or future Modeling/Sizes/...)
- Last fit result (χ², residuals, uncertainties)
- Metadata: file path, creation time

Sessions are stored in a module-level dict keyed by `session_id` (UUID string).

Lifetime: process lifetime. No persistence across restarts (per resolved
decision in [01](01-api-and-mcp-extensions.md)).

```python
# Sketch
_SESSIONS: dict[str, Session] = {}

class Session:
    session_id: str
    file_path: str
    data: dict  # {"Q": [...], "Intensity": [...], "Error": [...]}
    model_name: str | None
    model: object | None  # UnifiedFitModel | ... (per-tool union later)
    last_fit: FitResult | None
```

### Tool function pattern

Every public tool is a module-level function in `pyirena/api/control/` that:
- Takes plain JSON types (strings, numbers, lists, dicts)
- Returns a JSON-serializable dict
- Never raises across the API boundary; errors become `{"error": "...", "suggestion": "..."}`
- Looks up its session by `session_id`; returns an error dict if not found

This mirrors the read-only API's conventions exactly.

---

## Tool catalog for Phase 1 (Unified Fit)

Grouped by the categories from [01](01-api-and-mcp-extensions.md).
**Signatures are the contract — review carefully.** Return shapes are
illustrative; finalized during implementation.

### Category A — Session lifecycle
```python
open_dataset(file_path: str) -> {
    "session_id": str,
    "summary": {"file": str, "n_points": int, "q_range": [float, float],
                "intensity_range": [float, float]}
}

list_open_sessions() -> {"sessions": [{"session_id": str, "file": str, "model": str|None}, ...]}

close_session(session_id: str) -> {"ok": bool}

get_session_summary(session_id: str) -> {"file": str, "model": str|None,
                                          "n_points": int, "q_range": [...],
                                          "has_fit": bool, "chi_squared": float|None}
```

### Category B — Model selection and inspection
Phase 1 supports only `"unified_fit"`. Other tools added in later phases.

```python
list_available_models() -> {"models": ["unified_fit"]}  # grows over time

select_model(session_id: str, model_name: str = "unified_fit",
             nlevels: int = 1) -> {"ok": bool, "model": "unified_fit",
                                    "nlevels": int, "parameters": [...]}

get_model_parameters(session_id: str) -> {
    "parameters": [
        {"name": "Rg_1", "value": 10.0, "fixed": False,
         "lo": 1.0, "hi": 1000.0, "units": "Å",
         "description": "Radius of gyration of level 1"},
        ...
    ]
}

get_model_description(session_id: str) -> {
    "model": "unified_fit",
    "summary": "The Unified fit model describes hierarchical structures using ...",
    "parameters": {"Rg_1": "...", "G_1": "...", ...},
    "tips": ["Start with one level", "Add levels only if residuals show ..."]
}
```

Parameter naming convention: `<param>_<level_index>` for level parameters
(e.g. `Rg_1`, `G_2`). Model-wide parameters (background, scale) use plain
names. **Confirm with I-2 audit before locking.**

### Category C — Parameter control
```python
set_parameter_value(session_id, param_name: str, value: float) -> {"ok": bool, "value": float}

set_parameter_bounds(session_id, param_name: str,
                     lo: float, hi: float) -> {"ok": bool, "lo": float, "hi": float}

fix_parameter(session_id, param_name: str) -> {"ok": bool}

free_parameter(session_id, param_name: str) -> {"ok": bool}

fix_all_except(session_id, free_list: list[str]) -> {"ok": bool, "fixed_count": int, "free_count": int}

reset_parameters_to_defaults(session_id) -> {"ok": bool}
```

### Category C-prime — Unified-specific level operations
Per resolved decision: Unified Fit = one model with `nlevels` parameter.

```python
add_unified_level(session_id, position: int = -1) -> {"ok": bool, "nlevels": int, "new_parameters": [...]}
# position=-1 means append; otherwise insert at that level index (1-based)

remove_unified_level(session_id, level: int) -> {"ok": bool, "nlevels": int}
```

### Category D — Fit execution
```python
run_fit(session_id, max_iter: int|None = None,
        tolerance: float|None = None) -> {
    "success": bool,
    "chi_squared": float,
    "reduced_chi_squared": float,
    "iterations": int,
    "message": str,        # e.g. "Converged" / "Max iterations reached"
    "parameters_updated": [{"name": str, "value": float, "uncertainty": float|None}, ...]
}
```

Synchronous; no cancellation in Phase 1 (resolved).

### Category E — Quality assessment
```python
get_chi_squared(session_id) -> {"chi_squared": float, "reduced_chi_squared": float}

get_residuals(session_id) -> {
    "residuals": [...],            # decimated to max_points
    "summary": {"rms": float, "max_abs": float, "mean": float}
}

get_fit_image(session_id, width: int = 1024, height: int = 768,
              format: str = "png") -> {
    "image_base64": str,
    "format": "png",
    "width": int, "height": int
}
# Image: log-log data + model overlay + residuals subplot. Matches typical
# pyirena fit display so the AI sees what the user would see.

get_residuals_image(session_id, width: int = 1024, height: int = 768) -> {...}

get_parameter_uncertainties(session_id) -> {"uncertainties": {param_name: float, ...}}
```

### Category F — Persistence
```python
save_fit(session_id, output_path: str|None = None) -> {"ok": bool, "file_path": str}
# If output_path is None, save back to the original file. Writes NXcanSAS
# fit result block per existing pyirena conventions.

export_fit_report(session_id, format: str = "json") -> {
    "format": "json" | "markdown",
    "content": str         # JSON string or markdown text
}
```

---

## Implementation milestones

Each milestone is independently reviewable / committable. Each ends with
a working state demonstrable by a small script.

### M1 — Foundation: sessions + read-side tools (small-medium)
- `pyirena/api/control/` skeleton + `errors.py` shape
- `Session` class + registry
- `open_dataset`, `list_open_sessions`, `close_session`, `get_session_summary`
- `list_available_models`, `select_model`, `get_model_parameters`,
  `get_model_description`
- `get_fit_image`, `get_residuals_image` (works before any fit — shows raw data + initial model)

**Demo script:** open dataset, select Unified Fit, read parameters, save plot to disk.
**Unblocks:** Subproject 3 (advisor) MVP can start in parallel with M2.

### M2 — Parameter control + first fit (medium)
- All Category C tools (`set_*`, `fix_*`, `free_*`, `reset_*`)
- `fix_all_except`
- `run_fit` (synchronous, returns full result)
- `get_chi_squared`, `get_residuals`, `get_parameter_uncertainties`

**Demo script:** open, select model, set starting values, fix all but
scale + background, run fit, check χ², free Rg, run again. The
success-criteria script from [01](01-api-and-mcp-extensions.md) works.

### M3 — Level operations (small)
- `add_unified_level`, `remove_unified_level`
- Parameter renumbering: when you remove level 1, what happens to level
  2's name? **Decision needed in implementation:** renumber and notify, or
  stable IDs?

**Demo script:** fit with 1 level, residuals show extra structure → add
level 2, free its params, fit again.

### M4 — Persistence (small)
- `save_fit` (NXcanSAS write — reuse existing path from batch.py?)
- `export_fit_report` (JSON and markdown)

### M5 — JSON schemas + MCP wrapper (small)
- `pyirena/api/control/schemas.py` with one schema per tool, in Anthropic
  tool-use format
- New MCP tools in `pyirena/mcp/server.py` calling the control API
- Update [docs/ai_tools_reference.md](../../docs/ai_tools_reference.md)
  with the new tools

**End-of-M5 demo:** the success-criteria script runs, AND a Claude agent
session via MCP can fit a dataset by itself given the system prompt.

---

## Testing strategy

### Unit tests
- One test per tool function in `pyirena/tests/test_api_control_*.py`
- Use a small fixture HDF5 file (existing testData/ likely has candidates)
- Test happy path + at least one error path (bad session_id, bad param name)

### Integration test
- One end-to-end test that runs the success-criteria script of Phase 1
- Lives in `pyirena/tests/test_phase1_unified_fit.py`

### Schema validation
- Test that every schema in `schemas.py` is valid Anthropic tool-use JSON
  (probably via `anthropic.types.tool_param.ToolParam` or similar)

### Manual validation
- Run an actual Claude agent session against the MCP server with a real
  USAXS dataset; verify the agent can produce a reasonable fit
- This is the real test — automated tests can't catch "agent can't figure
  out what to do with these tools"

---

## Cross-cutting decisions for Phase 1

| Decision | Choice |
|----------|--------|
| Sandboxing | Reuse `PYIRENA_DATA_ROOT` env-var pattern from read-only API |
| Array decimation | Reuse `PYIRENA_MAX_ARRAY_POINTS` for `get_residuals` |
| Numpy types in returns | Cast to Python floats/lists everywhere (no `np.float64`) |
| Error shape | `{"error": "...", "suggestion": "...", "code": "..."}` (string code for programmatic handling) |
| Parameter naming | `<param>_<level>` for level params (confirm in I-2) |
| Image format | PNG, base64-encoded in JSON response, default 1024×768 |
| Logging | Use existing pyirena logging where available; no print() calls in control/ |

---

## What's explicitly OUT of Phase 1

- Other fitting tools (Modeling, Size distribution, Simple fits, WAXS) — separate phases
- Utility tools (Scattering contrast, Data merge, Data manipulation) — separate phases
- Async fits / fit cancellation — deferred per [01](01-api-and-mcp-extensions.md)
- Session persistence across restarts — deferred per [01](01-api-and-mcp-extensions.md)
- Strategy DSL / playbooks (lives in Subproject 2)
- Audit trail format (lives in Subproject 2)
- GUI changes (lives in Subproject 3)
- Multi-LLM provider abstraction (lives in Subprojects 2 and 3)

---

## Open issues to resolve before / during implementation

1. **(I-1)** State persistence interaction — see investigations above
2. **(I-2)** Full parameter list and naming convention — see investigations above
3. **(I-3)** Plot rendering path — see investigations above
4. **Numerical reproducibility** — re-running `run_fit` twice with identical
   inputs must produce identical results. Audit any random-seed usage in
   `UnifiedFitModel.fit()`.
5. **Level renumbering on removal** — when removing level 1, do parameter
   names `Rg_2`, `G_2`, etc. become `Rg_1`, `G_1`? Or do IDs stay stable
   with gaps? Affects how the AI tracks history of "I just set Rg_2 to 50".
   Decide before M3.
6. **NXcanSAS save path** — does the existing `batch.py` write path
   compose, or do we need a new write helper? Investigate in M4.

---

## Rough effort estimate

This is a guess until I-1/I-2/I-3 are done. Order of magnitude:

| Milestone | Effort (focused-day equivalents) |
|-----------|----------------------------------|
| Investigations (I-1, I-2, I-3) | 1-2 days |
| M1 — sessions + read-side | 2-3 days |
| M2 — parameter control + fit | 3-4 days |
| M3 — level operations | 1-2 days |
| M4 — persistence | 1-2 days |
| M5 — schemas + MCP + docs | 2-3 days |
| End-to-end testing + iteration | 2-3 days |
| **Total** | **~12-19 days of focused work** |

Calendar time depends on interruption rate. This is *Unified Fit only* —
each subsequent tool will be a similar (likely smaller) effort.

---

## Review checklist for this plan

Before starting M1, confirm:

- [ ] Tool list covers everything needed for a complete AI-driven fit
- [ ] Tool names and signatures feel natural for an LLM to call
- [ ] Milestone breakdown is the right grain
- [ ] Investigations I-1/I-2/I-3 are the right first steps
- [ ] No critical capability is missing (e.g., q-range setting?)
- [ ] Effort estimate is acceptable for what you're getting
