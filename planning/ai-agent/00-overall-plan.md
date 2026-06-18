# 00 — Overall Plan: AI-Driven Fitting for pyirena

**Status:** Draft
**Last updated:** 2026-06-09

---

## Vision

Enable AI agents to autonomously fit SAXS/USAXS data using pyirena's existing
analysis machinery. Move beyond the current read-only AI integration
(observation only) toward bidirectional control where the AI configures
models, runs fits, evaluates results, and iterates.

Three deliverables, sequenced so each builds on the previous, but each
shippable on its own:

1. **API & MCP extensions** — the shared foundation
2. **Standalone AI app** — the main product (autonomous folder/batch fitting)
3. **In-GUI AI advisor** — a low-effort win that gives existing pyirena users
   immediate AI value

---

## Scope of pyirena tools supported

Not every pyirena capability will be wrapped for AI use. The supported set
prioritizes tools that have well-defined fitting workflows or are useful
helpers AI can invoke for the user's convenience.

### In scope — analysis tools (AI-driven fitting)
- **Unified Fit** — first target; used as the test case for all infrastructure
- **Modeling** — parameter-rich; expect the largest per-tool wrapping effort
- **Size distribution**
- **Simple fits**
- **WAXS**

### In scope — workflow / utility tools (AI invokes for user convenience)
- **Scattering contrast calculator** — AI can compute contrast when a fit
  requires it, rather than asking the user to leave the app
- **Data merge** — useful for bio users; AI can prepare combined datasets
  before fitting
- **Data manipulation** — useful for bio users; AI can clean / trim / scale
  data before fitting

### Out of scope
- **saxsMorph** — not an analysis technique
- **Fractals** — not really an analysis technique
- Real-time beamline control, AI-driven data reduction, custom ML training

### Order of implementation
Will be decided after the per-tool wrapping cost is understood from the
Unified Fit implementation (the test case). Modeling is expected to be the
largest single effort.

---

## Why split into three projects

A single monolithic "add AI to pyirena" project would be too large to ship
and too coupled to evolve. The split has three benefits:

- **Foundation reuse**: subprojects 2 and 3 both depend on subproject 1, so
  the API/MCP work pays off twice
- **Independent release cadence**: the in-GUI advisor can ship while the
  standalone app is still in design
- **Risk isolation**: AI ecosystem moves fast; the standalone app may need
  rewrites as LLMs evolve, but pyirena's own GUI and API stay stable

## Primary user workflows

### Workflow A — Autonomous batch fitting (Subproject 2)
> "I have a folder of 50 USAXS scans from last night's beam run. Fit them all
> with a Unified Fit + size distribution, save results, flag anything weird."

User opens the standalone AI app, points it at a folder, selects/configures
a fitting strategy, presses go. AI fits each dataset, produces results,
flags outliers, generates a report. Optionally runs live as data arrives.

### Workflow B — Interactive AI assistance during manual fitting (Subproject 3)
> "I'm fitting this dataset in pyirena, the fit looks off, what should I try?"

User clicks "Ask AI advisor" in the existing pyirena GUI. Current
parameters + a screenshot of the fit get sent to an LLM. AI returns
suggestions in plain language ("the high-Q region suggests adding a third
level"; "try fixing background at the measured value before refining").
User reads, decides, applies manually.

### Workflow C — Live instrument-side fitting (variant of A)
> "Data is arriving from the beamline; fit each scan as it lands."

Same as A but with a folder watcher. Probably handled as a mode of the
standalone app, not a separate project.

---

## Cross-cutting concerns (apply to all three subprojects)

### Multi-LLM support
Even if Anthropic is the launch target, design behind a thin provider
abstraction. Different labs/institutions will require OpenAI, Azure OpenAI,
or local models (Ollama / llama.cpp).

### Audit trail / reproducibility
Every AI-driven fit must produce a transcript: prompts, tool calls,
arguments, intermediate results. Scientists must be able to answer "how did
you get this fit?" Save as JSON alongside results.

### Cost transparency
Show token usage and approximate cost per session/dataset. Users care.

### Custom instructions per user/lab
System prompts that encode local domain expertise should be a first-class
config feature, not hard-coded. Different materials, different conventions.

### API key management
Use OS keyring on macOS/Windows; env vars as fallback. Never plaintext config.

### Human-in-the-loop checkpoints
Agent should be able to pause for user confirmation on destructive or
ambiguous actions ("I want to fix radius at 5 nm — OK?"). Required for trust
in the autonomous app; nice-to-have in the advisor.

---

## Recommended sequencing

```
                ┌────────────────────────────────┐
                │ Subproject 1: API/MCP          │
                │ extensions (foundation)        │
                └─────────────┬──────────────────┘
                              │
              ┌───────────────┴──────────────────┐
              │                                  │
              ▼                                  ▼
┌──────────────────────────────┐   ┌──────────────────────────────┐
│ Subproject 3: In-GUI AI      │   │ Subproject 2: Standalone AI  │
│ advisor (small, fast)        │   │ app (large, main product)    │
│                              │   │                              │
│ Ships first as a quick win   │   │ Ships in phases after the    │
│ once basic param-read tools  │   │ full control surface and     │
│ exist.                       │   │ agent loop are proven.       │
└──────────────────────────────┘   └──────────────────────────────┘
```

### Phase ordering across projects

| Phase | What happens |
|-------|--------------|
| 1 | Design API surface for control tools (Subproject 1, design only) |
| 2 | Build minimum control tools needed for the advisor (read fit, read parameters, capture image) — Subproject 1 ships partial |
| 3 | Build in-GUI advisor (Subproject 3) on top of those tools |
| 4 | Build full control surface (set/fix/free param, run fit, get residuals, etc.) — Subproject 1 complete |
| 5 | Build headless agent prototype (Subproject 2 phase 1) using Subproject 1 |
| 6 | Build standalone app GUI (Subproject 2 phase 2) |
| 7 | Distribution polish (Subproject 2 phase 3): multi-LLM, audit, config, packaging |
| 8 | Optionally expose the new control tools via MCP too (Subproject 1 phase 4) |

---

## Decisions made so far

| Decision | Rationale |
|----------|-----------|
| Pyirena will be extended with proper API-level control tools | The existing pattern (`api/` → `services/pyirena_tools.py` → `_TOOL_FUNCS` + JSON schema) is already proven; this work is the foundation everything else needs |
| The standalone AI app will be a **separate package**, not a subpackage of pyirena | Keeps pyirena's dependency footprint clean; allows independent release cadence; AI users opt in to LLM SDKs and config complexity |
| MCP support for the new control tools is *post-hoc and cheap* once API tools exist | Don't let MCP design constrain the API design |

## Decisions still open

| Question | Options | Notes |
|----------|---------|-------|
| GUI framework for the standalone app | Gradio / Chainlit / Streamlit / Qt | See [02-standalone-ai-app.md](02-standalone-ai-app.md) — leaning Gradio |
| Name of the standalone package | `pyirena-ai`, `irena-copilot`, `saxs-agent`, ... | Bikeshed later |
| First LLM provider | Anthropic / OpenAI / both | Anthropic likely first given existing experimentation, but design for multi-provider from day one |
| Where audit trail lives | JSON sidecar / SQLite / NXcanSAS extension | Defer until Subproject 2 phase 1 |

---

## Out of scope for this initiative

- Replacing pyirena's manual fitting GUI with an AI-first interface
- AI-driven data *reduction* (1D/2D reduction is upstream of pyirena)
- Real-time beamline control (we read data; we don't drive the instrument)
- Training custom models on SAXS data (we use general-purpose LLMs as agents)
- Cloud-hosted SaaS deployment (everything runs locally on user's machine)

---

## Open risks

- **API granularity** — too coarse and the AI can't fit well; too fine and
  the JSON schema explosion overwhelms the LLM's context. Will require
  iteration with real fits before locking down.
- **Model drift** — what works with Claude Opus 4.7 today may need
  reprompting when newer models arrive. Mitigate by keeping system prompts
  in version-controlled files, not hard-coded.
- **Distribution complexity** — adding a second installable package doubles
  the support surface area. Mitigate with clear "you only need pyirena-ai if
  you want X" messaging.
- **User trust** — autonomous fitting is a scary phrase for scientists.
  Audit trails and human-in-the-loop checkpoints are not optional.
