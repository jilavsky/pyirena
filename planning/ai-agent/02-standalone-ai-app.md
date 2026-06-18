# 02 — Standalone AI App (Main Product)

**Status:** Draft
**Last updated:** 2026-06-09
**Depends on:** [01 API & MCP extensions](01-api-and-mcp-extensions.md) — needs the full control surface
**Related:** [03 In-GUI AI advisor](03-ai-advisor-in-gui.md) — shares LLM provider layer

---

## Goal

A separate installable package (working name: **`pyirena-ai`**) that imports
pyirena and lets users have an LLM **autonomously fit datasets and folders**
of SAXS/USAXS data. The main product of this initiative.

Distinct from pyirena itself: pyirena is a focused scientific library; this
package is the AI-driven workflow consumer.

The set of pyirena tools this app drives is defined in
[00-overall-plan.md — Scope of pyirena tools supported](00-overall-plan.md#scope-of-pyirena-tools-supported).
Unified Fit is the test case; the app inherits whatever support
[01](01-api-and-mcp-extensions.md) ships.

> **Note:** detailed decisions for this subproject (GUI framework, package
> name, audit format, etc.) will be revisited once subprojects 1 and 3 are
> further along and we have real experience with the agent loop. The plan
> below is a placeholder structure.

---

## Primary user workflows

### Workflow A — Folder batch fitting
> "Fit all 50 USAXS scans in this folder with Unified Fit + size distribution."

User opens the app, drags in a folder, selects a strategy (or accepts the
default), presses go. AI fits each dataset, produces results, flags
outliers, generates a summary report.

### Workflow B — Live instrument-side fitting
> "Watch this folder; fit each new scan as it arrives."

Same as A but with a folder watcher. Files appear → AI fits → results
update. Useful at beamlines during data collection.

### Workflow C — Interactive AI co-pilot for one dataset
> "Help me fit this tricky dataset; I want to watch your reasoning."

User uploads a single dataset; conversation-style interface. AI shows its
tool calls, results, decisions. User can interject, correct, override.
Bridges into Workflow A by saving the conversation as a reusable strategy.

---

## Architecture

```
pyirena-ai/
├── pyproject.toml             # depends on pyirena, anthropic, gradio (or chosen GUI), ...
├── pyirena_ai/
│   ├── __init__.py
│   ├── core/
│   │   ├── agent.py           # tool-use loop (LLM-agnostic)
│   │   ├── tools.py           # JSON schemas + dispatch → pyirena.api
│   │   ├── strategy.py        # System prompts, fitting playbooks
│   │   ├── session.py         # Conversation + audit trail
│   │   └── watcher.py         # Folder-watch for live mode
│   ├── llm/
│   │   ├── base.py            # Provider interface
│   │   ├── anthropic.py
│   │   ├── openai.py
│   │   └── local.py           # Ollama / llama.cpp (later phase)
│   ├── config/
│   │   ├── settings.py        # API keys (via keyring), model, custom instructions
│   │   └── strategies/        # Built-in fitting playbooks (.md files)
│   ├── gui/
│   │   ├── app.py             # Gradio / Chainlit / Streamlit entry point
│   │   ├── plots.py           # Render pyirena results
│   │   ├── conversation.py    # Chat / agent-step display
│   │   └── components.py
│   └── cli/
│       └── main.py            # Headless CLI for batch / scripting
└── tests/
```

### Layers
1. **Provider layer (`llm/`)** — thin interface over Anthropic, OpenAI, local.
   Handles tool-use loop primitives (send message, receive tool calls, send
   tool results). The rest of the app is LLM-agnostic.
2. **Agent layer (`core/`)** — orchestrates the loop. Owns the conversation,
   dispatches tool calls to `pyirena.api`, manages sessions, writes audit
   trail.
3. **Tool layer (`core/tools.py`)** — registers each `pyirena.api` function
   as an LLM-callable tool with JSON schema. The single bridge between the
   AI world and the pyirena world.
4. **GUI / CLI layer** — two front-ends sharing the same agent core.
   Headless CLI for batch and scripting; GUI for interactive use.

---

## GUI framework decision

| Framework | Strengths for this app | Weaknesses |
|-----------|------------------------|------------|
| **Gradio** (recommended) | Built for AI/ML demos; native chat, streaming, tool calls, file upload, image display; deploys like Streamlit | Less flexible than Streamlit for complex non-chat layouts |
| **Chainlit** | Purpose-built for LLM chat; agent step visualization (every tool call shown); great for "watch AI fit" UX | Smaller community; chat-only shape |
| **Streamlit** | Familiar to user (uses elsewhere); large ecosystem | Rerun-on-interaction model fights against streaming chat / agent loops; more plumbing |
| **PySide/PyQt** | True desktop app; tighter integration with pyirena's existing Qt GUI possible | Slow to develop chat UIs; packaging painful (PyInstaller, code signing) |

**Leaning: Gradio.** Native chat + agent step visualization + image
rendering + simple deploy, with enough flexibility to add result tables and
folder-progress views. Falls back to Streamlit if Gradio's layout
constraints become limiting.

Final decision deferred until after building a 1-day spike in both.

---

## Phasing

### Phase 1 — Headless agent prototype (small)
**Goal:** prove the agent loop with no GUI distraction.

- CLI: `pyirena-ai fit data.h5 --model unified --strategy default.md`
- One LLM provider (Anthropic), one model (Unified Fit), one strategy
- Produces: fitted HDF5 + JSON transcript of all tool calls
- No multi-LLM, no GUI, no folder mode

Validates: agent loop works, tool granularity is right, system prompts make
sense. Useful immediately for power users.

### Phase 2 — Minimal GUI (medium)
**Goal:** chosen framework + basic interactive use.

- File upload → chat with AI → see plots/tables → export
- Display current parameters live as AI changes them
- Show each tool call as an "agent step" visible to user
- Cancel / pause / approve checkpoints

Validates: GUI framework choice, UX patterns.

### Phase 3 — Batch + folder mode (medium)
**Goal:** the primary product workflow.

- Folder picker; AI fits each file
- Progress display: which file is being fit, how many done, how many failed
- Final summary report with parameter trends across the folder
- Folder watcher mode for live instrument-side fitting

### Phase 4 — Distribution polish (medium)
**Goal:** ship-ready for 10-100 users/month.

- Multi-LLM (OpenAI + Anthropic; local optional)
- API key management via OS keyring
- Custom instructions per lab (loaded from `~/.pyirena-ai/instructions.md`)
- Cost transparency (tokens / dollars per session)
- Audit trail format finalized (JSON sidecar alongside HDF5)
- Conda + pip packaging
- Installation and configuration docs

### Phase 5 — Strategy library (small, ongoing)
**Goal:** ship reusable fitting playbooks.

- Built-in strategies for common cases ("Unified Fit, 2 levels, mass
  fractal", "Sizes distribution, MaxEnt, log-spacing")
- User-saveable strategies from successful conversations
- Strategy = system prompt fragment + initial tool-call sequence

---

## Cross-cutting requirements

(Repeating items from `00-overall-plan.md` that apply specifically here.)

- **Audit trail**: every fit produces a JSON sidecar listing each tool call,
  arguments, results, χ² evolution. Stored next to the HDF5.
- **Multi-LLM from day one of phase 4**: don't lock to Anthropic-only.
- **Custom instructions**: lab-specific system prompts as a first-class
  config field.
- **Human-in-the-loop checkpoints**: agent can pause for confirmation on
  destructive actions or when it's uncertain.
- **Cost transparency**: tokens used + estimated cost per session/file.

---

## Open questions

| Question | Notes |
|----------|-------|
| Package name | `pyirena-ai`, `irena-copilot`, `saxs-agent`, ... |
| GUI framework final choice | Spike Gradio + Chainlit before committing |
| How to handle long-running fits | Streaming progress vs blocking? Likely OK to block initially since most pyirena fits are seconds-to-minutes |
| Audit trail format | JSON sidecar likely; could extend NXcanSAS schema instead |
| Where to store user strategies | `~/.pyirena-ai/strategies/` directory of markdown files |
| Folder watcher backend | `watchdog` library is the standard choice |
| How to handle fit failures in batch mode | Skip + log? Retry with different starting conditions? Ask AI to diagnose? |
| Cost guardrails | Per-session token cap? Per-user monthly budget warning? |
| Distribution channel | pip + conda? PyInstaller bundle for non-developers? |

---

## Out of scope

- Replacing pyirena's manual GUI
- Cloud-hosted SaaS (everything runs locally)
- Training custom models
- Multi-user / shared workspace features
- Real-time instrument control beyond reading the data folder

---

## Success criteria

1. User can install via `pip install pyirena-ai`
2. User can run `pyirena-ai fit folder/` and get reasonable fits for all files
3. User can open the GUI, drag a folder, watch the AI fit each file, and
   trust the results enough to use them
4. Every fit has an audit trail explaining how the AI arrived at the result
5. Cost per dataset is transparent and bounded
6. Works with at least Anthropic and OpenAI providers
