# AI Agent Initiative — Planning Folder

Internal planning artifact for AI-driven fitting in the pyIrena ecosystem —
not user-facing documentation.

**Status: one subproject left.** The initiative was originally three
subprojects sharing one foundation. Two of them shipped; their plans have been
deleted (git history has them) because a shipped plan is worse than no plan —
it describes intent that the code has since overtaken. What survives here is
the part that has *not* been built, plus the decisions and open questions that
outlive the plans.

| Subproject | Status |
|---|---|
| **1 — API & MCP control surface** | **Shipped** (1.1.0b5–b7). `pyirena/api/control/` covers all five fitting tools, exposed through `pyirena/mcp/dispatch.py`. |
| **2 — Standalone AI app** | **Not started.** The one open piece — see [02-standalone-ai-app.md](02-standalone-ai-app.md). |
| **3 — In-GUI AI advisor** | **Shipped.** `pyirena/gui/ai_advisor.py`. |
| Fit-quality metrics (a dependency of all three) | **Shipped.** `core/fit_metrics.py`, `gui/quality_display.py`, `io/nxcansas_fit_quality.py`; documented in [docs/fit_quality_metrics.md](../../docs/fit_quality_metrics.md). |

## What shipped

- **Control surface** for Unified Fit, Size Distribution, Simple Fits,
  Modeling and WAXS Peak Fit under `pyirena/api/control/`: open a dataset,
  configure a model, set/fix/free parameters, run the fit, read residuals and
  quality metrics, save to NXcanSAS.
- **Utility tools** the agent can invoke for the user's convenience:
  scattering contrast (`pyirena.api.calculators`), data merge
  (`merge_datasets` / `match_merge_files`) and data manipulation
  (`average_data` / `subtract_data` / `divide_data` / `scale_data` /
  `trim_data` / `rebin_data`).
- **Robust fit-quality metrics** — σ-scale-independent diagnostics, so an agent
  can tell a mis-scaled uncertainty from a genuine misfit.
- **In-GUI advisor** — screenshot + parameters → LLM → plain-language advice.
- **Setup state embedded in result files** (`_pyirena_config`), so an agent run
  can be reopened and continued interactively in the GUI.

## What is open

- **Subproject 2, the standalone agent app**, is untouched — the piece that
  would actually close the loop, and the one whose design assumptions are
  oldest. Treat [02-standalone-ai-app.md](02-standalone-ai-app.md) as a
  starting point to argue with, not a spec.
- **Whether the control surface is *sufficient*** for autonomous fitting has
  never been tested end to end by an agent working unattended. The tools exist;
  the agentic loop around them does not.
- **`mcp` is pinned to `<2`.** Migrating to 2.x is its own project and should
  be decided before anything new is built on top.
- **The remaining tools** (SAXS Morph, Fractals, Contrast, Merge, Manipulation)
  have no *control* surface. Decide per tool on merit rather than for
  completeness — SAXS Morph and Fractals are visualization, not analysis
  techniques, and were ruled out of scope at the start.

## Decisions that still stand

| Decision | Rationale |
|---|---|
| The standalone AI app is a **separate package**, not a pyIrena subpackage | Keeps pyIrena's dependency footprint clean; independent release cadence; AI users opt into LLM SDKs and config complexity |
| MCP exposure is **post-hoc and cheap** once API tools exist | Don't let MCP design constrain the API design — this held up in practice |
| Control tools live at the **api layer**, not in the GUI | Same surface serves the advisor, the standalone app, scripting and MCP |

Still open: GUI framework for the standalone app (leaning Gradio), package
name, first LLM provider, and where the audit trail lives (JSON sidecar /
SQLite / NXcanSAS extension).

## Cross-cutting requirements — apply to anything built here

- **Multi-LLM from day one.** Even with Anthropic as the launch target, design
  behind a thin provider abstraction; labs will require OpenAI, Azure or local
  models.
- **Audit trail.** Every AI-driven fit produces a transcript — prompts, tool
  calls, arguments, intermediate results. A scientist must be able to answer
  "how did you get this fit?"
- **Cost transparency.** Token usage and approximate cost per session.
- **Custom instructions per user/lab** as a first-class config feature, in
  version-controlled files, not hard-coded system prompts.
- **API keys** in the OS keyring, env vars as fallback, never plaintext config.
- **Human-in-the-loop checkpoints** — the agent pauses for confirmation on
  destructive or ambiguous actions. Not optional for trust.

## Explicitly out of scope

Replacing the manual fitting GUI with an AI-first interface; AI-driven data
*reduction* (upstream of pyIrena); real-time beamline control; training custom
models on SAXS data; cloud-hosted SaaS.

## Open risks

- **API granularity** — too coarse and the agent cannot fit well; too fine and
  the schema explosion overwhelms its context. Needs iteration against real
  fits.
- **Model drift** — prompts that work today need revisiting as models change.
  Keep them in version-controlled files.
- **Distribution complexity** — a second installable package doubles the
  support surface. Mitigate with clear "you only need this if…" messaging.
- **User trust** — "autonomous fitting" is a scary phrase for scientists. Audit
  trails and checkpoints are the answer.

## Related project docs

- [docs/ai_integration.md](../../docs/ai_integration.md) — the MCP server
- [docs/ai_tools_reference.md](../../docs/ai_tools_reference.md) — MCP tool catalog
- [pyirena/api/README.md](../../pyirena/api/README.md) — the API layer MCP wraps
- [docs/fit_quality_metrics.md](../../docs/fit_quality_metrics.md) — the metrics an agent judges fits by
