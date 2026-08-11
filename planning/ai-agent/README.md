# AI Agent Initiative — Planning Folder

This folder captures the planning work for adding **AI-driven fitting** to the
pyirena ecosystem. It is an internal planning artifact (not user-facing
documentation), intended to be edited iteratively as decisions are refined
and tracked.

## Context

Today, pyirena exposes its analysis *results* to AI clients through an MCP
server (read-only "output" tools). AI can read fit results, summarize
folders, tabulate parameters, plot data. This works well.

The next ambition is much larger: **let an AI agent actually run the fits**,
not just read their output. The AI would configure models, set parameter
starting values and bounds, run fits, evaluate quality, iterate. This
requires giving the AI *control* tools (in addition to the existing
*observation* tools) and building infrastructure around the agentic loop.

After discussion, the initiative is split into three subprojects that share
a common foundation (the API/MCP extensions) but ship independently:

## Subprojects

| # | Plan | Summary | Status |
|---|------|---------|--------|
| 0 | [Overall plan](00-overall-plan.md) | Strategic vision, sequencing, cross-cutting concerns | Partly delivered — see status note below |
| 1 | [API & MCP extensions](01-api-and-mcp-extensions.md) | Extend `pyirena/api/` with control tools (set/fix/free parameters, run fits, get residuals). Foundation for everything else. | **Shipped** (1.1.0b5–b7) — `pyirena/api/control/` covers all five fitting tools |
| 2 | [Standalone AI app](02-standalone-ai-app.md) | New separate package (`pyirena-ai` or similar) that imports pyirena and uses an LLM to autonomously fit datasets and folders. | **Not started** — the open subproject |
| 3 | [In-GUI AI advisor](03-ai-advisor-in-gui.md) | Small additional panel in the existing pyirena GUI: grab current fit screenshot + parameters, send to LLM, display advice. Low-effort, high-value. | **Shipped** — `pyirena/gui/ai_advisor.py` |
| 1a | [Phase 1 detailed plan: Unified Fit control](phase-1-unified-fit-control.md) | Concrete implementation plan for Subproject 1 with Unified Fit as the test case. Tool catalog, milestones, investigations. | **Shipped** — then extended to Sizes, Simple Fits, Modeling and WAXS |

## Status as of 1.1.0b7 (2026-08-10)

This folder is **kept** because the initiative is not finished, but much of it
has landed and the plans have not been rewritten to match.  Read them as the
original intent, not as the current state.

**Landed**

- The control surface (subproject 1) exists for all five fitting tools —
  Unified Fit, Size Distribution, Simple Fits, Modeling and WAXS Peak Fit —
  under `pyirena/api/control/`, exposed as ~119 MCP tools.
- Fit-quality metrics (plans 04/05) shipped: `core/fit_metrics.py`,
  `gui/quality_display.py`, `io/nxcansas_fit_quality.py`.
- The in-GUI advisor (subproject 3) shipped as `gui/ai_advisor.py`.
- Setup state is embedded in result files (`_pyirena_config`) so an agent run
  can be reopened and continued in the GUI — six tools, see
  `docs/HDF5_NxcanSAS_structure.md`.

**Open — worth a fresh look before more work**

- **Subproject 2, the standalone agent app**, is untouched.  It is the piece
  that would actually close the loop, and it is also the piece whose design
  assumptions are oldest.
- Whether the control surface is *sufficient* for autonomous fitting has not
  been tested end to end by an agent working unattended; the tools exist, the
  agentic loop around them does not.
- `mcp` is pinned to `<2`.  Migrating to 2.x is its own project and should be
  decided before building anything new on top.
- The remaining tools (SAXS Morph, Fractals, Contrast, Merge, Manipulation)
  have no control surface, and it is not obvious that they should — decide per
  tool rather than for completeness.

## Reading order

If you're new to this, read `00-overall-plan.md` first — it explains how the
three subprojects fit together and which one to tackle first. Then dive into
the individual plans.

## Status conventions

- **Draft** — initial brain dump, may have gaps and open questions
- **Reviewed** — discussed and refined, no major changes pending
- **Approved** — ready to begin implementation
- **In progress** — implementation under way, link to branch/PR
- **Done** — shipped; the plan can be moved to `planning/archive/`

## Related project docs

- [docs/ai_integration.md](../../docs/ai_integration.md) — current MCP server (read-only)
- [docs/ai_tools_reference.md](../../docs/ai_tools_reference.md) — existing MCP tool catalog
- [pyirena/api/README.md](../../pyirena/api/README.md) — API layer that MCP wraps
