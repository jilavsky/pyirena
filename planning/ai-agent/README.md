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
| 0 | [Overall plan](00-overall-plan.md) | Strategic vision, sequencing, cross-cutting concerns | Draft |
| 1 | [API & MCP extensions](01-api-and-mcp-extensions.md) | Extend `pyirena/api/` with control tools (set/fix/free parameters, run fits, get residuals). Foundation for everything else. | Draft |
| 2 | [Standalone AI app](02-standalone-ai-app.md) | New separate package (`pyirena-ai` or similar) that imports pyirena and uses an LLM to autonomously fit datasets and folders. | Draft |
| 3 | [In-GUI AI advisor](03-ai-advisor-in-gui.md) | Small additional panel in the existing pyirena GUI: grab current fit screenshot + parameters, send to LLM, display advice. Low-effort, high-value. | Draft |
| 1a | [Phase 1 detailed plan: Unified Fit control](phase-1-unified-fit-control.md) | Concrete implementation plan for Subproject 1 with Unified Fit as the test case. Tool catalog, milestones, investigations. | Draft — for review |

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
