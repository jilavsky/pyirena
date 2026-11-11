# Unfinished work — half-built things still in the tree

Internal planning artifact. First swept 29-09-2026; re-swept 11-11-2026 after
the 2026-11 review (issue #38) was merged. Finished items are removed — the
CHANGELOG records what shipped.

This is **not** a wish list (that is `planning/feature-candidates/`) and not
`PLAN.md` (open refactors and what was *decided against*). This file is the
third category: **work that got most of the way and was left, where the tree
and the documentation disagree** — a reader cannot tell these from finished
work without running them.

Each item says what exists, what is missing, and what "finished" means. Sized
**S** ≈ one session, **M** ≈ two to three.

---

**Nothing open (11-11-2026).** The 2026-11 leftovers were fixed on
`fix/unfinished-items`; dQ (#30) was deferred as a feature and moved to
`planning/feature-candidates/`. Add the next half-built thing here.

---

## Not on this list, and why

| | Where it is recorded |
|---|---|
| GUI fit orchestration (`FitPanelBase`), WAXS stored-results background bug | `planning/fit-runner/README.md` — the remaining step of the fit-runner work |
| Large-file refactors, ruff `UP`/`B`/format, reader consolidation, long functions | `PLAN.md` — open, no bug |
| dQ / Q-resolution smearing (#30, deferred 11-11-2026), Batch for Data Manipulation (#19), trim/delete (#36), new models, p(r), series workflow | `planning/feature-candidates/README.md` |
| No batch/setup/control for Fractals, Merge, Manipulation | `test_tool_registration.py:ABSENT_REASON` — each with a written reason |
| ZMQ deployment (#32), async jobs (#33) | `planning/zmq-service/README.md` |
| Carbon model follow-ups (#34) | `planning/carbon-fitting/03-open-questions.md` "Left open on purpose" |
