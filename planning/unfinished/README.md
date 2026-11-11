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

## 1. `dq` is stored and never used — issue #30 — **decide first**

`api/control/schemas.py` — *"Q resolution. Stored for provenance; not yet
used in fitting."* Sessions accept `dq`, the io layer round-trips it, and no
fit reads it. Only slit smearing is applied (`core/smearing.py`, which already
has the matrix machinery a resolution kernel would need).

**Finished means:** either implement pinhole/Q-resolution smearing in the fits
(the top feature candidate — see `planning/feature-candidates/`), or record
"decided against" in `PLAN.md`. Today it is neither, so every caller has to
discover that it does nothing.

## 2. Carbon model: `weighting='log'` is the same as `'relative'` — **S**

`CarbonFitModel.weighting` advertises `'auto' | 'sigma' | 'relative' | 'log'`,
and the GUI/API accept `'log'`, but `_sigma()` (`core/carbon_fit.py`) has no
branch for it: anything that is not `'sigma'` gets σ = |I|, i.e. relative
weighting.

**Finished means:** either implement log-space residuals (residual =
log I_model − log I_data, which is not just a different σ, so it touches the
residual function, not only `_sigma`), or remove `'log'` from the choices and
map old configs that carry it to `'relative'` in `from_dict()`.

---

## Not on this list, and why

| | Where it is recorded |
|---|---|
| GUI fit orchestration (`FitPanelBase`), WAXS stored-results background bug | `planning/fit-runner/README.md` — the remaining step of the fit-runner work |
| Large-file refactors, ruff `UP`/`B`/format, reader consolidation, long functions | `PLAN.md` — open, no bug |
| Batch for Data Manipulation (#19), trim/delete (#36), new models, p(r), series workflow | `planning/feature-candidates/README.md` |
| No batch/setup/control for Fractals, Merge, Manipulation | `test_tool_registration.py:ABSENT_REASON` — each with a written reason |
| ZMQ deployment (#32), async jobs (#33) | `planning/zmq-service/README.md` |
| Carbon model follow-ups (#34) | `planning/carbon-fitting/03-open-questions.md` "Left open on purpose" |
