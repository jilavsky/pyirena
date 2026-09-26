# Config serialisation: one vocabulary, not three

Internal planning artifact. Written 25-09-2026 against pyIrena 1.1.1, prompted
by a question while building the ZMQ service: *why does Unified Fit — the
simplest tool to configure — have the most trouble being scripted?*

The short answer is that it has nothing to do with Unified Fit being complex.
It is the oldest tool, its parameter vocabulary was fixed before `to_dict()`
existed, and that vocabulary is now written into every saved HDF5 file and
documented as the public batch config format. A migration was started and
finished for Modeling but not for Unified Fit.

**Status: not started. Independent of the ZMQ work; do not put it on that
branch.**

---

## 1. What is actually wrong

Two separate problems that look like one.

### Problem A — `to_dict()` is incomplete for some tools

Measured by comparing each model's scalar attributes against its `to_dict()`
output (fit results and data arrays correctly excluded):

| Tool | Settings `to_dict()` drops |
|---|---|
| Modeling | none |
| Carbon model | none |
| Unified Fit | `use_analytic_jacobian` |
| Simple Fits | `use_analytic_jacobian` |
| WAXS Peak Fit | `use_analytic_jacobian` |
| **Size Distribution** | **`use_slit_smearing`, `slit_length`**, plus `use_analytic_jacobian` |

`use_analytic_jacobian` looks like a harmless performance switch. In Simple
Fits it is not — see Problem C. That is the argument for the contract test in
Step 2: "probably harmless to omit" is not a judgement worth making tool by
tool, by hand.

**The Sizes entry is a bug, and it is a scientific one.** Demonstrated:

```python
>>> m = SizesDistribution(); m.use_slit_smearing = True; m.slit_length = 0.0217
>>> SizesDistribution.from_dict(m.to_dict()).use_slit_smearing
False
```

A slit-smeared size-distribution setup, saved and reloaded, silently comes
back as a pinhole fit. It does not error; it returns a different, plausible
answer. Unified Fit round-trips both fields correctly, which is what makes
this an inconsistency rather than a design choice.

### Problem C — `SimpleFitModel.from_dict()` produces an object that cannot fit

Found 26-09-2026 by running the new `testData/Scripting/` fixtures. Not a
dialect problem; the same root cause.

```python
>>> m = SimpleFitModel.from_dict({"model": "Guinier", "params": {}})
>>> m.fit(q, I, dI)
AttributeError: 'SimpleFitModel' object has no attribute 'use_analytic_jacobian'
```

`from_dict()` builds the object with `cls.__new__(cls)` and assigns fields one
by one, bypassing `__init__`. `use_analytic_jacobian` is set in `__init__`,
absent from `to_dict()`, and never assigned in `from_dict()` — so every model
rebuilt from a dict is missing it, and `fit()` reads it unconditionally.

**`pyirena.batch.fit_simple_from_config()` therefore fails for every real
model.** It is the documented scripting entry point for Simple Fits.

The test suite misses it because the only tests driving that path
(`test_invariant.py`) use the `Invariant` model, which returns before reaching
the Jacobian branch.

Fix: assign it in `from_dict()` with a `True` default, and add it to
`to_dict()`. Two lines. The interesting part is not the fix but that
`cls.__new__` + hand-assignment is a pattern that silently drops any field
added later — worth checking the other tools use the constructor instead
(they do; Simple Fits is the only one that loses an attribute outright).

### Problem B — Unified Fit has two vocabularies

Everyone else has one.

| | Core dialect (`to_dict`) | Panel dialect (GUI export, batch config, `_pyirena_config`) |
|---|---|---|
| A parameter | `"G": 1000.0` | `"G": {"value": 1000.0, "fit": true, "low_limit": 1e8, "high_limit": 1e12}` |
| Cutoff | `RgCO` | `RgCutoff` |
| Correlations | `correlations` | `correlated` |
| Link B | `link_B` | `estimate_B` |
| Link RgCO | `link_RGCO` | `link_rgco` |

`pyirena/batch/unified.py:_state_to_model` raises `AttributeError: 'float'
object has no attribute 'get'` when handed a core dict. `UnifiedLevel.
from_panel_params` translates one way (panel → core); there is no inverse, so
`api/control/unified_fit.py:_session_to_gui_state` hand-builds the panel dict
when it needs to save a setup.

For every other tool the two are the same thing. Verified for Modeling
against a real exported file (`testData/Core-shell-tests/pyirena_config.json`):
the population keys are identical to `to_dict()`, and the GUI simply omits
five top-level keys it has no widget for.

## 2. Why Unified Fit specifically

Not complexity — age, plus an unfinished migration.

`to_dict()` / `from_dict()` were added to `core/unified.py` and
`core/modeling.py` in the **same commit** (`dbffaf1`, "U6 and U10 fixes",
10-08-2026). That commit also rewrote the panels to use them —
`gui/modeling_panel.py` lost 231 lines, `gui/unified_fit.py` only 132, and
Unified kept four `from_panel_params` call sites. Modeling's panel came out
speaking the core dict; Unified's did not.

It could not, and that part was a reasonable decision at the time. By then the
panel vocabulary was already:

- **inside every saved file** — `io/setup_config.py:write_setup_config` embeds
  the panel state as `_pyirena_config` in the results group, and *Load Setup
  from File…* reads it back;
- **in every user's state file** — `state/state_manager.py` stores the Unified
  defaults with `RgCutoff`, `correlated`, `estimate_B`;
- **a documented public format** — `docs/batch_api.md` specifies it as *the*
  batch config, so users have hand-written files in it.

Changing it needed a compatibility bridge. The bridge was never built, so both
vocabularies stayed, and `from_panel_params` became the seam between them.

Modeling looks better here not because it was designed better, but because it
was young enough to be changed outright.

## 3. What "synchronised" should mean

Three rules, in priority order.

1. **`to_dict()` / `from_dict()` is the one canonical dialect.** Complete —
   every setting that changes a number is in it — and the only dialect
   pyIrena *writes* from now on.
2. **The panel dialect becomes a legacy input format.** Read forever, never
   written again. A user's existing config file and every saved `.h5` keep
   working; that is non-negotiable and is what backwards compatibility means
   here. Forward compatibility (an *older* pyIrena reading a *new* file) has
   never been promised and is not in scope.
3. **A setup is model config *plus* session settings.** The fitted Q range,
   slit smearing and background sub-ranges are not model fields and should not
   be forced into `to_dict()`. They need a defined envelope instead, the same
   shape for every tool:

   ```json
   {"tool": "unified_fit",
    "model": { ... to_dict() ... },
    "fit_q_range": {"q_min": 0.001, "q_max": 0.1},
    "data": {"is_slit_smeared": false, "slit_length": 0.0}}
   ```

   This is what `analyze()` (ZMQ Phase 3) should take, what *Export
   Parameters* should write, and what `pyirena.batch` should accept alongside
   the legacy format. Today each tool improvises: Sizes carries
   `cursor_q_min/max` in its batch section, WAXS carries `q_min`/`q_max`,
   Unified carries `cursor_left`/`cursor_right`. Three spellings of one idea.

## 4. Plan

Four steps, each shippable alone and in this order. Sizes are **S** ≈ one
session, **M** ≈ two to three.

### Step 1 — Fix the two serialisation bugs (S, do this first and separately)

1. **Simple Fits crash (Problem C).** Assign `use_analytic_jacobian` in
   `SimpleFitModel.from_dict()` and add it to `to_dict()`. This one is
   urgent: a documented scripting entry point raises `AttributeError` today.
2. **Sizes slit loss (Problem A).** Add `use_slit_smearing` and `slit_length`
   to `SizesDistribution.to_dict()` and read them in `from_dict()`. Additive,
   so old files still load.

Both need a round-trip case in `test_core_serialization.py`, and both are
independent of everything below. Worth their own small branch — one is a
crash, the other silently changes results.

### Step 2 — A completeness contract (S)

Two tests, in the spirit of `test_tool_registration.py` and
`test_gui_state_contract.py`:

1. For every fitting tool, every scalar/bool model attribute is either present
   in `to_dict()` or named in an explicit `_NOT_SERIALISED` allowlist on the
   class, with a comment saying why.
2. For every fitting tool, `from_dict(to_dict())` produces an object that can
   actually **run a fit** — not merely one that compares equal. Problem C
   would have been caught the day it was introduced by nothing more than
   that.

Together these turn "someone happened to run a fixture" into "the build
notices".

### Step 3 — The setup envelope (M)

Define the `{tool, model, fit_q_range, data}` envelope of §3.3 in one place —
`core/setup.py` or `api/setup.py` — with `build_setup(model, session)` and
`apply_setup(setup) -> (model, settings)`. Then:

- `export_results` returns it instead of a bare `config`;
- `analyze()` (ZMQ Phase 3) takes it;
- *Export Parameters* writes it;
- `pyirena.batch` accepts it **in addition to** the legacy per-tool sections.

This removes the Sizes and WAXS "missing key" gaps by putting those settings
where they belong rather than by pushing them into the model.

### Step 4 — Retire the Unified panel dialect as an output (M)

The actual migration, and the only step that touches the GUI.

1. Add `UnifiedLevel.to_panel_params()`, the inverse of `from_panel_params`,
   and make `api/control/unified_fit.py:_session_to_gui_state` use it instead
   of hand-building the dict. One translator, both directions, one place.
2. Make every *reader* accept both dialects, detected by shape — a parameter
   that is a dict rather than a number is the panel dialect. That is
   `batch/unified.py:_state_to_model`, the setup loader, and `analyze`.
   Route through `from_panel_params`; do not write a second translator.
3. Switch the Unified panel's `_collect_state()` to the core dialect, and the
   StateManager defaults with it. `_pyirena_config` then contains the core
   dialect for new files; old files still load through the reader from (2).
4. Update `docs/batch_api.md` to document the core dialect as current and the
   panel dialect as accepted-legacy, with the mapping table from §1.

**Do not** attempt to rename the core fields to the panel's names instead.
`RgCO`/`correlations`/`link_B` are the names in the maths and in the
docstrings, and the panel names are the historical accident.

## 5. Risks

- **Numerical regression is the only unacceptable outcome.** None of this
  changes maths, but Step 4 changes what the GUI restores from a setup, and a
  setup restored wrong produces a different fit. Every step needs a
  before/after on `validationData/` (`run_validation_report.py`) showing the
  numbers unchanged, not just a green test suite.
- **Step 4 touches `gui/unified_fit.py`**, the largest file in the package
  (4442 lines) and the one users exercise most. It should be its own branch
  and its own release note.
- **Fractals reuses Unified levels** (`core/fractals.py`,
  `gui/fractals_panel.py` both speak the panel dialect) — it must be migrated
  in the same step or explicitly left on the legacy reader.
- Cost of doing nothing is low and steady: one extra translator on the
  Unified path, a documented quirk in the ZMQ service, and the next tool that
  reads a config has to learn both. Nothing breaks. This is cleanup, not a
  fire — except Step 1, which is a bug.

## 6. What this is not

Not a rewrite of the state system, not a change to NXcanSAS layout, not a
change to any result group. `to_dict()` already exists and already works for
four tools; this is finishing a migration that is two-thirds done.
