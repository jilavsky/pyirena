# Unfinished work — the beta clean-up list

Internal planning artifact. Swept 29-09-2026 on `feature/zmq-service`
(pyIrena 1.1.1, Carbon model included), while reviewing the config dialects.

This is **not** a wish list and not the same thing as `PLAN.md`. `PLAN.md`
records what is open *and deliberately decided against*; `test_tool_registration.py`
records per-surface gaps *with a written reason*. Both of those are healthy.
This file is the third category: **work that was started, got most of the way,
and was left — where the half-built thing is still in the tree and nothing
says so.** A reader today cannot tell these apart from finished work without
running them.

Each item says what exists, what is missing, and what "finished" means. Sized
**S** ≈ one session, **M** ≈ two to three.

Ordered by what should close before the beta.

---

## 1. Config dialects — four silent-wrong-answer bugs — **M, blocks beta**

Fully worked out in **`planning/config-dialects/README.md`** (reviewed and
re-measured 29-09-2026). Not repeated here. The short version:

- `export_results()` → `analyze()` drops every Unified Fit per-parameter fit
  flag and every bound;
- a Modeling setup read back from HDF5 loses the parameters of every
  non-size-distribution population;
- a Size Distribution config replays a spheroid as a sphere, and a
  slit-smeared setup as a pinhole one;
- of the three config envelopes pyIrena *writes*, `build_setup()` accepts one.

All four return a plausible number rather than an error. **Step 0 of that
plan is the beta blocker**; Steps 1–4 are cleanup and can follow.

---

## 2. Monte-Carlo uncertainties are missing from the agent surface — **M**

The largest capability gap between the batch path and the MCP/ZMQ path, and
the one an agent notices first.

**What exists.** Every one of the six `pyirena.batch` fit functions takes
`with_uncertainty: bool` and `n_mc_runs: int`:
`batch/unified.py` (`_mc_uncertainty_unified`), `batch/sizes.py`
(`_mc_uncertainty_sizes`), `batch/simple.py`, `batch/modeling.py`,
`batch/carbon_fit.py`, `batch/waxs.py`.

**What is missing.** In `pyirena/api/control/`, only Carbon exposes it:

| control function | uncertainty parameter |
|---|---|
| `carbon_fit.run_carbon_fit(session_id, weighting, **n_mc_runs**)` | yes |
| `unified_fit.run_fit(session_id, max_iter, tolerance, random_seed, walk_limits)` | **no** |
| `sizes.run_sizes_fit(session_id, random_seed)` | **no** |
| `simple_fits.run_simple_fit(session_id, no_limits)` | **no** |
| `modeling.run_modeling_fit(session_id, fit_method)` | **no** |
| `waxs_peakfit.run_waxs_fit(session_id, weight_mode)` | **no** |

`api/control/unified_fit.py:get_parameter_uncertainties()` is a hard-coded
Phase-1 placeholder that returns `{"available": False, "note": "... not
computed in Phase 1 ..."}`, advertised as such in
`api/control/schemas.py:696` and in `api/control/export.py:26`. It exists only
for Unified Fit; the other five tools have no equivalent at all.

So a scientist running `pyirena.batch` gets error bars and the same scientist
driving the ZMQ service does not. `analyze()` already has the right instinct
here — it takes an uncertainty opt-in and reports the run count back, because
MC is the one setting that outruns a 60 s client timeout — so the design
question is answered; only the wiring is missing.

**Finished means:** `with_uncertainty` / `n_mc_runs` on all six `run_*_fit`
control functions, reusing the batch helpers rather than reimplementing them;
`get_parameter_uncertainties` generalised to all six and returning real
numbers; the schemas and `docs/ai_tools_reference.md` updated. Keep the
default off — the timeout reasoning in `planning/zmq-service/README.md` §13
still holds.

**Related, smaller:** `batch/waxs.py:29-38` accepts `with_uncertainty` and
`n_mc_runs` "for API consistency" and **does not use them** — WAXS
uncertainty comes from the `curve_fit` covariance. A signature that accepts a
parameter and ignores it is worse than one that omits it. Either return the
covariance-derived values under the same key, or drop the two arguments.

---

## 3. SAXS Morph: the GUI has no Fit, and half of Phase 4 is in the tree — **M**

**Done 01-10-2026** (issue #29). Resolved the second way: SAXS Morph is a
visualisation tool, so the fit was removed rather than wired. `Engine.fit`,
`calculate_uncertainty_mc`, `MAX_FIT_VOXEL_SIZE`, `voxel_size_fit`, the
`fit_*`/`*_limits`/`no_limits`/`n_mc_runs` config fields, both dead GUI
worker classes and the never-instantiated `ParamRow` are all gone. The
`(saxs_morph, viewer)` registration gap is now a written decision.

The clearest "stopped mid-way and never came back" in the package.

**What exists.** The maths is real and finished: `core/saxs_morph.py` fits via
`least_squares` / `minimize` (around line 1300), caps the fit voxel grid at
`MAX_FIT_VOXEL_SIZE = 256`, and `batch/saxs_morph.py:fit_saxs_morph` drives it
from a config. The tool is fully registered everywhere it should be
(`test_tool_registration.py`).

**What is missing.** Three things, all in `gui/saxs_morph_panel.py`:

1. The module docstring still says *"Fit/Cancel/MC uncertainty buttons are
   constructed but disabled until Phase 4 wires the QThread workers."*
   Measured: those buttons are **not constructed at all** — there is no
   `btn_fit`, `btn_mc` or `btn_cancel`. The panel offers *Graph Model* and the
   two background pre-fits (`btn_pl_fit`, `btn_bg_fit`) and nothing else. The
   docstring is stale in the direction that overstates completeness.
2. `_FitWorker` (line 654) and `_MCWorker` (line 691) are **defined and never
   instantiated** anywhere in the file. Phase 4's workers were written; the
   buttons that would start them were not. (`modeling_panel.py` has its own
   identically-named classes and does use them — do not be misled by a grep.)
3. `voxel_size_fit` is a dead knob. It is a field on `SAXSMorphConfig`, it is
   in `StateManager.DEFAULT_STATE['saxs_morph']`, and `batch/saxs_morph.py:145`
   reads it — but the panel's `_collect_state()` never writes it, and
   `saxs_morph_panel.py:1385` forces `voxel_size_fit=int(st['voxel_size_render'])`
   with the comment *"same — fit loop is not used"*. So the batch path honours
   a setting the GUI cannot express and never saves.

**Finished means:** either wire Phase 4 (buttons → the two workers that are
already written, plus the `voxel_size_fit` control), or delete the two dead
worker classes, remove `voxel_size_fit` from the GUI-facing state, and rewrite
the docstring to say the GUI is model-evaluation-only and fitting is a batch
operation. **Deciding which is the work**; both are small once decided. Do not
leave it as it is — the current state reads as a bug in the GUI.

**Also here:** `("saxs_morph", "viewer")` is the one entry in
`test_tool_registration.py:ABSENT_REASON` that says *"Open gap rather than a
decision"*. It is already honestly labelled, so it is not urgent — but if
SAXS Morph is not going to get trend plots, turn it into a decision.

---

## 4. Sizes `unc_n_runs` is collected, saved, exported and ignored — **S**

`unc_n_runs` is a real spin box (`gui/sizes_panel.py:1729`), it drives the GUI
uncertainty run (`sizes_panel.py:2478`), it is written into the panel state
(`:2729`), restored from it (`:2786`), and it is present in the real exported
config `testData/Scripting/SizeDis.json`.

Nothing on the replay path reads it. `core/tool_config.py:sizes_model_from_config`
does not mention it, so neither `batch.fit_sizes` nor `analyze()` honours it —
both take their run count from their own `n_mc_runs` argument. A user who sets
20 uncertainty runs in the GUI, exports, and replays gets whatever the caller's
default is.

**Finished means:** `ToolSetup` carries it (it is session/run state, not model
state, so it belongs beside `fit_q_min`, not in `to_dict()`), and
`batch.fit_sizes` / `analyze` use it as the default when the caller does not
override. One field, three call sites. Overlaps item 1 — do it with that work.

---

## 5. `dq` is stored and never used — **S, or make it a decision**

`api/control/schemas.py:83` — *"Q resolution. Stored for provenance; not yet
used in fitting."* The sessions accept `dq`, the io layer round-trips it, and
no fit reads it. For a USAXS package that is a real capability (resolution
smearing beyond the slit-length model), not a loose end — but it has been
"stored for provenance" long enough that it should either be planned or
recorded in `PLAN.md` as decided-against. Right now it is neither, so every
caller has to discover it does nothing.

---

## 6. Stale "Phase 1" language in shipped surfaces — **S**

**Done 01-10-2026** (issue #31). `Q_indices` was resolved by writing it
as `[0]` rather than deleting the line.

Cheap, and worth doing before a beta because these strings are user-visible.

- `api/control/unified_fit.py:1987` and `schemas.py:699` — *"not available in
  Phase 1"*, *"A future phase will add ..."*. Users of the MCP/ZMQ surface see
  this; it means nothing to them. Say what is true: uncertainties are not
  computed for this tool yet. (Superseded if item 2 lands.)
- `gui/saxs_morph_panel.py:8-16` — Phase 2/3/4 narrative in the module
  docstring, now wrong (item 3).
- `gui/ai_advisor.py:9` — *"Phase 1 supports: Anthropic ... local
  OpenAI-compatible"*. Harmless but the same habit; if there is no Phase 2,
  drop the word.
- `io/hdf5.py:815` and `:873` — `# TODO not sure what this means` on a
  commented-out `Q_indices` attribute, duplicated. Either work out what
  NXcanSAS wants there or delete both lines.

Internal `planning/` documents should keep their phase language — that is
what they are for. This item is only about strings that ship.

---

## Not on this list, and why

| | Where it is recorded |
|---|---|
| Batch path for Data Manipulation | `PLAN.md` — wanted, tool needs revising first, own GitHub issue |
| MCP 2.x migration | `PLAN.md` — pinned `<2`, its own project |
| Large-file refactors (`unified_fit.py`, `modeling_panel.py`, …) | `PLAN.md` — no bug, no user impact |
| No batch/setup/control for Fractals, Merge, Manipulation | `test_tool_registration.py:ABSENT_REASON` — each with a written reason |
| No core `to_dict()` for SAXS Morph, Contrast, Merge, Manipulation, Fractals | `AGENTS.md` §2 — acknowledged, those tools serialise in io/gui |
| Cross-machine ZMQ run, deployment mechanics | `planning/zmq-service/README.md` §11 — has a written procedure and a checklist |

These are all *known and labelled*. The five above are the ones where the tree
and the documentation currently disagree.
