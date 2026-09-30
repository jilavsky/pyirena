# Config serialisation: one vocabulary, not four

Internal planning artifact. Written 25-09-2026 against pyIrena 1.1.1, prompted
by a question while building the ZMQ service: *why does Unified Fit — the
simplest tool to configure — have the most trouble being scripted?*

**Reviewed and re-measured 29-09-2026** on `feature/zmq-service` (Carbon model
included), against the real exported configs and result files in
`testData/Scripting/`. Every claim below that says *measured* was reproduced by
running the code, not read off the source. The review changed the picture in
three ways, and all three are in §2:

- the problem is **four** dialects and **three** envelopes, not two dialects;
- **Modeling has a second dialect too**, in HDF5 only, and it silently drops
  every non-size-distribution population's parameters;
- the **`export_results` → `analyze` round trip is lossy today**, for Unified
  Fit and for Sizes, in the ZMQ service as shipped.

**Status: done, 30-09-2026**, on `feature/config-dialects` — all five steps,
against issues #23-#27. Two decisions the plan left open were taken: Step 3
was done in full rather than deferred, and Fractals was given the core
`to_dict()`/`from_dict()` pair rather than left on the legacy reader (§6).

Three things came out differently from the plan, each noted at the point it
applies below:

- **Step 0 fix 1** was not the four-line patch to `flatten_level_config`'s
  bare-number branch that §5 describes. Routing a core-shaped level through
  the panel translator loses `K`, `mass_fractal` and `RgCO`'s flag and bounds
  by construction, because the panel vocabulary has no room for them. Levels
  are instead dispatched by shape — `unified_level_from_config` — which is
  Step 3.2 arriving early, and is what Step 3.2 then reused.
- **Modeling's `_collect_state` still writes the type blocks** (§5 Step 3.3
  said to stop). They are how switching a population's type in the GUI keeps
  the settings of the type switched away from; removing them to satisfy a
  consistency rule would delete a working feature. It now writes the active
  type's fields flat *as well*, which is what the headless readers and a
  human reading the attribute in HDFView actually need.
- **Sizes has no fit flags or bounds to invert**, so Step 2.3's round trip is
  vacuous for it; its two Step 0 losses are pinned by Step 2.1 (the slit
  fields were in `vars()` and not in `to_dict()`) and by the value round trip
  instead. Noted in the test's tool table so the gap is not mistaken for
  coverage.

---

## 1. The inventory

Two questions decide everything: **what shape is one tool's settings blob**
(the *dialect*), and **what wrapper is it written inside** (the *envelope*).
They are independent, and both have drifted.

### 1.1 Four dialects

| # | Dialect | Produced by | Consumed by |
|---|---|---|---|
| **1** | **Core** — `model.to_dict()` | `export_results()["config"]`; Modeling's *Export Parameters* | `from_dict()`; `build_setup()` for 4 of 6 tools |
| **2** | **Panel** — `<panel>._collect_state()` | every panel's StateManager save, every HDF5 `_pyirena_config`, and *Export Parameters* for 5 of 6 tools | `_apply_state()`; `build_setup()` |
| **3** | **Report** — `get_<tool>_config()` in `api/control/` | the MCP / ZMQ read surface | nothing — it is not replayable, by design |
| **4** | **Igor/Irena legacy** — panel key names frozen into files users already have | history | `UnifiedLevel.from_panel_params` |

Dialect 3 is fine and should stay: an agent asking "what is this model set to"
wants `{"name": "Rg", "value": 250, "fit": true}` rows, not a serialisation
format. It is listed because it is a third thing called "config" in the code
and it gets mistaken for the other two.

### 1.2 Three envelopes written, one accepted — and a fourth that is a bug

Measured, `pyirena/core/tool_config.py:build_setup`:

| Envelope | Written by | `build_setup()` accepts it? |
|---|---|---|
| `{"_pyirena_config": {...}, "<tool>": {...}}` | *Export Parameters* (JSON sidecar) | **yes** |
| `{"_pyirena_config": {...}, "state": {...}}` | `io/setup_config.py:write_setup_config` → **every HDF5 result file** | **no** — `ConfigError: Config has no '<tool>' section` |
| `{ok, tool, config, fit_q_range, data, ...}` | `api/control/export.py:export_results` → **every ZMQ/MCP reply** | **no** — same error |
| `{"tool": ..., "model": ...}` | nobody (the §4 proposal) | `section_for`'s docstring claims yes; **measured: it returns the envelope as if it were the section, and builds an all-defaults model with no error** |

So the two envelopes a remote caller actually holds — the one inside the
result file and the one the service just handed back — are the two that cannot
be fed back in. That is the whole of the ZMQ "config awkwardness"; it is not
about Unified Fit.

The last row is a live bug in `section_for`. The branch
`isinstance(config.get("model"), dict) and detected is None` is unreachable
whenever the envelope carries `"tool"` (because `detect_tool` then returns
non-None), and when it *is* reached it returns `config`, not `config["model"]`.
Measured:

```python
>>> build_setup({"model": {"num_levels": 3, "levels": [...], "background": 0.42}},
...             "unified_fit").model
# num_levels 1, background 0.0, one level at Rg=10.0 — every value silently defaulted
```

### 1.3 Per-tool dialect map

Measured against `testData/Scripting/` — one real exported config and one real
result file per tool.

| Tool | *Export Parameters* sidecar | HDF5 `_pyirena_config.state` | `export_results()["config"]` | Same? |
|---|---|---|---|---|
| Unified Fit | panel | panel | core | **no** |
| Size Distribution | panel | panel | core | **no** |
| Simple Fits | panel | panel | core | near — `param_limits`/`param_fixed` vs `limits` |
| WAXS Peak Fit | panel ≡ core + Q range | same | core | yes |
| **Modeling** | **core** | **panel (nested)** | core | **no — and the two disagree** |
| Carbon model | panel (`model` nested + view state) | same | core | yes, modulo the nesting |

Modeling is the one this review changed. §12 of the ZMQ plan recorded Modeling
as "one dialect ✔", verified against
`testData/Core-shell-tests/pyirena_config.json`. That verification was sound
for the channel it looked at: the file is an *Export Parameters* sidecar, it
does contain a `unified_level` population, and its fields are flat and correct
— because that path serialises with `_pop_to_dict`, the core dialect. What was
never compared is Modeling's **other** writer, `_collect_state()`, which is
what goes into HDF5 and into StateManager. The two disagree. See §2.2.

### 1.4 Four spellings of "the Q range I fitted"

| Tool | Key | Where it lives |
|---|---|---|
| Unified Fit | `cursor_left` / `cursor_right` | config section |
| Size Distribution | `cursor_q_min` / `cursor_q_max` | config section |
| Simple Fits, WAXS, SAXS Morph | `q_min` / `q_max` | config section |
| **Modeling, Carbon** | `q_min` / `q_max` | **inside `to_dict()`** |
| Fractals | `q_range` (a nested `{q_min, q_max, n_points}` object) | config section |

Note the last two rows. §4 of the original plan said the fitted Q range "is not
a model field and should not be forced into `to_dict()`" — but for Modeling and
Carbon it already *is* one, and that is why those two are the only tools whose
core config replays with the right Q range. The envelope design has to pick one
of the two and say so, rather than leave half the tools on each side.

### 1.5 Five encodings of "a parameter with a fit flag and bounds"

| Encoding | Used by |
|---|---|
| `X: 1000.0`, `fit_X: true`, `X_limits: [lo, hi]` | Unified **core**, Carbon core, Modeling `UnifiedLevelPopulation` |
| `X: {value, fit, low_limit, high_limit}` | Unified **panel**, Modeling's `uf`/`peak`/`gp` blocks |
| `X: {value, fit, lo, hi}` | WAXS peaks and background |
| `params: {}`, `limits: {}`, `param_fixed: {}` | Simple Fits |
| `dist_params: {}`, `dist_params_fit: {}`, `dist_params_limits: {}` | Modeling size-dist populations |

Five encodings for one idea, and `low_limit`/`high_limit` vs `lo`/`hi` is two
spellings of a single one of them. Simple Fits goes further and inverts the
sense of the flag: `param_fixed[name] = not checked`, where everyone else
stores `fit_X = checked`. An inverted boolean under a different name is
exactly the kind of thing that reads as correct in review and produces a wrong
fit.

---

## 2. What is actually broken

Ordered by consequence. Everything here was reproduced on
`feature/zmq-service`.

### 2.1 `export_results()` → `analyze()` silently drops Unified Fit's fit flags and bounds — **worst**

This is the documented ZMQ workflow: set the fit up once in the GUI, export it,
replay it on every new measurement. Measured:

```python
m = UnifiedFitModel(num_levels=1)
# user pins Rg and G, frees P, and narrows the bounds
lv.fit_Rg = lv.fit_G = False; lv.fit_P = True; lv.Rg_limits = (100.0, 400.0)

back = build_setup({"unified_fit": m.to_dict()}).model.levels[0]
```

| field | exported | replayed |
|---|---|---|
| `fit_Rg` | False | **True** |
| `fit_G` | False | **True** |
| `fit_P` | True | **False** |
| `fit_B` | False | **True** |
| `fit_ETA` | True | **False** |
| `fit_PACK` | True | **False** |
| `Rg_limits` | (100, 400) | **(0.1, 1e6)** |
| `G_limits` | (1, 9999) | **(1e-10, 1e10)** |

Every value is right; every *decision about* the values is lost. The replayed
fit floats parameters the scientist pinned and drops the bounds that kept it
physical. It does not error and it returns a plausible number.

The cause is in `tool_config.flatten_level_config`. It handles a bare number
(`if not isinstance(entry, dict): flat[name] = float(entry); continue`) so the
core dialect *looks* supported — but the `continue` skips the three lines that
would carry `fit_X`, `X_low` and `X_high`, and the core dialect's own
`fit_Rg` / `Rg_limits` keys are never read at all. It is a half-written
core-dialect reader, which is worse than none: `batch/unified.py:_state_to_model`
at least raised `AttributeError`.

### 2.2 Modeling's HDF5 setup drops every non-size-distribution population — **new**

`ModelingPanel._collect_state()` serialises populations with
`PopulationTab.to_full_dict()`, which keeps the size-distribution fields flat
and nests the other five types under `uf` / `peak` / `gp` / `mf` / `sf2` — so
that switching a population's type in the GUI does not lose the settings of the
type you switched away from. A good reason, and invisible until something other
than the panel reads it.

`core/modeling.py:population_from_dict()` reads the fields flat. Measured:

```python
panel_pop = {"pop_type": "unified_level", "enabled": True,
             "uf": {"Rg": 250.0, "G": 5000.0, "P": 3.2, ...}}
pop = population_from_dict(panel_pop)
# UnifiedLevelPopulation -> Rg = 10.0, G = 1.0, P = 4.0   (the dataclass defaults)
```

Confirmed present on disk in `testData/Scripting/Modeling_PP15.h5`
(`entry/modeling_results`, `pop[0]` carries all five nested blocks).

Reach: anything that reads a Modeling setup out of an HDF5 file and replays it —
`build_setup`, `batch.fit_modeling`, `analyze`. It does not reach *Load Setup
from File…*, because the panel reads its own nesting back correctly. So the GUI
looks fine and every headless path is wrong, for five of the six population
types. Today this is partly masked by §1.2: the HDF5 envelope is rejected before
the populations are ever read. Fixing the envelope without fixing this would
turn a clear error into a silent wrong answer — **they must be fixed together.**

### 2.3 Size Distribution loses the slit settings and the aspect ratio

Two separate losses in one tool.

**Slit smearing** — the original Problem A, still present. Measured:

```python
>>> m = SizesDistribution(); m.use_slit_smearing = True; m.slit_length = 0.0217
>>> SizesDistribution.from_dict(m.to_dict()).use_slit_smearing
False
```

`to_dict()` omits both fields. A slit-smeared setup, saved and reloaded, comes
back as a pinhole fit — a different, plausible answer. Unified Fit, Simple Fits
and Modeling all round-trip both fields, which is what makes this an
inconsistency rather than a design choice. The panel state omits them too, so
the HDF5 setup and the exported config lose them as well.

**Aspect ratio** — not previously recorded. The core dialect writes
`shape_params: {"aspect_ratio": 3.0}`; `tool_config.sizes_model_from_config`
reads the panel's `aspect_ratio` key first and only falls back to
`shape_params` when the shape is not `spheroid` — which is the one case where
the aspect ratio matters. Measured: a spheroid config at aspect ratio 3.0
replays at **1.0**, i.e. as a sphere. Same failure shape: no error, different
physics.

`montecarlo_n_repetitions` is also forced to 1 on replay. That one is
deliberate and commented; it is listed only so it is not re-reported.

### 2.4 `use_analytic_jacobian` is dropped by three tools

| Tool | has the attribute | in `to_dict()` | survives `from_dict()` |
|---|---|---|---|
| Unified Fit | yes | no | yes (constructor restores the default) |
| WAXS Peak Fit | yes | no | yes |
| Simple Fits | yes | **yes** — fixed 26-09-2026 | yes |
| Size Distribution | **no such attribute** | — | — |

Correction to the original table: Sizes never had this field. Unified and WAXS
lose a user's choice on round trip but do not crash, because their `from_dict`
uses the constructor. Simple Fits did crash, because its `from_dict` uses
`cls.__new__` and hand-assigns; that is fixed. `cls.__new__` plus
hand-assignment remains the pattern worth banning — Simple Fits is the only
tool still using it.

### 2.5 Smaller things, worth a line each

- **Sizes `unc_n_runs`** — a GUI control, persisted in the state and present in
  the real exported config (`testData/Scripting/SizeDis.json`), read by
  *nothing* on the replay path. The user's uncertainty-run count is silently
  ignored by both `batch.fit_sizes` and `analyze`.
- **Carbon `last_folder`** — an absolute path from the author's machine, written
  into the exported config and into every result file. Harmless but wrong in a
  document meant to be shared; view state does not belong in a portable config.
  (`testData/Scripting/carbon_CE_950.json` carries one.)
- **Simple Fits `param_limits` vs core `limits`** — one rename, already bridged
  in `simple_model_from_config`, listed for the mapping table.

### 2.6 Fixed, kept for the record

- **Simple Fits `from_dict()` produced an object that could not fit** (Problem
  C) — fixed 26-09-2026 on `feature/zmq-service`.
- **`save_to_nexus=False` wrote to the data file anyway** (Problem D) — fixed
  26-09-2026.

Both were the same shape as everything in §2: a field or parameter that exists
in one place and is quietly dropped in another. Four found in a week, in one
small area, none caught by the suite.

This family has been fixed before. A round of the same bug — Unified's
`link_B`, Simple Fits' `param_fixed`, Modeling's three newer population types —
was closed on main on 24-06-2026, and `core/tool_config.py` exists *because* of
it: one shared config→model path instead of one per caller. Sharing the builder
was necessary and has not been sufficient. §2.1 is the new builder dropping the
same class of flag, and doing it more quietly than the code it replaced. The
missing half is a test, not another refactor — §5 Step 2.3.

---

## 3. Why Unified Fit, and why Modeling too

Not complexity — age, plus a migration that was started and not finished.

`to_dict()` / `from_dict()` were added to `core/unified.py` and
`core/modeling.py` in the **same commit** (`dbffaf1`, "U6 and U10 fixes",
10-08-2026). That commit also rewrote the panels to use them —
`gui/modeling_panel.py` lost 231 lines, `gui/unified_fit.py` only 132, and
Unified kept four `from_panel_params` call sites. Modeling's panel came out
speaking the core dict for *Export Parameters*; Unified's did not.

Unified could not, and that was a reasonable call at the time. By then the panel
vocabulary was already:

- **inside every saved file** — `io/setup_config.py:write_setup_config` embeds
  the panel state as `_pyirena_config`, and *Load Setup from File…* reads it
  back;
- **in every user's state file** — `state/state_manager.py` stores the Unified
  defaults with `cursor_left`, `RgCutoff`, `correlated`, `estimate_B`;
- **a documented public format** — `docs/batch_api.md` §"Per-level fields"
  specifies it as *the* batch config, so users have hand-written files in it.

Changing it needed a compatibility bridge. The bridge was never built, so both
vocabularies stayed and `from_panel_params` became the seam.

What the 29-09 review adds: **Modeling only got half-way too.** Its *Export
Parameters* path speaks the core dialect, but `_collect_state()` — the path into
HDF5 and into StateManager — still speaks a panel dialect of its own, with the
type-block nesting. Modeling looked finished because the two writers were never
compared with each other. So this is not "Unified Fit is the odd one out"; it is
"one migration, six tools, four of them done".

---

## 4. What "synchronised" should mean

Four rules, in priority order. Rules 1–3 are unchanged in intent from the
original plan; rule 2 is sharpened by §1.2 and rule 4 is new.

1. **`to_dict()` / `from_dict()` is the one canonical dialect.** Complete —
   every setting that changes a number is in it — and the only dialect pyIrena
   *writes* from now on.

2. **One envelope, written everywhere, read everywhere.** Three writers and one
   reader is the actual defect. The envelope already almost exists:
   `export_results` emits `{tool, config, fit_q_range, data, ...}`, which is the
   §4.3 shape with `config` where the proposal said `model`. Adopt what is
   already written rather than inventing a fourth:

   ```json
   {"tool": "unified_fit",
    "config":      { ... to_dict() ... },
    "fit_q_range": {"q_min": 0.001, "q_max": 0.1},
    "data":        {"is_slit_smeared": false, "slit_length": 0.0}}
   ```

   `build_setup()` must accept this, the HDF5 `{_pyirena_config, state}` form,
   and the sidecar `{_pyirena_config, <tool>}` form. Three readers, one internal
   representation, one writer going forward.

3. **The panel dialect becomes a legacy input format.** Read forever, never
   written again. A user's existing config file and every saved `.h5` keep
   working; that is non-negotiable. Forward compatibility (an *older* pyIrena
   reading a *new* file) has never been promised and is not in scope.

4. **Decide where the fitted Q range lives, once.** Modeling and Carbon keep it
   in `to_dict()`; the other four keep it beside. Both work; having both is what
   costs. **Recommendation: keep it in the envelope's `fit_q_range`, and have
   Modeling and Carbon continue to carry it in `to_dict()` as well** — removing
   it from their `to_dict()` would break files and gains nothing, and the
   envelope is the single place a reader has to look. `build_setup` prefers the
   envelope and falls back to the model dict. Write this down in
   `docs/batch_api.md`; it is the rule someone will otherwise re-derive wrongly.

### 4.1 The HDF5 files have to stay readable by hand

A constraint, not a nice-to-have: users open result files in HDFView, in Igor,
and with `h5dump`, and read `_pyirena_config` as text. That rules out two
tempting simplifications and adds one requirement.

- **Keep it one JSON string attribute on the results group.** Not a binary blob,
  not a sub-group per control. It is already this and should stay.
- **Keep the key names the ones in the GUI and in the maths.** `Rg`, `G`, `P`,
  `RgCO`, `phi` — a scientist reading the attribute should recognise the Irena
  names. This is the argument *against* renaming core fields to the panel's
  names, not just for it: `RgCO`, `correlations` and `link_B` are the names in
  the docstrings and in Irena.
- **New requirement: the embedded config must be self-describing.** Add a
  one-line `"_note"` to the envelope header saying what it is and how to replay
  it, e.g. `"pyIrena setup; replay with pyirena.batch or analyze()"`, plus the
  `tool` and `pyirena_version` that are already there. Costs one key; saves
  every user who finds the attribute and has to guess.
- Adopting the core dialect *improves* hand-readability for Unified Fit:
  `"G": 31300.0, "fit_G": true, "G_limits": [6260, 156000]` reads at least as
  well as the nested form, and matches what the other five tools already show.

---

## 5. Plan

Five steps. Step 0 is new and is the only one that should go near the beta;
1–4 are cleanup and can wait. **S** ≈ one session, **M** ≈ two to three.

### Step 0 — The four silent-wrong-answer bugs (S–M) — *do this first, own branch*

Not cleanup. Each one returns a plausible wrong number rather than an error,
and three of them are on the ZMQ path as shipped.

1. **Unified core-dialect replay loses fit flags and bounds (§2.1).** In
   `flatten_level_config`, read the core dialect's own `fit_X` and `X_limits`
   keys in the bare-number branch instead of `continue`-ing past them.
2. **Modeling's nested population blocks (§2.2).** Teach
   `population_from_dict` to look inside `uf`/`peak`/`gp`/`mf`/`sf2` for the
   matching `pop_type` before falling back to the flat keys. One function; the
   panel keeps writing what it writes.
3. **Sizes aspect ratio (§2.3).** Read `shape_params` when `aspect_ratio` is
   absent, for every shape.
4. **Sizes slit smearing (§2.3).** Add `use_slit_smearing` and `slit_length` to
   `SizesDistribution.to_dict()` and read them in `from_dict()`; add them to the
   panel's `_collect_state()`. Additive, so old files still load.

Each needs a round-trip case in `test_core_serialization.py` asserting values
*and* flags *and* bounds. Fix 2 must land before or with Step 1, per §2.2.

### Step 1 — One envelope (S)

Make `build_setup()` accept all three envelopes of §1.2, and fix the
unreachable-and-wrong `{tool, model}` branch in `section_for` while there.
This is the step that makes "read the config out of a result file and re-run
it" work at all — the thing the ZMQ service most obviously ought to do.

Then: `export_results` keeps its shape, *Export Parameters* keeps writing the
sidecar, and the HDF5 attribute keeps its `{_pyirena_config, state}` wrapper.
Nothing users hold changes. Only the reader widens.

### Step 2 — A completeness contract (S)

Three tests, in the spirit of `test_tool_registration.py` and
`test_gui_state_contract.py`. A loop over all six tools, not one test per tool.

1. Every scalar/bool model attribute is either in `to_dict()` or named in an
   explicit `_NOT_SERIALISED` allowlist on the class, with a comment saying why.
2. `from_dict(to_dict())` produces an object that can **run a fit** — not merely
   one that compares equal. This is what would have caught Problem C.
3. **New, and the one that matters most:** for every tool, a model with *every
   fit flag inverted from its default and every bound narrowed* survives
   `to_dict() → build_setup() → to_dict()` unchanged. Every bug in §2.1–2.3
   fails this test and passes a values-only comparison. Drive it from the real
   fixtures in `testData/Scripting/`, which cover all six tools.

### Step 3 — Retire the panel dialect as an output (M)

The actual migration, and the only step that touches the GUI.

1. Add `UnifiedLevel.to_panel_params()`, the inverse of `from_panel_params`, and
   make `api/control/unified_fit.py:_session_to_gui_state` use it instead of
   hand-building the dict. One translator, both directions, one place.
2. Make every *reader* accept both dialects, detected by shape — a parameter
   that is a dict rather than a number is the panel dialect. That is
   `batch/unified.py:_state_to_model`, the setup loader, and `analyze`. Route
   through `from_panel_params`; do not write a second translator.
3. Switch the Unified panel's `_collect_state()` to the core dialect, and the
   StateManager defaults with it. Do the same for Modeling's population blocks
   (§2.2) — write the core dialect, keep reading the nesting.
4. Drop `last_folder` and the other view-only keys from what *Export Parameters*
   writes (§2.5). Keep them in StateManager, where they belong.

### Step 4 — Documentation and the mapping table (S)

Update `docs/batch_api.md`: the core dialect is current, the panel dialect is
accepted-legacy, and the mapping tables of §1.4 and §1.5 belong there rather
than only here. Add the §4.1 note about hand-readable HDF5 configs to
`docs/HDF5_NxcanSAS_structure.md`.

**Do not** rename the core fields to the panel's names instead.
`RgCO`/`correlations`/`link_B` are the names in the maths, in the docstrings and
in Irena; the panel names are the historical accident. §4.1 is the user-facing
argument for the same conclusion.

---

## 6. Risks

- **Numerical regression is the only unacceptable outcome.** None of this
  changes maths, but Steps 0 and 3 change what gets restored from a setup, and a
  setup restored wrong produces a different fit. Every step needs a before/after
  on `validationData/` (`run_validation_report.py`) showing the numbers
  unchanged, not just a green suite.
- **Step 0 fix 2 changes results for existing Modeling setups — correctly.** A
  user replaying a saved non-size-dist Modeling setup has been getting default
  parameters; after the fix they get theirs. That is the point, and it needs a
  release note, because someone's pipeline output will move.
- **Step 3 touches `gui/unified_fit.py`**, the largest file in the package
  (4,442 lines) and the one users exercise most. Its own branch, its own release
  note.
- **Fractals reuses Unified levels** (`core/fractals.py`, `gui/fractals_panel.py`
  both speak the panel dialect, and Fractals has no core `to_dict()` at all) —
  migrate it in Step 3 or leave it explicitly on the legacy reader.
- Cost of doing nothing: §2.1–2.3 stay, and they are wrong answers rather than
  errors. That is the part that is not "cleanup, not a fire".

---

## 7. What this is not

Not a rewrite of the state system, not a change to NXcanSAS layout, not a change
to any result group, and not a change to the `get_<tool>_config()` report
dialect. `to_dict()` already exists and already works for four tools; this is
finishing a migration that is two-thirds done — plus four bugs that the
migration's unfinished half has been hiding.
