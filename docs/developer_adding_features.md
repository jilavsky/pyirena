# Developer Guide: Adding Features to pyIrena

How to add a new analysis feature (a model, a calculation method, or a whole
tool) so it works everywhere users expect: GUI, JSON scripting, batch/CLI,
HDF5 persistence, MCP/AI tools, Data Selector (graphing/report/tabulation),
HDF5 Data Explorer (trend plots, Igor export), and documentation.

pyIrena has grown a number of integration surfaces; forgetting one produces a
feature that "works in the GUI" but silently vanishes from reports, batch
runs, or Igor exports.  This guide is the map.  It was written while adding
the **Invariant** method to Simple Fits (July 2026) — that commit is a good
reference diff for the full wiring of a new Simple Fits entry.

---

## Contents

1. [Architecture in one paragraph](#architecture-in-one-paragraph)
2. [Decision: what kind of feature is it?](#decision-what-kind-of-feature-is-it)
3. [The master checklist](#the-master-checklist)
4. [Wiring points in detail](#wiring-points-in-detail)
5. [Case study: the Invariant](#case-study-the-invariant)
6. [Conventions](#conventions)
7. [Testing requirements](#testing-requirements)

---

## Architecture in one paragraph

Every tool follows the same layering: a **core** module
(`pyirena/core/<tool>.py`) holds all math and a serialisable model object
(`to_dict`/`from_dict`); a **GUI panel** (`pyirena/gui/<tool>_panel.py`) is a
thin Qt shell over the core object whose full control state round-trips
through `StateManager` (`_collect_state`/`load_state`); an **io** module
(`pyirena/io/nxcansas_<tool>.py`) saves/loads results into the NXcanSAS HDF5
file under `entry/<tool>_results` and embeds the GUI state as
`_pyirena_config` for setup restore; a **batch** function
(`pyirena/batch/<tool>.py`) runs the same core headlessly from a dict or a
JSON config section; the **api/mcp** layer (`pyirena/api/results.py`,
`pyirena/mcp/server.py`) exposes read access to the AI tools; and the two
browsers — **Data Selector** (`gui/data_selector/`) and **HDF5 Data
Explorer** (`gui/hdf5viewer/`) — consume the saved HDF5 groups for graphing,
reports, CSV tabulation, trend plots, and Igor-experiment export
(`io/h5xp_extractor.py` + `io/igor_names.py`).

The golden rule: **all state lives in the core model object and flows dict →
JSON → HDF5 unchanged.**  If you add a control that only exists as a Qt
widget, it will not survive scripting, batch, or setup restore.

## Decision: what kind of feature is it?

| Feature type | Effort | Start at section |
|---|---|---|
| **New Simple Fits model** (formula + params) | Small — the registry does most wiring | [4.1](#41-core) |
| **Calculation method** (no least-squares, e.g. Invariant) | Medium — registry + `calculation: True` flag + routing | [Case study](#case-study-the-invariant) |
| **New option/control on an existing tool** | Small–medium — model attr + state + JSON + HDF5 attr | [4.2](#42-gui-panel) |
| **New form/structure factor** | See `developer_adding_form_factors.md` / `developer_adding_structure_factors.md` | — |
| **Entirely new tool** | Large — every row of the checklist, new files throughout | all |

## The master checklist

Work top to bottom; check each box.  Items marked (auto) come for free if you
follow the registry/dict conventions — verify, don't implement.

```
CORE
[ ] 1. Math in pyirena/core/<tool>.py — pure numpy/scipy, no Qt imports
[ ] 2. Registered: MODEL_REGISTRY entry (Simple Fits) or model-object attribute
[ ] 3. New state serialised in to_dict()/from_dict() (defaults for old files!)
[ ] 4. Exported from pyirena/__init__.py if part of the public API

GUI
[ ] 5. Controls in the panel; visibility handled on model/option change
[ ] 6. State: _collect_state() + load_state() include every new key
       (this feeds StateManager, JSON export AND the HDF5-embedded setup)
[ ] 7. Plots updated + cleared on model change (no stale curves)
[ ] 8. Help button / tooltips updated

PERSISTENCE (HDF5, NXcanSAS)
[ ] 9. io/nxcansas_<tool>.py: save new scalars (params/, derived/, attrs)
       and arrays (datasets with units attr); loader reads them back
[ ] 10. Setup round-trip: "Store in File" → "Load Setup from File…" restores
        every new control (uses the state dict from step 6 — verify only)
[ ] 11. io/schema.py: add _scalar()/_plot() entries so generic consumers
        (HDF5 viewer scalar browser, merge tool) know about the values

SCRIPTING / BATCH / CLI
[ ] 12. JSON: "Save/Load params to JSON" round-trips (auto via step 6)
[ ] 13. batch/<tool>.py: headless function honours the new option
        (from_dict conventions make this automatic; check special routing)
[ ] 14. fit_pyirena / pipeline.py picks it up from the config section (auto)

MCP / AI API
[ ] 15. api/results.py read_<tool>: new fields returned (auto if generic
        params/derived loops; add explicit fields otherwise)
[ ] 16. api/aggregate.py: parameter reachable by tabulate/trend tools
        (fallback probes params/<name> and derived/<name> — auto)
[ ] 17. mcp/server.py: only needs changes for genuinely new tool functions

DATA SELECTOR (Data Browser)
[ ] 18. Create Report (data_selector/report.py): new values in the Markdown
        report (auto if stored under params/ or derived/)
[ ] 19. Create Graph (data_selector/results_windows.py): stored results
        render sensibly (special-case if the tool has no model curve)
[ ] 20. Tabulate Results: columns appear (auto via params/derived)
[ ] 21. <Tool> (script) button path works (auto via step 13)

HDF5 DATA EXPLORER
[ ] 22. plot_controls.py: new item names in the "collect" dropdown for
        trend plots (parameter vs file/temperature/time)
[ ] 23. pyirena_readers.py _collect_<tool>: new values readable
        (params/ auto; derived/ needs the fallback added July 2026)
[ ] 24. Igor export: io/igor_names.py wave-name entries (Y wave and
        RESULT_X_WAVE pairing) + io/h5xp_extractor.py writes the waves
        and wave-note parameters

DOCS & TESTS
[ ] 25. docs/<tool>_gui.md: user documentation + scripting example
[ ] 26. CHANGELOG.md entry
[ ] 27. Tests: math vs analytic ground truth; registry/serialisation
        round-trip; HDF5 save/load round-trip; batch-config path
[ ] 28. Run the full test suite, not just the new file
```

## Wiring points in detail

### 4.1 Core

- **Simple Fits models** live in `pyirena/core/simple_fits.py`.  Add a private
  formula function `_my_model(q, ParamA, ParamB) -> np.ndarray` and one
  `MODEL_REGISTRY` entry:

  ```python
  'My Model': {
      'params': [('ParamA', 1.0, 1e-30, None),    # (name, default, lo, hi)
                 ('ParamB', 50.0, 0.1, 10_000.0)],
      'formula': _my_model,
      'linearization': None,     # or 'guinier' / 'porod' family key
      'complex_bg': True,        # allow B·Q⁻ᴾ + flat to be added
      # 'calculation': True,     # only for no-fit methods (see case study)
  }
  ```

  The registry automatically propagates the model into: the GUI model combo,
  parameter grid, state persistence, JSON configs, HDF5 params/, batch
  fitting, MC uncertainty, the report and tabulation tools, and the existing
  parametrised tests in `tests/test_simple_fits_models.py` (add your model to
  the IDENTIFIABLE / DEGENERATE / CALCULATION list there).

- Derived quantities (computed from fitted params) go in
  `SimpleFitModel._compute_derived()`; they are saved under `derived/` and
  surface everywhere automatically.

- **Boolean/enum options** are attributes on the model object (like
  `use_complex_bg`, `invariant_porod_tail`).  They MUST be added to
  `to_dict()` and `from_dict()` with a default that keeps old JSON/HDF5 files
  loading.

### 4.2 GUI panel

`gui/simple_fits_panel.py` (other tools are analogous):

- Add widgets in `_create_control_panel()`; wrap option groups that only
  apply to one model in a container `QWidget` and toggle its visibility from
  a `_update_<feature>_ui()` helper called in **both** `_on_model_changed()`
  and `load_state()` (the two entry points into a model switch).
- Every new control needs three touch points: a change-handler syncing widget
  → model attribute, a line in `_collect_state()`, and a restore (with
  `blockSignals(True)`) in `load_state()`.  `_collect_state()` is used for
  StateManager persistence, JSON export, **and** the setup embedded in HDF5 —
  one dict, three consumers.
- Clear your plot items in `clear_fit()`/model change so nothing stale
  survives a model switch.

### 4.3 HDF5 (NXcanSAS)

`io/nxcansas_simple_fits.py` layout under `entry/simple_fit_results`:
group attrs (model, success, options…), scalar datasets (`chi_squared`),
`params/`, `params_std/`, `derived/` sub-groups (one scalar dataset each),
array datasets with a `units` attr, `fit_quality/`, and the embedded
`_pyirena_config` setup JSON.  New scalars: prefer `params/` (inputs) and
`derived/` (results) — generic consumers pick those up.  New arrays: add the
dataset with units in the saver, read it back in the loader, and document it
in the module docstring.

`io/schema.py` is the generic description used by the HDF5 viewer and the
merge tool — add `_scalar()` entries (path, units, label) and `_plot()`
entries for new arrays.

### 4.4 Batch / JSON

`batch/simple.py::fit_simple()` builds the model with
`SimpleFitModel.from_dict(config)` — so a correctly serialised feature is
scriptable with zero batch changes.  Watch for: GUI-state-only keys that
`from_dict` does not understand are stripped in `fit_simple_from_config()`
(`q_min`/`q_max`, `param_fixed`, `no_limits`, `schema_version`) — if you add
a state key that is *not* a model attribute, strip it there too.

### 4.5 MCP / AI API

`api/results.py` readers return `params`/`params_std`/`derived` dicts
verbatim — new values flow through to `pyirena_read_simple_fit` MCP tool
automatically.  `api/aggregate.py::_fallback_lookup()` resolves parameter
names against `params/<name>`, `derived/<name>`, top-level datasets and
attrs, so `pyirena_tabulate_parameter` / `pyirena_plot_parameter_trend` reach
new values without schema entries (schema entries still add units/labels).

### 4.6 Data Selector

- `report.py`: iterates `params` and `derived` generically — verify only.
- `results_windows.py`: plots `Q`/`I_model` from the stored group.  If your
  feature has no meaningful model curve, add a special case (see the
  Invariant block there).
- Tabulate/CSV: generic over params/derived.

### 4.7 HDF5 Data Explorer + Igor export

- `gui/hdf5viewer/plot_controls.py`: the per-tool "collect item" dropdown
  lists are **hardcoded** — add your new parameter/derived names.
- `gui/hdf5viewer/pyirena_readers.py::_collect_simple_fit()`: reads
  `params/` then falls back to `derived/`.
- Igor export (`io/h5xp_extractor.py` + `io/igor_names.py`): every exported
  curve needs a Y-wave name and an X-wave pairing.  Simple Fits model waves
  follow `SimFit<Name>I` → `SimFit<Name>Q` (a regex handles pyirena-only
  names ending in `I`; other names need a `RESULT_X_WAVE` entry).  Scalars go
  into the Igor wave note via the `params` dict passed to
  `write_result_wave()` — include `derived/` values there.

## Case study: the Invariant

The Invariant (added July 2026) is the template for a **calculation method**
— a registry entry with `'calculation': True`:

- `core/simple_fits.py`: standalone `calculate_invariant()` (pure function,
  analytic units test against 2π²Δρ²φ(1−φ)), a placeholder formula returning
  zeros (so "Graph model" shows just the background), and
  `SimpleFitModel.fit()` routing to `_run_calculation()` which returns a
  fit-shaped result dict: `chi2`/`dof`/`residuals` are `None`, results in
  `derived`, running-integral arrays in `extra_arrays` (saved as datasets by
  the io layer), plus a `warning` string (saved as a group attr).
- Because the result dict keeps the standard shape, the HDF5 saver, batch,
  report, tabulation, MCP and Igor export all worked with only additive
  changes (extra arrays, derived fallback, dropdown items, wave names).
- GUI: `_update_invariant_ui()` relabels the Fit button, hides the "Fit?"
  column, disables MC uncertainty, and shows the Porod-tail checkbox; the
  running integral is drawn on a second ViewBox tied to the right axis
  (log data must be log10-ed manually there — pyqtgraph does not transform
  items in a bare ViewBox).
- Guards added for calculation models: skip fit-quality metrics in the saver
  (`chi2 is None`), skip MC uncertainty in batch.

## Conventions

- Parameter names are case-sensitive everywhere; pick once, never rename
  (they are HDF5 dataset names and JSON keys).
- Units: Q in Å⁻¹, I in cm⁻¹ (absolute), contrast Δρ² in 10²⁰ cm⁻⁴,
  invariant in cm⁻⁴.  Always set the `units` attr on new datasets.
- Igor parity: when porting from `IR3_*.ipf`, replicate the algorithm
  (cite file + function in the docstring), then add improvements as
  **off-by-default options** so numbers match Igor unless the user opts in.
- No Qt imports in `core/` or `io/`.  No `print()` — use `logging`.
- Old-file compatibility: every new dict key and HDF5 attr needs a default
  on read.

## Testing requirements

Minimum for a new model/method (see `tests/test_invariant.py` as template):

1. **Analytic ground truth** — synthesise data where the answer is known in
   closed form and assert recovery within tolerance.
2. **Registry/serialisation** — `to_dict`/`from_dict` round-trip including
   new attributes; add the model to the categorised lists in
   `test_simple_fits_models.py`.
3. **HDF5 round-trip** — save → load → compare scalars and arrays.
4. **Batch config path** — `fit_simple(file, config_dict)` end-to-end on a
   temp NXcanSAS file, including re-loading the stored result.
5. Run the **whole** suite: `python -m pytest pyirena/tests` — registry-
   parametrised tests in other files will pick up your model.
