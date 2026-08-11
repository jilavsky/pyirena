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

STANDARD UX CONTRACT (see below — every new GUI surface)
[ ] 25. Every table: attach_table_copy() (+ enable_table_sorting() where the
        row order is not meaningful); CSV export for multi-row results
[ ] 26. Every plot: attach_plot_export() — or make_sas_plot(parent_widget=…)
[ ] 27. Every file list: shared filter (gui/file_filter.py), shared sort modes,
        remembered folder

DOCS & TESTS
[ ] 28. docs/<tool>_gui.md: user documentation + scripting example
[ ] 29. CHANGELOG.md entry
[ ] 30. Tests: math vs analytic ground truth; registry/serialisation
        round-trip; HDF5 save/load round-trip; batch-config path
[ ] 31. Run the full test suite, not just the new file
```

## Serialising a tool's state

The core model owns the state; everything else asks it for a dict.  Two
methods, whatever the tool:

```python
def to_dict(self) -> dict:      # settings, not results — no data arrays, no χ²
def from_dict(cls, d: dict):    # a default for *every* field
```

`from_dict` supplying a default for every field is not politeness, it is the
backwards-compatibility guarantee: a file written before a field existed has to
open with that field at its default.  `pyirena/tests/test_core_serialization.py`
tests both directions for every tool that has the pair, through
`json.dumps`/`loads` — which is where tuples silently become lists.  Bounds are
tuples in the model and lists in the file; restore them, because the fitter
unpacks them.

For a dataclass-based tool, inherit `_SerialisableDataclass`
(`core/modeling.py`) rather than writing the two methods out: it walks
`dataclasses.fields()`, so a new field is serialised the moment it is declared.
Modeling's six population types share one implementation this way, and
`population_from_dict` dispatches on `pop_type` through the
`POPULATION_CLASSES` registry — an unknown type from a newer pyIrena loads as a
size distribution with a warning instead of raising.

**The panel's key names are not the model's.**  `RgCutoff`, `correlated`,
`estimate_B`, `Rg_low` are baked into saved setups and batch config files and
cannot be renamed.  Translate in exactly one place —
`UnifiedLevel.from_panel_params()` is the pattern — and never a second time in
a panel or a batch module.

## Qt imports

One import point, no exceptions:

```python
from pyirena.gui._qt import QLabel, QTimer, Qt, Signal
```

If a class is missing, add it to `pyirena/gui/_qt.py` (and to its `__all__`) —
do **not** write a local `try: from PySide6 … except ImportError: …` block, not
even inside a function, and not in a test.  Those blocks are how the project
ended up with three shims and a stale PyQt5 branch that would have failed on a
PyQt6-only install.  `pyirena/tests/test_gui_qt_contract.py` enforces this: it
fails on any direct binding import in the package, on a second `_qt.py`
appearing in a subpackage, on a name imported from `_qt.py` that it does not
define, and on a name defined there but missing from `__all__`.

Lazy imports are still fine when you want them — `from pyirena.gui._qt import
QMenu` inside a method is one line and costs nothing after the first call.

## The standard UX contract

Users come from Igor Pro, Origin, Excel and Matlab.  Features they consider
part of "a table" or "a graph" are not optional extras, and when they live in
nobody's spec they get skipped — twice now (the filename filter, then table
clipboard copy).  Shared helpers exist so each line below is one function call.

**Every table** (`QTableWidget`)

```python
from pyirena.gui.table_utils import attach_table_copy, enable_table_sorting

attach_table_copy(my_table, on_save_csv=self._save_csv)   # Ctrl+C, Ctrl+Shift+C,
                                                          # right-click menu
enable_table_sorting(my_table)      # only if row order carries no meaning
```

Build numeric cells with `make_numeric_item()` so columns sort numerically, and
fill the table inside `with populating(my_table):` so an active sort cannot
scramble half-built rows.  Write CSV with `rows_to_csv_text()` /
`save_rows_as_csv()` rather than `",".join(...)` — quoting and precision are
then identical everywhere.  `pyirena/tests/test_gui_table_contract.py` fails
the build if a new module constructs a table without the helper.

**Every plot** (`pg.PlotItem`)

```python
from pyirena.gui.plot_export import attach_plot_export

attach_plot_export(my_plot, self, 'my_tool_iq',
                   window=self.graphics_layout,   # whole-window image
                   folder=self._export_folder)    # dialogs open next to data
```

Plots built with `make_sas_plot(..., parent_widget=...)` get this for free.
The menu is clipboard copy, image (PNG/JPEG/SVG), whole-window image, curve
CSV and Igor ITX — never add your own image-export action, and never import
`ImageExporter` outside `plot_export.py`; save dialogs must default to
`export_folder()`, not `Path.home()`.  `pyirena/tests/test_gui_plot_contract.py`
enforces all of this.

Two things the exporters cannot infer, so the panel must say them:

- **Name every curve** (`name=...`).  Unnamed curves are exported only as a
  fallback, labelled from the Y axis, and a plot mixing named and unnamed
  curves exports the named ones only — an unnamed residual trace next to a
  named one would silently vanish.
- **Call `tag_curve_uncertainty(item, dI)`** wherever you draw your own error
  bars.  Bars are NaN-separated line segments that no exporter can read back,
  so without this the `dY` column disappears from CSV and ITX while the bars
  stay visible on screen.  `plot_iq_data` already does it for you.

**Every file list**

```python
from pyirena.core.file_sorting import SORT_LABELS, SORT_TOOLTIP, sort_names
from pyirena.core.file_types import FILE_TYPES, files_in_folder
from pyirena.gui.file_filter import FILTER_PLACEHOLDER, FILTER_TOOLTIP

self.type_combo.addItems(FILE_TYPES)
self.sort_combo.addItems(SORT_LABELS)          # never a local label list
self.sort_combo.setToolTip(SORT_TOOLTIP)
...
files = files_in_folder(self.current_folder, self.type_combo.currentText())
files = sort_names(files, self.sort_combo.currentIndex())
```

plus the shared filter box and a remembered last folder via `StateManager`.
Never re-implement a `_sort_key_*` function, never hard-code an extension list,
and never write your own `os.listdir` loop — `files_in_folder()` (or
`files_with_extensions()`, when your dropdown is not `FILE_TYPES`) already
decides how sub-directories, letter case and an unreadable folder are handled,
and the four browsers disagreed about all three before it was shared.
`pyirena/tests/test_gui_browser_contract.py` fails the build on each of these.

There is deliberately **no** shared `FileBrowserWidget`. The logic is shared —
filtering, sorting, the file-type table, listing, drag-and-drop — while the
widget assembly is not, because the four browsers differ in ways a common
widget would have to be configured around anyway: the Data Explorer is a tree
with lazy sub-folder expansion, Data Merge shows two linked instances, Data
Manipulation adds a context menu, and the Data Selector alone offers text files
and the convert-on-load path. Duplicated assembly with no logic in it does not
drift; duplicated logic does.

A file list should also **accept dropped files**, in one call:

```python
from pyirena.gui.file_drop import drop_hint, enable_file_drop, first_folder

enable_file_drop(self, self.open_dropped_paths)      # whole panel, not just the list
self.file_list.setToolTip(drop_hint("a folder or data files"))
```

`enable_file_drop` installs an event filter (no subclassing), rejects the drag
while it is still over the window if nothing dropped is openable, expands a
dropped folder one level, and hands you absolute paths already de-duplicated
and sorted.  Your handler switches folder with `first_folder(paths)` and selects
the rest; a single-dataset target uses `paths[0]`.  Never parse the mime data
yourself — Finder sends URLs, some apps send text, and `paths_from_mime()`
already reads both.

**Every top-level window** — remember where it was:

```python
from pyirena.gui.window_state import install_window_state

install_window_state(self, 'my_tool', splitters={'main': self.main_splitter})
```

and one line at the launch site, so the tool answers Shift-click:

```python
reset_window_if_shift(self.my_tool_window)   # before .show()
self.my_tool_window.show()
```

One call restores size, position and splitter sizes, and saves them on close.
Do not call `resize()`/`move()` from saved numbers yourself: `window_state`
owns the policy for a screen layout that has changed since the geometry was
saved (see `resolve_geometry`), and that is the part that goes wrong.  Geometry
lives in its own `window_geometry.json`, *not* in the tool's `StateManager`
section — panels each hold their own `StateManager` and a `save()` rewrites the
whole file, so geometry written there is clobbered by the next panel to close.
Register only real windows.  A widget that ends up inside a splitter is a
pane, and setting a pane's geometry fights the layout (the symptom is the
*neighbouring* pane changing width); `window_state` checks `isWindow()` on
every apply, so a mistake is a no-op rather than a puzzle.  Pane widths are
recorded and restored at the first layout, because a splitter has no width in
the constructor — do not expect `install_window_state` to have touched the
splitter by the time it returns.

Panels are constructed once and reused, so the reset gesture cannot live in
the constructor: `install_window_state` records the window's coded default
before applying anything saved, and `reset_window()` restores *that*.  Adding a
tool without the launch-site line is caught by
`pyirena/tests/test_window_state.py`, which reads the Data Selector's
launchers.  Set `GEOMETRY_PERSISTENCE_ENABLED = False` in that module to switch
the whole feature off.

**Every fit panel** — results reachable as text, not only as HDF5:

```python
from pyirena.gui.report_buttons import make_report_buttons

layout.addLayout(make_report_buttons(
    self, self.results_for_report, tool_key='my_tool_results',
    default_stem='my_tool', status_setter=lambda m: self.status_label.setText(m),
    image_widget_provider=lambda: self.graph_window.graphics_layout,
))
```

`image_widget_provider` enables Ctrl/⌘-click → report plus embedded graph
figure; return the `GraphicsLayoutWidget` (the plots alone), not the window,
or the figure will contain tab bars and status lines.

`results_for_report()` returns the dict shape
`pyirena.io.load_<tool>_results()` produces — the *saved* key names, not the
fit object's internal ones — and `pyirena.core.reporting` renders it.  Add the
tool's section there rather than formatting text in the panel, so the panel,
the Data Selector report and the MCP `export_fit_report` cannot drift apart.

**Every panel** — state methods use one naming convention, so a helper can call
them and a reader of one panel knows the next:

| Method | Visibility | Does |
|---|---|---|
| `save_state()` | public | collect + write to `StateManager` (+ `save()`) |
| `load_state()` | public | read from `StateManager` + apply |
| `_collect_state()` | private | return the state dict |
| `_apply_state(state)` | private | apply a state dict to the widgets |

The private halves are optional for small panels, the public pair is not —
`pyirena/tests/test_gui_state_contract.py` fails the build if a panel is
missing it or brings back one of the retired names (`_save_state`,
`_load_state`, `_restore_state`, `_get_current_state`).  An **embedded
component** whose parent owns the persistence lifecycle (the Diffraction Lines
tab inside WAXS) instead exposes `collect_state()` / `apply_state()` publicly,
because the parent is the caller.

`UnifiedFitPanel` additionally keeps `get_current_state()` / `apply_state()` as
public aliases: those names are quoted in `pyirena/io/setup_config.py` and the
api/control layer as the shape of the embedded `_pyirena_config` setup state,
and agent scripts call them.

Also on every panel: params-to-JSON, and tooltips on buttons.

**Params-to-JSON — one deliberate split.**  Analysis tools append their section
to a shared `pyirena_config.json` that drives the batch pipeline.  **Data Merge**
and **Data Manipulation** write their *own* config file instead, and that is
intentional, not an oversight to be tidied away: both produce **new data files**
that a later pipeline stage consumes, so their config has to live next to that
stage, in a different folder and a different pipeline step from the analysis
config.  Unifying them would break instrument data-reduction pipelines.
**Fractals** is a visualization tool with no batch use case and no JSON by
design.  **Scattering Contrast** persists compounds through its own HDF5 library
(`io/contrast_io.py`).

**Every result** — reachable via HDF5 *and* `pyirena.api`.

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
- **Prefit replay** (a lesson learned mid-feature): the invariant depends
  critically on background subtraction, and exported B/P/flat values rarely
  transfer between files.  The fix was to make the GUI's
  "Fit … btwn cursors" buttons *record their Q windows* into a
  `bg_prefit` dict on the model (serialised like everything else), plus a
  `prefit_background(q, I)` method that both the GUI Calculate and
  `batch.fit_simple()` invoke with the FULL data before calculating.
  General principle: **if a scripted run cannot reproduce an interactive
  step, persist the step's inputs (ranges), not just its outputs
  (values).**

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
