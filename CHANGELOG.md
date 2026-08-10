# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Modeling is now agent-drivable** (feature parity review item U9, issue #13).
  `pyirena/api/control/modeling.py` adds 18 tools — population management (add,
  remove, enable, list), per-parameter value/fit/bounds, the non-numeric options
  (distribution, form factor, structure factor, peak type, correlations), the
  background, the Q range, the fit, results, a per-population fit image and save
  to NXcanSAS — exposed over MCP as `pyirena_ctrl_modeling_*`.
  - **One flat parameter namespace.** A population's parameters live in plain
    attributes and in three nested dicts; the surface flattens them to dotted
    names (`dist.mean_size`, `ff.sld_core`, `sf.eta`, `scale`) and lists only
    the ones active for the current distribution, form factor and structure
    factor — the same set the fitter packs. Switching to a core-shell form
    factor adds its SLD and thickness parameters with sensible defaults and
    removes `contrast`, whose role the SLDs take over.
  - The fit image draws **each population separately** alongside the total, so
    a population that has collapsed to nothing is visible at a glance.
  - The MCP surface grows from 83 to 101 tools.
  - Only WAXS Peak Fit remains without a control surface.
- **Simple Fits is now agent-drivable** (feature parity review item U9, issue
  #13). Interactive control sessions existed for Unified Fit and Size
  Distribution only, so an AI agent could read Simple Fits results but never
  produce them. `pyirena/api/control/simple_fits.py` adds 15 tools — model
  listing and selection, per-parameter value/bounds/fix/free, the complex
  background, the fit, results, a fit image, a linearization image, and save to
  NXcanSAS — exposed over MCP as `pyirena_ctrl_simple_*`. The session and
  Q-range tools are shared with the existing surfaces, so a conversation can
  open a dataset once and try Unified Fit, Sizes and Simple Fits on it.
  Modeling and WAXS Peak Fit remain to be wired.
  - The linearization tool returns the slope, intercept and R² alongside the
    image, which gives an agent a numeric handle on "is this model valid over
    this Q range?" rather than only a picture.
  - The MCP surface grows from 68 to 83 tools; the schema/callable/MCP parity
    tests and their locked counts are updated with it.
- **Data Merge gained the full sort dropdown** the other three file browsers
  have (feature parity review item U3, issue #13). It sorted silently by order
  number with no way to reorder a temperature or time series. Both dataset
  columns now offer all ten modes — filename, temperature, time, order number,
  pressure, each ascending and descending — defaulting to order number ↑ so
  the view opens exactly as before, and each column remembers its own choice
  between sessions (`data_merge` state schema 2).

- **"Copy results" and "Save report…" on the fit panels** (feature parity review
  item A5, issue #13). Igor Irena wrote every fit to the notebook; in pyIrena a
  fit panel's only exits were save-to-HDF5 (then tabulate in the Data Selector)
  or reading numbers off the widgets one at a time — while the api/control layer
  had `export_fit_report` for AI agents and not for people. **Unified Fit, Size
  Distribution, Simple Fits, Modeling and WAXS Peak Fit** now each have two
  buttons that put the current results on the clipboard, or into a `.md` file,
  as Markdown: parameters with their Monte-Carlo uncertainties, fit quality
  (including the robust metrics), the data summary and the setup. WAXS reports
  the peak rows on screen, so it works before a fit as well as after one.
- **Ctrl-click (⌘-click on macOS) "Save report…" saves the report *and* the
  graph**: a PNG of the panel's plots is written next to the `.md` with the
  same name and embedded in it as a figure, replacing the save-report,
  save-graph, find-both, insert-the-image sequence. The link is relative, so
  moving the pair together keeps the figure working, and the result renders in
  GitHub, VS Code, Jupyter, Obsidian and pandoc — `pandoc report.md -o
  report.docx` yields a Word document with the graph in place. The figure is
  the plots alone (the graphics layout is captured, not the window), so no tab
  bar or status line appears in it. A plain click still writes text only, and a
  graph that cannot be captured costs the figure, not the report.
- `pyirena.core.modeling.result_to_report_dict()` flattens a live
  `ModelingResult` into the saved-results dict shape, so the Modeling panel can
  report the fit it is holding without a save-and-reload round-trip.
  Populations are flattened with `dataclasses.asdict`, so a new population type
  or parameter reaches the report with no change to the converter.
- The text comes from **one builder shared by every consumer**
  (`pyirena/core/reporting.py`), so a value cannot be formatted one way in the
  panel and another way in the report: the panel buttons, the Data Selector's
  *Create Report*, and the MCP `export_fit_report` all render the same sections
  from the same dict shape (`load_<tool>_results()`, which is also what
  `pyirena.api.results` returns).

- **Copy any graph to the clipboard, and save it as PNG, SVG or CSV** (feature
  parity review items U2 / A3 / A4 / A8, issue #13). Every plot's right-click
  menu now offers the same five actions, from one shared implementation in
  `pyirena/gui/plot_export.py`:
  - **Copy graph to clipboard** — the graph as an image, ready to paste into
    PowerPoint, Word, email or an electronic notebook. This is Igor's
    Edit→Copy, which had no equivalent anywhere in pyIrena.
  - **Save graph as image…** — PNG (now the default), JPEG or SVG, selected by
    the file filter or the extension typed. JPEG compression visibly smears the
    thin lines and small text of a log-log SAS plot; it stays available for
    users who need it.
  - **Save whole window as image…** — every stacked panel (data + residuals +
    distribution) in one image, alongside the single-plot export.
  - **Save curve data as CSV…** — the plotted curves as text for Excel, Origin,
    Matlab or pandas: one `X`/`Y` column pair per curve, plus `dY` where
    uncertainties exist, padded when curves have different lengths. Previously
    the only text export was Igor `.itx`.
  - Curves the panel never labelled are exported too, named from the Y axis.
    Residual panels, the Simple Fits linearization, the Contrast (Δρ)² plot and
    the collected-values scatter all plot without a legend name, and every
    export used to refuse with "No named data curves found to export" while the
    curve was plainly visible. `ScatterPlotItem` and bare `PlotCurveItem`
    curves are collected as well — the linearization plot is built from those
    and could not be exported at all. Error-bar segments are still never
    mistaken for data.
  - Uncertainties now reach the export from panels that draw their own data
    scatter (Unified Fit including the Porod tab, Sizes, WAXS). Only
    `plot_iq_data` recorded them before, so a CSV or ITX from those panels
    silently lacked the `dY` column although error bars were on screen. Panels
    record them with the new `tag_curve_uncertainty()`.
  - **Save as Igor Pro ITX…** — now available on every plot rather than most of
    them, and the file reproduces how the curve looks in pyIrena: scatter data
    imports as markers (`mode=3, marker=19`), model curves as lines. Igor draws
    every imported wave as a line by default, so exported data points used to
    arrive as a zig-zag line indistinguishable from a model curve.
    Uncertainties import as Igor error bars rather than as an extra trace.
- **Export dialogs remember where you last saved**, falling back to the panel's
  data folder before anything has been exported, instead of defaulting to the
  home directory (five dialogs) or the process working directory (the HDF5
  Viewer's). Where the user last chose to save takes priority over the data
  folder — a panel that reset to its data folder every time looks like the
  memory is broken. The folder persists across restarts under a new `exports`
  section of the state file.
- `pyirena/tests/test_gui_plot_contract.py` — a source-level guard that fails
  the build if a module drives pyqtgraph's exporters directly, adds its own
  JPEG/PNG menu entry, or points a save dialog at `Path.home()`. The plot half
  of the standard UX contract is now written down in
  `docs/developer_adding_features.md`.

### Changed

- **Panel state methods follow one naming convention** (feature parity review
  item U5, issue #13). The same two operations had grown six names —
  `save_state` / `_save_state` / `load_state` / `_load_state` /
  `_restore_state` / `_get_current_state` — so reading one panel taught you
  nothing about the next and no shared helper could call them. Every panel now
  exposes public `save_state()` / `load_state()` with private
  `_collect_state()` / `_apply_state(state)`; embedded components driven by a
  parent (the Diffraction Lines tab inside WAXS) expose `collect_state()` /
  `apply_state()`. Pure rename, no behaviour change. `UnifiedFitPanel` keeps
  `get_current_state()` / `apply_state()` as public aliases because those names
  are documented as the `_pyirena_config` setup-state shape and agent scripts
  call them. WAXS Peak Fit gained the `load_state()` it never had (it applied
  saved state inline in `__init__`).
- **The form-factor parameter registry moved to `pyirena/core/form_factors.py`**
  (`FORM_FACTOR_PARAMS`, `FORM_FACTOR_PARAM_DEFAULTS`,
  `CONTRAST_FREE_FORM_FACTORS`) from `gui/modeling_panel.py`. Which parameters
  a core-shell shape needs, and their defaults, sat next to the widgets that
  drew them, so the api/control layer — which may not import `gui` — could not
  tell an agent what to set. The GUI now reads the same registry, so a new form
  factor is declared once, beside its G-matrix builder.
- **The data-file type table is shared** in `pyirena/core/file_types.py`
  (`FILE_TYPES`, `FILE_TYPE_EXTS`, `files_in_folder()`) instead of being copied
  verbatim into Data Manipulation and Data Merge along with their own
  `os.listdir` loops. Qt-free, so a batch script enumerates a folder with the
  same extensions and exclusions the GUI shows. An unreadable folder or an
  unknown type now lists nothing rather than raising. This is the tractable
  part of the review's "one file-browser widget" (U4); the widget itself needs
  the Data Selector's browser extracted from its panel first and remains open.
- **The filename sort keys are shared** in `pyirena/core/file_sorting.py`
  instead of existing as three verbatim copies (Data Selector, Data Explorer,
  Data Manipulation) plus a lone order-number copy in Data Merge. A new sort
  mode or a regex fix had to be applied in four places and would be missed in
  some — the failure mode that produced the filter incident. The module is
  Qt-free, so batch and api code can order input files exactly as the GUI
  displays them. All four dropdowns now share one label list (Data Selector's
  double-spaced labels are gone) and one tooltip documenting the recognised
  patterns — two of the four had no tooltip at all.
- **The Markdown report builder moved to `pyirena/core/reporting.py`** from
  `gui/data_selector/report.py`, which now re-exports it. `api/` and `core/`
  may not import from `gui/`, so the builder had to move for the panels and the
  control layer to share it; the generated text is byte-for-byte unchanged
  (verified against the pre-move output). `gui/fmt_utils.py` moved to
  `core/fmt_utils.py` for the same reason and also re-exports. The control
  layer's own smaller Unified-Fit-only Markdown report is gone; its
  agent-specific fit-flag and bounds table is appended to the shared report
  instead.
- **Eight parallel plot-export implementations collapsed into one.**
  `sas_plot.py`, `data_selector/plot_utils.py`, `unified_fit.py` (twice),
  `waxs_peakfit_panel.py`, `sizes_panel.py`, `modeling_panel.py`,
  `contrast_panel.py` and `simple_fits_panel.py` each had their own "add export
  to a plot" code, disagreeing on menu wording ("Save graph as JPEG…" vs "Save
  as JPEG…"), on what was captured (the plot vs the whole window), and on the
  default folder. They now all call `attach_plot_export()`. Curve collection is
  shared with the ITX exporter, so CSV and ITX always agree on what is on the
  plot. `save_itx_from_plot` and `_itx_folder_cmds` moved to the new module and
  are re-exported from `sas_plot` for existing callers.
- Plots that were missing exports gained them: **Simple Fits** had JPEG but no
  ITX and none on its linearization plot; **Modeling** and **Contrast** had
  JPEG but no ITX; the Unified Fit **residual and Porod** plots, the Sizes
  **residuals** plot and the Modeling **residuals** plot had no export menu at
  all. The HDF5 Viewer's graph window keeps its specialised NXcanSAS-aware
  PNG/HDF5/ITX exporters and gains clipboard copy.

- **Clipboard copy, column sorting and CSV export on every table** (feature
  parity review items U1 / A1 / A2 / A10, issue #13). Of the nine tables in
  the GUI, only two supported copying — with two different bespoke
  implementations — so collected values, similarity results and the isotope
  table were dead ends: a handful of numbers required a CSV file round-trip,
  and the similarity results could not be exported at all. All tables now
  share one implementation, `pyirena/gui/table_utils.py`, modelled on
  `file_filter.py`: one module, one behaviour, one tooltip, one test file.
  - **Ctrl+C** copies the selection as tab-separated text (pastes cleanly into
    Excel, Igor Pro and Origin); **Ctrl+Shift+C** includes the column headers;
    right-click offers Copy / Copy with Column Headers / Copy Whole Table and,
    where the panel supports it, *Save as CSV…*. With nothing selected, copy
    falls back to the whole table instead of silently doing nothing.
    Non-contiguous ctrl-click selections copy as a compact block of just the
    selected rows and columns.
  - **Click-a-header sorting** on the Data Selector *Tabulate Results* table,
    the HDF5 Viewer *Collect* and *Multi-Collect* windows, and the Data
    Manipulation similarity results. Numeric columns sort numerically —
    9 < 10 < 100, not "10" < "9" — with blanks and placeholders ("—", "(ref)")
    always last. Turning sorting on does not reorder anything by itself: rows
    stay in the order the panel supplied until a header is clicked, and (on
    Qt ≥ 6.1) a third click returns to that natural order. Sorting is
    deliberately *not* enabled on the Contrast results table (rows are grouped
    under section headers) or the Multi-Collect item list (its row order
    defines the output columns).
  - **Save CSV… for the Data Manipulation similarity results**, which
    previously had no export of any kind — filename, p-value, longest run,
    number of points and accepted/rejected at full precision.
  - The Contrast isotope table and the Multi-Collect item list gained copy
    support; a whole-table copy of the isotope table includes the isotope
    picked in each drop-down.
- `pyirena/tests/test_gui_table_contract.py` — a source-level guard that fails
  if a new module constructs a `QTableWidget` without calling
  `attach_table_copy()`, or hand-rolls clipboard code again. The standard UX
  contract it enforces is now written down in `docs/developer_adding_features.md`.

### Fixed

- **The GUI test modules errored instead of skipping without Qt.** The new
  table, plot-export and report tests used `pytest.importorskip`, which pytest
  ≥ 8.2 treats as a *broken* module when the import raises a plain ImportError —
  and `pyirena.gui._qt` raises one carrying an installation hint. The plain
  (no-GUI) CI job reported collection errors rather than skips; they now skip
  cleanly, like the older GUI tests.
- **Descending sorts put files without the pattern first.** Sorting a folder by
  Temperature ↓ opened with every unmatched file — logs, notes, a stray average
  — above the hottest measurement, because reversing the list also reversed the
  "no pattern" sentinel. Unmatched files now sort last in both directions, as
  the tooltip has always claimed, and keep their relative order among
  themselves. Affects all four file browsers.
- **An unknown WAXS peak shape aborted the whole report.** The area of a peak
  is recomputed for files written before areas were stored; a shape this build
  does not recognise raised out of `peak_area()` and no report was produced at
  all. It now costs one table cell ("N/A"). A peak area with no uncertainty
  shows "—" rather than "± 0", which read as a measured zero.
- **Igor wave names could collide.** Names were truncated to Igor's 31-character
  limit *after* the `_01`/`_02` index was appended, so two curves with long,
  similar labels produced the same wave name and Igor silently overwrote the
  first. The label is now truncated to leave room for the index.
- **The Simple Fits linearization exported no uncertainties.** `linearize()`
  propagates dY into linearized space, but the panel discarded it; the CSV and
  ITX now carry it (as an Igor error wave), and the grey out-of-range points
  are labelled *Data outside fit range* so they are not mistaken for one.
- **CSV export from the Collect and Multi-Collect windows corrupted rows whose
  labels contained a comma.** Both windows built their CSV with
  `",".join(...)`, so an item label such as `Rg, A` silently split into two
  columns. All CSV writing now goes through `rows_to_csv_text()`, which uses
  the stdlib `csv` writer (proper quoting) and formats floats with `%.10g` —
  the same precision the ITX exporters use.
- **Contrast results table: Ctrl+C did not work.** The shortcut was declared on
  a `QAction` created inside the context-menu handler, so it never fired while
  the menu was closed. It is now a real shortcut on the table.

## [1.1.0b6] - 2026-08-07

### Added

- **Cell/row/column selection and clipboard copy in Data Selector's Tabulate
  Results table.** The table only supported selecting whole rows, and had no
  clipboard support at all — copying a selection silently did nothing, so
  pasting elsewhere produced whatever had been on the clipboard before.
  Selection is now cell-based, so individual cells, whole rows or columns
  (via header click), and non-contiguous multi-selections (ctrl-click) are
  all supported. **Ctrl+C** copies the selection as tab-separated text;
  **Ctrl+Shift+C**, or right-click → *Copy with Column Headers*, includes the
  header row for the selected columns. Non-contiguous selections copy as a
  dense grid over just the selected rows and columns, so picking two
  non-adjacent columns pastes as a clean two-column block rather than
  everything in between.

### Changed

- **Minimum Python raised to 3.10.** `requires-python` said `>=3.9` while
  `environment.yml` already required `>=3.10,<3.14`, so the package advertised
  a floor no shipped or tested environment used — and `CLAUDE.md` told
  contributors to write 3.9-compatible code because of it. The floor is now
  `>=3.10` consistently across `pyproject.toml`, `conda/meta.yaml`, the CI
  matrix (now 3.10/3.11/3.13) and the README badge; the 3.9 classifier is
  dropped. `match` and PEP 604 `X | None` are allowed in new code, though
  existing modules keep `Optional[X]` — follow the file you are editing.
  No runtime behaviour changes; this only stops advertising support for a
  version that was never tested.
- **Ruff now selects `E`, `F`, `W`, `I`** (was the `E`/`F` defaults), with
  `E501` added to the existing ignore list — line length stays a formatting
  target, not a lint error. This sorts imports repo-wide and strips trailing
  whitespace across 131 files; the change is import order and whitespace only,
  no logic touched, and the full suite passes unchanged. The `UP` and `B` rule
  sets from the house standard are deliberately still off — enabling them
  rewrites ~1050 type annotations and surfaces ~100 findings that need hand
  review, which belongs in its own pass.
- Ruff `target-version` bumped `py39` → `py310`; the vestigial `[tool.black]`
  section is removed (black was never a dev dependency here, and `ruff format`
  supersedes it).
- `requirements.txt` declared `numpy>=1.22.0` while `pyproject.toml` requires
  `numpy>=2.0` for `numpy.trapezoid`; installing from the former produced an
  environment where `pyirena.core` fails on import. Floors now agree.
- `IMPROVEMENT_PLAN.md` renamed to `PLAN.md` to match the naming used across
  the other packages.

## [1.1.0b5] - 2026-08-07

### Fixed

- **Invariant: "Refit background from saved ranges" ignored the parameter
  "Fit?" checkboxes.** With complex background enabled and prefit Q windows
  saved, the replay refit `BG_B`, `BG_P` and `BG_flat` unconditionally — only
  `BG_P`'s checkbox was consulted, and only to choose between fitting P or
  holding it. A user who unchecked B/P (for example setting `BG_B` = 0 to drop
  the low-Q power-law term for part of a sequence) while still needing the flat
  background refit had their held values silently overwritten on every
  Calculate. The checkboxes now gate the replay: a background parameter with
  "Fit?" unchecked is left alone, `BG_B` held skips the power-law refit
  entirely, and `BG_flat` held skips the flat refit. `SimpleFitModel.
  prefit_background()` takes a `fixed_params` argument for this, and
  `batch.fit_simple()` passes the user's fixed parameters through, so
  scripted and batch runs behave the same as the GUI.

## [1.1.0b4] - 2026-08-06

### Added

- **Parallel Monte-Carlo uncertainty in the Modeling tool.** Each MC pass is an
  independent refit of noise-perturbed data, so passes now run across worker
  processes instead of one after another. This was impractical before for the
  models that need it most: a complex form factor (core-shell,
  core-shell-shell) or several populations can take minutes per pass, making a
  20-pass estimate a coffee break. A new **cores** spin box next to
  *Modeling → Passes:* controls it — **auto** (the default: all cores but two),
  **1** for the previous serial behaviour, or an explicit count. The same
  setting is available to scripts as `fit_modeling(..., mc_workers=N)` and as
  the `"mc_workers"` key in the exported JSON config. See
  [MC Uncertainty](docs/modeling_gui.md#mc-uncertainty).
  - Uncertainties are unchanged by the worker count: all the noise is drawn in
    the parent process before any pass starts, so a given run produces the same
    numbers serially and in parallel.
  - Short runs stay serial automatically. The first pass is timed and the pool
    is only started when the remaining passes are projected to take more than a
    few seconds, so simple models never pay for worker startup.
  - If the host cannot start worker processes the passes fall back to serial
    with a warning rather than failing.
- **Cancel for Modeling MC runs.** The **Calc. Uncertainty (MC)** button becomes
  **Cancel MC** while a run is in progress. Cancelling does not interrupt passes
  already in flight, so it takes effect within roughly one pass, and the
  uncertainties from the passes that did finish are still reported along with
  how many were used.

- **Regular expressions in every file Filter box.** The documentation promised
  grep-like filtering, but only the HDF5 Viewer actually interpreted the Filter
  text as a regex — Data Selector, Data Manipulation and Data Merge did a plain
  case-insensitive substring test, so patterns like `60C|100C`, `0[12]min`,
  `^sample`, `\.h5$` or `^(?!.*bkg)` silently matched nothing. All four
  browsers now share one implementation (`pyirena.gui.file_filter`) using full
  Python regular expressions matched anywhere in the name, case-insensitively.
  A plain fragment such as `Rg50` behaves exactly as before, and an incomplete
  or invalid pattern (common while typing) falls back to a substring match
  rather than emptying the list.

- **`pyirena-doctor` — installation troubleshooter.** A new console command
  that reports the running interpreter, every required and optional
  dependency, and whether the `pyirena-gui` launcher uses the same Python as
  your `pip`. Each dependency is classified as installed, missing, or
  *installed but failing to load* — the three states that the previous error
  messages collapsed into one. Users can paste its output into a bug report.
  See [Installation → Troubleshooting](docs/installation.md#troubleshooting).

### Fixed

- **"GUI dependencies not installed" when they demonstrably were.** `ImportError`
  covers both a genuinely absent package and one that is present but whose
  binaries will not load (wrong CPU architecture, missing system libraries,
  shiboken6/PySide6 version skew). pyIrena reported the second as the first,
  telling users to reinstall something they already had. Qt import failures are
  now diagnosed: the message says which of the two happened, where the package
  lives, the original loader error, the interpreter that failed, and — for
  recognised errors such as `incompatible architecture` or `libGL.so.1` — the
  likely cause and fix.
- **Wrong-environment installs are now detected.** When the `pyirena-gui`
  launcher script points at a different Python than the one running, both paths
  are printed along with the `python -m pip` / `python -m pyirena.gui.launch`
  form that keeps them in sync. This is the most common cause of the report
  above.
- **GUI startup errors are no longer flattened into a dependency message.** An
  `ImportError` from an internal module is now identified as a probable bug
  (with an issue-tracker pointer) rather than blamed on missing packages, and
  the full traceback is always written to `~/.pyirena/logs/gui.log`.
  `PYIRENA_DEBUG=1` also prints it to the terminal.
- `pyirena-viewer` imported Qt directly instead of going through
  `pyirena.gui._qt`, so it produced a bare `ModuleNotFoundError` instead of the
  diagnosed message.
- A Modeling test (`test_export_includes_selected_fit_method`) failed instead
  of skipping in environments without Qt, because `pytest.importorskip` treats
  a re-raised plain `ImportError` as a broken module rather than a missing one.

### Changed

- The **Filter** placeholder and tooltip in Data Selector, Data Manipulation,
  Data Merge and the HDF5 Viewer now come from one shared string and document
  the regex syntax inline, replacing the previous "text filter…" /
  "Enter text to filter files..." hints.

## [1.1.0b3] - 2026-08-03

### Added

- **Q-unit selector for text-file import.** Text files (`.dat`/`.txt`/`.csv`)
  carry no unit metadata, so pyIrena assumed Q was always in 1/Å. A new
  *Data Selector → Configure → Text File Options → Q unit in text files*
  setting (default 1/Å, also offering 1/nm, 1/pm, 1/µm, 1/mm) lets you tell
  pyIrena what unit your files actually use; Q (and dQ) is converted to 1/Å
  on import. Applied everywhere text files are read: Data Selector, Data
  Merge, Data Manipulation, the `pyirena.batch` API (`fit_sizes`,
  `fit_unified`, `fit_waxs`, `merge_data`, etc.), and
  `pyirena.plotting.plot_saxs`. The assumed unit is recorded in the
  converted HDF5 sibling's provenance, and changing the setting for a file
  you've already converted automatically invalidates the stale cached
  sibling rather than silently reusing it. See
  [Q units](docs/data_import_and_cleaning.md#q-units).
- **CSV file support in Data Selector, Data Merge, and Data Manipulation.**
  Comma-separated `.csv` files are now recognized alongside `.dat`/`.txt`
  everywhere those text formats were already supported:
  - Data Selector's "Text Files" / "All Supported Files" filters, and the
    load / plot / report / ASCII-export paths (auto-converted to a cleaned
    NXcanSAS `.h5` sibling on first use).
  - Data Merge's and Data Manipulation's file-type dropdown (now
    **Text (.dat/.txt/.csv)**) and reference-file loader.
  The underlying text reader (`readTextFile`) now auto-detects a comma
  delimiter from the first data-bearing line, so headerless or
  header-labeled CSV exports (e.g. `Q,I,dI` columns) parse the same way as
  whitespace-separated files.
- **Data Selector — GitHub-based update notification.** Replaces the retired
  Igor Pro APS/ANL-server version check (which also is not applicable to
  pyIrena). On startup, and at most once a week thereafter, the Data Selector
  reads GitHub's public `releases/latest` endpoint (pre-releases/betas are
  excluded automatically) and shows an info banner with a link to the GitHub
  releases page when a newer stable version is available.
  Stdlib-only (`urllib.request`, no new dependency), fails completely
  silently offline/on any error, and never blocks startup — the network call
  runs on a background `QThread` (`UpdateCheckWorker`). New
  `pyirena/version_check.py`; opt-out via a new **"Check for new pyIrena
  releases on startup"** checkbox in Data Selector's Configure… dialog
  (`check_for_updates`, default on). See `docs/gui_quickstart.md`.

## [1.1.0b2] - 2026-07-24

Beta 2 is a cleanup / hardening release on top of the 1.1.0b1 slit-smearing
beta. It closes gaps found in an independent code review — the control/MCP
write surface, packaging metadata, and CI — with **no change to the fitting
science or numerical results**.

### Security

- **Control / MCP file access is now confined to `PYIRENA_DATA_ROOT`.** The
  read-only `pyirena.api` discovery/data tools already resolved user paths
  through `resolve_safe*`, but the mutating control tools
  (`open_dataset`, `save_fit`, `save_sizes_fit` — all exposed as MCP tools)
  built `Path(...)` directly, so setting `PYIRENA_DATA_ROOT` did **not** confine
  the part of the API that can read *and* write HDF5 files. All three now
  resolve inputs with `resolve_safe_file()` and write targets with
  `resolve_safe(..., must_exist=False)`, rejecting absolute paths and `..`
  traversal outside the root with a `PATH_NOT_ALLOWED` error dict (never an
  exception across the API boundary). The stale "information disclosure only"
  note in `api/_paths.py` was corrected.

### Fixed

- **Control API saved slit-smeared fits as if they were pinhole fits.** A
  control-/MCP-driven Unified Fit or Size Distribution on slit-smeared data
  stored `slit_length = 0` and no ideal (pinhole) model curve, silently losing
  the slit-smearing provenance that the batch and GUI paths record.
  `save_fit()` now writes the smeared model, the ideal `intensity_model_ideal`
  curve, and `slit_length`; `save_sizes_fit()` now writes `slit_length`,
  `data_is_slit_smeared`, and `intensity_model_ideal`. The `open_dataset`
  JSON schema also gained the `use_slit_smeared` property (the function always
  accepted it, but schema-driven clients could not request the `_SMR` dataset).
- **Control API `output_path` produced an incomplete data file.** Saving a fit
  to a *new* `output_path` passed the path straight to the result writer, which
  created a results-only HDF5 (no reduced-data SASdata group) that
  `readGenericNXcanSAS()` could not reopen as reduced data. Both save functions
  now seed a new target from the source file (`copy_and_strip_results`, reduced
  data + metadata, stale results stripped) before appending — the original file
  is never modified.
- **Declared NumPy floor now matches the code (NumPy ≥ 2.0).** The core calls
  `numpy.trapezoid` (a NumPy 2.0 API) in ~30 places, but `pyproject.toml`
  declared `numpy>=1.22` and the conda recipe `numpy>=1.20`, so an install on
  NumPy 1.x would crash. The floors were raised to `numpy>=2.0` in both.
- **Slit-smearing engine now validates its inputs** (`core/smearing.py`). When
  smearing is active (slit length > 0), `build_smearing_matrix`, `smear_curve`,
  and `build_extended_q` raise a clear `ValueError` on non-finite, unsorted, or
  too-short `q`, on `n_l < 2`, and on a q/intensity length mismatch — instead of
  returning silently wrong output or a cryptic `IndexError`. The `slit_length
  <= 0` no-op contract is preserved (inputs pass through untouched), so the
  already-sanitised fit paths are unaffected.
- **Import Igor Experiment: Irena/Nika "Use QRS Names" data not recognised.**
  Igor experiments (`.h5xp`/`.pxp`) whose reduced 1-D data followed the common
  Irena/Nika "Use QRS Names" convention — waves named `Q_<folder>` /
  `R_<folder>` / `S_<folder>` (uppercase prefix, data-folder name as suffix;
  R = intensity, S = error) — imported as **zero** samples: every folder was
  silently skipped. The importer's wave-name picker tables only knew pyirena's
  own lowercase `q_<folder>` output and a few other conventions. Added the
  uppercase QRS triple to all three techniques (USAXS/SAXS/WAXS) in both
  `WAVE_PICKERS_H5XP` and `WAVE_PICKERS` in `pyirena/io/pxp_to_nexus.py`. A real
  desktop-SAXS (Xenocs, ANSTO) file that previously imported 0/44 samples now
  imports 44/44. *Note: a follow-up will make Igor wave-name matching fully
  case-insensitive, since Igor names are case-insensitive by design.*
- **"Teubner-Strey" model name typo.** The Simple Fits model was misspelled
  "Treubner-Strey" (the model is named for M. Teubner & R. Strey, J. Chem.
  Phys. 87, 1987) in the registry key, internal function name, derived-value
  docs, and the Igor-compatibility wave name (`SimFitTreubnerStreyI` →
  `SimFitTeubnerStreyI`). Corrected throughout `pyirena/core/simple_fits.py`,
  `pyirena/io/igor_names.py`, docs, and README. Old saved GUI state,
  NXcanSAS result files, and JSON exports that stored the misspelled name
  keep loading correctly — `SimpleFitModel.set_model()`,
  `SimpleFitModel.from_dict()`, and the Simple Fits panel's `load_state()`
  all transparently map the legacy spelling to the corrected one via a new
  `_resolve_model_name()` helper, so no existing files need to be migrated.

### Changed / Internal

- **Packaging metadata cleaned up.** Added a real `plotting` extra
  (`pip install pyirena[plotting]`, referenced by the docs but previously
  undefined); `import pyirena` no longer eagerly imports matplotlib (the
  `plot_saxs` convenience export is now lazy via module `__getattr__`);
  migrated to the SPDX `license = "MIT"` form (+ `license-files`) ahead of the
  setuptools deprecation; the conda recipe now matches `pyproject.toml`
  (Python ≥ 3.9, NumPy ≥ 2.0, SciPy ≥ 1.8, and the previously-missing `igor2`
  dependency) with a documented `sha256` placeholder; and the test suite is no
  longer packaged into the built wheel/sdist (developers and the conda recipe
  build from the source archive).
- **CI is green and bounded.** Fixed the 12 outstanding `ruff` findings so the
  lint job passes; added `pytest-timeout` with a 300 s per-test watchdog so one
  hung test can no longer stall a CI job; and fixed the long-hanging
  `test_export_includes_selected_fit_method` — it blocked forever on an
  unmocked modal `QMessageBox` overwrite confirmation, now auto-accepted.
- **New contract tests.** The full MCP surface (68 tools) and all 51 control
  JSON schemas are now checked for structural validity and schema↔signature
  parity (this would have caught the `use_slit_smeared` schema gap); new
  behaviour tests cover control-API path confinement, `output_path`
  completeness, and slit-smearing provenance on save; and the slit-smearing
  input validation is covered by unit tests.
- **Docs / repo hygiene.** `docs/distribution.md` no longer tells maintainers to
  edit a version string in `pyirena/__init__.py` (the version is single-sourced
  from `pyproject.toml`); `publish.yml` now verifies the release tag equals the
  `pyproject.toml` version before publishing to PyPI; and the tracked
  `scratch_sizes_diagnosis/` development directory (already gitignored) was
  removed from version control.

## [1.1.0b1] - 2026-07-22

### Added

- **Slit smearing across all fitting tools (USAXS/Matilda).** pyIrena can fit
  slit-smeared USAXS data directly by smearing the **model** to match the data
  (Lake infinite-slit, `I_sm(q)=(1/SL)∫₀^SL I(√(q²+l²))dl`) — the data are never
  modified and fitted parameters are always ideal-space (pinhole-equivalent).
  Core engine `pyirena/core/smearing.py` (`SlitSmearer`, `smear_model`,
  `smear_curve`, a fixed sparse operator `W` so a fit loop costs one matvec per
  iteration). Slit-smeared data are detected from NXcanSAS (`Q@resolutions`
  containing `dQl` + a scalar `dQl`); files with both a desmeared and a
  slit-smeared copy (Matilda) show a **"Slit smeared data"** checkbox to select
  which to load. Wired through **Unified Fit** (all levels + background + local
  Guinier/Porod cursor fits + "Show selected level" overlay; invariant computed
  from the ideal model), **Size Distribution** (G matrix smeared once, recovered
  distribution is ideal-space), **Simple Fits** (each analytic model smeared;
  Invariant disabled on smeared data with a message), and **Modeling** (total +
  per-population curves smeared). **Fractals** smears its comparison overlay.
  GUI control shared via `pyirena/gui/slit_smearing_ui.py::SlitSmearingMixin`
  (used by every fitting panel, incl. Unified Fit). Scripting/MCP contract:
  `load_slit_smeared: true` (+ optional `slit_length`) in a tool's JSON block,
  enforced with a hard error on files lacking `dQl`. Saved results record
  `slit_length` / `data_is_slit_smeared` and an ideal (`*_ideal`) model curve
  alongside the smeared one. See `docs/slit_smearing.md`.
- **Data Explorer — "Show all attributes" checkbox.** The HDF5 tree browser
  previously only surfaced a curated set of attributes (`NX_class`, `units`,
  `analysis_type`, ...), silently hiding others such as `data_is_slit_smeared`
  / `slit_length` on `unified_fit_results`. A new checkbox above the tree
  toggles between the curated set and every attribute on a node; internal
  bookkeeping (`_pyirena_config` and any other `_`-prefixed attribute) stays
  hidden either way. Toggling preserves the tree's current expand state.
- **Data Merge — slit-smearing provenance.** Merging a slit-smeared USAXS curve
  with a pinhole SAXS curve produces a slit-smeared output: the merged file
  gets a `dQl` dataset (so downstream tools auto-detect it) and the merge
  provenance records `slit_length_ds1/ds2` and `slit_length_merged`. Two inputs
  with different nonzero slit lengths warn (the larger is kept). Optimization is
  unchanged — the slit length sits at/below the SAXS Qmin, negligible in the
  overlap.

### Changed

- **Data Manipulation — slit-smearing safety.** Subtract/divide refuse to mix a
  slit-smeared curve with a pinhole one (or two different slit lengths); the
  check now lives in the core engine (`DataManipulation.check_slit_compatible`),
  so batch scripting inherits it. Manipulation/merge outputs drop any stale
  `_SMR` twin entry copied from the source and clear an orphaned `dQl`, so a
  later slit-smeared load can't return an inconsistent curve.
- **Create Report** notes when a tool's saved results used slit smearing
  (Unified/Sizes/Simple/Modeling), reading `slit_length` from the HDF5 file.

### Fixed

- **Clear error for non-scattering / empty data files, instead of a different
  cryptic crash per tool.** A file whose Q/Intensity arrays are empty (sample
  didn't scatter, aborted measurement, corrupted file) used to load
  "successfully" and then crash deep inside each tool with an unrelated numpy
  error (`len() of unsized object`, `too many indices for array: array is
  0-dimensional`, mismatched broadcast shapes, ...). `readGenericNXcanSAS` /
  `readSimpleHDF5` (`pyirena/io/hdf5.py`) now raise a single
  `NoScatteringDataError` right at load time with a message that says what's
  actually wrong, shown identically by every tool (Unified Fit, Size
  Distribution, Modeling, Simple Fits, SAXS Morph, Data Selector "Create
  Graph"). The batch API (`_load_data` and all `fit_*` functions) already
  catches load errors per file and returns `None`, so unattended/scripted
  batch runs now skip a bad file with one clear log line instead of risking a
  crash — verified across `fit_simple`, `fit_sizes`, and the Data Selector's
  batch-script worker. "Create Graph" now also reports which files were
  skipped and why instead of silently plotting nothing for them.
- **Simple Fits — displayed model is now slit smeared.** "Graph model" and
  auto-graph previously drew the ideal (sharp) curve even with smearing on, so
  Sphere/Spheroid/Teubner-Strey oscillations looked unsmeared; `compute()` now
  smears, matching the fit. Guinier/Porod linearization is labelled best-effort
  (ideal-space) for smeared data (no closed-form linearization exists).
- **Background prefits are ideal-space under smearing** (Size Distribution and
  Simple Fits): the prefit power-law/flat is now smeared before comparison, so
  the returned B/P/flat are not double-smeared by the main fit.
- **Modeling** saves the ideal (pinhole) model curve (`model_I_ideal`) alongside
  the smeared one, and no longer reuses a wrong-length cached G matrix when
  producing that ideal curve.
- **Fractals** no longer silently shows an unsmeared overlay when smearing the
  comparison curve fails — it logs and flags the discrepancy.
- **Unified Fit — large slit-smeared fit slowdown fixed (~50× on affected
  fits).** ETA/PACK are now only free fit parameters when a level's
  *correlations* are enabled. With correlations off they have no effect on the
  model, so fitting them made the least-squares problem rank-deficient and the
  solver thrashed for thousands of no-op iterations — cheap for pinhole data but
  badly amplified by slit smearing (each iteration evaluates the model on the
  extended grid). A real 2-level USAXS fit dropped from ~2.5 s to ~0.05 s with
  identical χ² and parameters. Also speeds up pinhole fits and fixes a latent
  correctness issue (fitting parameters that cannot change the fit).
- **Size Distribution — power-law/complex background now displayed and
  subtracted slit-smeared.** The background preview drew and subtracted the
  *ideal* (pinhole) background on slit-smeared data, so the curve sat below the
  data and looked unfitted and `I−bg` was wrong; it now smears the background
  for both display and subtraction, matching the fit.
- **Size Distribution / Modeling panels — layout no longer forced too wide.**
  The shared "Slit smeared data" control row is now compact (short labels;
  status text on its own wrapped line), so it no longer stretches narrow control
  panels and hides widgets. The Size Distribution control panel is now
  user-widenable via the splitter (420 px minimum instead of a hard-fixed
  width).
- **Size Distribution — fitted model + background now displayed slit-smeared.**
  The fit was correct, but the result plot added the *ideal* (pinhole)
  background to the (smeared) model scattering and drew the ideal "Complex bg"
  curve, so the red model+background curve visibly misfit the smeared data. All
  three display paths (fit result, "Graph model", background preview) now show
  the smeared background — and the "Graph model" preview smears the scattering
  too — matching the data.
- **Data Merge — select the slit-smeared copy.** When a merge input file carries
  both a desmeared and a slit-smeared (USAXS) copy, a **"Use slit-smeared copy"**
  checkbox appears under that dataset's file list; checking it reloads the
  slit-smeared data so the merged output is written slit smeared (`dQl`).
  Previously only the desmeared `@default` entry could be loaded, so the
  provenance plumbing added earlier had no way to be triggered from the GUI.

## [1.0.1] - 2026-07-15

### Added

- **Simple Fits — Invariant (calculation model).** New `Invariant` entry in the
  Simple Fits model selector calculates the Porod invariant
  Q* = ∫q²·I(q) dq = 2π²·Δρ²·φ(1−φ) and, with contrast and absolute-scale data,
  the two-phase volume fraction φ. Port of Igor Irena `IR3J_CalculateInvariant`
  with the same numerics: data are extended to Q=0, background is removed with
  the standard complex background (B·Q⁻ᴾ + flat, incl. the prefit-between-cursor
  buttons), and the invariant is read where the running integral flattens
  (QmaxUsed, smoothed-derivative plateau detection). Optional (off by default,
  beyond Igor) Porod tail extension Kp/Qmax compensates high-Q truncation. The
  running integral is drawn on a right axis of the I(Q) graph and stored in the
  HDF5 file (`Q_integral`/`running_integral`); results land in `derived/`
  (`Invariant` [cm⁻⁴], `VolumeFraction`, `PhiOneMinusPhi`, `QmaxUsed`,
  optional `PorodTail`/`PorodKp`). Wired through GUI, JSON scripting, batch
  (`fit_simple` with `model: "Invariant"`), MCP read tools, Data Selector
  graph/report/tabulate, HDF5 Data Explorer trend plots, and Igor experiment
  export (`SimFitInvariantI`, `SimFitInvariantIntegral` waves). Registry
  entries now support a `calculation: True` flag for no-fit methods; χ²,
  residuals, linearization and MC uncertainty are disabled/None for these.
- **Simple Fits — background prefit replay for the Invariant.** The
  "Fit B/P btwn cursors" / "Fit Flat btwn cursors" buttons now remember the
  Q windows they were used on (`bg_prefit` in the model/state/JSON/HDF5
  setup). A new "Refit background from saved ranges" checkbox (Invariant
  options) re-runs those prefits on the full data before every Calculate —
  identically in the GUI and in scripted/batch runs — so batch invariant
  results re-determine the background per file instead of trusting exported
  B/P/flat values. The BG_* "Fit?" checkboxes remain visible in Invariant
  mode (BG_P's controls fit-both vs. hold-P prefit). Also fixed: right-axis
  label now survives model switching, and the residuals panel is hidden for
  calculation models.
- **Developer guide: adding features.** New
  `docs/developer_adding_features.md` — the master checklist of every wiring
  point (core, GUI, HDF5, JSON/batch, MCP, Data Selector, Data Explorer, Igor
  export, docs, tests) a new feature must touch, with the Invariant as a
  worked case study.

## [1.0.0] - 2026-07-15

### Added

- **Size Distribution — "Set r min/max from Q range" button.** The Size Grid group
  in the Sizes panel gains a button that fills both `r min` and `r max` from the
  current cursor Q-range via `r ≈ π/Q`, rounded outward to tidy 1-2-5 values so
  the grid comfortably brackets the resolvable range. Backed by a shared helper
  `pyirena.core.form_factors.r_bounds_from_q_range()` used by the GUI and the
  recommendation tool alike.

### Changed

- **Size Distribution — suitability & auto-setup recognise broad distributions.**
  `recommend_sizes_setup` (exposed as the `suggest_sizes_setup` control tool, used
  by pyirena-ai) no longer rejects the everyday case of a **broad size
  distribution on a low-Q power law + high-Q flat background** (precipitates in
  metals; pores in rocks/minerals/solids). Previously it required a Guinier knee
  and ≤2 detected "levels", so these routine datasets were flagged *not suitable*
  and the agent refused to fit them. Now:
  - Suitability is based on whether a Q-band exists where the particle signal is
    clearly above the fitted background (I(Q) ≳ 2× background, relaxing the ratio
    when the signal is weak), not on Guinier knees or level counts.
  - The inversion Q-range is chosen from that signal-to-background band (trying the
    full complex background first, then flat-only), which also keeps the noisy,
    background-dominated high-Q tail out of the inversion.
  - The high-Q flat background is detected even when the feature detector labels it
    a "guinier_plateau"; the estimated `flat_background` level is returned.
  - Recommended `r_min`/`r_max` now come from `r_bounds_from_q_range` (tidy,
    slightly-padded bounds) instead of raw π/Q.
  - `recommended` gains `flat_background`; `features` gains `signal_to_bg_ratio`
    and `inversion_span_decades`. "Multiple knees / several levels" is now an
    advisory warning, not a hard *unsuitable* verdict.

### Fixed

- **Size Distribution — Monte Carlo distribution shape.** The Monte Carlo
  (McSAS) method reported a volume distribution shifted to roughly 2× too large
  in radius. Its solution vector `x = A·count` is already a per-bin *volume
  fraction* (the G-matrix columns are intensity per unit volume fraction), but
  the post-processing multiplied it by `V(r) = (4/3)πr³` again on the mistaken
  assumption that it was a number distribution. That spurious r³ re-weighting is
  removed, so Monte Carlo now flows through the same normalisation as MaxEnt,
  TNNLS, and Regularization. χ² and volume fraction were already correct (they
  are computed from the un-weighted solution), so only the displayed/saved
  distribution shape is affected. A ground-truth inversion test confirms the
  recovered Rg and mean radius now match the input distribution.

- **Size Distribution — Monte Carlo confined to the resolvable size band.**
  Monte Carlo contributions are now restricted to `r ∈ [π/Q_max, π/Q_min]` (the
  standard McSAS size-range rule). Previously, radius bins outside the range the
  data can constrain were left free, so Monte Carlo dropped stray contributions
  into very large-r bins — which, being r²-weighted, wildly inflated (and
  destabilised, run-to-run) the reported **Rg** — and into sub-resolution
  small-r bins, producing a spurious spike that soaked up high-Q noise. On a real
  USAXS dataset this brought Monte Carlo's Rg and peak radius into agreement with
  MaxEnt and Regularization (previously Rg varied from ~2000 to ~19000 Å between
  runs and the density peak sat at the smallest bin). The deterministic methods
  already suppress these bins via their priors; this gives Monte Carlo the same
  discipline.

- **Size Distribution — Regularization robustness when χ² = M is unachievable.**
  When the discrepancy target χ² = M cannot be reached (fit-window error bars too
  small, or the model cannot describe the data at the window edges — typically
  noisy, background-dominated high-Q points), the fallback previously selected
  the *least-smoothed* (minimum-χ²) solution. That could collapse into a single
  huge spike at the smallest radius bin and made the result extremely sensitive
  to the exact high-Q cut-off. The fallback now selects the *smoothest* solution
  whose χ² is within a small factor (`regularization_fallback_factor`, default
  1.05) of the achievable minimum χ², and logs a warning suggesting the user trim
  the high-Q end of the inversion window. Results within a well-chosen fit window
  are unchanged.

### Documentation

- **Size Distribution methods guide** (`docs/sizes_methods.md`): updated the
  Regularization fallback description and added a *"Choosing the fit (inversion)
  window"* section explaining why background-dominated high-Q points destabilise
  the inversion and how to pick the Q-range; a *"Reading the distribution plot on
  a log radius axis"* note explaining the dV/dr-vs-dV/dln(r) visual effect (why a
  distribution can look concentrated at tiny r yet report a large Rg); and a note
  on the Monte Carlo resolvable size band and setting `R max ≈ π/Q_min`.

## [0.10.0] - 2026-07-10

### Added

- **Data Explorer — Collect Values: Export to ITX.** The Collect Values graph now
  has a "Save ITX" toolbar button and a "Save ITX (Igor Pro)…" right-click context
  menu entry. Exports the scatter plot (X, Y, and optional error-bar waves) as an
  Igor Pro Text (.itx) file with display commands, axis labels, and a legend.

- **Export to Igor — derived parameters in wave notes.** The h5xp export now
  includes additional derived / calculated parameters in the wave notes of result
  waves:
  - **Unified Fit**: per-level `Sv` (surface-to-volume ratio, m²/cm³) and
    `Invariant` (scattering invariant, cm⁻⁴) are now written as `Sv_L{n}` and
    `Invariant_L{n}` in the wave note of `UnifiedFitIntensity` (when present in
    the HDF5 file, i.e. the GUI calculated them before saving).
  - **Size Distribution**: `specific_surface_invA` (Å⁻¹) is now written in the
    wave note of `SizesFitIntensity`, `SizesVolumeDistribution`, and related
    distribution waves.

### Fixed

- **Order-number sort is now robust against arbitrary filename suffixes.** The
  "Order number" sort in Data Explorer, Data Selector, Data Merge, Data
  Manipulation, and the HDF5 Viewer file tree now scans `_`-separated segments
  right-to-left and picks the first bare integer (digits only, no letters).
  This handles all pyirena-generated suffixes (`_merged`, `_mrg`, `_scaled`,
  `_trimmed`, `_rebinned`, `_avg`, `_sub`, `_div`) as well as Irena-style `_mrg`
  and unit-bearing tokens like `_10min` or `_5C`, without requiring an explicit
  allowlist.

### Internal / Code quality

- **Local Guinier / power-law fits unified into `core.unified`.** The inline
  `curve_fit` Guinier and Porod/power-law models in `gui/unified_fit.py` and the
  `fit_local_guinier` / `fit_local_power_law` tools in
  `api/control/unified_fit.py` now delegate to a single implementation
  (`core.unified.fit_local_guinier`, `core.unified.fit_local_power_law`). Numbers
  are pinned by new regression tests (`tests/test_local_fits.py`) so both the GUI
  and API paths reproduce their previous results exactly. As part of this the GUI
  power-law fit now drops non-positive-intensity points before fitting (matching
  the API), a minor correctness improvement in the over-subtracted-background edge
  case.
- **`sphere_amplitude` deduplicated.** `core/modeling.py` and `core/unified.py`
  shared identical Born-Green sphere-amplitude code; there is now one
  module-level `core.unified.sphere_amplitude`, reused by both.
- **`pyirena/gui/_qt.py` single Qt import point.** ~30 GUI modules previously
  repeated a `try: from PySide6 … except ImportError: from PyQt6 …` block; they
  now import Qt names from one shim. Removes the duplication and all of the
  associated unused-import lint noise.
- **Logging follow-through.** 137 silent `except …: pass` handlers in `gui/` now
  emit `log.debug(…, exc_info=True)` (file-only DEBUG, no console noise) so
  swallowed exceptions are observable. Diagnostic `print()` calls in `gui/` and
  `core/` were converted to module loggers (intentional CLI echoes and the
  dependency-missing message in `gui/launch.py` were kept).
- **Lint clean.** `ruff check pyirena` now reports zero findings (was 63:
  45 F401 + 18 F841); each unused local was individually reviewed.
- **CI.** Added a `test-gui` job that installs the `[gui]` extra and runs the
  suite headless (`QT_QPA_PLATFORM=offscreen`) so the optional-dependency tests
  (periodictable / xraydb / Dans_Diffraction) execute in CI.
- **Version single-sourced.** `pyirena.__version__` is now read from installed
  package metadata (`importlib.metadata`) instead of being duplicated alongside
  `pyproject.toml`.
- **Packaging / repo hygiene.** Development scratch directories
  (`codeFragments/`, `IgorCodeFragments/`, `planning/`) and the root
  `pyIrena_icon.png` are excluded from the sdist (verified with `python -m
  build`). Older changelog entries (0.7.2 and earlier) moved to
  `docs/CHANGELOG_archive.md`.

## [0.9.9] — 2026-07-05

### Fixed

- **ITX export now includes error bars.** Both `save_itx_from_plot` (all analysis
  tool windows) and `save_itx` (HDF5 Viewer) now emit a `Yerr` wave and an
  `ErrorBars` command for each data curve that has measurement uncertainties.
  Model and fit curves without uncertainties are unaffected. The `dI` array is
  stored on the scatter item by `plot_iq_data` and retrieved at export time, so
  the exporter is safe for data with or without error bars.

- **Modeling & Size Distribution: MC uncertainties now reach the graph label.**
  - **Modeling**: the Monte-Carlo result is no longer shown in a transient
    blocking dialog; the ±1σ parameter list is written into the graph-window
    status area (green report), matching the Unified Fit tool. **Results to
    Graph** now shows each fitted parameter's MC uncertainty in the annotation.
    The lookup is group-agnostic, so it works for every population type:
    size distribution (incl. `volume_fraction`, propagated from the fitted
    `scale`), Unified Fit level (G/Rg/B/P), Guinier-Porod (G/Rg1/s1/P/Rg2/s2),
    diffraction peak (position/amplitude/width), and mass/surface fractals.
    Purely-derived quantities with no single fitted parameter (e.g. the Unified
    `invariant`) remain without a ± by design.
  - **Size Distribution**: the MC scalar uncertainties (Rg, Vf, peak r) are now
    retained after the MC run so **Results to Graph** re-displays them in the
    plot annotation; they are cleared on a new fit or data load so stale values
    never carry over.

### Added

- **Complex-background prefit-between-cursors helpers (Simple Fits + Modeling).**
  The power-law background prefactor spans many decades and is painful to set by
  hand, so both tools now offer small "Fit … btwn cursors" helper buttons that
  prefit the background over the cursor-selected Q window (mirroring the Unified
  Fit tool's local-fit pattern). No full fit is run — results are written as
  starting values.
  - **Simple Fits** (shown only when "Complex background" is active): **Fit B/P
    btwn cursors** fits the power-law `B·Q⁻ᴾ`; if P's "Fit?" box is unchecked,
    only B is fit at the current (model-guided) P. **Fit Flat btwn cursors**
    estimates the flat term from the median residual over the window.
  - **Modeling** (Unified-Fit-Level population): **Fit B/P btwn cursors**
    (enabled only when the active population tab is a Unified Fit Level) fits
    that population's B/P, respecting its P "Fit?" box. **Fit Flat btwn cursors**
    sets the global flat Background to the median intensity over the cursor
    window (place cursors on a flat high-Q region).
  - New reusable core helper `fit_power_law_bg_fixed_p()` in `saxs_morph.py`
    (closed-form median estimate of B at a fixed P), alongside the existing
    `fit_power_law_bg()` / `fit_flat_bg()`.

### Changed

- **Modeling: parameter limits auto-update when a value is changed.** Scrubbing
  or typing a new value in any population-tab parameter field (including Scale and
  Contrast) now immediately updates that parameter's lo/hi limits to a 0.2×…5×
  bracket around the new value, clamped to the parameter's hard physical bounds.
  Matches the Unified Fit tool's long-standing behaviour. The "Fix limits?" button
  continues to reset all parameters at once.

- **Simple Fits: complex-background symbols now use B / P / flat.** The checkbox
  now reads `B·Q⁻ᴾ + flat` (was `A·Q⁻ⁿ + flat`) and the power-law prefactor is
  labelled **B** (was `BG_G`), matching the Unified Fit convention (B = prefactor,
  P = exponent). Internal state keys were renamed `BG_G → BG_B` (P and flat keys
  unchanged); the parameter grid shows friendly B/P/flat labels. Old saved states
  containing the previous `BG_G` key load without error — the stale key is simply
  not restored (re-enter B if needed).

### Fixed

- **Simple Fits: bottom linearization graph now tracks the selected model,
  bounds its axes, and hides when unavailable.** The linearization panel
  (Guinier ln(I) vs Q², Porod I·Q⁴ vs Q⁴, …) was only redrawn after a *Fit*, so
  switching the model or pressing "Graph model" left it stuck on the previous
  model's transform (usually Guinier). It now refreshes on Fit, on "Graph
  model", and on model change. Both axes are bounded to the cursor-selected Q
  range (converted to the transform's X units) with a robust 1st–99th-percentile
  Y window computed from the in-range points, so tiny high-Q intensities (huge
  negative ln values) no longer blow up the auto-scale and squash the linear
  region out of view. Models with no linearization (Power Law, Sphere,
  Debye-Bueche, Treubner-Strey, Hermans, etc.) now hide the bottom panel
  entirely, giving the I(Q) plot more vertical space; the panel reappears when a
  linearizable model (Guinier family or Porod) is selected.

- **Modeling: "Fix limits?" now updates distribution, scale, and contrast
  limits (previously a no-op for narrowed fields).** The button is meant to set
  every fit limit to a 0.2×…5× bracket around the current value, but it clamped
  that bracket against each parameter's *current editable limits*
  (`new_lo = max(cur_lo, bracket_lo)`, `new_hi = min(cur_hi, bracket_hi)`).
  Because `max`/`min` can only tighten, the button did nothing whenever the
  current limits were already narrower than the bracket — e.g. a Gauss
  distribution with `mean_size` limits 200–400 stayed 200–400. Form-factor
  parameters *appeared* to work only because their default limits are typically
  much wider than the bracket. Distribution parameters (`mean_size`, `width`,
  etc.) were the most visibly affected. `fix_limits()` now clamps the bracket to
  each parameter's *hard physical bounds* (from the per-type default-limit
  tables, captured at row construction) instead of the current fields, so the
  button always resets the bracket — while still keeping naturally-restricted
  params (fractal dimensions, power-law exponents, deviations) within valid
  ranges. Also added editable low/high limit fields for **Scale** and
  **Contrast** in the Physical Parameters section (previously fixed at
  hardcoded defaults with no GUI control); these are now shown/hidden by
  "No limits?" and driven by "Fix limits?" like all other parameters.

- **Unified Fit: `success` flag is now consistent between GUI and scripting.**
  The GUI showed "Fit completed successfully!" (green bar) for any fit that
  returned a finite result — including ones with large chi² on complex or noisy
  data — but `batch.fit_unified()` / scripting reported the same fit as
  `'success': False` because it passed scipy's internal convergence status
  through directly. scipy's `success` flag means "a tolerance was met by the
  optimizer" — it is independent of chi² magnitude and can be `False` even on a
  fully-converged restart loop. The `fit()` method in `unified.py` now sets
  `success = np.isfinite(chi_squared)` (fit ran to completion and returned a
  valid result), matching the GUI's effective rule. A fit with chi²=5000 /
  reduced chi²=111 is now reported as `success=True` everywhere.

- **Unified Fit and Modeling: fit now converges in a single press.**
  Least-squares fitting previously terminated far short of the true minimum
  and only crept toward the optimum each time the user re-pressed **Fit**
  (often needing several presses) — a problem especially in scripts, where
  Fit cannot be re-pressed. Three combined fixes:
  1. **Parameter scaling.** `scipy.optimize.least_squares` (TRF) was called
     without `x_scale`, so its trust region and convergence tests operated on
     the raw parameter vector. With parameters spanning many orders of
     magnitude (`G` ~ 10⁴, `B` ~ 10⁻¹⁰, `Rg` ~ 10¹, background ~ 10⁻²), no
     single trust-region step could be meaningful for both large and tiny
     parameters at once. Both engines now pass `x_scale='jac'` (auto-rescale
     each parameter by its Jacobian-column norm each iteration — the scipy
     equivalent of Igor Pro's per-parameter fit-step / epsilon on
     log-dependent parameters).
  2. **Tight convergence tolerances.** With `x_scale='jac'` the `xtol`/`ftol`
     tests run in *scaled* space, where scipy's loose defaults (and the
     Modeling engine's former `1e-5`) fired a spurious "converged" on the
     first small step while still far from the minimum. Tightened to `1e-12`.
  3. **Internal restart loop.** As a safety net, each engine now re-seeds the
     solver from its own result until χ² stops improving — the automated
     equivalent of pressing Fit a few times — so scripts get the fully-settled
     result on the first call (typically 1–3 restarts, ~50–70 evaluations).
     In the Modeling engine every restart calls `_residuals`, so the GUI's
     "Cancel Fit" stays responsive across the loop.

  A single **Fit** press (or one scripted `fit()` call) now reaches the
  minimum; a genuinely bad starting point in the wrong basin still needs the
  existing **global (differential-evolution)** fit option. Affects
  `pyirena/core/unified.py` (standalone Unified Fit tool) and
  `pyirena/core/modeling.py` (Modeling tool's unified-level and other
  local-fit populations).

## [0.9.7] — 2026-07-02

### Added

- **Automatic text-file import and cleaning.** ASCII SAS files (`.dat`,
  `.txt`) are now automatically cleaned and converted to a full NXcanSAS
  HDF5 sibling (`<stem>.h5`) on first use. All fitting tools, result saving,
  and viewers then work on the HDF5 file — text-file awareness is confined to
  a single import layer (`pyirena/io/text_import.py`).
  - **Cleaning rules (silent, recorded in HDF5 provenance):** points with
    `Q ≤ 0` removed (occasional Q=0 beamstop rows); points with `I ≤ 0`
    removed (beamstop zeros invisible on log-scale but fatal for numerical
    fits); surviving points with missing or zero uncertainty have `E` replaced
    by `I × error_fraction` (default 5%, configurable in Data Selector →
    Configure).
  - **Naming and caching:** converted file is placed next to the original
    as `mydata.dat → mydata.h5`. Reused on subsequent calls via mtime cache.
    Collision guard: if `<stem>.h5` already exists and was not created by
    pyirena, falls back silently to `<stem>_NX.h5`.
  - **All consumers updated:** Data Selector (plotting, all 6 tool launchers,
    report, ASCII export), `pyirena.batch` (headless fits also get cleaned
    data and always save results to a valid file), `plot_saxs.py`, Data Merge,
    Data Manipulation panels.
  - **Documentation:** new `docs/data_import_and_cleaning.md` covering the
    workflow, cleaning rules, naming convention, batch API, and low-level API.
  - **22 unit tests** covering `clean_sas_arrays`, `ensure_nxcansas_sibling`,
    mtime cache, collision guard, and regression that the produced file
    contains both `sasdata` and can receive fit results.

- **Shared in-panel data loader — all 6 tool panels.** Every tool panel
  (Unified Fit, Sizes, Simple Fits, WAXS Peak Fit, Modeling, SAXS Morph) now
  has a uniform `Open…` button at the top of its left panel so data can be
  loaded directly within the tool, without going through the Data Selector.
  - New module `pyirena/gui/data_loading.py` provides: `load_data_file`
    (text → clean/convert → HDF5, with multi-dataset picker for HDF5 files),
    `read_nxcansas_with_picker`, `prompt_dataset_choice`, and the
    `DataFileLoaderRow` widget (filename field + `Open…` button, emits
    `data_loaded` signal).
  - Dialog filter includes text files (`.dat`, `.txt`) in all tools — text
    files are automatically cleaned and converted; Modeling and SAXS Morph
    previously accepted HDF5 only.
  - `DataSelectorPanel` now calls the shared functions (thin wrappers); no
    behavior change for Data Selector users.
  - Last-used folder is shared across all tools via the `data_selector/
    last_folder` state key.

### Changed

- **Sizes: "Fit All" button renamed and improved.** Button renamed to
  "Fit Cmplx. Bckg. & Sizes" (wider, 180px) to clearly reflect that it runs
  power-law fit → background fit → size distribution fit sequentially on the
  **current loaded data**, not on all selected files (use the Data Selector's
  "Size distribution (script)" button for batch fitting). Updated tooltip and
  docstring to match.

- **GUI layout consistency across all tool panels.** The top-of-panel layout
  is now uniform: bold tool title on the left, red `? Help` button on the
  right, data-file loader row below.
  - WAXS Peak Fit: "No limits?" checkbox moved from title row to the right
    end of the Q-range display row.
  - Modeling: "No limits?" moved from its own row into the Q-range row;
    label shortened (tooltip carries the explanation).
  - Unified Fit: block-style title (with background fill) replaced by plain
    bold title + Help button on right; "Identify Features" button enlarged to
    26 px / 12 pt for readability.
  - Simple Fits: "Simple Fits" bold title row with Help button added above
    the data loader; Help removed from the Model selector row.
  - Sizes: title background block removed to match the other panels;
    "Identify Features" button enlarged to 26 px / 12 pt (matching Unified).

## [0.9.5] — 2026-06-30

### Added

- **Sizes: "Identify Features…" button.** The Size Distribution panel now has
  the same Feature Identifier as the Unified Fit panel (top-right, next to
  Help). It segments the I(Q) log-log slope profile, overlays the segments and
  Guinier knees on the graph, and — below the segment list — shows the
  **size-distribution recommendation** (suggested radius grid, inversion
  Q-range, low-Q power-law and high-Q flat-background windows) with a
  suitability verdict and warnings. The GUI and the AI `suggest_sizes_setup`
  tool share one core function (`pyirena.core.sizes.recommend_sizes_setup`), so
  the displayed recommendation matches what the assistant receives.
  Visualisation only — it never modifies the fit. (`SizesFeatureIdentifierDialog`
  in `pyirena/gui/sizes_feature_identifier.py`; the base
  `FeatureIdentifierDialog` gained overridable title/help/summary hooks.)
- **AI control tools for Size Distribution (MCP).** The `pyirena.api.control`
  surface and the `pyirena-mcp` server now let an AI agent drive a particle
  size-distribution fit end-to-end, alongside the existing Unified Fit control
  tools. New module `pyirena/api/control/sizes.py` exposes 16 functions
  (MCP-prefixed `pyirena_ctrl_sizes_*`): `select_sizes_model` (MaxEnt /
  Regularization / TNNLS / Monte Carlo), `set_size_grid`, `set_shape`
  (sphere/spheroid, contrast, aspect ratio), `set_method`,
  `set_error_handling` (error scaling or fractional errors), the complex
  background workflow (`set_background`, `fit_power_law_background`,
  `fit_flat_background`, `get_background_preview_image`), `run_sizes_fit`,
  `get_sizes_distribution`, `get_sizes_results`, `get_sizes_fit_image`,
  `suggest_sizes_setup` (data-driven suitability check + auto-setup), and
  `save_sizes_fit`. The session lifecycle and Q-range tools are shared with
  Unified Fit (`set_fit_q_range` defines the inversion window);
  `list_available_models` now returns `["unified_fit", "sizes"]`. JSON
  tool-schemas added for all new tools. Monte-Carlo per-bin uncertainty is not
  exposed yet. Verified to produce results identical to `batch.fit_sizes` for
  the same settings. Docs: see
  [ai_tools_reference.md](docs/ai_tools_reference.md) and
  [ai_integration.md](docs/ai_integration.md).
- **Sizes: fractional error option.** A new "Fractional error" checkbox in the
  Error Scaling box lets you ignore the uncertainties from the data file and
  generate them as `error = |I| × fraction` (default `0.03` = 3%). Useful when
  collected uncertainties are unreliable (e.g. after merging subsets, or when
  error estimation failed). Mutually exclusive with the error-scale field; the
  setting is persisted in the GUI state, honored by the batch API, and written
  to / read from the NXcanSAS file (`fractional_error`,
  `fractional_error_value`).
- **Sizes: surface-area distribution.** The volume distribution is now also
  converted to a surface-area distribution `S(r) = sv(r)·P_V(r)` and its
  cumulative (running integral), saved into the NXcanSAS file as
  `surface_dist` and `cumul_surf_dist`, plus the total **specific surface
  area** (`specific_surface`, [Å⁻¹]) as a browseable scalar. The surface-to-
  volume ratio is shape-aware: `3/r` for spheres and a closed-form `C(AR)/r`
  for spheroids (`pyirena.core.sizes.surface_distribution`).
- **Data Explorer: surface-area presets.** New "Size Dist. surf. S(r)" and
  "Size Dist. cumul. surf." checkboxes in the 1D Graph tab, and
  `specific_surface` added to the Size Distribution "Collect Values" items.
- **Data Explorer: curve offsets (waterfall / stacked view).** A new **Offset…**
  toolbar button and **Offset curves…** right-click action open a dialog that
  separates overlapping curves. Each axis has an independent offset type —
  additive (`v + off`, for linear axes) or multiplicative (`v × off`, a constant
  visual shift on log axes) — defaulting to match the current axis scale.
  Includes auto-stagger (curve *i* → `i × inc` additive, `inc`ⁱ multiplicative)
  and a per-curve fine-tuning table. Offsets are display-only: the raw data is
  never modified, so CSV/HDF5/ITX/matplotlib exports stay pristine.

### Changed

- **Data Explorer: one graph per quantity.** Selecting several distribution
  presets no longer dumps differential and cumulative curves onto a single Y
  axis. Each quantity (I(Q), vol P(r), num N(r), surf S(r), and the three
  cumulatives) now opens its own graph with the correct axis labels and title;
  the *same* quantity from multiple files is grouped into one graph so files
  can be compared. "Add to active graph" sends only I(Q) curves to the active
  graph; other quantities open their own.

### Fixed

- **Sizes: number distribution now shape-aware.** `number_dist` /
  `cumul_num_dist` previously assumed a sphere volume `(4/3)πr³` even for
  spheroid fits. They now use the true particle volume (`×AR` for spheroids)
  via `pyirena.core.sizes.number_distribution`, making them consistent with
  the new surface-area distribution. Sphere fits are unchanged; older files
  written before this fix need a re-store to update spheroid number
  distributions.
- **Sizes: "Save as Igor Pro ITX" on the distribution graph.** The size
  distribution is drawn as a stepMode bar (a `PlotCurveItem` in log-x space),
  which the ITX exporter skipped — it reported "No named data curves found to
  export." The exporter now honours an explicit `_itx_export` payload, so the
  distribution exports with linear radii and `P(r)` values.
- **Data Explorer: legend now populates.** The legend was created *after* the
  first curve was plotted, so pyqtgraph never registered it and the "Legend"
  button appeared to do nothing. The legend is now rebuilt explicitly from each
  curve's label, identifying both the quantity and the source file.
- **Sizes: "Fit B?" / "Fit P?" checkboxes are now persisted.** Their state is
  saved to the state file and exported config (and restored on restart),
  instead of resetting to unchecked every session. `batch.fit_sizes` now acts
  on these flags, running the power-law (and flat-background) pre-fit before
  the size fit — so a script can fit the background first, then the model,
  matching the GUI "Fit All" sequence.

## [0.9.4] — 2026-06-27

### Added

- **About dialog reads version from `pyirena.__version__`** — no longer a
  manually maintained string that lagged behind releases.
- **Data Merge batch: live progress reporting.** Status bar shows
  `Merging N/total: filename …` during a batch run so the user knows the
  tool is active when merging hundreds of pairs.
- **Data Merge batch: detailed failure summary.** When pairs are skipped or
  fail, a warning dialog lists every affected pair with its specific reason
  (load failure, no valid data points, no Q overlap, optimizer error, etc.)
  instead of a silent stop.
- **Sizes: "Set Q from cursors" button** enlarged and placed inline to the
  right of the Q min / Q max fields (spanning both rows) in both the
  Power-Law and Flat Background sections. "Q range for fit:" label centred.
- **Sizes: number of bins max raised to 501.** 501 bins across 5 decades
  gives exactly 100 bins per decade for log-spaced grids. Entering a round
  multiple of 100 (100, 200, 300, 400, 500) is silently rounded up by 1 to
  maintain the decade-aligned count.

### Fixed

- **Unified Fit: fitting stopped prematurely.** `max_nfev` (maximum function
  evaluations passed to `scipy.least_squares`) raised from 1 000 → 5 000,
  reducing early convergence on complex multi-level fits.
- **Unified Fit GUI: control panel width was fixed.** Size policy changed
  from `Fixed` to `Preferred`; `setMaximumWidth(400)` removed so the
  splitter can be dragged wider.
- **Unified Fit GUI: Copy/Swap button stretched full width.** Now capped at
  120 px (same as Graph Unified), matching its secondary importance.
- **Batch scripting (`fit_unified`, `fit_sizes`, `fit_simple`, `fit_waxs`,
  `fit_modeling`): `_pyirena_config` missing from HDF5 output.** The
  IO writers already accepted `setup_state` but the batch callers never
  passed it. "Load Setup from File" in the GUI therefore failed with
  *"No … setup is stored"* on files produced by scripts. Each batch function
  now builds and passes the setup state; `fit_unified` additionally updates
  parameter values to the fitted result so the GUI starts from the solution.
- **HDF5 `num_levels` / `level_number` attributes unreadable in external
  viewers** (e.g. Igor Pro). Stored as Python `int` (h5py writes these as
  object-typed scalars on some versions); changed to `numpy.int32`.
- **Save Params to JSON — confusing double dialog.** All five panels
  (Unified Fit, Sizes, Simple Fits, WAXS Peak Fit, Modeling) showed the
  macOS native *"Replace file?"* dialog followed by a pyirena *"Overwrite
  section?"* dialog. The OS dialog was misleading (we only update one JSON
  section, not the whole file). Both replaced by a single clear message:
  *"Only the [Tool] section will be updated — all other tool settings in
  this file are preserved."* Implemented via `DontConfirmOverwrite |
  DontUseNativeDialog` on the file picker.
- **Modeling: Save Params to JSON replaced the whole file.** Previously the
  Modeling export wrote `{'_pyirena_config': …, 'modeling': …}` from
  scratch, destroying any Sizes or Unified Fit sections in the same config
  file. Now loads the existing file and updates only the `modeling` key,
  consistent with all other tools.
- **WAXS Peak Fit: Save Params to JSON had no section-exists check.** Would
  silently overwrite an existing `waxs_peakfit` section without asking.
  Consistent check and dialog added.
- **Data Merge batch: silent fail on bad data pairs.** Two root causes: (1)
  `_load_file` showed a blocking modal `QMessageBox` inside the loop, which
  on some Qt/macOS combinations caused the outer method to return early
  without the final summary. Fixed with `quiet=True` batch mode — errors
  printed to console, not shown as dialogs. (2) The per-pair `try/except
  Exception` missed non-`Exception` subclasses (certain C-extension errors).
  Changed to `except BaseException` with explicit re-raise of
  `KeyboardInterrupt` / `SystemExit`.
- **Data Merge batch: Q overlap check silently passed for all-NaN arrays.**
  `q.max()` returns `nan` for arrays containing only non-finite values;
  `nan >= anything` is `False`, so the guard never triggered. Now uses
  `q[np.isfinite(q)].max()` to guarantee a finite comparison value.
- **Data Merge: scale fitting bounds too narrow.** Occasional hardware
  miscalibration produces intensity ratios outside the old 0.01–100 window,
  causing the optimizer to be clamped at the boundary. Bounds widened to
  0.001–1 000 throughout (initial clip, fallback median, `_wls_bg_scale2`).
- **Igor import (`extract_h5xp_to_nexus`): USAXS data not found in
  Igor-exported h5xp files.** `WAVE_PICKERS_H5XP["USAXS"]` only knew
  pyirena-produced wave names (`q_<folder>`, `Q`/`R`/`S`). Igor's own h5xp
  export retains the original USAXS pipeline names (`DSM_Qvec`/`DSM_Int`/
  `DSM_Error` for desmeared data, `SMR_*` for slit-smeared). Both naming
  conventions added, matching the existing pxp path. `R_Qvec`/`R_Int`/
  `R_Error` also added for SAXS/WAXS h5xp to keep parity.
- **Data Merge panel: folder paths not restored after restart.** State was
  correctly written on close but `launch_data_merge` unconditionally called
  `set_folder(1, current_folder)` on every open, overwriting the restored
  DS1 path. DS1 pre-populate now only runs when no saved state exists.
  `hideEvent` added so Cmd+W on macOS also triggers `save_state`.
- **Sizes: power-law fit error not visible in control panel.** When the user
  clicked "Fit P/B" without checking "Fit B?" or "Fit P?", the error
  appeared only in the graph window's status area. Now also shown in the
  control-panel status label in bold red.
- **HDF5 Viewer — Collect tab: Level/Peak selector too narrow.** Spinner was
  fixed at 60 px; text content clipped. Now uses `setMinimumWidth(160)` to
  match the Item combo directly above it.

## [0.9.3] — 2026-06-26

### Added

- **Modeling: Global fit (Differential Evolution → local polish).** New "Fit"
  method selector (right of the Background field) offering **Standard (local)**
  — the unchanged default — and **Global (DE→local)**. The global method runs
  `scipy.optimize.differential_evolution` to locate the correct basin of a
  multimodal χ² surface, then polishes with TRF least-squares. Intended for
  monodisperse **core-shell** and **core-shell-shell** spheres whose sharp Bessel
  oscillations trap local fitters in the wrong minimum. Parameters spanning many
  decades are searched in log₁₀ space internally so the global search samples
  small and large values evenly. Global requires finite bounds (disabled in
  "No limits?" mode); pairs naturally with **Fix limits?**. Cancellation works
  during the global stage; Monte-Carlo uncertainty always uses the fast local
  refinement. Threaded through GUI, session state (schema 2→3), NXcanSAS
  setup save/load, and the `fit_modeling` batch API (`"fit_method": "global"`).
- **Modeling: Parallel global fits (`cores`).** Spinbox beside the method
  selector sets worker processes for the Global (DE) search (`de_workers`;
  default 1 = serial). Higher values evaluate the DE population in parallel —
  e.g. a core-shell global fit drops from ~60 s to ~17 s on 6 cores — with
  identical results. Pins workers to single-threaded BLAS; cancellation via
  per-generation callback; automatic serial fallback on any multiprocessing
  failure. Threaded through GUI, session state (schema 3→4), JSON export, and
  the `fit_modeling` batch API (`"de_workers"` key).
- **Modeling: Core-shell-shell sphere form factor** (`css_sphere_by_core`) —
  distribution over core radius; both shell thicknesses are fixed parameters.
- **Modeling panel improvements:** Autoupdate (150 ms debounce, off by default);
  Show individual population curve; Fix limits? button; Background row moved
  below Population tabs; 4-slot equal-width button layout — all mirroring the
  Unified Fit panel.

### Fixed

- **Modeling: JSON export (`Save params to JSON`) dropped `fit_method`**, so
  headless `fit_modeling` batch runs silently fell back to Standard even when
  Global was selected in the GUI.
- **Modeling: Create Report and Tabulate Results CSV dropped form-factor and
  structure-factor parameters** (SLDs, shell thicknesses, etc.) — they were
  stored in the HDF5 file but missing from both text outputs. Both now enumerate
  dist/ff/sf parameters dynamically so new form factors are never silently omitted.

## [0.9.2] — 2026-06-21

### Added

- **Unified Fit level display reworked** for clarity; Save State buttons removed
  (session state is always auto-saved).

### Fixed

- **Batch scripting dropped GUI fit/link flags** across Unified Fit, Sizes, and
  WAXS Peak Fit (`fit_method`, linked-Rg flags, etc.).
- Two long-standing test failures cleared.

## [0.9.1]

### Fixed

- **Data Explorer and Simple Fits crash on Python 3.9** with `TypeError` on
  startup. `hdf5viewer/__init__.py` and `simple_fits_panel.py` used `X | Y`
  union type annotations (Python 3.10+ syntax) without
  `from __future__ import annotations`. Fixed by adding the future import.

## [0.9.0]

### Fixed

- **Unified Fit: SAS Background field shows unparseable string on first launch.**
  On a fresh install, loading state set the Background field via `eng_fmt()`,
  which formats small values like `1e-6` as `'1×10^-6'` (display-only notation).
  `float()` cannot parse this string, causing a `ValueError` on every Graph Model
  or Fit action. Fixed by using `eng_fmt_edit()` (parseable `e`-notation) for all
  three places that write to the editable background field: state restore, post-fit
  result display, and undo/restore.

- **Data Selector squashed buttons on high-DPI / scaled Windows displays.**
  First-time users on 1080p screens with 125–150% Windows display scaling
  reported the right-column buttons (View, Analysis Tools, Data Processing)
  appearing vertically compressed and unreadable — several wrote in believing
  the application was corrupted. The fix was to manually drag the window
  taller; not at all obvious.

  Root cause: the right column needed ~514 logical pixels of button content
  but the window's 600 px minimum height left only ~86 px for chrome
  (title row, folder row, type/sort row, menu bar). On 150% scaling the
  usable screen height drops below the layout's natural requirement, and Qt
  silently compresses the buttons rather than clipping them — producing
  the squashed look.

  Three coordinated changes resolve this:
  1. **Right column wrapped in a `QScrollArea`** — a vertical scrollbar
     appears only when the column doesn't fit; on large displays it is
     invisible and the layout looks identical to before.
  2. **File Type / Sort row moved into the left column** above the Filter
     row. The right column now starts one row higher, gaining ~26 logical
     pixels of usable height before scrolling kicks in. The combos still
     refresh the file list the same way.
  3. **Title row shrunk by ~25%** (font-size 24→18 px, padding 10→4 px),
     reclaiming another ~18 px of vertical space.

  Window minimum height lowered from 600 to 500 px to match — the file
  list's own 400 px minimum still dominates. No other tool windows were
  changed.

## [0.8.5] — 2026-06-17

### Changed

- **Feature Identifier algorithm rewritten as change-point detection.**
  The v0.8.4 variance-based stability check had two failure modes flagged
  by user testing:
  - **Sample15**: a 0.4-decade Guinier plateau at the very low-Q end was
    invisible because the stability window straddled the transition zone,
    leaving only 0.04-decade "stable" runs that the min-width filter
    discarded.
  - **Sample25**: three distinct power-law regions (P ≈ 2.4, 2.8, 4.1)
    connected by smooth slope drifts over ~2 decades were lumped into a
    single segment with average slope ≈ −3.1, because the variance check
    cannot distinguish "constant slope" from "slope changing slowly".

  The root cause is the same: variance-based stability answers the wrong
  question.  The right question is "does the mean slope here differ from
  the mean slope nearby?" — change-point detection.

  New algorithm:
  - **Change-point statistic**: at each candidate boundary point, compute
    `|mean_left − mean_right|` of the slope profile over a configurable
    window.  Local maxima exceeding a threshold are change-points;
    segments are the intervals between them.
  - **Two-pass refinement**: a loose first pass finds major boundaries;
    a tighter second pass re-scans any segment wider than
    `wide_region_decades` (default 1.0) to detect hidden sub-structure,
    catching the sample25 smooth-drift case.
  - **Edge-aware width filter**: segments touching the data extremes use
    a looser `edge_min_segment_decades` (default 0.05) than interior
    segments (`min_segment_decades` default 0.10), preserving narrow
    low-Q Guinier plateaus and high-Q backgrounds.

  Validated against the 31 hand-labelled ground-truth samples in
  `testData/StructureIdentificationExamples/`: **100% of human-marked
  PLS / GP / Background ranges are now matched** by a detected segment
  (was 94.8% in v0.8.4 with looser thresholds).

  Config schema changed:
  - Removed: `stability_window`, `stability_std_max`,
    `merge_max_gap_decades`.
  - Added: `change_window_1`, `change_threshold_1`, `change_window_2`,
    `change_threshold_2`, `wide_region_decades`,
    `edge_min_segment_decades`.
  - Other classification / merge / knee thresholds retained with updated
    defaults tuned to the new algorithm.

  GUI dialog state version bumped to 2; old saved state under v0.8.4
  field names is silently ignored on first open (defaults are restored).

  Files: `pyirena/core/feature_detect.py`,
  `pyirena/api/control/unified_fit.py`,
  `pyirena/mcp/server.py` (unchanged signature),
  `pyirena/gui/feature_identifier.py` (advanced-params field list +
  schema version bump),
  `pyirena/tests/test_feature_detect.py` (25 tests, including parametrised
  ground-truth fidelity tests for sample15 and sample25),
  `docs/feature_identifier.md` (rewritten parameter explanations + new
  algorithm section + tuning guidance).

## [0.8.4] — 2026-06-17

### Changed

- **Feature Identifier rewritten as power-law segmentation.**  The v0.8.3
  threshold-based plateau detector was insufficient — it produced one or
  zero features on many real SAXS curves because steep Porod slopes
  (P ≈ 2-4) fell between the SAXS plateau threshold and the USAXS
  power-law threshold and were silently ignored.  The new algorithm:
  - Segments the **entire I(Q) curve** into contiguous power-law-slope
    regions using a sliding-window stability test on the smoothed
    d(log I)/d(log Q) profile, with adjacent-segment merging by
    mean-slope similarity.
  - Classifies each segment as ``background`` (small |slope| AND
    touches the high-Q end), ``guinier_plateau`` (small |slope|
    elsewhere), or ``power_law``.
  - Derives ``guinier_knees`` between adjacent segments whose slopes
    differ by ≥ 0.5.
  - Clips data to Q ≤ 0.6 Å⁻¹ by default (SAS-approximation limit;
    user-overridable).
  - Validated against 31 hand-labelled SAXS samples in
    ``testData/StructureIdentificationExamples/``: 95% of human-marked
    Power-Law-Slope / Guinier-Plateau / Background ranges are matched
    by a detected segment.
  - Result schema changed: ``segments`` + ``guinier_knees`` replace the
    v0.8.3 ``plateaus`` / ``peaks`` / ``power_law_regions`` lists.
    No SAXS/USAXS presets — the algorithm is range-independent.
  - Updated: ``pyirena.core.feature_detect``,
    ``pyirena.api.control.unified_fit.detect_features``,
    ``pyirena.mcp.server.pyirena_ctrl_detect_features``,
    ``pyirena.gui.feature_identifier`` (drops preset combo; renders
    one overlay per segment plus knee markers).
  - 9 of 21 unit tests are parametrised ground-truth fidelity checks.

## [0.8.3] — 2026-06-13

### Added

- **Feature Identifier add-on for the Unified Fit panel and AI agent.**
  Both interactive users and the LLM agent driving the Unified Fit tool now
  have a way to ask "what does this curve actually contain?" before choosing
  the number of levels and Q-windows.
  - **Core** (`pyirena/core/feature_detect.py`): slope-profile detector that
    operates on `d(log I)/d(log Q)` in log(Q) space, with reflection +
    linear-extrapolation boundary handling to avoid edge bias. Classifies
    regions as Guinier plateaus, structure-factor peaks, or power-law
    sections, and proposes initial Guinier Q-windows per detected feature.
    All thresholds are expressed in log decades so behaviour is independent
    of point count. Two presets — `saxs_preset()` (≤2 decades) and
    `usaxs_preset()` (>2.5 decades, relaxed slope thresholds) — with
    `auto(q)` picking by data extent.
  - **GUI** (`pyirena/gui/feature_identifier.py`): non-modal
    `FeatureIdentifierDialog` opened from a new "Identify Features…" button
    in the Unified Fit top control row. Draws plateau / power-law regions
    as semi-transparent overlays and peaks as vertical lines on the main
    I(Q) graph. Visualisation only — never modifies level parameters.
  - **MCP / AI agent** (`pyirena/api/control/unified_fit.py`,
    `pyirena/mcp/server.py`): new `detect_features(session_id, preset, …)`
    tool with the same JSON output the dialog renders, so the agent can
    base level-count and Q-window decisions on the curve's slope profile
    instead of guessing.
  - 11 new unit tests in `pyirena/tests/test_feature_detect.py`.

### Fixed

- **`set_parameter_value` now auto-expands bounds when the assigned value
  falls outside them** (`pyirena/api/control/unified_fit.py`). Previously,
  setting G or Rg to a sentinel value (e.g. G=0 to "remove" a level) would
  leave the bounds unchanged, causing scipy to raise "Initial guess is outside
  of provided bounds" at `run_fit` time — even when that parameter was fixed
  and not actually passed to the optimizer. The fix silently widens the lower
  or upper bound to contain the new value, so AI agents can set extreme
  sentinel values without needing to also call `set_parameter_bounds`.

## [0.8.2] — 2026-05-25

### Fixed

- **`.pxp` importer: many small fixes from end-to-end testing on real
  legacy USAXS experiments**. Aggregating one release because each
  fix on its own was small and they were diagnosed in a single session:

  1. **macOS native-alignment bug** (`fix 61333a3`). The igor2 packed
     record-header struct was being built with native alignment on
     macOS Python builds, inflating the 8-byte logical header to 16
     bytes and immediately misreading `numDataBytes` as a giant
     garbage int. Force explicit `<`/`>` byte order to disable
     alignment padding — header is now reliably 8 bytes everywhere.
  2. **OOM-safety on malformed records** (`fix 3cfd756`). Added a
     256 MB sanity ceiling on `numDataBytes` plus a file-size check,
     so a corrupt or unhandled record type can't make us preallocate
     gigabytes of garbage.
  3. **Igor 8/10 v7 wave format** (`fix 63a38a5`). `igor2 0.5.x`
     rejects binary-wave version 7 (Igor 8's long-name format) even
     though its numeric data layout is identical to v5. The loader
     now patches the leading version field from 7→5 and retries —
     recovering hundreds of waves per file in real experiments.
  4. **Per-sample (root-level) folder layout** (`fix 35eec97`). Some
     experiments organise data by sample rather than by technique —
     one folder per sample at root level, with USAXS/SAXS/WAXS waves
     all inside the same folder. The walker now recognises this and
     yields up to three output files per such folder, one per
     technique whose wave triple is present. Infrastructure folders
     (`Packages/`, `SavedSampleSets/`) are skipped explicitly so
     they don't pollute the summary.
  5. **Igor 8 long-name records that we can't decode are now
     surfaced explicitly**. Igor Pro 8 introduced new packed-record
     types (26 and 33) used when folder or wave names exceed 31
     characters. `igor2 0.5.x` skips these as unknown, so any
     sample whose folder name is > 31 chars is silently invisible
     to the importer. The loader now counts these markers and
     reports `n_igor8_longname_markers` on `ExtractionResult`; the
     CLI prints a `*** WARNING ***` block, the batch API surfaces
     it via the result dict, and the GUI shows a Warning-icon
     dialog with the recommendation to re-save the experiment as
     `.h5xp` (which has no such limit and which pyirena reads
     perfectly). See `docs/igor_pxp_import.md` for details.

  Net effect: a real 22 MB legacy experiment went from "crashes with
  MemoryError" to "imports cleanly, with a clear warning if any
  samples are missing due to Igor 8 long names". A 125 MB
  time-series experiment with ~700 samples went from "imports 22
  samples" to "imports every sample whose folder name fits, with
  warning that the rest need .h5xp re-save".

  +4 new unit tests (23 total).

## [0.8.1] — 2026-05-25

### Added

- **Import Igor Pro `.h5xp` packed experiments** (Wavemetrics' HDF5
  packed-experiment format), completing the Igor import story. The
  GUI button, batch API, and CLI now all accept `.h5xp` alongside
  the existing `.pxp`. Format is auto-detected from the file extension.

  New entry points:
  - `pyirena.batch.igor_to_nexus("file.pxp"|"file.h5xp", …)` — single
    function for both formats; the old `pxp_to_nexus` is kept as a
    deprecated alias for back-compat.
  - `pyirena.io.pxp_to_nexus.extract_igor_experiment(...)` — dispatcher
    used by the GUI and CLI; also exposes the format-specific
    `extract_h5xp_to_nexus()` for explicit calls.
  - GUI file dialog filter now includes `*.h5xp`.

  Implementation notes:
  - h5xp tree is walked with `h5py.File(...).visit()` over the
    `/Packed Data/` subtree; no special parser needed (HDF5 is
    well-defined). The `/Packed Data/Results/` group is intentionally
    skipped — its per-parameter collected-value waves don't match the
    per-sample-folder shape this importer expects.
  - h5xp wave notes use `key:value;` (colon) while pxp uses
    `key=value;` (equals); the parser auto-detects the separator
    per-note so the same metadata extraction code handles both.
  - h5xp wave names use either the literal triple `Q`/`R`/`S`/`dQ` or
    the suffixed `q_<folder>`/`r_<folder>`/`s_<folder>`/`dq_<folder>`
    convention emitted by `pyirena.io.h5xp_writer.write_iq_data`. Both
    patterns are recognised via the new `<folder>` substitution token in
    `WAVE_PICKERS_H5XP`.
  - 9 new tests use `h5xp_writer` to synthesise fixtures, so they
    don't depend on any external data files and run on every CI build.

  Files: `pyirena/io/pxp_to_nexus.py` (added `_load_h5xp_filesystem`,
  `_H5Wave` adapter, `extract_h5xp_to_nexus`, `extract_igor_experiment`
  dispatcher, `WAVE_PICKERS_H5XP`); `pyirena/batch.py`
  (added `igor_to_nexus`, kept `pxp_to_nexus` as alias);
  `pyirena/gui/data_selector.py` (extended file dialog filter);
  `pyirena/tests/test_pxp_to_nexus.py` (+9 tests).

## [0.8.0] — 2026-05-25

### Added

- **Import Igor Pro packed experiments (`.pxp`) → per-sample NeXus files.**
  Users with legacy Igor data can now bring it into pyIrena for analysis
  without re-reducing from raw detector files. The importer walks the
  experiment's USAXS / SAXS / WAXS folder hierarchy, picks the
  standard wave triples (`DSM_Qvec`/`DSM_Int`/`DSM_Error`(+`DSM_dQ`)
  for USAXS, `R_Qvec`/`R_Int`/`R_Error` for SAXS and WAXS), and
  writes one NXcanSAS `.h5` per sample into a sibling `<pxp>_data/`
  folder organised as `USAXS/`, `SAXS/`, `WAXS/`.

  Rich wave-note metadata produced by the APS USAXS pipeline (sample
  name, thickness, temperature, wavelength, plus the full instrument
  state) is parsed from the `NXSampleStart`/`End`, `NXInstrumentStart`/
  `End`, `NXMetadataStart`/`End` sentinel blocks and lands in the
  canonical NXcanSAS locations (`entry/sample/`, `entry/instrument/`,
  `entry/notes/`). The Q resolution column `Qdev` is written when
  the source wave has `DSM_dQ` (per the NXcanSAS optional-array spec).

  Entry points:
  - GUI: **Data Processing & Reference → Import Igor Experiment…**
    (purple button in the Data Selector). A modal dialog lets users
    pick output folder, technique subset, and overwrite policy; on
    success the output folder loads automatically as the current
    data folder.
  - Headless: `pyirena.batch.pxp_to_nexus("legacy.pxp", techniques=["USAXS"])`
  - CLI: `python -m pyirena.io.pxp_to_nexus legacy.pxp -v`

  New dependency: `igor2 >= 0.5.13` (pure Python, numpy-only).

  Implementation notes:
  - Files: `pyirena/io/pxp_to_nexus.py` (reader + writer engine),
    `pyirena/io/nxcansas_unified.py` (extended `create_nxcansas_file`
    with `dq=` and `metadata=` parameters), `pyirena/batch.py`
    (`pxp_to_nexus` entry point), `pyirena/gui/data_selector.py`
    (button + `_IgorImportDialog`), `pyirena/tests/test_pxp_to_nexus.py`
    (11 tests).
  - The wave-name picker (`WAVE_PICKERS`) and folder-name classifier
    (`TECHNIQUE_FOLDERS`) are data-driven dicts at the top of
    `pxp_to_nexus.py`; users with non-standard pipelines can extend
    them without touching extractor logic.
  - The reader is defensive: a single corrupt wave record (real Igor
    files occasionally have them) is skipped and reported in the
    summary rather than aborting the whole import.


---

Older releases (0.7.2 and earlier) are archived in
[docs/CHANGELOG_archive.md](docs/CHANGELOG_archive.md).
