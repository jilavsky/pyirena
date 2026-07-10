# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
