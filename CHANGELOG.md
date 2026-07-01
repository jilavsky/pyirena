# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- **Sizes: "Fit All" button renamed and improved.** Button renamed to
  "Fit Cmplx. Bckg. & Sizes" (wider, 180px) to clearly reflect that it runs
  power-law fit → background fit → size distribution fit sequentially on the
  **current loaded data**, not on all selected files (use the Data Selector's
  "Size distribution (script)" button for batch fitting). Updated tooltip and
  docstring to match.

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

## [0.7.2] — 2026-05-20

### Fixed

- **`api.tabulate_parameter` silently returned `null` for model-specific
  Simple-Fit parameters.** The schema enumerates only the most common
  Simple-Fits parameters (I0, Rg, B, P, BG_G, BG_P), but actual fits
  often use model-specific names (e.g. Porod `Kp`, Debye-Bueche `Lc`).
  `_extract_scalar` did a strict spec lookup and returned `None` for
  anything not on the list, even though `read_simple_fit()` returned
  the same value correctly. Reported by an AI consumer hitting the
  pyirena MCP server.

  Fix: when no schema spec matches, `_extract_scalar` now falls back
  to probing `params/<leaf>` (with sibling `params_std/<leaf>`),
  `<leaf>` as a top-level dataset, and `<leaf>` as a group attribute,
  where *leaf* = parameter name with optional `param_` prefix stripped.
  Both `tabulate_parameter(tool="simple_fits", parameter="Kp")` and
  `parameter="param_Kp"` now resolve correctly for any model-specific
  parameter, without requiring `TOOL_REGISTRY` updates. A truly absent
  parameter still returns `value=None` rather than raising — matches
  the rest of the api's missing-data convention. New regression tests
  in `pyirena/tests/api/test_aggregate.py`.

## [0.7.1] — 2026-05-19

### Fixed

- **Conda + pip Qt mix breaks GUI on Windows.** `environment.yml` previously
  installed `pyside6` and `pyqtgraph` from conda-forge, but `pyproject.toml`
  (and `pip install pyirena[gui]`) pulls PySide6 from PyPI. Users who ran
  both — e.g. `conda env create -f environment.yml` followed by
  `pip install -e ".[gui]"` to pick up a new dep — ended up with conda's
  PySide6 and pip's PySide6 coexisting. On Windows this produces
  `ImportError: DLL load failed while importing QtCore: The specified
  procedure could not be found` because `shiboken6` from one version
  pairs with `Qt6*.dll` from the other. The launcher's broad
  `except ImportError` then misreports it as "Neither PySide6 nor PyQt6
  found", masking the real cause.

  `environment.yml` now installs Python + the heavy scientific stack
  (numpy/scipy/h5py/matplotlib) via conda and the **entire Qt6/GUI stack
  via pip** (`pip: - -e ".[gui]"`). One installer owns Qt = no mix
  possible. Python is also capped at `<3.14` because Python 3.14 wheels
  for several compiled bindings (PySide6 included) are still flaky.

  Existing broken envs: `conda env remove -n pyirena -y` then
  `conda env create -f environment.yml`.

### Changed

- `docs/installation.md` and `README.md` updated to describe the new
  install architecture, warn against `conda install pyside6` inside the
  pyirena env, and add a dedicated troubleshooting entry for the
  "DLL load failed" error.

## [0.7.0] — 2026-05-19

### Added

- **`pyirena.api` — AI-friendly facade over the HDF5 readers.** A new top-level
  package that wraps the existing `pyirena.io.nxcansas_*` loaders and the
  TOOL_REGISTRY schema in a stable, documented surface intended for non-GUI
  consumers (control agents, automation scripts, LLM tools). Every function
  returns a JSON-serializable dict; numpy arrays are decimated to a bounded
  length and NaN/inf are replaced with None so output is strict-JSON-safe.
  Categories:
  - **Discovery** — `list_files`, `summarize_folder`, `inspect_file`.
  - **Reading** — `read_reduced_data`, `read_metadata`, plus
    `read_<tool>()` for each of the 9 supported analyses (simple_fits,
    unified_fit, size_distribution, modeling, saxs_morph, waxs_peakfit,
    fractals, data_merge, data_manipulation).
  - **Aggregation** — `tabulate_parameter`, `summarize_sample`
    (cross-file parameter trends + per-sample condensation).
  - **Plotting** — headless matplotlib `plot_iq` and `plot_parameter_trend`,
    each returning the saved PNG path and optional base64.
  - **Safety** — optional `PYIRENA_DATA_ROOT` env var restricts all file
    access to a configured subtree; `PYIRENA_MAX_ARRAY_POINTS` (default 500)
    tunes the array decimation cap; `PYIRENA_PLOT_CACHE` selects where PNGs
    are written.
  - **Zero new core dependencies.** Uses stdlib `dataclasses` only.

- **`pyirena.mcp` — Model Context Protocol server.** Optional package
  (`pip install pyirena[mcp]`) that exposes the `pyirena.api` surface to
  any MCP-compatible client — Claude Desktop, Claude Code, AnythingLLM,
  custom AI agents. Stdio transport. New entry point `pyirena-mcp`. Plot
  tools return inline images so the LLM sees what the user sees.

  Example Claude Desktop config:
  ```json
  {"mcpServers": {"pyirena": {"command": "pyirena-mcp",
                              "env": {"PYIRENA_DATA_ROOT": "/data"}}}}
  ```

- **Documentation:**
  [docs/ai_integration.md](docs/ai_integration.md) — install &
  configuration walk-through with Claude Desktop / Claude Code /
  AnythingLLM configs, env-var reference, verification steps,
  troubleshooting table, and security notes.
  [docs/ai_tools_reference.md](docs/ai_tools_reference.md) — tool-by-tool
  reference written for AI agents (usable as system-prompt material) and
  for humans authoring agent prompts: discovery workflow, recommended
  call patterns, example interactions, things to avoid.

## [0.6.5] — 2026-05-19

### Added

- **Modeling tool: three new population types** (Modeling tool grows from 3 to 6 types).
  All three types support fitting with bounds, MC uncertainty, HDF5/NXcanSAS output,
  JSON state persistence, and headless batch use via `fit_modeling()`.
  - **Guinier-Porod Level** (`guinier_porod`) — piecewise analytical model for non-spherical
    or non-dilute scatterers (Hammouda 2010, *J. Appl. Cryst.* **43**, 716–719).
    Parameters: G, Rg1, s1, P, Rg2, s2, RgCO; optional Born-Green correlations (ETA, PACK).
    Default Rg2 = 1e10 gives single-level behaviour; a finite Rg2 activates the low-Q regime.
  - **Mass Fractal** (`mass_fractal`) — fractal aggregate scattering (Teixeira 1988,
    *J. Appl. Cryst.* **21**, 781–785).  Primary particles are spheroids with aspect ratio
    β (default 1 = sphere).  β ≠ 1 uses the orientation-averaged spheroid form factor
    (50-point Gauss-Legendre quadrature over particle orientation), eliminating the
    unphysical Bessel-function oscillations of the monodisperse sphere.
    Parameters: Phi, Radius, Beta, Dv, Ksi, Eta, Contrast.
  - **Surface Fractal** (`surface_fractal`) — scattering from a surface with fractal
    roughness (Teixeira 1988).  Parameters: Surface, Ds (2–3), Ksi, Contrast;
    optional smooth Porod transition to Q⁻⁴ above Qc via error-function blending.

- **Similarity-based outlier rejection in the Data Manipulation Average tab.**
  Detects radiation-damaged frames before averaging using the CorMap test
  (Franke et al., *Nature Methods* 12, 419–422, 2015).
  - Pure NumPy/Python implementation in `pyirena/core/similarity.py` — no
    silx or freesas dependency.
  - Two reference modes: **First frame** (compare each frame vs. frame 0;
    frame 0 is always accepted — standard bioSAXS approach) and **Majority
    vote** (compare each frame vs. the median of all frames; more robust when
    the first frame itself may be damaged).
  - P-value threshold is user-adjustable (default 0.01).
  - Results table with colour-coded rows (green = accepted, red = rejected)
    appears after clicking **Check Similarity**; **Auto-reject** button
    removes outliers from the selection in one click.
  - Similarity settings (method, reference, p-value) are persisted in state.
  - Extensible registry (`SIMILARITY_METHODS` dict): add a new method by
    registering a function — the GUI dropdown updates automatically.

- **`average_data()` headless API gains similarity filtering.**
  New keyword arguments: `similarity_check=False`, `similarity_p_min=0.01`,
  `similarity_method='cormap'`, `similarity_reference='first'`.  When
  `similarity_check=True`, frames failing the p-value test are discarded
  before averaging and reported in the return dict (`'rejected'` list).

- **Per-peak Area reported in the WAXS Peak Fit tool.**  Integral under each
  peak profile (∫I_peak(q)dq) computed in closed form for all four peak shapes
  (Gauss, Lorentz, Pseudo-Voigt, LogNormal); linearised 1-σ uncertainty propagated
  from fitted-parameter uncertainties.
  - Saved as `area` + `area_std` in HDF5; older files recomputed transparently on load.
  - Included in Data Selector report, Tabulate Results CSV, and ASCII export.
  - **Data Explorer → Collect Values** pulldown lists `Area` and `Area_err` for WAXS Peak Fit.
  - **Data Explorer → Export to Igor** writes a `peak_Area_P{n}` wave for every peak.

- **d-spacing top axis (d = 2π/Q [Å]) on all WAXS linear-scale plots.**
  New `DSpacingAxisItem` in `sas_plot.py`; `make_sas_plot()` gains a `d_spacing_axis=True`
  keyword.  Active on the Data Merge WAXS plot and the WAXS Peak Fit main plot.

- **Draggable power-law slope guide lines in the Unified Fit I(Q) graph.**
  Right-click the graph → *Add slope guide line* submenu.  Presets: n = −2, −3,
  −4 (Porod), −5, −6; *Custom slope…* accepts any value in [−8, −0.5].
  Lines carry a label, cycle through a colour palette, drag vertically, and can be
  individually removed or cleared in bulk.  Implemented as `SlopeLine` in `sas_plot.py`.

- **Clearer Save/Load parameter button labels across all tools.**
  "Export Parameters" → **"Save params to JSON"**; "Import Parameters" →
  **"Load params from JSON"** in Simple Fits, Sizes, Unified Fit, Modeling, and
  WAXS Peak Fit.

- **Tooltips added to all main action buttons across all GUI panels.**
  Covers Unified Fit, Simple Fits, Sizes, Modeling, WAXS Peak Fit, Data Merge,
  Data Manipulation, and Data Selector.

### Fixed

- **Modeling: "No limits?" checkbox now hides Min/Max columns for Guinier-Porod,
  Mass Fractal, and Surface Fractal populations.**  `set_no_limits()` was missing
  the five new row-store dicts; the fitting engine already applied unconstrained
  Nelder-Mead for all types — this was a GUI-only gap.

- **Unified Fit top axis now shows R = π/Q numerical labels.**
  `setStyle(showValues=False)` on the `RadiusAxisItem` was removed; behaviour
  now matches Simple Fits.

- **Unified Fit: B limits correctly hidden when "No limits?" is checked after
  unchecking "Estimate B from G/Rg/P".**  `UnifiedFitLevelWidget` now tracks
  the limits-visibility state in `_limits_visible` and consults it in
  `_on_estimate_b_changed`.

## [0.6.4]

### Added

- **Igor Pro h5xp export pipeline (`pyirena/io/h5xp_extractor.py`, `h5xp_writer.py`).**
  `batch_extract_to_h5xp()` reads one or more pyirena NXcanSAS HDF5 files and
  writes a single Igor Pro packed-experiment file (`.h5xp`) containing:
  - I(Q) data waves (`Q_0`, `R_0`, `S_0`, `dQ_0`) with rich wave-note metadata
    (sample name, energy, wavelength, temperature, time, pressure from filename tokens).
  - Analysis result waves for all supported tools (Unified Fit, Size Distribution,
    WAXS Peak Fit, Simple Fits, Modeling, SAXS Morph, Fractals) with canonical
    Igor wave names and `_0` run-index suffix on both Y and X waves.
  - Per-tool scalar Results table (one column per file) for Rg, G, B, P, χ², etc.
  - Filename-derived metadata waves (Temperature_C, Time_min, Time_sec, Pressure_psi)
    written only when at least one file contains the token.
  - `Packed Notebooks/pyIrena_ExportNotes` plain-text notebook listing all processed
    files, with correct `H5T_CSET_UTF8` encoding on dataset and all attributes so Igor
    Pro recognises and reopens it.
  - `Recreation/Recreation Procedures` macro so Igor reopens the notebook on next load.

- **"Export to Igor" tab in the Data Explorer** (`pyirena/gui/hdf5viewer/export_to_igor_tab.py`).
  New tab inside the Data Explorer window lets users select files, choose which
  analysis tools to include, pick an output path, and run the export in a background
  thread with a live log display.

### Changed

- **Data Selector GUI reorganised into three `QGroupBox` sections.**
  - *View & Export*: checkboxes + four output buttons (Create Graph, Create Report,
    Tabulate Results, Export to ASCII) + Data Explorer on its own half-width row.
    Export to ASCII recoloured blue; Data Explorer recoloured red to distinguish it
    from the checkbox-driven output tools.
  - *Analysis Tools*: all six GUI + script button pairs.
  - *Data Processing & Reference*: Data Merge, Data Manipulation, Scattering
    Contrast Calculator, Fractals.
  - **Manage Config** moved below the file list alongside Configure…

- **"HDF5 Viewer" renamed "Data Explorer"** throughout: button, window title,
  panel header, menu action, docstrings, README, and documentation pages.

## [0.6.3]

### Added

- **`pyirena/io/schema.py` — machine-readable HDF5 schema registry.**
  New module `TOOL_REGISTRY` provides a single importable source of truth
  for every pyirena result group: HDF5 path, `@analysis_type` tag,
  available plot curves (X/Y names, units, plot type), extractable scalar
  quantities (including per-subgroup paths for level_N, pop_N, peak_N),
  and whether the group is stripped by Data Merge / Data Manipulation.
  Includes convenience functions `tool_group_present()`, `available_tools()`,
  `tool_key_for_group()`, `tool_key_for_analysis_type()`, plus lookup
  tables `GROUP_TO_TOOL`, `ANALYSIS_TYPE_TO_TOOL`, `STRIPPED_GROUPS`.
  This registry is the foundation for the upcoming HDF5 Extractor and
  Plotter tools.

- **`HDF5_Structure_Reference.md`: three missing SAXS Morph items documented.**
  Section §3.6 now covers `rg_A` (Rg of scattering structure, Å),
  `q_max_model_A` (Q upper bound for model, Å⁻¹), and the optional
  `morphology_metrics/` sub-group with all 10 Tier A/B porous-media
  descriptors (connectivity, percolation, Euler number, EDT pore-size
  percentiles).

### Fixed

- **Strip-list bug: Simple Fit results now correctly stripped on merge/manipulation.**
  `nxcansas_data_merge.py` and `nxcansas_data_manipulation.py` listed
  `entry/simple_fits_results` (with trailing `s`) in their strip lists,
  but the actual group path is `entry/simple_fit_results`.  The mismatch
  silently caused Simple Fit results to survive Data Merge and Data
  Manipulation operations instead of being stripped.  Both strip lists
  and the corresponding note in `HDF5_Structure_Reference.md` §3.8 are
  now corrected.

## [0.6.2]

### Added

- **Diffraction Lines tab: per-CIF lattice / d-spacing scale (`d×`).**
  New per-row spinbox (0.9000 – 1.1000, default 1.0000) lets the user
  scale every reflection's d-spacing for that CIF independently — useful
  when the real material is an alloy whose lattice parameter differs
  from the tabulated nominal element (e.g. Al-based solid-solution
  whose `a` is 0.5 % larger than pure Al → set `d× = 1.0050`).  All
  peak Q's are divided by the factor (Q = 2π/d).  **Exact for cubic
  systems** (Al, Cu, Ni, FCC/BCC steels …); first-order accurate for
  non-cubic.  Per-CIF, so multiple phases on the same graph each get
  their own correction.  Composes with the existing global ΔL distance
  correction (d× applied first, ΔL second).  Right-click row →
  **Reset d× to 1.0**.  Persisted in state (`diffraction_lines`
  schema_version 2; old states default to 1.0).

### Fixed

- **WAXS Peak Fitting: Q0 pre-search now respects the per-peak "Fit?"
  flag.**  Previously, even with `Q0 → Fit?` unchecked for a peak, the
  Q0 pre-search step (cross-correlation global shift and per-peak
  brute-force scan) would still move that peak's starting Q0 before
  the scipy fit, so locked Q0s appeared to drift after fitting.  Both
  pre-search steps now skip any peak whose `Q0['fit']` is False —
  locked Q0 starting values are honoured exactly.  Affects both the
  GUI panel and the `fit_waxs_peaks()` batch API.

### Added

- **Diffraction Lines tab: distance correction (ΔL).**  New ΔL (mm) and
  L_cal (mm) controls compensate for sample-to-detector distance
  miscalibration without touching the data — overlay sticks (and HKL
  labels) shift to where the peaks actually appear in mis-calibrated
  reduced data.  Uses the exact non-linear arctan/sin geometry, which is
  important at WAXS angles where a uniform Q-shift or Q-scale would be
  wrong.  L_cal auto-loaded from NXcanSAS metadata when available
  (mirrors the existing wavelength Auto behaviour); manual override
  always available.  Persisted in state.

## [0.6.0] — 2026-05-06

Major feature release.  Two brand-new visualisation tools (Fractals,
SAXS Morph), comprehensive 3-D viewer enhancements, and a morphology-
metrics module.  No breaking changes to existing tools.

### Highlights

- **NEW TOOL: SAXS Morph (3D voxelgram).**  Generates a 3-D voxelgram
  of a two-phase porous structure from experimental I(Q) using the
  Gaussian Random Fields (GRF) method (Berk / Roberts / Levitz),
  computes the model I(Q) from that voxelgram, and shows it alongside
  a 3-D PyVista isosurface and 2-D slice viewer.  Full GUI panel with
  background pre-fits, save/append HDF5, save-config-to-JSON for
  batch use, and Data Selector integration (button + checkbox +
  tabulate + report + standalone viewer).
- **NEW TOOL: Fractals (mass-fractal aggregate visualization).**
  Grows random fractal aggregates by Monte-Carlo random walk on a
  simple cubic lattice (port of Irena's IR3A_*).  Computes fractal
  parameters (Z, dmin, c, df, Rg primary, Rg aggregate) and back-
  calculates I(Q) by two paths: closed-form Beaucage Unified-fit
  intensity (always) and Shape2SAS-style point-cloud Debye sum (on
  demand).  Three growth modes (Grow One, Grow Many queue, Optimizer
  bisection).  GUI-only, no batch entry.
- **3-D viewer enhancements (applies to all 3-D viewers).**  Rich
  right-click menu: View ▶ (XY/XZ/YZ/Isometric, Perspective ↔
  Orthographic), Lighting ▶ (Headlight / 3-point Lightkit, EDL,
  SSAO, cast shadows), explicit bounding box, white background, axis
  units, default emerald-green colour.  VTK shader-error spam
  silenced at module load.
- **NEW MODULE: morphology metrics for the minority phase of binary
  voxelgrams** (connectivity, percolation, Euler χ, pore-size
  percentiles).  Used by SAXS Morph; saved to HDF5 + tabulated +
  reported in Data Selector.

For the detailed per-commit change log of this release, see
git log v0.5.8..v0.6.0 .

### Added

- **saxsMorph: morphology metrics for the minority phase.**  New
  `pyirena.core.morphology` module computes connectivity / topology /
  pore-size descriptors of the binary voxelgram, applied to the
  minority (rarer-by-volume) phase by default.  Tier-A connectivity:
  number of independent clusters, open vs closed porosity, percolation
  flags along X / Y / Z.  Tier-B topology: Euler-Poincaré characteristic
  χ, plus pore-radius percentiles (Q1 / median / Q3) extracted from the
  Euclidean distance transform.  Compute is well under one second at
  N≤256, uses only `scipy.ndimage` + `skimage.measure.euler_number`
  (graceful degradation if skimage is missing — Euler reports 0).
- **saxsMorph: morphology panel** to the right of the existing fit
  result block in the GUI.  Shows the 11 metrics in a compact table
  with a small "Single realisation — depends on RNG seed and voxel
  pitch" footnote so users don't over-interpret single-shot values.
- **saxsMorph results saved + reloaded with morphology metrics.**
  `nxcansas_saxs_morph.save_saxs_morph_results` writes a
  `morphology_metrics` sub-group under the SAXS Morph group in the
  HDF5 file; `load_saxs_morph_results` reconstructs a
  `MorphologyMetrics` dataclass.  Round-trip preserves all 12 fields.
- **Data Selector: Tabulate Results extended for saxsMorph.**  Check
  the "3D saxsMorph" checkbox and click Tabulate to add 19 SM_*
  columns to the table for each selected file with stored saxsMorph
  results: chi², φ, contrast, Rg, S/V, voxel size / box / pitch, plus
  all 11 morphology metrics.  CSV-exportable like the existing
  per-tool tables.
- **Data Selector: Create Report extended for saxsMorph.**  The same
  checkbox now also drives a "## SAXS Morph Results" section in the
  per-file Markdown report, including a "### Morphology of minority
  phase" sub-section with the 11 metrics in a single Markdown table.

### Changed

- **2D slice viewer: full-white panel background + complete bounding
  box.**  After we switched the 2D slice axis labels to black for
  better contrast against the white data area, those labels became
  invisible against pyqtgraph's default dark-grey panel chrome.  The
  ImageView widget background is now also white (set via stylesheet),
  so the entire viewer reads as a clean white canvas with black
  axes and labels everywhere.  Top and right axes are also shown
  now (without tick values) so the plot has a complete black
  bounding box matching the 3D viewer's outline.

- **3D viewer: default isosurface color is now emerald green
  (`#2ecc71`).**  Was dark grey.  Green is the convention used by
  molecular-visualisation tools (PyMOL, VMD) for the default highlight
  / first object, reads well against the white background under any
  of the lighting modes (default headlight, 3-point lightkit, EDL,
  SSAO, shadows), and stays distinct from the data-curve red /
  reference-fit colours used elsewhere in the I(Q) panel.  Users
  who want a different colour can still pick one via the right-click
  menu (`Pick isosurface color…`).

### Fixed

- **3D viewer: silence VTK shader-error spam on macOS.**  After
  enabling SSAO / cast-shadows / occasionally EDL on macOS (where
  Apple has deprecated OpenGL), VTK's modern shader templates miss
  substitutions on Apple's translation layer and the driver rejects
  the post-processing program.  The viewer recovers (falls back to
  the previous valid shader) but the console used to fill with red
  `vtkShaderProgram` / `vtkOpenGLPolyDataMapper` "Could not set
  shader program" / "attempt to add attribute without a bound program"
  walls per render frame.  We now call
  `vtk.vtkObject.GlobalWarningDisplayOff()` once at module import to
  mute VTK's global warning display.  Trade-off: VTK warnings that
  don't bubble through Python exceptions are silenced — acceptable
  for an end-user tool.  Real exceptions still surface normally.

### Added

- **3D viewer: rich right-click menu for visual tuning.**  Applies to
  ALL `Voxel3DViewer` instances — saxsMorph panel, Fractals panel, and
  the standalone `VoxelViewerWindow` opened from Data Selector.  New
  options on the right-click menu:
    - **View ▶** — quick-set standard orientations
      (XY top, XZ front, YZ side, Isometric) and **Perspective ↔
      Orthographic** projection toggle.
    - **Lighting ▶** — switch between **Default headlight** and
      **3-point light kit** (key/fill/back rig — gives surfaces a
      stronger sense of form).  Three independent post-processing
      toggles, each composable on top of any lighting mode:
        - **Eye-dome lighting (EDL)** — darkens silhouette edges so
          features pop visually; massive readability boost for fractal
          aggregates.
        - **SSAO (depth shading in crevices)** — screen-space ambient
          occlusion; darkens concavities to give a true 3-D-cavity
          feel to porous saxsMorph volumes.
        - **Cast shadows** — real cast shadows from positional lights
          (most visible with the 3-point light kit).
  All toggles default to OFF so the initial render matches the
  previous behaviour; users opt in.  Each handler is wrapped in
  `try/except` so an unsupported feature on an older VTK falls back
  silently rather than crashing the viewer.

### Added

- **saxsMorph: "Save config to JSON…" button.**  Persists the current
  saxsMorph parameters (voxel grid, box, input-mode, φ / contrast,
  background pre-fit Q ranges + values, smoothing σ, RNG seed, fit Q
  range) as a `saxs_morph` section inside a `pyirena_config.json`
  file.  This is the file the Data Selector "SAXS Morph (script)"
  button reads to batch-process multiple files with the same
  parameters.  Behaviour matches the canonical Unified-Fit / Modeling
  pattern: existing `_pyirena_config` files have only the saxs_morph
  section replaced (other tools' sections preserved); new files get a
  fresh envelope; non-pyIrena .json files are rejected to prevent
  accidental overwrite.
- **Data Selector ConfigManagerDialog: now lists saxs_morph.**
  Added `'saxs_morph'` to the dialog's `_KNOWN_TOOLS`, so the user
  can prune the saxs_morph section from a `pyirena_config.json` by
  unchecking its checkbox in **Manage Config…** — same as for any
  other batch-capable tool.

### Changed

- **Data Selector "Create Graph" with Fractals / 3D-saxsMorph
  checkboxes now opens a standalone read-only viewer**, not the full
  tool.  Previously checking either box and clicking Create Graph
  launched the full SaxsMorphPanel / FractalsGraphWindow with the
  selected file pre-loaded — useful but heavyweight (the panel re-runs
  the engine).  The new behaviour:
    - Reads the saved voxelgram (saxsMorph) or aggregate positions
      (Fractals, re-voxelised with the chunky-sphere display geometry)
      directly from the NeXus file.
    - Displays it in a new `VoxelViewerWindow` that holds the same
      `Slice2DViewer` + `Voxel3DViewer` pair the parent tools use.
    - When multiple files / aggregates are selected, a dropdown at
      the top lets the user switch between them inside the same window.
  No engine, no fitting, no Calculate button — pure visualisation of
  saved data.  Tabulate, Create Report and Export ASCII intentionally
  ignore these checkboxes (no meaningful per-file table or curve to
  emit for these visualisation tools).

### Added

- **saxsMorph: vertical Q_box marker on the I(Q) plot.**  Sits at
  `Q_box = 2π / box_size_A` and indicates the LOW-Q model-validity
  limit: below this Q the simulation box is smaller than the structure
  the data is asking the model to represent (no Guinier plateau can
  appear because the box truncates the autocorrelation).  The marker
  updates immediately when the user types a new box size.  Together
  with the existing `Q_nyq` (high-Q voxel-resolution marker, now
  augmented by the analytical Porod extension), users see the full
  valid Q range at a glance.
- **saxsMorph: warning before Calculate when Q_min < 3 × Q_box.**
  The simulation box of side L can only represent structures with a
  longest scale of ~L/2; the lowest Q the model can faithfully
  reproduce is therefore ~Q_box = 2π/L.  When the cursor-set Q_min is
  within 3× of Q_box, the panel pops a Yes/No dialog explaining the
  ratio and recommending the user either move the LEFT cursor to
  higher Q or increase the box size.  Calculate proceeds only on Yes.
- **Fractals: red warning when a loaded Unified-fit gives c < 1.**
  Connectivity dimension `c` must be ≥ 1 for any branched / fractal
  aggregate (chains have c=1; branched fractals c≈1.1–1.4).  When the
  Igor-style derivation from a NeXus Unified-fit pair yields c < 1
  the data simply cannot be modelled as a mass-fractal aggregate
  (Z, dmin, df all become non-physical).  The targets box at the top
  of the panel turns red with a "⚠ NOT a valid fractal" banner, and
  the compact summary above the Active Aggregate Parameters widget
  also flips to red so the warning stays visible no matter where the
  user has scrolled.  Matches the equivalent Igor-Pro warning.

### Added

- **saxsMorph: model now extends above the voxel-Nyquist Q with a clean
  analytical Porod tail K/Q⁴.**  `voxelgram_to_iq` gains a
  `high_q_mode` parameter with three settings:
    - `'porod'` (new default) — extract Porod prefactor `K = median(I·Q⁴)`
      from the FFT in `[0.5·Q_nyq, 0.9·Q_nyq]` (where the FFT is still
      reliable but Porod scattering dominates), then continue as `K/Q⁴`
      above `0.7·Q_nyq`.  No surface-area calculation needed; the
      transition is continuous with the FFT result.
    - `'truncate'` — old behaviour (zero above `Q_nyq`, leaving the
      user-supplied Porod + flat background to take over).
    - `'raw'` — leave the FFT result alone (aliasing-dominated above
      `Q_nyq`); for diagnostic use only.
  The deprecated `truncate_above_nyquist` keyword is honoured for
  back-compat (True → `'truncate'`, False → `'raw'`).
- **saxsMorph: `SaxsMorphResult` gains `porod_K_struct` and
  `specific_surface_area_inv_A`** (the Porod prefactor in I_struct units
  and the derived specific surface area S/V in Å⁻¹).  The GUI result
  table now shows S/V both in Å⁻¹ and cm⁻¹ for a direct physical
  comparison.  The "Q_max_model" row is renamed "Q_nyq (vox.)" with the
  in-plot marker label "Q_nyq (Porod ext. above)" so users see the
  transition rather than expecting a hard cutoff.
- New helpers `extract_porod_prefactor(I_struct, q, Q_nyq)` and
  `specific_surface_area_from_K(K)` exposed at the `pyirena.core.saxs_morph`
  level for callers that want to compute these from a precomputed
  I_struct / Result without re-running the engine.

### Fixed

- **Fractals: spurious high-Q "intensity rises above plateau" with large
  primary particles (e.g. Rg ≥ 100 Å)** caused by histogram bin-width
  vs Q_max coupling.  Each bin contributes `p_bin · sinc(Q · r_centre)`
  to the Debye sum, accurate only while `Q · bin_width ≪ 1`.  With
  default `n_bins = 200` and a 13 000 Å aggregate (Rg = 200 Å), bin
  width was 70 Å → the sinc-at-bin-centre approximation broke down at
  Q ≈ 0.01 and produced large pseudo-random fluctuations (sometimes
  back above the I(0) plateau) at higher Q.  `intensity_montecarlo` now
  picks `n_bins` adaptively from `r_max · Q_max · 10` (clipped to
  [200, 200 000]) so `bin_width ≤ 0.1 / Q_max` everywhere — Q_max ·
  bin_width ≤ 0.1 throughout the requested Q range.  Verified across
  Rg = 10..500 Å: no above-plateau values, plateau holds at ~1.0,
  high-Q tail decays cleanly through Porod into the numerical-noise
  floor.  `MCIntensityWorker.configure(n_bins=None)` (default) selects
  auto mode; pass an explicit integer to override.

### Changed

- **Fractals: `intensity_montecarlo` rewritten as a Shape2SAS-style
  point-cloud Debye sum** (replaces the previous voxel-based random-pair
  Monte-Carlo).  The old method rasterised each primary sphere into
  thousands of voxels then sampled ~10⁷ random pairs from ~10¹⁰
  possible — leaving the high-Q tail dominated by Poisson shot noise
  rather than the expected smooth Q⁻⁴ Porod envelope.  The new method,
  modelled on Andreas H. Larsen's Shape2SAS
  (github.com/andreashlarsen/Shape2SAS):
    1. Generates `n_points_per_sphere` (default 50) uniformly-distributed
       points inside every primary sphere via cube-root inverse-CDF for r
       and Gaussian normalisation for direction.
    2. Computes ALL N(N−1)/2 unique pair distances in a triangular loop,
       histogramming each row's distances on the fly (memory stays in MB,
       not GB) — no random sampling, no shot noise.
    3. Applies the Debye sum  I(Q) = Σ p(r) · sinc(Q·r) / Σ p(r),
       normalised so I(0) = 1.
  New `polydispersity` parameter (default 0.10) gives each primary sphere
  its own R drawn from `N(R_mean, polydispersity · R_mean)` so the
  identical-sphere form-factor zeros decohere — without it, monodisperse
  fringes hide the Porod envelope.  Verified slopes for typical
  Z=150/df=1.9 aggregate:
    - Fractal regime (Q in 0.005–0.05): slope = −2.22 (expected −df ≈ −1.92).
    - Porod region (Q in 0.15–0.4): slope = **−4.58** (expected −4).
  Function signature changed from `(oversample, sphere_voxel_radius,
  max_pairs, time_budget_s)` to `(n_points_per_sphere, n_bins,
  polydispersity, seed)`; `mc_q_max` now reports the typical
  inter-point Q (informational only — `intensity_montecarlo` no longer
  NaN-truncates).  `MCIntensityWorker.configure` updated to match.
- **Fractals: `voxelize` is now used only for 3-D / 2-D display**
  (the Igor-style "fat sphere" rendering).  All MC scattering goes
  through the new point-cloud path; the heavy 80–340 MB voxelgrams the
  scattering used to need are gone.

### Added

- **Fractals: documentation file** [`docs/fractals_gui.md`](docs/fractals_gui.md)
  covering panel layout, workflow, target derivation formulas, growth
  parameter reference, MC vs analytical I(Q), HDF5 layout, and Igor
  compatibility notes.
- **Data Selector: top-row checkboxes for "3D saxsMorph" and "Fractals"**.
  Clicking **Create Graph** with these checked opens the saxsMorph panel
  (loaded with the first selected file's I(Q)) or the Fractals window
  (pre-populated with every aggregate found in the selected files).
  Files without stored results are silently skipped.
- **Fractals: auto-load NeXus reference file from saved state.**  When
  the Fractals tool is opened and a previously-loaded NeXus file path
  is in the application state, the file is now opened automatically
  (data + Unified-fit + targets) without the user needing to re-pick
  it.  If the file has been moved or deleted, the path is silently
  cleared.
- **Fractals: target Z, dmin, c, df derived using Igor formulas.**
  Previously `z_target` was just `G2` (the Guinier prefactor of level 2,
  ~3e6 for typical aggregates) and `dmin`/`c` targets were missing
  altogether.  Now matches Irena's `IR3A_*` exactly:
    `df    = P2`
    `dmin  = B2 · Rg2^P2 / (Γ(P2/2) · G2)`
    `c     = P2 / dmin`
    `Z     = G2 / G1 + 1`
  Verified: a Unified-fit pair that gives Z=40 in Igor now also gives
  Z=40 in pyirena (previously gave 3 000 000).
- **Fractals: target summary duplicated next to the Active Aggregate
  Parameters widget**, so the comparison stays visible without scrolling
  back to the top of the panel.  Each comparable parameter row now also
  shows `(target: X)` inline next to the actual value, colour-coded
  green/orange/red by relative agreement.

### Changed

- **Data Selector: Fractals button moved to share a row with Scattering
  Contrast Calculator** (Support Tools section) instead of taking a
  full row of its own.
- **Data Selector: SAXS Morph buttons renamed** to "3D saxsMorph (GUI)"
  and "3D saxsMorph (script)" for consistency with the new "Fractals"
  3D-tool labelling.
- **saxsMorph 3D viewer: dark grey isosurface on white background**
  (was light blue) so the structure stands out cleanly with visible
  edges.
- **saxsMorph 2D slice viewer: white = void / minority phase, dark
  grey = solid / majority phase** (was black/white inverted), matching
  the 3D viewer convention.  Black axes and tick labels for clear
  contrast with the white canvas.
- **saxsMorph 2D + 3D viewers: axis units now `[A]` instead of `(Å)`**.
  VTK's default font does not include the U+00C5 (Å) glyph, so the 3D
  axes were previously rendering as `X()`.  Plain ASCII `[A]` works
  in every Qt font and VTK build.

### Fixed

- **Fractals: MC I(Q) low-Q Guinier knee was shifted to lower Q vs the
  analytical Unified curve (aggregate Rg approximation bug).**
  `rg_aggregate` was computed with the Alex McGlasson approximation
  `Rg_agg = Rg_primary · Z^((1/c−1)/(dmin−df))`, which systematically
  underestimates the true aggregate Rg by 30–100 % depending on Z and
  df.  As a result, the analytical Unified curve placed its low-Q Guinier
  knee at too-high Q relative to the MC curve (factor of ~1.6 for Z=80,
  larger for smaller Z).  Replaced with the direct measurement from the
  grown particle positions via the parallel-axis (Steiner) theorem:
    `Rg_agg = sqrt(Rg_centers² + Rg_primary²)`
  where `Rg_centers = sqrt( mean( ||p_i − centroid||² ) ) × D` and D is
  the physical lattice spacing (`primary_diameter`).  The residual Q
  offset after the fix is ~10–15 %, which is expected: the MC I(Q) has a
  fractal power-law in the intermediate regime that shifts the half-fall
  Q slightly lower than the Unified model's mathematical Guinier knee.

- **Fractals: Monte-Carlo Guinier knees were a factor of 2 too low in Q
  (radius vs diameter geometry bug).**  In the voxelization, lattice
  spacing 1 = `primary_diameter` and oversample = 10 voxels per lattice
  unit, so 1 voxel = `primary_radius / 5`.  The sphere kernel radius was
  `sphere_voxel_radius=10`, which made the physical sphere radius
  `10 × R/5 = 2R = primary_diameter` — i.e. each particle in the MC
  voxelization was exactly twice the correct primary radius.  This
  inflated the Rg of the solid voxel cloud and pushed both the primary
  and aggregate Guinier knees of the MC curve to half their true Q.
  Corrected default to `sphere_voxel_radius=5` so the physical sphere
  radius equals the primary radius and edge-neighbor particles are
  exactly tangent.  Smoke-test: aggregate-knee Q is now 0.01047 1/Å for
  Rg_agg = 100 Å (expected 1/Rg = 0.01) — was previously ~0.005.

- **Fractals: high-Q MC intensity was unphysical noise.**  Above
  `Q_voxel = π / pitch_A` the discrete voxel grid cannot resolve sphere
  surfaces, so the MC curve there is dominated by quantization noise
  rather than the expected Porod tail.  `intensity_montecarlo` now
  returns NaN for Q > Q_voxel; the panel filters those out and `clip`s
  the curve at Q_voxel cleanly.  New helper `mc_q_max(aggregate)`
  exposes the cutoff for callers.

### Changed

- **Fractals: # points spinbox max raised 2000 → 20000.**  MC scattering
  computation is fast (≤1 s for 20000 points) so a fine Q grid is cheap.
- **Fractals: Monte-Carlo curve color changed from orange to blue
  (`#2980b9`).**  The previous orange was hard to distinguish from the
  red "Unified fit (loaded)" curve.
- **Fractals: voxelization defaults are now `oversample=20`,
  `sphere_voxel_radius=10`.**  The first attempt at the
  radius/diameter fix dropped the sphere kernel from 10→5 voxels with
  oversample=10; that produced the correct physical sphere radius R but
  rendered each sphere with only 5 voxels — flying-edges iso-surfaces
  could barely show the 1-voxel-wide neck where edge-neighbor spheres
  meet, so the 3D view looked disconnected.  The new defaults preserve
  the invariant `sphere_voxel_radius / oversample = 1/2` (so physical R
  is still correct) but use a 2× finer voxel grid, giving a 10-voxel
  kernel that renders smoothly with a clearly visible neck at edge-
  neighbor tangent points.  Costs: 8× more memory for the voxelgram
  (~80 MB at Z=80, ~340 MB at Z=500) and ~8× slower MC sampling
  (still sub-second to a few seconds).  Q_voxel cutoff doubles
  (e.g. 1.21 → 2.43 1/Å for primary_diameter = 26 Å), so the MC curve
  stays valid further into Porod.
- **Fractals: model is now scaled to data over the fractal regime
  Q ∈ [0.5π/Rg_agg, 1.5π/Rg_primary]** (matches Igor's
  `IR3A_Calculate1DIntensity`) instead of the full visible Q range.
  Previously the model sat too high relative to data because the data's
  full integral is inflated by sample-level low-Q power-law and high-Q
  flat background that the single-aggregate model never tries to
  reproduce.  Falls back to the visible-Q overlap when no usable Rg pair
  exists.
- **Fractals: I(Q) plot now shows a legend** ("Data", "Unified fit
  (loaded)", "Aggregate Unified (analytical)", "Aggregate Monte-Carlo").
  The legend is rebuilt on every plot refresh so stale entries from prior
  aggregates do not accumulate.

(Earlier in this Unreleased cycle:)

- **Fractals: Monte-Carlo I(Q) had a spurious low-Q rise.**  The Glatter-
  Kratky transform was applied with an extra `r²` factor.  The pair-distance
  histogram from random voxel-pair sampling already IS the Glatter PDD
  `p(r) = 4π·r²·γ(r)` — the spherical-shell `r²` weighting is built in by
  construction.  Multiplying by another `r²` over-weighted large pair
  separations and distorted the low-Q shape (the Guinier plateau expected
  for a finite particle below `Q ≈ 1/Rg_aggregate` was wrong).  Corrected
  to `I(Q) ∝ ∫ p(r) · sinc(Qr) dr`.
- **Fractals: changing "# points" after Grow made all curves vanish from
  the I(Q) plot.**  When the Monte-Carlo worker returned with a different
  Q-grid length than the previously stored `i_unified`, the bookkeeping
  arrays (`agg.q`, `agg.i_unified`, `agg.i_montecarlo`) ended up with
  mismatched shapes.  The plot's mask `(q > 0) & (I > 0)` then raised a
  silent broadcast error.  `_on_compute_mc` now syncs `agg.q` and
  `agg.i_unified` to the new Q-grid *before* starting the MC worker, and
  `_on_mc_finished` re-evaluates `i_unified` on the returned grid as a
  belt-and-braces measure.  `_refresh_plot` also includes a defensive
  shape check that drops mismatched arrays rather than crashing.

### Added

- **Fractals tool — new Support Tool for mass-fractal aggregate visualization.**
  Grows random fractal aggregates by Monte-Carlo random walk on a simple
  cubic lattice (port of Irena's `IR3A_MakeAgg` / `IR3_3DModels.ipf`),
  computes their fractal parameters (Z, dmin, c, df, R, p, s, Rg primary,
  Rg aggregate, true sticking probability), and back-calculates I(Q) by
  two paths: a fast closed-form two-level Beaucage Unified-fit intensity
  (always computed) and a slow Monte-Carlo PDF-based intensity (on demand).
  - Three growth modes: Grow One, Grow Many (queued batch), and an
    Optimizer ("Find Best Growth") that bisects sticking probability to
    match a target dmin and c.
  - Background QThread queue keeps the UI responsive while aggregates
    grow; users can inspect already-completed aggregates in the in-session
    list while new ones are still growing.
  - Optionally loads a NeXus file containing Unified-fit results: when two
    consecutive levels have `RgCutoff_high ≈ Rg_low`, the pair is treated
    as a fractal representation and target Rg primary, Rg aggregate, and
    df are displayed for visual comparison.
  - GUI mirrors SAXS Morph: scrollable left panel + right vertical
    splitter (top: log-log I(Q), bottom: 2D slice + 3D PyVista isosurface
    of the voxelized sphere structure).  Reuses `Voxel3DViewer`,
    `Slice2DViewer`, `make_popout_button`, and the standard `sas_plot`
    helpers — no duplication.
  - Save selected aggregate to a NeXus file as
    `entry/fractals_results/aggregate_{N}` (NXprocess group with
    positions, neighbor list, computed parameters, input parameters, and
    optional intensity).  Multiple aggregates auto-increment.  Loader
    supports reading them back into the session list.
  - GUI-only — no batch / scripting / headless API (matches the tool's
    visualization-focused workflow).
  - New files: `pyirena/core/fractals.py`,
    `pyirena/io/nxcansas_fractals.py`, `pyirena/gui/fractals_workers.py`,
    `pyirena/gui/fractals_panel.py`.  New `fractals` block in
    `state_manager` (schema_version 1).  New "Fractals" button in
    `data_selector` Support Tools section, immediately below the
    Scattering Contrast Calculator.

## [0.5.8]

### Fixed

- **SAXS Morph: Berk inversion missing — model I(Q) was systematically
  too low by 5–10× in the data Q range.**  In Berk's GRF method, if the
  underlying Gaussian field has autocorrelation `g(r)` and is thresholded
  at level α, the binary indicator's autocorrelation is `T(g, α)` —
  **not** `g` itself, where

      T(g, α) = (1/2π) · ∫₀^g exp(-α²/(1+t)) / √(1−t²) dt

  For φ=0.3 (α≈0.52), `T(g, α) ≈ 0.16·g` at small g, so feeding
  `gamma_data(r)` directly as the field's autocorrelation (the previous
  behaviour) produced an indicator whose autocorrelation was ~6× lower
  than the data, and a model I(Q) correspondingly suppressed.
  - New `berk_lut(alfa, n)` builds a forward LUT of T(g, α) using the
    `t = sin(θ)` substitution to remove the integrable singularity at
    g = ±1; `berk_invert(target_T, alfa, lut)` solves T(g, α) = target_T
    by 1-D interpolation (T is monotonic in g).
  - `compute_voxelgram` now: (1) rescales the normalised γ to its
    physical absolute scale `γ_phys = γ_norm · φ(1−φ)`, (2) inverts
    via `berk_invert` to obtain g(r), (3) computes the spectral function
    F(k) as the Fourier transform of g(r) — not γ(r) — so that the
    field synthesised by FFT-filtering white noise has the correct g,
    and after thresholding gives the desired indicator γ.
  - Also: the spectral-function k-grid is now sized to cover the full
    3D FFT range (`k_max = √3 · π / pitch`), so high-k content is no
    longer silently zeroed during the 3D resampling.
  - Four new tests verify Berk LUT properties: T(0)=0, T(g→1) ≈ φ(1−φ),
    invertibility, and monotonicity.

### Changed

- **SAXS Morph: Result block re-laid out as 2-column HTML table** with
  font bumped from 9pt monospace to 11pt sans.  All 14 result fields fit
  comfortably; the parameters block is much easier to scan.
- **SAXS Morph: 2D and 3D viewers now show the same boundary** when
  smoothing is on.  Previously the 2D slice showed the raw binary
  voxelgram (sharp, voxel-scale aliasing) while the 3D viewer showed a
  Gaussian-smoothed isosurface (much smoother) — the two views looked
  inconsistent.  Now the smoothed scalar field is computed once; the 3D
  viewer renders its isosurface at 0.5 and the 2D slice shows the same
  smoothed field thresholded at 0.5.  Both views display the same
  microstructure boundary as a binary black/white image.

## [0.5.7]

### Changed

- **SAXS Morph: removed "Cube side (fit)" from the GUI** — it had no
  effect. The two-combo design came from an abandoned idea of running an
  optimiser loop at a small voxel size (fit) then rendering at a larger
  one (render). Since there is no iterative fitting in the current
  workflow (Calculate 3D always runs once at the chosen resolution),
  both controls pointed at the same step and the "fit" one was silently
  ignored. The Voxel grid box now has a single "Cube side" combo with
  an improved tooltip listing expected compute times and RAM usage.
  The underlying `SaxsMorphConfig.voxel_size_fit` field is kept for
  backward compatibility with the deprecated `engine.fit()` API and is
  automatically set equal to `voxel_size_render` by the GUI.

## [0.5.6]

### Fixed

- **SAXS Morph: beach ball / frozen GUI during Calculate 3D.** The
  `QTimer.singleShot(0, ...)` deferral still ran on the main thread.
  Replaced with a proper `_CalcWorker(QThread)` that runs
  `compute_voxelgram()` on a background thread, keeping the GUI fully
  responsive (no macOS spinning beach ball).
- **SAXS Morph: model I(Q) systematically low (~2x) when Gaussian
  smoothing was enabled.** The engine was computing I(Q) from the
  *smoothed* float32 voxelgram, which has lower variance than the
  original binary field (smoothing suppresses high-frequency content).
  Fix: I(Q) is always computed from the **binary uint8 voxelgram**
  (correct physics).  Gaussian smoothing (`smooth_sigma`) is applied
  **separately** at display time — only to the 3D PyVista viewer, where
  it softens isosurface aliasing.  The 2D slice viewer and the HDF5
  result always hold the binary cube.  `result.voxelgram` is always
  `uint8` regardless of `smooth_sigma`.
- **SAXS Morph: 2D slice viewer grayscale was confusing.** The slice
  viewer now always displays the binary microstructure (black = void,
  white = solid) with a 2-colour LUT.  Only the 3D isosurface viewer
  uses the smoothed float field.
- **SAXS Morph: Save button label misleading.** Renamed from
  "Save Result to HDF5…" to **"Save/Append to HDF5…"** with a tooltip
  explaining that all other result groups in the HDF5 file are
  preserved and only `entry/saxs_morph_results` is replaced.  In
  batch/script mode (`fit_saxs_morph()`) the save is always automatic
  — no dialog shown.

## [0.5.5]

### Added

- **SAXS Morph: invariant extrapolation outside the data Q range** —
  matches the Igor reference's `IR3T_ExtendDataCalcParams` behaviour
  and removes the residual ~1.6× discrepancy with Igor noted in 0.5.2.
  - Low-Q: `I(Q→0) ≈ <I[0:5]>` (constant), contributing
    `I0 · q_min³ / 3`.
  - High-Q: `I(Q→∞) ≈ K·Q⁻⁴` Porod tail with
    `K = <I·Q⁴ over last 5 points>`, contributing closed-form
    `K / q_max` to the invariant integral.
  - New helper `compute_invariant_extrapolated()`; both
    `derive_contrast_from_invariant()` and
    `derive_phi_from_invariant()` now use it.
- **SAXS Morph: post-threshold Gaussian smoothing of the voxelgram** —
  mimics Igor's `ImageFilter/N=5 gauss3d`, knocking down the per-voxel
  numerical noise that comes out of the GRF realisation.
  - New `SaxsMorphConfig.smooth_sigma` (in voxels; default 1.0).
  - When > 0, voxelgram is `float32` in `[0, 1]` (3D viewer's
    isosurface at 0.5 still works; 2D slice viewer renders as a
    smooth grayscale gradient via a 256-level LUT).
  - When 0, voxelgram stays `uint8` binary (legacy behaviour).
  - Exposed in the GUI as **Smoothing σ [vox]** in the Voxel grid box.
- **SAXS Morph: snappy Calculate 3D button + elapsed-time report.**
  - The yellow "Calculating..." status banner now appears immediately
    on click (was delayed by the synchronous compute starting before
    Qt could repaint).  Fix uses `QTimer.singleShot(0, ...)` to defer
    the heavy compute to the next event-loop tick.
  - Re-entrancy guard (`_calculating` flag) prevents queued OS-level
    clicks from re-triggering the compute after `setEnabled(False)`.
  - Status banner now reports elapsed time, e.g.
    `Calculate 3D done: χ² = 12.34, φ_actual = 0.31, voxel = 256³ [took 18.7 s]`.

## [0.5.4]

### Fixed

- **SAXS Morph: dead 3D viewer on second open.** The 0.5.3 shutdown fix
  nulled `Voxel3DViewer.plotter` on close, but the Data Selector cached
  the panel object — so reopening it gave back the same panel with a
  dead plotter (`AttributeError: 'NoneType' object has no attribute
  'add_mesh'` on the next Calculate 3D).
  - `SaxsMorphPanel` now sets `WA_DeleteOnClose`, so the panel is fully
    destroyed when closed (rather than just hidden) — the embedded
    PyVista QtInteractor goes away with it.
  - `DataSelectorPanel.launch_saxs_morph()` connects to the panel's
    `destroyed` signal to clear `self.saxs_morph_window`, so the next
    launch creates a fresh panel with a live VTK render window.
  - State / parameters persist across recreations because they're
    saved to `~/.pyirena/state.json`.  The data file is reloaded from
    Data Selector's current selection.
- **CHANGELOG correction**: the 0.5.3 note suggesting
  `pip install "PySide6==6.6.*"` is wrong on Python 3.13+ — only
  PySide6 6.10.1+ ships wheels for those Python versions.  The 0.5.3
  shutdown fix already makes the panel safe to close on PySide6 6.10.x;
  the WA_DeleteOnClose addition above closes the loop on reopening.

## [0.5.3]

### Fixed

- **SAXS Morph: macOS segfault when closing the panel.** Closing the
  SAXS Morph window killed the entire pyirena process with
  `EXC_BAD_ACCESS` in `vtkCocoaRenderWindow::Render` (PySide6 6.10.x +
  PyVista on macOS).  Root cause: Qt does not propagate `closeEvent` to
  child widgets, so `Voxel3DViewer.closeEvent` was never called and the
  VTK render window stayed alive long enough to receive a late
  "backing layer changed" signal from AppKit, which then tried to
  `Render()` onto a freed `NSView`.
  - Added `Voxel3DViewer.shutdown()` that explicitly calls
    `render_window.Finalize()` before `plotter.close()` and nulls out
    `self.plotter` so any late signal sees `None` and bails out.
  - Added `SaxsMorphPanel.closeEvent()` that calls `shutdown()` on the
    embedded 3D viewer first, and also stops any in-flight worker
    threads (`_FitWorker`, `_MCWorker`) before allowing Qt to destroy
    the widget tree.
  - `_PopoutDialog.closeEvent()` now calls `shutdown()` on its widget
    if the original layout has been destroyed (e.g. main panel closed
    while the popout was open).
- **Note**: PySide6 6.10.x is still flagged as incompatible in the
  `gui` extra (`!=6.10.*`) but pip will not downgrade an
  already-installed PySide6.  If you continue to see VTK shutdown
  warnings, pin manually: `pip install "PySide6==6.6.*"`.

## [0.5.2]

### Fixed

- **SAXS Morph: critical Å↔cm unit conversion bug.** Two related errors in
  the math gave nonsensical contrast values and grossly inflated model
  intensity:
  - `derive_contrast_from_invariant()` and `derive_phi_from_invariant()`
    integrated `Q²·I·dQ` numerically with Q in Å⁻¹ and I in cm⁻¹, giving
    a result in Å⁻³·cm⁻¹ — but treated it as if it were in cm⁻⁴ (or
    10²⁰ cm⁻⁴). Missing factor: 10⁴. For a typical sample where Igor
    reports Δρ² = 34, pyirena 0.5.1 returned 0.005559 (6116× too low,
    or 1.6× off after this fix once you account for Igor's Q*
    extrapolations beyond the data window).
  - `voxelgram_to_iq()` multiplied the per-bin `<|FFT|²>/N³` by the box
    volume `(N·pitch)³` instead of just `pitch³`, double-counting N³.
    It also omitted the Å³→cm³ conversion (10⁻²⁴) needed when contrast
    is in 10²⁰ cm⁻⁴ and the result must be in cm⁻¹. Combined error:
    `N³ × 10⁴` too high in the structure factor — for N=256, ~10¹¹×
    too high. With contrast self-derived from the (also buggy) invariant,
    the two 10⁴ factors cancelled but the spurious N³ remained, leaving
    the model curve ~N³ ≈ 10⁷–10⁸ above the data.
  - Both formulas now use explicit `× 1e4` and `× 1e-4` factors with
    inline derivations in the docstrings/comments, and a new self-
    consistency test (`test_voxelgram_iq_is_in_per_cm_when_contrast_is_e20cm4`)
    verifies that integrating the model `Q²·I·dQ` reproduces
    `2π²·φ(1−φ)·Δρ²` to within a numerical-aperture factor of ~1.5
    (a stray 10⁴ or N³ would push it off by orders of magnitude).
  - 28 tests pass (was 26).

## [0.5.1]

### Changed

- **SAXS Morph: workflow rework to match the Igor Pro reference.** The
  initial 0.5.0 release treated φ, contrast, B/P, and the flat background
  as iteratively-fittable parameters via least_squares — but in the
  Berk/Roberts/Levitz GRF method these are not free parameters. The
  voxelgram is computed deterministically from φ (or contrast) plus a
  spectral function derived from the data; the only true fits are the
  two background pre-fits.
  - Background tabs renamed to **Power-law Bckg** and **Flat Bckg**.
    Each tab now has its own Q range fields, a **Set from cursors**
    button, and a **Fit *Bckg** button that runs the pre-fit over that
    range. The Power-law fit is a log-log linear regression; the flat
    fit is a median over the high-Q window after the power-law has been
    subtracted.
  - Two-phase parameters box now has an **Input mode** combo with three
    options: *Input φ → derive Δρ²*, *Input Δρ² → derive φ*, *Use both
    as-is*. The non-input field is greyed out and auto-updated.
  - **Fit / Cancel / MC uncertainty / Revert** buttons removed from
    the GUI; the main action is now a single big green **Calculate 3D**
    button that runs `compute_voxelgram()` synchronously at the render
    resolution.
  - Headless `pyirena.batch.fit_saxs_morph()` rewritten to follow the
    same three-step workflow: Power-law pre-fit → Flat pre-fit →
    Calculate 3D. Each pre-fit is skipped if its Q range is missing.
  - State schema bumped to v2 with new `input_mode`, `power_law_q_min/max`,
    `background_q_min/max` fields. Old v1 states migrate automatically
    (`link_phi_contrast=True` ⇒ `input_mode='phi'`).
  - New engine helpers `fit_power_law_bg()`, `fit_flat_bg()`, and
    `derive_phi_from_invariant()`. Documented `Engine.fit()` and
    `calculate_uncertainty_mc()` as deprecated (kept for back-compat).
  - 11 new unit tests (26 total) covering the helpers, the φ ↔ contrast
    round-trip, and the three input-mode resolutions.
- Documentation rewrite of `docs/saxs_morph_gui.md` to describe the new
  workflow, panel layout, parameter table, and scripting recipes.

### Fixed

- **gui3d**: vtk version pin `<9.5` excluded the only currently-available
  PyPI wheels (vtk 9.6.x). Upper bounds dropped on `pyvista` and `vtk`
  so pip can resolve to the latest available release on any platform.

## [0.5.0]

### Added

- **SAXS Morph (3D voxelgram) tool** (Issue #5) — generates a 3D voxelgram of
  a two-phase porous structure from experimental I(Q) using the Gaussian
  Random Fields method (Berk 1991, Roberts 1997, Levitz 2007), then computes
  the model I(Q) from the voxelgram and refines volume fraction, contrast,
  and Power-law + Flat background by fitting the model back to the data.
  Includes:
  - `pyirena.core.saxs_morph.SaxsMorphEngine` — pure-numpy/scipy engine with
    `compute_voxelgram()`, `fit()` (least_squares + Nelder-Mead), and
    `calculate_uncertainty_mc()`. The fit loop is hard-clamped to ≤256³
    voxels for memory safety; the final post-convergence voxelgram uses the
    user-selected render resolution (up to 512³).
  - GUI panel `pyirena-saxsmorph` with:
    - left controls: Voxel grid (fit/render combos + box size + RNG seed),
      Two-phase parameters (φ + contrast + invariant link), Power-law / Flat
      background tabs, action row with Graph Model / Fit / Cancel / MC
      uncertainty / Revert.
    - top right: log-log I(Q) plot with two cursors and three traces
      (data, data−background, red model).
    - bottom right: 2D slice viewer (axis combo + slider) and interactive
      3D PyVista viewer (binary isosurface via flying_edges, right-click
      menu for color / outline / screenshot). Each viewer has a Pop out ⤢
      button that reparents it to a resizable QDialog.
  - Headless `pyirena.batch.fit_saxs_morph()` for batch-fitting.
  - HDF5 storage at `entry/saxs_morph_results/` with the binary uint8
    voxelgram stored gzip-compressed and chunked `(N, N, 1)` for cheap 2D
    slice loads. Stores the RNG seed for reproducibility.
- **Data Selector**: new SAXS Morph (GUI) and SAXS Morph (script) buttons
  on row 8 (purple `#8e44ad` / `#6c3483`); new SAXS Morph (3D) entry in the
  Models menu.
- **New optional dependency group `gui3d`** (PyVista 0.45-0.48 + pyvistaqt
  ≥ 0.11 + VTK 9.3-9.4). Install with `pip install pyirena[gui3d]`. The
  rest of pyirena still works without VTK; the SAXS Morph 3D viewer pane
  shows an install hint when PyVista is missing, while the 2D slice viewer
  and I(Q) plot still function.
- **New script entry point** `pyirena-saxsmorph`.

### Changed

- **`gui` extra**: pinned `PySide6 != 6.7.*, != 6.10.*` to avoid the known
  pyvistaqt incompatibility in those releases (also matters for
  `_SafeInfiniteLine` cursors used by existing tools).
- **Version**: bumped to 0.5.0 to mark the addition of the new 3D-tools
  category.

## [0.4.8]

### Added

- **Data Selector: Export to ASCII** — new orange "Export to ASCII" button next
  to "Tabulate Results".  For each selected HDF5 file, writes plain-text `.dat`
  files into an `ascii_export/` subfolder next to the source file.  Each `.dat`
  carries up to 25 lines of `# key = value` metadata (sample, instrument,
  wavelength, energy, Kfactor, OmegaFactor, Q range, units, model parameters).
  Filenames per source HDF5 (gated by checkboxes, see "Checkbox semantics"):
    - `{stem}.dat` — primary data, columns Q I dI
    - `{stem}_unif.dat` — Unified Fit, Q I_model I_data dI
    - `{stem}_simp.dat` — Simple Fits, Q I_model I_data dI
    - `{stem}_waxs.dat` — WAXS Peak Fit, Q I_fit I_data dI
    - `{stem}_sdQI.dat` — Sizes I(Q) curve, Q I_model I_data dI
    - `{stem}_sdSD.dat` — Sizes distribution, r volume_dist number_dist [std]
    - `{stem}_modQI.dat` — Modeling total, Q I_model I_data dI
    - `{stem}_modP1.dat`, `_modP2.dat`, … — Modeling per population:
        size_dist → r volume_dist number_dist;
        diffraction_peak → Q I_peak;
        unified_level populations are skipped (no distribution to plot).
  USAXS files: only the desmeared sasdata group is exported.  Files with only
  a slit-smeared `_SMR` variant are silently skipped (typically indicates
  reduction never produced desmeared data — a problem with the file itself).
- **Data Selector: checkbox semantics for export** — the "Data" checkbox now
  gates whether the primary `{stem}.dat` is written; uncheck to export only
  model curves without re-writing data.  Each fit-result checkbox gates its
  own model `.dat` files when "Also write model curves" is enabled in
  Configure.
- **Configure dialog: ASCII Export Options** — four new settings under a new
  "ASCII Export Options" group: column delimiter (Space default / Comma),
  significant figures (7 single-precision-safe default / 12 double-precision),
  write metadata header (on by default), and also write model curves
  (on by default).  Defaults match old Fortran-tolerant tools.  Adds a
  reusable `'choice'` field type to the Configure dialog spec system.
- **Simple Fits: Debye Polymer Chain model** — single-chain Debye form factor
  `I(Q) = Scale · 2(exp(−x) − 1 + x) / x²` where `x = Q²·Rg²`.  Parameters:
  Scale (forward intensity) and Rg (radius of gyration).  Supports complex background.
- **Documentation: HDF5 Structure Reference** — new `docs/HDF5_Structure_Reference.md`
  documents every HDF5 group structure pyirena reads and writes (raw NXcanSAS
  data, all 7 result groups, sample/instrument/metadata sub-trees, desmeared
  vs SMR detection, backward-compatibility notes, Igor Pro and Python loading
  recipes, glossary of acronyms).  Suitable for an AI assistant writing Igor Pro
  importer code.

### Changed

- **Simple Fits: complex background parameter rename** — `BG_A` → `BG_G` and
  `BG_n` → `BG_P`, matching the Unified Fit amplitude/exponent naming convention
  (G for Guinier-type prefactor, P for power-law exponent).  Existing HDF5 results
  with the old names are not automatically migrated (use-complex-bg fits will need
  to re-export parameters from the GUI after upgrading).
- **Data Selector: button layout** — Tabulate Results moved from full-width
  to half-width to accommodate Export to ASCII alongside it.

### Fixed

- **Modeling: Unified Fit Level G=0 enforcement** — when the user enters G=0
  on a Unified Fit Level population in Modeling, Rg is now automatically set
  to 1e10 and both "Fit?" checkboxes are unchecked.  This mirrors the
  long-standing behavior of the standalone Unified Fit tool and keeps the
  Guinier term `G·exp(-q²Rg²/3)` numerically inert; previously the fitter
  could oscillate on Rg when the Guinier amplitude was zero.
- **Simple Fits: panel resizes with window** — the graphics widget and the
  control/graph splitter now expand vertically when the window is enlarged;
  previously they stayed at preferred size while the window grew empty
  margins.  Other tools were unaffected because their graph widgets sit
  inside QTabWidgets which expand by default.

## [0.4.7]

### Fixed

#### Unified Fit: stop Porod-tab frame from "breathing" on parameter scrub

`init_plots()` was calling `porod_layout.clear()` on every auto-update
(every parameter scroll tick).  `clear()` empties the
`QGraphicsScene`, collapsing its bounding rect, which fires
`updateGeometry()` on the parent widget — the tab frame briefly
shrinks.  Then `addPlot()` grows the scene back, firing a second
`updateGeometry()` — the frame expands again.  Two layout passes per
tick = visible breathing.  The I-Q tab has the same `clear()` call
but two stacked plots give the layout a stable anchor; the single-plot
Porod layout exposes the jitter.

`init_plots()` now accepts a `rebuild_porod` flag (default `True`).
When `False`, the Porod layout widget is left intact and only the
`PlotItem`'s data items are cleared via `porod_plot.clear()`
(pyqtgraph's `removeItem` loop, which also strips legend entries).
The `QGraphicsScene` bounding rect never changes, no `updateGeometry`
fires, and the frame is completely stable.

`graph_unified()` (the auto-update path triggered by parameter scrubs)
now calls `init_plots(rebuild_porod=False)`.  All other callers
(`set_data`, after-fit, state-restore) still use the default full
rebuild.

#### Unified Fit: eliminate Porod-tab scale flicker during parameter scrub

When parameters were scrubbed with the mouse wheel, the Porod plot
flashed an auto-ranged "everything visible" view for a split-second
before snapping back to the user's saved range.  Root cause:
`init_plots()` left auto-range enabled on the Porod view-box, so each
`.plot()` call inside `plot_data_porod` (scatter and the wide error
bars) triggered a transient `sigRangeChanged` → repaint at the
auto-computed extent before the explicit `setYRange()` at the end of
the method took effect.  The I-Q tab uses the same pattern but is
imperceptible because raw I has a much tighter log-decade range than
I·Q⁴ — the Porod transform amplifies high-Q dramatically and makes
the flash visible.  Fix: call `viewBox.disableAutoRange()` before
adding any items in `plot_data_porod`, so the view stays put through
the `.plot()` calls and only repaints once at the final explicit
range.

#### Unified Fit: cap tick-label density on all axes (~12 max)

When zoomed into ~3 decades, pyqtgraph's log axis labelled every
sub-decade value (1, 2, 3, …, 9 within each decade), producing ~30
cramped, unreadable labels.  A new `_LimitedAxisItem` subclass in
`pyirena/gui/sas_plot.py` overrides `tickValues()` to keep the major
decade tier intact and stride-thin the finer minor tiers so the total
label count stays ≤ 12.  Applied to the left+bottom axes of the I-Q
main plot, the residuals plot, and the new Porod plot.  The top
`RadiusAxisItem` already had its own thinning (max 8) and is unchanged.

#### Unified Fit: Porod tab now preserves user zoom across auto-updates

The new Porod tab inherited the same auto-range bug the I-Q tab had
before its zoom-preservation fix: every parameter scroll triggered
`init_plots()` + `plot_data_porod()`, which re-applied the percentile-
based default Y-range and snapped any user zoom back to the full data
extent.  `graph_unified()` now saves the Porod view-box range before
the rebuild (mirroring the existing main + residuals logic) and
restores it afterward, so a zoomed Porod view stays put while the user
scrubs parameters with the mouse wheel.

### Added

#### Unified Fit: Porod-presentation tab (I·Q⁴ vs Q)

The Unified Fit graph window is now a tabbed view.  The first tab
("I vs Q") is unchanged — main I-Q log-log plot, residuals strip, and
draggable cursors.  The second tab ("Porod (I·Q⁴)") shows the same
data and unified-fit model in Porod presentation: log(I·Q⁴) vs log(Q)
in a single full-height plot.  Local fits and Porod fits drawn via
"Display local fits" appear on both tabs.  This presentation makes
it easier to judge how many levels are needed and where they should be
placed.  No cursors on the Porod tab — it is purely for visual
inspection.  Right-click on the Porod plot offers the same JPEG and
Igor ITX export actions.

## [0.4.6]

### Fixed

#### Unified Fit: zoom preserved across recalculations

When "Update automatically?" is checked, every parameter change
triggered a full `init_plots()` call which re-enabled auto-range,
snapping the view back to the full data extent.  The graph window now
saves the pyqtgraph view-box ranges before each plot rebuild and
restores them afterward, so a user-zoomed region is preserved across
automatic recalculations.  Auto-range continues to apply on the first
plot or whenever the user explicitly zooms out to the full view.

#### Unified Fit: scroll-wheel parameter changes no longer hang the GUI

With "Update automatically?" enabled, each scroll-wheel notch on a
parameter field emitted an `editingFinished` signal and immediately
called `graph_unified()`.  Rapid scrolling therefore queued dozens of
full recalculations that were executed one after another, freezing the
GUI for several seconds after the user stopped scrolling.

Auto-updates now use a **throttle** rather than a debounce: the first
change arms a 150 ms single-shot timer that is *not* restarted by
subsequent events.  When the timer fires, `graph_unified()` reads the
current (latest) field values, then re-arms once if further events
arrived during the blocking call.  The maximum lag before the graph
updates is ~150 ms + one recalculation, and the GUI stays responsive
throughout because the timer is never reset mid-scroll.

#### Unified Fit: out-of-range cursors now clamp to data boundaries

Previously, if either cursor was positioned outside the loaded data
range (common when fitting a different dataset with a narrower Q
range), both cursors were reset to the central 20–80% of the log Q
span.  Cursors are now clamped individually: cursor A snaps to Q_min
if it is below the data, cursor B snaps to Q_max if it is above.  The
other cursor is left untouched.  This preserves the user's intent
(e.g. "include all high-Q data") rather than discarding both
positions.  The default placement on first load is also widened from
20–80% to 10–90% of the log Q range.

## [0.4.5] - 2026-04-25

### Fixed

- **Startup crash on macOS (and any pyqtgraph build with compiled
  subpackages).** The `sendDragEvent` safeguard introduced in 0.4.4 used a
  module import path (`pyqtgraph.GraphicsScene.GraphicsScene`) that works on
  Windows but fails on macOS conda where pyqtgraph compiles that subpackage,
  raising `AttributeError`/`ImportError` before the GUI could open. The import
  now tries both paths with a full fallback chain and, if the class cannot be
  located at all, skips the patch silently rather than aborting startup.

## [0.4.4] - 2026-04-24

### Fixed

- **Unified Fit: GUI no longer crashes when moving cursors after parameter
  changes.** Toggling controls such as the number of levels or "Display local
  fits?" rebuilds the plot, which was destroying the cursor C++ objects while
  pyqtgraph's `GraphicsScene` still held references to them via
  `lastHoverEvent`. The next mouse drag then raised
  `RuntimeError: Internal C++ object (_SafeInfiniteLine) already deleted` and,
  on macOS, took the whole window down. Cursors are now explicitly detached
  before each plot rebuild, and a defensive guard inside the shared SAS plot
  module swallows any remaining stale references rather than crashing.
- **Unified Fit: Sv/Invariant calculation no longer raises
  `ZeroDivisionError` when fitted P = 3.** The Porod-tail integrand is
  `B/Q` at P = 3, whose closed-form `(3 − P)` denominator is singular; the
  tail term is now skipped in that degenerate case so the calculation
  completes instead of erroring out.

## [0.4.3] - 2026-04-22

### Added

#### Diffraction Lines overlay in WAXS Peak Fit (issue #4)

The WAXS Peak Fit window now has two tabs in the left control area:

1. **WAXS Peak Fit** — the existing peak-fitting controls, unchanged.
2. **Diffraction Lines** — a new tab for overlaying theoretical powder
   diffraction stick patterns on the experimental I(Q) curve, to help
   identify crystallographic phases.

Workflow:

- **Import CIF** files via a file-picker that remembers the last folder used.
- **AMCSD** and **COD** buttons open the free crystallography databases in
  your browser so you can download CIFs.
- Each imported CIF appears in a list with: visibility checkbox,
  colour swatch (click to recolour), phase name (auto-detected from the CIF
  formula), per-phase scale factor (× auto-scale to the data peak), and a
  toggle to show Miller-index `(hkl)` labels above each stick.
- **Right-click** a CIF row or click the red `×` button to remove that phase.
  Reset to Defaults clears the entire CIF list.
- **Wavelength** is auto-detected from the loaded NXcanSAS file
  (`/entry/instrument/wavelength`) when "Auto from file" is enabled, and
  can be overridden manually (default 1.5406 Å, Cu Kα).
- The CIF list, wavelength, colours, and last-folder are persisted in
  `state.json` and restored on next launch.

The theoretical patterns are computed by the
[Dans-Diffraction](https://github.com/DanPorter/Dans_Diffraction) library
(Apache-2.0). It is added as a GUI optional-dependency; install with
`pip install pyirena[gui]` or `pip install Dans-Diffraction`.

New files: `pyirena/core/diffraction_lines.py`,
`pyirena/gui/diffraction_lines_panel.py`.

## [0.4.2] - 2026-04-21

### Added

#### Igor Pro ITX export on all graph right-click menus
Every interactive graph in pyIrena now offers **"Save as Igor Pro ITX…"** in its
ViewBox right-click menu, alongside the existing "Save graph as JPEG…" action.

- New public helper `save_itx_from_plot()` in `pyirena/gui/sas_plot.py` collects
  all named `PlotDataItem` objects from a plot (data scatter, model curves) and
  writes them as Igor Pro waves with display, log-axis, color, axis-label, and
  legend commands.  Unnamed error-bar segments are skipped automatically.
- Axis labels and title are auto-extracted from the plot when not supplied.
- Correct log-mode handling: `getOriginalDataset()` is used (not `getData()`) so
  that the physical linear values (Q in Å⁻¹, I in cm⁻¹, radius in Å, …) are
  written to the ITX waves, and the `ModifyGraph log=1` commands applied by Igor
  Pro provide the log scale — no double-log artefact.
- Works correctly for mixed-mode plots (e.g. size distribution: log-x / linear-y).
- Covered tools: Data Selector raw viewer, all fit-result windows, Simple Fits,
  Unified Fit, Size Distribution, Modeling, WAXS Peak Fit, Data Manipulation,
  Data Explorer (already had ITX; menu entry unchanged).

#### Legend text color
All legends across every tool now use **black text** (`labelTextColor='k'`), which
is more readable on the white background used throughout pyIrena.  The canonical
value `SASPlotStyle.LEGEND_TEXT_COLOR = 'k'` is defined in `sas_plot.py`.

### Fixed

- **Size Distribution — shape mismatch with negative intensities**: fitting
  data that contained any non-positive or non-finite I values (typical of
  Data Merge output where flat-background subtraction can drive a few
  end-of-DS1 points negative) crashed with
  `operands could not be broadcast together with shapes (323,) (326,)`.
  - `SizesDistribution.fit()` now exposes the actual q / I_data / err it used
    after its internal `(q > 0) & (I > 0) & finite` mask, via the new
    `result['q']`, `result['I_data']`, `result['err']` keys.
  - `SizesFitPanel.run_fit()` and `store_results_to_file()` use these fit-side
    arrays for the complex-background curve, model overlay, residuals plot,
    and saved HDF5 datasets — guaranteeing matching shapes.
  - Backwards-compatible: if `result['q']` is absent (legacy fit objects), the
    panel falls back to its own cursor-range arrays.

- **Data Merge / Data Manipulation — strip non-positive intensities at save**:
  `save_merged_data()` and `save_manipulated_data()` now drop any points with
  non-positive or non-finite Q / I (with matching dI / dQ entries) before
  writing the HDF5 file, mirroring Igor Pro's behavior.  Prints
  `[data_merge] Stripped N non-positive/non-finite point(s)…` (or
  `[data_manipulation]`) when stripping occurs.  Refuses to save if fewer
  than 2 points remain.

- **Size Distribution — MaxEnt `sky_background` auto-correction** (closes #3):
  two-layer self-correction prevents divergence when the starting sky-background
  value is out of range.
  - Layer 1: if χ² > 100 × M after the first run, retries with sky / 100, / 1000,
    / 10 000 until convergence is restored.
  - Layer 2: after convergence, if sky > 5 % of max(distribution) the value is
    recalibrated to 1 % of max and the fit re-runs once.
  - The corrected value is written back to `maxent_sky_background`, reported via
    `result['sky_note']`, shown in the GUI sky-background field, and flagged with
    an orange status message.  Headless/batch runs receive a `log.info` note.

- **Size Distribution batch defaults** (`batch.py`): default `maxent_max_iter` and
  `tnnls_max_iter` reduced from 1000 to 300 when no saved state is present,
  matching the typical GUI default and avoiding unexpectedly long batch runs.

## Released as first beta

## [0.3.2] - 2026-04-05 — First public beta

First beta release published to PyPI. Install with `pip install pyirena[gui]`.

### Added

#### Cylinder form factors for Modeling tool
Two new form factor entries for disk-like and rod-like particles:

- **Cylinder (Aspect Ratio)** (`cylinder_ar`): half-length L = AR·R, where R is the
  radius from the distribution. Disk-like for AR < 1, rod-like for AR > 1.
- **Cylinder (Length)** (`cylinder_length`): fixed total height H [Å], independent of
  radius. Volume scales as R² (not R³).
- Derived quantities (Rg, specific surface) use cylinder geometry.
- Orientationally-averaged form factor via 50-point Gauss-Legendre quadrature.

#### Core-shell form factors for Modeling tool
Five core-shell form factor variants with three polydispersity modes, following the
Igor Pro implementation:

- **Core-Shell Sphere** — by core R, by shell thickness, or by total R.
- **Core-Shell Spheroid** — by core R or by total R (with shared aspect ratio).
- SLD parameters (sld_core, sld_shell, sld_solvent in 10⁻⁶ Å⁻²) replace the scalar
  contrast parameter, which is automatically locked to 1.0 and hidden in the GUI.
- SLD-contrast-weighted Rg; specific surface based on outer radius (Porod law).
- Volume convention: `scale` and `volume_fraction` are based on total particle volume
  (core + shell together), not core-only. Documented in code and user guide.

#### Modeling GUI — standard buttons
Added standard buttons matching other tools (Simple Fits, Unified Fit, Sizes):

- **Results to graph** (#81c784 green) — annotate I(Q) plot with fitted parameter values.
- **Save State** (#3498db blue) — save current parameters to state file.
- **Import Parameters** (lightgreen) — import from pyIrena JSON config file.
- **Reset to Defaults** (#e67e22 orange) — reset all populations to default values.
- **Store in File** / **Export Parameters** renamed and recoloured to match standard.
- Auto-save state on window close (`closeEvent`).

### Fixed

- **Modeling "Open..." button**: loading a new file now clears all previous data, model
  curves, distributions, residuals, and annotations from the graph. Previously old items
  persisted because `graph_model()` and `_on_fit_complete()` bypassed the data tracking.
- **Modeling `load_file()`**: replaced nonexistent `load_nxcansas` import with the correct
  `readGenericNXcanSAS` function from `pyirena.io.hdf5`.
- **Core-shell G-matrix units**: corrected from `F²×1e-16` to `F²/V_total×1e-4`. The old
  formula double-counted the SLD unit conversion and missed the division by particle
  volume, causing intensities ~10⁶ times too low.

## [0.2.1] - 2026-03-20

### Added

#### Modeling tool — parametric size-distribution fitting
New analysis tool for forward-modelling small-angle scattering data using
parametric size distributions, Unified Fit levels, and diffraction peaks.

- **5 distribution functions**: Gaussian, LogNormal (3-parameter shifted), LSW,
  Schulz-Zimm (Gamma), Ardell — Igor-style CDF-inversion radius grid.
- **3 population types** combinable in up to 5 simultaneous populations:
  - `SizeDistPopulation`: distribution × form factor (sphere, spheroid) × optional
    structure factor (interferences / hard-sphere Percus-Yevick).
  - `UnifiedLevelPopulation`: Beaucage Unified level G·exp(−q²Rg²/3) + B·Q*⁻ᴾ
    with optional Born-Green correlations.
  - `DiffractionPeakPopulation`: Gaussian, Lorentzian, or pseudo-Voigt peak at Q₀.
- **Engine** (`pyirena/core/modeling.py`): `ModelingEngine` with G-matrix caching,
  `scipy.optimize.least_squares` (TRF) or Nelder-Mead fitting, and MC uncertainty.
- **HDF5 I/O** (`pyirena/io/nxcansas_modeling.py`): save/load for all population types.
- **GUI panel** (`pyirena/gui/modeling_panel.py`): multi-tab population editor,
  I(Q) log-log plot with per-population overlays, distribution preview, Export/Import
  Parameters, and MC uncertainty estimation.
- **Batch API** (`pyirena/batch.py`): `fit_modeling()` headless fitting function;
  registered in `fit_pyirena()` under the `modeling` config key.
- **Data Selector** integration: Modeling (GUI/script) buttons at row 4.

#### Data Merge tool — SAXS/WAXS merging
New tool for merging two SAS datasets (e.g. SAXS + WAXS) onto a common Q scale.

- **Engine** (`pyirena/core/data_merge.py`): Nelder-Mead optimisation of scale factor,
  flat background (DS1), and optional Q-shift; log-log linear interpolation in the
  overlap region.
- **HDF5 I/O** (`pyirena/io/nxcansas_data_merge.py`): copies DS1 NXcanSAS file and
  replaces Q/I/Idev/Qdev with merged data; appends `entry/data_merge_results`.
- **GUI panel** (`pyirena/gui/data_merge_panel.py`): dual dataset loader, SAXS/WAXS
  plot mode toggle, cursor-driven overlap range, optimisation controls, and batch mode.
- **Batch API**: `merge_data()` headless merge function.
- **CLI entry point**: `pyirena-datamerge` console script.

#### Scattering Contrast Calculator
New tool for computing X-ray and neutron scattering length densities and contrast.

- Compound-level SLD calculation from chemical formula, density, and X-ray energy.
- Supports neutron SLDs with isotope substitution.
- Phase-pair contrast (ΔρX)² and (ΔρN)² tables.
- **GUI panel** (`pyirena/gui/contrast_panel.py`) with compound library, interactive
  crosshair energy cursor, and JPEG export.
- **CLI entry point**: `pyirena-contrast` console script.

### Fixed

- **`fit_modeling` batch function**: previously only deserialized `SizeDistPopulation`
  from config; now correctly rebuilds `UnifiedLevelPopulation` and
  `DiffractionPeakPopulation` via `pop_type` dispatch.
- **`fit_pyirena`**: `modeling` config key was missing from `_TOOL_REGISTRY`; added
  so `fit_pyirena()` automatically runs `fit_modeling()` when a `modeling` section is
  present in the config file.
- **`pyirena/__init__.py`**: `__version__` was not updated from 0.1.1; now kept in
  sync with `pyproject.toml`.

## [0.1.2] - 2026-02-22

### Added

#### Simple Fits tool — 13 single-model analytical fits
New analysis tool ported from Igor Pro `IR3_SimpleFits.ipf` and
`IR3_SystemSpecificModels.ipf`.  Provides direct analytical model fitting with
linearization plots, Monte Carlo uncertainty estimation, and full HDF5 I/O.

**Core (`pyirena/core/simple_fits.py`)**
- `SimpleFitModel` dataclass: holds model name, parameters, limits, and complex-background
  flag; provides `fit(q, I, dI)`, `to_dict()`, `from_dict()`.
- `MODEL_REGISTRY` dict mapping every model name to its parameter definitions, formula
  function, linearization type, and complex-background flag.
- 13 supported models (parameter names in parentheses):

  | Model | Parameters | Linearization |
  |-------|-----------|---------------|
  | Guinier | I0, Rg | ln(I) vs Q² |
  | Guinier Rod | I0, Rc | ln(QI) vs Q² |
  | Guinier Sheet | I0, Rg | ln(Q²I) vs Q² |
  | Porod | Kp, Background | IQ⁴ vs Q⁴ |
  | Power Law | P, Exponent, Background | — |
  | Sphere | Scale, R | — |
  | Spheroid | Scale, R, Beta | — |
  | Debye-Bueche | Prefactor, Eta, CorrLength | — |
  | Treubner-Strey | Prefactor, A, C1, C2 | — |
  | Benedetti-Ciccariello | SolidSLD, VoidSLD, LayerSLD, Sp, t | — |
  | Hermans | B, s, d1, d2, sigma1, sigma2 | — |
  | Hybrid Hermans | Hermans params + G2, Rg2, G3, Rg3, B3, P3 | — |
  | Unified Born Green | G1, Rg1, B1, P1, G2, Rg2, B2, P2, eta, ksi | — |

- Optional complex background `BG_A·Q⁻ⁿ + BG_flat` on all models except Porod and
  Power Law (which have an explicit flat Background parameter).
- Sphere uses Gauss-Legendre quadrature orientational average for efficiency.
- Treubner-Strey derives correlation length ξ and repeat distance d from fit params.
- Guinier Sheet derives layer thickness from Rg.
- `fit()` uses `scipy.optimize.curve_fit`; returns chi², reduced chi², DOF,
  `params_std` (from covariance diagonal), `residuals`, and model-specific `derived`.

**HDF5 I/O (`pyirena/io/nxcansas_simple_fits.py`)**
- `save_simple_fit_results(filepath, result, model_obj, intensity_data, intensity_error)`:
  writes/overwrites `entry/simple_fit_results` NXprocess group.
- `load_simple_fit_results(filepath)`: loads all fit scalars, arrays, params, params_std,
  and derived quantities.
- `print_simple_fit_results(result)`: formatted console summary.

**GUI panel (`pyirena/gui/simple_fits_panel.py`)**
- Three-panel layout: I(Q) log-log + model overlay (with cursor-selectable Q range),
  normalised residuals, linearization plot (Guinier/Porod families only).
- Cursor-driven Q range selection linked to the I(Q) plot.
- Per-parameter widgets: value, lower/upper limits, "Fit?" checkbox; "No Limits" toggle.
- Complex background checkbox: adds BG_A, BG_n, BG_flat sub-parameters.
- **Fit** button: runs `SimpleFitModel.fit()`, displays results with ± uncertainties.
- **Calculate Uncertainty** button: Monte Carlo loop (configurable N runs); updates ± labels.
- **Store in File** button: saves results to HDF5.
- **Export / Import Parameters** buttons: read/write `simple_fits` section in a shared
  `pyirena_config.json` (interoperable with Unified Fit and Sizes config files).
- Linearization panel shows transformed data + linear fit + intercept annotation for
  Guinier/Porod models; shows "No linearization available" for other models.
- Highlighted Q-range band on linearization scatter plot.
- State fully persisted via `StateManager`; stale parameters removed on model change.

**Batch API (`pyirena/batch.py`)**
- `fit_simple(data_file, config, with_uncertainty, n_mc_runs, q_min, q_max)`: headless
  fitting from a `SimpleFitModel` or plain dict config; saves results to HDF5.
- `fit_simple_from_config(data_file, config_file, save_to_nexus, ...)`: wrapper that
  reads the `simple_fits` section of a `pyirena_config.json` (normalises GUI state-dict
  keys, extracts q_min/q_max) and calls `fit_simple()`.
- `fit_pyirena()` now dispatches `simple_fits` config sections to `fit_simple_from_config()`.
- `fit_simple` exported from `pyirena.__init__`.

**Data Selector integration (`pyirena/gui/data_selector.py`)**
- **"Simple Fits"** checkbox added alongside Data / Unified Fit / Size Dist. checkboxes.
- *Create Graph*: opens `SimpleFitResultsWindow` (I(Q) + model + residuals, matching
  colour scheme to the dataset; model line uses `color.darker(280)` for clear contrast).
- *Create Report*: loads `simple_fit_results` from HDF5; adds "## Simple Fits" section
  with model name, chi², reduced chi², DOF, Q range, all parameters ± std, derived quantities.
- *Tabulate Results*: adds `SF_model`, `SF_chi2`, `SF_reduced_chi2`, `SF_dof`,
  `SF_q_min`, `SF_q_max`, `SF_use_complex_bg`, `SF_<param>`, `SF_<param>_std`, and
  `SF_derived_<name>` columns (dynamic — adapts to the parameters of the fitted model).
- **"Simple Fits (script)"** batch button now correctly calls `fit_simple_from_config()`
  and saves results to the HDF5 file (was falling through to fit_sizes with no save).

#### Public API
- **`load_result(filepath, analysis)`** (`pyirena/io/results.py`) — importable directly
  as `from pyirena import load_result`.  Pass a file path and an analysis name
  (`'unified_fit'`, `'size_distribution'`, or `'simple_fits'`) to retrieve a fully
  documented result dict.  Returns a safe empty structure (`found=False`) if results
  are absent — no exception raised.
- `SUPPORTED_ANALYSES` tuple exported from `pyirena`.

#### Plot improvements (all panels)
- **Phantom-point fix** (I(Q) graph, Unified Fit, Size Distribution): error bar tops
  clipped at `min(I·1000, P99(I)·1000)` — prevents cosmic rays or WAXS peaks at
  extreme intensities from driving the y-axis to 10³⁸.
- **ViewBox hard limits**: y-axis locked to ±3 decades beyond the 99th-percentile data
  range; x-axis locked to the nearest full decade beyond the data Q range ±1 extra
  decade.  Prevents accidental zoom-out to empty space.
- **x-axis zoom constraint**: applied to Simple Fits, Size Distribution, and Unified Fit
  panels.  Example: data 0.003–0.8 Å⁻¹ → zoom limits 0.0001–10 Å⁻¹.

### Fixed
- **Data Selector sort order**: the saved "Sort" pulldown state is now applied to the
  file list immediately on startup (was restored to the combo but the list remained
  alphabetical until the user manually changed the setting).
- **Report generation** (Size Distribution): section now reads scalar parameters from
  the flat dict returned by `load_sizes_results()` (no longer looks for a non-existent
  nested `'params'` key).  Previously all parameters showed as `nan`/`unknown`.
- **Report generation** (Size Distribution): now includes all stored parameters:
  contrast, aspect ratio, log spacing, background, error scale, power law B/P, Q range,
  n_iterations, plus method-specific sub-tables for MaxEnt, Regularization, TNNLS, and
  Monte Carlo.

## [0.1.1] - 2026-02-19

### Added

#### Size Distribution — Complex Background Model
- `SizesDistribution.compute_complex_background(q)`: evaluates `B·q⁻ᴾ + background`
  for arbitrary q arrays.
- `SizesDistribution.fit_power_law(q, I, q_min, q_max, fit_B, fit_P)`: fits the
  power-law amplitude B and/or exponent P to a user-selected Q range using
  `scipy.optimize.curve_fit`.
- `SizesDistribution.fit_background_term(q, I, q_min, q_max)`: estimates the flat
  background by averaging `I − B·q⁻ᴾ` in a selected Q range.
- New class attributes `power_law_B` and `power_law_P`; `fit()` now subtracts the
  full complex background before inverting, and returns `n_data` in the result dict.

#### Size Distribution GUI (`sizes_panel.py`)
- **Background tab**: new panel tab with "Power-Law Term B·q⁻ᴾ" and "Flat Background"
  groups; each has its own Q-range fields, "Set Q from cursors" button, and individual
  fit buttons ("Fit P/B", "Fit Background").
- **"Fit Sizes" / "Fit All"** buttons: "Fit All" runs power-law fit → background fit
  → size distribution fit sequentially in one click.
- **Corrected-data overlay**: I(Q) − complex background shown as blue triangles on the
  main graph, limited to the cursor Q range; complex background shown as a dashed grey
  line.
- **"Calculate Uncertainty (MC)"** button: runs 10 Monte-Carlo fits on
  Gaussian-perturbed data; reports per-bin mean/std of the distribution as error bars,
  and propagates Rg, Vf, and peak-r uncertainties with ± notation in the Results box.
- **"Export Parameters" / "Import Parameters"** buttons: save and load all Sizes
  parameters to/from a `pyirena_config.json` file sharing the same `_pyirena_config`
  envelope used by Unified Fit (files are interoperable between tools).
- **Layout improvements**:
  - Q min and Q max shown on one row in the Q Range group.
  - Number-of-bins spinbox and "Logarithmic spacing" checkbox share one row.
  - MaxEnt, Regularization, and TNNLS sub-controls each condensed to a single row.
  - Results box condensed to two rows (χ²/Vf and Rg/Peak r).
- **Finer mouse-wheel steps** for Error scale: `ScrubbableLineEdit` now accepts a
  `step_factor` parameter (default 0.1); Error scale uses 0.02 for precise control.
- **Fit curve** rendered with `width=4` pen and plotted on top of all other items.
- **Error scale live preview**: changing the Error scale field immediately redraws the
  error bars without requiring a new fit.

#### Batch API (`batch.py`)
- `fit_sizes(data_file, config_file, save_to_nexus=True)`: headless size-distribution
  fitting function analogous to `fit_unified()`.  Reads a `'sizes'` group from a
  pyIrena config JSON, applies the saved cursor Q range, runs the fit, optionally saves
  results to NXcanSAS HDF5, and returns a structured result dict.
- `fit_pyirena()` now automatically dispatches to `fit_sizes()` when a config file
  contains a `'sizes'` group.
- `fit_sizes` added to the public API (`pyirena.__init__`).

#### Data Selector (`data_selector.py`)
- **"Size Dist." checkbox** added next to the existing "Data" and "Unified Fit"
  checkboxes.
  - *Create Graph*: if checked, opens the Size Distribution panel with the selected
    file's data.
  - *Create Report*: if checked, loads stored size-distribution results from HDF5 and
    includes a "## Size Distribution" section (chi², Vf, Rg, peak r, method, shape,
    residuals stats) in the Markdown report.

#### State management (`state_manager.py`)
- `DEFAULT_STATE["sizes"]` extended with `power_law_B`, `power_law_P`,
  `power_law_q_min`, `power_law_q_max`, `background_q_min`, `background_q_max`.
- Schema migration from version 1 → 2 resets `n_bins` (50 → 200), `log_spacing`
  (False → True), and initialises `error_scale` on old saved states.

### Changed
- `SizesDistribution.fit()` now subtracts the complex background `B·q⁻ᴾ + background`
  (previously only the flat background was removed).
- Default `n_bins` changed from 50 to 200 and `log_spacing` from False to True
  (better resolution for typical SAXS size distributions); old saved states are
  migrated automatically via `schema_version`.
- Sizes panel minimum width increased to 420 px to accommodate the wider layout.

## [0.1.0] - 2024-02-13

### Added
- Initial release of pyIrena
- Unified Fit model for small-angle scattering analysis
- Support for multi-level hierarchical structures
- Mass fractal mode
- Correlation function (Born-Green approximation)
- Parameter linking capabilities
- NXcanSAS HDF5 file support
- Basic plotting and analysis utilities
- Comprehensive documentation and examples
