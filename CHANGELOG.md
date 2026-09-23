# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Carbon model** — a new analysis tool that fits the whole measured range of
  a disordered carbonaceous material (USAXS → SAXS → WAXS, often five decades
  in Q) as one model: grain Porod scattering with optional surface roughness,
  micropore scattering, and turbostratic diffraction peaks. Following
  Saurel et al., *Energy Storage Materials* **21** (2019) 162–173 and its 2020
  corrigendum.

  All three contributions are refined together, because the contrast that
  scales them is computed from the fitted WAXS peak positions and the porosity
  — fitting the regions in sequence converges somewhere else.

  Beyond the formulas, the tool wires the materials-science layer: a chemical
  formula and the fitted lattice spacings give the structural density, the
  porosity gives the sample density, and each density gives an SLD and a
  contrast — yielding about twenty derived quantities (BET-comparable specific
  surface areas, pore and wall widths, stack height, layers per stack,
  d-spacings) rather than just fit coefficients.

  Micropores are described either as dilute pores with optional fractal
  aggregation or by the Teubner-Strey two-phase model (a mode toggle, never
  summed); diffraction peaks can optionally carry the crumpled-layer envelope
  of eq. (16)–(17), whose geometry can be linked to the SAXS region's.

  Wired into every surface: GUI panel, Data Browser launcher and results
  graph, `entry/carbon_fit_results` in NXcanSAS, `fit_carbon` batch path,
  `pyirena.api.read_carbon_fit`, the `carbon` MCP control category (16 tools),
  HDF5 Data Explorer trend plots, Igor export, and
  `docs/carbon_fit_gui.md`.
- **Carbon model: "Revert back".** One level of undo, snapshotted at the start
  of every fit, matching the button Unified Fit has. A full-range model with
  fifteen-plus free parameters will sometimes run to a wrong solution, and
  rebuilding the starting point by hand is the expensive part. Restores the
  whole model — values, bounds, Fit? flags and the peak list — but not the Q
  range, which belongs to the user's cursors rather than to the fit.
- **True Voigt peak shape in WAXS Peak Fit.** The real Lorentzian⊗Gaussian
  convolution via the Faddeeva function, not the existing linear pseudo-Voigt
  mix. Added for the Carbon model, where crystallite size (Gaussian) and layer
  curvature (Lorentzian) broaden independently, but available to WAXS Peak Fit
  too: choose `Voigt` and the Lorentzian component appears as `FWHM_L`
  alongside the Gaussian `FWHM`. Wired through the panel, the HDF5 schema,
  reporting and both export paths.

- **Tool-registration contract test** (`pyirena/tests/test_tool_registration.py`).
  Adding an analysis tool means touching ~45 files, about a dozen of which are a
  key in a hand-maintained registry — and missing one fails silently. The test
  holds the canonical tool table (one row per tool, one column per surface,
  either the registering key or a written reason for its absence) and checks all
  of them. Adding a tool is now: add the row, run the test, work the failures.
  Documented in `docs/developer_adding_features.md` § "Adding a whole new tool".
- **`docs/module_map.md`** — one line per module across all ~120 modules, so a
  shared helper can be found without grepping 113,000 lines.
- **`AGENTS.md`** — the repository orientation file, previously `CLAUDE.md`,
  under the name every coding agent looks for. `CLAUDE.md` now points to it.

### Fixed

- **Carbon model: a fitted fractal dimension could be silently unfittable.**
  `teixeira_structure_factor` clamps D into (1.001, 2.999) because Γ(D−1)
  diverges at one end and the formula degenerates at the other — and outside
  the clamp the function is *flat*, so the finite-difference gradient is
  exactly zero. A D parked there burns the whole evaluation budget without
  moving, which is indistinguishable from a parameter that was never wired up.
  Widening the bounds to the physical (1, 3), which is the natural thing to
  type, put D in exactly that dead zone. Fit bounds are now narrowed into the
  differentiable range automatically and the fit reports that it did so.
- **Carbon model: the fit did not scale its parameters.** They span five orders
  of magnitude (a specific surface area ~1e4 cm²/cm³ beside a fractal dimension
  ~2.5), so an unscaled trust region is sized by the largest and the small ones
  barely move. With `x_scale='jac'` a far-from-solution fit drops from ~3000
  evaluations to ~560, and a case with a large Σ/r ratio now recovers the
  fractal dimension it previously missed by 5 % (χ² better by 14×).
- **Carbon model: a parameter pinned at a bound is now reported.** It is the
  one failure that looks like success — the fit runs, a value is reported, and
  it came from the limit. The panel highlights the field, the status line and
  the report say so, and the agent API returns `pinned_parameters`. A parameter
  resting on a zero floor is reported separately, because "widen the bound" is
  the wrong advice for it.
- **Carbon model: the WAXS zoom panel went blank after a fit.** The cached
  curves shared one Q array between the measured data (full range) and the
  model (the fit range only), so a fit over a restricted range paired a
  284-point model with a 400-point data array; the redraw raised part-way
  through, leaving an empty panel — and the exception escaped the fit itself.
  Data and model are now cached on separate grids and the redraw can no longer
  propagate.
- **Carbon model: the WAXS zoom panel was almost unusable to frame.** It took
  its Q window from the peak *starting guesses* (graphite positions, before a
  fit possibly nowhere near the sample) and left Y to autoscale over the whole
  curve, whose low-Q end is orders of magnitude above the diffraction peaks —
  so the peaks sat flat on the baseline and the only way to see them was typing
  limits into the axis dialog. It now opens on Q ≥ 1 Å⁻¹ out to the highest
  measured Q with the Y axis scaled to that window, re-frames after a fit, and
  has a ⟲ button to restore the view by hand.
- **Carbon model: the filename field read "(no file selected)" when the panel
  was opened from the Data Browser.** `set_data` did not call the loader row's
  `set_filename`, so the field only ever tracked this panel's own Open… button.
- **Cursors were fiddly to grab.** pyqtgraph derives an `InfiniteLine`'s hit
  area from its pen width, giving pyIrena's 2-pixel cursors a band about four
  pixels wide — accurate, and hard to hit, especially where two cursors sit
  close together on a log axis. `_SafeInfiniteLine` now widens the bounding
  rectangle (which is what Qt hit-tests against) to an 18-pixel band without
  drawing a thicker line. Applies to every tool that uses the shared cursors.
- **Carbon model: peak parameters scrubbed far too fast.** The shared
  wheel-scrub step is 10 % of the value's leading decade, which moves a peak
  centre at Q₀ ≈ 1.9 Å⁻¹ by 0.1 Å⁻¹ per notch — about a third of a carbon
  peak's width, so it jumps straight past the target. Peak positions and both
  Voigt widths now scrub ten times finer (0.01 Å⁻¹ and 0.001 Å⁻¹ per notch).
  The amplitude keeps the coarse step, being the one peak parameter you want to
  move in large relative jumps.
- **Carbon model: the fit progress callback was repainting too often.** It
  pumped the Qt event loop every tenth residual evaluation; on a
  badly-conditioned model that is thousands of repaints and made the fit
  several times slower than the arithmetic it was reporting on. Now throttled
  to ten updates a second. Fit tolerances also relaxed from 1e-10 to 1e-8 —
  precision no SAS measurement carries, and worth roughly 2× on a degenerate
  model for results that agree to four significant figures.
- **`set_cursor_q_range()` did not place the cursors where it was asked to.**
  On the path where the cursors do not exist yet it created them with
  `make_cursors()`, which insets them by 10 % of the log span, and then
  returned without moving them — so a restored Q range came back 10 % narrow
  at both ends. Harmless for a single-region fit; for the Carbon model, 10 % of
  five decades is half a decade, and the top half-decade holds the (100)
  reflection whose position the contrast calculation reads.
- **`fit_pyirena` silently skipped `saxs_morph` config sections.**
  `batch.saxs_morph.fit_saxs_morph()` existed with the right signature but was
  never registered in the pipeline's tool registry, so a `saxs_morph` block in
  a `pyirena_config.json` did nothing. Found by the new registration test.
- **`test_window_state.py` was intermittently failing on developer machines.**
  Its live-window tests went through the real `restore_window_state`, which
  asks `screen_rects()` where the monitors are and declines to place a window
  it cannot see — so the tests depended on the ambient display configuration.
  When `screen_rects()` came back empty (macOS under load, display sleep, fast
  user switching) every restore became a no-op and a window asserted at
  (70, 55) arrived at (0, 0). The module now pins the screen layout with an
  autouse fixture; the placement policy itself is pure and was already tested
  with explicit screen rectangles, so no coverage is lost.
- **SAXS Morph and Fractals results leaked into derived data files.**
  `PYIRENA_RESULT_GROUPS` — the list of result groups stripped when a tool
  seeds a new output file from a source file — was maintained by hand and had
  fallen two tools behind, so Data Merge, Data Manipulation and the control
  API's `output_path` carried stale SAXS Morph and Fractals results into the
  new file. The list is now derived from `io/schema.py::TOOL_REGISTRY`.

- Fixed a `QSplitter.setSizes()` quirk that could make the Modeling and
  Unified Fit control panel balloon to a large, unusable fraction of the
  window width when the saved pane-width state was stale or corrupted
  (values much smaller than the window's actual width). Saved splitter sizes
  are now rescaled to the panel's current width before being applied.

### Changed

- Unified Fit: the "Fit" button is now noticeably larger than "Fix limits?"
  so the primary action is easier to hit; "Fix limits?" is slightly shorter.
- Data Selector: right-clicking a file in the file list now offers "Show
  file in Finder/Explorer..." to reveal it in the OS file browser.

## [1.1.1] - 2026-09-12

A maintenance release on top of 1.1.0. It completes the AI/scripting surface
for the three tools that had stopped at the GUI (Scattering Contrast, Data
Manipulation, Data Merge), keeps the MCP tool count inside provider-side
caps, and fixes a data-corrupting misread of slit-smeared NXcanSAS files.
Saved-file compatibility is unchanged; files written by 1.0.x and 1.1.0 load
as before.

### Added

- Data Manipulation and Data Merge are now reachable from scripts and AI
  agents. A new `pyirena.api.data_ops` group adds `average_data`,
  `subtract_data`, `divide_data`, `scale_data`, `trim_data`, `rebin_data`,
  `merge_datasets` and the read-only `match_merge_files` pairing helper,
  exposed over MCP as the dispatcher's `data` category (again no new
  top-level MCP tools — the count stays at 26). Previously an agent asked
  to average or subtract data had no tool for it: only the two read-only
  provenance readers existed, and those return `{"found": false}` on a raw
  file, so the agent would retry and then sweep the whole tool list.
  Output naming is unchanged from the GUI and batch layers — a sibling
  `_manip` / `_merged` folder plus the per-operation filename suffix
  (`_avg`, `_sub`, `_div`, `_scaled`, `_trimmed`, `_rebinned`, `_merged`).
  These wrap `core` + `io` directly rather than `pyirena.batch`, which
  attaches a stdout log handler that would corrupt MCP's stdio transport
  and returns a bare `None` on every failure. The api layer adds guards
  core does not have: mismatched slit smearing is refused for averaging
  (not only for subtract/divide), an empty trim window is an error rather
  than a silent zero-length result, `auto_scale` without both Q bounds is
  rejected instead of silently ignored, swapped merge inputs are caught,
  and every call reports `n_points_in` / `n_points_written` /
  `n_dropped_nonpositive` so an over-subtraction is visible.

- Scattering Contrast is now reachable from scripts and AI agents. A new
  stateless `pyirena.api.calculators` group wraps
  `pyirena/core/scattering_contrast.py`: `calc_contrast` (X-ray and neutron
  contrast between two compounds from formulas + densities, optionally
  anomalous at a given energy), `calc_compound` (one material's SLDs),
  `calc_contrast_energy_scan` (contrast vs energy, for anomalous SAXS
  planning), `lookup_element`, and read-only access to the saved compound
  library (`list_compound_library`, `load_compound`). Previously the tool
  stopped at the core/io/gui layers, so `pyirena-mcp` — a pure wrapper over
  `pyirena.api` — could not see it and agents derived contrasts by hand.
  `calc_contrast` returns (Δρ)² in 10²⁰ cm⁻⁴, the same units the Sizes and
  Modeling `contrast` parameter expects, so results feed straight into a fit.
- `pyirena-mcp` exposes these as a sixth dispatcher category,
  `calculators`, alongside the five fitting categories. The dispatcher
  (`pyirena/mcp/dispatch.py`) now serves several schema registries rather
  than only `pyirena.api.control`, so the group adds **no** top-level MCP
  tools — the registered count stays at 26, and future calculators are free.
- New `contrast` extra (`periodictable`, `xraydb`), referenced from `gui`,
  `mcp` and `all`. `pip install pyirena[mcp]` now gets working calculators
  without pulling in Qt; previously these two packages were reachable only
  through the `gui` extra. `pyirena-doctor` reports them against the new
  extra.

### Changed

- `pyirena-mcp` collapsed ~90 `pyirena_ctrl_*` control tools (Sizes, Simple
  Fits, Modeling, WAXS Peak Fit, and most of Unified Fit) into a fixed
  4-tool dispatcher (`pyirena_list_categories`, `pyirena_list_tools`,
  `pyirena_describe_tool`, `pyirena_call`), reusing the existing
  `pyirena.api.control.schemas.TOOL_SCHEMA_BY_NAME` registry
  (`pyirena/mcp/dispatch.py`). The server's registered tool count drops
  from ~110 to ~26. This was needed because a large combined MCP tool
  count from clients with several active servers can exceed provider-side
  caps on the number of tools in a single request (observed via a gateway
  proxy's 128-tool limit). Session-lifecycle tools (`pyirena_ctrl_
  open_dataset`, `list_open_sessions`, `close_session`,
  `get_session_summary`) are unchanged and still top-level. See
  `docs/ai_tools_reference.md` for the new calling convention.

### Fixed

- **Slit-smeared data with no per-point Q resolution was read incorrectly,
  crashing or corrupting every manipulation.** NXcanSAS lets such data
  declare `Q@resolutions='dQl'` — the scalar slit length as the only
  resolution contribution — and `readGenericNXcanSAS` took the first
  `resolutions` token as the per-point resolution dataset, so `dQ` came
  back as a 0-d scalar. Every consumer that indexes `dQ` then raised
  `IndexError: invalid index to scalar variable` (trim, rebin, average,
  subtract, divide — in the GUI panels as well as batch), while `scale`
  quietly copied the scalar through and the saver wrote the slit length out
  as a 0-d `Qdev`, corrupting the file's resolution metadata in a way that
  read back as a scalar again next time. The file was standards-correct;
  the reader was not. Fixed at three levels: `readGenericNXcanSAS` never
  selects `dQl` as the per-point resolution and validates the shape of
  whatever it does select; `DataManipulation` and `DataMerge` treat a `dQ`
  that is not parallel to `Q` as absent; and `create_nxcansas_file` ignores
  a non-conforming `dq` with a warning instead of writing a malformed
  `Qdev`. Files already carrying a 0-d `Qdev` now read cleanly, and shed it
  the next time any pyirena operation rewrites them. No resolution
  information is lost — for slit-smeared data the resolution is the slit
  length, carried separately as `slit_length` / `dQl`.
- `pyirena.batch.average_data` did not pass `slit_length` to the saver, so
  averaging a slit-smeared series silently produced a pinhole file with no
  `dQl` — corrupting any later desmearing or smeared fit. `manipulate_data`
  already did this correctly.
- Windows: pyqtgraph could fail with `DLL load failed while importing QtCore`
  in an environment that carries both PySide6 and PyQt6. pyqtgraph tries
  PyQt6 *before* PySide6, so a half-installed PyQt6 broke it even though
  PySide6 loaded fine. `pyirena/gui/_qt.py` now pins `PYQTGRAPH_QT_LIB` to
  the binding pyIrena actually imported.
- `pyirena-doctor` reported a broken PyQt6 as `ok (version unknown)` — it
  only imported the (nearly empty) top-level package. It now imports
  `<binding>.QtCore`, warns when both bindings are installed, and no longer
  blames the Visual C++ redistributable for what is a mixed-bindings
  collision.
- Declared `six` as a core runtime dependency. PyIrena's public HDF5 reader
  uses it for metadata handling, so clean installations can now import the
  scattering-data API without an undeclared dependency.

## [1.1.0] - 2026-09-08

First stable release of the 1.1 line. The code is identical to `1.1.0b12`;
this entry summarises what 1.1.0 brings over **1.0.1**. Per-change detail,
including every fix, lives in the `1.1.0b1` – `1.1.0b12` entries below.

Saved-file compatibility is preserved: NXcanSAS result files written by 1.0.x
load unchanged, and `from_dict()` supplies a default for every field added
during the 1.1 cycle.

### Highlights

- **Slit smearing across all fitting tools** (b1). Unified Fit, Sizes,
  Modeling, Simple Fits and WAXS Peak Fit can fit slit-smeared USAXS/Matilda
  data directly, with the smearing geometry carried through result files and
  merge provenance.
- **Fits are substantially faster, with the same numbers.** Analytic Jacobians
  replaced finite differences in the four gradient-based tools — 1.7–2.1× on
  Unified Fit, 1.4–2.2× on Modeling, up to 2× on Simple Fits, ~1.8× on large
  WAXS fits (b11) — on top of a ~3× vectorisation of the Modeling hot path
  (b9) and parallel Monte-Carlo uncertainty (b4). An exact gradient also fixes
  fits finite differences could not do at all, such as a Debye chain over a
  USAXS Q range.
- **pyIrena is agent-drivable.** `pyirena/api/control/` exposes interactive
  fitting sessions for all five fitting tools, wrapped by an MCP server
  (b5–b7), with the control/MCP write surface confined to `PYIRENA_DATA_ROOT`
  (b2).
- **A validation suite with exactly known answers** (`validationData/`, b10):
  synthetic data generated from known parameters, refitted end-to-end into
  `VALIDATION_RESULTS.md`, with **Igor Pro Irena results tabulated alongside**
  for 98 of the quantities. Current status: **195 of 195 scored comparisons
  within tolerance**; median deviation from truth 0.185 % for pyIrena versus
  0.308 % for Igor Irena, and 0.040 % between the two packages.
- **A uniform UX contract across every panel** (b7): clipboard copy, column
  sorting and CSV export on every table; graph copy/PNG/SVG/CSV export
  everywhere; "Copy results" and "Save report…" on the fit panels; regex file
  filters; shared filename sorting; drag-and-drop file opening; windows that
  reopen where you left them.
- **Readable on every desktop** (b10): pyIrena ships its own light theme, so
  controls no longer become unreadable under a dark OS appearance.
- **Core model objects serialise themselves** (b7). `to_dict`/`from_dict` on
  Unified Fit, Sizes, Modeling, Simple Fits and WAXS Peak Fit make panel
  state, batch configs and HDF5 contents one shape rather than three.
- **A stable, Qt-independent scattering I/O boundary** for companion tools
  (b12): `discover_scattering()` / `load_scattering()` plus the minimal
  `pyirena[qtplot]` extra.
- **Tested on Linux, macOS and Windows** (b10), with Windows-specific
  encoding, path-separator and BOM bugs fixed.
- **Installation troubleshooting**: `pyirena-doctor` diagnoses wrong-environment
  and missing-extra installs instead of reporting "GUI dependencies not
  installed" (b4).

### Changed (compatibility notes since 1.0.1)

- **Minimum Python is now 3.10** (was 3.9); tested on 3.10, 3.11 and 3.13.
- **NumPy floor raised to ≥ 2.0**, matching the `numpy.trapezoid` calls the
  core already made.
- **`pyirena[mcp]` pins `mcp<2`** — pyIrena has not been migrated to the
  mcp 2.x API.
- **The Simple Fits model formerly misspelled "Treubner-Strey" is now
  "Teubner-Strey"**; old result files still load.


## [1.1.0b12] - 2026-09-04

### Added

- Added a stable, Qt-independent scattering I/O boundary for companion tools:
  `ScatteringLocation`, `ScatteringRecord`, `discover_scattering()`, and
  `load_scattering()`. Discovery includes all NXcanSAS entries, slit-smeared
  variants, conventional simple-HDF5 Q/I groups, and read-only text imports.
- Exposed `create_h5xp()` and `write_iq_data()` through `pyirena.io` for
  supported Igor data export without importing an implementation module.
- Added the minimal `pyirena[qtplot]` extra containing only PySide6 and
  PyQtGraph for plotting companion applications such as Bernardyn.

### Fixed

- Copied HDF5 dataset attributes before closing their source file, so units and
  metadata returned from `readGenericNXcanSAS()` remain usable by callers.

## [1.1.0b11] - 2026-09-04

### Changed

- **The four gradient-based fitting tools now use analytic Jacobians.** Unified
  Fit, WAXS Peak Fit, Simple Fits and Modeling hand
  `scipy.optimize.least_squares` / `curve_fit` a closed-form Jacobian instead of
  letting it approximate one by finite differences. The finite-difference path
  spent `N_params + 1` model evaluations per iteration purely probing gradients;
  on a two-level Unified fit that is 418 of 478 model evaluations (87 %), which
  the analytic path removes outright — 60 evaluations for the same
  χ² = 3.80967e4. Wall-clock gain is smaller than the evaluation-count gain,
  because scipy's own trust-region linear algebra is unchanged: measured
  **1.7–2.1× on Unified Fit**, **1.4–2.2× on Modeling** with unified-level or
  diffraction-peak populations, **1.0–2.0× on Simple Fits** depending on model,
  and **~1.8× on large many-peak WAXS Peak Fit** (see the tolerance note
  below).

  Each tool keeps the analytic path on by default and falls back to finite
  differences automatically — both for sub-models where a closed form is
  impractical (LogNormal WAXS peaks; the size-distribution shape parameters in
  Modeling, whose radius grid comes from a numerical CDF inversion) and for any
  unexpected failure, so a fit can never fail because of the Jacobian. Modeling
  is a hybrid: closed-form columns for unified-level, diffraction-peak and
  background parameters, finite differences for the rest, assembled into one
  matrix. The toggle is an implementation detail and is not serialised. Size
  Distribution is unaffected — it solves a regularised linear inverse problem
  and already uses exact gradients.

  Every derivative is verified against high-accuracy central differences, and
  the full `validationData/` suite refits identically: **195/195 checks still
  PASS with no status change**. Two caveats on "results are unchanged":
  Debye-Bueche's `Prefactor` and `Eta` shift by a few percent because only their
  product `Prefactor·Eta²` is determined by the data (that product is conserved
  to 7 significant figures, and `CorrLength` is unchanged); and a Modeling fit
  containing a size-distribution population can land in a different local
  minimum, because an exact gradient traverses a multimodal landscape by a
  different route — neither path is reliably better.

  An exact gradient also fixes fits that finite differences could not do at all.
  A Debye polymer chain fitted over a USAXS range (q from 1e-4 Å⁻¹) previously
  stalled at its starting Rg — the forward model's small-`q²Rg²` cancellation
  noise swamps scipy's default `sqrt(eps)` difference step, so the gradient was
  pure noise. It now converges: Rg = 10.05 against a true 10.0, reduced
  χ² = 1.01 instead of 205.

- **WAXS Peak Fit convergence tolerances stay at 1e-5 — a deliberate decision.**
  The analytic-Jacobian work initially tightened `ftol`/`xtol`/`gtol` to 1e-8 on
  the grounds that exact gradients had made iterations cheap. That is reverted.
  WAXS data with real counting statistics does not support convergence criteria
  far below its own uncertainties — 1e-8 is well past the point where the fit is
  describing noise rather than structure — and Igor Irena has always fitted
  these peaks with comparably loose settings without trouble. Keeping 1e-5 is
  the scientifically honest setting, not a performance compromise.

  Tightening also turned out to be actively unsafe at scale, because tolerance
  and iteration budget are coupled and only one was changed. With
  `maxfev = 10_000` unchanged, a large many-peak fit at 1e-8 exhausts its
  evaluation budget before meeting the criterion; `curve_fit` then raises and
  the handler returns the user's **starting** parameters with
  `success = False`. Measured on 15 Gaussians / 4000 points / 47 free
  parameters: 1e-5 converges in 27 s, while 1e-8 spent 268 s and then failed
  outright. The comment in `fit()` records this so the tolerances are not
  tightened again without also raising `maxfev`.

  With tolerances equal on both paths, the analytic Jacobian's benefit in WAXS
  is a clean speed-up on the same answer: the 15-peak fit above takes 28.7 s
  against 52.7 s on finite differences, and a 20-Gaussian / 8000-point fit
  takes 23.4 s against 61.2 s — both converging to the same reduced χ² as the
  finite-difference path (459.55 on the 20-peak fit, matching to five
  significant figures).

### Fixed

- **`ruff check pyirena` failed CI.** Five lint errors had accumulated in the
  GUI layer as fallout from the b10 theme work: four unsorted import blocks
  (`data_selector/panel.py`, `modeling_panel.py`, `saxs_morph_panel.py`,
  `unified_fit.py`, all from `pyirena.gui.theme` imports being appended rather
  than merged in order) and one unused import (`READONLY_FIELD_CSS` in
  `simple_fits_panel.py`). All are import-ordering only, with no change in
  behaviour. Verified clean under both ruff 0.15 and 0.16 — the workflow
  installs ruff unpinned, so it picks up whatever is current.

## [1.1.0b10] - 2026-09-04

Usability and correctness release, driven by what users hit in practice.

The largest fix is that pyIrena was unreadable for anyone whose operating
system is set to a dark colour scheme: it now ships its own theme instead of
inheriting the desktop's. The fit Q range can be typed again, as in Irena, with
one identical control in Simple Fits, Modeling and Unified Fit. Two parameter
controls that invited meaningless fits — Modeling's Contrast *and* Scale both
fittable, Simple Fits' bounds on a Contrast the Invariant never fits — are
gone.

Underneath, CI now runs the test suite on macOS (arm64) and Windows as well as
Linux, which immediately turned up four platform-specific defects (three in the
tests, one real and shipping), and a new `validationData/` set of 27 synthetic
datasets with exactly known parameters lets pyIrena's mathematics be checked
independently — and compared, quantity by quantity, against Igor Pro Irena.

### Added

- **`pyirena/gui/theme.py` — pyIrena now ships its own light theme.**
  `apply_theme(app)` installs the Fusion style, an explicit `QPalette` and a
  baseline stylesheet on the `QApplication`, so every panel renders the same
  on every platform whatever the desktop's light/dark setting (see *Fixed*
  below for what this cures). It also exports the semantic colour tokens and
  the `accent_button_css` / `chip_button_css` / `soft_button_css` /
  `readonly_field_css` / `status_css` helpers panels should use instead of
  hardcoding hex values — every helper emits a background **and** a text
  colour. `PYIRENA_NATIVE_THEME=1` falls back to the platform style.
- **`pyirena/gui/q_range_ui.py` — one editable "Q range for fit" control,
  shared by every tool.** `QRangeFields` restores Irena's behaviour of
  *typing* the fit limits as well as dragging the graph cursors, and is now
  embedded in Simple Fits, Modeling and Unified Fit so the three cannot drift
  apart. Typed values are validated (positive, distinct), swapped when entered
  the wrong way round, clamped to the loaded data's Q range, and then pushed
  to the cursors, which remain the single source of truth.
- `pyirena/tests/test_gui_theme_contract.py`, `test_gui_q_range_fields.py`
  and `test_gui_fit_flag_rules.py` — regression tests for all of the above.
  The theme contract fails the build if any new inline stylesheet sets a
  background without a text colour.

- **`validationData/` — synthetic data with exactly known parameters, for
  validating pyIrena and for comparing it with Igor Pro Irena.** 27 datasets
  covering Size Distribution, Unified Fit (including slit smearing), Modeling
  (size distributions, Unified levels, diffraction peaks, hard-sphere structure
  factor, mass and surface fractals, Guinier-Porod), Simple Fits, WAXS Peak Fit,
  Data Merge and Data Manipulation. Each is written as a 3-column ASCII file
  with the ground truth in its header, a noise-free `_ideal.dat` companion, and
  an NXcanSAS HDF5 file carrying the exact curve and a JSON of the generating
  parameters under `entry/ground_truth/`.

  The intensities are computed by `validationData/_models.py`, which implements
  every model from the published literature and **does not import pyIrena**, so
  recovering the parameters tests pyIrena's mathematics rather than its
  self-consistency. `generate_validation_data.py` regenerates the files
  deterministically (per-dataset seed derived from the dataset name) and writes
  the `README.md` manifest and `ground_truth.json`.
  `run_validation_report.py` fits every file and writes
  `VALIDATION_RESULTS.md` / `.csv`, with empty Irena columns for a
  side-by-side comparison against the Igor package.
- `pyirena/tests/test_validation_data.py` — asserts that pyIrena's forward
  models agree with the independent implementations to machine precision, that
  every generated file loads, and that a fast subset of the datasets recovers
  its generating parameters.
- `pyirena/tests/test_modeling_structure_factor.py` — regression tests for the
  structure-factor fitting bug below.
- **Igor Pro Irena results are now part of the validation report.**
  `validationData/irena_values.csv` holds the values obtained with the Igor
  package, keyed by dataset and quantity, and `run_validation_report.py` reads
  it to fill the `Irena` and `Irena dev %` columns, compute the agreement
  statistics, and append the commentary in `irena_notes.md`. Because those
  numbers are data in the repository rather than text inside a document,
  regenerating the report never loses them. Over the 98 quantities both
  packages currently report, the median deviation from the known truth is
  0.19 % for pyIrena and 0.31 % for Irena, and the median difference between
  the two packages is 0.04 %.
- `validationData/fill_irena_deviations.py --export-csv` writes a hand-edited
  table's Irena column back into `irena_values.csv`.
- `validationData/md_to_docx.py` renders a results table as a landscape Word
  document (real tables, built-in heading styles) for a manuscript draft.

### Fixed

- **Controls were unreadable when the operating system was set to a dark
  colour scheme.** Reported from a beamline user's machine: the Modeling
  panel's "Fit B/P btwn cursors" and "Fit Flat btwn cursors" buttons appeared
  as blank whitish boxes the user never realised were buttons, and the
  *Population type* / *Distribution* pull-downs could not be read. The cause
  was structural rather than local — several hundred inline stylesheets set a
  `background-color` written for a light scheme without naming a text colour,
  so Qt took the text colour from the *system* palette and painted near-white
  on near-white; controls with no inline style at all were left entirely to
  the platform's dark rendering. Rather than chase every OS × theme
  combination, pyIrena now applies its own palette at startup
  (`pyirena.gui.theme.apply_theme`, wired into all eight GUI entry points),
  which also matches the plots — they are drawn on a hard white background
  throughout. The stylesheets that set a background with no text colour have
  been converted to theme helpers, and a test now blocks new ones.
- **Modeling → size-distribution populations let *Contrast* and *Scale* be
  fitted at the same time.** The two enter the model only as their product, so
  freeing both leaves the least-squares problem with a flat direction: the
  solver wanders and the covariance matrix is singular, making the reported
  uncertainties meaningless. Checking either "Fit" box now clears the other
  (last one ticked wins); holding both fixed is still allowed, and a setup
  loaded with both flags set is normalised to fitting *Scale*.
- **Simple Fits → Invariant showed lo/hi limit fields for *Contrast*.** The
  Invariant is a direct calculation with no least-squares step, so its
  Contrast is never fitted and the bounds did nothing but prompt users to ask
  what they were for. Bound fields (and the lo/hi column headers) now follow
  the "Fit?" box: they are shown only for parameters that can actually be
  fitted, so the background `B`, `P` and flat terms keep theirs.

- **`pyirena-mcp` gave the wrong advice when mcp 2.x was installed.** The
  import guard in `pyirena/mcp/server.py` caught every `ImportError` and
  answered with "install with: pip install pyirena[mcp]". mcp 2.0 replaced
  `mcp.server.fastmcp` with a stub that raises `ModuleNotFoundError` — an
  `ImportError` subclass — so users who already had mcp installed, just the
  wrong major version, were told to install a package they already had, and
  the SDK's own migration hint was swallowed. The guard now reports the
  installed mcp version and the actual fix (`pip install 'mcp>=1.0.0,<2.0'`,
  or upgrade pyirena, which has pinned `mcp<2` since 1.1.0b7). Reached only
  by installs that predate that pin or that upgrade `mcp` by hand — the extra
  itself is unchanged.

- **Modeling: structure-factor parameters were never actually fitted.** The
  hard-sphere radius and volume fraction (and the Interferences `eta`/`pack`)
  were packed into the fit vector under the key group `'sf'`, which the
  *surface fractal* population already used for its own attributes.
  `_unpack_params` matched the surface-fractal branch first, so each optimiser
  step wrote those values onto the population as stray attributes instead of
  into `pop.sf_params`. The model never saw them change: they contributed
  nothing to chi-squared and `fit()` returned them at exactly their starting
  values, with no warning. Structure-factor parameters now use their own key
  group `'sfp'`. Found while building `validationData/`; on the hard-sphere
  validation dataset the fit went from reduced chi-squared 397 (all four
  structure-factor and distribution parameters wrong by 8-40 %) to 1.15 with
  every parameter within 0.2 % of truth.

- **A UTF-8 BOM silently dropped the first data point of a text file.** This is
  not Windows-specific — it affected every platform. `readTextFile` opened the
  file with the platform default encoding, so a byte-order mark from a Windows
  editor arrived as a `\ufeff` glued onto the first Q value. That failed the
  `float()` probe used to find where the header ends, so the line was
  classified as another header row and skipped: a four-point file loaded as
  three points, with no warning. The file is now read as `utf-8-sig`, which
  consumes the BOM.
- **Text data files whose header was not UTF-8 failed to load on Windows.**
  Same root cause from the other direction: an older instrument header carrying
  a cp1252/ISO-8859-1 `Å` or degree sign is not valid UTF-8, and the platform
  default on Windows is cp1252, so the same file loaded on Linux and raised
  `UnicodeDecodeError` on Windows. `readTextFile` now tries UTF-8 and falls
  back to latin-1, which cannot fail — only header text is affected either way,
  since the numeric columns are ASCII in all of these encodings. The fix had to
  cover two reads, not one: `np.loadtxt` defaults to `encoding=None`, meaning
  the platform default, so it would have raised on the second pass over a file
  the first pass had just recovered. Both now work from a single decoded read.
- **Dropped files arrived with the wrong path separators on Windows.**
  `paths_from_mime` returned whatever Qt handed it, and `QUrl.toLocalFile()`
  gives `C:/Users/...` with forward slashes. The plain-text branch was worse: it
  stripped seven characters off `file://` unconditionally, turning Explorer's
  `file:///C:/Users/x` into `/C:/Users/x`, which is not a path at all. A dropped
  file therefore compared unequal to the same file found by a folder scan.
  Real file URLs now go through `urllib.request.url2pathname` (which knows
  `/C:/x` means `C:\x`, and un-escapes `%20`), the malformed `file://C:\x`
  spelling that some applications send is taken literally, and every result is
  normalised to native separators.
- **`export_fit_results` crashed on Windows.** It writes `Å` and `cm⁻¹` into a
  plain text file without specifying an encoding, so it raised
  `UnicodeEncodeError` under cp1252. It has no test, which is why CI did not
  flag it; found by inspection while fixing the above.

### Changed

- **Tests run on three operating systems.** `test` became a deliberately sparse
  matrix rather than the full cartesian product: Linux sweeps Python
  3.10/3.11/3.13, and macOS-arm64 and Windows each get 3.12. Version bugs need
  the Python sweep; platform bugs need the platform, and the extra six
  combinations of a full matrix buy almost nothing. `test-gui` runs on all
  three under `QT_QPA_PLATFORM=offscreen`, with the `apt-get` step gated to
  Linux. Both jobs carry a 45-minute timeout.
- **Three test-only Windows failures fixed.** Six `Path.read_text()` calls
  reading UTF-8 content — a source file, and `.itx` exports containing `cm⁻¹` —
  decoded as cp1252 and raised; the ITX *writer* was already correct. Four more
  in `test_report_buttons` were one `χ²` away from the same failure and got the
  same treatment. Separately,
  `test_run_fit_walks_to_distant_optimum_by_default` called
  `add_unified_level()` on top of the level `select_model()` already creates,
  leaving a second, entirely default level free in the fit — a duplicate power
  law the data cannot distinguish from the first, whose G and B drift onto
  their `1e-10`/`1e-20` lower bounds. Whether they land *on* a bound or just
  above it is numerical noise that differs by platform, which is why only
  Windows reported them as pinned. The spurious level is gone and the
  one-level assumption is now asserted rather than assumed.

## [1.1.0b9] - 2026-08-21

Performance release: Modeling fits run about three times faster with unchanged
results. Also fixes a class of silently wrong Modeling fit, in which a
population that could not be evaluated was dropped from the model without
failing the fit.

### Fixed

- **Modeling: a population that could not be evaluated was silently dropped,
  and the fit reported success without it.** `total_intensity` caught every
  per-population exception, issued a `RuntimeWarning`, and left that population
  out of the sum. During a fit that happens identically at every trial, so the
  optimiser converged happily on the remaining terms and returned a result
  whose reduced chi-squared belonged to a different model than the one asked
  for — observed for real at reduced chi-squared 62 where the intact model
  gives 0.51. Fitting is now strict: `fit()` evaluates the model once at the
  starting parameters and raises `PopulationEvaluationError` if it cannot,
  and a population that fails only at some trial point makes that trial
  maximally bad (so the optimiser walks away from it) and is counted in the
  fit warnings instead of vanishing. Display paths are unchanged — they still
  warn and skip, so a half-typed parameter cannot take a GUI panel down.
- **Modeling: the core-shell G matrix lost precision in the Guinier regime.**
  Introduced by the vectorisation below and caught before release. Factoring
  the r³ out of Δρ·V(r)·f_sph(qr) leaves `sin(x) − x·cos(x)`, which cancels
  catastrophically as x → 0: ~5 significant digits by qr = 1e-5 and 32 % error
  by qr = 1e-7, a regime small particles at USAXS Q genuinely reach. The
  vectorised builders now switch to the exact Taylor expansion below the same
  qr = 1e-3 threshold `_sphere_amplitude` uses, restoring full precision (the
  branch is skipped entirely, at no cost, when no element is small).
- **MCP: every image tool failed to return its picture.** `pyirena_plot_iq`,
  `pyirena_plot_parameter_trend` and all eight `pyirena_ctrl_*_image` tools
  rendered and saved their PNG, then died on the way back to the client with
  `Unable to serialize unknown type: ... fastmcp.utilities.types.Image`.
  Cause: FastMCP ≥ 1.10 derives a *structured output* JSON schema from a
  tool's return annotation and serialises the return value against it, and
  the `Image` object is not JSON-serialisable — so the annotation
  `-> list[Any]` silently opted every image tool into a code path that could
  never succeed. Images belong in the unstructured content array; image tools
  are now registered through `_image_tool()`, which disables structured
  output (and falls back to the plain decorator on mcp < 1.10, where the
  feature does not exist). Claude Desktop / claude.ai reported this as "the
  plot tool has a serialization bug but it did save the PNG to disk".

### Added

- **MCP: image tools now report the PNG's path as well as the image.** Each
  returns a text block (`"<label>\nPNG saved to: /abs/path.png"`) followed by
  the inline `image/png` content block. Previously the control-API image
  tools returned base64 only, so a client that does not render inline images
  — AnythingLLM in some modes, plain agent loops — had no way to reach the
  plot and told users no file existed. The control API dicts gained a
  matching `image_path` key alongside `image_base64`.
- **`PYIRENA_PLOT_CACHE` now also directs control-API fit images.** It
  previously only affected `plot_iq` / `plot_parameter_trend`; fit images
  were hardwired to `<tempdir>/pyirena-ctrl` (still the default). Point it at
  a browsable folder when the user will want to open the PNGs themselves.

### Changed

- **Modeling fits are ~3x faster, with identical results.** A user-supplied
  core-shell model (1406 Q-points, 199 radius bins, 10 free parameters) went
  from 50 s to 16 s for a global fit and 39 s to 12 s for 50 Monte-Carlo
  passes, on the same 18-core machine. Every fitted parameter and every MC
  uncertainty is unchanged; the full test suite is unchanged. Four
  independent causes, all in the per-iteration hot path:
  - `schulz_zimm_cdf` and `gauss_cdf` called `scipy.stats.gamma.cdf` /
    `norm.cdf`, whose `rv_continuous` wrapper costs ~14 µs per *scalar* call.
    `generate_radius_grid` makes a few hundred of those per fit iteration
    while walking the CDF out to the distribution tails. They now call
    `scipy.special.gammainc` / `ndtr` — the functions scipy itself dispatches
    to — for bit-identical results at ~1/33 the cost.
  - `total_intensity` rebuilt each population's radius grid a second time
    after `calculate_pop_intensity` had already built it. The grid is now
    returned by the same call, and memoised on the distribution parameters
    that define it, so a Jacobian step that leaves them alone reuses it.
  - The core-shell sphere G matrix (all three of `cs_sphere_by_core`,
    `_by_shell` and `_by_total`) was assembled with a Python loop over radius
    bins. It is now one vectorised expression; where the shell thickness is
    constant across bins, half the sin/cos work is recovered by angle
    addition.
  - Each population's I(Q) is memoised during a fit on the full set of
    parameters it depends on, so a step in the background or in another
    population skips it entirely.
- **The remaining two per-bin loops are vectorised too**, on the same identity,
  through a shared `_shell_step_amplitude` helper: the core-shell-shell sphere
  over both axes at once, and the core-shell spheroid over its 50
  Gauss-Legendre orientations rather than its (typically 200+) radius bins, so
  each iteration is one full Q x bin evaluation and there are far fewer of
  them. Both are ~1.8x faster and agree with the loops they replace to 4e-14
  relative-to-peak. This is what a slow core-shell-shell or core-shell-spheroid
  fit was waiting on; the spheroid is still much the most expensive form factor
  because of the orientational average.
- **The fit-time G-matrix cache now keys on the Q grid as well.** It keyed
  only on the radius grid and form-factor parameters, which is why
  `total_intensity_ideal` had to avoid it: with slit smearing on, a fit
  evaluates the model on both the data grid and the smearer's extended grid.
  Mismatched lengths raised; two same-length grids would have silently reused
  the wrong matrix.
- **Documented where parallel fitting stops paying off** (`docs/modeling_gui.md`).
  Both the global-fit and MC-uncertainty `cores` settings flatten out around
  8–10 workers — model evaluation is bound by memory bandwidth over the G
  matrix, not by CPU — so values above ~10 cost resources without saving time.
- **PNG rendering for the control API lives in one place**,
  `pyirena/api/control/_images.py` (`render_png`, `image_cache_dir`),
  replacing five near-identical private copies in `unified_fit.py`,
  `sizes.py`, `simple_fits.py`, `modeling.py` and `waxs_peakfit.py`.

## [1.1.0b8] - 2026-08-13

### Fixed

- **Unified Fit: "press Fit repeatedly" caused by the moving limits window.**
  The panel keeps each parameter's fit limits centred on its current value
  (B: value×0.2 … ×5; see `fix_limits`). When the optimum lies far outside
  that window — changing P by ~1 moves the matching B by orders of magnitude
  in the USAXS range — a single bounded fit converged correctly *within its
  box* but ended with B pinned at a limit, and each further Fit press only
  moved B another factor of 5 (this was never a tolerance, restart-loop, or
  step-size problem). The Fit button now uses
  `UnifiedFitModel.fit_with_limit_walking()`: while a fitted parameter ends
  pinned at a panel auto-limit, the limits are recentred on the fitted values
  and the fit rerun (up to 4 walks), so one press reaches the interior
  optimum. The auto-limits policy itself now lives once, in
  `panel_auto_limits()` (core), used by both the panel and the walker. Batch
  fits with explicit config limits are unchanged — those limits are treated
  as deliberate constraints. The agent-control `run_fit` uses the same
  walking by default (`walk_limits=True`, opt out with False) and now
  returns `pinned_parameters` and `warnings` so a bound-stopped fit is
  visible instead of silently wrong.
- **Modeling & Unified Fit: micro-stepping local fits (`diff_step=1e-3`).**
  `least_squares` now uses a finite-difference Jacobian step of 0.1 % of each
  parameter's value — the scipy analogue of Igor's per-parameter epsilon,
  which Igor Irena needed for the same reason. With scipy's default step
  (√eps ≈ 1.5e-8 relative), the locally exact Jacobian made TRF advance in
  micro-steps on stiff, strongly correlated surfaces (a mean size of 5×10⁴
  fitted alongside peak widths of 4×10⁻³; B and P in the USAXS power-law
  region), exhausting the ~2400-evaluation budget every round without meeting
  any tolerance. With the 0.1 % secant step the same fits converge in ~30–90
  evaluations to an equal or better χ² (measured 20–35× fewer evaluations on
  two user datasets; a pathological 10-minute fit now finishes in ~10 s).
  The Unified Fit restart loop also gets the relative (1e-4·χ²) stop
  threshold described below.
- **Modeling: runaway 10-minute fits on misspecified models.** When a model
  structurally could not reach the data (e.g. a strong low-Q upturn with
  `scale` pinned at its limit), every local-fit round exhausted its evaluation
  budget in a nearly flat χ² valley and the internal restart loop ran all five
  rounds — χ² improvements of a few counts on χ² ~ 10⁶ passed the old
  `1e-8·χ²` restart test. The threshold is now `1e-4·χ²` (relative), which
  stops such fits after two rounds while leaving genuine "press Fit again"
  improvements (orders of magnitude larger) unaffected. Reproduced with
  `testData/bad1.h5`.

### Added

- **Modeling: post-fit diagnostics** on `ModelingResult.fit_warnings`, shown
  in the GUI status line, logged by `batch.fit_modeling`, and returned as
  `warnings` by the agent-control `modeling_run_fit`. Warns when a fitted
  parameter ends pinned at a fit limit, when the fit stopped at its evaluation
  budget without converging, when a lognormal `sdeviation` exceeds 3
  (unphysically broad — usually mimicking a power law), and when reduced χ²
  is above 10³ (model cannot reach the data).

## [1.1.0b7] - 2026-08-10

A consolidation release, implementing the feature-parity review in issue #13.
The theme is that behaviour a user thinks of as belonging to "a table", "a
graph" or "a file browser" is now implemented once and available everywhere,
rather than per panel: tables copy and sort, graphs export the same five ways,
browsers accept dropped files, windows reopen where you left them, and all five
fitting tools write Markdown reports and can be driven by an agent.  Nine
duplicated implementations were removed along the way, three of which had
already drifted apart.  See `PLAN.md` for the full disposition, including what
was deliberately *not* done.

### Added

- **File browsers share their remaining logic** (feature parity review item
  U4, issue #13 — scoped down deliberately, see below).
  - The Data Selector was the last browser keeping its own extension lists and
    its own `os.listdir` loop; it now uses `core.file_types`, so all four
    browsers answer the awkward questions the same way (skip sub-directories,
    match case-insensitively, treat an unreadable folder as empty rather than
    raising — the Data Selector used to show an error message instead).
    New `files_with_extensions()` serves a browser whose dropdown is not
    `FILE_TYPES`; `files_in_folder()` now delegates to it, so there is one
    listing implementation.
  - The post-drop "switch folder, select what was dropped" step, written three
    times during the drag-and-drop work, is now `select_dropped_in_list()`.
  - `test_gui_browser_contract.py` grew four checks: every browser accepts
    drops, none walks a folder by hand, none hard-codes extensions, and the
    list browsers share the drop-selection step.
  - **A shared `FileBrowserWidget` was considered and deliberately not built.**
    After U1–U3 and the drag-and-drop work, what remains duplicated across the
    four browsers is widget assembly with no logic in it — and they differ in
    ways a common widget would have to be configured around regardless: the
    Data Explorer is a tree with lazy sub-folder expansion, Data Merge shows
    two linked instances, Data Manipulation adds a context menu, and the Data
    Selector alone offers text files and the convert-on-load path. The reason
    is recorded in `docs/developer_adding_features.md` so the question does not
    get re-opened from scratch.
- **Core model objects serialise themselves — `to_dict`/`from_dict`**
  (feature parity review item U6, issue #13).  Added to Unified Fit, Modeling
  and WAXS Peak Fit, joining Sizes and Simple Fits; the callers that used to
  build those dicts by hand now go through the model.
  - **Unified Fit.**  `UnifiedLevel`/`UnifiedFitModel` gained the pair, plus
    `UnifiedLevel.from_panel_params()` — the one place that translates the
    panel's historical vocabulary (`RgCutoff` → `RgCO`, `correlated` →
    `correlations`, `estimate_B` → `link_B`, `Rg_low`/`Rg_high` → `Rg_limits`),
    which is baked into saved setups and batch configs and so cannot be
    renamed.  It replaced **eight** hand-written conversions: six in the panel
    (graph, both fit branches, undo, both Monte-Carlo branches) and two in
    `batch/unified.py`.  The flags `with_limits` / `with_links` /
    `with_fit_flags` keep each call site's exact behaviour — notably the graph
    and undo paths, which must *not* enable `link_B`, since that recomputes B
    from G, Rg and P and would change the drawn curve.
  - **Modeling.**  A `_SerialisableDataclass` mixin walks
    `dataclasses.fields()`, so all six population types share one
    implementation and a new parameter is serialised the moment it is declared.
    `population_from_dict()` dispatches through a `POPULATION_CLASSES`
    registry; an unrecognised `pop_type` (a file from a newer pyIrena) loads as
    a size distribution with a warning instead of raising.  This removed three
    near-identical six-branch deserialisers — in `batch/modeling.py`, in
    `gui/modeling_panel.py`, and a fourth that was already dead — and the
    panel's hand-written serialiser, which was verified field-by-field to
    produce exactly the same dict for all six types before being deleted.
  - **WAXS Peak Fit.**  `WAXSPeakFitModel.to_dict`/`from_dict` cover background
    shape, background parameters and peaks; an unknown background shape falls
    back with a warning, and `to_dict` deep-copies so a caller cannot edit the
    live peak list by accident.
  - Two historical quirks were preserved on purpose rather than tidied away: a
    size-distribution population whose config or state omits `enabled` still
    loads *disabled* (batch and GUI both did this; every other population type
    defaults to enabled), and the panel's size-dist restore still falls back to
    the log-normal defaults for a missing `dist_params`.
- **SAXS Morph result files now carry the panel setup** (feature parity review
  item U10, issue #13).  `save_saxs_morph_results(..., setup_state=...)` embeds
  the `_pyirena_config` attribute the other five fitting tools already wrote,
  the panel gained a **Load Setup from File…** button, and batch runs embed
  their config section too — so a batch-produced result opens in the GUI with
  every control set.  The physics scalars were already stored; what was missing
  was the input mode, the cursor Q range and the background pre-fit windows.
  - The per-tool policy is now written down in
    `docs/HDF5_NxcanSAS_structure.md`: six tools embed the setup, Data Merge
    and Data Manipulation record **NXprocess provenance** instead (they run
    inside reduction pipelines driven by their own JSON, where a GUI session
    never existed), and Fractals stores every growth parameter as explicit
    datasets that `load_fractal_aggregate()` rebuilds.  A test fails the build
    if a tool is missing from that table.
- **One Qt import point at last** (feature parity review item U7, issue #13).
  The invariant said every Qt import goes through `pyirena/gui/_qt.py`; in
  practice there were three shims — `gui/_qt.py`, `gui/data_selector/_qt.py` and
  an inline one in `slit_smearing_ui.py` — plus 22 local
  `try: PySide6 / except: PyQt6` blocks scattered through ten modules.
  - `gui/data_selector/_qt.py` is **removed** and its five importers now use
    `pyirena.gui._qt`; every local block is gone; `QMimeData` was the only name
    the shim was missing.
  - Three of those blocks had a **stale `PyQt5` third branch** (`unified_fit`,
    `modeling_panel`, `ai_advisor`), which on a PyQt6-only install would have
    raised `ModuleNotFoundError` from inside a paint or a message box rather
    than falling back.
  - `pyirena/tests/test_gui_qt_contract.py` keeps it fixed: it fails the build
    on a direct binding import anywhere in the package (tests included), on a
    second `_qt.py` in a subpackage, on a name imported from the shim that it
    does not define, and on a name defined there but absent from `__all__`.
    It reads the source, so it runs on a machine with no Qt installed.
  - No behaviour change intended; the one visible edit is that Unified Fit's
    cursor paint code now uses the shared `Qt` enum instead of a local alias
    that shadowed it as `QtCore`.
- **Drag a file onto pyIrena to open it** (feature parity review item A6,
  issue #13).  New `pyirena/gui/file_drop.py`; wired into the Data Selector,
  every fitting panel (through the shared `DataFileLoaderRow`, so Unified Fit,
  Sizes, Simple Fits, Modeling, WAXS and SAXS Morph all gained it at once), the
  Data Explorer file tree, and the Data Merge and Data Manipulation file lists.
  - Dropping a **folder** offers the data files inside it, one level deep — a
    measurement directory can be dragged in whole.
  - The drag is **refused while it is still over the window** if nothing in it
    is openable, so the cursor answers before the user lets go.  Files pyIrena
    cannot read are dropped from the list rather than failing later.
  - Both drop payloads are read: Finder/Explorer send file URLs, some
    applications send plain text paths.
  - Text files still go through the clean-and-convert path, so a dropped
    `.dat` behaves exactly like one opened from the dialog.
  - Installed through an event filter, so no existing widget had to be
    subclassed; `collect_dropped_paths()` is pure and covered by 23 tests.
- **Windows reopen where you left them** (feature parity review item A9,
  issue #13).  New `pyirena/gui/window_state.py` saves size, position and
  splitter sizes for every pyIrena window — so the width you gave the left
  control panel survives a restart, not just the window frame.
  - **The failure mode this is really about is a changed display setup.**  A
    position saved on a monitor that is no longer attached would reopen the
    window somewhere it can be neither seen nor dragged back.  Saved geometry
    is checked against the screens that exist *now*: a window that would be
    unreachable is moved onto the nearest real screen and shrunk only if it no
    longer fits, and if it cannot be placed sensibly the tool falls back to its
    default size.  The whole policy is the pure function `resolve_geometry()`,
    tested against removed monitors, negative screen origins,
    larger-than-screen windows and a grid of off-screen positions.
  - **Per-tool reset, the Irena gesture:** hold **Shift** while clicking a
    tool's button and that tool forgets its position and opens centred at its
    default size, control-panel width included.  Because panels are built once
    and reused, `install_window_state()` records each window's coded default
    *before* applying anything saved, so the gesture works on the tenth launch
    as well as the first; a panel's separate graph window is reset with it.
    `pyirena/tests/test_window_state.py` reads the launcher source and fails
    the build if a new tool is wired up without it.
  - Pane widths are recorded and applied **at the first layout**, not in the
    constructor: a panel sets its splitter while that splitter still has no
    width, so anything read or written before Qt's first layout pass is a
    placeholder.  They are stored as *fractions*, since a reset re-centres the
    window at a different size than the splitter currently has.
  - Geometry is only ever applied to real windows (`is_top_level`).  Unified
    Fit's and WAXS's `…GraphWindow` classes are the right-hand *pane* of their
    panel, and setting a pane's geometry fought the layout — which surfaced as
    the control panel coming back much too wide after a reset.  Both are no
    longer registered, and the guard makes the same mistake harmless in future.
  - **Global escape hatch:** hold **Shift** while launching `pyirena-gui`, or
    set `PYIRENA_RESET_WINDOWS=1`, to discard every saved geometry and open all
    windows at their defaults.
  - Minimised and full-screen windows are not saved.
  - Geometry is stored in its own `window_geometry.json` next to the pyIrena
    state file, deliberately: each panel holds its own `StateManager` and
    `save()` rewrites the whole file from that copy, so geometry kept there was
    clobbered by whichever panel closed last.
  - `GEOMETRY_PERSISTENCE_ENABLED = False` in `window_state.py` disables the
    feature outright, and the placement tunables (minimum visible area,
    title-bar grab height, edge margin) are named constants at the top of the
    module for fine-tuning.
- **WAXS Peak Fit is now agent-drivable — U9 is complete** (feature parity
  review item U9, issue #13). `pyirena/api/control/waxs_peakfit.py` adds 18
  tools: background choice (adaptive SNIP and friends, or fitted polynomials),
  peak finding, add/remove/list peaks, per-peak shape and parameters with fit
  flags and bounds, the fit with its three weighting modes, results, an image
  and save to NXcanSAS — exposed over MCP as `pyirena_ctrl_waxs_*`.
  - **`find_waxs_peaks` is the entry point**: it runs pyIrena's peak detector
    over the fit Q range and creates peaks with position, amplitude and width
    already close, so an agent starts from the data rather than from guesses.
    `prominence_frac` trades spurious detections against missed shoulders.
  - Every peak reports its **integrated area** (derived from the fitted shape,
    with a propagated uncertainty) — usually the quantity a WAXS question is
    actually about.
  - The fit image draws each peak separately over the background, which is how
    a peak that has drifted onto its neighbour becomes visible.
  - The MCP surface grows from 101 to 119 tools. **All five fitting tools —
    Unified Fit, Sizes, Simple Fits, Modeling and WAXS Peak Fit — now have a
    control surface.**
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

- **`pip install pyirena[mcp]` broke against the newly released mcp 2.0.**
  The 2.x SDK removed `mcp.server.fastmcp`, which `pyirena/mcp/server.py` is
  built on, so a fresh install picked up 2.x and failed at import with a
  misleading "install with pyirena[mcp]" message. The extra is pinned to
  `mcp>=1.0.0,<2.0` until the server is migrated to the 2.x API.
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

---

Entries for **1.0.1 and earlier** live in
[docs/CHANGELOG_archive.md](docs/CHANGELOG_archive.md).
