# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
  HDF5 Viewer (already had ITX; menu entry unchanged).

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
