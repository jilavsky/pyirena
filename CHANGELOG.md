# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
