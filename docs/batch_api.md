# pyIrena Batch Fitting API

`pyirena.batch` provides a headless (no GUI required) fitting API for scripting,
automation, and integration with other analysis packages.  It reads the same
configuration files produced by the GUI's **Export Parameters** button, so there
is no need to write configuration code by hand.

---

## Contents

1. [Quick start](#quick-start)
2. [Configuration file format](#configuration-file-format)
3. [API reference — `fit_unified`](#fit_unified)
4. [API reference — `fit_sizes`](#fit_sizes)
5. [API reference — `fit_simple`](#fit_simple)
6. [API reference — `fit_waxs`](#fit_waxs)
7. [API reference — `fit_modeling`](#fit_modeling)
8. [API reference — `fit_saxs_morph`](#fit_saxs_morph)
9. [API reference — `merge_data`](#merge_data)
10. [API reference — `average_data`](#average_data)
11. [API reference — `igor_to_nexus`](#igor_to_nexus)
12. [API reference — `fit_pyirena`](#fit_pyirena)
13. [Return structures](#return-structures)
14. [Error handling](#error-handling)
15. [Batch processing patterns](#batch-processing-patterns)
16. [Extending with new tools](#extending-with-new-tools)

---

## Quick start

```python
from pyirena.batch import fit_unified, fit_sizes, fit_simple, fit_waxs, fit_modeling, fit_saxs_morph, fit_pyirena

# Fit one file with Unified Fit, save results to NXcanSAS
result = fit_unified("sample.h5", "pyirena_config.json")
if result and result['success']:
    p = result['parameters']
    print(f"chi² = {p['chi_squared']:.4g}")
    print(f"Level 1 Rg = {p['levels'][0]['Rg']:.2f} Å")

# Fit a size distribution (reads 'sizes' section from config)
result = fit_sizes("sample.h5", "pyirena_config.json")
if result and result['success']:
    print(f"Volume fraction = {result['volume_fraction']:.4g}")

# Fit a simple analytical model (no config file needed)
result = fit_simple("sample.h5",
                    config={'model': 'Guinier', 'params': {'I0': 1.0, 'Rg': 50.0}})
if result and result['success']:
    print(f"Rg = {result['params']['Rg']:.2f} ± {result['params_std']['Rg']:.2f} Å")

# Fit WAXS diffraction peaks (reads 'waxs_peakfit' section from config)
result = fit_waxs("waxs_data.h5", "pyirena_config.json")
if result and result['success']:
    for pk in result['peaks']:
        print(f"Q0={pk['Q0']['value']:.4f}  FWHM={pk['FWHM']['value']:.4f}")

# Fit a parametric size-distribution / modeling model (reads 'modeling' section from config)
result = fit_modeling("sample.h5", "pyirena_config.json")
if result and result['success']:
    print(f"χ²/dof = {result['result'].reduced_chi_squared:.4g}")

# Generate a 3D voxelgram from SAXS data (reads 'saxs_morph' section from config)
result = fit_saxs_morph("sample.h5", "pyirena_config.json")
if result and result['success']:
    r = result['result']
    print(f"χ²/dof = {r.reduced_chi_squared:.4g},  φ_actual = {r.phi_actual:.4f}")
    voxelgram = r.voxelgram   # uint8 ndarray, shape (N, N, N), values 0 or 1

# Run ALL configured tools on one file (one function call)
results = fit_pyirena("sample.h5", "pyirena_config.json")

# Fit many files, collect summaries
from pathlib import Path
data_files = sorted(Path("data/").glob("*.h5"))
config = "pyirena_config.json"

summaries = []
for f in data_files:
    r = fit_unified(f, config, save_to_nexus=True)
    if r and r['success']:
        summaries.append({
            'file': f.name,
            'Rg_1': r['parameters']['levels'][0]['Rg'],
            'chi2': r['parameters']['chi_squared'],
        })
```

---

## Configuration file format

Configuration files are plain JSON, written by the GUI's **Export Parameters**
button (or by hand).  The structure is designed to grow as new tools are added:
each tool occupies its own top-level key alongside the required `_pyirena_config`
header.

### Annotated example

```json
{
  "_pyirena_config": {
    "file_type": "pyIrena Configuration File",
    "version": "0.1.0",
    "created": "2026-02-17T10:30:00",
    "modified": "2026-02-17T14:22:05",
    "written_by": "pyIrena 0.1.0"
  },

  "unified_fit": {
    "num_levels": 2,
    "cursor_left": 0.003,
    "cursor_right": 0.45,
    "no_limits": false,
    "update_auto": false,
    "display_local": false,

    "background": {
      "value": 1e-6,
      "fit": false
    },

    "levels": [
      {
        "level": 1,
        "G":    { "value": 1e10, "fit": true,  "low_limit": 1e8,  "high_limit": 1e12 },
        "Rg":   { "value": 100,  "fit": true,  "low_limit": 10,   "high_limit": 1000 },
        "B":    { "value": 1e6,  "fit": true,  "low_limit": 1e4,  "high_limit": 1e8  },
        "P":    { "value": 4.0,  "fit": false, "low_limit": 0,    "high_limit": 6    },
        "ETA":  { "value": 0,    "fit": false, "low_limit": 0.1,  "high_limit": 1e6  },
        "PACK": { "value": 0,    "fit": false, "low_limit": 0,    "high_limit": 16   },
        "RgCutoff": 0.0,
        "correlated": false,
        "estimate_B": false,
        "link_rgco": false
      },
      {
        "level": 2,
        "G":    { "value": 1e8,  "fit": true,  "low_limit": 1e6,  "high_limit": 1e10 },
        "Rg":   { "value": 10,   "fit": true,  "low_limit": 1,    "high_limit": 100  },
        "B":    { "value": 1e4,  "fit": true,  "low_limit": 1e2,  "high_limit": 1e6  },
        "P":    { "value": 3.5,  "fit": false, "low_limit": 0,    "high_limit": 6    },
        "ETA":  { "value": 0,    "fit": false, "low_limit": 0.1,  "high_limit": 1e6  },
        "PACK": { "value": 0,    "fit": false, "low_limit": 0,    "high_limit": 16   },
        "RgCutoff": 100.0,
        "correlated": false,
        "estimate_B": false,
        "link_rgco": true
      }
    ]
  }
}
```

### `_pyirena_config` header fields

| Field | Type | Description |
|-------|------|-------------|
| `file_type` | str | Always `"pyIrena Configuration File"` — used to validate the file |
| `version` | str | pyIrena version that originally created the file |
| `written_by` | str | `"pyIrena <version>"` — version that last wrote to this file |
| `created` | str | ISO-8601 timestamp of file creation |
| `modified` | str | ISO-8601 timestamp of last write |

### `unified_fit` fields

| Field | Type | Description |
|-------|------|-------------|
| `num_levels` | int | Number of active structural levels (1–5) |
| `cursor_left` | float | Lower Q bound used for fitting (Å⁻¹) |
| `cursor_right` | float | Upper Q bound used for fitting (Å⁻¹) |
| `no_limits` | bool | If true, ignore user limits and use wide defaults |
| `background.value` | float | Flat background value (cm⁻¹) |
| `background.fit` | bool | Whether to fit the background |
| `levels` | list | One entry per level (see below) |

### Per-level fields

Each entry in `levels` describes one structural level.

| Field | Type | Description |
|-------|------|-------------|
| `level` | int | Level number (1-based, informational only) |
| `G.value` | float | Guinier prefactor (cm⁻¹) |
| `G.fit` | bool | Fit G? |
| `G.low_limit` | float | Lower bound for G during fitting |
| `G.high_limit` | float | Upper bound for G during fitting |
| `Rg.value` | float | Radius of gyration (Å) |
| `Rg.fit` | bool | Fit Rg? |
| `Rg.low_limit` / `Rg.high_limit` | float | Bounds for Rg |
| `B.value` | float | Porod/power-law prefactor |
| `B.fit` | bool | Fit B? |
| `P.value` | float | Power-law slope |
| `P.fit` | bool | Fit P? |
| `ETA.value` | float | Correlation distance (Å) — used when `correlated=true` |
| `ETA.fit` | bool | Fit ETA? |
| `PACK.value` | float | Packing factor — used when `correlated=true` |
| `PACK.fit` | bool | Fit PACK? |
| `RgCutoff` | float | Cutoff radius linking to smaller level (Å); always 0 for level 1 |
| `correlated` | bool | Enable Born-Green correlation function for this level |
| `link_rgco` | bool | Automatically link RgCutoff to Rg of the level below |
| `estimate_B` | bool | Estimate B from mass fractal assumption |

---

## `fit_unified`

```python
from pyirena.batch import fit_unified

result = fit_unified(
    data_file,          # str or Path — input SAS data file
    config_file,        # str or Path — pyIrena JSON config file
    save_to_nexus=True  # bool — write results to NXcanSAS HDF5 file
)
```

### What it does, step by step

1. **Loads the config file** — validates the `_pyirena_config` header and reads
   the `unified_fit` section.
2. **Loads data** — detects format by file extension:
   - `.dat` / `.txt` → column-delimited text (Q, Intensity, Error)
   - `.h5` / `.hdf5` → NXcanSAS HDF5 via `readGenericNXcanSAS`
3. **Applies Q range** — masks data to `[cursor_left, cursor_right]` from config.
   If either cursor is absent, the full Q range is used.
4. **Builds `UnifiedFitModel`** — constructs `UnifiedLevel` objects from config
   values, fit flags, and limits.  Respects `no_limits` flag.
5. **Runs the fit** — calls `model.fit(q, intensity, error)`.
6. **Saves results** *(if `save_to_nexus=True`)* — writes to NXcanSAS HDF5:
   - NXcanSAS input → results appended to the same file
   - Text input → new file `<stem>_NX.h5` created alongside the data
7. **Returns** a result dict (see [Return structures](#return-structures)).

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_file` | `str` or `Path` | required | Path to SAS data file |
| `config_file` | `str` or `Path` | required | Path to pyIrena JSON config |
| `save_to_nexus` | `bool` | `True` | Save fit results to NXcanSAS HDF5 |

### Returns

`dict` on success (see [Return structures](#return-structures)), or `None` if a
fatal error occurs before fitting (file not found, invalid config, data unreadable).
A failed fit (optimizer did not converge) still returns a dict with `success=False`.

---

## `fit_sizes`

```python
from pyirena.batch import fit_sizes

result = fit_sizes(
    data_file,           # str or Path — input SAS data file
    config_file,         # str or Path — pyIrena JSON config with 'sizes' section
    save_to_nexus=True,  # bool — write results to NXcanSAS HDF5 file
    with_uncertainty=False,
    n_mc_runs=10,
)
```

### What it does

Reads the `sizes` section from the config file and runs the size-distribution
inversion (same four methods as the GUI: Regularization, MaxEnt, TNNLS, Monte Carlo).
Saves the volume size distribution, number size distribution, cumulative
distributions (volume and number), residuals, and all scalar parameters to
`entry/sizes_results` in the HDF5 file.

### Returns

`dict` with keys `success`, `message`, `input_file`, `config_file`, `output_file`,
and all the scalar values produced by the inversion (same as `load_sizes_results()`).

---

## `fit_simple`

```python
from pyirena.batch import fit_simple

result = fit_simple(
    data_file,            # str or Path — input SAS data file
    config,               # dict or SimpleFitModel
    with_uncertainty=False,
    n_mc_runs=50,
    q_min=None,           # float or None — lower Q bound (Å⁻¹)
    q_max=None,           # float or None — upper Q bound (Å⁻¹)
    verbose=True,
)
```

### What it does

Fits a single analytical model to one SAS file and, for NXcanSAS HDF5 inputs,
saves the results to `entry/simple_fit_results`.  The `config` argument accepts
either a `SimpleFitModel` instance or a plain dict:

```python
# Minimal dict config
result = fit_simple("sample.h5",
                    {'model': 'Guinier', 'params': {'I0': 1.0, 'Rg': 50.0}})

# With bounds and complex background
result = fit_simple("sample.h5", {
    'model':          'Sphere',
    'params':         {'Scale': 1e6, 'R': 100.0},
    'limits':         {'Scale': [1e3, 1e9], 'R': [10.0, 1000.0]},
    'use_complex_bg': True,
})

# With Monte Carlo uncertainty
result = fit_simple("sample.h5",
                    {'model': 'Guinier'},
                    with_uncertainty=True,
                    n_mc_runs=100,
                    q_min=0.005, q_max=0.1)
```

### Supported models

| Model | Key parameters |
|-------|---------------|
| `Guinier` | I0, Rg |
| `Guinier Rod` | I0, Rc |
| `Guinier Sheet` | I0, Rg |
| `Porod` | Kp, Background |
| `Power Law` | P, Exponent, Background |
| `Sphere` | Scale, R |
| `Spheroid` | Scale, R, Beta |
| `Debye-Bueche` | Prefactor, Eta, CorrLength |
| `Teubner-Strey` | Prefactor, A, C1, C2 |
| `Benedetti-Ciccariello` | SolidSLD, VoidSLD, LayerSLD, Sp, t |
| `Hermans` | B, s, d1, d2, sigma1, sigma2 |
| `Hybrid Hermans` | Hermans params + G2, Rg2, G3, Rg3, B3, P3 |
| `Unified Born Green` | G1, Rg1, B1, P1, G2, Rg2, B2, P2, eta, ksi |

Use `from pyirena import MODEL_NAMES` for the full list at runtime.

### Returns

`dict` with keys `success`, `model`, `params`, `params_std`, `I_model`, `q`,
`residuals`, `chi2`, `reduced_chi2`, `dof`, `derived`.  On failure: `success=False`
and `error` message.  Returns `None` only if the data file cannot be loaded.

---

## `fit_waxs`

```python
from pyirena.batch import fit_waxs

result = fit_waxs(
    data_file,          # str or Path — input WAXS data file (HDF5 or text)
    config_file,        # str or Path — pyIrena JSON config with 'waxs_peakfit' section
    save_to_nexus=True, # bool — write results to entry/waxs_peakfit_results in HDF5
)
```

`fit_waxs` is a short alias for `fit_waxs_peaks_from_config`.  It reads the
`waxs_peakfit` section from a `pyirena_config.json` produced by the **WAXS Peak
Fit** GUI panel's **Export Parameters** button, fits all configured peaks plus
background, and (if `save_to_nexus=True`) writes results back to the HDF5 file.

### What it does, step by step

1. **Loads the config file** — reads `'waxs_peakfit'` section; validates the
   `_pyirena_config` header.
2. **Loads data** — same format detection as `fit_unified` (`.h5`/`.hdf5` via
   NXcanSAS, `.dat`/`.txt` via text reader).
3. **Applies Q range** — masks data to `[q_min, q_max]` stored in the config
   (set by the GUI cursors at Export time).
4. **Auto-finds peaks** *(if config has `peak_find` section and no explicit
   peaks)*  — Savitzky-Golay background subtraction + `scipy.signal.find_peaks`.
5. **Fits peaks + background** — `scipy.optimize.curve_fit` with the peak
   shapes and background polynomial specified in the config; bounds are applied
   unless `no_limits=true`.
6. **Saves results** *(if `save_to_nexus=True`)* — writes
   `entry/waxs_peakfit_results` (NXprocess group) to the HDF5 file, including
   per-peak parameter arrays, fit curve, background curve, and residuals.
7. **Returns** a result dict (see below).

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_file` | `str` or `Path` | required | Path to WAXS data file |
| `config_file` | `str` or `Path` | required | Path to pyIrena JSON config |
| `save_to_nexus` | `bool` | `True` | Save fit results to NXcanSAS HDF5 |

### Returns

`dict` with keys:

| Key | Type | Description |
|-----|------|-------------|
| `success` | bool | True if `curve_fit` converged |
| `message` | str | Human-readable status |
| `n_peaks` | int | Number of peaks fitted |
| `bg_shape` | str | Background shape used (e.g. `"Cubic"`) |
| `chi2` | float | χ² (sum of squared weighted residuals) |
| `reduced_chi2` | float | χ² / degrees of freedom |
| `dof` | int | Degrees of freedom |
| `q_min` | float | Lower Q bound used (Å⁻¹) |
| `q_max` | float | Upper Q bound used (Å⁻¹) |
| `q` | ndarray | Q array in fit range (Å⁻¹) |
| `I_fit` | ndarray | Total model intensity over `q` |
| `I_bg` | ndarray | Background-only curve over `q` |
| `residuals` | ndarray | Normalized residuals `(I − fit) / σ` |
| `bg_params` | dict | Fitted background coefficients with uncertainties |
| `peaks` | list | Per-peak dicts: `shape`, `Q0`, `A`, `FWHM`, (`eta`) each with `value` and `std` |

Returns `None` only on fatal pre-fit errors (file not found, data unreadable).

### Peak profile shapes

| Shape | Parameters | Notes |
|-------|-----------|-------|
| `Gauss` | A, Q0, FWHM | Symmetric Gaussian; A = peak height |
| `Lorentz` | A, Q0, FWHM | Symmetric Lorentzian; A = peak height |
| `Pseudo-Voigt` | A, Q0, FWHM, eta | η·Lorentz + (1−η)·Gauss; eta ∈ [0, 1] |
| `LogNormal` | A, Q0, FWHM | Asymmetric; mode at Q0, width from FWHM |

### Background shapes

| Shape | Parameters |
|-------|-----------|
| `Constant` | bg0 |
| `Linear` | bg0, bg1 |
| `Cubic` | bg0, bg1, bg2, bg3 |
| `5th Polynomial` | bg0 … bg5 |

### Example

```python
from pyirena.batch import fit_waxs

# Fit one file
result = fit_waxs("NaHCO3_sample.h5", "pyirena_config.json")
if result and result['success']:
    print(f"Fitted {result['n_peaks']} peak(s),  "
          f"reduced-χ² = {result['reduced_chi2']:.4g}")
    for i, pk in enumerate(result['peaks']):
        q0   = pk['Q0']['value']
        fwhm = pk['FWHM']['value']
        A    = pk['A']['value']
        print(f"  Peak {i+1}: Q0={q0:.4f}  FWHM={fwhm:.4f}  A={A:.4g}")

# Batch over many files
from pathlib import Path
data_files = sorted(Path("data/").glob("*_waxs.h5"))
for f in data_files:
    r = fit_waxs(f, "pyirena_config.json")
    if r and r['success']:
        print(f"{f.name}: {r['n_peaks']} peaks, χ²r={r['reduced_chi2']:.3g}")
```

---

## `fit_modeling`

```python
from pyirena.batch import fit_modeling

result = fit_modeling(
    data_file,              # str or Path — input SAS data file
    config_file,            # str or Path — pyIrena JSON config with 'modeling' section
    save_to_nexus=True,     # bool — write results to NXcanSAS HDF5 file
    with_uncertainty=False, # bool — run MC uncertainty estimation after the main fit
    n_mc_runs=10,           # int — number of MC runs (used when with_uncertainty=True)
)
```

### What it does, step by step

1. **Loads the config file** — validates the `_pyirena_config` header and reads
   the `modeling` section (written by the **Modeling** GUI panel's **Export Parameters**
   button).
2. **Loads data** — same format detection as `fit_unified` (`.h5`/`.hdf5` via
   NXcanSAS, `.dat`/`.txt` via text reader).
3. **Builds populations** — deserializes each population dict from the config into
   the appropriate dataclass based on the `pop_type` field.
4. **Runs the fit** — `ModelingEngine.fit()` using `scipy.optimize.least_squares`
   (TRF with bounds) or Nelder-Mead when `no_limits=True`. If the config's
   `modeling` section sets `"fit_method": "global"`, a `differential_evolution`
   global search runs first to locate the basin, followed by a local TRF polish
   (recommended for multimodal monodisperse core-shell / core-shell-shell fits).
   `fit_method` is ignored — forced to local — when `no_limits` is true, since
   the global search needs finite bounds. The `"de_workers"` key (int, default
   `1`) sets the number of worker processes for the global search; `>1` or `-1`
   (all cores) parallelizes the DE population evaluation, with automatic
   fallback to serial if the host cannot start workers.
5. **MC uncertainty** *(if `with_uncertainty=True`)* — re-fits `n_mc_runs` noise-perturbed
   copies of the data and accumulates per-parameter standard deviations (always
   using the fast local refinement, even when `fit_method` is `global`).
6. **Saves results** *(if `save_to_nexus=True` and input is HDF5)* — writes
   `entry/modeling_results` to the HDF5 file.
7. **Returns** a result dict (see below).

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_file` | `str` or `Path` | required | Path to SAS data file |
| `config_file` | `str` or `Path` | required | Path to pyIrena JSON config containing a `modeling` section |
| `save_to_nexus` | `bool` | `True` | Save fit results to NXcanSAS HDF5; ignored for text-file inputs |
| `with_uncertainty` | `bool` | `False` | Run Monte Carlo uncertainty estimation |
| `n_mc_runs` | `int` | `10` | Number of MC passes for uncertainty estimation |

### Population types

The Modeling tool supports six population types, each stored in the config with
a `pop_type` field:

| `pop_type` | Dataclass | Description |
|------------|-----------|-------------|
| `size_dist` | `SizeDistPopulation` | Parametric size distribution (Gauss, LogNormal, LSW, Schulz-Zimm, Ardell) convolved with a form factor and optional structure factor |
| `unified_level` | `UnifiedLevelPopulation` | Beaucage Unified Fit level: G·exp(−q²Rg²/3) + B·Q*⁻ᴾ + optional Born-Green correlations |
| `diffraction_peak` | `DiffractionPeakPopulation` | Gaussian, Lorentzian, or pseudo-Voigt peak at Q₀ |
| `guinier_porod` | `GuinierPorodPopulation` | Piecewise Guinier-Porod model (Hammouda 2010); 6 shape params G, Rg1, s1, P, Rg2, s2 + optional RgCO and Born-Green correlations |
| `mass_fractal` | `MassFractalPopulation` | Mass fractal aggregate (Teixeira 1988); sphere primary particles; params Phi, Radius, Dv, Ksi, Eta, Contrast |
| `surface_fractal` | `SurfaceFractalPopulation` | Surface fractal (Teixeira 1988); params Surface, Ds, Ksi, Contrast + optional Porod transition at Qc |

Up to 10 populations of any type can be combined in a single fit.

### Returns

```python
{
    'success':     bool,            # True if the optimizer converged
    'message':     str,             # human-readable status line, e.g.
                                    # "Modeling fit complete — χ²/dof=1.23, 2 active population(s)"
    'result':      ModelingResult,  # fully populated result object; see below
    'output_file': Path | None,     # HDF5 file written, or None
}
```

`result` is a `ModelingResult` dataclass with fields:

| Field | Type | Description |
|-------|------|-------------|
| `success` | bool | Whether the optimizer converged |
| `message` | str | Status message from the engine |
| `chi_squared` | float | χ² (weighted sum of squared residuals) |
| `reduced_chi_squared` | float | χ² / degrees of freedom |
| `n_active_pops` | int | Number of enabled populations |
| `pop_indices` | list[int] | 0-based indices of enabled populations |
| `params_std` | dict | Per-parameter standard deviations (populated after MC uncertainty) |
| `q` | ndarray | Q array used for fitting (Å⁻¹) |
| `I_total` | ndarray | Total model intensity over `q` |
| `I_per_pop` | list[ndarray] | Per-population intensity contributions |

Returns `None` only on fatal pre-fit errors (file not found, config unreadable,
data load failure).  A non-converging optimizer returns a dict with `success=False`.

### Example

```python
from pyirena.batch import fit_modeling

# Single file
result = fit_modeling("sample.h5", "pyirena_config.json")
if result and result['success']:
    r = result['result']
    print(f"χ²/dof = {r.reduced_chi_squared:.4g}")
    print(f"Active populations: {r.n_active_pops}")

# With MC uncertainty
result = fit_modeling(
    "sample.h5", "pyirena_config.json",
    with_uncertainty=True, n_mc_runs=50,
)

# Batch over many files
from pathlib import Path
for f in sorted(Path("data/").glob("*.h5")):
    r = fit_modeling(f, "pyirena_config.json")
    status = "OK" if r and r['success'] else "FAIL"
    chi2 = f"{r['result'].reduced_chi_squared:.4g}" if r and r['success'] else "—"
    print(f"{status}  {f.name}  χ²/dof={chi2}")
```

---

## `fit_saxs_morph`

```python
from pyirena.batch import fit_saxs_morph

result = fit_saxs_morph(
    data_file,              # str or Path — input SAS data file
    config_file,            # str or Path — pyIrena JSON config with 'saxs_morph' section
    save_to_nexus=True,     # bool — write voxelgram + metadata to NXcanSAS HDF5 file
    with_uncertainty=False, # bool — reserved for API symmetry; ignored (no fittable parameters)
    n_mc_runs=10,           # int — reserved for API symmetry; ignored
)
```

SAXS Morph reconstructs a three-dimensional two-phase voxelgram whose simulated
I(Q) matches your measured SAXS curve.  The algorithm uses the measured scattering
curve to derive the spectral density function F(k), then draws a random Gaussian
field in reciprocal space, band-limits it at the data's q_max to suppress ringing,
thresholds it at the volume fraction that produces the requested φ, and evaluates
χ² against the experimental data.  No iterative parameter fitting takes place;
the volume fraction and contrast (or scattering length density difference) are set
by the user — this is a direct inversion, not an optimisation.

See the [SAXS Morph GUI documentation](saxs_morph_gui.md) for background theory
and a description of all GUI controls.

### What it does, step by step

1. **Loads the config file** — validates the `_pyirena_config` header and reads
   the `saxs_morph` section (written by the **SAXS Morph** panel's **Export
   Parameters** button).
2. **Loads data** — same format detection as `fit_unified` (`.h5`/`.hdf5` via
   NXcanSAS, `.dat`/`.txt` via text reader).
3. **Applies Q range** — masks data to `[q_min, q_max]` from config.
4. **Pre-fits background** *(if windows are configured)*:
   - **Power-law pre-fit** over `[power_law_q_min, power_law_q_max]` → resolves
     `B` and `P` in I_bg(Q) = B·Q⁻ᴾ.
   - **Flat-background pre-fit** over `[background_q_min, background_q_max]` →
     resolves the flat background level.
5. **Subtracts background** from the data; derives contrast Δρ² from the Porod
   invariant if `input_mode='phi'`, or derives φ if `input_mode='contrast'`.
6. **Computes voxelgram** at fit resolution (`voxel_size_fit`) — builds the
   autocorrelation γ(r), Fourier-transforms to the spectral density F(k),
   band-limits at q_max, generates a Gaussian random field, and thresholds to
   obtain a binary {0, 1} voxelgram.
7. **Re-renders at render resolution** (`voxel_size_render`) using the same RNG
   seed for the final high-resolution voxelgram.
8. **Evaluates χ²** — simulates I(Q) from the voxelgram and compares to data.
9. **Computes topology metrics** — Euler characteristic, genus, connectivity, and
   specific surface area S/V.
10. **Saves results** *(if `save_to_nexus=True` and input is HDF5)* — writes the
    voxelgram array and all scalar results to `entry/saxs_morph_results` in the HDF5 file.
11. **Returns** a result dict (see below).

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_file` | `str` or `Path` | required | Path to SAS data file |
| `config_file` | `str` or `Path` | required | Path to pyIrena JSON config containing a `saxs_morph` section |
| `save_to_nexus` | `bool` | `True` | Save voxelgram and metrics to NXcanSAS HDF5; ignored for text-file inputs |
| `with_uncertainty` | `bool` | `False` | Reserved; ignored (SAXS Morph has no fittable parameters) |
| `n_mc_runs` | `int` | `10` | Reserved; ignored |

### `saxs_morph` config section fields

| Field | Type | Description |
|-------|------|-------------|
| `q_min` | float or null | Lower Q bound of the fitting window (Å⁻¹) |
| `q_max` | float or null | Upper Q bound of the fitting window (Å⁻¹) |
| `power_law_q_min` | float or null | Lower Q bound for power-law background pre-fit (Å⁻¹); null = skip |
| `power_law_q_max` | float or null | Upper Q bound for power-law background pre-fit (Å⁻¹); null = skip |
| `background_q_min` | float or null | Lower Q bound for flat-background pre-fit (Å⁻¹); null = skip |
| `background_q_max` | float or null | Upper Q bound for flat-background pre-fit (Å⁻¹); null = skip |
| `voxel_size_fit` | int | Cube side length (voxels) used during the F(k) → voxelgram step; smaller = faster (e.g. `64` or `128`) |
| `voxel_size_render` | int | Cube side length (voxels) for the final high-resolution voxelgram (e.g. `256` or `512`) |
| `box_size_A` | float | Physical edge length of the voxel cube (Å); sets the real-space scale |
| `input_mode` | str | `"phi"` — φ given, Δρ² derived; `"contrast"` — Δρ² given, φ derived; `"both"` — both given (no invariant derivation) |
| `volume_fraction` | float | Target solid-phase volume fraction φ ∈ (0, 1) |
| `contrast` | float | Scattering contrast Δρ² (cm⁻⁴) |
| `power_law_B` | float | Power-law prefactor B (cm⁻¹·Å⁻ᴾ); 0 = no power-law background |
| `power_law_P` | float | Power-law exponent P; used with `power_law_B` |
| `background` | float | Flat background level (cm⁻¹); 0 = no flat background |
| `smooth_sigma` | float | Gaussian smoothing σ applied after thresholding (voxels); 0 = disabled |
| `rng_seed` | int or null | RNG seed for reproducible voxelgram generation; null = random |

### Returns

`dict` on success, `None` on fatal error (file not found, config missing, data unreadable).

```python
{
    'success':     bool,             # True if voxelgram was computed without error
    'message':     str,              # e.g. "SAXS Morph complete — χ²/dof=1.07, φ=0.312, pitch=3.91 Å"
    'result':      SaxsMorphResult,  # complete result object (see below)
    'output_file': Path | None,      # HDF5 file written, or None
}
```

`result` is a `SaxsMorphResult` dataclass with the following key attributes:

| Attribute | Type | Description |
|-----------|------|-------------|
| `chi_squared` | float | χ² (weighted sum of squared residuals) |
| `reduced_chi_squared` | float | χ² / degrees of freedom |
| `dof` | int | Degrees of freedom |
| `voxelgram` | ndarray | `uint8`, shape `(N, N, N)`, values 0 (void) or 1 (solid) |
| `voxel_size` | int | Cube side length N (voxels) |
| `box_size_A` | float | Physical box size (Å) |
| `voxel_pitch_A` | float | Physical voxel edge length (Å) = `box_size_A / voxel_size` |
| `phi_actual` | float | Realized solid-phase volume fraction after thresholding |
| `rg_A` | float | Radius of gyration estimated from the model I(Q) (Å) |
| `q_max_model_A` | float | Voxel Nyquist limit π / voxel_pitch_A (Å⁻¹) |
| `porod_K_struct` | float | Porod prefactor extracted from the structural I(Q) |
| `specific_surface_area_inv_A` | float | Specific surface area S/V (Å⁻¹) |
| `morphology_metrics` | `MorphologyMetrics` or `None` | Topology metrics: Euler characteristic, genus, connected components |
| `data_q` | ndarray | Q array over the fit window (Å⁻¹) |
| `data_I` | ndarray | Measured intensity over the fit window (cm⁻¹) |
| `data_I_corr` | ndarray | Background-subtracted intensity (cm⁻¹) |
| `model_I` | ndarray | Model intensity over the fit window (cm⁻¹) |
| `gamma_r` | ndarray | Autocorrelation function γ(r) |
| `r_grid` | ndarray | Real-space distance grid for γ(r) (Å) |
| `spectral_F` | ndarray | Spectral density F(k) used to generate the random field |
| `spectral_k` | ndarray | Wavenumber grid for F(k) (Å⁻¹) |
| `config` | `SaxsMorphConfig` | Full configuration used for this run |

### Example

```python
from pyirena.batch import fit_saxs_morph

# Single file — uses all settings from the exported config
result = fit_saxs_morph("sample.h5", "pyirena_config.json")
if result and result['success']:
    r = result['result']
    print(f"χ²/dof          = {r.reduced_chi_squared:.4g}")
    print(f"φ_actual        = {r.phi_actual:.4f}")
    print(f"voxel pitch     = {r.voxel_pitch_A:.2f} Å")
    print(f"Rg              = {r.rg_A:.1f} Å")
    print(f"S/V             = {r.specific_surface_area_inv_A:.4g} Å⁻¹")
    print(f"voxelgram shape = {r.voxelgram.shape}")   # (N, N, N)

# In-memory only — no file written
result = fit_saxs_morph("sample.h5", "pyirena_config.json", save_to_nexus=False)
voxelgram = result['result'].voxelgram   # uint8 ndarray ready for VTK, numpy, etc.

# Access topology metrics
m = result['result'].morphology_metrics
if m is not None:
    print(f"Euler number = {m.euler_number},  genus = {m.genus}")

# Visualise I(Q) fit with matplotlib
import matplotlib.pyplot as plt
r = result['result']
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
ax = axes[0]
ax.loglog(r.data_q, r.data_I, 'k.', ms=3, label='data')
ax.loglog(r.data_q, r.model_I, 'r-', lw=1.5, label='model')
ax.set_xlabel('Q (Å⁻¹)'); ax.set_ylabel('I(Q) (cm⁻¹)'); ax.legend()
ax = axes[1]
ax.imshow(r.voxelgram[r.voxel_size // 2], cmap='gray', origin='lower')
ax.set_title(f"Central slice  φ={r.phi_actual:.3f}")
plt.tight_layout(); plt.show()

# Batch over many files
from pathlib import Path
for f in sorted(Path("data/").glob("*.h5")):
    r = fit_saxs_morph(f, "pyirena_config.json")
    if r and r['success']:
        res = r['result']
        print(f"{f.name}: χ²/dof={res.reduced_chi_squared:.3g}  φ={res.phi_actual:.4f}")
```

### Lower-level engine API

For scripting workflows that bypass the config file entirely:

```python
from pyirena.core.saxs_morph import SaxsMorphConfig, SaxsMorphEngine
import numpy as np

cfg = SaxsMorphConfig(
    q_min=0.001, q_max=0.30,
    power_law_q_min=0.001, power_law_q_max=0.005,
    background_q_min=0.25,  background_q_max=0.30,
    voxel_size_fit=128,
    voxel_size_render=256,
    box_size_A=1000.0,
    input_mode='phi',        # derive contrast from invariant
    volume_fraction=0.30,
    smooth_sigma=1.0,
    rng_seed=42,
)

engine = SaxsMorphEngine()
result = engine.compute_voxelgram(cfg, q, I, dI)   # numpy arrays
print(f"χ²/dof = {result.reduced_chi_squared:.4g}")
print(f"φ_actual = {result.phi_actual:.4f}")
```

---

## `merge_data`

```python
from pyirena.batch import merge_data

result = merge_data(
    file1,                  # str or Path — DS1 file (lower-Q, absolute scale)
    file2,                  # str or Path — DS2 file (higher-Q)
    config_file=None,       # str or Path or None — JSON merge config file
    save_to_nexus=True,     # bool — write merged file to disk
    output_folder=None,     # str or Path or None — where to write the file
    verbose=True,           # bool — print progress to stdout
)
```

Unlike the fitting functions, `merge_data` uses its own **merge config file**
(produced by the GUI's **Save JSON Config** button) rather than the shared
`pyirena_config.json` format.  See the [Data Merge tool documentation](data_merge_gui.md)
for the full config file reference.

### What it does, step by step

1. **Loads data** from both files — detects format by extension:
   - `.h5` / `.hdf5` / `.hdf` → NXcanSAS or generic HDF5
   - `.dat` / `.txt` → column-delimited text (Q, I, dI)
2. **Reads the merge config** *(if provided)* — sets Q overlap range and
   optimisation options (scale, Q shift, split mode, etc.).
   If no config is given, the overlap is auto-detected as the central 80 % of
   the Q intersection and all default options are used.
3. **Optimises** — Nelder-Mead minimisation of weighted χ² in the overlap
   region to find the best scale factor, constant background, and optional
   Q shift.
4. **Merges** — assembles the combined Q/I/dI/dQ arrays, trimming DS1 at the
   right cursor and DS2 at the left cursor (or splitting hard at the left
   cursor when `split_at_left_cursor=True`).
5. **Saves** *(if `save_to_nexus=True`)* — copies DS1 (if NXcanSAS) to the
   output folder, replaces the Q/I/Idev/Qdev arrays with the merged data, and
   appends an `entry/data_merge_results` provenance group.  For non-NXcanSAS
   DS1 files a fresh NXcanSAS file is created.
6. **Returns** a result dict (see below).

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `file1` | `str` or `Path` | required | DS1 file (lower-Q, absolute intensity scale cm⁻¹) |
| `file2` | `str` or `Path` | required | DS2 file (higher-Q, any scale) |
| `config_file` | `str`, `Path`, or `None` | `None` | Merge JSON config (see [Config file format](data_merge_gui.md#config-file-format)).  All keys optional; missing keys use defaults |
| `save_to_nexus` | `bool` | `True` | Write merged data to a NXcanSAS HDF5 file |
| `output_folder` | `str`, `Path`, or `None` | `None` | Output directory.  If `None`, a folder `{DS1_folder}_merged/` is created next to DS1's parent |
| `verbose` | `bool` | `True` | Print scale / BG / χ² results to stdout |

### Returns

`dict` on success, `None` on fatal error (file not found, no Q overlap, etc.).
A non-converged optimisation still returns a dict with `success=False`.

```python
{
    'q':           np.ndarray,        # merged Q array (Å⁻¹), sorted ascending
    'I':           np.ndarray,        # merged intensity (cm⁻¹)
    'dI':          np.ndarray,        # merged uncertainty (cm⁻¹)
    'dQ':          np.ndarray | None, # merged Q-resolution, or None
    'scale':       float,             # optimised (or fixed) scale factor
    'q_shift':     float,             # optimised (or fixed) Q shift (Å⁻¹)
    'background':  float,             # optimised background subtracted from DS1 (cm⁻¹)
    'chi_squared': float,             # weighted χ² in the overlap region
    'success':     bool,              # True if the optimiser converged
    'output_path': str | None,        # full path to the written file, or None
}
```

### Examples

```python
from pyirena.batch import merge_data

# Merge two files — auto-detect overlap, save to default folder
result = merge_data("saxs/sample_001.h5", "waxs/sample_001.h5")
if result and result['success']:
    print(f"scale={result['scale']:.4g}  BG={result['background']:.4g}")
    print(f"Saved to: {result['output_path']}")

# Use a saved config file
result = merge_data(
    "saxs/sample_001.h5",
    "waxs/sample_001.h5",
    config_file="merge_config.json",
    output_folder="merged/",
)

# In-memory only — no file written
result = merge_data(
    "saxs/sample_001.h5",
    "waxs/sample_001.h5",
    config_file="merge_config.json",
    save_to_nexus=False,
)
q, I, dI = result['q'], result['I'], result['dI']

# Batch over matched pairs
from pathlib import Path
saxs_files = sorted(Path("saxs/").glob("*.h5"))
waxs_files = sorted(Path("waxs/").glob("*.h5"))

for f1, f2 in zip(saxs_files, waxs_files):
    r = merge_data(f1, f2, config_file="merge_config.json", output_folder="merged/")
    status = "OK" if r and r['success'] else "FAIL"
    print(f"{status}  {f1.name} + {f2.name}")
```

### Command-line script

For folder-level batch processing from the shell, use the included helper script:

```bash
python scripts/batch_merge.py \
    --folder1 saxs/ \
    --folder2 waxs/ \
    --config merge_config.json \
    --output merged/
```

Run `python scripts/batch_merge.py --help` for full option reference.

---

## `average_data`

```python
from pyirena.batch import average_data

result = average_data(
    data_files,                         # list[str | Path] — files to average
    output_folder=None,                 # str | Path | None — where to save
    verbose=True,                       # bool — print progress to stdout
    similarity_check=False,             # bool — run similarity filter first
    similarity_p_min=0.01,              # float — rejection threshold (0–1)
    similarity_method='cormap',         # str — test algorithm (see below)
    similarity_reference='first',       # str — 'first' or 'majority'
    similarity_normalize_scale=True,    # bool — normalise amplitude before comparing
)
```

`average_data` does **not** use a JSON config file.  All parameters are
passed directly.  The GUI saves its similarity settings to pyIrena's internal
state file (auto-restored on next launch) but those settings are not read
by the batch function — pipeline scripts should hard-code or compute the
values they want, or read them from a user-defined config dict.

### What it does, step by step

1. **Loads** each file in `data_files` — accepts `.h5`/`.hdf5`/`.hdf`
   (NXcanSAS or generic HDF5) and `.dat`/`.txt` (column-delimited text).
   Files that fail to load are skipped with a warning.
2. **Similarity filter** *(if `similarity_check=True`)* — runs the CorMap
   test (or another registered method) on all loaded datasets and marks
   frames with `p_value < similarity_p_min` as rejected.  Rejected frames
   are logged (when `verbose=True`) and excluded before averaging.
   If fewer than 2 frames remain, the function returns `None`.
3. **Averages** — interpolates every dataset onto the first accepted frame's
   Q grid (log-log), then computes:
   ```
   I_avg   = mean(I_i)
   dI_avg  = max( sqrt(sum(dI_i^2)) / N,  std(I_i) )  per Q point
   ```
   Only Q points valid in **all** datasets contribute to the result.
4. **Saves** to a NXcanSAS HDF5 file in `output_folder` (created if absent).
   The file includes an `entry/data_manipulation_results` NXprocess provenance
   group recording the operation, number of datasets, and timestamp.
5. **Returns** a result dict (see below).

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_files` | `list[str \| Path]` | required | Ordered list of files to average.  Order matters only when `similarity_reference='first'`. |
| `output_folder` | `str`, `Path`, or `None` | `None` | Output directory.  If `None`, a folder `{first_file_folder}_manip/` is created next to the first input file. |
| `verbose` | `bool` | `True` | Print per-frame similarity results and final summary to stdout. |
| `similarity_check` | `bool` | `False` | Enable similarity filtering.  When `True`, frames with p < `similarity_p_min` are excluded before averaging. |
| `similarity_p_min` | `float` | `0.01` | P-value rejection threshold.  Frames with p below this are discarded.  Typical range: 0.001 (strict) – 0.05 (loose). |
| `similarity_method` | `str` | `'cormap'` | Test algorithm.  Currently only `'cormap'` (Franke et al., 2015) is registered.  Additional methods can be added to `pyirena/core/similarity.py`. |
| `similarity_reference` | `str` | `'first'` | `'first'` — compare each frame vs. frame #1, which is always kept.  Standard bioSAXS approach when the first exposure is known clean.  `'majority'` — compare each frame vs. the **median** of all frames; more robust when the first frame may also be damaged or when damage is not monotonically increasing. |
| `similarity_normalize_scale` | `bool` | `True` | Recommended: `True`.  Before comparing, rescales each frame to match the reference amplitude using the geometric-median intensity ratio.  This removes flux-drift and absorption differences so that CorMap detects *shape* changes only — the relevant signature of radiation damage.  Set `False` only when data are on a perfectly calibrated absolute scale and you want amplitude changes to count as damage. |

### Returns

`dict` on success, `None` on fatal error (no files loaded, fewer than 2
frames survive similarity filtering, save error, etc.).

```python
{
    'success':    bool,              # True when result was saved
    'operation':  'avg',
    'output_file': Path,            # full path to the written NXcanSAS file
    'message':    str,              # human-readable summary line
    'n_datasets': int,              # number of frames that were averaged
    'rejected':   list[tuple],      # [(filename, p_value), ...] — frames
                                    # discarded by similarity filter; [] when
                                    # similarity_check=False
}
```

### Examples

#### Simple average, no filtering

```python
from pathlib import Path
from pyirena.batch import average_data

files = sorted(Path("frames/").glob("sample_*.h5"))
result = average_data(files, output_folder="averaged/")
if result:
    print(f"Averaged {result['n_datasets']} frames → {result['output_file'].name}")
```

#### Automatic radiation-damage rejection

```python
from pathlib import Path
from pyirena.batch import average_data

files = sorted(Path("frames/").glob("sample_*.h5"))

result = average_data(
    files,
    output_folder="averaged/",
    similarity_check=True,
    similarity_p_min=0.01,
    similarity_reference='first',   # first frame known clean
    similarity_normalize_scale=True,
    verbose=True,                   # prints [OK] / [REJECT] per frame
)

if result:
    print(f"Averaged {result['n_datasets']} of {len(files)} frames")
    for fname, p in result['rejected']:
        print(f"  REJECTED  {fname}  (p={p:.4f})")
```

#### Majority-vote reference (robust when first frame may be damaged)

```python
result = average_data(
    files,
    output_folder="averaged/",
    similarity_check=True,
    similarity_p_min=0.01,
    similarity_reference='majority',  # compare each vs. median of all
    similarity_normalize_scale=True,
)
```

#### Pipeline: process many sample folders

```python
from pathlib import Path
from pyirena.batch import average_data

# Each sub-folder contains sequential exposures for one sample
data_root = Path("raw/")
out_root  = Path("averaged/")

# Similarity settings determined empirically on test data
SIM_PARAMS = dict(
    similarity_check=True,
    similarity_p_min=0.01,
    similarity_reference='first',
    similarity_normalize_scale=True,
)

for sample_dir in sorted(data_root.iterdir()):
    if not sample_dir.is_dir():
        continue
    frames = sorted(sample_dir.glob("*.h5"))
    if len(frames) < 2:
        continue

    out_dir = out_root / sample_dir.name
    result = average_data(frames, output_folder=out_dir, **SIM_PARAMS)

    if result is None:
        print(f"FAILED  {sample_dir.name}")
    else:
        n_kept = result['n_datasets']
        n_total = len(frames)
        n_rej = len(result['rejected'])
        print(f"OK  {sample_dir.name}: {n_kept}/{n_total} frames "
              f"({n_rej} rejected)")
```

#### In-memory usage (no file written, just get the averaged arrays)

```python
from pyirena.batch import _load_data
from pyirena.core.data_manipulation import DataManipulation
from pyirena.core.similarity import check_similarity
from pathlib import Path

files = sorted(Path("frames/").glob("*.h5"))

# Load data
datasets, names = [], []
for fp in files:
    d = _load_data(fp)
    if d is None:
        continue
    datasets.append((d['Q'], d['Intensity'], d.get('Error', d['Intensity']*0.05), d.get('dQ')))
    names.append(fp.name)

# Check similarity
results = check_similarity(datasets, filenames=names,
                           reference='first', p_min=0.01)
accepted = [datasets[r.idx] for r in results if r.accepted]

# Average accepted frames
avg = DataManipulation.average(accepted, reference_index=0)
q, I, dI = avg.q, avg.I, avg.dI
```

### Choosing a p-value threshold

The p-value is the probability that two truly identical curves would
produce a run of same-sign residuals at least as long as the observed one
by random chance alone.  Low p → curves differ → likely damaged.

| Threshold | Effect |
|-----------|--------|
| 0.001 | Strict: only reject severe outliers (large, obvious damage) |
| 0.010 | Moderate: good default starting point |
| 0.050 | Loose: rejects borderline frames; cleaner average but fewer frames used |

**Recommended calibration workflow:**

1. Run with `verbose=True` on a clean sample (known undamaged exposures).
   Note the lowest p-value observed.  Set `similarity_p_min` just below that.
2. Run on a sample known to contain damaged frames.  Raise the threshold
   until the damaged frames appear in `result['rejected']`.
3. Fix the threshold in your pipeline configuration.

---

## `igor_to_nexus`  *(formerly `pxp_to_nexus`)*

```python
from pyirena.batch import igor_to_nexus

result = igor_to_nexus(
    igor_file,                   # str or Path — input .pxp OR .h5xp
    output_folder=None,          # str/Path/None — default: <stem>_data
    techniques=None,             # list[str]/None — e.g. ["USAXS", "SAXS"]
    overwrite=False,             # bool — append _2,_3,… if False
    verbose=False,               # bool — print per-file summary
)
```

Import an Igor Pro packed experiment and export each reduced USAXS /
SAXS / WAXS sample as a stand-alone NXcanSAS HDF5 file. Accepts both
`.pxp` (legacy binary) and `.h5xp` (modern HDF5) — the format is
auto-detected from the file extension. This is the headless counterpart
of the GUI's **Import Igor Experiment…** button.

> `pxp_to_nexus(...)` remains as a deprecated alias for back-compat.

### What it does

1. Walks the experiment's folder tree (`root:USAXS:...`, `root:SAXS:...`,
   `root:WAXS:...`).
2. For each sample folder, picks the documented wave triple
   (`DSM_*` for USAXS, `R_*` for SAXS/WAXS).
3. Parses the wave note — including APS USAXS's `NXSampleStart`/`End`
   and `NXInstrumentStart`/`End` sentinel blocks — into NXcanSAS
   sample/instrument metadata.
4. Writes one `.h5` per sample into
   `<output_folder>/<technique>/<sample>.h5`.

### Returns

```python
{
    'success': True,
    'output_folder': str,
    'n_written': int,
    'n_skipped': int,
    'n_errors':  int,
    'n_unparseable_records': int,
    'files': [
        {
            'source':    str,   # Igor folder path, e.g. "root:SAXS:Sample01"
            'output':    str,   # absolute path to the .h5 file (or None)
            'technique': str,   # "USAXS" | "SAXS" | "WAXS"
            'n_points':  int,
            'status':    str,   # "ok" | "skipped" | "error"
            'message':   str,   # explanation for skipped/error rows
        },
        ...
    ],
}
```

Returns `None` only if the input `.pxp` cannot be found.

### Example

```python
from pyirena.batch import igor_to_nexus, fit_pyirena
from pathlib import Path

# 1. Bring legacy Igor data into NeXus.
imported = igor_to_nexus("2020_beamtime.pxp", techniques=["USAXS"])
print(f"Got {imported['n_written']} USAXS files")

# 2. Fit each one with the standard pyirena config.
out_dir = Path(imported['output_folder']) / "USAXS"
for h5 in sorted(out_dir.glob("*.h5")):
    fit_pyirena(str(h5), "pyirena_config.json")
```

See the [Import Igor Experiment guide](igor_pxp_import.md) for the
wave-name and folder-name conventions, and how to extend them for
non-standard reduction pipelines.

---

## `fit_pyirena`

```python
from pyirena.batch import fit_pyirena

results = fit_pyirena(
    data_file,              # str or Path — input SAS data file
    config_file,            # str or Path — pyIrena JSON config file
    save_to_nexus=True,     # bool — passed to each tool's fitting function
    with_uncertainty=False,
    n_mc_runs=10,
    tools=None,             # list[str] or None — limit which tools run (see below)
)
```

### What it does

Reads the config file, discovers which tool sections are present, and dispatches
to the appropriate fitting function for each.  Recognised tool sections:

| Config key | Function called |
|------------|----------------|
| `unified_fit` | `fit_unified()` |
| `sizes` | `fit_sizes()` |
| `simple_fits` | `fit_simple_from_config()` |
| `waxs_peakfit` | `fit_waxs()` |
| `modeling` | `fit_modeling()` |
| `saxs_morph` | `fit_saxs_morph()` |

Unknown sections in the config file are silently skipped.  This means a config
file created today will still work correctly when new tools are added to pyIrena,
and vice-versa.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_file` | `str` or `Path` | required | Input SAS data file |
| `config_file` | `str` or `Path` | required | pyIrena JSON config file |
| `save_to_nexus` | `bool` | `True` | Passed to each tool's fitting function |
| `with_uncertainty` | `bool` | `False` | Run MC uncertainty estimation |
| `n_mc_runs` | `int` | `10` | Number of MC passes (when `with_uncertainty=True`) |
| `tools` | `list[str]` or `None` | `None` | If given, only the listed tool names are run, even if other sections exist in the config. When `None` (default), every recognised section is executed — identical to behaviour before this parameter was added, so existing callers are unaffected. |

**Tip — limiting tools at the scripting level:**

```python
# Run only unified_fit and sizes, skip simple_fits / modeling even if present
results = fit_pyirena("sample.h5", "pyirena_config.json",
                      tools=["unified_fit", "sizes"])
```

**Tip — removing stale sections from the config file:**

If you have accumulated tool sections in `pyirena_config.json` that you no longer
want to run, use the **Manage Config...** button in the Data Selector utility row.
It opens a dialog that lists all tool sections as checkboxes; uncheck a section and
click Save to remove it from the file permanently.

### Returns

`dict` or `None`.  Returns `None` only if the config file itself cannot be loaded.

```python
{
    'input_file':  Path,               # data file that was processed
    'config_file': Path,               # config file used
    'tools_run':   ['unified_fit'],    # list of tool names that were executed
    'results': {
        'unified_fit': { ... },        # fit_unified() return dict, or None if it failed
        # future tools will appear here automatically
    }
}
```

---

## Return structures

### `fit_unified` return dict

```python
{
    # --- Status ---
    'success':     bool,    # True if the optimizer converged
    'tool':        str,     # always 'unified_fit'
    'message':     str,     # human-readable one-liner, e.g.
                            # "Unified Fit: 2 level(s), chi²=1.23, success=True"

    # --- File paths ---
    'input_file':  Path,    # data file that was fitted
    'config_file': Path,    # config file used
    'output_file': Path,    # NXcanSAS output file, or None if save_to_nexus=False
                            # or if the save failed (non-fatal)

    # --- Model object ---
    'model': UnifiedFitModel,   # fully configured & fitted model object;
                                # call model.calculate_intensity(q) for any Q array

    # --- Raw fit result from model.fit() ---
    'fit_result': {
        'success':             bool,
        'message':             str,
        'chi_squared':         float,
        'reduced_chi_squared': float,
        'n_iterations':        int,
        'levels':              list,    # list of UnifiedLevel objects (fitted values)
        'background':          float,
        'fit_intensity':       np.ndarray,
        'residuals':           np.ndarray,
    },

    # --- Structured parameter summary (most useful for downstream code) ---
    'parameters': {
        'num_levels':          int,
        'background':          float,   # fitted background (cm⁻¹)
        'chi_squared':         float,
        'reduced_chi_squared': float,
        'n_iterations':        int,
        'levels': [
            {
                'level':      int,      # 1-based
                'G':          float,    # Guinier prefactor (cm⁻¹)
                'Rg':         float,    # radius of gyration (Å)
                'B':          float,    # Porod prefactor
                'P':          float,    # power-law slope
                'RgCutoff':   float,    # cutoff radius (Å); 0 for level 1
                'ETA':        float,    # correlation distance (Å)
                'PACK':       float,    # packing factor
                'correlated': bool,     # Born-Green correlation active
            },
            # ... one entry per fitted level
        ],
    },

    # --- Data arrays ---
    'data': {
        'Q':               np.ndarray,  # full Q array from file (Å⁻¹)
        'Intensity':       np.ndarray,  # full experimental intensity (cm⁻¹)
        'Error':           np.ndarray | None,  # experimental uncertainty
        'Q_fit':           np.ndarray,  # Q sub-range actually used for fitting
                                        # (after cursor masking; may equal Q)
        'intensity_model': np.ndarray,  # model evaluated on full Q
        'residuals':       np.ndarray,  # (I_data - I_model) / Error (or / I_data)
                                        # evaluated on full Q
    },
}
```

### Accessing results

```python
result = fit_unified("sample.h5", "pyirena_config.json")

# Guard against total failure
if result is None:
    print("Fatal error — see messages above")
    sys.exit(1)

# Guard against non-convergence
if not result['success']:
    print(f"Fit did not converge: {result['message']}")

# Scalar parameters
p = result['parameters']
print(p['chi_squared'])
print(p['background'])

# Per-level parameters
for lv in p['levels']:
    print(f"Level {lv['level']}: Rg={lv['Rg']:.1f} Å, G={lv['G']:.3g}")

# Plot with matplotlib
import matplotlib.pyplot as plt
d = result['data']
plt.loglog(d['Q'], d['Intensity'], 'k.', label='data')
plt.loglog(d['Q'], d['intensity_model'], 'r-', label='fit')
plt.legend(); plt.show()

# Further model evaluation at arbitrary Q
import numpy as np
q_dense = np.logspace(-3, 0, 500)
i_dense = result['model'].calculate_intensity(q_dense)
```

---

## Error handling

The batch functions are designed to be safe inside automated pipelines.
No exception will propagate out of `fit_unified` or `fit_pyirena`.

| Situation | Return value | Printed message |
|-----------|-------------|-----------------|
| Config file not found / unreadable | `None` | Yes |
| File is not a pyIrena config | `None` | Yes |
| `unified_fit` section missing from config | `None` | Yes |
| Data file not found | `None` | Yes |
| Data file unreadable | `None` | Yes |
| Model construction error | `None` | Yes (with traceback) |
| Optimizer did not converge | `dict` with `success=False` | Yes |
| NXcanSAS save failed | `dict` with `output_file=None` | Yes (non-fatal) |
| Any other unexpected exception | `None` | Yes (with traceback) |

### Recommended pattern for automated processing

```python
from pathlib import Path
from pyirena.batch import fit_unified

config = "pyirena_config.json"
data_files = sorted(Path("data/").glob("*.h5"))

results = []
failed = []

for f in data_files:
    r = fit_unified(f, config, save_to_nexus=True)
    if r is None:
        failed.append({'file': f.name, 'reason': 'fatal error'})
    elif not r['success']:
        failed.append({'file': f.name, 'reason': r['message']})
    else:
        results.append(r)

print(f"Fitted {len(results)}/{len(data_files)} files successfully.")
if failed:
    print("Failed files:")
    for item in failed:
        print(f"  {item['file']}: {item['reason']}")
```

---

## Batch processing patterns

### Collect results into a table

```python
import pandas as pd
from pyirena.batch import fit_unified

config = "pyirena_config.json"
rows = []

for f in data_files:
    r = fit_unified(f, config)
    if not r or not r['success']:
        continue
    p = r['parameters']
    row = {
        'file':    r['input_file'].name,
        'chi2':    p['chi_squared'],
        'rchi2':   p['reduced_chi_squared'],
        'bg':      p['background'],
    }
    for lv in p['levels']:
        n = lv['level']
        row[f'Rg_{n}']  = lv['Rg']
        row[f'G_{n}']   = lv['G']
        row[f'B_{n}']   = lv['B']
        row[f'P_{n}']   = lv['P']
    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv("fit_results.csv", index=False)
```

### Use without saving to NXcanSAS (library mode)

```python
# Call from another package that manages its own data storage
result = fit_unified(data_file, config_file, save_to_nexus=False)

if result and result['success']:
    model      = result['model']
    parameters = result['parameters']
    q_model    = result['data']['Q']
    i_model    = result['data']['intensity_model']
    # hand off to your own code ...
```

### Parallel batch fitting

```python
from concurrent.futures import ProcessPoolExecutor, as_completed
from pyirena.batch import fit_unified

def _fit_one(args):
    return fit_unified(*args)

jobs = [(f, "pyirena_config.json") for f in data_files]
results = []

with ProcessPoolExecutor() as pool:
    futures = {pool.submit(_fit_one, job): job[0] for job in jobs}
    for future in as_completed(futures):
        r = future.result()
        if r and r['success']:
            results.append(r)
```

---

## Extending with new tools

When a new analysis tool is added to pyIrena, three things happen:

1. A fitting function `fit_<toolname>(data_file, config_file, save_to_nexus)` is
   added, following the same conventions as `fit_unified`.

2. The GUI's **Export Parameters** writes the tool's config into a new top-level
   key in the JSON file (e.g. `"size_distribution": { ... }`).

3. The tool is registered in `fit_pyirena`'s `_TOOL_REGISTRY` dict in
   [pyirena/batch.py](../pyirena/batch.py):

```python
_TOOL_REGISTRY = {
    'unified_fit':  lambda: fit_unified(data_file, config_file, save_to_nexus, ...),
    'sizes':        lambda: fit_sizes(data_file, config_file, save_to_nexus, ...),
    'simple_fits':  lambda: fit_simple_from_config(data_file, config_file, save_to_nexus, ...),
    'waxs_peakfit': lambda: fit_waxs(data_file, config_file, save_to_nexus, ...),
    'modeling':     lambda: fit_modeling(data_file, config_file, save_to_nexus, ...),
    'saxs_morph':   lambda: fit_saxs_morph(data_file, config_file, save_to_nexus, ...),
    # add new tools here
}
```

That is all.  Existing config files that do not contain the new key are unaffected.
Config files that do contain it will automatically have the new tool run when
`fit_pyirena` is called.
