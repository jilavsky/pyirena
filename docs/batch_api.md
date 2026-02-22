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
6. [API reference — `fit_pyirena`](#fit_pyirena)
7. [Return structures](#return-structures)
8. [Error handling](#error-handling)
9. [Batch processing patterns](#batch-processing-patterns)
10. [Extending with new tools](#extending-with-new-tools)

---

## Quick start

```python
from pyirena.batch import fit_unified, fit_sizes, fit_simple, fit_pyirena

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
Saves the full size distribution, residuals, and all scalar parameters to
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
| `Treubner-Strey` | Prefactor, A, C1, C2 |
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

## `fit_pyirena`

```python
from pyirena.batch import fit_pyirena

results = fit_pyirena(
    data_file,          # str or Path — input SAS data file
    config_file,        # str or Path — pyIrena JSON config file
    save_to_nexus=True, # bool — passed to each tool's fitting function
    with_uncertainty=False,
    n_mc_runs=10,
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

Unknown sections in the config file are silently skipped.  This means a config
file created today will still work correctly when new tools are added to pyIrena,
and vice-versa.

### Parameters

Same as `fit_unified`.

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
    'unified_fit': lambda: fit_unified(data_file, config_file, save_to_nexus, ...),
    'sizes':       lambda: fit_sizes(data_file, config_file, save_to_nexus, ...),
    'simple_fits': lambda: fit_simple_from_config(data_file, config_file, save_to_nexus, ...),
    # add new tools here
}
```

That is all.  Existing config files that do not contain the new key are unaffected.
Config files that do contain it will automatically have the new tool run when
`fit_pyirena` is called.
