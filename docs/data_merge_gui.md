# Data Merge Tool

The **Data Merge** tool combines two SAS datasets measured on different
instruments or Q ranges (e.g. USAXS + SAXS, SAXS + WAXS/diffraction) into a
single consistent dataset.  It finds the best-fit scale factor, constant
background, and optional Q-shift by optimising a weighted chi-squared in the
user-defined overlap region.

The tool is available three ways:

| Mode | How to launch |
|------|--------------|
| Interactive GUI | `pyirena-datamerge` (no arguments), or **Data Merge** button in the main hub |
| Headless CLI | `pyirena-datamerge --file1 DS1.h5 --file2 DS2.h5 [options]` |
| Python API | `from pyirena.batch import merge_data` |

---

## Contents

1. [Concepts](#concepts)
2. [GUI walkthrough](#gui-walkthrough)
3. [Config file format](#config-file-format)
4. [CLI reference](#cli-reference)
5. [Python API reference — `merge_data`](#python-api-reference)
6. [Output file format](#output-file-format)
7. [Optimisation details](#optimisation-details)

---

## Concepts

### DS1 and DS2

- **DS1** — the lower-Q dataset, assumed to be on an **absolute intensity
  scale** (cm⁻¹).  DS1 is never rescaled when Scale dataset = DS2.
- **DS2** — the higher-Q dataset (arbitrary or calibrated scale).

### Overlap region

The two cursors **A** (left) and **B** (right) on the plot define the Q range
used for optimisation.  In the output merged file:

| Checkbox state | DS1 contribution | DS2 contribution |
|----------------|-----------------|-----------------|
| **Split at left cursor** ON | Q < cursor A | Q ≥ cursor A |
| **Split at left cursor** OFF (default) | Q ≤ cursor B | Q ≥ cursor A |

When *Split at left cursor* is OFF the overlap region [A, B] contains
interleaved points from both datasets.

### What is optimised?

The optimiser minimises a weighted chi-squared in the overlap region:

```
χ² = Σ  [(I1_interp(q) − BG − I2(q)·scale)² / σ²(q)]
```

where:

- `I1_interp` — DS1 intensity log-log-interpolated onto DS2's Q grid
- `BG` — constant background **subtracted from DS1** before comparison
- `scale` — multiplicative factor applied to the selected dataset
- `σ` — per-point uncertainty (at least 5 % of the local intensity)

Three parameters may be free or fixed:

| Parameter | Fit checkbox | Fixed value |
|-----------|-------------|-------------|
| Scale | **Fit** checkbox next to Scale | Type value in the box when Fit is unchecked |
| Q shift (additive) | **Fit** checkbox next to Q shift | Type value in the box when Fit is unchecked |
| Background | Always optimised | — |

---

## GUI walkthrough

### 1 — Load datasets

- Click **Select Folder…** in the DS1 column and choose the folder containing
  your lower-Q files.  Repeat for DS2.
- Select a **Type** (HDF5 NXcanSAS, HDF5 Generic, or Text .dat/.txt).
- Use the **Filter** field to narrow the file list.
- **Double-click** a file to load it and display it on the plot.

### 2 — Set the overlap region

Drag cursors **A** and **B** on the plot to bracket the Q region where both
datasets overlap.  The current Q min / Q max values are shown in the
**Overlap Q range** box below the controls.

### 3 — Configure optimisation

| Control | Effect |
|---------|--------|
| **Mode** | SAXS (log-log plot) or WAXS / diffraction (lin-lin) |
| **Scale: DS1 / DS2** | Which dataset the scale factor is applied to |
| **Fit** (Scale) | Optimise scale; uncheck to fix it at the typed value |
| **Q shift: None / DS1 / DS2** | Which dataset an additive Q shift is applied to |
| **Fit** (Q shift) | Optimise Q shift; uncheck to fix it at the typed value |
| **BG (DS1)** | Constant background (always optimised, result shown read-only) |
| **Split at left cursor** | Hard split mode (see [Concepts](#concepts)) |
| **Method** | Interpolation method (currently log-log linear interpolation) |

### 4 — Run optimisation

Click **Optimize Merge**.  The result fields update immediately:

- **Scale** — fitted (or fixed) scale factor
- **Q shift** — fitted (or fixed) Q shift (Å⁻¹)
- **BG (DS1)** — fitted constant background (cm⁻¹)
- Status bar shows χ², overlap point count, and convergence status

The merged curve is overlaid on the plot in green.

### 5 — Save

1. Click **Create Default Folder** (creates `<DS1_folder>_merged/` next to DS1)
   or **Select Existing Folder…** to choose an output location.
2. Click **Save Merged Data**.  The output file is a NXcanSAS HDF5 copy of DS1
   with the Q/I/Idev/Qdev arrays replaced by the merged data and a
   `data_merge_results` provenance group appended.

### 6 — Batch run

Enable **Match files** to automatically pair DS1 and DS2 files that share the
same prefix (text before the first `_`) and the same trailing number.
Dataset **Filter** fields are applied first — matching only considers files
that pass their respective filter.  Changing a filter while Match mode is
active re-runs the matching automatically.
Then click **Batch Run** to process all pairs sequentially.

Alternatively, select specific files in both lists using Ctrl+click or
Shift+click.  The batch run then uses the selected pairs in list order
(first selected DS1 with first selected DS2, etc.).

### Saving and loading config

- **Save JSON Config…** — saves current Q range and optimisation settings to a
  JSON file.  Load this file in a later session with **Load JSON Config…** or
  pass it to the CLI / Python API.

---

## Config file format

The JSON config file stores the overlap Q range and all optimisation settings.
It is written by **Save JSON Config…** and read by **Load JSON Config…**,
`pyirena-datamerge --config`, and `merge_data(config_file=...)`.

### Annotated example

```json
{
  "version": "1.0",
  "q_overlap_min": 0.08,
  "q_overlap_max": 0.25,
  "fit_scale": true,
  "scale_dataset": 2,
  "fixed_scale_value": 1.0,
  "fit_qshift": false,
  "fixed_qshift_value": 0.0,
  "qshift_dataset": 0,
  "split_at_left_cursor": false
}
```

### Field reference

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `q_overlap_min` | float | auto-detect | Left cursor position (Å⁻¹) |
| `q_overlap_max` | float | auto-detect | Right cursor position (Å⁻¹) |
| `fit_scale` | bool | `true` | Optimise scale factor |
| `scale_dataset` | int | `2` | Dataset to scale: `1` = DS1, `2` = DS2 |
| `fixed_scale_value` | float | `1.0` | Scale used when `fit_scale` is `false` |
| `fit_qshift` | bool | `false` | Optimise additive Q shift |
| `fixed_qshift_value` | float | `0.0` | Q shift (Å⁻¹) used when `fit_qshift` is `false` |
| `qshift_dataset` | int | `0` | Dataset to shift: `0` = none, `1` = DS1, `2` = DS2 |
| `split_at_left_cursor` | bool | `false` | Hard split at left cursor |

If `q_overlap_min` / `q_overlap_max` are absent, `merge_data()` auto-detects
the overlap as the central 80 % of the Q intersection.

---

## CLI reference

```
pyirena-datamerge [options]
```

### GUI mode (no `--file1` / `--file2`)

```bash
# Open the GUI with no pre-filled folders
pyirena-datamerge

# Open the GUI with folders pre-filled
pyirena-datamerge --folder1 /data/saxs --folder2 /data/waxs

# Open GUI with folders and match mode enabled
pyirena-datamerge --folder1 /data/saxs --folder2 /data/waxs --match
```

### Headless mode (`--file1` + `--file2`)

When both `--file1` and `--file2` are supplied, the tool runs without a GUI
and exits with code 0 on success or 1 on failure.

```bash
# Merge two files, auto-detect overlap, save to default folder
pyirena-datamerge --file1 saxs/sample_001.h5 --file2 waxs/sample_001.h5

# Use a saved config file for Q range and settings
pyirena-datamerge --file1 saxs/sample_001.h5 --file2 waxs/sample_001.h5 \
                  --config merge_config.json

# Save to a specific output folder
pyirena-datamerge --file1 saxs/sample_001.h5 --file2 waxs/sample_001.h5 \
                  --config merge_config.json --output merged/
```

### All options

| Option | Description |
|--------|-------------|
| `--folder1 DIR` | Pre-populate DS1 folder in the GUI |
| `--folder2 DIR` | Pre-populate DS2 folder in the GUI |
| `--file1 FILE` | DS1 file — triggers headless mode |
| `--file2 FILE` | DS2 file — triggers headless mode |
| `--config JSON` | JSON config file (Q range + optimisation settings) |
| `--output DIR` | Output folder for headless mode |
| `--match` | Enable file-matching mode when opening the GUI |

---

## Python API reference

```python
from pyirena.batch import merge_data
```

### Signature

```python
result = merge_data(
    file1,                  # str or Path — DS1 file (lower Q, absolute scale)
    file2,                  # str or Path — DS2 file (higher Q)
    config_file=None,       # str or Path or None — JSON config file
    save_to_nexus=True,     # bool — write merged file to disk
    output_folder=None,     # str or Path or None — where to write the file
    verbose=True,           # bool — print progress to stdout
)
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `file1` | `str` or `Path` | required | DS1 file.  Accepts `.h5`/`.hdf5`/`.hdf` (NXcanSAS or generic HDF5) or `.dat`/`.txt` (text) |
| `file2` | `str` or `Path` | required | DS2 file.  Same formats as `file1` |
| `config_file` | `str`, `Path`, or `None` | `None` | JSON config file (see [Config file format](#config-file-format)).  All keys are optional; missing keys use defaults |
| `save_to_nexus` | `bool` | `True` | Write merged data to a NXcanSAS HDF5 file |
| `output_folder` | `str`, `Path`, or `None` | `None` | Output directory.  If `None`, a folder named `{DS1_folder}_merged` is created next to DS1's parent directory |
| `verbose` | `bool` | `True` | Print progress and results to stdout |

### Returns

`dict` on success, `None` on fatal error.

```python
{
    'q':           np.ndarray,   # merged Q array (Å⁻¹), sorted ascending
    'I':           np.ndarray,   # merged intensity (cm⁻¹)
    'dI':          np.ndarray,   # merged uncertainty (cm⁻¹)
    'dQ':          np.ndarray or None,  # merged Q-resolution, or None
    'scale':       float,        # optimised (or fixed) scale factor
    'q_shift':     float,        # optimised (or fixed) Q shift (Å⁻¹)
    'background':  float,        # optimised background subtracted from DS1 (cm⁻¹)
    'chi_squared': float,        # weighted χ² in the overlap region
    'success':     bool,         # True if the optimizer converged
    'output_path': str or None,  # full path to the written file, or None
}
```

Returns `None` if either file cannot be loaded or if there is no Q overlap.
A failed optimisation (not converged) still returns a dict with `success=False`.

### Quick-start examples

```python
from pyirena.batch import merge_data

# Merge two files with default settings (auto-detect overlap)
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

# In-memory only (no file written)
result = merge_data(
    "saxs/sample_001.h5",
    "waxs/sample_001.h5",
    config_file="merge_config.json",
    save_to_nexus=False,
)
import numpy as np
q, I, dI = result['q'], result['I'], result['dI']

# Batch over matched file pairs
from pathlib import Path

saxs_files = sorted(Path("saxs/").glob("*.h5"))
waxs_files = sorted(Path("waxs/").glob("*.h5"))
config = "merge_config.json"

for f1, f2 in zip(saxs_files, waxs_files):
    r = merge_data(f1, f2, config_file=config, output_folder="merged/")
    status = "OK" if r and r['success'] else "FAIL"
    print(f"{status}  {f1.name} + {f2.name}")
```

### Using the engine directly (advanced)

For fine-grained control you can call the engine and merge functions directly:

```python
import numpy as np
from pyirena.core.data_merge import DataMerge, MergeConfig

engine = DataMerge()

config = MergeConfig(
    q_overlap_min=0.08,
    q_overlap_max=0.25,
    fit_scale=True,
    scale_dataset=2,     # scale DS2
    fit_qshift=False,
    fixed_qshift_value=0.0,
    split_at_left_cursor=False,
)

result = engine.optimize(q1, I1, dI1, q2, I2, dI2, config)
print(f"scale={result.scale:.4g}  BG={result.background:.4g}  χ²={result.chi_squared:.4g}")

q_m, I_m, dI_m, dQ_m = engine.merge(
    q1, I1, dI1, dQ1,
    q2, I2, dI2, dQ2,
    result, config,
)
```

---

## Output file format

Merged data is saved as a NXcanSAS HDF5 file.

### When DS1 is NXcanSAS

The DS1 file is copied to the output folder with `_merged` appended to the
stem.  The Q/I/Idev/Qdev arrays in the first `SASdata` group are replaced with
the merged arrays.  All other metadata (sample name, instrument, etc.) is
preserved.

### When DS1 is text or generic HDF5

A new NXcanSAS file `<DS1_stem>_merged.h5` is created from scratch.

### Provenance group

A group `entry/data_merge_results` (class `NXprocess`) is added to every
output file, recording the full merge provenance:

| Dataset / attribute | Content |
|--------------------|---------|
| `ds1_file` | Full path to DS1 input file |
| `ds2_file` | Full path to DS2 input file |
| `scale` | Scale factor applied |
| `q_shift` | Q shift applied (Å⁻¹) |
| `background` | Background subtracted from DS1 (cm⁻¹) |
| `chi_squared` | χ² in the overlap region |
| `n_overlap_points` | Number of DS2 points in the overlap |
| `q_overlap_min` | Left cursor position (Å⁻¹) |
| `q_overlap_max` | Right cursor position (Å⁻¹) |
| `scale_dataset` | 1 or 2 |
| `fit_scale` | 1 (true) or 0 (false) |
| `qshift_dataset` | 0, 1, or 2 |
| `fit_qshift` | 1 (true) or 0 (false) |
| `split_at_left_cursor` | 1 (true) or 0 (false) |
| `timestamp` (attr) | ISO-8601 timestamp of the merge |

---

## Optimisation details

### Algorithm

Nelder-Mead simplex (scipy) minimising a weighted chi-squared in the overlap
region.  DS1 is log-log-linearly interpolated onto DS2's Q grid.

### Free parameters

| Slot | Free when | Bound |
|------|----------|-------|
| Background (BG) | Always | ±max(I1 in overlap) |
| Scale | `fit_scale=True` | [0.01, 100] |
| Q shift | `fit_qshift=True` and `qshift_dataset≠0` | [−0.1, 0.1] Å⁻¹ |

### Degeneracy handling

For `scale_dataset=2` the objective has a degenerate valley where
`BG = I1 − I2·scale`, making large-BG / small-scale solutions equally valid
from a chi-squared standpoint.  This is resolved by a background regularisation
term `(BG / I1_median)²` added to the objective, which prefers smaller
backgrounds when data chi-squared is otherwise equal.

### Initial guess

The initial scale is estimated as the median-intensity ratio of the two
datasets in the overlap region (`median(I1) / median(I2)` for `scale_dataset=2`).
Background always starts at 0.
