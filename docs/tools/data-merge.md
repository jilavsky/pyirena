(data-merge)=
# Data Merge Tool

:::{admonition} At a glance
:class: tip

**Launch:** `pyirena-gui` → **Tools → Data Merge**, or the **Data Merge** button in the main
hub. Standalone: `pyirena-datamerge` (no arguments); headless:
`pyirena-datamerge --file1 DS1.h5 --file2 DS2.h5 [options]`. \
**Python:** `from pyirena.batch import merge_data`; `pyirena.api.merge_datasets`. \
**MCP:** `merge_datasets` and `match_merge_files` through `pyirena_call` (dispatcher category
`"data"`); `pyirena_read_merge_provenance` reads the provenance of a merged file. \
**Igor equivalent:** {doc}`Irena Data manipulation tools <igor:Irena/DataManipulation>`
:::

The **Data Merge** tool combines two SAS datasets measured on different instruments or Q ranges
(e.g. USAXS + SAXS, SAXS + WAXS/diffraction) into a single consistent dataset. It finds the
best-fit scale factor, constant background, and optional Q-shift by optimising a weighted
chi-squared in the user-defined overlap region.

(data-merge-window)=
## The window

```{figure} img/data-merge/data-merge-01-window.png
:alt: Data Merge window after Optimize Merge on merge_usaxs.h5 and merge_saxs.h5: the DS1 and DS2 file columns on the left, the merge controls and overlap Q range above the log-log plot with cursors A and B in the middle, and the output-folder and save buttons on the right. Numbered circles mark the parts described in the legend.

The Data Merge window after **Optimize Merge** on `validationData/merge/merge_usaxs.h5` (DS1) and
`merge_saxs.h5` (DS2).
```

1. **Dataset 1** column (lower Q, absolute scale): **Select Folder…**, **Type**, **Sort**,
   **Filter**, **Strip** and the file list.
2. **Dataset 2** column (higher Q), the same controls.
3. **Mode** (display) and **DS2 on right axis**.
4. **Scale**, **Q shift** and **BG (DS1)**, each with its **Fit** box and value field.
5. **Split at left cursor**, **Optimize Merge** and the method menu.
6. **Overlap Q range**: Q min / Q max of cursors A and B.
7. **Plot**: DS1, DS2, the merged curve and cursors A and B.
8. **Output folder**: **Create Default Folder**, **Select Existing Folder…**.
9. **Save Merged Data**, **Batch Run**, **Save JSON Config…**, **Load JSON Config…**.
10. **Match files**.
11. **Status bar**: χ², overlap point count and convergence.

(data-merge-workflow)=
## Typical workflow

The steps use `validationData/merge/merge_usaxs.h5` (DS1, the low-Q branch of a two-level Unified
curve on the correct scale) and `merge_saxs.h5` (DS2, the high-Q branch of the same curve,
deliberately mis-scaled by 1.15). With **Scale: DS2** the merge should recover
1/1.15 = 0.870.

1. **Load datasets.**
   - Click **Select Folder…** in the DS1 column and choose the folder containing your lower-Q
     files. Repeat for DS2.
   - Select a **Type** (HDF5 NXcanSAS, HDF5 Generic, or Text .dat/.txt/.csv). Text files are
     assumed to be in the Q unit set in the Data Selector's *Configure → Text File Options → Q
     unit in text files* (default 1/Å) — see
     [Q units](../guides/data-import-and-cleaning.md#q-units).
   - Use the **Sort** dropdown and the **Filter** field to order and narrow each file list (see
     [Controls](#data-merge-file-columns)).
   - **Double-click** a file to load it and display it on the plot.

   ```{figure} img/data-merge/data-merge-02-loaded.png
   :alt: Log-log plot with the USAXS branch (blue, low Q) and the SAXS branch (red, high Q) overlapping between about 0.005 and 0.05 inverse angstrom, cursors A at 0.01 and B at 0.04.

   Both datasets loaded, cursors A and B at Q = 0.01 and 0.04 Å⁻¹.
   ```

2. **Set the overlap region.** Drag cursors **A** and **B** on the plot to bracket the Q region
   where both datasets overlap. The current Q min / Q max values are shown in the **Overlap Q
   range** box below the controls.
3. **Configure optimisation.** Choose which dataset is scaled and which parameters are fitted;
   see [Controls](#data-merge-controls). For the example: **Scale: DS2**, Fit on, **BG (DS1)**
   Fit on (default).
4. **Run optimisation.** Click **Optimize Merge**. The result fields update immediately:
   - **Scale** — fitted (or fixed) scale factor
   - **Q shift** — fitted (or fixed) Q shift (Å⁻¹)
   - **BG (DS1)** — fitted constant background (cm⁻¹)
   - Status bar shows χ², overlap point count, and convergence status

   The merged curve is overlaid on the plot in green. On the example files the scale comes out at
   0.872 (expected 0.870) and BG (DS1) at 0.0016 cm⁻¹.
5. **Adjust by hand** if needed. Optimize does not always converge to something sensible — poorly
   subtracted data at low Q, or beam-stop scattering left in the high-Q segment, can leave it with
   nothing good to find. When that happens, set the three numbers yourself and watch the merged
   curve redraw; see [Adjusting by hand](#data-merge-manual).
6. **Save.**
   1. Click **Create Default Folder** (creates `<DS1_folder>_merged/` next to DS1) or **Select
      Existing Folder…** to choose an output location.
   2. Click **Save Merged Data**. The output file is a NXcanSAS HDF5 copy of DS1 with the
      Q/I/Idev/Qdev arrays replaced by the merged data and a `data_merge_results` provenance group
      appended.

   Many file pairs: see [Batch run](#data-merge-batch).

% SCREENSHOT data-merge-03-porod-real-data: Porod (log I·Q⁴ vs log Q) display mode on a real USAXS + SAXS
% pair from testData/ (Jan picks one) where a residual background or beam-stop scattering shows as
% a step or slope in the overlap; both curves on the same flat line after a good merge.

(data-merge-controls)=
## Controls

(data-merge-file-columns)=
### Dataset columns

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Select Folder…** | Choose the folder containing the files for this dataset | | |
| **Type** | HDF5 NXcanSAS, HDF5 Generic, or Text .dat/.txt/.csv | HDF5 Nexus | |
| **Sort** | Order the file list — filename, temperature (`_500C`), time (`_10min`), pressure (`_100PSI`) or order number (trailing `_03`), ascending or descending | Order number ↑ | Each column remembers its own choice between sessions |
| **Filter** | Narrow the file list with a **regular expression** matched anywhere in the name, case-insensitively | blank | |
| **Strip** | Regex removed from the start of each filename before matching | blank | See [Batch run](#data-merge-batch) |
| File list | **Double-click** a file to load it and display it on the plot | | Ctrl/Shift-click selects pairs for **Batch Run** |
| **Use slit-smeared copy** | Merge the slit-smeared copy instead of the desmeared one | off | Appears only when the file has both; see [Slit smearing](#data-merge-slit-smearing) |

This panel previously sorted by order number with no way to change it; the default is unchanged,
and each column remembers its own choice between sessions. Files whose name lacks the pattern sort
last, so they never split an ordered series. Filter examples: `60C|100C`, or `^(?!.*bkg)` to
exclude backgrounds (grep semantics). An invalid pattern falls back to a plain substring match.

### Merge controls

| Control | Effect | Default | Tips |
|---------|--------|---|---|
| **Mode** | How the data is drawn — see [Display modes](#data-merge-display-modes) | SAXS (log-log) | The merge always runs on the original Q and I |
| **DS2 on right axis** | Give DS2 its own Y scale — see [Display modes](#data-merge-display-modes) | off | To find a feature, not to judge the match |
| **Scale: DS1 / DS2** | Which dataset the scale factor is applied to | DS1 | |
| **Fit** (Scale) | Whether **Optimize** is free to vary the scale | on | |
| **Q shift: None / DS1 / DS2** | Which dataset an additive Q shift is applied to | None | |
| **Fit** (Q shift) | Whether **Optimize** is free to vary the Q shift | off | |
| **BG (DS1)** | Constant background subtracted from DS1 | 0 | |
| **Fit** (BG) | Whether **Optimize** is free to vary the background (default: on) | on | |
| **Split at left cursor** | Hard split mode (see [Concepts](#data-merge-concepts)) | off | |
| **Optimize Merge** | Optimise the free parameters in the overlap region; results are previewed on the graph | | Use **Save Merged Data** to write the file |
| **Method** | Interpolation method (currently log-log linear interpolation) | Log-log linear interpolation | |

% CHECK(Jan): the panel's Scale menu starts at DS1, while the config default `scale_dataset` is 2
% (DS2) and Concepts says DS1 is the absolute-scale reference.

The three value fields are **always editable**, whether or not their **Fit** box is ticked. A
ticked **Fit** box means only that **Optimize** may change that number — so running Optimize again
will overwrite a value you set by hand, which is what asking to optimise it means.

At least one of **Scale** / **Background** must stay checked — a Q shift alone has nothing to
optimize against; the panel warns and refuses to run otherwise.

(data-merge-manual)=
### Adjusting by hand

| Gesture | Step |
|---------|------|
| Type a value and press Enter | exactly what you typed |
| Mouse wheel over the field | 1% of the current value |
| **↑** / **↓** with the field focused | 1% |
| **▲** / **▼** buttons beside the field | 1% |
| **Shift** + any of the above | 10% — coarse |
| **Ctrl/Cmd** + any of the above | 0.1% — fine |

Steps are a percentage of the current value, so the same gesture is useful for a scale near 1 and
a background near 1e-5. A background sitting at exactly zero is the exception — a percentage of
zero is zero — so it steps by a small absolute amount derived from the data instead.

The preview redraws about every 150 ms while you scroll, and it is computed the same way the saved
file is: **what you see is what Save writes**, trimmed according to **Split at left cursor**
exactly as an optimised merge would be. The status bar shows `Manual —` while you are tuning, so it
is always clear whether the numbers on screen came from Optimize or from you.

Switch to the **Porod** display mode while doing this. It is far easier to see a few-percent
mismatch as a step in a flat line than as a slight thickening of a log-log curve.

(data-merge-display-modes)=
### Display modes

The **Mode** menu changes only what is drawn. **The merge always runs on the original Q and I**, in
every mode, and the overlap cursors keep their physical Q positions when you switch.

| Mode | Plot | Use it when |
|------|------|-------------|
| **SAXS (log-log)** | I vs Q, both axes log | The default; the whole Q range at once |
| **WAXS / diffraction (lin-lin)** | I vs Q, both axes linear | Peaks, where log scaling flattens the thing you care about |
| **Porod (log I·Q⁴ vs log Q)** | I·Q⁴ vs Q, both axes log | Judging whether the two curves really agree |

**Porod** is the one to reach for when a merge looks fine on a log-log plot and the result is
still wrong. Multiplying by Q⁴ flattens a Q⁻⁴ tail and, on a log y axis, stretches the intensity
scale over just the range the overlap occupies. A residual background or leftover beam-stop
scattering then shows up as a slope or a step in the overlap region — a difference of a few percent
that a plain log I vs log Q plot hides. Judge the overlap by whether the two curves lie on the
*same* line, not by whether each looks plausible on its own.

The axes are log-log, as in the Porod tab of Unified Fit. The horizontal axis is plain Q (not Q⁴,
which crushes everything below the overlap against the left edge and was tried and rejected).

**Changing an axis from the right-click menu.** *Plot options → Transforms* can flip either axis
between linear and log in any mode. The overlap cursors follow: they stay on the same Q, and the
Q range handed to Optimize and Save is the one you see in the Q min / Q max boxes.

**Lost or invisible cursors.** If a cursor is off the graph or Q min / Q max show nonsense, press
**Reset cursors** (next to the Q boxes). They return to the middle 80% of the Q range the two
datasets share. A saved overlap range that does not make sense for the loaded data is ignored in
the same way.

Note that **Simple Fits plots its Porod linearization differently** — I·Q⁴ against *Q⁴* — and that
is correct there, not an inconsistency. For `I = Kp·Q⁻⁴ + Background`, I·Q⁴ = Kp + Background·Q⁴
is a straight line against Q⁴ whose intercept is Kp and whose slope is the background, which is
what you want when fitting those two numbers. Here you are not fitting anything; you are comparing
two curves, so readable Q spacing wins.

**DS2 on right axis** gives DS2 its own Y scale, so two datasets decades apart in intensity both
fill the plot. Use it to *find* a feature. Do not use it to judge the match: with two scales the
curves no longer overlap in any physical sense, and they will look aligned whatever the scale
factor is. Unticked (the default) the right-hand axis simply mirrors the left, so you can read a
value near the right edge without tracking back across the plot.

(data-merge-results)=
## Results and outputs

| Control | What it does |
|---|---|
| **Create Default Folder** | Create `<DS1 folder>_merged` next to the DS1 folder and use it as the output folder |
| **Select Existing Folder…** | Choose an existing folder for saving merged files |
| **Save Merged Data** | Write the merged dataset to an NXcanSAS HDF5 file in the output folder |
| **Batch Run** | Merge all matched (or selected) pairs; see [Batch run](#data-merge-batch) |
| **Save JSON Config…** | Saves current Q range and optimisation settings to a JSON file. Load this file in a later session with **Load JSON Config…** or pass it to the CLI / Python API |

% CHECK(Jan): the Save JSON Config… tooltip says it also saves the folders; the page says "Q range
% and optimisation settings".

Merged data is saved as a NXcanSAS HDF5 file.

**When DS1 is NXcanSAS:** the DS1 file is copied to the output folder with `_merged` appended to
the stem. The Q/I/Idev/Qdev arrays in the first `SASdata` group are replaced with the merged
arrays. All other metadata (sample name, instrument, etc.) is preserved.

**When DS1 is text or generic HDF5:** a new NXcanSAS file `<DS1_stem>_merged.h5` is created from
scratch.

**Provenance group.** A group `entry/data_merge_results` (class `NXprocess`) is added to every
output file, recording the full merge provenance (see also {ref}`hdf5-results-groups`):

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

The merged file is ordinary reduced data: open it in any fitting tool. An AI agent reads the
provenance with `pyirena_read_merge_provenance(path)`.

(data-merge-scripting)=
## From Python or AIDA

::::{tab-set}
:::{tab-item} Python
Run from the repository root (outputs go to `merged/` in the current folder):

```python
import json

from pyirena.batch import merge_data

# Merge two files with default settings (auto-detect overlap)
result = merge_data("validationData/merge/merge_usaxs.h5",
                    "validationData/merge/merge_saxs.h5",
                    output_folder="merged/")
if result and result['success']:
    print(f"scale={result['scale']:.4g}  BG={result['background']:.4g}")   # scale ≈ 1/1.15
    print(f"Saved to: {result['output_path']}")

# Use a saved config file (written here by hand; the GUI's Save JSON Config… writes one)
with open("merge_config.json", "w") as f:
    json.dump({"version": "1.0", "q_overlap_min": 0.01, "q_overlap_max": 0.04,
               "fit_scale": True, "scale_dataset": 2, "fit_background": True}, f, indent=2)
result = merge_data(
    "validationData/merge/merge_usaxs.h5",
    "validationData/merge/merge_saxs.h5",
    config_file="merge_config.json",
    output_folder="merged/",
)

# In-memory only (no file written)
result = merge_data(
    "validationData/merge/merge_usaxs.h5",
    "validationData/merge/merge_saxs.h5",
    config_file="merge_config.json",
    save_to_nexus=False,
)
import numpy as np
q, I, dI = result['q'], result['I'], result['dI']

# Batch over matched file pairs
from pathlib import Path

usaxs_files = sorted(Path("validationData/merge/").glob("merge_usaxs*.h5"))
saxs_files = sorted(Path("validationData/merge/").glob("merge_saxs*.h5"))
config = "merge_config.json"

for f1, f2 in zip(usaxs_files, saxs_files):
    r = merge_data(f1, f2, config_file=config, output_folder="merged/")
    status = "OK" if r and r['success'] else "FAIL"
    print(f"{status}  {f1.name} + {f2.name}")
```

The full signature, parameters and return value are under
[Python API reference](#data-merge-python-api).
:::
:::{tab-item} AIDA / MCP
Merging is available headlessly through `pyirena.api` and over MCP as the dispatcher's `data`
category, using the same engine and the same output naming as the panel.

```python
from pyirena import api

# pair two folders, then merge each pair
# (with real folders: pairs = api.match_merge_files("/data/usaxs", "/data/saxs")["pairs"])
pairs = [{"file1": "validationData/merge/merge_usaxs.h5",
          "file2": "validationData/merge/merge_saxs.h5"}]
for pair in pairs:
    print(api.merge_datasets(pair["file1"], pair["file2"], output_folder="merged")["scale"])
```

`file1` must be the lower-Q dataset; passing them the other way round returns a
`SWAPPED_INPUTS` error naming the correct order rather than silently rescaling the wrong curve.
Omitting `q_overlap_min` / `q_overlap_max` auto-detects the overlap the same way batch mode does.

The returned dict carries the fitted `scale`, `background`, `chi_squared`, `n_overlap_points` and
any `slit_warning`. As in the panel, a background is fitted by default (`fit_background=True`) —
pass `fit_background=False` for data with no background assumption; at least one of `fit_scale` /
`fit_background` must stay True — see [Optimisation details](#data-merge-optimisation).

See [ai_tools_reference.md](../reference/ai-tools-reference.md) for the agent-facing tool
descriptions.
:::
::::

(data-merge-method)=
## Method notes

(data-merge-concepts)=
### DS1 and DS2

- **DS1** — the lower-Q dataset, assumed to be on an **absolute intensity scale** (cm⁻¹). DS1 is
  never rescaled when Scale dataset = DS2.
- **DS2** — the higher-Q dataset (arbitrary or calibrated scale).

### Overlap region

The two cursors **A** (left) and **B** (right) on the plot define the Q range used for
optimisation. In the output merged file:

| Checkbox state | DS1 contribution | DS2 contribution |
|----------------|-----------------|-----------------|
| **Split at left cursor** ON | Q < cursor A | Q ≥ cursor A |
| **Split at left cursor** OFF (default) | Q ≤ cursor B | Q ≥ cursor A |

When *Split at left cursor* is OFF the overlap region [A, B] contains interleaved points from both
datasets.

% CHECK(Jan): the Split at left cursor tooltip says "Default (unchecked): include all data from
% both datasets", not DS1 only up to cursor B and DS2 only from cursor A.

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
| Background | **Fit** checkbox next to BG (DS1); checked by default | Type value in the box when Fit is unchecked |

(data-merge-optimisation)=
### Optimisation details

**Algorithm.** Nelder-Mead simplex (scipy) minimising a weighted chi-squared in the overlap region.
DS1 is log-log-linearly interpolated onto DS2's Q grid.

% CHECK(Jan): with Scale = DS2 and no Q shift the engine reports "Analytical WLS" (closed form;
% the api docstring says the same), not Nelder-Mead.

**Free parameters.**

| Slot | Free when | Bound |
|------|----------|-------|
| Background (BG) | `fit_background=True` (default) | ±max(I1 in overlap) |
| Scale | `fit_scale=True` | [0.01, 100] |
| Q shift | `fit_qshift=True` and `qshift_dataset≠0` | [−0.1, 0.1] Å⁻¹ |

At least one of `fit_scale` / `fit_background` must be True — otherwise there is nothing to
optimize against.

**Degeneracy handling.** For `scale_dataset=2` the objective has a degenerate valley where
`BG = I1 − I2·scale`, making large-BG / small-scale solutions equally valid from a chi-squared
standpoint. This is resolved by a background regularisation term `(BG / I1_median)²` added to the
objective, which prefers smaller backgrounds when data chi-squared is otherwise equal.

**Initial guess.** The initial scale is estimated as the median-intensity ratio of the two datasets
in the overlap region (`median(I1) / median(I2)` for `scale_dataset=2`). Background starts at 0
when fitted, or at the fixed value when not.

(data-merge-slit-smearing)=
### Slit smearing (USAXS)

When an input file carries **both** a desmeared and a slit-smeared (USAXS) copy (Matilda), a **"Use
slit-smeared copy"** checkbox appears under that dataset's file list. It loads the desmeared
`@default` entry by default; check it to merge the slit-smeared data instead.

Merging a slit-smeared curve (typically the low-Q USAXS dataset) with a pinhole one produces a
**slit-smeared output**: the merged file gets a `dQl` dataset so downstream tools auto-detect it,
and the provenance records the input and merged slit lengths (`slit_length_ds1/ds2`,
`slit_length_merged`). Two inputs with *different* nonzero slit lengths warn (the larger is kept).
The optimisation itself is unchanged — the slit length sits at or below the SAXS Qmin, so its
effect in the overlap region is negligible. See **[Slit smearing](../methods/slit-smearing.md)**
for the full reference.

(data-merge-batch)=
## Batch run

Enable **Match files** to automatically pair DS1 and DS2 files that share the same prefix (text
before the first `_`) and the same trailing number. Dataset **Filter** fields are applied first —
matching only considers files that pass their respective filter. Changing a filter while Match mode
is active re-runs the matching automatically. Then click **Batch Run** to process all pairs
sequentially.

If the two datasets glue an instrument code onto the sample name instead of keeping the base name
identical — e.g. `SmySample_0001.dat` (SAXS) vs. `WmySample_0001.dat` (WAXS) — use the
per-dataset **Strip** field to remove that instrument code before matching. It takes a regex,
applied once to the start of each filename, e.g. `^S` for DS1 and `^W` for DS2; both then key on
`mySample_0001` and pair correctly. Leave both **Strip** fields blank for the default behavior.
Like **Filter**, changing **Strip** while Match mode is active re-runs the matching automatically.

Alternatively, select specific files in both lists using Ctrl+click or Shift+click. The batch run
then uses the selected pairs in list order (first selected DS1 with first selected DS2, etc.).

(data-merge-config)=
## Config file format

The JSON config file stores the overlap Q range and all optimisation settings. It is written by
**Save JSON Config…** and read by **Load JSON Config…**, `pyirena-datamerge --config`, and
`merge_data(config_file=...)`.

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
  "fit_background": true,
  "fixed_background_value": 0.0,
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
| `fit_background` | bool | `true` | Optimise constant background subtracted from DS1 |
| `fixed_background_value` | float | `0.0` | Background (cm⁻¹) used when `fit_background` is `false` |
| `split_at_left_cursor` | bool | `false` | Hard split at left cursor |

At least one of `fit_scale` / `fit_background` must be `true` — a Q shift alone has nothing to
optimize against.

If `q_overlap_min` / `q_overlap_max` are absent, `merge_data()` auto-detects the overlap as the
central 80 % of the Q intersection.

(data-merge-cli)=
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

When both `--file1` and `--file2` are supplied, the tool runs without a GUI and exits with code 0
on success or 1 on failure.

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

(data-merge-python-api)=
## Python API reference

```python
from pyirena.batch import merge_data
```

### Signature

```text
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
| `file1` | `str` or `Path` | required | DS1 file. Accepts `.h5`/`.hdf5`/`.hdf` (NXcanSAS or generic HDF5) or `.dat`/`.txt` (text) |
| `file2` | `str` or `Path` | required | DS2 file. Same formats as `file1` |
| `config_file` | `str`, `Path`, or `None` | `None` | JSON config file (see [Config file format](#data-merge-config)). All keys are optional; missing keys use defaults |
| `save_to_nexus` | `bool` | `True` | Write merged data to a NXcanSAS HDF5 file |
| `output_folder` | `str`, `Path`, or `None` | `None` | Output directory. If `None`, a folder named `{DS1_folder}_merged` is created next to DS1's parent directory |
| `verbose` | `bool` | `True` | Print progress and results to stdout |

### Returns

`dict` on success, `None` on fatal error.

```text
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

Returns `None` if either file cannot be loaded or if there is no Q overlap. A failed optimisation
(not converged) still returns a dict with `success=False`.

### Using the engine directly (advanced)

For fine-grained control you can call the engine and merge functions directly:

```python
import numpy as np
from pyirena.core.data_merge import DataMerge, MergeConfig
from pyirena.io import discover_scattering, load_scattering

ds1 = load_scattering(discover_scattering("validationData/merge/merge_usaxs.h5")[0])
ds2 = load_scattering(discover_scattering("validationData/merge/merge_saxs.h5")[0])
q1, I1, dI1, dQ1 = ds1.q, ds1.intensity, ds1.uncertainty, ds1.dq
q2, I2, dI2, dQ2 = ds2.q, ds2.intensity, ds2.uncertainty, ds2.dq

engine = DataMerge()

config = MergeConfig(
    q_overlap_min=0.01,
    q_overlap_max=0.04,
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
