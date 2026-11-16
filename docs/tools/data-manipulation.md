(data-manipulation)=
# Data Manipulation Tool

:::{admonition} At a glance
:class: tip

**Launch:** `pyirena-gui` → **Tools > Data Manipulation** menu, or the **Data Manipulation**
button in the main hub. Standalone: `pyirena-datamanip` (no arguments). \
**Python:** `from pyirena.batch import manipulate_data, average_data`; `pyirena.api.scale_data`,
`average_data`, `subtract_data`, … \
**MCP:** `average_data`, `subtract_data`, `divide_data`, `scale_data`, `trim_data`, `rebin_data`
through `pyirena_call` (dispatcher category `"data"`); `pyirena_read_manipulation_provenance`
reads the provenance of an output file. \
**Igor equivalent:** {doc}`Irena Data manipulation tools <igor:Irena/DataManipulation>`
:::

The **Data Manipulation** tool provides common post-reduction operations on SAS datasets: scaling,
trimming, rebinning, averaging, buffer subtraction, and division. It handles datasets on different
Q grids by interpolating in log-log space and propagates uncertainties through every operation.

(data-manipulation-window)=
## The window

```{figure} img/data-manipulation/data-manipulation-01-window.png
:alt: Data Manipulation window on the Scale tab with manipulate_input.h5 loaded, Scale I 0.4 and Background 0.12: the file browser on the left, the operation tabs and the log-log plot with the input (blue) and the scaled result (green) in the middle, and the output buttons on the right. Numbered circles mark the parts described in the legend.

The Data Manipulation window on the **Scale** tab, `validationData/manipulation/manipulate_input.h5`
scaled by 0.4 with a background of 0.12 subtracted.
```

1. **Select Folder…**, **Type**, **Sort** and **Filter**: the file browser.
2. **File list**: double-click to plot; select several for Average and batch runs; right-click
   for **Set as buffer** / **Set as denominator**.
3. **Operation tabs**: Scale, Trim, Rebin, Average, Subtract, Divide.
4. **Tab-specific controls** (here Scale I, Background, Scale uncert.).
5. **Plot**: input data and the result (green).
6. **Output folder**: **Create Default Folder**, **Select Existing Folder…**.
7. **Apply & Save** and **Batch All Selected**.
8. **Status bar**.

(data-manipulation-workflow)=
## Typical workflow

The steps use `validationData/manipulation/manipulate_input.h5`, which is
`manipulate_reference.h5` multiplied by 2.5 with 0.3 cm⁻¹ added. Undoing that with **Scale I** = 0.4
and **Background** = 0.12 must give back the reference exactly.

1. Open the tool from the pyIrena hub or run `pyirena-datamanip`.
2. Click **Select Folder** and choose a folder containing your data files.
3. Select the appropriate file type (HDF5 Nexus, HDF5 Generic, or Text).
4. Double-click a file to load it into the graph.
5. Choose an operation tab (Scale, Trim, Rebin, Average, Subtract, Divide).
6. Adjust parameters --- results update automatically. For the example: **Scale** tab, Scale I =
   0.4, Background = 0.12; the green result is the reference curve (it agrees with
   `manipulate_reference.h5` to about 1e-13).

   ```{figure} img/data-manipulation/data-manipulation-02-scale-result.png
   :alt: Log-log plot of manipulate_input.h5 (blue, levelling off at 0.3 at high Q) and the scaled result (green), which is a factor 2.5 lower and keeps falling at high Q because the 0.3 offset is removed.

   Scale tab: input (blue) and result (green).
   ```

   % CHECK(Jan): the plot's y range follows the input data, so the green result is cut off below
   % I ≈ 0.25 cm⁻¹ at Q > 0.04 Å⁻¹.

7. Set an output folder, then click **Apply & Save**.

The same file pair checks **Subtract**: set `manipulate_reference.h5` as the buffer, double-click
`manipulate_input.h5` and set **Buffer scale** to 2.5. What is left is the 0.3 cm⁻¹ offset.

```{figure} img/data-manipulation/data-manipulation-03-subtract.png
:alt: Subtract tab with sample manipulate_input.h5 (blue) and buffer manipulate_reference.h5 (red), buffer scale 2.5; the green result is a horizontal line at 0.3 inverse centimetres.

Subtract tab: sample (blue) − 2.5 × buffer (red) = a flat 0.3 cm⁻¹ (green).
```

% SCREENSHOT data-manipulation-04-similarity: Average tab on a real series of sequential exposures
% from testData/ (Jan picks one with radiation damage) after Check Similarity: coloured datasets
% on the graph and the results table with green accepted and red rejected rows.

(data-manipulation-controls)=
## Controls

### File browser and graph

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Select Folder** | Choose the folder containing your data files | | |
| **Type** | File format: HDF5 Nexus (NXcanSAS), HDF5 Generic, or Text (.dat/.txt/.csv) | HDF5 Nexus | Text files are assumed to be in the Q unit set in the Data Selector's *Configure → Text File Options → Q unit in text files* (default 1/Å) — see [Q units](../guides/data-import-and-cleaning.md#q-units) |
| **Sort** | 10 sort options: filename, temperature, time, order number, pressure (ascending or descending) | Order number ↑ | |
| **Filter** | Type text to filter filenames | blank | Regular expression, see below |
| File list | **Double-click** a file to load and plot it | | Right-click: **Set as buffer**, **Set as denominator** |

The Filter text is a **regular expression** matched anywhere in the name, case-insensitively (grep
semantics): `60C` , `60C|100C` , `0[12]min` , `^sample` , `\.h5$` , `^(?!.*bkg)` . An invalid
pattern falls back to a plain substring match.

Graph (centre):

- Log-log I(Q) plot with the top axis showing R = pi/Q.
- Input data displayed as scatter points; result overlaid as a thick green line.
- Right-click the plot for the standard export menu (clipboard copy, image, curve CSV, Igor ITX —
  see [Data Selector](data-selector.md)) and, on the Average tab, dataset removal.

### Operations

All operations update automatically when you change parameters. There is no Preview button ---
results recalculate within ~500 ms of any change.

(data-manipulation-scale)=
#### Scale + Background

Multiply intensity by a constant and/or subtract a flat background.

```
I_out = scale * I_in - background
dI_out = scale_uncertainty * dI_in
```

| Control | Meaning | Default | Tips |
|---------|------------|---|---|
| **Scale I** | Multiplicative factor (mouse-wheel adjustable, step 0.01) | 1.0 | |
| **Background** | Constant subtracted after scaling (mouse-wheel, step 0.001) | 0.0 | |
| **Scale uncert.** | Factor for uncertainty; leave empty to use Scale I | empty (= Scale I) | |

**Output suffix:** `_scaled`

(data-manipulation-trim)=
#### Trim Q range

Remove data points outside a Q range defined by two draggable cursors.

- Drag the red and blue vertical cursor lines on the plot.
- Q min and Q max update in real time.
- Only data between the cursors is kept.

**Output suffix:** `_trimmed`

(data-manipulation-rebin)=
#### Rebin

Rebin data onto a new Q grid.

| Control | Meaning | Default | Tips |
|---------|------------|---|---|
| **Grid** | Log-spaced, Linear, or From reference file | Log-spaced | |
| **Points** | Number of points in the new grid (10--10000, default 200) | 200 | |
| **Q min / Q max** | Range for the new grid (auto-filled from data) | from data | |
| **Load Reference** | Select "From reference file" first, then click to load a data file whose Q values become the new grid | | |

Rebinning uses log-log interpolation with relative-uncertainty propagation.

**Output suffix:** `_rebinned`

(data-manipulation-average)=
#### Average

Average multiple datasets to improve statistics. Common in any experiment where sequential
exposures are collected and radiation-damaged frames must be identified and discarded before
averaging.

**Workflow:**

1. Select multiple files in the file list (Ctrl+click or Shift+click).
2. All selected datasets are plotted with distinct colours and averaged automatically.
3. Optionally run the **Similarity analysis** to flag damaged frames automatically (see below).
4. **Right-click on the graph** for manual removal:
   - **Remove dataset** --- submenu with colour-coded entries matching the plot.
   - **Remove all after** --- remove a dataset and everything measured after it.
5. The average recomputes automatically after each removal.

The average is computed on the first selected file's Q grid. Other datasets are log-log
interpolated onto that grid before averaging.

**Uncertainty:** the larger of propagated error and sample spread:

```
dI_avg = max( sqrt(sum(dI_i^2)) / N,  std(I_values) )
```

This is robust: propagated error dominates for consistent data, while the spread correctly
reflects systematic drift (e.g. radiation damage).

**Output suffix:** `_avg`

(data-manipulation-similarity)=
##### Similarity analysis — automatic radiation damage detection

The **Similarity analysis** group inside the Average tab uses the **CorMap test** (Franke et al.,
*Nature Methods* 12, 419–422, 2015) to compute a p-value for each frame. A low p-value means the
frame differs from the reference beyond what random noise alone can explain — the hallmark of
radiation damage or a transient contamination event (e.g. a bubble in the beam path).

| Control | Meaning | Default | Tips |
|---------|-------------|---|---|
| **Method** | Test algorithm. Currently *CorMap (Franke 2015)*. Additional methods can be registered in `pyirena/core/similarity.py`. | CorMap (Franke 2015) | |
| **Reference** | *First frame* — each dataset is compared to dataset #1, which is always accepted. Standard bioSAXS approach. Use when the first exposure is known to be clean. *Majority vote* — each dataset is compared to the **median** of all datasets; more robust when the first frame may itself be damaged or when damage is not monotonically increasing. | First frame | |
| **P-value threshold** | Frames with p-value below this are flagged as *Rejected* (default 0.01). Lower threshold = stricter (fewer rejections); higher = looser (more rejections but cleaner average). | 0.01 | |
| **Normalize scale** | Recommended: **on**. Rescales each frame to match the reference amplitude before comparing, removing flux-drift and absorption differences so that CorMap detects *shape* changes only. Turn off only when data are on an identical absolute scale and you want scale changes to count as damage. | on | |
| **Check Similarity** | Runs the test and shows a colour-coded results table | | |
| **Auto-reject N below threshold** | Deselects all rejected frames and re-computes the average | | |
| **Save CSV…** | Writes every row of the results table | | |

% CHECK(Jan): "Lower threshold = stricter (fewer rejections); higher = looser (more rejections)" —
% with p < threshold rejected, a lower threshold rejects fewer frames, which is usually called
% looser; the threshold guide below ("0.001 (strict)") uses the same wording.

**Check Similarity** runs the test and shows a colour-coded results table: green rows are
accepted, red rows are rejected. The status bar reports `N results for M files` — if N < M, some
files could not be loaded.

**Auto-reject N below threshold** deselects all rejected frames in one click and re-computes the
average. Manual right-click removal on the graph still works independently.

**Getting the numbers out.** The results table is no longer view-only:

- **Ctrl+C** copies the selected cells as tab-separated text (Ctrl+Shift+C includes the column
  headers); with nothing selected the whole table is copied. Right-click for the same options.
- **Save CSV…** writes every row — filename, p-value, longest run, number of points,
  accepted/rejected — at full precision.
- Click a column header to sort; sorting the p-value column numerically brings the worst frames to
  the top. The reference frame ("— (ref)") sorts last.

> **How to choose a p-value threshold**
> Run the tool on a clean test dataset (frames you know are undamaged) and note the lowest p-value
> you see. Set the threshold just below that value so good frames are accepted. Then test on a
> dataset known to contain damaged frames and raise the threshold until they are rejected.
> Typical starting values: 0.001 (strict) or 0.010 (moderate).

(data-manipulation-subtract)=
#### Subtract (buffer)

Subtract buffer/solvent scattering from a sample. Common in bioSAXS.

```
I_out = I_sample - buffer_scale * I_buffer
dI_out = sqrt( dI_sample^2 + (buffer_scale * dI_buffer)^2 )
```

**Workflow:**

1. Right-click a file in the list and choose **Set as buffer**.
2. Double-click a sample file. Both sample (blue) and buffer (red) are plotted.
3. Adjust **Buffer scale** with the mouse wheel (step 0.001) or type a value.
4. Enable **Auto-scale** to match buffer to sample intensity over a cursor Q range.

| Control | Meaning | Default | Tips |
|---------|------------|---|---|
| **Buffer scale** | Multiplicative factor for buffer (mouse-wheel adjustable) | 1.0 | |
| **Auto-scale** | Compute scale from integral intensity ratio in cursor range | off | |
| **Q min / Q max** | Read-only; drag cursors when Auto-scale is enabled | | Shown only with Auto-scale |

Negative intensity after subtraction is preserved (physically meaningful over-subtraction), not
clamped to zero.

**Batch mode:** with multiple sample files selected, Batch All applies the same buffer subtraction
to each.

**Output suffix:** `_sub`

(data-manipulation-divide)=
#### Divide (structure factor)

Divide datasets to extract structure factors. For concentrated vs. dilute samples:
S(Q) = I_concentrated / I_dilute.

```
D = denom_scale * I_denominator - denom_background
I_out = I_numerator / D
dI_out = |I_out| * sqrt( (dI_num/I_num)^2 + (denom_scale * dI_den / D)^2 )
```

**Workflow:**

1. Right-click a file and choose **Set as denominator**.
2. Double-click the numerator file.
3. Adjust **Denom. scale** and **Denom. bg** as needed.

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Denom. scale** | `denom_scale` in the equation above | 1.0 | |
| **Denom. bg** | `denom_background` in the equation above | 0.0 | |

**Output suffix:** `_div`

(data-manipulation-results)=
## Results and outputs

| Control | What it does |
|---|---|
| **Create Default Folder** | Creates `{data_folder}_manip` next to your data |
| **Select Existing Folder** | Choose any folder |
| **Apply & Save** | Save the currently previewed result to an HDF5 file |
| **Batch All Selected** | Apply the current operation to all selected files |

Results are saved as NXcanSAS HDF5 files in the chosen output folder.

| Source type | Strategy |
|-------------|---------|
| NXcanSAS input | Copy source file, strip pyirena result groups, replace Q/I/Idev/Qdev |
| Text input | Create fresh NXcanSAS file |

All output files include an `entry/data_manipulation_results` NXprocess group recording the
operation name, parameters, source file, and timestamp for full provenance (see also
{ref}`hdf5-results-groups`). An AI agent reads it with
`pyirena_read_manipulation_provenance(path)`.

**Filename convention:** `{input_stem}_{suffix}{ext}`

(data-manipulation-scripting)=
## From Python or AIDA

::::{tab-set}
:::{tab-item} Python
Run from the repository root (outputs go to `output/` in the current folder).

**`manipulate_data()`**

```python
import numpy as np
from pyirena import manipulate_data
from pyirena.io import discover_scattering, load_scattering

result = manipulate_data(
    data_file="validationData/manipulation/manipulate_input.h5",
    operation="scale",            # 'scale', 'trim', 'rebin', 'subtract', 'divide'
    config={"scale_I": 0.4, "background": 0.12},
    buffer_file=None,             # required for 'subtract' and 'divide'
    output_folder="./output",
)
print(result['message'])

# The input is 2.5 × reference + 0.3, so this undoes it exactly:
out = load_scattering(discover_scattering(result['output_file'])[0])
ref = load_scattering(discover_scattering(
    "validationData/manipulation/manipulate_reference.h5")[0])
print(np.max(np.abs(out.intensity / ref.intensity - 1)))   # ~1e-13
```

**`average_data()`**

```python
from pyirena.batch import average_data, manipulate_data

# Two "frames" of the same curve: the reference, and the input scaled back onto it.
scaled = manipulate_data("validationData/manipulation/manipulate_input.h5", "scale",
                         config={"scale_I": 0.4, "background": 0.12},
                         output_folder="./output")["output_file"]
result = average_data(
    data_files=["validationData/manipulation/manipulate_reference.h5", scaled],
    output_folder="./output",
    # optional similarity filtering:
    similarity_check=True,
    similarity_p_min=0.01,
    similarity_reference='first',      # 'first' or 'majority'
    similarity_normalize_scale=True,
)
print(result['message'], result['rejected'])
```

See the [full `average_data` API reference](../guides/batch-api.md#average_data) in
[Batch fitting API](../guides/batch-api.md) for all parameters, return value structure, and
pipeline examples.

`manipulate_data()` returns a dict with keys `success`, `operation`, `output_file`, and `message`,
or `None` on fatal errors. `average_data()` returns the same plus `rejected` (list of
`(filename, p_value)` pairs for frames discarded by the similarity filter).
:::
:::{tab-item} AIDA / MCP
Every operation on this panel is also available headlessly through `pyirena.api`, and over MCP as
the dispatcher's `data` category — so an AI agent can average or subtract without the GUI. Same
core engine, same output naming as described under [Results and outputs](#data-manipulation-results).

```python
from pyirena import api

ref = "validationData/manipulation/manipulate_reference.h5"
inp = "validationData/manipulation/manipulate_input.h5"   # = 2.5 × ref + 0.3
out = "output"
scaled = api.scale_data(inp, scale_I=0.4, background=0.12, output_folder=out)
api.average_data([ref, scaled["output_path"]], similarity_check=True, output_folder=out)
api.subtract_data(inp, ref, buffer_scale=2.5, output_folder=out)   # leaves 0.3
api.divide_data(inp, ref, output_folder=out)
api.trim_data(inp, q_min=0.01, q_max=0.3, output_folder=out)
r = api.rebin_data(inp, mode="log", n_points=200, output_folder=out)
print(r["n_points_in"], r["n_points_written"])
```

Unlike `pyirena.batch`, these return a dict describing what happened rather than `None` on
failure, including `n_points_in` / `n_points_written` / `n_dropped_nonpositive` so a caller can see
when points were stripped (an over-subtraction, most often). Errors come back as
`{"error", "suggestion", "code"}`.

The api layer also adds two guards the panel relies on the user for: it refuses to average
datasets with mismatched slit smearing, and it rejects `auto_scale` without both Q bounds instead
of silently ignoring it.

See [ai_tools_reference.md](../reference/ai-tools-reference.md) for the agent-facing tool
descriptions.
:::
::::

(data-manipulation-method)=
## Method notes

(data-manipulation-interpolation)=
### Interpolation and uncertainty propagation

When combining datasets on different Q grids (average, subtract, divide), one dataset's Q values
are used as the reference grid and the others are interpolated onto it.

**Log-log interpolation.** Intensity is interpolated in log10(I) vs log10(Q) space (linear
interpolation of logarithms), then transformed back. This is appropriate because SAS intensity
typically follows power-law behaviour. Out-of-range points return NaN and are excluded from the
result.

**Uncertainty propagation through interpolation.** The **relative uncertainty** (dI/I) is
interpolated linearly in log10(Q) space, then multiplied by the interpolated intensity. This
preserves the signal-to-noise character of counting-statistics data, where relative uncertainty
varies smoothly with Q.

(data-manipulation-slit-smearing)=
### Slit smearing (USAXS)

Subtract and divide **refuse to mix** a slit-smeared curve with a pinhole one (or two different
slit lengths), since that would silently produce a wrong result; the guard lives in the core engine
so batch scripting inherits it. Manipulation outputs also drop any stale `_SMR` (slit-smeared twin)
entry copied from the source and clear an orphaned `dQl`, so a later slit-smeared load can't return
an inconsistent curve. See **[Slit smearing](../methods/slit-smearing.md)** for the full reference.

% CHECK(Jan): this says the guard for subtract/divide lives in the core engine; the AIDA/MCP text
% says the api layer adds the slit-smearing guard for averaging "the panel relies on the user for".

## References

- Franke, D., Jeffries, C. M. & Svergun, D. I. (2015). Correlation Map, a goodness-of-fit test
  for one-dimensional X-ray scattering spectra. *Nature Methods* **12**, 419–422.
