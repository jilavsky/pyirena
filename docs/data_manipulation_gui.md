# Data Manipulation Tool

The **Data Manipulation** tool provides common post-reduction operations on
SAS datasets: scaling, trimming, rebinning, averaging, buffer subtraction,
and division.  It handles datasets on different Q grids by interpolating in
log-log space and propagates uncertainties through every operation.

The tool is available three ways:

| Mode | How to launch |
|------|--------------|
| Interactive GUI | `pyirena-datamanip` (no arguments), or **Data Manipulation** button in the main hub |
| Python API | `from pyirena.batch import manipulate_data, average_data` |
| Within pyIrena GUI | **Tools > Data Manipulation** menu |

---

## Contents

1. [Quick start](#quick-start)
2. [GUI walkthrough](#gui-walkthrough)
3. [Operations](#operations)
   - [Scale + Background](#scale--background)
   - [Trim Q range](#trim-q-range)
   - [Rebin](#rebin)
   - [Average](#average)
   - [Subtract (buffer)](#subtract-buffer)
   - [Divide (structure factor)](#divide-structure-factor)
4. [Output files](#output-files)
5. [Python API reference](#python-api-reference)
6. [Interpolation and uncertainty propagation](#interpolation-and-uncertainty-propagation)

---

## Quick start

1. Open the tool from the pyIrena hub or run `pyirena-datamanip`.
2. Click **Select Folder** and choose a folder containing your data files.
3. Select the appropriate file type (HDF5 Nexus, HDF5 Generic, or Text).
4. Double-click a file to load it into the graph.
5. Choose an operation tab (Scale, Trim, Rebin, Average, Subtract, Divide).
6. Adjust parameters --- results update automatically.
7. Set an output folder, then click **Apply & Save**.

---

## GUI walkthrough

```
+----------------------------------------------------------------------+
| pyIrena -- Data Manipulation                             [? Help]    |
+---------------+--------------------------------------+---------------+
| Data Files    | [Scale][Trim][Rebin][Avg][Sub][Div]  | Output        |
|               |                                      |               |
| [Select Fldr] | (tab-specific controls)              | Folder:       |
| folder_name   |                                      | (not set)     |
| Type:[HDF5  ] | +---- I(Q) Plot ----+                | [Create Def]  |
| Sort:[Order ] | |  o input data     |                | [Select Fldr] |
| Filter:[    ] | |  -- result curve  |                |               |
|               | |  |  |  cursors    |                | [Apply & Save]|
| file_001.h5   | +-------------------+                | [Batch All]   |
| file_002.h5   |                                      |               |
+---------------+--------------------------------------+---------------+
| Status: Ready                                                        |
+----------------------------------------------------------------------+
```

### File browser (left panel)

- **Select Folder** --- choose the folder containing your data files.
- **Type** --- file format: HDF5 Nexus (NXcanSAS), HDF5 Generic, or Text (.dat/.txt).
- **Sort** --- 10 sort options: filename, temperature, time, order number, pressure
  (ascending or descending).
- **Filter** --- type text to filter filenames (case-insensitive substring match).
- **Double-click** a file to load and plot it.

### Graph (centre)

- Log-log I(Q) plot with the top axis showing R = pi/Q.
- Input data displayed as scatter points; result overlaid as a thick green line.
- Right-click the plot for JPEG export and (on the Average tab) dataset removal.

### Output panel (right)

- **Create Default Folder** --- creates `{data_folder}_manip` next to your data.
- **Select Existing Folder** --- choose any folder.
- **Apply & Save** --- save the currently previewed result to an HDF5 file.
- **Batch All Selected** --- apply the current operation to all selected files.

---

## Operations

All operations update automatically when you change parameters.  There is no
Preview button --- results recalculate within ~500 ms of any change.

### Scale + Background

Multiply intensity by a constant and/or subtract a flat background.

```
I_out = scale * I_in - background
dI_out = scale_uncertainty * dI_in
```

| Control | Description |
|---------|------------|
| **Scale I** | Multiplicative factor (mouse-wheel adjustable, step 0.01) |
| **Background** | Constant subtracted after scaling (mouse-wheel, step 0.001) |
| **Scale uncert.** | Factor for uncertainty; leave empty to use Scale I |

**Output suffix:** `_scaled`

### Trim Q range

Remove data points outside a Q range defined by two draggable cursors.

- Drag the red and blue vertical cursor lines on the plot.
- Q min and Q max update in real time.
- Only data between the cursors is kept.

**Output suffix:** `_trimmed`

### Rebin

Rebin data onto a new Q grid.

| Control | Description |
|---------|------------|
| **Grid** | Log-spaced, Linear, or From reference file |
| **Points** | Number of points in the new grid (10--10000, default 200) |
| **Q min / Q max** | Range for the new grid (auto-filled from data) |
| **Load Reference** | Select "From reference file" first, then click to load a data file whose Q values become the new grid |

Rebinning uses log-log interpolation with relative-uncertainty propagation.

**Output suffix:** `_rebinned`

### Average

Average multiple datasets to improve statistics.  Common in any experiment
where sequential exposures are collected and radiation-damaged frames must
be identified and discarded before averaging.

**Workflow:**

1. Select multiple files in the file list (Ctrl+click or Shift+click).
2. All selected datasets are plotted with distinct colours and averaged
   automatically.
3. Optionally run the **Similarity analysis** to flag damaged frames
   automatically (see below).
4. **Right-click on the graph** for manual removal:
   - **Remove dataset** --- submenu with colour-coded entries matching the plot.
   - **Remove all after** --- remove a dataset and everything measured after it.
5. The average recomputes automatically after each removal.

The average is computed on the first selected file's Q grid.  Other datasets
are log-log interpolated onto that grid before averaging.

**Uncertainty:** the larger of propagated error and sample spread:

```
dI_avg = max( sqrt(sum(dI_i^2)) / N,  std(I_values) )
```

This is robust: propagated error dominates for consistent data, while the
spread correctly reflects systematic drift (e.g. radiation damage).

**Output suffix:** `_avg`

#### Similarity analysis — automatic radiation damage detection

The **Similarity analysis** group inside the Average tab uses the
**CorMap test** (Franke et al., *Nature Methods* 12, 419–422, 2015) to
compute a p-value for each frame.  A low p-value means the frame differs
from the reference beyond what random noise alone can explain — the
hallmark of radiation damage or a transient contamination event (e.g. a
bubble in the beam path).

| Control | Description |
|---------|-------------|
| **Method** | Test algorithm.  Currently *CorMap (Franke 2015)*.  Additional methods can be registered in `pyirena/core/similarity.py`. |
| **Reference** | *First frame* — each dataset is compared to dataset #1, which is always accepted.  Standard bioSAXS approach.  Use when the first exposure is known to be clean.  *Majority vote* — each dataset is compared to the **median** of all datasets; more robust when the first frame may itself be damaged or when damage is not monotonically increasing. |
| **P-value threshold** | Frames with p-value below this are flagged as *Rejected* (default 0.01).  Lower threshold = stricter (fewer rejections); higher = looser (more rejections but cleaner average). |
| **Normalize scale** | Recommended: **on**.  Rescales each frame to match the reference amplitude before comparing, removing flux-drift and absorption differences so that CorMap detects *shape* changes only.  Turn off only when data are on an identical absolute scale and you want scale changes to count as damage. |

**Check Similarity** runs the test and shows a colour-coded results table:
green rows are accepted, red rows are rejected.
The status bar reports `N results for M files` — if N < M, some files
could not be loaded.

**Auto-reject N below threshold** deselects all rejected frames in one
click and re-computes the average.  Manual right-click removal on the
graph still works independently.

> **How to choose a p-value threshold**
> Run the tool on a clean test dataset (frames you know are undamaged) and
> note the lowest p-value you see.  Set the threshold just below that value
> so good frames are accepted.  Then test on a dataset known to contain
> damaged frames and raise the threshold until they are rejected.  Typical
> starting values: 0.001 (strict) or 0.010 (moderate).

### Subtract (buffer)

Subtract buffer/solvent scattering from a sample.  Common in bioSAXS.

```
I_out = I_sample - buffer_scale * I_buffer
dI_out = sqrt( dI_sample^2 + (buffer_scale * dI_buffer)^2 )
```

**Workflow:**

1. Right-click a file in the list and choose **Set as buffer**.
2. Double-click a sample file.  Both sample (blue) and buffer (red) are plotted.
3. Adjust **Buffer scale** with the mouse wheel (step 0.001) or type a value.
4. Enable **Auto-scale** to match buffer to sample intensity over a cursor Q range.

| Control | Description |
|---------|------------|
| **Buffer scale** | Multiplicative factor for buffer (mouse-wheel adjustable) |
| **Auto-scale** | Compute scale from integral intensity ratio in cursor range |
| **Q min / Q max** | Read-only; drag cursors when Auto-scale is enabled |

Negative intensity after subtraction is preserved (physically meaningful
over-subtraction), not clamped to zero.

**Batch mode:** with multiple sample files selected, Batch All applies the
same buffer subtraction to each.

**Output suffix:** `_sub`

### Divide (structure factor)

Divide datasets to extract structure factors.  For concentrated vs. dilute
samples: S(Q) = I_concentrated / I_dilute.

```
D = denom_scale * I_denominator - denom_background
I_out = I_numerator / D
dI_out = |I_out| * sqrt( (dI_num/I_num)^2 + (denom_scale * dI_den / D)^2 )
```

**Workflow:**

1. Right-click a file and choose **Set as denominator**.
2. Double-click the numerator file.
3. Adjust **Denom. scale** and **Denom. bg** as needed.

**Output suffix:** `_div`

---

## Output files

Results are saved as NXcanSAS HDF5 files in the chosen output folder.

| Source type | Strategy |
|-------------|---------|
| NXcanSAS input | Copy source file, strip pyirena result groups, replace Q/I/Idev/Qdev |
| Text input | Create fresh NXcanSAS file |

All output files include an `entry/data_manipulation_results` NXprocess group
recording the operation name, parameters, source file, and timestamp for
full provenance.

**Filename convention:** `{input_stem}_{suffix}{ext}`

---

## Python API reference

### `manipulate_data()`

```python
from pyirena import manipulate_data

result = manipulate_data(
    data_file="sample.h5",
    operation="scale",            # 'scale', 'trim', 'rebin', 'subtract', 'divide'
    config={"scale_I": 2.0, "background": 0.5},
    buffer_file=None,             # required for 'subtract' and 'divide'
    output_folder="./output",
)
```

### `average_data()`

```python
from pyirena.batch import average_data

result = average_data(
    data_files=["frame_001.h5", "frame_002.h5", "frame_003.h5"],
    output_folder="./output",
    # optional similarity filtering:
    similarity_check=True,
    similarity_p_min=0.01,
    similarity_reference='first',      # 'first' or 'majority'
    similarity_normalize_scale=True,
)
```

See the [full `average_data` API reference](batch_api.md#average_data) in
`batch_api.md` for all parameters, return value structure, and pipeline
examples.

`manipulate_data()` returns a dict with keys `success`, `operation`,
`output_file`, and `message`, or `None` on fatal errors.  `average_data()`
returns the same plus `rejected` (list of `(filename, p_value)` pairs for
frames discarded by the similarity filter).

---

## Slit smearing (USAXS)

Subtract and divide **refuse to mix** a slit-smeared curve with a pinhole one
(or two different slit lengths), since that would silently produce a wrong
result; the guard lives in the core engine so batch scripting inherits it.
Manipulation outputs also drop any stale `_SMR` (slit-smeared twin) entry copied
from the source and clear an orphaned `dQl`, so a later slit-smeared load can't
return an inconsistent curve. See **[Slit smearing](slit_smearing.md)** for the
full reference.

## Interpolation and uncertainty propagation

When combining datasets on different Q grids (average, subtract, divide),
one dataset's Q values are used as the reference grid and the others are
interpolated onto it.

### Log-log interpolation

Intensity is interpolated in log10(I) vs log10(Q) space (linear
interpolation of logarithms), then transformed back.  This is appropriate
because SAS intensity typically follows power-law behaviour.  Out-of-range
points return NaN and are excluded from the result.

### Uncertainty propagation through interpolation

The **relative uncertainty** (dI/I) is interpolated linearly in log10(Q)
space, then multiplied by the interpolated intensity.  This preserves the
signal-to-noise character of counting-statistics data, where relative
uncertainty varies smoothly with Q.
