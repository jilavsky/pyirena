(data-explorer)=
# Data Explorer

:::{admonition} At a glance
:class: tip

**Launch:** `pyirena-gui` → **Tools → Data Explorer**, or the **Data Explorer** button (the
viewer opens pre-loaded with the folder that is currently active in the Data Selector).
Standalone: `pyirena-viewer` (opens with an empty folder; use **Select Folder** in the left
panel to load a directory). \
**Python:** `pyirena.api` read functions (`inspect_file`, `read_reduced_data`,
`read_unified_fit`, …, `tabulate_parameter`). \
**MCP:** `pyirena_inspect_file`, `pyirena_read_reduced_data`, `pyirena_read_<tool>`,
`pyirena_tabulate_parameter`, `pyirena_plot_iq`, `pyirena_plot_parameter_trend`. \
**Igor equivalent:** {doc}`Irena HDF5 Browser <igor:Irena/HDF5Browser>`
:::

> **First-time setup:** if `pyirena-viewer` is not found after upgrading, run
> `pip install -e .` from the repository root to register the new entry point.

The **Data Explorer** lets you browse HDF5 files produced by
pyirena (and compatible software), plot 1D data such as I(Q) curves and fit
models, and collect scalar parameters (e.g., Rg, peak Q₀, χ²) across a series
of files into a table and scatter plot.

The viewer supports results from all pyirena analysis tools: **Unified Fit**
(including Guinier analysis, Porod's law exponents, and fractal/agglomerate
models), **Size Distribution** (volume/number distributions, cumulative distributions, and I(Q) models),
**WAXS Peak Fit** (diffraction peak fitting), **Simple Fits** (Guinier, Porod,
Debye-Bueche, Teubner-Strey, and other models), and **Data Merge**
(merged SAXS/USAXS datasets).

% CHECK(Jan): "Data Merge (merged SAXS/USAXS datasets)" — there is no Data Merge preset or
% Collect type in the current panel; merged files are plotted as NXcanSAS data.

(data-explorer-window)=
## The window

```{figure} img/data-explorer/data-explorer-01-window.png
:alt: Data Explorer with the folder "unified" open and three files selected; the centre panel shows the HDF5 tree of unified_two_level.h5 with the unified_fit_results group expanded; the right panel shows the 1D Graph tab with the NXcanSAS and Unified Fit model presets checked. Numbered circles mark the parts described in the legend.

The Data Explorer on a folder of three Unified Fit results (copies of
`validationData/unified/`, each fitted and stored).
```

The window has three panels — **File Select** (left, ≈260 px), **Data Select** (centre,
≈320 px) and **Plot Controls** (right, flexible width) — separated by draggable splitters.

1. **Select Folder** and the name of the folder in use.
2. **Sort** order of the file list.
3. **Filter** (regular expressions allowed) and the file count.
4. **File list**: click, Ctrl+click or Shift+click to select files.
5. **HDF5 tree** of the first selected file; right-click for plot and collect actions.
6. **Value field**: the value of the selected scalar dataset or attribute.
7. **Plot Controls tabs**: 1D Graph, Collect Values, Multi-Collect, Export to Igor.
8. **Pyirena presets**: what the next graph shows.
9. **New Graph** and **Add to active graph**.
10. **Status bar**.

(data-explorer-workflow)=
## Typical workflow

The steps use a folder named `unified` holding copies of the three files in
`validationData/unified/`, each with a Unified Fit stored in it (fit them in the
[Unified Fit](unified-fit.md) tool and click **Store in File**).

1. **Open the folder.** Click **Select Folder** and choose `unified`. The three files appear in
   the list.
2. **Select the files.** Click the first file and Shift+click the last; the status bar reports
   "3 files selected", and the centre panel shows the tree of the first one.
3. **Browse a file.** Expand `entry` and `unified_fit_results` in the centre panel; click a
   scalar such as `chi_squared` to see its value in the field below the tree.
4. **Plot all curves.** On the **1D Graph** tab keep **NXcanSAS (I vs Q)** and **Unified Fit
   model** checked, choose **All selected files** under Source and click **New Graph**.

   ```{figure} img/data-explorer/data-explorer-02-graph.png
   :alt: Graph window with three I(Q) curves from the folder unified and the Unified Fit model of each, on log-log axes, with the X Log, Y Log, Labels, Legend, Offset, Save JPEG and Save CSV toolbar.

   **New Graph** with *All selected files*: data and Unified Fit model for each file.
   ```

5. **Collect a parameter.** On the **Collect Values** tab choose Type *Unified Fit*, Item
   *Rg*, Level 1, X axis *File order (1, 2, 3…)*, and click **Collect from all selected
   files**.

   ```{figure} img/data-explorer/data-explorer-03-collect.png
   :alt: Collect window with a table of the level-1 Rg of the three files (about 120, 200 and 200 Å) and a scatter plot of Rg against file order.

   The level-1 Rg of the three files, as table and scatter plot.
   ```

6. **Take the numbers out.** Select cells and press Ctrl+C, or use **Save CSV** / **Save ITX**
   in the Collect window (see [Results and outputs](#data-explorer-results)).

(data-explorer-controls)=
## Controls

(data-explorer-file-select)=
### Left panel — File Select

**Selecting a folder.** Click **Select Folder** to open a folder browser. The viewer scans the
chosen directory for HDF5 files (`.h5`, `.hdf5`, `.hdf`, `.nxs`, `.h5xp`) and lists them in
the tree. Subdirectories that contain at least one HDF5 file appear as
**bold-italic** folder nodes; they expand lazily on first click.

The folder is remembered across sessions.

**Sort order.** Use the **Sort** combo box to order files by:

| Option | Description |
|--------|-------------|
| Filename A→Z / Z→A | Alphabetical |
| Temperature ↑ / ↓ | Numeric value of `_XXC_` token in filename |
| Time ↑ / ↓ | Numeric value of `_XXmin_` token in filename |
| Order number ↑ / ↓ | Trailing numeric index in filename (default) |
| Pressure ↑ / ↓ | Numeric value of `_XXPSI_` token in filename |

The same dropdown, patterns and tooltip appear in **all four** file browsers
(Data Selector, Data Explorer, Data Manipulation, Data Merge). Matching is
case-insensitive, and files whose name does not carry the pattern — logs,
notes, a stray average — sort **last** in both directions, so they never break
up an ordered series.

**Filtering files.** Type in the **Filter** box to restrict which files are shown. The filter
supports **regular expressions** (regex / grep syntax):

| Filter text | Matches |
|-------------|---------|
| `60C` | any filename containing "60C" |
| `60C\|0[12]min` | filenames with "60C" **or** "01min" / "02min" |
| `^sample` | filenames starting with "sample" |
| `\.h5$` | filenames ending with ".h5" |
| Plain text | plain substring match (case-insensitive) |

If the filter text is not valid regex it falls back to plain substring matching.

**Selecting files.**

- Click a file to select it (shown in the centre panel and used by Plot Controls).
- **Ctrl+click** / **Shift+click** to select multiple files.
- All selected files are available for "plot all" and "collect values" operations.

The file count is shown below the filter box.

(data-explorer-data-select)=
### Centre panel — Data Select

**File name display.** The name of the currently browsed file appears at the top in **bold
black** text (full path shown in tooltip).

**HDF5 tree.** The tree shows the internal structure of the first selected file. Groups
(HDF5 groups) appear as **bold** items with an expand triangle; datasets appear
as normal items showing their shape and dtype in the second column. Scalar
datasets also show their value inline (e.g., `scalar  [float64] = 1.234`).

Key HDF5 attributes (NX_class, canSAS_class, signal, units, …) appear as
indented italic `@attribute` sub-items.

Children are loaded **lazily** — the tree only reads from disk when you expand
a group, keeping large files responsive.

**Scalar value display.** Clicking any item shows its value in the read-only field at the
bottom of the panel:

- **Scalar dataset** — e.g., `chi_squared = 0.00347`
- **Attribute** — e.g., `@NX_class = NXentry`

Array datasets do not show a value (use Plot Controls to plot them).

**Right-click on a dataset**

| Action | Effect |
|--------|--------|
| Add as Y axis | Sets as Y data in Tab 1 Custom |
| Add as X axis | Sets as X data in Tab 1 Custom |
| Add as Y error | Sets as Y error in Tab 1 Custom |
| Add as X error | Sets as X error in Tab 1 Custom |
| Collect value across selected files | Pre-fills the Collect Values tab with this HDF5 path |
| Add to Multi-Collect list | Appends this HDF5 path to the item list in Tab 3 (Multi-Collect) |
| Set as X-axis metadata path | Pre-fills the X-axis metadata path in Tab 2 |

**Right-click on a group**

| Action | When available |
|--------|---------------|
| Plot NXcanSAS data (I vs Q) | Group has `@canSAS_class = SASentry` |
| Plot Unified Fit model | Group is `unified_fit_results` |
| Plot Size Distribution model | Group is `sizes_results` |
| Plot WAXS Peak Fit model | Group is `waxs_peakfit_results` |
| Plot Simple Fit model | Group is `simple_fit_results` |

These actions read the data immediately and open a new Graph Window.

% CHECK(Jan): the group menu also offers "Plot Modeling model" (not in the table).

(data-explorer-graph-tab)=
### Right panel, Tab 1 — 1D Graph

**Pyirena presets.** Check one or more boxes to include the corresponding data in the next
graph:

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **NXcanSAS (I vs Q)** | Reduced scattering data + uncertainties | on | |
| **Unified Fit model** | Fit model curve; data overlay if available | on | |
| **Size Dist. model (I vs Q)** | Model I(Q) from size distribution | off | |
| **Size Dist. vol. P(r)** | Volume size distribution vs r (separate graph) | off | |
| **Size Dist. num. N(r)** | Number size distribution vs r (separate graph) | off | |
| **Size Dist. cumul. vol.** | Cumulative volume distribution vs r (separate graph) | off | |
| **Size Dist. cumul. num.** | Cumulative number distribution vs r (separate graph) | off | |
| **WAXS model** | Total fit curve | off | Called "WAXS Peak Fit" in the old text |
| **Simple Fit model** | Model curve from simple fit | off | |

Size Distribution distribution plots (P(r), N(r), cumulative) always open in a
**separate** Graph Window because their axes are incompatible with I(Q) data.

% CHECK(Jan): the panel also has "Size Dist. surf. S(r)", "Size Dist. cumul. surf.",
% "Modeling model", "Modeling vol. P(r)" and "Modeling num. P(r)", not described here.

**Custom data (from HDF5 browser).** Right-click datasets in the centre panel to populate
X / Y / Y-error / X-error paths. Use **Clear** to reset.

**Source**

- **First selected file** — plot data from the file shown in the centre panel only.
- **All selected files** — loop over all files selected in File Select and add
  one curve per file to the same graph, labelled by filename.

| Button | Action |
|--------|--------|
| **New Graph** | Open a new Graph Window with the current curves |
| **Add to active graph** | Add curves to the most recently used Graph Window |

(data-explorer-collect-tab)=
### Right panel, Tab 2 — Collect Values

Collects a **single scalar** from each selected file and assembles a table and
scatter plot.

**What to collect.** Choose a **Type** and **Item** from the combo boxes:

| Type | Item examples |
|------|--------------|
| Unified Fit | Rg, G, B, P, ETA, PACK (per level); background; chi2. Covers Guinier knee parameters and Porod law exponents. |
| WAXS Peak Fit | Q0, A, FWHM (per peak); chi2. Suitable for diffraction peak fitting results. |
| Size Distribution | chi2, volume_fraction. Covers volume-weighted size distribution models. |
| Simple Fits | any fitted parameter by name; chi2. Covers Guinier, Porod, Debye-Bueche, Teubner-Strey, and other models. |
| Custom HDF5 path | any dataset or `path@attribute` in the file |

% CHECK(Jan): the Type list also has "Modeling" and "Carbon model".

When MC uncertainty has been calculated, the standard deviation is collected
automatically and shown as error bars.

For **Unified Fit** and **WAXS Peak Fit**, use the **Level** / **Peak**
spin box to select which level or peak to read from.

**X axis**

| Option | X value |
|--------|---------|
| File order (1, 2, 3…) | Sequential integer position |
| Filename sort key | Numeric value extracted from filename (temperature, time, etc.) |
| HDF5 metadata path | Scalar read from that HDF5 path in each file |

Click **Collect from all selected files** to open a **Collect Window** showing
the table and scatter plot.

(data-explorer-multi-collect-tab)=
### Right panel, Tab 3 — Multi-Collect

Collects **multiple HDF5 items simultaneously** from each selected file and
assembles them into a single wide table. Unlike Tab 2 (one scalar per file),
Multi-Collect lets you extract 8 or more items per file in one pass — useful
for ion-chamber readings, detector counts, or any set of related quantities.

**Item list.** Each row in the item list specifies one quantity to collect:

| Column | Description |
|--------|-------------|
| Path | HDF5 path inside the file (e.g., `entry/monitor/counts`) |
| Label | Column header in the output table |
| Reduction | How to reduce an array to a scalar: **sum**, **mean**, **min**, or **max** |

**Adding items:**

- Type a path and label in the text fields below the list and click **Add item**,
  OR
- Right-click any dataset in the centre panel and choose **Add to Multi-Collect
  list** — the path is inserted automatically and the viewer switches to this tab.

Use **Remove selected** to delete highlighted rows; **Clear all** removes the
entire list.

% CHECK(Jan): the buttons are labelled "+ Add", "Remove sel." and "Clear all", and the
% collect button "Collect from all selected files" (not "Multi-Collect from all selected files").

Both scalars and arrays are supported. String datasets are read as-is.
Arrays are reduced to a single number using the chosen Reduction function
(`sum` is the default and is appropriate for raw detector counts).

The same X-axis options as Tab 2 are available (file order, filename sort key, HDF5 metadata
path).

Click **Multi-Collect from all selected files** to open a **Multi-Collect
Window** showing the wide table (one column per item).

% CHECK(Jan): a fourth tab, "Export to Igor", exists in the panel and is not described on this
% page (Igor export is described on the Data Selector page).

(data-explorer-graph-window)=
### Graph Window

Each "New Graph" opens a floating Graph Window. Multiple windows can be open
simultaneously. The most recently activated window is the **active graph** for
"Add to active graph".

**Toolbar**

| Button | Action |
|--------|--------|
| X log / X lin | Toggle X axis between log₁₀ and linear scale |
| Y log / Y lin | Toggle Y axis between log₁₀ and linear scale |
| Labels… | Open dialog to set graph title, X label, Y label |
| Legend | Show / hide the curve legend |
| Offset… | Offset curves into a waterfall / stacked view (see below) |
| Save JPEG | Save graph image at 1600 px width |
| Save CSV | Save all curves as a comma-separated text file |

**Right-click menu (on graph area)**

| Action | Description |
|--------|-------------|
| Set X range… | Enter explicit min/max for the X axis |
| Set Y range… | Enter explicit min/max for the Y axis |
| Curve styles… | Per-curve color, line width, and symbol picker |
| Offset curves… | Offset curves into a waterfall / stacked view (see below) |
| Remove error bars | Remove all error bar overlays from the graph |
| Add slope guide line ▸ | Add a draggable power-law guide line at a preset slope (n = −2, −3, −4 Porod, −5, −6) or a **Custom slope…** value; lines are draggable and right-clicking one lets you change or remove it |
| Remove all slope lines | Remove every slope guide line from the graph |
| Save PNG… | Save as PNG image |
| Save HDF5… | Export all curves as NXdata groups in an HDF5 file |
| Save ITX (Igor Pro)… | Export in Igor Pro text wave format |
| Open as matplotlib figure… | Transfer curves to a matplotlib window |

Slope guide lines are most useful on log–log I(Q) plots for comparing the
local power-law slope of the data against reference exponents (e.g., −4 for
Porod scattering).

**Adding curves.** Use **Add to active graph** in the Plot Controls to append more curves to
an existing window without opening a new one. The active window is highlighted
when clicked.

(data-explorer-offsets)=
#### Offsetting curves (waterfall / stacked view)

Open the **Offset…** dialog from the toolbar button or the **Offset curves…**
right-click action to separate overlapping curves into a stacked view. Because
the meaning of an offset differs by axis type, each axis has an independent
**offset type**:

| Offset type | Effect | Best for |
|-------------|--------|----------|
| Additive (+) | `value + offset` | linear axes |
| Multiplicative (×) | `value × offset` (a constant *visual* shift on a log axis) | log axes |

When the dialog opens, the offset type defaults to match the current axis scale
(log → multiplicative, linear → additive); you can override it per axis.

The dialog provides:

- **Auto-stagger** — enter a single increment and click **Apply** to spread all
  curves automatically: curve *i* gets `i × increment` (additive) or
  `increment` ⁱ (multiplicative). The first curve stays in place.
- **Per-curve offsets** — a table to fine-tune the X and Y offset of any
  individual curve.
- **Clear offsets** — reset everything to no offset.

Offsets are **display-only**: they change only how curves are drawn. The status
line shows `[offsets applied]` while any offset is active, and all exports
(CSV, HDF5, ITX, matplotlib) always write the original, unshifted data.

(data-explorer-results)=
## Results and outputs

The Data Explorer never changes the files it reads. What it produces are graphs, tables and
the files you export from them.

### Collect Window

Opened by the **Collect** button in Tab 2. Shows:

- **Table** — one row per file: filename, X value, collected Y value, ±error.
- **Scatter plot** — X vs Y with error bars (where available).

The table supports **clipboard copy and column sorting** — see
[Working with tables](#data-explorer-tables) below. For a few numbers there is
no need to go through a CSV file: select the cells and press Ctrl+C.

Right-clicking the **plot** gives the standard graph export menu (clipboard
copy, PNG/JPEG/SVG image, curve CSV) plus this window's own *Save collected
values as ITX*, which writes the collected X/Y/error waves with the right axis
labels.

| Button | Action |
|--------|--------|
| Save CSV | Save table as a comma-separated file (default: window title + `.csv` in CWD) |
| Save JPEG | Save the full window (table + plot) as a PNG or JPEG image |

% CHECK(Jan): the Collect window toolbar also has a "Save ITX" button.

### Multi-Collect Window

Opened by the **Multi-Collect** button in Tab 3. Shows a wide table with one
row per file and one column per collected item:

```
File | X value | Item 1 | Item 2 | … | Item N
```

The status bar reports how many files yielded at least one non-empty value and
how many items were collected.

The table supports **clipboard copy and column sorting** — see
[Working with tables](#data-explorer-tables) below.

| Button | Action |
|--------|--------|
| Save CSV | Save table as a comma-separated file (default: window title + `.csv` in CWD) |
| Save JPEG | Save the full window as a PNG or JPEG image |

(data-explorer-tables)=
### Working with tables

Every table in pyIrena behaves the same way (Collect, Multi-Collect, the
Multi-Collect item list, Data Selector's Tabulate Results, Data Manipulation's
similarity results, and the Scattering Contrast tables):

| Action | How |
|--------|-----|
| Copy selected cells | **Ctrl+C** — tab-separated text, pastes straight into Excel, Igor Pro or Origin |
| Copy with column headers | **Ctrl+Shift+C**, or right-click → *Copy with Column Headers* |
| Copy the whole table | Right-click → *Copy Whole Table* (also what Ctrl+C does when nothing is selected) |
| Save as CSV | Right-click → *Save as CSV…*, or the toolbar/panel button where one exists |
| Sort | Click a column header; click again to reverse. Numeric columns sort numerically (9 before 10 before 100), and blank or "—" cells sort last |

Selections can be single cells, rectangular blocks, whole rows or columns
(click the header), or non-contiguous ctrl-clicked cells; a non-contiguous
selection is copied as a compact block of just the selected rows and columns.

Sorting is deliberately **not** offered where the row order carries meaning —
the Contrast results table (rows are grouped under section headers) and the
Multi-Collect item list (its order defines the output columns).

### State persistence

The viewer saves and restores these settings automatically across sessions:

| Setting | Saved |
|---------|-------|
| Last folder | Yes |
| Sort order | Yes |
| Filter text | Yes |

(data-explorer-scripting)=
## From Python or AIDA

::::{tab-set}
:::{tab-item} Python
The `pyirena.api` read functions return the same information as plain dictionaries. Run from
the repository root:

```python
from pyirena import api

path = "validationData/unified/unified_two_level.h5"
info = api.inspect_file(path)
print(info["name"], info["analyses_present"])    # [] : no fit stored in this file yet

data = api.read_reduced_data(path)               # Q, I, dI as lists
print(data["n_points"], data["q_min"], data["q_max"], data["I_units"])
```

For a folder of stored results, `api.tabulate_parameter(folder, "unified_fit", "Rg",
subgroup_index=1)` returns one row per file — the scripted counterpart of **Collect Values** —
and `api.read_unified_fit(path)` (and the other `read_*` functions) one file's results.
:::
:::{tab-item} AIDA / MCP
Ask in plain language, for example *"Tabulate Rg of level 1 for every file in this folder and
plot the trend."* The agent uses `pyirena_summarize_folder`, `pyirena_tabulate_parameter`
and `pyirena_plot_parameter_trend`; `pyirena_inspect_file` and `pyirena_read_reduced_data`
look inside one file. See the
[AI tools reference](../reference/ai-tools-reference.md).
:::
::::

(data-explorer-method)=
## Method notes

The Data Explorer computes nothing itself: it reads what the analysis tools stored. The layout
of those results groups is described in {ref}`hdf5-results-groups`.

(data-explorer-supported)=
## Supported files and analyses

| Extension | Format |
|-----------|--------|
| `.h5`, `.hdf5`, `.hdf` | Generic HDF5 |
| `.nxs` | NeXus / NXcanSAS |

The viewer can browse **any** HDF5 file. Pyirena-specific plot actions
(NXcanSAS, Unified Fit, etc.) require the pyirena result groups to be present
in the file.

% CHECK(Jan): this table omits `.h5xp`, which the file list does scan (see "Selecting a folder").

The following analysis methods are supported for result viewing and parameter
collection via the HDF5 Viewer:

| Analysis | Keyword / technique |
|----------|-------------------|
| Unified Fit | Guinier analysis, Porod's law (power-law exponent), fractal dimension, Beaucage model |
| Size Distribution | size distribution, volume/number P(r), cumulative distributions, particle size |
| WAXS Peak Fit | diffraction peak fitting, d-spacing, FWHM |
| Simple Fits | Guinier, Porod, Debye-Bueche (two-phase model), Teubner-Strey (microemulsion / lamellar), power law, Lorentzian, Gaussian |
| Data Merge | SAXS/USAXS merge, scale matching, background subtraction |

## References

None; the file layout is in {ref}`hdf5-results-groups`.
