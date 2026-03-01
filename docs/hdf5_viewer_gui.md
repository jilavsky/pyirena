# HDF5 Viewer / Data Extractor — GUI Guide

The **HDF5 Viewer / Data Extractor** lets you browse HDF5 files produced by
pyirena (and compatible software), plot 1D data such as I(Q) curves and fit
models, and collect scalar parameters (e.g., Rg, peak Q₀, χ²) across a series
of files into a table and scatter plot.

---

## Opening the tool

### From the command line (standalone)

```bash
pyirena-viewer
```

If you installed pyirena with `pip install -e .` (or a regular install) the
`pyirena-viewer` command is available.  The viewer opens with an empty folder;
use **Select Folder** in the left panel to load a directory.

> **First-time setup:** if `pyirena-viewer` is not found after upgrading, run
> `pip install -e .` from the repository root to register the new entry point.

### From the Data Selector

Click **HDF5 Viewer** in the tool-launch row of the Data Selector.  The viewer
opens pre-loaded with the folder that is currently active in the Data Selector.

---

## Window layout

```
┌──────────────────────────────────────────────────────────────────────────┐
│  pyIrena HDF5 Viewer / Data Extractor                                    │
├────────────────────┬───────────────────────┬─────────────────────────────┤
│  File Select       │  Data Select          │  Plot Controls              │
│  (left panel)      │  (centre panel)       │  (right panel)              │
│  ≈260 px           │  ≈320 px              │  flex width                 │
└────────────────────┴───────────────────────┴─────────────────────────────┘
```

The three panels are separated by draggable splitters.

---

## Left panel — File Select

### Selecting a folder

Click **Select Folder** to open a folder browser.  The viewer scans the chosen
directory for HDF5 files (`.h5`, `.hdf5`, `.hdf`, `.nxs`) and lists them in
the tree.  Subdirectories that contain at least one HDF5 file appear as
**bold-italic** folder nodes; they expand lazily on first click.

The folder is remembered across sessions.

### Sort order

Use the **Sort** combo box to order files by:

| Option | Description |
|--------|-------------|
| Filename A→Z / Z→A | Alphabetical |
| Temperature ↑ / ↓ | Numeric value of `_XXC_` token in filename |
| Time ↑ / ↓ | Numeric value of `_XXmin_` token in filename |
| Order number ↑ / ↓ | Trailing numeric index in filename (default) |
| Pressure ↑ / ↓ | Numeric value of `_XXPSI_` token in filename |

### Filtering files

Type in the **Filter** box to restrict which files are shown.  The filter
supports **regular expressions** (regex / grep syntax):

| Filter text | Matches |
|-------------|---------|
| `60C` | any filename containing "60C" |
| `60C\|0[12]min` | filenames with "60C" **or** "01min" / "02min" |
| `^sample` | filenames starting with "sample" |
| `\.h5$` | filenames ending with ".h5" |
| Plain text | plain substring match (case-insensitive) |

If the filter text is not valid regex it falls back to plain substring matching.

### Selecting files

- Click a file to select it (shown in the centre panel and used by Plot Controls).
- **Ctrl+click** / **Shift+click** to select multiple files.
- All selected files are available for "plot all" and "collect values" operations.

The file count is shown below the filter box.

---

## Centre panel — Data Select

### File name display

The name of the currently browsed file appears at the top in **bold black**
text (full path shown in tooltip).

### HDF5 tree

The tree shows the internal structure of the first selected file.  Groups
(HDF5 groups) appear as **bold** items with an expand triangle; datasets appear
as normal items showing their shape and dtype in the second column.  Scalar
datasets also show their value inline (e.g., `scalar  [float64] = 1.234`).

Key HDF5 attributes (NX_class, canSAS_class, signal, units, …) appear as
indented italic `@attribute` sub-items.

Children are loaded **lazily** — the tree only reads from disk when you expand
a group, keeping large files responsive.

### Scalar value display

Clicking any item shows its value in the read-only field at the bottom of the
panel:

- **Scalar dataset** — e.g., `chi_squared = 0.00347`
- **Attribute** — e.g., `@NX_class = NXentry`

Array datasets do not show a value (use Plot Controls to plot them).

### Right-click on a dataset

| Action | Effect |
|--------|--------|
| Add as Y axis | Sets as Y data in Tab 1 Custom |
| Add as X axis | Sets as X data in Tab 1 Custom |
| Add as Y error | Sets as Y error in Tab 1 Custom |
| Add as X error | Sets as X error in Tab 1 Custom |
| Collect value across selected files | Pre-fills the Collect Values tab with this HDF5 path |
| Set as X-axis metadata path | Pre-fills the X-axis metadata path in Tab 2 |

### Right-click on a group

| Action | When available |
|--------|---------------|
| Plot NXcanSAS data (I vs Q) | Group has `@canSAS_class = SASentry` |
| Plot Unified Fit model | Group is `unified_fit_results` |
| Plot Size Distribution model | Group is `sizes_results` |
| Plot WAXS Peak Fit model | Group is `waxs_peakfit_results` |
| Plot Simple Fit model | Group is `simple_fit_results` |

These actions read the data immediately and open a new Graph Window.

---

## Right panel — Plot Controls

### Tab 1 — 1D Graph

#### Pyirena presets

Check one or more boxes to include the corresponding data in the next graph:

| Checkbox | Data plotted |
|----------|-------------|
| NXcanSAS (I vs Q) | Reduced scattering data + uncertainties |
| Unified Fit model | Fit model curve; data overlay if available |
| Size Distribution I(Q) | Model I(Q) from size distribution |
| Size Distribution P(r) | Pair-distance distribution (separate graph) |
| WAXS Peak Fit | Total fit curve |
| Simple Fit model | Model curve from simple fit |

Size Distribution P(r) always opens in a **separate** Graph Window because its
axes are incompatible with I(Q) data.

#### Custom data (from HDF5 browser)

Right-click datasets in the centre panel to populate X / Y / Y-error / X-error
paths.  Use **Clear** to reset.

#### Source

- **First selected file** — plot data from the file shown in the centre panel only.
- **All selected files** — loop over all files selected in File Select and add
  one curve per file to the same graph, labelled by filename.

#### Buttons

| Button | Action |
|--------|--------|
| New Graph | Open a new Graph Window with the current curves |
| Add to active graph | Add curves to the most recently used Graph Window |

---

### Tab 2 — Collect Values

Collects a **single scalar** from each selected file and assembles a table and
scatter plot.

#### What to collect

Choose a **Type** and **Item** from the combo boxes:

| Type | Item examples |
|------|--------------|
| Unified Fit | Rg, G, B, P, ETA, PACK (per level); background; chi2 |
| WAXS Peak Fit | Q0, A, FWHM (per peak); chi2 |
| Size Distribution | chi2, volume_fraction |
| Simple Fits | any fitted parameter by name; chi2 |
| Custom HDF5 path | any dataset or `path@attribute` in the file |

When MC uncertainty has been calculated, the standard deviation is collected
automatically and shown as error bars.

For **Unified Fit** and **WAXS Peak Fit**, use the **Level** / **Peak**
spin box to select which level or peak to read from.

#### X axis

| Option | X value |
|--------|---------|
| File order (1, 2, 3…) | Sequential integer position |
| Filename sort key | Numeric value extracted from filename (temperature, time, etc.) |
| HDF5 metadata path | Scalar read from that HDF5 path in each file |

#### Collect button

Click **Collect from all selected files** to open a **Collect Window** showing
the table and scatter plot.

---

## Graph Window

Each "New Graph" opens a floating Graph Window.  Multiple windows can be open
simultaneously.  The most recently activated window is the **active graph** for
"Add to active graph".

### Toolbar

| Button | Action |
|--------|--------|
| X log / X lin | Toggle X axis between log₁₀ and linear scale |
| Y log / Y lin | Toggle Y axis between log₁₀ and linear scale |
| Labels… | Open dialog to set graph title, X label, Y label |
| Legend | Show / hide the curve legend |
| Save JPEG | Save graph image at 1600 px width |
| Save CSV | Save all curves as a comma-separated text file |

### Right-click menu (on graph area)

| Action | Description |
|--------|-------------|
| Set X range… | Enter explicit min/max for the X axis |
| Set Y range… | Enter explicit min/max for the Y axis |
| Curve styles… | Per-curve color, line width, and symbol picker |
| Remove error bars | Remove all error bar overlays from the graph |
| Save PNG… | Save as PNG image |
| Save HDF5… | Export all curves as NXdata groups in an HDF5 file |
| Save ITX (Igor Pro)… | Export in Igor Pro text wave format |
| Open as matplotlib figure… | Transfer curves to a matplotlib window |

### Adding curves

Use **Add to active graph** in the Plot Controls to append more curves to an
existing window without opening a new one.  The active window is highlighted
when clicked.

---

## Collect Window

Opened by the **Collect** button in Tab 2.  Shows:

- **Table** — one row per file: filename, X value, collected Y value, ±error.
- **Scatter plot** — X vs Y with error bars (where available).

### Toolbar

| Button | Action |
|--------|--------|
| Save CSV | Save table as a comma-separated file (default: window title + `.csv` in CWD) |
| Save JPEG | Save the full window (table + plot) as a JPEG image |

---

## State persistence

The viewer saves and restores these settings automatically across sessions:

| Setting | Saved |
|---------|-------|
| Last folder | Yes |
| Sort order | Yes |
| Filter text | Yes |

---

## Supported file types

| Extension | Format |
|-----------|--------|
| `.h5`, `.hdf5`, `.hdf` | Generic HDF5 |
| `.nxs` | NeXus / NXcanSAS |

The viewer can browse **any** HDF5 file.  Pyirena-specific plot actions
(NXcanSAS, Unified Fit, etc.) require the pyirena result groups to be present
in the file.
