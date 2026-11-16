(data-selector)=
# Data Selector

:::{admonition} At a glance
:class: tip

**Launch:** `pyirena-gui`, or `python -m pyirena.gui.launch`. It is the main pyIrena window:
every other tool is opened from it. \
**Python:** `pyirena.api.list_files`, `summarize_folder`, `plot_iq`; the **(script)** buttons
run `pyirena.batch` (`fit_pyirena(data_file, config_file)`). \
**MCP:** `pyirena_list_files`, `pyirena_summarize_folder`, `pyirena_plot_iq`,
`pyirena_tabulate_parameter`. \
**Igor equivalent:** none, new in pyIrena.
:::

Graphical user interface for pyIrena data analysis. The main GUI provides:

1. **Folder Selection** - Browse and select data directories
2. **File Type Filtering** - Choose between:
   - HDF5 files (.hdf, .h5, .hdf5)
   - Text files (.txt, .dat, .csv) — automatically converted to cleaned NXcanSAS
     HDF5 on first use (see [Text File Import and Cleaning](../guides/data-import-and-cleaning.md))
   - All supported files
3. **File List** - View and select files from the folder
4. **Text Filter** - Grep-like filtering of file names, with full regular
   expression support
5. **Multi-selection** - Select multiple files (Ctrl/Cmd + click)
6. **Graphing** - Plot selected files with:
   - Log-log scale
   - Multiple datasets overlaid
   - Legend and labels

This page also describes behaviour shared by every pyIrena window:
[getting graphs and tables out](#getting-graphs-and-tables-out),
[opening data by dragging it in](#data-selector-drag-and-drop) and
[windows that remember where they were](#data-selector-window-geometry).

(data-selector-window)=
## The window

```{figure} img/data-selector/data-selector-01-window.png
:alt: pyIrena Data Selector with the folder /tmp/pyirena-docs/unified open and three files selected; on the right the Show in graph / reports checkboxes, the Create Graph, Create Report, Tabulate Results, Export to ASCII and Data Explorer buttons, the Analysis Tools buttons and the Data Processing & Reference buttons. Numbered circles mark the parts described in the legend.

The Data Selector on a folder of three stored Unified Fits.
```

1. **Select Data Folder** and **Refresh**; the folder in use is shown beside them.
2. **File Type**: HDF5, text or all supported files.
3. **Sort** order of the file list.
4. **Filter** (regular expressions allowed).
5. **File list**: click, Ctrl/Cmd+click, Shift+click; double-click plots.
6. **Configure…**, **Manage Config…**, **Locate Logs…**.
7. **Show in graph / reports**: which stored results the buttons below use.
8. **Create Graph**, **Create Report**, **Tabulate Results**, **Export to ASCII**.
9. **Data Explorer**.
10. **Analysis Tools**: each tool's window (GUI) or a batch fit of all selected files (script).
11. **Data Processing & Reference**: Data Merge, Data Manipulation, Scattering Contrast
    Calculator, Fractals, Import Igor Experiment….
12. **Status line**.

(data-selector-workflow)=
## Typical workflow

The steps use a folder named `unified` holding copies of the three files in
`validationData/unified/`, each with a Unified Fit stored in it.

1. **Launch the GUI**
   ```bash
   pyirena-gui
   ```

2. **Select Data Folder**
   - Click "Select Data Folder" button
   - Navigate to your data directory (here `unified`)
   - Click "Select Folder"

3. **Choose File Type**
   - Use the dropdown to select file type
   - Default: HDF5 files

4. **Filter Files (Optional)**
   - Type in the filter box to search
   - Works like grep: the text is a full regular expression, matched anywhere
     in the file name, case-insensitively
   - Plain fragments such as `60C` still work; invalid patterns fall back to a
     substring match

5. **Select Files**
   - Single click to select one file
   - Ctrl/Cmd + click for multiple files
   - Double-click to plot immediately

6. **Create Graph**
   - Click "Create Graph" button
   - Graph opens in new window
   - Multiple files plotted together

   With **Unified Fit** checked under *Show in graph / reports*, the graph shows the stored
   fits with their residuals:

   ```{figure} img/data-selector/data-selector-02-graph.png
   :alt: Unified Fit Results graph for the three files: data and fitted model of each file with chi-squared in the legend, and rescaled residuals below.

   **Create Graph** with *Data* and *Unified Fit* checked.
   ```

7. **Tabulate the results.** Click **Tabulate Results** for a table of the stored parameters
   of every selected file; **Save as CSV…** writes it out.

   ```{figure} img/data-selector/data-selector-03-tabulate.png
   :alt: Tabulated results window: one row per file with UF_chi2, UF_background, UF_num_levels, L1_G and further columns, and a Save as CSV button.

   **Tabulate Results** for the three files.
   ```

(data-selector-controls)=
## Controls

(data-selector-file-controls)=
### Folder and files

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Select Data Folder** | Browse to a folder containing NXcanSAS/HDF5 data files | | The last folder is reopened next time |
| **Refresh** | Re-scan the current folder and update the file list | | Use after adding or removing files from the folder |
| **File Type** | HDF5 Files (.hdf, .h5, .hdf5), Text Files (.txt, .dat, .csv), All Supported Files | HDF5 Files | Text files are converted to cleaned NXcanSAS on first use |
| **Sort** | Order of the file list (same choices as the Data Explorer, see {ref}`data-explorer-file-select`) | Filename A→Z | |
| **Filter** | Grep-like filtering of file names, with full regular expression support | empty | |
| **Configure…** | Configure data loading options | | |
| **Manage Config…** | Inspect and remove tool sections from `pyirena_config.json` | | Use this to prevent unwanted tools from running during batch fitting |
| **Locate Logs…** | Open the pyirena log folder (`~/.pyirena/logs`) | | Attach the newest `gui.log` when reporting a problem |

(data-selector-results-controls)=
### Graphs, reports and tables

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Show in graph / reports** checkboxes | **Data**, **Unified Fit**, **Size Dist.**, **Simple Fits**, **WAXS Peaks**, **Carbon model**, **Modeling**, **3D saxsMorph**, **Fractals**: which stored results the graph, report and table use | Data | Files without the chosen results are skipped |
| **Create Graph** | Plot I(Q) for all selected files in a new graph window | | Double-clicking a file does the same |
| **Create Report** | Generate a Markdown report (.md) summarising the Unified Fit results for each selected HDF5 file | | Files without stored fit results are skipped |
| **Tabulate Results** | Build a table of fit results for selected files and display it | | Results included depend on the checked checkboxes; save as CSV |
| **Export to ASCII** | Write space-separated `.dat` files (Q, I, dI) into an `ascii_export` subfolder next to each selected HDF5 file | | Model curves are added per checked result type when enabled in **Configure…** |
| **Data Explorer** | Open the [Data Explorer](data-explorer.md) to browse, plot, and export data from HDF5 files | | Including export to Igor Pro h5xp format |

(data-selector-tool-buttons)=
### Tool buttons

| Control | Meaning | Default | Tips |
|---|---|---|---|
| ***Tool* (GUI)** | Open the tool's panel for the first selected file | | Hold **Shift** while clicking to reset the tool's window position (see below) |
| ***Tool* (script)** | Batch-fit all selected files with the tool using `pyirena_config.json`; results are saved into each file's NXcanSAS record | | See [Batch API](../guides/batch-api.md) for the config file |
| **Data Merge**, **Data Manipulation**, **Scattering Contrast Calculator**, **Fractals** | Open that tool | | |
| **Import Igor Experiment…** | Open an Igor Pro packed experiment (.pxp or .h5xp) and export each reduced USAXS/SAXS/WAXS sample into a NeXus (.h5) file | | See [Igor experiment import](../guides/igor-pxp-import.md) |

The same tools are in the **Models** and **Tools** menus.

(data-selector-shortcuts)=
### Keyboard shortcuts

Every window:

- **Cmd + W** (macOS) / **Ctrl + W** (Windows, Linux): Close this window
- **F1**: Open this tool's documentation — same as the **? Help** button
- **Cmd/Ctrl + C**: Copy selected table cells (**+ Shift** to include headers)

In the file list:

- **Double-click** on file: Plot immediately
- **Ctrl/Cmd + A**: Select all visible files
- **Ctrl/Cmd + Click**: Add/remove from selection

(data-selector-results)=
## Results and outputs

The Data Selector reads the results the tools stored in each file; **(script)** buttons write
new results into the selected files. Graphs, reports, tables and ASCII files are the outputs
it makes itself, as described in the tables above and in the next section.

(getting-graphs-and-tables-out)=
## Getting graphs and tables out of pyIrena

Every graph and every table behaves the same way, in every panel.

**Graphs — right-click anywhere on the plot:**

| Action | What it does |
|--------|--------------|
| Copy graph to clipboard | The graph as an image, ready to paste into PowerPoint, Word, email or an electronic notebook — the equivalent of Igor's Edit→Copy |
| Save graph as image… | PNG (default), JPEG or SVG — pick the format in the file dialog or just type the extension.  PNG is recommended: JPEG compression smears the thin lines and small text of a log-log plot |
| Save whole window as image… | All stacked panels at once (data + residuals + distribution) |
| Save curve data as CSV… | The plotted curves as text — one `X`/`Y` (and `dY` where errors exist) column pair per curve, for Excel, Origin, Matlab or pandas |
| Save as Igor Pro ITX… | Igor Pro Text waves with display, log-axis, colour and legend commands |

% SCREENSHOT data-selector-04-export-menu: right-click context menu on a graph (e.g. the Create
% Graph window) showing Copy graph to clipboard, Save graph as image…, Save whole window as image…,
% Save curve data as CSV…, Save as Igor Pro ITX… (a popup menu does not grab offscreen).

Exported curves are always in **linear (physical) units**, even when the plot
is logarithmic, and error bars are exported as an uncertainty column/wave
rather than as line segments — in Igor they arrive as error bars, not as an
extra trace.  ITX files also reproduce how each curve is drawn: data points
import as markers, model curves as lines.  Curves without a legend name — residuals, the
Simple Fits linearization, the Contrast (Δρ)² curve — are exported too, labelled
from the Y axis.

Save dialogs open in the folder you last exported to; before you have exported
anything they start in your data folder.  The choice persists across restarts.

**Fit results as text — two buttons on the panel:**

| Button | What it does |
|--------|--------------|
| Copy results | The current fit as Markdown on the clipboard: parameters with their uncertainties, fit quality, setup and data summary.  Paste into an e-mail, logbook or manuscript |
| Save report… | The same text as a `.md` file |
| **Ctrl-click** (⌘-click on macOS) Save report… | Saves the report **and** the graph: a PNG of the plots is written next to the `.md` with the same name, and embedded in it as a figure |

The embedded figure uses a *relative* link, so moving the `.md` and its `.png`
together — into a shared drive, a git repository, an Obsidian vault — keeps the
figure working.  The Markdown renders with its figure in GitHub, VS Code,
Jupyter, Obsidian and pandoc (so `pandoc report.md -o report.docx` gives you a
Word document with the graph in place).

Available on Unified Fit, Size Distribution, Simple Fits, Modeling and WAXS
Peak Fit.  The text is built by the same code as the Data Selector's
*Create Report* and the AI agents' `export_fit_report`, so a number is never
formatted one way in the panel and another way in the report.

**Tables — Ctrl+C, or right-click:**

| Action | What it does |
|--------|--------------|
| Ctrl+C | Copy selected cells as tab-separated text (pastes straight into Excel, Igor, Origin) |
| Ctrl+Shift+C | The same, with the column headers |
| Right-click → Copy Whole Table | Everything, headers included (also what Ctrl+C does with nothing selected) |
| Right-click → Save as CSV… | The table as a CSV file, where the panel supports it |
| Click a column header | Sort; numeric columns sort numerically and blanks sort last.  Not offered where row order carries meaning |

(data-selector-scripting)=
## From Python or AIDA

::::{tab-set}
:::{tab-item} Python
The folder listing and the I(Q) graph without the window. Run from the repository root:

```python
from pyirena import api

files = api.list_files("validationData/unified", sort="name_asc")
for f in files:
    print(f["name"], f["analyses"])               # analyses stored in each file

print(api.summarize_folder("validationData/unified")["n_files"])
api.plot_iq([f["path"] for f in files], output_path="unified_iq.png")
```

The **(script)** buttons call `pyirena.batch.fit_pyirena(data_file, config_file)` on each
selected file; see [Batch API](../guides/batch-api.md).
:::
:::{tab-item} AIDA / MCP
Ask in plain language, for example *"What is in this folder, and which files already have a
Unified Fit?"* The agent calls `pyirena_summarize_folder` and `pyirena_list_files`;
`pyirena_plot_iq` draws the curves and `pyirena_tabulate_parameter` builds the table. See the
[AI tools reference](../reference/ai-tools-reference.md).
:::
::::

(data-selector-method)=
## Method notes

### HDF5 Files (NXcanSAS)

The GUI uses `pyirena.io.hdf5.readGenericNXcanSAS()` to load:
- Q: Scattering vector [1/Å]
- Intensity: [cm⁻¹]
- Error: Uncertainties (optional)
- dQ: Q resolution (optional)

### Text Files (Future)

Support for tab/space delimited text files will be added:
- Column 1: Q
- Column 2: Intensity
- Column 3: Error (optional)

% CHECK(Jan): text files are already supported (File Type "Text Files", converted to NXcanSAS on
% first use, see guides/data-import-and-cleaning); this "Future" section is stale.

## References

None.

(data-selector-drag-and-drop)=
## Opening data by dragging it in

Drag files from Finder or Explorer onto a pyIrena window and they open — no
Select Folder, no hunting through a list.  It works on:

| Drop it on | What happens |
|------------|--------------|
| The Data Selector | Switches to the dropped file's folder and selects the dropped files |
| Any fitting panel (Unified Fit, Sizes, Simple Fits, Modeling, WAXS, SAXS Morph) | Loads the file, exactly as the panel's *Open…* button would |
| Data Explorer's file tree | Browses the dropped folder and selects the dropped files |
| Data Merge / Data Manipulation file lists | Points that dataset at the dropped folder and selects the files |

Dropping a **folder** works too: its data files are offered, one level deep, so
a measurement directory can be dragged in whole.  Files pyIrena cannot open are
ignored, and the drag is refused while it is still over the window if nothing in
it is openable — the cursor tells you before you let go.  Text files
(`.dat`/`.txt`/`.csv`) are converted to cleaned NXcanSAS on the way in, the same
as when opened from the dialog.  A fitting panel holds one dataset, so dropping
several files there loads the first.

(data-selector-window-geometry)=
## Windows remember where they were

Each pyIrena window reopens at the size and position you left it, and panels
with a draggable divider also remember how wide you made the controls.  This is
stored separately from the tool settings, in `window_geometry.json` beside the
pyIrena state file.

Screens change — you undock a laptop, unplug a second monitor, change
resolution — and a window remembered on a screen that no longer exists would
otherwise reopen somewhere you cannot see *or* drag back.  pyIrena checks the
saved position against the screens that exist now: if the window would be
unreachable, it is moved onto the nearest real screen (and shrunk only if it no
longer fits); if it cannot be placed sensibly at all, the tool simply opens at
its default size.

**To reset one tool:** hold **Shift** while clicking that tool's button in the
Data Selector.  That tool forgets where it was and opens centred at its default
size, with the control panel back at its original width — the same gesture as
Irena in Igor.  It works every time, not only the first time a tool is opened in
a session, and it resets the tool's graph window along with its panel.

**To reset everything:** hold **Shift** while launching pyIrena (`pyirena-gui`),
or start it with the environment variable `PYIRENA_RESET_WINDOWS=1`, which is the
reliable route on a system where the modifier key is intercepted:

```bash
PYIRENA_RESET_WINDOWS=1 pyirena-gui
```

(Shift-clicking a *file* still extends the selection as usual — the reset
gesture is on the tool buttons only.)

Minimised and full-screen windows are not saved — restoring a window that opens
minimised would be its own trap.

(data-selector-installation)=
## Installation and running

Install pyIrena with GUI support:

```bash
pip install pyirena[gui]
```

This will install:
- PySide6 (Qt6 for Python)
- matplotlib (for plotting)

**Method 1: Command Line Entry Point**

```bash
pyirena-gui
```

**Method 2: Python Module**

```bash
python -m pyirena.gui.launch
```

**Method 3: From Python Script**

% The block below is not run by the snippet check: it opens the window and waits for it to close.

```{code-block} python
from pyirena.gui.data_selector import main

main()
```

**Example test data.** Create a test data directory:

```bash
mkdir testData
# Copy your HDF5 files to testData/
```

The GUI will automatically detect files with supported extensions.

**Requirements**

- Python 3.10+
- PySide6
- matplotlib
- numpy
- h5py

(data-selector-troubleshooting)=
## Troubleshooting

### "No module named PySide6"

Install GUI dependencies:
```bash
pip install pyirena[gui]
```

### "Error loading file"

Check that:
- File is valid HDF5/NXcanSAS format
- File is not corrupted
- You have read permissions

### Graph window doesn't appear

- Check console for error messages
- Ensure matplotlib Qt backend is working
- Try: `pip install --upgrade matplotlib PySide6`

% CHECK(Jan): the Data Selector graphs are pyqtgraph, not matplotlib; "matplotlib Qt backend"
% advice may be stale.

(data-selector-development)=
## Development

### Adding New Features

The GUI is modular:
- `data_selector.py` - Main selector panel
- `__init__.py` - Qt backend detection
- `launch.py` - Entry point

To add new panels:
1. Create new module in `pyirena/gui/`
2. Import in `__init__.py`
3. Add launcher if needed

### Testing

```bash
# Run the GUI in development mode
cd pyirena/gui
python data_selector.py
```

% CHECK(Jan): stale paths — the panel is now pyirena/gui/data_selector/panel.py, Qt detection is
% pyirena/gui/_qt.py, and pyirena/gui/data_selector.py no longer exists (so the Testing command
% fails). developer/adding-features.md is the current checklist.

## Future Features

- [ ] Text file support
- [ ] Unified fit parameter panel
- [ ] Interactive fitting
- [ ] Export graphs (PNG, PDF, SVG)
- [ ] Data processing tools
- [ ] Batch processing
- [ ] Session saving/loading

% CHECK(Jan): most of these exist now (text import, Unified Fit panel, interactive fitting,
% graph export, data processing tools, batch processing, setup save/restore).

## Contributing

See [CONTRIBUTING.md](https://github.com/jilavsky/pyirena/blob/main/CONTRIBUTING.md) for guidelines.

## License

MIT License - See [LICENSE](https://github.com/jilavsky/pyirena/blob/main/LICENSE)
