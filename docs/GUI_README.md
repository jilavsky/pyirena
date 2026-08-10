# pyIrena GUI Module

Graphical user interface for pyIrena data analysis.

## Installation

Install pyIrena with GUI support:

```bash
pip install pyirena[gui]
```

This will install:
- PySide6 (Qt6 for Python)
- matplotlib (for plotting)

## Running the GUI

### Method 1: Command Line Entry Point

```bash
pyirena-gui
```

### Method 2: Python Module

```bash
python -m pyirena.gui.launch
```

### Method 3: From Python Script

```python
from pyirena.gui.data_selector import main

main()
```

## Features

### Data Selector Panel

The main GUI provides:

1. **Folder Selection** - Browse and select data directories
2. **File Type Filtering** - Choose between:
   - HDF5 files (.hdf, .h5, .hdf5)
   - Text files (.txt, .dat, .csv) — automatically converted to cleaned NXcanSAS
     HDF5 on first use (see [Text File Import and Cleaning](data_import_and_cleaning.md))
   - All supported files
3. **File List** - View and select files from the folder
4. **Text Filter** - Grep-like filtering of file names, with full regular
   expression support
5. **Multi-selection** - Select multiple files (Ctrl/Cmd + click)
6. **Graphing** - Plot selected files with:
   - Log-log scale
   - Multiple datasets overlaid
   - Legend and labels

## Usage

### Basic Workflow

1. **Launch the GUI**
   ```bash
   pyirena-gui
   ```

2. **Select Data Folder**
   - Click "Select Data Folder" button
   - Navigate to your data directory
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

### Keyboard Shortcuts

- **Double-click** on file: Plot immediately
- **Ctrl/Cmd + A**: Select all visible files
- **Ctrl/Cmd + Click**: Add/remove from selection

### Getting graphs and tables out of pyIrena

Every graph and every table behaves the same way, in every panel.

**Graphs — right-click anywhere on the plot:**

| Action | What it does |
|--------|--------------|
| Copy graph to clipboard | The graph as an image, ready to paste into PowerPoint, Word, email or an electronic notebook — the equivalent of Igor's Edit→Copy |
| Save graph as image… | PNG (default), JPEG or SVG — pick the format in the file dialog or just type the extension.  PNG is recommended: JPEG compression smears the thin lines and small text of a log-log plot |
| Save whole window as image… | All stacked panels at once (data + residuals + distribution) |
| Save curve data as CSV… | The plotted curves as text — one `X`/`Y` (and `dY` where errors exist) column pair per curve, for Excel, Origin, Matlab or pandas |
| Save as Igor Pro ITX… | Igor Pro Text waves with display, log-axis, colour and legend commands |

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

### Opening data by dragging it in

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

### Windows remember where they were

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

## Data Format Support

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

## Example Test Data

Create a test data directory:

```bash
mkdir testData
# Copy your HDF5 files to testData/
```

The GUI will automatically detect files with supported extensions.

## Requirements

- Python 3.10+
- PySide6
- matplotlib
- numpy
- h5py

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

## Future Features

- [ ] Text file support
- [ ] Unified fit parameter panel
- [ ] Interactive fitting
- [ ] Export graphs (PNG, PDF, SVG)
- [ ] Data processing tools
- [ ] Batch processing
- [ ] Session saving/loading

## Screenshots

(Add screenshots here when available)

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for guidelines.

## License

MIT License - See [LICENSE](../../LICENSE)
