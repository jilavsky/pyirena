# pyIrena GUI Quick Start Guide

## Installation

### Step 1: Install pyIrena with GUI support

```bash
# Install from source with GUI dependencies
pip install -e ".[gui]"
```

This installs:
- pyirena (core package)
- PySide6 (Qt6 for Python)
- matplotlib (for plotting)
- All required dependencies (numpy, scipy, h5py)

### Step 2: Verify Installation

```bash
# Check that pyirena-gui command is available
pyirena-gui --help 2>&1 || python -m pyirena.gui.launch
```

## Running the GUI

### Method 1: Command-line launcher

```bash
pyirena-gui
```

### Method 2: Python module

```bash
python -m pyirena.gui.launch
```

### Method 3: From Python script

```python
from pyirena.gui.data_selector import main
main()
```

## Using the Data Selector

### Interface Layout

```
┌─────────────────────────────────────────────────────────┐
│                       pyIrena                           │
├─────────────────────────────────────────────────────────┤
│  [Select Data Folder]  /path/to/data                   │
│                                                         │
│  File Type: [HDF5 Files (.hdf, .h5, .hdf5)  ▼]        │
│                                                         │
│  Filter: [Enter text to filter files...        ]       │
│                                                         │
│  ┌──────────────────────────────┐                      │
│  │ ☑ complexUnified.h5         │  [Create Graph]      │
│  │ ☐ fractal_aggregate.h5      │                      │
│  │ ☐ hierarchical_structure.h5 │                      │
│  │ ☑ spheres_Rg50.h5           │                      │
│  │ ☐ spheres_Rg100.h5          │                      │
│  │ ☐ surface_fractal.h5        │                      │
│  │                              │                      │
│  └──────────────────────────────┘                      │
│                                                         │
│  Status: Showing 6 of 6 files                          │
└─────────────────────────────────────────────────────────┘
```

### Step-by-Step Workflow

#### 1. Select Data Folder

- Click **"Select Data Folder"** button
- Navigate to your data directory (e.g., `testData/`)
- Click "Select Folder" or "Choose"

The folder path appears next to the button, and files are automatically loaded.

#### 2. Choose File Type (Optional)

Use the dropdown to filter by file type:
- **HDF5 Files** (.hdf, .h5, .hdf5) - Default
- **Text Files** (.txt, .dat) - For ASCII data
- **All Supported Files** - Show all

#### 3. Filter Files (Optional)

- Type in the **Filter** box to search for files
- Works like grep (case-insensitive substring match)
- Examples:
  - `sphere` → shows only files containing "sphere"
  - `Rg50` → shows spheres_Rg50.h5
  - `fractal` → shows fractal files

#### 4. Select Files to Plot

**Single selection:**
- Click on a file name

**Multiple selection:**
- Hold Ctrl (Windows/Linux) or Cmd (Mac)
- Click multiple files
- Or: Click first file, hold Shift, click last file (range selection)

**Quick plot:**
- Double-click any file to plot immediately

#### 5. Create Graph

- Click **"Create Graph"** button
- A new window opens with the plot
- Multiple selected files are overlaid on the same graph

### Graph Window Features

The graph displays:
- **Log-log scale** for Q vs Intensity
- **Multiple datasets** with different colors
- **Legend** showing file names
- **Grid** for easy reading
- **Axis labels**: Q (Å⁻¹) and Intensity (cm⁻¹)

## Test Data

The package includes synthetic test data in `testData/`:

| File | Description | Features |
|------|-------------|----------|
| `spheres_Rg50.h5` | Small spherical particles | Rg ≈ 50 Å, Porod tail |
| `spheres_Rg100.h5` | Large spherical particles | Rg ≈ 100 Å |
| `fractal_aggregate.h5` | Mass fractal | Power law slope ≈ 2.5 |
| `surface_fractal.h5` | Surface fractal | Power law slope ≈ 3.5 |
| `hierarchical_structure.h5` | Multi-level | Two size scales |
| `complexUnified.h5` | Complex 3-level | Primary + aggregates + clusters |

### Regenerate Test Data

```bash
python create_test_data.py
```

## Features

### Current Features ✅

- [x] Folder browsing
- [x] File type filtering (HDF5/Text)
- [x] File list with scroll
- [x] Text-based filtering (grep-like)
- [x] Multiple file selection
- [x] Double-click to plot
- [x] Graph in separate window
- [x] Log-log plotting
- [x] Multiple datasets overlay
- [x] NXcanSAS HDF5 support

### Planned Features 🚧

- [ ] Text file support (.txt, .dat)
- [ ] Save/export graphs (PNG, PDF)
- [ ] Zoom and pan tools
- [ ] Data table view
- [ ] Unified fit parameter panel
- [ ] Interactive fitting
- [ ] Residuals plot
- [ ] Batch processing

## Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| **Double-click** | Plot selected file immediately |
| **Ctrl/Cmd + A** | Select all visible files |
| **Ctrl/Cmd + Click** | Add/remove file from selection |
| **Shift + Click** | Select range of files |

## Supported File Formats

### HDF5 Files (NXcanSAS)

Expected structure:
```
/entry1/data1/
    Q           - Scattering vector [1/Å]
    I           - Intensity [cm⁻¹]
    Idev        - Error/uncertainty (optional)
    Qdev        - Q resolution (optional)
```

The GUI uses `pyirena.io.hdf5.readGenericNXcanSAS()` for reading.

### Text Files (Coming Soon)

Expected format:
```
# Comment lines start with #
Q           Intensity   Error
0.001       1000.0      10.0
0.002       950.0       9.5
0.003       900.0       9.0
...
```

## Troubleshooting

### Issue: "No module named PySide6"

**Solution:**
```bash
pip install pyirena[gui]
# or
pip install PySide6 matplotlib
```

### Issue: GUI doesn't start

**Solution 1:** Check Python version
```bash
python --version  # Should be 3.8+
```

**Solution 2:** Try running directly
```bash
python -m pyirena.gui.data_selector
```

**Solution 3:** Check error messages
```bash
python -c "from pyirena.gui.data_selector import main; main()"
```

### Issue: "Error loading file"

**Possible causes:**
- File is not valid NXcanSAS format
- File is corrupted
- Missing read permissions

**Solution:** Check file with h5py
```python
import h5py
with h5py.File('testData/complexUnified.h5', 'r') as f:
    print(list(f.keys()))
    print(list(f['entry1']['data1'].keys()))
```

### Issue: Graph window is blank

**Solution 1:** Update matplotlib
```bash
pip install --upgrade matplotlib
```

**Solution 2:** Try different backend
```python
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg'
```

### Issue: Can't select multiple files

**Solution:** Use Ctrl/Cmd key while clicking

## Tips

1. **Performance:** For many files (>100), use the filter to narrow the list
2. **Quick plotting:** Double-click is fastest for single files
3. **Comparison:** Select multiple files to overlay and compare
4. **File organization:** Name files systematically for easy filtering
5. **Test first:** Use the provided test data to learn the interface

## Advanced Usage

### Creating Custom Test Data

```python
import numpy as np
import h5py

# Generate synthetic data
q = np.logspace(-3, 0, 100)
I = 1000 * q**(-4) + 0.01

# Save as HDF5
with h5py.File('my_data.h5', 'w') as f:
    entry = f.create_group('entry1')
    data = entry.create_group('data1')
    data.create_dataset('Q', data=q)
    data.create_dataset('I', data=I)
```

### Programmatic Access

```python
from pyirena.gui.data_selector import DataSelectorPanel
from PySide6.QtWidgets import QApplication

app = QApplication([])
window = DataSelectorPanel()
window.show()
app.exec()
```

## Next Steps

After getting familiar with the data selector:

1. **Analyze your data** - Load your experimental files
2. **Try unified fitting** - Use `pyirena.core.unified`
3. **Explore examples** - See `pyirena/examples/`
4. **Read documentation** - Check the full user guide

## Getting Help

- **GitHub Issues:** https://github.com/jilavsky/pyirena/issues
- **Documentation:** https://github.com/jilavsky/pyirena#readme
- **Email:** ilavsky@aps.anl.gov

## Screenshots

(Screenshots will be added here)

---

**Enjoy using pyIrena GUI!** 🚀
