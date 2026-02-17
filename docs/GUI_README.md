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
   - Text files (.txt, .dat)
   - All supported files
3. **File List** - View and select files from the folder
4. **Text Filter** - Grep-like filtering of file names
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
   - Works like grep (case-insensitive)

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

- Python 3.8+
- PySide6 or PyQt6
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
