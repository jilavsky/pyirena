# Testing the pyIrena GUI

## Quick Installation & Test

### 1. Install with GUI support

```bash
# From the project directory
pip install -e ".[gui]"
```

**What this installs:**
- pyirena package
- PySide6 (Qt6 for Python)
- matplotlib (for plotting)
- numpy, scipy, h5py (dependencies)

### 2. Generate test data

```bash
python create_test_data.py
```

**This creates 6 HDF5 test files in `testData/`:**
- `complexUnified.h5` - Your requested test file ⭐
- `spheres_Rg50.h5` - Small particles
- `spheres_Rg100.h5` - Large particles
- `fractal_aggregate.h5` - Mass fractal
- `surface_fractal.h5` - Surface fractal
- `hierarchical_structure.h5` - Two-level

### 2b. Generate text (.dat) files (Optional)

```bash
python convert_to_dat.py
```

**This creates corresponding .dat files from the HDF5 files:**
- `complexUnified.dat`
- `spheres_Rg50.dat`
- `spheres_Rg100.dat`
- `fractal_aggregate.dat`
- `surface_fractal.dat`
- `hierarchical_structure.dat`

### 3. Launch the GUI

```bash
pyirena-gui
```

**Alternative methods:**
```bash
# Using Python module
python -m pyirena.gui.launch

# Direct execution
python pyirena/gui/data_selector.py
```

### 4. Use the GUI

**Step-by-step test:**

1. **Click "Select Data Folder"**
   - Navigate to the `testData` folder
   - Click "Select Folder"
   - ✅ Files should appear in the list

2. **Verify file list shows (HDF5 files by default):**
   ```
   ☐ complexUnified.h5
   ☐ fractal_aggregate.h5
   ☐ hierarchical_structure.h5
   ☐ ProperNxcanSASNexus.h5
   ☐ spheres_Rg100.h5
   ☐ spheres_Rg50.h5
   ☐ surface_fractal.h5
   ```

3. **Test filtering:**
   - Type "complex" in filter box
   - ✅ Only `complexUnified.h5` should show
   - Clear filter (delete text)
   - ✅ All files should reappear

4. **Plot a single file:**
   - Double-click `complexUnified.h5`
   - ✅ Graph window should open
   - ✅ See log-log plot of Intensity vs Q

5. **Plot multiple files:**
   - Hold Ctrl/Cmd and click multiple files
   - Click "Create Graph" button
   - ✅ All selected files plotted together
   - ✅ Different colors, with legend

## What the Graph Should Show

### complexUnified.h5

Expected features in the plot:
- **High Q region** (Q > 0.1): Power law decay
- **Mid Q region** (0.01 < Q < 0.1): Multiple features
- **Low Q region** (Q < 0.01): Guinier plateau
- **Three distinct levels** visible as changes in slope

**Note:** The test data files use a simple HDF5 structure. The GUI automatically detects this and uses a simple reader as a fallback when full NXcanSAS structure is not found.

### Comparison Plot

Try plotting all files together:
1. Select all files (Ctrl/Cmd + A)
2. Click "Create Graph"
3. You should see:
   - Different Rg values (different Guinier regions)
   - Different power laws (different slopes)
   - Clear separation between samples

## Verification Checklist

### GUI Layout ✅
- [ ] Title "pyIrena" at top
- [ ] "Select Data Folder" button visible
- [ ] "Refresh" button next to folder selection
- [ ] File type dropdown (HDF5/Text/All)
- [ ] Filter input box below file type
- [ ] File listbox (wide, 15+ rows visible)
- [ ] "Create Graph" button to the right (enabled only when files selected)
- [ ] Status bar at bottom
- [ ] Window is ~2x width of listbox

### Functionality ✅
- [ ] Folder selection works
- [ ] Folder selection remembers last folder
- [ ] Refresh button reloads folder contents
- [ ] Files load and display
- [ ] File type filter works (try switching types)
- [ ] Text filter works (type partial names)
- [ ] Single selection works (click)
- [ ] Multiple selection works (Ctrl/Cmd + click)
- [ ] "Create Graph" button enables when files selected
- [ ] "Create Graph" button disables when no files selected
- [ ] Double-click plots immediately
- [ ] "Create Graph" button works
- [ ] Graph opens in new window
- [ ] Graph shows log-log scale
- [ ] Multiple files overlay correctly
- [ ] Legend shows filenames
- [ ] Status bar updates
- [ ] Both NXcanSAS and simple HDF5 files load correctly

### Data Display ✅
- [ ] X-axis labeled "Q (Å⁻¹)"
- [ ] Y-axis labeled "Intensity (cm⁻¹)"
- [ ] Both axes are logarithmic
- [ ] Grid visible
- [ ] Data points connected with lines
- [ ] Markers visible on lines
- [ ] Colors distinguish different files

## Troubleshooting

### Problem: "No module named PySide6"

**Solution:**
```bash
pip install PySide6
# or
pip install -e ".[gui]"
```

### Problem: GUI doesn't start

**Check 1:** Python version
```bash
python --version  # Must be 3.8+
```

**Check 2:** Import test
```python
python -c "from pyirena.gui.data_selector import main; print('OK')"
```

**Check 3:** Dependencies
```bash
pip install -e ".[gui]" --force-reinstall
```

### Problem: Graph window is blank

**Solution 1:** Update matplotlib
```bash
pip install --upgrade matplotlib
```

**Solution 2:** Check backend
```python
import matplotlib
print(matplotlib.get_backend())  # Should be Qt5Agg or similar
```

### Problem: Can't load HDF5 files

**Check file structure:**
```python
import h5py
with h5py.File('testData/complexUnified.h5', 'r') as f:
    print(f.keys())
    # Should show: ['entry1']
```

**Regenerate test data:**
```bash
rm testData/*.h5
python create_test_data.py
```

### Problem: "Create Graph" button is disabled

**Cause:** No files selected

**Solution:** Click at least one file in the list

## Advanced Testing

### Test with Real Data

If you have your own HDF5 files:

1. Copy them to a new folder
2. Select that folder in the GUI
3. Verify they load correctly
4. Try plotting

### Performance Test

Create many files:

```python
# Modify create_test_data.py to generate 100 files
for i in range(100):
    create_nxcansas_file(f'testData/test_{i:03d}.h5', q, I, error)
```

Test:
- Loading time
- Filter responsiveness
- Plotting speed

### File Type Test

1. Switch to "Text Files (.txt, .dat)" in dropdown
2. Verify .dat files appear:
   ```
   ☐ complexUnified.dat
   ☐ fractal_aggregate.dat
   ☐ hierarchical_structure.dat
   ☐ spheres_Rg100.dat
   ☐ spheres_Rg50.dat
   ☐ surface_fractal.dat
   ```
3. Verify no .h5 files show
4. Try plotting a .dat file (should work identically to .h5)
5. Switch to "All Supported Files"
6. Verify both .h5 and .dat files appear together

## Expected Performance

| Operation | Expected Time |
|-----------|--------------|
| Launch GUI | < 2 seconds |
| Select folder | < 0.5 seconds |
| Load 10 files | < 0.1 seconds |
| Filter files | Instant |
| Plot 1 file | < 1 second |
| Plot 5 files | < 2 seconds |

## Success Criteria

The GUI test is successful if:

1. ✅ GUI launches without errors
2. ✅ Folder selection works
3. ✅ Test files load correctly
4. ✅ Filter works as expected
5. ✅ Plots display correctly
6. ✅ Multiple files can be overlaid
7. ✅ Graph window shows proper labels and scales
8. ✅ complexUnified.h5 plots showing expected features

## Next Steps After Testing

Once GUI works:

1. **Try your own data** - Load real experimental files
2. **Customize appearance** - Modify styles in data_selector.py
3. **Add features** - Extend the GUI with new panels
4. **Report issues** - Open GitHub issues for bugs
5. **Request features** - Suggest improvements

## Getting Help

If you encounter issues:

1. Check console output for error messages
2. Review [GUI_QUICKSTART.md](GUI_QUICKSTART.md)
3. Check [pyirena/gui/README.md](pyirena/gui/README.md)
4. Open issue: https://github.com/jilavsky/pyirena/issues
5. Email: ilavsky@aps.anl.gov

## Reporting Bugs

Include:
- Python version (`python --version`)
- OS (Windows/Mac/Linux)
- Installation method
- Error messages (full traceback)
- Steps to reproduce

---

**Happy Testing! 🧪**
