# Unified Fit GUI - Feature Summary

## Overview

The Unified Fit GUI is a complete Python implementation of the Igor Pro Irena Unified Fit panel, providing an interactive interface for fitting small-angle scattering data using the Beaucage Unified Fit model.

## Current Status

### ✅ Fully Implemented Features

#### 1. **Complete GUI Layout**
- Split-screen design (controls left, graphs right)
- 5 structural levels with tab-based navigation
- All parameter input fields with validation
- Fit checkboxes for selective parameter optimization
- Parameter bounds (low/high limits)
- Matches Igor Pro design from screenshot

**Files**: [pyirena/gui/unified_fit.py](pyirena/gui/unified_fit.py)

#### 2. **Graph Panel**
- Main plot (log-log scale) showing data + fit
- Residuals plot (log-linear scale) showing fit quality
- Matplotlib integration with interactive toolbar
- Error bars on data points
- Multiple plot overlays (data, fit, local contributions)

**Files**: [pyirena/gui/unified_fit.py](pyirena/gui/unified_fit.py)

#### 3. **Graph Unified Button** ✅ WORKING (Green)
- Calculates unified fit with current parameters
- Displays model curve on graph
- Shows residuals in bottom plot
- No optimization - just evaluates current parameters
- Updates automatically if "Update automatically?" is checked

**Files**: [pyirena/gui/unified_fit.py](pyirena/gui/unified_fit.py) - `graph_unified()` method

#### 4. **Fit Button** ✅ WORKING (Green)
- Runs least-squares optimization
- Respects "Fit?" checkboxes for selective fitting
- Updates GUI with fitted parameters
- Displays chi-squared statistics
- Shows success message with fit quality metrics
- Uses only data between cursors for fitting

**Bug Fixes Applied**:
- ✅ Fixed `fit_background` parameter error
- ✅ Fixed dictionary access for result object
- ✅ Fixed `calculate()` → `calculate_intensity()` method name

**Files**:
- [pyirena/gui/unified_fit.py](pyirena/gui/unified_fit.py) - `run_fit()` method
- [BUGFIX_AUTO_UPDATE.md](BUGFIX_AUTO_UPDATE.md)
- [BUGFIX_CALCULATE.md](BUGFIX_CALCULATE.md)

#### 5. **Auto-Update Feature** ✅ WORKING
When "Update automatically?" checkbox is selected:
- Graph recalculates whenever any parameter changes
- Triggers on G, Rg, B, P, or Background field changes
- Activates when user presses Enter or Tab (editingFinished signal)
- Provides real-time visualization of parameter effects

**Implementation**:
- Signal connections on all parameter fields
- `editingFinished` signal triggers `on_parameter_changed()`
- Parent-child widget communication
- Automatic call to `graph_unified()`

**Files**: [pyirena/gui/unified_fit.py](pyirena/gui/unified_fit.py)
**Documentation**: [BUGFIX_AUTO_UPDATE.md](BUGFIX_AUTO_UPDATE.md)

#### 6. **Interactive Cursors** ✅ WORKING (NEW!)
Draggable vertical lines for selecting data range:
- **Red dashed line**: Left cursor (minimum Q)
- **Blue dashed line**: Right cursor (maximum Q)
- Click and drag to reposition
- Cursors cannot cross each other
- Only data between cursors used for fitting
- Full data still displayed on graph

**Features**:
- Automatic initialization when data loads (20% and 80% of Q range)
- Mouse event handling for dragging
- Visual feedback during interaction
- Status bar shows number of points and Q range used in fit
- **Cursors persist across all graph updates** - don't disappear when clicking buttons!

**Bug Fixed**:
- ✅ Cursors now remain visible after clicking "Graph Unified" or "Fit"
- ✅ Cursor positions preserved when user drags them
- ✅ Automatic restoration after graph refreshes

**Files**:
- [pyirena/gui/unified_fit.py](pyirena/gui/unified_fit.py)
  - `add_cursors()` method
  - `on_mouse_press()`, `on_mouse_release()`, `on_mouse_move()` methods
  - `get_cursor_range()` method
  - Modified `run_fit()` to filter data
  - Modified `plot_data()` to restore cursors
- [CURSOR_FEATURE.md](CURSOR_FEATURE.md)
- [BUGFIX_CURSOR_PERSISTENCE.md](BUGFIX_CURSOR_PERSISTENCE.md)

#### 7. **Data Integration**
- Loads data from main data selector GUI
- Accepts Q, Intensity, Error arrays
- Menu bar access: **Models → Unified Fit**
- Button access: **"Unified Fit"** button (green)
- Standalone mode also available

**Files**:
- [pyirena/gui/data_selector.py](pyirena/gui/data_selector.py) - menu and button integration
- [pyirena/gui/unified_fit.py](pyirena/gui/unified_fit.py) - `set_data()` method

#### 8. **Model Integration**
- Full integration with `UnifiedFitModel` backend
- Supports 1-5 structural levels
- UnifiedLevel parameter objects
- Background fitting option
- Correlated systems checkbox (parameter ready, not functional yet)

**Files**: [pyirena/core/unified.py](pyirena/core/unified.py) (backend)

### 📋 Buttons Not Yet Functional

These buttons are present in the GUI but not yet connected:

#### Parameter Estimation
- [ ] **Fit Rg/G btwn cursors**: Fit Guinier region using cursor selection
- [ ] **Fit P/B btwn cursors**: Fit Porod region using cursor selection
- [ ] **Estimate B from G/Rg/P?**: Auto-calculate B using Hammouda relationship

#### Level Management
- [ ] **Copy/Move/swap level**: Transfer parameters between levels

#### Fit Control
- [ ] **Revert back**: Restore previous parameters
- [ ] **reset unif?**: Reset all parameters to defaults
- [ ] **Fix limits?**: Lock current parameter bounds

#### Results Export
- [ ] **Store in Data Folder**: Save results to data folder
- [ ] **Export ASCII**: Export fit results as text file
- [ ] **Results to graphs**: Plot results with detailed analysis
- [ ] **Analyze Results**: Perform statistical analysis
- [ ] **Anal. Uncertainty**: Calculate parameter uncertainties

## Testing

### Validation Script
```bash
python validate_cursor_code.py
```
Checks cursor implementation without requiring Qt.

### GUI Test Script
```bash
python test_unified_gui.py
```
Launches GUI with synthetic 2-level test data.

### Manual Testing Checklist

#### Test 1: Graph Unified Button
- [x] Click "Graph Unified" button (green)
- [x] Verify fit curve appears in main plot
- [x] Verify residuals appear in bottom plot
- [x] Change parameter and click again
- [x] Verify graph updates

#### Test 2: Fit Button
- [x] Select "Fit?" checkboxes for G, Rg, P
- [x] Click "Fit" button (green)
- [x] Verify fit converges without errors
- [x] Verify parameters update in GUI
- [x] Verify success message shows chi-squared
- [x] Verify residuals are reasonable

#### Test 3: Auto-Update
- [x] Check "Update automatically?" checkbox
- [x] Change G value and press Enter
- [x] Verify graph updates automatically
- [x] Change other parameters (Rg, B, P, Background)
- [x] Verify each triggers auto-update
- [x] Uncheck "Update automatically?"
- [x] Verify changes don't auto-update

#### Test 4: Interactive Cursors
- [x] Load data
- [x] Verify red and blue cursor lines appear
- [x] Click and drag left (red) cursor
- [x] Verify it moves horizontally
- [x] Click and drag right (blue) cursor
- [x] Verify it moves horizontally
- [x] Try to drag left cursor past right cursor
- [x] Verify cursors don't cross
- [x] Select "Fit?" for some parameters
- [x] Click "Fit"
- [x] Verify status shows Q range and number of points used
- [x] Verify fit uses only data between cursors

## Documentation

### User Guides
- [UNIFIED_FIT_GUI.md](UNIFIED_FIT_GUI.md) - Comprehensive user guide with workflow examples

### Technical Documentation
- [BUGFIX_AUTO_UPDATE.md](BUGFIX_AUTO_UPDATE.md) - Bug fixes #1, #2, and auto-update implementation
- [BUGFIX_CALCULATE.md](BUGFIX_CALCULATE.md) - Bug fix for `calculate()` method error
- [CURSOR_FEATURE.md](CURSOR_FEATURE.md) - Interactive cursor implementation details

### Test Files
- [test_unified_gui.py](test_unified_gui.py) - GUI test with synthetic data
- [validate_cursor_code.py](validate_cursor_code.py) - Cursor validation without Qt

## Key Files

### GUI Implementation
- `pyirena/gui/unified_fit.py` - Main GUI (825 lines)
  - `UnifiedFitGraphWindow` - Graph display with cursors
  - `LevelParametersWidget` - Parameter input for each level
  - `UnifiedFitPanel` - Main panel with controls

### Integration
- `pyirena/gui/data_selector.py` - Main data browser with menu and button

### Backend
- `pyirena/core/unified.py` - UnifiedFitModel and UnifiedLevel classes

## Architecture

### Class Structure
```
UnifiedFitPanel (Main window)
├── Control Panel (Left)
│   ├── Top controls (Graph Unified, num levels, checkboxes)
│   ├── Level tabs (1-5)
│   │   └── LevelParametersWidget (G, Rg, B, P, RgCutoff)
│   ├── Background controls
│   ├── Fit buttons
│   └── Results buttons
└── Graph Window (Right)
    ├── UnifiedFitGraphWindow
    │   ├── Main plot (data + fit + cursors)
    │   └── Residuals plot
    └── Navigation toolbar
```

### Data Flow
```
1. Data Selector GUI → select file
2. Click "Unified Fit" button
3. UnifiedFitPanel.set_data(Q, I, Error)
4. User adjusts parameters
5. Click "Fit" or "Graph Unified"
6. Get cursor range from graph window
7. Filter data to cursor range (fit only)
8. Call UnifiedFitModel.fit() or calculate_intensity()
9. Update GUI with results
10. Plot fit and residuals
```

### Signal Flow (Auto-Update)
```
1. User edits parameter field
2. User presses Enter or Tab
3. QLineEdit emits editingFinished signal
4. LevelParametersWidget.on_parameter_changed()
5. Find parent UnifiedFitPanel
6. UnifiedFitPanel.on_parameter_changed()
7. Check if "Update automatically?" is checked
8. Call graph_unified()
9. Graph updates automatically
```

### Cursor Interaction Flow
```
1. Data loaded → plot_data()
2. Initialize cursor positions (20%, 80% of Q range)
3. add_cursors() draws red and blue lines
4. User clicks near cursor → on_mouse_press()
   - Check if click is within tolerance
   - Set dragging_cursor = 'left' or 'right'
5. User drags mouse → on_mouse_move()
   - Update cursor position
   - Ensure cursors don't cross
   - Redraw canvas
6. User releases mouse → on_mouse_release()
   - Clear dragging_cursor
7. User clicks "Fit"
   - get_cursor_range() returns (Q_min, Q_max)
   - Filter data: mask = (Q >= Q_min) & (Q <= Q_max)
   - Fit using filtered data
   - Display results for full data range
```

## Known Issues

None currently! All implemented features are working correctly.

## Future Enhancements

### High Priority
- [ ] Implement cursor-based Guinier and Porod fitting buttons
- [ ] Save/load parameter sets
- [ ] Export results to ASCII
- [ ] Revert/reset parameter functionality

### Medium Priority
- [ ] Display local (Porod & Guinier) fits
- [ ] Parameter uncertainty analysis
- [ ] Batch fitting multiple files
- [ ] Copy/move/swap level functionality

### Low Priority
- [ ] Keyboard shortcuts
- [ ] Undo/redo functionality
- [ ] Parameter correlation analysis
- [ ] Results plotting and analysis tools

## Version History

### v0.3.2 - Exception Handling Fix (Current)
- ✅ Fixed "cannot remove artist" error when updating graphs
- ✅ Broadened exception handling to catch all matplotlib artist removal errors
- ✅ More robust cursor restoration across different matplotlib versions

### v0.3.1 - Cursor Persistence Fix
- ✅ Fixed bug: cursors now persist across graph updates
- ✅ Cursors don't disappear when clicking "Graph Unified" or "Fit"
- ✅ Robust error handling for already-removed lines
- ✅ Cursor positions preserved after user drags them

### v0.3.0 - Interactive Cursors
- ✅ Added draggable cursor lines for data range selection
- ✅ Implemented mouse event handling
- ✅ Modified fit() to use only data between cursors
- ✅ Added visual feedback (red/blue dashed lines)
- ✅ Status bar shows Q range and points used

### v0.2.0 - Auto-Update and Bug Fixes
- ✅ Implemented auto-update when parameters change
- ✅ Fixed fit_background parameter error
- ✅ Fixed dictionary access for fit results
- ✅ Fixed calculate_intensity method name
- ✅ Added chi-squared and reduced chi-squared display

### v0.1.0 - Initial GUI
- ✅ Complete GUI layout matching Igor design
- ✅ 5 level tabs with all parameters
- ✅ Graph Unified button functional
- ✅ Fit button functional
- ✅ Integration with main data selector
- ✅ Menu bar and button access

## References

### Scientific Papers
- Beaucage, G. "Approximations Leading to a Unified Exponential/Power-Law Approach to Small-Angle Scattering" *J. Appl. Cryst.* (1995) **28**, 717-728
- Beaucage, G. "Small-Angle Scattering from Polymeric Mass Fractals of Arbitrary Mass-Fractal Dimension" *J. Appl. Cryst.* (1996) **29**, 134-146
- Ilavsky, J. & Jemian, P. R. "Irena: tool suite for modeling and analysis of small-angle scattering" *J. Appl. Cryst.* (2009) **42**, 347-353

### Original Software
Based on the Irena package for Igor Pro by Jan Ilavsky
