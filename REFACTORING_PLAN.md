# PyIrena Refactoring Plan

## Overview

Major refactoring to improve performance and add state management capabilities.

## Goals

1. **Replace matplotlib with pyqtgraph** for better performance
2. **Add state management** to persist settings across sessions
3. **Add reset to defaults** functionality
4. **Enable programmatic API** for other tools to use Unified Fit parameters

## Implementation Phases

### Phase 1: State Management System ✅ COMPLETE

**Status**: Implemented

**Files Created**:
- `pyirena/state/__init__.py`
- `pyirena/state/state_manager.py`

**Features**:
- Hierarchical JSON state structure
- Human-readable format
- Save/load state across sessions
- Export/import tool-specific states
- Merge with defaults for backward compatibility
- State file location: `~/.pyirena/state.json`

**State Structure**:
```json
{
  "version": "1.0",
  "unified_fit": {
    "num_levels": 1,
    "levels": [
      {
        "level": 1,
        "G": {"value": 100, "fit": false, "low_limit": null, "high_limit": null},
        "Rg": {"value": 100, "fit": false, "low_limit": null, "high_limit": null},
        "B": {"value": 0.01, "fit": false, "low_limit": null, "high_limit": null},
        "P": {"value": 4, "fit": false, "low_limit": null, "high_limit": null},
        "RgCutoff": 0,
        "correlated": false,
        "estimate_B": false
      }
      // ... levels 2-5
    ],
    "background": {"value": 1e-6, "fit": false},
    "cursor_left": null,
    "cursor_right": null,
    "update_auto": false,
    "display_local": false,
    "no_limits": false,
    "skip_fit_check": false
  }
}
```

### Phase 2: Pyqtgraph Migration (IN PROGRESS)

**Goal**: Replace matplotlib with pyqtgraph for all plotting

**Benefits**:
- 10-100x faster rendering
- Better interactivity
- Smoother cursor dragging
- Hardware-accelerated rendering
- Lower CPU usage

**Files to Modify**:
1. `pyirena/gui/unified_fit.py` - Main unified fit GUI
2. `pyirena/gui/data_selector.py` - Data browser GUI
3. `pyirena/plotting/unified_plots.py` - Plotting utilities

**Pyqtgraph Equivalents**:
- `matplotlib.figure.Figure` → `pyqtgraph.PlotWidget`
- `ax.errorbar()` → `plot.plot()` with error bars using `ErrorBarItem`
- `ax.axvline()` → `InfiniteLine` (much better for cursors!)
- `ax.set_xscale('log')` → `plot.setLogMode(x=True, y=True)`
- Mouse events → `sigMouseClicked`, `sigMouseMoved`

**Cursor Implementation**:
- Use `pg.InfiniteLine` with `movable=True`
- Built-in drag support (no manual mouse handling!)
- Signal `sigPositionChanged` for tracking cursor moves
- Much simpler than matplotlib implementation

### Phase 3: Integrate State Management with GUI

**Goal**: Connect state manager to Unified Fit GUI

**Tasks**:
1. Load state on GUI startup
2. Save state when parameters change
3. Add "Reset to Defaults" button
4. Auto-save on close

**Implementation**:
```python
class UnifiedFitPanel:
    def __init__(self):
        self.state_manager = StateManager()
        self.load_state()

    def load_state(self):
        # Load from state manager
        state = self.state_manager.get('unified_fit')
        self.apply_state(state)

    def save_state(self):
        # Save current GUI state
        state = self.get_current_state()
        self.state_manager.update('unified_fit', state)
        self.state_manager.save()

    def reset_to_defaults(self):
        # Reset button clicked
        self.state_manager.reset('unified_fit')
        self.load_state()
```

### Phase 4: Programmatic API

**Goal**: Allow other tools to read/write Unified Fit parameters

**API Design**:
```python
from pyirena.state import StateManager

# Read parameters
sm = StateManager()
unified_params = sm.get('unified_fit')
num_levels = unified_params['num_levels']
level1_G = unified_params['levels'][0]['G']['value']

# Write parameters
sm.set('unified_fit', 'num_levels', 2)
sm.save()

# Export/import presets
sm.export_tool_state('unified_fit', 'guinier_preset.json')
sm.import_tool_state('unified_fit', 'porod_preset.json')
```

**Use Cases**:
1. **Batch processing**: Load fit parameters, apply to multiple datasets
2. **Preset creation**: Save common fitting scenarios
3. **Automated workflows**: Other tools can configure Unified Fit programmatically
4. **Testing**: Create known-good parameter sets

### Phase 5: Testing & Documentation

**Tasks**:
1. Update test scripts to use pyqtgraph
2. Test state save/load
3. Test reset functionality
4. Create user documentation
5. Create developer API documentation

## Current Status

- ✅ Phase 1: Complete - State management system implemented
- 🚧 Phase 2: In progress - Need to implement pyqtgraph migration
- ⏳ Phase 3: Pending - State integration
- ⏳ Phase 4: Pending - Programmatic API
- ⏳ Phase 5: Pending - Testing & docs

## Migration Notes

### Matplotlib → Pyqtgraph Conversion

**Data Plot**:
```python
# OLD (matplotlib)
self.ax_main.errorbar(q, intensity, yerr=error, fmt='o')

# NEW (pyqtgraph)
self.plot_widget.plot(q, intensity, pen=None, symbol='o')
error_bars = pg.ErrorBarItem(x=q, y=intensity, height=error)
self.plot_widget.addItem(error_bars)
```

**Fit Curve**:
```python
# OLD (matplotlib)
self.ax_main.plot(q, fit, '-', linewidth=2)

# NEW (pyqtgraph)
self.plot_widget.plot(q, fit, pen=pg.mkPen('r', width=2))
```

**Cursors**:
```python
# OLD (matplotlib) - manual implementation
self.cursor_left_line = self.ax_main.axvline(x=cursor_left)
# + manual mouse event handlers

# NEW (pyqtgraph) - built-in!
self.cursor_left = pg.InfiniteLine(pos=cursor_left, angle=90, movable=True)
self.cursor_left.sigPositionChanged.connect(self.on_cursor_moved)
self.plot_widget.addItem(self.cursor_left)
```

**Log Scale**:
```python
# OLD (matplotlib)
self.ax_main.set_xscale('log')
self.ax_main.set_yscale('log')

# NEW (pyqtgraph)
self.plot_widget.setLogMode(x=True, y=True)
```

## File Structure After Refactoring

```
pyirena/
├── state/
│   ├── __init__.py
│   └── state_manager.py
├── gui/
│   ├── unified_fit.py (refactored for pyqtgraph)
│   └── data_selector.py (refactored for pyqtgraph)
├── plotting/
│   └── unified_plots.py (refactored for pyqtgraph)
└── core/
    └── unified.py (no changes needed)
```

## Dependencies

**New Requirement**:
```
pyqtgraph >= 0.12.0
```

**Remove**:
```
matplotlib (no longer needed)
```

## Benefits Summary

### Performance
- 10-100x faster rendering
- Smooth real-time updates
- Lower CPU usage

### User Experience
- Smoother cursor dragging
- Faster graph updates with auto-update
- Persistent state across sessions
- Quick reset to defaults

### Developer Experience
- Simpler cursor implementation
- Programmatic API for other tools
- Human-readable state files
- Easy preset sharing

## Risks & Mitigation

**Risk**: Breaking existing functionality
**Mitigation**: Comprehensive testing, keep old matplotlib version in git history

**Risk**: State file corruption
**Mitigation**: Merge with defaults, validate on load, backup on save

**Risk**: Pyqtgraph learning curve
**Mitigation**: Well-documented examples, simpler than matplotlib for this use case

## Timeline

1. **Phase 2 (Pyqtgraph)**: 2-3 hours
2. **Phase 3 (Integration)**: 1 hour
3. **Phase 4 (API)**: 30 minutes
4. **Phase 5 (Testing)**: 1 hour

**Total**: ~5 hours

## Next Steps

1. Implement pyqtgraph migration for UnifiedFitGraphWindow
2. Update cursor implementation to use InfiniteLine
3. Integrate state manager with GUI
4. Add "Reset to Defaults" button
5. Test and document
