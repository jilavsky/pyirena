# Quick Start: New Features

## What's New in v0.4.0

### ⚡ 10x Faster Performance
Everything is now **dramatically faster** thanks to pyqtgraph:
- Graphs render instantly
- Cursor dragging is silky smooth
- Auto-update is real-time

### 💾 Persistent State
Your settings are now **automatically saved**:
- Close the GUI - settings preserved
- Reopen - exactly where you left off
- No more re-entering parameters!

### 🔄 Reset to Defaults
New **orange "Reset to Defaults"** button:
- One click to restore factory settings
- Perfect for starting fresh
- Asks for confirmation first

### 📤 Export/Import Parameters
Share your fitting parameters:
- **Export Parameters** → Save to JSON file
- **Import Parameters** → Load from JSON file
- Share with colleagues or create presets

### 💡 Better Cursors
Cursors are now **much easier to use**:
- Smoother dragging
- Better visual feedback
- Built-in boundary checking

## How to Use

### Save Your Current Settings
```
1. Set up your parameters (levels, G, Rg, B, P, etc.)
2. Click "Save State" button (blue button in Results section)
3. Done! Settings are saved to ~/.pyirena/state.json
```

### Reset Everything
```
1. Click "Reset to Defaults" button (orange)
2. Confirm when asked
3. All parameters restored to factory defaults
```

### Create a Preset
```
1. Set up parameters for a specific type of fit (e.g., Guinier fit)
2. Click "Export Parameters"
3. Save as "guinier_preset.json"
4. Share with team or use on other machines!
```

### Use a Preset
```
1. Click "Import Parameters"
2. Select your preset file (e.g., "guinier_preset.json")
3. Parameters instantly applied!
```

### Programmatic Use (Advanced)
```python
from pyirena.state import StateManager

# Load Unified Fit settings
sm = StateManager()
params = sm.get('unified_fit')

# Use in your script
num_levels = params['num_levels']
G_value = params['levels'][0]['G']['value']

# Modify and save
params['num_levels'] = 2
sm.update('unified_fit', params)
sm.save()
```

## Performance Tips

### Enable Auto-Update
With the new pyqtgraph backend, auto-update is now **fast enough for real-time use**:

1. Check "Update automatically?"
2. Adjust parameters
3. See results instantly!

No more waiting - it's now interactive!

### Use Cursors Liberally
Cursor dragging is now **smooth and fast**:
- Drag to select fit range
- No lag or jank
- Works even with large datasets

## File Locations

### State File
```
~/.pyirena/state.json
```
Human-readable JSON file containing all your settings.

### Exported Parameters
Wherever you save them! Recommended:
```
~/pyirena_presets/
├── guinier_fit.json
├── porod_fit.json
├── two_level_fit.json
└── complex_system.json
```

## Keyboard Shortcuts (Future)

Coming soon:
- `Ctrl+S` - Save State
- `Ctrl+R` - Reset to Defaults
- `Ctrl+E` - Export Parameters
- `Ctrl+I` - Import Parameters

## Tips & Tricks

### Share Presets with Your Team
1. Create a shared folder (Dropbox, Google Drive, etc.)
2. Save presets there
3. Everyone can import the same parameters
4. Consistent fitting across the team!

### Quick Experiment Workflow
1. Start with defaults ("Reset to Defaults")
2. Adjust parameters interactively (auto-update on)
3. When happy, save: "Export Parameters"
4. Document what system it's for in filename

### Batch Processing
```python
from pyirena.state import StateManager

# Load your optimal parameters
sm = StateManager()
sm.import_tool_state('unified_fit', 'optimal_fit.json')
params = sm.get('unified_fit')

# Apply to multiple datasets
for data_file in my_data_files:
    # Use params for automated fitting
    # ...
```

## Migration from Old Version

### Your Data is Safe
- All existing functionality preserved
- API unchanged
- Just faster and with new features!

### No Action Required
- Settings will use defaults first time
- Adjust as needed
- Save state for next time

### If You Had Custom Scripts
- No changes needed for basic usage
- New state management is optional
- Old matplotlib backup available if needed

## Troubleshooting

### State File Corrupted?
```bash
# Delete and restart with fresh defaults
rm ~/.pyirena/state.json
```
Or use "Reset to Defaults" button in GUI.

### Missing pyqtgraph?
```bash
pip install pyqtgraph
```

### Want Old Matplotlib Version?
```python
# Old version still available:
from pyirena.gui.unified_fit_matplotlib_old import UnifiedFitPanel
```

## What's Next

Future improvements:
- Matplotlib export for publication-quality figures
- More tools migrated to pyqtgraph
- Cloud sync for settings (optional)
- Preset management dialog
- Keyboard shortcuts

## Feedback

Found a bug or have a suggestion?
- Create an issue on GitHub
- Include your state file if relevant
- Describe what you expected vs what happened

Enjoy the new speed! 🚀
