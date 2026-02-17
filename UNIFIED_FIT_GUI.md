# Unified Fit GUI Guide

## Overview

The Unified Fit GUI provides an interactive interface for fitting small-angle scattering data using the Unified Fit model (Beaucage model). The GUI closely follows the design of the Igor Pro Irena package.

## Launching the Unified Fit GUI

There are two ways to launch the Unified Fit panel:

### Method 1: From the Menu Bar
1. Launch the main pyIrena GUI: `pyirena-gui`
2. Select one or more data files
3. Click **Models → Unified Fit** from the menu bar

### Method 2: Using the Button
1. Launch the main pyIrena GUI: `pyirena-gui`
2. Select one or more data files
3. Click the **"Unified Fit"** button on the right panel (green button below "Create Graph")

### Method 3: Standalone
```bash
python -m pyirena.gui.unified_fit
```

## GUI Layout

The Unified Fit panel uses a split-screen layout:

```
┌─────────────────────────────────────────────────────────────┐
│  Left Panel: Controls     │  Right Panel: Graph             │
│                           │                                 │
│  - Model parameters       │  ┌──────────────────────────┐  │
│  - Level tabs (1-5)       │  │   Main Plot (data+fit)   │  │
│  - Fit controls           │  │                          │  │
│  - Results buttons        │  └──────────────────────────┘  │
│                           │  ┌──────────────────────────┐  │
│                           │  │   Residuals Plot         │  │
│                           │  └──────────────────────────┘  │
└─────────────────────────────────────────────────────────────┘
```

## Left Panel: Controls

### Top Controls

- **Graph Unified**: Calculate and plot the model with current parameters (no fitting)
- **Number of levels**: Select 1-5 structural levels
- **Update automatically?**: Auto-recalculate when parameters change
- **Display local (Porod & Guinier) fits?**: Show individual level contributions
- **No limits?**: Disable parameter bounds during fitting

### Level Tabs

Each level (1-5) has its own tab with identical controls:

#### Guinier Region Parameters
- **G**: Guinier prefactor (low-Q amplitude) [cm⁻¹]
  - Value input field
  - "Fit?" checkbox to enable fitting this parameter
  - Low and High limit fields (optional bounds)

- **Rg**: Radius of gyration [Å]
  - Value input field
  - "Fit?" checkbox
  - Low and High limit fields

- **Fit Rg/G btwn cursors**: Fit G and Rg using data between cursor positions

#### Porod Region Parameters
- **Estimate B from G/Rg/P?**: Auto-calculate B using Hammouda relationship

- **B**: Porod constant [cm⁻¹ Å⁻ᴾ]
  - Value input field
  - "Fit?" checkbox
  - Low and High limit fields

- **P**: Power law slope [dimensionless]
  - Value input field
  - "Fit?" checkbox
  - Low and High limit fields
  - Warning: "Level may not be physically feasible" (shown if P is unusual)

- **Fit P/B btwn cursors**: Fit P and B using data between cursor positions

#### Additional Level Controls
- **RgCutoff**: Transition Q value between levels [Å⁻¹]
- **pi B/Q [m2/cm3]**: Calculated Porod invariant (read-only)
- **Is this correlated system?**: Enable correlation function for this level
- **Copy/Move/swap level**: Manipulate level parameters

### Bottom Controls

#### Background
- **SAS Background**: Constant background value
- **Fit Bckg?**: Enable background fitting
- **Skip Fit Check?**: Skip convergence checks

#### Fitting Controls
- **Fit**: Run least-squares optimization
- **Revert back**: Restore previous parameters
- **reset unif?**: Reset all parameters to defaults
- **Fix limits?**: Lock current parameter bounds
- **Store local (Porod & Guinier) fits?**: Save individual level fits

#### Results Section
- **Store in Data Folder**: Save results to data folder
- **Export ASCII**: Export fit results as text file
- **Results to graphs**: Plot results with detailed analysis
- **Analyze Results**: Perform statistical analysis
- **Anal. Uncertainty**: Calculate parameter uncertainties
- **Ext. warnings?**: Show extended warning messages

## Right Panel: Graphs

### Main Plot (Top)
- Log-log scale
- Shows experimental data (points with error bars)
- Shows unified fit curve (solid line)
- Shows individual level contributions (if enabled)
- Interactive zoom and pan (using matplotlib toolbar)

### Residuals Plot (Bottom)
- Log-linear scale
- Shows (Data - Fit) / Error
- Horizontal line at y=0 for reference
- Helps assess fit quality

## Unified Fit Model

The model combines multiple structural levels:

For each level i:
```
I_i(q) = G_i × exp(-q² Rg_i² / 3)
         + exp(-q² Rg_{i-1}² / 3) × B_i × {[erf(k×q×Rg_i / √6)]³ / q}^P_i
```

Where:
- **Guinier term**: `G × exp(-q² Rg² / 3)` describes low-Q scattering
- **Porod term**: `B × q^(-P)` describes high-Q scattering
- **Transition**: Error function smoothly connects the regions

### Parameters for Each Level

| Parameter | Symbol | Units | Typical Range | Physical Meaning |
|-----------|--------|-------|---------------|------------------|
| Guinier prefactor | G | cm⁻¹ | 1 - 10¹⁰ | Low-Q scattering amplitude |
| Radius of gyration | Rg | Å | 1 - 10,000 | Size of scattering object |
| Power law slope | P | - | 1 - 6 | Surface/mass fractal dimension |
| Porod constant | B | cm⁻¹Å⁻ᴾ | 10⁻²⁰ - 10¹⁰ | High-Q scattering amplitude |
| RgCutoff | RgCO | Å | 0 - 10,000 | Transition Q between levels |

### Typical P Values

- **P = 1**: Rods (1D objects)
- **P = 2**: Platelets (2D objects)
- **P = 3**: Rough interface (Porod law)
- **P = 4**: Sharp interface (Porod law)
- **P > 4**: Mass fractals
- **P < 3**: Surface fractals

## Workflow Example

### 1. Load Data
- In main GUI, select data file(s)
- Click "Unified Fit" button
- Data appears in the graph

### 2. Set Number of Levels
- Start with 1 level
- Increase if residuals show structure
- Typically use 1-3 levels

### 3. Initial Parameter Estimation

For **Level 1** (largest structures):
1. Look at low-Q plateau → estimate G
2. Look at Guinier region → estimate Rg (where I drops to ~0.37 × G)
3. Look at high-Q slope → estimate P
4. Check "Estimate B from G/Rg/P?" or estimate from high-Q intensity

### 4. Fitting Strategy

**Approach 1: Sequential fitting**
1. Fit only G and Rg first (check only these boxes)
2. Click "Fit"
3. Then enable B and P fitting
4. Click "Fit" again

**Approach 2: Cursor-based fitting**
1. Use "Fit Rg/G btwn cursors" for low-Q region
2. Use "Fit P/B btwn cursors" for high-Q region

**Approach 3: Full fit**
1. Check all "Fit?" boxes
2. Click "Fit"
3. May need good initial guesses

### 5. Add More Levels if Needed

If residuals show structure:
1. Increase "Number of levels"
2. Switch to new level tab
3. Estimate parameters for smaller structures (higher Q)
4. Fit again

### 6. Save Results
- Click "Store in Data Folder" to save
- Click "Export ASCII" for text export
- Click "Results to graphs" for detailed plots

## Tips and Best Practices

### Parameter Estimation
- **G**: Y-intercept on log-log plot
- **Rg**: `Rg ≈ 1 / (Q at inflection point)
- **P**: Slope of high-Q region
- **B**: Intensity at high Q × Q^P

### Fitting Tips
1. **Start simple**: Begin with 1 level, add more only if needed
2. **Fix parameters**: Uncheck "Fit?" for well-known values
3. **Use bounds**: Set realistic low/high limits
4. **Check residuals**: Should be random, centered at 0
5. **Physical meaning**: Ensure Rg values make sense for your sample

### Common Issues

**Problem**: Fit doesn't converge
- **Solution**: Better initial guesses, reduce number of fitted parameters, add bounds

**Problem**: Unphysical P values
- **Solution**: Check data quality, consider different structural model, add constraints

**Problem**: Large residuals at specific Q
- **Solution**: May need additional level, check data quality in that region

## Keyboard Shortcuts

*(To be implemented)*

| Shortcut | Action |
|----------|--------|
| **Ctrl+F** | Run fit |
| **Ctrl+R** | Revert parameters |
| **Ctrl+G** | Graph unified model |
| **Ctrl+S** | Save results |

## Implementation Status

### ✅ Implemented
- [x] GUI layout matching Igor design
- [x] Level tabs (1-5 levels)
- [x] All parameter input fields
- [x] Fit checkboxes
- [x] Parameter bounds (low/high limits)
- [x] Graph panel with main plot + residuals
- [x] Integration with main data selector
- [x] Menu bar access
- [x] Button access from main GUI
- [x] Data loading and display

### 🚧 In Progress
- [ ] Connect to UnifiedFitModel backend
- [ ] Implement fitting algorithm
- [ ] Cursor-based fitting
- [ ] Auto-calculate B from G/Rg/P
- [ ] Display local fits
- [ ] Parameter uncertainty analysis

### 📋 Planned
- [ ] Save/load parameter sets
- [ ] Results export (ASCII, graphs)
- [ ] Batch fitting multiple files
- [ ] Parameter correlation analysis
- [ ] Copy/move/swap levels
- [ ] Undo/redo functionality
- [ ] Real-time parameter updates

## Testing

To test the GUI:

```bash
# Test unified fit GUI standalone
python -m pyirena.gui.unified_fit

# Test from main GUI
pyirena-gui
# Then: select file → click "Unified Fit" button
```

## References

- Beaucage, G. "Approximations Leading to a Unified Exponential/Power-Law Approach to Small-Angle Scattering" *J. Appl. Cryst.* (1995) **28**, 717-728
- Beaucage, G. "Small-Angle Scattering from Polymeric Mass Fractals of Arbitrary Mass-Fractal Dimension" *J. Appl. Cryst.* (1996) **29**, 134-146
- Ilavsky, J. & Jemian, P. R. "Irena: tool suite for modeling and analysis of small-angle scattering" *J. Appl. Cryst.* (2009) **42**, 347-353

## Support

For issues or questions:
- GitHub Issues: https://github.com/yourusername/pyirena/issues
- Based on original Irena package by Jan Ilavsky
