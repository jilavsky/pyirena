# Unified Fit GUI Guide

## Overview

The Unified Fit GUI provides an interactive interface for fitting small-angle scattering data using the Unified Fit model (Beaucage model). The GUI closely follows the design of the Igor Pro Irena package.

## Launching the Unified Fit GUI

### Method 1: From the Main GUI (recommended)
1. Launch the main pyIrena GUI: `pyirena-gui`
2. Select one or more data files
3. Click **Models → Unified Fit** from the menu bar, or click the **"Unified Fit"** button in the right panel

### Method 2: Standalone
```bash
python -m pyirena.gui.unified_fit
```

## GUI Layout

The Unified Fit panel uses a split-screen layout: controls on the left, graphs on the right.

```
┌──────────────────────────────────────────────────────────────────┐
│  Left Panel: Controls          │  Right Panel: Graphs            │
│                                │                                 │
│  - Top controls (Graph, Fit)   │  [I vs Q] [Porod (I·Q⁴)]  ← tabs
│  - Level tabs (1-5)            │                                 │
│  - Background & fit controls   │  ┌───────────────────────────┐  │
│  - Results buttons             │  │  Main Plot (data + fit)   │  │
│                                │  │  + draggable cursors      │  │
│                                │  └───────────────────────────┘  │
│                                │  ┌───────────────────────────┐  │
│                                │  │  Residuals                │  │
│                                │  └───────────────────────────┘  │
│                                │  (or switch to Porod tab)       │
└──────────────────────────────────────────────────────────────────┘
```

## Left Panel: Controls

### Top Controls

- **Graph Unified**: Calculate and plot the model with current parameters (no fitting)
- **Number of levels**: Select 1–5 structural levels
- **Update automatically?**: Auto-recalculate whenever a parameter changes (throttled to ~150 ms)
- **Display local (Porod & Guinier) fits?**: Overlay individual-level Guinier and power-law contributions on the graph
- **No limits?**: Disable parameter bounds during fitting

### Level Tabs (1–5)

Each level has its own tab. Only levels up to the selected "Number of levels" are active.

#### Guinier Region
- **G** — Guinier prefactor [cm⁻¹]: low-Q amplitude
- **Rg** — radius of gyration [Å]: characteristic size
- **Fit Rg/G btwn cursors**: local Guinier fit to data between the A/B cursors

#### Porod Region
- **Estimate B from G/Rg/P?**: auto-calculate B using the Hammouda relationship
- **B** — Porod constant [cm⁻¹ Å⁻ᴾ]: high-Q amplitude
- **P** — power-law slope (1–6)
- **Fit P/B btwn cursors**: local power-law fit to data between the cursors

#### Each parameter row has:
- Value field (supports mouse-wheel scrubbing: plain = ±1 %, Shift = ±10 %, Ctrl = ±0.1 %)
- **Fit?** checkbox — uncheck to hold parameter fixed during fitting
- **Low / High** limit fields (hidden when "No limits?" is checked)

#### Additional Level Controls
- **RgCutoff** — Beaucage transition radius between levels [Å]
- **Is this correlated system?** — enable Born-Green structure-factor correlations (ETA, PACK)
- **Sv [m²/cm³] / Invariant** — read-only surface-area and scattering-invariant estimates
- **Copy / Move / Swap level** — reorganize level parameters across tabs

### Background
- **SAS Background** — constant incoherent background [cm⁻¹]
- **Fit Bckg?** — include background in fitting

### Fitting Controls
- **Fit** — run least-squares optimisation (scipy TRF with bounds, or Nelder-Mead when "No limits?" is checked)
- **Revert back** — restore parameters to state before last fit
- **Reset Unified** — reset all parameters to defaults
- **Fix limits?** — lock current low/high values
- **Store local fits** — persist Guinier/Porod level fits for the "Display local fits" overlay

### Results Section
- **Store in Data Folder** — write fit result to the HDF5 source file
- **Export ASCII** — write Q, I(Q), fit, and residuals as a text file
- **Results to graphs** — annotate the main graph with level parameters
- **Passes / Calc. Uncertainty (MC)** — Monte Carlo uncertainty analysis: run N forward-model samples with noise to estimate parameter confidence intervals

## Right Panel: Graphs

The graph window has two tabs you can switch between at any time. Both update together whenever the model is recalculated.

### Tab 1 — "I vs Q"

**Main plot (top, ~80 % of height)**

- Log-log scale: Intensity [cm⁻¹] vs Q [Å⁻¹]
- Data points (blue) with error bars
- Unified Fit model curve (red)
- Individual level Guinier (green dashed) and Porod (green dotted) contributions, if "Display local fits" is checked
- SAS background level (brown dashed horizontal line), if enabled
- Secondary top axis showing feature radius R = π/Q [Å]
- Two draggable vertical cursors (**A** red, **B** blue) — define the Q range for "Fit btwn cursors" and for the global fit

**Residuals plot (bottom, ~20 % of height)**

- Log(Q) vs (Data − Fit) / σ
- Dashed horizontal line at 0
- A flat, random residuals trace indicates a good fit

**Graph interaction**

- Scroll wheel: zoom in/out
- Click-drag: pan
- Right-click → "Save graph as JPEG…" or "Save as Igor Pro ITX…"
- Zoom is preserved across auto-updates; only resets when new data is loaded

**Cursors**

- Drag cursor A (red) or cursor B (blue) to bracket the Q range of interest
- Cursor positions feed into "Fit Rg/G btwn cursors", "Fit P/B btwn cursors", and the global Fit Q range
- Cursors clamp to data boundaries when a new dataset is loaded

### Tab 2 — "Porod (I·Q⁴)"

- Log-log scale: I·Q⁴ [cm⁻¹·Å⁻⁴] vs Q [Å⁻¹]
- Same data (blue) and model (red) as Tab 1, transformed into Porod presentation
- Local Guinier and Porod level fits also shown (if enabled), scaled by Q⁴
- **No cursors** — this tab is for visual inspection only
- The Porod plot makes it much easier to judge how many levels are needed and where each level's Q range lies: a flat region indicates a well-described Porod scatterer; shoulders or bumps indicate additional levels
- Right-click → save as JPEG or ITX
- Zoom is preserved independently of Tab 1

## Unified Fit Model

The model sums contributions from N structural levels:

```
I(q) = Σᵢ [ Gᵢ exp(−q² Rgᵢ²/3)
           + exp(−q² Rgᵢ₋₁²/3) · Bᵢ · ([erf(q Rgᵢ/√6)]³/q)^Pᵢ ]
       + Background
```

Each level describes one population of scatterers at a characteristic size Rgᵢ.

### Parameter Table

| Parameter | Symbol | Units | Typical Range | Physical meaning |
|-----------|--------|-------|---------------|-----------------|
| Guinier prefactor | G | cm⁻¹ | 1–10¹⁰ | Low-Q scattering amplitude |
| Radius of gyration | Rg | Å | 1–10 000 | Characteristic size |
| Power-law slope | P | — | 1–6 | Surface / mass fractal |
| Porod constant | B | cm⁻¹Å⁻ᴾ | 10⁻²⁰–10¹⁰ | High-Q amplitude |
| Cutoff radius | RgCO | Å | 0–10 000 | Level-to-level transition |

### Typical P Values

| P | Physical interpretation |
|---|------------------------|
| 1 | Rods (1-D objects) |
| 2 | Platelets (2-D objects) |
| 3 | Rough surface (surface fractal) |
| 4 | Smooth / sharp interface (Porod law) |
| 5–6 | Mass fractals |

## Workflow

### 1. Load Data
Select a data file in the main GUI and open Unified Fit. Data appears immediately in the graph.

### 2. Judge the Number of Levels

Switch to the **Porod (I·Q⁴)** tab. Flat plateaus separated by transitions indicate distinct levels. Count the plateaus to estimate how many levels to use.

### 3. Initial Parameter Estimation

Start with **1 level** (increase later if residuals show structure).

For each level, working from largest structure (Level 1) to smallest:

- **G** — read the low-Q plateau height on the I-Q plot
- **Rg** — position the cursors in the Guinier region and click **Fit Rg/G btwn cursors**
- **P** — position the cursors in the power-law region and click **Fit P/B btwn cursors**
- **B** — use **Estimate B from G/Rg/P?** to get a starting value, or read from Fit P/B result

### 4. Fitting Strategy

**Sequential** (most reliable):
1. Check "Fit?" only for G and Rg of Level 1 → **Fit**
2. Then also enable B and P → **Fit**
3. Enable background if needed → **Fit**

**Cursor-based** (quick local estimates before global fit):
1. Bracket each region with cursors A and B
2. Click **Fit Rg/G btwn cursors** or **Fit P/B btwn cursors**
3. Then run the global **Fit**

**Full simultaneous fit**:
1. Set reasonable initial values for all parameters
2. Enable all "Fit?" checkboxes
3. Click **Fit** — may need tight bounds if starting far from solution

### 5. Add More Levels if Needed

If the residuals show systematic structure (humps, dips):
1. Increase "Number of levels"
2. Switch to the Porod tab to identify where the new level sits in Q
3. Set initial parameters for the new level
4. Fit

### 6. Assess Quality

- Residuals should be flat and random around 0
- χ² (shown in status bar after fit) should be close to 1 for well-scaled data
- Physical checks: Rg values should be ordered (Level 1 largest), P in 1–6, G positive

### 7. Save Results

- **Store in Data Folder** — writes to the HDF5 source file
- **Export ASCII** — saves Q, I, fit, residuals as a text file
- **Calc. Uncertainty (MC)** — run Monte Carlo uncertainty before saving if needed

## Tips and Best Practices

### Reading the Porod Tab
The Porod presentation (I·Q⁴ vs Q) is the fastest way to assess structural complexity:
- Flat regions = well-described Porod scattering at that scale
- Peaks or bumps = additional levels needed
- The slope of the tail in Porod space relates to P − 4

### Mouse-Wheel Scrubbing
All parameter fields support the mouse wheel with modifier keys:
- **No modifier** — ±1 % of current value per notch
- **Shift** — ±10 % (coarse)
- **Ctrl / Cmd** — ±0.1 % (fine)

Watch both graphs update in real time when "Update automatically?" is checked.

### Cursor Placement
- Place cursors to bracket a single structural feature
- For Level 1 Guinier: A at low-Q flat region onset, B just past the Guinier knee
- For Level 1 Porod: A just past the Guinier region, B before the next structural feature

### Common Issues

**Fit doesn't converge**
- Improve initial parameter guesses
- Reduce the number of simultaneously fitted parameters
- Tighten the bounds (set realistic low/high limits)

**Unphysical P values (< 1 or > 6)**
- Check if another level is needed in that Q range
- Add a constraint via the limit fields

**Large residuals at specific Q**
- Zoom in to that region on the I-Q tab to inspect the data
- Consider whether a Bragg peak or crystalline feature requires a different model
- An additional level may help

**Fit snaps to wrong region**
- Use cursor-based local fitting first to guide the global fit
- Fix the well-determined parameters and only fit the uncertain ones

## References

- Beaucage, G. "Approximations Leading to a Unified Exponential/Power-Law Approach to Small-Angle Scattering" *J. Appl. Cryst.* (1995) **28**, 717-728
- Beaucage, G. "Small-Angle Scattering from Polymeric Mass Fractals of Arbitrary Mass-Fractal Dimension" *J. Appl. Cryst.* (1996) **29**, 134-146
- Ilavsky, J. & Jemian, P. R. "Irena: tool suite for modeling and analysis of small-angle scattering" *J. Appl. Cryst.* (2009) **42**, 347-353

## Support

For issues or questions: https://github.com/jilavsky/pyirena/issues
