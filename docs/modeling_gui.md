# Modeling Tool

The Modeling tool performs parametric forward-modeling of small-angle scattering data.
It sums contributions from up to 10 independent **populations**, each of which can be a
Size Distribution, a Unified Fit Level (Beaucage equation), or a Diffraction Peak.
The total model I(Q) is fitted to experimental data using non-linear least squares.

---

## Contents

1. [Overview](#overview)
2. [Launching the tool](#launching-the-tool)
3. [Panel layout](#panel-layout)
4. [Population types](#population-types)
   - [Size Distribution](#41-size-distribution)
   - [Unified Fit Level](#42-unified-fit-level)
   - [Diffraction Peak](#43-diffraction-peak)
5. [Distributions](#distributions)
6. [Form Factors](#form-factors)
7. [Structure Factors](#structure-factors)
8. [Population label and color coding](#population-label-and-color-coding)
9. [Fitting controls](#fitting-controls)
10. [MC Uncertainty](#mc-uncertainty)
11. [Derived results](#derived-results)
12. [Revert and state persistence](#revert-and-state-persistence)
13. [HDF5 output format](#hdf5-output-format)
14. [Integration with other tools](#integration-with-other-tools)
15. [Scripting API](#scripting-api)

---

## Overview

Unlike the Sizes tool (which inverts the scattering matrix to recover a size distribution),
the Modeling tool works **forward**: you define a parametric model and fit its parameters
to the data. This gives you full control over the physical model and meaningful uncertainties
for each parameter.

Typical use cases:
- Multi-population systems (e.g., small primary particles + large agglomerates)
- Systems where a power-law or Guinier term is needed alongside a size distribution
- Fitting diffraction peaks on top of a smooth scattering background
- Mixed morphologies (e.g., spheres + cylinders modeled as two populations)

---

## Launching the tool

**From the Data Selector:**
- Select one or more HDF5 files in the file list
- Click **Modeling (GUI)** (orange button, row 4)

The tool opens with the last-used file pre-loaded and all settings restored from the
previous session.

**From the command line:**

```bash
python -m pyirena.gui.modeling_panel path/to/file.h5
```

**Headless / batch fitting** — see [Scripting API](#scripting-api).

---

## Panel layout

```
┌────────────────────────────────────┬──────────────────────────────────────────┐
│  Modeling Input          [? Help]  │  I(Q) graph                              │
│  ──────────────────────────────    │                                          │
│  Data file: [____________] [Open…] │   log I(Q) vs log Q                      │
│                                    │   data (grey) + total model (black)      │
│  ┌──────────────────────────────┐  │   per-population curves (colored)        │
│  │ P1 │ P2 │ P3 │ … │ P10     │  │                                          │
│  │                              │  │  [Q cursor left]  [Q cursor right]      │
│  │  Population type: [combo]    │  │                                          │
│  │  ─────────────────────────── │  ├──────────────────────────────────────────┤
│  │  [Size Dist / UF / Peak      │  │  Derived results (per population)        │
│  │   parameters]                │  │                                          │
│  └──────────────────────────────┘  ├──────────────────────────────────────────┤
│                                    │  Status bar                              │
│  Background: [____] [✓ Fit]        └──────────────────────────────────────────┘
│  [✓ No limits?]
│  Qmin: [____]  Qmax: [____]
│  MC passes: [10]
│
│  [Graph Model]  [Fit]  [Calc. Uncertainty (MC)]
│  [Revert]       [Save JSON]  [Load JSON]
│  [Save HDF5]
└────────────────────────────────────┘
```

---

## Population types

Each population tab has a **Population type** combo at the top with three options:

| Type | Description |
|------|-------------|
| **Size Distribution** | Parametric size distribution convolved with a form factor G-matrix |
| **Unified Fit Level** | Single Beaucage level: Guinier + power-law (+ optional Born-Green correlations) |
| **Diffraction Peak** | Gaussian, Lorentzian, or pseudo-Voigt peak at a chosen Q₀ |

Switching population type preserves the parameters of all three types independently —
switching back restores the values from before.

### 4.1 Size Distribution

The Size Distribution population uses a **parametric distribution** (see [Distributions](#distributions))
to describe particle sizes, a **form factor** to compute scattering per unit volume fraction,
and an optional **structure factor** for inter-particle correlations.

**Parameters:**

| Parameter | Description |
|-----------|-------------|
| Distribution | Choice of size distribution shape |
| Distribution params | Shape-specific parameters (see Distributions table) |
| Scale | Overall scale factor (proportional to volume fraction × (1 − volume fraction) for hard-sphere-like systems) |
| Contrast | (Δρ)² in units of 10²⁰ cm⁻⁴ |
| Form factor | Particle shape (see Form Factors table) |
| FF params | Extra form factor parameters (e.g. aspect ratio for spheroids) |
| Structure factor | Inter-particle structure factor |
| SF params | Structure factor parameters |
| Use number dist. | Fit a number-weighted distribution instead of volume-weighted |

**Derived quantities** (shown in the Derived panel after Graph Model or Fit):

| Quantity | Description |
|----------|-------------|
| Volume fraction | Total particle volume fraction |
| Vol. mean radius | Volume-weighted mean radius [Å] |
| Num. mean radius | Number-weighted mean radius [Å] |
| Specific surface | Specific surface area [Å⁻¹] |

### 4.2 Unified Fit Level

Implements a single level of the Beaucage Unified model:

```
I(Q) = G · exp(−Q²Rg²/3) + B · Q*⁻ᴾ · exp(−Q²RgCO²/3)

where  Q* = Q / [erf(K·Q·Rg/√6)]³
```

**Parameters:**

| Parameter | Description | Fitted by default |
|-----------|-------------|-------------------|
| G | Guinier pre-factor [cm⁻¹] | Yes |
| Rg | Radius of gyration [Å] | Yes |
| B | Power-law pre-factor | Yes |
| P | Power-law exponent | No |
| RgCO | Low-Q cutoff radius of gyration (0 = no cutoff) | No |

**Correlations (Born-Green):** When enabled, the Unified Fit Level intensity is divided
by a Born-Green structure factor S(Q) = 1 / (1 + PACK × F(Q, ETA)):

| Parameter | Description |
|-----------|-------------|
| ETA [Å] | Correlation length |
| PACK | Packing factor |

### 4.3 Diffraction Peak

Models a single scattering peak at position Q₀.

**Peak shape:**

| Shape | Formula |
|-------|---------|
| Gaussian | A · exp[−(Q−Q₀)² / (2σ²)] |
| Lorentzian | A / [1 + ((Q−Q₀)/σ)²] |
| Pseudo-Voigt | A · [η·Lorentzian + (1−η)·Gaussian] |

**Parameters:**

| Parameter | Description | Fitted by default |
|-----------|-------------|-------------------|
| Position Q₀ [Å⁻¹] | Peak centre | Yes |
| Amplitude [cm⁻¹] | Peak height | Yes |
| Width σ [Å⁻¹] | Gaussian/Lorentzian width | Yes |
| η (mixing) | Voigt mixing ratio (0=Gaussian, 1=Lorentzian) | No |

---

## Distributions

Available for Size Distribution populations:

| Distribution | Parameters | Shape |
|-------------|-----------|-------|
| **Log-Normal** | Median [Å], σ (width), Shift [Å] | Asymmetric, right-tailed |
| **Gaussian** | Mean [Å], σ (width) | Symmetric bell curve |
| **LSW** | Location [Å] | Fixed shape (Lifshitz-Slyozov-Wagner coarsening) |
| **Schulz-Zimm** | Mean [Å], Z (sharpness) | Gamma-like, narrow for large Z |
| **Ardell** | Location [Å], m (≥ 2) | Asymmetric, hard cutoff at m/(m−1) × location |

**Note on Ardell distribution:** The location parameter is **not** the mean or mode.
The distribution has a hard upper cutoff at `r_max = m/(m−1) × location` and peaks
at approximately 1.35× location for typical m values. The formula follows
Ardell (1966) Acta Metall. 14, 1573.

---

## Form Factors

| Form Factor | Key | Extra parameters | Description |
|------------|-----|-----------------|-------------|
| Sphere | `sphere` | — | Solid sphere, F = 3[sin(Qr)−Qr·cos(Qr)]/(Qr)³ |
| Spheroid | `spheroid` | Aspect ratio (AR) | Oblate (AR<1) or prolate (AR>1) ellipsoid of revolution; orientation-averaged via 50-point Gauss-Legendre quadrature |

The G-matrix element is `G[i,j] = V(r_j) × F²(Q_i, r_j) × (Δρ)² × 10⁻⁴` in cm⁻¹.

---

## Structure Factors

| Structure Factor | Key | Parameters | Description |
|----------------|-----|-----------|-------------|
| None | `none` | — | S(Q) = 1 |
| Interferences (Born-Green) | `interferences` | ETA [Å], PACK | S(Q) = 1/(1 + PACK·F(Q,ETA)) |
| Hard Sphere (Percus-Yevick) | `hard_sphere` | Radius [Å], Vol. frac. | Analytical hard-sphere S(Q) |

---

## Population label and color coding

Each population tab can have an optional **label** (text field in the tab header area).
The label appears in the tab as `P1: label_text` (truncated to 14 characters).

Population colors are fixed per index:

| Pop | Color |
|-----|-------|
| P1  | Blue (#2980b9) |
| P2  | Red (#e74c3c) |
| P3  | Green (#27ae60) |
| P4  | Purple (#8e44ad) |
| P5  | Orange (#e67e22) |
| P6–P10 | Additional colors |

---

## Fitting controls

| Control | Description |
|---------|-------------|
| **Background** | Flat additive background [cm⁻¹]; check **Fit** to optimize it |
| **No limits?** | When checked, fitting is unconstrained (Nelder-Mead); unchecked uses TRF with bounds |
| **Qmin / Qmax** | Q range for fitting (drag the vertical cursor lines in the graph) |
| **Graph Model** | Compute and display the forward model without fitting |
| **Fit** | Run least-squares optimization and update all parameters |

**Fit parameters:** Each population has "Fit?" checkboxes per parameter. Unchecked
parameters are held fixed during optimization. This lets you fix physically known
quantities (e.g. P=4 for Porod scattering) while fitting others.

**Fitting algorithm:**
- With limits: `scipy.optimize.least_squares` with Trust Region Reflective (TRF)
- No limits: `scipy.optimize.minimize` with Nelder-Mead

---

## MC Uncertainty

Click **Calc. Uncertainty (MC)** to estimate parameter uncertainties by Monte Carlo:

1. The last fit result is used as the starting point
2. `N` refit runs with Gaussian noise added to I(Q) (scaled by dI)
3. Standard deviations of the fitted parameters are reported in the Derived panel

Set the number of passes with the **Passes:** spin box (1–500, default 10).
More passes give better uncertainty estimates at the cost of longer compute time.

---

## Derived results

The Derived Results panel (below the graph) shows post-fit quantities:

| Quantity | Populations |
|----------|-------------|
| Volume fraction | Size Distribution |
| Vol. mean radius [Å] | Size Distribution |
| Num. mean radius [Å] | Size Distribution |
| Specific surface [Å⁻¹] | Size Distribution |
| G, Rg, B, P | Unified Fit Level |
| Position, Amplitude, Width | Diffraction Peak |

---

## Revert and state persistence

- **Revert** restores all parameters to the state they were in just before the last **Fit** run.
  This lets you undo a bad fit without manually resetting values.
- Parameters, population settings, background, and Q-range cursors are automatically
  saved to the pyIrena state file and restored the next time the tool is opened
  (even if the application is closed).
- **Save JSON / Load JSON** exports/imports all population parameters to a portable
  JSON file for sharing or archiving.

---

## HDF5 output format

Modeling results are saved under `entry/modeling_results/` in the same HDF5 file as the
experimental data (NXcanSAS format). The schema:

```
entry/modeling_results/
├── chi_squared              scalar
├── reduced_chi_squared      scalar
├── dof                      scalar (int)
├── background               scalar [cm⁻¹]
├── q_min, q_max             scalars [Å⁻¹]
├── model_q                  array [Å⁻¹]
├── model_I                  array [cm⁻¹]   (total model)
├── background_err           scalar (if MC was run)
│
├── pop_01/                  (first population)
│   ├── @pop_type            attribute: 'size_dist' | 'unified_level' | 'diffraction_peak'
│   ├── @label               attribute: user-defined label string
│   ├── @enabled             attribute: bool
│   ├── @population_index    attribute: 1
│   ├── model_I              array [cm⁻¹]   (this population's contribution)
│   │
│   │   [unified_level only]
│   ├── G, Rg, B, P, RgCO, ETA, PACK   scalars
│   │
│   │   [diffraction_peak only]
│   ├── @peak_type           attribute
│   ├── position, amplitude, width, eta_voigt   scalars
│   │
│   │   [size_dist only]
│   ├── @dist_type, @form_factor, @structure_factor   attributes
│   ├── scale, contrast      scalars
│   ├── dist_params/         sub-group (distribution parameters)
│   ├── ff_params/           sub-group (form factor parameters)
│   ├── sf_params/           sub-group (structure factor parameters)
│   ├── radius_grid          array [Å]
│   ├── volume_dist          array
│   └── number_dist          array
│
├── pop_02/, pop_03/, …      (additional populations)
```

---

## Integration with other tools

### HDF5 Viewer

Files containing Modeling results are automatically detected. In the **1D Graph** tab:
- Check **Modeling model** to overlay the total I(Q) model curve
- Right-click the `modeling_results` group in the browser tree → **Plot Modeling model**

In the **Collect Values** tab:
- Set Type = **Modeling**
- Choose Population (1–10) and Item (chi2, background, Rg, G, B, P, position, vol_fraction, …)
- Collect from all selected files into a single graph

### Data Selector

In the "Show in graph/reports:" section, check **Modeling** to include modeling results in:
- **Create Graph** — overlays total model I(Q) curves for all selected files
- **Create Report** — adds a Modeling section to the Markdown report, with chi², background,
  and per-population parameters (type-specific)
- **Tabulate Results** — adds `MOD_` columns to the CSV table (chi², background, Q range,
  per-population type/label/parameters)

---

## Scripting API

```python
from pyirena.batch import fit_modeling
from pathlib import Path

result = fit_modeling(
    data_file='path/to/data.h5',
    config_file='path/to/config.json',   # JSON exported from Save JSON
    save_to_nexus=True,                  # save results back into data.h5
    with_uncertainty=True,               # run MC uncertainty estimation
    n_mc_runs=50,
)

if result['success']:
    mod_result = result['result']
    print(f"chi² = {mod_result.chi_squared:.4f}")
    for pop in mod_result.populations:
        print(f"  {pop.pop_type}: {pop}")
```

The config JSON is the same format produced by **Save JSON** in the GUI and can be
generated programmatically by constructing a `ModelingConfig` object:

```python
from pyirena.core.modeling import ModelingConfig, SizeDistPopulation

config = ModelingConfig(
    populations=[
        SizeDistPopulation(
            enabled=True,
            dist_type='lognormal',
            dist_params={'median': 100.0, 'sigma': 0.3, 'shift': 0.0},
            scale=0.001,
            contrast=1.0,
            form_factor='sphere',
        )
    ],
    background=0.0,
    fit_background=True,
    q_min=0.01,
    q_max=0.3,
)
```
