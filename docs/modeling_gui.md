# Modeling Tool

The Modeling tool performs parametric forward-modeling of small-angle scattering data.
It sums contributions from up to 10 independent **populations**, each of which can be a
Size Distribution, a Unified Fit Level (Beaucage equation), a Diffraction Peak,
a Guinier-Porod Level, a Mass Fractal, or a Surface Fractal.
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
   - [Guinier-Porod Level](#44-guinier-porod-level)
   - [Mass Fractal](#45-mass-fractal)
   - [Surface Fractal](#46-surface-fractal)
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

Each population tab has a **Population type** combo at the top with six options:

| Type | Description |
|------|-------------|
| **Size Distribution** | Parametric size distribution convolved with a form factor G-matrix |
| **Unified Fit Level** | Single Beaucage level: Guinier + power-law (+ optional Born-Green correlations) |
| **Diffraction Peak** | Gaussian, Lorentzian, or pseudo-Voigt peak at a chosen Q₀ |
| **Guinier-Porod Level** | Piecewise Guinier-Porod model for non-spherical / non-dilute systems |
| **Mass Fractal** | Fractal aggregate scattering (sphere primary particles + fractal S(Q)) |
| **Surface Fractal** | Scattering from a surface with fractal roughness |

Switching population type preserves the parameters of all six types independently —
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

### 4.4 Guinier-Porod Level

Implements the Guinier-Porod piecewise scattering model
(Hammouda 2010, *J. Appl. Cryst.* **43**, 716–719).
Suitable for non-spherical or elongated scatterers, and for systems where the low-Q
slope deviates from the Guinier-regime flat behaviour.

**Formula:** Three regimes joined at crossover Q values:
```
Q1 = sqrt((P − s1)(3 − s1)/2) / Rg1
D  = G · exp(−(P−s1)/2) · ((3−s1)(P−s1)/2)^((P−s1)/2) / Rg1^(P−s1)
Q2 = sqrt((1 − s2) / (2Rg2²/(3−s2) − 2Rg1²/(3−s1)))  [0 if non-finite]
G2 = G · exp(−Q2²·(Rg1²/(3−s1) − Rg2²/(3−s2))) · Q2^(s2−s1)

I(Q) = G2/Q^s2 · exp(−Q²Rg2²/(3−s2))   Q < Q2  (active only if Q2 > 0)
I(Q) = G/Q^s1 · exp(−Q²Rg1²/(3−s1))    Q2 ≤ Q < Q1
I(Q) = D/Q^P                             Q ≥ Q1
```
Setting Rg2 = 1e10 (default) collapses Q2 to zero — single-level behaviour.

Optional: `I *= exp(−RgCO²·Q²/3)` (RgCO > 0).
Optional Born-Green correlations: `I /= (1 + PACK · F(Q, ETA))`.

**Parameters:**

| Parameter | Default | Limits | Fitted | Description |
|-----------|---------|--------|--------|-------------|
| G [cm⁻¹] | 1.0 | 1e-10 … 1e10 | Yes | Guinier amplitude |
| Rg1 [Å] | 10.0 | 0.1 … 1e6 | Yes | Radius of gyration (high-Q level) |
| Slope s1 | 0.0 | 0 … 3 | No | Low-Q slope (0=sphere, 1=rod, 2=lamella) |
| Power P | 4.0 | 0 … 6 | No | Porod exponent at high Q |
| Rg2 [Å] | 1e10 | 0.1 … 1e12 | No | Rg of second (lower-Q) level; 1e10=inactive |
| Slope s2 | 0.0 | 0 … 3 | No | Low-Q slope of second level |
| RgCO [Å] | 0.0 | 0 … 1e6 | No | Cutoff Rg (0 = no cutoff) |
| ETA [Å] | 10.0 | 0.1 … 1e6 | No | Born-Green correlation length |
| PACK | 0.0 | 0 … 16 | No | Born-Green packing factor |

### 4.5 Mass Fractal

Implements the Teixeira (1988) mass fractal aggregate scattering model
(*J. Appl. Cryst.* **21**, 781–785).
Primary particles are spheroids with aspect ratio β; the fractal aggregate structure is described analytically.
β = 1 gives the monodisperse sphere form factor (Bessel function oscillations in the Porod region);
β ≠ 1 uses the orientation-averaged spheroid form factor which eliminates those oscillations.
A modest β ≈ 0.5–2 is often sufficient to obtain a physically smooth I(Q).

**Formula:**
```
V  = (4/3)·π·Radius³
P(Q) = (3·(sin(qR) − qR·cos(qR)) / (qR)³)²        sphere form factor

Bracket = Eta · 8 · (Ksi / (2·Radius))^Dv

I(Q) = Phi · Contrast · 1e-4 · V
       · [Bracket · sin((Dv−1)·atan(Q·Ksi)) / ((Dv−1)·Q·Ksi·(1+(Q·Ksi)²)^((Dv−1)/2))
          + (1−Eta)²]
       · P(Q)
```
The factor `1e-4` is the Igor-convention contrast unit conversion
(contrast in Å⁻⁴ units → cm⁻⁴ scale). Eta (volume filling factor) is typically 0.3–0.8.

**Parameters:**

| Parameter | Default | Limits | Fitted | Description |
|-----------|---------|--------|--------|-------------|
| Phi | 0.001 | 1e-8 … 1 | Yes | Volume fraction of primary particles |
| Radius [Å] | 50.0 | 0.1 … 1e6 | Yes | Primary particle equatorial radius |
| Aspect ratio β | 1.0 | 0.01 … 100 | No | Spheroid aspect ratio (1 = sphere, < 1 = oblate, > 1 = prolate) |
| Fractal dim. Dv | 2.5 | 1 … 3 | Yes | Mass fractal dimension |
| Ksi [Å] | 500.0 | 1 … 1e7 | Yes | Fractal correlation length (aggregate size) |
| Eta | 0.5 | 0.3 … 0.8 | No | Volume filling factor within the aggregate |
| Contrast | 1.0 | 0 … 1e10 | No | Scattering contrast (Δρ)² |

### 4.6 Surface Fractal

Implements the Teixeira (1988) surface fractal scattering model
(*J. Appl. Cryst.* **21**, 781–785).
Describes scattering from a surface with fractal roughness (2 ≤ Ds ≤ 3).

**Formula:**
```
I(Q) = π · Contrast · 1e20 · Ksi⁴ · 1e-32 · Surface · Γ(5−Ds)
       · sin((3−Ds)·atan(Q·Ksi))
       / ((1+(Q·Ksi)²)^((5−Ds)/2) · Q·Ksi)
```
Limiting slopes: `I(Q) ~ Q^(Ds−6)` for Q·Ksi >> 1 (fractal regime, slope −3 to −4).

**Optional Porod transition:** Above Qc, smoothly blends to `A·Q⁻⁴` using an error-function
step, with continuity condition `A = I(Qc)·Qc⁴`.

**Parameters:**

| Parameter | Default | Limits | Fitted | Description |
|-----------|---------|--------|--------|-------------|
| Surface [cm⁻¹] | 1e4 | 1 … 1e12 | Yes | Surface area per unit volume |
| Fractal dim. Ds | 2.5 | 2 … 3 | Yes | Surface fractal dimension |
| Ksi [Å] | 500.0 | 1 … 1e7 | Yes | Correlation length (upper cutoff of fractal regime) |
| Contrast | 1.0 | 0 … 1e10 | No | Scattering contrast |
| Qc [Å⁻¹] | 0.1 | 0.001 … 10 | No | Crossover to Porod Q⁻⁴ (enabled by checkbox) |
| QcWidth | 0.1 | 0.01 … 1 | No | Width of erf blending as fraction of Qc |

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
| Spheroid | `spheroid` | Aspect ratio AR (dimensionless) | Oblate (AR<1) or prolate (AR>1) ellipsoid of revolution; orientation-averaged via 50-point Gauss-Legendre quadrature |
| Cylinder (Aspect Ratio) | `cylinder_ar` | Aspect ratio AR = L/R (dimensionless) | Finite cylinder; "size" R is the radius, half-length L = AR·R. Disk-like for AR≪1, rod-like for AR≫1. Orientation-averaged. |
| Cylinder (Length) | `cylinder_length` | Length H [Å] (total height) | Finite cylinder; "size" R is the radius, total height H is a fixed (or fitted) parameter independent of R. V = πR²·H scales as R². |

The G-matrix element is `G[i,j] = V(r_j) × F²(Q_i, r_j) × (Δρ)² × 10⁻⁴` in cm⁻¹.

### Cylinder form factor details

Both cylinder entries use the same orientationally-averaged form factor:

```
⟨F²(Q)⟩ = ∫₀¹ [2J₁(QR√(1−u²)) / (QR√(1−u²))]² · [sin(QL·u) / (QL·u)]² du
```

where u = cos(α), α is the angle between Q and the cylinder axis, R is the equatorial
radius, and L is the half-length. The integral is evaluated by 50-point Gauss-Legendre
quadrature (same as spheroid).

**Derived quantities for cylinder populations:**

| Quantity | Cylinder formula |
|----------|-----------------|
| Rg | √(R²/2 + L²/3), volume-weighted over distribution |
| Specific surface Sv | 2/L + 2/R (= 2(R+L)/(RL)), volume-weighted |

For `cylinder_ar`: L = AR·R (varies with R), so Rg = R·√(1/2 + AR²/3).  
For `cylinder_length`: L = H/2 is constant, so Rg² = ⟨R²⟩/2 + H²/12.

**Parameterisation guide:**

- **Disk-like particles** (platelets, lamellae): use `cylinder_ar` with AR < 0.5, or `cylinder_length` with H ≪ ⟨R⟩.
- **Rod-like particles** (fibres, elongated): use `cylinder_ar` with AR > 2, or `cylinder_length` with H ≫ ⟨R⟩.
- Use `cylinder_ar` when the aspect ratio is physically meaningful and expected to be correlated with particle size.
- Use `cylinder_length` when the thickness/length is a fixed structural parameter (e.g., bilayer thickness, layer period).

### Core-Shell form factor details

Core-shell form factors embed the scattering length densities (SLDs) directly into the form factor calculation. The **Contrast** field is automatically locked to 1.0 when a core-shell form factor is selected — do not change it; contrast is encoded by the SLD difference parameters.

**Available core-shell entries:**

| Key | Combo label | "Size" axis | Extra FF parameters |
|-----|-------------|-------------|---------------------|
| `cs_sphere_by_core` | Core-Shell Sphere (by core R) | R_core | sld_core, sld_shell, sld_solvent, t_shell |
| `cs_sphere_by_shell` | Core-Shell Sphere (by shell t) | t_shell | sld_core, sld_shell, sld_solvent, r_core_fixed |
| `cs_sphere_by_total` | Core-Shell Sphere (by total R) | R_total | sld_core, sld_shell, sld_solvent, t_shell |
| `cs_spheroid_by_core` | Core-Shell Spheroid (by core R) | R_core | sld_core, sld_shell, sld_solvent, t_shell, aspect_ratio |
| `cs_spheroid_by_total` | Core-Shell Spheroid (by total R) | R_total | sld_core, sld_shell, sld_solvent, t_shell, aspect_ratio |

**SLD parameters** are in units of **10⁻⁶ Å⁻²** (standard SLD units). Typical values:

| Material | SLD (10⁻⁶ Å⁻²) | Notes |
|----------|----------------|-------|
| H₂O | 9.46 | X-ray solvent default |
| D₂O | 6.36 | Neutron solvent default |
| Silica SiO₂ | 18.85 | X-ray |
| Gold Au | 121 | X-ray |
| Polystyrene | 9.5 | X-ray |

**Scattering amplitude:**
```
F_cs(Q, R_core, R_total) =
    (ρ_core − ρ_shell)   · V_core  · f_sph(Q·R_core)  +
    (ρ_shell − ρ_solvent) · V_total · f_sph(Q·R_total)
```
where `f_sph(x) = 3(sin x − x·cos x)/x³` (normalised sphere amplitude) and V = (4/3)πR³.

**Polydispersity mode guide:**
- **by_core** — polydisperse core radius (most common; shell thickness uniform). R_total = R_core + t_shell.
- **by_shell** — polydisperse shell thickness (core radius fixed). R_total = r_core_fixed + t_shell.
- **by_total** — polydisperse outer radius (shell thickness uniform). R_core = R_total − t_shell.

**Derived quantities for core-shell populations:**

| Quantity | Description |
|----------|-------------|
| vol_mean_r | Volume-weighted mean of the distributed dimension (R_core, t_shell, or R_total) |
| r_total_mean | Volume-weighted mean total (outer) radius [Å] |
| volume_fraction | Derived from `scale` parameter (based on total particle volume) |
| Rg | SLD-contrast-weighted radius of gyration (matches Guinier plot) |
| specific_surface | Based on outer surface: Sv = 3/R_total (Porod law, outer interface) |

> **Volume convention (important):** The `scale` parameter and `volume_fraction`
> derived quantity are based on the **total particle volume** V_total = (4/3)π R_total³
> (core + shell together), **not** the core volume alone.  This means: if a sample
> contains 1% by volume of core-shell particles, `scale ≈ Vf·(1−Vf) ≈ 0.01`.
> The same convention applies to the number↔volume distribution conversion
> when "Number distribution" is selected: the volume used is V_total, not V_core.
>
> The G-matrix is normalised as G = |F_cs|²/V_total × 1e-4, so that
> `I(Q) = G × vol_dist × dr` gives the correct cm⁻¹ intensity when vol_dist
> integrates to the total-particle volume fraction.

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

### Data Explorer

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

## Fit Quality

Residuals are displayed **rescaled** as r' = r / σ(robust): the normalized
residual divided by a robust (MAD-based) estimate of the actual noise scale, so
the scatter is judged against the data's own noise floor rather than the (often
mis-scaled) reported σ. The fit status line also reports **σ-scale** (how many ×
the actual scatter exceeds the reported σ), the realistic reduced-χ² floor, and
the largest fractional misfit **max|(I−M)/I|** — a σ-independent gross-misfit
backstop.

See the **[Fit Quality Metrics guide](fit_quality_metrics.md)** for full
interpretation of these numbers.

