# Simple Fits Tool

The Simple Fits tool provides direct analytical model fitting for common small-angle
scattering functions.  It is ported from the Igor Pro `IR3_SimpleFits.ipf` and
`IR3_SystemSpecificModels.ipf` modules of the Irena package.

---

## Contents

1. [Overview](#overview)
2. [Launching the tool](#launching-the-tool)
3. [Panel layout](#panel-layout)
4. [Models](#models)
5. [Background options](#background-options)
6. [Q range selection](#q-range-selection)
7. [Fitting workflow](#fitting-workflow)
8. [Linearization plots](#linearization-plots)
9. [Monte Carlo uncertainty](#monte-carlo-uncertainty)
10. [Storing and loading results](#storing-and-loading-results)
11. [Export / Import Parameters (batch scripting)](#export--import-parameters)
12. [Data Selector integration](#data-selector-integration)
13. [Scripting API](#scripting-api)

---

## Overview

Simple Fits fits a single analytical scattering model to I(Q) data over a
user-selected Q range.  It is most useful when the scattering is well described
by one function — for example, a Guinier region at low Q, a Porod tail at high Q,
or scattering from monodisperse spherical particles.

For complex multi-level structures use the **Unified Fit** tool instead.
For model-independent particle size distributions use the **Size Distribution** tool.

---

## Launching the tool

From the Data Selector, select one file and click **"Simple Fits (GUI)"**.  The
panel opens with that file's data loaded and plotted.

To process many files without the GUI, use **"Simple Fits (script)"** — this reads
the `simple_fits` section from `pyirena_config.json` and fits all selected files,
saving results to each HDF5 file automatically.

---

## Panel layout

```
┌─────────────────────────────────────────────────────────────────────┐
│ Control panel (left ~300 px)      │  Plots (right)                  │
│                                   │                                  │
│  Model: [Guinier            ▼]    │  ┌─────────────── 50% ────────┐ │
│                                   │  │  I(Q) log-log               │ │
│  Q range                          │  │  data (symbols) + model     │ │
│  From [0.003] To [0.08] Å⁻¹      │  │  (line) + complex BG        │ │
│  [Set Q from cursors]             │  │  (dashed)                   │ │
│                                   │  └─────────────────────────────┘ │
│  □ Complex background             │  ┌─────────────── 10% ────────┐ │
│  □ No limits                      │  │  Residuals r' (rescaled)    │ │
│                                   │  └─────────────────────────────┘ │
│  ─── Parameters ───               │  ┌─────────────── 40% ────────┐ │
│  I0   [1.0]  lo [0] hi [—]  [✓] │  │  Linearization              │ │
│  Rg   [50.]  lo[0.01] hi[500][✓] │  │  (Guinier/Porod families)   │ │
│                                   │  │  or "not available" label   │ │
│  N MC runs: [50]                  │  └─────────────────────────────┘ │
│                                   │                                  │
│  [Fit]  [Calculate Uncertainty]   │                                  │
│  [Store in File]  [Revert]        │                                  │
│  [Export parameters]              │                                  │
│  [Import parameters]              │                                  │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Models

### Guinier family

| Model | Formula | Linearization |
|-------|---------|--------------|
| **Guinier** | `I(Q) = I0·exp(−Q²Rg²/3)` | ln I vs Q² → slope = −Rg²/3 |
| **Guinier Rod** | `I(Q) = I0·exp(−Q²Rc²/2)/Q` | ln(QI) vs Q² → slope = −Rc²/2 |
| **Guinier Sheet** | `I(Q) = I0·exp(−Q²Rg²)/Q²` | ln(Q²I) vs Q² → slope = −Rg² |

Guinier Sheet also derives thickness `t = 2√3·Rg`.

Valid Q range for Guinier approximation: `Q·Rg < 1.3`.

### Porod / Power Law

| Model | Formula | Parameters |
|-------|---------|-----------|
| **Porod** | `I(Q) = Kp·Q⁻⁴ + Background` | Kp, Background |
| **Power Law** | `I(Q) = P·Q⁻ⁿ + Background` | P (prefactor), Exponent n, Background |

These models have an explicit flat Background parameter.  The complex background
option is not available for these two models.

### Debye Polymer Chain

| Model | Formula |
|-------|---------|
| **Debye Polymer Chain** | `I(Q) = Scale · 2(exp(−x) − 1 + x) / x²`  where `x = Q²·Rg²` |

The Debye function (P. Debye, 1947) describes single-chain scattering from an ideal
Gaussian polymer.  It is exact for a random-walk chain with no excluded volume.
Key properties:

- `P(Q) → 1` as `Q → 0` (forward-scattering limit)
- At high Q, `P(Q) ∝ Q⁻²` — characteristic of a chain with fractal dimension 2

Parameters: **Scale** (forward-scattering intensity at Q=0), **Rg** (radius of
gyration of the chain).  Complex background is supported.

### Sphere and Spheroid

| Model | Formula |
|-------|---------|
| **Sphere** | `I(Q) = Scale·|FF(QR)|²`;  `FF(x) = 3(sin x − x cos x)/x³` |
| **Spheroid** | orientational average of sphere FF with semi-axis ratio β (oblate β<1, prolate β>1) |

The spheroid orientational average is computed with Gauss-Legendre quadrature
(50-point rule) for efficiency.

Both models support complex background.  A flat background can be added via the
complex background option (set B = 0, flat = desired value).

### Correlation models and Invariant

| Model | Formula / description |
|-------|-----------------------|
| **Debye-Bueche** | `I(Q) = Prefactor·Eta²·ξ³/(1+Q²ξ²)²` — two-phase random media |
| **Treubner-Strey** | `I(Q) = Prefactor/(A + C1·Q² + C2·Q⁴)` — micro-emulsion / lamellar |
| **Benedetti-Ciccariello** | Correlated two-phase model with interfacial layer; derives from SLD contrast |
| **Hermans** | Paracrystalline lamellar model; complex exponentials for d1, d2, σ1, σ2 |
| **Hybrid Hermans** | Hermans + one Guinier term + one Unified Fit level |
| **Unified Born Green** | Two-level Unified model with Born-Green structure factor S(Q); gives correlation length ξ and inter-particle distance |

Treubner-Strey derives correlation length ξ = (A/C2)^¼/√2 and repeat distance d.
These appear in the "Derived" section of the results.

### Invariant (calculation — no fitting)

The **Invariant** entry in the model selector is not a fit: it directly
integrates the data to obtain the Porod invariant

```
Q* = ∫ q²·I(q) dq  =  2π² · Δρ² · φ(1−φ)
```

which is exact for any two-phase system with sharp interfaces, independent of
morphology.  This makes it a very robust way to obtain volume fractions for
polydisperse single-"particle"-type systems (porosity in rocks, precipitates,
voids) — no model, no fitting, just calibrated data and a contrast value.

Workflow:

1. Load data and select **Invariant** as the model.
2. If a background must be removed, check **Complex background** and use the
   **Fit B/P btwn cursors** / **Fit Flat btwn cursors** prefit buttons on a
   background-dominated Q region — the background `B·Q⁻ᴾ + flat` is
   subtracted from the data before integration.  The Q window used by each
   prefit button is **remembered** (shown below the Invariant options and
   stored with the setup).  Check **"Refit background from saved ranges"**
   to re-run these prefits automatically before every Calculate — in the
   GUI *and* in scripted/batch runs.  This is essential for scripting:
   the invariant is only as good as the background subtraction, so batch
   runs re-determine B/P/flat per file from the saved windows instead of
   trusting exported values.  The BG_P "Fit?" checkbox controls whether
   the power-law prefit fits both B and P or holds P and fits B only.
3. Enter **Contrast** Δρ² in units of 10²⁰ cm⁻⁴ (default 100; use the
   Scattering Contrast calculator if needed).
4. Position the cursors over the Q range to integrate and press
   **Calculate Invariant**.

The graph shows the background-corrected data (green), the background curve,
and the running integral `∫q²I dq` on the **right axis** (red) — analogous to
the integrated-volume-fraction curve in the Size Distribution tool.  The
invariant is read at the point where the running integral flattens (smoothed
derivative crosses zero); this Q value is reported as **QmaxUsed** and guards
against a slightly over-subtracted background at high Q.

Details and options:

- **Extrapolation to Q = 0** is always applied (the first integrand value is
  repeated at Q = 0, following Igor Irena).  Its contribution is ~I₁q₁³ —
  negligible for any reasonable q_min — so there is no user control for it.
- **Extend by Porod tail (Kp/Qmax)** (optional, off by default = Igor
  behaviour): adds the analytic high-Q tail `∫_{Qmax}^{∞} Kp·q⁻² dq = Kp/Qmax`
  to compensate for measurement truncation.  Kp is estimated as the median of
  `(I−bg)·q⁴` over the last half-decade before QmaxUsed.
- **Volume fraction**: with data on absolute scale [cm⁻¹] and correct
  contrast, φ is obtained by solving `φ(1−φ) = Q*/(2π²Δρ²)` (lower root,
  φ < 0.5 — the invariant cannot distinguish φ from 1−φ).  If
  `φ(1−φ) > 0.249` the value is unphysical: φ is set to NaN and a warning
  points at calibration/contrast problems.

Results stored in `derived/`: `Invariant` [cm⁻⁴], `VolumeFraction`,
`PhiOneMinusPhi`, `QmaxUsed` [Å⁻¹], and — when the tail option is on —
`PorodTail` [cm⁻⁴] and `PorodKp` [cm⁻¹Å⁻⁴].  The running integral is stored
as `Q_integral` / `running_integral` datasets and is exported to Igor
experiments as the `SimFitInvariantIntegral` wave.

**Caveats**: requires absolute intensity calibration and desmeared data
(slit-smeared USAXS must be desmeared first); q² amplifies background errors
at high Q, so background subtraction matters; systems with more than two
phases or diffuse interfaces bias the result.  χ², residuals, linearization
and Monte Carlo uncertainty do not apply and are disabled for this model.

---

## Background options

### Models with explicit Background (Porod, Power Law)

Background is a regular fitted parameter — enter a starting value and bounds like
any other parameter.

### Models with complex background (all others)

Check **"Complex background"** to add a background term of the form:

```
BG(Q) = B · Q^(−P) + flat
```

Three additional parameters appear in the grid, labelled **B**, **P**, **flat**
(stored internally as `BG_B`, `BG_P`, `BG_flat`):

| Parameter | Meaning | Typical use |
|-----------|---------|-------------|
| `B` | Power-law prefactor | Parasitic scattering |
| `P` | Power-law exponent | Usually 4 (Porod) |
| `flat` | Flat (incoherent) background | Incoherent scattering |

The naming mirrors the Unified Fit convention (B for prefactor, P for exponent).

To add only a flat background: set `B = 0` (fixed) and fit `flat`.
To model a Porod-law background: fix `P = 4` and fit `B`.

**Prefit helpers.**  When complex background is active, two small buttons appear
below the parameter grid:

- **Fit B/P btwn cursors** — fits the power-law `B·Q⁻ᴾ` to the data between the
  two cursors and writes the result as starting values.  If P's "Fit?" box is
  unchecked, only B is fit at the current (model-guided) P.
- **Fit Flat btwn cursors** — estimates the flat term from the median residual
  over the cursor window (using the current B and P).

Place the cursors over the low-Q power-law tail for **Fit B/P**, and over a
high-Q region where the data has flattened for **Fit Flat**.

---

## Q range selection

Two methods:

**Manual entry:** type Q min and Q max directly into the fields.

**Cursor-driven:** drag the two vertical cursor lines on the I(Q) plot to bracket
the region of interest, then click **"Set Q from cursors"**.  The Q min/max fields
update automatically.

Only data points inside the Q range are used for fitting.  The model is evaluated
and plotted over the full data Q range so you can see extrapolation quality.

---

## Fitting workflow

1. Select a **model** from the dropdown.
2. Set the **Q range** (see above).
3. Enter **starting parameter values** and optionally adjust lower/upper **limits**.
   - Uncheck "Fit?" to hold a parameter fixed at its current value.
   - Check "No limits" to ignore all bounds during fitting (uses very wide defaults).
4. Click **Fit**.  The model curve appears on the I(Q) plot; residuals update; the
   linearization panel updates (if applicable).  Each parameter field shows `± σ`
   (uncertainty from the covariance matrix).
5. If convergence is poor:
   - Try better starting values (the linearization plot can guide initial Rg/Rc estimates).
   - Narrow the Q range to exclude non-model contributions.
   - Try "No limits" if the optimiser is hitting a bound unexpectedly.

---

## Linearization plots

| Model | X axis | Y axis | Slope encodes |
|-------|--------|--------|--------------|
| Guinier | Q² | ln I | −Rg²/3 |
| Guinier Rod | Q² | ln(QI) | −Rc²/2 |
| Guinier Sheet | Q² | ln(Q²I) | −Rg² |
| Porod | Q⁴ | IQ⁴ | — (extrapolates to Kp at Q=0) |

The fitted Q-range points are highlighted; points outside the range are shown in grey.
A linear fit line is overlaid with the intercept labelled.

For all other models the linearization panel shows "No linearization available for
\<model\>".

---

## Monte Carlo uncertainty

The covariance-matrix uncertainties from `scipy.optimize.curve_fit` can underestimate
real parameter uncertainty when the model is non-linear or the data have correlated
errors.  The Monte Carlo option provides an independent estimate:

1. Set **N MC runs** (default 50; increase for publication-quality estimates).
2. Click **Calculate Uncertainty**.  The tool runs N independent fits on
   Gaussian-perturbed copies of the data (perturbation amplitude = dI).
3. The standard deviation of each parameter across the N fits updates the ± labels.

MC uncertainties reflect data noise propagation more faithfully than the covariance
matrix for models with strong parameter correlations (e.g., Spheroid Scale and R).

---

## Storing and loading results

### Store in File (HDF5 only)

Click **"Store in File"** to write results into `entry/simple_fit_results` in the
NXcanSAS HDF5 file.  The group contains:

```
entry/simple_fit_results/
    attrs:  model, success, chi_squared, reduced_chi_squared, dof,
            q_min, q_max, use_complex_bg, timestamp, program
    Q             [array, Å⁻¹]
    I_model       [array]
    residuals     [array]
    intensity_data  [array]   (copy of fitted I data)
    intensity_error [array]
    params/         one scalar per parameter
    params_std/     one scalar per parameter (MC or covariance uncertainty)
    derived/        model-specific computed values (xi, d, thickness, …)
```

### Loading stored results

```python
from pyirena.io.nxcansas_simple_fits import load_simple_fit_results
result = load_simple_fit_results("sample.h5")
print(result['model'], result['chi_squared'])
print(result['params'])        # dict of parameter values
print(result['params_std'])    # dict of uncertainties
print(result['derived'])       # model-specific computed values
```

Or via the top-level convenience function:

```python
from pyirena import load_result
result = load_result("sample.h5", "simple_fits")
```

---

## Export / Import Parameters

**Export Parameters** writes the current model configuration — model name, parameter
values, limits, Q range, complex-background flag — into a `pyirena_config.json` file
under the `simple_fits` key.  Multiple tool configurations can coexist in the same
JSON file.

**Import Parameters** reads the `simple_fits` section back and restores all values.

These files drive batch processing (see [Scripting API](#scripting-api)).

---

## Data Selector integration

With **"Simple Fits"** checked in the Data Selector:

| Button | Behaviour |
|--------|-----------|
| **Create Graph** | Opens a two-panel window (I(Q) + residuals) showing stored fit results from all selected HDF5 files.  Each file gets a unique colour; the model line uses a darker shade for contrast. |
| **Create Report** | Adds a "## Simple Fits" section to each file's Markdown report: model, chi², reduced chi², DOF, Q range, all parameters ± std, derived quantities. |
| **Tabulate Results** | Adds `SF_model`, `SF_chi2`, `SF_reduced_chi2`, `SF_dof`, `SF_q_min`, `SF_q_max`, `SF_use_complex_bg`, `SF_<param>`, `SF_<param>_std`, and `SF_derived_<name>` columns to the CSV table. |
| **Simple Fits (script)** | Reads `pyirena_config.json`, fits all selected files, saves results. |

---

## Scripting API

```python
from pyirena.batch import fit_simple, fit_pyirena
from pyirena import MODEL_NAMES

# List all available model names
print(MODEL_NAMES)

# Fit one file — minimal config
result = fit_simple("sample.h5",
                    config={'model': 'Guinier'})   # uses registry defaults

# Fit with explicit starting values and bounds
result = fit_simple("sample.h5", config={
    'model':  'Guinier',
    'params': {'I0': 5.0, 'Rg': 80.0},
    'limits': {'I0': [1e-3, 1e6], 'Rg': [1.0, 500.0]},
}, q_min=0.005, q_max=0.08)

if result and result['success']:
    print(f"Rg = {result['params']['Rg']:.2f} ± "
          f"{result['params_std']['Rg']:.2f} Å")
    print(f"Reduced chi² = {result['reduced_chi2']:.3f}")

# With Monte Carlo uncertainty
result = fit_simple("sample.h5",
                    config={'model': 'Sphere',
                            'params': {'Scale': 1e6, 'R': 100.0}},
                    with_uncertainty=True, n_mc_runs=100)

# Invariant (calculation — no fitting): volume fraction from Q*
result = fit_simple("sample.h5", config={
    'model': 'Invariant',
    'use_complex_bg': True,                  # subtract B·Q⁻ᴾ + flat first
    'params': {'Contrast': 100.0,            # Δρ² in 10²⁰ cm⁻⁴
               'BG_B': 0.0, 'BG_P': 4.0, 'BG_flat': 0.05},
    'invariant_porod_tail': False,           # True → add Kp/Qmax tail
    # Recommended for batch: refit the background per file from saved
    # Q windows (recorded by the GUI prefit buttons) before calculating.
    'bg_prefit': {
        'enabled': True,
        'power_law': {'use': True, 'q_min': 1e-4, 'q_max': 8e-4,
                      'fit_P': True},        # low-Q upturn window
        'flat':      {'use': True, 'q_min': 0.4, 'q_max': 1.0},  # high-Q
    },
}, q_min=0.001, q_max=0.3)
if result and result['success']:
    d = result['derived']
    print(f"Q* = {d['Invariant']:.4g} cm⁻⁴,  φ = {d['VolumeFraction']:.4g},"
          f"  QmaxUsed = {d['QmaxUsed']:.4g} Å⁻¹")

# Batch: use pyirena_config.json (exported from GUI)
from pathlib import Path
config = "pyirena_config.json"
for f in sorted(Path("data/").glob("*.h5")):
    r = fit_simple(f, config=__import__('json').load(open(config))['simple_fits'])
    if r and r['success']:
        print(f.name, r['params'])

# Or use fit_pyirena to run ALL tools at once
results = fit_pyirena("sample.h5", "pyirena_config.json")
sf = results['results'].get('simple_fits')
if sf and sf['success']:
    print(sf['params'])
```

For the full batch API reference, see [batch_api.md](batch_api.md).

## Slit smearing (USAXS)

When slit-smeared USAXS data are loaded (an NXcanSAS `dQl`, e.g. from Matilda),
a **"Slit smeared"** control row appears and each analytic model is smeared
before comparison, so the displayed "Graph model" curve matches the data
(e.g. Sphere/Spheroid oscillations are correctly smoothed) and fitted parameters
remain ideal-space. Guinier/Porod **linearization** has no closed form for
smeared data, so it is drawn best-effort (ideal-space) and labelled as such. The
**Invariant** integrates the data directly, which is invalid for smeared data,
so it is **disabled with a message** — load the desmeared dataset, or use Unified
Fit's model-based invariant. See **[Slit smearing](slit_smearing.md)** for the
full reference.

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

