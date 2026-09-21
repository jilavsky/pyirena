# Carbon model — GUI Guide

The **Carbon model** tool fits the *whole* measured range of a disordered
carbonaceous material — USAXS through SAXS into WAXS, often five decades in Q
— as one model with three physically distinct contributions:

```
I(Q) = I_Porod(Q) + I_mp(Q) + I_waxs(Q) + background
```

| Contribution | What it is | Where it dominates |
|---|---|---|
| `I_Porod` | The outer surface of the powder grains, plus an optional nanometre-scale surface roughness | Lowest Q |
| `I_mp` | The micropore structure — either dilute pores with optional fractal aggregation, or a Teubner-Strey two-phase description | Middle of the range |
| `I_waxs` | Turbostratic stacking: Voigt diffraction peaks with a Debye-Waller factor and the powder-orientation average | Highest Q |

Everything is refined **together**. That is not a convenience: the contrast
that scales the background and the micropore term is computed from the fitted
WAXS peak positions and the porosity, so the regions are coupled through the
physics. There is deliberately no per-region fit button.

The tool's real value beyond the formulas is the **materials-science
convenience layer**. Enter a chemical formula once and the fitted peak
positions give the lattice spacings, the spacings give the structural density,
the porosity gives the sample density, and each density gives a scattering
length density and a contrast — feeding back into the SAXS and background
terms. Out the other side come about twenty quantities a carbon paper actually
quotes: BET-comparable specific surface areas, pore and wall widths, stack
height, layers per stack, d-spacings.

**Source**: D. Saurel, J. Segalini, M. Jauregui, A. Pendashteh, B. Daffos,
P. Simon, M. Casas-Cabanas, "A SAXS outlook on disordered carbonaceous
materials for electrochemical energy storage", *Energy Storage Materials*
**21** (2019) 162–173, `doi:10.1016/j.ensm.2019.05.007`, with the corrigendum
*Energy Storage Materials* **28** (2020) 418,
`doi:10.1016/j.ensm.2020.03.013`. Equation numbers below are the paper's own.

**Units**: Å and Å⁻¹ throughout, as everywhere else in pyIrena. The source
paper works in nm — its tabulated lengths are ×10 and its Q values ×0.1 before
they can be compared with anything here. Intensity is cm⁻¹, contrast
10²⁰ cm⁻⁴, surface area cm²/cm³ (the m²/g number is derived), ⟨δz²⟩ Å².

---

## Opening the tool

In the Data Browser:

1. Select one HDF5 data file.
2. Click **Carbon model (GUI)**, or use **Models → Carbon model**.

Shift-clicking the button forgets the window's remembered position and opens it
at its default size.

---

## Panel layout

A **left control panel** with five tabs, and a **right graph area**.

### Graph area

| Panel | Content |
|---|---|
| Top | Log-log I(Q): data, total fit (red), and the three components dashed — grain Porod (blue), micropores (orange), diffraction (green) |
| Middle | Weighted residuals, x-axis linked to the main plot |
| Bottom (optional) | The WAXS region on a **linear** Q *and* intensity axis |

Two vertical cursors set the fit range. Unlike the other pyIrena tools they
start on the data's own edges rather than 10 % inside: 10 % of five decades is
half a decade, and the top half-decade is where the (100) reflection lives.

**Show components** overlays the three contributions. Leave it on — seeing
which component owns which decade is most of what makes a full-range fit
debuggable, and it is how you catch the characteristic failure of this model,
the grain Porod term creeping up under the micropore region.

**WAXS zoom** adds the linear panel. A diffraction peak's *shape* is what the
two Voigt widths are fitted to, and that shape is unreadable on a log-log plot
spanning five decades.

**Auto-update** redraws the model whenever a control changes. Turn it off for
very large datasets.

Every plot carries the standard export menu (clipboard, PNG/JPEG/SVG, whole
window, curve CSV, Igor ITX).

### Parameter rows

Every fittable quantity, in every tab, is one row:

| Column | Meaning |
|---|---|
| **Fit?** | Refine this parameter |
| Parameter | Name and units |
| Value | Editable; **scroll over the field** to nudge it |
| lo / hi | Fitting bounds |
| ± std | 1-σ uncertainty after a fit |

A greyed-out row is a parameter the current model shape does not use — a
linked geometry parameter, or the branch of the SAXS-region selector that is
not chosen.

---

## Tab: Background

The grain Porod term, main text eq. (3):

```
I_Porod(Q) = 2π(Δρ)²·[ S_macro·Q⁻ⁿ + S_rough·f_rough(Q, R_rough) ]

f_rough(Q, R) = (2/9)·R⁴ / [ 1 + (1/5)(QR)² + (2/9)(QR)⁴ ]
```

| Parameter | Units | Notes |
|---|---|---|
| `S_macro` | cm²/cm³ | Specific surface area of the grain surface |
| `porod_exponent` | — | Fixed at 4 by default |
| `S_rough` | cm²/cm³ | Extra area from nanoscale roughness |
| `R_rough` | Å | Roughness correlation length |
| `flat_background` | cm⁻¹ | Instrumental background |

**Add nanoscale surface roughness** switches on the second term. It is
invisible below Q ≈ 1/R_rough and above it simply adds its surface area to the
Porod law, so at high Q the bracket becomes `(S_macro + S_rough)·Q⁻⁴` — the
paper's eq. (4). The derived **particle surface area** `S_part` is the sum, and
that is the number to compare with a BET nitrogen measurement.

The **Porod exponent** is exposed because real grain surfaces are rarely
perfectly sharp, but it is fixed at 4 by default and should usually stay there.
Away from exactly 4 the prefactor stops being a surface area, and the tool
reports every derived surface area as blank rather than as a wrong number.
Freeing the exponent and `S_rough` together is over-fitting.

---

## Tab: SAXS region

The micropores. Two **mutually exclusive** models — they are alternative
descriptions of the same pore population, and the paper never sums them.
Selecting one hides the other's parameters.

### Pores + fractal aggregation (eq. 5–6)

```
I_mp(Q) = I₀ · S_fractal(Q; D, Σ, r) · P_globule(Q; r, k)
I₀ = φ·(Δρ)²·(4/3)πr³
```

| Parameter | Units | Notes |
|---|---|---|
| `phi` | — | Micropore volume fraction φ. Also sets the sample density, and through it the grain contrast |
| `pore_radius` | Å | Pore radius r. Micropores in carbons are typically 3–10 Å |
| `globule_k` | — | Shape factor: 1 for monodisperse spheres, larger for pores of less well-defined shape |
| `fractal_D` | — | Mass fractal dimension of the aggregate |
| `fractal_sigma` | Å | Cutoff Σ above which the fractal ordering is lost |

`P_globule` is Beaucage's unified level for a compact object written out for a
sphere, and `S_fractal` is the Teixeira (1988) structure factor. Leave **Pores
aggregate into a mass fractal** off for dilute, uncorrelated pores; `D` and `Σ`
then disappear.

Derived: `S_mp = 3φk/r`, the micropore specific surface area.

### Teubner-Strey (eq. 8)

```
I_mp(Q) = I₀ / (1 + C₁·Q² + C₂·Q⁴)
```

Use this when the pores show a broad correlation peak rather than
particle-like scattering — the paper's glassy, hard, activated and CDC900
carbons. This is exactly the Simple Fits `Teubner-Strey` model, evaluated by
the same function.

Derived: the correlation length ξ, the repeat distance d, the amphiphilicity
`f_a = C₁/(2√C₂)`, the porosity φ inverted from I₀, and — following the 2020
corrigendum — the **average pore width** `w_P = ξ/(1−φ)` and **average wall
width** `w_C = ξ/φ`. These two replace the original main-text pore radius
`r = √(5·C₁)`, which the corrigendum restricts to `f_a ≈ 0.4`; that estimate is
still reported, but only when `f_a` is actually close to 0.4.

`f_a` reads as: negative → short-range order, peak-like; ≈ 0.4 → equivalent to
the globule model above; 1 → fully disordered (Debye-Bueche).

---

## Tab: WAXS region

```
I_waxs(Q) = envelope(Q) · exp(−Q²⟨δz²⟩/3) · Σ_peaks Voigt_i(Q) / Q²
```

Each peak is a **true Voigt** — the actual Lorentzian⊗Gaussian convolution, not
a linear pseudo-Voigt mix — because the two widths mean different things:

| Per-peak parameter | Units | Physical origin |
|---|---|---|
| `Q0` | Å⁻¹ | Peak centre; d = 2π/Q₀ |
| `K` | — | Amplitude, before the 1/Q² and Debye-Waller factors |
| `FWHM_G` | Å⁻¹ | Gaussian — finite crystallite size (Annex 3 eq. A3.32) |
| `FWHM_L` | Å⁻¹ | Lorentzian — layer bending and curvature, distortions of the second kind (Annex 3 eq. A3.38–41) |

The (002), (100) and (004) reflections are there by default at their graphite
positions, with (004) disabled because it is often too weak to see. **+ Add
peak** adds more; the **✕** removes one.

> The **label** matters. The Material tab looks peaks up by label, so the
> `002` and `100` names are what feed the density calculation. Renaming them
> breaks the chain — the Material tab then falls back to the manual density
> rather than reporting a wrong one.

**⟨δz²⟩ (stacking disorder)** is shared by every peak, not set per peak. It is
a property of the material — random fluctuations of the interlayer spacing,
distortions of the *first* kind — and it attenuates intensity without touching
peak width or position. Sharing it is what makes it identifiable: it sets the
(004)/(002) intensity ratio. Fitted per peak it would be degenerate with that
peak's amplitude.

**Apply 1/Q² orientation average** is the powder average of a locally
one-dimensional stacking correlation (Annex 3 eq. A3.25). Leave it on unless
you know why not: a fit without it looks perfectly fine and reports a wrong
amplitude, which is exactly why it is an explicit switch here rather than
something buried in the peak shape.

### Layer shape: crumpled layers (eq. 16–17)

For sp² carbons whose layers bend rather than forming flat nanocrystallites,
the summed peaks are multiplied by a crumpled-layer envelope — a Teixeira
fractal structure factor times a unified discoid form factor:

| Parameter | Units | Notes |
|---|---|---|
| `R_layer` | Å | Radius over which a layer stays flat before it bends |
| `fractal_D` | — | Crumpling fractal dimension |
| `fractal_sigma` | Å | Crumpling cutoff Σ |

This is the *same kind* of geometry the SAXS region may be describing, so each
of the three can be **linked** to its SAXS counterpart instead of fitted twice.
A linked parameter takes the SAXS value, leaves the fit vector, and greys out —
so the panel says the same thing the optimiser sees. Linking `R` to the pore
radius is a strong physical assumption (a disc radius is not a pore radius) and
is off by default.

---

## Tab: Material

The composition → density → contrast chain, and the place to break it when the
data cannot supply a stage.

```
ρ_struc  = ρ_graphite · (d₀₀₂,graphite/d₀₀₂) · (d₁₀₀,graphite/d₁₀₀)²
ρ_sample = (1 − φ) · ρ_struc
```

Two contrasts come out, and they are different:

- **Grain vs. vacuum** uses `ρ_sample` — the grain *including* its pores. The
  Background tab uses this one.
- **Pore vs. matrix** uses `ρ_struc` — the pore-free skeleton. The SAXS region
  uses this one.

The **Computed now** box shows the whole chain live, updating as the peaks move.

| Control | When to change it |
|---|---|
| Chemical formula | `C` by default; `C0.95N0.05` for an N-doped carbon, etc. Parsed by the Scattering Contrast engine |
| Structural density: *From fitted peaks* / *Manual* | Manual when there is no usable (100) peak |
| Porosity for density: *From SAXS region* / *Manual* | The Teubner-Strey branch derives φ from I₀, which depends on the very contrast it would feed — that branch needs a manual estimate here |
| Contrast: *From formula + density* / *Manual* | Manual bypasses the chain entirely |
| Graphite reference values | Only if you have a better reference |

> **On the graphite (100) reference.** The default is 2.1315 Å, the (100)
> d-spacing (a·√3/2 with a = 2.4612 Å). The source paper prints 0.246 nm, which
> is the *a* lattice parameter, not the d-spacing. Since the formula uses the
> ratio `d₁₀₀,graphite/d₁₀₀` and the denominator is read off a fitted peak as
> 2π/Q_c, the numerator has to be a d-spacing too — using 2.46 Å would scale
> the density wrongly by (√3/2)² ≈ 0.75. Keep this in mind when comparing
> against the paper's tables.
>
> Note also that both spacing ratios have the *graphite* value on top. Density
> goes as 1/(a²c), so a carbon whose lattice is swollen relative to graphite —
> which every disordered carbon is — is correspondingly *less* dense.

---

## Tab: Results

Every derived quantity, live, in the units a carbon paper quotes. The table is
copyable (Ctrl+C, right-click → Copy) and exports to CSV.

| Group | Quantities |
|---|---|
| Grain morphology | `S_part_m2_g` (BET-comparable), and the macroscopic / roughness split |
| Micropores, fractal branch | φ, pore radius r, I₀, `S_mp_m2_g` |
| Micropores, Teubner-Strey | ξ, d, f_a, w_P, w_C, `S_mp_m2_g` |
| Stacking | `L_c` (stack height), layers per stack, `L_a` (layer extent), ⟨δz²⟩ |
| Crumpling | D, Σ, R |
| Material | ρ_struc, ρ_sample, both SLDs, both contrasts, d₀₀₂, d₁₀₀ |
| Per peak | d-spacing, total Voigt FWHM, coherence length, peak height |

Quantities the current model shape does not define are omitted rather than
shown as N/A, so the table reflects what was actually fitted.

`L_c` comes from the Gaussian width alone (Scherrer, K = 0.9), because the
Lorentzian part is curvature, not finite size. `L_a` uses Warren's K = 1.84 for
a two-dimensional hk band.

---

## Fitting

| Button | Does |
|---|---|
| **Graph model** | Evaluate at the current parameters without fitting |
| **Fit all** | Refine every ticked parameter across the whole Q range at once |
| **Stop** | Abort; the parameters go back to where the fit started |
| **Calc. Uncertainty (MC)** | Re-fit noise-perturbed copies of the data |

Monte-Carlo uncertainties are more honest than the covariance estimate here,
because the parameters are always correlated — contrast couples the regions.

A practical order of work:

1. Load the data and press **Graph model**. The default parameters will be
   wrong, but the three dashed components show you which decade each one is
   trying to explain.
2. Get the WAXS peaks roughly in place first, with the **WAXS zoom** on. Their
   positions set the density and therefore the contrast for everything else.
3. Set a plausible porosity and pore radius so the middle of the range is in
   the right order of magnitude.
4. Fit everything. If it will not converge, untick the parameters you are least
   sure of, fit, then free them again.
5. Read the answers off the **Results** tab, not off the parameter values.

---

## Saving and scripting

| Button | Does |
|---|---|
| **Store in File** | Writes `entry/carbon_fit_results` into the loaded HDF5 file, with the full panel state embedded so the file reopens exactly as you left it |
| **Load Setup from File…** | Restores every control from a saved result file |
| **Save / Load params to JSON** | Writes a `carbon_fit` section into a shared `pyirena_config.json` for batch and scripted runs |
| **Copy report / Save report** | Markdown summary; Ctrl/⌘-click also writes the graph beside it |

### Batch

```python
from pyirena.batch import fit_carbon

result = fit_carbon("carbon_powder.h5", "pyirena_config.json")
if result["success"]:
    print(result["derived"]["S_part_m2_g"], "m2/g")
    print(result["derived"]["L_c"], "A stack height")
```

Or across a folder, through the shared pipeline:

```python
from pathlib import Path
from pyirena.batch import fit_pyirena

for path in sorted(Path("data").glob("*.h5")):
    fit_pyirena(path, "pyirena_config.json", tools=["carbon_fit"])
```

The config section *is* the model's `to_dict()`, so anything the panel can set
is scriptable with no separate translation layer.

### Reading results back

```python
import pyirena.api as api

res = api.read_carbon_fit("carbon_powder.h5")
print(res["derived"]["w_pore"], res["derived"]["w_carbon"])
```

### Driving it from an agent

The tool has a full control surface under the MCP category `carbon`.
Parameters are addressed by the model's own dotted keys rather than through one
setter per control:

```python
import pyirena.api.control as ctrl

sid = ctrl.open_dataset("/data/hard_carbon.h5")["session_id"]
ctrl.select_carbon_model(sid)
ctrl.configure_carbon_model(sid, saxs_mode="teubner_strey", use_roughness=True)
ctrl.list_carbon_parameters(sid)          # which keys exist now
ctrl.set_carbon_parameter(sid, "peak.002.Q0", 1.80)
ctrl.run_carbon_fit(sid)
ctrl.get_carbon_results(sid)              # derived quantities
ctrl.save_carbon_fit(sid)
```

Which keys exist depends on the model's shape, so `list_carbon_parameters`
after every `configure_carbon_model` is the discovery step, not an optional
one.

---

## Trend plots and Igor export

The HDF5 Data Explorer's **Collect** panel offers **Carbon model** as a type,
with both fitted parameters and derived quantities in one list — for a trend
across a series both are just "a number per file".

Igor export writes four waves per file: `CarbonModelI`, `CarbonPorodI`,
`CarbonMicroporeI` and `CarbonWAXSI`, all sharing a wave note that carries every
fitted parameter *and* every derived quantity.

---

## Notes and limitations

- The two SAXS-region models are alternatives, never summed. If neither
  describes the data, the sample may not be a two-phase microporous carbon.
- A free Porod exponent invalidates every derived surface area; the tool blanks
  them rather than reporting a number in the wrong units.
- The Teubner-Strey porosity is derived from I₀, which depends on the contrast,
  which depends on the porosity. The Material tab's manual porosity breaks that
  loop; there is no way to have it both ways.
- Validation is against analytic ground truth (synthesise from known
  parameters, refit, recover them) rather than against Igor Irena, which has no
  equivalent tool. See `pyirena/tests/test_carbon_fit.py`.
