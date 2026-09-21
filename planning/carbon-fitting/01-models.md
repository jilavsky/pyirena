# Models — extracted from Saurel et al. 2019/2020

All equation numbers below are the source paper's own numbering. "Main text
eq. (N)" refers to Saurel et al. 2019; "Annex 3 eq. (A3.N)" refers to the
Supplementary Information (`mmc2`); "Corrigendum eq." refers to the 2020
erratum.

## 0. Unit convention — DECIDED: Ångström everywhere

**The paper works in nm.** Q is in nm⁻¹, all lengths (R, Σ, ξ, d, r, w_P, w_C)
are in nm, and SLD/contrast are in 10¹⁰ cm⁻² / 10²⁰ cm⁻⁴ (their own SI
tables).

**pyIrena uses Å and Å⁻¹, and so does this tool — in the core math, in the
HDF5 file, in the JSON config and on the panel. nm appears nowhere.** Decided
2026-09-21 (Jan): same units for Q, sizes and contrast as the rest of
pyIrena; no nm/Å display split, no unit toggle. Every formula transcribed
from the paper below therefore needs 1 nm⁻¹ = 0.1 Å⁻¹ and 1 nm = 10 Å applied
to the paper's own tabulated values when comparing against them.

Units in force throughout:

| Quantity | Unit |
|---|---|
| Q | Å⁻¹ |
| lengths (r, R, Σ, ξ, d, w_P, w_C, R_rough) | Å |
| intensity | cm⁻¹ (absolute) |
| SLD | 10¹⁰ cm⁻² |
| contrast (Δρ)² | 10²⁰ cm⁻⁴ |
| specific surface area S | cm²/cm³ ≡ cm⁻¹ (reported *also* as m²/g) |
| ⟨δz²⟩ | Å² |
| density | g/cm³ |

Note on surface area: the paper's `S_macro`/`S_rough`/`S_mp` are per **gram**
(cm²/g) because their intensity is per gram. pyIrena's intensity is per
**volume** (cm⁻¹), so the fitted S values are cm²/cm³ and the m²/g number
comparable with BET is a derived quantity, `S[cm⁻¹] / ρ_sample[g/cm³] × 1e-4`.

## 1. Overall model — three additive components

Main text eq. (1):

```
I(Q) = I_Porod(Q) + I_mp(Q) + I_waxs(Q)
```

Three regions of Q dominate each term (main text Fig. 2, bottom panel):
I_Porod — macroscopic particle/grain surface, dominant at low Q, Q⁻⁴-like;
I_mp — micropore structure, dominant at intermediate Q; I_waxs — atomic-scale
turbostratic stacking, dominant at high Q (diffraction peaks). All three are
fit **simultaneously** over the full Q range — this is the point of the tool.

This maps directly onto the three tabs Jan described: **Complex background**
= I_Porod, **SAXS region** = I_mp, **WAXS region** = I_waxs.

---

## 2. Complex background — I_Porod (particle/grain morphology)

Main text eq. (3) (§3.2 "Particles morphological model"):

```
I_Porod(Q) = 2π(ΔSLD)² [ S_macro·Q⁻⁴  +  S_rough · f_rough(Q, R_rough) ]

f_rough(Q, R) = (2/9)·R⁴ / [ 1 + (1/5)(QR)² + (2/9)(QR)⁴ ]
```

- `S_macro` — macroscopic specific surface area of the powder grains
  (cm²/g if intensity is in cm²/g). Governs the Q⁻⁴ Porod slope visible at
  the very lowest Q.
- `S_rough` — additional specific surface area contributed by nanometer-scale
  surface roughness, only visible once Q ≳ 1/R_rough.
- `R_rough` — characteristic roughness length (nm in paper; convert to Å).
- `ΔSLD` — contrast between the powder grain and vacuum (`SLD_sample`,
  computed from structural density — see §5).

`f_rough` is the empirical "algebraic globule" envelope (Annex 3 eq.
A3.8/A3.9). **Coefficients verified against Jan's copy of the source PDF
(2026-09-21): the middle term is (1/5)(QR)², not (1/3)(QR)².** The earlier
planning pass mis-read it from a low-resolution page render; the value in the
box above is the paper's. Limits: `f_rough(0, R) = (2/9)R⁴`, and at
Q ≫ 1/R the leading term is `(2/9)R⁴ / [(2/9)(QR)⁴] = Q⁻⁴`, so the whole
bracket → `(S_macro + S_rough)·Q⁻⁴`, which is the paper's own eq. (4). Both
limits are asserted in `pyirena/tests/test_carbon_fit.py`.

**Generalization Jan asked for:** "low-Q power law slope + flat background"
implies the exponent on the macro term need not be fixed at −4 (real
interfaces are rarely perfectly sharp). Recommend keeping `S_macro·Q⁻⁴` as
the physically-motivated default but exposing the exponent as a fittable
parameter (defaulting to 4, fixed by default) — this is exactly the existing
`use_complex_bg` mechanism in Simple Fits (`Prefactor·Q⁻Exponent + Background`,
see 02-implementation-plan.md §2) and should be reused rather than
reimplemented.

**Derived:** `S_part = S_macro + S_rough` (total specific surface area,
directly comparable to BET N₂ adsorption SSA — this cross-check is the
paper's own validation method, Table 1 / Fig. 8b).

---

## 3. SAXS region — I_mp (micropore models)

Two alternative sub-models, selected per sample in the paper (a mode toggle
in the GUI):

### 3a. Dilute pores/particles + optional fractal aggregation

Main text eq. (5)–(6) (§3.3):

```
I_mp(Q) = I0 · ( 1 + 4π (D/(4πR^D)) Σ^D · [D·Γ(D−1) / (1+(QΣ)²)^((D−1)/2)]
                 · sin((D−1)·atan(QΣ)) / (QΣ) )
          · ( exp(−(Qr)²/5) + [erf(Qr/√10)]^12 · (9/2) · (k/(Qr)^4) )

I0 = φ·(ΔSLD)²·V_P¹      (paper's eq. 6; V_P¹ = (4/3)π r³, single-pore volume)
```

- First bracket is Teixeira's fractal-aggregate structure factor S(Q)
  (Annex 3 eq. A3.17) with fractal dimension `D` and cutoff length `Σ`
  (correlation length above which the fractal ordering is lost) — **this is
  already implemented** as `MassFractalPopulation`/`_mass_fractal_intensity`
  in `pyirena/core/modeling.py` (Teixeira 1988), see 02-implementation-plan.md.
  Set `D=3` (or `S≈1`) to disable aggregation and recover the dilute limit.
- Second bracket is the "algebraic globule" unified form factor for a pore of
  radius `r` (Annex 3 eq. A3.8/A3.9, Beaucage unified fit) — pore is
  isolated, aspect ratio ≈ 1. `k` is an adjustable prefactor (k=1 for
  monodisperse spheres, k>1 for globules of undefined/broader shape — see
  paper text after eq. 6). For slit-shaped pores, the paper notes the disk
  form factor (Annex 3 eq. A3.14, "unified discoid") can be substituted —
  **flag this as a v2 alternative pore shape**, not needed for v1.
- `φ` — pore volume fraction. `ΔSLD` — contrast between pore (SLD=0, vacuum)
  and carbon matrix.

**Derived** (main text eq. 13–15, Annex 3 eq. A3.20/A3.7):
```
φ  = I0 · ρ_struc / (8π(ΔSLD)²) · ξ⁻³ (1+(2πξ/d)²)²      [only for TS branch below — see 3b]
r  = sqrt(5·C1)                                            [algebraic-globule branch: r from Guinier term]
S_mp = 9k / (2r²) · I0 / (2π(ΔSLD)²)                       [specific surface area of micropores]
```
(Careful: the paper's φ/S_mp derived-quantity formulas are written for the
Teubner-Strey branch in the main text; for the fractal/globule branch the
natural derived quantities are simply `φ` and `S_mp` from `I0`, `r`, `ΔSLD`
directly via the Porod-consistency relation, `S_mp = 3φ/r` for spheres —
verify against Annex 3 A3.6/A3.7 when implementing, don't just transcribe
eq. 13-15 verbatim for this branch.)

### 3b. Teubner-Strey (semi-empirical two-phase / bi-continuous)

Main text eq. (8)–(12) (§3.3, used when pores show a broad correlation peak
rather than discrete particle-like scattering — e.g. GC, HC, AC, CDC900):

```
I_mp(Q) = I0 / (1 + C1·Q² + C2·Q⁴)

I0 = 8π·φ(1−φ)·(ΔSLD)² · ξ³ / (1+(2πξ/d)²)²
d  = 2π · [ ½C2^(−1/2) − C1/(4C2) ]^(−1/2)
ξ  = [ ½C2^(−1/2) + C1/(4C2) ]^(−1/2)
f_a = C1 / (2√C2)                       (disorder/"amphiphilicity" parameter)
```

- **Already implemented, verbatim**, as pyIrena's `Teubner-Strey` Simple
  Fits model: `I = Prefactor / (A + C1·Q² + C2·Q⁴)` with `Prefactor↔I0`,
  `A↔1`. It already computes `CorrLength` (↔ξ) and `RepeatDist` (↔d) as
  derived quantities (`core/simple_fits.py::_compute_derived`). **Missing**:
  `f_a`, and the materials-science-layer derived quantities below (φ, r,
  S_mp, w_P, w_C) which need `ΔSLD`/`ρ_struc` as extra inputs Simple Fits
  doesn't carry.
- `f_a` interpretation (Annex 3 A4.23, main text after eq. 12): `f_a<0` →
  short-range order, diffraction-peak-like; `f_a≈0.4` → equivalent to the
  globule form factor (3a); `f_a=1` → fully disordered, reduces to the
  Debye-Bueche model (also already in pyIrena's Simple Fits registry).

**Corrected derived quantities — corrigendum (2020), supersedes main-text
eq. (13)/(14):**
```
φ  from main text eq. (13), Annex 3 eq. (A3.20):
   φ = I0 · ρ_struc / (8π(ΔSLD)²) · ξ⁻³ · (1+(2πξ/d)²)²

w_P = ξ / (1 − φ)     (average pore width — CORRECTED formula, general, any f_a)
w_C = ξ / φ           (average carbon-matrix wall width — new in corrigendum)
S_mp = I0 / (C2 · 2π(ΔSLD)²)     (main text eq. 15, Annex 3 eq. A3.7 Porod limit)
```
The corrigendum explicitly says the **original** main-text method for
average pore radius (`r = sqrt(5·C1)`, eq. 14 — this is literally
`sqrt(5·C1)` reused from the algebraic-globule branch 3a) is **only accurate
when f_a≈0.4** (i.e., `C2 ≈ 2.8·C1²`) and in the low-pore-concentration
limit for spheroidal pores. `w_P`/`w_C` above are the general replacement,
derived from Babinet's principle applied to the Teubner-Strey correlation
length, and are valid for any `f_a`/φ. **Implement w_P/w_C as the primary
derived quantities; keep `r=sqrt(5·C1)` as a secondary "spheroid-limit
estimate" shown only when f_a is close to 0.4, with a tooltip explaining the
restriction.** This is exactly the kind of "convenience wiring" Jan flagged
as the hard/valuable part.

---

## 4. WAXS region — I_waxs (turbostratic diffraction peaks)

Main text eq. (2) (§3.1) is the general form; eq. (16)–(17) is the extension
used when layers are curved/crumpled (CPC, and generally any sp² carbon
showing peak broadening from curvature rather than finite crystallite size).

### 4a. Base model — one lamellar stacking peak

Main text eq. (2), Annex 3 §A3.6 (eq. A3.23–A3.43):

```
I_waxs(Q) = K · [P_1D(Q)/Q²] · [L(Q_c,w_L) ⊗ G(Q_c,w_G)] · exp(−Q²⟨δz²⟩/3)

Q_c = 2π/d          (peak center; d = layer spacing, e.g. d002)
```

- `K` — intensity scaling.
- `P_1D(Q)` — form factor of the layer/lamella (Annex 3 treats this as the
  atomic electron-density form factor of the graphene sheet for a real
  (002)-type peak — in practice, for a fitting tool, `P_1D` is usually
  absorbed into `K` as a slowly-varying prefactor near the peak, i.e. treat
  `K/Q²` as the fittable intensity-scale term rather than trying to model
  `P_1D` atomistically; **decide explicitly in implementation, don't silently
  drop the `/Q²`** — the `1/Q²` **orientation-averaging factor (powder
  average of a locally-1D/lamellar structure) is the physically important
  and easy-to-forget piece**, see the paper's extended discussion in §3.1
  and Annex 3 eq. A3.25.
- `L ⊗ G` — the peak shape is a **true Voigt profile** (convolution of
  Lorentzian width `w_L` and Gaussian width `w_G`), not a pseudo-Voigt
  linear mix. Physical origin: `w_L` (Lorentzian) comes from **distortions
  of the second kind** — layer bending/curvature (Annex 3 eq. A3.38–A3.41,
  Vonk); `w_G` (Gaussian) approximates **finite crystallite size** broadening
  (the true shape is squared-sinc, Annex 3 eq. A3.32, but a Gaussian is an
  excellent approximation to its upper 2/3 — Annex 3 Fig. A3.4 shows the
  ~1:2 Lorentzian:Gaussian FWHM ratio that works well across the whole
  SAXS+WAXS range, better than either pure shape alone). **pyIrena's existing
  `WAXSPeakFitModel`/`waxs_peakfit.py` only has Gauss/Lorentz/Pseudo-Voigt
  (linear mix) — a true Voigt (convolution) is new physics to add**, ideally
  as a fourth peak shape there so both tools share it (see
  02-implementation-plan.md §3).
- `exp(−Q²⟨δz²⟩/3)` — **distortions of the first kind** (Debye-Waller-like):
  random local fluctuations of interlayer spacing (or, in the general
  Debye-Waller sense, atomic thermal motion). Affects **intensity only**,
  not peak width or position (Annex 3, "Distortions of the first kind").
  New physics — not present anywhere in pyIrena yet.

**Practical note on peak count**: real carbon WAXS patterns show (002),
(100), and sometimes (004) — Fig. 2 inset in the main paper. The tool should
support an **arbitrary, user-managed list of peaks** (like `WAXSPeakFitModel`
already does for its own peak table) rather than hardcoding exactly one
Bragg reflection, even though the paper's worked examples mostly fit just the
(002).

### 4b. Extension for crumpled/curved layers — shared with SAXS region

Main text eq. (16)–(17) (§3.4, used for CPC soft carbon and generally
whenever layers are bent/crumpled rather than forming flat nanocrystallites):

```
I_waxs(Q) = K · [S_3D(Q)/Q²] · P_1D(Q) · L⊗G · exp(−Q²⟨δz²⟩/3)

S_3D(Q)/Q² = ( 1 + 4π(D/4πR^D)Σ^D · [DΓ(D−1)/(1+(QΣ)²)^((D-1)/2)] · sin((D−1)atan(QΣ))/(QΣ) )
             · ( exp(−(Qr)²/6) + [erf(1.06·QR/√12)]^6 · 2/(QR)² )   -- Annex 3 eq. A3.47 form
```

This is **structurally identical** to the SAXS-region fractal+globule model
in §3a above (same Teixeira fractal-SF term, same unified-fit-style
form-factor term for a "disk" of transition radius `R`, cutoff `Σ`, fractal
dimension `D`) — meaning **`R`, `Σ`, `D` are physically the same crumpling
geometry whether you're describing the low-Q pore/surface signature or the
high-Q peak broadening**. Table 2 in the paper (columns "R (nm)", "Σ (nm)",
"D") lists these per-sample for the WAXS fit; comparing to Table 1's Porod
roughness (`R_rough`) shows they're related but not forced identical in the
paper's own fits. **Recommend an optional "link geometry" checkbox that ties
the SAXS-region and WAXS-region R/Σ/D together as one set of fit parameters
when the user believes both signals come from the same crumpled-layer
geometry**, defaulting off (independent fits) — see
[03-open-questions.md](03-open-questions.md).

This crumpled-layer extension is the most complex single piece of physics in
the whole tool. **Recommend it be v2** (flagged in the implementation plan):
ship v1 with the simple Voigt-peak WAXS model (4a) plus the independent SAXS
fractal/globule or Teubner-Strey model, get that fitting correctly and
tested, then add the shared crumpled-geometry option as an enhancement.

---

## 5. Density / SLD / contrast conversion layer

This is the "wiring the conversions for density and other things" Jan
specifically called out as the hard/valuable part. Source: SI (`mmc3`)
eq. (S1)–(S2), **corrected** by the 2020 corrigendum.

```
ρ_struc = ρ_graphite · (d002_graphite/d002) · (d100/d100_graphite)²        [CORRECTED — corrigendum eq. S1]
          (original main-text/SI eq. S1 had this inverted — do not use the
          pre-corrigendum form)

ρ_sample = (1 − φ) · ρ_struc          [SI eq. S2; φ = micropore volume fraction from SAXS fit]

SLD_C      = f(ρ_struc, "C")          via NIST/Chantler tables — pyIrena already
                                       has this: pyirena.core.scattering_contrast
SLD_sample = f(ρ_sample, "C")         same function, different density input

ΔSLD (Porod, grain-vs-vacuum)   = SLD_sample − 0   (grain vs. surrounding vacuum/air)
ΔSLD (micropore, pore-vs-matrix) = SLD_struc  − 0   (pore is vacuum, contrast is against ρ_struc,
                                                       i.e. the carbon matrix's own SLD)
```

with `ρ_graphite = 2.26 g/cm³`, `d002_graphite = 0.3354 nm`,
`d100_graphite = 0.246 nm` (standard graphite reference values — confirm
exact digits against a citable source before hardcoding, e.g. Franklin 1951,
already a reference in the paper's bibliography).

`d002` and `d100` are read directly off the fitted WAXS peaks (§4): `d002 =
2π/Q_c` for the (002) peak, `d100 = 2π/Q_c` for the (100) peak. **This is the
cross-tab wiring**: the WAXS tab's fitted peak positions feed the density
calculation, which feeds the SAXS tab's contrast, which feeds `φ`, `r`,
`S_mp`, `w_P`, `w_C`. All of it should update live/on-demand rather than
requiring the user to re-type numbers between tabs.

**Reuse**: `pyirena/core/scattering_contrast.py` (`compute_compound`,
`compute_contrast`) already does formula→SLD/contrast for X-rays given a
density, using the NIST/Chantler anomalous-scattering tables. This is the
correct engine to call for `SLD_C`/`SLD_sample`/`ΔSLD` — **do not
reimplement SLD physics**; add a thin wrapper that takes `ρ_struc`/`ρ_sample`
and a fixed formula (default `"C"`, editable for e.g. N-doped or
Si-containing carbons) and returns the two SLDs and contrast(s) needed.

---

## 6. Full parameter inventory (v1 scope)

| Tab | Parameter | Symbol | Units (paper / pyIrena) | Source eq. |
|---|---|---|---|---|
| Background | Macro surface area | S_macro | cm²/g | eq. 3 |
| Background | Porod exponent | n (default 4, fixed) | — | generalization |
| Background | Roughness surface area | S_rough | cm²/g | eq. 3 |
| Background | Roughness length | R_rough | nm / Å | eq. 3 |
| SAXS (fractal branch) | Pore volume fraction | φ | — | eq. 6 |
| SAXS (fractal branch) | Pore radius | r | nm / Å | eq. 6 |
| SAXS (fractal branch) | Globule shape factor | k | — | eq. 6 |
| SAXS (fractal branch) | Fractal dimension | D | — | eq. 5 |
| SAXS (fractal branch) | Fractal cutoff length | Σ | nm / Å | eq. 5 |
| SAXS (Teubner-Strey branch) | Prefactor | I0 | cm²/g | eq. 8 |
| SAXS (TS branch) | Curvature coeff. | C1 | nm² / Å² | eq. 8 |
| SAXS (TS branch) | Curvature coeff. | C2 | nm⁴ / Å⁴ | eq. 8 |
| WAXS (per peak) | Peak center / d-spacing | Q_c (or d) | nm⁻¹ (nm) / Å⁻¹ (Å) | eq. 2 |
| WAXS (per peak) | Lorentzian FWHM | w_L | nm⁻¹ / Å⁻¹ | A3.38–41 |
| WAXS (per peak) | Gaussian FWHM | w_G | nm⁻¹ / Å⁻¹ | A3.32 |
| WAXS (per peak) | Intensity scale | K | — | eq. 2 |
| WAXS (per peak) | Mean-square disorder | ⟨δz²⟩ | nm² / Å² | eq. 2 |
| WAXS (crumpled, v2) | Transition radius | R | nm / Å | eq. 17 |
| WAXS (crumpled, v2) | Fractal cutoff | Σ | nm / Å | eq. 17 (shared w/ SAXS?) |
| WAXS (crumpled, v2) | Fractal dimension | D | — | eq. 17 (shared w/ SAXS?) |
| Material | Chemical formula | (default "C") | — | §5 |
| Material | Structural density | ρ_struc | g/cm³ | S1 (corrected) |

That's 15 core fit parameters in v1 (18 if the Teubner-Strey branch's own 3
are counted separately from the fractal branch's 5 — the two SAXS branches
are mutually exclusive, not additive), consistent with Jan's estimate of
"10-15 parameters" once background+peaks are included, before the v2
crumpled-layer extension adds 2-3 more (shared) geometric parameters.

**Derived quantities** (computed, not fit): `S_part` (=S_macro+S_rough),
`d` and `ξ` and `f_a` (TS branch), `φ`, `r` or `w_P`/`w_C`, `S_mp`,
`ρ_struc`, `ρ_sample`, `SLD_C`, `SLD_sample`, `ΔSLD` (×2, Porod and
micropore), and per-WAXS-peak `d-spacing` from `Q_c`. ~15 derived values —
this is the "convenience" layer.

---

## References

- D. Saurel, J. Segalini, M. Jauregui, A. Pendashteh, B. Daffos, P. Simon,
  M. Casas-Cabanas, "A SAXS outlook on disordered carbonaceous materials for
  electrochemical energy storage," *Energy Storage Materials* **21** (2019)
  162–173. `doi:10.1016/j.ensm.2019.05.007`
- Corrigendum, *Energy Storage Materials* **28** (2020) 418.
  `doi:10.1016/j.ensm.2020.03.013`
- Key cited methods reused inside the paper's own models: Beaucage (unified
  fit, `J. Appl. Crystallogr. 1995`), Teixeira (fractal SAS, `J. Appl.
  Crystallogr. 1988`), Vonk (distorted lamellae, `J. Appl. Crystallogr.
  1978`), Teubner & Strey (`J. Chem. Phys. 1987`), Schubert et al. (`J. Chem.
  Phys. 1994`, the `f_a` amphiphilicity factor), Debye & Bueche (`J. Appl.
  Phys. 1949`).
