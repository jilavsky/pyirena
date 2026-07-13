# SAXS-Morph (3D two-phase solid): cross-implementation method comparison

**Purpose.** Evaluate whether the four implementations of the clipped-Gaussian-random-field
("leveled-wave") two-phase reconstruction implement the *same method*, and identify any
functional or approximation differences that are wrong or that constitute genuine improvements.
The concern is logic/method, not coding style.

**Implementations compared**

| Tag | What it is | Location |
|---|---|---|
| **Java** | Original SAXSMorph (V2.2), ~25 yr old, single-core | `java source/` (`Analyser.class`, `Morphology.class`, `MorphVoxel.class`) |
| **Fortran** | Roberts/Quintanilla reference GRF generator | `Fortran 2D two phase code Quintanilla/` (`gen3F.f`, `grfsubs.f`, `genco.f`, `dens_modele.f`) |
| **Igor** | Irena port | `SAXS_IgorCode/.../IR3_3DTwoPhaseSolid.ipf` |
| **Python** | pyIrena engine | `pyirena/core/saxs_morph.py` (+ `morphology.py`) |

The Java folder shipped only compiled `.class` files (no `.java`). No JDK/decompiler was
available in the sandbox, so the Java algorithm was recovered from the constant pool
(method names, field names, call order, referenced library calls). Data flow *inside* Java
methods is therefore inferred from method-call ordering + physics, not read line-by-line; the
one place this matters is flagged below.

---

## 1. The canonical method (reference)

From Quintanilla, Chen, Reidy & Allen, *Modelling Simul. Mater. Sci. Eng.* **15** (2007) S337
(the paper the Igor code cites step-by-step), plus Berk *Phys. Rev.* and Roberts *Phys. Rev. E*
**51** 4141 / **56** 3203 (PDFs in the Fortran folder). The data-driven "one-cut GRF" recipe is:

1. **χ(r)** — phase (Debye) autocorrelation from I(q):
   χ(r) = (1/2π²Vη²) ∫ [I(q)−I\*] q² sinc(qr) dq. Note χ(0) = φ(1−φ), χ(∞) → 0.
2. **α** — threshold for the target volume fraction φ: α = √2·erf⁻¹(1−2φ).
3. **G(r)** — the *field* autocorrelation of a zero-mean, unit-variance Gaussian field,
   obtained by **inverting the clipping (Berk) relation** χ(r) = T(G(r), α), where
   T(G,α) = (1/2π) ∫₀^G exp(−α²/(1+t)) / √(1−t²) dt. G is a *different* function from χ.
4. **G̃(k)** — the spectral density = FT of **G** (not χ):
   G̃(k) = (2/π) ∫ r² G(r) sinc(kr) dr, which must be ≥ 0.
5. **Field** — synthesize a Gaussian field with power spectrum G̃(k). Two equivalent routes:
   (a) FFT-filter white noise by **√G̃(k)**; (b) sum many equal-amplitude cosines
   √(2/N)·Σ cos(kₙ·r+φₙ) with directions isotropic, phases random, and magnitudes kₙ
   sampled ∝ G̃(k). Route (b) is the "leveled wave"; it is mathematically the √-weighting of (a).
6. **Threshold** the field at α → binary voxelgram of volume fraction φ.
7. **Model I(q)** back from the structure for comparison to data.

The two load-bearing subtleties are **(3) you must transform G, the *clipping-inverted* field
correlation — not χ itself** and **(5) the field is filtered by the *square root* of the spectral
density.** Both the Fortran reference and the Python port honor them; the Igor port violates both.

---

## 2. Step-by-step implementation matrix

| Step | Java (orig) | Fortran (ref) | Igor | Python (pyIrena) |
|---|---|---|---|---|
| 1. χ(r) from I(q) | `calcGammaOfR`, sinc transform | n/a (analytic G assumed) | `IR3T_ConvertIntToDACF` ✓ | `debye_autocorr` ✓ |
| 2. α = √2 erf⁻¹(1−2φ) | `erfinv` ✓ | from φ ✓ | `inverseErf` ✓ | `alfa_threshold` ✓ |
| 3. Berk invert χ→G | `calcGOfR` (`erfinv`/`calc_gr`) ✓ | n/a | computes `AutoCorfnctGr` ✓ **but not used for the field** | `berk_invert` ✓ **used** |
| 4. spectral density | FT of g(r) *(inferred)* | analytic `density(k)` | **FT of χ (`PhaseAutocorFnct`)** ✗ | **FT of g_r** ✓ |
| 5. field synthesis | cosine sum, k∝f(k), amp √(2/N) ✓ | FFT, a,b = **√σ**·dev ✓ | FFT, × `SpectralFk` — **no √** ✗ | FFT, × **√F** ✓ |
| 6. threshold at α·σ | `tval<α` ✓ | downstream ✓ | `greater(field, α·V_sdev)` ✓ | `field > α·σ` ✓ |
| 7. model I(q) | Debye back-transform | — | analytic from `TheorAutoCorrFnct` (**not** from the voxelgram) | **FFT of the actual voxelgram** ✓ |

Legend: ✓ matches canonical, ✗ deviates. "Inferred" = recovered from Java method ordering, not read directly.

---

## 3. Findings

### Finding 1 — Igor builds the field from the wrong correlation function (the main bug)

Igor **does** compute the Berk-inverted field autocorrelation G(r) (`AutoCorfnctGr`, step 3,
lines ~1256-1267) — correctly, using the same clipping integrand as everyone else. But at step 4
it builds the spectral density from the **phase** autocorrelation χ, not from G:

```igor
SpectralFk = IR3T_Formula4_SpectralFnct(Kvalues[p], PhaseAutocorFnct, Radii)   // line 1285
```

`PhaseAutocorFnct` is χ (the DACF). `AutoCorfnctGr` (= G) is computed and then used **only** for a
separate "theoretical intensity" display curve (lines ~1334-1347); it never reaches the generator.
So the entire Berk-inversion step is effectively dead code for the morphology. The generated field
therefore carries (approximately) the covariance of χ, and after thresholding its two-phase
correlation is T(χ,α) ≠ χ — i.e. the built structure does **not** reproduce the input χ(r).

This is precisely the systematic error the Python source warns about in its own comments
(`saxs_morph.py` lines ~1046-1056: "Without this step the model I(Q) is systematically too low by
a factor T′(0,α)… which can be 5-10×"). Python avoids it by transforming `g_r` (the berk_invert
output), matching the canonical recipe and the Java structure.

### Finding 2 — Igor omits the square root in the FFT generator (compounding bug)

Canonical generation filters white noise by **√(spectral density)** (Fortran `gencoeffs`:
`a = sqrt(sigma)·dev`; Python: `sqrt_F = np.sqrt(spectral_F)` then `noise_k * spectrum_3d`).
Igor's active generator multiplies the noise FFT by `SpectralFk` itself:

```igor
Gamma3D = cmplx(SpectralFkLoc[...], 0)               // line 1423  (no sqrt)
GaussNoise3DFFT = GaussNoise3DFFT * Gamma3D           // line 1425
```

so the field's power spectrum ∝ `SpectralFk²`. Combined with Finding 1, the Igor field covariance
≈ autocorrelation of χ (i.e. χ⋆χ) — doubly distorted from the target G. The dominant length scale
still survives, which is why Igor produces visually plausible structures that nonetheless
"give different values than SAXSMorph" (the author's own note). Notably, Igor's *other*,
commented-out generator (`IR3T_UseFFTtoGenerateMatrix`, line 1567) **did** take the square root
(`powC(...,1/2)`), with the author's margin note "even for 3D this needs to be sqrt, so why is it
for 1D not sqrt???" — the correct instinct sits right next to the shipped, incorrect path.

**Net:** the real methodological problems are in **Igor**, not pyIrena. Python corrects both
deviations and reproduces the Fortran reference recipe and the Java structure. On method grounds,
**pyIrena is the more faithful implementation of the published algorithm.**

### Finding 3 — Python verifies against the *actual* voxelgram (an improvement)

Python computes model I(q) by FFT of the generated binary cube (`voxelgram_to_iq`), so residuals
against data test what was actually built — discretization, finite box, and clipping included.
Igor instead computes a *theoretical* intensity from the analytic forward-clip of G, which can look
correct even when the voxelgram is wrong (as in Findings 1-2). For a stochastic method where
"you can't verify everything," measuring the realized structure is the better choice. The cost is a
noisier curve, which Python manages with log-binning + smoothing and a Porod tail (Finding 4).

### Finding 4 — Python's deliberate approximations (not bugs, but know they're there)

These are Python-specific regularizations with no counterpart in the reference; each trades a little
fidelity for stability:

- **Hann taper of F(k)** between 0.7·q_max and q_max, zeroing above q_max (lines ~1087-1099).
  Kills single-voxel speckle, but also removes real detail finer than ~π/q_max and slightly softens
  the sharpest interfaces / lowers the Porod amplitude.
- **`high_q_mode='porod'`** replaces the FFT tail above ~0.7·Q_nyq with an analytic K/Q⁴ (K = median
  of I·Q⁴). Clean and robust for a sharp two-phase interface; **wrong if the sample genuinely
  departs from Porod** (diffuse interface, surface-fractal) because it *imposes* Q⁻⁴.
- **Log-binned + Gaussian-smoothed** spherical average (~10% Q resolution). Good variance reduction,
  but can wash out sharp interference fringes — relevant because the current test data is a
  monodisperse-sphere profile whose ringing is exactly what smoothing attenuates.
- **Display smoothing kept out of the physics** (`smooth_sigma` applied only for the 3D view; I(q)
  uses the binary cube — lines ~1100-1110). This is *more* correct than Igor, which applies
  `ImageFilter/N=5 gauss3d` to the actual matrix (line 1446) before stats.
- **Fit voxel clamp** `MAX_FIT_VOXEL_SIZE=256`, render up to 512 — pure performance, no physics.

### Finding 5 — What is consistent everywhere (reassurance)

- α = √2·erf⁻¹(1−2φ) is identical in all four.
- The Berk clipping integrand exp(−α²/(1+t))/√(1−t²) is identical in Igor (`IR3T_JanCalcOfRInt`)
  and Python (`berk_lut`, after the t=sinθ substitution that removes the ±1 singularity). Same physics,
  and Python's LUT+interp is the same trick Igor now uses (`Fg_LUT`).
- Self-normalizing threshold at α·σ(field) (instead of enforcing G(0)=1) — Igor (`V_sdev`) and
  Python (`field.std()`) both do this; equivalent and robust.
- Which phase is labelled "1" and the threshold sign differ between codes but are immaterial (you
  choose which phase is which afterward) — consistent with the author's own notes.
- Java/Igor's cosine-sum path carries a "×10" units kludge (Igor line 2028, Å↔nm) that Python avoids
  by working in consistent physical k-units with explicit 2π factors — cleaner in Python.

---

## 4. Recommended actions

**To make Igor correct (and match pyIrena):** two one-line changes, then re-verify.

1. Step 4 — transform G, not χ:
   `SpectralFk = IR3T_Formula4_SpectralFnct(Kvalues[p], AutoCorfnctGr, Radii)` (line 1285).
   This finally makes the Berk inversion matter.
2. Step 5 — take the square root in `IR3T_MakeGRF`:
   `Gamma3D = cmplx(sqrt(SpectralFkLoc[...]), 0)` (line 1423) — matching the Fortran and the author's
   own commented note.
3. Verify: autocorrelate the thresholded matrix (`IR3T_Autocorelate3D` already exists) and confirm it
   reproduces T(G,α) = input χ. Before the fix it will not; after, it should.

**To validate pyIrena (a stochastic method, so verify the deterministic core against closed forms):**

1. **Analytic round-trip.** Feed a known analytic field correlation (e.g. Debye G(r)=exp(−r/a),
   whose I(q) ∝ 1/(1+(qa)²)²) straight into `spectral_function`→`generate_voxelgram`, bypassing the
   data→G step. Then (a) autocorrelate the voxelgram and confirm it returns T(G,α); (b) confirm
   `voxelgram_to_iq` returns the analytic Debye I(q). This isolates and proves steps 4-7 with no data.
2. **Cross-check against the Fortran.** Use the exact reference model and params from `GRinF.dat`
   (Roberts model (e): ξ=53.744, r_c=53.743, d=249.146; box T=1500, M=128) as the spectral density in
   Python's generator and compare field statistics, χ, and I(q). Cleanest apples-to-apples test of the
   generator core against the accepted reference.
3. **Sphere-profile sanity.** On the current sphere test data, temporarily disable the F(k) taper,
   the Porod extension, and the spectrum smoothing, and confirm the model-I(q) fringe *positions*
   line up with the data — i.e. that the smoothing hides no length-scale error.
4. **Invariant/contrast unit factor.** `derive_contrast_from_invariant` uses
   Q*=2π²φ(1−φ)Δρ² with a 1e4 unit factor the code says was wrong in pyirena ≤ 0.5.1. Run one dataset
   through Igor and Python and confirm Δρ² agrees numerically.

---

## 5. Bottom line

All four are the same clipped-GRF ("leveled-wave") method; the differences are in *how faithfully*
each executes steps 3-5. The Fortran is the reference generator; Java is structurally faithful
(cosine-sum sampling ∝ the spectrum of G). **The Igor port has two real, compounding bugs** — it
transforms χ instead of the Berk-inverted G (Finding 1) and drops the √ on the spectral density
(Finding 2) — so its morphology is only qualitatively right, matching the author's long-standing
"different from SAXSMorph" observations. **pyIrena fixes both and is the most faithful implementation
of the published algorithm**; its remaining caveats are deliberate regularizations (taper, Porod
tail, smoothing), not errors, but they should be disclosed when quoting fine-scale features or
non-Porod tails. Nothing found suggests pyIrena is producing physically wrong structures; the two
concrete corrections above belong in Igor.
