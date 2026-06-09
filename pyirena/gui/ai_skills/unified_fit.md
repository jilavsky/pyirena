# Unified Fit — Expert Fitting Guidance

## Model structure
The Unified Fit model (Beaucage, 1995/1996) sums structural levels from
smallest to largest scale. Level 1 is the **small-scale** structure;
higher-numbered levels describe progressively larger features. Each level
has a Guinier term (G, Rg) and a power-law term (B, P).

## Parameter interpretation

**Rg** — radius of gyration in Å. The Guinier knee in the data (where log-log
slope begins to flatten) is roughly at Q ≈ 1/Rg. If the model knee is
displaced from the data knee, Rg is wrong.

**G** — Guinier prefactor. Sets the intensity of this level at Q → 0. If G
is too small, the level contributes negligibly; if too large, it dominates
inappropriately. For a single-level fit, G ≈ I(Q→0).

**P** — power-law exponent (slope in log-log space at high Q for this level).
- P = 4: smooth surface (Porod law)
- 3 < P < 4: rough surface (fractal surface)
- P = 3: mass fractal surface
- P < 3: mass fractal interior

Start with P = 4 fixed; free it only if residuals show a systematic slope
mismatch in the power-law region of the affected level.

**B** — power-law prefactor. Controls the absolute intensity of the power-law
tail. For smooth surfaces (P ≈ 4), B can be linked to G and Rg via the
Porod invariant: B ≈ (1.62/Rg)⁴ × G. Enabling estimate_B automates this.

**RgCO** (RgCutoff) — high-Q exponential roll-off radius. Physically, this is
the lower length-scale limit of the power-law regime — i.e., the primary
particle size when fitting an aggregate. Link RgCO to the next (smaller)
level's Rg when fitting hierarchical structures. Usually at most two levels
are linked together, reflecting primary particles (lower level) and aggregate
size (higher level). Optionally, for particles with two main dimensions, the
lower level represents the smaller dimension (for example, the radius of an
elongated cylinder) and the higher level represents its length. For disk-like
particles the lower level represents thickness (the smaller dimension) and
the higher level represents the disk radius (the larger dimension).

**ETA, PACK** — correlation length and packing factor for the Born-Green
liquid-like ordering correction. Only meaningful for concentrated,
interacting systems. Leave at zero unless there is a clear interference
peak.

## Staged fitting strategy
1. Fix everything except background; fit background to correct the baseline
   using high-Q data only, where the flat incoherent background dominates.
2. Free Rg and G for the smallest-features level (lowest level number); keep
   P fixed at 4 and fit to the Q range where the Guinier region of those
   features dominates. Enabling the "Estimate B from G, Rg, P" option reduces
   the number of free parameters and is recommended at this stage.
3. Evaluate residuals. If a systematic slope mismatch remains in the power-law
   region, free P for that level. If needed and physically appropriate, also
   free B.
4. Add a higher-numbered level only if residuals still show systematic
   structure at lower Q values (larger-scale features) than the current levels
   capture. Lower level numbers represent smaller features.
5. When adding a level, start with Rg ~ 5–10× the previous level's Rg
   (hierarchy principle). Free only Rg and G of the new level first.
6. Do not free B and Rg simultaneously without a good starting estimate for
   B; the correlation is strong and the fit may diverge.
7. The highest-numbered level (largest features) often has only a power-law
   slope because its Guinier region falls outside the measured Q range (below
   the minimum Q). In this case set G = 0 and Rg = 10¹⁰ for that level,
   leaving only the power-law slope. This is the normal way to fit a level
   when the data show a low-Q power-law slope without an obvious Guinier knee.

## Residual pattern recognition
- **Horizontal band around zero, random scatter** → good fit.
- **S-shaped (negative then positive with increasing Q)** → Rg is too large;
  reduce Rg or tighten its upper bound.
- **Inverted S-shape (positive then negative)** → Rg is too small.
- **Monotone slope in residuals** → P is wrong; free P.
- **Bump or dip near a specific Q** → a structural feature at that scale is
  not modelled; consider adding a level or enabling correlations.
- **High residuals only at very low Q** → beam-stop artifact or missing
  large-scale level; restrict Q range or add a Guinier level for the large-scale
  structure.
- **High residuals only at very high Q** → background is too low or there is
  an incoherent baseline; increase background or restrict Q range.
- **High-frequency noise on residuals** → the uncertainty estimates for the
  measured data are too small. This is a relatively common problem. Do not
  focus on the absolute value of χ²; look for systematic low-frequency
  variations instead. High-frequency residual structure cannot be described by
  the Unified Fit model and is either measurement noise, underestimated
  uncertainties, or diffraction peaks or other features in the data.

## Common mistakes
- **Fitting with too many free parameters at once** — parameters are strongly
  correlated. Always use staged fitting; freeing everything simultaneously
  rarely converges to a physically meaningful solution.
- **Rg of one level near Rg of another** — levels overlap in Q space and
  compete; the fit is degenerate. Ensure Rg values are separated by at least
  a factor of 3.
- **B too small or zero for a free-P level** — B and P are correlated; if B
  is near its lower bound, P will drift to compensate. Reset B to a
  physically reasonable value and refit.
- **G of a level fixed at zero** — that level contributes no Guinier term and
  is effectively a power-law-only level, which is physically unusual. Only
  do this intentionally for the highest-numbered level when its Guinier
  region is outside the measured Q range (see item 7 in Staged fitting).

## Q-range advice
Always check whether the full Q range is appropriate before fitting:
- Exclude low-Q points dominated by beam-stop shadow (they drive the
  background up artificially).
- Exclude high-Q points where the signal-to-noise ratio is poor
  (they pull P to unphysical values).
- State explicitly if you are advising the user to change the Q range.

## Reporting
When the fit converges, report: Rg and G for each level, P and whether it
is consistent with the expected morphology, reduced χ², and any caveats
about parameter correlations or parameters hitting their bounds.
