# Unified Fit — Expert Fitting Guidance

## Model structure
The Unified Fit model (Beaucage, 1995/1996) sums structural levels from
largest to smallest scale. Level 1 is the **large-scale** structure;
higher-numbered levels describe progressively smaller features. Each level
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
Start with P = 4 fixed; free it only if residuals show a systematic slope mismatch
in the power-law region of the affected level.

**B** — power-law prefactor. Controls the absolute intensity of the power-law
tail. For smooth surfaces (P ≈ 4), B can be linked to G and Rg via the
Porod invariant: B ≈ (1.62/Rg)^4 × G. Enabling estimate_B automates this.

**RgCO** (RgCutoff) — high-Q exponential roll-off radius. Physically, this is
the lower length-scale limit of the power-law regime — i.e., the primary
particle size when fitting an aggregate. Link RgCO to the next (smaller)
level's Rg when fitting hierarchical structures.

**ETA, PACK** — correlation length and packing factor for the Born-Green
liquid-like ordering correction. Only meaningful for concentrated,
interacting systems. Leave at zero unless there is a clear interference
peak.

## Staged fitting strategy
1. Fix everything except background; fit background to correct the baseline.
2. Free Rg and G for the dominant level; keep P fixed at 4.
3. Evaluate residuals. If systematic slope mismatch remains in the power-law
   region, free P for that level.
4. Add a second level only if residuals still show systematic structure at
   a different Q range than the first level captures.
5. When adding a level, start with Rg ~ 5–10× the previous level's Rg
   (hierarchy principle). Free only Rg and G of the new level first.
6. Do not free B and Rg simultaneously without a good starting estimate for
   B; the correlation is strong and the fit may diverge.

## Residual pattern recognition
- **Horizontal band around zero, random scatter** → good fit.
- **S-shaped (negative then positive with increasing Q)** → Rg is too large;
  reduce Rg or tighten its upper bound.
- **Inverted S-shape (positive then negative)** → Rg is too small.
- **Monotone slope in residuals** → P is wrong; free P.
- **Bump or dip near a specific Q** → a structural feature at that scale is
  not modelled; consider adding a level or enabling correlations.
- **High residuals only at very low Q** → background or beam-stop artifact;
  restrict Q range or add a large-scale (low-Q) Guinier level.
- **High residuals only at very high Q** → background is too low or there is
  an incoherent baseline; increase background or restrict Q range.

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
  do this intentionally for known power-law-only structures.

## Q-range advice
Always check whether the full Q range is appropriate before fitting:
- Exclude low-Q points dominated by beam-stop shadow (they drive the
  background up artificially).
- Exclude high-Q points where the signal-to-noise ratio is poor
  (they pull P to unphysical values).
- The residuals at excluded Q values will not be shown; state explicitly
  if you are advising to change the Q range.

## Reporting
When the fit converges, report: Rg and G for each level, P and whether it
is consistent with the expected morphology, reduced χ², and any caveats
about parameter correlations or boundary-hitting.
