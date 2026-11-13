# Size Distribution

The Size Distribution tool recovers a particle size distribution P(r) by one of four
inversion methods: MaxEnt, Tikhonov regularization, TNNLS/IPG and Monte Carlo. How the
methods work, and when to use which, is described in
[Size Distribution fitting methods](../methods/size-distribution-methods.md).

## Shared Parameters (all methods)

### Radius grid

| Parameter | Default | Meaning |
|-----------|---------|---------|
| **R min** | 10 Å | Smallest particle radius modelled. Should be consistent with the maximum Q in your data: `R_min ≈ π/Q_max`. |
| **R max** | 1000 Å | Largest particle radius modelled. Should be consistent with the minimum Q: `R_max ≈ π/Q_min`. |
| **N bins** | 200 | Number of radius bins. More bins give finer resolution but make the problem more underdetermined. 50–200 is typical. |
| **Log spacing** | Yes | Logarithmically spaced bins give equal resolution per decade — strongly recommended for SAS data spanning multiple decades of r. |

### Background

| Parameter | Default | Meaning |
|-----------|---------|---------|
| **Background** | 0 | Flat incoherent background [cm⁻¹] subtracted before fitting. |
| **B (power-law amplitude)** | 0 | Amplitude of the low-Q power-law contribution B·Q⁻ᴾ. |
| **P (power-law exponent)** | 4 | Exponent of the low-Q power law (4 = Porod). |

### Other

| Parameter | Default | Meaning |
|-----------|---------|---------|
| **Contrast (Δρ)²** | 1.0 | Scattering length density contrast in units of 10²⁰ cm⁻⁴. Scales the absolute amplitude of P(r): volume fraction = ∫P(r)dr / contrast. Getting this right is essential for a physically meaningful volume fraction. |
| **Error scale** | 1.0 | Multiplies all measurement errors before fitting. Values > 1 make the fit trust the data less (more regularisation); < 1 does the opposite. Use only if your error bars are known to be systematically under- or over-estimated. |

### Choosing the fit (inversion) window

The Q-range you invert over matters as much as the method. All four methods
assume every point in the window is *signal from the particles* (after the
flat + power-law background is subtracted). Two rules of thumb:

- **Low-Q end** should reach small enough Q to see the largest particles
  (`Q_min ≲ π/R_max`). Truncating it too high loses the big sizes.
- **High-Q end** should stop where the corrected signal still sits clearly
  **above the background**. USAXS/SAXS data typically fall into a flat,
  background-dominated, noisy plateau at high Q; the particle model cannot
  describe that region, so including it does not add information — it injects
  noise that the inversion tries to fit with spurious small particles.

Including background-dominated high-Q points is the most common cause of a
distribution that looks wrong despite a plausible χ². It shows up most sharply
with **Regularization**, where it can drive the χ² = M target out of reach and
trigger the fallback path (a warning is logged). The symptom is a distribution
that changes a lot when you move the high-Q cut-off by a little. If you see
that, **pull the high-Q end back** to where I(Q) departs visibly from the
background curve, and re-fit. Realistic error bars help too: overly tight
errors (e.g. a small fixed fractional error) at high Q make the discrepancy
principle unachievable and push every method toward over-fitting.

### Reading the distribution plot on a log radius axis

The distribution is plotted as **D(r) = dV/dr** (volume fraction per Å). When the
radius axis is **logarithmic**, be careful interpreting where the "bulk" of the
material is: equal *visual area* under D(r) does **not** correspond to equal
volume, because the linear bin width shrinks toward small r. A small amount of
volume in the smallest bins can therefore appear as a tall spike at r ≈ R_min
even when it represents a negligible fraction of the total — while the
physically dominant large-r volume looks compressed and small.

If you want the plot's visual weight to match actual volume on a log axis, read
it as **dV/dln(r) = r·D(r)** instead (volume per decade of radius). The reported
**Rg** and **volume fraction** are integrals of the true distribution and are not
affected by this plotting choice — so if a distribution *looks* concentrated at
tiny r yet reports a large Rg, that is the log-axis effect, not an inconsistency.
This is most visible with Monte Carlo, which (having no smoothness prior) is the
method most prone to leaving a little volume in the smallest bins.

---

## Fit Quality

Residuals are displayed **rescaled** as r' = r / σ(robust): the normalized
residual divided by a robust (MAD-based) estimate of the actual noise scale, so
the scatter is judged against the data's own noise floor rather than the (often
mis-scaled) reported σ. The fit status line also reports **σ-scale** (how many ×
the actual scatter exceeds the reported σ), the realistic reduced-χ² floor, and
the largest fractional misfit **max|(I−M)/I|** — a σ-independent gross-misfit
backstop.

See the **[Fit Quality Metrics guide](../methods/fit-quality-metrics.md)** for full
interpretation of these numbers.

## Slit smearing (USAXS)

When slit-smeared USAXS data are loaded (an NXcanSAS `dQl`, e.g. from Matilda),
the G matrix is smeared once (`G_smeared = W · G(q_ext)`), so every inversion
method above inherits smearing unchanged and the recovered distribution is
ideal-space (pinhole-equivalent). The power-law/complex background is smeared
too — the background preview shows the *smeared* curve so it overlays the data.
See **[Slit smearing](../methods/slit-smearing.md)** for the full reference and the
scripting contract.
