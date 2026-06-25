# Modeling — Expert Fitting Guidance

## Model structure
The Modeling tool builds the scattered intensity as a sum of independent
**populations** (P1–P10) plus a flat background. Each population is one of
several types — most commonly a **size-distribution population**: a
distribution of sizes (lognormal, Gaussian, Schulz-Zimm, etc.) convolved with
a **form factor** (sphere, spheroid, cylinder, core-shell sphere,
core-shell-shell sphere, …) and optionally multiplied by a **structure
factor** (hard-sphere, etc.). Other population types include Unified Level,
Guinier-Porod, mass/surface fractal, and diffraction peak.

Each population carries a `scale` (volume fraction term), a `contrast` (or
explicit SLDs for core-shell shapes), distribution parameters, and form-factor
parameters. Every numeric parameter has a "Fit?" checkbox and min/max limits.

## Parameter interpretation
- **mean_size / radius** — the characteristic size of the distribution [Å].
  The form-factor oscillations and the high-Q falloff are controlled by this.
- **width / sdeviation** — distribution breadth. Small widths give sharp
  form-factor (Bessel) oscillations; larger widths smear them into a smooth
  curve. Polydispersity strongly affects how easy the fit is (see below).
- **scale** — overall intensity of the population; spans many decades.
- **SLD core / shell / solvent** (core-shell, core-shell-shell) — scattering
  length densities in 10⁻⁶ Å⁻². Contrast steps between adjacent layers drive
  the relative amplitude of each shell's contribution.
- **t_shell / t_shell1 / t_shell2** — fixed shell thicknesses [Å]. For
  core-shell-shell the distribution is over the **core radius only**; both
  shell thicknesses are fixed parameters.
- **Background** — flat additive level [cm⁻¹]; usually fit.

## Choosing the fit method (Standard vs Global)

The Modeling panel has a **Fit method** selector (right of the Background
field): **Standard (local)** and **Global (DE→local)**.

- **Standard (local)** is the default — a fast gradient/Gauss-Newton search
  (Trust Region Reflective, or Nelder-Mead in "No limits" mode). It converges
  to the nearest minimum of the starting guess. Use it for routine, smooth,
  well-behaved models — especially **polydisperse** size distributions where
  the χ² surface is smooth.

- **Global (DE→local)** runs a differential-evolution (genetic/evolutionary)
  global search to find the right basin, then polishes locally. It is much
  slower but robust against local minima.

**When to recommend Global.** Monodisperse and near-monodisperse form factors
— particularly **core-shell** and **core-shell-shell spheres** — produce sharp
Bessel-function oscillations and a highly multimodal χ² surface: one local
minimum per oscillation lobe. A local fit started from the wrong radius gets
trapped in the wrong lobe and produces a visibly poor fit (model oscillations
out of phase with the data, large χ²/dof). The classic symptom is a Standard
fit that "sticks" — barely moving the radius and leaving an obviously bad
curve. In that situation, advise switching to Global.

**Workflow advice for Global fits:**
1. Set physically reasonable min/max limits on the fitted parameters. Global
   needs finite bounds — it is disabled in "No limits" mode.
2. Press **Fix limits?** to bracket each parameter around its current value
   (≈0.2×…5×). Tighter, sensible bounds make the global search converge much
   faster and more reliably.
3. Run Fit. Expect it to take noticeably longer than Standard; this is normal.
4. The global result is automatically polished with a local least-squares step,
   so reported parameters and uncertainties are at the true minimum.

If a Global fit is still poor, the most common causes are (a) limits too tight
and excluding the true value, (b) the wrong form factor / population type for
the data, or (c) too many simultaneously free parameters — fix the ones that
are physically known (e.g. SLDs, a known shell thickness) and re-run.

## General fitting strategy
- Start from a Graph Model preview to get the model into the right ballpark
  before fitting; a good starting guess helps even the global search.
- Fit the dominant population first; add and fit additional populations once
  the main feature is captured.
- Fix physically known quantities (SLDs, known thicknesses, P=4 for sharp
  interfaces) rather than letting everything float.
- Judge fit quality from χ²/dof and the residuals panel, not just the overlay.
