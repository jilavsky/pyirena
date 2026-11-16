(size-distribution)=
# Size Distribution

:::{admonition} At a glance
:class: tip

**Launch:** `pyirena-gui` → **Models → Size Distribution**, or the **Size Distribution (GUI)**
button (it opens the first selected file). \
**Python:** `pyirena.api.control` (`open_dataset`, `select_sizes_model`, `run_sizes_fit`, …) or
`pyirena.batch.fit_sizes(data_file, config_file)`. \
**MCP:** `pyirena_ctrl_open_dataset`, then the category `"sizes"` tools through `pyirena_call`;
`pyirena_read_size_distribution` reads stored results. \
**Igor equivalent:** {doc}`Irena Size Distribution <igor:Irena/SizeDistribution>`
:::

The Size Distribution tool recovers a particle size distribution P(r) by one of four
inversion methods: MaxEnt, Tikhonov regularization, TNNLS/IPG and Monte Carlo. How the
methods work, and when to use which, is described in
[Size Distribution fitting methods](../methods/size-distribution-methods.md).

(size-distribution-window)=
## The window

```{figure} img/size-distribution/size-distribution-01-window.png
:alt: Size Distribution window after a regularization fit of sizes_sphere_lognormal.h5. Controls on the left in the Sizes tab; on the right the I versus Q graph with data, fitted model and cursors, the residuals, and the size distribution P(r) peaked near 100 Å. Numbered circles mark the parts described in the legend.

The Size Distribution window after a regularization fit of
`validationData/sizes/sizes_sphere_lognormal.h5` (contrast 100, grid 10–1000 Å).
```

1. **Data file**: the loaded file; **Open…** picks another (or drop a file on the panel).
2. **Slit smeared** and slit length **SL**: see [Slit smearing](#size-distribution-slit-smearing).
3. **Sizes / Background tabs**: the inversion setup, and the background terms.
4. **Particle Shape**: particle type (sphere, spheroid with aspect ratio) and contrast (Δρ)².
5. **Q Range (cursors)**: the inversion window set by the cursors, and the matching R max / R min.
6. **Size Grid**: r min, r max, number of bins, logarithmic spacing.
7. **Error Scaling**: error scale, or generated fractional errors.
8. **Fitting Method** and its parameters.
9. **Graph Model**, **Fit Sizes** and **Fit Cmplx. Bckg. & Sizes**.
10. **Results**: χ² (with its target), volume fraction Vf, Rg and peak r.
11. **Result buttons**: annotate the graph, revert, store, save and load setups, reports.
12. **Passes** and **Calc. Uncertainty (MC)**: Monte Carlo uncertainty.
13. **I(Q) plot**: data, model, background and the cursors that set the inversion window.
14. **Residuals plot** (rescaled, see [Fit quality](#size-distribution-fit-quality)).
15. **Distribution plot**: P(r) against r.
16. **Message bar**: the result of the last fit, with the fit-quality metrics.

(size-distribution-workflow)=
## Typical workflow

The steps use `validationData/sizes/sizes_sphere_flat_background.h5`: log-normal spheres
(R_median = 120 Å, volume fraction 0.006, contrast 100 ×10²⁰ cm⁻⁴) on a flat background of
0.02 cm⁻¹ that dominates above Q ≈ 0.15 Å⁻¹ (see `validationData/README.md`).

1. **Load data** and open Size Distribution from the main GUI.
2. **Set the particle shape and contrast** on the **Sizes** tab: sphere, contrast 100.
3. **Set the size grid**: r min 10 Å, r max 1000 Å, 200 bins, logarithmic spacing (or
   **Set r min/max from Q range**).
4. **Fit the background** on the **Background** tab: Q min 0.2, Q max 0.35 under **Flat
   Background (constant)**, then **Fit Background**. The corrected data I − (B·q⁻ᴾ + background)
   are drawn in blue.

   ```{figure} img/size-distribution/size-distribution-02-background-tab.png
   :alt: Background tab after Fit Background with Q min 0.2 and Q max 0.35; the flat background line and the background-corrected data in blue on the I versus Q graph.

   **Fit Background** at Q = 0.2–0.35 Å⁻¹ gives 0.0223 cm⁻¹; the corrected data fall away above
   Q ≈ 0.1 Å⁻¹.
   ```

5. **Choose the inversion window**: drag cursor B back to where the corrected signal still sits
   clearly above the background (here Q ≈ 0.15 Å⁻¹); see
   [Choosing the fit window](#size-distribution-fit-window).
6. **Choose the method** and click **Fit Sizes** (or **Fit Cmplx. Bckg. & Sizes** to refit the
   background terms first).
7. **Check the result**: flat residuals, and Vf, Rg and peak r in **Results**. **Results to
   graph** writes them on the I(Q) plot.

   ```{figure} img/size-distribution/size-distribution-03-fit-result.png
   :alt: Graph area after a regularization fit on Q 0.0008 to 0.15 with Results to graph: model through the corrected data, residuals, and a single peak in P(r) near 114 Å.

   Regularization on Q = 0.0008–0.15 Å⁻¹: Vf = 0.0060, peak r ≈ 114 Å.
   ```

8. **Estimate uncertainties** with **Calc. Uncertainty (MC)**, then **Store in File**.

% SCREENSHOT size-distribution-04-real-data: a regularization or MaxEnt fit on real USAXS data
% from testData/ (Jan picks the sample), with power-law and flat background, Results to graph.

(size-distribution-controls)=
## Controls

### Shared parameters (all methods)

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **r min [Å]** | Smallest particle radius modelled. Should be consistent with the maximum Q in your data: `R_min ≈ π/Q_max`. | 10 Å | **÷10 / ×10** scale the value |
| **r max [Å]** | Largest particle radius modelled. Should be consistent with the minimum Q: `R_max ≈ π/Q_min`. | 1000 Å | |
| **Set r min/max from Q range** | Set both from the cursor Q-range using r ≈ π/Q, rounded outward | | A too-tight grid piles the distribution against an edge |
| **Number of bins** | Number of radius bins. More bins give finer resolution but make the problem more underdetermined. 50–200 is typical. | 200 | |
| **Logarithmic spacing** | Logarithmically spaced bins give equal resolution per decade — strongly recommended for SAS data spanning multiple decades of r. | on | |
| **Particle type** | Form factor of the particles | sphere | **Aspect ratio** appears for spheroid (default 1.5) |
| **Contrast (Δρ²)** | Scattering length density contrast in units of 10²⁰ cm⁻⁴. Scales the absolute amplitude of P(r): volume fraction = ∫P(r)dr / contrast. Getting this right is essential for a physically meaningful volume fraction. | 1.0 | |
| **Error scale** | Multiplies all measurement errors before fitting. Values > 1 make the fit trust the data less (more regularisation); < 1 does the opposite. Use only if your error bars are known to be systematically under- or over-estimated. | 1.0 | Increase if MaxEnt cannot find a solution within the error bars |
| **Fractional error** / **Fraction** | Ignore the file's uncertainties and use error = \|I\| × fraction | off / 0.03 | When on, **Error scale** is ignored |

% CHECK(Jan): "volume fraction = ∫P(r)dr / contrast" (Contrast row) — validationData/README.md
% normalises P_V(r) so that ∫P_V(r)dr is the volume fraction itself.

### Method parameters

| Control | Method | Default | Details |
|---|---|---|---|
| **Method** | Inversion method | regularization | [Methods page](../methods/size-distribution-methods.md) |
| **Sky bg**, **Max iter** | MaxEnt | 1e-06, 1000 | {ref}`sizes-method-maxent` |
| **χ² (evalue)**, **Min ratio** | Regularization | 1.0, 1e-04 | {ref}`sizes-method-regularization` |
| **Approach**, **Max iter** | TNNLS | 0.95, 1000 | {ref}`sizes-method-tnnls` |
| **Convergence**, **Max iter** | Monte Carlo | 1.0, 100000 | {ref}`sizes-method-montecarlo` |

### Background (Background tab)

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Background** | Flat incoherent background [cm⁻¹] subtracted before fitting. | 0 | |
| **B** (power-law amplitude) | Amplitude of the low-Q power-law contribution B·Q⁻ᴾ. | 0 | **Fit B?** includes it in **Fit P/B** |
| **P** (power-law exponent) | Exponent of the low-Q power law (4 = Porod). | 4 | **Fit P?** includes it in **Fit P/B** |
| **Q min / Q max**, **Set Q from cursors** | Q range for **Fit P/B** and for **Fit Background** | | |
| **Fit P/B** | Fit B·Q⁻ᴾ to the data within its Q range; sets the starting values for the size fit | | |
| **Fit Background** | Fit a flat constant to the data within its Q range | | |

### Fitting

| Control | Meaning | Tips |
|---|---|---|
| **Graph Model** | Compute and display the current model without fitting | Preview before fitting |
| **Fit Sizes** | Fit the size distribution on the cursor Q range | |
| **Fit Cmplx. Bckg. & Sizes** | Fit the power law (if Fit B?/Fit P? is checked), then the flat background (if its Q range is set), then the sizes, in one click | |
| **Revert back** | Restore the parameters from before the last fit | |
| **Reset to Defaults** | Reset all parameters to factory defaults | Cannot be undone; **Revert back** undoes only the last fit |

(size-distribution-results)=
## Results and outputs

| Control | What it does |
|---|---|
| **χ², Vf, Rg, Peak r** | Fit result: χ² and its target, volume fraction, radius of gyration, peak radius |
| **Results to graph** | Annotate the I(Q) plot with Rg, peak r, volume fraction and χ², near Q ≈ 2π/Rg |
| **Store in File** | Save parameters and model curves into the source NXcanSAS file |
| **Load Setup from File…** | Restore every Sizes control (shape, r-range, method, bounds, …) from a stored file |
| **Save params to JSON / Load params from JSON** | Save the current parameters to a pyIrena JSON file, and restore them |
| **Copy results / Save report…** | The fit summary as Markdown, on the clipboard or in a `.md` file |
| **Passes / Calc. Uncertainty (MC)** | Repeat the fit on noise-perturbed data to estimate parameter uncertainties |

**Store in File** appends the results to the source file as `entry/sizes_results`, with the
complete panel setup embedded as `_pyirena_config`; see {ref}`hdf5-results-groups`. To reuse a
stored fit, **Load Setup from File…** restores the controls (also from files written by the
pyirena-ai agent), the [Data Selector](data-selector.md) graphs and tabulates results across
files, and `pyirena_read_size_distribution(path)` returns them to an AI agent.

(size-distribution-scripting)=
## From Python or AIDA

::::{tab-set}
:::{tab-item} Python
The control API drives the same model as the window. Run from the repository root:

```python
from pyirena.api.control import (
    open_dataset, select_sizes_model, set_shape, set_size_grid,
    fit_flat_background, set_fit_q_range, run_sizes_fit, save_sizes_fit,
)

sid = open_dataset("validationData/sizes/sizes_sphere_flat_background.h5")["session_id"]
select_sizes_model(sid, method="regularization")
set_shape(sid, shape="sphere", contrast=100)            # (Δρ)² in 10²⁰ cm⁻⁴
set_size_grid(sid, r_min=10, r_max=1000, n_bins=200, log_spacing=True)
fit_flat_background(sid, q_min=0.2, q_max=0.35)        # the flat high-Q plateau
set_fit_q_range(sid, q_min=0.0008, q_max=0.15)         # invert only above the background
fit = run_sizes_fit(sid)
print(fit["volume_fraction"], fit["rg"], fit["peak_r"])   # about 0.0060, 126, 114
save_sizes_fit(sid, output_path="sizes_flat_background_fit.h5")   # input file stays untouched
```

For many files with one setup, export the panel's setup with **Save params to JSON** and use
`pyirena.batch.fit_sizes(data_file, config_file)`; see [Batch API](../guides/batch-api.md).
:::
:::{tab-item} AIDA / MCP
Ask in plain language, for example *"Invert validationData/sizes/sizes_sphere_flat_background.h5
for a sphere size distribution, contrast 100, after fitting the flat background above
Q = 0.2."* The agent calls `pyirena_ctrl_open_dataset`, then `pyirena_call` with
`select_sizes_model`, `set_shape`, `set_size_grid`, `fit_flat_background`, `run_sizes_fit` and
`save_sizes_fit`. See the
[AI tools reference](../reference/ai-tools-reference.md#control-tools-reference--size-distribution-dispatcher-category-sizes).
:::
::::

(size-distribution-method)=
## Method notes

The four inversion methods, their parameters and when to use each are described in
[Size Distribution fitting methods](../methods/size-distribution-methods.md).

(size-distribution-fit-window)=
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

(size-distribution-log-axis)=
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

% CHECK(Jan): the distribution plot's axis is labelled "P(r) (a.u.)", while this section says
% it shows D(r) = dV/dr in volume fraction per Å.

(size-distribution-fit-quality)=
### Fit quality

Residuals are displayed **rescaled** as r' = r / σ(robust): the normalized
residual divided by a robust (MAD-based) estimate of the actual noise scale, so
the scatter is judged against the data's own noise floor rather than the (often
mis-scaled) reported σ. The fit status line also reports **σ-scale** (how many ×
the actual scatter exceeds the reported σ), the realistic reduced-χ² floor, and
the largest fractional misfit **max|(I−M)/I|** — a σ-independent gross-misfit
backstop.

See the **[Fit Quality Metrics guide](../methods/fit-quality-metrics.md)** for full
interpretation of these numbers.

(size-distribution-slit-smearing)=
### Slit smearing (USAXS)

When slit-smeared USAXS data are loaded (an NXcanSAS `dQl`, e.g. from Matilda),
the G matrix is smeared once (`G_smeared = W · G(q_ext)`), so every inversion
method above inherits smearing unchanged and the recovered distribution is
ideal-space (pinhole-equivalent). The power-law/complex background is smeared
too — the background preview shows the *smeared* curve so it overlays the data.
See **[Slit smearing](../methods/slit-smearing.md)** for the full reference and the
scripting contract.

## References

The references for each method are listed in
[Size Distribution fitting methods](../methods/size-distribution-methods.md#sizes-methods-references).
