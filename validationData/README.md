# Synthetic validation data for pyIrena and Irena

This folder holds **synthetic small-angle scattering data with exactly known
parameters**.  Its purpose is to make the correctness of pyIrena checkable, and
to let the same files be analysed with the Igor Pro *Irena* package so that the
two implementations can be compared parameter by parameter on identical input.

Every curve was produced by `generate_validation_data.py`, which computes the
scattering from independent implementations of the published model equations
(`_models.py`) and **does not import pyIrena**.  Recovering the tabulated
parameters is therefore a test of the analysis mathematics, not a
self-consistency check of a package inverting its own forward model.

## What is in each dataset

| file | content |
|---|---|
| `<name>.dat` | 3-column ASCII `Q  I  dI` with the noisy "measurement". The complete ground truth is repeated in the `#` header. |
| `<name>_ideal.dat` | 2-column ASCII `Q  I` — the exact, noise-free model. Use this to compare *forward models* between packages without any fitting. |
| `<name>.h5` | NXcanSAS HDF5 of the noisy data (`entry/sasdata`), with `entry/ground_truth/I_ideal` and `entry/ground_truth/parameters_json` alongside. |

`ground_truth.json` carries the same information machine-readably; it is what
`run_validation_report.py` reads when it fits every file and tabulates the
recovered parameters.

## Conventions

* `Q` in Å⁻¹, `I` in cm⁻¹ (absolute), radii and `Rg` in Å.
* `contrast` is (Δρ)² in units of 10²⁰ cm⁻⁴, matching Irena and pyIrena.
* Size distributions are quoted as **volume** distributions `P_V(r)`, normalised
  so that `∫P_V(r)dr` is the stated volume fraction.
* A log-normal `P_V(r)` is written with its **median** `R_median` and log-space
  standard deviation `sigma_log`; this is exactly Irena's `LogNormal` with
  `min_size = 0`, `mean_size = R_median`, `sdeviation = sigma_log`.

## Noise

`σ(Q) = sqrt(I/κ + (0.01·I)²)`, then `I_obs = I + σ·N(0,1)`, with `κ` chosen per
dataset so the counting term is 6 % at the 5th percentile of the intensity.
The random seed is `SHA-256(dataset name)[:8]`, so the files regenerate
byte-for-byte and adding a dataset never perturbs the others.

## Reproducing

```bash
python3 validationData/generate_validation_data.py    # regenerate the data
python3 validationData/run_validation_report.py       # fit it all, write the tables
```

## Tolerances

Recovered parameters are expected to agree with the ground truth to the
tolerances listed in `run_validation_report.py`.  Two systematic effects are
known and deliberately not hidden:

1. **Grid truncation.** Irena and pyIrena integrate a modelled size
   distribution over a finite radius grid that omits 1 % of the cumulative
   distribution at each tail, so a recovered volume fraction or `scale` may sit
   ~1–2 % below the generated value. The generator integrates over ±6σ.
2. **Regularised inversions.** The Size Distribution tool solves an
   ill-conditioned inverse problem; the *shape* it returns is smoothed relative
   to the generating distribution. Integral quantities (volume fraction, mean
   and RMS radius) are the meaningful comparison, not bin-by-bin values.

---

## Datasets

| # | dataset | tool | model | Q range (1/Å) | points |
|---|---|---|---|---|---|
| 1 | `sizes/sizes_sphere_lognormal` | Size Distribution | Single log-normal volume distribution of solid spheres | 0.0008 – 0.35 | 350 |
| 2 | `sizes/sizes_sphere_bimodal` | Size Distribution | Bimodal log-normal volume distribution of solid spheres | 0.0005 – 0.4 | 400 |
| 3 | `sizes/sizes_sphere_flat_background` | Size Distribution | Log-normal spheres on a flat background | 0.0008 – 0.35 | 350 |
| 4 | `sizes/sizes_spheroid_lognormal` | Size Distribution | Log-normal volume distribution of oblate spheroids (AR = 0.3) | 0.0005 – 0.3 | 320 |
| 5 | `unified/unified_one_level` | Unified Fit | Single Beaucage Unified level (Guinier + Porod), flat background | 0.0005 – 0.3 | 350 |
| 6 | `unified/unified_two_level` | Unified Fit | Two-level Beaucage Unified model with level-2 low-Q cut-off | 0.0001 – 0.3 | 450 |
| 7 | `unified/unified_one_level_slitsmeared` | Unified Fit (slit smearing) | Single Unified level, infinite-slit smeared | 0.0005 – 0.3 | 350 |
| 8 | `modeling/modeling_sizedist_plus_unified` | Modeling | Population 1: log-normal spheres.  Population 2: Unified level. | 5e-05 – 0.5 | 480 |
| 9 | `modeling/modeling_sizedist_plus_peak` | Modeling | Population 1: log-normal spheres.  Population 2: Gaussian diffraction peak. | 0.003 – 0.5 | 450 |
| 10 | `modeling/modeling_hardsphere` | Modeling | Log-normal spheres with a Percus-Yevick hard-sphere structure factor | 0.003 – 0.5 | 450 |
| 11 | `modeling/modeling_mass_fractal` | Modeling | Mass-fractal aggregate of spherical primary particles (Teixeira 1988) | 0.0003 – 0.3 | 400 |
| 12 | `modeling/modeling_surface_fractal` | Modeling | Surface-fractal scattering (Teixeira 1988) | 0.0005 – 0.3 | 380 |
| 13 | `modeling/modeling_guinier_porod` | Modeling | Guinier-Porod model, single level with dimensionality s1 = 1 (rod-like) | 0.002 – 0.4 | 400 |
| 14 | `simple_fits/simple_guinier` | Simple Fits | Guinier | 0.002 – 0.03 | 200 |
| 15 | `simple_fits/simple_guinier_rod` | Simple Fits | Guinier Rod | 0.005 – 0.08 | 200 |
| 16 | `simple_fits/simple_guinier_sheet` | Simple Fits | Guinier Sheet | 0.01 – 0.1 | 200 |
| 17 | `simple_fits/simple_porod` | Simple Fits | Porod | 0.05 – 0.5 | 200 |
| 18 | `simple_fits/simple_power_law` | Simple Fits | Power Law | 0.01 – 0.4 | 250 |
| 19 | `simple_fits/simple_sphere` | Simple Fits | Sphere (monodisperse) + flat background | 0.002 – 0.15 | 500 |
| 20 | `simple_fits/simple_debye_chain` | Simple Fits | Debye Polymer Chain | 0.003 – 0.3 | 300 |
| 21 | `simple_fits/simple_debye_bueche` | Simple Fits | Debye-Bueche | 0.002 – 0.2 | 300 |
| 22 | `waxs/waxs_three_gauss` | WAXS Peak Fit | Three Gaussian peaks on a linear background | 1.2 – 4 | 700 |
| 23 | `waxs/waxs_pseudovoigt` | WAXS Peak Fit | Two pseudo-Voigt peaks on a cubic-polynomial background | 1.5 – 4 | 650 |
| 24 | `merge/merge_usaxs` | Data Merge | Low-Q (USAXS-like) branch of a two-level Unified curve, correct scale | 0.0001 – 0.05 | 300 |
| 25 | `merge/merge_saxs` | Data Merge | High-Q (SAXS-like) branch of the SAME curve, deliberately mis-scaled by 1.15 | 0.005 – 0.5 | 320 |
| 26 | `manipulation/manipulate_reference` | Data Manipulation | Reference curve (single Unified level) | 0.001 – 0.3 | 300 |
| 27 | `manipulation/manipulate_input` | Data Manipulation | The reference curve scaled and offset by known constants | 0.001 – 0.3 | 300 |

## Ground truth, dataset by dataset

### 1. `sizes_sphere_lognormal`

**Tool:** Size Distribution  
**Model:** Single log-normal volume distribution of solid spheres  
**Equation:** `I(Q) = 1e-4 * contrast * INT V(r) F_sph^2(Qr) P_V(r) dr`  
**Q:** 0.0008 – 0.35 Å⁻¹, 350 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `distribution` | log-normal (volume weighted) | - |  |
| `R_median` | 100 | A | median of P_V(r) |
| `sigma_log` | 0.25 | - | log-space standard deviation |
| `volume_fraction` | 0.01 | - | INT P_V(r) dr |
| `contrast` | 100 | 1e20 cm^-4 | (delta-rho)^2 |
| `background` | 0 | 1/cm |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `R_mode` | 93.9413 | A | median*exp(-sigma^2) |
| `R_mean_vol` | 103.174 | A | median*exp(sigma^2/2) |
| `R_rms` | 106.449 | A | reported by pyIrena as 'Rg' |

Noise: κ = 63358.3, seed = 3027962536, median relative error 1.00 %, max 11.08 %.

### 2. `sizes_sphere_bimodal`

**Tool:** Size Distribution  
**Model:** Bimodal log-normal volume distribution of solid spheres  
**Equation:** `I(Q) = 1e-4 * contrast * INT V(r) F_sph^2(Qr) [P_1(r)+P_2(r)] dr`  
**Q:** 0.0005 – 0.4 Å⁻¹, 400 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `R_median_1` | 40 | A |  |
| `sigma_log_1` | 0.18 | - |  |
| `volume_fraction_1` | 0.004 | - |  |
| `R_median_2` | 250 | A |  |
| `sigma_log_2` | 0.22 | - |  |
| `volume_fraction_2` | 0.01 | - |  |
| `contrast` | 100 | 1e20 cm^-4 |  |
| `background` | 0 | 1/cm |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `R_mode_1` | 38.7248 | A |  |
| `R_mode_2` | 238.188 | A |  |
| `volume_fraction_total` | 0.014 | - |  |

Noise: κ = 68857.2, seed = 3407787525, median relative error 1.00 %, max 11.78 %.

### 3. `sizes_sphere_flat_background`

**Tool:** Size Distribution  
**Model:** Log-normal spheres on a flat background  
**Equation:** `I(Q) = 1e-4*contrast*INT V(r) F_sph^2(Qr) P_V(r) dr + flat`  
**Q:** 0.0008 – 0.35 Å⁻¹, 350 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `R_median` | 120 | A |  |
| `sigma_log` | 0.22 | - |  |
| `volume_fraction` | 0.006 | - |  |
| `contrast` | 100 | 1e20 cm^-4 |  |
| `background_flat` | 0.02 | 1/cm | dominates above Q ~ 0.15 1/A; fit it there before inverting |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `R_mode` | 114.33 | A |  |
| `R_mean_vol` | 122.939 | A |  |
| `R_rms` | 125.951 | A |  |

Noise: κ = 12526.5, seed = 1916699605, median relative error 1.00 %, max 6.30 %.

### 4. `sizes_spheroid_lognormal`

**Tool:** Size Distribution  
**Model:** Log-normal volume distribution of oblate spheroids (AR = 0.3)  
**Equation:** `I(Q) = 1e-4*contrast*INT V(r,AR) <F_ell^2(Q,r,AR)>_orient P_V(r) dr`  
**Q:** 0.0005 – 0.3 Å⁻¹, 320 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `shape` | spheroid | - |  |
| `aspect_ratio` | 0.3 | - | semi-axes (r, r, AR*r) |
| `R_median` | 150 | A |  |
| `sigma_log` | 0.2 | - |  |
| `volume_fraction` | 0.008 | - |  |
| `contrast` | 50 | 1e20 cm^-4 |  |
| `background` | 0 | 1/cm |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `R_mode` | 144.118 | A |  |
| `R_rms` | 156.122 | A |  |

Noise: κ = 61853.2, seed = 312993427, median relative error 1.00 %, max 11.44 %.

### 5. `unified_one_level`

**Tool:** Unified Fit  
**Model:** Single Beaucage Unified level (Guinier + Porod), flat background  
**Equation:** `I(Q) = G exp(-Q^2 Rg^2/3) + B/Q*^P + bg,  Q* = Q/erf(K Q Rg/sqrt6)^3,  K = 1 for P>3`  
**Q:** 0.0005 – 0.3 Å⁻¹, 350 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `G` | 100 | 1/cm | Guinier prefactor |
| `Rg` | 200 | A |  |
| `P` | 4 | - | power-law slope |
| `B` | 3.045e-07 | cm^-1 A^-P | = G exp(-P/2)(3P/2)^(P/2)/Rg^P (smooth join) |
| `background` | 0.01 | 1/cm |  |
| `K` | 1 | - |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `I(0)` | 100.01 | 1/cm |  |
| `Q_rollover` | 0.0141421 | 1/A |  |

Noise: κ = 27407.2, seed = 1260015518, median relative error 1.01 %, max 6.11 %.

### 6. `unified_two_level`

**Tool:** Unified Fit  
**Model:** Two-level Beaucage Unified model with level-2 low-Q cut-off  
**Equation:** `I(Q) = SUM_i [G_i exp(-Q^2 Rg_i^2/3) + B_i/Q*_i^P_i exp(-Q^2 RgCO_i^2/3)] + bg`  
**Q:** 0.0001 – 0.3 Å⁻¹, 450 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `G_1` | 8 | 1/cm |  |
| `Rg_1` | 120 | A |  |
| `P_1` | 4 | - |  |
| `B_1` | 1.88e-07 | cm^-1 A^-P |  |
| `RgCO_1` | 0 | A |  |
| `G_2` | 4000 | 1/cm |  |
| `Rg_2` | 1200 | A |  |
| `P_2` | 3.2 | - |  |
| `B_2` | 1.393e-06 | cm^-1 A^-P |  |
| `RgCO_2` | 120 | A | linked to Rg of level 1 |
| `background` | 0.02 | 1/cm |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `I(0)` | 4008.02 | 1/cm |  |

Noise: κ = 13809.4, seed = 286680045, median relative error 1.01 %, max 6.10 %.

### 7. `unified_one_level_slitsmeared`

**Tool:** Unified Fit (slit smearing)  
**Model:** Single Unified level, infinite-slit smeared  
**Equation:** `I_smr(Q) = (1/SL) INT_0^SL I_ideal(sqrt(Q^2+l^2)) dl`  
**Q:** 0.0005 – 0.3 Å⁻¹, 350 points, log-spaced  
**Slit length dQl:** 0.03 Å⁻¹  

| parameter | value | unit | note |
|---|---|---|---|
| `G` | 100 | 1/cm |  |
| `Rg` | 200 | A |  |
| `P` | 4 | - |  |
| `B` | 3.045e-07 | cm^-1 A^-P |  |
| `background` | 0.01 | 1/cm |  |
| `slit_length_dQl` | 0.03 | 1/A | written as entry/sasdata/dQl in the HDF5 |

Noise: κ = 27411.7, seed = 2570568956, median relative error 1.04 %, max 6.11 %.

### 8. `modeling_sizedist_plus_unified`

**Tool:** Modeling  
**Model:** Population 1: log-normal spheres.  Population 2: Unified level.  
**Equation:** `I(Q) = I_sizedist(Q) + I_unified(Q) + bg`  
**Q:** 5e-05 – 0.5 Å⁻¹, 480 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `pop1_type` | size_dist / lognormal / sphere | - |  |
| `pop1_min_size` | 0 | A |  |
| `pop1_mean_size` | 30 | A | log-normal median |
| `pop1_sdeviation` | 0.3 | - |  |
| `pop1_scale` | 0.005 | - | INT P_V(r) dr |
| `pop1_contrast` | 100 | 1e20 cm^-4 |  |
| `pop2_type` | unified_level | - |  |
| `pop2_G` | 20000 | 1/cm |  |
| `pop2_Rg` | 4000 | A |  |
| `pop2_P` | 3.5 | - |  |
| `pop2_B` | 1.56345e-08 | cm^-1 A^-P |  |
| `background` | 0.01 | 1/cm |  |

Noise: κ = 20749.9, seed = 1656842697, median relative error 1.02 %, max 6.84 %.

### 9. `modeling_sizedist_plus_peak`

**Tool:** Modeling  
**Model:** Population 1: log-normal spheres.  Population 2: Gaussian diffraction peak.  
**Equation:** `I(Q) = I_sizedist(Q) + A exp(-(Q-Q0)^2/(2 w^2)) + bg`  
**Q:** 0.003 – 0.5 Å⁻¹, 450 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `pop1_type` | size_dist / lognormal / sphere | - |  |
| `pop1_min_size` | 0 | A |  |
| `pop1_mean_size` | 60 | A |  |
| `pop1_sdeviation` | 0.25 | - |  |
| `pop1_scale` | 0.004 | - |  |
| `pop1_contrast` | 80 | 1e20 cm^-4 |  |
| `pop2_type` | diffraction_peak / gaussian | - |  |
| `pop2_amplitude` | 5 | 1/cm |  |
| `pop2_position` | 0.15 | 1/A |  |
| `pop2_width` | 0.012 | 1/A | Gaussian sigma, NOT FWHM |
| `background` | 0.005 | 1/cm |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `peak_FWHM` | 0.0282578 | 1/A |  |

Noise: κ = 50847.9, seed = 1155364236, median relative error 1.01 %, max 6.25 %.

### 10. `modeling_hardsphere`

**Tool:** Modeling  
**Model:** Log-normal spheres with a Percus-Yevick hard-sphere structure factor  
**Equation:** `I(Q) = I_sizedist(Q) * S_PY(Q; R_HS, phi_HS) + bg`  
**Q:** 0.003 – 0.5 Å⁻¹, 450 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `pop1_type` | size_dist / lognormal / sphere | - |  |
| `pop1_min_size` | 0 | A |  |
| `pop1_mean_size` | 50 | A |  |
| `pop1_sdeviation` | 0.12 | - |  |
| `pop1_scale` | 0.2 | - |  |
| `pop1_contrast` | 100 | 1e20 cm^-4 |  |
| `structure_factor` | hard_sphere | - |  |
| `sf_radius` | 50 | A |  |
| `sf_volume_fraction` | 0.25 | - |  |
| `background` | 0.01 | 1/cm |  |

Noise: κ = 6306.94, seed = 1762171967, median relative error 1.01 %, max 8.50 %.

### 11. `modeling_mass_fractal`

**Tool:** Modeling  
**Model:** Mass-fractal aggregate of spherical primary particles (Teixeira 1988)  
**Equation:** `I(Q) = phi*contrast*1e-4*V*[bracket*S_f(Q)+(1-eta)^2]*F_sph^2(QR) + bg`  
**Q:** 0.0003 – 0.3 Å⁻¹, 400 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `pop_type` | mass_fractal | - |  |
| `Phi` | 0.005 | - |  |
| `Radius` | 40 | A |  |
| `Dv` | 2.4 | - | mass fractal dimension |
| `Ksi` | 800 | A | aggregate correlation length |
| `Eta` | 0.5 | - | packing term |
| `Beta` | 1 | - | aspect ratio (sphere) |
| `Contrast` | 100 | 1e20 cm^-4 |  |
| `background` | 0.01 | 1/cm |  |

Noise: κ = 24414.7, seed = 2358675931, median relative error 1.00 %, max 6.47 %.

### 12. `modeling_surface_fractal`

**Tool:** Modeling  
**Model:** Surface-fractal scattering (Teixeira 1988)  
**Equation:** `I(Q) = pi*contrast*Ksi^4*S*Gamma(5-Ds)*sin[(3-Ds)atan(Q Ksi)]/[(1+Q^2Ksi^2)^((5-Ds)/2) Q Ksi] + bg`  
**Q:** 0.0005 – 0.3 Å⁻¹, 380 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `pop_type` | surface_fractal | - |  |
| `Surface` | 20000 | 1/cm |  |
| `Ds` | 2.4 | - | surface fractal dimension |
| `Ksi` | 600 | A |  |
| `Contrast` | 100 | 1e20 cm^-4 |  |
| `background` | 0.01 | 1/cm |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `Porod_slope_high_Q` | 3.6 | - |  |

Noise: κ = 8527.71, seed = 3936854831, median relative error 1.00 %, max 8.33 %.

### 13. `modeling_guinier_porod`

**Tool:** Modeling  
**Model:** Guinier-Porod model, single level with dimensionality s1 = 1 (rod-like)  
**Equation:** `I(Q<Q1) = G Q^-s1 exp(-Q^2 Rg1^2/(3-s1)); I(Q>=Q1) = D Q^-P`  
**Q:** 0.002 – 0.4 Å⁻¹, 400 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `pop_type` | guinier_porod | - |  |
| `G` | 150 | 1/cm |  |
| `Rg1` | 90 | A |  |
| `s1` | 1 | - | 0 = globular, 1 = rod, 2 = lamella |
| `P` | 3.6 | - |  |
| `Rg2` | 1e+10 | A | collapsed (single level) |
| `s2` | 0 | - |  |
| `background` | 0.005 | 1/cm |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `Q1` | 0.0179161 | 1/A |  |

Noise: κ = 3172.66, seed = 1241706239, median relative error 1.00 %, max 9.31 %.

### 14. `simple_guinier`

**Tool:** Simple Fits  
**Model:** Guinier  
**Equation:** `I(Q) = I0 exp(-Q^2 Rg^2/3)`  
**Q:** 0.002 – 0.03 Å⁻¹, 200 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `I0` | 250 | 1/cm |  |
| `Rg` | 45 | A |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `Q_max*Rg` | 1.35 | - | Guinier validity limit ~1.3 |

Noise: κ = 1.76605, seed = 192184972, median relative error 4.96 %, max 6.53 %.

### 15. `simple_guinier_rod`

**Tool:** Simple Fits  
**Model:** Guinier Rod  
**Equation:** `I(Q) = I0 exp(-Q^2 Rc^2/2)/Q`  
**Q:** 0.005 – 0.08 Å⁻¹, 200 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `I0` | 5 | cm^-1 A^-1 |  |
| `Rc` | 25 | A |  |

Noise: κ = 17.614, seed = 538831288, median relative error 1.89 %, max 8.25 %.

### 16. `simple_guinier_sheet`

**Tool:** Simple Fits  
**Model:** Guinier Sheet  
**Equation:** `I(Q) = I0 exp(-Q^2 Rg^2)/Q^2`  
**Q:** 0.01 – 0.1 Å⁻¹, 200 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `I0` | 0.02 | cm^-1 A^-2 |  |
| `Rg` | 15 | A |  |

Noise: κ = 658.905, seed = 2739133858, median relative error 1.40 %, max 8.54 %.

### 17. `simple_porod`

**Tool:** Simple Fits  
**Model:** Porod  
**Equation:** `I(Q) = Kp Q^-4 + bg`  
**Q:** 0.05 – 0.5 Å⁻¹, 200 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `Kp` | 2.5e-06 | cm^-1 A^-4 |  |
| `Background` | 0.05 | 1/cm |  |

Noise: κ = 5548.52, seed = 2585338660, median relative error 5.86 %, max 6.08 %.

### 18. `simple_power_law`

**Tool:** Simple Fits  
**Model:** Power Law  
**Equation:** `I(Q) = Prefactor Q^-Exponent + bg`  
**Q:** 0.01 – 0.4 Å⁻¹, 250 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `Prefactor` | 0.0012 | cm^-1 A^-n |  |
| `Exponent` | 3.2 | - |  |
| `Background` | 0.02 | 1/cm |  |

Noise: κ = 4580.16, seed = 4090444554, median relative error 1.12 %, max 7.24 %.

### 19. `simple_sphere`

**Tool:** Simple Fits  
**Model:** Sphere (monodisperse) + flat background  
**Equation:** `I(Q) = Scale |3(sin x - x cos x)/x^3|^2 + BG_flat, x = QR`  
**Q:** 0.002 – 0.15 Å⁻¹, 500 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `Scale` | 8 | 1/cm |  |
| `R` | 120 | A |  |
| `BG_flat` | 0.01 | 1/cm | fit with the complex background, flat term only |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `Q_first_minimum` | 0.037445 | 1/A |  |

Noise: κ = 27045, seed = 3077917839, median relative error 1.06 %, max 6.16 %.

### 20. `simple_debye_chain`

**Tool:** Simple Fits  
**Model:** Debye Polymer Chain  
**Equation:** `I(Q) = Scale 2(exp(-x)-1+x)/x^2, x = Q^2 Rg^2`  
**Q:** 0.003 – 0.3 Å⁻¹, 300 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `Scale` | 15 | 1/cm |  |
| `Rg` | 60 | A |  |

Noise: κ = 1902.13, seed = 1210075857, median relative error 1.34 %, max 7.61 %.

### 21. `simple_debye_bueche`

**Tool:** Simple Fits  
**Model:** Debye-Bueche  
**Equation:** `I(Q) = Prefactor Eta^2 xi^3/(1+Q^2 xi^2)^2`  
**Q:** 0.002 – 0.2 Å⁻¹, 300 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `Prefactor` | 1 | cm^-1 A^-3 |  |
| `Eta` | 0.05 | - |  |
| `CorrLength` | 80 | A |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `I(0)` | 1280 | 1/cm | only the product Prefactor*Eta^2 is determined by a fit |

Noise: κ = 5731.79, seed = 1541744988, median relative error 1.01 %, max 9.54 %.

### 22. `waxs_three_gauss`

**Tool:** WAXS Peak Fit  
**Model:** Three Gaussian peaks on a linear background  
**Equation:** `I(Q) = bg0 + bg1 Q + SUM_i A_i exp[-4 ln2 (Q-Q0_i)^2/FWHM_i^2]`  
**Q:** 1.2 – 4 Å⁻¹, 700 points, linear-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `background_shape` | Linear | - |  |
| `bg0` | 8 | 1/cm |  |
| `bg1` | 1.5 | cm^-1 A |  |
| `peak1_shape` | Gauss | - |  |
| `peak1_A` | 100 | 1/cm | peak height at Q0 |
| `peak1_Q0` | 1.85 | 1/A |  |
| `peak1_FWHM` | 0.09 | 1/A |  |
| `peak2_shape` | Gauss | - |  |
| `peak2_A` | 60 | 1/cm | peak height at Q0 |
| `peak2_Q0` | 2.55 | 1/A |  |
| `peak2_FWHM` | 0.12 | 1/A |  |
| `peak3_shape` | Gauss | - |  |
| `peak3_A` | 35 | 1/cm | peak height at Q0 |
| `peak3_Q0` | 3.05 | 1/A |  |
| `peak3_FWHM` | 0.15 | 1/A |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `peak1_area` | 9.5802 | cm^-1 A^-1 | A*FWHM*sqrt(pi/(4 ln2)) |
| `peak1_d_spacing` | 3.39632 | A |  |
| `peak2_area` | 7.66416 | cm^-1 A^-1 | A*FWHM*sqrt(pi/(4 ln2)) |
| `peak2_d_spacing` | 2.46399 | A |  |
| `peak3_area` | 5.58845 | cm^-1 A^-1 | A*FWHM*sqrt(pi/(4 ln2)) |
| `peak3_d_spacing` | 2.06006 | A |  |

Noise: κ = 27.75, seed = 1368163545, median relative error 5.32 %, max 6.15 %.

### 23. `waxs_pseudovoigt`

**Tool:** WAXS Peak Fit  
**Model:** Two pseudo-Voigt peaks on a cubic-polynomial background  
**Equation:** `I(Q) = SUM_k bg_k Q^k + SUM_i A_i[eta L_i(Q) + (1-eta) G_i(Q)]`  
**Q:** 1.5 – 4 Å⁻¹, 650 points, linear-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `background_shape` | Cubic | - |  |
| `bg0` | 10 | cm^-1 A^0 |  |
| `bg1` | -2 | cm^-1 A^1 |  |
| `bg2` | 0.8 | cm^-1 A^2 |  |
| `bg3` | -0.05 | cm^-1 A^3 |  |
| `peak1_shape` | Pseudo-Voigt | - |  |
| `peak1_A` | 200 | 1/cm |  |
| `peak1_Q0` | 2.1 | 1/A |  |
| `peak1_FWHM` | 0.07 | 1/A |  |
| `peak1_eta` | 0.4 | - |  |
| `peak2_shape` | Pseudo-Voigt | - |  |
| `peak2_A` | 90 | 1/cm |  |
| `peak2_Q0` | 2.95 | 1/A |  |
| `peak2_FWHM` | 0.11 | 1/A |  |
| `peak2_eta` | 0.4 | - |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `peak1_area` | 17.738 | cm^-1 A^-1 |  |
| `peak1_d_spacing` | 2.99199 | A |  |
| `peak2_area` | 12.5433 | cm^-1 A^-1 |  |
| `peak2_d_spacing` | 2.12989 | A |  |

Noise: κ = 30.3851, seed = 4129793739, median relative error 5.48 %, max 6.14 %.

### 24. `merge_usaxs`

**Tool:** Data Merge  
**Model:** Low-Q (USAXS-like) branch of a two-level Unified curve, correct scale  
**Equation:** `I(Q) = I_unified_2level(Q)`  
**Q:** 0.0001 – 0.05 Å⁻¹, 300 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `scale_factor` | 1 | - | reference branch |
| `G1` | 8 | - |  |
| `Rg1` | 120 | - |  |
| `P1` | 4 | - |  |
| `B1` | 1.88e-07 | - |  |
| `G2` | 4000 | - |  |
| `Rg2` | 1200 | - |  |
| `P2` | 3.2 | - |  |
| `B2` | 1.393e-06 | - |  |
| `RgCO2` | 120 | - |  |
| `bg` | 0.02 | - |  |

Noise: κ = 2240.84, seed = 183252891, median relative error 1.00 %, max 9.51 %.

### 25. `merge_saxs`

**Tool:** Data Merge  
**Model:** High-Q (SAXS-like) branch of the SAME curve, deliberately mis-scaled by 1.15  
**Equation:** `I(Q) = 1.15 * I_unified_2level(Q)`  
**Q:** 0.005 – 0.5 Å⁻¹, 320 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `scale_factor` | 1.15 | - | merge should recover 1/1.15 = 0.869565 for this branch |
| `G1` | 8 | - |  |
| `Rg1` | 120 | - |  |
| `P1` | 4 | - |  |
| `B1` | 1.88e-07 | - |  |
| `G2` | 4000 | - |  |
| `Rg2` | 1200 | - |  |
| `P2` | 3.2 | - |  |
| `B2` | 1.393e-06 | - |  |
| `RgCO2` | 120 | - |  |
| `bg` | 0.02 | - |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `expected_recovered_scale` | 0.869565 | - |  |

Noise: κ = 12072.7, seed = 652841110, median relative error 3.93 %, max 6.08 %.

### 26. `manipulate_reference`

**Tool:** Data Manipulation  
**Model:** Reference curve (single Unified level)  
**Equation:** `I(Q) = G exp(-Q^2Rg^2/3) + B/Q*^P`  
**Q:** 0.001 – 0.3 Å⁻¹, 300 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `G` | 100 | 1/cm |  |
| `Rg` | 200 | A |  |
| `P` | 4 | - |  |
| `B` | 3.045e-07 | cm^-1 A^-P |  |

Noise: κ = 2.36109e+06, seed = 3138783688, median relative error 1.00 %, max 10.66 %.

### 27. `manipulate_input`

**Tool:** Data Manipulation  
**Model:** The reference curve scaled and offset by known constants  
**Equation:** `I_input(Q) = a * I_reference(Q) + b  (applied to the NOISY reference, so the inverse is exact)`  
**Q:** 0.001 – 0.3 Å⁻¹, 300 points, log-spaced  

| parameter | value | unit | note |
|---|---|---|---|
| `a_multiply` | 2.5 | - |  |
| `b_add` | 0.3 | 1/cm |  |

Derived (not fitted):

| quantity | value | unit | note |
|---|---|---|---|
| `inverse_operation` | subtract 0.3 then divide by 2.5 | - | must reproduce manipulate_reference to machine precision |

Noise: κ = nan, seed = 3138783688, median relative error 0.96 %, max 1.03 %.

