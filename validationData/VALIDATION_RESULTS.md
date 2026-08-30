# pyIrena validation results

Every row compares a parameter recovered by **pyIrena** with the **exact value
used to synthesise the data**.  The data files, and the generator that produced
them, are in this folder; see `README.md` for the full ground truth and
`ground_truth.json` for the machine-readable version.

The **Irena** columns are deliberately empty.  Analyse the same files in the
Igor Pro Irena package using the settings quoted for each dataset, enter the
values, and the table becomes a direct implementation-to-implementation
comparison against a common, exactly known reference.

Deviation is `100 (fitted - true) / true`.  `Tol` is the tolerance the pyIrena
regression test enforces; a blank tolerance marks a quantity that is reported
for information but is not expected to be individually determined by the data
(a correlated prefactor, a polynomial coefficient, a reduced chi-squared).

Starting values for every fit were offset from the truth (typically 0.6x) so
that the tables demonstrate convergence rather than assume it.


Generated 2026-08-30 by `validationData/run_validation_report.py`.

**Summary: 195 of 195 scored comparisons within tolerance** (225 rows in total).

## Size Distribution

### `sizes_sphere_lognormal`

*pyIrena settings — use the same in Irena:* maxent, sphere, r = 20-320 A in 80 linear bins, contrast = 100e20 cm^-4

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| volume fraction | - | 0.01 | 0.00995975 | -0.402 | 6 | yes |  |  |
| R at distribution peak | A | 93.9413 | 95.9494 | +2.138 | 10 | yes |  |  |
| volume-weighted mean R | A | 103.174 | 103.201 | +0.026 | 6 | yes |  |  |
| RMS radius (pyIrena 'Rg') | A | 106.449 | 106.544 | +0.089 | 6 | yes |  |  |
| reduced chi^2 | - | 1 | 0.999466 | -0.053 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.14132 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 1.18765 |  | 15 | yes |  |  |

> regularised inversions broaden the distribution, so the modal radius has a wider tolerance than the integral moments

> fitted model compared point-by-point with <name>_ideal.dat

### `sizes_sphere_bimodal`

*pyIrena settings — use the same in Irena:* maxent, sphere, r = 10-600 A in 120 log bins, contrast = 100e20 cm^-4

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| mode 1: R at peak | A | 38.7248 | 39.6002 | +2.261 | 10 | yes |  |  |
| mode 1: volume fraction | - | 0.004 | 0.00400567 | +0.142 | 10 | yes |  |  |
| mode 2: R at peak | A | 238.188 | 245.272 | +2.974 | 10 | yes |  |  |
| mode 2: volume fraction | - | 0.01 | 0.00992221 | -0.778 | 10 | yes |  |  |
| total volume fraction | - | 0.014 | 0.0139301 | -0.499 | 6 | yes |  |  |
| reduced chi^2 | - | 1 | 0.999984 | -0.002 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.493856 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 7.62201 |  | 15 | yes |  |  |

> population separated at r = 100 A

> fitted model compared point-by-point with <name>_ideal.dat

### `sizes_sphere_flat_background`

*pyIrena settings — use the same in Irena:* averaged over Q = 0.3-0.35 1/A, the same recipe Irena uses; the residual particle signal there biases it slightly high

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| background_flat | 1/cm | 0.02 | 0.0203018 | +1.509 | 15 | yes |  |  |
| volume fraction | - | 0.006 | 0.00599547 | -0.076 | 6 | yes |  |  |
| R at distribution peak | A | 114.33 | 111.392 | -2.570 | 10 | yes |  |  |
| volume-weighted mean R | A | 122.939 | 122.879 | -0.049 | 6 | yes |  |  |
| RMS radius (pyIrena 'Rg') | A | 125.951 | 125.933 | -0.014 | 6 | yes |  |  |
| reduced chi^2 | - | 1 | 1.03303 | +3.303 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.118391 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 2.32986 |  | 15 | yes |  |  |

> regularised inversions broaden the distribution, so the modal radius has a wider tolerance than the integral moments

> fitted model compared point-by-point with <name>_ideal.dat

### `sizes_spheroid_lognormal`

*pyIrena settings — use the same in Irena:* maxent, spheroid, r = 30-450 A in 80 linear bins, contrast = 50e20 cm^-4

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| volume fraction | - | 0.008 | 0.00798119 | -0.235 | 6 | yes |  |  |
| R at distribution peak | A | 144.118 | 152.278 | +5.662 | 10 | yes |  |  |
| RMS radius (pyIrena 'Rg') | A | 156.122 | 156.277 | +0.100 | 6 | yes |  |  |
| reduced chi^2 | - | 1 | 1.00062 | +0.062 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.0889834 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 1.10375 |  | 15 | yes |  |  |

> regularised inversions broaden the distribution, so the modal radius has a wider tolerance than the integral moments

> fitted model compared point-by-point with <name>_ideal.dat

## Unified Fit

### `unified_one_level`

*pyIrena settings — use the same in Irena:* 1 level(s), fit G/Rg/B/P + flat background, TRF, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| G | 1/cm | 100 | 99.9499 | -0.050 | 3 | yes |  |  |
| Rg | A | 200 | 199.844 | -0.078 | 2 | yes |  |  |
| P | - | 4 | 4.00924 | +0.231 | 2 | yes |  |  |
| B | cm^-1 A^-P | 3.04500e-07 | 2.95156e-07 | -3.069 | 25 | yes |  |  |
| background | 1/cm | 0.01 | 0.01006 | +0.600 | 15 | yes |  |  |
| reduced chi^2 | - | 1 | 0.979006 | -2.099 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.0832303 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.590335 |  | 12 | yes |  |  |

> B and P are strongly correlated (B carries units cm^-1 A^-P), so B has a wide tolerance whenever P is a free parameter

> 1.0 expected for a correct model and correct error bars

> fitted model compared point-by-point with <name>_ideal.dat

### `unified_two_level`

*pyIrena settings — use the same in Irena:* 2 level(s), fit G/Rg/B/P + flat background, TRF, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| G_1 | 1/cm | 8 | 7.96184 | -0.477 | 3 | yes |  |  |
| Rg_1 | A | 120 | 119.829 | -0.143 | 2 | yes |  |  |
| P_1 | - | 4 | 3.93526 | -1.618 | 2 | yes |  |  |
| B_1 | cm^-1 A^-P | 1.88000e-07 | 2.32120e-07 | +23.468 | 25 | yes |  |  |
| G_2 | 1/cm | 4000 | 4001.17 | +0.029 | 3 | yes |  |  |
| Rg_2 | A | 1200 | 1202 | +0.167 | 2 | yes |  |  |
| P_2 | - | 3.2 | 3.18946 | -0.329 | 2 | yes |  |  |
| B_2 | cm^-1 A^-P | 1.39300e-06 | 1.47606e-06 | +5.963 | 25 | yes |  |  |
| background | 1/cm | 0.02 | 0.0201605 | +0.802 | 15 | yes |  |  |
| reduced chi^2 | - | 1 | 1.06538 | +6.538 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.153929 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 1.7138 |  | 12 | yes |  |  |

> B and P are strongly correlated (B carries units cm^-1 A^-P), so B has a wide tolerance whenever P is a free parameter

> 1.0 expected for a correct model and correct error bars

> fitted model compared point-by-point with <name>_ideal.dat

## Unified Fit (slit smearing)

### `unified_one_level_slitsmeared`

*pyIrena settings — use the same in Irena:* 1 level(s), fit G/Rg/B/P + flat background, TRF, start = 0.60x truth, slit smearing dQl=0.03

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| G | 1/cm | 100 | 99.8555 | -0.144 | 3 | yes |  |  |
| Rg | A | 200 | 200.098 | +0.049 | 2 | yes |  |  |
| P | - | 4 | 4.00443 | +0.111 | 2 | yes |  |  |
| B | cm^-1 A^-P | 3.04500e-07 | 3.00025e-07 | -1.470 | 25 | yes |  |  |
| background | 1/cm | 0.01 | 0.00996062 | -0.394 | 15 | yes |  |  |
| reduced chi^2 | - | 1 | 0.930125 | -6.988 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.181991 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.412501 |  | 12 | yes |  |  |

> B and P are strongly correlated (B carries units cm^-1 A^-P), so B has a wide tolerance whenever P is a free parameter

> 1.0 expected for a correct model and correct error bars

> fitted model compared point-by-point with <name>_ideal.dat

## Modeling

### `modeling_sizedist_plus_unified`

*pyIrena settings — use the same in Irena:* populations: size_dist, unified_level; local (TRF) fit, size-distribution grid 200 bins, contrast held fixed, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| pop1 mean_size (log-normal median) | A | 30 | 30.0419 | +0.140 | 4 | yes |  |  |
| pop1 sdeviation | - | 0.3 | 0.316207 | +5.402 | 6 | yes |  |  |
| pop1 scale | - | 0.005 | 0.00500931 | +0.186 | 6 | yes |  |  |
| pop2 G | 1/cm | 20000 | 19996 | -0.020 | 4 | yes |  |  |
| pop2 Rg | A | 4000 | 3998.8 | -0.030 | 4 | yes |  |  |
| pop2 P | - | 3.5 | 3.49239 | -0.217 | 3 | yes |  |  |
| pop2 B | cm^-1 A^-P | 1.56345e-08 | 1.64675e-08 | +5.328 | 25 | yes |  |  |
| background | 1/cm | 0.01 | 0.00994733 | -0.527 | 30 | yes |  |  |
| reduced chi^2 | - | 1 | 1.03921 | +3.921 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.213212 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 2.16489 |  | 15 | yes |  |  |

> a small flat background is a weakly determined nuisance parameter

> fitted model compared point-by-point with <name>_ideal.dat

### `modeling_sizedist_plus_peak`

*pyIrena settings — use the same in Irena:* populations: size_dist, diffraction_peak; local (TRF) fit, size-distribution grid 200 bins, contrast held fixed, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| pop1 mean_size (log-normal median) | A | 60 | 60.1816 | +0.303 | 4 | yes |  |  |
| pop1 sdeviation | - | 0.25 | 0.260352 | +4.141 | 6 | yes |  |  |
| pop1 scale | - | 0.004 | 0.00400364 | +0.091 | 6 | yes |  |  |
| pop2 amplitude | 1/cm | 5 | 4.99142 | -0.172 | 5 | yes |  |  |
| pop2 position | 1/A | 0.15 | 0.149981 | -0.013 | 1 | yes |  |  |
| pop2 width (sigma) | 1/A | 0.012 | 0.0120011 | +0.009 | 6 | yes |  |  |
| background | 1/cm | 0.005 | 0.00501291 | +0.258 | 30 | yes |  |  |
| reduced chi^2 | - | 1 | 1.18594 | +18.594 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.335475 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 2.45529 |  | 15 | yes |  |  |

> a small flat background is a weakly determined nuisance parameter

> fitted model compared point-by-point with <name>_ideal.dat

### `modeling_hardsphere`

*pyIrena settings — use the same in Irena:* populations: size_dist; local (TRF) fit, size-distribution grid 200 bins, contrast held fixed, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| pop1 mean_size (log-normal median) | A | 50 | 49.983 | -0.034 | 4 | yes |  |  |
| pop1 sdeviation | - | 0.12 | 0.123946 | +3.288 | 6 | yes |  |  |
| pop1 scale | - | 0.2 | 0.200044 | +0.022 | 6 | yes |  |  |
| hard-sphere radius | A | 50 | 49.9432 | -0.114 | 6 | yes |  |  |
| hard-sphere volume fraction | - | 0.25 | 0.249643 | -0.143 | 10 | yes |  |  |
| background | 1/cm | 0.01 | 0.0093383 | -6.617 | 30 | yes |  |  |
| reduced chi^2 | - | 1 | 1.15067 | +15.067 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.103617 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 3.05998 |  | 15 | yes |  |  |

> a small flat background is a weakly determined nuisance parameter

> fitted model compared point-by-point with <name>_ideal.dat

### `modeling_mass_fractal`

*pyIrena settings — use the same in Irena:* populations: mass_fractal; local (TRF) fit, size-distribution grid 200 bins, contrast held fixed, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| Phi | - | 0.005 | 0.004999 | -0.020 | 6 | yes |  |  |
| Radius | A | 40 | 40.019 | +0.047 | 4 | yes |  |  |
| Dv | - | 2.4 | 2.39923 | -0.032 | 3 | yes |  |  |
| Ksi | A | 800 | 800.352 | +0.044 | 8 | yes |  |  |
| background | 1/cm | 0.01 | 0.0102148 | +2.148 | 30 | yes |  |  |
| reduced chi^2 | - | 1 | 1.09383 | +9.383 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.0725479 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 2.27086 |  | 15 | yes |  |  |

> a small flat background is a weakly determined nuisance parameter

> fitted model compared point-by-point with <name>_ideal.dat

### `modeling_surface_fractal`

*pyIrena settings — use the same in Irena:* populations: surface_fractal; local (TRF) fit, size-distribution grid 200 bins, contrast held fixed, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| Surface | 1/cm | 20000 | 20016.1 | +0.081 | 6 | yes |  |  |
| Ds | - | 2.4 | 2.39989 | -0.004 | 3 | yes |  |  |
| Ksi | A | 600 | 599.633 | -0.061 | 8 | yes |  |  |
| background | 1/cm | 0.01 | 0.00967284 | -3.272 | 30 | yes |  |  |
| reduced chi^2 | - | 1 | 1.03326 | +3.326 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.0466197 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 1.89966 |  | 15 | yes |  |  |

> a small flat background is a weakly determined nuisance parameter

> fitted model compared point-by-point with <name>_ideal.dat

### `modeling_guinier_porod`

*pyIrena settings — use the same in Irena:* populations: guinier_porod; local (TRF) fit, size-distribution grid 200 bins, contrast held fixed, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| G | 1/cm | 150 | 149.306 | -0.463 | 5 | yes |  |  |
| Rg1 | A | 90 | 89.8951 | -0.117 | 4 | yes |  |  |
| s1 | - | 1 | 1.00091 | +0.091 | 5 | yes |  |  |
| P | - | 3.6 | 3.60163 | +0.045 | 3 | yes |  |  |
| background | 1/cm | 0.005 | 0.00631306 | +26.261 | 30 | yes |  |  |
| reduced chi^2 | - | 1 | 1.02814 | +2.814 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.0695124 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 3.25451 |  | 15 | yes |  |  |

> a small flat background is a weakly determined nuisance parameter

> fitted model compared point-by-point with <name>_ideal.dat

## Simple Fits

### `simple_guinier`

*pyIrena settings — use the same in Irena:* model 'Guinier', all listed parameters free, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| I0 | 1/cm | 250 | 250.686 | +0.275 | 3 | yes |  |  |
| Rg | A | 45 | 45.9727 | +2.162 | 3 | yes |  |  |
| reduced chi^2 | - | 1 | 0.871892 | -12.811 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.240806 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 2.35248 |  | 12 | yes |  |  |

> fitted model compared point-by-point with <name>_ideal.dat

### `simple_guinier_rod`

*pyIrena settings — use the same in Irena:* model 'Guinier Rod', all listed parameters free, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| I0 | cm^-1 A^-1 | 5 | 5.0042 | +0.084 | 3 | yes |  |  |
| Rc | A | 25 | 25.072 | +0.288 | 3 | yes |  |  |
| reduced chi^2 | - | 1 | 1.04912 | +4.912 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.0742042 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 1.06439 |  | 12 | yes |  |  |

> fitted model compared point-by-point with <name>_ideal.dat

### `simple_guinier_sheet`

*pyIrena settings — use the same in Irena:* model 'Guinier Sheet', all listed parameters free, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| I0 | cm^-1 A^-2 | 0.02 | 0.0200201 | +0.100 | 3 | yes |  |  |
| Rg | A | 15 | 15.0031 | +0.021 | 3 | yes |  |  |
| reduced chi^2 | - | 1 | 0.988868 | -1.113 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.0910625 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.0995459 |  | 12 | yes |  |  |

> fitted model compared point-by-point with <name>_ideal.dat

### `simple_porod`

*pyIrena settings — use the same in Irena:* model 'Porod', all listed parameters free, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| Kp | cm^-1 A^-4 | 2.50000e-06 | 2.51216e-06 | +0.486 | 3 | yes |  |  |
| Background | 1/cm | 0.05 | 0.0499857 | -0.029 | 3 | yes |  |  |
| reduced chi^2 | - | 1 | 1.06544 | +6.544 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.0280001 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.42908 |  | 12 | yes |  |  |

> fitted model compared point-by-point with <name>_ideal.dat

### `simple_power_law`

*pyIrena settings — use the same in Irena:* model 'Power Law', all listed parameters free, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| Prefactor | cm^-1 A^-n | 0.0012 | 0.0012047 | +0.391 | 3 | yes |  |  |
| Exponent | - | 3.2 | 3.19894 | -0.033 | 1 | yes |  |  |
| Background | 1/cm | 0.02 | 0.0195288 | -2.356 | 3 | yes |  |  |
| reduced chi^2 | - | 1 | 0.933632 | -6.637 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.079124 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.952092 |  | 12 | yes |  |  |

> fitted model compared point-by-point with <name>_ideal.dat

### `simple_sphere`

*pyIrena settings — use the same in Irena:* model 'Sphere', all listed parameters free, start = 0.60x truth, complex background (flat term only)

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| Scale | 1/cm | 8 | 7.99862 | -0.017 | 3 | yes |  |  |
| R | A | 120 | 119.965 | -0.029 | 3 | yes |  |  |
| BG_flat | 1/cm | 0.01 | 0.0100646 | +0.646 | 3 | yes |  |  |
| reduced chi^2 | - | 1 | 1.01504 | +1.504 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.0346539 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.983124 |  | 12 | yes |  |  |

> fitted model compared point-by-point with <name>_ideal.dat

### `simple_debye_chain`

*pyIrena settings — use the same in Irena:* model 'Debye Polymer Chain', all listed parameters free, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| Scale | 1/cm | 15 | 14.99 | -0.067 | 3 | yes |  |  |
| Rg | A | 60 | 60.0224 | +0.037 | 3 | yes |  |  |
| reduced chi^2 | - | 1 | 1.08716 | +8.716 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.11394 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.14104 |  | 12 | yes |  |  |

> fitted model compared point-by-point with <name>_ideal.dat

### `simple_debye_bueche`

*pyIrena settings — use the same in Irena:* model 'Debye-Bueche', all listed parameters free, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| Prefactor | cm^-1 A^-3 | 1 | 0.817249 | -18.275 |  | - |  |  |
| Eta | - | 0.05 | 0.0553398 | +10.680 |  | - |  |  |
| CorrLength | A | 80 | 79.951 | -0.061 | 3 | yes |  |  |
| Prefactor * Eta^2 | cm^-1 A^-3 | 0.0025 | 0.00250282 | +0.113 | 3 | yes |  |  |
| reduced chi^2 | - | 1 | 0.967369 | -3.263 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.105291 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.173263 |  | 12 | yes |  |  |

> degenerate with the other prefactor - only the product is determined by the data

> the combination the data actually determine

> fitted model compared point-by-point with <name>_ideal.dat

## WAXS Peak Fit

### `waxs_three_gauss`

*pyIrena settings — use the same in Irena:* 3 x Gauss peaks, Linear background, all parameters free; A and background start at 0.60x truth, Q0 within 1 %, FWHM 1.30x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| bg0 | - | 8 | 7.87062 | -1.617 | 10 | yes |  |  |
| bg1 | - | 1.5 | 1.54147 | +2.765 | 10 | yes |  |  |
| background at Q = 1.2 | 1/cm | 9.8 | 9.72039 | -0.812 | 5 | yes |  |  |
| background at Q = 1.9 | 1/cm | 10.85 | 10.7994 | -0.466 | 5 | yes |  |  |
| background at Q = 2.6 | 1/cm | 11.9 | 11.8784 | -0.181 | 5 | yes |  |  |
| background at Q = 3.3 | 1/cm | 12.95 | 12.9575 | +0.058 | 5 | yes |  |  |
| background at Q = 4 | 1/cm | 14 | 14.0365 | +0.261 | 5 | yes |  |  |
| peak1_A | 1/cm | 100 | 99.5683 | -0.432 | 3 | yes |  |  |
| peak1_Q0 | 1/A | 1.85 | 1.85011 | +0.006 | 0.5 | yes |  |  |
| peak1_FWHM | 1/A | 0.09 | 0.0902364 | +0.263 | 3 | yes |  |  |
| peak1_area | cm^-1 A^-1 | 9.5802 | 9.5639 | -0.170 | 4 | yes |  |  |
| peak2_A | 1/cm | 60 | 59.7818 | -0.364 | 3 | yes |  |  |
| peak2_Q0 | 1/A | 2.55 | 2.55047 | +0.018 | 0.5 | yes |  |  |
| peak2_FWHM | 1/A | 0.12 | 0.120372 | +0.310 | 3 | yes |  |  |
| peak2_area | cm^-1 A^-1 | 7.66416 | 7.65996 | -0.055 | 4 | yes |  |  |
| peak3_A | 1/cm | 35 | 34.9834 | -0.047 | 3 | yes |  |  |
| peak3_Q0 | 1/A | 3.05 | 3.04892 | -0.035 | 0.5 | yes |  |  |
| peak3_FWHM | 1/A | 0.15 | 0.150923 | +0.615 | 3 | yes |  |  |
| peak3_area | cm^-1 A^-1 | 5.58845 | 5.62018 | +0.568 | 4 | yes |  |  |
| reduced chi^2 | - | 1 | 1.03093 | +3.093 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.327711 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 1.84068 |  | 12 | yes |  |  |

> polynomial background coefficients are strongly correlated and individually ill-determined; compare the background CURVE below

> background curve, the quantity that is determined

> analytic integral of the fitted peak

> fitted model compared point-by-point with <name>_ideal.dat

### `waxs_pseudovoigt`

*pyIrena settings — use the same in Irena:* 2 x Pseudo-Voigt peaks, Cubic background, all parameters free; A and background start at 0.60x truth, Q0 within 1 %, FWHM 1.30x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| bg0 | - | 10 | 11.9979 | +19.979 |  | - |  |  |
| bg1 | - | -2 | -4.51983 | +125.992 |  | - |  |  |
| bg2 | - | 0.8 | 1.84822 | +131.027 |  | - |  |  |
| bg3 | - | -0.05 | -0.189004 | +278.007 |  | - |  |  |
| background at Q = 1.5 | 1/cm | 8.63125 | 8.7387 | +1.245 | 5 | yes |  |  |
| background at Q = 2.12 | 1/cm | 8.88271 | 8.92544 | +0.481 | 5 | yes |  |  |
| background at Q = 2.75 | 1/cm | 9.51016 | 9.61476 | +1.100 | 5 | yes |  |  |
| background at Q = 3.38 | 1/cm | 10.4403 | 10.5298 | +0.857 | 5 | yes |  |  |
| background at Q = 4 | 1/cm | 11.6 | 11.3937 | -1.778 | 5 | yes |  |  |
| peak1_A | 1/cm | 200 | 197.484 | -1.258 | 3 | yes |  |  |
| peak1_Q0 | 1/A | 2.1 | 2.09964 | -0.017 | 0.5 | yes |  |  |
| peak1_FWHM | 1/A | 0.07 | 0.0703344 | +0.478 | 3 | yes |  |  |
| peak1_eta | - | 0.4 | 0.406171 | +1.543 | 15 | yes |  |  |
| peak1_area | cm^-1 A^-1 | 17.738 | 17.6419 | -0.542 | 4 | yes |  |  |
| peak2_A | 1/cm | 90 | 89.99 | -0.011 | 3 | yes |  |  |
| peak2_Q0 | 1/A | 2.95 | 2.95028 | +0.010 | 0.5 | yes |  |  |
| peak2_FWHM | 1/A | 0.11 | 0.109408 | -0.538 | 3 | yes |  |  |
| peak2_eta | - | 0.4 | 0.397555 | -0.611 | 15 | yes |  |  |
| peak2_area | cm^-1 A^-1 | 12.5433 | 12.4622 | -0.646 | 4 | yes |  |  |
| reduced chi^2 | - | 1 | 1.02501 | +2.501 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.676722 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 1.79415 |  | 12 | yes |  |  |

> polynomial background coefficients are strongly correlated and individually ill-determined; compare the background CURVE below

> background curve, the quantity that is determined

> analytic integral of the fitted peak

> fitted model compared point-by-point with <name>_ideal.dat

## Data Merge

### `merge_usaxs + merge_saxs`

*pyIrena settings — use the same in Irena:* DS1 = merge_usaxs (absolute), DS2 = merge_saxs (scaled); overlap Q = 0.006-0.04 1/A, log-log interpolation, scale DS2, no Q shift

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| recovered scale for DS2 | - | 0.869565 | 0.871304 | +0.200 | 1 | yes |  |  |
| fitted background offset | 1/cm | 0 | 0.0017514 |  |  | - |  |  |

> DS2 was generated 1.15x too high, so the correct scale is 1/1.15

> should be ~0: the two branches differ by a pure multiplicative factor

## Data Manipulation

### `manipulate_input`

*pyIrena settings — use the same in Irena:* Scale + Background: I_out = (1/2.5)*I - 0.12, the exact inverse of the generating operation

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| max abs. deviation from manipulate_reference | % | 0 | 3.15555e-11 |  | 1.00000e-08 | yes |  |  |

> the inverse operation must reproduce the reference to machine precision

## Scattering Contrast

### `contrast: SiO2 (2.2 g/cm3)`

*pyIrena settings — use the same in Irena:* free-electron X-ray SLD and bound-coherent neutron SLD; reference values from IUPAC 2021 atomic weights and Sears (1992) scattering lengths

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| formula weight | g/mol | 60.083 | 60.083 | +0.000 | 0.2 | yes |  |  |
| X-ray SLD | 1e10 cm^-2 | 18.6412 | 18.6412 | +0.000 | 0.5 | yes |  |  |
| neutron SLD | 1e10 cm^-2 | 3.47411 | 3.47477 | +0.019 | 1 | yes |  |  |

> amorphous silica

### `contrast: H2O (1 g/cm3)`

*pyIrena settings — use the same in Irena:* free-electron X-ray SLD and bound-coherent neutron SLD; reference values from IUPAC 2021 atomic weights and Sears (1992) scattering lengths

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| formula weight | g/mol | 18.015 | 18.015 | +0.000 | 0.2 | yes |  |  |
| X-ray SLD | 1e10 cm^-2 | 9.41995 | 9.41995 | +0.000 | 0.5 | yes |  |  |
| neutron SLD | 1e10 cm^-2 | -0.559927 | -0.560963 | +0.185 | 1 | yes |  |  |

> light water

### `contrast: D2O (1.107 g/cm3)`

*pyIrena settings — use the same in Irena:* free-electron X-ray SLD and bound-coherent neutron SLD; reference values from IUPAC 2021 atomic weights and Sears (1992) scattering lengths

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| formula weight | g/mol | 20.027 | 20.0272 | +0.001 | 0.2 | yes |  |  |
| X-ray SLD | 1e10 cm^-2 | 9.38025 | 9.38015 | -0.001 | 0.5 | yes |  |  |
| neutron SLD | 1e10 cm^-2 | 6.37291 | 6.37115 | -0.028 | 1 | yes |  |  |

> heavy water

### `contrast: Al2O3 (3.97 g/cm3)`

*pyIrena settings — use the same in Irena:* free-electron X-ray SLD and bound-coherent neutron SLD; reference values from IUPAC 2021 atomic weights and Sears (1992) scattering lengths

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| formula weight | g/mol | 101.961 | 101.96 | -0.001 | 0.2 | yes |  |  |
| X-ray SLD | 1e10 cm^-2 | 33.0376 | 33.0379 | +0.001 | 0.5 | yes |  |  |
| neutron SLD | 1e10 cm^-2 | 5.69953 | 5.70007 | +0.010 | 1 | yes |  |  |

> corundum

### `contrast: Fe (7.874 g/cm3)`

*pyIrena settings — use the same in Irena:* free-electron X-ray SLD and bound-coherent neutron SLD; reference values from IUPAC 2021 atomic weights and Sears (1992) scattering lengths

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| formula weight | g/mol | 55.845 | 55.845 | +0.000 | 0.2 | yes |  |  |
| X-ray SLD | 1e10 cm^-2 | 62.211 | 62.211 | +0.000 | 0.5 | yes |  |  |
| neutron SLD | 1e10 cm^-2 | 8.02405 | 8.02405 | +0.000 | 1 | yes |  |  |

> alpha iron

### `contrast: C8H8 (1.05 g/cm3)`

*pyIrena settings — use the same in Irena:* free-electron X-ray SLD and bound-coherent neutron SLD; reference values from IUPAC 2021 atomic weights and Sears (1992) scattering lengths

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| formula weight | g/mol | 104.152 | 104.152 | +0.000 | 0.2 | yes |  |  |
| X-ray SLD | 1e10 cm^-2 | 9.58059 | 9.58059 | +0.000 | 0.5 | yes |  |  |
| neutron SLD | 1e10 cm^-2 | 1.41191 | 1.41157 | -0.024 | 1 | yes |  |  |

> polystyrene

### `contrast: SiO2 vs H2O`

*pyIrena settings — use the same in Irena:* free-electron X-ray SLD and bound-coherent neutron SLD; reference values from IUPAC 2021 atomic weights and Sears (1992) scattering lengths

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| X-ray contrast (delta-rho)^2 | 1e20 cm^-4 | 85.0324 | 85.0324 | +0.000 | 1 | yes |  |  |
| neutron contrast (delta-rho)^2 | 1e20 cm^-4 | 16.2734 | 16.2871 | +0.084 | 2 | yes |  |  |

### `contrast: SiO2 vs D2O`

*pyIrena settings — use the same in Irena:* free-electron X-ray SLD and bound-coherent neutron SLD; reference values from IUPAC 2021 atomic weights and Sears (1992) scattering lengths

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| X-ray contrast (delta-rho)^2 | 1e20 cm^-4 | 85.7661 | 85.7679 | +0.002 | 1 | yes |  |  |
| neutron contrast (delta-rho)^2 | 1e20 cm^-4 | 8.40308 | 8.38902 | -0.167 | 2 | yes |  |  |

### `contrast: Fe vs H2O`

*pyIrena settings — use the same in Irena:* free-electron X-ray SLD and bound-coherent neutron SLD; reference values from IUPAC 2021 atomic weights and Sears (1992) scattering lengths

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| X-ray contrast (delta-rho)^2 | 1e20 cm^-4 | 2786.89 | 2786.89 | +0.000 | 1 | yes |  |  |
| neutron contrast (delta-rho)^2 | 1e20 cm^-4 | 73.6847 | 73.7025 | +0.024 | 2 | yes |  |  |

