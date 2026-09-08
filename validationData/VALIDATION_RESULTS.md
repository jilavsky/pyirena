# pyIrena validation results

Every row compares a parameter recovered by **pyIrena** with the **exact value
used to synthesise the data**.  The data files, and the generator that produced
them, are in this folder; see `README.md` for the full ground truth and
`ground_truth.json` for the machine-readable version.

The **Irena** column holds the corresponding value obtained with the Igor Pro
Irena package, analysing the same file with the settings quoted for each
dataset.  Where it is filled, the table is a three-way comparison: both
implementations measured against one common, exactly known reference.  A blank
Irena cell means that quantity has not been analysed in Irena or that Irena does
not report it; `*` marks a quantity Irena defines differently, so the numbers
are not directly comparable.  These values are read from `irena_values.csv`.

Both deviations are `100 (fitted - true) / true`.  `Tol` is the tolerance the pyIrena
regression test enforces; a blank tolerance marks a quantity that is reported
for information but is not expected to be individually determined by the data
(a correlated prefactor, a polynomial coefficient, a reduced chi-squared).

Starting values for every fit were offset from the truth (typically 0.6x) so
that the tables demonstrate convergence rather than assume it.


Generated 2026-09-08 by `validationData/run_validation_report.py`.

**Summary: 195 of 195 scored comparisons within tolerance** (225 rows in total).

Over the 98 quantities both packages report, the median deviation from the known truth is **0.19 % for pyIrena and 0.31 % for Irena**, and the median difference **between the two packages is 0.04 %**.

## Size Distribution

### `sizes_sphere_lognormal`

*pyIrena settings — use the same in Irena:* maxent, sphere, r = 20-320 A in 80 linear bins, contrast = 100e20 cm^-4

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| volume fraction | - | 0.01 | 0.00995975 | -0.402 | 6 | yes | 0.00996809 | -0.319 |
| R at distribution peak | A | 93.9413 | 95.9494 | +2.138 | 10 | yes | 95.95 | +2.138 |
| volume-weighted mean R | A | 103.174 | 103.201 | +0.026 | 6 | yes | 103.19 | +0.015 |
| RMS radius (pyIrena 'Rg') | A | 106.449 | 106.544 | +0.089 | 6 | yes | 106.52 | +0.066 |
| reduced chi^2 | - | 1 | 0.999466 | -0.053 |  | - | 1.0001 | +0.010 |
| model vs exact curve: median abs. deviation | % | 0 | 0.14132 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 1.18765 |  | 15 | yes |  |  |

> regularised inversions broaden the distribution, so the modal radius has a wider tolerance than the integral moments

> fitted model compared point-by-point with <name>_ideal.dat

### `sizes_sphere_bimodal`

*pyIrena settings — use the same in Irena:* maxent, sphere, r = 10-600 A in 120 log bins, contrast = 100e20 cm^-4

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| mode 1: R at peak | A | 38.7248 | 39.6002 | +2.261 | 10 | yes | 36.9 | -4.712 |
| mode 1: volume fraction | - | 0.004 | 0.00400567 | +0.142 | 10 | yes | 0.00400055 | +0.014 |
| mode 2: R at peak | A | 238.188 | 245.272 | +2.974 | 10 | yes | 245.27 | +2.973 |
| mode 2: volume fraction | - | 0.01 | 0.00992221 | -0.778 | 10 | yes | 0.0098972 | -1.028 |
| total volume fraction | - | 0.014 | 0.0139301 | -0.499 | 6 | yes | 0.0139006 | -0.710 |
| reduced chi^2 | - | 1 | 0.999984 | -0.002 |  | - | 1 | +0.000 |
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
| volume fraction | - | 0.008 | 0.00798119 | -0.235 | 6 | yes | 0.0079583 | -0.521 |
| R at distribution peak | A | 144.118 | 152.278 | +5.662 | 10 | yes | 147.5 | +2.346 |
| RMS radius (pyIrena 'Rg') | A | 156.122 | 156.277 | +0.100 | 6 | yes | 156.43 | +0.198 |
| reduced chi^2 | - | 1 | 1.00062 | +0.062 |  | - | 1 | +0.000 |
| model vs exact curve: median abs. deviation | % | 0 | 0.0889834 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 1.10375 |  | 15 | yes |  |  |

> regularised inversions broaden the distribution, so the modal radius has a wider tolerance than the integral moments

> fitted model compared point-by-point with <name>_ideal.dat

## Unified Fit

### `unified_one_level`

*pyIrena settings — use the same in Irena:* 1 level(s), fit G/Rg/B/P + flat background, TRF, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| G | 1/cm | 100 | 99.9498 | -0.050 | 3 | yes | 99.95 | -0.050 |
| Rg | A | 200 | 199.844 | -0.078 | 2 | yes | 199.8 | -0.100 |
| P | - | 4 | 4.00923 | +0.231 | 2 | yes | 4.009 | +0.225 |
| B | cm^-1 A^-P | 3.04500e-07 | 2.95171e-07 | -3.064 | 25 | yes | 2.95e-7 | -3.120 |
| background | 1/cm | 0.01 | 0.01006 | +0.600 | 15 | yes | 0.01006 | +0.600 |
| reduced chi^2 | - | 1 | 0.979006 | -2.099 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.0828329 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.589873 |  | 12 | yes |  |  |

> B and P are strongly correlated (B carries units cm^-1 A^-P), so B has a wide tolerance whenever P is a free parameter

> 1.0 expected for a correct model and correct error bars

> fitted model compared point-by-point with <name>_ideal.dat

### `unified_two_level`

*pyIrena settings — use the same in Irena:* 2 level(s), fit G/Rg/B/P + flat background, TRF, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| G_1 | 1/cm | 8 | 7.96186 | -0.477 | 3 | yes | 7.9617 | -0.479 |
| Rg_1 | A | 120 | 119.829 | -0.143 | 2 | yes | 119.83 | -0.142 |
| P_1 | - | 4 | 3.93526 | -1.619 | 2 | yes | 3.9354 | -1.615 |
| B_1 | cm^-1 A^-P | 1.88000e-07 | 2.32125e-07 | +23.471 | 25 | yes | 2.3206e-7 | +23.436 |
| G_2 | 1/cm | 4000 | 4001.17 | +0.029 | 3 | yes | 4001. | +0.025 |
| Rg_2 | A | 1200 | 1202 | +0.167 | 2 | yes | 1202 | +0.167 |
| P_2 | - | 3.2 | 3.18946 | -0.329 | 2 | yes | 3.1894 | -0.331 |
| B_2 | cm^-1 A^-P | 1.39300e-06 | 1.47601e-06 | +5.959 | 25 | yes | 1.4765e-6 | +5.994 |
| background | 1/cm | 0.02 | 0.0201604 | +0.802 | 15 | yes | 0.02016 | +0.800 |
| reduced chi^2 | - | 1 | 1.06538 | +6.538 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.153888 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 1.71387 |  | 12 | yes |  |  |

> B and P are strongly correlated (B carries units cm^-1 A^-P), so B has a wide tolerance whenever P is a free parameter

> 1.0 expected for a correct model and correct error bars

> fitted model compared point-by-point with <name>_ideal.dat

## Unified Fit (slit smearing)

### `unified_one_level_slitsmeared`

*pyIrena settings — use the same in Irena:* 1 level(s), fit G/Rg/B/P + flat background, TRF, start = 0.60x truth, slit smearing dQl=0.03

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| G | 1/cm | 100 | 99.8558 | -0.144 | 3 | yes |  |  |
| Rg | A | 200 | 200.098 | +0.049 | 2 | yes |  |  |
| P | - | 4 | 4.00444 | +0.111 | 2 | yes |  |  |
| B | cm^-1 A^-P | 3.04500e-07 | 3.00013e-07 | -1.474 | 25 | yes |  |  |
| background | 1/cm | 0.01 | 0.00996066 | -0.393 | 15 | yes |  |  |
| reduced chi^2 | - | 1 | 0.930125 | -6.988 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.181945 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.412444 |  | 12 | yes |  |  |

> B and P are strongly correlated (B carries units cm^-1 A^-P), so B has a wide tolerance whenever P is a free parameter

> 1.0 expected for a correct model and correct error bars

> fitted model compared point-by-point with <name>_ideal.dat

## Modeling

### `modeling_sizedist_plus_unified`

*pyIrena settings — use the same in Irena:* populations: size_dist (contrast 100e20 cm^-4), unified_level; local (TRF) fit, size-distribution grid 200 bins, contrast held fixed at the value quoted above, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| pop1 mean_size (log-normal median) | A | 30 | 30.0419 | +0.140 | 4 | yes | 31.34 | +4.467 |
| pop1 sdeviation | - | 0.3 | 0.316207 | +5.402 | 6 | yes | 0.322 | +7.333 |
| pop1 scale | - | 0.005 | 0.00500931 | +0.186 | 6 | yes | 0.004947 | -1.060 |
| pop2 G | 1/cm | 20000 | 19996 | -0.020 | 4 | yes | 19999.7 | -0.001 |
| pop2 Rg | A | 4000 | 3998.8 | -0.030 | 4 | yes | 4000.4 | +0.010 |
| pop2 P | - | 3.5 | 3.49238 | -0.218 | 3 | yes | 3.503 | +0.086 |
| pop2 B | cm^-1 A^-P | 1.56345e-08 | 1.64682e-08 | +5.333 | 25 | yes | 1.551e-8 | -0.796 |
| background | 1/cm | 0.01 | 0.00994733 | -0.527 | 30 | yes | 0.00998 | -0.200 |
| reduced chi^2 | - | 1 | 1.03921 | +3.921 |  | - | 1.0535 | +5.350 |
| model vs exact curve: median abs. deviation | % | 0 | 0.213211 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 2.16483 |  | 15 | yes |  |  |

> a small flat background is a weakly determined nuisance parameter

> fitted model compared point-by-point with <name>_ideal.dat

### `modeling_sizedist_plus_peak`

*pyIrena settings — use the same in Irena:* populations: size_dist (contrast 80e20 cm^-4), diffraction_peak; local (TRF) fit, size-distribution grid 200 bins, contrast held fixed at the value quoted above, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| pop1 mean_size (log-normal median) | A | 60 | 60.1816 | +0.303 | 4 | yes | 59.119 | -1.468 |
| pop1 sdeviation | - | 0.25 | 0.260352 | +4.141 | 6 | yes | 0.261 | +4.400 |
| pop1 scale | - | 0.004 | 0.00400364 | +0.091 | 6 | yes | 0.003955 | -1.125 |
| pop2 amplitude | 1/cm | 5 | 4.99144 | -0.171 | 5 | yes | 4.992 | -0.160 |
| pop2 position | 1/A | 0.15 | 0.149981 | -0.013 | 1 | yes | 0.1494 | -0.400 |
| pop2 width (sigma) | 1/A | 0.012 | 0.0120011 | +0.009 | 6 | yes | 0.01413 | +17.750 |
| background | 1/cm | 0.005 | 0.00501291 | +0.258 | 30 | yes | 0.005008 | +0.160 |
| reduced chi^2 | - | 1 | 1.18594 | +18.594 |  | - | 1.131 | +13.100 |
| model vs exact curve: median abs. deviation | % | 0 | 0.335414 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 2.45532 |  | 15 | yes |  |  |

> a small flat background is a weakly determined nuisance parameter

> fitted model compared point-by-point with <name>_ideal.dat

### `modeling_hardsphere`

*pyIrena settings — use the same in Irena:* populations: size_dist (contrast 100e20 cm^-4); local (TRF) fit, size-distribution grid 200 bins, contrast held fixed at the value quoted above, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| pop1 mean_size (log-normal median) | A | 50 | 49.983 | -0.034 | 4 | yes | 49.222 | -1.556 |
| pop1 sdeviation | - | 0.12 | 0.123946 | +3.288 | 6 | yes | 0.1248 | +4.000 |
| pop1 scale | - | 0.2 | 0.200044 | +0.022 | 6 | yes | 0.19763 | -1.185 |
| hard-sphere radius | A | 50 | 49.9432 | -0.114 | 6 | yes | 49.9542 | -0.092 |
| hard-sphere volume fraction | - | 0.25 | 0.249643 | -0.143 | 10 | yes | 0.2497 | -0.120 |
| background | 1/cm | 0.01 | 0.0093383 | -6.617 | 30 | yes | 0.00936 | -6.400 |
| reduced chi^2 | - | 1 | 1.15067 | +15.067 |  | - | 1.1246 | +12.460 |
| model vs exact curve: median abs. deviation | % | 0 | 0.103617 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 3.05998 |  | 15 | yes |  |  |

> a small flat background is a weakly determined nuisance parameter

> fitted model compared point-by-point with <name>_ideal.dat

### `modeling_mass_fractal`

*pyIrena settings — use the same in Irena:* populations: mass_fractal (contrast 100e20 cm^-4); local (TRF) fit, size-distribution grid 200 bins, contrast held fixed at the value quoted above, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| Phi | - | 0.005 | 0.004999 | -0.020 | 6 | yes |  |  |
| Radius | A | 40 | 40.019 | +0.047 | 4 | yes |  |  |
| Dv | - | 2.4 | 2.39923 | -0.032 | 3 | yes |  |  |
| Ksi | A | 800 | 800.352 | +0.044 | 8 | yes |  |  |
| background | 1/cm | 0.01 | 0.0102148 | +2.148 | 30 | yes |  |  |
| reduced chi^2 | - | 1 | 1.09383 | +9.383 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.072548 |  | 3 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 2.27086 |  | 15 | yes |  |  |

> a small flat background is a weakly determined nuisance parameter

> fitted model compared point-by-point with <name>_ideal.dat

### `modeling_surface_fractal`

*pyIrena settings — use the same in Irena:* populations: surface_fractal (contrast 100e20 cm^-4); local (TRF) fit, size-distribution grid 200 bins, contrast held fixed at the value quoted above, start = 0.60x truth

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

*pyIrena settings — use the same in Irena:* populations: guinier_porod; local (TRF) fit, size-distribution grid 200 bins, contrast held fixed at the value quoted above, start = 0.60x truth

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
| I0 | 1/cm | 250 | 250.686 | +0.275 | 3 | yes | 250.77 | +0.308 |
| Rg | A | 45 | 45.9727 | +2.162 | 3 | yes | 46.065 | +2.367 |
| reduced chi^2 | - | 1 | 0.871892 | -12.811 |  | - | 0.88 | -12.000 |
| model vs exact curve: median abs. deviation | % | 0 | 0.240806 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 2.35248 |  | 12 | yes |  |  |

> fitted model compared point-by-point with <name>_ideal.dat

### `simple_guinier_rod`

*pyIrena settings — use the same in Irena:* model 'Guinier Rod', all listed parameters free, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| I0 | cm^-1 A^-1 | 5 | 5.0042 | +0.084 | 3 | yes | 5.0044 | +0.088 |
| Rc | A | 25 | 25.072 | +0.288 | 3 | yes | 25.069 | +0.276 |
| reduced chi^2 | - | 1 | 1.04912 | +4.912 |  | - | 1.05 | +5.000 |
| model vs exact curve: median abs. deviation | % | 0 | 0.0742042 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 1.06439 |  | 12 | yes |  |  |

> fitted model compared point-by-point with <name>_ideal.dat

### `simple_guinier_sheet`

*pyIrena settings — use the same in Irena:* model 'Guinier Sheet', all listed parameters free, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| I0 | cm^-1 A^-2 | 0.02 | 0.0200201 | +0.100 | 3 | yes | 0.02002 | +0.100 |
| Rg | A | 15 | 15.0031 | +0.021 | 3 | yes | 15 | +0.000 |
| reduced chi^2 | - | 1 | 0.988868 | -1.113 |  | - | 0.99 | -1.000 |
| model vs exact curve: median abs. deviation | % | 0 | 0.0910625 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.0995459 |  | 12 | yes |  |  |

> fitted model compared point-by-point with <name>_ideal.dat

### `simple_porod`

*pyIrena settings — use the same in Irena:* model 'Porod', all listed parameters free, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| Kp | cm^-1 A^-4 | 2.50000e-06 | 2.51216e-06 | +0.486 | 3 | yes | 2.5137e-06 | +0.548 |
| Background | 1/cm | 0.05 | 0.0499857 | -0.029 | 3 | yes | 0.049994 | -0.012 |
| reduced chi^2 | - | 1 | 1.06544 | +6.544 |  | - | 1.07 | +7.000 |
| model vs exact curve: median abs. deviation | % | 0 | 0.0280001 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.42908 |  | 12 | yes |  |  |

> fitted model compared point-by-point with <name>_ideal.dat

### `simple_power_law`

*pyIrena settings — use the same in Irena:* model 'Power Law', all listed parameters free, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| Prefactor | cm^-1 A^-n | 0.0012 | 0.0012047 | +0.391 | 3 | yes | 0.0012048 | +0.400 |
| Exponent | - | 3.2 | 3.19894 | -0.033 | 1 | yes | 3.1989 | -0.034 |
| Background | 1/cm | 0.02 | 0.0195288 | -2.356 | 3 | yes | 0.019507 | -2.465 |
| reduced chi^2 | - | 1 | 0.933632 | -6.637 |  | - | 0.94 | -6.000 |
| model vs exact curve: median abs. deviation | % | 0 | 0.079124 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.952092 |  | 12 | yes |  |  |

> fitted model compared point-by-point with <name>_ideal.dat

### `simple_sphere`

*pyIrena settings — use the same in Irena:* model 'Sphere', all listed parameters free, start = 0.60x truth, complex background (flat term only)

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| Scale | 1/cm | 8 | 7.99862 | -0.017 | 3 | yes | 7.9986 | -0.018 |
| R | A | 120 | 119.965 | -0.029 | 3 | yes | 119.96 | -0.033 |
| BG_flat | 1/cm | 0.01 | 0.0100646 | +0.646 | 3 | yes | 0.010049 | +0.490 |
| reduced chi^2 | - | 1 | 1.01504 | +1.504 |  | - | 1.002 | +0.200 |
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
| model vs exact curve: max abs. deviation | % | 0 | 0.141041 |  | 12 | yes |  |  |

> fitted model compared point-by-point with <name>_ideal.dat

### `simple_debye_bueche`

*pyIrena settings — use the same in Irena:* model 'Debye-Bueche', all listed parameters free, start = 0.60x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| Prefactor | cm^-1 A^-3 | 1 | 0.782174 | -21.783 |  | - |  |  |
| Eta | - | 0.05 | 0.056567 | +13.134 |  | - |  |  |
| CorrLength | A | 80 | 79.951 | -0.061 | 3 | yes |  |  |
| Prefactor * Eta^2 | cm^-1 A^-3 | 0.0025 | 0.00250282 | +0.113 | 3 | yes |  |  |
| reduced chi^2 | - | 1 | 0.967369 | -3.263 |  | - |  |  |
| model vs exact curve: median abs. deviation | % | 0 | 0.10529 |  | 2 | yes |  |  |
| model vs exact curve: max abs. deviation | % | 0 | 0.173262 |  | 12 | yes |  |  |

> degenerate with the other prefactor - only the product is determined by the data

> the combination the data actually determine

> fitted model compared point-by-point with <name>_ideal.dat

## WAXS Peak Fit

### `waxs_three_gauss`

*pyIrena settings — use the same in Irena:* 3 x Gauss peaks, Linear background, all parameters free; A and background start at 0.60x truth, Q0 within 1 %, FWHM 1.30x truth

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| bg0 | - | 8 | 7.87062 | -1.617 | 10 | yes | * |  |
| bg1 | - | 1.5 | 1.54147 | +2.765 | 10 | yes | * |  |
| background at Q = 1.2 | 1/cm | 9.8 | 9.72039 | -0.812 | 5 | yes |  |  |
| background at Q = 1.9 | 1/cm | 10.85 | 10.7994 | -0.466 | 5 | yes |  |  |
| background at Q = 2.6 | 1/cm | 11.9 | 11.8784 | -0.181 | 5 | yes |  |  |
| background at Q = 3.3 | 1/cm | 12.95 | 12.9575 | +0.058 | 5 | yes |  |  |
| background at Q = 4 | 1/cm | 14 | 14.0365 | +0.261 | 5 | yes |  |  |
| peak1_A | 1/cm | 100 | 99.5683 | -0.432 | 3 | yes | 99.4159 | -0.584 |
| peak1_Q0 | 1/A | 1.85 | 1.85011 | +0.006 | 0.5 | yes | 1.85019 | +0.010 |
| peak1_FWHM | 1/A | 0.09 | 0.0902364 | +0.263 | 3 | yes | 0.09144 | +1.600 |
| peak1_area | cm^-1 A^-1 | 9.5802 | 9.5639 | -0.170 | 4 | yes |  |  |
| peak2_A | 1/cm | 60 | 59.7818 | -0.364 | 3 | yes | 59.7985 | -0.336 |
| peak2_Q0 | 1/A | 2.55 | 2.55047 | +0.018 | 0.5 | yes | 2.55028 | +0.011 |
| peak2_FWHM | 1/A | 0.12 | 0.120372 | +0.310 | 3 | yes | 0.122968 | +2.473 |
| peak2_area | cm^-1 A^-1 | 7.66416 | 7.65996 | -0.055 | 4 | yes |  |  |
| peak3_A | 1/cm | 35 | 34.9834 | -0.047 | 3 | yes | 34.9566 | -0.124 |
| peak3_Q0 | 1/A | 3.05 | 3.04892 | -0.035 | 0.5 | yes | 3.04925 | -0.025 |
| peak3_FWHM | 1/A | 0.15 | 0.150923 | +0.615 | 3 | yes | 0.1555942 | +3.729 |
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
| bg0 | - | 10 | 11.9978 | +19.978 |  | - | * |  |
| bg1 | - | -2 | -4.51982 | +125.991 |  | - | * |  |
| bg2 | - | 0.8 | 1.84821 | +131.027 |  | - | * |  |
| bg3 | - | -0.05 | -0.189003 | +278.006 |  | - | * |  |
| background at Q = 1.5 | 1/cm | 8.63125 | 8.7387 | +1.245 | 5 | yes |  |  |
| background at Q = 2.12 | 1/cm | 8.88271 | 8.92544 | +0.481 | 5 | yes |  |  |
| background at Q = 2.75 | 1/cm | 9.51016 | 9.61476 | +1.100 | 5 | yes |  |  |
| background at Q = 3.38 | 1/cm | 10.4403 | 10.5298 | +0.857 | 5 | yes |  |  |
| background at Q = 4 | 1/cm | 11.6 | 11.3937 | -1.778 | 5 | yes |  |  |
| peak1_A | 1/cm | 200 | 197.484 | -1.258 | 3 | yes |  |  |
| peak1_Q0 | 1/A | 2.1 | 2.09964 | -0.017 | 0.5 | yes | 2.09957 | -0.020 |
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
| formula weight | g/mol | 60.083 | 60.083 | +0.000 | 0.2 | yes | 60.08 | -0.005 |
| X-ray SLD | 1e10 cm^-2 | 18.6412 | 18.6412 | +0.000 | 0.5 | yes | 18.64 | -0.007 |
| neutron SLD | 1e10 cm^-2 | 3.47411 | 3.47477 | +0.019 | 1 | yes | 3.474 | -0.003 |

> amorphous silica

### `contrast: H2O (1 g/cm3)`

*pyIrena settings — use the same in Irena:* free-electron X-ray SLD and bound-coherent neutron SLD; reference values from IUPAC 2021 atomic weights and Sears (1992) scattering lengths

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| formula weight | g/mol | 18.015 | 18.015 | +0.000 | 0.2 | yes | 18.02 | +0.028 |
| X-ray SLD | 1e10 cm^-2 | 9.41995 | 9.41995 | +0.000 | 0.5 | yes | 9.42 | +0.001 |
| neutron SLD | 1e10 cm^-2 | -0.559927 | -0.560963 | +0.185 | 1 | yes | -0.5606 | +0.120 |

> light water

### `contrast: D2O (1.107 g/cm3)`

*pyIrena settings — use the same in Irena:* free-electron X-ray SLD and bound-coherent neutron SLD; reference values from IUPAC 2021 atomic weights and Sears (1992) scattering lengths

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| formula weight | g/mol | 20.027 | 20.0272 | +0.001 | 0.2 | yes | 20.03 | +0.015 |
| X-ray SLD | 1e10 cm^-2 | 9.38025 | 9.38015 | -0.001 | 0.5 | yes | 9.38 | -0.003 |
| neutron SLD | 1e10 cm^-2 | 6.37291 | 6.37115 | -0.028 | 1 | yes | 6.375 | +0.033 |

> heavy water

### `contrast: Al2O3 (3.97 g/cm3)`

*pyIrena settings — use the same in Irena:* free-electron X-ray SLD and bound-coherent neutron SLD; reference values from IUPAC 2021 atomic weights and Sears (1992) scattering lengths

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| formula weight | g/mol | 101.961 | 101.96 | -0.001 | 0.2 | yes | 102.0 | +0.038 |
| X-ray SLD | 1e10 cm^-2 | 33.0376 | 33.0379 | +0.001 | 0.5 | yes | 33.04 | +0.007 |
| neutron SLD | 1e10 cm^-2 | 5.69953 | 5.70007 | +0.010 | 1 | yes | 5.698 | -0.027 |

> corundum

### `contrast: Fe (7.874 g/cm3)`

*pyIrena settings — use the same in Irena:* free-electron X-ray SLD and bound-coherent neutron SLD; reference values from IUPAC 2021 atomic weights and Sears (1992) scattering lengths

| quantity | unit | true | pyIrena | dev % | tol % | in tol? | **Irena** | **Irena dev %** |
|---|---|---|---|---|---|---|---|---|
| formula weight | g/mol | 55.845 | 55.845 | +0.000 | 0.2 | yes | 55.85 | +0.009 |
| X-ray SLD | 1e10 cm^-2 | 62.211 | 62.211 | +0.000 | 0.5 | yes | 62.21 | -0.002 |
| neutron SLD | 1e10 cm^-2 | 8.02405 | 8.02405 | +0.000 | 1 | yes | 8.1 | +0.946 |

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

---

## Notes on the pyIrena-Irena comparison

### Quantities Irena defines differently

**Modeling diffraction-peak width.**  pyIrena reports the Gaussian standard
deviation sigma; Irena reports the half-width at half-maximum.  This is the only
row in the table where the two packages differ by more than 8 %, and it is
entirely a difference in the reported quantity.  On
`modeling_sizedist_plus_peak` the ratio of the two reported widths is
0.01413 / 0.012 = **1.17750**, against sqrt(2 ln 2) = **1.17741** for a Gaussian
— agreement to five digits.  Converted to sigma the Irena width is 0.012001,
i.e. **+0.008 %** from the true value rather than the +17.8 % the raw numbers
suggest.

**WAXS polynomial background coefficients** (`bg0`, `bg1`, … marked `*`).  The
two packages parameterise the peak-fit background differently, so the
coefficients cannot be compared term by term.  The quantity the data actually
determine is the background *curve*, which is tabulated at five Q values in each
WAXS section and agrees with the truth to better than 2 % in pyIrena.

### Getting a like-for-like comparison

Three settings have to match before the parameters mean the same thing in both
packages.  Each was found the hard way while this table was assembled, and each
produces a large apparent disagreement while leaving the fit quality untouched
— so a good reduced chi-squared on both sides is not on its own evidence that
the two analyses are comparable.

1. **The contrast.** Each dataset's contrast is stated in the `#` header of its
   `.dat` file and in the settings line of its section here; it is not the same
   for every dataset (`modeling_sizedist_plus_peak` uses 80, its neighbours
   100).  Since the fitted scale absorbs the contrast, using the wrong one
   shifts the recovered volume fraction and nothing else.
2. **Cut-off radii in the Unified model.** `unified_two_level` links `RgCO` on
   level 2 to `Rg` of level 1.  Fitting the same data without that cut-off gives
   a model that describes the curve just as well (0.11 % median deviation from
   the exact curve, against 0.14 % with it) but with visibly different `G_1`,
   `B_1`, `B_2` and `P` values.
3. **Width conventions**, as above.

### The remaining rows that deviate by more than a percent or two

Every one of these is a quantity the data determine only weakly.  They are worth
listing because in each case the *pattern* of the disagreement is informative.

**Power-law prefactors.**  `B` carries units cm^-1 A^-P, so it is strongly
correlated with `P`: a 1.6 % shift in `P_1` moves `B_1` by more than 20 %.  In
`unified_two_level` both packages are off the true `B_1` by the same +23.4 %,
and off the true `B_2` by the same +6.0 % — they agree with **each other** to
0.03 % on both.  Two independent implementations landing on the same value of an
ill-determined parameter is a stronger result than either agreeing with the
truth would be.  The `pop2 B` of `modeling_sizedist_plus_unified` is the one
prefactor where the packages differ from each other (+5.3 % against -0.8 %), for
the same reason: their fitted `P` differ by 0.3 %.

**Modal radii from the Size Distribution tool.**  `sizes_sphere_bimodal` mode 1
(+2.3 % against -4.7 %) and `sizes_spheroid_lognormal` (+5.7 % against +2.3 %)
are the two largest.  A regularised inversion returns a smoothed distribution, so
where its peak sits depends on the regularisation, which the two packages do not
tune identically.  The integral moments, which do not depend on the smoothing,
agree to a few tenths of a percent throughout.

**The log-normal shape in `modeling_sizedist_plus_unified`.**  Both packages put
`sdeviation` high (+5.4 % and +7.3 %) and Irena also puts the median high
(+4.5 %).  This population's Guinier region sits on the Porod tail of the Unified
level, so its width and median trade against that level's `B` and `P`; the two
packages are biased in the same direction, which is what a shared degeneracy
looks like rather than an implementation difference.

**Reduced chi-squared** is a property of a fit, not a parameter of the model, and
the two packages weight and count degrees of freedom slightly differently.  The
values are listed for information; a difference of a few percent there carries no
meaning.

### One value read against the table

**`simple_porod`, `Kp`** was entered as `2.5137`, three decades from the true
`2.5e-06`, and has been read as `2.5137e-06` (deviation +0.548 %).  Two entries
in `unified_two_level` (`B_1`, `B_2`) arrived as `2.3206-7` and `1.4765-6`, with
the exponent's `e` lost in transcription, and have been rewritten as `2.3206e-7`
and `1.4765e-6`.

### Provenance of the numbers

The pyIrena column was produced by `validationData/run_validation_report.py`
against pyIrena **1.1.0b9**.  One row depends on the version: the hard-sphere
structure-factor parameters in the Modeling tool were not actually varied during
fitting before 1.1.0b9 (they were packed into the fit vector under a key group
the surface-fractal population already used, so each optimiser step wrote them to
the wrong place).  On `modeling_hardsphere` that defect gave reduced chi-squared
397 with the structure-factor and distribution parameters wrong by 8-40 %; the
row in this table is the fixed behaviour.

### Reproducing this table

```bash
python3 validationData/generate_validation_data.py    # rebuild the data
python3 validationData/run_validation_report.py       # fit it all, write the tables
```

The pyIrena columns come from that second command.  The Irena columns are read
from `validationData/irena_values.csv`, which holds the values obtained with the
Igor package, keyed by dataset and quantity; a quantity with no entry there is
left blank.  Adding more Irena results means adding rows to that CSV and
re-running the report — the numbers are data in the repository, not text inside
a document, so regenerating the report never loses them.

`validationData/fill_irena_deviations.py` fills the `Irena dev %` column of a
hand-edited copy of this table, and with `--export-csv` writes that copy's Irena
column back into `irena_values.csv`.

This prose section is maintained by hand in `validationData/irena_notes.md` and
is appended to the report verbatim.

