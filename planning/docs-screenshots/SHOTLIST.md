# Documentation screenshots: shot list

One row per screenshot in `docs/tools/`. The ID is `<page>-<nn>-<what>`; the image lives at
`docs/tools/img/<page>/<id>.png`. Pages are in the order of the Help buttons
(`pyirena/_docs.py`), then pages without a Help button.

- **automated: yes**: `python scripts/make_doc_screenshots.py --only <id>` makes it (see
  `docs/developer/testing.md`, "Updating documentation screenshots").
- **automated: no**: needs judgment (real data, a tricky state, an OS dialog). The page holds an
  invisible `% SCREENSHOT <id>: …` comment until the image exists; then replace the comment with
  a `figure` that has alt text and a caption, and set the status to `done`.
- **status**: `done` (image committed and referenced), `todo` (placeholder only).

`pyirena/tests/test_doc_screenshots.py` checks that this list and the pages agree.

Column order: ID · page · dataset · state to set up · must show · automated? · status.

## unified-fit

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
| unified-fit-01-window | unified-fit | validationData/unified/unified_two_level.h5 | 2 levels, rough start values, RgCO of level 2 linked, Fit; callouts 1–13 | whole window, fitted model, flat residuals, quality message | yes | done |
| unified-fit-02-porod-tab | unified-fit | validationData/unified/unified_two_level.h5 | data loaded, Porod (I·Q⁴) tab selected | two plateaus/transitions | yes | done |
| unified-fit-03-guinier-fit | unified-fit | validationData/unified/unified_two_level.h5 | 2 levels, cursors Q 0.0002–0.0012, level 2 tab, Fit Rg/G btwn cursors | green local Guinier curve, updated G/Rg, message | yes | done |
| unified-fit-04-fit-result | unified-fit | validationData/unified/unified_two_level.h5 | as 01, then Results to graphs | graph area with parameter boxes and residuals | yes | done |
| unified-fit-05-slit-smeared | unified-fit | validationData/unified/unified_one_level_slitsmeared.h5 | 1 level, Fit (slit smearing switches on from dQl) | Slit smeared checked, SL 0.03, fit through smeared data | yes | done |
| unified-fit-06-real-data | unified-fit | testData/ (Jan picks a sample) | good multi-level fit, Results to graphs | a real multi-level curve and its fit | no | todo |

## size-distribution

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
| size-distribution-01-window | size-distribution | validationData/sizes/sizes_sphere_lognormal.h5 | contrast 100, grid 10–1000 Å, regularization, Fit Sizes; callouts 1–16 | whole window, model on data, P(r) peak near 100 Å, results fields | yes | done |
| size-distribution-02-background-tab | size-distribution | validationData/sizes/sizes_sphere_flat_background.h5 | Background tab, flat background Q 0.2–0.35, Fit Background | background line and blue corrected data | yes | done |
| size-distribution-03-fit-result | size-distribution | validationData/sizes/sizes_sphere_flat_background.h5 | as 02, cursors Q 0.0008–0.15, Fit Sizes, Results to graph | graph area: annotation, residuals, P(r) peak near 114 Å | yes | done |
| size-distribution-04-real-data | size-distribution | testData/ (Jan picks a sample) | regularization or MaxEnt fit with power-law and flat background, Results to graph | a real distribution and its fit | no | todo |

## modeling

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
| modeling-01-window | modeling | validationData/modeling/modeling_sizedist_plus_unified.h5 | window 1500×1000; P1 log-normal spheres, P2 Unified level, rough starts, background fitted, Standard Fit; callouts 1–15 | whole window, model on data, residuals, P1 distribution, Fit Results | yes | done |
| modeling-02-population-tab | modeling | validationData/modeling/modeling_sizedist_plus_unified.h5 | as 01 before fitting, P1 tab, Graph Model | P1 settings and the starting model off the data | yes | done |
| modeling-03-fit-result | modeling | validationData/modeling/modeling_sizedist_plus_unified.h5 | as 01, Show individual? on with P2 tab, Results to graph | graph area: P2 curve, derived-values box, residuals, distribution | yes | done |
| modeling-04-real-data | modeling | testData/ (Jan picks a sample) | multi-population fit on real data, Show individual? on | a real multi-population fit | no | todo |

## simple-fits

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
| simple-fits-01-window | simple-fits | validationData/simple_fits/simple_guinier.h5 | window 1200×900; Guinier, I0 100, Rg 30 start, default cursors, Fit; callouts 1–14 | whole window, model on data, residuals, Guinier linearization | yes | done |
| simple-fits-02-sphere | simple-fits | validationData/simple_fits/simple_sphere.h5 | Sphere, complex background with B = 0 and P fixed, flat fitted, Fit | graph area: model through the form-factor oscillations | yes | done |
| simple-fits-03-invariant | simple-fits | validationData/sizes/sizes_sphere_flat_background.h5 | Invariant, contrast 100, complex background B = 0 fixed, Fit Flat btwn cursors at Q 0.2–0.35, cursors 0.0008–0.35, Calculate Invariant | corrected data, background line, running integral, φ in status line | yes | done |
| simple-fits-04-real-data | simple-fits | testData/ (Jan picks a sample) | Guinier or Porod fit on real data | linearization with the fitted range highlighted | no | todo |

## fractals

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
| fractals-01-window | fractals | validationData/modeling/modeling_mass_fractal.h5 (copy with a stored, RgCO-linked two-level Unified fit) | NeXus loaded, Grow One (Z 250, seed 12345); PyVista off (offscreen); callouts 1–11 | whole window, data + Unified + green analytical curve, 2D slice, parameters | yes | done |
| fractals-02-montecarlo-iq | fractals | as 01 | as 01, then Compute Monte Carlo I(Q) (MC seed pinned) | I(Q) plot with blue dashed MC curve | yes | done |
| fractals-03-optimizer | fractals | as 01 | Z 100, seed 12345, Optimizer: dmin 1.3, c 1.4, 2 iterations, Find Best Growth | controls panel: six trials, ★ BEST selected, green values vs targets | yes | done |
| fractals-04-3d-viewer | fractals | any grown aggregate (Z 250, seed 12345) | real display, PyVista installed | 3D isosurface in the right pane | no | todo |
| fractals-05-targets | fractals | file with a Unified fractal level pair | Open NeXus…; blocked until the panel finds the pair (RgCutoff/RgCO key) | target summary + colour-coded parameters | no | todo |

## waxs-peak-fit

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
| waxs-peak-fit-01-window | waxs-peak-fit | validationData/waxs/waxs_three_gauss.h5 | Background Linear, Find Peaks, remove peaks < 20 % of tallest, Fit; callouts 1–11 | whole window, 3 fitted Gaussians, background, residuals, quality line | yes | done |
| waxs-peak-fit-02-find-peaks | waxs-peak-fit | validationData/waxs/waxs_three_gauss.h5 | Background Linear, Find Peaks (defaults) | detected peaks incl. false detections | yes | done |
| waxs-peak-fit-03-fit-result | waxs-peak-fit | validationData/waxs/waxs_three_gauss.h5 | as 01, then Results to graphs | graph area with annotation | yes | done |
| waxs-peak-fit-04-diffraction-lines | waxs-peak-fit | validationData/waxs/waxs_three_gauss.h5 + testData/Si.cif | as 01, Diffraction Lines tab, Si.cif added with hkl | CIF row, sticks with hkl labels | yes | done |
| waxs-peak-fit-05-phase-id | waxs-peak-fit | testData/ real powder pattern + matching CIFs (Jan picks) | two phases imported, sticks aligned with peaks | a real phase identification | no | todo |

## carbon-fit

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|

## saxs-morph

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
| saxs-morph-01-window | saxs-morph | validationData/simple_fits/simple_debye_bueche.h5 | box 4000 Å, seed 42, φ 0.3 (derive Δρ²), cursors Q 0.005–0.2, Calculate 3D; 3D pane covered (PyVista off offscreen); callouts 1–12 | whole window, model on data, slice, result block | yes | done |
| saxs-morph-02-calculated | saxs-morph | validationData/simple_fits/simple_debye_bueche.h5 | as 01; I(Q) plot only | red model following the data, Q_box and Q_nyq markers | yes | done |
| saxs-morph-03-slice | saxs-morph | validationData/simple_fits/simple_debye_bueche.h5 | as 01; 2D slice viewer only | bicontinuous slice, slice label | yes | done |
| saxs-morph-04-3d-view | saxs-morph | validationData/simple_fits/simple_debye_bueche.h5 | as 01 in the real GUI with PyVista; rotate the isosurface | 3D bicontinuous isosurface with bounding box | no | todo |

## data-merge

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
| data-merge-01-window | data-merge | validationData/merge/merge_usaxs.h5 + merge_saxs.h5 | DS1/DS2 folders USAXS/ and SAXS/ (copies), cursors Q 0.01–0.04, Scale: DS2, Optimize Merge; callouts 1–11 | whole window, merged curve, scale ≈ 0.872 in the field and status | yes | done |
| data-merge-02-loaded | data-merge | validationData/merge/merge_usaxs.h5 + merge_saxs.h5 | both loaded, cursors Q 0.01–0.04, before Optimize | both branches and the overlap cursors | yes | done |
| data-merge-03-porod-real-data | data-merge | testData/ (Jan picks a real USAXS + SAXS pair) | Porod display mode, before and/or after a good merge | mismatch as a step/slope in the flat I·Q⁴ line, then one flat line | no | todo |

## data-manipulation

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
| data-manipulation-01-window | data-manipulation | validationData/manipulation/manipulate_input.h5 (+ manipulate_reference.h5 in the same folder) | folder "manipulation" (copies), Scale tab, input double-clicked, Scale I 0.4, Background 0.12; callouts 1–8 | whole window, input and green result | yes | done |
| data-manipulation-02-scale-result | data-manipulation | validationData/manipulation/manipulate_input.h5 | as 01; graph only | input and the scaled result (= reference) | yes | done |
| data-manipulation-03-subtract | data-manipulation | validationData/manipulation/manipulate_input.h5 + manipulate_reference.h5 | Subtract tab, reference set as buffer, input double-clicked, Buffer scale 2.5 | sample, buffer and a flat 0.3 cm⁻¹ result | yes | done |
| data-manipulation-04-similarity | data-manipulation | testData/ (Jan picks a real exposure series with radiation damage) | Average tab, frames selected, Check Similarity | coloured frames and the accepted/rejected results table | no | todo |

## scattering-contrast

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
| scattering-contrast-01-window | scattering-contrast | none: magnetite (Fe3O4, 5.17) in water (H2O, 1.0) | both compounds set, Parse, Calculate at 12 keV; window 1120×1140; callouts 1–11 | whole window, filled results table, status line | yes | done |
| scattering-contrast-02-energy-scan | scattering-contrast | none: magnetite in water | Energy Scan 6.5–8.0 keV, 300 points; graph window 900×760 | contrast dip and absorption jump at the Fe K edge | yes | done |
| scattering-contrast-03-import-dialog | scattering-contrast | a compounds file with several entries (Jan prepares) | Import from File… on that file | the Load Selected into Slot / Add All to Local Library dialog | no | todo |

## data-explorer

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
| data-explorer-01-window | data-explorer | scratch folder `unified`: copies of validationData/unified/*.h5, each Unified-fitted and saved | folder open, all 3 files selected, unified_two_level.h5 tree with unified_fit_results expanded, 1D Graph tab; callouts 1–10 | three panels, results group in the tree, presets | yes | done |
| data-explorer-02-graph | data-explorer | as 01 | NXcanSAS + Unified Fit model, All selected files, New Graph; graph 760×560 | three curves with models, toolbar | yes | done |
| data-explorer-03-collect | data-explorer | as 01 | Collect Values: Unified Fit, Rg, level 1, file order; Collect | table and scatter plot of Rg | yes | done |

## data-selector

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
| data-selector-01-window | data-selector | /tmp/pyirena-docs/unified: copies of validationData/unified/*.h5, each Unified-fitted and saved | folder loaded, all files selected, Data + Unified Fit checked, update check off; window 1000×720; callouts 1–12 | whole window, no update banner | yes | done |
| data-selector-02-graph | data-selector | as 01 | Create Graph; results window 800×600 | data, fits and residuals of three files | yes | done |
| data-selector-03-tabulate | data-selector | as 01 | Tabulate Results; window 900×420 | one row per file, Save as CSV… | yes | done |
| data-selector-04-export-menu | data-selector | any graph (e.g. 02) | right-click on the plot | the five export actions of the context menu | no | todo |

## feature-identifier

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
| feature-identifier-01-window | feature-identifier | validationData/unified/unified_two_level.h5 | loaded in Unified Fit, Identify Features…, Show advanced params, Detect segments; dialog 440×760; callouts 1–5 | parameters, buttons, result summary | yes | done |
| feature-identifier-02-markers | feature-identifier | validationData/unified/unified_two_level.h5 | as 01 without advanced params; grab the Unified Fit graph | GP/PLS regions and GK bands with labels | yes | done |

## diffraction-lines

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
