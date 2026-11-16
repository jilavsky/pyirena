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

## simple-fits

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|

## fractals

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|

## waxs-peak-fit

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|

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

## data-explorer

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|

## data-selector

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|

## feature-identifier

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|

## diffraction-lines

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|
