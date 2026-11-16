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

## data-merge

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|

## data-manipulation

| ID | page | dataset | state to set up | must show | automated? | status |
|---|---|---|---|---|---|---|

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
