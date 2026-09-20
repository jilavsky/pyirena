# Carbon Fitting — Planning Folder

Internal planning artifact for a new pyIrena analysis tool: full-range
SAXS+WAXS fitting of disordered carbonaceous materials (particle Porod
scattering + micropore SAXS models + turbostratic-stacking WAXS diffraction
peaks), following Saurel et al., *Energy Storage Materials* **21** (2019)
162–173 (`10.1016/j.ensm.2019.05.007`) and its corrigendum,
*Energy Storage Materials* **28** (2020) 418 (`10.1016/j.ensm.2020.03.013`).

**Status: planning only.** No code has been written. This folder is the
output of a planning session (2026-09-20) that read the source paper plus its
two supplementary-information PDFs, matched the required models against what
already exists in pyIrena, and produced a step-by-step wiring plan against
`docs/developer_adding_features.md`. Implementation happens later, on a
separate branch, per file at a time, starting from [02-implementation-plan.md](02-implementation-plan.md).

## Contents

| File | What's in it |
|---|---|
| [01-models.md](01-models.md) | The physics: every model, its parameters, derived quantities, and which source-paper equation it comes from. Read this before writing any core math. |
| [02-implementation-plan.md](02-implementation-plan.md) | The wiring: what's new vs. what already exists in pyIrena to reuse, proposed architecture, the full registry checklist from `docs/developer_adding_features.md` with proposed keys, phased build order, testing strategy. |
| [03-open-questions.md](03-open-questions.md) | Decisions for Jan to make before or during implementation — tool name, scope of v1 vs v2, a couple of equations that need re-verification against the original PDF at higher zoom, the nm↔Å unit-convention issue. |

## Source documents (not checked into the repo)

Three PDFs were read for this plan, supplied by Jan:

1. Saurel et al. 2019 main text — "A SAXS outlook on disordered carbonaceous
   materials for electrochemical energy storage" (`Energy Storage Materials
   21, 162–173`). Has the full working three-component model
   (`I = I_Porod + I_mp + I_waxs`), with per-sample refined parameters in
   Tables 1–3, and the 2020 corrigendum appended at the end (correct `w_P`,
   `w_C` formulas and a fixed sign in the density equation).
2. Supplementary Information, Annexes 1–3 (`mmc2`, 19 pp.) — Annex 1
   (intensity calibration, **not needed** — Jan does this at the instrument),
   Annex 2 (general SAXS theory, **not needed**, background only), Annex 3
   (**the model catalogue**: sphere/spheroid/globule form factors, disk/lamellar
   form factors, Teixeira fractal structure factor, Teubner-Strey, and the
   full lamellar-stacking diffraction-peak derivation with Voigt broadening —
   this is also where the WAXS peak model actually lives, not in `mmc3`).
3. Supplementary Information (`mmc3`, 4 pp.) — density/SLD method (`ρ_struc`
   formula, later corrected) and PXRD figures/tables. No additional model
   equations beyond what's in Annex 3.

Jan has an old Igor Pro implementation of these models but judged it not
worth digging out — the models are simple formulas; the main plan below
takes them straight from the paper.

## One-paragraph summary of the ask

Add a new pyIrena tool, GUI-similar to Modeling/Size Distribution, that fits
full-range SAXS→WAXS carbon powder data as the sum of three tabs: **Complex
background** (macroscopic-particle Porod scattering, optionally with a
surface-roughness term), **SAXS region** (micropore scattering — either a
dilute particle/pore form factor with an optional fractal-aggregate structure
factor, or the Teubner-Strey two-phase model), and **WAXS region** (one or
more turbostratic diffraction peaks, e.g. (002)/(100)/(004), with proper
Voigt broadening, orientation-averaging, and Debye-Waller-type intensity
correction). The tool's main value-add beyond the raw formulas is the
*materials-science convenience layer*: density → structural density → SLD →
contrast → micropore volume fraction → pore/wall widths, wired through so the
user enters a chemical formula and density once and gets ~15 physically
meaningful derived quantities, not just fit coefficients.
