# Manuals to do

Everything in the manuals that needs Jan's eyes, in one list for the walk-through.
Generated 16-11-2026 by `python scripts/make_manuals_todo.py` from the
`% CHECK(Jan): …` comments in `docs/` (invisible on the site) and the `automated: no` rows
of `planning/docs-screenshots/SHOTLIST.md`. Resolve an item by fixing the page and deleting
its comment (or taking the screenshot), then re-run the script.

**62 text checks, 13 screenshots.**

## Screenshots to take

Save as `docs/tools/img/<page>/<id>.png`, replace the `% SCREENSHOT <id>` comment in the
page with a `figure` (alt text and caption), and set the row in `SHOTLIST.md` to `done`.

| ID | dataset | state to set up | must show |
|---|---|---|---|
| unified-fit-06-real-data | testData/ (Jan picks a sample) | good multi-level fit, Results to graphs | a real multi-level curve and its fit |
| size-distribution-04-real-data | testData/ (Jan picks a sample) | regularization or MaxEnt fit with power-law and flat background, Results to graph | a real distribution and its fit |
| modeling-04-real-data | testData/ (Jan picks a sample) | multi-population fit on real data, Show individual? on | a real multi-population fit |
| simple-fits-04-real-data | testData/ (Jan picks a sample) | Guinier or Porod fit on real data | linearization with the fitted range highlighted |
| fractals-04-3d-viewer | any grown aggregate (Z 250, seed 12345) | real display, PyVista installed | 3D isosurface in the right pane |
| fractals-05-targets | real mass-fractal data with a good linked two-level Unified fit (Jan picks) | Open NeXus…, Grow One | target summary + colour-coded parameters |
| waxs-peak-fit-05-phase-id | testData/ real powder pattern + matching CIFs (Jan picks) | two phases imported, sticks aligned with peaks | a real phase identification |
| carbon-fit-04-real-carbon | real hard carbon USAXS+SAXS+WAXS (Jan picks) | full-range fit, Show components, WAXS zoom | a real carbon fit |
| saxs-morph-04-3d-view | validationData/simple_fits/simple_debye_bueche.h5 | as 01 in the real GUI with PyVista; rotate the isosurface | 3D bicontinuous isosurface with bounding box |
| data-merge-03-porod-real-data | testData/ (Jan picks a real USAXS + SAXS pair) | Porod display mode, before and/or after a good merge | mismatch as a step/slope in the flat I·Q⁴ line, then one flat line |
| data-manipulation-04-similarity | testData/ (Jan picks a real exposure series with radiation damage) | Average tab, frames selected, Check Similarity | coloured frames and the accepted/rejected results table |
| scattering-contrast-03-import-dialog | a compounds file with several entries (Jan prepares) | Import from File… on that file | the Load Selected into Slot / Add All to Local Library dialog |
| data-selector-04-export-menu | any graph (e.g. 02) | right-click on the plot | the five export actions of the context menu |

## Text to check, by page

### tools/unified-fit.md

- [ ] [line 132](../docs/tools/unified-fit.md#L132) — "Store local fits" was not found in the current panel; "Show selected level?" is listed where the old page had "Display local (Porod & Guinier) fits?".
- [ ] [line 228](../docs/tools/unified-fit.md#L228) — "Export ASCII" was not found in the current panel.

### tools/size-distribution.md

- [ ] [line 106](../docs/tools/size-distribution.md#L106) — "volume fraction = ∫P(r)dr / contrast" (Contrast row) — validationData/README.md normalises P_V(r) so that ∫P_V(r)dr is the volume fraction itself.
- [ ] [line 246](../docs/tools/size-distribution.md#L246) — the distribution plot's axis is labelled "P(r) (a.u.)", while this section says it shows D(r) = dV/dr in volume fraction per Å.

### tools/modeling.md

- [ ] [line 177](../docs/tools/modeling.md#L177) — no "Derived Results panel" below the graph was found in the current window; derived values appear in the "Results to graph" annotation and the reports.
- [ ] [line 242](../docs/tools/modeling.md#L242) — @pop_type lists three types; Guinier-Porod, Mass Fractal and Surface Fractal populations exist too.
- [ ] [line 439](../docs/tools/modeling.md#L439) — the reference titles above were added from the citations in the text (Hammouda 2010, Teixeira 1988, Ardell 1966); please confirm.
- [ ] [line 517](../docs/tools/modeling.md#L517) — "a core-shell global fit that takes ~60 s serially drops to ~17 s on 6 cores" (Fit method section) and this table (130 s serial, 21 s on 6 cores) quote different timings.
- [ ] [line 557](../docs/tools/modeling.md#L557) — the second table above was a continuation of the first one in the old page (the quote block split it); it is shown as its own table now.
- [ ] [line 737](../docs/tools/modeling.md#L737) — the window labels the Log-Normal parameters Min Size, Mean Size and Sdeviation (core keys min_size, mean_size, sdeviation), not Median, σ and Shift.

### tools/simple-fits.md

- [ ] [line 212](../docs/tools/simple-fits.md#L212) — the old page also showed `from pyirena import load_result; load_result("sample.h5", "simple_fits")` as a top-level convenience function. It raises ValueError: load_result supports only 'unified_fit' and 'size_distribution'. Removed from the runnable example; restore once load_result supports Simple Fits.
- [ ] [line 229](../docs/tools/simple-fits.md#L229) — the buttons are labelled "Save params to JSON" and "Load params from JSON" in the current window.
- [ ] [line 340](../docs/tools/simple-fits.md#L340) — the old Scripting API block used placeholder files (sample.h5, data/, pyirena_config.json). It now runs on validationData copies; the Sphere example adds the flat complex background the dataset needs, and the Invariant example uses sizes_sphere_flat_background.h5 (true φ = 0.006) with a flat-only background prefit. The loop over data/*.h5 and the fit_pyirena call, which need a config file, are described in the sentence below the code instead.
- [ ] [line 434](../docs/tools/simple-fits.md#L434) — in the Sphere shot the linearization panel is hidden rather than showing the "No linearization available" text.
- [ ] [line 452](../docs/tools/simple-fits.md#L452) — the window's controls are "Passes" (default 10) and "Calc. Uncertainty (MC)"; 50 is the default of fit_simple(n_mc_runs=...) in scripts.
- [ ] [line 490](../docs/tools/simple-fits.md#L490) — the old page had no References section; the two entries above come from the text. Add the papers for Teubner-Strey, Debye-Bueche, Benedetti-Ciccariello and Hermans.

### tools/fractals.md

- [ ] [line 24](../docs/tools/fractals.md#L24) — "GUI-only; no batch / scripting / headless API" — pyirena.core.fractals exposes grow_aggregate, intensity_unified and optimize_growth (used in "From Python" below), and pyirena_read_fractals reads saved aggregates.
- [ ] [line 88](../docs/tools/fractals.md#L88) — which dataset and Unified fit make a good fractal reference? A linked two-level fit of validationData/modeling/modeling_mass_fractal.h5 gives c = 0.67 < 1 (free B) or c = 1.00 with level 1's B at its lower bound (mass-fractal option), so the panel reports "not a valid fractal". Target detection itself works since 16-11-2026 (RgCutoff fix).
- [ ] [line 412](../docs/tools/fractals.md#L412) — the Monte-Carlo steps above describe a voxel-pair sampling method; the code (core/fractals.py intensity_montecarlo) now places uniform points in every primary sphere (optional 10 % polydispersity) and histograms ALL pair distances (Shape2SAS-style).

### tools/waxs-peak-fit.md

- [ ] [line 22](../docs/tools/waxs-peak-fit.md#L22) — the peak-shape list (Gaussian, Lorentzian, Pseudo-Voigt, log-normal) omits "Voigt", which the Shape combo and core PEAK_SHAPES offer.
- [ ] [line 97](../docs/tools/waxs-peak-fit.md#L97) — steps 10 and 12 name "Revert" and "Export Parameters"; the buttons are now "Revert back" and "Save params to JSON".
- [ ] [line 118](../docs/tools/waxs-peak-fit.md#L118) — the cursors are drawn blue (Q_min) and red (Q_max) dashed, and the Q range is shown as a "Fit Q range: 1.2 – 4 Å⁻¹" label, not boxes.
- [ ] [line 191](../docs/tools/waxs-peak-fit.md#L191) — there is no "Add Peak Manually" button in the current panel (code comment: add and remove are available via right-click on the graph). The Shape combo also offers "Voigt".
- [ ] [line 348](../docs/tools/waxs-peak-fit.md#L348) — "Save State", "Export Parameters" and "Import Parameters" are not in the current panel; the JSON buttons are "Save params to JSON" / "Load params from JSON".
- [ ] [line 461](../docs/tools/waxs-peak-fit.md#L461) — the two batch examples above are shown as text, not Python: they need a data/pyirena_config.json and *_waxs.h5 files that do not exist in the repository, so they cannot run as written (phase-5 rule). The second one also uses Path without importing it.

### tools/carbon-fit.md

- [ ] [line 61](../docs/tools/carbon-fit.md#L61) — "Data Browser" — the main window is called the Data Selector elsewhere in the docs.
- [ ] [line 481](../docs/tools/carbon-fit.md#L481) — the four examples below were ```python blocks; they reference carbon_powder.h5, data/*.h5 and /data/hard_carbon.h5, which do not exist, so they cannot run as written (phase-5 rule) and are shown as text.

### tools/saxs-morph.md

- [ ] [line 15](../docs/tools/saxs-morph.md#L15) — Igor equivalent marked uncertain (?) in irena-pyirena-map.md; the page's own Igor reference is IR3_3DTwoPhaseSolid.ipf, which suggests TwoPhaseSolid is right.
- [ ] [line 26](../docs/tools/saxs-morph.md#L26) — "refines parameters ... by fitting back to the data" conflicts with Method notes ("not an iterative fit"; the only fits are the two background pre-fits).
- [ ] [line 112](../docs/tools/saxs-morph.md#L112) — step 6 says "black data − background trace"; in the screenshot the data − background trace is green.
- [ ] [line 118](../docs/tools/saxs-morph.md#L118) — the button is now labelled "Save/Append to HDF5…".
- [ ] [line 172](../docs/tools/saxs-morph.md#L172) — the current panel has no Fit / limits rows for B, P or background (the code says they are not fittable in this workflow).
- [ ] [line 191](../docs/tools/saxs-morph.md#L191) — the current panel runs Calculate 3D on a background thread (_CalcWorker); the GUI stays responsive and only the button is disabled. Same in Limitations below.
- [ ] [line 378](../docs/tools/saxs-morph.md#L378) — "low-Q window where the data is dominated by the Porod tail" (also workflow step 2: "low-Q window ... dominated by the Q⁻⁴ Porod slope") — a Porod tail is normally a high-Q feature; the low-Q power law is usually a different (larger-scale) contribution.
- [ ] [line 396](../docs/tools/saxs-morph.md#L396) — the Cube side tooltip gives different numbers (e.g. 256³ ~15 s, 125 MB RAM; 512³ ~120 s, 1 GB). The 256³ Calculate 3D of the example took about 4 s.
- [ ] [line 419](../docs/tools/saxs-morph.md#L419) — "Berk, N. F. (1991). Phys. Rev. E 51, 4141" — Phys. Rev. E volume 51 is 1995 (Berk's 1991 paper is Phys. Rev. A 44, 5069); year or journal looks inconsistent.

### tools/data-merge.md

- [ ] [line 137](../docs/tools/data-merge.md#L137) — the panel's Scale menu starts at DS1, while the config default `scale_dataset` is 2 (DS2) and Concepts says DS1 is the absolute-scale reference.
- [ ] [line 218](../docs/tools/data-merge.md#L218) — the Save JSON Config… tooltip says it also saves the folders; the page says "Q range and optimisation settings".
- [ ] [line 363](../docs/tools/data-merge.md#L363) — the Split at left cursor tooltip says "Default (unchecked): include all data from both datasets", not DS1 only up to cursor B and DS2 only from cursor A.
- [ ] [line 395](../docs/tools/data-merge.md#L395) — with Scale = DS2 and no Q shift the engine reports "Analytical WLS" (closed form; the api docstring says the same), not Nelder-Mead.

### tools/data-manipulation.md

- [ ] [line 63](../docs/tools/data-manipulation.md#L63) — the plot's y range follows the input data, so the green result is cut off below I ≈ 0.25 cm⁻¹ at Q > 0.04 Å⁻¹.
- [ ] [line 204](../docs/tools/data-manipulation.md#L204) — "Lower threshold = stricter (fewer rejections); higher = looser (more rejections)" — with p < threshold rejected, a lower threshold rejects fewer frames, which is usually called looser; the threshold guide below ("0.001 (strict)") uses the same wording.
- [ ] [line 431](../docs/tools/data-manipulation.md#L431) — this says the guard for subtract/divide lives in the core engine; the AIDA/MCP text says the api layer adds the slit-smearing guard for averaging "the panel relies on the user for".

### tools/scattering-contrast.md

- [ ] [line 232](../docs/tools/scattering-contrast.md#L232) — the linear-absorption axis in the graph window is labelled "(kcm⁻¹)" (pyqtgraph SI prefix on cm⁻¹), while this table says μ [cm⁻¹].

### tools/data-explorer.md

- [ ] [line 33](../docs/tools/data-explorer.md#L33) — "Data Merge (merged SAXS/USAXS datasets)" — there is no Data Merge preset or Collect type in the current panel; merged files are plotted as NXcanSAS data.
- [ ] [line 194](../docs/tools/data-explorer.md#L194) — the group menu also offers "Plot Modeling model" (not in the table).
- [ ] [line 217](../docs/tools/data-explorer.md#L217) — the panel also has "Size Dist. surf. S(r)", "Size Dist. cumul. surf.", "Modeling model", "Modeling vol. P(r)" and "Modeling num. P(r)", not described here.
- [ ] [line 250](../docs/tools/data-explorer.md#L250) — the Type list also has "Modeling" and "Carbon model".
- [ ] [line 295](../docs/tools/data-explorer.md#L295) — the buttons are labelled "+ Add", "Remove sel." and "Clear all", and the collect button "Collect from all selected files" (not "Multi-Collect from all selected files").
- [ ] [line 308](../docs/tools/data-explorer.md#L308) — a fourth tab, "Export to Igor", exists in the panel and is not described on this page (Igor export is described on the Data Selector page).
- [ ] [line 410](../docs/tools/data-explorer.md#L410) — the Collect window toolbar also has a "Save ITX" button.
- [ ] [line 515](../docs/tools/data-explorer.md#L515) — this table omits `.h5xp`, which the file list does scan (see "Selecting a folder").

### tools/data-selector.md

- [ ] [line 285](../docs/tools/data-selector.md#L285) — text files are already supported (File Type "Text Files", converted to NXcanSAS on first use, see guides/data-import-and-cleaning); this "Future" section is stale.
- [ ] [line 424](../docs/tools/data-selector.md#L424) — the Data Selector graphs are pyqtgraph, not matplotlib; "matplotlib Qt backend" advice may be stale.
- [ ] [line 450](../docs/tools/data-selector.md#L450) — stale paths — the panel is now pyirena/gui/data_selector/panel.py, Qt detection is pyirena/gui/_qt.py, and pyirena/gui/data_selector.py no longer exists (so the Testing command fails). developer/adding-features.md is the current checklist.
- [ ] [line 464](../docs/tools/data-selector.md#L464) — most of these exist now (text import, Unified Fit panel, interactive fitting, graph export, data processing tools, batch processing, setup save/restore).

### tools/feature-identifier.md

- [ ] [line 74](../docs/tools/feature-identifier.md#L74) — on this two-level synthetic curve the detector reports 5 segments and "Suggested Unified Fit levels: 5", and labels the flat high-Q background as "PLS P=0.6" rather than background (visible in both shots). Step 4's "you may need fewer" applies; worth a look.
- [ ] [line 77](../docs/tools/feature-identifier.md#L77) — step 3 says the summary is on the "right side"; it is below the buttons.
- [ ] [line 126](../docs/tools/feature-identifier.md#L126) — "Min segment width" default is 0.2 in the code (core FeatureDetectConfig and the dialog), not 0.10.
- [ ] [line 171](../docs/tools/feature-identifier.md#L171) — the knee label reads "GK ΔP=X.X" in the current GUI.
- [ ] [line 376](../docs/tools/feature-identifier.md#L376) — the knee keys are now `P_low_q`, `P_high_q`, `delta_P` (core and api).
- [ ] [line 382](../docs/tools/feature-identifier.md#L382) — there is no `pyirena_ctrl_detect_features` MCP tool any more; it is reached as `pyirena_call("detect_features", …)` through the dispatcher.
