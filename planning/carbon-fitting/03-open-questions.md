# Open questions — decide before or at the start of implementation

1. **Tool name/key.** Proposed `carbon_fit` (see 02-implementation-plan.md
   §0). Alternatives: `carbon`, `waxs_saxs_carbon`, something matching how
   Jan refers to it verbally with beamline users. Cheap to change now,
   expensive once it's baked into HDF5 group names / JSON config keys /
   Igor wave names (per AGENTS.md: "these are baked into saved files and
   config sections and cannot be renamed").

2. **Units shown to the user: nm⁻¹ or Å⁻¹?** The source paper and this
   sub-field (carbon SAXS/WAXS) work natively in nm. pyIrena's internal
   convention is Å throughout. Recommend: core math and HDF5 storage in Å
   (consistent with the rest of pyIrena, and with cross-tool parameter
   sharing e.g. `S_mp` compared against other tools' surface-area outputs),
   but worth asking whether the **panel's displayed** numbers (Q_c, d, r, ξ,
   R_rough, Σ) should be shown in nm for readability, with a small
   unit-toggle or just a fixed nm display — carbon/porous-material
   practitioners will expect nm-scale numbers, not "18 Å⁻¹". Precedent:
   check whether any other pyIrena panel already does an nm-for-display /
   Å-internally split before inventing a new convention.

3. **SAXS-region mode: fractal+globule vs. Teubner-Strey — one tool, toggle,
   or two?** The paper uses both, sample-dependent. Plan assumes a single
   dropdown/toggle within the SAXS tab (like Simple Fits' model combo). Confirm
   that's the right UX rather than, say, letting both be summed simultaneously
   (the paper never does this — they're alternative descriptions of the same
   micropore population, not additive contributions) — this affects the core
   model's shape (`mode: str` vs. a list of SAXS sub-components).

4. **Should WAXS peaks and SAXS-region crumpling share `R`/`Σ`/`D`
   (v2 feature) be plumbed as day-one architecture** (i.e., build the core
   model with the linkage hook even if the UI checkbox ships later), **or
   deferred entirely** until v2 is actually scheduled? Building the hook
   early avoids a breaking core-model change later, but adds complexity to
   v1. Recommend deciding based on how confident Jan is that v2 will
   actually happen soon vs. being aspirational.

5. **Re-verify the Porod-roughness formula (eq. 3) against a clean, zoomed
   read of the original PDF page.** The planning pass extracted it from a
   moderate-resolution page render and reconstructed the exact bracket power
   analytically (by requiring the paper's own stated Q→∞ limit, eq. 4, to
   hold) rather than trusting the OCR literally — see 01-models.md §2 for
   the reasoning. High confidence in the reconstructed form, but worth a
   30-second visual double-check on Jan's own copy of the PDF (page 4 of the
   main text, labeled page "165", equation (3)) before committing it to
   code and a test suite.

6. **How much of the "already exists" reuse should be literal function
   calls vs. a deliberate carbon-tool-local reimplementation for
   independence/stability?** E.g., calling into `modeling.py`'s
   `_mass_fractal_intensity` directly ties the new tool's correctness to
   Modeling's internals not changing incompatibly; a thin private copy in
   `carbon_fit.py` is more isolated but reintroduces the duplication
   AGENTS.md warns against. Precedent in the codebase (e.g. how
   `fit_metrics.py` or `core/reporting.py` are shared) suggests **shared
   module, not duplication** is the house style — recommend that — but
   flagging because it's a real design choice with real tradeoffs for a
   scientific-correctness-critical codebase per AGENTS.md's stated
   priorities.

7. **Agent-drivability priority.** Jan asked for "new scripting API as well
   as MCP wiring" explicitly — treating this as equal priority to the GUI
   from the start (registry #9/#10 in the implementation plan), rather than
   the more common pattern in this codebase of building control-surface
   support as a later pass once the GUI is stable. Confirm that's the
   intent (build API/control alongside the GUI, phase 5 of §7 in the
   implementation plan) rather than strictly after it.
