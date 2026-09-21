# Decisions — settled 2026-09-21 (Jan)

The seven open questions from the planning pass are answered. Recorded here
so they are not re-litigated; the answers are already folded into
[01-models.md](01-models.md) and [02-implementation-plan.md](02-implementation-plan.md).

Implementation branch: `feature/carbon-model`.

---

**1. Tool name/key — `carbon_fit` internally, "Carbon model" to users.**

`carbon_fit` is the module/schema/HDF5/JSON key everywhere it is not forced
otherwise; the MCP dispatch category is the short noun `carbon`, matching
`unified`/`sizes`/`simple`/`modeling`/`waxs`. Every user-visible string —
panel title, Data Browser checkbox, menu entry, report heading — reads
**"Carbon model"**, never "carbon_fit".

**2. Units — Ångström, everywhere, no exceptions.**

Q in Å⁻¹, lengths in Å, contrast in 10²⁰ cm⁻², exactly as the rest of
pyIrena. No nm anywhere: not in the core, not in HDF5, not on the panel, no
unit toggle. The source paper's nm values are converted when comparing
against its tables, and only there. See 01-models.md §0 for the full unit
table. *Rationale: no scientific impact, and a second unit convention inside
one package is a permanent source of error.*

**3. SAXS region — one tool, one tab, a mode toggle.**

Single Carbon model tool with one SAXS-region tab carrying a model selector
(`fractal+globule` | `teubner_strey`). Selecting a mode shows only that
mode's parameters. The two are alternative descriptions of the same
micropore population and are never summed. More SAXS-region models are
expected later — the selector is built as a registry so adding one is a
dict entry, not a new branch through the panel.

**4. v2 (linked crumpled-layer geometry) — build it now.**

Not deferred. The WAXS section carries the crumpled-layer envelope
(Teixeira S(Q) × Beaucage discoid, eq. 16–17) as a per-section mode, and the
linkage to the SAXS region's fractal parameters is a real, fittable link
from day one rather than a hook. Cheaper now than a breaking core-model
change later.

**5. Porod-roughness formula (eq. 3) — corrected.**

The middle bracket term is `(1/5)(QR)²`, not `(1/3)(QR)²`. The planning pass
mis-read it off a low-resolution render. 01-models.md §2 now carries the
correct form and the question is closed.

**6. Reuse, not reimplementation.**

The value of this tool is the GUI, the WAXS→density→contrast→SAXS wiring,
and simultaneous fitting of the whole Q range — not the maths, which is
mostly already in pyIrena. Shared code wherever a shared function exists:
`simple_fits._teubner_strey`, `scattering_contrast.compute_compound`,
`waxs_peakfit`'s peak machinery (extended with a true Voigt, which the WAXS
Peak Fit tool gains too), `fit_metrics`, `reporting`. Where the paper's
parameterisation genuinely differs from an existing Irena one — the Teixeira
structure factor is the case in point, Irena's `_mass_fractal_intensity` uses
a different (η-based, no Γ) normalisation — the paper's form is written as a
small, separately tested pure function rather than bent onto the existing
signature. Correctness outranks reuse; reuse outranks convenience.

**7. Priority — GUI first, then scripting, then MCP; but all three are on the list.**

Build order: core → io → **GUI** (debugging and verification happen there) →
batch/scripting → api/control → MCP. Scripting and MCP are committed
deliverables, not optional extras, so the NXcanSAS result schema and the JSON
config shape are designed up front, in the core's `to_dict()`, before the
panel is written — that is what stops the tool being coded into a corner.

---

## Left open on purpose

- **Number of graphs.** The full USAXS+SAXS+WAXS range is ~5 decades in Q and
  a single log-log plot will crowd the diffraction peaks. Start with one
  full-range plot plus residuals; add a second linear-Q WAXS-zoom plot once
  there is real data on screen and the need is obvious. Usability decides.
- **Additional SAXS-region models** beyond fractal+globule and Teubner-Strey.
  The selector is a registry so this stays cheap.
- **Digitised regression fixtures from the paper's Tables 1–3** — a stretch
  goal for `validationData/`, not a blocker.
