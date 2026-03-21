# Scattering Contrast Calculator

The **Scattering Contrast Calculator** computes X-ray and neutron scattering length
densities (SLDs) and the contrast between two user-defined compounds (a phase and a
matrix/solvent).  It is a **support / experiment-planning tool** — it does not analyse
measured data, but helps you predict whether a measurement will produce detectable
signal and at which X-ray energy anomalous contrast matching or enhancement is useful.

This tool is ported from the Igor Pro `IR1K_ScattContrast.ipf` module of the Irena package.

---

## Contents

1. [Launching the tool](#launching-the-tool)
2. [Panel layout](#panel-layout)
3. [Compound definition](#compound-definition)
4. [Isotope selection (neutrons)](#isotope-selection-neutrons)
5. [Calculation parameters](#calculation-parameters)
6. [Calculate button — results table](#calculate-button--results-table)
7. [Energy scan and graph window](#energy-scan-and-graph-window)
8. [Compound library](#compound-library)
9. [File-based export / import of compounds](#file-based-export--import-of-compounds)
10. [Exporting results](#exporting-results)
11. [Physics reference](#physics-reference)

---

## Launching the tool

| Method | Command |
|--------|---------|
| From the Data Selector | Click **Scattering Contrast (GUI)** in the *Support Tools* section |
| CLI entry point | `pyirena-contrast` |

No loaded dataset is required — the tool is always available.

---

## Panel layout

```
┌─────────────────────────────────────────────┬───────────────────────────────┐
│  Compound 1 (Phase)            Compound 2 (Matrix / Solvent)  [? Help]      │
│  Name / Formula / Mode / Density            same                            │
│  Isotope Selection table                    same                            │
│  Saved Compounds Library                    same                            │
├─────────────────────────────────────────────┤                               │
│  Calculation Parameters                     │  Results table                │
│  Energy Scan Range                          │  (read-only, right-click copy)│
│  [Calculate]  [Energy Scan]                 │                               │
│                                             │  [Export CSV]  [Save HDF5]    │
└─────────────────────────────────────────────┴───────────────────────────────┘
```

The left scroll panel defines both compounds and controls all calculations.  The right
panel shows the results table.  A separate floating **Graph Window** (opened by the
Energy Scan button or the *Show Graphs* button) shows energy-dependent plots.

---

## Compound definition

Each compound is defined by four fields:

| Field | Description |
|-------|-------------|
| **Name** | Human-readable label (appears in the results table header) |
| **Formula** | Chemical formula or composition string (see modes below) |
| **Mode** | How the formula is interpreted |
| **Density** | Bulk density in g/cm³ |

Click **Set Vacuum** to set the compound to empty space (zero SLD) — useful when you
only care about the SLD of a single material.

### Composition modes

#### Atomic formula  *(e.g. H2O, Fe3O4, SiO2, Ca5(PO4)3OH)*

Enter a standard chemical formula.  Element symbols are capitalised; numbers give the
count of each atom per formula unit.  Parentheses are supported.  The formula is parsed
by the `periodictable` library.

Examples:
```
H2O          → water
SiO2         → silica
Fe3O4        → magnetite
Ca5(PO4)3OH  → hydroxyapatite
```

#### Wt-fractions of elements  *(e.g. Au0.35Ag0.65)*

Numbers immediately after element symbols are interpreted as **mass fractions** (not
counts).  All fractions must sum to 1.  The tool converts them to molar ratios
internally.

Examples:
```
Au0.35Ag0.65          → 35 wt% gold, 65 wt% silver alloy
Pb0.5Sn0.3Cu0.2       → Pb/Sn/Cu alloy by weight
```

#### Wt-fractions of compounds  *(e.g. Y2O3:0.10 ZrO2:0.90)*

Tokens are `CompoundFormula:fraction` separated by spaces.  Each token contributes
proportionally to the composite.  Useful for oxide mixtures, glasses, and ceramics.

Examples:
```
Y2O3:0.10 ZrO2:0.90     → 10 wt% yttria-stabilised zirconia
SiO2:0.72 Al2O3:0.28    → silica-alumina mix
```

### Parsing the formula

Click **Parse ▸** (next to the Formula field) to populate the **Isotope Selection**
table from the current formula without running a full calculation.  This is useful to
set isotope choices before clicking Calculate.

---

## Isotope selection (neutrons)

The **Neutron Isotope Selection** table (collapsible) lists every element found in the
formula.  Each element has a drop-down with:

- **natural (default)** — uses the natural-abundance coherent scattering length b_c
- Named isotopes — e.g. `2  [b_c = −3.739 fm]` for deuterium

Selecting an isotope affects only neutron calculations (SLD, contrast, neutron b).
X-ray calculations always use the natural atomic number Z.

**Common use case:** switch H → ²H (deuterium) to model D₂O solvent or deuterated
polymers.

---

## Calculation parameters

| Parameter | Description |
|-----------|-------------|
| **X-ray energy** | Energy for anomalous X-ray calculations (keV).  Does not affect free-electron SLD. |
| **Sample thickness** | Thickness in mm used to compute transmission T = exp(−μ·d). |
| **Vol. frac. compound 1** | Volume fraction of Compound 1 in the mixed sample.  Used to compute the *combined* sample transmission. |

---

## Calculate button — results table

Click **Calculate** to compute SLDs, contrast, and anomalous corrections at the
specified X-ray energy.  The right-hand results table fills with:

### Molecular properties

| Row | Units | Notes |
|-----|-------|-------|
| Molecular weight | g/mol | Sum of atomic masses for one formula unit |
| Weight / formula unit | g | M / N_A |
| Formula units / cm³ | cm⁻³ | ρ N_A / M |
| Electrons / formula unit | — | Σ n_i Z_i |
| Electrons / cm³ | cm⁻³ | electrons per formula unit × formula units/cm³ |
| Volume / formula unit | cm³ | M / (ρ N_A) |

### X-ray SLD (free electron)

Computed from the number of electrons per unit volume using the classical electron
radius r_e = 2.818 fm.  Independent of energy.

| Row | Units |
|-----|-------|
| X-ray SLD | 10¹⁰ cm⁻² |
| X-ray SLD / gram | 10¹⁰ cm/g |

### Neutron SLD

Computed from coherent scattering lengths b_c (using your isotope selections).

| Row | Units |
|-----|-------|
| Total neutron b | cm (per formula unit) |
| Neutron SLD | 10¹⁰ cm⁻² |
| Neutron SLD / gram | 10¹⁰ cm/g |

### Contrast

| Row | Units | Notes |
|-----|-------|-------|
| X-ray contrast (Δρ)² | 10²⁰ cm⁻⁴ | Free-electron X-ray, energy-independent |
| Neutron contrast (Δρ)² | 10²⁰ cm⁻⁴ | |
| X-ray / Neutron ratio | — | Useful for SAXS vs SANS comparison |

### Anomalous X-ray (Chantler tables)

These rows use `xraydb` Chantler tables at the specified energy.

| Row | Units | Notes |
|-----|-------|-------|
| X-ray SLD (anomalous) | 10¹⁰ cm⁻² | One value per compound |
| Linear absorption μ | cm⁻¹ | One value per compound |
| Transmission (compound) | — | T = exp(−μ · d) for the pure compound |
| Sample transmission | — | Combined T using the volume fraction setting |
| X-ray contrast (anomalous) (Δρ)² | 10²⁰ cm⁻⁴ | Using anomalous SLDs |

### Copying values from the results table

- Click any cell to select it.
- Right-click → **Copy cell** (or **Ctrl+C**) to copy to clipboard.
- Multi-cell selection with Ctrl+click or Shift+click; Ctrl+C copies all selected
  cells as tab-separated text.

---

## Energy scan and graph window

Set the **Energy Scan Range** (start, end, number of points) and click **Energy Scan**
to compute anomalous contrast and transmission across the full range.

A floating **Graph Window** opens automatically showing three plots:

| Plot | Content |
|------|---------|
| X-ray Contrast (Δρ)² | Anomalous contrast vs energy [keV] |
| Linear absorption μ | μ [cm⁻¹] for Compound 1 (red) and Compound 2 (blue) |
| Transmission | T for Compound 1, Compound 2, and Sample |

**Crosshair cursor** — move the mouse over any plot to see the energy and Y value in
the readout at the bottom of the window.  All three plots share the same X axis; the
vertical crosshair line moves across all plots simultaneously.

**Log Y axis** — the top checkbox toggles logarithmic Y scale on the contrast plot.

**JPEG export** — right-click on any plot → *Save graph as JPEG…*

Click **Show Graphs** at any time to bring the graph window back to the front.

---

## Compound library

Each compound editor has a **Saved Compounds Library** section backed by a single HDF5
file at `~/.pyirena/contrast_compounds.h5`.

| Button | Action |
|--------|--------|
| **↻** | Refresh the drop-down from the library file |
| **Load** | Load the selected compound into this slot |
| **Save** | Save the current compound definition to the library (prompts for name) |
| **Del** | Delete the selected compound from the library |

The library file is shared between both compound slots and persists between sessions.

---

## File-based export / import of compounds

These buttons let you share compound definitions between computers or with colleagues.

### Export to File…

Saves the current compound definition to a portable HDF5 file you choose.

- If the file **does not exist** it is created.
- If the file **already contains compounds**, a dialog asks whether to
  **Append** the new compound to the file or **Replace** the entire file.
- The status bar reports how many compounds are in the file after saving.

You can build up a shared library by exporting multiple compounds to the same file.

### Import from File…

Loads a compound from a shared HDF5 file into this slot.

- If the file contains **one compound**, it is loaded immediately (no dialog).
- If the file contains **multiple compounds**, a dialog lets you either:
  - **Load Selected into Slot** — loads the highlighted compound into this compound slot.
  - **Add All to Local Library** — copies every compound from the file into the local
    library (`~/.pyirena/contrast_compounds.h5`) without loading into a slot.
- After loading into a slot, you are offered the option to also add the compound to
  the local library.

---

## Exporting results

| Button | Output |
|--------|--------|
| **Export Results CSV** | Compound properties and contrast values in a readable CSV table |
| **Export Scan CSV** | Energy scan arrays (E, SLD, μ, T) as a CSV spreadsheet |
| **Save Scan HDF5** | Energy scan arrays plus compound definitions in an HDF5 file |

---

## Physics reference

### X-ray SLD — free electron approximation

```
r_e = 2.8179 × 10⁻¹³ cm        (classical electron radius)
N_A = 6.0221 × 10²³ mol⁻¹

N_e = Σ nᵢ Zᵢ                   electrons per formula unit
n   = ρ N_A / M                  formula units / cm³

ρ_xray = N_e × n × r_e          [cm⁻²]  → divide by 10¹⁰ → [10¹⁰ cm⁻²]
```

Energy-independent; uses atomic number Z (not anomalous factors).

### X-ray SLD — anomalous (Chantler)

```
f1ᵢ(E) = xraydb.f1_chantler(element, E_eV)   (includes Z)
ρ_anom = n × r_e × Σ nᵢ f1ᵢ(E)
```

### Neutron SLD

```
bᵢ = periodictable[element].neutron.b_c       [fm = 10⁻¹³ cm]
     (isotope: e.g. periodictable.H[2].neutron.b_c)
B_total = Σ nᵢ bᵢ × 10⁻¹³                    [cm/formula unit]
V = M / (ρ N_A)                               [cm³/formula unit]

ρ_neutron = B_total / V                       [cm⁻²]
```

### Contrast

```
(Δρ)²_xray   = (ρ₁_xray  − ρ₂_xray )²       [10²⁰ cm⁻⁴]
(Δρ)²_neut   = (ρ₁_neut  − ρ₂_neut )²       [10²⁰ cm⁻⁴]
(Δρ)²_anom   = (ρ₁_anom  − ρ₂_anom )²       [10²⁰ cm⁻⁴]
```

### Linear absorption and transmission

```
μ_mass,i = xraydb.mu_chantler(element, E_eV)   [cm²/g]
w_i      = nᵢ Mᵢ / M                           weight fraction
μ_linear = ρ × Σ w_i μ_mass,i                  [cm⁻¹]
T        = exp(−μ_linear × d_cm)               [dimensionless]
```

`d_cm` = sample thickness in cm (convert from the *Sample thickness* field in mm).
