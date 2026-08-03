"""
pyirena/io/igor_names.py — Cross-reference between pyirena result names and
Igor Pro / Irena wave-type names (AllCurrentlyAllowedTypes).

This module answers: "given a pyirena tool and dataset, what Igor wave name
should be used when writing to an h5xp file?"

Structure
---------
TOOL_CROSS_REF maps each pyirena tool key to a dict with:
  igor_tool   : str  — Igor AllKnownToolsResults display name
  waves       : list[dict]  — plottable wave mappings (one dict per wave)

Each wave dict has:
  igor_name   : str  — Igor AllCurrentlyAllowedTypes name, with optional
                       "{n}" placeholder for 1-based sub-group index
  pyirena_x   : str  — HDF5 dataset path for X axis (relative to tool group)
  pyirena_y   : str  — HDF5 dataset path for Y axis (relative to tool group)
  per_subgroup: bool — True if {n} must be substituted (1-based)
  subgroup_prefix: str — HDF5 sub-group prefix (e.g. "level_", "pop_", "peak_")
  subgroup_fmt: str  — Python format spec for index, e.g. "d" for plain int,
                       "02d" for zero-padded 2-digit
  note        : str  — human-readable comment

IGOR_ONLY lists Igor wave types that have no pyirena equivalent (read-only
context: these can appear in h5xp files from Igor but pyirena cannot produce
them from its own results).

SIMPLE_FIT_MODEL_WAVE maps each pyirena SimpleFits model name (key in
MODEL_REGISTRY) to the Igor AllCurrentlyAllowedTypes wave name.  New pyirena
models that have no Igor counterpart use the "SimFit<Name>I" convention.

Notes
-----
- WAXS peak sub-groups: pyirena uses zero-padded 2-digit names ("peak_01"),
  but Igor peak numbers are 1-based integers (SADModelIntPeak1).
- Modeling pop sub-groups: pyirena uses zero-padded 2-digit names ("pop_01"),
  Igor population numbers are 1-based integers (ModelingVolDist_Pop1).
- Unified Fit per-level component waves (UniLocalLevel{n}Unified, Pwrlaw,
  Guinier) are Igor-only: pyirena stores only the total model intensity.
- saxs_morph has no Igor Irena equivalent; a descriptive name is used.
"""

from __future__ import annotations

import re

# ---------------------------------------------------------------------------
# Simple Fits model name → Igor AllCurrentlyAllowedTypes wave name
# ---------------------------------------------------------------------------
# Models present in Igor's list use the canonical Igor name.
# pyirena-only models follow the "SimFit<CamelCase>I" convention so that
# Igor can still load the waves (they just won't appear in a preset menu).

SIMPLE_FIT_MODEL_WAVE: dict[str, str] = {
    "Guinier":              "SimFitGuinierI",
    "Guinier Rod":          "SimFitGuinierRI",
    "Guinier Sheet":        "SimFitGuinierSII",
    "Porod":                "SimFitPorodI",
    "Power Law":            "SimFitPwrLawI",
    "Sphere":               "SimFitSphereI",
    "Spheroid":             "SimFitSpheroidI",
    # pyirena-only models (no Igor AllCurrentlyAllowedTypes entry):
    "Debye Polymer Chain":  "SimFitDebyePolymerI",
    "Debye-Bueche":         "SimFitDebyeBuecheI",
    "Teubner-Strey":        "SimFitTeubnerStreyI",
    "Benedetti-Ciccariello": "SimFitBenedettiI",
    "Hermans":              "SimFitHermansI",
    "Hybrid Hermans":       "SimFitHybridHermansI",
    "Unified Born Green":   "SimFitUnifiedBornGreenI",
    # Invariant is a calculation, not a fit; the "model" wave holds the
    # background curve.  The running integral is exported separately as
    # SimFitInvariantIntegral (vs SimFitInvariantIntegralQ).
    "Invariant":            "SimFitInvariantI",
}


# ---------------------------------------------------------------------------
# Main cross-reference table
# ---------------------------------------------------------------------------

def _wave(igor_name: str, x: str, y: str, note: str = "",
          per_subgroup: bool = False, prefix: str = "",
          fmt: str = "d") -> dict:
    return {
        "igor_name":       igor_name,
        "pyirena_x":       x,
        "pyirena_y":       y,
        "per_subgroup":    per_subgroup,
        "subgroup_prefix": prefix,
        "subgroup_fmt":    fmt,
        "note":            note,
    }


TOOL_CROSS_REF: dict[str, dict] = {

    # ── Unified Fit ───────────────────────────────────────────────────────
    "unified_fit": {
        "igor_tool": "Unified Fit",
        "waves": [
            _wave(
                "UnifiedFitIntensity",
                "Q", "intensity_model",
                note="Total Unified Fit model I(Q)",
            ),
            # Note: Igor also has UniLocalLevel{n}Unified / Pwrlaw / Guinier
            # (per-level components).  pyirena does not write these separately,
            # so they cannot be exported.  Marked as IGOR_ONLY below.
        ],
    },

    # ── Simple Fits ───────────────────────────────────────────────────────
    "simple_fits": {
        "igor_tool": "Simple Fits",
        "waves": [
            # Wave name is model-dependent; use SIMPLE_FIT_MODEL_WAVE at runtime.
            # The placeholder "<model>" is resolved by the writer.
            _wave(
                "<model>",
                "Q", "I_model",
                note="Simple Fit model I(Q); look up Igor name via SIMPLE_FIT_MODEL_WAVE",
            ),
        ],
    },

    # ── Size Distribution ─────────────────────────────────────────────────
    "size_distribution": {
        "igor_tool": "Size Distribution",
        "waves": [
            _wave("SizesFitIntensity",       "Q",      "intensity_model",
                  note="Size distribution model I(Q)"),
            _wave("SizesVolumeDistribution", "r_grid", "distribution",
                  note="Volume-weighted size distribution P(r)"),
            _wave("SizesNumberDistribution", "r_grid", "number_dist",
                  note="Number-weighted size distribution N(r)"),
            _wave("CumulativeSizeDist",      "r_grid", "cumul_vol_dist",
                  note="Cumulative volume distribution"),
        ],
    },

    # ── WAXS Peak Fit ─────────────────────────────────────────────────────
    "waxs_peakfit": {
        "igor_tool": "Small-angle diffraction",
        "waves": [
            _wave("SADModelIntensity", "Q", "I_fit",
                  note="Total WAXS peak fit I(Q)"),
            # Per-peak model: peak_{n:02d}/Q_peak + I_peak → SADModelIntPeak{n}
            _wave(
                "SADModelIntPeak{n}",
                "peak_{n:02d}/Q_peak", "peak_{n:02d}/I_peak",
                note="Per-peak model I(Q); n is 1-based",
                per_subgroup=True, prefix="peak_", fmt="02d",
            ),
        ],
    },

    # ── Modeling (parametric forward) ─────────────────────────────────────
    "modeling": {
        "igor_tool": "Modeling",
        "waves": [
            _wave("ModelingIntensity",
                  "model_q", "model_I",
                  note="Total modeling model I(Q)"),
            # Per-population (size_dist type only; n is 1-based):
            _wave(
                "ModelingVolDist_Pop{n}",
                "pop_{n:02d}/radius_grid", "pop_{n:02d}/volume_dist",
                note="Per-population volume distribution; size_dist populations only",
                per_subgroup=True, prefix="pop_", fmt="02d",
            ),
            _wave(
                "ModelingNumDist_Pop{n}",
                "pop_{n:02d}/radius_grid", "pop_{n:02d}/number_dist",
                note="Per-population number distribution; size_dist populations only",
                per_subgroup=True, prefix="pop_", fmt="02d",
            ),
            # Per-population model I(Q); Igor names it IntensityModelLSQF2pop{n}
            # but pyirena uses pop_{n:02d}/model_I:
            _wave(
                "IntensityModelLSQF2pop{n}",
                "model_q", "pop_{n:02d}/model_I",
                note="Per-population model I(Q) stored on model_q grid",
                per_subgroup=True, prefix="pop_", fmt="02d",
            ),
        ],
    },

    # ── SAXS Morph ────────────────────────────────────────────────────────
    "saxs_morph": {
        "igor_tool": None,   # no Igor Irena equivalent
        "waves": [
            _wave("SAXSMorphModelI",    "model_q",    "model_I",
                  note="SAXS Morph model I(Q); pyirena-only wave name"),
            _wave("SAXSMorphDataI",     "data_q",     "data_I",
                  note="SAXS Morph input data I(Q)"),
            _wave("SAXSMorphGammaR",    "r_grid",     "gamma_r",
                  note="Autocorrelation γ(r)"),
            _wave("SAXSMorphSpectralF", "spectral_k", "spectral_F",
                  note="Spectral density F(k)"),
        ],
    },

    # ── Fractals ──────────────────────────────────────────────────────────
    "fractals": {
        "igor_tool": "Fractals",
        "waves": [
            _wave("FractFitIntensity",
                  "intensity/Q", "intensity/I_unified",
                  note="Total fractal aggregate analytical I(Q)"),
            # Per-aggregate mass fractal component (n is 1-based):
            _wave(
                "Mass{n}FractFitInt",
                "aggregate_{n}/intensity/Q",
                "aggregate_{n}/intensity/mass_fractal",
                note="Per-aggregate mass fractal contribution",
                per_subgroup=True, prefix="aggregate_", fmt="d",
            ),
            _wave(
                "Surf{n}FractFitInt",
                "aggregate_{n}/intensity/Q",
                "aggregate_{n}/intensity/surface_fractal",
                note="Per-aggregate surface fractal contribution",
                per_subgroup=True, prefix="aggregate_", fmt="d",
            ),
        ],
    },

    # ── Data Merge / Data Manipulation ────────────────────────────────────
    # These are provenance tools; no model arrays → no wave export.
    "data_merge": {
        "igor_tool": None,
        "waves": [],
    },
    "data_manipulation": {
        "igor_tool": None,
        "waves": [],
    },
}


# ---------------------------------------------------------------------------
# Igor-only wave types (no pyirena equivalent)
# ---------------------------------------------------------------------------
# These appear in Igor experiments but cannot be generated from pyirena results.
# Provided here for documentation and for future round-trip reading.

IGOR_ONLY: dict[str, list[str]] = {
    "Unified Fit per-level components": [
        "UniLocalLevel{n}Unified",
        "UniLocalLevel{n}Pwrlaw",
        "UniLocalLevel{n}Guinier",
    ],
    "Analytical models": [
        "DebyeBuecheModelInt",   # Igor Analytical Models tool (not Simple Fits)
        "AnalyticalModelInt",
        "SysSpecModelInt",
    ],
    "PDDF": [
        "PDDFIntensity",
        "PDDFDistFunction",
        "PDDFChiSquared",
        "PDDFGammaFunction",
        "SADUnifiedIntensity",
    ],
    "Reflectivity": [
        "ReflModel",
        "SLDProfile",
    ],
    "Guinier-Porod": [
        "GuinierPorodFitIntensity",
        "GuinierPorodIntLevel{n}",  # n = 0-5
    ],
    "Evaluate Size Distribution": [
        "UnifSizeDistVolumeDist",
        "UnifSizeDistNumberDist",
        "CumulativeSfcArea",
        "MIPVolume",
    ],
    "1D Correlation (Simple Fits)": [
        "Corr1DK",
        "Corr1DGammaA",
        "Corr1DGammaI",
    ],
    "LSQF2 Modeling (legacy Igor)": [
        "IntensityModelLSQF2",
        "NumberDistModelLSQF2",
        "VolumeDistModelLSQF2",
        "IntensityModelLSQF2pop{n}",   # n = 0-6
        "NumberDistModelLSQF2pop{n}",
        "VolumeDistModelLSQF2pop{n}",
    ],
}


# ---------------------------------------------------------------------------
# Helper: resolve igor_name template at runtime
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Y-wave → X-wave lookup (ResultsDataTypesLookup from Igor/Irena)
# ---------------------------------------------------------------------------
# Static map: Y-wave name → canonical Igor X-wave name.
# Per-n entries (Pop{n}, Peak{n}, …) are handled by _RESULT_X_PATTERNS below.

RESULT_X_WAVE: dict[str, str] = {
    # Size Distribution
    "SizesFitIntensity":       "SizesFitQvector",
    "SizesVolumeDistribution": "SizesDistDiameter",
    "SizesNumberDistribution": "SizesDistDiameter",
    "CumulativeSizeDist":      "CumulativeDistDiameters",
    # Unified Fit
    "UnifiedFitIntensity":     "UnifiedFitQvector",
    # Modeling (total)
    "ModelingIntensity":       "ModelingQvector",
    # Fractals (total)
    "FractFitIntensity":       "FractFitQvector",
    # WAXS Peak Fit (total)
    "SADModelIntensity":       "SADModelQ",
    # Simple Fits — Igor AllCurrentlyAllowedTypes entries
    "SimFitGuinierI":          "SimFitGuinierQ",
    "SimFitGuinierRI":         "SimFitGuinierRQ",
    "SimFitGuinierSII":        "SimFitGuinierSQ",
    "SimFitSphereI":           "SimFitSphereQ",
    "SimFitSpheroidI":         "SimFitSpheroidQ",
    "SimFitPorodI":            "SimFitPorodQ",
    "SimFitPwrLawI":           "SimFitPwrLawQ",
    # Invariant running integral (pyirena-only; SimFitInvariantI is handled
    # by the SimFit<Name>I → SimFit<Name>Q regex below)
    "SimFitInvariantIntegral": "SimFitInvariantIntegralQ",
    # SAXS Morph (pyirena-only wave names)
    "SAXSMorphModelI":         "SAXSMorphModelQ",
    "SAXSMorphDataI":          "SAXSMorphDataQ",
    "SAXSMorphGammaR":         "SAXSMorphR",
    "SAXSMorphSpectralF":      "SAXSMorphK",
}

# Regex patterns for per-n Y-waves.  Each tuple is (compiled pattern, replacement).
# ``re.Match.expand(replacement)`` is used, so ``\1`` backreferences work.
_RESULT_X_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(r"^ModelingVolDist_Pop(\d+)$"),    r"ModelingDia_Pop\1"),
    (re.compile(r"^ModelingNumDist_Pop(\d+)$"),    r"ModelingDia_Pop\1"),
    (re.compile(r"^IntensityModelLSQF2pop(\d+)$"), r"ModelingQvector"),
    (re.compile(r"^SADModelIntPeak(\d+)$"),         r"SADModelQPeak\1"),
    (re.compile(r"^Mass(\d+)FractFitInt$"),         r"Mass\1FractFitQvec"),
    (re.compile(r"^Surf(\d+)FractFitInt$"),         r"Surf\1FractFitQvec"),
    # pyirena-only Simple Fits not in the static table: SimFit<Name>I → SimFit<Name>Q
    (re.compile(r"^(SimFit\w+?)I$"),                r"\1Q"),
]


def result_x_wave_name(y_name: str) -> str:
    """Return the canonical Igor X-wave name paired with *y_name*.

    Checks the static :data:`RESULT_X_WAVE` dict first, then applies regex
    patterns for per-n waves (e.g. ``ModelingVolDist_Pop3``).  Falls back to
    ``y_name + "_X"`` for unknown wave types.
    """
    if y_name in RESULT_X_WAVE:
        return RESULT_X_WAVE[y_name]
    for pattern, replacement in _RESULT_X_PATTERNS:
        m = pattern.match(y_name)
        if m:
            return m.expand(replacement)
    return y_name + "_X"


def resolve_igor_name(igor_name_template: str, n: int) -> str:
    """Substitute {n} in an Igor wave name template with the given 1-based index."""
    return igor_name_template.replace("{n}", str(n))


def resolve_pyirena_path(path_template: str, n: int, fmt: str = "d") -> str:
    """Substitute {n:fmt} in a pyirena HDF5 path template."""
    formatted_n = format(n, fmt)
    return path_template.replace(f"{{n:{fmt}}}", formatted_n).replace("{n}", str(n))
