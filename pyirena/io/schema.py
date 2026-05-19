"""
pyirena/io/schema.py — Machine-readable HDF5 result-group schema registry.

This module is the single authoritative Python source of truth for:
  - which HDF5 group each tool writes
  - what arrays are available to plot
  - what scalar quantities are extractable across a series of files
  - whether Data Merge / Data Manipulation strips the group

It complements (does not replace) docs/HDF5_Structure_Reference.md.
That document is the human-readable reference; this module is what
code should import to avoid hard-coding group paths, dataset names, or
unit strings in multiple places.

No heavy dependencies — only builtins.  Import freely.

Usage examples
--------------
Check whether a tool's group is present in a file::

    import h5py
    from pyirena.io.schema import TOOL_REGISTRY, tool_group_present

    with h5py.File("data.h5", "r") as f:
        if tool_group_present(f, "unified_fit"):
            ...

Enumerate all available tools in a file::

    from pyirena.io.schema import available_tools

    with h5py.File("data.h5", "r") as f:
        for key, schema in available_tools(f):
            print(key, schema["label"])

Get plottable curves for a tool::

    from pyirena.io.schema import TOOL_REGISTRY

    for plot in TOOL_REGISTRY["size_distribution"]["plots"]:
        print(plot["label"], plot["x"], "vs", plot["y"])

Get extractable scalars (for building a cross-file parameter series)::

    for sc in TOOL_REGISTRY["unified_fit"]["scalars"]:
        print(sc["key"], sc["path"], sc["units"])
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import h5py

# ---------------------------------------------------------------------------
# Schema version — bump when the registry structure itself changes (rare)
# ---------------------------------------------------------------------------
SCHEMA_VERSION = "1.0"


# ---------------------------------------------------------------------------
# Internal helpers for constructing per-subgroup scalar entries
# ---------------------------------------------------------------------------

def _per_subgroup(key: str, rel_path: str, units: str, label: str,
                  subgroup_prefix: str) -> dict:
    """Build a ScalarSpec entry for a quantity that repeats per sub-group."""
    return {
        "key":              key,
        "path":             rel_path,    # template: {n} = sub-group index (1-based)
        "units":            units,
        "label":            label,
        "per_subgroup":     True,
        "subgroup_prefix":  subgroup_prefix,
    }


def _scalar(key: str, path: str, units: str, label: str) -> dict:
    """Build a ScalarSpec entry for a top-level (once-per-file) scalar."""
    return {
        "key":              key,
        "path":             path,
        "units":            units,
        "label":            label,
        "per_subgroup":     False,
        "subgroup_prefix":  "",
    }


def _plot(x: str, y: str, x_units: str, y_units: str,
          x_label: str, y_label: str, label: str,
          plot_type: str = "custom") -> dict:
    """Build a PlotSpec entry."""
    return {
        "x":        x,
        "y":        y,
        "x_units":  x_units,
        "y_units":  y_units,
        "x_label":  x_label,
        "y_label":  y_label,
        "label":    label,
        "plot_type": plot_type,   # 'iq' | 'distribution' | 'spectral' |
                                  # 'correlation' | 'residuals' | 'custom'
    }


# ---------------------------------------------------------------------------
# Sub-group descriptors
# ---------------------------------------------------------------------------

_LEVEL_SUBGROUP = {
    "prefix":        "level_",
    "subgroup_fmt":  "d",       # level_1, level_2, …
    "attr_key":      None,      # sub-groups are uniform, no type dispatch
    "label":         "Level",
}

_POP_SUBGROUP = {
    "prefix":        "pop_",
    "subgroup_fmt":  "02d",     # pop_01, pop_02, …
    "attr_key":      "pop_type",
    "label":         "Population",
}

_PEAK_SUBGROUP = {
    "prefix":        "peak_",
    "subgroup_fmt":  "02d",     # peak_01, peak_02, …
    "attr_key":      "shape",
    "label":         "Peak",
}

_AGGREGATE_SUBGROUP = {
    "prefix":        "aggregate_",
    "subgroup_fmt":  "d",       # aggregate_1, aggregate_2, …
    "attr_key":      None,
    "label":         "Aggregate",
}


# ---------------------------------------------------------------------------
# TOOL_REGISTRY
# ---------------------------------------------------------------------------
# Keys are short stable identifiers used throughout the codebase.
# Do NOT change them — they will appear in state files and configs.
#
# Each entry is a dict with:
#   group          str  — HDF5 path from file root (e.g. "entry/unified_fit_results")
#   analysis_type  str  — value of @analysis_type attribute on the HDF5 group
#                          (empty string if the tool does not write this attribute)
#   label          str  — human-readable tool name (for UI labels)
#   nx_class       str  — value of @NX_class on the group
#   strips_on_merge bool — True: Data Merge / Data Manipulation strip this group
#                          False: group survives merge/manipulation
#   plots          list[PlotSpec]   — available Y-vs-X curves
#   scalars        list[ScalarSpec] — scalar quantities extractable across files
#   sub_groups     dict | None      — repeating numbered sub-groups, or None

TOOL_REGISTRY: dict[str, dict] = {

    # ── Unified Fit (Beaucage) ────────────────────────────────────────────
    "unified_fit": {
        "group":          "entry/unified_fit_results",
        "analysis_type":  "Unified Fit",
        "label":          "Unified Fit",
        "nx_class":       "NXprocess",
        "strips_on_merge": True,
        "plots": [
            _plot("Q", "intensity_data",  "1/angstrom", "1/cm",
                  "Q (Å⁻¹)", "I (cm⁻¹)", "Data",    plot_type="iq"),
            _plot("Q", "intensity_model", "1/angstrom", "1/cm",
                  "Q (Å⁻¹)", "I (cm⁻¹)", "Model",   plot_type="iq"),
            _plot("Q", "residuals",       "1/angstrom", "",
                  "Q (Å⁻¹)", "Residuals", "Residuals", plot_type="residuals"),
        ],
        "scalars": [
            _scalar("background",  "background",  "1/cm", "Background"),
            _scalar("chi_squared", "chi_squared", "",     "χ²"),
            # Per-level quantities — path uses {n} as 1-based index placeholder
            _per_subgroup("G",        "level_{n}/G",        "",    "G (Guinier prefactor)",  "level_"),
            _per_subgroup("Rg",       "level_{n}/Rg",       "Å",   "Rg",                     "level_"),
            _per_subgroup("B",        "level_{n}/B",        "",    "B (Porod constant)",      "level_"),
            _per_subgroup("P",        "level_{n}/P",        "",    "P (Porod exponent)",      "level_"),
            _per_subgroup("ETA",      "level_{n}/ETA",      "Å",   "ETA (corr. hole size)",  "level_"),
            _per_subgroup("PACK",     "level_{n}/PACK",     "",    "PACK (packing factor)",   "level_"),
            _per_subgroup("RgCutoff", "level_{n}/RgCutoff", "Å",   "RgCutoff",               "level_"),
        ],
        "sub_groups": _LEVEL_SUBGROUP,
    },

    # ── Simple Fits ───────────────────────────────────────────────────────
    "simple_fits": {
        "group":          "entry/simple_fit_results",
        "analysis_type":  "",          # this tool does not write @analysis_type
        "label":          "Simple Fits",
        "nx_class":       "NXprocess",
        "strips_on_merge": True,
        "plots": [
            _plot("Q", "intensity_data",  "1/angstrom", "1/cm",
                  "Q (Å⁻¹)", "I (cm⁻¹)", "Data",    plot_type="iq"),
            _plot("Q", "I_model",         "1/angstrom", "arb",
                  "Q (Å⁻¹)", "I (arb)",  "Model",   plot_type="iq"),
            _plot("Q", "residuals",       "1/angstrom", "",
                  "Q (Å⁻¹)", "Residuals", "Residuals", plot_type="residuals"),
        ],
        "scalars": [
            _scalar("chi_squared",         "chi_squared",         "",    "χ²"),
            _scalar("reduced_chi_squared", "reduced_chi_squared", "",    "Reduced χ²"),
            # params/* contents are model-dependent; enumerated at runtime
            # from the file's params/ sub-group.  Common ones listed here:
            _scalar("param_I0",   "params/I0",   "",   "I₀"),
            _scalar("param_Rg",   "params/Rg",   "Å",  "Rg"),
            _scalar("param_B",    "params/B",    "",   "B (Porod constant)"),
            _scalar("param_P",    "params/P",    "",   "P (Porod exponent)"),
            _scalar("param_BG_G", "params/BG_G", "1/cm", "Background (Guinier)"),
            _scalar("param_BG_P", "params/BG_P", "1/cm", "Background (Porod)"),
        ],
        "sub_groups": None,
        # Note: model name is in group attribute @model; params/* keys vary per model.
        # Consumers should enumerate params/* from the file when reading scalars.
    },

    # ── Size Distribution ─────────────────────────────────────────────────
    "size_distribution": {
        "group":          "entry/sizes_results",
        "analysis_type":  "",          # not written by this tool
        "label":          "Size Distribution",
        "nx_class":       "NXprocess",
        "strips_on_merge": True,
        "plots": [
            _plot("Q",      "intensity_data",  "1/angstrom", "1/cm",
                  "Q (Å⁻¹)", "I (cm⁻¹)", "Data",           plot_type="iq"),
            _plot("Q",      "intensity_model", "1/angstrom", "1/cm",
                  "Q (Å⁻¹)", "I (cm⁻¹)", "Model",          plot_type="iq"),
            _plot("Q",      "residuals",       "1/angstrom", "",
                  "Q (Å⁻¹)", "Residuals", "Residuals",      plot_type="residuals"),
            _plot("r_grid", "distribution",    "angstrom", "volume_fraction/angstrom",
                  "r (Å)", "P(r) (vol. frac./Å)", "Volume dist.", plot_type="distribution"),
            _plot("r_grid", "number_dist",     "angstrom", "1/angstrom",
                  "r (Å)", "N(r) (Å⁻¹)", "Number dist.",    plot_type="distribution"),
            _plot("r_grid", "cumul_vol_dist",  "angstrom", "volume_fraction",
                  "r (Å)", "Cumul. vol. frac.", "Cumul. volume", plot_type="distribution"),
        ],
        "scalars": [
            _scalar("chi_squared",      "chi_squared",      "",    "χ²"),
            _scalar("volume_fraction",  "volume_fraction",  "",    "Volume fraction"),
            _scalar("rg",               "rg",               "Å",   "Rg"),
            _scalar("q_power",          "q_power",          "",    "Q power"),
        ],
        "sub_groups": None,
    },

    # ── WAXS Peak Fit ─────────────────────────────────────────────────────
    "waxs_peakfit": {
        "group":          "entry/waxs_peakfit_results",
        "analysis_type":  "",          # not written by this tool
        "label":          "WAXS Peak Fit",
        "nx_class":       "NXprocess",
        "strips_on_merge": True,
        "plots": [
            _plot("Q", "intensity_data", "1/angstrom", "arb",
                  "Q (Å⁻¹)", "I (arb)", "Data",     plot_type="iq"),
            _plot("Q", "I_fit",          "1/angstrom", "arb",
                  "Q (Å⁻¹)", "I (arb)", "Total fit", plot_type="iq"),
            _plot("Q", "I_bg",           "1/angstrom", "arb",
                  "Q (Å⁻¹)", "I (arb)", "Background", plot_type="iq"),
            _plot("Q", "residuals",      "1/angstrom", "",
                  "Q (Å⁻¹)", "Residuals", "Residuals", plot_type="residuals"),
        ],
        "scalars": [
            _scalar("chi_squared",         "chi_squared",         "",       "χ²"),
            _scalar("reduced_chi_squared", "reduced_chi_squared", "",       "Reduced χ²"),
            # Per-peak quantities — path uses {n} (zero-padded 2-digit, 1-based)
            _per_subgroup("peak_Q0",   "peak_{n}/params/Q0",   "1/angstrom", "Peak Q₀ (Å⁻¹)",  "peak_"),
            _per_subgroup("peak_A",    "peak_{n}/params/A",    "arb",        "Peak amplitude",  "peak_"),
            _per_subgroup("peak_FWHM", "peak_{n}/params/FWHM", "1/angstrom", "Peak FWHM (Å⁻¹)", "peak_"),
            _per_subgroup("peak_eta",  "peak_{n}/params/eta",  "",           "Peak η (Voigt)",  "peak_"),
            # Derived: integral under the peak profile, ∫I_peak(q)dq.
            # NaN for HDF5 files written before this dataset was added.
            _per_subgroup("peak_Area", "peak_{n}/area",        "arb/angstrom", "Peak area (∫I dq)", "peak_"),
        ],
        "sub_groups": _PEAK_SUBGROUP,
        # Note: peak sub-group names are zero-padded, e.g. "peak_01", "peak_02".
        # {n} in path templates must be formatted with f"{n:02d}".
    },

    # ── Modeling (parametric forward) ─────────────────────────────────────
    "modeling": {
        "group":          "entry/modeling_results",
        "analysis_type":  "Modeling",
        "label":          "Modeling",
        "nx_class":       "NXprocess",
        "strips_on_merge": True,
        "plots": [
            # Modeling does not copy raw data arrays; only model output is stored.
            _plot("model_q", "model_I", "1/angstrom", "1/cm",
                  "Q (Å⁻¹)", "I (cm⁻¹)", "Total model", plot_type="iq"),
            # Per-population model_I is also available but requires sub-group traversal.
        ],
        "scalars": [
            _scalar("chi_squared",         "chi_squared",         "",    "χ²"),
            _scalar("reduced_chi_squared", "reduced_chi_squared", "",    "Reduced χ²"),
            _scalar("background",          "background",          "1/cm", "Background"),
            # Per-population scalars depend on @pop_type; listed for unified_level:
            _per_subgroup("pop_G",    "pop_{n}/G",    "",   "G (pop. Guinier prefactor)", "pop_"),
            _per_subgroup("pop_Rg",   "pop_{n}/Rg",   "Å",  "Rg (pop.)",                  "pop_"),
            _per_subgroup("pop_B",    "pop_{n}/B",    "",   "B (pop. Porod constant)",     "pop_"),
            _per_subgroup("pop_P",    "pop_{n}/P",    "",   "P (pop. Porod exponent)",     "pop_"),
        ],
        "sub_groups": _POP_SUBGROUP,
        # Note: pop sub-group names are zero-padded, e.g. "pop_01", "pop_02".
        # Scalar paths for size_dist and diffraction_peak populations differ —
        # consumers must dispatch on pop_N/@pop_type.
    },

    # ── SAXS Morph ────────────────────────────────────────────────────────
    "saxs_morph": {
        "group":          "entry/saxs_morph_results",
        "analysis_type":  "SAXS Morph",
        "label":          "SAXS Morph",
        "nx_class":       "NXprocess",
        "strips_on_merge": False,   # deliberately preserved across merge/manipulation
        "plots": [
            _plot("data_q",    "data_I",    "1/angstrom", "1/cm",
                  "Q (Å⁻¹)", "I (cm⁻¹)", "Data",                    plot_type="iq"),
            _plot("model_q",   "model_I",   "1/angstrom", "1/cm",
                  "Q (Å⁻¹)", "I (cm⁻¹)", "Model",                   plot_type="iq"),
            _plot("data_q",    "data_I_corr", "1/angstrom", "1/cm",
                  "Q (Å⁻¹)", "I (cm⁻¹)", "Data − background",      plot_type="iq"),
            _plot("r_grid",    "gamma_r",   "angstrom", "",
                  "r (Å)", "γ(r)", "Autocorrelation γ(r)",           plot_type="correlation"),
            _plot("spectral_k", "spectral_F", "1/angstrom", "angstrom**3",
                  "k (Å⁻¹)", "F(k) (Å³)", "Spectral density F(k)",  plot_type="spectral"),
        ],
        "scalars": [
            _scalar("chi_squared",         "chi_squared",         "",    "χ²"),
            _scalar("reduced_chi_squared", "reduced_chi_squared", "",    "Reduced χ²"),
            _scalar("volume_fraction",     "volume_fraction",     "",    "Volume fraction (φ)"),
            _scalar("contrast",            "contrast",            "1e20 cm^-4", "Contrast (Δρ)²"),
            _scalar("phi_actual",          "phi_actual",          "",    "Realised φ of voxelgram"),
            _scalar("rg_A",                "rg_A",                "Å",   "Rg"),
            _scalar("background",          "background",          "1/cm", "Background"),
            _scalar("power_law_B",         "power_law_B",         "",    "Power law B"),
            _scalar("power_law_P",         "power_law_P",         "",    "Power law P"),
            # morphology_metrics sub-group scalars (optional; present only when computed)
            _scalar("mm_n_clusters",             "morphology_metrics/n_clusters",             "",  "# clusters"),
            _scalar("mm_open_porosity",          "morphology_metrics/open_porosity_fraction", "",  "Open porosity"),
            _scalar("mm_closed_porosity",        "morphology_metrics/closed_porosity_fraction","", "Closed porosity"),
            _scalar("mm_euler_number",           "morphology_metrics/euler_number",           "",  "Euler number χ"),
            _scalar("mm_pore_size_median_A",     "morphology_metrics/pore_size_median_A",     "Å", "Median pore size"),
            _scalar("mm_pore_size_q25_A",        "morphology_metrics/pore_size_q25_A",        "Å", "Pore size Q25"),
            _scalar("mm_pore_size_q75_A",        "morphology_metrics/pore_size_q75_A",        "Å", "Pore size Q75"),
        ],
        "sub_groups": None,
    },

    # ── Fractals ──────────────────────────────────────────────────────────
    "fractals": {
        "group":          "entry/fractals_results",
        "analysis_type":  "Mass Fractal Aggregate",
        "label":          "Fractals",
        "nx_class":       "NXcollection",     # parent is NXcollection, not NXprocess
        "strips_on_merge": False,              # explicitly preserved across merge/manipulation
        "plots": [
            # Per-aggregate: intensity/Q vs intensity/I_unified or I_montecarlo
            _plot("intensity/Q", "intensity/I_unified",    "1/angstrom", "1/cm",
                  "Q (Å⁻¹)", "I (cm⁻¹)", "I(Q) analytical",  plot_type="iq"),
            _plot("intensity/Q", "intensity/I_montecarlo", "1/angstrom", "1/cm",
                  "Q (Å⁻¹)", "I (cm⁻¹)", "I(Q) Monte Carlo", plot_type="iq"),
        ],
        "scalars": [
            # Per-aggregate scalars — sub-group names are "aggregate_1", "aggregate_2", ...
            _per_subgroup("agg_Z",    "aggregate_{n}/parameters/Z",             "",  "Z (# primary spheres)", "aggregate_"),
            _per_subgroup("agg_df",   "aggregate_{n}/parameters/df",            "",  "df (mass-fractal dim.)", "aggregate_"),
            _per_subgroup("agg_dmin", "aggregate_{n}/parameters/dmin",          "",  "dmin",                   "aggregate_"),
            _per_subgroup("agg_c",    "aggregate_{n}/parameters/c",             "",  "c (connectivity dim.)",  "aggregate_"),
            _per_subgroup("agg_Rg_agg", "aggregate_{n}/parameters/RgAggregate", "Å", "Rg aggregate",           "aggregate_"),
            _per_subgroup("agg_Rg_prim","aggregate_{n}/parameters/RgPrimary",   "Å", "Rg primary sphere",      "aggregate_"),
        ],
        "sub_groups": _AGGREGATE_SUBGROUP,
    },

    # ── Data Merge (provenance only) ──────────────────────────────────────
    "data_merge": {
        "group":          "entry/data_merge_results",
        "analysis_type":  "",
        "label":          "Data Merge",
        "nx_class":       "NXprocess",
        "strips_on_merge": True,
        "plots":   [],    # no plottable arrays stored (merged data is in sasdata/)
        "scalars": [
            _scalar("scale",         "scale",         "",  "Merge scale factor"),
            _scalar("q_shift",       "q_shift",       "1/angstrom", "Q shift"),
            _scalar("background",    "background",    "1/cm", "Background"),
            _scalar("chi_squared",   "chi_squared",   "",  "χ² of merge"),
            _scalar("q_overlap_min", "q_overlap_min", "1/angstrom", "Q overlap min"),
            _scalar("q_overlap_max", "q_overlap_max", "1/angstrom", "Q overlap max"),
        ],
        "sub_groups": None,
    },

    # ── Data Manipulation (provenance only) ───────────────────────────────
    "data_manipulation": {
        "group":          "entry/data_manipulation_results",
        "analysis_type":  "Data Manipulation",
        "label":          "Data Manipulation",
        "nx_class":       "NXprocess",
        "strips_on_merge": True,
        "plots":   [],    # no plottable arrays; manipulated data is in sasdata/
        "scalars": [],    # operation-specific scalars; enumerated at runtime
        "sub_groups": None,
    },
}


# ---------------------------------------------------------------------------
# Convenience lookup tables (derived from TOOL_REGISTRY — do not edit)
# ---------------------------------------------------------------------------

#: Map HDF5 group path → tool key.
GROUP_TO_TOOL: dict[str, str] = {
    schema["group"]: key
    for key, schema in TOOL_REGISTRY.items()
}

#: Map @analysis_type attribute value → tool key (only for tools that write it).
ANALYSIS_TYPE_TO_TOOL: dict[str, str] = {
    schema["analysis_type"]: key
    for key, schema in TOOL_REGISTRY.items()
    if schema["analysis_type"]
}

#: Set of group paths that Data Merge / Data Manipulation strip.
STRIPPED_GROUPS: frozenset[str] = frozenset(
    schema["group"]
    for schema in TOOL_REGISTRY.values()
    if schema["strips_on_merge"]
)


# ---------------------------------------------------------------------------
# Query helpers
# ---------------------------------------------------------------------------

def tool_group_present(f: "h5py.File", tool_key: str) -> bool:
    """Return True if *tool_key*'s result group exists in the open HDF5 file *f*."""
    schema = TOOL_REGISTRY.get(tool_key)
    if schema is None:
        return False
    return schema["group"] in f


def available_tools(f: "h5py.File") -> list[tuple[str, dict]]:
    """Return list of (tool_key, schema) for every result group present in *f*.

    *f* must be an open h5py.File.  Results are in TOOL_REGISTRY insertion order.
    """
    return [
        (key, schema)
        for key, schema in TOOL_REGISTRY.items()
        if schema["group"] in f
    ]


def tool_key_for_group(group_path: str) -> str | None:
    """Return the tool key for an HDF5 group path, or None if unrecognised."""
    return GROUP_TO_TOOL.get(group_path)


def tool_key_for_analysis_type(analysis_type: str) -> str | None:
    """Return the tool key for a @analysis_type attribute value, or None."""
    return ANALYSIS_TYPE_TO_TOOL.get(analysis_type)
