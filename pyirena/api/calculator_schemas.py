"""JSON schemas for pyirena.api.calculators in Anthropic tool-use format.

Mirrors pyirena/api/control/schemas.py, but kept as a separate registry:
calculators are stateless (no session_id) and are not part of the control
surface, so they must not inflate TOOL_SCHEMA_BY_NAME.

pyirena/mcp/dispatch.py merges this registry with the control one to expose
the "calculators" category over MCP.

Usage
-----
from pyirena.api.calculator_schemas import CALCULATOR_TOOL_SCHEMAS
# pass to client.messages.create(tools=CALCULATOR_TOOL_SCHEMAS, ...)
"""
from __future__ import annotations

from pyirena.api.calculators import COMPOSITION_MODES

_MODE_PROPERTY = {
    "type": "string",
    "enum": list(COMPOSITION_MODES),
    "description": (
        "How to read the formula. 'atomic_ratio' is standard chemical "
        "notation ('Ti2O3'). 'weight_fraction_elements' reads the numbers as "
        "weight fractions summing to 1 ('Au0.35Ag0.65'). "
        "'weight_fraction_compounds' takes space-separated 'formula:fraction' "
        "tokens ('Y2O3:0.10 ZrO2:0.90'). In the two weight-fraction modes the "
        "per-formula-unit quantities are basis-dependent; the SLDs and "
        "contrasts are still correct."
    ),
    "default": "atomic_ratio",
}

_ISOTOPES_PROPERTY_DESC = (
    "Optional isotope overrides as {element_symbol: mass_number_string}, "
    "e.g. {\"H\": \"2\"} for deuterium. Affects the molar mass and the "
    "neutron scattering length, not the electron count. Use lookup_element "
    "to see which isotopes have neutron data."
)


CALCULATOR_TOOL_SCHEMAS: list[dict] = [

    # -----------------------------------------------------------------------
    # Scattering contrast
    # -----------------------------------------------------------------------
    {
        "name": "calc_contrast",
        "description": (
            "Compute the X-ray and neutron scattering contrast between two "
            "compounds from their chemical formulas and mass densities. No "
            "dataset needed. The returned 'xray_contrast' is (delta-rho)^2 in "
            "10^20 cm^-4 — exactly the units and convention the 'contrast' "
            "parameter of a Sizes set_shape() or a Modeling population "
            "expects, so the result can be used directly in a fit. Leave one "
            "formula empty (or its density 0) to get the contrast against "
            "vacuum. Pass energy_keV to also get the anomalous "
            "(Chantler-corrected) contrast, absorption and transmission, "
            "which you need near an absorption edge where the free-electron "
            "approximation fails."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "formula_1": {
                    "type": "string",
                    "description": (
                        "Chemical formula of the first compound, e.g. 'TiO'. "
                        "Empty string means vacuum."
                    ),
                },
                "density_1": {
                    "type": "number",
                    "description": "Mass density of the first compound in g/cm^3.",
                },
                "formula_2": {
                    "type": "string",
                    "description": (
                        "Chemical formula of the second compound, e.g. "
                        "'Ti2O3'. Empty string means vacuum."
                    ),
                },
                "density_2": {
                    "type": "number",
                    "description": "Mass density of the second compound in g/cm^3.",
                },
                "name_1": {
                    "type": "string",
                    "description": "Optional display name for compound 1.",
                    "default": "",
                },
                "name_2": {
                    "type": "string",
                    "description": "Optional display name for compound 2.",
                    "default": "",
                },
                "mode": _MODE_PROPERTY,
                "isotopes_1": {
                    "type": ["object", "null"],
                    "description": _ISOTOPES_PROPERTY_DESC,
                },
                "isotopes_2": {
                    "type": ["object", "null"],
                    "description": _ISOTOPES_PROPERTY_DESC,
                },
                "energy_keV": {
                    "type": ["number", "null"],
                    "description": (
                        "X-ray energy in keV. Omit for the free-electron "
                        "approximation (fine away from edges); supply it for "
                        "anomalous contrast, absorption and transmission."
                    ),
                },
                "thickness_mm": {
                    "type": "number",
                    "description": (
                        "Sample thickness in mm, used only for the "
                        "transmission calculation."
                    ),
                    "default": 1.0,
                },
                "vol_frac_1": {
                    "type": "number",
                    "minimum": 0.0,
                    "maximum": 1.0,
                    "description": (
                        "Volume fraction of compound 1 in the sample, used "
                        "only to combine the two absorptions into "
                        "'transmission_sample'."
                    ),
                    "default": 0.01,
                },
            },
            "required": ["formula_1", "density_1", "formula_2", "density_2"],
        },
    },

    {
        "name": "calc_compound",
        "description": (
            "Compute the X-ray and neutron scattering length densities of a "
            "single compound from its chemical formula and mass density, plus "
            "molar mass, electron density and molar volume. Use this when you "
            "need one material's SLD rather than a contrast between two. "
            "SLDs are in 10^10 cm^-2. Pass energy_keV to add the anomalous "
            "SLD, the linear absorption coefficient and the transmission."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "formula": {
                    "type": "string",
                    "description": (
                        "Chemical formula, e.g. 'SiO2'. Empty string means "
                        "vacuum."
                    ),
                },
                "density": {
                    "type": "number",
                    "description": "Mass density in g/cm^3.",
                },
                "mode": _MODE_PROPERTY,
                "name": {
                    "type": "string",
                    "description": "Optional display name. Defaults to the formula.",
                    "default": "",
                },
                "isotopes": {
                    "type": ["object", "null"],
                    "description": _ISOTOPES_PROPERTY_DESC,
                },
                "energy_keV": {
                    "type": ["number", "null"],
                    "description": (
                        "X-ray energy in keV. Omit for the free-electron "
                        "approximation; supply it for the anomalous SLD, "
                        "absorption and transmission."
                    ),
                },
                "thickness_mm": {
                    "type": "number",
                    "description": "Sample thickness in mm, for the transmission.",
                    "default": 1.0,
                },
            },
            "required": ["formula", "density"],
        },
    },

    {
        "name": "calc_contrast_energy_scan",
        "description": (
            "Scan the anomalous X-ray contrast, absorption and transmission "
            "between two compounds across an energy range. Use this to plan "
            "an anomalous SAXS experiment: it answers 'which energy maximises "
            "the contrast' and shows where an absorption edge sits. Returns "
            "the arrays plus a 'best' summary giving the energy of maximum "
            "absolute contrast — always check 'transmission_sample' there, "
            "since the highest contrast is often at an unusably absorbing "
            "energy. Cost is linear in n_points; start with the default 200 "
            "over a wide range, then rescan a narrow range around the edge."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "formula_1": {
                    "type": "string",
                    "description": "Chemical formula of the first compound.",
                },
                "density_1": {
                    "type": "number",
                    "description": "Mass density of the first compound in g/cm^3.",
                },
                "formula_2": {
                    "type": "string",
                    "description": "Chemical formula of the second compound.",
                },
                "density_2": {
                    "type": "number",
                    "description": "Mass density of the second compound in g/cm^3.",
                },
                "e_start_keV": {
                    "type": "number",
                    "description": "Start of the energy range in keV.",
                },
                "e_end_keV": {
                    "type": "number",
                    "description": (
                        "End of the energy range in keV; must be greater than "
                        "e_start_keV."
                    ),
                },
                "n_points": {
                    "type": "integer",
                    "minimum": 2,
                    "maximum": 2000,
                    "description": "Number of energies to evaluate.",
                    "default": 200,
                },
                "name_1": {
                    "type": "string",
                    "description": "Optional display name for compound 1.",
                    "default": "",
                },
                "name_2": {
                    "type": "string",
                    "description": "Optional display name for compound 2.",
                    "default": "",
                },
                "mode": _MODE_PROPERTY,
                "isotopes_1": {
                    "type": ["object", "null"],
                    "description": _ISOTOPES_PROPERTY_DESC,
                },
                "isotopes_2": {
                    "type": ["object", "null"],
                    "description": _ISOTOPES_PROPERTY_DESC,
                },
                "thickness_mm": {
                    "type": "number",
                    "description": "Sample thickness in mm, for the transmission.",
                    "default": 1.0,
                },
                "vol_frac_1": {
                    "type": "number",
                    "minimum": 0.0,
                    "maximum": 1.0,
                    "description": (
                        "Volume fraction of compound 1, for the combined "
                        "sample transmission."
                    ),
                    "default": 0.01,
                },
                "max_points": {
                    "type": ["integer", "null"],
                    "description": (
                        "Decimation cap for the returned arrays (default 500). "
                        "The 'best' summary is computed before decimation, so "
                        "lowering this does not degrade it."
                    ),
                },
            },
            "required": [
                "formula_1", "density_1", "formula_2", "density_2",
                "e_start_keV", "e_end_keV",
            ],
        },
    },

    {
        "name": "lookup_element",
        "description": (
            "Look up an element's atomic number, atomic mass, neutron "
            "coherent scattering length b_c (in fm) and the list of isotopes "
            "that have neutron data. Use it to sanity-check a formula, or to "
            "find the isotope label to pass in an isotopes override for a "
            "neutron contrast-variation calculation."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "symbol": {
                    "type": "string",
                    "description": (
                        "Element symbol with a capitalised first letter, "
                        "e.g. 'Ti'."
                    ),
                },
            },
            "required": ["symbol"],
        },
    },

    {
        "name": "list_compound_library",
        "description": (
            "List the compounds the user has saved in their Scattering "
            "Contrast compound library. Read-only. Use it to reuse a "
            "material the user already defined in the GUI instead of asking "
            "them for the formula and density again. The library starts "
            "empty, so an empty list is normal."
        ),
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },

    {
        "name": "load_compound",
        "description": (
            "Load one saved compound definition from the user's Scattering "
            "Contrast library. Read-only — saving and deleting are not "
            "exposed. The returned formula_str, density, composition_mode and "
            "isotope_overrides can be passed straight into calc_compound or "
            "calc_contrast."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "name": {
                    "type": "string",
                    "description": (
                        "Compound name as returned by list_compound_library."
                    ),
                },
            },
            "required": ["name"],
        },
    },
]

# Convenience: look up a schema by name
CALCULATOR_SCHEMA_BY_NAME: dict[str, dict] = {
    t["name"]: t for t in CALCULATOR_TOOL_SCHEMAS
}
