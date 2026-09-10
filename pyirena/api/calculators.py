"""pyirena.api.calculators — stateless support calculators for AI/scripting.

Unlike the rest of ``pyirena.api``, these functions do not read a measured
dataset: they answer *experiment-planning* questions from first principles.
The first group wraps :mod:`pyirena.core.scattering_contrast` — the same
engine behind the Scattering Contrast GUI panel — so an agent can compute a
scattering contrast and feed it straight into a Sizes or Modeling fit
(``set_shape(..., contrast=...)`` wants exactly the ``xray_contrast`` value
returned here, in 10^20 cm^-4).

Two deliberate deviations from the rest of the stateless api layer:

1. **No path arguments.** Nothing here takes a user-supplied path, so none
   of these functions go through ``pyirena.api._paths`` / the
   ``PYIRENA_DATA_ROOT`` sandbox. The only file touched is the user's own
   compound library at a fixed location in the home directory, and only
   for reading.
2. **Errors are returned, not raised.** ``discovery``/``data``/``results``
   raise on bad input, but calculators are reached through the MCP
   dispatcher alongside the control tools, and their failure modes (a typo
   in a chemical formula, an unknown element, a missing optional
   dependency) are all agent-recoverable. So they return the control
   layer's ``{"error", "suggestion", "code"}`` shape instead. The shape is
   rebuilt locally rather than imported from ``pyirena.api.control.errors``
   so that ``import pyirena.api`` does not pull in the control package —
   the same choice ``pyirena/mcp/dispatch.py`` makes.

Units follow the house convention: densities g/cm^3, energies keV,
thicknesses mm, scattering length densities 10^10 cm^-2, contrasts
10^20 cm^-4, linear absorption cm^-1.

Requires the optional ``pyirena[contrast]`` extra (``periodictable`` for
masses/Z/neutron b_c, ``xraydb`` additionally for anomalous quantities).
Both are imported lazily by the core, so importing this module never fails.
"""
from __future__ import annotations

from typing import Any, Optional

import numpy as np

from pyirena.api.schemas import array_to_list

# The three composition modes accepted by parse_formula(). The core only
# validates these with a bare `else: raise ValueError`, and the tuple itself
# lived only in the GUI panel (_MODE_KEYS), so it is published here.
COMPOSITION_MODES: tuple[str, ...] = (
    "atomic_ratio",
    "weight_fraction_elements",
    "weight_fraction_compounds",
)

# Guard rail for calc_contrast_energy_scan: the core calls xraydb's
# f1_chantler/mu_chantler once per (element, energy) in a Python loop, so
# cost is linear in n_points and a careless request is slow, not just large.
_MAX_SCAN_POINTS = 2000

# Fields of CompoundProperties worth returning. Excludes the two private
# fields (_element_counts / _isotope_overrides), which are core-internal
# plumbing for the anomalous calculations and carry no meaning to a caller.
_COMPOUND_FIELDS = (
    "name",
    "formula_str",
    "composition_mode",
    "density",
    "mol_weight",
    "weight_1mol",
    "n_mol_per_cm3",
    "n_electrons_per_mol",
    "n_electrons_per_cm3",
    "volume_1mol",
    "xray_sld",
    "xray_sld_per_gram",
    "neutron_total_b",
    "neutron_sld",
    "neutron_sld_per_gram",
)

# Quantities that are basis-dependent (and so not chemically meaningful) in
# the two weight-fraction modes, where parse_formula returns moles per gram
# of mixture rather than per formula unit.
_BASIS_DEPENDENT_FIELDS = (
    "mol_weight",
    "weight_1mol",
    "n_mol_per_cm3",
    "n_electrons_per_mol",
    "volume_1mol",
)


# ---------------------------------------------------------------------------
# Error helpers — same shape as pyirena.api.control.errors.make_error
# ---------------------------------------------------------------------------

def _error(message: str, suggestion: str = "", code: str = "ERROR") -> dict:
    return {"error": message, "suggestion": suggestion, "code": code}


def _missing_dependency(exc: ImportError) -> dict:
    return _error(
        f"An optional dependency for scattering-contrast calculations is "
        f"missing: {exc}.",
        suggestion="Install it with: pip install 'pyirena[contrast]'",
        code="MISSING_DEPENDENCY",
    )


def _bad_mode(mode: str) -> dict:
    return _error(
        f"Unknown composition mode '{mode}'.",
        suggestion=f"Use one of: {', '.join(COMPOSITION_MODES)}.",
        code="BAD_MODE",
    )


# Per-mode syntax reminder, used in error suggestions and tool descriptions.
_MODE_SYNTAX = {
    "atomic_ratio": "standard chemical notation, e.g. 'TiO', 'Ti2O3', 'SiO2'",
    "weight_fraction_elements": (
        "element symbols followed by weight fractions summing to 1, "
        "e.g. 'Au0.35Ag0.65' for 35 wt% Au / 65 wt% Ag"
    ),
    "weight_fraction_compounds": (
        "space-separated 'formula:fraction' tokens summing to 1, "
        "e.g. 'Y2O3:0.10 ZrO2:0.90' (percentages such as '10%' also work)"
    ),
}


def _bad_formula(formula: str, mode: str, exc: Exception) -> dict:
    return _error(
        f"Could not parse formula '{formula}' in mode '{mode}': {exc}",
        suggestion=f"For mode '{mode}' use {_MODE_SYNTAX[mode]}.",
        code="BAD_FORMULA",
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _compound_to_dict(comp: Any, mode: str) -> dict:
    """Flatten a core CompoundProperties into a JSON-safe dict."""
    out: dict[str, Any] = {}
    for name in _COMPOUND_FIELDS:
        value = getattr(comp, name, None)
        if isinstance(value, (float, int)) and not isinstance(value, bool):
            value = float(value)
            if not np.isfinite(value):
                value = None
        out[name] = value
    if mode in ("weight_fraction_elements", "weight_fraction_compounds"):
        out["basis_warning"] = (
            "In weight-fraction modes the composition is normalised per gram "
            "of mixture, so " + ", ".join(_BASIS_DEPENDENT_FIELDS) + " are "
            "basis-dependent and not per-formula-unit values. The scattering "
            "length densities and contrasts are unaffected."
        )
    return out


def _make_compound(
    formula: str,
    density: float,
    mode: str,
    name: str,
    isotopes: Optional[dict],
):
    """Build a CompoundProperties, or return (None, error_dict).

    Replicates the GUI panel's vacuum rule: an empty formula or a zero
    density means vacuum. A *fresh* zeroed instance is built rather than
    handing out the core's mutable module-level VACUUM singleton.
    """
    from pyirena.core.scattering_contrast import CompoundProperties, compute_compound

    if mode not in COMPOSITION_MODES:
        return None, _bad_mode(mode)

    formula = (formula or "").strip()
    if not formula or formula.lower() == "vacuum" or float(density) == 0.0:
        return CompoundProperties(
            name=name or "vacuum",
            formula_str="vacuum",
            composition_mode=mode,
            density=0.0,
        ), None

    if float(density) < 0.0:
        return None, _error(
            f"Density must be positive, got {density} g/cm^3.",
            suggestion="Pass the mass density in g/cm^3 (e.g. 4.95 for TiO).",
            code="BAD_DENSITY",
        )

    try:
        comp = compute_compound(
            formula_str=formula,
            density=float(density),
            mode=mode,
            isotope_overrides=isotopes or {},
            name=name or formula,
        )
    except ValueError as exc:
        return None, _bad_formula(formula, mode, exc)
    return comp, None


def _contrast_to_dict(res: Any, anomalous: bool) -> dict:
    """Flatten a core ContrastResult into a JSON-safe dict."""
    def _f(value):
        value = float(value)
        return value if np.isfinite(value) else None

    out = {
        "xray_contrast": _f(res.xray_contrast),
        "neutron_contrast": _f(res.neutron_contrast),
        "ratio_xray_neutron": _f(res.ratio_xn),
        "units": {
            "xray_contrast": "10^20 cm^-4",
            "neutron_contrast": "10^20 cm^-4",
            "sld": "10^10 cm^-2",
        },
    }
    if anomalous:
        out.update({
            "xray_sld_anom_1": _f(res.xray_sld_anom_1),
            "xray_sld_anom_2": _f(res.xray_sld_anom_2),
            "xray_contrast_anom": _f(res.xray_contrast_anom),
            "mu_1": _f(res.mu_1),
            "mu_2": _f(res.mu_2),
            "transmission_1": _f(res.transmission_1),
            "transmission_2": _f(res.transmission_2),
            "transmission_sample": _f(res.transmission_sample),
        })
        out["units"].update({
            "mu": "cm^-1",
            "transmission": "dimensionless (0-1)",
        })
    return out


# ---------------------------------------------------------------------------
# Public calculators
# ---------------------------------------------------------------------------

def calc_compound(
    formula: str,
    density: float,
    mode: str = "atomic_ratio",
    name: str = "",
    isotopes: Optional[dict] = None,
    energy_keV: Optional[float] = None,
    thickness_mm: float = 1.0,
) -> dict:
    """Compute scattering length densities for a single compound.

    Parameters
    ----------
    formula : str
        Chemical formula, e.g. "SiO2". An empty string (or a zero density)
        is treated as vacuum and yields all-zero SLDs.
    density : float
        Mass density in g/cm^3.
    mode : str
        Composition mode; one of COMPOSITION_MODES. In the two
        weight-fraction modes the numbers are fractions summing to 1, and
        the per-formula-unit quantities become basis-dependent (a
        "basis_warning" key is added in that case).
    name : str
        Display name. Defaults to the formula.
    isotopes : dict, optional
        Isotope overrides as {element_symbol: mass_number_string},
        e.g. {"H": "2"} for deuterium. Affects molar mass and the neutron
        scattering length, not the electron count.
    energy_keV : float, optional
        When given, adds the anomalous (Chantler-corrected) X-ray SLD, the
        linear absorption coefficient and the transmission at this energy.
    thickness_mm : float
        Sample thickness in mm, used only for the transmission.

    Returns
    -------
    dict
        Compound properties (see _COMPOUND_FIELDS), plus an "anomalous"
        sub-dict when energy_keV is given. On failure, a dict with an
        "error" key.
    """
    try:
        comp, err = _make_compound(formula, density, mode, name, isotopes)
    except ImportError as exc:
        return _missing_dependency(exc)
    if err is not None:
        return err

    out = _compound_to_dict(comp, mode)
    out["units"] = {
        "density": "g/cm^3",
        "xray_sld": "10^10 cm^-2",
        "neutron_sld": "10^10 cm^-2",
        "mol_weight": "g/mol",
    }

    if energy_keV is not None:
        try:
            from pyirena.core.scattering_contrast import compute_anomalous
            anom = compute_anomalous(comp, float(energy_keV), float(thickness_mm))
        except ImportError as exc:
            return _missing_dependency(exc)
        out["anomalous"] = {
            "energy_keV": float(anom.energy_keV),
            "xray_sld_anom": float(anom.xray_sld_anom),
            "mu_linear": float(anom.mu_linear),
            "transmission": float(anom.transmission),
            "thickness_mm": float(thickness_mm),
            "units": {
                "xray_sld_anom": "10^10 cm^-2",
                "mu_linear": "cm^-1",
                "transmission": "dimensionless (0-1)",
            },
        }
    return out


def calc_contrast(
    formula_1: str,
    density_1: float,
    formula_2: str,
    density_2: float,
    name_1: str = "",
    name_2: str = "",
    mode: str = "atomic_ratio",
    isotopes_1: Optional[dict] = None,
    isotopes_2: Optional[dict] = None,
    energy_keV: Optional[float] = None,
    thickness_mm: float = 1.0,
    vol_frac_1: float = 0.01,
) -> dict:
    """Compute the scattering contrast between two compounds.

    This is the number a Sizes or Modeling fit needs: the returned
    "xray_contrast" is (delta-rho)^2 in 10^20 cm^-4, the same units and
    convention as the ``contrast`` parameter of ``set_shape`` and of a
    Modeling population.

    Use an empty formula (or zero density) for either compound to represent
    vacuum — that gives the contrast of the other compound against vacuum.

    Parameters
    ----------
    formula_1, formula_2 : str
        Chemical formulas, e.g. "TiO" and "Ti2O3".
    density_1, density_2 : float
        Mass densities in g/cm^3.
    name_1, name_2 : str
        Optional display names.
    mode : str
        Composition mode applied to both formulas; one of COMPOSITION_MODES.
    isotopes_1, isotopes_2 : dict, optional
        Isotope overrides per compound, e.g. {"H": "2"} for deuteration.
        Relevant to the neutron contrast.
    energy_keV : float, optional
        When given, also returns the anomalous (Chantler-corrected) X-ray
        contrast, the linear absorption of each compound, and the sample
        transmission. Needed near an absorption edge, where the
        free-electron approximation is wrong.
    thickness_mm : float
        Sample thickness in mm, for the transmission calculation.
    vol_frac_1 : float
        Volume fraction of compound 1 in the sample, used only to combine
        the two absorptions into "transmission_sample".

    Returns
    -------
    dict
        "xray_contrast" and "neutron_contrast" (10^20 cm^-4), their ratio,
        the two compounds under "compound_1"/"compound_2", and the
        anomalous quantities when energy_keV is given. On failure, a dict
        with an "error" key.
    """
    try:
        comp1, err = _make_compound(formula_1, density_1, mode, name_1, isotopes_1)
        if err is not None:
            return err
        comp2, err = _make_compound(formula_2, density_2, mode, name_2, isotopes_2)
        if err is not None:
            return err

        from pyirena.core.scattering_contrast import (
            compute_contrast,
            compute_contrast_anomalous,
        )
        if energy_keV is None:
            res = compute_contrast(comp1, comp2)
        else:
            res = compute_contrast_anomalous(
                comp1, comp2,
                energy_keV=float(energy_keV),
                thickness_mm=float(thickness_mm),
                vol_frac_comp1=float(vol_frac_1),
            )
    except ImportError as exc:
        return _missing_dependency(exc)

    out = _contrast_to_dict(res, anomalous=energy_keV is not None)
    out["compound_1"] = _compound_to_dict(comp1, mode)
    out["compound_2"] = _compound_to_dict(comp2, mode)
    if energy_keV is not None:
        out["energy_keV"] = float(energy_keV)
        out["thickness_mm"] = float(thickness_mm)
        out["vol_frac_1"] = float(vol_frac_1)
    return out


def calc_contrast_energy_scan(
    formula_1: str,
    density_1: float,
    formula_2: str,
    density_2: float,
    e_start_keV: float,
    e_end_keV: float,
    n_points: int = 200,
    name_1: str = "",
    name_2: str = "",
    mode: str = "atomic_ratio",
    isotopes_1: Optional[dict] = None,
    isotopes_2: Optional[dict] = None,
    thickness_mm: float = 1.0,
    vol_frac_1: float = 0.01,
    max_points: Optional[int] = None,
) -> dict:
    """Scan anomalous X-ray contrast, absorption and transmission vs energy.

    Answers "which energy maximises the contrast" and "where is the edge"
    — the anomalous-SAXS experiment-planning question. Returns the arrays
    plus a "best" summary giving the energy of maximum absolute contrast.

    Parameters
    ----------
    formula_1, formula_2 : str
        Chemical formulas.
    density_1, density_2 : float
        Mass densities in g/cm^3.
    e_start_keV, e_end_keV : float
        Energy range in keV; e_end_keV must be greater than e_start_keV.
    n_points : int
        Number of energies to evaluate (default 200, maximum 2000). Cost is
        linear in this: each point does a table lookup per element.
    name_1, name_2, mode, isotopes_1, isotopes_2 : see calc_contrast
    thickness_mm : float
        Sample thickness in mm, for transmission.
    vol_frac_1 : float
        Volume fraction of compound 1, for the combined sample transmission.
    max_points : int, optional
        Decimation cap applied to the returned arrays (default:
        PYIRENA_MAX_ARRAY_POINTS, 500). The "best" summary is computed on
        the full-resolution scan before decimation.

    Returns
    -------
    dict
        "energy" plus "xray_contrast_anom", "xray_sld_anom_1/2", "mu_1/2"
        and "transmission_1/2/sample" as lists, a "best" summary, and
        "units". On failure, a dict with an "error" key.
    """
    if float(e_end_keV) <= float(e_start_keV):
        return _error(
            f"Energy range is empty or reversed: e_start_keV={e_start_keV}, "
            f"e_end_keV={e_end_keV}.",
            suggestion="Pass e_end_keV greater than e_start_keV, both in keV.",
            code="BAD_ENERGY_RANGE",
        )
    n_points = int(n_points)
    if n_points < 2:
        return _error(
            f"n_points must be at least 2, got {n_points}.",
            suggestion="Use n_points=200 for a typical scan.",
            code="BAD_N_POINTS",
        )
    if n_points > _MAX_SCAN_POINTS:
        return _error(
            f"n_points={n_points} exceeds the maximum of {_MAX_SCAN_POINTS}.",
            suggestion=(
                f"Use n_points<={_MAX_SCAN_POINTS}; 200 is usually enough to "
                "locate an edge, then rescan a narrow range around it."
            ),
            code="BAD_N_POINTS",
        )

    try:
        comp1, err = _make_compound(formula_1, density_1, mode, name_1, isotopes_1)
        if err is not None:
            return err
        comp2, err = _make_compound(formula_2, density_2, mode, name_2, isotopes_2)
        if err is not None:
            return err

        from pyirena.core.scattering_contrast import compute_anomalous_scan
        scan = compute_anomalous_scan(
            comp1, comp2,
            e_start_keV=float(e_start_keV),
            e_end_keV=float(e_end_keV),
            n_points=n_points,
            thickness_mm=float(thickness_mm),
            vol_frac_comp1=float(vol_frac_1),
        )
    except ImportError as exc:
        return _missing_dependency(exc)

    out: dict[str, Any] = {
        "compound_1": _compound_to_dict(comp1, mode),
        "compound_2": _compound_to_dict(comp2, mode),
        "e_start_keV": float(e_start_keV),
        "e_end_keV": float(e_end_keV),
        "n_points": n_points,
        "thickness_mm": float(thickness_mm),
        "vol_frac_1": float(vol_frac_1),
    }
    for key, arr in scan.items():
        out[key] = array_to_list(arr, max_points=max_points)

    # Summary on the full-resolution scan, before decimation.
    energy = np.asarray(scan["energy"], dtype=float)
    contrast = np.asarray(scan["xray_contrast_anom"], dtype=float)
    finite = np.isfinite(contrast)
    if finite.any():
        idx = int(np.nanargmax(np.abs(np.where(finite, contrast, np.nan))))
        out["best"] = {
            "energy_keV": float(energy[idx]),
            "xray_contrast_anom": float(contrast[idx]),
            "transmission_sample": float(
                np.asarray(scan["transmission_sample"], dtype=float)[idx]
            ),
            "note": (
                "Energy of maximum |anomalous X-ray contrast| over the scanned "
                "range. Check transmission_sample is workable at this energy."
            ),
        }
    out["units"] = {
        "energy": "keV",
        "xray_contrast_anom": "10^20 cm^-4",
        "xray_sld_anom_1": "10^10 cm^-2",
        "xray_sld_anom_2": "10^10 cm^-2",
        "mu_1": "cm^-1",
        "mu_2": "cm^-1",
        "transmission_1": "dimensionless (0-1)",
        "transmission_2": "dimensionless (0-1)",
        "transmission_sample": "dimensionless (0-1)",
    }
    return out


def lookup_element(symbol: str) -> dict:
    """Look up an element's Z, atomic mass, neutron b_c and isotopes.

    Useful for sanity-checking a formula, or for picking an isotope for a
    neutron contrast-variation calculation before calling calc_contrast
    with an isotope override.

    Parameters
    ----------
    symbol : str
        Element symbol, e.g. "Ti". Case-sensitive in the standard way
        ("Ti", not "TI").

    Returns
    -------
    dict
        "symbol", "Z", "mass" (g/mol), "neutron_b_c" (fm, may be None) and
        "isotopes" — a list of {"label", "neutron_b_c"} where "label" is
        either "natural" or a mass number usable as an isotope override.
        On failure, a dict with an "error" key.
    """
    symbol = (symbol or "").strip()
    if not symbol:
        return _error(
            "No element symbol given.",
            suggestion="Pass an element symbol such as 'Ti'.",
            code="BAD_ELEMENT",
        )
    try:
        import periodictable as pt
    except ImportError as exc:
        return _missing_dependency(exc)

    # Validate first: the core's get_element_info() leaves `el` unbound in
    # its except branch and raises NameError on an unknown symbol.
    try:
        pt.elements.symbol(symbol)
    except Exception:
        return _error(
            f"Unknown element symbol '{symbol}'.",
            suggestion=(
                "Use a standard symbol with a capitalised first letter, "
                "e.g. 'Ti', 'Si', 'O'."
            ),
            code="BAD_ELEMENT",
        )

    from pyirena.core.scattering_contrast import (
        get_element_info,
        get_isotopes_for_element,
    )
    info = get_element_info(symbol)
    b_c = info.get("neutron_b_c")
    return {
        "symbol": info.get("symbol"),
        "Z": int(info.get("Z") or 0),
        "mass": float(info.get("mass") or 0.0),
        "neutron_b_c": None if b_c is None else float(b_c),
        "isotopes": [
            {"label": label, "neutron_b_c": float(value)}
            for label, value in get_isotopes_for_element(symbol)
        ],
        "units": {"mass": "g/mol", "neutron_b_c": "fm"},
    }


def list_compound_library() -> dict:
    """List the compounds saved in the user's contrast compound library.

    Read-only. The library is the one written by the Scattering Contrast
    GUI panel; it lives in the user's home directory and starts empty.

    Returns
    -------
    dict
        "compounds" (a list of names, possibly empty) and "library_path".
    """
    from pyirena.io.contrast_io import DEFAULT_LIBRARY_PATH, list_compounds_in_library

    return {
        "compounds": list(list_compounds_in_library()),
        "library_path": str(DEFAULT_LIBRARY_PATH),
    }


def load_compound(name: str) -> dict:
    """Load one saved compound definition from the user's library.

    Read-only — saving and deleting are deliberately not exposed. The
    returned "formula_str", "density", "composition_mode" and
    "isotope_overrides" can be passed straight into calc_compound or
    calc_contrast.

    Parameters
    ----------
    name : str
        Compound name as returned by list_compound_library().

    Returns
    -------
    dict
        The stored definition plus any cached computed properties, or a
        dict with an "error" key if the library or the compound is missing.
    """
    from pyirena.io.contrast_io import (
        DEFAULT_LIBRARY_PATH,
        list_compounds_in_library,
        load_compound_from_library,
    )

    try:
        stored = load_compound_from_library(name)
    except FileNotFoundError:
        return _error(
            f"No compound library found at {DEFAULT_LIBRARY_PATH}.",
            suggestion=(
                "Save a compound from the Scattering Contrast GUI panel "
                "first, or pass the formula and density directly to "
                "calc_contrast."
            ),
            code="NO_LIBRARY",
        )
    except KeyError:
        available = list_compounds_in_library()
        return _error(
            f"No compound named '{name}' in the library.",
            suggestion=(
                f"Available: {', '.join(available)}."
                if available
                else "The library is empty."
            ),
            code="UNKNOWN_COMPOUND",
        )

    out = {k: v for k, v in stored.items()}
    out["library_path"] = str(DEFAULT_LIBRARY_PATH)
    return out
