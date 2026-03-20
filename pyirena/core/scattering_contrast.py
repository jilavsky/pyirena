"""
Scattering Contrast Calculator — core engine.

Computes X-ray and neutron scattering length densities (SLDs) and contrast
for two user-defined compounds.  No GUI/IO dependencies.

Physics:
  X-ray SLD (free electron):
      rho = (density * N_A / M_mol) * Z_total * r_e     [cm⁻²]
  X-ray SLD (anomalous, Chantler tables via xraydb):
      rho_anom = (density * N_A / M_mol) * f1_total * r_e
      f1_total = Σ n_i * f1_chantler(element_i, E_eV)   (includes Z)
  Neutron SLD:
      rho_n = B_total / V_molecule
      B_total = Σ n_i * b_c_i * 1e-13  [cm]   (b_c in fm from periodictable)
      V_molecule = M_mol / (density * N_A)      [cm³/molecule]

Units returned:
  SLD in 10¹⁰ cm⁻²   (common SAXS/SANS convention)
  SLD per gram in 10¹⁰ cm/g
  contrast (Δρ)² in 10²⁰ cm⁻⁴
  mu (linear absorption) in cm⁻¹
  transmission dimensionless [0, 1]
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

# Physical constants
_R_E = 2.8179403227e-13   # classical electron radius [cm]
_N_A = 6.02214076e23      # Avogadro's number [mol⁻¹]
_FM_TO_CM = 1.0e-13       # 1 fm = 1e-13 cm


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class CompoundProperties:
    """All computed properties for one compound (or vacuum)."""
    name: str = ""
    formula_str: str = ""
    composition_mode: str = "atomic_ratio"
    density: float = 0.0          # g/cm³

    # Molecular / formula-unit properties
    mol_weight: float = 0.0       # g/mol (per formula unit)
    weight_1mol: float = 0.0      # g per molecule/formula-unit
    n_mol_per_cm3: float = 0.0    # formula-units / cm³
    n_electrons_per_mol: float = 0.0   # electrons per formula-unit
    n_electrons_per_cm3: float = 0.0   # electrons / cm³
    volume_1mol: float = 0.0      # cm³ per formula-unit

    # X-ray (free electron approximation)
    xray_sld: float = 0.0         # 10¹⁰ cm⁻²
    xray_sld_per_gram: float = 0.0  # 10¹⁰ cm/g

    # Neutron
    neutron_total_b: float = 0.0  # cm per formula-unit
    neutron_sld: float = 0.0      # 10¹⁰ cm⁻²
    neutron_sld_per_gram: float = 0.0  # 10¹⁰ cm/g

    # Parsed element counts (stored for anomalous calculations)
    _element_counts: Dict[str, float] = field(default_factory=dict, repr=False)
    _isotope_overrides: Dict[str, str] = field(default_factory=dict, repr=False)


@dataclass
class AnomalousResult:
    """Anomalous X-ray results at a single energy for one compound."""
    energy_keV: float = 0.0
    xray_sld_anom: float = 0.0   # 10¹⁰ cm⁻²  (Chantler-corrected)
    mu_linear: float = 0.0       # cm⁻¹
    transmission: float = 1.0    # dimensionless


@dataclass
class ContrastResult:
    """Contrast between two compounds."""
    xray_contrast: float = 0.0       # 10²⁰ cm⁻⁴ (free electron)
    neutron_contrast: float = 0.0    # 10²⁰ cm⁻⁴
    ratio_xn: float = float('nan')   # X-ray / neutron contrast ratio
    # Anomalous (filled by compute_contrast_anomalous)
    xray_sld_anom_1: float = 0.0     # comp1 anomalous SLD
    xray_sld_anom_2: float = 0.0     # comp2 anomalous SLD
    xray_contrast_anom: float = 0.0  # 10²⁰ cm⁻⁴ (Chantler-corrected)
    mu_1: float = 0.0                # comp1 linear absorption [cm⁻¹]
    mu_2: float = 0.0                # comp2 linear absorption [cm⁻¹]
    transmission_1: float = 1.0      # comp1 transmission
    transmission_2: float = 1.0      # comp2 transmission
    transmission_sample: float = 1.0 # combined sample transmission


# Vacuum singleton (all zeros)
VACUUM = CompoundProperties(name="vacuum", formula_str="vacuum",
                             composition_mode="atomic_ratio", density=0.0)


# ---------------------------------------------------------------------------
# Formula parsing helpers
# ---------------------------------------------------------------------------

def get_isotopes_for_element(element_symbol: str) -> List[Tuple[str, float]]:
    """
    Return list of (label, b_c_fm) for all isotopes of an element that have
    neutron scattering data.  First entry is 'natural'.
    """
    import periodictable as pt
    try:
        el = pt.elements.symbol(element_symbol)
    except Exception:
        return []
    result = []
    # Natural
    if el.neutron is not None and el.neutron.b_c is not None:
        result.append(('natural', el.neutron.b_c))
    # Specific isotopes
    for iso in el:
        try:
            if iso.neutron is not None and iso.neutron.b_c is not None:
                result.append((str(iso.isotope), iso.neutron.b_c))
        except Exception:
            pass
    return result


def parse_formula(formula_str: str, mode: str = "atomic_ratio") -> Dict[str, float]:
    """
    Parse a compound formula string into a dict {element_symbol: molar_count}.

    Parameters
    ----------
    formula_str : str
        The compound specification.
        - atomic_ratio: standard formula, e.g. 'H2O', 'Fe3O4', 'SiO2'
        - weight_fraction_elements: numbers are weight fractions of elements,
          e.g. 'Au0.35Ag0.65'  (35 wt% Au, 65 wt% Ag)
        - weight_fraction_compounds: space-separated 'formula:frac' tokens,
          e.g. 'Y2O3:0.10 ZrO2:0.90'  (10 wt% Y2O3, 90 wt% ZrO2)
          Fractions can also be in percent, e.g. 'Y2O3:10% ZrO2:90%'

    mode : str
        One of 'atomic_ratio', 'weight_fraction_elements',
        'weight_fraction_compounds'.

    Returns
    -------
    dict[str, float]
        Molar counts per formula unit (may be un-normalized).

    Raises
    ------
    ValueError
        If the formula cannot be parsed or weights don't sum to ~1.
    """
    import periodictable as pt

    formula_str = formula_str.strip()
    if not formula_str:
        raise ValueError("Formula string is empty.")

    if mode == "atomic_ratio":
        try:
            f = pt.formula(formula_str)
        except Exception as e:
            raise ValueError(f"Cannot parse formula '{formula_str}': {e}")
        return {str(el): float(cnt) for el, cnt in f.atoms.items()}

    elif mode == "weight_fraction_elements":
        # Numbers in the formula are weight fractions, not atom counts.
        # Parse the raw formula to get the "atoms" dict, then convert.
        try:
            f = pt.formula(formula_str)
        except Exception as e:
            raise ValueError(f"Cannot parse formula '{formula_str}': {e}")
        raw = {str(el): float(cnt) for el, cnt in f.atoms.items()}
        # Validate that fractions sum to ~1
        total = sum(raw.values())
        if not (0.99 <= total <= 1.01):
            raise ValueError(
                f"Weight fractions sum to {total:.4f}, expected ~1.0. "
                "Use values that sum to 1 (e.g. 'Au0.35Ag0.65')."
            )
        # Convert wt_frac → molar count: n_i = wt_frac_i / M_i
        element_counts: Dict[str, float] = {}
        for sym, wf in raw.items():
            try:
                M_i = pt.elements.symbol(sym).mass
            except Exception:
                raise ValueError(f"Unknown element '{sym}'.")
            element_counts[sym] = wf / M_i
        return element_counts

    elif mode == "weight_fraction_compounds":
        return _parse_weight_fraction_compounds(formula_str)

    else:
        raise ValueError(f"Unknown composition mode: '{mode}'")


def _parse_weight_fraction_compounds(formula_str: str) -> Dict[str, float]:
    """
    Parse 'Compound1:frac1 Compound2:frac2 ...' into element molar counts.

    Fractions can be decimals (0.10) or percentages (10%).
    """
    import periodictable as pt

    tokens = formula_str.split()
    if not tokens:
        raise ValueError("No compound:fraction tokens found.")

    element_counts: Dict[str, float] = {}
    total_frac = 0.0

    for token in tokens:
        token = token.strip()
        if not token:
            continue
        # Split at last ':' to handle formulas like 'Ca3(PO4)2:0.5'
        if ':' not in token:
            raise ValueError(
                f"Token '{token}' missing ':fraction'. "
                "Expected format: 'formula:fraction', e.g. 'Y2O3:0.10 ZrO2:0.90'"
            )
        colon_idx = token.rfind(':')
        comp_str = token[:colon_idx].strip()
        frac_str = token[colon_idx + 1:].strip()

        # Parse fraction (accept 0.10 or 10%)
        frac_str = frac_str.rstrip('%').strip()
        try:
            frac = float(frac_str)
        except ValueError:
            raise ValueError(f"Cannot parse fraction '{frac_str}' in token '{token}'.")
        if frac_str.endswith('%') or (frac > 1.0 and not token.endswith('%')):
            # If fraction looks like a percent value (>1), convert
            if frac > 1.0:
                frac = frac / 100.0
        # Re-check after potential %-normalization
        # (if user wrote '10%' we already stripped %, so frac=10 → 0.10)
        # Actually we need to re-check the original frac_str before stripping
        orig_frac_str = token[colon_idx + 1:].strip()
        if '%' in orig_frac_str:
            frac = float(orig_frac_str.rstrip('%')) / 100.0

        total_frac += frac

        # Parse compound formula
        try:
            f = pt.formula(comp_str)
        except Exception as e:
            raise ValueError(f"Cannot parse compound formula '{comp_str}': {e}")

        M_compound = f.mass  # g/mol of this compound

        # Contribution: frac [g compound / g total] / M_compound [g/mol]
        # = moles of compound per gram of mixture
        moles_per_gram = frac / M_compound
        for el, cnt in f.atoms.items():
            sym = str(el)
            element_counts[sym] = element_counts.get(sym, 0.0) + cnt * moles_per_gram

    # Validate total fraction
    if not (0.98 <= total_frac <= 1.02):
        raise ValueError(
            f"Weight fractions sum to {total_frac:.4f}, expected ~1.0."
        )

    return element_counts


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

def compute_compound(
    formula_str: str,
    density: float,
    mode: str = "atomic_ratio",
    isotope_overrides: Optional[Dict[str, str]] = None,
    name: str = "",
) -> CompoundProperties:
    """
    Compute all scattering properties for one compound.

    Parameters
    ----------
    formula_str : str   Chemical formula / composition string.
    density : float     Density in g/cm³.
    mode : str          'atomic_ratio', 'weight_fraction_elements', or
                        'weight_fraction_compounds'.
    isotope_overrides : dict, optional
        Map element symbol → isotope mass number string, e.g. {'H': '2'}.
        Used only for neutron b_c lookup.
    name : str          Optional display name for the compound.

    Returns
    -------
    CompoundProperties
    """
    import periodictable as pt

    if isotope_overrides is None:
        isotope_overrides = {}

    element_counts = parse_formula(formula_str, mode)

    # --- Molecular weight ---
    # Use isotope mass if an override is specified (e.g. H→2 means use D mass).
    M_mol = 0.0
    n_electrons_per_mol = 0.0
    for sym, cnt in element_counts.items():
        el = pt.elements.symbol(sym)
        if sym in isotope_overrides:
            try:
                iso_num = int(isotope_overrides[sym])
                mass_i = el[iso_num].mass
            except Exception:
                mass_i = el.mass
        else:
            mass_i = el.mass
        M_mol += cnt * mass_i
        n_electrons_per_mol += cnt * el.number   # Z = number of electrons

    if M_mol <= 0:
        raise ValueError("Computed molecular weight is zero — check formula.")

    weight_1mol = M_mol / _N_A                    # g per formula-unit
    n_mol_per_cm3 = density * _N_A / M_mol if density > 0 else 0.0
    n_electrons_per_cm3 = n_electrons_per_mol * n_mol_per_cm3
    volume_1mol = 1.0 / n_mol_per_cm3 if n_mol_per_cm3 > 0 else 0.0

    # --- X-ray SLD (free electron approximation) ---
    xray_sld_cm2 = n_electrons_per_cm3 * _R_E        # [cm⁻²]
    xray_sld = xray_sld_cm2 / 1e10                   # [10¹⁰ cm⁻²]
    xray_sld_per_gram = xray_sld / density if density > 0 else 0.0

    # --- Neutron SLD ---
    neutron_total_b_fm = 0.0
    for sym, cnt in element_counts.items():
        bc = _get_neutron_b_c(sym, isotope_overrides)
        if bc is not None:
            neutron_total_b_fm += cnt * bc
        # Elements without data contribute 0 (with implicit warning)

    neutron_total_b_cm = neutron_total_b_fm * _FM_TO_CM   # [cm per formula-unit]
    neutron_sld_cm2 = neutron_total_b_cm / volume_1mol if volume_1mol > 0 else 0.0
    neutron_sld = neutron_sld_cm2 / 1e10                  # [10¹⁰ cm⁻²]
    neutron_sld_per_gram = neutron_sld / density if density > 0 else 0.0

    return CompoundProperties(
        name=name or formula_str,
        formula_str=formula_str,
        composition_mode=mode,
        density=density,
        mol_weight=M_mol,
        weight_1mol=weight_1mol,
        n_mol_per_cm3=n_mol_per_cm3,
        n_electrons_per_mol=n_electrons_per_mol,
        n_electrons_per_cm3=n_electrons_per_cm3,
        volume_1mol=volume_1mol,
        xray_sld=xray_sld,
        xray_sld_per_gram=xray_sld_per_gram,
        neutron_total_b=neutron_total_b_cm,
        neutron_sld=neutron_sld,
        neutron_sld_per_gram=neutron_sld_per_gram,
        _element_counts=element_counts,
        _isotope_overrides=isotope_overrides,
    )


def _get_neutron_b_c(element_symbol: str,
                     isotope_overrides: Dict[str, str]) -> Optional[float]:
    """Return coherent neutron scattering length b_c in fm for an element."""
    import periodictable as pt
    el = pt.elements.symbol(element_symbol)
    if element_symbol in isotope_overrides:
        iso_str = isotope_overrides[element_symbol]
        try:
            iso_num = int(iso_str)
            iso = el[iso_num]
            if iso.neutron is not None and iso.neutron.b_c is not None:
                return iso.neutron.b_c
        except (ValueError, KeyError, AttributeError):
            pass
    # Fall back to natural element
    if el.neutron is not None and el.neutron.b_c is not None:
        return el.neutron.b_c
    return None


# ---------------------------------------------------------------------------
# Contrast calculation
# ---------------------------------------------------------------------------

def compute_contrast(
    comp1: CompoundProperties,
    comp2: CompoundProperties,
) -> ContrastResult:
    """
    Compute free-electron X-ray and neutron contrast between two compounds.

    Returns ContrastResult (anomalous fields left at 0.0 — call
    compute_contrast_anomalous to fill them).
    """
    drho_xray = comp1.xray_sld - comp2.xray_sld       # 10¹⁰ cm⁻²
    drho_neut = comp1.neutron_sld - comp2.neutron_sld  # 10¹⁰ cm⁻²

    xray_contrast = drho_xray ** 2    # 10²⁰ cm⁻⁴
    neut_contrast = drho_neut ** 2    # 10²⁰ cm⁻⁴

    if neut_contrast > 0:
        ratio = xray_contrast / neut_contrast
    else:
        ratio = float('nan')

    return ContrastResult(
        xray_contrast=xray_contrast,
        neutron_contrast=neut_contrast,
        ratio_xn=ratio,
    )


# ---------------------------------------------------------------------------
# Anomalous X-ray calculations
# ---------------------------------------------------------------------------

def compute_anomalous(
    comp: CompoundProperties,
    energy_keV: float,
    thickness_mm: float = 1.0,
) -> AnomalousResult:
    """
    Compute anomalous X-ray SLD and linear absorption at one energy.

    Parameters
    ----------
    comp : CompoundProperties   Pre-computed compound (needs _element_counts).
    energy_keV : float          X-ray energy in keV.
    thickness_mm : float        Sample thickness in mm (for transmission).

    Returns
    -------
    AnomalousResult
    """
    import xraydb

    energy_eV = energy_keV * 1000.0
    element_counts = comp._element_counts

    if not element_counts or comp.density <= 0 or comp.mol_weight <= 0:
        return AnomalousResult(energy_keV=energy_keV)

    # Anomalous SLD: sum f1_i * count_i  (f1 includes Z in xraydb)
    f1_total = 0.0
    for sym, cnt in element_counts.items():
        try:
            f1 = xraydb.f1_chantler(sym, energy_eV)
            z = xraydb.atomic_number(sym)
            f1_total += cnt * (z + f1)   # total real scattering factor
        except Exception:
            # Fall back to Z (free electron) if element not in Chantler tables
            try:
                z = xraydb.atomic_number(sym)
                f1_total += cnt * z
            except Exception:
                pass

    n_mol_per_cm3 = comp.density * _N_A / comp.mol_weight
    xray_sld_anom_cm2 = n_mol_per_cm3 * f1_total * _R_E   # [cm⁻²]
    xray_sld_anom = xray_sld_anom_cm2 / 1e10               # [10¹⁰ cm⁻²]

    # Linear absorption: mu = density * Σ(w_i * mu_mass_i)
    mu_linear = 0.0
    for sym, cnt in element_counts.items():
        try:
            mu_mass = xraydb.mu_chantler(sym, energy_eV)   # [cm²/g]
        except Exception:
            continue
        import periodictable as pt
        try:
            M_i = pt.elements.symbol(sym).mass
        except Exception:
            continue
        w_i = cnt * M_i / comp.mol_weight   # weight fraction of element i
        mu_linear += comp.density * w_i * mu_mass   # [cm⁻¹]

    thickness_cm = thickness_mm / 10.0
    transmission = float(np.exp(-mu_linear * thickness_cm))

    return AnomalousResult(
        energy_keV=energy_keV,
        xray_sld_anom=xray_sld_anom,
        mu_linear=mu_linear,
        transmission=transmission,
    )


def compute_contrast_anomalous(
    comp1: CompoundProperties,
    comp2: CompoundProperties,
    energy_keV: float,
    thickness_mm: float = 1.0,
    vol_frac_comp1: float = 0.01,
) -> ContrastResult:
    """
    Compute full contrast including anomalous X-ray correction.

    Parameters
    ----------
    comp1, comp2 : CompoundProperties
    energy_keV : float   X-ray energy in keV.
    thickness_mm : float Sample thickness in mm.
    vol_frac_comp1 : float  Volume fraction of comp1 in sample (for combined T).

    Returns
    -------
    ContrastResult with all fields populated.
    """
    base = compute_contrast(comp1, comp2)

    anom1 = compute_anomalous(comp1, energy_keV, thickness_mm)
    anom2 = compute_anomalous(comp2, energy_keV, thickness_mm)

    drho_anom = anom1.xray_sld_anom - anom2.xray_sld_anom
    xray_contrast_anom = drho_anom ** 2

    # Combined sample transmission: T_sample = T1^vf * T2^(1-vf)
    # (assuming Beer-Lambert law for mixed phases)
    vf2 = 1.0 - vol_frac_comp1
    if anom1.mu_linear > 0 or anom2.mu_linear > 0:
        mu_combined = (vol_frac_comp1 * anom1.mu_linear +
                       vf2 * anom2.mu_linear)
        t_cm = thickness_mm / 10.0
        transmission_sample = float(np.exp(-mu_combined * t_cm))
    else:
        transmission_sample = 1.0

    base.xray_sld_anom_1 = anom1.xray_sld_anom
    base.xray_sld_anom_2 = anom2.xray_sld_anom
    base.xray_contrast_anom = xray_contrast_anom
    base.mu_1 = anom1.mu_linear
    base.mu_2 = anom2.mu_linear
    base.transmission_1 = anom1.transmission
    base.transmission_2 = anom2.transmission
    base.transmission_sample = transmission_sample

    return base


def compute_anomalous_scan(
    comp1: CompoundProperties,
    comp2: CompoundProperties,
    e_start_keV: float,
    e_end_keV: float,
    n_points: int = 500,
    thickness_mm: float = 1.0,
    vol_frac_comp1: float = 0.01,
) -> Dict[str, np.ndarray]:
    """
    Compute anomalous X-ray quantities over an energy range.

    Returns
    -------
    dict with keys:
      'energy'            [keV]  shape (n_points,)
      'xray_contrast_anom'       [10²⁰ cm⁻⁴]
      'xray_sld_anom_1'          [10¹⁰ cm⁻²]
      'xray_sld_anom_2'          [10¹⁰ cm⁻²]
      'mu_1'                     [cm⁻¹]
      'mu_2'                     [cm⁻¹]
      'transmission_1'
      'transmission_2'
      'transmission_sample'
    """
    import xraydb
    import periodictable as pt

    energies = np.linspace(e_start_keV, e_end_keV, n_points)
    energies_eV = energies * 1000.0

    def _sld_scan(comp: CompoundProperties) -> Tuple[np.ndarray, np.ndarray]:
        """Return (xray_sld_anom, mu_linear) arrays over energies_eV."""
        if not comp._element_counts or comp.density <= 0 or comp.mol_weight <= 0:
            return np.zeros(n_points), np.zeros(n_points)

        n_mol = comp.density * _N_A / comp.mol_weight

        # Pre-fetch per-element arrays
        elem_data = []
        for sym, cnt in comp._element_counts.items():
            try:
                f1_arr = np.array(
                    [xraydb.f1_chantler(sym, e) for e in energies_eV]
                )
                z = xraydb.atomic_number(sym)
                f_total_arr = z + f1_arr  # includes Z
            except Exception:
                try:
                    z = xraydb.atomic_number(sym)
                    f_total_arr = np.full(n_points, float(z))
                except Exception:
                    f_total_arr = np.zeros(n_points)

            try:
                mu_arr = np.array(
                    [xraydb.mu_chantler(sym, e) for e in energies_eV]
                )
            except Exception:
                mu_arr = np.zeros(n_points)

            try:
                M_i = pt.elements.symbol(sym).mass
            except Exception:
                M_i = 0.0

            w_i = cnt * M_i / comp.mol_weight
            elem_data.append((cnt, f_total_arr, mu_arr, w_i))

        sld_arr = np.zeros(n_points)
        mu_arr_total = np.zeros(n_points)
        for cnt, f_arr, mu_arr, w_i in elem_data:
            sld_arr += cnt * f_arr
            mu_arr_total += comp.density * w_i * mu_arr

        sld_arr = n_mol * sld_arr * _R_E / 1e10   # [10¹⁰ cm⁻²]
        return sld_arr, mu_arr_total

    sld1, mu1 = _sld_scan(comp1)
    sld2, mu2 = _sld_scan(comp2)

    contrast_anom = (sld1 - sld2) ** 2

    t_cm = thickness_mm / 10.0
    T1 = np.exp(-mu1 * t_cm)
    T2 = np.exp(-mu2 * t_cm)

    vf2 = 1.0 - vol_frac_comp1
    mu_sample = vol_frac_comp1 * mu1 + vf2 * mu2
    T_sample = np.exp(-mu_sample * t_cm)

    return {
        'energy': energies,
        'xray_contrast_anom': contrast_anom,
        'xray_sld_anom_1': sld1,
        'xray_sld_anom_2': sld2,
        'mu_1': mu1,
        'mu_2': mu2,
        'transmission_1': T1,
        'transmission_2': T2,
        'transmission_sample': T_sample,
    }


# ---------------------------------------------------------------------------
# Utility: element Z lookup (for display)
# ---------------------------------------------------------------------------

def get_element_info(symbol: str) -> Dict[str, object]:
    """Return basic element info dict for display purposes."""
    import periodictable as pt
    import xraydb
    try:
        el = pt.elements.symbol(symbol)
        z = el.number
        mass = el.mass
    except Exception:
        z, mass = 0, 0.0
    result = {'symbol': symbol, 'Z': z, 'mass': mass, 'neutron_b_c': None}
    if el.neutron is not None:
        result['neutron_b_c'] = el.neutron.b_c
    return result
