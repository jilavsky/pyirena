"""
Scattering Contrast Calculator — I/O module.

Provides:
  - HDF5 compound library (save, load, list, delete)
  - CSV export of computed properties
  - HDF5 export of energy-scan arrays
"""

from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import h5py
import numpy as np

DEFAULT_LIBRARY_PATH = Path.home() / ".pyirena" / "contrast_compounds.h5"

# HDF5 group path for the compound library
_LIBRARY_ROOT = "compounds"


# ---------------------------------------------------------------------------
# Compound library — HDF5
# ---------------------------------------------------------------------------

def list_compounds_in_library(library_path: Optional[Path] = None) -> List[str]:
    """Return sorted list of compound names stored in the library file."""
    library_path = Path(library_path or DEFAULT_LIBRARY_PATH)
    if not library_path.exists():
        return []
    try:
        with h5py.File(library_path, 'r') as f:
            if _LIBRARY_ROOT not in f:
                return []
            return sorted(f[_LIBRARY_ROOT].keys())
    except Exception:
        return []


def save_compound_to_library(
    compound_dict: Dict[str, Any],
    library_path: Optional[Path] = None,
    name: Optional[str] = None,
) -> str:
    """
    Save a compound definition to the HDF5 library.

    Parameters
    ----------
    compound_dict : dict
        Must contain keys: 'formula_str', 'composition_mode', 'density',
        'isotope_overrides' (dict), and optionally 'name'.
        May also contain computed properties from CompoundProperties.
    library_path : Path, optional
        Path to the HDF5 library file.  Defaults to DEFAULT_LIBRARY_PATH.
    name : str, optional
        Compound name (used as HDF5 group key).  Falls back to 'name' in
        compound_dict, then to formula_str.

    Returns
    -------
    str   The name used to store the compound.
    """
    library_path = Path(library_path or DEFAULT_LIBRARY_PATH)
    library_path.parent.mkdir(parents=True, exist_ok=True)

    compound_name = (
        name
        or compound_dict.get('name')
        or compound_dict.get('formula_str', 'unknown')
    )
    # Sanitise for HDF5 group name (no slashes etc.)
    compound_name = compound_name.replace('/', '_').replace('\\', '_').strip()
    if not compound_name:
        compound_name = "unnamed"

    mode = 'a' if library_path.exists() else 'w'
    with h5py.File(library_path, mode) as f:
        root = f.require_group(_LIBRARY_ROOT)

        # Delete existing group with same name so we can overwrite
        if compound_name in root:
            del root[compound_name]

        grp = root.require_group(compound_name)

        # Metadata attributes
        grp.attrs['name'] = compound_name
        grp.attrs['formula_str'] = compound_dict.get('formula_str', '')
        grp.attrs['composition_mode'] = compound_dict.get('composition_mode', 'atomic_ratio')
        grp.attrs['density'] = float(compound_dict.get('density', 0.0))
        grp.attrs['saved_at'] = datetime.now().isoformat()
        # Isotope overrides stored as JSON string
        grp.attrs['isotope_overrides'] = json.dumps(
            compound_dict.get('isotope_overrides', {})
        )

        # Computed properties as datasets
        _float_fields = [
            'mol_weight', 'weight_1mol', 'n_mol_per_cm3',
            'n_electrons_per_mol', 'n_electrons_per_cm3', 'volume_1mol',
            'xray_sld', 'xray_sld_per_gram',
            'neutron_total_b', 'neutron_sld', 'neutron_sld_per_gram',
        ]
        computed = grp.require_group('computed')
        for field in _float_fields:
            val = compound_dict.get(field)
            if val is not None:
                computed.create_dataset(field, data=float(val))

    return compound_name


def load_compound_from_library(
    name: str,
    library_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """
    Load a compound from the HDF5 library.

    Returns
    -------
    dict with all stored keys (suitable for passing to compute_compound or
    restoring CompoundProperties).
    """
    library_path = Path(library_path or DEFAULT_LIBRARY_PATH)
    if not library_path.exists():
        raise FileNotFoundError(f"Library not found: {library_path}")

    with h5py.File(library_path, 'r') as f:
        if _LIBRARY_ROOT not in f or name not in f[_LIBRARY_ROOT]:
            raise KeyError(f"Compound '{name}' not found in library.")
        grp = f[f'{_LIBRARY_ROOT}/{name}']

        result: Dict[str, Any] = {
            'name': grp.attrs.get('name', name),
            'formula_str': grp.attrs.get('formula_str', ''),
            'composition_mode': grp.attrs.get('composition_mode', 'atomic_ratio'),
            'density': float(grp.attrs.get('density', 0.0)),
            'saved_at': grp.attrs.get('saved_at', ''),
        }
        # Decode isotope overrides
        iso_json = grp.attrs.get('isotope_overrides', '{}')
        try:
            result['isotope_overrides'] = json.loads(iso_json)
        except Exception:
            result['isotope_overrides'] = {}

        # Computed properties
        if 'computed' in grp:
            cgrp = grp['computed']
            for key in cgrp.keys():
                result[key] = float(cgrp[key][()])

    return result


def delete_compound_from_library(
    name: str,
    library_path: Optional[Path] = None,
) -> None:
    """Remove a compound group from the HDF5 library."""
    library_path = Path(library_path or DEFAULT_LIBRARY_PATH)
    if not library_path.exists():
        return
    with h5py.File(library_path, 'a') as f:
        if _LIBRARY_ROOT in f and name in f[_LIBRARY_ROOT]:
            del f[f'{_LIBRARY_ROOT}/{name}']


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

def export_results_csv(results: Dict[str, Any], filepath: Path) -> None:
    """
    Export computed properties for both compounds to a CSV file.

    Parameters
    ----------
    results : dict
        Expected keys (all optional):
          'comp1': dict of CompoundProperties fields
          'comp2': dict of CompoundProperties fields
          'contrast': dict of ContrastResult fields
          'anomalous': dict of AnomalousResult fields per compound
    filepath : Path
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    c1 = results.get('comp1', {})
    c2 = results.get('comp2', {})
    ct = results.get('contrast', {})

    rows = [
        ('Property', 'Units', 'Compound 1', 'Compound 2'),
        ('Name', '', c1.get('name', ''), c2.get('name', '')),
        ('Formula', '', c1.get('formula_str', ''), c2.get('formula_str', '')),
        ('Density', 'g/cm³', _fmt(c1.get('density')), _fmt(c2.get('density'))),
        ('', '', '', ''),
        ('--- Basic Properties ---', '', '', ''),
        ('Molecular weight', 'g/mol', _fmt(c1.get('mol_weight')), _fmt(c2.get('mol_weight'))),
        ('Weight of 1 molecule', 'g', _fmte(c1.get('weight_1mol')), _fmte(c2.get('weight_1mol'))),
        ('Number of molecules/cm³', '1/cm³', _fmte(c1.get('n_mol_per_cm3')), _fmte(c2.get('n_mol_per_cm3'))),
        ('Electrons per formula unit', '', _fmt(c1.get('n_electrons_per_mol')), _fmt(c2.get('n_electrons_per_mol'))),
        ('Electrons per cm³', '1/cm³', _fmte(c1.get('n_electrons_per_cm3')), _fmte(c2.get('n_electrons_per_cm3'))),
        ('Volume of 1 formula unit', 'cm³', _fmte(c1.get('volume_1mol')), _fmte(c2.get('volume_1mol'))),
        ('', '', '', ''),
        ('--- X-ray (free electron) ---', '', '', ''),
        ('X-ray SLD (rho)', '10^10 cm^-2', _fmt4(c1.get('xray_sld')), _fmt4(c2.get('xray_sld'))),
        ('X-ray SLD per gram', '10^10 cm/g', _fmt4(c1.get('xray_sld_per_gram')), _fmt4(c2.get('xray_sld_per_gram'))),
        ('', '', '', ''),
        ('--- Neutron ---', '', '', ''),
        ('Total neutron b (molecule)', 'cm', _fmte(c1.get('neutron_total_b')), _fmte(c2.get('neutron_total_b'))),
        ('Neutron SLD (rho)', '10^10 cm^-2', _fmt4(c1.get('neutron_sld')), _fmt4(c2.get('neutron_sld'))),
        ('Neutron SLD per gram', '10^10 cm/g', _fmt4(c1.get('neutron_sld_per_gram')), _fmt4(c2.get('neutron_sld_per_gram'))),
        ('', '', '', ''),
        ('--- Contrast ---', '', '', ''),
        ('X-ray contrast (delta-rho)^2', '10^20 cm^-4', _fmt4(ct.get('xray_contrast')), ''),
        ('Neutron contrast (delta-rho)^2', '10^20 cm^-4', _fmt4(ct.get('neutron_contrast')), ''),
        ('X-ray / Neutron contrast ratio', '', _fmt4(ct.get('ratio_xn')), ''),
    ]

    # Anomalous section
    anom = results.get('anomalous', {})
    if anom:
        rows += [
            ('', '', '', ''),
            ('--- Anomalous X-ray (Chantler) ---', '', '', ''),
            ('Energy', 'keV', _fmt(anom.get('energy_keV')), ''),
            ('Anomalous X-ray SLD', '10^10 cm^-2',
             _fmt4(ct.get('xray_sld_anom_1')), _fmt4(ct.get('xray_sld_anom_2'))),
            ('Linear absorption coefficient', 'cm^-1',
             _fmt4(ct.get('mu_1')), _fmt4(ct.get('mu_2'))),
            ('Transmission', '',
             _fmt4(ct.get('transmission_1')), _fmt4(ct.get('transmission_2'))),
            ('Sample transmission', '', _fmt4(ct.get('transmission_sample')), ''),
            ('Anomalous X-ray contrast (delta-rho)^2', '10^20 cm^-4',
             _fmt4(ct.get('xray_contrast_anom')), ''),
        ]

    with open(filepath, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerows(rows)


def _fmt(val) -> str:
    if val is None:
        return ''
    try:
        return f'{float(val):.6g}'
    except Exception:
        return str(val)


def _fmt4(val) -> str:
    if val is None:
        return ''
    try:
        return f'{float(val):.4g}'
    except Exception:
        return str(val)


def _fmte(val) -> str:
    if val is None:
        return ''
    try:
        return f'{float(val):.4e}'
    except Exception:
        return str(val)


# ---------------------------------------------------------------------------
# HDF5 energy-scan export
# ---------------------------------------------------------------------------

def save_scan_to_hdf5(
    scan_data: Dict[str, np.ndarray],
    filepath: Path,
    group_name: str = "contrast_scan",
    metadata: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Save energy-scan arrays to HDF5.

    Parameters
    ----------
    scan_data : dict
        Keys: 'energy', 'xray_contrast_anom', 'xray_sld_anom_1',
              'xray_sld_anom_2', 'mu_1', 'mu_2',
              'transmission_1', 'transmission_2', 'transmission_sample'.
    filepath : Path
    group_name : str   HDF5 group name.
    metadata : dict, optional  Stored as group attributes.
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    mode = 'a' if filepath.exists() else 'w'
    with h5py.File(filepath, mode) as f:
        if group_name in f:
            del f[group_name]
        grp = f.require_group(group_name)
        grp.attrs['saved_at'] = datetime.now().isoformat()
        grp.attrs['description'] = 'Scattering contrast energy scan (pyirena)'
        if metadata:
            for k, v in metadata.items():
                try:
                    grp.attrs[k] = v
                except Exception:
                    grp.attrs[k] = str(v)

        units = {
            'energy': 'keV',
            'xray_contrast_anom': '10^20 cm^-4',
            'xray_sld_anom_1': '10^10 cm^-2',
            'xray_sld_anom_2': '10^10 cm^-2',
            'mu_1': 'cm^-1',
            'mu_2': 'cm^-1',
            'transmission_1': '',
            'transmission_2': '',
            'transmission_sample': '',
        }

        for key, arr in scan_data.items():
            if arr is not None:
                ds = grp.create_dataset(key, data=np.asarray(arr, dtype='f8'),
                                        compression='gzip')
                if key in units:
                    ds.attrs['units'] = units[key]


def export_scan_csv(
    scan_data: Dict[str, np.ndarray],
    filepath: Path,
) -> None:
    """Export energy-scan data to CSV."""
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    keys_ordered = [
        'energy', 'xray_sld_anom_1', 'xray_sld_anom_2', 'xray_contrast_anom',
        'mu_1', 'mu_2', 'transmission_1', 'transmission_2', 'transmission_sample',
    ]
    header = [
        'Energy [keV]',
        'Xray SLD comp1 [10^10 cm^-2]',
        'Xray SLD comp2 [10^10 cm^-2]',
        'Xray Contrast (anom) [10^20 cm^-4]',
        'Mu comp1 [cm^-1]',
        'Mu comp2 [cm^-1]',
        'Transmission comp1',
        'Transmission comp2',
        'Sample Transmission',
    ]

    available = [k for k in keys_ordered if k in scan_data]
    avail_header = [header[keys_ordered.index(k)] for k in available]

    n = len(scan_data.get('energy', []))
    with open(filepath, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerow(avail_header)
        for i in range(n):
            row = [f'{scan_data[k][i]:.6g}' for k in available]
            writer.writerow(row)
