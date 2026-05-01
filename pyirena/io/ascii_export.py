"""
ASCII export for SAS data files.

Exports HDF5 (NXcanSAS) datasets and any pyirena fit-result groups they
contain to plain-text .dat files, suitable for legacy data-analysis tools
(old Igor Pro, Fortran codes, gnuplot, …).

Output convention
-----------------
- Default delimiter: single space (LoadWave/J-friendly in Igor Pro;
  also accepted by awk, np.loadtxt, gnuplot).  Comma optional.
- Default precision: 7 significant figures (single-precision-safe;
  old Fortran code parsed double-precision exponent strings poorly).
  12 sig figs available for tools that need full double-precision.
- Header lines start with '# ', followed by 'key = value [unit]'.
  Capped at 25 lines.  Last header line is always '# Columns = ...'.
- For USAXS files, the desmeared sasdata group is preferred; if only
  the slit-smeared variant (group name ending in _SMR) is present,
  the export refuses (raises) — old tools cannot interpret the
  resolution columns correctly.
- Negative or NaN values are rendered as the literal token `nan`
  (space delimiter) or as an empty field (comma delimiter).
- All ASCII-only header text — Unicode characters in source data
  (Å, χ², ±, …) are rewritten to ASCII (A, chi^2, +/-).

Public entry point
------------------
``export_dataset_to_ascii(h5_path, out_dir, **opts)`` returns a manifest
dict listing the files written and any skipped suffixes.  Designed for
use by the Data Selector GUI; safe to call programmatically.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from pyirena.io.hdf5 import find_matching_groups


# Tools whose result groups can be exported as a 4-column model file.
# Tuple = (group path, loader function name, header function, curves function)
# Loaders are imported lazily so that a missing optional dependency in one
# tool does not break the whole module.
_MODEL_SUFFIX_GROUPS = {
    "unif": "entry/unified_fit_results",
    "simp": "entry/simple_fit_results",
    "mod":  "entry/modeling_results",
    "sd":   "entry/sizes_results",
    "waxs": "entry/waxs_peakfit_results",
}

_MAX_HEADER_LINES = 25


# ──────────────────────────────────────────────────────────────────────────────
# Formatting helpers
# ──────────────────────────────────────────────────────────────────────────────

def _ascii_safe(s) -> str:
    """Replace non-ASCII characters that frequently appear in pyirena strings."""
    if s is None:
        return ""
    if isinstance(s, bytes):
        try:
            s = s.decode("utf-8")
        except Exception:
            s = s.decode("ascii", errors="replace")
    text = str(s)
    return (text
            .replace("Å^-1", "1/A")
            .replace("Å⁻¹", "1/A")
            .replace("Å", "A")
            .replace("χ²", "chi^2")
            .replace("χ2", "chi^2")
            .replace("±", "+/-")
            .replace("²", "^2")
            .replace("³", "^3")
            .replace("⁻¹", "^-1")
            .replace("μ", "u"))


def _format_value(v, precision: int, delimiter: str) -> str:
    """Format a single numeric value for column output."""
    if v is None:
        return "" if delimiter == "," else "nan"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return _ascii_safe(v)
    if np.isnan(f) or np.isinf(f):
        return "" if delimiter == "," else "nan"
    return f"{f:.{precision}g}"


def _format_scalar_for_header(v, precision: int = 7) -> str:
    """Format a scalar for the metadata header (always lowercase 'g').

    Default precision matches the user's default column precision so that
    values shown in the header resolve to the same accuracy.
    """
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    try:
        f = float(v)
        if np.isnan(f) or np.isinf(f):
            return ""
        return f"{f:.{precision}g}"
    except (TypeError, ValueError):
        return _ascii_safe(v)


# ──────────────────────────────────────────────────────────────────────────────
# Primary data loading (USAXS desmeared / SAXS / WAXS)
# ──────────────────────────────────────────────────────────────────────────────

def _attr_value(group, key, default=None):
    """Read attribute, decoding bytes if necessary."""
    if key not in group.attrs:
        return default
    v = group.attrs[key]
    if isinstance(v, bytes):
        return v.decode("utf-8", errors="replace")
    return v


def _dataset_value(h5file, path):
    """Read dataset value at path or return None.  Decodes bytes."""
    if path not in h5file:
        return None
    obj = h5file[path]
    if not isinstance(obj, h5py.Dataset):
        return None
    v = obj[()]
    if isinstance(v, bytes):
        return v.decode("utf-8", errors="replace")
    if hasattr(v, "size") and v.size == 1 and not isinstance(v, str):
        try:
            scalar = v.item()
            if isinstance(scalar, bytes):
                return scalar.decode("utf-8", errors="replace")
            return scalar
        except Exception:
            return v
    return v


def _pick_desmeared_entry(f) -> Optional[str]:
    """
    Return the path of the desmeared NXcanSAS subentry (the one whose group
    name does NOT end in '_SMR').  Returns None if no SAS entry exists.

    For files containing only an _SMR variant, returns None — caller decides
    whether to raise or proceed.
    """
    sas_entries = find_matching_groups(
        f,
        required_attributes={"canSAS_class": "SASentry", "NX_class": "NXsubentry"},
        required_items={"definition": "NXcanSAS"},
    )
    if not sas_entries:
        return None

    # Prefer entries that do not end in '_SMR'
    desmeared = [e for e in sas_entries if not e.split("/")[-1].endswith("_SMR")]
    if desmeared:
        # Honour entry/@default if it points to one of the desmeared candidates
        if "entry" in f:
            default_name = _attr_value(f["entry"], "default", default=None)
            if default_name:
                preferred = f"entry/{default_name}"
                if preferred in desmeared:
                    return preferred
        return desmeared[0]
    return None


def _read_sasdata(f, entry_path: str) -> dict:
    """
    Read Q / I / Idev / Qdev arrays from the sasdata group inside *entry_path*.

    Honours NXcanSAS attribute conventions (signal, I_axes, uncertainties,
    resolutions) so that non-standard array names still resolve.
    """
    sasdata_path = f"{entry_path}/sasdata"
    if sasdata_path not in f:
        # Some hand-built files have the data arrays directly under entry/
        sasdata_path = entry_path
    sd = f[sasdata_path]

    signal = _attr_value(sd, "signal", default="I")
    i_axes = _attr_value(sd, "I_axes", default="Q")

    I = _dataset_value(f, f"{sasdata_path}/{signal}")
    Q = _dataset_value(f, f"{sasdata_path}/{i_axes}")

    # Intensity attrs (units, uncertainties pointer, instrument-specific extras)
    int_attrs = {}
    if f"{sasdata_path}/{signal}" in f:
        for k, v in f[f"{sasdata_path}/{signal}"].attrs.items():
            int_attrs[k] = v.decode("utf-8", errors="replace") if isinstance(v, bytes) else v

    units_I = int_attrs.get("units")

    err_key = int_attrs.get("uncertainties")
    Idev = _dataset_value(f, f"{sasdata_path}/{err_key}") if err_key else None

    # Q resolution (Qdev for desmeared, dQw + dQl for SMR)
    q_attrs = {}
    if f"{sasdata_path}/{i_axes}" in f:
        for k, v in f[f"{sasdata_path}/{i_axes}"].attrs.items():
            q_attrs[k] = v.decode("utf-8", errors="replace") if isinstance(v, bytes) else v
    res_key = q_attrs.get("resolutions")
    Qdev = _dataset_value(f, f"{sasdata_path}/{res_key}") if res_key else None

    return {
        "Q": np.asarray(Q, dtype=float) if Q is not None else None,
        "I": np.asarray(I, dtype=float) if I is not None else None,
        "Idev": np.asarray(Idev, dtype=float) if Idev is not None else None,
        "Qdev": np.asarray(Qdev, dtype=float) if Qdev is not None else None,
        "units_I": units_I,
        "units_Q": q_attrs.get("units"),
        "Kfactor": int_attrs.get("Kfactor"),
        "OmegaFactor": int_attrs.get("OmegaFactor"),
        "blankname": int_attrs.get("blankname"),
        "thickness": int_attrs.get("thickness"),
        "label": int_attrs.get("label"),
    }


def _read_metadata(f, entry_path: str) -> dict:
    """
    Read sample / instrument / metadata into a flat dict for the header.

    Tolerant of missing groups — pyirena-created files often have only
    sasdata, no instrument/sample sub-trees.
    """
    md: dict = {}

    # Sample
    if "entry/sample" in f and isinstance(f["entry/sample"], h5py.Group):
        s = f["entry/sample"]
        md["sample_name"] = _dataset_value(f, "entry/sample/name")
        md["sample_title"] = _dataset_value(f, "entry/sample/title")
        md["sample_thickness"] = _dataset_value(f, "entry/sample/thickness")
    # Sub-entry title (per-sample)
    md["sample_title"] = md.get("sample_title") or _dataset_value(f, f"{entry_path}/title")

    # Instrument
    md["instrument_name"] = _dataset_value(f, "entry/instrument/name")
    md["wavelength_A"]  = _dataset_value(f, "entry/instrument/monochromator/wavelength")
    md["energy_keV"]    = _dataset_value(f, "entry/instrument/monochromator/energy")
    if md["wavelength_A"] is None:
        md["wavelength_A"] = _dataset_value(f, "entry/instrument/wavelength")
    if md["energy_keV"] is None:
        md["energy_keV"] = _dataset_value(f, "entry/instrument/energy")

    # USAXS metadata
    if md["energy_keV"] is None:
        md["energy_keV"] = _dataset_value(f, "entry/metadata/DCM_energy")

    # Root file attrs
    md["instrument_root"] = _attr_value(f, "instrument", default=None)
    md["file_time"]       = _attr_value(f, "file_time", default=None)

    return md


def _load_primary_data(h5_path: Path) -> dict:
    """
    Open *h5_path*, locate the desmeared sasdata group, return a dict with
    Q, I, dI plus header metadata.  Raises on hard failure.

    Strategy
    --------
    1. Look for a proper NXcanSAS SASentry group.  Prefer the desmeared
       variant (group name without ``_SMR`` suffix) when both exist.
    2. Fall back to "simple" HDF5 layouts with Q/I/Idev under one of a
       few well-known paths (entry1/data1, entry/data, …) — the same
       paths that ``readSimpleHDF5`` in ``hdf5.py`` recognises.
    """
    h5_path = Path(h5_path)
    with h5py.File(h5_path, "r") as f:
        entry_path = _pick_desmeared_entry(f)
        if entry_path is not None:
            sd = _read_sasdata(f, entry_path)
            if sd["Q"] is None or sd["I"] is None:
                raise ValueError(f"sasdata at {entry_path} missing Q or I")
            md = _read_metadata(f, entry_path)
            md["source_group"] = entry_path
            return {**sd, **md}

        # No NXcanSAS group — check if only an _SMR variant is present
        sas_entries = find_matching_groups(
            f,
            required_attributes={"canSAS_class": "SASentry", "NX_class": "NXsubentry"},
            required_items={"definition": "NXcanSAS"},
        )
        if sas_entries:
            raise ValueError("only slit-smeared (_SMR) data present; export refuses")

        # Fallback: simple HDF5 layouts (entry1/data1/Q,I,Idev etc.)
        for base in ("entry1/data1", "entry/data", "data", ""):
            q_path = f"{base}/Q".strip("/")
            i_path = f"{base}/I".strip("/")
            if q_path in f and i_path in f:
                Q  = np.asarray(f[q_path][()], dtype=float)
                I  = np.asarray(f[i_path][()], dtype=float)
                dI = None
                for err_name in ("Idev", "Error", "error", "I_error"):
                    p = f"{base}/{err_name}".strip("/")
                    if p in f:
                        dI = np.asarray(f[p][()], dtype=float)
                        break
                return {
                    "Q": Q, "I": I, "Idev": dI, "Qdev": None,
                    "units_I": None, "units_Q": None,
                    "Kfactor": None, "OmegaFactor": None,
                    "blankname": None, "thickness": None, "label": None,
                    "sample_name": None, "sample_title": None,
                    "sample_thickness": None, "instrument_name": None,
                    "instrument_root": _attr_value(f, "instrument", default=None),
                    "wavelength_A": None, "energy_keV": None,
                    "file_time": _attr_value(f, "file_time", default=None),
                    "source_group": (base or "/"),
                }

        raise ValueError("No NXcanSAS entry found")


# ──────────────────────────────────────────────────────────────────────────────
# Header formatting
# ──────────────────────────────────────────────────────────────────────────────

def _format_data_header(
    h5_path: Path,
    primary: dict,
    extra_lines: Optional[list] = None,
    columns_label: str = "Q  I  dI",
    notes: Optional[str] = None,
) -> list:
    """Build the metadata header for {stem}.dat or {stem}_*.dat files."""
    extra_lines = extra_lines or []
    lines = ["# pyirena ASCII export"]
    lines.append(f"# Source file = {_ascii_safe(h5_path.name)}")
    lines.append(f"# Source group = {_ascii_safe(primary.get('source_group', ''))}  (desmeared)")
    lines.append(f"# Export time = {datetime.now().isoformat(timespec='seconds')}")

    sample_name = primary.get("sample_name") or ""
    if sample_name:
        lines.append(f"# Sample name = {_ascii_safe(sample_name)}")
    sample_title = primary.get("sample_title") or ""
    if sample_title and sample_title != sample_name:
        lines.append(f"# Sample title = {_ascii_safe(sample_title)}")

    thickness = primary.get("sample_thickness") or primary.get("thickness")
    if thickness is not None:
        lines.append(f"# Thickness mm = {_format_scalar_for_header(thickness)}")
    blankname = primary.get("blankname")
    if blankname:
        lines.append(f"# Blank name = {_ascii_safe(blankname)}")

    instr = primary.get("instrument_name") or primary.get("instrument_root") or "unknown"
    lines.append(f"# Instrument = {_ascii_safe(instr)}")

    if primary.get("energy_keV") is not None:
        lines.append(f"# Energy keV = {_format_scalar_for_header(primary['energy_keV'])}")
    if primary.get("wavelength_A") is not None:
        lines.append(f"# Wavelength A = {_format_scalar_for_header(primary['wavelength_A'])}")
    if primary.get("Kfactor") is not None:
        lines.append(f"# Kfactor = {_format_scalar_for_header(primary['Kfactor'])}")
    if primary.get("OmegaFactor") is not None:
        lines.append(f"# OmegaFactor = {_format_scalar_for_header(primary['OmegaFactor'])}")

    Q = primary.get("Q")
    if Q is not None and len(Q):
        lines.append(f"# Q points = {len(Q)}")
        q_pos = Q[Q > 0]
        q_min = float(q_pos.min()) if q_pos.size else float(Q.min())
        q_max = float(Q.max())
        lines.append(
            f"# Q range 1/A = {_format_scalar_for_header(q_min)} to "
            f"{_format_scalar_for_header(q_max)}"
        )

    units_I = primary.get("units_I") or "1/cm"
    units_Q = primary.get("units_Q") or "1/A"
    lines.append(f"# Units I = {_ascii_safe(units_I)}   Units Q = {_ascii_safe(units_Q)}   "
                 f"Units dI = {_ascii_safe(units_I)}")

    if notes:
        lines.append(f"# Notes = {_ascii_safe(notes)}")

    # Append any model-specific lines (truncate if needed)
    for line in extra_lines:
        if not line.startswith("#"):
            line = "# " + line
        lines.append(_ascii_safe(line))

    # Cap at _MAX_HEADER_LINES - 1 (reserve last for Columns =)
    if len(lines) > _MAX_HEADER_LINES - 1:
        lines = lines[:_MAX_HEADER_LINES - 2]
        lines.append("# (additional fields truncated; see HDF5 file for full metadata)")

    lines.append(f"# Columns = {columns_label}")
    return lines


# ──────────────────────────────────────────────────────────────────────────────
# Model header builders — one per tool acronym
# ──────────────────────────────────────────────────────────────────────────────

def _hdr_unif(uf: dict) -> list:
    """Header lines for Unified Fit results."""
    lines = ["# Model = Unified Fit"]
    n_levels = int(uf.get("num_levels", len(uf.get("levels", []))))
    bg = uf.get("background", None)
    chi2 = uf.get("chi_squared", None)
    parts = [f"Levels = {n_levels}"]
    if bg is not None:
        parts.append(f"Background 1/cm = {_format_scalar_for_header(bg)}")
    if chi2 is not None:
        parts.append(f"chi^2 = {_format_scalar_for_header(chi2)}")
    lines.append("# " + "   ".join(parts))

    levels = uf.get("levels", []) or []
    for i, lev in enumerate(levels, start=1):
        parts = []
        for key in ("G", "Rg", "B", "P", "ETA", "PACK", "RgCutoff"):
            if key in lev:
                err = lev.get(f"{key}_err")
                if err is not None and float(err) > 0:
                    parts.append(f"{key} = {_format_scalar_for_header(lev[key])} +/- "
                                 f"{_format_scalar_for_header(err)}")
                else:
                    parts.append(f"{key} = {_format_scalar_for_header(lev[key])}")
        if "correlated" in lev:
            parts.append(f"correlated = {lev['correlated']}")
        lines.append(f"# L{i}: " + "   ".join(parts))
    return lines


def _hdr_simp(sf: dict) -> list:
    """Header lines for Simple Fits results."""
    lines = ["# Model = Simple Fit (" + _ascii_safe(sf.get("model", "")) + ")"]
    chi2 = sf.get("chi_squared")
    rchi2 = sf.get("reduced_chi_squared")
    parts = []
    if chi2 is not None:
        parts.append(f"chi^2 = {_format_scalar_for_header(chi2)}")
    if rchi2 is not None:
        parts.append(f"reduced chi^2 = {_format_scalar_for_header(rchi2)}")
    if sf.get("dof") is not None:
        parts.append(f"dof = {int(sf['dof'])}")
    if parts:
        lines.append("# " + "   ".join(parts))

    params = sf.get("params", {}) or {}
    stds = sf.get("params_std", {}) or {}
    for name, val in params.items():
        err = stds.get(name)
        if err is not None and np.isfinite(err) and err > 0:
            lines.append(f"# {name} = {_format_scalar_for_header(val)} +/- "
                         f"{_format_scalar_for_header(err)}")
        else:
            lines.append(f"# {name} = {_format_scalar_for_header(val)}")

    derived = sf.get("derived", {}) or {}
    for name, val in derived.items():
        lines.append(f"# derived_{name} = {_format_scalar_for_header(val)}")
    return lines


def _hdr_mod(mod: dict) -> list:
    """Header lines for Modeling results."""
    lines = ["# Model = Modeling (parametric forward)"]
    parts = []
    if mod.get("chi_squared") is not None:
        parts.append(f"chi^2 = {_format_scalar_for_header(mod['chi_squared'])}")
    if mod.get("reduced_chi_squared") is not None:
        parts.append(f"reduced chi^2 = {_format_scalar_for_header(mod['reduced_chi_squared'])}")
    if mod.get("background") is not None:
        parts.append(f"background = {_format_scalar_for_header(mod['background'])}")
    if parts:
        lines.append("# " + "   ".join(parts))

    pops = mod.get("populations", []) or []
    lines.append(f"# Populations = {len(pops)}")
    for k, p in enumerate(pops, start=1):
        ptype = p.get("pop_type", "size_dist")
        label = _ascii_safe(p.get("label", ""))
        if ptype == "unified_level":
            extras = []
            for key in ("G", "Rg", "B", "P"):
                if p.get(key) is not None:
                    extras.append(f"{key}={_format_scalar_for_header(p[key])}")
            lines.append(f"# P{k} [unified_level] {label}: " + " ".join(extras))
        elif ptype == "diffraction_peak":
            lines.append(
                f"# P{k} [{_ascii_safe(p.get('peak_type', 'gaussian'))}] {label}: "
                f"Q0={_format_scalar_for_header(p.get('position'))} "
                f"A={_format_scalar_for_header(p.get('amplitude'))} "
                f"FWHM={_format_scalar_for_header(p.get('width'))}"
            )
        else:
            ff = _ascii_safe(p.get("form_factor", ""))
            sf = _ascii_safe(p.get("structure_factor", ""))
            dt = _ascii_safe(p.get("dist_type", ""))
            scale = p.get("scale")
            contrast = p.get("contrast")
            extras = [f"dist={dt}", f"ff={ff}", f"sf={sf}"]
            if scale is not None:
                extras.append(f"scale={_format_scalar_for_header(scale)}")
            if contrast is not None:
                extras.append(f"contrast={_format_scalar_for_header(contrast)}")
            lines.append(f"# P{k} [size_dist] {label}: " + " ".join(extras))
    return lines


def _hdr_sd(sz: dict) -> list:
    """Header lines for Size Distribution results."""
    lines = ["# Model = Size Distribution"]
    parts = []
    if sz.get("method") is not None:
        parts.append(f"method = {_ascii_safe(sz['method'])}")
    if sz.get("shape") is not None:
        parts.append(f"shape = {_ascii_safe(sz['shape'])}")
    if sz.get("contrast") is not None:
        parts.append(f"contrast = {_format_scalar_for_header(sz['contrast'])}")
    if parts:
        lines.append("# " + "   ".join(parts))

    fit_parts = []
    if sz.get("chi_squared") is not None:
        fit_parts.append(f"chi^2 = {_format_scalar_for_header(sz['chi_squared'])}")
    if sz.get("volume_fraction") is not None:
        fit_parts.append(f"volume fraction = {_format_scalar_for_header(sz['volume_fraction'])}")
    if sz.get("rg") is not None:
        fit_parts.append(f"Rg A = {_format_scalar_for_header(sz['rg'])}")
    if sz.get("n_iterations") is not None:
        fit_parts.append(f"iterations = {int(sz['n_iterations'])}")
    if fit_parts:
        lines.append("# " + "   ".join(fit_parts))

    grid_parts = []
    if sz.get("r_min") is not None and sz.get("r_max") is not None:
        grid_parts.append(
            f"r range A = {_format_scalar_for_header(sz['r_min'])} to "
            f"{_format_scalar_for_header(sz['r_max'])}"
        )
    if sz.get("n_bins") is not None:
        grid_parts.append(f"bins = {int(sz['n_bins'])}")
    if sz.get("background") is not None:
        grid_parts.append(f"background = {_format_scalar_for_header(sz['background'])}")
    if grid_parts:
        lines.append("# " + "   ".join(grid_parts))

    lines.append("# Note: distribution P(r) and N(r) are stored in the HDF5 file only")
    return lines


def _hdr_waxs(wp: dict) -> list:
    """Header lines for WAXS Peak Fit results."""
    lines = ["# Model = WAXS Peak Fit"]
    bg_shape = _ascii_safe(wp.get("bg_shape", ""))
    n_peaks = int(wp.get("n_peaks", 0))
    parts = [f"peaks = {n_peaks}", f"bg_shape = {bg_shape}"]
    if wp.get("chi_squared") is not None:
        parts.append(f"chi^2 = {_format_scalar_for_header(wp['chi_squared'])}")
    lines.append("# " + "   ".join(parts))

    # Background coefficients
    bg = wp.get("bg_params", {}) or {}
    bg_std = wp.get("bg_params_std", {}) or {}
    if bg:
        bg_strs = []
        for name, info in bg.items():
            val = info.get("value") if isinstance(info, dict) else info
            err = bg_std.get(name)
            if err is not None and np.isfinite(err) and err > 0:
                bg_strs.append(f"{name}={_format_scalar_for_header(val)}+/-"
                               f"{_format_scalar_for_header(err)}")
            else:
                bg_strs.append(f"{name}={_format_scalar_for_header(val)}")
        lines.append("# Background: " + " ".join(bg_strs))

    # Per-peak parameters (one line each)
    peaks = wp.get("peaks", []) or []
    peaks_std = wp.get("peaks_std", []) or [{} for _ in peaks]
    for i, peak in enumerate(peaks, start=1):
        shape = _ascii_safe(peak.get("shape", "Gauss"))
        pstd = peaks_std[i - 1] if i - 1 < len(peaks_std) else {}
        params_strs = []
        for pn in ("Q0", "A", "FWHM", "eta"):
            if pn in peak:
                pd = peak[pn]
                val = pd.get("value") if isinstance(pd, dict) else pd
                err = pstd.get(pn)
                if err is not None and np.isfinite(err) and err > 0:
                    params_strs.append(f"{pn}={_format_scalar_for_header(val)}+/-"
                                       f"{_format_scalar_for_header(err)}")
                else:
                    params_strs.append(f"{pn}={_format_scalar_for_header(val)}")
        lines.append(f"# Peak {i} [{shape}]: " + " ".join(params_strs))
    return lines


# ──────────────────────────────────────────────────────────────────────────────
# Model curves extraction — return (Q, I_model, I_data, dI) for 4-col file
# ──────────────────────────────────────────────────────────────────────────────

def _curves_unif(uf: dict, primary: dict):
    Q = uf.get("Q")
    I_model = uf.get("intensity_model")
    I_data = uf.get("intensity_data")
    dI = uf.get("intensity_error")
    if dI is None:
        dI = _resample_dI(primary, Q)
    return Q, I_model, I_data, dI


def _curves_simp(sf: dict, primary: dict):
    Q = sf.get("Q")
    I_model = sf.get("I_model")
    I_data = sf.get("intensity_data")
    dI = sf.get("intensity_error")
    if I_data is None:
        I_data = _resample_I(primary, Q)
    if dI is None:
        dI = _resample_dI(primary, Q)
    return Q, I_model, I_data, dI


def _curves_mod(mod: dict, primary: dict):
    Q = mod.get("model_q")
    I_model = mod.get("model_I")
    I_data = _resample_I(primary, Q)
    dI = _resample_dI(primary, Q)
    return Q, I_model, I_data, dI


def _curves_sd(sz: dict, primary: dict):
    Q = sz.get("Q")
    I_model = sz.get("intensity_model")
    I_data = sz.get("intensity_data")
    dI = sz.get("intensity_error")
    if dI is None:
        dI = _resample_dI(primary, Q)
    return Q, I_model, I_data, dI


def _curves_waxs(wp: dict, primary: dict):
    Q = wp.get("Q")
    I_model = wp.get("I_fit")
    I_data = wp.get("intensity_data")
    dI = wp.get("intensity_error")
    if I_data is None:
        I_data = _resample_I(primary, Q)
    if dI is None:
        dI = _resample_dI(primary, Q)
    return Q, I_model, I_data, dI


def _resample_I(primary: dict, Q_target):
    """Log-log linear interpolation of primary I onto Q_target."""
    if Q_target is None:
        return None
    Q_src = primary.get("Q")
    I_src = primary.get("I")
    if Q_src is None or I_src is None:
        return None
    return _interp_loglog(Q_target, Q_src, I_src)


def _resample_dI(primary: dict, Q_target):
    """Log-log linear interpolation of primary dI onto Q_target."""
    if Q_target is None:
        return None
    Q_src = primary.get("Q")
    dI_src = primary.get("Idev")
    if Q_src is None or dI_src is None:
        return None
    return _interp_loglog(Q_target, Q_src, dI_src)


def _interp_loglog(q_new, q_src, y_src):
    """Linear interpolation in log10 space; out-of-range returns NaN."""
    q_new = np.asarray(q_new, dtype=float)
    q_src = np.asarray(q_src, dtype=float)
    y_src = np.asarray(y_src, dtype=float)
    mask = (q_src > 0) & (y_src > 0) & np.isfinite(q_src) & np.isfinite(y_src)
    if mask.sum() < 2:
        return np.full_like(q_new, np.nan)
    lq_src = np.log10(q_src[mask])
    ly_src = np.log10(y_src[mask])
    out = np.full_like(q_new, np.nan, dtype=float)
    pos = q_new > 0
    if not np.any(pos):
        return out
    lq_new = np.log10(q_new[pos])
    ly_new = np.interp(lq_new, lq_src, ly_src, left=np.nan, right=np.nan)
    out[pos] = 10.0 ** ly_new
    return out


# ──────────────────────────────────────────────────────────────────────────────
# Lazy loader dispatch
# ──────────────────────────────────────────────────────────────────────────────

def _load_model(suffix: str, h5_path: Path):
    """
    Load a result group, returning the loader's dict or None if missing.

    Loader functions are imported lazily and may raise either KeyError /
    ValueError (group missing) or other exceptions (e.g. UnicodeEncodeError
    from a print() call inside the loader on Windows cp1252 consoles).  We
    catch broadly so a transient loader bug never aborts the export.

    Returns None for any failure — the caller distinguishes "no result"
    from "format error" by re-checking the HDF5 group existence directly
    when needed; for our purposes both are equivalent (skip silently).
    """
    try:
        if suffix == "unif":
            from pyirena.io.nxcansas_unified import load_unified_fit_results
            return load_unified_fit_results(h5_path)
        if suffix == "simp":
            from pyirena.io.nxcansas_simple_fits import load_simple_fit_results
            return load_simple_fit_results(h5_path)
        if suffix == "mod":
            from pyirena.io.nxcansas_modeling import load_modeling_results
            return load_modeling_results(h5_path)
        if suffix == "sd":
            from pyirena.io.nxcansas_sizes import load_sizes_results
            return load_sizes_results(h5_path)
        if suffix == "waxs":
            from pyirena.io.nxcansas_waxs_peakfit import load_waxs_peakfit_results
            return load_waxs_peakfit_results(h5_path)
    except Exception:
        return None
    return None


_HDR_DISPATCH = {
    "unif": _hdr_unif, "simp": _hdr_simp, "mod": _hdr_mod,
    "sd": _hdr_sd, "waxs": _hdr_waxs,
}
_CURVES_DISPATCH = {
    "unif": _curves_unif, "simp": _curves_simp, "mod": _curves_mod,
    "sd": _curves_sd, "waxs": _curves_waxs,
}


# ──────────────────────────────────────────────────────────────────────────────
# File writing
# ──────────────────────────────────────────────────────────────────────────────

def _write_dat(path: Path, columns: list, header_lines: list,
               delimiter: str, precision: int) -> None:
    """Write a .dat file with given header (already prefixed with '# ') and columns."""
    n_rows = len(columns[0])
    n_cols = len(columns)
    sep = delimiter
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        for line in header_lines:
            fh.write(line + "\n")
        for r in range(n_rows):
            row_strs = [_format_value(columns[c][r], precision, delimiter)
                        for c in range(n_cols)]
            fh.write(sep.join(row_strs) + "\n")


# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────

def export_dataset_to_ascii(
    h5_path,
    out_dir,
    *,
    delimiter: str = " ",
    precision: int = 7,
    include_header: bool = True,
    include_models: bool = True,
    model_flags: Optional[dict] = None,
    error_fraction: float = 0.05,
) -> dict:
    """
    Export one HDF5 NXcanSAS file to ASCII.

    Parameters
    ----------
    h5_path : Path or str
        Source HDF5 file.
    out_dir : Path or str
        Output directory; created if missing.
    delimiter : {' ', ','}
        Column separator.  Default single space.
    precision : int
        Significant figures for numeric columns and parameters in the header.
        Default 7 (single-precision-safe).
    include_header : bool
        When False, the .dat files contain only data lines (no '#' header).
    include_models : bool
        When False, only the primary {stem}.dat file is written; per-model
        .dat files are skipped even if model_flags asks for them.
    model_flags : dict or None
        Map of suffix → bool.  Suffixes: 'unif', 'simp', 'mod', 'sd', 'waxs'.
        When None, all five are False (only primary data is exported).
    error_fraction : float
        When the source HDF5 has no Idev dataset, dI is generated as
        intensity * error_fraction.

    Returns
    -------
    dict
        ``{'data': Path|None, 'models': [Path...], 'skipped': [(suffix, reason)...]}``
    """
    h5_path = Path(h5_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if model_flags is None:
        model_flags = {}

    manifest = {"data": None, "models": [], "skipped": []}

    # 1. Primary data
    primary = _load_primary_data(h5_path)
    Q  = primary["Q"]
    I  = primary["I"]
    dI = primary.get("Idev")
    notes = None
    if dI is None:
        dI = I * float(error_fraction)
        notes = (f"dI generated as I * {error_fraction:.4f} "
                 f"(no measurement uncertainty in source)")

    stem = h5_path.stem
    data_path = out_dir / f"{stem}.dat"

    if include_header:
        header = _format_data_header(h5_path, primary,
                                     extra_lines=None,
                                     columns_label="Q  I  dI",
                                     notes=notes)
    else:
        header = []
    _write_dat(data_path, [Q, I, dI], header, delimiter, precision)
    manifest["data"] = data_path

    # 2. Model results (one .dat per enabled checkbox that has data)
    if include_models:
        # Pre-check which result groups exist so we can distinguish
        # "missing" (silent skip, expected) from "format error" (worth noting).
        with h5py.File(h5_path, "r") as f:
            present = {
                suffix: (path in f)
                for suffix, path in _MODEL_SUFFIX_GROUPS.items()
            }
        for suffix in ("unif", "simp", "mod", "sd", "waxs"):
            if not model_flags.get(suffix, False):
                continue
            if not present.get(suffix, False):
                # Group not in the file — silent skip (expected case)
                continue
            loaded = _load_model(suffix, h5_path)
            if loaded is None:
                manifest["skipped"].append((suffix, "loader error"))
                continue
            try:
                hdr_extra = _HDR_DISPATCH[suffix](loaded)
                Qm, Im_model, Im_data, dIm = _CURVES_DISPATCH[suffix](loaded, primary)
            except Exception as e:
                manifest["skipped"].append((suffix, f"format error: {e}"))
                continue
            if Qm is None or Im_model is None:
                manifest["skipped"].append((suffix, "missing Q/I_model arrays"))
                continue
            # If I_data or dI couldn't be reconstructed, fill with NaN
            if Im_data is None:
                Im_data = np.full_like(np.asarray(Qm, dtype=float), np.nan)
            if dIm is None:
                dIm = np.full_like(np.asarray(Qm, dtype=float), np.nan)
            if include_header:
                hdr_lines = _format_data_header(
                    h5_path, primary,
                    extra_lines=hdr_extra,
                    columns_label="Q  I_model  I_data  dI",
                    notes=notes,
                )
            else:
                hdr_lines = []
            mp = out_dir / f"{stem}_{suffix}.dat"
            _write_dat(mp, [Qm, Im_model, Im_data, dIm],
                       hdr_lines, delimiter, precision)
            manifest["models"].append(mp)

    return manifest
