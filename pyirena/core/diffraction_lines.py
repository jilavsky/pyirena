"""
diffraction_lines.py — Compute theoretical powder diffraction stick patterns.

Wraps Dans_Diffraction so the rest of the codebase has a single, stable entry
point. Returns peak (Q, I, hkl) — no profile, no peak shapes. Used by the WAXS
Peak Fit GUI to overlay theoretical lines on experimental I(Q) data for
crystallographic phase identification.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np


@dataclass
class DiffractionPattern:
    """Theoretical stick pattern for one CIF / phase."""

    name: str                          # display name (formula or filename)
    cif_path: str                      # source file
    wavelength_a: float                # Å, used to compute pattern
    q: np.ndarray                      # Å⁻¹, peak positions (sorted ascending)
    intensity: np.ndarray              # relative intensity (max normalised to 1)
    hkl: List[Tuple[int, int, int]]    # Miller indices, one per peak


def _formula_from_crystal(xtl) -> str:
    """Best-effort short formula string for the CIF."""
    try:
        return str(xtl.Properties.molname()).strip()
    except Exception:
        try:
            return str(xtl.name).strip()
        except Exception:
            return "phase"


def compute_pattern(
    cif_path: str | Path,
    wavelength_a: float,
    q_min: float = 0.1,
    q_max: float = 10.0,
    intensity_threshold: float = 1e-3,
) -> DiffractionPattern:
    """Compute the powder diffraction stick pattern from a CIF file.

    Parameters
    ----------
    cif_path
        Path to the CIF file.
    wavelength_a
        X-ray wavelength in Ångström. Q = 4π sin(θ)/λ.
    q_min, q_max
        Q-range (Å⁻¹) to retain. Reflections outside are discarded.
    intensity_threshold
        Reflections below this fraction of the maximum intensity are dropped
        (default 0.001 = 0.1%).

    Returns
    -------
    DiffractionPattern
        Sorted by Q ascending. Intensities normalised so max == 1.

    Raises
    ------
    ImportError
        If Dans_Diffraction is not installed.
    ValueError
        If the CIF cannot be parsed or yields no reflections in the range.
    """
    try:
        import Dans_Diffraction as dif
    except ImportError as exc:
        raise ImportError(
            "Dans_Diffraction is required for theoretical diffraction patterns. "
            "Install with: pip install Dans-Diffraction"
        ) from exc

    cif_path = Path(cif_path)
    if not cif_path.is_file():
        raise FileNotFoundError(f"CIF not found: {cif_path}")

    xtl = dif.Crystal(str(cif_path))
    name = _formula_from_crystal(xtl)

    # output=False suppresses Dans_Diffraction's stdout dump of settings
    xtl.Scatter.setup_scatter(
        scattering_type="x-ray",
        wavelength_a=float(wavelength_a),
        output=False,
    )

    # powder() returns (profile_q, profile_I, reflections); we only use the
    # reflections (sticks). The throwaway profile is cheap (~16k floats).
    # NOTE: do not pass pixels=1 — that path returns NaN intensities.
    _, _, refl = xtl.Scatter.powder(units="q")

    refl = np.asarray(refl, dtype=float)
    if refl.size == 0:
        raise ValueError(f"No reflections returned for {cif_path.name}")

    # Columns: h, k, l, q, intensity
    h = refl[:, 0].astype(int)
    k = refl[:, 1].astype(int)
    l = refl[:, 2].astype(int)
    q = refl[:, 3]
    intens = refl[:, 4]

    # Keep only physical (positive Q) reflections inside requested range
    mask = (q >= q_min) & (q <= q_max) & (intens > 0)
    if not np.any(mask):
        raise ValueError(
            f"No reflections for {cif_path.name} in Q=[{q_min:.3f},{q_max:.3f}] Å⁻¹ "
            f"at λ={wavelength_a:.4f} Å"
        )

    h, k, l, q, intens = h[mask], k[mask], l[mask], q[mask], intens[mask]

    imax = float(np.max(intens))
    if imax <= 0:
        raise ValueError(f"All-zero intensities for {cif_path.name}")
    intens = intens / imax

    # Drop tiny reflections after normalisation
    keep = intens >= intensity_threshold
    h, k, l, q, intens = h[keep], k[keep], l[keep], q[keep], intens[keep]

    # Sort by Q
    order = np.argsort(q)
    q = q[order]
    intens = intens[order]
    hkl = [(int(h[i]), int(k[i]), int(l[i])) for i in order]

    return DiffractionPattern(
        name=name,
        cif_path=str(cif_path),
        wavelength_a=float(wavelength_a),
        q=q,
        intensity=intens,
        hkl=hkl,
    )


def hkl_label(hkl: Tuple[int, int, int]) -> str:
    """Format Miller indices for plot labels: negatives shown with overbar fallback."""
    def _one(i: int) -> str:
        return f"{i}" if i >= 0 else f"-{abs(i)}"
    return f"({_one(hkl[0])}{_one(hkl[1])}{_one(hkl[2])})"
