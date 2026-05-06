"""
Morphology metrics for binary 3-D voxelgrams.

Computes connectivity / topology / pore-size descriptors that are
standard in the porous-media / micro-CT analysis literature, applied to
the **minority phase** of a two-phase voxelgram (the structure of
interest — pores in a solid matrix, or solid inclusions in a void
matrix, whichever is rarer by volume).

All metrics use only `scipy.ndimage` and (optionally for one item)
`skimage.measure`, both of which are common dependencies.  Compute time
for a 256³ voxelgram is well under one second on a typical laptop.

Reusable from both the SAXS-Morph engine and (in principle) the Fractals
voxel display, though it is currently called only by SAXS-Morph because
that's where the structure analysis adds genuine new information vs the
direct lattice we already know in Fractals.

References
----------
* MicroCT analysis of connectivity in porous structures, R. Soc.
  Interface 17 20190833 (2020) doi:10.1098/rsif.2019.0833
* PoreSpy: A Python Toolkit for Quantitative Analysis of Porous Media
  Images.  J. Open Source Softw. 4(37) 1296 (2019)

Public API
----------
MorphologyMetrics      — dataclass holding all computed scalars
compute_morphology_metrics(voxelgram, voxel_pitch_A) -> MorphologyMetrics
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Optional

import numpy as np
from scipy import ndimage as _ndi


# ---------------------------------------------------------------------------
# Dataclass
# ---------------------------------------------------------------------------

@dataclass
class MorphologyMetrics:
    """Topology / connectivity / pore-size metrics of a binary voxelgram.

    All metrics describe the **minority phase** (the rarer-by-volume
    phase) — the "structure of interest" for users analysing a porous
    sample.

    Note on resolution: every metric depends on the chosen voxel
    pitch.  A coarser grid loses fine pores → fewer clusters, simpler
    topology, smaller measured pore radii.  Always report alongside
    the voxel pitch to allow others to reproduce.

    Note on stochasticity: SAXS-Morph generates one stochastic
    realisation of a structure compatible with the data invariant.
    Different RNG seeds give different topology metrics on the same
    sample; metrics quoted here are *single-realisation* values.  A
    later release may add seed-averaged "expected" metrics.
    """

    # Identity / context
    minority_phase_value: int           # 0 or 1: which phase value was analysed
    minority_volume_fraction: float     # = phi or (1-phi), whichever is smaller
    voxel_pitch_A: float                # echoed for self-contained interpretation

    # ── Connectivity (Tier A) ─────────────────────────────────────────
    n_clusters: int                     # number of disconnected components (6-conn)
    open_porosity_fraction: float       # voxels in clusters touching ANY box face
    closed_porosity_fraction: float     # = 1 − open_porosity_fraction
    percolating_x: bool                 # spans box from −X face to +X face
    percolating_y: bool                 #   ditto Y
    percolating_z: bool                 #   ditto Z

    # ── Topology (Tier B) ─────────────────────────────────────────────
    euler_number: int                   # χ = N_components − N_handles + N_cavities
                                        # χ < 0: spongy with many handles
                                        # χ > 0: many isolated objects

    # ── Pore size from EDT (Tier B) ───────────────────────────────────
    # These are PERCENTILES of the maximum-inscribed-sphere radius at
    # every minority-phase voxel.  Comparable in spirit (not exactly
    # identical to) mercury-intrusion porosimetry results.
    pore_size_median_A: float           # 50th percentile, Å
    pore_size_q25_A: float              # 25th percentile, Å
    pore_size_q75_A: float              # 75th percentile, Å

    @classmethod
    def empty(cls) -> "MorphologyMetrics":
        """Return a sentinel instance with NaN / zero / False values.

        Used as a placeholder when the input voxelgram is empty or
        when a metric cannot be computed.
        """
        return cls(
            minority_phase_value=-1,
            minority_volume_fraction=float("nan"),
            voxel_pitch_A=float("nan"),
            n_clusters=0,
            open_porosity_fraction=float("nan"),
            closed_porosity_fraction=float("nan"),
            percolating_x=False, percolating_y=False, percolating_z=False,
            euler_number=0,
            pore_size_median_A=float("nan"),
            pore_size_q25_A=float("nan"),
            pore_size_q75_A=float("nan"),
        )

    def as_dict(self) -> dict:
        """Return a flat dict of name → primitive value (for HDF5 / JSON)."""
        d = asdict(self)
        # Cast booleans to int for clean HDF5 storage; downstream readers
        # can still treat as truthy.
        for k, v in list(d.items()):
            if isinstance(v, bool):
                d[k] = int(v)
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "MorphologyMetrics":
        """Reconstruct from the dict produced by `as_dict()`."""
        bool_keys = {"percolating_x", "percolating_y", "percolating_z"}
        kwargs = {}
        for k, v in d.items():
            if k in bool_keys:
                kwargs[k] = bool(int(v))
            else:
                kwargs[k] = v
        return cls(**kwargs)


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

def _identify_minority_phase(voxelgram: np.ndarray) -> int:
    """Return which phase value (0 or 1) is the minority by voxel count."""
    n_total = int(voxelgram.size)
    n_one = int(np.sum(voxelgram > 0))
    return 1 if n_one <= n_total - n_one else 0


def _percolates(labels: np.ndarray, axis: int) -> bool:
    """True when at least one labelled cluster touches BOTH opposite faces
    along `axis` (i.e. there is a continuous path through the box)."""
    if axis == 0:
        neg = set(np.unique(labels[0, :, :]).tolist())
        pos = set(np.unique(labels[-1, :, :]).tolist())
    elif axis == 1:
        neg = set(np.unique(labels[:, 0, :]).tolist())
        pos = set(np.unique(labels[:, -1, :]).tolist())
    else:
        neg = set(np.unique(labels[:, :, 0]).tolist())
        pos = set(np.unique(labels[:, :, -1]).tolist())
    neg.discard(0)
    pos.discard(0)
    return bool(neg & pos)


def _open_closed_voxel_counts(labels: np.ndarray) -> tuple[int, int]:
    """Return `(n_open_voxels, n_closed_voxels)` for the labelled image.

    "Open" = belonging to any cluster that touches at least one box face.
    "Closed" = isolated, fully interior cluster.  Background (label 0)
    is excluded.
    """
    face_labels: set[int] = set()
    face_labels.update(np.unique(labels[0, :, :]).tolist())
    face_labels.update(np.unique(labels[-1, :, :]).tolist())
    face_labels.update(np.unique(labels[:, 0, :]).tolist())
    face_labels.update(np.unique(labels[:, -1, :]).tolist())
    face_labels.update(np.unique(labels[:, :, 0]).tolist())
    face_labels.update(np.unique(labels[:, :, -1]).tolist())
    face_labels.discard(0)

    n_total = int(np.sum(labels > 0))
    if not face_labels:
        return 0, n_total

    open_mask = np.isin(labels, list(face_labels))
    n_open = int(np.sum(open_mask))
    n_closed = n_total - n_open
    return n_open, n_closed


# ---------------------------------------------------------------------------
# Public function
# ---------------------------------------------------------------------------

def compute_morphology_metrics(
    voxelgram: np.ndarray,
    voxel_pitch_A: float,
    minority_phase: Optional[int] = None,
) -> MorphologyMetrics:
    """Compute connectivity / topology / pore-size metrics for the
    minority phase of a binary voxelgram.

    Parameters
    ----------
    voxelgram : (N, N, N) array, ideally `uint8` with values in `{0, 1}`.
        Floats are accepted and binarised at `> 0.5`.
    voxel_pitch_A : float
        Edge length of one voxel in Å (used to convert pore-size
        percentiles to physical units).
    minority_phase : int or None
        Force which phase value (0 or 1) to analyse.  None → auto-pick
        whichever has the smaller volume fraction.

    Returns
    -------
    MorphologyMetrics
        Dataclass with the 12 + identity/context fields described in
        `MorphologyMetrics`'s docstring.
    """
    vox_bool = (np.asarray(voxelgram) > 0.5)
    n_total = int(vox_bool.size)
    if n_total == 0:
        return MorphologyMetrics.empty()

    n_one = int(np.sum(vox_bool))
    phi_one = n_one / n_total
    if minority_phase is None:
        minority_phase = 1 if phi_one <= 0.5 else 0

    if minority_phase == 1:
        minority_mask = vox_bool
        minority_phi = phi_one
    else:
        minority_mask = ~vox_bool
        minority_phi = 1.0 - phi_one

    # Empty-phase guard
    if int(np.sum(minority_mask)) == 0:
        m = MorphologyMetrics.empty()
        m.minority_phase_value = int(minority_phase)
        m.minority_volume_fraction = float(minority_phi)
        m.voxel_pitch_A = float(voxel_pitch_A)
        return m

    # ── Connected components (6-connectivity) ─────────────────────────
    # 6-conn (face-only neighbours) is the physically right choice for
    # "can fluid flow through this network" — diagonal-only contacts at
    # corners / edges are not real connections for transport.
    structure = _ndi.generate_binary_structure(3, 1)
    labels, n_clusters = _ndi.label(minority_mask, structure=structure)

    # ── Open vs closed porosity ───────────────────────────────────────
    n_open, n_closed = _open_closed_voxel_counts(labels)
    n_minority = n_open + n_closed
    open_frac = (n_open / n_minority) if n_minority > 0 else 0.0
    closed_frac = 1.0 - open_frac

    # ── Percolation per axis ──────────────────────────────────────────
    perc_x = _percolates(labels, 0)
    perc_y = _percolates(labels, 1)
    perc_z = _percolates(labels, 2)

    # ── Euler-Poincaré characteristic ─────────────────────────────────
    # χ = N_components − N_handles + N_cavities.  Single integer that
    # captures topological complexity: a sphere has χ=2, a torus χ=0,
    # a sponge with many handles is large negative.  Reported per box,
    # not per unit volume.
    try:
        from skimage.measure import euler_number as _euler_number
        euler = int(_euler_number(minority_mask, connectivity=1))
    except Exception:
        # skimage missing or older signature — graceful zero so the rest
        # of the metric set still ships.
        euler = 0

    # ── Pore-size percentiles from Euclidean distance transform ───────
    # For each minority-phase voxel, the EDT value = distance to the
    # nearest non-minority voxel = radius of the maximum inscribed
    # sphere centred at that voxel.  The distribution of these radii is
    # a quick "pore size summary" — biased toward pore CENTRES being
    # over-represented vs PoreSpy's `local_thickness`, but
    # qualitatively the same answer for an order-of-magnitude check.
    try:
        edt = _ndi.distance_transform_edt(minority_mask)
        radii_A = edt[minority_mask] * float(voxel_pitch_A)
        if radii_A.size > 0:
            pore_q25 = float(np.percentile(radii_A, 25))
            pore_med = float(np.percentile(radii_A, 50))
            pore_q75 = float(np.percentile(radii_A, 75))
        else:
            pore_q25 = pore_med = pore_q75 = float("nan")
    except Exception:
        pore_q25 = pore_med = pore_q75 = float("nan")

    return MorphologyMetrics(
        minority_phase_value=int(minority_phase),
        minority_volume_fraction=float(minority_phi),
        voxel_pitch_A=float(voxel_pitch_A),
        n_clusters=int(n_clusters),
        open_porosity_fraction=float(open_frac),
        closed_porosity_fraction=float(closed_frac),
        percolating_x=bool(perc_x),
        percolating_y=bool(perc_y),
        percolating_z=bool(perc_z),
        euler_number=int(euler),
        pore_size_median_A=float(pore_med),
        pore_size_q25_A=float(pore_q25),
        pore_size_q75_A=float(pore_q75),
    )
