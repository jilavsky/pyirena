"""
Mass-fractal aggregate generator and analyzer.

This module ports the fractal-aggregate algorithms from Irena's
`IR3_3DModels.ipf` and the relevant pieces of `IR3_3DSupportFunctions.ipf`
to Python.  Pure-numpy / scipy implementation; no Qt imports.  Threadsafe
(callable from a QThread).

Public API
----------
GrowthConfig         : input parameters for `grow_aggregate`.
FractalParams        : computed fractal descriptors of an aggregate.
FractalAggregate     : grown aggregate (positions + neighbors + params).
OptimizerConfig      : targets for `optimize_growth`.
grow_aggregate(...)         : random-walk growth on a simple cubic lattice.
compute_fractal_params(...) : path enumeration → df, c, dmin, R, p, s, Rg.
intensity_unified(...)      : closed-form two-level Unified-fit intensity.
intensity_montecarlo(...)   : PDF-based Monte-Carlo intensity (slow).
voxelize(...)               : convert aggregate positions → 3D sphere voxels.
optimize_growth(...)        : bisection over sticking probability.

Algorithm fidelity
------------------
Several details are reproduced verbatim from Irena (deviation will produce
df/c/dmin values that no longer compare against published Irena results):

* Path-walk junctions with > 3 outgoing neighbors take only the **first 3**
  (matches `IR3A_MT_NextPathStep`).
* Sticking-probability piecewise table — see `_sticking_prob_for_chcnt`.
* Neighbor distance thresholds 1.1, 1.05·√2, 1.05·√3 (lenient by design).
* Voxelization: 10× oversampled cubic grid, sphere kernel of 10 voxels,
  threshold 0.5 (matches `IR3T_ConvertToVoxelGram` + `IR3T_CreateSpheresStructure`).
"""

from __future__ import annotations

import math
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime
from typing import Callable, Optional

import numpy as np

from pyirena.core.modeling import UnifiedLevelPopulation, _unified_level_intensity


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class GrowthConfig:
    """Input parameters for `grow_aggregate`."""

    z: int = 250                       # Degree of aggregation (number of particles)
    sticking_prob: float = 75.0        # Base sticking probability, percent (1-100)
    num_test_paths: int = 2500         # Max unique paths per endpoint during evaluation
    rg_primary: float = 10.0           # Rg of the primary sphere [Å]
    allowed_near_dist: int = 3         # 1=axis-only, 2=face diag, 3=body diag
    attraction: str = "Neutral"        # 'Neutral' / 'Attractive' / 'Repulsive' / 'Not allowed'
    seed: int = 0                      # 0 → random; non-zero → reproducible


@dataclass
class FractalParams:
    """Descriptive parameters computed from a grown aggregate."""

    z: int                             # Actual degree of aggregation
    dmin: float                        # Minimum dimension
    c: float                           # Connectivity dimension
    df: float                          # Mass-fractal dimension
    R_dimensionless: float             # Weighted end-to-end distance (lattice units)
    p: float                           # Weighted average path length
    s: float                           # exp(ln(z)/dmin)
    rg_primary: float                  # Primary-particle Rg [Å] (echoed from config)
    rg_aggregate: float                # Aggregate Rg [Å] (Alex McGlasson formula)
    primary_diameter: float            # 2·√(5/3)·Rg_primary [Å]
    true_sticking_prob: float          # 100·z / attempt_value
    num_endpoints: int                 # Particles with neighbor_count < 2
    num_paths_used: int                # Total unique paths enumerated


@dataclass
class FractalAggregate:
    """A grown aggregate plus its computed parameters and (optionally) I(Q)."""

    positions: np.ndarray              # (Z, 3) int32 — lattice coordinates
    neighbor_list: np.ndarray          # (Z, 26) int32 — neighbor indices, -1 padded
    neighbor_count: np.ndarray         # (Z,) uint8  — outgoing neighbor count
    params: FractalParams
    config: GrowthConfig
    attempt_value: int                 # Total random-walk steps taken
    q: Optional[np.ndarray] = None
    i_unified: Optional[np.ndarray] = None
    i_montecarlo: Optional[np.ndarray] = None
    uuid: str = field(default_factory=lambda: uuid.uuid4().hex[:8])
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    label: str = ""                    # Optional UI label (e.g. "opt-iter-3")


@dataclass
class OptimizerConfig:
    """Targets and growth template for `optimize_growth`."""

    target_dmin: float = 2.0
    target_c: float = 1.2
    tolerance: float = 0.05
    max_iter: int = 10
    z: int = 250
    num_test_paths: int = 2500
    rg_primary: float = 10.0
    allowed_near_dist: int = 3
    attraction: str = "Neutral"
    seed: int = 0


# ---------------------------------------------------------------------------
# Lattice geometry helpers
# ---------------------------------------------------------------------------

# 26 cubic-lattice neighbor offsets (all (dx,dy,dz) with dx,dy,dz ∈ {-1,0,1} except (0,0,0))
_ALL_OFFSETS = np.array(
    [(dx, dy, dz)
     for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)
     if not (dx == 0 and dy == 0 and dz == 0)],
    dtype=np.int32,
)


def _neighbor_distance_threshold(allowed_near_dist: int) -> float:
    """Distance below which a particle is considered 'in contact'."""
    if allowed_near_dist == 1:
        return 1.1                       # Axis-only (a = 1)
    if allowed_near_dist == 2:
        return 1.05 * math.sqrt(2.0)     # +face diagonals (a·√2)
    if allowed_near_dist == 3:
        return 1.05 * math.sqrt(3.0)     # +body diagonals (a·√3)
    raise ValueError(f"allowed_near_dist must be 1/2/3, got {allowed_near_dist}")


def _sticking_prob_for_chcnt(base_sp: float, chcnt: int, attraction: str) -> float:
    """Piecewise sticking-probability table from `IR3A_MakeAgg`.

    chcnt is the number of existing aggregate particles within the neighbor
    distance threshold of the moving particle.  The base SP applies for
    chcnt=1; multi-contact attachment is modulated by the `attraction` mode.
    Returns a percentage in [0, 100].
    """
    if chcnt <= 0:
        return 0.0
    if chcnt == 1:
        return float(base_sp)
    a = attraction.lower()
    if a.startswith("neutral"):
        return float(base_sp)
    if a.startswith("attractive"):
        if chcnt == 2:
            return (base_sp + 100.0) * 0.5
        return (base_sp + 300.0) * 0.25
    if a.startswith("repulsive"):
        if chcnt == 2:
            return (base_sp + 10.0) * 0.5
        return (base_sp + 30.0) * 0.25
    if a.startswith("not"):
        return 0.0
    return float(base_sp)


# ---------------------------------------------------------------------------
# Growth: random-walk MC on a simple cubic lattice
# ---------------------------------------------------------------------------

# Internal safety limit; matches Igor's bail-out.
_MAX_ATTEMPTS = 100_000_000


def grow_aggregate(
    config: GrowthConfig,
    progress_cb: Optional[Callable[[int, int], None]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> FractalAggregate:
    """Grow a mass-fractal aggregate by random-walk MC.

    Ports `IR3A_MakeAgg`.  Particles are placed one at a time on a random
    face of an expanding cubic box, undergo a random ±1 step in x/y/z until
    they are within the neighbor threshold of any existing particle, then
    stick with a probability that depends on neighbor count and the chosen
    `attraction` mode.

    Parameters
    ----------
    config : GrowthConfig
    progress_cb : callable(current_z, total_z) or None
        Optional progress callback invoked after each successful stick.
    cancel_check : callable() -> bool or None
        Polled inside the inner walk loop; if True, raises `RuntimeError`.

    Returns
    -------
    FractalAggregate with the computed `params` already populated.
    """
    z = int(config.z)
    if z < 2:
        raise ValueError(f"z must be ≥ 2, got {z}")

    rng = np.random.default_rng(config.seed if config.seed != 0 else None)
    near_thresh = _neighbor_distance_threshold(config.allowed_near_dist)
    near_thresh_sq = near_thresh ** 2

    positions = np.zeros((z, 3), dtype=np.int32)
    neighbor_list = np.full((z, 26), -1, dtype=np.int32)
    neighbor_count = np.zeros(z, dtype=np.uint8)

    # First particle at origin
    aggct = 1

    attempt_value = 0
    cancel_poll_every = 256        # Throttle cancel/check cost

    while aggct < z:
        # Box edge length (must contain all existing particles + 5 voxel margin)
        farthest = int(np.max(np.abs(positions[:aggct]))) if aggct > 0 else 0
        gl = max(2 * farthest + 10, 12)
        half = gl // 2

        # Choose a random launch wall (1..6) — for very small aggregates use corners (7..8)
        if aggct < 6:
            wall = 1 + int(rng.integers(0, 8))   # 1..8
        else:
            wall = 1 + int(rng.integers(0, 6))   # 1..6
        if wall == 1:
            px, py, pz = -half + 1, int(rng.integers(-half, half + 1)), int(rng.integers(-half, half + 1))
        elif wall == 2:
            px, py, pz = +half - 1, int(rng.integers(-half, half + 1)), int(rng.integers(-half, half + 1))
        elif wall == 3:
            px, py, pz = int(rng.integers(-half, half + 1)), -half + 1, int(rng.integers(-half, half + 1))
        elif wall == 4:
            px, py, pz = int(rng.integers(-half, half + 1)), +half - 1, int(rng.integers(-half, half + 1))
        elif wall == 5:
            px, py, pz = int(rng.integers(-half, half + 1)), int(rng.integers(-half, half + 1)), -half + 1
        elif wall == 6:
            px, py, pz = int(rng.integers(-half, half + 1)), int(rng.integers(-half, half + 1)), +half - 1
        else:
            # Corner launch (wall 7 or 8) — used for very small aggregates
            sx, sy, sz = (1 if rng.integers(0, 2) else -1,
                          1 if rng.integers(0, 2) else -1,
                          1 if rng.integers(0, 2) else -1)
            px, py, pz = sx * (half - 1), sy * (half - 1), sz * (half - 1)

        stuck = False
        local_attempts = 0
        existing = positions[:aggct]            # view into already-placed particles

        while not stuck:
            # ── Cancel poll ────────────────────────────────────────────────
            if cancel_check is not None and (local_attempts & (cancel_poll_every - 1)) == 0:
                if cancel_check():
                    raise RuntimeError("Growth cancelled by user.")

            # ── One random step in x, y, or z ──────────────────────────────
            dim = int(rng.integers(0, 3))
            step = 1 if rng.integers(0, 2) else -1
            if dim == 0:
                px += step
            elif dim == 1:
                py += step
            else:
                pz += step
            attempt_value += 1
            local_attempts += 1

            # ── Boundary check: particle that wanders outside box restarts ─
            if abs(px) > half or abs(py) > half or abs(pz) > half:
                # Re-launch on a fresh wall
                if aggct < 6:
                    wall = 1 + int(rng.integers(0, 8))
                else:
                    wall = 1 + int(rng.integers(0, 6))
                if wall == 1:
                    px, py, pz = -half + 1, int(rng.integers(-half, half + 1)), int(rng.integers(-half, half + 1))
                elif wall == 2:
                    px, py, pz = +half - 1, int(rng.integers(-half, half + 1)), int(rng.integers(-half, half + 1))
                elif wall == 3:
                    px, py, pz = int(rng.integers(-half, half + 1)), -half + 1, int(rng.integers(-half, half + 1))
                elif wall == 4:
                    px, py, pz = int(rng.integers(-half, half + 1)), +half - 1, int(rng.integers(-half, half + 1))
                elif wall == 5:
                    px, py, pz = int(rng.integers(-half, half + 1)), int(rng.integers(-half, half + 1)), -half + 1
                elif wall == 6:
                    px, py, pz = int(rng.integers(-half, half + 1)), int(rng.integers(-half, half + 1)), +half - 1
                else:
                    sx, sy, sz = (1 if rng.integers(0, 2) else -1,
                                  1 if rng.integers(0, 2) else -1,
                                  1 if rng.integers(0, 2) else -1)
                    px, py, pz = sx * (half - 1), sy * (half - 1), sz * (half - 1)

            # ── Distance to all existing particles (squared) ───────────────
            dx = existing[:, 0] - px
            dy = existing[:, 1] - py
            dz = existing[:, 2] - pz
            dsq = dx * dx + dy * dy + dz * dz

            # Hard exclusion: never overlap a fully occupied site
            if np.any(dsq == 0):
                # Random walk would have to pass *through* an occupied site;
                # bump in a random direction to avoid being trapped.
                continue

            within = dsq <= near_thresh_sq
            chcnt = int(np.sum(within))
            if chcnt == 0:
                # Safety circuit-break: if we're wandering forever, bail.
                if attempt_value > _MAX_ATTEMPTS:
                    raise RuntimeError(
                        f"Growth aborted: exceeded {_MAX_ATTEMPTS} steps "
                        f"after placing {aggct}/{z} particles."
                    )
                continue

            # ── Sticking decision ──────────────────────────────────────────
            sp_loc = _sticking_prob_for_chcnt(
                config.sticking_prob, chcnt, config.attraction,
            )
            if sp_loc > 0 and rng.uniform(0.0, 100.0) < sp_loc:
                # Stick at (px, py, pz) — record position and update neighbor lists.
                positions[aggct, 0] = px
                positions[aggct, 1] = py
                positions[aggct, 2] = pz

                # Append the bond to both this particle's and each neighbor's list.
                neighbor_indices = np.where(within)[0]
                for ni in neighbor_indices:
                    if neighbor_count[aggct] < 26:
                        neighbor_list[aggct, neighbor_count[aggct]] = ni
                        neighbor_count[aggct] += 1
                    if neighbor_count[ni] < 26:
                        neighbor_list[ni, neighbor_count[ni]] = aggct
                        neighbor_count[ni] += 1

                aggct += 1
                stuck = True
                if progress_cb is not None:
                    progress_cb(aggct, z)
            else:
                # Did not stick — keep walking, but bail if absurdly stuck.
                if attempt_value > _MAX_ATTEMPTS:
                    raise RuntimeError(
                        f"Growth aborted: exceeded {_MAX_ATTEMPTS} steps "
                        f"after placing {aggct}/{z} particles."
                    )

    # Compute parameters before returning
    params = compute_fractal_params(
        positions, neighbor_list, neighbor_count,
        rg_primary=config.rg_primary,
        num_test_paths=config.num_test_paths,
        attempt_value=attempt_value,
        cancel_check=cancel_check,
    )
    return FractalAggregate(
        positions=positions,
        neighbor_list=neighbor_list,
        neighbor_count=neighbor_count,
        params=params,
        config=config,
        attempt_value=attempt_value,
    )


# ---------------------------------------------------------------------------
# Parameter evaluation: path enumeration + endpoint-pair distance statistics
# ---------------------------------------------------------------------------

def _walk_unique_paths(
    start: int,
    neighbor_list: np.ndarray,
    neighbor_count: np.ndarray,
    max_paths: int,
) -> list[list[int]]:
    """Enumerate up to `max_paths` unique paths from `start` to dead-ends.

    Iterative DFS that mirrors `IR3A_MT_NextPathStep`:

      - Visited set prevents loops.
      - At dead-ends (1 outgoing neighbor that is the prior point), the path
        is recorded.
      - At junctions with 2 or 3 outgoing neighbors, all are explored.
      - At junctions with > 3 outgoing neighbors, only the first 3 are
        explored (Igor compatibility — required to reproduce historical
        df/c/dmin values).
    """
    paths: list[list[int]] = []
    # Stack of (current_path_list, prior_point_or_-1)
    stack: list[tuple[list[int], int]] = [([start], -1)]
    while stack and len(paths) < max_paths:
        path, prior = stack.pop()
        current = path[-1]

        # Outgoing neighbors = neighbor_list[current] minus prior
        nc = int(neighbor_count[current])
        if nc == 0:
            # Isolated particle (unlikely once growth completes, but guard anyway)
            paths.append(path)
            continue
        outgoing = []
        visited = set(path)
        for k in range(nc):
            n = int(neighbor_list[current, k])
            if n < 0 or n == prior:
                continue
            if n in visited:
                continue
            outgoing.append(n)

        if not outgoing:
            # Dead end
            paths.append(path)
            continue

        if len(outgoing) > 3:
            # Igor takes only the first 3 — preserve that bias
            outgoing = outgoing[:3]

        if len(outgoing) == 1:
            # Linear segment — extend in place to avoid stack growth
            n = outgoing[0]
            stack.append((path + [n], current))
        else:
            for n in outgoing:
                stack.append((path + [n], current))

    return paths


def compute_fractal_params(
    positions: np.ndarray,
    neighbor_list: np.ndarray,
    neighbor_count: np.ndarray,
    rg_primary: float,
    num_test_paths: int = 2500,
    attempt_value: int = 0,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> FractalParams:
    """Compute df, c, dmin, R, p, s, Rg of an aggregate from its positions.

    Ports `IR3A_EvaluateAggregateUsingMT` + `IR3A_CalculateParametersMT`.
    """
    z = int(positions.shape[0])
    if z < 4:
        # Too small to compute meaningful fractal dimensions
        primary_diameter = 2.0 * math.sqrt(5.0 / 3.0) * float(rg_primary)
        return FractalParams(
            z=z, dmin=float("nan"), c=float("nan"), df=float("nan"),
            R_dimensionless=float("nan"), p=float("nan"), s=float("nan"),
            rg_primary=float(rg_primary), rg_aggregate=float("nan"),
            primary_diameter=primary_diameter,
            true_sticking_prob=(100.0 * z / attempt_value) if attempt_value > 0 else float("nan"),
            num_endpoints=0, num_paths_used=0,
        )

    # ── Endpoints (chains-end particles, neighbor_count < 2) ───────────────
    endpoint_idx = np.where(neighbor_count < 2)[0]
    num_endpoints = int(endpoint_idx.size)
    if num_endpoints < 2:
        # Pathological compact case — return NaNs but supply primary diameter
        primary_diameter = 2.0 * math.sqrt(5.0 / 3.0) * float(rg_primary)
        return FractalParams(
            z=z, dmin=float("nan"), c=float("nan"), df=float("nan"),
            R_dimensionless=float("nan"), p=float("nan"), s=float("nan"),
            rg_primary=float(rg_primary), rg_aggregate=float("nan"),
            primary_diameter=primary_diameter,
            true_sticking_prob=(100.0 * z / attempt_value) if attempt_value > 0 else float("nan"),
            num_endpoints=num_endpoints, num_paths_used=0,
        )

    # ── Path enumeration: parallel across endpoints ────────────────────────
    def _walk_one(start: int) -> list[list[int]]:
        if cancel_check is not None and cancel_check():
            return []
        return _walk_unique_paths(int(start), neighbor_list, neighbor_count,
                                  int(num_test_paths))

    sum_l = 0.0
    sum_l2 = 0.0
    total_paths = 0
    # Use threads — _walk_unique_paths is pure-Python and contention is minimal
    # compared with the CPU cost of building the path lists.
    n_workers = min(8, max(1, num_endpoints))
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        for paths in pool.map(_walk_one, endpoint_idx):
            for p in paths:
                L = len(p)
                if L > 2:
                    sum_l += L
                    sum_l2 += L * L
                    total_paths += 1

    if total_paths == 0 or sum_l == 0:
        # No usable paths — return NaNs
        primary_diameter = 2.0 * math.sqrt(5.0 / 3.0) * float(rg_primary)
        return FractalParams(
            z=z, dmin=float("nan"), c=float("nan"), df=float("nan"),
            R_dimensionless=float("nan"), p=float("nan"), s=float("nan"),
            rg_primary=float(rg_primary), rg_aggregate=float("nan"),
            primary_diameter=primary_diameter,
            true_sticking_prob=(100.0 * z / attempt_value) if attempt_value > 0 else float("nan"),
            num_endpoints=num_endpoints, num_paths_used=0,
        )

    p = sum_l2 / sum_l         # Weighted average path length

    # ── Endpoint-pair distance statistics ──────────────────────────────────
    end_pos = positions[endpoint_idx].astype(np.float64)
    # All pairwise distances
    diff = end_pos[:, None, :] - end_pos[None, :, :]
    d_sq = np.sum(diff * diff, axis=-1)
    iu = np.triu_indices(num_endpoints, k=1)
    d = np.sqrt(d_sq[iu])
    if d.size == 0 or np.sum(d) == 0:
        R = float("nan")
    else:
        R = float(np.sum(d * d) / np.sum(d))   # Σd² / Σd

    # ── Fractal dimensions ─────────────────────────────────────────────────
    if R > 1.0 and p > 1.0 and z > 1:
        df = math.log(z) / math.log(R)
        c = math.log(z) / math.log(p)
    else:
        df = float("nan")
        c = float("nan")
    if math.isfinite(df) and math.isfinite(c) and c != 0:
        dmin = df / c
    else:
        dmin = float("nan")
    if math.isfinite(dmin) and dmin > 0 and z > 1:
        s = math.exp(math.log(z) / dmin)
    else:
        s = float("nan")

    # ── Physical scaling ──────────────────────────────────────────────────
    primary_diameter = 2.0 * math.sqrt(5.0 / 3.0) * float(rg_primary)
    if math.isfinite(c) and math.isfinite(dmin) and math.isfinite(df) \
            and (dmin - df) != 0:
        rg_aggregate = float(rg_primary) * (z ** ((1.0 / c - 1.0) / (dmin - df)))
    else:
        rg_aggregate = float("nan")
    true_sticking_prob = (100.0 * z / attempt_value) if attempt_value > 0 else float("nan")

    return FractalParams(
        z=z,
        dmin=float(dmin),
        c=float(c),
        df=float(df),
        R_dimensionless=float(R),
        p=float(p),
        s=float(s),
        rg_primary=float(rg_primary),
        rg_aggregate=float(rg_aggregate),
        primary_diameter=float(primary_diameter),
        true_sticking_prob=float(true_sticking_prob),
        num_endpoints=int(num_endpoints),
        num_paths_used=int(total_paths),
    )


# ---------------------------------------------------------------------------
# Intensity (fast): two-level Unified-fit closed form
# ---------------------------------------------------------------------------

def intensity_unified(params: FractalParams, q: np.ndarray) -> np.ndarray:
    """Two-level Beaucage Unified-fit intensity from fractal parameters.

    Ports `IR3A_Calculate1DIntensity` (analytical / fast path).

    Level 1 — primary sphere:
        G = 1, Rg = Rg_primary, P = 4, B = 4π / Rg_primary⁴
    Level 2 — aggregate (mass fractal):
        G = z, Rg = Rg_aggregate, P = df,
        B = (G·P / Rg^P) · Γ(P/2),  RgCutoff = Rg_primary
    """
    q = np.asarray(q, dtype=np.float64)
    if not (math.isfinite(params.rg_primary) and params.rg_primary > 0):
        return np.zeros_like(q)

    rg1 = float(params.rg_primary)
    pop1 = UnifiedLevelPopulation(
        G=1.0, Rg=rg1,
        B=4.0 * math.pi / (rg1 ** 4),
        P=4.0,
        RgCO=0.0,
        correlations=False,
    )
    I1 = _unified_level_intensity(q, pop1)

    if (math.isfinite(params.rg_aggregate) and params.rg_aggregate > 0
            and math.isfinite(params.df) and params.df > 0):
        rg2 = float(params.rg_aggregate)
        df = float(params.df)
        z = float(params.z)
        # B for the aggregate level (matches Igor's IR3A_Calculate1DIntensity)
        try:
            B2 = (z * df / (rg2 ** df)) * math.gamma(df * 0.5)
        except (OverflowError, ValueError):
            B2 = 0.0
        pop2 = UnifiedLevelPopulation(
            G=z, Rg=rg2,
            B=B2, P=df,
            RgCO=rg1,                  # Rg-cutoff at primary level
            correlations=False,
        )
        I2 = _unified_level_intensity(q, pop2)
    else:
        I2 = np.zeros_like(q)

    return I1 + I2


# ---------------------------------------------------------------------------
# Voxelization: positions → 3D sphere voxel cube
# ---------------------------------------------------------------------------

def voxelize(
    positions: np.ndarray,
    oversample: int = 10,
    sphere_voxel_radius: int = 10,
) -> tuple[np.ndarray, float]:
    """Convert (N, 3) lattice positions to a binary voxelgram of touching spheres.

    Ports `IR3T_ConvertToVoxelGram` + `IR3T_CreateSpheresStructure`.

    Returns
    -------
    voxelgram : (N, N, N) uint8
    pitch_lattice : float
        Voxel pitch in **lattice units** (= 1 / oversample).  Multiply by
        `primary_diameter` to convert to Å.
    """
    pos = np.asarray(positions, dtype=np.int32)
    if pos.size == 0:
        return np.zeros((1, 1, 1), dtype=np.uint8), 1.0 / oversample

    max_coord = int(np.max(np.abs(pos))) + 1
    box_size = 2 * max_coord * oversample + 4 * sphere_voxel_radius
    # Round up to even to keep FFT happy
    box_size = int(2 * math.ceil(box_size / 2))

    center = box_size // 2
    voxelgram_seeds = np.zeros((box_size, box_size, box_size), dtype=np.float32)

    # Place a single voxel at each particle's oversampled lattice position.
    for x, y, z in pos:
        ix = center + int(x) * oversample
        iy = center + int(y) * oversample
        iz = center + int(z) * oversample
        if 0 <= ix < box_size and 0 <= iy < box_size and 0 <= iz < box_size:
            voxelgram_seeds[ix, iy, iz] = 1.0

    # Build a sharp sphere kernel, then convolve via FFT.
    r = int(sphere_voxel_radius)
    k_size = 2 * r + 1
    kernel = np.zeros((k_size, k_size, k_size), dtype=np.float32)
    grid = np.arange(k_size) - r
    dx, dy, dz = np.meshgrid(grid, grid, grid, indexing="ij")
    dist_sq = dx * dx + dy * dy + dz * dz
    kernel[dist_sq <= r * r] = 1.0

    # FFT convolution (faster than direct for these sizes)
    from scipy.signal import fftconvolve
    convolved = fftconvolve(voxelgram_seeds, kernel, mode="same")
    voxelgram = (convolved > 0.5).astype(np.uint8)
    return voxelgram, 1.0 / float(oversample)


# ---------------------------------------------------------------------------
# Intensity (slow): Monte-Carlo PDF
# ---------------------------------------------------------------------------

def intensity_montecarlo(
    aggregate: FractalAggregate,
    q: np.ndarray,
    oversample: int = 10,
    sphere_voxel_radius: int = 10,
    max_pairs: int = 10_000_000,
    time_budget_s: float = 20.0,
    progress_cb: Optional[Callable[[float], None]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> np.ndarray:
    """Monte-Carlo PDF-based scattering intensity from the 3D voxelized model.

    Ports `IR3A_Model1DIntensity` + `IR3T_CreatePDFIntensity` +
    `IR3T_CalcIntensityPDF`.  Voxelizes the aggregate, samples random pairs
    of solid voxels for distances, builds a PDF, then Fourier-sine-transforms
    it onto the supplied Q grid.

    Parameters
    ----------
    aggregate : FractalAggregate
    q : (Nq,) Q values [Å⁻¹]
    oversample : voxels per lattice-distance unit (10 matches Igor)
    sphere_voxel_radius : sphere kernel radius in voxels (10 matches Igor)
    max_pairs : hard cap on sampled distance pairs
    time_budget_s : soft cap on sampling duration
    progress_cb : callable(percent) for UI feedback
    cancel_check : callable() -> bool

    Returns
    -------
    I_q : (Nq,) intensity, **not** scaled to data — caller may invariant-scale.
    """
    q = np.asarray(q, dtype=np.float64)

    # 1) Build voxelgram (binary uint8)
    voxelgram, pitch_lattice = voxelize(
        aggregate.positions, oversample=oversample,
        sphere_voxel_radius=sphere_voxel_radius,
    )

    # 2) Convert voxel pitch into Å:
    #    1 lattice unit  ≡  primary_diameter Å
    #    1 voxel         ≡  primary_diameter / oversample Å
    pitch_A = aggregate.params.primary_diameter / float(oversample)

    # 3) Coordinates of all solid voxels
    solid = np.argwhere(voxelgram > 0).astype(np.float32)
    n_solid = solid.shape[0]
    if n_solid < 2:
        return np.zeros_like(q)

    # 4) Random-pair sampling for distances
    rng = np.random.default_rng()
    distances: list[np.ndarray] = []
    pairs_drawn = 0
    chunk = 100_000
    t_start = time.perf_counter()
    while pairs_drawn < max_pairs:
        if cancel_check is not None and cancel_check():
            raise RuntimeError("Monte-Carlo intensity cancelled by user.")
        if (time.perf_counter() - t_start) > time_budget_s:
            break

        # Draw `chunk` random pairs of distinct indices
        i = rng.integers(0, n_solid, size=chunk)
        j = rng.integers(0, n_solid, size=chunk)
        same = (i == j)
        if np.any(same):
            i = i[~same]; j = j[~same]
        d = np.linalg.norm(solid[i] - solid[j], axis=1) * pitch_A
        distances.append(d.astype(np.float32))
        pairs_drawn += d.size
        if progress_cb is not None:
            elapsed = time.perf_counter() - t_start
            pct = min(99, int(100 * max(pairs_drawn / max_pairs,
                                         elapsed / time_budget_s)))
            progress_cb(pct)

    if not distances:
        return np.zeros_like(q)

    all_d = np.concatenate(distances)
    if all_d.size == 0:
        return np.zeros_like(q)

    # 5) Histogram → PDF
    d_max = float(np.max(all_d))
    bin_width = max(0.5 * pitch_A, d_max / 1024.0)
    n_bins = max(32, int(math.ceil(d_max / bin_width)))
    hist, edges = np.histogram(all_d, bins=n_bins, range=(0.0, d_max))
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Normalize PDF to unit area
    area = np.trapezoid(hist.astype(np.float64), centers)
    if area > 0:
        pdf = hist.astype(np.float64) / area
    else:
        pdf = hist.astype(np.float64)

    # 6) Glatter-Kratky sine transform: I(Q) = 4π · Σ p(r) · r² · sinc(Qr) Δr
    delta_r = bin_width
    qr = np.outer(q, centers)             # (Nq, n_bins)
    sinc_term = np.where(qr > 0,
                          np.sin(qr) / np.where(qr == 0, 1.0, qr),
                          1.0)
    integrand = pdf * (centers ** 2) * sinc_term
    I_q = 4.0 * math.pi * np.sum(integrand, axis=1) * delta_r
    if progress_cb is not None:
        progress_cb(100)
    return I_q


# ---------------------------------------------------------------------------
# Optimizer: bisection over sticking probability
# ---------------------------------------------------------------------------

def optimize_growth(
    opt_config: OptimizerConfig,
    on_iter_complete: Optional[Callable[[int, FractalAggregate, float], None]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> tuple[FractalAggregate, list[FractalAggregate]]:
    """Bisection search over `sticking_prob` to match target dmin and c.

    Objective = (dmin − target_dmin)² + (c − target_c)².  Iteratively grows
    aggregates at low/mid/high sticking probabilities, narrows the bracket
    around the best, and stops at `max_iter` or when the objective drops
    below `tolerance²`.

    Returns
    -------
    (best_aggregate, all_aggregates_attempted)
    """
    sp_lo = 10.0
    sp_hi = 90.0
    target_dmin = float(opt_config.target_dmin)
    target_c = float(opt_config.target_c)

    all_attempts: list[FractalAggregate] = []
    best: Optional[FractalAggregate] = None
    best_obj = math.inf

    def _grow_and_score(sp: float, iter_label: str) -> tuple[FractalAggregate, float]:
        cfg = GrowthConfig(
            z=opt_config.z,
            sticking_prob=sp,
            num_test_paths=opt_config.num_test_paths,
            rg_primary=opt_config.rg_primary,
            allowed_near_dist=opt_config.allowed_near_dist,
            attraction=opt_config.attraction,
            seed=opt_config.seed,
        )
        agg = grow_aggregate(cfg, cancel_check=cancel_check)
        agg.label = iter_label
        if math.isfinite(agg.params.dmin) and math.isfinite(agg.params.c):
            obj = ((agg.params.dmin - target_dmin) ** 2
                   + (agg.params.c - target_c) ** 2)
        else:
            obj = math.inf
        return agg, obj

    for it in range(int(opt_config.max_iter)):
        if cancel_check is not None and cancel_check():
            break
        sp_mid = 0.5 * (sp_lo + sp_hi)
        candidates = []
        for sp in (sp_lo, sp_mid, sp_hi):
            agg, obj = _grow_and_score(sp, f"opt-iter-{it+1}-sp{sp:.0f}")
            all_attempts.append(agg)
            candidates.append((sp, agg, obj))
            if on_iter_complete is not None:
                on_iter_complete(it + 1, agg, obj)
            if obj < best_obj:
                best_obj = obj
                best = agg

        # Narrow bracket around the best of the three
        candidates.sort(key=lambda t: t[2])
        winner_sp = candidates[0][0]
        if winner_sp == sp_lo:
            sp_hi = sp_mid
        elif winner_sp == sp_hi:
            sp_lo = sp_mid
        else:
            # Centered: tighten symmetrically around mid
            half = 0.25 * (sp_hi - sp_lo)
            sp_lo = sp_mid - half
            sp_hi = sp_mid + half

        if best_obj < (opt_config.tolerance ** 2):
            break

    if best is None:
        # No usable aggregate produced — surface the most recent one
        best = all_attempts[-1] if all_attempts else _grow_and_score(50.0, "opt-fallback")[0]
    return best, all_attempts
