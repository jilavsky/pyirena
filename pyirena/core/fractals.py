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
    rg_aggregate: float                # Aggregate Rg [Å] from particle positions (parallel-axis)
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

    # Aggregate Rg directly from particle positions (parallel-axis theorem):
    #   Rg_total² = Rg_centers² + Rg_sphere²
    # where Rg_sphere = Rg_primary (for a solid sphere of radius R with
    # Rg_sphere = R·√(3/5) = Rg_primary when R = primary_radius).
    # Rg_centers is the rms distance of the Z particle centres from their
    # centroid, in Ångströms.  The lattice unit = primary_diameter, so
    # multiplying lattice-unit Rg by primary_diameter converts to Å.
    #
    # This replaces the earlier Alex McGlasson approximation
    #   Rg_agg = Rg_primary · Z^((1/c - 1)/(dmin - df))
    # which underestimates the true aggregate Rg by 30–100 % depending
    # on Z and df, causing the analytical Unified I(Q) to place its
    # low-Q Guinier knee at too-high Q relative to the MC curve.
    pos_f = positions.astype(np.float64)
    CM = np.mean(pos_f, axis=0)
    Rg_centers_lattice = float(np.sqrt(np.mean(np.sum((pos_f - CM) ** 2, axis=1))))
    Rg_centers_A = Rg_centers_lattice * primary_diameter
    rg_aggregate = math.sqrt(Rg_centers_A ** 2 + float(rg_primary) ** 2)

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
    oversample: int = 20,
    sphere_voxel_radius: int = 10,
) -> tuple[np.ndarray, float]:
    """Convert (N, 3) lattice positions to a binary voxelgram of touching spheres.

    Ports `IR3T_ConvertToVoxelGram` + `IR3T_CreateSpheresStructure`.

    Geometry — the key invariant:

      Physical sphere radius = sphere_voxel_radius · pitch_A
                             = sphere_voxel_radius · (primary_diameter / oversample)

      For this to equal the correct primary radius R = primary_diameter / 2:

                   sphere_voxel_radius / oversample = 1 / 2

      Default (oversample=20, sphere_voxel_radius=10) satisfies the
      invariant AND uses a "fat" 10-voxel kernel that renders smoothly
      and produces a visible neck at edge-neighbor tangent points.
      An equivalent lighter setting is (oversample=10,
      sphere_voxel_radius=5) — same physics, 8× less memory, but the
      5-voxel kernel renders thin tangent connections that VTK's
      flying-edges iso-surface barely shows.

      Lattice spacing in voxels = `oversample`; edge-neighbor centers
      are `oversample` voxels apart.  With sphere_voxel_radius =
      oversample/2, edge neighbors are exactly tangent.  Face- and
      body-diagonal neighbors are at √2 and √3 times the lattice
      spacing and have visible gaps — that is the true geometry for
      physical primary radius R.

    Avoid setting sphere_voxel_radius = oversample (Irena's classic
    choice): the resulting physical sphere radius is the lattice
    spacing = 2R, which inflates the Rg of the solid voxel cloud and
    shifts BOTH MC Guinier knees by a factor of 2.

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
# Intensity (slow): point-cloud Debye sum (Shape2SAS-style)
# ---------------------------------------------------------------------------
#
# Method change vs the previous voxel-MC implementation:
#
# Previously each primary sphere was rasterized into ~thousands of solid
# voxels on a 20× oversampled grid, then RANDOM PAIRS of solid voxels
# were sampled to estimate the pair-distance distribution.  Two flaws
# accumulated:
#   1. Random sampling left only ~10⁷ of the ~10¹⁰ possible pairs as
#      the histogram input → severe shot noise per bin → noisy I(Q),
#      especially at high Q where the Porod amplitude is small relative
#      to the noise floor.
#   2. Stair-stepped voxel surfaces have the wrong fine-scale geometry,
#      destroying the expected Q⁻⁴ Porod tail above the voxel-Nyquist Q.
#
# The new method follows Andreas H. Larsen's Shape2SAS
# (github.com/andreashlarsen/Shape2SAS) exactly:
#   1. Generate ~50 uniformly-distributed points inside each primary
#      sphere (cube-root inverse CDF for r, Gaussian normalize for
#      direction).
#   2. Compute ALL N(N−1)/2 unique pair distances in a triangular loop,
#      histogramming each row's distances on the fly to keep memory
#      bounded.
#   3. Apply the Debye sum  I(Q) = Σ p(r_i) · sinc(Q·r_i) / Σ p(r_i)
#      over all bins, normalized so I(0) = 1.
#
# Outcome: deterministic (no MC noise), smooth Q⁻⁴ Porod tail, faster
# and lower memory than the voxel approach for typical aggregates.

def _uniform_points_in_sphere(
    n: int, radius: float, rng: np.random.Generator,
) -> np.ndarray:
    """Generate `n` uniformly-distributed points inside a sphere of given
    radius centered at the origin.

    Uses inverse-CDF for the radial coordinate (r = R · u^(1/3) gives
    uniform 3-D volume density) and Gaussian normalization for an
    isotropic direction.
    """
    u = rng.uniform(0.0, 1.0, size=n)
    r = radius * np.cbrt(u)
    g = rng.standard_normal((n, 3))
    g /= np.linalg.norm(g, axis=1, keepdims=True)
    return (r[:, None] * g).astype(np.float64)


def _point_cloud_from_aggregate(
    positions: np.ndarray, primary_diameter: float,
    n_points_per_sphere: int, seed: Optional[int],
    polydispersity: float = 0.0,
) -> np.ndarray:
    """Build the (Z·n_points_per_sphere, 3) point cloud in physical Å.

    Each lattice center contributes `n_points_per_sphere` uniformly-
    distributed points inside a sphere of physical radius
    R = primary_diameter / 2, translated to the lattice center's
    physical position (lattice spacing = primary_diameter Å).

    `polydispersity` (relative σ, e.g. 0.10 = 10 % size variation):
    each primary sphere gets its own R drawn from a Gaussian around
    the mean radius.  Without polydispersity all primary spheres are
    identical, and their form-factor zeros align coherently → strong
    fringes in I(Q) above Q ≈ π/(2R).  A small polydispersity
    decoheres those fringes and reveals the smooth Porod envelope.
    """
    rng = np.random.default_rng(seed)
    R_mean = 0.5 * float(primary_diameter)
    Z = positions.shape[0]
    n = int(n_points_per_sphere)
    cloud = np.empty((Z * n, 3), dtype=np.float64)

    if polydispersity > 0:
        # One R per primary sphere, drawn from N(R_mean, polydispersity·R_mean).
        # Clip to avoid pathological tiny / huge spheres.
        sigma = float(polydispersity) * R_mean
        R_per_sphere = rng.normal(R_mean, sigma, size=Z)
        np.clip(R_per_sphere, 0.1 * R_mean, 3.0 * R_mean, out=R_per_sphere)
    else:
        R_per_sphere = np.full(Z, R_mean)

    for i in range(Z):
        center = positions[i].astype(np.float64) * float(primary_diameter)
        cloud[i * n:(i + 1) * n] = center + _uniform_points_in_sphere(
            n, float(R_per_sphere[i]), rng,
        )
    return cloud


def _pair_distance_histogram(
    points: np.ndarray, n_bins: int, r_max: float,
    progress_cb: Optional[Callable[[float], None]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> np.ndarray:
    """Compute the histogram of all N(N−1)/2 unique pair distances.

    Uses a triangular loop (i, then j > i) and accumulates the histogram
    incrementally so the full distance array is never materialized.  For
    N = 12_500 (Z=250 × n=50) this is 78 M distances; the loop runs in
    a few seconds and uses only ~MB of RAM (one row of distances at a
    time).
    """
    N = int(points.shape[0])
    hist = np.zeros(n_bins, dtype=np.int64)
    if N < 2 or r_max <= 0:
        return hist
    inv_bw = float(n_bins) / float(r_max)
    poll_every = 1024
    for i in range(N - 1):
        if cancel_check is not None and (i & (poll_every - 1)) == 0 and cancel_check():
            raise RuntimeError("Monte-Carlo intensity cancelled by user.")
        diff = points[i] - points[i + 1:]
        d = np.sqrt(np.sum(diff * diff, axis=1))
        # Direct bin-index → bincount is faster than np.histogram for
        # known fixed range (avoids np.histogram's edge-search overhead).
        idx = (d * inv_bw).astype(np.int64)
        np.clip(idx, 0, n_bins - 1, out=idx)
        hist += np.bincount(idx, minlength=n_bins)[:n_bins]
        if progress_cb is not None and (i & 255) == 0 and N > 1:
            progress_cb(int(99 * i / (N - 1)))
    return hist


def _debye_sum(q: np.ndarray, r: np.ndarray, p: np.ndarray) -> np.ndarray:
    """I(Q) = (1 / I0) · Σ p(r_i) · sinc(Q·r_i),  sinc(x) = sin(x)/x.

    Normalised so I(Q=0) = 1.
    """
    qr = np.outer(q, r)
    # sin(x)/x with the L'Hôpital limit 1 at x = 0
    sinc_term = np.where(
        qr > 1e-12,
        np.sin(qr) / np.where(qr > 0, qr, 1.0),
        1.0,
    )
    I0 = float(p.sum())
    if I0 <= 0:
        return np.zeros_like(q)
    return (p[None, :] * sinc_term).sum(axis=1) / I0


def intensity_montecarlo(
    aggregate: FractalAggregate,
    q: np.ndarray,
    n_points_per_sphere: int = 50,
    n_bins: Optional[int] = None,
    polydispersity: float = 0.10,
    seed: Optional[int] = None,
    progress_cb: Optional[Callable[[float], None]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> np.ndarray:
    """Debye-sum scattering intensity via a uniform point-cloud
    representation of the aggregate (Shape2SAS-style).

    Algorithm:
      1. Place `n_points_per_sphere` uniformly-distributed points inside
         every primary sphere of radius R = primary_diameter / 2,
         translated to its lattice position (lattice spacing = D).
         If `polydispersity > 0` each primary sphere's R is drawn from
         a Gaussian around the mean — see below.
      2. Compute ALL N(N−1)/2 pair distances deterministically; histogram
         them into linear r-bins.  By default `n_bins` is chosen so that
         the bin width is small compared to the smallest sinc-oscillation
         period at the user's `Q_max` (see notes below) — required when
         the aggregate spans thousands of Ångströms.
      3. Apply the Debye sum  I(Q) = Σ p(r) · sinc(Q·r) / Σ p(r)
         over every bin, normalised so I(0) = 1.

    Bin-width / Q_max coupling
    --------------------------
    Each histogram bin contributes `p_bin · sinc(Q · r_centre)` to I(Q).
    This is accurate only while sinc varies little within the bin —
    i.e. while `Q · bin_width ≪ 1`.  When `Q · bin_width > 1`, sinc
    oscillates *within* one bin and the centre-evaluation gives a
    pseudo-random value of order ±1/(Q·r) instead of the smooth average.
    Result: spurious fluctuations and even negative I(Q) at high Q.

    For tiny primaries (Rg ≈ 10 Å) the aggregate is small (~hundreds of
    Å) so 200 bins gives bin_width ~ 1 Å — fine up to Q ≈ 1.  For large
    primaries (Rg ≈ 200 Å) the aggregate grows to ~10⁴ Å so 200 bins
    gives bin_width ~ 70 Å and the breakdown starts at Q ≈ 0.01.

    Auto-mode (default): `n_bins` is chosen so `bin_width ≤ 0.1 / Q_max`,
    which keeps Q·bin_width ≤ 0.1 throughout the requested Q range.
    Capped at 200 000 bins for safety (still O(MB) memory and O(few %)
    extra Debye-sum cost vs the histogram itself).  Pass an explicit
    `n_bins` to override (e.g. for testing).

    Parameters
    ----------
    aggregate : FractalAggregate
    q : (Nq,) Q values [Å⁻¹]
    n_points_per_sphere : int (default 50)
        Number of uniform points per primary sphere.  Higher → finer
        high-Q resolution at O(N²) cost.
    n_bins : int or None (default None → auto)
        Number of r-bins in the pair-distance histogram.  None →
        adaptive (`r_max · Q_max · 10`, clipped to [200, 200 000]).
        Override only for testing or to force a coarser histogram.
    polydispersity : float (default 0.10)
        Relative size variation of the primary spheres (σ_R / R_mean).
        Without it (= 0) all primaries are identical and their form-
        factor zeros align coherently, producing strong fringes
        starting near Q ≈ 4.5/R.  A small polydispersity (≈ 0.05–0.20)
        decoheres those fringes so the smooth Porod Q⁻⁴ envelope
        becomes visible — exactly as in real polydisperse SAXS samples.
        Set 0 to see the strict monodisperse-aggregate fringes.
    seed : int or None
        RNG seed for reproducible point clouds (None → fresh random).
    progress_cb : callable(percent) for UI feedback.
    cancel_check : callable() → bool.

    Returns
    -------
    I_q : (Nq,) intensity, normalised so I(0) = 1.  No NaN truncation —
        the Debye sum is mathematically valid at all Q; deviations from
        the smooth Porod envelope above
        Q ≈ π / mean_inter_point_distance reflect the discreteness of
        the point cloud.  Use a larger `n_points_per_sphere` if you
        need a wider clean-Porod range.
    """
    q = np.asarray(q, dtype=np.float64)

    # 1) Uniform point cloud (with optional per-sphere polydispersity)
    points = _point_cloud_from_aggregate(
        aggregate.positions,
        aggregate.params.primary_diameter,
        n_points_per_sphere=n_points_per_sphere,
        seed=seed,
        polydispersity=float(polydispersity),
    )
    if points.shape[0] < 2:
        return np.zeros_like(q)

    # 2) Determine r_max from bounding-box diagonal
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    r_max = float(np.sqrt(np.sum((maxs - mins) ** 2))) * 1.05
    if not (r_max > 0):
        return np.zeros_like(q)

    # 3) Pick n_bins adaptively if not specified.  Goal: bin_width small
    #    enough that Q_max · bin_width ≤ 0.1, so the sinc-at-bin-centre
    #    approximation stays accurate over the full requested Q range.
    #    See "Bin-width / Q_max coupling" in the docstring.
    if n_bins is None:
        q_pos = q[q > 0]
        Q_max = float(np.max(q_pos)) if q_pos.size > 0 else 1.0
        # bin_width target = 0.1 / Q_max  →  n_bins = r_max / bin_width = 10 · r_max · Q_max
        n_bins_auto = int(math.ceil(10.0 * r_max * Q_max))
        n_bins = max(200, min(n_bins_auto, 200_000))
    else:
        n_bins = int(n_bins)

    # 4) Histogram all pair distances (deterministic, on-the-fly)
    hist = _pair_distance_histogram(
        points, n_bins=n_bins, r_max=r_max,
        progress_cb=progress_cb, cancel_check=cancel_check,
    )
    if hist.sum() == 0:
        return np.zeros_like(q)

    edges = np.linspace(0.0, r_max, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # 5) Debye sum
    I_q = _debye_sum(q, centers, hist.astype(np.float64))

    if progress_cb is not None:
        progress_cb(100)
    return I_q


def mc_q_max(aggregate: FractalAggregate, n_points_per_sphere: int = 50) -> float:
    """Estimated upper Q above which the discrete-point-cloud Debye sum
    starts deviating from true continuous-particle scattering.

    Estimate: Q_max ≈ π / d_typical, where d_typical = 2R / n^(1/3) is
    the mean inter-point distance for `n` uniform points in a sphere of
    radius R = primary_diameter / 2.  For the default n=50, R≈12.9 Å this
    gives Q_max ≈ 0.45 Å⁻¹.  Increase `n_points_per_sphere` to push it
    higher.

    No longer used as a hard truncation by `intensity_montecarlo`
    (returned for informational use only — callers may display a
    vertical guide at this Q).
    """
    R = 0.5 * float(aggregate.params.primary_diameter)
    n = max(int(n_points_per_sphere), 1)
    d_typical = 2.0 * R / (float(n) ** (1.0 / 3.0))
    if d_typical <= 0:
        return float("inf")
    return math.pi / d_typical


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
