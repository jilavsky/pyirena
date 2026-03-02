"""
data_merge.py — Core engine for merging two SAS datasets.

No GUI dependencies. Handles optimization (scale, Q-shift, background) in the
overlap region using log-log interpolation, then concatenates the two datasets.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
from scipy.optimize import minimize


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class MergeConfig:
    """All optimization configuration for a single merge operation.

    Parameters
    ----------
    q_overlap_min, q_overlap_max : float
        Linear Q bounds (Å⁻¹) of the overlap/optimization region.
        Set by the user via cursors.
    fit_scale : bool
        Whether to optimize a multiplicative scale factor.
    scale_dataset : int
        1 = scale DS1, 2 = scale DS2.  DS1 is assumed absolute-scale, so
        default is 2 (scale DS2 to match DS1).
    fit_qshift : bool
        Whether to optimize an additive Q shift.
    qshift_dataset : int
        0 = no shift, 1 = shift DS1, 2 = shift DS2.
    method : str
        Key into DataMerge.METHODS.  Currently only 'interpolation'.
    split_at_left_cursor : bool
        If True, final merged array takes DS1 only for Q < q_overlap_min and
        DS2 for Q >= q_overlap_min (no overlap region in output).
        Default False: all data from both datasets is included (interleaved in
        the overlap region).
    """
    q_overlap_min: float
    q_overlap_max: float
    fit_scale: bool = True
    scale_dataset: int = 2
    fit_qshift: bool = False
    qshift_dataset: int = 0
    method: str = 'interpolation'
    split_at_left_cursor: bool = False


@dataclass
class MergeResult:
    """Output of DataMerge.optimize().

    Parameters
    ----------
    scale : float
        Optimal scale factor applied to the selected dataset (1.0 if fit_scale
        is False).
    q_shift : float
        Optimal additive Q shift in Å⁻¹ (0.0 if fit_qshift is False).
    background : float
        Optimal constant background subtracted from DS1 before comparison.
    success : bool
        Whether the optimizer converged.
    chi_squared : float
        Sum of squared weighted residuals at the optimum.
    n_overlap_points : int
        Number of DS2 overlap-region points used in the objective.
    message : str
        Optimizer status message.
    """
    scale: float = 1.0
    q_shift: float = 0.0
    background: float = 0.0
    success: bool = False
    chi_squared: float = float('nan')
    n_overlap_points: int = 0
    message: str = ''


# ---------------------------------------------------------------------------
# Main engine
# ---------------------------------------------------------------------------

class DataMerge:
    """Core merge engine — no GUI dependencies."""

    #: Registry of available merge methods.  Values are human-readable labels.
    METHODS: dict[str, str] = {
        'interpolation': 'Log-log linear interpolation',
    }

    # ------------------------------------------------------------------ #
    #  Public API                                                          #
    # ------------------------------------------------------------------ #

    def optimize(
        self,
        q1: np.ndarray, I1: np.ndarray, dI1: np.ndarray,
        q2: np.ndarray, I2: np.ndarray, dI2: np.ndarray,
        config: MergeConfig,
    ) -> MergeResult:
        """Find optimal scale, Q-shift, and background in the overlap region.

        All arrays must be in physical (linear) units.  Negative or zero
        values in the overlap region are masked out automatically.

        Parameters
        ----------
        q1, I1, dI1 : arrays
            Dataset 1 (lower-Q, assumed absolute scale).
        q2, I2, dI2 : arrays
            Dataset 2 (higher-Q, to be brought onto DS1 scale).
        config : MergeConfig
            Optimization settings.

        Returns
        -------
        MergeResult
            Populated with the optimizer output.
        """
        # --- initial guess for scale ---
        mask1_init = (
            (q1 >= config.q_overlap_min) & (q1 <= config.q_overlap_max)
            & (I1 > 0) & np.isfinite(I1)
        )
        mask2_init = (
            (q2 >= config.q_overlap_min) & (q2 <= config.q_overlap_max)
            & (I2 > 0) & np.isfinite(I2)
        )
        if mask1_init.sum() < 2 or mask2_init.sum() < 1:
            return MergeResult(
                success=False,
                message=(
                    f"Insufficient overlap-region data: {mask1_init.sum()} DS1 points, "
                    f"{mask2_init.sum()} DS2 points in [{config.q_overlap_min:.4g}, "
                    f"{config.q_overlap_max:.4g}]"
                ),
            )

        med1 = float(np.median(I1[mask1_init]))
        med2 = float(np.median(I2[mask2_init]))
        if config.fit_scale and med2 > 0:
            scale_init = (med1 / med2) if config.scale_dataset == 2 else (med2 / med1)
        else:
            scale_init = 1.0

        # --- build free-parameter vector and bounds ---
        # Always 3 slots: [background, scale, q_shift]
        # Unused slots are fixed at their initial values via a closure.
        p0 = [0.0, scale_init, 0.0]

        max_bg = float(np.max(np.abs(I1[mask1_init]))) * 0.5
        bounds = [
            (-max_bg, max_bg),  # background
            (0.01, 100.0),      # scale
            (-0.1, 0.1),        # q_shift
        ]

        # Freeze slots that are not fitted
        fixed_bg = not True            # background is always fit
        fixed_scale = not config.fit_scale
        fixed_qshift = not config.fit_qshift

        def _objective_wrapper(free_params: np.ndarray) -> float:
            params = [p0[0], p0[1], p0[2]]
            idx = 0
            if not fixed_bg:
                params[0] = free_params[idx]; idx += 1
            if not fixed_scale:
                params[1] = free_params[idx]; idx += 1
            if not fixed_qshift:
                params[2] = free_params[idx]; idx += 1
            return self._objective(params, q1, I1, dI1, q2, I2, dI2, config)

        # Build initial guess and bounds for the free subset
        free_p0 = []
        free_bounds = []
        if not fixed_bg:
            free_p0.append(p0[0]); free_bounds.append(bounds[0])
        if not fixed_scale:
            free_p0.append(p0[1]); free_bounds.append(bounds[1])
        if not fixed_qshift:
            free_p0.append(p0[2]); free_bounds.append(bounds[2])

        if not free_p0:
            # Nothing to optimise — just evaluate at defaults
            chi2 = self._objective(p0, q1, I1, dI1, q2, I2, dI2, config)
            n_pts = int(mask2_init.sum())
            return MergeResult(
                scale=1.0, q_shift=0.0, background=0.0,
                success=True, chi_squared=chi2,
                n_overlap_points=n_pts,
                message="No free parameters — nothing to optimise.",
            )

        opt = minimize(
            _objective_wrapper,
            x0=np.array(free_p0),
            method='Nelder-Mead',
            bounds=free_bounds,
            options={'xatol': 1e-6, 'fatol': 1e-8, 'maxiter': 5000},
        )

        # Reconstruct full params vector from result
        params_final = [p0[0], p0[1], p0[2]]
        idx = 0
        if not fixed_bg:
            params_final[0] = opt.x[idx]; idx += 1
        if not fixed_scale:
            params_final[1] = opt.x[idx]; idx += 1
        if not fixed_qshift:
            params_final[2] = opt.x[idx]; idx += 1

        # Count overlap points in optimised coordinates
        q1_adj, I1_adj, q2_adj, I2_adj = self._apply_params(
            q1, I1, q2, I2,
            background=params_final[0],
            scale=params_final[1],
            q_shift=params_final[2],
            config=config,
        )
        mask2_final = (
            (q2_adj >= config.q_overlap_min) & (q2_adj <= config.q_overlap_max)
            & (I2_adj > 0) & np.isfinite(I2_adj)
        )

        return MergeResult(
            scale=params_final[1],
            q_shift=params_final[2],
            background=params_final[0],
            success=opt.success or opt.fun < 1e10,
            chi_squared=float(opt.fun),
            n_overlap_points=int(mask2_final.sum()),
            message=opt.message,
        )

    def merge(
        self,
        q1: np.ndarray, I1: np.ndarray, dI1: np.ndarray,
        dQ1: Optional[np.ndarray],
        q2: np.ndarray, I2: np.ndarray, dI2: np.ndarray,
        dQ2: Optional[np.ndarray],
        result: MergeResult,
        config: MergeConfig,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]:
        """Apply optimised parameters and produce merged (sorted) arrays.

        Returns
        -------
        (q_merged, I_merged, dI_merged, dQ_merged)
            dQ_merged is None if both dQ1 and dQ2 are None.
            All arrays are sorted ascending by Q.
        """
        q1_adj, I1_adj, q2_adj, I2_adj = self._apply_params(
            q1, I1, q2, I2,
            background=result.background,
            scale=result.scale,
            q_shift=result.q_shift,
            config=config,
        )

        # Propagate dI for the scaled dataset
        dI1_adj = dI1.copy()
        dI2_adj = dI2.copy()
        if config.fit_scale and config.scale_dataset == 1:
            dI1_adj = dI1 * result.scale
        elif config.fit_scale and config.scale_dataset == 2:
            dI2_adj = dI2 * result.scale

        if config.split_at_left_cursor:
            # DS1 only below overlap start, DS2 from overlap start onward
            mask1 = q1_adj < config.q_overlap_min
            mask2 = q2_adj >= config.q_overlap_min
        else:
            # Include all data from both (interleaved in overlap region)
            mask1 = np.ones(len(q1_adj), dtype=bool)
            mask2 = np.ones(len(q2_adj), dtype=bool)

        q_out = np.concatenate([q1_adj[mask1], q2_adj[mask2]])
        I_out = np.concatenate([I1_adj[mask1], I2_adj[mask2]])
        dI_out = np.concatenate([dI1_adj[mask1], dI2_adj[mask2]])

        if dQ1 is not None or dQ2 is not None:
            dQ1_eff = dQ1 if dQ1 is not None else np.zeros(len(q1))
            dQ2_eff = dQ2 if dQ2 is not None else np.zeros(len(q2))
            dQ_out = np.concatenate([dQ1_eff[mask1], dQ2_eff[mask2]])
        else:
            dQ_out = None

        order = np.argsort(q_out)
        if dQ_out is not None:
            return q_out[order], I_out[order], dI_out[order], dQ_out[order]
        return q_out[order], I_out[order], dI_out[order], None

    def match_files(
        self,
        filenames1: List[str],
        filenames2: List[str],
    ) -> List[Tuple[str, str]]:
        """Match files by (prefix-before-first-underscore, last-integer-in-stem).

        Files that share both parts of the key are considered a pair.
        Unmatched files are silently skipped.  If multiple files in one folder
        share the same key, a warning is printed and the first is used.

        Parameters
        ----------
        filenames1, filenames2 : list of str
            Base filenames (not full paths) from each folder.

        Returns
        -------
        list of (name1, name2) tuples sorted by the common key.
        """
        def _extract_key(fname: str) -> Optional[Tuple[str, str]]:
            stem = Path(fname).stem
            prefix = stem.split('_')[0]
            numbers = re.findall(r'\d+', stem)
            if not numbers:
                return None
            return (prefix, numbers[-1])

        dict1: dict[Tuple[str, str], str] = {}
        for f in filenames1:
            k = _extract_key(f)
            if k is None:
                continue
            if k in dict1:
                print(f"[data_merge] Warning: duplicate key {k!r} in folder1: "
                      f"{dict1[k]!r} and {f!r}. Using first.")
            else:
                dict1[k] = f

        dict2: dict[Tuple[str, str], str] = {}
        for f in filenames2:
            k = _extract_key(f)
            if k is None:
                continue
            if k in dict2:
                print(f"[data_merge] Warning: duplicate key {k!r} in folder2: "
                      f"{dict2[k]!r} and {f!r}. Using first.")
            else:
                dict2[k] = f

        common_keys = sorted(set(dict1.keys()) & set(dict2.keys()))
        return [(dict1[k], dict2[k]) for k in common_keys]

    # ------------------------------------------------------------------ #
    #  Private helpers                                                     #
    # ------------------------------------------------------------------ #

    def _objective(
        self,
        params: list,
        q1: np.ndarray, I1: np.ndarray, dI1: np.ndarray,
        q2: np.ndarray, I2: np.ndarray, dI2: np.ndarray,
        config: MergeConfig,
    ) -> float:
        """Weighted chi-squared in the overlap region.

        params = [background, scale, q_shift]
        """
        background, scale, q_shift = params

        q1_adj, I1_adj, q2_adj, I2_adj = self._apply_params(
            q1, I1, q2, I2, background, scale, q_shift, config
        )

        # DS1 points in the overlap (source for interpolation)
        mask1 = (
            (q1_adj >= config.q_overlap_min) & (q1_adj <= config.q_overlap_max)
            & (I1_adj > 0) & np.isfinite(I1_adj) & np.isfinite(q1_adj)
        )

        if mask1.sum() < 2:
            return 1e30

        # DS2 points in the overlap, clipped to DS1's q range to avoid extrapolation
        q1_ov_min = float(q1_adj[mask1].min())
        q1_ov_max = float(q1_adj[mask1].max())
        mask2 = (
            (q2_adj >= config.q_overlap_min) & (q2_adj <= config.q_overlap_max)
            & (q2_adj >= q1_ov_min) & (q2_adj <= q1_ov_max)
            & (I2_adj > 0) & np.isfinite(I2_adj) & np.isfinite(q2_adj)
        )

        if mask2.sum() < 1:
            return 1e30

        try:
            I1_interp = self._interp_loglog(
                q2_adj[mask2], q1_adj[mask1], I1_adj[mask1]
            )
        except Exception:
            return 1e30

        # Use DS2 uncertainty for weighting; fall back to 5% fractional if zero/missing
        dI2_ov = dI2[mask2]
        weights = np.where(dI2_ov > 0, dI2_ov, I2_adj[mask2] * 0.05)
        weights = np.where(np.isfinite(weights) & (weights > 0), weights, I2_adj[mask2] * 0.05)

        residuals = (I1_interp - I2_adj[mask2]) / weights
        return float(np.sum(residuals ** 2))

    def _apply_params(
        self,
        q1: np.ndarray, I1: np.ndarray,
        q2: np.ndarray, I2: np.ndarray,
        background: float,
        scale: float,
        q_shift: float,
        config: MergeConfig,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return (q1_adj, I1_adj, q2_adj, I2_adj) after all corrections."""
        # Background always subtracted from DS1
        I1_adj = I1 - background

        # Scale applied to the chosen dataset
        if config.fit_scale:
            if config.scale_dataset == 1:
                I1_adj = I1_adj * scale
            else:
                I2 = I2 * scale

        # Q shift applied to the chosen dataset
        if config.fit_qshift and config.qshift_dataset == 1:
            q1 = q1 + q_shift
        elif config.fit_qshift and config.qshift_dataset == 2:
            q2 = q2 + q_shift

        return q1, I1_adj, q2, I2

    @staticmethod
    def _interp_loglog(
        q_target: np.ndarray,
        q_known: np.ndarray,
        I_known: np.ndarray,
    ) -> np.ndarray:
        """Interpolate I at q_target using log-log linear interpolation.

        q_known and I_known must be positive and finite (caller is responsible).
        Out-of-range q_target values return NaN (no extrapolation).

        Returns I values in linear units.
        """
        log_I = np.interp(
            np.log10(q_target),
            np.log10(q_known),
            np.log10(I_known),
            left=float('nan'),
            right=float('nan'),
        )
        return 10.0 ** log_I
