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
    fixed_scale_value : float
        Scale factor applied when *fit_scale* is False.  Set by the user;
        allows manual scaling without optimization.  Default 1.0 (no scaling).
    fit_qshift : bool
        Whether to optimize an additive Q shift.
    fixed_qshift_value : float
        Q shift (Å⁻¹) applied when *fit_qshift* is False.  Allows the user to
        apply a known Q offset without running optimization.  Default 0.0.
    qshift_dataset : int
        0 = no shift, 1 = shift DS1, 2 = shift DS2.
    method : str
        Key into DataMerge.METHODS.  Currently only 'interpolation'.
    split_at_left_cursor : bool
        If True, final merged array takes DS1 only for Q < q_overlap_min and
        DS2 for Q >= q_overlap_min (hard split at left cursor, no overlap in output).
        If False (default), DS1 is trimmed at the right cursor and DS2 at the left
        cursor; the overlap region [q_overlap_min, q_overlap_max] contains
        interleaved points from both datasets.
    """
    q_overlap_min: float
    q_overlap_max: float
    fit_scale: bool = True
    scale_dataset: int = 2
    fixed_scale_value: float = 1.0
    fit_qshift: bool = False
    fixed_qshift_value: float = 0.0
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

        # --- initial guess for scale and background ---
        med1 = float(np.median(I1[mask1_init]))
        med2 = float(np.median(I2[mask2_init]))

        # Scale initial guess: simple median ratio.
        # The log-log-interpolation regression used previously was UNRELIABLE
        # when a Q-shift exists: interpolating DS1 at DS2's Q positions without
        # shift correction misaligns diffraction peaks, causing the regression to
        # return the wrong slope (matching the degenerate BG≈I1, scale≈0 solution).
        # The median ratio is not sensitive to peak positions and always gives a
        # reasonable starting point.
        if config.fit_scale:
            if config.scale_dataset == 2 and med2 > 0:
                scale_init = float(np.clip(med1 / med2, 0.01, 100.0))
            elif config.scale_dataset == 1 and med1 > 0:
                scale_init = float(np.clip(med2 / med1, 0.01, 100.0))
            else:
                scale_init = 1.0
        else:
            scale_init = config.fixed_scale_value

        # Background always starts at 0 — the BG regularisation in
        # _objective_wrapper keeps it small and avoids the scale–BG valley.
        bg_init = 0.0

        # --- build free-parameter vector and bounds ---
        # Always 3 slots: [background, scale, q_shift]
        # Unused slots are fixed at their initial values via a closure.
        p0 = [bg_init, scale_init, 0.0 if config.fit_qshift else config.fixed_qshift_value]

        # Background bound: allow up to the full max of DS1 in overlap.
        # Using 50% of max (the old value) prevents finding backgrounds that
        # exceed half the DS1 signal — a common situation (e.g. instrument
        # background ≈ 80% of DS1 in the tail).
        max_bg = float(np.max(I1[mask1_init]))
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
            # Hard-clip to bounds — Nelder-Mead does not enforce bounds reliably;
            # returning 1e30 for out-of-range values acts as a hard wall.
            if params[0] < bounds[0][0] or params[0] > bounds[0][1]:
                return 1e30
            if params[1] < bounds[1][0] or params[1] > bounds[1][1]:
                return 1e30
            if params[2] < bounds[2][0] or params[2] > bounds[2][1]:
                return 1e30
            chi2 = self._objective(params, q1, I1, dI1, q2, I2, dI2, config)
            # BG regularisation to break the scale–BG valley degeneracy.
            #
            # For scale_dataset=2 the chi-squared is exactly zero along the
            # line (BG = I1_interp − I2*scale), so BOTH (BG≈0, scale≈1) and
            # (BG≈I1, scale≈0) are valid chi-squared minima.  Nelder-Mead finds
            # whichever it hits first, often the wrong one.
            #
            # Adding (BG / med1)² costs ~0 at BG=0 and ~1 chi²-unit at BG=I1_median.
            # This breaks the tie in favour of smaller background without
            # preventing the optimizer from finding genuinely large backgrounds
            # (where the data chi² improvement >> 1 chi²-unit).
            if med1 > 0:
                chi2 += (params[0] / med1) ** 2
            return chi2

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

        # Propagate dI for the scaled dataset (scale always applied, whether
        # fitted or user-specified via fixed_scale_value)
        dI1_adj = dI1.copy()
        dI2_adj = dI2.copy()
        if config.scale_dataset == 1:
            dI1_adj = dI1 * result.scale
        else:
            dI2_adj = dI2 * result.scale

        if config.split_at_left_cursor:
            # Hard split at left cursor: DS1 only below, DS2 from left cursor onward
            mask1 = q1_adj < config.q_overlap_min
            mask2 = q2_adj >= config.q_overlap_min
        else:
            # DS1 trimmed at right cursor, DS2 trimmed at left cursor.
            # The overlap region [left, right] is covered by both (interleaved).
            mask1 = q1_adj <= config.q_overlap_max
            mask2 = q2_adj >= config.q_overlap_min

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

        # Use DS2 uncertainty for weighting, but enforce a minimum relative weight
        # of 5% of I2_adj.  Without this floor, the objective has a degenerate
        # minimum at (BG ≈ I1, scale ≈ 0): dI2 is fixed (not scaled with I2_adj),
        # so residuals (I1_interp - 0) / dI2 are bounded while scale → lower-limit.
        # The relative floor ensures scale=0 is penalised (weights → 0 → division
        # by zero → inf), steering the optimizer away from the trivial solution.
        dI2_ov = dI2[mask2]
        rel_floor = I2_adj[mask2] * 0.05
        abs_weights = np.where(
            (dI2_ov > 0) & np.isfinite(dI2_ov), dI2_ov, rel_floor
        )
        weights = np.maximum(abs_weights, rel_floor)
        weights = np.where(np.isfinite(weights) & (weights > 0), weights, rel_floor)

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

        # Scale always applied to the chosen dataset, whether the scale was
        # fitted (fit_scale=True) or fixed by the user (fit_scale=False).
        # When fit_scale=False the optimizer passes fixed_scale_value as the
        # scale argument, so the result is the same code path.
        if config.scale_dataset == 1:
            I1_adj = I1_adj * scale
        else:
            I2 = I2 * scale

        # Q shift applied to the chosen dataset (whether fitted or user-specified)
        if config.qshift_dataset == 1:
            q1 = q1 + q_shift
        elif config.qshift_dataset == 2:
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
