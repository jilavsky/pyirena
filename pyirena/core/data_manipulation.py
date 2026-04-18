"""Data Manipulation engine — headless operations on SAS datasets.

Provides stateless methods for scaling, trimming, rebinning, averaging,
subtracting, and dividing I(Q) datasets.  All methods work on raw numpy
arrays with no GUI dependencies.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np


# ======================================================================
# Configuration dataclasses
# ======================================================================

@dataclass
class ScaleConfig:
    """Parameters for the Scale + Background operation."""
    scale_I: float = 1.0
    background: float = 0.0              # subtracted *after* scaling
    scale_uncertainty: Optional[float] = None  # None → use scale_I for dI

@dataclass
class TrimConfig:
    """Parameters for Q-range trimming."""
    q_min: float = 0.0
    q_max: float = float('inf')

@dataclass
class RebinConfig:
    """Parameters for rebinning to a new Q grid."""
    mode: str = 'log'                    # 'log' | 'linear' | 'reference'
    n_points: int = 200
    q_min: Optional[float] = None        # None → use data min
    q_max: Optional[float] = None        # None → use data max
    reference_q: Optional[np.ndarray] = None  # for mode='reference'

@dataclass
class SubtractConfig:
    """Parameters for buffer/solvent subtraction."""
    buffer_scale: float = 1.0
    auto_scale: bool = False
    auto_q_min: Optional[float] = None
    auto_q_max: Optional[float] = None

@dataclass
class DivideConfig:
    """Parameters for dataset division (structure factor extraction)."""
    denominator_scale: float = 1.0
    denominator_background: float = 0.0


# ======================================================================
# Result dataclass
# ======================================================================

@dataclass
class ManipResult:
    """Result of a data manipulation operation."""
    q: np.ndarray
    I: np.ndarray
    dI: np.ndarray
    dQ: Optional[np.ndarray] = None
    operation: str = ''
    metadata: dict = field(default_factory=dict)


# ======================================================================
# Engine
# ======================================================================

class DataManipulation:
    """Core data manipulation engine — no GUI dependencies.

    All public methods are static and return a :class:`ManipResult`.
    """

    # ------------------------------------------------------------------
    #  Scale + Background
    # ------------------------------------------------------------------

    @staticmethod
    def scale(
        q: np.ndarray,
        I: np.ndarray,
        dI: np.ndarray,
        dQ: Optional[np.ndarray],
        config: ScaleConfig,
    ) -> ManipResult:
        """Apply ``I_out = scale_I * I - background``.

        Uncertainty is scaled by *scale_uncertainty* (or *scale_I* if not set).
        """
        s = config.scale_I
        s_unc = config.scale_uncertainty if config.scale_uncertainty is not None else s
        I_out = s * I - config.background
        dI_out = np.abs(s_unc) * dI
        return ManipResult(
            q=q.copy(),
            I=I_out,
            dI=dI_out,
            dQ=dQ.copy() if dQ is not None else None,
            operation='scaled',
            metadata={
                'scale_I': s,
                'background': config.background,
                'scale_uncertainty': s_unc,
            },
        )

    # ------------------------------------------------------------------
    #  Trim Q range
    # ------------------------------------------------------------------

    @staticmethod
    def trim(
        q: np.ndarray,
        I: np.ndarray,
        dI: np.ndarray,
        dQ: Optional[np.ndarray],
        config: TrimConfig,
    ) -> ManipResult:
        """Keep only points where ``q_min <= Q <= q_max``."""
        mask = (q >= config.q_min) & (q <= config.q_max)
        return ManipResult(
            q=q[mask].copy(),
            I=I[mask].copy(),
            dI=dI[mask].copy(),
            dQ=dQ[mask].copy() if dQ is not None else None,
            operation='trimmed',
            metadata={'q_min': config.q_min, 'q_max': config.q_max},
        )

    # ------------------------------------------------------------------
    #  Rebin
    # ------------------------------------------------------------------

    @staticmethod
    def rebin(
        q: np.ndarray,
        I: np.ndarray,
        dI: np.ndarray,
        dQ: Optional[np.ndarray],
        config: RebinConfig,
    ) -> ManipResult:
        """Rebin data onto a new Q grid via log-log interpolation.

        Modes
        -----
        ``'log'``  : ``np.geomspace(q_min, q_max, n_points)``
        ``'linear'``: ``np.linspace(q_min, q_max, n_points)``
        ``'reference'``: use ``config.reference_q`` directly
        """
        if config.mode == 'reference':
            if config.reference_q is None:
                raise ValueError("reference_q must be provided for mode='reference'")
            q_new = config.reference_q.copy()
        else:
            q_lo = config.q_min if config.q_min is not None else q.min()
            q_hi = config.q_max if config.q_max is not None else q.max()
            if config.mode == 'log':
                q_new = np.geomspace(q_lo, q_hi, config.n_points)
            else:  # linear
                q_new = np.linspace(q_lo, q_hi, config.n_points)

        I_new = DataManipulation._interp_loglog(q_new, q, I)
        dI_new = DataManipulation._interp_uncertainty_loglog(q_new, q, I, dI)

        # Interpolate dQ if available
        dQ_new: Optional[np.ndarray] = None
        if dQ is not None:
            dQ_new = np.interp(
                np.log10(q_new), np.log10(q), dQ,
                left=float('nan'), right=float('nan'),
            )

        # Drop NaN points (outside source Q range)
        valid = np.isfinite(I_new) & np.isfinite(dI_new)
        return ManipResult(
            q=q_new[valid],
            I=I_new[valid],
            dI=dI_new[valid],
            dQ=dQ_new[valid] if dQ_new is not None else None,
            operation='rebinned',
            metadata={
                'mode': config.mode,
                'n_points': config.n_points,
                'q_min': float(q_new[valid].min()) if valid.any() else None,
                'q_max': float(q_new[valid].max()) if valid.any() else None,
            },
        )

    # ------------------------------------------------------------------
    #  Average
    # ------------------------------------------------------------------

    @staticmethod
    def average(
        datasets: list,
        reference_index: int = 0,
    ) -> ManipResult:
        """Average multiple datasets on a common Q grid.

        Parameters
        ----------
        datasets : list of (q, I, dI, dQ) tuples
            Each element is ``(np.ndarray, np.ndarray, np.ndarray,
            np.ndarray | None)``.
        reference_index : int
            Index of the dataset whose Q grid is used as the common grid.

        Returns
        -------
        ManipResult
            Uncertainty is ``max(propagated_error, spread)`` per Q point.
        """
        if not datasets:
            raise ValueError("At least one dataset is required")

        n = len(datasets)
        if n == 1:
            q, I, dI, dQ = datasets[0]
            return ManipResult(
                q=q.copy(), I=I.copy(), dI=dI.copy(),
                dQ=dQ.copy() if dQ is not None else None,
                operation='avg',
                metadata={'n_datasets': 1},
            )

        q_ref = datasets[reference_index][0]
        n_pts = len(q_ref)

        # Interpolate all datasets onto reference Q grid
        I_stack = np.empty((n, n_pts))
        dI_stack = np.empty((n, n_pts))
        for i, (qi, Ii, dIi, _) in enumerate(datasets):
            if i == reference_index:
                I_stack[i] = Ii
                dI_stack[i] = dIi
            else:
                I_stack[i] = DataManipulation._interp_loglog(q_ref, qi, Ii)
                dI_stack[i] = DataManipulation._interp_uncertainty_loglog(
                    q_ref, qi, Ii, dIi,
                )

        # Points valid in ALL datasets
        valid_all = np.all(np.isfinite(I_stack), axis=0) & np.all(
            np.isfinite(dI_stack), axis=0
        )

        I_avg = np.nanmean(I_stack[:, valid_all], axis=0)

        # Propagated uncertainty: sqrt(sum(dI^2)) / N
        dI_prop = np.sqrt(np.nansum(dI_stack[:, valid_all] ** 2, axis=0)) / n

        # Spread (sample std dev of values)
        dI_spread = np.nanstd(I_stack[:, valid_all], axis=0, ddof=0)

        # Use whichever is larger
        dI_avg = np.maximum(dI_prop, dI_spread)

        # dQ from reference dataset
        dQ_ref = datasets[reference_index][3]
        dQ_out = dQ_ref[valid_all].copy() if dQ_ref is not None else None

        return ManipResult(
            q=q_ref[valid_all].copy(),
            I=I_avg,
            dI=dI_avg,
            dQ=dQ_out,
            operation='avg',
            metadata={'n_datasets': n},
        )

    # ------------------------------------------------------------------
    #  Subtract (sample - buffer)
    # ------------------------------------------------------------------

    @staticmethod
    def subtract(
        q_sample: np.ndarray,
        I_sample: np.ndarray,
        dI_sample: np.ndarray,
        dQ_sample: Optional[np.ndarray],
        q_buffer: np.ndarray,
        I_buffer: np.ndarray,
        dI_buffer: np.ndarray,
        config: SubtractConfig,
    ) -> ManipResult:
        """``I_out = I_sample - buffer_scale * I_buffer``.

        Buffer is interpolated onto the sample's Q grid.
        If *auto_scale* is True, *buffer_scale* is computed from
        :meth:`auto_scale_buffer` over the cursor Q range first.
        """
        scale = config.buffer_scale
        if config.auto_scale and config.auto_q_min is not None and config.auto_q_max is not None:
            scale = DataManipulation.auto_scale_buffer(
                q_sample, I_sample, q_buffer, I_buffer,
                config.auto_q_min, config.auto_q_max,
            )

        # Interpolate buffer onto sample Q grid
        I_buf_interp = DataManipulation._interp_loglog(q_sample, q_buffer, I_buffer)
        dI_buf_interp = DataManipulation._interp_uncertainty_loglog(
            q_sample, q_buffer, I_buffer, dI_buffer,
        )

        # Only keep points where interpolation succeeded
        valid = np.isfinite(I_buf_interp) & np.isfinite(dI_buf_interp)

        q_out = q_sample[valid].copy()
        I_out = I_sample[valid] - scale * I_buf_interp[valid]
        dI_out = np.sqrt(dI_sample[valid] ** 2 + (scale * dI_buf_interp[valid]) ** 2)
        dQ_out = dQ_sample[valid].copy() if dQ_sample is not None else None

        return ManipResult(
            q=q_out,
            I=I_out,
            dI=dI_out,
            dQ=dQ_out,
            operation='sub',
            metadata={'buffer_scale': float(scale), 'auto_scale': config.auto_scale},
        )

    # ------------------------------------------------------------------
    #  Divide (numerator / denominator)
    # ------------------------------------------------------------------

    @staticmethod
    def divide(
        q_num: np.ndarray,
        I_num: np.ndarray,
        dI_num: np.ndarray,
        dQ_num: Optional[np.ndarray],
        q_den: np.ndarray,
        I_den: np.ndarray,
        dI_den: np.ndarray,
        config: DivideConfig,
    ) -> ManipResult:
        r"""Compute ``I_out = I_num / (scale * I_den - background)``.

        Denominator is interpolated onto the numerator's Q grid.
        """
        s = config.denominator_scale
        bg = config.denominator_background

        # Interpolate denominator onto numerator Q grid
        I_den_interp = DataManipulation._interp_loglog(q_num, q_den, I_den)
        dI_den_interp = DataManipulation._interp_uncertainty_loglog(
            q_num, q_den, I_den, dI_den,
        )

        valid = np.isfinite(I_den_interp) & np.isfinite(dI_den_interp)
        q_out = q_num[valid].copy()
        I_n = I_num[valid]
        dI_n = dI_num[valid]
        I_d = s * I_den_interp[valid] - bg
        dI_d = s * dI_den_interp[valid]

        # Avoid division by zero
        nonzero = I_d != 0.0
        I_out = np.full_like(I_n, np.nan)
        dI_out = np.full_like(I_n, np.nan)

        I_out[nonzero] = I_n[nonzero] / I_d[nonzero]
        # Relative error propagation: σ(A/B) = |A/B| * sqrt((σA/A)² + (σB/B)²)
        with np.errstate(divide='ignore', invalid='ignore'):
            rel_n = np.where(I_n[nonzero] != 0, dI_n[nonzero] / np.abs(I_n[nonzero]), 0.0)
            rel_d = np.where(I_d[nonzero] != 0, dI_d[nonzero] / np.abs(I_d[nonzero]), 0.0)
        dI_out[nonzero] = np.abs(I_out[nonzero]) * np.sqrt(rel_n ** 2 + rel_d ** 2)

        dQ_out = dQ_num[valid].copy() if dQ_num is not None else None

        return ManipResult(
            q=q_out,
            I=I_out,
            dI=dI_out,
            dQ=dQ_out,
            operation='div',
            metadata={
                'denominator_scale': s,
                'denominator_background': bg,
            },
        )

    # ------------------------------------------------------------------
    #  Auto-scale buffer
    # ------------------------------------------------------------------

    @staticmethod
    def auto_scale_buffer(
        q_sample: np.ndarray,
        I_sample: np.ndarray,
        q_buffer: np.ndarray,
        I_buffer: np.ndarray,
        q_min: float,
        q_max: float,
    ) -> float:
        """Compute buffer scale factor by matching integral intensity.

        Scale = integral(I_sample) / integral(I_buffer) over [q_min, q_max].
        Both datasets are evaluated on the sample's Q grid points that
        fall within the given range.
        """
        mask_s = (q_sample >= q_min) & (q_sample <= q_max)
        q_region = q_sample[mask_s]
        I_s_region = I_sample[mask_s]

        if len(q_region) < 2:
            return 1.0

        I_b_interp = DataManipulation._interp_loglog(q_region, q_buffer, I_buffer)
        valid = np.isfinite(I_b_interp)
        if valid.sum() < 2:
            return 1.0

        integral_s = np.trapezoid(I_s_region[valid], q_region[valid])
        integral_b = np.trapezoid(I_b_interp[valid], q_region[valid])

        if integral_b == 0.0:
            return 1.0

        return float(integral_s / integral_b)

    # ------------------------------------------------------------------
    #  Interpolation helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _interp_loglog(
        q_target: np.ndarray,
        q_known: np.ndarray,
        I_known: np.ndarray,
    ) -> np.ndarray:
        """Interpolate I at *q_target* using log-log linear interpolation.

        *q_known* and *I_known* must be positive and finite.
        Out-of-range *q_target* values return NaN (no extrapolation).
        Returns I values in linear units.
        """
        # Mask out non-positive values for log safety
        pos = (q_known > 0) & (I_known > 0) & np.isfinite(q_known) & np.isfinite(I_known)
        if pos.sum() < 2:
            return np.full_like(q_target, np.nan, dtype=float)
        log_I = np.interp(
            np.log10(q_target),
            np.log10(q_known[pos]),
            np.log10(I_known[pos]),
            left=float('nan'),
            right=float('nan'),
        )
        return 10.0 ** log_I

    @staticmethod
    def _interp_uncertainty_loglog(
        q_target: np.ndarray,
        q_known: np.ndarray,
        I_known: np.ndarray,
        dI_known: np.ndarray,
    ) -> np.ndarray:
        """Propagate uncertainty through log-log interpolation.

        Interpolates the *relative* uncertainty (dI/I) linearly in log10(Q),
        then multiplies by the interpolated I.  This preserves the
        signal-to-noise character of counting-statistics data.
        """
        pos = (q_known > 0) & (I_known > 0) & np.isfinite(q_known) & np.isfinite(I_known)
        if pos.sum() < 2:
            return np.full_like(q_target, np.nan, dtype=float)

        rel_unc = np.where(I_known[pos] > 0, dI_known[pos] / I_known[pos], 0.0)
        rel_target = np.interp(
            np.log10(q_target),
            np.log10(q_known[pos]),
            rel_unc,
            left=float('nan'),
            right=float('nan'),
        )
        I_target = DataManipulation._interp_loglog(q_target, q_known, I_known)
        return np.abs(rel_target * I_target)
