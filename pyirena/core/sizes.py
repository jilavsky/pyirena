"""
Particle size distribution fitting for SAS data.

Four inversion methods are provided, all sharing the same G-matrix framework:

  * Maximum Entropy (MaxEnt)  — entropy-regularised inversion; produces the
    smoothest distribution consistent with the data.
  * Tikhonov Regularization   — 2nd-derivative smoothness penalty; binary
    search on the regularisation parameter α.
  * TNNLS / IPG               — Interior-Point Gradient non-negative least
    squares; no explicit smoothness penalty but enforces x ≥ 0 strictly.
  * McSAS (Monte Carlo)       — Monte Carlo replacement algorithm; draws N_c
    random contributions and optimises them one at a time.  A single run is
    returned by ``fit()``; per-bin uncertainties come from the GUI's
    "Calculate Uncertainty" function (data perturbation, same as other methods).

The first three methods are faithful ports of the Igor Pro IR1R_sizes.ipf
package (Jan Ilavsky, Argonne National Laboratory).  McSAS follows the
algorithm of Bressler et al., J. Appl. Cryst. 48 (2015).

Usage example
-------------
>>> from pyirena.core.sizes import SizesDistribution
>>> import numpy as np
>>> s = SizesDistribution()
>>> s.r_min, s.r_max, s.n_bins = 50, 800, 60
>>> s.contrast = 1.0          # (Δρ)² in 10^20 cm^-4
>>> s.method = 'regularization'
>>> result = s.fit(q, intensity, error)
>>> print(result['chi_squared'], result['volume_fraction'])
"""

from __future__ import annotations

import logging
import math
from typing import Optional

import numpy as np
from scipy.optimize import nnls, curve_fit

from pyirena.core.form_factors import (
    build_g_matrix,
    make_r_grid,
    bin_widths,
)

log = logging.getLogger(__name__)


class SizesDistribution:
    """
    Particle size distribution model for SAS data.

    Attributes
    ----------
    r_min, r_max : float
        Radius range [Å] for the size grid.
    n_bins : int
        Number of radius bins.
    log_spacing : bool
        If True, bin centres are log-spaced.
    shape : str
        Particle shape: ``'sphere'`` or ``'spheroid'``.
    contrast : float
        Scattering contrast (Δρ)² in units of 10^20 cm^-4.
    shape_params : dict
        Extra shape parameters forwarded to the form factor.
        For spheroid: ``{'aspect_ratio': 1.5}``.
    background : float
        Constant background subtracted from the data before fitting [cm^-1].
    method : str
        Fitting method: ``'maxent'``, ``'regularization'``, or ``'tnnls'``.
    MaxEnt parameters
    -----------------
    maxent_sky_background : float
        Initial uniform model value (positive; acts as entropy prior).
    maxent_stability : float
        Convergence tolerance scale: tol = stability × √(2M).
    maxent_max_iter : int
        Maximum number of iterations.

    Regularization parameters
    -------------------------
    regularization_evalue : float
        Chi² tolerance scale: tol = evalue × √(2M).
    regularization_min_ratio : float
        Minimum distribution value = min_ratio × max(distribution).

    TNNLS parameters
    ----------------
    tnnls_approach_param : float
        Step safety factor (< 1); prevents x going negative.
    tnnls_max_iter : int
        Maximum number of iterations.

    McSAS parameters
    ----------------
    mcsas_n_repetitions : int
        Number of independent MC repetitions (default 10).  The mean and std
        of the per-bin volume fractions across repetitions are returned.
    mcsas_convergence : float
        Convergence criterion: stop the MC loop when χ²/M ≤ this value
        (default 1.0 = one standard deviation per point on average).
    mcsas_max_iter : int
        Maximum number of replacement iterations per repetition (default 100000).

    Results (populated after ``fit()``)
    ------------------------------------
    distribution : np.ndarray or None
        Normalised size distribution P(r)  [volume fraction / Å].
    distribution_std : np.ndarray or None
        Per-bin standard deviation of P(r) across McSAS repetitions.
        None for deterministic methods (MaxEnt, Regularization, TNNLS).
    r_grid : np.ndarray or None
        Radius bin centres [Å].
    model_intensity : np.ndarray or None
        Fitted scattering intensity [cm^-1].
    residuals : np.ndarray or None
        Normalised residuals (I_data - I_model) / error.
    chi_squared : float or None
        Sum of squared normalised residuals.
    volume_fraction : float or None
        Total volume fraction  ∫ P(r) dr.
    rg : float or None
        Volume-weighted radius of gyration [Å].
    n_iterations : int or None
        Number of iterations used by the fitting method.
    """

    # ── Grid ──────────────────────────────────────────────────────────────────
    r_min: float = 10.0
    r_max: float = 1000.0
    n_bins: int = 200
    log_spacing: bool = True

    # ── Shape ─────────────────────────────────────────────────────────────────
    shape: str = 'sphere'
    contrast: float = 1.0
    shape_params: dict

    # ── Background ────────────────────────────────────────────────────────────
    background: float = 0.0

    # ── Method ────────────────────────────────────────────────────────────────
    method: str = 'regularization'

    # ── MaxEnt parameters ─────────────────────────────────────────────────────
    maxent_sky_background: float = 1e-6
    maxent_stability: float = 0.01
    maxent_max_iter: int = 1000

    # ── Regularization parameters ─────────────────────────────────────────────
    regularization_evalue: float = 1.0
    regularization_min_ratio: float = 1e-4

    # ── TNNLS parameters ──────────────────────────────────────────────────────
    tnnls_approach_param: float = 0.95
    tnnls_max_iter: int = 1000

    # ── McSAS parameters ──────────────────────────────────────────────────────
    mcsas_n_repetitions: int = 10
    mcsas_convergence: float = 1.0
    mcsas_max_iter: int = 100000

    # ── Error scaling ─────────────────────────────────────────────────────────
    error_scale: float = 1.0   # Multiply measurement errors by this factor before fitting

    # ── Power-law background ───────────────────────────────────────────────────
    power_law_B: float = 0.0   # Amplitude of B·q^(-P) term  [same units as I]
    power_law_P: float = 4.0   # Exponent of power-law term

    # ── Results ───────────────────────────────────────────────────────────────
    distribution: Optional[np.ndarray]
    distribution_std: Optional[np.ndarray]
    r_grid: Optional[np.ndarray]
    model_intensity: Optional[np.ndarray]
    residuals: Optional[np.ndarray]
    chi_squared: Optional[float]
    volume_fraction: Optional[float]
    rg: Optional[float]
    n_iterations: Optional[int]

    def __init__(self):
        self.shape_params = {}
        self.distribution = None
        self.distribution_std = None
        self.r_grid = None
        self.model_intensity = None
        self.residuals = None
        self.chi_squared = None
        self.volume_fraction = None
        self.rg = None
        self.n_iterations = None

    # ── Public API ────────────────────────────────────────────────────────────

    def fit(
        self,
        q: np.ndarray,
        intensity: np.ndarray,
        error: Optional[np.ndarray] = None,
    ) -> dict:
        """
        Fit a size distribution to the experimental SAS data.

        Steps:
          1. Subtract background.
          2. Generate synthetic errors if none provided (I × 0.05).
          3. Build G matrix from the current shape settings.
          4. Apply Q-weighting if ``q_power > 0`` and no user errors.
          5. Run the selected fitting method.
          6. Post-process: normalise by bin width, compute chi², Rg, Vf.

        Args:
            q:         Q vector [Å^-1], 1-D.
            intensity: Measured I(Q) [cm^-1], 1-D.
            error:     Measurement uncertainty [cm^-1], 1-D.
                       If None, generated as ``intensity × 0.05``.

        Returns:
            dict with keys:
              ``success``, ``message``,
              ``distribution``, ``r_grid``, ``model_intensity``,
              ``residuals``, ``chi_squared``, ``volume_fraction``,
              ``rg``, ``n_iterations``.
        """
        q = np.asarray(q, dtype=float)
        I = np.asarray(intensity, dtype=float)

        # Remove NaN / non-finite points
        mask = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
        if error is not None:
            err = np.asarray(error, dtype=float)
            mask &= np.isfinite(err) & (err > 0)
        else:
            err = None

        q, I = q[mask], I[mask]
        err = err[mask] if err is not None else None

        if len(q) == 0:
            return self._fail("No valid data points after cleaning.")

        # Subtract complex background (B·q^(-P) + flat)
        I = I - self.compute_complex_background(q)

        # Generate errors if absent
        no_user_errors = err is None
        if no_user_errors:
            err = I * 0.05
            err[err <= 0] = np.abs(I[err <= 0]) * 0.05 + 1e-20

        # Apply user-specified error scaling (default 1.0 = no change)
        if self.error_scale != 1.0:
            err = err * float(self.error_scale)
            err = np.maximum(err, 1e-300)

        # Build G matrix
        r_grid = make_r_grid(self.r_min, self.r_max, self.n_bins, self.log_spacing)
        G = self._build_g_matrix(q, r_grid)

        # Fit
        try:
            method = self.method.lower()
            x_raw_std = None
            if method == 'maxent':
                x_raw, n_iter = self._fit_maxent(G, I, err)
            elif method == 'regularization':
                x_raw, n_iter = self._fit_regularization(G, I, err)
            elif method in ('tnnls', 'ipg', 'nnls'):
                x_raw, n_iter = self._fit_tnnls(G, I, err)
            elif method == 'mcsas':
                x_raw, n_iter, x_raw_std = self._fit_mcsas(G, I, err, r_grid)
            else:
                return self._fail(f"Unknown method '{self.method}'.")
        except Exception as exc:
            log.exception("Size distribution fit failed")
            return self._fail(str(exc))

        # Post-process
        result = self._post_process(x_raw, r_grid, G, I, err)

        # McSAS: the MC fit uses the standard G matrix (numerically stable, 6-decade
        # dynamic range) so x_raw = A × count reflects the number distribution.
        # Multiply by V(r) = (4/3)πr³ and renormalise to convert to volume distribution,
        # matching the output of MaxEnt, Regularisation, and TNNLS for consistent display.
        if method == 'mcsas':
            _V_r = (4.0 / 3.0) * np.pi * r_grid ** 3
            dist_vol = result['distribution'] * _V_r
            _vf_orig = result['volume_fraction']
            _vf_new = float(np.trapezoid(dist_vol, r_grid))
            _mcsas_scale = (_vf_orig / _vf_new) if _vf_new > 0 else 1.0
            dist_vol *= _mcsas_scale
            result['distribution'] = dist_vol
            result['volume_fraction'] = float(np.trapezoid(dist_vol, r_grid))
            if result['volume_fraction'] > 0:
                result['rg'] = float(np.sqrt(
                    np.trapezoid(r_grid ** 2 * dist_vol, r_grid) / result['volume_fraction']
                ))
            else:
                result['rg'] = 0.0

        # McSAS: add per-bin uncertainty from spread across repetitions
        if x_raw_std is not None:
            dw_safe = np.maximum(bin_widths(r_grid), 1e-300)
            dist_std = x_raw_std / dw_safe
            if method == 'mcsas':
                dist_std = dist_std * _V_r * _mcsas_scale
            result['distribution_std'] = dist_std
        else:
            result['distribution_std'] = None

        result['n_iterations'] = n_iter
        result['n_data'] = len(I)   # number of Q points used; chi² target ≈ this value
        result['success'] = True
        result['message'] = (
            f"Fit converged in {n_iter} iterations; "
            f"χ² = {result['chi_squared']:.4f}"
        )

        # Store results on self
        for k in ('distribution', 'distribution_std', 'r_grid', 'model_intensity',
                  'residuals', 'chi_squared', 'volume_fraction', 'rg', 'n_iterations'):
            setattr(self, k, result.get(k))
        self.r_grid = r_grid

        return result

    # ── G-matrix helper ───────────────────────────────────────────────────────

    def _build_g_matrix(self, q: np.ndarray, r_grid: np.ndarray) -> np.ndarray:
        return build_g_matrix(q, r_grid, self.shape, self.contrast, **self.shape_params)

    # ── Complex background ─────────────────────────────────────────────────────

    def compute_complex_background(self, q: np.ndarray) -> np.ndarray:
        """Return B·q^(-P) + background for each q value.

        When ``power_law_B`` is zero only the flat ``background`` is returned.
        """
        q = np.asarray(q, dtype=float)
        bg = np.full(len(q), self.background, dtype=float)
        if self.power_law_B != 0.0:
            with np.errstate(divide='ignore', invalid='ignore'):
                pl = self.power_law_B * np.power(q, -self.power_law_P)
            bg += np.where(np.isfinite(pl), pl, 0.0)
        return bg

    def fit_power_law(
        self,
        q: np.ndarray,
        intensity: np.ndarray,
        q_min: float,
        q_max: float,
        fit_B: bool = True,
        fit_P: bool = True,
    ) -> dict:
        """Fit B·q^(-P) to data in [q_min, q_max].

        ``fit_B`` / ``fit_P`` select which parameters are free (at least one
        must be True).  Updates ``self.power_law_B`` and/or ``self.power_law_P``
        on success.

        Returns dict with keys ``'success'``, ``'B'``, ``'P'``, ``'message'``.
        """
        q = np.asarray(q, dtype=float)
        I = np.asarray(intensity, dtype=float)

        mask = ((q >= q_min) & (q <= q_max)
                & np.isfinite(q) & np.isfinite(I)
                & (q > 0) & (I > 0))
        qf, If = q[mask], I[mask]

        if len(qf) < 2:
            return {
                'success': False,
                'message': 'Too few data points in Q range for power-law fit.',
                'B': self.power_law_B, 'P': self.power_law_P,
            }

        try:
            if fit_B and fit_P:
                def model(q_, B, P):
                    return B * np.power(q_, -P)
                p0 = [max(self.power_law_B, float(If.mean())), self.power_law_P]
                popt, _ = curve_fit(model, qf, If, p0=p0,
                                    bounds=([0.0, 0.1], [np.inf, 12.0]),
                                    maxfev=10000)
                B_new, P_new = float(popt[0]), float(popt[1])

            elif fit_B:
                P_fixed = self.power_law_P
                def model(q_, B):
                    return B * np.power(q_, -P_fixed)
                p0 = [max(self.power_law_B, float(If.mean()))]
                popt, _ = curve_fit(model, qf, If, p0=p0,
                                    bounds=([0.0], [np.inf]), maxfev=10000)
                B_new, P_new = float(popt[0]), self.power_law_P

            elif fit_P:
                B_fixed = self.power_law_B
                if B_fixed == 0.0:
                    return {
                        'success': False,
                        'message': 'B is zero — set a non-zero B before fitting P alone.',
                        'B': self.power_law_B, 'P': self.power_law_P,
                    }
                def model(q_, P):
                    return B_fixed * np.power(q_, -P)
                p0 = [self.power_law_P]
                popt, _ = curve_fit(model, qf, If, p0=p0,
                                    bounds=([0.1], [12.0]), maxfev=10000)
                B_new, P_new = self.power_law_B, float(popt[0])

            else:
                return {
                    'success': False,
                    'message': 'No parameters selected for fitting (check Fit B or Fit P).',
                    'B': self.power_law_B, 'P': self.power_law_P,
                }

            self.power_law_B = B_new
            self.power_law_P = P_new
            return {
                'success': True,
                'message': f'Power-law fit: B = {B_new:.4g}, P = {P_new:.4g}',
                'B': B_new, 'P': P_new,
            }
        except Exception as exc:
            return {
                'success': False,
                'message': f'Power-law fit failed: {exc}',
                'B': self.power_law_B, 'P': self.power_law_P,
            }

    def fit_background_term(
        self,
        q: np.ndarray,
        intensity: np.ndarray,
        q_min: float,
        q_max: float,
    ) -> dict:
        """Fit flat background by averaging I − B·q^(-P) in [q_min, q_max].

        Updates ``self.background`` on success.
        Returns dict with keys ``'success'``, ``'background'``, ``'message'``.
        """
        q = np.asarray(q, dtype=float)
        I = np.asarray(intensity, dtype=float)

        mask = ((q >= q_min) & (q <= q_max)
                & np.isfinite(q) & np.isfinite(I)
                & (q > 0))
        qf, If = q[mask], I[mask]

        if len(qf) < 1:
            return {
                'success': False,
                'message': 'No data points in Q range for background fit.',
                'background': self.background,
            }

        with np.errstate(divide='ignore', invalid='ignore'):
            if self.power_law_B != 0.0:
                pl = self.power_law_B * np.power(qf, -self.power_law_P)
                pl = np.where(np.isfinite(pl), pl, 0.0)
            else:
                pl = np.zeros(len(qf))

        bg = float(np.mean(If - pl))
        self.background = bg
        return {
            'success': True,
            'message': f'Background fit: {bg:.4g} cm⁻¹',
            'background': bg,
        }

    # ── Maximum Entropy ───────────────────────────────────────────────────────

    def _fit_maxent(
        self, G: np.ndarray, I: np.ndarray, err: np.ndarray
    ) -> tuple[np.ndarray, int]:
        """
        Maximum Entropy inversion — faithful port of IR1R_MaximumEntropy (Igor Pro).

        Maximises the Skilling-Bryan entropy S = Σ(x − sky − x·ln(x/sky))
        subject to χ² = M (number of data points).

        Uses three model-space search directions with x-weighted inner products,
        and IR1R_Move bisection to find the optimal step at each iteration.
        """
        M = len(I)
        N = G.shape[1]
        sky = float(self.maxent_sky_background)
        tolerance = self.maxent_stability * math.sqrt(2.0 * M)
        tst_lim = 0.05
        chizer = float(M)   # chi² target = number of data points

        # Initialise model at sky background level
        x = np.full(N, sky, dtype=float)
        err2 = err ** 2

        # Working arrays (model space N×3, data space M×3)
        xi  = np.zeros((N, 3))
        eta = np.zeros((M, 3))
        c1 = np.zeros(3)
        s1 = np.zeros(3)
        c2 = np.zeros((3, 3))
        s2 = np.zeros((3, 3))

        test = 0.0
        chi2_cur = 0.0
        n_iter = 0

        for n_iter in range(self.maxent_max_iter):
            # ── Forward model and chi² ────────────────────────────────────────
            ox = G @ x                              # [M] predicted
            ascratch = (ox - I) / err              # [M] (Gx − data) / err
            ox_scaled = 2.0 * ascratch / err       # [M] 2·(Gx − data)/err²
            chi2_cur = float(np.dot(ascratch, ascratch))

            if chi2_cur < 1e-300:
                break

            # chi² gradient in model space (UPHILL: points toward increasing χ²)
            cgrad = G.T @ ox_scaled                # [N]  G^T · (2·(Gx−I)/err²)

            # ── Entropy gradient (Skilling-Bryan): ∂S/∂xᵢ = −ln(xᵢ/sky) ──────
            # Igor scales by 1/(sky·e) so the gradient is zero at x = sky
            x_safe = np.maximum(x, sky * 1e-10)
            sgrad = -np.log(x_safe / sky) / (sky * math.e)   # [N]

            # ── x-weighted norms (natural gradient metric) ────────────────────
            snorm = math.sqrt(max(0.0, float(np.sum(x * sgrad ** 2))))
            cnorm = math.sqrt(max(0.0, float(np.sum(x * cgrad ** 2))))
            tnorm = float(np.sum(x * sgrad * cgrad))

            if cnorm < 1e-300:
                break

            # ── Test statistic: alignment of entropy & chi² gradients ─────────
            # Converged when test → 0 (gradients parallel → saddle point)
            if n_iter > 0:
                denom = snorm * cnorm
                test = math.sqrt(max(0.0, 0.5 * (1.0 - tnorm / denom))) if denom > 1e-300 else 1.0

            # ── Parameters a, b for search directions ─────────────────────────
            # iter 0: a=1, b=1/cnorm  (pure chi² direction + un-scaled entropy)
            # iter>0: scale to make xi[:,1] ≈ orthogonal to xi[:,0] in x-metric
            a = 1.0
            b = 1.0 / cnorm
            if n_iter > 0 and test > 1e-10:
                a = 0.5 / (snorm * test)
                b = 0.5 / (cnorm * test)

            # ── Search directions in model space ──────────────────────────────
            xi[:, 0] = x * cgrad / cnorm
            xi[:, 1] = x * (a * sgrad - b * cgrad)

            # Forward-map xi[:,0,1] to data space
            eta[:, 0] = G @ xi[:, 0]
            eta[:, 1] = G @ xi[:, 1]

            # Third direction: TrOpus(eta[:,1] / err²), then x-weighted normalise
            xiscratch = G.T @ (eta[:, 1] / err2)
            a_val = float(np.sum(x * xiscratch ** 2))
            if a_val > 0.0:
                xi[:, 2] = x * xiscratch / math.sqrt(a_val)
            else:
                xi[:, 2] = 0.0
            eta[:, 2] = G @ xi[:, 2]

            # ── Build 3×3 matrices c1, s1, c2, s2 ────────────────────────────
            for i in range(3):
                s1[i] = float(np.dot(xi[:, i], sgrad))
                c1[i] = float(np.dot(xi[:, i], cgrad))
            c1 /= chi2_cur          # normalise c1 by χ²

            s2[:] = 0.0
            c2[:] = 0.0
            for k in range(3):
                for l in range(k + 1):
                    s2[k, l] = -float(np.sum(xi[:, k] * xi[:, l] / x_safe))
                    c2[k, l] =  float(np.sum(eta[:, k] * eta[:, l] / err2))
            s2 /= sky
            c2 *= 2.0 / chi2_cur
            # Symmetrise upper triangle
            for k in range(3):
                for l in range(k):
                    s2[l, k] = s2[k, l]
                    c2[l, k] = c2[k, l]

            # ── Compute optimal step beta ─────────────────────────────────────
            if n_iter == 0:
                # First iteration: simple chi² gradient step, no bisection
                beta = np.zeros(3)
                if abs(c2[0, 0]) > 1e-300:
                    beta[0] = -0.5 * c1[0] / c2[0, 0]
            else:
                beta = self._maxent_move(
                    c1, s1, c2, s2, chi2_cur, chizer, float(x.sum()), sky
                )

            # ── Update model ──────────────────────────────────────────────────
            for i in range(N):
                df = beta[0]*xi[i, 0] + beta[1]*xi[i, 1] + beta[2]*xi[i, 2]
                if df < -x[i]:
                    df = 0.001 * sky - x[i]   # floor patch: keep x ≥ 0.001·sky
                x[i] += df

            # ── Recompute chi² on updated model ──────────────────────────────
            resid2 = (I - G @ x) / err
            chi2_new = float(np.dot(resid2, resid2))

            if n_iter % 10 == 0 or n_iter < 5:
                log.debug(
                    "MaxEnt iter %d: chi2=%.4g  target=%.4g  test=%.4g",
                    n_iter, chi2_new, chizer, test,
                )

            # ── Convergence check ─────────────────────────────────────────────
            if abs(chi2_new - chizer) < tolerance and test < tst_lim:
                log.debug("MaxEnt converged at iter %d, chi2=%.4g", n_iter, chi2_new)
                break

        log.debug(
            "MaxEnt finished: %d iters, chi2=%.4g, target=%.4g",
            n_iter, float(np.sum(((I - G @ x) / err) ** 2)), chizer,
        )
        return x, n_iter

    def _maxent_move(
        self,
        c1: np.ndarray, s1: np.ndarray,
        c2: np.ndarray, s2: np.ndarray,
        chi2: float, chizer: float,
        fSum: float, sky: float,
    ) -> np.ndarray:
        """
        IR1R_Move: find optimal step beta by bisection on ax ∈ [0, 1].

        At ax=0: pure chi² step.  At ax=1: pure entropy step.
        Bisects to find the ax that brings chi²_predicted to the target chizer.
        Then applies IR1R_Dist step-magnitude limiter.
        """
        MAX_LOOP = 500
        PASSES = 0.001      # bisection convergence tolerance

        def chi_now(ax: float):
            """
            IR1R_ChiNow: solve (bx·c2 − ax·s2)·beta = −(bx·c1 − ax·s1)
            and return (ChiNow, beta) where ChiNow = χ²_new / χ²_current.
            """
            bx = 1.0 - ax
            A = bx * c2 - ax * s2
            b_vec = -(bx * c1 - ax * s1)
            try:
                beta_loc = np.linalg.solve(A, b_vec)
            except np.linalg.LinAlgError:
                beta_loc, _, _, _ = np.linalg.lstsq(A, b_vec, rcond=None)
            # ChiNow = 1 + beta^T (c1 + 0.5 c2 beta)
            w = float(np.dot(beta_loc, c1 + 0.5 * (c2 @ beta_loc)))
            return 1.0 + w, beta_loc

        # Initial bracket
        a1, a2 = 0.0, 1.0
        cmin, beta_final = chi_now(a1)

        # Choose bisection target based on whether pure chi² step overshoots
        if cmin * chi2 > chizer:
            ctarg = 0.5 * (1.0 + cmin)
        else:
            ctarg = chizer / chi2

        f1 = cmin - ctarg
        cn2, _ = chi_now(a2)
        f2 = cn2 - ctarg

        # Bisection
        for _ in range(MAX_LOOP):
            anew = 0.5 * (a1 + a2)
            fx_val, beta_new = chi_now(anew)
            fx = fx_val - ctarg
            if f1 * fx > 0:
                a1, f1 = anew, fx
            if f2 * fx > 0:
                a2, f2 = anew, fx
            beta_final = beta_new
            if abs(fx) < PASSES:
                break

        # IR1R_Dist: w = −beta^T s2 beta  (step magnitude in entropy space)
        # s2 is negative-definite so −s2 is positive-definite → w ≥ 0
        w = float(-np.dot(beta_final, s2 @ beta_final))

        # Step limiter: scale down if step is too large relative to model total
        limit = 0.1 * fSum / sky if sky > 0 else 1e300
        if w > limit:
            beta_final = beta_final * math.sqrt(limit / w)

        return beta_final

    # ── Tikhonov Regularization ───────────────────────────────────────────────

    def _fit_regularization(
        self, G: np.ndarray, I: np.ndarray, err: np.ndarray
    ) -> tuple[np.ndarray, int]:
        """
        Tikhonov regularization (port of IR1R_DoInternalRegularization).

        Minimises  ||W^{1/2}(I - G·x)||² + α ||L·x||²  subject to  x ≥ 0,

        where L is the 2nd-order finite-difference matrix (so LᵀL ≈ H, the
        4th-order smoothness operator) and α is found by binary search in
        log₁₀ space to achieve  χ² ≈ M.

        Each trial uses ``scipy.optimize.nnls`` (Lawson-Hanson) on the
        augmented system [W^{1/2}G; √α·L] to correctly enforce non-negativity.
        This avoids the oscillation artefacts that arise from unconstrained
        lstsq followed by clamping to zero, which can produce pathologically
        large chi² at small alpha for ill-conditioned G matrices.

        If the chi² target is genuinely unachievable (minimum chi² > M), the
        method falls back to the L-curve elbow (minimum chi² alpha).
        """
        M, N = G.shape
        chi_target = float(M)
        tol = self.regularization_evalue * np.sqrt(2.0 * M)
        min_ratio = self.regularization_min_ratio

        # ── 2nd-order FD matrix L: (N-2) × N with rows [1, -2, 1]  ──────────
        L = self._make_l_matrix(N)
        nL = L.shape[0]

        # Error weights: W^{1/2} = diag(1/err)
        W12 = 1.0 / err
        Gw = G * W12[:, np.newaxis]   # [M, N]
        Iw = I * W12                   # [M]
        zeros_L = np.zeros(nL)

        def _solve(log_alpha: float) -> tuple[np.ndarray, float]:
            """Solve non-negative regularised LS and return (x, chi²)."""
            alpha = 10.0 ** float(np.clip(log_alpha, -300, 300))
            sa = np.sqrt(alpha)
            A_aug = np.vstack([Gw, sa * L])
            b_aug = np.concatenate([Iw, zeros_L])
            x, _ = nnls(A_aug, b_aug, maxiter=10 * N)
            resid = (I - G @ x) / err
            c2 = float(np.dot(resid, resid))
            return x, c2 if np.isfinite(c2) else 1e300

        # ── Establish bracket: lo → chi² ≤ target, hi → chi² ≥ target ────────
        lo, hi = -5.0, 15.0
        _, c2_lo = _solve(lo)
        _, c2_hi = _solve(hi)

        # Shift lo toward smaller alpha until chi² drops to or below target.
        # With NNLS, chi²(alpha) is monotone: decreasing as alpha → 0.
        for _ in range(30):
            if c2_lo <= chi_target:
                break
            lo -= 2.0
            _, c2_lo = _solve(lo)

        # Shift hi toward larger alpha until chi² exceeds target.
        for _ in range(20):
            if c2_hi >= chi_target:
                break
            hi += 2.0
            _, c2_hi = _solve(hi)

        # ── Fallback: chi² target not achievable ──────────────────────────────
        if c2_lo > chi_target:
            # The minimum achievable chi² exceeds chi_target.
            # Find the L-curve elbow (minimum chi²) by a coarse scan and
            # return that solution instead.
            log.warning(
                "Regularization: chi² target %.1f not achievable "
                "(minimum chi² %.1f); using L-curve elbow.", chi_target, c2_lo
            )
            scan_alphas = np.arange(-5.0, 16.0, 1.0)
            scan_c2 = [_solve(la)[1] for la in scan_alphas]
            best_log_alpha = float(scan_alphas[int(np.argmin(scan_c2))])
            x, _ = _solve(best_log_alpha)
            if x.max() > 0:
                x = np.maximum(x, min_ratio * x.max())
            return x, 0

        # ── Binary search on log₁₀(α) ────────────────────────────────────────
        n_iter = 0
        for n_iter in range(1, 100):
            mid = 0.5 * (lo + hi)
            x, c2 = _solve(mid)

            if abs(c2 - chi_target) < tol:
                log.debug(f"Regularization converged at iteration {n_iter}, log_α={mid:.2f}")
                break

            # Large α → more smoothing → higher χ²
            if c2 > chi_target:
                hi = mid
            else:
                lo = mid

        # ── Final solution with floor applied ─────────────────────────────────
        x, _ = _solve(0.5 * (lo + hi))
        if x.max() > 0:
            x = np.maximum(x, min_ratio * x.max())

        return x, n_iter

    @staticmethod
    def _make_l_matrix(N: int) -> np.ndarray:
        """
        2nd-order finite-difference matrix L of shape (max(N-2,1), N).

        Each interior row has the pattern [1, -2, 1], so LᵀL approximates
        the 4th-order smoothness operator H.  Using L in the augmented system
        [W^{1/2}G; √α·L] is numerically much stabler than building H = LᵀL
        explicitly.
        """
        n_rows = max(N - 2, 1)
        L = np.zeros((n_rows, N), dtype=float)
        for i in range(min(N - 2, n_rows)):
            L[i, i]     =  1.0
            L[i, i + 1] = -2.0
            L[i, i + 2] =  1.0
        return L

    @staticmethod
    def _make_h_matrix(N: int) -> np.ndarray:
        """
        4th-order finite-difference smoothness matrix H [N×N].

        In the interior: row pattern [1, -4, 6, -4, 1] (discrete 4th derivative).
        Boundary rows use lower-order differences to stay within bounds.
        Matches Igor IR1R_MakeHmatrix.
        """
        H = np.zeros((N, N), dtype=float)

        if N < 3:
            return H

        for i in range(N):
            if i == 0:
                H[i, 0] = 1.0
                if N > 1:
                    H[i, 1] = -2.0
                if N > 2:
                    H[i, 2] = 1.0
            elif i == 1:
                H[i, 0] = -2.0
                H[i, 1] = 5.0
                if N > 2:
                    H[i, 2] = -4.0
                if N > 3:
                    H[i, 3] = 1.0
            elif i == N - 2:
                if N > 3:
                    H[i, N - 4] = 1.0
                H[i, N - 3] = -4.0
                H[i, N - 2] = 5.0
                H[i, N - 1] = -2.0
            elif i == N - 1:
                if N > 2:
                    H[i, N - 3] = 1.0
                H[i, N - 2] = -2.0
                H[i, N - 1] = 1.0
            else:
                H[i, i - 2] = 1.0
                H[i, i - 1] = -4.0
                H[i, i    ] = 6.0
                H[i, i + 1] = -4.0
                H[i, i + 2] = 1.0

        return H

    # ── TNNLS / IPG ───────────────────────────────────────────────────────────

    def _fit_tnnls(
        self,
        G: np.ndarray,
        I: np.ndarray,
        err: np.ndarray,
    ) -> tuple[np.ndarray, int]:
        """
        Interior-Point Gradient (IPG) non-negative least squares
        (port of IR1R_TNNLS in Igor Pro).

        Minimises  ||A·x - b||²  subject to  x ≥ 0,
        where A = G/err (and optionally Q-weighted) and b = I/err.

        The preconditioned gradient direction is scaled component-wise by
        x / (AᵀA·x) to stay in the positive orthant.
        """
        M, N = G.shape
        approach = self.tnnls_approach_param

        # Scale by errors
        A = G / err[:, np.newaxis]   # [M, N]
        b = I / err                  # [M]

        # Precompute
        AtA = A.T @ A    # [N, N]
        Atb = A.T @ b    # [N]

        # Initialise near zero (positive)
        x = np.full(N, 1e-32, dtype=float)

        n_iter = 0
        for n_iter in range(1, self.tnnls_max_iter + 1):
            Atax = AtA @ x              # [N]

            # Gradient of ||A·x - b||²:  2(Aᵀ A x - Aᵀ b)
            grad = Atax - Atb           # [N]

            # Diagonal preconditioner:  D_k = x / (AᵀA·x)
            denom = np.maximum(Atax, 1e-300)
            D = x / denom               # [N]

            # Search direction:  P = -D · grad  (element-wise)
            P = -D * grad               # [N]

            # Optimal step size (exact line search for quadratic objective)
            AtP = AtA @ P               # [N]
            denom_alpha = float(P @ AtP)
            if abs(denom_alpha) < 1e-300:
                break

            alpha_star = -float(P @ grad) / denom_alpha
            alpha_star = max(alpha_star, 0.0)

            # Step limiting: prevent any component from going negative
            for i in range(N):
                if P[i] < 0.0:
                    alpha_max = approach * (-x[i] / P[i])
                    if alpha_star > alpha_max:
                        alpha_star = alpha_max

            x = x + alpha_star * P
            x = np.maximum(x, 1e-300)  # numerical floor

            # Convergence check: reduced chi²
            resid = I - G @ x
            c2 = float(np.dot(resid / err, resid / err))
            if c2 / M <= 1.0:
                log.debug(f"TNNLS converged at iteration {n_iter}")
                break

        return x, n_iter

    # ── McSAS / Monte-Carlo ────────────────────────────────────────────────────

    def _fit_mcsas(
        self,
        G: np.ndarray,
        I: np.ndarray,
        err: np.ndarray,
        r_grid: np.ndarray,
        n_reps: Optional[int] = None,
    ) -> tuple[np.ndarray, int, Optional[np.ndarray]]:
        """
        Monte Carlo Size Analysis (McSAS) inversion.

        Uses the same G matrix and r_grid as the other methods.  N_c = n_bins
        contributions are randomly assigned to r_grid bins and then optimised
        by single-contribution replacement to minimise χ².

        The forward model is the standard intensity model:

            I(q) = A × Σ_k count[k] × G[:,k]

        where G[:,k] = V(r_k) × F²_norm(q, r_k) × contrast × 1e-4  [cm⁻¹ per
        unit volume fraction].  The optimal scale factor A (least-squares, ≥ 0)
        equals  Vf / N_c  for a monodisperse sample, so:

            x_raw[k] = A × count[k]   [volume fraction per bin]

        This is identical to the x_raw produced by MaxEnt, Regularization, and
        TNNLS, so all four methods feed into the same ``_post_process`` pipeline
        and produce consistent volume distributions.

        Note: an earlier implementation used volume-weighted G columns
        (G_V[:,k] = G[:,k]×V(r_k)) to make the number-density interpretation
        explicit, but that caused a 12-decade dynamic-range problem with log-
        spaced grids, leading to ringing artefacts.  The standard G model has
        only a 6-decade range (r³ instead of r⁶) and is numerically stable.

        Parameters
        ----------
        G : (M, N) array
            Pre-built G matrix from ``_build_g_matrix``.
        I : (M,) array
            Background-subtracted, error-scaled experimental intensity.
        err : (M,) array
            Measurement uncertainties after scaling.
        r_grid : (N,) array
            Radius bin centres [Å].
        n_reps : int, optional
            Number of independent MC repetitions.  Defaults to
            ``self.mcsas_n_repetitions``.  Pass ``1`` for a single-run fit
            (the ``distribution_std`` return value will be ``None``).

        Returns
        -------
        x_raw_mean : (N,) array
            Mean volume-fraction per bin across repetitions.
        n_iter_avg : int
            Average number of MC iterations per repetition.
        x_raw_std : (N,) array or None
            Std of volume-fraction per bin across repetitions.
            ``None`` when *n_reps* == 1 (no meaningful spread from a single run).
        """
        M, N = G.shape
        N_c = N   # one contribution per r_grid bin; they can cluster during MC
        n_rep = int(n_reps if n_reps is not None else self.mcsas_n_repetitions)
        convergence = float(self.mcsas_convergence)
        max_iter = int(self.mcsas_max_iter)

        rng = np.random.default_rng()
        inv_err2 = 1.0 / (err ** 2)               # (M,)

        def _scale_factor(g_sum: np.ndarray) -> float:
            """Optimal A = (g·I/σ²)/(g·g/σ²), clamped ≥ 0."""
            num = float(np.dot(g_sum * inv_err2, I))
            den = float(np.dot(g_sum * inv_err2, g_sum))
            return max(num / den, 0.0) if den > 1e-300 else 0.0

        def _chi2(g_sum: np.ndarray, A: float) -> float:
            resid = (I - A * g_sum) / err
            return float(np.dot(resid, resid))

        x_raw_reps = np.zeros((n_rep, N), dtype=float)
        total_iters = 0

        for rep in range(n_rep):
            # ── Initialise: assign each of N_c contributions to a random bin ──
            contrib_bins = rng.integers(0, N, size=N_c)      # (N_c,) bin indices
            counts = np.bincount(contrib_bins, minlength=N).astype(float)  # (N,)

            g_sum = G @ counts                               # (M,)
            A = _scale_factor(g_sum)
            chi2_val = _chi2(g_sum, A)

            # ── MC replacement loop ────────────────────────────────────────────
            rep_iters = 0
            for rep_iters in range(1, max_iter + 1):
                j = int(rng.integers(N_c))
                k_old = int(contrib_bins[j])
                k_new = int(rng.integers(N))

                if k_new == k_old:
                    continue

                # Incremental update (O(M) not O(M×N))
                g_sum_trial = g_sum + G[:, k_new] - G[:, k_old]
                A_trial = _scale_factor(g_sum_trial)
                chi2_trial = _chi2(g_sum_trial, A_trial)

                if chi2_trial <= chi2_val:
                    contrib_bins[j] = k_new
                    counts[k_old] -= 1.0
                    counts[k_new] += 1.0
                    g_sum = g_sum_trial
                    A = A_trial
                    chi2_val = chi2_trial

                if chi2_val / M <= convergence:
                    break

            total_iters += rep_iters
            log.debug(
                "McSAS rep %d/%d: %d iters, chi2=%.4g (target %.4g)",
                rep + 1, n_rep, rep_iters, chi2_val, convergence * M,
            )

            # x_raw[k] = A × count[k]  →  volume fraction per bin (same as other methods)
            x_raw_reps[rep] = A * counts

        x_raw_mean = np.mean(x_raw_reps, axis=0)
        x_raw_std  = np.std(x_raw_reps, axis=0) if n_rep > 1 else None
        n_iter_avg = total_iters // max(n_rep, 1)

        return x_raw_mean, n_iter_avg, x_raw_std

    # ── Post-processing ───────────────────────────────────────────────────────

    def _post_process(
        self,
        x_raw: np.ndarray,
        r_grid: np.ndarray,
        G: np.ndarray,
        I: np.ndarray,
        err: np.ndarray,
    ) -> dict:
        """
        Convert raw solution vector to size distribution and compute diagnostics.

        Normalisation: divide each bin value by the bin width so that
          ∫ P(r) dr = total volume fraction (discrete trapz).
        """
        dw = bin_widths(r_grid)

        # Normalise to density  [volume fraction / Å]
        # Guard against zero bin widths
        dw_safe = np.maximum(dw, 1e-300)
        distribution = x_raw / dw_safe

        I_model = G @ x_raw
        residuals = (I - I_model) / err
        chi_squared = float(np.dot(residuals, residuals))

        vol_fraction = float(np.trapezoid(distribution, r_grid))

        # Volume-weighted Rg:  Rg = √(∫r²P(r)dr / ∫P(r)dr)
        if vol_fraction > 0:
            rg = float(np.sqrt(np.trapezoid(r_grid ** 2 * distribution, r_grid) / vol_fraction))
        else:
            rg = 0.0

        return {
            'distribution':  distribution,
            'r_grid':        r_grid,
            'model_intensity': I_model,
            'residuals':     residuals,
            'chi_squared':   chi_squared,
            'volume_fraction': vol_fraction,
            'rg':            rg,
        }

    # ── Utilities ─────────────────────────────────────────────────────────────

    def calculate_volume_fraction(self) -> float:
        """Integrate distribution over the radius grid (trapezoidal rule)."""
        if self.distribution is None or self.r_grid is None:
            return 0.0
        return float(np.trapezoid(self.distribution, self.r_grid))

    def calculate_rg(self) -> float:
        """
        Volume-weighted radius of gyration [Å].

        Rg = √( ∫r² P(r) dr / ∫ P(r) dr )
        """
        if self.distribution is None or self.r_grid is None:
            return 0.0
        vf = self.calculate_volume_fraction()
        if vf <= 0:
            return 0.0
        return float(
            np.sqrt(np.trapezoid(self.r_grid ** 2 * self.distribution, self.r_grid) / vf)
        )

    def to_dict(self) -> dict:
        """Serialise all parameters (not results) to a plain dict."""
        return {
            'r_min':                    self.r_min,
            'r_max':                    self.r_max,
            'n_bins':                   self.n_bins,
            'log_spacing':              self.log_spacing,
            'shape':                    self.shape,
            'contrast':                 self.contrast,
            'shape_params':             dict(self.shape_params),
            'background':               self.background,
            'method':                   self.method,
            'maxent_sky_background':    self.maxent_sky_background,
            'maxent_stability':         self.maxent_stability,
            'maxent_max_iter':          self.maxent_max_iter,
            'regularization_evalue':    self.regularization_evalue,
            'regularization_min_ratio': self.regularization_min_ratio,
            'tnnls_approach_param':     self.tnnls_approach_param,
            'tnnls_max_iter':           self.tnnls_max_iter,
            'mcsas_n_repetitions':      self.mcsas_n_repetitions,
            'mcsas_convergence':        self.mcsas_convergence,
            'mcsas_max_iter':           self.mcsas_max_iter,
            'error_scale':              self.error_scale,
            'power_law_B':              self.power_law_B,
            'power_law_P':              self.power_law_P,
        }

    @classmethod
    def from_dict(cls, d: dict) -> 'SizesDistribution':
        """Reconstruct a ``SizesDistribution`` from a serialised dict."""
        obj = cls()
        for key, val in d.items():
            if hasattr(obj, key):
                setattr(obj, key, val)
        return obj

    # ── Private helpers ───────────────────────────────────────────────────────

    @staticmethod
    def _fail(message: str) -> dict:
        return {
            'success':           False,
            'message':           message,
            'distribution':      None,
            'distribution_std':  None,
            'r_grid':            None,
            'model_intensity':   None,
            'residuals':         None,
            'chi_squared':       None,
            'volume_fraction':   None,
            'rg':                None,
            'n_iterations':      None,
        }
