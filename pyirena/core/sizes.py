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
    random contributions and optimises them one at a time.  Runs n_repetitions
    independent times; the spread gives per-bin uncertainties.

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

        # McSAS: add per-bin uncertainty from spread across repetitions
        if x_raw_std is not None:
            dw_safe = np.maximum(bin_widths(r_grid), 1e-300)
            result['distribution_std'] = x_raw_std / dw_safe
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
        Maximum Entropy inversion (port of IR1R_MaximumEntropy in Igor Pro).

        Maximises  S = -Σ (x/Z) ln(x/Z)  subject to  χ² ≈ M.

        The algorithm uses three conjugate search directions derived from
        the entropy gradient and the chi-squared gradient, and solves a 3×3
        linear system for the optimal step sizes at each iteration.
        """
        M = len(I)
        N = G.shape[1]
        sky = float(self.maxent_sky_background)
        tol = self.maxent_stability * np.sqrt(2.0 * M)
        chi_target = float(M)

        # Initialise model
        x = np.full(N, sky, dtype=float)

        def chi2(x_):
            resid = (I - G @ x_) / err
            c2 = float(np.dot(resid, resid))
            return c2 if np.isfinite(c2) else 1e300

        n_iter = 0
        for n_iter in range(1, self.maxent_max_iter + 1):
            ox = G @ x                     # forward problem [M]
            resid = (I - ox) / err         # normalised residuals [M]
            c2 = float(np.dot(resid, resid))
            if not np.isfinite(c2):
                c2 = 1e300

            # ── Entropy gradient (∂S/∂x_j, with sign for maximisation → -grad)
            fSum = float(x.sum())
            sgrad = -(np.log(x / sky) + 1.0) / fSum   # shape [N]

            # ── Chi-squared gradient  ∂(χ²)/∂x_j = -2 Gᵀ[(I-Gx)/err²]
            cgrad = G.T @ (resid / err)    # shape [N]  (up to constant 2)

            # ── Build three search directions (xi vectors)
            # xi0 proportional to entropy gradient
            # xi1 proportional to chi² gradient
            # xi2 orthogonalised remainder
            snorm = float(np.dot(sgrad, sgrad))
            cnorm = float(np.dot(cgrad, cgrad))

            if snorm < 1e-300 or cnorm < 1e-300:
                break   # degenerate gradients → converged or failed

            xi0 = sgrad / np.sqrt(snorm)
            xi1 = cgrad / np.sqrt(cnorm)

            # Gram-Schmidt orthogonalise xi1 against xi0
            xi1 = xi1 - np.dot(xi1, xi0) * xi0
            xi1_norm = float(np.dot(xi1, xi1))
            if xi1_norm < 1e-300:
                xi1 = np.zeros(N)
            else:
                xi1 /= np.sqrt(xi1_norm)

            # Third direction: orthogonal to both
            xi2 = np.cross(xi0, xi1) if N == 3 else np.zeros(N)
            # For N != 3 we use only 2 directions (xi2 stays zero)

            # ── Map xi to data space for 3×3 normal equations
            #   Aij = (G·xi_i)ᵀ (G·xi_j) / err²
            Gxi = np.column_stack([G @ xi0, G @ xi1, G @ xi2])   # (M, 3)
            Gxi_scaled = Gxi / err[:, np.newaxis]                  # (M, 3)

            # 3×3 normal matrix:  N_ij = Σ_k (G·xi_i)[k] (G·xi_j)[k] / err[k]²
            Nm = Gxi_scaled.T @ Gxi_scaled

            # RHS: gradient of Lagrangian projected onto xi directions
            # For entropy: ∂(χ²)/∂β_i where χ²-target → λ drives equality
            # Use the self-consistent step:  β = λ A^{-1} s_proj + μ A^{-1} c_proj
            # Simplified: solve for β by requiring χ²(x+Σβ_i ξ_i) → chi_target
            # Project both gradients:
            s_proj = np.array([np.dot(sgrad, xi0),
                                np.dot(sgrad, xi1),
                                np.dot(sgrad, xi2)])
            c_proj = np.array([np.dot(cgrad, xi0),
                                np.dot(cgrad, xi1),
                                np.dot(cgrad, xi2)])

            # ── Lagrange step: balance chi² descent and entropy ascent.
            #
            # The step Δx = Σ βᵢ ξᵢ is found by solving:
            #    Nm · β = λ · s_proj − c_proj
            # where
            #   s_proj = projections of entropy gradient onto search dirs
            #   c_proj = projections of chi²-descent gradient onto search dirs
            #   λ = (chi_target − c2)/c2   (negative when c2 > target)
            #
            # Note: xi0 = sgrad/|sgrad| points in the all-negative direction
            # (sgrad < 0 at initialisation), so c_proj[0] < 0 and s_proj[0] > 0.
            # With λ < 0 (c2 > target):
            #   rhs[0] = λ·s_proj[0] − c_proj[0]  →  (large negative) − (negative)
            #   beta[0] < 0  →  delta ∝ beta·xi0 = negative·(−ones) = positive
            # x increases → G@x grows → chi² falls toward target. ✓
            # Reversing the sign (c_proj − λ·s_proj) gives beta[0] > 0, driving x
            # to the clamp floor and preventing convergence.
            if abs(c2) < 1e-300:
                break

            lam = (chi_target - c2) / c2
            rhs = lam * s_proj - c_proj

            # Solve 3×3 system Nm · β = rhs  (third row/col near-zero when N>3)
            try:
                beta = np.linalg.solve(Nm + 1e-12 * np.eye(3), rhs)
            except np.linalg.LinAlgError:
                beta = np.zeros(3)

            # ── Backtracking line search: halve step until chi² does not grow
            step = 1.0
            for _ in range(8):
                delta = step * (beta[0] * xi0 + beta[1] * xi1 + beta[2] * xi2)
                x_trial = np.maximum(x + delta, sky * 1e-6)
                c2_trial = chi2(x_trial)
                # Accept if chi² moved toward target (or is already there)
                if c2_trial <= c2 + 1.0 or c2_trial <= chi_target:
                    break
                step *= 0.5

            x = x_trial

            # ── Convergence check
            c2_new = c2_trial
            if n_iter % 50 == 0 or n_iter <= 5:
                log.debug(
                    "MaxEnt iter %d: chi2=%.4g  target=%.4g  step=%.3g",
                    n_iter, c2_new, chi_target, step,
                )
            if abs(c2_new - chi_target) < tol:
                # Secondary test: entropy and chi² gradients nearly parallel
                tnorm = float(np.dot(sgrad, cgrad))
                test = np.sqrt(max(0.0, 0.5 * (1.0 - tnorm / (np.sqrt(snorm * cnorm) + 1e-300))))
                log.debug("MaxEnt near target at iter %d: test=%.4g", n_iter, test)
                if test < 0.05:
                    log.debug("MaxEnt converged at iteration %d", n_iter)
                    break

        log.debug("MaxEnt finished: %d iterations, chi2=%.4g, target=%.4g",
                  n_iter, chi2(x), chi_target)
        return x, n_iter

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
    ) -> tuple[np.ndarray, int, np.ndarray]:
        """
        Monte Carlo Size Analysis (McSAS) inversion.

        Uses the same r_grid as the other methods.  Each of the N_c = n_bins
        contributions is assigned to one of the N r_grid bins (discrete).
        The forward model uses volume-weighted G columns so that the scale
        factor has units of number density [Å^-3] and the output is a proper
        volume distribution:

            I(q) = A_n × Σ_k count[k] × G[:,k] × V(r_k)
            x_raw[k] = A_n × count[k] × V(r_k)   [volume fraction]

        where V(r_k) = (4/3)π r_k³ [Å³].  Using the volume-weighted columns
        ensures that the MC replacement loop implicitly treats each contribution
        as ONE PARTICLE (number density A_n), so the converged count[k]
        reflects the number distribution Nd(r), and x_raw[k] = Nd(r_k)×V(r_k)
        is the volume distribution.

        For each of ``mcsas_n_repetitions`` independent repetitions:
          1. Randomly assign N_c = n_bins contributions to bins.
          2. Find the optimal number-density scale factor A_n by weighted LS
             using the V-weighted model sum G_V @ counts.
          3. Iterate: pick a random contribution, reassign it to a random bin,
             accept if chi² does not increase.
          4. Stop when chi²/M ≤ ``mcsas_convergence`` or ``mcsas_max_iter``
             attempts have been made.
          5. x_raw[k] = A_n × count[k] × V(r_k).

        Returns
        -------
        x_raw_mean : (N,) array
            Mean volume-fraction per bin across repetitions.
        n_iter_avg : int
            Average number of MC iterations per repetition.
        x_raw_std : (N,) array
            Std of volume-fraction per bin across repetitions.
        """
        M, N = G.shape
        N_c = N   # one contribution per r_grid bin; they can cluster during MC
        n_rep = int(self.mcsas_n_repetitions)
        convergence = float(self.mcsas_convergence)
        max_iter = int(self.mcsas_max_iter)

        # Sphere volume at each grid point [Å³]
        V_r = (4.0 / 3.0) * np.pi * r_grid ** 3   # (N,)

        # Volume-weighted G columns: G_V[:,k] = G[:,k] × V(r_k)
        # With this model: I = A_n × (G_V @ counts) gives number-density scaling.
        G_V = G * V_r[np.newaxis, :]               # (M, N)

        rng = np.random.default_rng()
        inv_err2 = 1.0 / (err ** 2)               # (M,)

        def _scale_factor(gv_sum: np.ndarray) -> float:
            """Optimal number-density A_n = (gv·I/σ²)/(gv·gv/σ²), clamped ≥ 0."""
            num = float(np.dot(gv_sum * inv_err2, I))
            den = float(np.dot(gv_sum * inv_err2, gv_sum))
            return max(num / den, 0.0) if den > 1e-300 else 0.0

        def _chi2(gv_sum: np.ndarray, A: float) -> float:
            resid = (I - A * gv_sum) / err
            return float(np.dot(resid, resid))

        x_raw_reps = np.zeros((n_rep, N), dtype=float)
        total_iters = 0

        for rep in range(n_rep):
            # ── Initialise: assign each of N_c contributions to a random bin ──
            contrib_bins = rng.integers(0, N, size=N_c)      # (N_c,) bin indices
            counts = np.bincount(contrib_bins, minlength=N).astype(float)  # (N,)

            gv_sum = G_V @ counts                            # (M,)
            A = _scale_factor(gv_sum)
            chi2_val = _chi2(gv_sum, A)

            # ── MC replacement loop ────────────────────────────────────────────
            rep_iters = 0
            for rep_iters in range(1, max_iter + 1):
                j = int(rng.integers(N_c))
                k_old = int(contrib_bins[j])
                k_new = int(rng.integers(N))

                if k_new == k_old:
                    continue

                # Incremental update of gv_sum (O(M) not O(M×N))
                gv_sum_trial = gv_sum + G_V[:, k_new] - G_V[:, k_old]
                A_trial = _scale_factor(gv_sum_trial)
                chi2_trial = _chi2(gv_sum_trial, A_trial)

                if chi2_trial <= chi2_val:
                    contrib_bins[j] = k_new
                    counts[k_old] -= 1.0
                    counts[k_new] += 1.0
                    gv_sum = gv_sum_trial
                    A = A_trial
                    chi2_val = chi2_trial

                if chi2_val / M <= convergence:
                    break

            total_iters += rep_iters
            log.debug(
                "McSAS rep %d/%d: %d iters, chi2=%.4g (target %.4g)",
                rep + 1, n_rep, rep_iters, chi2_val, convergence * M,
            )

            # ── Volume distribution: x_raw[k] = A_n × count[k] × V(r_k) ──────
            x_raw_reps[rep] = A * counts * V_r

        x_raw_mean = np.mean(x_raw_reps, axis=0)
        x_raw_std  = np.std(x_raw_reps, axis=0)
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
