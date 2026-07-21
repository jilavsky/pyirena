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

    Monte Carlo parameters
    ----------------------
    montecarlo_n_repetitions : int
        Number of independent MC repetitions (default 10).  The mean and std
        of the per-bin volume fractions across repetitions are returned.
    montecarlo_convergence : float
        Convergence criterion: stop the MC loop when χ²/M ≤ this value
        (default 1.0 = one standard deviation per point on average).
    montecarlo_max_iter : int
        Maximum number of replacement iterations per repetition (default 100000).

    Results (populated after ``fit()``)
    ------------------------------------
    distribution : np.ndarray or None
        Normalised size distribution P(r)  [volume fraction / Å].
    distribution_std : np.ndarray or None
        Per-bin standard deviation of P(r) across Monte Carlo repetitions.
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
    # When the chi²=M discrepancy target is unreachable (error bars too small or
    # model/data mismatch at the fit-window edges), fall back to the smoothest
    # solution whose chi² is within this factor of the achievable minimum chi².
    # 1.05 → within 5 %.  Larger → smoother/more conservative fallback.
    regularization_fallback_factor: float = 1.05

    # ── TNNLS parameters ──────────────────────────────────────────────────────
    tnnls_approach_param: float = 0.95
    tnnls_max_iter: int = 1000

    # ── Monte Carlo parameters ─────────────────────────────────────────────────
    montecarlo_n_repetitions: int = 10
    montecarlo_convergence: float = 1.0
    montecarlo_max_iter: int = 100000

    # ── Error scaling ─────────────────────────────────────────────────────────
    error_scale: float = 1.0   # Multiply measurement errors by this factor before fitting
    # When True, ignore the measured/file uncertainties entirely and use a
    # fractional error err = |I| × fractional_error_value (e.g. 0.03 = 3%).
    # Useful when collected uncertainties are unreliable (e.g. after merging
    # subsets or when error estimation failed).  Mutually exclusive with
    # error_scale (error_scale is not applied in this mode).
    fractional_error: bool = False
    fractional_error_value: float = 0.03

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

        # Slit smearing.  When active, the whole G matrix is smeared once
        # (G_sm = W @ G(q_ext)) so every inversion method (MaxEnt/Reg/TNNLS/MC)
        # inherits smearing unchanged and the recovered distribution is
        # ideal-space.  The power-law background is smeared too; a flat
        # background smears to itself.  See pyirena.core.smearing.
        self.use_slit_smearing = False
        self.slit_length = 0.0

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

        # Preserve the raw observed intensities (post-mask, pre-bg-subtraction)
        # so callers can plot/store data and model on consistent arrays even
        # when the input contained non-positive or non-finite points that the
        # mask above filtered out.
        I_observed = I.copy()

        # Slit smearing operator (None => pinhole no-op).  Built once on the
        # cleaned q grid and reused for the background and the G matrix.
        smearer = None
        if self.use_slit_smearing and self.slit_length > 0:
            from pyirena.core.smearing import SlitSmearer
            smearer = SlitSmearer(q, self.slit_length)

        # Subtract complex background (B·q^(-P) + flat).  When the data are slit
        # smeared the background is smeared too before subtraction (the flat
        # part is invariant; the power-law part is not).
        if smearer is not None:
            I = I - smearer.smear_model(self.compute_complex_background)
        else:
            I = I - self.compute_complex_background(q)

        # Fractional error: replace any file/measured uncertainties with a
        # user-chosen fraction of the observed intensity.  Based on the observed
        # (pre-background-subtraction) intensity so the bars match what callers
        # plot, consistent with the ASCII-import fractional-error path.
        if self.fractional_error:
            frac = float(self.fractional_error_value)
            err = np.abs(I_observed) * frac
            err = np.maximum(err, 1e-300)
        else:
            # Generate errors if absent
            no_user_errors = err is None
            if no_user_errors:
                err = I * 0.05
                err[err <= 0] = np.abs(I[err <= 0]) * 0.05 + 1e-20

            # Apply user-specified error scaling (default 1.0 = no change)
            if self.error_scale != 1.0:
                err = err * float(self.error_scale)
                err = np.maximum(err, 1e-300)

        # Build G matrix.  When smearing, evaluate G on the extended grid and
        # apply the smearing operator to every column once; the inversion
        # machinery then sees a smeared forward model and returns an
        # ideal-space distribution.
        r_grid = make_r_grid(self.r_min, self.r_max, self.n_bins, self.log_spacing)
        if smearer is not None:
            G_ideal = self._build_g_matrix(q, r_grid)
            G = smearer.smear_columns(self._build_g_matrix(smearer.q_ext, r_grid))
        else:
            G_ideal = None
            G = self._build_g_matrix(q, r_grid)

        # Fit
        _sky_note: str | None = None
        try:
            method = self.method.lower()
            x_raw_std = None
            if method == 'maxent':
                _sky_notes: list[str] = []
                _sky_orig = float(self.maxent_sky_background)
                x_raw, n_iter = self._fit_maxent(G, I, err)

                # Layer 1: catastrophic-failure retry.
                # chi² >> M almost always means sky_background is far too large;
                # progressively reduce by 100×/1000×/10000× until chi² is sane.
                _chi2 = float(np.sum(((G @ x_raw - I) / err) ** 2))
                if _chi2 > 100.0 * len(I):
                    for _factor in (100, 1000, 10000):
                        _sky_trial = _sky_orig / _factor
                        if _sky_trial < 1e-20:
                            break
                        self.maxent_sky_background = _sky_trial
                        x_raw, n_iter = self._fit_maxent(G, I, err)
                        _chi2 = float(np.sum(((G @ x_raw - I) / err) ** 2))
                        if _chi2 <= 100.0 * len(I):
                            _sky_notes.append(
                                f"sky_background reduced {_factor}× to {_sky_trial:.2e} for convergence"
                            )
                            break

                # Layer 2: post-fit calibration.
                # sky should be ~0.01 × max(distribution); if it is more than 5×
                # above that target the fit is likely pulled toward the sky prior.
                _dw_safe = np.maximum(bin_widths(r_grid), 1e-300)
                _pos = (x_raw / _dw_safe)
                _pos = _pos[_pos > 0]
                if len(_pos) > 0:
                    _dist_max = float(np.max(_pos))
                    _sky_now = float(self.maxent_sky_background)
                    _target = 0.01 * _dist_max
                    if _sky_now > 5.0 * _target and _target > 1e-20:
                        self.maxent_sky_background = _target
                        x_raw, n_iter = self._fit_maxent(G, I, err)
                        _sky_notes.append(
                            f"sky_background calibrated {_sky_now:.2e} → {_target:.2e}"
                        )

                if _sky_notes:
                    _sky_note = '; '.join(_sky_notes)

            elif method == 'regularization':
                x_raw, n_iter = self._fit_regularization(G, I, err)
            elif method in ('tnnls', 'ipg', 'nnls'):
                x_raw, n_iter = self._fit_tnnls(G, I, err)
            elif method == 'montecarlo':
                x_raw, n_iter, x_raw_std = self._fit_montecarlo(G, I, err, r_grid, q=q)
            else:
                return self._fail(f"Unknown method '{self.method}'.")
        except Exception as exc:
            log.exception("Size distribution fit failed")
            return self._fail(str(exc))

        # Post-process
        result = self._post_process(x_raw, r_grid, G, I, err)

        # Slit-smearing provenance.  model_intensity (from _post_process) is the
        # SMEARED distribution scattering when smearing is on; also expose the
        # ideal (pinhole) curve so downstream can save/plot both.  The
        # distribution itself is ideal-space by construction.
        result['slit_length'] = float(self.slit_length) if smearer is not None else 0.0
        result['data_is_slit_smeared'] = smearer is not None
        result['model_intensity_ideal'] = (G_ideal @ x_raw) if G_ideal is not None else None

        # Monte Carlo (McSAS): the MC fit uses the standard G matrix whose columns
        # G[:,k] = V(r_k)·F²_norm·contrast·1e-4 are the intensity per unit *volume
        # fraction* of bin k.  The optimal scale A therefore makes x_raw = A·count a
        # per-bin *volume fraction* already — identical in meaning to the x_raw
        # produced by MaxEnt, Regularization, and TNNLS.  So MC flows through the
        # same _post_process (x_raw / bin_width → volume distribution) with no extra
        # weighting.  A previous version multiplied by V(r)=(4/3)πr³ here on the
        # mistaken assumption that x_raw was a number distribution; that spurious r³
        # re-weighting shifted the reported distribution to ~2× larger radius while
        # χ² (computed from the correct x_raw in _post_process) still looked good.

        # Monte Carlo: add per-bin uncertainty from spread across repetitions
        if x_raw_std is not None:
            dw_safe = np.maximum(bin_widths(r_grid), 1e-300)
            result['distribution_std'] = x_raw_std / dw_safe
        else:
            result['distribution_std'] = None

        result['n_iterations'] = n_iter
        result['n_data'] = len(I)   # number of Q points used; chi² target ≈ this value
        result['sky_note'] = _sky_note

        # Expose the actual arrays the fit used (post-mask).  Callers should
        # use these for any post-fit plotting / storing so that array shapes
        # match `model_intensity` and `residuals`.  Without this, when input
        # contains negative or non-finite intensities, the mask drops points
        # and downstream code mixing the original input arrays with fit
        # outputs hits a shape mismatch (e.g. (323,) vs (326,)).
        result['q']      = q
        result['I_data'] = I_observed   # raw observed I, post-mask, pre-bg-subtraction
        result['err']    = err
        if _sky_note:
            log.info("MaxEnt sky_background auto-adjusted: %s", _sky_note)
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

        # Slit smearing: the data are slit smeared, so the power-law model must
        # be smeared before comparison — otherwise the prefit returns
        # smeared-space B/P which the main fit would smear a second time
        # (double-smearing).  A single SlitSmearer on qf builds W once and is
        # reused across curve_fit iterations.  The flat background is invariant
        # under smearing, so only the q^-P term needs it.
        _smearer = None
        if self.use_slit_smearing and self.slit_length > 0:
            from pyirena.core.smearing import SlitSmearer
            _smearer = SlitSmearer(qf, self.slit_length)

        def _pl(q_, B, P):
            if _smearer is not None and q_ is qf:
                return _smearer.smear_model(lambda qq: B * np.power(qq, -P))
            return B * np.power(q_, -P)

        try:
            if fit_B and fit_P:
                def model(q_, B, P):
                    return _pl(q_, B, P)
                p0 = [max(self.power_law_B, float(If.mean())), self.power_law_P]
                popt, _ = curve_fit(model, qf, If, p0=p0,
                                    bounds=([0.0, 0.1], [np.inf, 12.0]),
                                    maxfev=10000)
                B_new, P_new = float(popt[0]), float(popt[1])

            elif fit_B:
                P_fixed = self.power_law_P
                def model(q_, B):
                    return _pl(q_, B, P_fixed)
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
                    return _pl(q_, B_fixed, P)
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
                # Smear the subtracted power-law when the data are slit smeared,
                # so the flat background is estimated consistently with the
                # main fit (flat itself is invariant under smearing).
                if self.use_slit_smearing and self.slit_length > 0:
                    from pyirena.core.smearing import SlitSmearer
                    _sm = SlitSmearer(qf, self.slit_length)
                    pl = _sm.smear_model(
                        lambda qq: self.power_law_B * np.power(qq, -self.power_law_P))
                else:
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

        # ── Fallback: chi² = M target not achievable ──────────────────────────
        # The minimum achievable chi² (at the smallest alpha) already exceeds M.
        # This happens when the fit-window error bars are too small and/or the
        # model cannot describe the data at the window edges (typically the noisy,
        # background-dominated high-Q points): the discrepancy principle chi²=M is
        # then unreachable no matter how little smoothing is applied.
        #
        # In that regime chi²(alpha) is nearly flat over a wide alpha range (all
        # solutions fit the data equally poorly), yet the DISTRIBUTIONS differ
        # enormously — small alpha yields a spike (over-fitting the noise), large
        # alpha yields a smooth curve.  The previous fallback returned the
        # *minimum-chi²* (smallest-alpha) solution, i.e. the spikiest member of
        # that family — the origin of the "single huge spike at the lowest bin"
        # artifact and of the extreme sensitivity to the high-Q cut-off.
        #
        # Instead, apply the discrepancy principle *relative to the achievable
        # floor*: seek the SMOOTHEST solution (largest alpha) whose chi² is still
        # within a small factor of the minimum achievable chi².  This degrades
        # gracefully and is far less sensitive to exactly where the fit window
        # ends, while remaining faithful to the (genuinely poor) data misfit.
        if c2_lo > chi_target:
            chi_fallback = self.regularization_fallback_factor * c2_lo
            log.warning(
                "Regularization: chi²=M target %.1f not achievable "
                "(minimum chi² %.1f); falling back to the smoothest solution "
                "within %.0f%% of the achievable misfit (target chi² %.1f). "
                "This usually means the fit-window error bars are too small or "
                "the model cannot describe the high-Q edge — consider trimming "
                "the high-Q end of the inversion window.",
                float(M), c2_lo,
                100.0 * (self.regularization_fallback_factor - 1.0), chi_fallback,
            )
            # chi²(alpha) is nearly flat over a wide alpha range here, so an
            # equality search would stop at an arbitrary mid-plateau alpha and
            # still return a spiky solution.  Instead find the LARGEST alpha
            # (⇒ smoothest solution) whose chi² is still ≤ chi_fallback, i.e. the
            # upper edge of the acceptable-misfit plateau.
            a_lo, a_hi = lo, hi
            for _ in range(20):          # ensure a_hi brackets from above
                if _solve(a_hi)[1] > chi_fallback:
                    break
                a_hi += 2.0
            n_iter = 0
            for n_iter in range(1, 100):
                mid = 0.5 * (a_lo + a_hi)
                if _solve(mid)[1] <= chi_fallback:
                    a_lo = mid
                else:
                    a_hi = mid
                if a_hi - a_lo < 0.01:
                    break
            x, _ = _solve(a_lo)
            if x.max() > 0:
                x = np.maximum(x, min_ratio * x.max())
            return x, n_iter

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

    # ── Monte Carlo (McSAS algorithm) ──────────────────────────────────────────

    def _fit_montecarlo(
        self,
        G: np.ndarray,
        I: np.ndarray,
        err: np.ndarray,
        r_grid: np.ndarray,
        n_reps: Optional[int] = None,
        q: Optional[np.ndarray] = None,
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
            ``self.montecarlo_n_repetitions``.  Pass ``1`` for a single-run fit
            (the ``distribution_std`` return value will be ``None``).
        q : (M,) array, optional
            The fitted Q values.  When supplied, Monte Carlo contributions are
            confined to the resolvable size band  r ∈ [π/Q_max, π/Q_min]  (the
            standard McSAS rule).  This prevents stray contributions in
            unconstrained bins from inflating the r²-weighted Rg and from
            over-fitting high-Q noise as a sub-resolution spike.  If omitted, all
            bins are eligible (legacy behaviour).

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
        n_rep = int(n_reps if n_reps is not None else self.montecarlo_n_repetitions)
        convergence = float(self.montecarlo_convergence)
        max_iter = int(self.montecarlo_max_iter)

        # ── Restrict contributions to the RESOLVABLE size band ─────────────────
        # A radius bin is only constrained by the data when its characteristic
        # scattering feature falls inside the measured Q window.  The standard
        # McSAS size-range rule is  r ∈ [π/Q_max, π/Q_min].  Bins outside this band
        # are invisible to the fit, so unconstrained: left free, Monte Carlo drops
        # stray contributions into them — a few huge-r bins wildly inflate the
        # (r²-weighted) Rg, and sub-resolution small-r bins soak up high-Q noise as
        # a spurious spike.  We forbid contributions outside the band; every other
        # method suppresses these bins implicitly (entropy prior / smoothness), so
        # this simply gives Monte Carlo the same discipline.
        allowed_mask = np.ones(N, dtype=bool)
        if q is not None and len(q) > 0:
            qv = np.asarray(q, dtype=float)
            qv = qv[np.isfinite(qv) & (qv > 0)]
            if qv.size:
                r_hi = np.pi / qv.min()      # largest resolvable radius
                r_lo = np.pi / qv.max()      # smallest resolvable radius
                band = (r_grid >= r_lo) & (r_grid <= r_hi)
                if band.sum() >= 3:          # keep a sane number of active bins
                    allowed_mask = band
        allowed_idx = np.flatnonzero(allowed_mask)
        n_allowed = allowed_idx.size

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
            # ── Initialise: assign each of N_c contributions to a random bin
            #    drawn only from the resolvable band ─────────────────────────────
            contrib_bins = allowed_idx[rng.integers(0, n_allowed, size=N_c)]
            counts = np.bincount(contrib_bins, minlength=N).astype(float)  # (N,)

            g_sum = G @ counts                               # (M,)
            A = _scale_factor(g_sum)
            chi2_val = _chi2(g_sum, A)

            # ── MC replacement loop ────────────────────────────────────────────
            rep_iters = 0
            for rep_iters in range(1, max_iter + 1):
                j = int(rng.integers(N_c))
                k_old = int(contrib_bins[j])
                k_new = int(allowed_idx[rng.integers(n_allowed)])

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
                "Monte Carlo rep %d/%d: %d iters, chi2=%.4g (target %.4g)",
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
            'montecarlo_n_repetitions': self.montecarlo_n_repetitions,
            'montecarlo_convergence':   self.montecarlo_convergence,
            'montecarlo_max_iter':      self.montecarlo_max_iter,
            'error_scale':              self.error_scale,
            'fractional_error':         self.fractional_error,
            'fractional_error_value':   self.fractional_error_value,
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


# ──────────────────────────────────────────────────────────────────────────────
# Surface-area distribution (derived from the volume distribution)
# ──────────────────────────────────────────────────────────────────────────────

def _spheroid_sv_factor(aspect_ratio: float) -> float:
    """Dimensionless surface-to-volume factor C for a spheroid: ``sv(r) = C / r``.

    The spheroid has equatorial semi-axes ``a = b = r`` and polar semi-axis
    ``c = r·AR``.  Because the surface area scales as ``r²`` and the volume as
    ``r³``, the ratio ``A(r)/V(r) = C(AR)/r`` where C depends only on AR:

        C(AR) = 1.5 · g(AR) / AR

    with ``g(AR)`` the surface-area shape factor (g = 2 for a sphere, giving
    C = 3, consistent with the sphere relation ``sv = 3/r``).

      * Prolate (AR > 1):  e = √(1 − 1/AR²),  g = 1 + (AR/e)·arcsin(e)
      * Oblate  (AR < 1):  e = √(1 − AR²),    g = 1 + (AR²/e)·artanh(e)
      * Sphere  (AR ≈ 1):  C = 3

    Args:
        aspect_ratio: polar/equatorial axis ratio (AR).  AR ≤ 0 is treated as 1.

    Returns:
        The factor C such that ``sv(r) = C / r``.
    """
    AR = float(aspect_ratio)
    if AR <= 0.0 or abs(AR - 1.0) < 1e-6:
        return 3.0
    if AR > 1.0:                       # prolate
        e = math.sqrt(1.0 - 1.0 / AR ** 2)
        g = 1.0 + (AR / e) * math.asin(e)
    else:                             # oblate
        e = math.sqrt(1.0 - AR ** 2)
        g = 1.0 + (AR ** 2 / e) * math.atanh(e)
    return 1.5 * g / AR


def particle_volume(
    r_grid: np.ndarray,
    shape: str = 'sphere',
    aspect_ratio: float = 1.0,
) -> np.ndarray:
    """Per-particle volume V(r) for the given shape [Å³].

      * Sphere:   ``V = (4/3)π r³``
      * Spheroid: ``V = (4/3)π r³ · AR``  (semi-axes r, r, r·AR)

    Args:
        r_grid:       Radius bin centres [Å].
        shape:        Particle shape (``'sphere'`` or ``'spheroid'``).  Any
                      unrecognised shape falls back to the sphere volume.
        aspect_ratio: Spheroid aspect ratio (used only when shape=='spheroid').

    Returns:
        Array of particle volumes matching ``r_grid``.
    """
    r_grid = np.asarray(r_grid, dtype=float)
    V = (4.0 / 3.0) * np.pi * r_grid ** 3
    if str(shape).lower() == 'spheroid':
        AR = float(aspect_ratio)
        if AR > 0.0:
            V = V * AR
    return V


def number_distribution(
    distribution: np.ndarray,
    r_grid: np.ndarray,
    shape: str = 'sphere',
    aspect_ratio: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Derive the number distribution from a volume distribution.

    The number distribution is ``N(r) = P_V(r) / V(r)`` with the shape-aware
    particle volume from :func:`particle_volume`.  The cumulative number
    distribution is the running integral ``∫ N(r) dr``.

    Args:
        distribution: Volume size distribution P_V(r) [vol-fraction / Å].
        r_grid:       Radius bin centres [Å].
        shape:        Particle shape (``'sphere'`` or ``'spheroid'``).
        aspect_ratio: Spheroid aspect ratio (used only when shape=='spheroid').

    Returns:
        ``(number_dist, cumul_num_dist)`` arrays matching ``r_grid``.
    """
    distribution = np.asarray(distribution, dtype=float)
    r_grid = np.asarray(r_grid, dtype=float)

    V_r = particle_volume(r_grid, shape, aspect_ratio)
    number_dist = np.where(V_r > 0, distribution / V_r, 0.0)

    dr = np.diff(r_grid, prepend=r_grid[0])
    cumul_num_dist = np.cumsum(number_dist * dr)

    return number_dist, cumul_num_dist


def surface_distribution(
    distribution: np.ndarray,
    r_grid: np.ndarray,
    shape: str = 'sphere',
    aspect_ratio: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Derive the surface-area distribution from a volume distribution.

    The surface-area distribution is ``S(r) = sv(r) · P_V(r)`` where
    ``sv(r) = A(r)/V(r)`` is the particle surface-to-volume ratio:

      * Sphere:   ``sv(r) = 3 / r``
      * Spheroid: ``sv(r) = C(AR) / r``  (see :func:`_spheroid_sv_factor`)

    The cumulative surface distribution is the running integral
    ``∫ S(r) dr``; its final value is the total specific surface area
    (surface area per unit sample volume, [Å⁻¹]) — the same normalisation
    convention used by the cumulative volume/number distributions.

    Args:
        distribution: Volume size distribution P_V(r) [vol-fraction / Å].
        r_grid:       Radius bin centres [Å].
        shape:        Particle shape (``'sphere'`` or ``'spheroid'``).  Any
                      unrecognised shape falls back to the sphere relation.
        aspect_ratio: Spheroid aspect ratio (used only when shape=='spheroid').

    Returns:
        ``(surface_dist, cumul_surf_dist, specific_surface)`` where the arrays
        match ``r_grid`` in shape and ``specific_surface`` is a float.
    """
    distribution = np.asarray(distribution, dtype=float)
    r_grid = np.asarray(r_grid, dtype=float)

    r_safe = np.maximum(r_grid, 1e-30)
    if str(shape).lower() == 'spheroid':
        factor = _spheroid_sv_factor(aspect_ratio)
    else:
        factor = 3.0
    sv = factor / r_safe

    surface_dist = sv * distribution

    dr = np.diff(r_grid, prepend=r_grid[0])
    cumul_surf_dist = np.cumsum(surface_dist * dr)
    specific_surface = float(cumul_surf_dist[-1]) if cumul_surf_dist.size else 0.0

    return surface_dist, cumul_surf_dist, specific_surface


def recommend_sizes_setup(
    q: np.ndarray,
    intensity: np.ndarray,
    sigma: Optional[np.ndarray] = None,
    config=None,
) -> dict:
    """Data-driven suggestions for a size-distribution fit + a suitability check.

    Segments the I(Q) curve (via :func:`pyirena.core.feature_detect.detect_features`)
    and derives a recommended setup: radius grid, inversion Q-range, and the
    low-Q power-law / high-Q flat background windows.  Also flags whether the
    data looks like a viable single size-distribution candidate.

    This is the single source of truth shared by the AI control tool
    (``pyirena.api.control.sizes.suggest_sizes_setup``) and the GUI's
    Feature Identifier dialog, so both show identical recommendations.

    Parameters
    ----------
    q, intensity : np.ndarray
        The scattering curve.
    sigma : np.ndarray, optional
        Intensity uncertainties (improves the segmentation noise filter).
    config : FeatureDetectConfig, optional
        Override the segmentation parameters.  Defaults to
        ``FeatureDetectConfig()`` (the same defaults the AI tool uses).

    Returns
    -------
    dict with keys ``suitable`` (bool), ``recommended`` (dict),
    ``warnings`` (list[str]), and ``features`` (dict).
    """
    from pyirena.core.feature_detect import detect_features  # noqa: PLC0415

    q = np.asarray(q, dtype=float)
    I = np.asarray(intensity, dtype=float)
    sig = np.asarray(sigma, dtype=float) if sigma is not None else None

    fr = detect_features(q, I, sigma_I=sig, config=config)
    feat = fr.to_dict()
    segments = feat.get("segments", []) or []
    knees = feat.get("guinier_knees", []) or []

    valid = np.isfinite(q) & np.isfinite(I) & (q > 0)
    if np.any(valid):
        q_lo = float(np.nanmin(q[valid]))
        q_hi = float(np.nanmax(q[valid]))
    else:
        q_lo, q_hi = float("nan"), float("nan")

    warnings: list[str] = []

    # ── Low-Q power-law window (the steep upturn that is *background*, not a
    #    particle) ────────────────────────────────────────────────────────────
    # Segments are ordered high-Q -> low-Q; kind in {background, guinier_plateau,
    # power_law}.  A broad size distribution scatters as a smoothly rolling,
    # power-law-like region, so the detector often labels the whole particle
    # region "power_law" too — we only treat the *lowest-Q* steep segment as the
    # power-law background.
    pl_window = None
    power_law_segs = [seg for seg in segments if seg.get("kind") == "power_law"]
    if power_law_segs:
        low_seg = min(power_law_segs, key=lambda seg: seg.get("q_min", q_lo))
        pl_window = (float(low_seg["q_min"]), float(low_seg["q_max"]))

    # ── High-Q flat background window ──────────────────────────────────────────
    # detect_features may hand us an explicit background segment, or (common for
    # these data) label the high-Q flat tail a "guinier_plateau".  Accept the
    # highest-Q flat-ish segment either way; otherwise fall back to the top Q decade.
    bg_q_min = feat.get("background_q_min")
    background_q_min = float(bg_q_min) if bg_q_min is not None else None
    if background_q_min is None:
        flat_kinds = [s for s in segments if s.get("kind") in ("background", "guinier_plateau")]
        hi_flat = [s for s in flat_kinds if float(s.get("q_max", 0.0)) >= 0.5 * q_hi]
        if hi_flat:
            background_q_min = float(min(hi_flat, key=lambda s: s.get("q_min", q_hi))["q_min"])
    if background_q_min is None and np.isfinite(q_hi):
        background_q_min = q_hi / (10.0 ** 0.5)      # top ~half-decade as a fallback
    background_q_max = q_hi if background_q_min is not None else None

    # ── Estimate the background so the inversion range can be chosen where the
    #    particle signal is clearly above it (the physically meaningful
    #    criterion), not from Guinier knees. ────────────────────────────────────
    flat_bg = _estimate_flat_bg(q, I, valid, background_q_min, background_q_max)
    B_pl, P_pl = _estimate_power_law(q, I, valid, pl_window)
    flat_floor = flat_bg if flat_bg > 0 else 1e-300
    with np.errstate(over="ignore", invalid="ignore"):
        bg_complex = B_pl * np.power(np.where(q > 0, q, np.nan), -P_pl) + flat_bg
    bg_complex = np.where(np.isfinite(bg_complex) & (bg_complex > 0), bg_complex, flat_floor)
    bg_flat = np.full_like(q, flat_floor, dtype=float)

    # ── Inversion (particle) Q-range from the signal-to-background ratio ────────
    # Keep the widest contiguous Q-band where I(q) ≥ sb_ratio·background, relaxing
    # the ratio if the strict band is too short (weak but real signal).  Try the
    # full complex background first — it best separates a genuine low-Q upturn —
    # then fall back to the flat background only, which is robust when the
    # power-law estimate is unreliable (e.g. it absorbed the particle signal).
    inv_q_min = inv_q_max = None
    sb_used = None
    best_band = None          # widest band seen at any ratio (last-resort fallback)
    best_span = -1.0
    best_sb = None
    for bg_curve in (bg_complex, bg_flat):
        for sb_ratio in (2.0, 1.5, 1.2, 1.05):
            band = _widest_band(q, I, bg_curve, valid, sb_ratio)
            if band is None or band[1] <= band[0]:
                continue
            span_dec = float(np.log10(band[1] / band[0]))
            if span_dec > best_span:
                best_span, best_band, best_sb = span_dec, band, sb_ratio
            if span_dec >= 0.5:          # accept only a genuinely usable band
                inv_q_min, inv_q_max = band
                sb_used = sb_ratio
                break
        if inv_q_min is not None:
            break

    # Last resort: no ratio gave ≥0.5 decades — take the widest band found.
    if inv_q_min is None and best_band is not None:
        inv_q_min, inv_q_max = best_band
        sb_used = best_sb

    have_range = inv_q_min is not None and inv_q_max is not None and inv_q_min < inv_q_max
    if not have_range:
        inv_q_min = pl_window[1] if pl_window is not None else q_lo
        inv_q_max = background_q_min if background_q_min is not None else q_hi
        if not (inv_q_min < inv_q_max):
            inv_q_min, inv_q_max = q_lo, q_hi
        warnings.append(
            "Could not isolate a Q-band where the particle signal clearly exceeds "
            "the background; defaulting the inversion range — inspect the "
            "background preview and adjust the cursors."
        )

    inv_span_dec = float(np.log10(inv_q_max / inv_q_min)) if inv_q_min and inv_q_max else 0.0

    # ── Radius grid from the inversion Q-range (nice, slightly-padded bounds) ──
    from pyirena.core.form_factors import r_bounds_from_q_range  # noqa: PLC0415
    r_min, r_max = r_bounds_from_q_range(inv_q_min, inv_q_max, pad=True)

    # ── Suitability ────────────────────────────────────────────────────────────
    # A size distribution is suitable whenever there is a usable particle-signal
    # band between the low-Q power-law upturn and the high-Q flat background —
    # which is exactly the common "broad distribution + power-law + flat bg" case
    # (precipitates, pores in rocks/minerals).  A broad distribution has NO clean
    # Guinier knee and looks like several power-law "levels", so those are NOT
    # disqualifiers — they only add advisory notes.
    log_decades = float(feat.get("log_decades", 0.0) or 0.0)
    rec_nlevels = int(feat.get("recommended_nlevels", 0) or 0)

    suitable = have_range and inv_span_dec >= 0.3

    if not have_range:
        warnings.append(
            "No clear particle-scattering band was found above the background — a "
            "size distribution may not be meaningful for this curve."
        )
    elif inv_span_dec < 0.5:
        warnings.append(
            f"The particle-signal band spans only {inv_span_dec:.1f} decades of Q — "
            "the size range that can be resolved is narrow."
        )
    if len(knees) > 1:
        warnings.append(
            f"{len(knees)} Guinier knees detected. If these are genuinely distinct "
            "particle populations, Unified Fit is an alternative; a single broad "
            "size distribution over the recommended range is still valid and is the "
            "usual choice for precipitation / porosity data."
        )
    if sb_used is not None and sb_used < 2.0:
        warnings.append(
            f"Particle signal is weak (used a signal/background ratio of {sb_used:.2f} "
            "to define the range). Background subtraction quality matters more here — "
            "check the background preview."
        )

    return {
        "suitable": suitable,
        "recommended": {
            "r_min": round(float(r_min), 3) if np.isfinite(r_min) else None,
            "r_max": round(float(r_max), 3) if np.isfinite(r_max) else None,
            "inversion_q_min": inv_q_min,
            "inversion_q_max": inv_q_max,
            "power_law_q_min": pl_window[0] if pl_window else None,
            "power_law_q_max": pl_window[1] if pl_window else None,
            "background_q_min": background_q_min,
            "background_q_max": background_q_max,
            "flat_background": float(flat_bg) if np.isfinite(flat_bg) else None,
        },
        "warnings": warnings,
        "features": {
            "n_segments_found": feat.get("n_segments_found", 0),
            "recommended_nlevels": rec_nlevels,
            "log_decades": log_decades,
            "guinier_knees": knees,
            "signal_to_bg_ratio": sb_used,
            "inversion_span_decades": round(inv_span_dec, 2),
            "segments": segments,
        },
    }


def _estimate_flat_bg(q, I, valid, bg_q_min, bg_q_max) -> float:
    """Robust flat-background level = median I over the high-Q flat window."""
    if bg_q_min is None:
        return 0.0
    m = valid & (q >= bg_q_min)
    if bg_q_max is not None:
        m = m & (q <= bg_q_max)
    vals = I[m]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return 0.0
    return float(np.median(vals))


def _estimate_power_law(q, I, valid, pl_window) -> tuple[float, float]:
    """Fit  I = B·q^(-P)  (log-log line) over the low-Q power-law window.

    Returns (B, P); (0.0, 4.0) if it cannot be estimated (→ no power-law term).
    """
    if pl_window is None:
        return 0.0, 4.0
    m = valid & (q >= pl_window[0]) & (q <= pl_window[1]) & (I > 0)
    if int(np.count_nonzero(m)) < 3:
        return 0.0, 4.0
    lq, lI = np.log(q[m]), np.log(I[m])
    try:
        slope, intercept = np.polyfit(lq, lI, 1)
    except (np.linalg.LinAlgError, ValueError):
        return 0.0, 4.0
    P = float(-slope)
    B = float(np.exp(intercept))
    if not (np.isfinite(P) and np.isfinite(B)) or P <= 0:
        return 0.0, 4.0
    return B, P


def _widest_band(q, I, bg_curve, valid, sb_ratio):
    """Return (q_lo, q_hi) of the widest contiguous run where I ≥ sb_ratio·bg."""
    m = valid & np.isfinite(I) & np.isfinite(bg_curve) & (I >= sb_ratio * bg_curve)
    if not np.any(m):
        return None
    order = np.argsort(q)
    qs, ms = q[order], m[order]
    best_lo = best_hi = None
    best_len = 0
    i = 0
    n = len(qs)
    while i < n:
        if ms[i]:
            j = i
            while j + 1 < n and ms[j + 1]:
                j += 1
            if (j - i) >= best_len:
                best_len = j - i
                best_lo, best_hi = float(qs[i]), float(qs[j])
            i = j + 1
        else:
            i += 1
    if best_lo is None:
        return None
    return best_lo, best_hi
