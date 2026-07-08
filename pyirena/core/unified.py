"""
Unified Fit Model for Small-Angle Scattering Data Analysis

This module implements the Unified fit model based on work by Greg Beaucage
for analyzing hierarchical structures using SAXS/SANS/USAXS data.

References:
- http://www.eng.uc.edu/~gbeaucag/PDFPapers/Beaucage2.pdf
- http://www.eng.uc.edu/~gbeaucag/PDFPapers/Beaucage1.pdf
- http://www.eng.uc.edu/~gbeaucag/PDFPapers/ma970373t.pdf

The model combines multiple structural levels, each described by:
    I_i(q) = G_i * exp(-q² R_g,i² / 3) + exp(-q² R_g(i-1)² / 3) * B_i * {[erf(k*q*R_g,i / √6)]³ / q}^P_i

With optional correlations using Born-Green approximation.
"""

import numpy as np
from scipy.special import erf, gamma
from scipy.optimize import least_squares
from dataclasses import dataclass
from typing import List, Optional, Dict, Tuple
import warnings 


@dataclass
class UnifiedLevel:
    """
    Parameters for one structural level in the Unified fit model.

    Attributes:
        Rg: Radius of gyration [Angstroms]
        G: Guinier prefactor (low-q amplitude) [cm^-1]
        P: Power law slope [dimensionless]
        B: Power law prefactor (Porod constant) [cm^-1 Angstrom^-P]
        ETA: Correlation distance [Angstroms]
        PACK: Packing factor [dimensionless]
        RgCO: Cutoff radius [Angstroms]
        K: Correction factor (1.0 or 1.06) [dimensionless]
        correlations: Enable correlation function
        mass_fractal: Enable mass fractal mode
        link_RGCO: Link RgCO to previous level's Rg
        link_B: Estimate B from G, Rg, and P

    Fitting control:
        fit_Rg, fit_G, fit_P, fit_B, fit_ETA, fit_PACK, fit_RgCO: Enable fitting for each parameter
        Rg_limits, G_limits, P_limits, B_limits, etc.: (min, max) bounds for fitting
    """
    # Core parameters
    Rg: float = 10.0
    G: float = 1.0
    P: float = 4.0
    B: float = 1e-4
    ETA: float = 10.0
    PACK: float = 0.0
    RgCO: float = 0.0
    K: float = 1.0

    # Boolean flags
    correlations: bool = False
    mass_fractal: bool = False
    link_RGCO: bool = False
    link_B: bool = False

    # Fitting control
    fit_Rg: bool = True
    fit_G: bool = True
    fit_P: bool = False
    fit_B: bool = True
    fit_ETA: bool = False
    fit_PACK: bool = False
    fit_RgCO: bool = False

    # Parameter bounds (min, max)
    Rg_limits: Tuple[float, float] = (0.1, 1e6)
    G_limits: Tuple[float, float] = (1e-10, 1e10)
    P_limits: Tuple[float, float] = (0.0, 6.0)
    B_limits: Tuple[float, float] = (1e-20, 1e10)
    ETA_limits: Tuple[float, float] = (0.1, 1e6)
    PACK_limits: Tuple[float, float] = (0.0, 16.0)
    RgCO_limits: Tuple[float, float] = (0.0, 1e6)

    def auto_calculate_K(self):
        """Automatically set K based on P value (version 2.36 logic)."""
        self.K = 1.0 if self.P > 3 else 1.06

    def auto_calculate_B(self):
        """Calculate B from G, Rg, and P using Hammouda relationship."""
        if self.link_B and self.Rg > 0:
            self.B = self.G * np.exp(-self.P / 2.0) * ((3.0 * self.P / 2.0) ** (self.P / 2.0)) * (1.0 / self.Rg ** self.P)

    def auto_calculate_B_mass_fractal(self):
        """Calculate B for mass fractal mode."""
        if self.mass_fractal and self.Rg > 0:
            self.B = (self.G * self.P / self.Rg ** self.P) * np.exp(np.log(gamma(self.P / 2.0)))

    def check_physical_feasibility(self) -> bool:
        """
        Check if the level parameters represent a physically meaningful set.

        A level is feasible if the Guinier and power-law regions smoothly connect
        at the rollover Q point (defined by Hammouda's theory). The continuity is
        checked by comparing their values at the rollover: a large discontinuity
        indicates unphysical parameter combinations.

        Returns:
            True if physically feasible, False otherwise.
        """
        # Special case: removed level (G=0 and Rg very large)
        if self.G <= 0 and self.Rg > 1e9:
            return True

        # Rollover Q value (Hammouda definition)
        Q_rollover = 2.0 * (1.0 / self.Rg) * np.sqrt(self.P / 2.0)

        # Guinier value at rollover Q
        guinier_value = self.G * np.exp((-Q_rollover**2 * self.Rg**2) / 3.0)

        # Power law value at rollover Q
        power_law_value = self.B / (Q_rollover**self.P)

        # Relative difference (fraction of Guinier value)
        if guinier_value <= 0:
            return False
        difference = (power_law_value - guinier_value) / guinier_value

        # Valid range: parameters must produce smooth connection
        LOW_LIMIT = -0.633
        HIGH_LIMIT = 2.25

        return LOW_LIMIT <= difference <= HIGH_LIMIT


class UnifiedFitModel:
    """
    Complete Unified fit model supporting multiple structural levels.
    """

    def __init__(self, num_levels: int = 1):
        """
        Initialize the Unified fit model.

        Args:
            num_levels: Number of structural levels (1-5)
        """
        if not 1 <= num_levels <= 5:
            raise ValueError("Number of levels must be between 1 and 5")

        self.num_levels = num_levels
        self.levels: List[UnifiedLevel] = [UnifiedLevel() for _ in range(num_levels)]
        self.background: float = 0.0
        self.fit_background: bool = True
        self.background_limits: Tuple[float, float] = (0.0, 1e10)

        # Data storage
        self.q_data: Optional[np.ndarray] = None
        self.I_data: Optional[np.ndarray] = None
        self.error_data: Optional[np.ndarray] = None
        self.dq_data: Optional[np.ndarray] = None

        # Slit smearing parameters
        self.use_slit_smearing: bool = False
        self.slit_length: float = 0.0

        # Fit results
        self.fit_result = None
        self.fit_intensity: Optional[np.ndarray] = None
        self.chi_squared: Optional[float] = None
        self.reduced_chi_squared: Optional[float] = None

    def sphere_amplitude(self, q: np.ndarray, eta: float) -> np.ndarray:
        """
        Calculate sphere amplitude function for Born-Green correlation.

        f(q, eta) = 3 * [sin(q*eta) - q*eta*cos(q*eta)] / (q*eta)^3

        Args:
            q: Scattering vector [1/Angstrom]
            eta: Correlation distance [Angstroms]

        Returns:
            Sphere amplitude values
        """
        q_eta = q * eta
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Handle q*eta → 0 limit (result → 1)
            result = np.where(
                q_eta < 1e-10,
                1.0,
                3.0 * (np.sin(q_eta) - q_eta * np.cos(q_eta)) / (q_eta ** 3)
            )
        return result

    def calculate_level_intensity(self, q: np.ndarray, level_idx: int,
                                  prev_Rg: float = 0.0) -> np.ndarray:
        """
        Calculate intensity contribution from one structural level.

        Args:
            q: Scattering vector [1/Angstrom]
            level_idx: Level index (0-based)
            prev_Rg: Rg from previous level (for RgCO linking)

        Returns:
            Intensity from this level [cm^-1]
        """
        level = self.levels[level_idx]

        # Auto-calculate parameters if needed
        level.auto_calculate_K()

        # Handle RgCO linking
        RgCO = level.RgCO
        if level.link_RGCO and level_idx > 0:
            RgCO = prev_Rg

        # Handle mass fractal mode
        if level.mass_fractal:
            level.auto_calculate_B_mass_fractal()
        elif level.link_B:
            level.auto_calculate_B()

        # Calculate Q* (characteristic Q)
        k_factor = level.K
        sqrt_6 = np.sqrt(6.0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            erf_term = erf(k_factor * q * level.Rg / sqrt_6)
            # Avoid division by zero
            erf_term = np.maximum(erf_term, 1e-10)
            Q_star = q / (erf_term ** 3)

        # Guinier term
        guinier = level.G * np.exp(-q ** 2 * level.Rg ** 2 / 3.0)

        # Power law term with cutoff
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cutoff_exp = np.exp(-RgCO ** 2 * q ** 2 / 3.0)
            # Avoid division by zero
            Q_star_safe = np.maximum(Q_star, 1e-100)
            power_law = (level.B / (Q_star_safe ** level.P)) * cutoff_exp

        # Combine terms
        intensity = guinier + power_law

        # Apply correlation function if enabled
        if level.correlations and level.PACK > 0:
            sphere_amp = self.sphere_amplitude(q, level.ETA)
            intensity = intensity / (1.0 + level.PACK * sphere_amp)

        return intensity

    def calculate_intensity(self, q: np.ndarray) -> np.ndarray:
        """
        Calculate total intensity from all levels plus background.

        Args:
            q: Scattering vector [1/Angstrom]

        Returns:
            Total intensity [cm^-1]
        """
        intensity = np.zeros_like(q, dtype=float)

        prev_Rg = 0.0
        for i in range(self.num_levels):
            intensity += self.calculate_level_intensity(q, i, prev_Rg)
            prev_Rg = self.levels[i].Rg

        # Add background
        intensity += self.background

        return intensity

    def apply_slit_smearing(self, q: np.ndarray, intensity: np.ndarray) -> np.ndarray:
        """
        Apply slit smearing correction to the calculated intensity.

        This is a simplified implementation. For production use, implement
        proper slit smearing integration over the slit length.

        Args:
            q: Scattering vector [1/Angstrom]
            intensity: Calculated intensity [cm^-1]

        Returns:
            Slit-smeared intensity [cm^-1]
        """
        if not self.use_slit_smearing or self.slit_length <= 0:
            return intensity

        # Simplified slit smearing: convolve with rectangular slit function
        # For proper implementation, see Lake (1967) or Strobl (1970)
        warnings.warn("Slit smearing is simplified. Use proper integration for production.")

        # Check if q range is sufficient
        if q[-1] < 3 * self.slit_length:
            warnings.warn(f"Q range should extend to at least 3*slit_length = {3*self.slit_length:.4f}")

        return intensity  # Placeholder - implement proper smearing if needed

    def _pack_parameters(self) -> np.ndarray:
        """Pack all fittable parameters into a 1D array."""
        params = []

        for level in self.levels:
            if level.fit_Rg:
                params.append(level.Rg)
            if level.fit_G:
                params.append(level.G)
            if level.fit_P:
                params.append(level.P)
            if level.fit_B and not level.link_B and not level.mass_fractal:
                params.append(level.B)
            if level.fit_ETA:
                params.append(level.ETA)
            if level.fit_PACK:
                params.append(level.PACK)
            if level.fit_RgCO and not level.link_RGCO:
                params.append(level.RgCO)

        if self.fit_background:
            params.append(self.background)

        return np.array(params)

    def _unpack_parameters(self, params: np.ndarray):
        """Unpack 1D parameter array back into level objects."""
        idx = 0

        for level in self.levels:
            if level.fit_Rg:
                level.Rg = params[idx]
                idx += 1
            if level.fit_G:
                level.G = params[idx]
                idx += 1
            if level.fit_P:
                level.P = params[idx]
                idx += 1
            if level.fit_B and not level.link_B and not level.mass_fractal:
                level.B = params[idx]
                idx += 1
            if level.fit_ETA:
                level.ETA = params[idx]
                idx += 1
            if level.fit_PACK:
                level.PACK = params[idx]
                idx += 1
            if level.fit_RgCO and not level.link_RGCO:
                level.RgCO = params[idx]
                idx += 1

        if self.fit_background:
            self.background = params[idx]

    def _get_bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get parameter bounds for fitting."""
        lower = []
        upper = []

        for level in self.levels:
            if level.fit_Rg:
                lower.append(level.Rg_limits[0])
                upper.append(level.Rg_limits[1])
            if level.fit_G:
                lower.append(level.G_limits[0])
                upper.append(level.G_limits[1])
            if level.fit_P:
                lower.append(level.P_limits[0])
                upper.append(level.P_limits[1])
            if level.fit_B and not level.link_B and not level.mass_fractal:
                lower.append(level.B_limits[0])
                upper.append(level.B_limits[1])
            if level.fit_ETA:
                lower.append(level.ETA_limits[0])
                upper.append(level.ETA_limits[1])
            if level.fit_PACK:
                lower.append(level.PACK_limits[0])
                upper.append(level.PACK_limits[1])
            if level.fit_RgCO and not level.link_RGCO:
                lower.append(level.RgCO_limits[0])
                upper.append(level.RgCO_limits[1])

        if self.fit_background:
            lower.append(self.background_limits[0])
            upper.append(self.background_limits[1])

        return np.array(lower), np.array(upper)

    def _residuals(self, params: np.ndarray) -> np.ndarray:
        """Calculate weighted residuals for fitting."""
        self._unpack_parameters(params)

        model_intensity = self.calculate_intensity(self.q_data)

        if self.use_slit_smearing:
            model_intensity = self.apply_slit_smearing(self.q_data, model_intensity)

        # Calculate weighted residuals
        if self.error_data is not None and np.any(self.error_data > 0):
            weights = 1.0 / self.error_data
        else:
            weights = np.ones_like(self.I_data)

        residuals = (self.I_data - model_intensity) * weights

        return residuals

    def fit(self, q: np.ndarray, intensity: np.ndarray,
            error: Optional[np.ndarray] = None,
            method: str = 'trf',
            max_iterations: int = 500,
            verbose: int = 0) -> Dict:
        """
        Fit the Unified model to experimental data.

        Args:
            q: Scattering vector [1/Angstrom]
            intensity: Measured intensity [cm^-1]
            error: Uncertainty in intensity [cm^-1]
            method: Optimization method ('trf', 'dogbox', 'lm')
            max_iterations: Maximum number of iterations
            verbose: Verbosity level (0, 1, 2)

        Returns:
            Dictionary with fit results
        """
        self.q_data = np.asarray(q)
        self.I_data = np.asarray(intensity)
        self.error_data = np.asarray(error) if error is not None else None

        # Initial parameters
        p0 = self._pack_parameters()

        # Get bounds
        lower, upper = self._get_bounds()

        # Perform fit
        # x_scale='jac' auto-rescales each parameter by its Jacobian-column
        # norm every iteration. Unified-fit parameters span many decades
        # (G ~ 10³, B ~ 10⁻⁴, Rg ~ 10¹, P ~ 10⁰, background ~ 10⁻²); without
        # this the TRF trust region and xtol test — which act on the raw
        # parameter vector — cannot take a step that is meaningful for both a
        # large and a tiny parameter at once, so the fit terminates far short
        # of the minimum. Auto-scaling gives every parameter a natural unit
        # step (the scipy equivalent of Igor's per-parameter fit-step / epsilon
        # on log-dependent parameters). 'lm' does not support 'jac'-scaling, so
        # only pass it for the bounded methods.
        #
        # Tight ftol/xtol/gtol: with 'jac' scaling the convergence tests run in
        # SCALED space, where scipy's default 1e-8 is far too loose — the xtol
        # test fires on the first small step and the solver returns "converged"
        # while still far from the minimum. This was the cause of the
        # "press Fit 2–3 times to converge" behaviour. 1e-12 defers termination
        # to a genuine minimum.
        # Tolerance kwargs are only meaningful (and 'x_scale' only accepted) for
        # the bounded methods; 'lm' keeps scipy defaults.
        tol_kwargs: Dict = {}
        if method in ('trf', 'dogbox'):
            tol_kwargs = dict(x_scale='jac', ftol=1e-12, xtol=1e-12, gtol=1e-12)

        # Internal restart loop: re-seed the solver from its own result until
        # χ² stops improving. This is the automated equivalent of a user
        # pressing "Fit" a few times — a safety net for the rare stiff case
        # where a single least_squares call still terminates one step short,
        # so scripts (which cannot re-press) always get the fully-settled
        # result on the first call. Converges in 1–3 restarts (~50–70 total
        # function evaluations); a fully-converged first call exits after one.
        x_current = np.asarray(p0, dtype=float)
        prev_chi2 = np.inf
        max_restarts = 5
        result = least_squares(
            self._residuals, x_current,
            bounds=(lower, upper), method=method,
            max_nfev=max_iterations, verbose=verbose, **tol_kwargs,
        )
        for _ in range(max_restarts - 1):
            chi2 = float(np.sum(result.fun ** 2))
            if prev_chi2 - chi2 <= 1e-8 * max(chi2, 1.0):
                break
            prev_chi2 = chi2
            result = least_squares(
                self._residuals, result.x,
                bounds=(lower, upper), method=method,
                max_nfev=max_iterations, verbose=verbose, **tol_kwargs,
            )
        self.fit_result = result

        # Unpack final parameters
        self._unpack_parameters(self.fit_result.x)

        # Calculate final fit
        self.fit_intensity = self.calculate_intensity(self.q_data)
        if self.use_slit_smearing:
            self.fit_intensity = self.apply_slit_smearing(self.q_data, self.fit_intensity)

        # Calculate chi-squared
        residuals = self._residuals(self.fit_result.x)
        self.chi_squared = np.sum(residuals ** 2)

        # Degrees of freedom
        n_data = len(self.q_data)
        n_params = len(p0)
        dof = n_data - n_params

        self.reduced_chi_squared = self.chi_squared / dof if dof > 0 else np.inf

        # 'success' means the fit ran to completion and returned finite, valid
        # parameters — NOT scipy's internal convergence status. scipy marks a
        # result as success=False when the last restart hit max_nfev or a loose
        # tolerance, even when chi2 is fully converged and the parameters are
        # physically meaningful (e.g. chi2=5000 on noisy/complex data is a
        # valid, usable result that the GUI correctly accepts as "success").
        # Tying success to scipy's flag causes batch scripts to report the same
        # fit as "failed" that the GUI shows as "Fit completed successfully!".
        fit_succeeded = np.isfinite(self.chi_squared) and len(self.levels) > 0
        results = {
            'success': fit_succeeded,
            'message': self.fit_result.message,
            'chi_squared': self.chi_squared,
            'reduced_chi_squared': self.reduced_chi_squared,
            'n_iterations': self.fit_result.nfev,
            'levels': self.levels,
            'background': self.background,
            'fit_intensity': self.fit_intensity,
            'residuals': self.I_data - self.fit_intensity
        }

        return results

    def calculate_invariant(self, level_idx: int = 0,
                           q_min: float = 0.0,
                           q_max: Optional[float] = None) -> Dict[str, float]:
        """
        Calculate the scattering invariant Q for a given level.

        Q = ∫₀^∞ I(q) q² dq

        Note:
            This integrates over an explicit [q_min, q_max] range (default:
            the data range).  For the Igor-equivalent calculation used by the
            GUI and batch API — fixed Qmax = 2*pi/(Rg/10) with Simpson
            integration — use :func:`compute_invariant_sv` instead.

        Args:
            level_idx: Level index (0-based)
            q_min: Minimum q for integration [1/Angstrom]
            q_max: Maximum q for integration [1/Angstrom]

        Returns:
            Dictionary with invariant and related quantities
        """
        level = self.levels[level_idx]

        if q_max is None:
            if self.q_data is not None:
                q_max = self.q_data[-1]
            else:
                q_max = 1.0  # Default

        # Create q array for integration
        q = np.linspace(q_min, q_max, 1000)

        # Calculate intensity
        intensity = self.calculate_level_intensity(q, level_idx)

        # Integrate I(q) * q^2
        # At q=0 the power-law term blows up to +inf (the erf denominator is
        # clamped to a constant instead of vanishing like q), so the bin
        # evaluates to inf*0 = nan and poisons the integral.  Its true
        # contribution to I(q)*q^2 there is 0, so zero out any non-finite bin
        # (and silence the expected inf*0 warning from that single bin).
        with np.errstate(invalid='ignore', over='ignore'):
            integrand = intensity * q ** 2
        integrand = np.where(np.isfinite(integrand), integrand, 0.0)
        invariant = np.trapezoid(integrand, q)

        # Add analytic Porod tail beyond q_max when there is no cutoff.
        # Skip when P ~ 3 (the (3-P) denominator is singular there).
        if level.RgCO < 0.1:
            denom = 3.0 - abs(level.P)
            if abs(denom) > 1e-6:
                invariant += -level.B * q_max ** denom / denom

        # Convert from cm^-1 Angstrom^-3 to cm^-4
        invariant_cm4 = invariant * 1e24

        # Surface/volume ratio, only meaningful in the Porod regime (P ~ 4).
        # 1e4 converts pi*B/Q from A^-1 to m^2/cm^3 (matches Igor Irena).
        surf_to_vol = None
        if 3.95 <= level.P <= 4.05 and invariant > 0:
            surf_to_vol = 1e4 * np.pi * level.B / invariant  # m^2/cm^3

        return {
            'invariant': invariant,
            'invariant_cm4': invariant_cm4,
            'surface_to_volume': surf_to_vol,
            'q_min': q_min,
            'q_max': q_max
        }

    def get_parameter_summary(self) -> str:
        """Generate a formatted summary of all parameters."""
        summary = []
        summary.append("=" * 70)
        summary.append("UNIFIED FIT MODEL - PARAMETER SUMMARY")
        summary.append("=" * 70)

        for i, level in enumerate(self.levels, 1):
            summary.append(f"\nLevel {i}:")
            summary.append(f"  Rg     = {level.Rg:.4e} Å")
            summary.append(f"  G      = {level.G:.4e} cm⁻¹")
            summary.append(f"  P      = {level.P:.4f}")
            summary.append(f"  B      = {level.B:.4e} cm⁻¹ Å⁻ᴾ")
            summary.append(f"  RgCO   = {level.RgCO:.4e} Å")
            summary.append(f"  K      = {level.K:.4f}")

            if level.correlations:
                summary.append(f"  ETA    = {level.ETA:.4e} Å")
                summary.append(f"  PACK   = {level.PACK:.4f}")

            flags = []
            if level.mass_fractal:
                flags.append("MassFractal")
            if level.link_B:
                flags.append("LinkB")
            if level.link_RGCO:
                flags.append("LinkRGCO")
            if level.correlations:
                flags.append("Correlations")

            if flags:
                summary.append(f"  Flags: {', '.join(flags)}")

        summary.append(f"\nBackground: {self.background:.4e} cm⁻¹")

        if self.chi_squared is not None:
            summary.append("\nFit Quality:")
            summary.append(f"  χ² = {self.chi_squared:.4e}")
            summary.append(f"  Reduced χ² = {self.reduced_chi_squared:.4e}")

        summary.append("=" * 70)

        return "\n".join(summary)


def compute_invariant_sv(
    G: float,
    Rg: float,
    B: float,
    P: float,
    RgCO: float = 0.0,
    ETA: float = 0.0,
    PACK: float = 0.0,
    correlated: bool = False,
) -> Tuple[Optional[float], Optional[float]]:
    """
    Compute the scattering invariant and surface-to-volume ratio for one
    Unified Fit level.

    This is the Igor-faithful implementation (port of
    ``IR1A_SurfToVolCalcInvarVec`` / ``IR1A_UpdatePorodSfcandInvariant``)
    and the single source of truth used by both the GUI panel and the
    batch API, so both report identical numbers.

    Method
    ------
    - Evaluate the single-level unified intensity on 2000 points from
      Q = 0 to Qmax = 2*pi / (Rg/10).
    - Invariant = Simpson integral of I(Q)*Q^2, plus (when RgCO < 0.1 and
      P is not ~3) the analytic Porod tail  -B * Qmax^(3-P) / (3-P).
    - Sv = 1e4 * pi * B / invariant  [m^2/cm^3], reported only in the
      Porod regime 3.95 <= P <= 4.05.

    Returns
    -------
    (invariant_cm4, sv)
        invariant in cm^-4 (converted from cm^-1 A^-3 via 1e24), or None
        when not computable (Rg <= 0, or a negative invariant indicating a
        bad extrapolation).  sv is None outside the Porod regime.
    """
    from scipy.integrate import simpson

    if Rg <= 0:
        return None, None
    try:
        maxQ = 2 * np.pi / (Rg / 10)
        surf_q = np.linspace(0, maxQ, 2000)

        # Single-level unified intensity, Igor-style q -> 0 regularization
        K = 1.0 if P > 3 else 1.06
        q_safe = np.where(np.abs(surf_q) < 1e-10, 1e-10, surf_q)
        erf_cubed = erf(K * q_safe * Rg / np.sqrt(6)) ** 3
        erf_cubed = np.where(np.abs(erf_cubed) < 1e-10, 1e-10, erf_cubed)
        qstar = q_safe / erf_cubed
        qstar = np.where(np.isfinite(qstar), qstar, 1e-10)
        intensity = (G * np.exp(-surf_q ** 2 * Rg ** 2 / 3)
                     + (B / qstar ** P) * np.exp(-RgCO ** 2 * surf_q ** 2 / 3))
        intensity = np.where(np.isfinite(intensity), intensity, 0.0)
        if correlated and PACK > 0 and ETA > 0:
            qr = np.where(surf_q * ETA == 0, 1e-10, surf_q * ETA)
            sphere_amp = 3 * (np.sin(qr) - qr * np.cos(qr)) / qr ** 3
            intensity = intensity / (1 + PACK * sphere_amp)
        # Handle the Q=0 point (use the Q=dQ value)
        if intensity[0] == 0 or np.isnan(intensity[0]):
            intensity[0] = intensity[1]

        # Invariant = integral of I(Q)*Q^2 dQ
        invariant = simpson(intensity * surf_q ** 2, x=surf_q)

        # Analytic Porod tail beyond Qmax when there is no cutoff.
        # Skip when P ~ 3: the integrand becomes B/Q which integrates to
        # B*ln(Q), and the closed-form (3-P) denominator is singular.
        if RgCO < 0.1:
            denom = 3 - abs(P)
            if abs(denom) > 1e-6:
                invariant += -B * maxQ ** denom / denom

        # Negative invariant means the extrapolation is invalid
        if invariant <= 0:
            return None, None

        # Sv is only meaningful in the Porod regime (P ~ 4).
        # 1e4 converts pi*B/Q from A^-1 to m^2/cm^3.
        sv = 1e4 * np.pi * B / invariant if 3.95 <= P <= 4.05 else None
        return invariant * 1e24, sv   # invariant converted to cm^-4
    except Exception:
        return None, None


def load_data_from_nxcansas(file_path: str,
                            use_slit_smeared: bool = False) -> Dict[str, np.ndarray]:
    """
    Load data from NXcanSAS HDF5 file using hdf5code.readGenericNXcanSAS.

    Args:
        file_path: Full path to the HDF5 file
        use_slit_smeared: If True, look for slit-smeared data (SMR)

    Returns:
        Dictionary with Q, Intensity, Error, dQ arrays
    """
    import os
    from pyirena.io.hdf5 import readGenericNXcanSAS

    path, filename = os.path.split(file_path)
    data = readGenericNXcanSAS(path, filename)

    if data is None:
        raise ValueError(f"Could not read data from {file_path}")

    return {
        'Q': data['Q'],
        'Intensity': data['Intensity'],
        'Error': data['Error'],
        'dQ': data['dQ']
    }


# Example usage
if __name__ == "__main__":
    # Create a simple example with synthetic data
    print("Unified Fit Model - Example Usage")
    print("=" * 70)

    # Generate synthetic data (simple power law + Guinier)
    q = np.logspace(-3, 0, 100)  # 0.001 to 1.0 Å^-1

    # True parameters for level 1
    true_Rg = 50.0
    true_G = 1000.0
    true_P = 4.0
    true_B = 1e-3
    true_background = 0.01

    # Generate "measured" intensity
    model = UnifiedFitModel(num_levels=1)
    model.levels[0].Rg = true_Rg
    model.levels[0].G = true_G
    model.levels[0].P = true_P
    model.levels[0].B = true_B
    model.background = true_background

    I_true = model.calculate_intensity(q)

    # Add noise
    np.random.seed(42)
    noise = 0.05 * I_true * np.random.randn(len(q))
    I_measured = I_true + noise
    I_error = 0.05 * I_true

    # Fit with initial guess
    fit_model = UnifiedFitModel(num_levels=1)
    fit_model.levels[0].Rg = 40.0  # Initial guess
    fit_model.levels[0].G = 800.0
    fit_model.levels[0].P = 4.0
    fit_model.levels[0].B = 1e-4
    fit_model.levels[0].fit_P = False  # Fix P for this example
    fit_model.background = 0.0

    print("\nPerforming fit...")
    results = fit_model.fit(q, I_measured, I_error, verbose=1)

    print("\n" + fit_model.get_parameter_summary())

    print("\nTrue vs Fitted Parameters:")
    print(f"  Rg:  True = {true_Rg:.2f}, Fitted = {fit_model.levels[0].Rg:.2f}")
    print(f"  G:   True = {true_G:.2f}, Fitted = {fit_model.levels[0].G:.2f}")
    print(f"  B:   True = {true_B:.4e}, Fitted = {fit_model.levels[0].B:.4e}")
    print(f"  Bkg: True = {true_background:.4f}, Fitted = {fit_model.background:.4f}")

    # Calculate invariant
    inv_results = fit_model.calculate_invariant(0)
    print(f"\nInvariant Q = {inv_results['invariant']:.4e} cm⁻¹ Å⁻³")
    if inv_results['surface_to_volume'] is not None:
        print(f"Surface/Volume = {inv_results['surface_to_volume']:.4e} m²/cm³")
