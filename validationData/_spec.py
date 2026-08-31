"""Ground-truth specification of every synthetic validation dataset.

Each entry is built by a function returning ``(q, I_ideal, meta)``.  `meta`
carries the complete, human-readable ground truth that goes into README.md,
ground_truth.json and the header of every ASCII file.

Nothing here imports pyIrena: the intensities come from :mod:`_models`, which
implements every formula independently from the published literature.
"""

from __future__ import annotations

import numpy as np

import _models as M

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def logq(qmin, qmax, n):
    return np.logspace(np.log10(qmin), np.log10(qmax), int(n))


def linq(qmin, qmax, n):
    return np.linspace(qmin, qmax, int(n))


def P(name, value, unit, note=""):
    return {"name": name, "value": value, "unit": unit, "note": note}


def lognormal_moments(median, sigma, phi):
    """Closed-form moments of a log-normal *volume* distribution P_V(r).

    P_V(r) = phi * lognormal(r; median, sigma), so INTEGRAL P_V dr = phi.
      mode                = median * exp(-sigma^2)
      volume-weighted <r> = median * exp(sigma^2 / 2)
      RMS radius          = median * exp(sigma^2)
        (this is the quantity pyIrena's Size Distribution reports as "Rg":
         sqrt(INTEGRAL r^2 P dr / INTEGRAL P dr) -- a radius moment, not a
         true radius of gyration)
    """
    return {
        "mode": median * np.exp(-sigma ** 2),
        "mean_volume_weighted": median * np.exp(sigma ** 2 / 2.0),
        "rms_radius": median * np.exp(sigma ** 2),
        "volume_fraction": phi,
    }


def _fine_r_grid(median, sigma, shift=0.0, n=4000, n_sigma=6.0):
    """Dense integration grid spanning +/- n_sigma of a log-normal in log space."""
    lo = shift + median * np.exp(-n_sigma * sigma)
    hi = shift + median * np.exp(+n_sigma * sigma)
    return np.logspace(np.log10(lo), np.log10(hi), n)


# ===========================================================================
# Size Distribution
# ===========================================================================

def sizes_sphere_lognormal():
    median, sigma, phi, contrast = 100.0, 0.25, 0.010, 100.0
    q = logq(8e-4, 0.35, 350)
    r = _fine_r_grid(median, sigma)
    fv = M.lognormal_pdf(r, median, sigma)
    fv *= phi / np.trapezoid(fv, r)
    I = M.size_distribution_intensity(q, r, fv, contrast)
    mom = lognormal_moments(median, sigma, phi)
    meta = dict(
        tool="Size Distribution",
        model="Single log-normal volume distribution of solid spheres",
        equation="I(Q) = 1e-4 * contrast * INT V(r) F_sph^2(Qr) P_V(r) dr",
        parameters=[
            P("distribution", "log-normal (volume weighted)", "-"),
            P("R_median", median, "A", "median of P_V(r)"),
            P("sigma_log", sigma, "-", "log-space standard deviation"),
            P("volume_fraction", phi, "-", "INT P_V(r) dr"),
            P("contrast", contrast, "1e20 cm^-4", "(delta-rho)^2"),
            P("background", 0.0, "1/cm"),
        ],
        derived=[
            P("R_mode", mom["mode"], "A", "median*exp(-sigma^2)"),
            P("R_mean_vol", mom["mean_volume_weighted"], "A", "median*exp(sigma^2/2)"),
            P("R_rms", mom["rms_radius"], "A", "reported by pyIrena as 'Rg'"),
        ],
    )
    return q, I, meta


def sizes_sphere_bimodal():
    contrast = 100.0
    m1, s1, p1 = 40.0, 0.18, 0.0040
    m2, s2, p2 = 250.0, 0.22, 0.0100
    q = logq(5e-4, 0.4, 400)
    r = np.logspace(np.log10(5.0), np.log10(1500.0), 6000)
    f1 = M.lognormal_pdf(r, m1, s1); f1 *= p1 / np.trapezoid(f1, r)
    f2 = M.lognormal_pdf(r, m2, s2); f2 *= p2 / np.trapezoid(f2, r)
    I = M.size_distribution_intensity(q, r, f1 + f2, contrast)
    meta = dict(
        tool="Size Distribution",
        model="Bimodal log-normal volume distribution of solid spheres",
        equation="I(Q) = 1e-4 * contrast * INT V(r) F_sph^2(Qr) [P_1(r)+P_2(r)] dr",
        parameters=[
            P("R_median_1", m1, "A"), P("sigma_log_1", s1, "-"),
            P("volume_fraction_1", p1, "-"),
            P("R_median_2", m2, "A"), P("sigma_log_2", s2, "-"),
            P("volume_fraction_2", p2, "-"),
            P("contrast", contrast, "1e20 cm^-4"),
            P("background", 0.0, "1/cm"),
        ],
        derived=[
            P("R_mode_1", m1 * np.exp(-s1 ** 2), "A"),
            P("R_mode_2", m2 * np.exp(-s2 ** 2), "A"),
            P("volume_fraction_total", p1 + p2, "-"),
        ],
    )
    return q, I, meta


def sizes_sphere_flat_background():
    median, sigma, phi, contrast = 120.0, 0.22, 0.0060, 100.0
    bg_flat = 0.020
    q = logq(8e-4, 0.35, 350)
    r = _fine_r_grid(median, sigma)
    fv = M.lognormal_pdf(r, median, sigma)
    fv *= phi / np.trapezoid(fv, r)
    I = M.size_distribution_intensity(q, r, fv, contrast) + bg_flat
    mom = lognormal_moments(median, sigma, phi)
    meta = dict(
        tool="Size Distribution",
        model="Log-normal spheres on a flat background",
        equation=("I(Q) = 1e-4*contrast*INT V(r) F_sph^2(Qr) P_V(r) dr + flat"),
        parameters=[
            P("R_median", median, "A"), P("sigma_log", sigma, "-"),
            P("volume_fraction", phi, "-"),
            P("contrast", contrast, "1e20 cm^-4"),
            P("background_flat", bg_flat, "1/cm",
              "dominates above Q ~ 0.15 1/A; fit it there before inverting"),
        ],
        derived=[
            P("R_mode", mom["mode"], "A"),
            P("R_mean_vol", mom["mean_volume_weighted"], "A"),
            P("R_rms", mom["rms_radius"], "A"),
        ],
    )
    return q, I, meta


def sizes_spheroid_lognormal():
    median, sigma, phi, contrast, ar = 150.0, 0.20, 0.0080, 50.0, 0.30
    q = logq(5e-4, 0.30, 320)
    r = _fine_r_grid(median, sigma, n=1200)
    fv = M.lognormal_pdf(r, median, sigma)
    fv *= phi / np.trapezoid(fv, r)
    I = M.size_distribution_intensity(q, r, fv, contrast,
                                      shape="spheroid", aspect_ratio=ar)
    mom = lognormal_moments(median, sigma, phi)
    meta = dict(
        tool="Size Distribution",
        model="Log-normal volume distribution of oblate spheroids (AR = 0.3)",
        equation=("I(Q) = 1e-4*contrast*INT V(r,AR) <F_ell^2(Q,r,AR)>_orient"
                  " P_V(r) dr"),
        parameters=[
            P("shape", "spheroid", "-"),
            P("aspect_ratio", ar, "-", "semi-axes (r, r, AR*r)"),
            P("R_median", median, "A"), P("sigma_log", sigma, "-"),
            P("volume_fraction", phi, "-"),
            P("contrast", contrast, "1e20 cm^-4"),
            P("background", 0.0, "1/cm"),
        ],
        derived=[P("R_mode", mom["mode"], "A"),
                 P("R_rms", mom["rms_radius"], "A")],
    )
    return q, I, meta


# ===========================================================================
# Unified Fit
# ===========================================================================

def unified_one_level():
    G, Rg, P_, B, bg = 100.0, 200.0, 4.00, 3.045e-7, 0.010
    q = logq(5e-4, 0.30, 350)
    I = M.unified_level(q, G, Rg, B, P_) + bg
    meta = dict(
        tool="Unified Fit",
        model="Single Beaucage Unified level (Guinier + Porod), flat background",
        equation=("I(Q) = G exp(-Q^2 Rg^2/3) + B/Q*^P + bg,"
                  "  Q* = Q/erf(K Q Rg/sqrt6)^3,  K = 1 for P>3"),
        parameters=[
            P("G", G, "1/cm", "Guinier prefactor"),
            P("Rg", Rg, "A"),
            P("P", P_, "-", "power-law slope"),
            P("B", B, "cm^-1 A^-P", "= G exp(-P/2)(3P/2)^(P/2)/Rg^P (smooth join)"),
            P("background", bg, "1/cm"),
            P("K", 1.0, "-"),
        ],
        derived=[
            P("I(0)", G + bg, "1/cm"),
            P("Q_rollover", 2.0 / Rg * np.sqrt(P_ / 2.0), "1/A"),
        ],
    )
    return q, I, meta


def unified_two_level():
    # level 1 = small (high Q), level 2 = large (low Q), Irena numbering
    G1, Rg1, P1, B1 = 8.0, 120.0, 4.00, 1.880e-7
    G2, Rg2, P2, B2, RgCO2 = 4000.0, 1200.0, 3.20, 1.393e-6, 120.0
    bg = 0.020
    q = logq(1e-4, 0.30, 450)
    I = (M.unified_level(q, G1, Rg1, B1, P1)
         + M.unified_level(q, G2, Rg2, B2, P2, RgCO=RgCO2) + bg)
    meta = dict(
        tool="Unified Fit",
        model="Two-level Beaucage Unified model with level-2 low-Q cut-off",
        equation="I(Q) = SUM_i [G_i exp(-Q^2 Rg_i^2/3) + B_i/Q*_i^P_i exp(-Q^2 RgCO_i^2/3)] + bg",
        parameters=[
            P("G_1", G1, "1/cm"), P("Rg_1", Rg1, "A"),
            P("P_1", P1, "-"), P("B_1", B1, "cm^-1 A^-P"),
            P("RgCO_1", 0.0, "A"),
            P("G_2", G2, "1/cm"), P("Rg_2", Rg2, "A"),
            P("P_2", P2, "-"), P("B_2", B2, "cm^-1 A^-P"),
            P("RgCO_2", RgCO2, "A", "linked to Rg of level 1"),
            P("background", bg, "1/cm"),
        ],
        derived=[P("I(0)", G1 + G2 + bg, "1/cm")],
    )
    return q, I, meta


def unified_one_level_slitsmeared():
    G, Rg, P_, B, bg = 100.0, 200.0, 4.00, 3.045e-7, 0.010
    SL = 0.030
    q = logq(5e-4, 0.30, 350)

    def ideal(qq):
        return M.unified_level(qq, G, Rg, B, P_) + bg

    I = M.slit_smear(ideal, q, SL)
    meta = dict(
        tool="Unified Fit (slit smearing)",
        model="Single Unified level, infinite-slit smeared",
        equation="I_smr(Q) = (1/SL) INT_0^SL I_ideal(sqrt(Q^2+l^2)) dl",
        parameters=[
            P("G", G, "1/cm"), P("Rg", Rg, "A"), P("P", P_, "-"),
            P("B", B, "cm^-1 A^-P"), P("background", bg, "1/cm"),
            P("slit_length_dQl", SL, "1/A", "written as entry/sasdata/dQl in the HDF5"),
        ],
        derived=[],
        slit_length=SL,
    )
    return q, I, meta


# ===========================================================================
# Modeling
# ===========================================================================

def modeling_sizedist_plus_unified():
    median, sigma, scale, contrast = 30.0, 0.30, 0.0050, 100.0
    G, Rg, P_ = 20000.0, 4000.0, 3.50
    B = M.unified_B_from_G(G, Rg, P_)
    bg = 0.010
    q = logq(5e-5, 0.50, 480)
    r = _fine_r_grid(median, sigma)
    fv = M.lognormal_pdf(r, median, sigma)
    fv *= scale / np.trapezoid(fv, r)
    I = (M.size_distribution_intensity(q, r, fv, contrast)
         + M.unified_level(q, G, Rg, B, P_) + bg)
    meta = dict(
        tool="Modeling",
        model="Population 1: log-normal spheres.  Population 2: Unified level.",
        equation="I(Q) = I_sizedist(Q) + I_unified(Q) + bg",
        parameters=[
            P("pop1_type", "size_dist / lognormal / sphere", "-"),
            P("pop1_min_size", 0.0, "A"),
            P("pop1_mean_size", median, "A", "log-normal median"),
            P("pop1_sdeviation", sigma, "-"),
            P("pop1_scale", scale, "-", "INT P_V(r) dr"),
            P("pop1_contrast", contrast, "1e20 cm^-4"),
            P("pop2_type", "unified_level", "-"),
            P("pop2_G", G, "1/cm"), P("pop2_Rg", Rg, "A"),
            P("pop2_P", P_, "-"), P("pop2_B", B, "cm^-1 A^-P"),
            P("background", bg, "1/cm"),
        ],
        derived=[],
    )
    return q, I, meta


def modeling_sizedist_plus_peak():
    median, sigma, scale, contrast = 60.0, 0.25, 0.0040, 80.0
    amp, pos, width = 5.0, 0.150, 0.0120
    bg = 0.005
    q = logq(3e-3, 0.50, 450)
    r = _fine_r_grid(median, sigma)
    fv = M.lognormal_pdf(r, median, sigma)
    fv *= scale / np.trapezoid(fv, r)
    peak = amp * np.exp(-(q - pos) ** 2 / (2.0 * width ** 2))
    I = M.size_distribution_intensity(q, r, fv, contrast) + peak + bg
    meta = dict(
        tool="Modeling",
        model="Population 1: log-normal spheres.  Population 2: Gaussian diffraction peak.",
        equation="I(Q) = I_sizedist(Q) + A exp(-(Q-Q0)^2/(2 w^2)) + bg",
        parameters=[
            P("pop1_type", "size_dist / lognormal / sphere", "-"),
            P("pop1_min_size", 0.0, "A"),
            P("pop1_mean_size", median, "A"),
            P("pop1_sdeviation", sigma, "-"),
            P("pop1_scale", scale, "-"),
            P("pop1_contrast", contrast, "1e20 cm^-4"),
            P("pop2_type", "diffraction_peak / gaussian", "-"),
            P("pop2_amplitude", amp, "1/cm"),
            P("pop2_position", pos, "1/A"),
            P("pop2_width", width, "1/A", "Gaussian sigma, NOT FWHM"),
            P("background", bg, "1/cm"),
        ],
        derived=[P("peak_FWHM", 2.0 * np.sqrt(2.0 * np.log(2.0)) * width, "1/A")],
    )
    return q, I, meta


def modeling_hardsphere():
    median, sigma, scale, contrast = 50.0, 0.12, 0.200, 100.0
    hs_R, hs_phi = 50.0, 0.250
    bg = 0.010
    q = logq(3e-3, 0.50, 450)
    r = _fine_r_grid(median, sigma)
    fv = M.lognormal_pdf(r, median, sigma)
    fv *= scale / np.trapezoid(fv, r)
    I = (M.size_distribution_intensity(q, r, fv, contrast)
         * M.hard_sphere_sf(q, hs_R, hs_phi)) + bg
    meta = dict(
        tool="Modeling",
        model="Log-normal spheres with a Percus-Yevick hard-sphere structure factor",
        equation="I(Q) = I_sizedist(Q) * S_PY(Q; R_HS, phi_HS) + bg",
        parameters=[
            P("pop1_type", "size_dist / lognormal / sphere", "-"),
            P("pop1_min_size", 0.0, "A"),
            P("pop1_mean_size", median, "A"),
            P("pop1_sdeviation", sigma, "-"),
            P("pop1_scale", scale, "-"),
            P("pop1_contrast", contrast, "1e20 cm^-4"),
            P("structure_factor", "hard_sphere", "-"),
            P("sf_radius", hs_R, "A"),
            P("sf_volume_fraction", hs_phi, "-"),
            P("background", bg, "1/cm"),
        ],
        derived=[],
    )
    return q, I, meta


def modeling_mass_fractal():
    phi, R, Dv, ksi, eta, contrast = 0.0050, 40.0, 2.40, 800.0, 0.50, 100.0
    bg = 0.010
    q = logq(3e-4, 0.30, 400)
    I = M.mass_fractal(q, phi, R, Dv, ksi, eta, contrast) + bg
    meta = dict(
        tool="Modeling",
        model="Mass-fractal aggregate of spherical primary particles (Teixeira 1988)",
        equation=("I(Q) = phi*contrast*1e-4*V*[bracket*S_f(Q)+(1-eta)^2]*F_sph^2(QR) + bg"),
        parameters=[
            P("pop_type", "mass_fractal", "-"),
            P("Phi", phi, "-"), P("Radius", R, "A"),
            P("Dv", Dv, "-", "mass fractal dimension"),
            P("Ksi", ksi, "A", "aggregate correlation length"),
            P("Eta", eta, "-", "packing term"),
            P("Beta", 1.0, "-", "aspect ratio (sphere)"),
            P("Contrast", contrast, "1e20 cm^-4"),
            P("background", bg, "1/cm"),
        ],
        derived=[],
    )
    return q, I, meta


def modeling_surface_fractal():
    surface, Ds, ksi, contrast = 2.0e4, 2.40, 600.0, 100.0
    bg = 0.010
    q = logq(5e-4, 0.30, 380)
    I = M.surface_fractal(q, surface, Ds, ksi, contrast) + bg
    meta = dict(
        tool="Modeling",
        model="Surface-fractal scattering (Teixeira 1988)",
        equation=("I(Q) = pi*contrast*Ksi^4*S*Gamma(5-Ds)"
                  "*sin[(3-Ds)atan(Q Ksi)]/[(1+Q^2Ksi^2)^((5-Ds)/2) Q Ksi] + bg"),
        parameters=[
            P("pop_type", "surface_fractal", "-"),
            P("Surface", surface, "1/cm"),
            P("Ds", Ds, "-", "surface fractal dimension"),
            P("Ksi", ksi, "A"),
            P("Contrast", contrast, "1e20 cm^-4"),
            P("background", bg, "1/cm"),
        ],
        derived=[P("Porod_slope_high_Q", 6.0 - Ds, "-")],
    )
    return q, I, meta


def modeling_guinier_porod():
    G, Rg1, s1, P_ = 150.0, 90.0, 1.00, 3.60
    bg = 0.005
    q = logq(2e-3, 0.40, 400)
    I = M.guinier_porod(q, G, Rg1, s1, P_) + bg
    meta = dict(
        tool="Modeling",
        model="Guinier-Porod model, single level with dimensionality s1 = 1 (rod-like)",
        equation=("I(Q<Q1) = G Q^-s1 exp(-Q^2 Rg1^2/(3-s1)); I(Q>=Q1) = D Q^-P"),
        parameters=[
            P("pop_type", "guinier_porod", "-"),
            P("G", G, "1/cm"), P("Rg1", Rg1, "A"),
            P("s1", s1, "-", "0 = globular, 1 = rod, 2 = lamella"),
            P("P", P_, "-"),
            P("Rg2", 1e10, "A", "collapsed (single level)"),
            P("s2", 0.0, "-"),
            P("background", bg, "1/cm"),
        ],
        derived=[P("Q1", np.sqrt((P_ - s1) * (3.0 - s1) / 2.0) / Rg1, "1/A")],
    )
    return q, I, meta


# ===========================================================================
# Simple Fits
# ===========================================================================

def simple_guinier():
    I0, Rg = 250.0, 45.0
    q = logq(2e-3, 0.030, 200)
    meta = dict(tool="Simple Fits", model="Guinier",
                equation="I(Q) = I0 exp(-Q^2 Rg^2/3)",
                parameters=[P("I0", I0, "1/cm"), P("Rg", Rg, "A")],
                derived=[P("Q_max*Rg", 0.030 * Rg, "-", "Guinier validity limit ~1.3")])
    return q, M.guinier(q, I0, Rg), meta


def simple_guinier_rod():
    I0, Rc = 5.0, 25.0
    q = logq(5e-3, 0.080, 200)
    meta = dict(tool="Simple Fits", model="Guinier Rod",
                equation="I(Q) = I0 exp(-Q^2 Rc^2/2)/Q",
                parameters=[P("I0", I0, "cm^-1 A^-1"), P("Rc", Rc, "A")],
                derived=[])
    return q, M.guinier_rod(q, I0, Rc), meta


def simple_guinier_sheet():
    I0, Rg = 0.020, 15.0
    q = logq(1e-2, 0.100, 200)
    meta = dict(tool="Simple Fits", model="Guinier Sheet",
                equation="I(Q) = I0 exp(-Q^2 Rg^2)/Q^2",
                parameters=[P("I0", I0, "cm^-1 A^-2"), P("Rg", Rg, "A")],
                derived=[])
    return q, M.guinier_sheet(q, I0, Rg), meta


def simple_porod():
    Kp, bg = 2.5e-6, 0.050
    q = logq(0.050, 0.500, 200)
    meta = dict(tool="Simple Fits", model="Porod",
                equation="I(Q) = Kp Q^-4 + bg",
                parameters=[P("Kp", Kp, "cm^-1 A^-4"), P("Background", bg, "1/cm")],
                derived=[])
    return q, M.porod(q, Kp, bg), meta


def simple_power_law():
    A, n, bg = 1.2e-3, 3.20, 0.020
    q = logq(0.010, 0.400, 250)
    meta = dict(tool="Simple Fits", model="Power Law",
                equation="I(Q) = Prefactor Q^-Exponent + bg",
                parameters=[P("Prefactor", A, "cm^-1 A^-n"),
                            P("Exponent", n, "-"),
                            P("Background", bg, "1/cm")],
                derived=[])
    return q, M.power_law(q, A, n, bg), meta


def simple_sphere():
    scale, R, bg = 8.0, 120.0, 0.010
    q = logq(2e-3, 0.150, 500)
    meta = dict(tool="Simple Fits", model="Sphere (monodisperse) + flat background",
                equation="I(Q) = Scale |3(sin x - x cos x)/x^3|^2 + BG_flat, x = QR",
                parameters=[P("Scale", scale, "1/cm"), P("R", R, "A"),
                            P("BG_flat", bg, "1/cm",
                              "fit with the complex background, flat term only")],
                derived=[P("Q_first_minimum", 4.4934 / R, "1/A")])
    return q, M.sphere_model(q, scale, R) + bg, meta


def simple_debye_chain():
    scale, Rg = 15.0, 60.0
    q = logq(3e-3, 0.300, 300)
    meta = dict(tool="Simple Fits", model="Debye Polymer Chain",
                equation="I(Q) = Scale 2(exp(-x)-1+x)/x^2, x = Q^2 Rg^2",
                parameters=[P("Scale", scale, "1/cm"), P("Rg", Rg, "A")],
                derived=[])
    return q, M.debye_polymer_chain(q, scale, Rg), meta


def simple_debye_bueche():
    pref, eta, xi = 1.0, 0.050, 80.0
    q = logq(2e-3, 0.200, 300)
    meta = dict(tool="Simple Fits", model="Debye-Bueche",
                equation="I(Q) = Prefactor Eta^2 xi^3/(1+Q^2 xi^2)^2",
                parameters=[P("Prefactor", pref, "cm^-1 A^-3"),
                            P("Eta", eta, "-"),
                            P("CorrLength", xi, "A")],
                derived=[P("I(0)", pref * eta ** 2 * xi ** 3, "1/cm",
                           "only the product Prefactor*Eta^2 is determined by a fit")])
    return q, M.debye_bueche(q, pref, eta, xi), meta


# ===========================================================================
# WAXS Peak Fit
# ===========================================================================

def waxs_three_gauss():
    peaks = [(100.0, 1.850, 0.090), (60.0, 2.550, 0.120), (35.0, 3.050, 0.150)]
    bg0, bg1 = 8.0, 1.50
    q = linq(1.20, 4.00, 700)
    I = bg0 + bg1 * q
    for A, Q0, W in peaks:
        I = I + M.gauss_peak(q, A, Q0, W)
    params = [P("background_shape", "Linear", "-"),
              P("bg0", bg0, "1/cm"), P("bg1", bg1, "cm^-1 A")]
    derived = []
    for i, (A, Q0, W) in enumerate(peaks, 1):
        params += [P(f"peak{i}_shape", "Gauss", "-"),
                   P(f"peak{i}_A", A, "1/cm", "peak height at Q0"),
                   P(f"peak{i}_Q0", Q0, "1/A"),
                   P(f"peak{i}_FWHM", W, "1/A")]
        derived += [P(f"peak{i}_area", A * W * np.sqrt(np.pi / (4.0 * np.log(2.0))),
                      "cm^-1 A^-1", "A*FWHM*sqrt(pi/(4 ln2))"),
                    P(f"peak{i}_d_spacing", 2.0 * np.pi / Q0, "A")]
    meta = dict(tool="WAXS Peak Fit",
                model="Three Gaussian peaks on a linear background",
                equation="I(Q) = bg0 + bg1 Q + SUM_i A_i exp[-4 ln2 (Q-Q0_i)^2/FWHM_i^2]",
                parameters=params, derived=derived)
    return q, I, meta


def waxs_pseudovoigt():
    peaks = [(200.0, 2.100, 0.070, 0.40), (90.0, 2.950, 0.110, 0.40)]
    coeffs = [10.0, -2.00, 0.80, -0.050]
    q = linq(1.50, 4.00, 650)
    I = sum(c * q ** i for i, c in enumerate(coeffs))
    for A, Q0, W, eta in peaks:
        I = I + M.pseudo_voigt_peak(q, A, Q0, W, eta)
    params = [P("background_shape", "Cubic", "-")]
    params += [P(f"bg{i}", c, f"cm^-1 A^{i}") for i, c in enumerate(coeffs)]
    derived = []
    for i, (A, Q0, W, eta) in enumerate(peaks, 1):
        params += [P(f"peak{i}_shape", "Pseudo-Voigt", "-"),
                   P(f"peak{i}_A", A, "1/cm"), P(f"peak{i}_Q0", Q0, "1/A"),
                   P(f"peak{i}_FWHM", W, "1/A"), P(f"peak{i}_eta", eta, "-")]
        K = eta * np.pi / 2.0 + (1.0 - eta) * np.sqrt(np.pi / (4.0 * np.log(2.0)))
        derived += [P(f"peak{i}_area", A * W * K, "cm^-1 A^-1"),
                    P(f"peak{i}_d_spacing", 2.0 * np.pi / Q0, "A")]
    meta = dict(tool="WAXS Peak Fit",
                model="Two pseudo-Voigt peaks on a cubic-polynomial background",
                equation="I(Q) = SUM_k bg_k Q^k + SUM_i A_i[eta L_i(Q) + (1-eta) G_i(Q)]",
                parameters=params, derived=derived)
    return q, I, meta


# ===========================================================================
# Data Merge  (two overlapping curves from one true model)
# ===========================================================================

_MERGE_TRUTH = dict(G1=8.0, Rg1=120.0, P1=4.00, B1=1.880e-7,
                    G2=4000.0, Rg2=1200.0, P2=3.20, B2=1.393e-6,
                    RgCO2=120.0, bg=0.020)
_MERGE_SCALE = 1.150


def _merge_model(q):
    t = _MERGE_TRUTH
    return (M.unified_level(q, t["G1"], t["Rg1"], t["B1"], t["P1"])
            + M.unified_level(q, t["G2"], t["Rg2"], t["B2"], t["P2"], RgCO=t["RgCO2"])
            + t["bg"])


def merge_usaxs():
    q = logq(1.0e-4, 0.050, 300)
    meta = dict(tool="Data Merge",
                model="Low-Q (USAXS-like) branch of a two-level Unified curve, correct scale",
                equation="I(Q) = I_unified_2level(Q)",
                parameters=[P("scale_factor", 1.0, "-", "reference branch")]
                + [P(k, v, "-") for k, v in _MERGE_TRUTH.items()],
                derived=[])
    return q, _merge_model(q), meta


def merge_saxs():
    q = logq(0.0050, 0.500, 320)
    meta = dict(tool="Data Merge",
                model=("High-Q (SAXS-like) branch of the SAME curve, deliberately "
                       "mis-scaled by 1.15"),
                equation="I(Q) = 1.15 * I_unified_2level(Q)",
                parameters=[P("scale_factor", _MERGE_SCALE, "-",
                              "merge should recover 1/1.15 = 0.869565 for this branch")]
                + [P(k, v, "-") for k, v in _MERGE_TRUTH.items()],
                derived=[P("expected_recovered_scale", 1.0 / _MERGE_SCALE, "-")])
    return q, _MERGE_SCALE * _merge_model(q), meta


# ===========================================================================
# Data Manipulation
# ===========================================================================

def manipulate_reference():
    q = logq(1e-3, 0.30, 300)
    meta = dict(tool="Data Manipulation",
                model="Reference curve (single Unified level)",
                equation="I(Q) = G exp(-Q^2Rg^2/3) + B/Q*^P",
                parameters=[P("G", 100.0, "1/cm"), P("Rg", 200.0, "A"),
                            P("P", 4.0, "-"), P("B", 3.045e-7, "cm^-1 A^-P")],
                derived=[])
    return q, M.unified_level(q, 100.0, 200.0, 3.045e-7, 4.00), meta


def manipulate_input():
    q, I, _ = manipulate_reference()
    a, b = 2.50, 0.300
    meta = dict(
        derive_from=("manipulate_reference", a, b),tool="Data Manipulation",
                model="The reference curve scaled and offset by known constants",
                equation="I_input(Q) = a * I_reference(Q) + b",
                parameters=[P("a_multiply", a, "-"), P("b_add", b, "1/cm")],
                derived=[P("inverse_operation",
                           "subtract 0.3 then divide by 2.5", "-",
                           "must reproduce manipulate_reference to machine precision")])
    meta["equation"] += "  (applied to the NOISY reference, so the inverse is exact)"
    return q, a * I + b, meta
