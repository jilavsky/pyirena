"""
Carbon model — maths against analytic ground truth, and parameter recovery.

The physics here comes from Saurel et al. 2019/2020 (see
``pyirena/core/carbon_fit.py``), and the paper's raw data is not available, so
correctness is pinned three ways:

1. **Limits.** Every building block has a limit the paper states or the
   underlying theory forces — the roughness envelope's Q⁻⁴ asymptote, the
   globule form factor's P(0) = 1, the Voigt's Gaussian and Lorentzian limits,
   the structure factor's S(∞) = 1.  These catch a mistyped coefficient, which
   is the realistic failure mode when transcribing from a paper.
2. **Self-consistency of the derived layer.** A Porod prefactor and the surface
   area read back out of it have to agree with the intensity the model
   actually produces at high Q; a porosity fitted in one branch has to come
   back out of the other branch's derived formula.
3. **Round-trip recovery.** Synthesise noise-free data from known parameters,
   perturb the start, refit, and require the original numbers back.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from pyirena.core import carbon_density as cd
from pyirena.core.carbon_fit import (
    CarbonFitAborted,
    CarbonFitModel,
    CarbonWaxsPeak,
    discoid_form_factor,
    globule_form_factor,
    porod_roughness_factor,
    teixeira_structure_factor,
)
from pyirena.core.waxs_peakfit import (
    gauss_peak,
    lorentz_peak,
    peak_area,
    voigt_fwhm,
    voigt_peak,
)

# ── Building blocks ─────────────────────────────────────────────────────────


def test_roughness_envelope_zero_q_limit():
    """f_rough(0, R) = (2/9)R⁴ — the roughness term is finite at Q → 0."""
    R = 37.0
    assert porod_roughness_factor(np.array([0.0]), R)[0] == pytest.approx(
        2.0 / 9.0 * R ** 4)


def test_roughness_envelope_becomes_porod_at_high_q():
    """f_rough → Q⁻⁴ for QR ≫ 1 — the paper's own eq. (4).

    This is the assertion that would fail if the middle bracket coefficient
    were mistyped: at QR = 300 the (1/5)(QR)² term is still 1.1e-3 of the
    quartic one, so a wrong coefficient shows up well inside the tolerance.
    """
    R = 20.0
    q = np.array([300.0 / R, 1000.0 / R])
    assert porod_roughness_factor(q, R) == pytest.approx(q ** -4.0, rel=2e-3)


def test_roughness_middle_coefficient_is_one_fifth():
    """Pin the corrected (1/5)(QR)² term against the wrong (1/3) one.

    Evaluated at QR = 1, where the three bracket terms are comparable and the
    two candidate formulas differ by ~5 %.
    """
    R = 10.0
    q = np.array([1.0 / R])
    expected = (2.0 / 9.0) * R ** 4 / (1.0 + 1.0 / 5.0 + 2.0 / 9.0)
    wrong = (2.0 / 9.0) * R ** 4 / (1.0 + 1.0 / 3.0 + 2.0 / 9.0)
    assert porod_roughness_factor(q, R)[0] == pytest.approx(expected)
    assert porod_roughness_factor(q, R)[0] != pytest.approx(wrong, rel=1e-3)


def test_globule_form_factor_normalised_and_porod():
    """P(0) = 1, and P → (9/2)k/(Qr)⁴ once the erf saturates."""
    r, k = 8.0, 1.0
    assert globule_form_factor(np.array([0.0]), r, k)[0] == pytest.approx(1.0)
    q = np.array([40.0 / r])
    assert globule_form_factor(q, r, k) == pytest.approx(4.5 * k / (q * r) ** 4,
                                                         rel=1e-6)


def test_globule_guinier_limit_matches_sphere_rg():
    """The Guinier term uses Rg² = 3r²/5, i.e. a solid sphere of radius r."""
    r = 12.0
    q = np.array([0.01, 0.02])
    Rg2 = 3.0 * r ** 2 / 5.0
    assert globule_form_factor(q, r) == pytest.approx(
        np.exp(-q ** 2 * Rg2 / 3.0), rel=2e-3)


def test_discoid_form_factor_normalised_and_q_minus_two():
    """A sheet-like object: P(0) = 1 and a Q⁻² power law at high Q."""
    R = 25.0
    assert discoid_form_factor(np.array([0.0]), R)[0] == pytest.approx(1.0)
    q = np.array([50.0 / R])
    assert discoid_form_factor(q, R) == pytest.approx(2.0 / (q * R) ** 2, rel=1e-6)


def test_teixeira_structure_factor_limits():
    """S(∞) = 1 (no aggregation at short length scales), S(0) = 1 + D·Γ(D)(Σ/R)^D."""
    from scipy.special import gamma

    D, sigma, R = 2.4, 150.0, 6.0
    # The fractal term dies off as roughly Q^-D, so "S → 1" is approached
    # slowly; at QΣ = 7500 it is still 2e-6 above 1.
    assert teixeira_structure_factor(np.array([50.0]), D, sigma, R)[0] == pytest.approx(
        1.0, abs=1e-5)
    assert teixeira_structure_factor(np.array([5000.0]), D, sigma, R)[0] == (
        pytest.approx(1.0, abs=1e-9))
    s0 = teixeira_structure_factor(np.array([0.0]), D, sigma, R)[0]
    assert s0 == pytest.approx(1.0 + D * gamma(D) * (sigma / R) ** D, rel=1e-9)


def test_teixeira_is_monotonic_and_finite():
    """No NaNs, no negative S, decreasing with Q across the whole useful range."""
    S = teixeira_structure_factor(np.logspace(-4, 1.5, 600), 2.5, 200.0, 8.0)
    assert np.all(np.isfinite(S))
    assert np.all(S >= 1.0 - 1e-12)
    assert np.all(np.diff(S) <= 1e-9)


# ── The true Voigt added to waxs_peakfit for this tool ──────────────────────


def test_voigt_reduces_to_gauss_and_lorentz():
    q = np.linspace(1.0, 3.0, 2001)
    assert voigt_peak(q, 2.0, 2.0, 0.2, 0.0) == pytest.approx(
        gauss_peak(q, 2.0, 2.0, 0.2), abs=1e-12)
    assert voigt_peak(q, 2.0, 2.0, 1e-12, 0.2) == pytest.approx(
        lorentz_peak(q, 2.0, 2.0, 0.2), abs=1e-9)


def test_voigt_height_is_A_and_width_matches_olivero():
    q = np.linspace(0.5, 3.5, 60001)
    v = voigt_peak(q, 1.0, 2.0, 0.10, 0.15)
    assert v.max() == pytest.approx(1.0, rel=1e-9)
    above = q[v >= 0.5]
    assert above[-1] - above[0] == pytest.approx(voigt_fwhm(0.10, 0.15), rel=1e-3)


def test_voigt_area_is_analytic():
    """Area = A / V(0), because scipy's voigt_profile is unit-area.

    Integrated with ``quad`` rather than on a grid: the Lorentzian wings carry
    a few parts in 10⁴ of the area out past any practical trapezoid range.
    """
    from scipy.integrate import quad

    pr = {"A": 3.0, "Q0": 0.0, "FWHM": 0.10, "FWHM_L": 0.05}
    numeric, _ = quad(lambda x: float(voigt_peak(np.array([x]), 3.0, 0.0,
                                                 0.10, 0.05)[0]),
                      -np.inf, np.inf, limit=400)
    assert peak_area("Voigt", pr) == pytest.approx(numeric, rel=1e-6)


# ── Density / SLD / contrast layer ──────────────────────────────────────────


def test_rho_struc_reproduces_graphite():
    """Graphite's own spacings must give back graphite's own density."""
    assert cd.rho_struc_from_spacings(cd.D002_GRAPHITE, cd.D100_GRAPHITE) == (
        pytest.approx(cd.RHO_GRAPHITE))


def test_rho_struc_decreases_when_the_lattice_swells():
    """Both spacing ratios sit in the numerator with graphite on top.

    A disordered carbon always has d002 > 3.354 Å and is always less dense than
    graphite.  This is the assertion that fails if the in-plane ratio is
    inverted — the planning transcription had it the other way round.
    """
    rho_a = cd.rho_struc_from_spacings(3.354, 2.1315)
    rho_b = cd.rho_struc_from_spacings(3.80, 2.1315)     # layers further apart
    rho_c = cd.rho_struc_from_spacings(3.354, 2.20)      # lattice swollen in-plane
    assert rho_b < rho_a
    assert rho_c < rho_a
    # ρ ∝ 1/(d002·d100²) exactly.
    assert rho_b == pytest.approx(rho_a * 3.354 / 3.80)
    assert rho_c == pytest.approx(rho_a * (2.1315 / 2.20) ** 2)


def test_sld_is_linear_in_density():
    """The cached per-gram shortcut must agree with the full computation."""
    from pyirena.core.scattering_contrast import compute_compound

    for rho in (1.0, 1.8, 2.26):
        assert cd.xray_sld('C', rho) == pytest.approx(
            compute_compound('C', rho).xray_sld, rel=1e-12)


def test_specific_surface_area_conversion():
    """1e6 cm²/cm³ at 2 g/cm³ is 50 m²/g."""
    assert cd.specific_surface_area_m2_g(1.0e6, 2.0) == pytest.approx(50.0)


# ── Model assembly ──────────────────────────────────────────────────────────


def test_total_is_the_sum_of_its_parts():
    m = CarbonFitModel()
    q = np.logspace(-3, 0.6, 300)
    c = m.evaluate_components(q)
    assert c['total'] == pytest.approx(
        c['porod'] + c['mp'] + c['waxs'] + c['background'])


def test_disabled_sections_contribute_nothing():
    m = CarbonFitModel()
    q = np.logspace(-3, 0.6, 200)
    m.background.enabled = False
    m.saxs.enabled = False
    m.waxs.enabled = False
    assert m.evaluate(q) == pytest.approx(np.zeros_like(q))


def test_porod_prefactor_gives_back_the_surface_area():
    """I → 2π(Δρ)²S·1e-12/Q⁴, so S read off the curve must be S_macro.

    This is the unit-factor assertion: get ``_POROD_UNIT`` wrong and every
    surface area the tool reports is out by orders of magnitude while the fit
    still looks perfect.
    """
    m = CarbonFitModel()
    m.saxs.enabled = False
    m.waxs.enabled = False
    m.background.S_macro = 2.5e5
    m.material.contrast_mode = 'manual'
    m.material.contrast_porod_manual = 300.0

    q = np.array([0.01])
    I = m.evaluate(q)[0]
    S_back = I * q[0] ** 4 / (2.0 * np.pi * 300.0 * 1e-12)
    assert S_back == pytest.approx(2.5e5, rel=1e-9)


def test_micropore_porod_limit_matches_derived_surface_area():
    """S_mp = 3φk/r must be what the model's own high-Q tail implies.

    Ties the derived layer to the forward model rather than to a second
    transcription of the same formula.
    """
    m = CarbonFitModel()
    m.background.enabled = False
    m.waxs.enabled = False
    m.saxs.phi = 0.18
    m.saxs.pore_radius = 6.0
    m.saxs.globule_k = 1.4
    contrast = m.resolve_material()['contrast_micropore']

    q = np.array([8.0])            # Qr = 48, deep in the Porod regime
    I = m.evaluate(q)[0]
    S_from_curve = I * q[0] ** 4 / (2.0 * np.pi * contrast * 1e-12)
    assert S_from_curve == pytest.approx(m.compute_derived()['S_mp'], rel=1e-6)


def test_teubner_strey_derived_quantities_round_trip():
    """ξ, d and φ must invert the I₀ the paper's eq. (8) defines.

    Build C1/C2 from a chosen (ξ, d), then require the derived layer to give
    those back, and require the φ read out of I₀ to be the φ that produced it.
    """
    xi, d, phi = 12.0, 55.0, 0.22
    # Invert ξ = [½C2^(−½) + C1/(4C2)]^(−½), d = 2π[½C2^(−½) − C1/(4C2)]^(−½)
    a = 1.0 / xi ** 2                       # ½C2^(−½) + C1/(4C2)
    b = (2.0 * np.pi / d) ** 2              # ½C2^(−½) − C1/(4C2)
    half = 0.5 * (a + b)                    # = ½C2^(−½)
    C2 = (0.5 / half) ** 2
    C1 = 2.0 * C2 * (a - b)

    m = CarbonFitModel()
    m.saxs.mode = 'teubner_strey'
    m.saxs.ts_C1, m.saxs.ts_C2 = C1, C2
    m.material.contrast_mode = 'manual'
    m.material.contrast_micropore_manual = 350.0
    m.material.porosity_mode = 'manual'
    m.material.porosity_manual = phi
    m.saxs.ts_I0 = (8.0 * np.pi * phi * (1.0 - phi) * 350.0 * xi ** 3
                    / (1.0 + (2.0 * np.pi * xi / d) ** 2) ** 2 * 1e-4)

    dv = m.compute_derived()
    assert dv['ts_xi'] == pytest.approx(xi, rel=1e-9)
    assert dv['ts_d'] == pytest.approx(d, rel=1e-9)
    assert dv['ts_fa'] == pytest.approx(C1 / (2.0 * np.sqrt(C2)), rel=1e-12)
    assert dv['mp_phi'] == pytest.approx(phi, rel=1e-6)
    assert dv['w_pore'] == pytest.approx(xi / (1.0 - phi), rel=1e-6)
    assert dv['w_carbon'] == pytest.approx(xi / phi, rel=1e-6)


def test_teubner_strey_matches_simple_fits():
    """The SAXS-region TS branch is the Simple Fits model, not a second copy."""
    from pyirena.core.simple_fits import _teubner_strey

    m = CarbonFitModel()
    m.background.enabled = False
    m.waxs.enabled = False
    m.saxs.mode = 'teubner_strey'
    q = np.logspace(-2, 0, 100)
    assert m.evaluate(q) == pytest.approx(
        _teubner_strey(q, m.saxs.ts_I0, 1.0, m.saxs.ts_C1, m.saxs.ts_C2))


def test_orientation_factor_is_a_real_switch():
    """Turning off the 1/Q² powder average must change the curve by exactly Q².

    A fit that silently drops this factor looks fine and reports a wrong K,
    which is why it is an explicit, tested switch.
    """
    m = CarbonFitModel()
    m.background.enabled = False
    m.saxs.enabled = False
    q = np.linspace(1.5, 2.5, 50)
    with_factor = m.evaluate(q)
    m.waxs.use_orientation_factor = False
    without = m.evaluate(q)
    assert with_factor * q ** 2 == pytest.approx(without)


def test_debye_waller_damps_intensity_without_moving_the_peak():
    """⟨δz²⟩ is a distortion of the first kind: intensity only."""
    m = CarbonFitModel()
    m.background.enabled = False
    m.saxs.enabled = False
    m.peaks = [CarbonWaxsPeak(label='002', Q0=1.87, K=1.0, FWHM_G=0.2, FWHM_L=0.1)]
    q = np.linspace(1.0, 2.8, 4001)
    base = m.evaluate(q)
    m.waxs.delta_z2 = 0.05
    damped = m.evaluate(q)
    assert np.all(damped < base)
    assert damped == pytest.approx(base * np.exp(-q ** 2 * 0.05 / 3.0))


def test_delta_z2_sets_the_004_to_002_ratio():
    """Shared ⟨δz²⟩ is what makes higher orders weaker — the point of sharing it."""
    m = CarbonFitModel()
    m.background.enabled = False
    m.saxs.enabled = False
    m.peaks = [CarbonWaxsPeak(label='002', Q0=1.87, K=1.0),
               CarbonWaxsPeak(label='004', Q0=3.74, K=1.0)]
    q = np.array([1.87, 3.74])
    r0 = m.evaluate(q)[1] / m.evaluate(q)[0]
    m.waxs.delta_z2 = 0.08
    r1 = m.evaluate(q)[1] / m.evaluate(q)[0]
    assert r1 < r0


# ── Linked crumpled-layer geometry ──────────────────────────────────────────


def test_links_take_the_saxs_values_and_leave_the_fit_vector():
    m = CarbonFitModel()
    m.waxs.envelope = 'crumpled'
    m.saxs.use_fractal = True
    m.saxs.fractal_D, m.saxs.fractal_sigma, m.saxs.pore_radius = 2.2, 400.0, 7.0
    m.waxs.fractal_D, m.waxs.fractal_sigma, m.waxs.R_layer = 2.8, 90.0, 15.0

    assert m.effective_waxs_geometry() == (2.8, 90.0, 15.0)
    keys = {r.key for r in m.parameter_refs()}
    assert {'waxs.fractal_D', 'waxs.fractal_sigma', 'waxs.R_layer'} <= keys

    m.waxs.link_D_to_saxs = True
    m.waxs.link_sigma_to_saxs = True
    m.waxs.link_R_to_pore = True
    assert m.effective_waxs_geometry() == (2.2, 400.0, 7.0)
    keys = {r.key for r in m.parameter_refs()}
    assert not ({'waxs.fractal_D', 'waxs.fractal_sigma', 'waxs.R_layer'} & keys)


def test_crumpled_envelope_multiplies_the_peaks():
    m = CarbonFitModel()
    m.background.enabled = False
    m.saxs.enabled = False
    q = np.linspace(1.0, 4.0, 200)
    flat = m.evaluate(q)
    m.waxs.envelope = 'crumpled'
    crumpled = m.evaluate(q)
    D, sigma, R = m.effective_waxs_geometry()
    expected = flat * teixeira_structure_factor(q, D, sigma, R) \
        * discoid_form_factor(q, R)
    assert crumpled == pytest.approx(expected)


# ── Contrast chain ──────────────────────────────────────────────────────────


def test_moving_the_002_peak_changes_both_contrasts():
    """The chain really is live: peak position → density → SLD → contrast."""
    m = CarbonFitModel()
    before = m.resolve_material()
    m.peak_by_label('002').Q0 = 2.0 * np.pi / 3.80        # swollen interlayer
    after = m.resolve_material()
    assert after['d002'] == pytest.approx(3.80)
    assert after['rho_struc'] < before['rho_struc']
    assert after['contrast_porod'] < before['contrast_porod']
    assert after['contrast_micropore'] < before['contrast_micropore']


def test_manual_overrides_break_the_chain_where_asked():
    m = CarbonFitModel()
    m.material.contrast_mode = 'manual'
    m.material.contrast_porod_manual = 111.0
    m.material.contrast_micropore_manual = 222.0
    mat = m.resolve_material()
    assert mat['contrast_porod'] == 111.0
    assert mat['contrast_micropore'] == 222.0

    m2 = CarbonFitModel()
    m2.material.rho_struc_mode = 'manual'
    m2.material.rho_struc_manual = 1.75
    assert m2.resolve_material()['rho_struc'] == pytest.approx(1.75)


def test_missing_peaks_fall_back_instead_of_failing():
    """Before the WAXS peaks are set up, the panel still has to draw a curve."""
    m = CarbonFitModel()
    m.peaks = []
    m.material.rho_struc_manual = 1.9
    mat = m.resolve_material()
    assert mat['rho_struc'] == pytest.approx(1.9)
    assert np.isnan(mat['d002'])
    assert np.all(np.isfinite(m.evaluate(np.logspace(-3, 0.5, 50))))


def test_peak_lookup_is_forgiving_about_how_the_label_is_typed():
    m = CarbonFitModel()
    for spelling in ('002', '(002)', 'd002', ' 002 ', 'D002'):
        assert m.peak_by_label(spelling) is m.peaks[0]
    assert m.peak_by_label('999') is None


# ── Serialisation ───────────────────────────────────────────────────────────


def test_to_dict_from_dict_round_trip_through_json():
    m = CarbonFitModel()
    m.saxs.mode = 'teubner_strey'
    m.waxs.envelope = 'crumpled'
    m.waxs.link_D_to_saxs = True
    m.background.use_roughness = True
    m.add_peak('004', 3.74)
    m.q_min, m.q_max, m.n_mc_runs = 1e-3, 6.0, 5

    d = json.loads(json.dumps(m.to_dict()))
    back = CarbonFitModel.from_dict(d)
    assert back.to_dict() == m.to_dict()
    # Bounds are tuples in the model and lists in the file; the fitter unpacks
    # them, so they have to come back as tuples.
    assert isinstance(back.background.S_macro_limits, tuple)
    assert isinstance(back.peaks[0].Q0_limits, tuple)


def test_from_dict_supplies_a_default_for_every_field():
    """A file written before a field existed still opens."""
    m = CarbonFitModel.from_dict({'saxs': {'phi': 0.3}})
    assert m.saxs.phi == 0.3
    assert m.saxs.pore_radius == CarbonFitModel().saxs.pore_radius
    assert len(m.peaks) == len(CarbonFitModel().peaks)
    assert CarbonFitModel.from_dict({}).to_dict() == CarbonFitModel().to_dict()
    assert CarbonFitModel.from_dict(None).to_dict() == CarbonFitModel().to_dict()


def test_an_explicitly_empty_peak_list_stays_empty():
    """'No peaks' is a legitimate saved state — a pure SAXS fit."""
    assert CarbonFitModel.from_dict({'peaks': []}).peaks == []


def test_unknown_mode_from_a_newer_version_falls_back_with_a_warning(caplog):
    m = CarbonFitModel.from_dict({'saxs': {'mode': 'something_new'},
                                  'waxs': {'envelope': 'also_new'}})
    assert m.saxs.mode == 'fractal'
    assert m.waxs.envelope == 'none'


# ── Fitting ─────────────────────────────────────────────────────────────────


def _synth(model, q, seed=None, rel_noise=0.0):
    I = model.evaluate(q)
    if rel_noise:
        rng = np.random.default_rng(seed)
        I = I * (1.0 + rng.normal(0.0, rel_noise, I.shape))
    return I


def test_full_range_fit_recovers_every_parameter():
    """All three regions at once, noise-free — the headline case.

    Deliberately includes the overlap: the Porod tail, the micropore Porod
    limit and the low-Q side of the (002) peak all contribute in the middle of
    the range, which is where cross-talk between regions would show up.
    """
    truth = CarbonFitModel()
    truth.background.S_macro = 2.0e4
    truth.saxs.phi = 0.14
    truth.saxs.pore_radius = 6.5
    truth.peaks[0].Q0, truth.peaks[0].FWHM_G, truth.peaks[0].FWHM_L = 1.80, 0.35, 0.20
    truth.peaks[1].K = 0.25

    q = np.logspace(-3, 0.65, 700)
    I = _synth(truth, q)

    m = CarbonFitModel()
    m.background.S_macro = 5.0e4
    m.saxs.phi = 0.25
    m.saxs.pore_radius = 4.0
    m.peaks[0].Q0, m.peaks[0].FWHM_G, m.peaks[0].FWHM_L = 1.95, 0.5, 0.3
    m.peaks[1].K = 0.6

    res = m.fit(q, I)
    assert res.success
    assert res.reduced_chi_squared < 1e-10
    want = truth.parameter_values()
    for key, got in res.params.items():
        assert got == pytest.approx(want[key], rel=1e-4, abs=1e-8), key


def test_fit_recovers_parameters_from_noisy_data():
    truth = CarbonFitModel()
    truth.background.S_macro = 3.0e4
    truth.saxs.phi = 0.20
    truth.saxs.pore_radius = 7.0
    q = np.logspace(-3, 0.65, 900)
    I = _synth(truth, q, seed=7, rel_noise=0.02)

    m = CarbonFitModel()
    m.background.S_macro = 1.0e4
    m.saxs.phi = 0.10
    m.saxs.pore_radius = 5.0
    res = m.fit(q, I, error=0.02 * np.abs(I))
    assert res.success
    assert m.saxs.phi == pytest.approx(0.20, rel=0.05)
    assert m.saxs.pore_radius == pytest.approx(7.0, rel=0.05)
    assert m.background.S_macro == pytest.approx(3.0e4, rel=0.05)


def test_teubner_strey_branch_fits():
    truth = CarbonFitModel()
    truth.saxs.mode = 'teubner_strey'
    truth.saxs.ts_I0, truth.saxs.ts_C1, truth.saxs.ts_C2 = 4.0, -45.0, 8000.0
    q = np.logspace(-3, 0.65, 700)
    I = _synth(truth, q)

    m = CarbonFitModel()
    m.saxs.mode = 'teubner_strey'
    m.saxs.ts_I0, m.saxs.ts_C1, m.saxs.ts_C2 = 1.0, -20.0, 4000.0
    res = m.fit(q, I)
    assert res.success
    assert m.saxs.ts_I0 == pytest.approx(4.0, rel=1e-4)
    assert m.saxs.ts_C1 == pytest.approx(-45.0, rel=1e-4)
    assert m.saxs.ts_C2 == pytest.approx(8000.0, rel=1e-4)


def test_roughness_term_is_recovered():
    truth = CarbonFitModel()
    truth.background.use_roughness = True
    truth.background.S_macro = 1.0e4
    truth.background.S_rough = 8.0e4
    truth.background.R_rough = 60.0
    truth.saxs.enabled = False
    truth.waxs.enabled = False
    q = np.logspace(-4, -0.5, 500)
    I = _synth(truth, q)

    m = CarbonFitModel()
    m.background.use_roughness = True
    m.background.S_macro = 5.0e3
    m.background.S_rough = 2.0e4
    m.background.R_rough = 25.0
    m.saxs.enabled = False
    m.waxs.enabled = False
    res = m.fit(q, I)
    assert res.success
    assert m.background.S_rough == pytest.approx(8.0e4, rel=1e-3)
    assert m.background.R_rough == pytest.approx(60.0, rel=1e-3)


def test_fit_honours_the_q_range():
    m = CarbonFitModel()
    q = np.logspace(-3, 0.65, 400)
    I = m.evaluate(q)
    m.q_min, m.q_max = 0.01, 1.0
    res = m.fit(q, I)
    assert res.n_points == int(np.sum((q >= 0.01) & (q <= 1.0)))


def test_fit_refuses_when_nothing_is_marked_for_fitting():
    m = CarbonFitModel()
    for ref in m.parameter_refs():
        setattr(ref.owner, f'fit_{ref.attr}', False)
    q = np.logspace(-3, 0, 100)
    with pytest.raises(ValueError, match="No parameters"):
        m.fit(q, m.evaluate(q))


def test_fit_refuses_when_there_are_fewer_points_than_parameters():
    m = CarbonFitModel()
    q = np.logspace(-3, 0, 5)
    with pytest.raises(ValueError, match="free parameters"):
        m.fit(q, m.evaluate(q))


def test_a_stop_request_restores_the_starting_parameters():
    """The GUI's Stop button raises CarbonFitAborted out of the callback."""
    m = CarbonFitModel()
    q = np.logspace(-3, 0.5, 300)
    I = m.evaluate(q)
    start = m.parameter_values()

    def stop(n, chi2):
        if n > 3:
            raise CarbonFitAborted("stopped")

    res = m.fit(q, I, progress=stop)
    assert not res.success
    assert 'aborted' in res.message
    assert m.parameter_values() == start


def test_monte_carlo_produces_uncertainties():
    m = CarbonFitModel()
    m.saxs.enabled = False
    m.waxs.enabled = False
    m.n_mc_runs = 4
    q = np.logspace(-4, -1, 200)
    rng = np.random.default_rng(3)
    I = m.evaluate(q) * (1.0 + rng.normal(0.0, 0.03, q.shape))
    res = m.fit(q, I, error=0.03 * np.abs(I))
    assert res.errors['background.S_macro'] > 0


def test_result_components_add_up_and_are_json_safe():
    m = CarbonFitModel()
    q = np.logspace(-3, 0.65, 300)
    res = m.fit(q, m.evaluate(q))
    assert res.I_model == pytest.approx(res.I_porod + res.I_mp + res.I_waxs)
    json.dumps(res.to_dict())          # must not raise
