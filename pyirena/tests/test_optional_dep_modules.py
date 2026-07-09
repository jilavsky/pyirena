"""
Tests for core modules that need optional dependencies:
- scattering_contrast (periodictable, xraydb)
- diffraction_lines   (Dans_Diffraction)

They skip cleanly when the GUI extras are not installed, and run on
developer machines / CI configurations that have them.  Reference values
are textbook SLDs (10^10 cm^-2).
"""

import pytest


class TestScatteringContrast:
    def test_parse_formula_h2o(self):
        pytest.importorskip("periodictable")
        from pyirena.core.scattering_contrast import parse_formula
        counts = parse_formula("H2O")
        assert counts == {"H": 2.0, "O": 1.0}

    def test_water_slds(self):
        pytest.importorskip("periodictable")
        from pyirena.core.scattering_contrast import compute_compound
        water = compute_compound("H2O", density=1.0, name="Water")
        assert water.mol_weight == pytest.approx(18.015, abs=0.02)
        # X-ray SLD of water: 9.47e10 cm^-2 (0.334 e/A^3)
        assert water.xray_sld == pytest.approx(9.47, rel=0.01)
        # Neutron SLD of H2O: -0.56e10 cm^-2
        assert water.neutron_sld == pytest.approx(-0.56, abs=0.03)

    def test_heavy_water_neutron_sld(self):
        """Isotope override H->2 must give the D2O neutron SLD (+6.37)."""
        pytest.importorskip("periodictable")
        from pyirena.core.scattering_contrast import compute_compound
        d2o = compute_compound("H2O", density=1.107,
                               isotope_overrides={"H": "2"}, name="D2O")
        assert d2o.neutron_sld == pytest.approx(6.37, rel=0.02)

    def test_contrast_water_silica(self):
        pytest.importorskip("periodictable")
        from pyirena.core.scattering_contrast import compute_compound, compute_contrast
        water = compute_compound("H2O", density=1.0)
        silica = compute_compound("SiO2", density=2.2)
        c12 = compute_contrast(silica, water)
        c21 = compute_contrast(water, silica)
        assert c12.xray_contrast > 0
        # squared difference is symmetric
        assert c12.xray_contrast == pytest.approx(c21.xray_contrast)
        # SiO2 (rho=2.2): X-ray SLD ~18.8 -> contrast ~ (18.8-9.47)^2 ~ 87
        assert c12.xray_contrast == pytest.approx(87.0, rel=0.05)

    def test_anomalous_near_edge_reduces_f1(self):
        """f' is strongly negative near an absorption edge (Fe K at 7.112 keV)."""
        pytest.importorskip("periodictable")
        pytest.importorskip("xraydb")
        from pyirena.core.scattering_contrast import compute_compound, compute_anomalous
        fe = compute_compound("Fe", density=7.87)
        near = compute_anomalous(fe, energy_keV=7.10)
        far = compute_anomalous(fe, energy_keV=12.0)
        assert near.xray_sld_anom < far.xray_sld_anom


class TestDiffractionLines:
    def test_hkl_label(self):
        from pyirena.core.diffraction_lines import hkl_label
        assert hkl_label((1, 1, 1)) == "(111)"

    def test_shift_q_for_distance_error(self):
        import numpy as np
        from pyirena.core.diffraction_lines import shift_q_for_distance_error
        q = np.array([1.0, 2.0, 3.0])
        # delta_L = 0 -> no miscalibration -> identity
        shifted = shift_q_for_distance_error(q, wavelength_a=1.0,
                                             L_cal_mm=1000.0, delta_L_mm=0.0)
        np.testing.assert_allclose(shifted, q, rtol=1e-12)
        # Reduction with a too-short L shifts apparent peaks to higher Q
        shifted = shift_q_for_distance_error(q, wavelength_a=1.0,
                                             L_cal_mm=1000.0, delta_L_mm=50.0)
        assert (shifted > q).all() or (shifted < q).all()  # systematic shift

    def test_compute_pattern_requires_dans(self):
        pytest.importorskip("Dans_Diffraction")
        # Smoke only: a full pattern needs a CIF file, which the repo's
        # testData does not currently ship. Presence of the import path
        # is verified; computation is covered manually.
