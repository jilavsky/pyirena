"""Scattering-contrast inputs that used to fall back silently (issue #38).

An explicit isotope override that cannot be honoured must fail loudly: the old
fallback to the natural element changed the molecular weight and, for
neutrons, swapped D (b_c = +6.67 fm) for H (−3.74 fm) — the opposite contrast.
"""
import pytest

pytest.importorskip("periodictable")

from pyirena.core.scattering_contrast import compute_compound, get_element_info  # noqa: E402


def test_valid_deuterium_override_unchanged():
    d2o = compute_compound("H2O", 1.107, isotope_overrides={"H": "2"})
    h2o = compute_compound("H2O", 1.0)
    assert d2o.neutron_sld > 0 > h2o.neutron_sld


@pytest.mark.parametrize("iso", ["7", "x", ""])
def test_unknown_isotope_override_raises(iso):
    with pytest.raises(ValueError, match="isotope"):
        compute_compound("H2O", 1.0, isotope_overrides={"H": iso})


def test_get_element_info_unknown_symbol_raises_valueerror():
    with pytest.raises(ValueError, match="Zz"):
        get_element_info("Zz")
    assert get_element_info("Ti")["Z"] == 22
