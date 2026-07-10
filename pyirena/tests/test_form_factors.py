"""
Tests for pyirena.core.form_factors — pure functions, tested against
physical limits and known analytic values.
"""

import numpy as np
import pytest

from pyirena.core.form_factors import (
    sphere_ff,
    spheroid_ff,
    build_g_matrix,
    make_r_grid,
    bin_widths,
)


def _sphere_volume(r):
    return (4.0 / 3.0) * np.pi * r ** 3


class TestSphereFF:
    def test_q_zero_limit_is_volume(self):
        """At Q -> 0 the sphere FF (per unit vol-frac) equals V(r) [A^3]."""
        r = 50.0
        ff = sphere_ff(np.array([1e-8]), r)
        assert ff[0] == pytest.approx(_sphere_volume(r), rel=1e-6)

    def test_first_minimum_position(self):
        """First zero of the sphere amplitude is at q*r ~ 4.4934."""
        r = 50.0
        q = np.linspace(3.5, 5.5, 20001) / r
        ff = sphere_ff(q, r)
        q_min = q[np.argmin(ff)]
        assert q_min * r == pytest.approx(4.4934, abs=0.001)

    def test_monotonic_decrease_in_guinier_region(self):
        r = 50.0
        q = np.linspace(1e-4, 1.0 / r, 100)
        ff = sphere_ff(q, r)
        assert (np.diff(ff) < 0).all()

    def test_nonnegative_everywhere(self):
        q = np.logspace(-4, 0, 500)
        assert (sphere_ff(q, 123.0) >= 0).all()


class TestSpheroidFF:
    def test_aspect_ratio_one_equals_sphere(self):
        r = 40.0
        q = np.logspace(-3, -0.8, 100)
        np.testing.assert_allclose(
            spheroid_ff(q, r, aspect_ratio=1.0),
            sphere_ff(q, r),
            rtol=1e-4,
        )

    def test_q_zero_limit_is_spheroid_volume(self):
        """V_spheroid = (4/3) pi r^3 AR for semi-axes (r, r, r*AR)."""
        r, ar = 40.0, 2.0
        ff = spheroid_ff(np.array([1e-8]), r, aspect_ratio=ar)
        assert ff[0] == pytest.approx(_sphere_volume(r) * ar, rel=1e-4)

    def test_prolate_differs_from_sphere_at_high_q(self):
        r = 40.0
        q = np.array([0.05])
        assert spheroid_ff(q, r, aspect_ratio=3.0)[0] != pytest.approx(
            sphere_ff(q, r)[0], rel=0.01)


class TestBuildGMatrix:
    def test_shape_and_positivity(self):
        q = np.logspace(-3, -1, 30)
        r_grid = make_r_grid(10, 500, 20)
        G = build_g_matrix(q, r_grid, shape="sphere")
        assert G.shape == (30, 20)
        assert (G >= 0).all()
        assert np.isfinite(G).all()

    def test_sphere_column_matches_sphere_ff(self):
        q = np.logspace(-3, -1, 30)
        r_grid = np.array([50.0, 100.0])
        G = build_g_matrix(q, r_grid, shape="sphere", contrast=1.0)
        np.testing.assert_allclose(G[:, 0], sphere_ff(q, 50.0) * 1e-4, rtol=1e-10)

    def test_contrast_scaling(self):
        q = np.logspace(-3, -1, 10)
        r_grid = np.array([50.0])
        G1 = build_g_matrix(q, r_grid, shape="sphere", contrast=1.0)
        G2 = build_g_matrix(q, r_grid, shape="sphere", contrast=3.0)
        np.testing.assert_allclose(G2, 3.0 * G1, rtol=1e-12)

    def test_unknown_shape_raises(self):
        with pytest.raises(ValueError, match="Unknown shape"):
            build_g_matrix(np.array([0.01]), np.array([50.0]), shape="banana")

    def test_core_shell_zero_shell_reduces_to_sphere(self):
        """cs_sphere_by_core with t_shell=0 and shell SLD == solvent SLD
        must reduce to a plain sphere with contrast (sld_core-sld_solvent)^2."""
        q = np.logspace(-3, -1, 50)
        r_grid = np.array([60.0])
        sld_core, sld_solvent = 10.0, 2.0
        G_cs = build_g_matrix(
            q, r_grid, shape="cs_sphere_by_core",
            sld_core=sld_core, sld_shell=sld_solvent, sld_solvent=sld_solvent,
            t_shell=0.0,
        )
        contrast = (sld_core - sld_solvent) ** 2  # in 10^20 cm^-4 units
        G_sphere = build_g_matrix(q, r_grid, shape="sphere", contrast=contrast)
        np.testing.assert_allclose(G_cs, G_sphere, rtol=1e-6)


class TestGridHelpers:
    def test_make_r_grid_linear(self):
        g = make_r_grid(10, 100, 10)
        assert len(g) == 10
        assert g[0] == 10 and g[-1] == 100
        np.testing.assert_allclose(np.diff(g), np.diff(g)[0])

    def test_make_r_grid_log(self):
        g = make_r_grid(10, 1000, 21, log_spacing=True)
        assert g[0] == pytest.approx(10) and g[-1] == pytest.approx(1000)
        ratios = g[1:] / g[:-1]
        np.testing.assert_allclose(ratios, ratios[0])

    def test_bin_widths_sum_covers_range(self):
        g = np.linspace(10.0, 100.0, 10)
        w = bin_widths(g)
        assert len(w) == len(g)
        assert (w > 0).all()
        # Trapezoidal widths: interior = (r[i+1]-r[i-1])/2, ends = one step;
        # for a uniform grid every width equals the step.
        np.testing.assert_allclose(w, 10.0)
