"""
Tests for pyirena.core.form_factors — pure functions, tested against
physical limits and known analytic values.
"""

import numpy as np
import pytest

from pyirena.core.form_factors import (
    bin_widths,
    build_g_matrix,
    make_r_grid,
    sphere_ff,
    spheroid_ff,
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


class TestSpheroidQuadratureOrder:
    """The orientation average drops nodes near AR = 1; it must not drift.

    ``_orientation_nodes`` trades quadrature order against |AR - 1| so that a
    near-spherical spheroid does not pay for 50 orientations it cannot
    resolve. These tests pin the two things that buys us: the answer still
    matches the 50-node Igor reference, and AR = 1 costs nothing at all.
    """

    @staticmethod
    def _reference_g(q, r_grid, contrast, aspect_ratio):
        """G matrix with the full 50-node average — the pre-tier behaviour."""
        from pyirena.core.form_factors import (
            _GL_NODES,
            _GL_WEIGHTS,
            _sphere_amplitude,
        )
        AR = float(aspect_ratio)
        cos_t = 0.5 * (_GL_NODES + 1.0)
        weights = 0.5 * _GL_WEIGHTS
        G = np.empty((len(q), len(r_grid)))
        for j, r in enumerate(r_grid):
            V = (4.0 / 3.0) * np.pi * r ** 3 * AR
            r_eff = r * np.sqrt(1.0 + (AR ** 2 - 1.0) * cos_t ** 2)
            F = _sphere_amplitude(q[:, np.newaxis] * r_eff[np.newaxis, :])
            G[:, j] = V * (F ** 2 @ weights)
        return G * (contrast * 1e-4)

    # Spans every tier boundary, prolate and oblate, plus AR beyond the last
    # tier (which must keep all 50 nodes and so be bit-exact).
    @pytest.mark.parametrize("aspect_ratio", [
        1.0001, 0.9999, 1.02, 0.98, 1.05, 0.95, 1.1, 0.9, 1.2, 0.8,
        1.25, 0.75, 1.5, 0.5, 3.0, 5.0,
    ])
    def test_spheroid_quadrature_matches_50_nodes(self, aspect_ratio):
        """Reduced-order average agrees with the 50-node reference to 1e-8.

        The q and r ranges are deliberately wide (q*r up to ~5000) so the
        sphere form factor's minima are fully resolved — that is where a
        too-coarse orientation average shows up first.
        """
        q = np.logspace(-4, 0.0, 400)
        r_grid = np.concatenate([np.linspace(5, 300, 80),
                                 np.linspace(300, 5000, 80)])
        got = build_g_matrix(q, r_grid, shape="spheroid", contrast=100.0,
                             aspect_ratio=aspect_ratio)
        ref = self._reference_g(q, r_grid, 100.0, aspect_ratio)
        assert np.max(np.abs(got - ref)) / np.max(np.abs(ref)) < 1e-8

    def test_spheroid_ar_one_is_exactly_sphere(self):
        """AR = 1 must return the sphere G matrix itself, not an average of it.

        A spheroid of AR = 1 is a sphere, and the GUI offers both; picking
        "Spheroid" and pinning AR = 1 used to cost ~40x a "Sphere" population
        for a bit-identical answer.
        """
        q = np.logspace(-4, -0.5, 300)
        r_grid = np.linspace(5.0, 1000.0, 150)
        g_spheroid = build_g_matrix(q, r_grid, shape="spheroid",
                                    contrast=100.0, aspect_ratio=1.0)
        g_sphere = build_g_matrix(q, r_grid, shape="sphere", contrast=100.0)
        np.testing.assert_array_equal(g_spheroid, g_sphere)

    def test_spheroid_ff_ar_one_is_exactly_sphere(self):
        """Same short-circuit on the single-curve spheroid_ff."""
        q = np.logspace(-4, -0.5, 300)
        np.testing.assert_array_equal(spheroid_ff(q, 120.0, aspect_ratio=1.0),
                                      sphere_ff(q, 120.0))

    @pytest.mark.parametrize("aspect_ratio", [1.01, 1.05, 1.15, 1.3, 2.0])
    def test_core_shell_spheroid_keeps_all_50_nodes(self, aspect_ratio):
        """Core-shell spheroid is deliberately NOT reduced for AR != 1.

        Its amplitude carries two length scales (R_core and R_total), so the
        integrand oscillates faster than the plain spheroid's and the plain
        tiers are not accurate enough (~3e-7 at AR = 1.05). It passes
        ``tiers=()``, so results must be bit-for-bit unchanged.
        """
        from pyirena.core.form_factors import (
            _GL_NODES,
            _GL_WEIGHTS,
            _cs_spheroid_g_from_pairs,
            _shell_step_amplitude,
        )
        q = np.logspace(-4, 0.0, 300)
        r_c = np.linspace(10.0, 500.0, 120)
        r_t = r_c + 20.0
        AR = float(aspect_ratio)

        cos_a = 0.5 * (_GL_NODES + 1.0)
        weights = 0.5 * _GL_WEIGHTS
        stretch = np.sqrt(1.0 + (AR ** 2 - 1.0) * cos_a ** 2)
        acc = np.zeros((len(q), len(r_c)))
        for w_k, s_k in zip(weights, stretch):
            F = _shell_step_amplitude(
                q[:, np.newaxis] * s_k,
                (r_c[np.newaxis, :], r_t[np.newaxis, :]),
                (10.0 - 1.0, 1.0 - 9.46), vol_factor=AR)
            acc += w_k * F * F
        ref = acc / ((4.0 / 3.0) * np.pi * r_t[np.newaxis, :] ** 3 * AR) * 1e-4

        got = _cs_spheroid_g_from_pairs(q, r_c, r_t, 10.0, 1.0, 9.46, AR)
        np.testing.assert_array_equal(got, ref)
