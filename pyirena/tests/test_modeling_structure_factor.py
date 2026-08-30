"""Structure-factor parameters must actually be fitted by the Modeling engine.

Regression test for a bug found while building ``validationData/``: the
parameter-key group for a size-distribution population's structure-factor
parameters was ``'sf'``, the same tag the *surface fractal* population uses for
its own attributes.  ``_unpack_params`` tested the surface-fractal group first,
so every optimiser step wrote the hard-sphere radius and volume fraction onto
the population object as stray attributes instead of into ``pop.sf_params``.
The model therefore never saw them change: they were packed into the parameter
vector, contributed nothing to chi-squared, and came back from ``fit()`` at
exactly their starting values.

The structure-factor group is now ``'sfp'``.  These tests pin both halves of
that contract: the key group, and the recovery of hard-sphere parameters from
data generated with a known Percus-Yevick structure factor.
"""

import numpy as np

from pyirena.core.modeling import (
    ModelingConfig,
    ModelingEngine,
    SizeDistPopulation,
    _hard_sphere_sf,
)


def _hard_sphere_pop(radius=50.0, volume_fraction=0.25, mean_size=50.0,
                     sdeviation=0.12, scale=0.20):
    pop = SizeDistPopulation()
    pop.dist_type = 'lognormal'
    pop.dist_params = {'min_size': 0.0, 'mean_size': mean_size,
                       'sdeviation': sdeviation}
    pop.dist_params_fit = {'min_size': False, 'mean_size': True,
                           'sdeviation': True}
    pop.form_factor = 'sphere'
    pop.contrast = 100.0
    pop.fit_contrast = False
    pop.scale = scale
    pop.fit_scale = True
    pop.n_bins = 120
    pop.structure_factor = 'hard_sphere'
    pop.sf_params = dict(pop.sf_params)
    pop.sf_params['radius'] = radius
    pop.sf_params['volume_fraction'] = volume_fraction
    pop.sf_params_fit = dict(pop.sf_params_fit)
    pop.sf_params_fit['radius'] = True
    pop.sf_params_fit['volume_fraction'] = True
    return pop


def _config(pop, background=0.0):
    cfg = ModelingConfig()
    cfg.populations = [pop]
    cfg.background = background
    cfg.fit_background = False
    cfg.q_min, cfg.q_max = 1e-4, 10.0
    return cfg


def test_structure_factor_params_use_their_own_key_group():
    """Packed under 'sfp' — never 'sf', which belongs to the surface fractal."""
    eng = ModelingEngine()
    cfg = _config(_hard_sphere_pop())
    _, _, _, keys = eng._pack_params(cfg)
    sf_keys = [k for k in keys if k[0] == 'pop' and k[2] in ('sf', 'sfp')]
    assert {k[3] for k in sf_keys} == {'radius', 'volume_fraction'}
    assert all(k[2] == 'sfp' for k in sf_keys)


def test_unpack_writes_structure_factor_params_into_sf_params():
    """The round trip must land in pop.sf_params, not on stray attributes."""
    eng = ModelingEngine()
    cfg = _config(_hard_sphere_pop(radius=50.0, volume_fraction=0.25))
    x0, _, _, keys = eng._pack_params(cfg)
    x = list(x0)
    for i, k in enumerate(keys):
        if k[0] == 'pop' and k[2] == 'sfp' and k[3] == 'radius':
            x[i] = 77.0
        if k[0] == 'pop' and k[2] == 'sfp' and k[3] == 'volume_fraction':
            x[i] = 0.33
    eng._unpack_params(np.asarray(x, float), keys, cfg)
    assert cfg.populations[0].sf_params['radius'] == 77.0
    assert cfg.populations[0].sf_params['volume_fraction'] == 0.33


def test_structure_factor_changes_the_model_intensity():
    """A different hard-sphere volume fraction must change I(Q)."""
    eng = ModelingEngine()
    q = np.logspace(-2.5, -0.3, 200)
    i_low, *_ = eng.total_intensity(_config(_hard_sphere_pop(volume_fraction=0.05)),
                                    q, use_cache=False)
    i_high, *_ = eng.total_intensity(_config(_hard_sphere_pop(volume_fraction=0.40)),
                                     q, use_cache=False)
    assert np.max(np.abs(i_low - i_high)) / np.max(i_low) > 0.1


def test_hard_sphere_parameters_are_recovered_from_synthetic_data():
    """Fit data built with a known S(Q) and check the parameters come back."""
    from pyirena.core.form_factors import bin_widths, build_g_matrix

    q = np.logspace(np.log10(3e-3), np.log10(0.5), 300)
    r = np.logspace(np.log10(20.0), np.log10(120.0), 600)
    fv = np.exp(-np.log(r / 50.0) ** 2 / (2 * 0.12 ** 2)) / (r * 0.12 * np.sqrt(2 * np.pi))
    fv *= 0.20 / np.trapezoid(fv, r)
    G = build_g_matrix(q, r, 'sphere', 100.0)
    I = (G @ (fv * bin_widths(r))) * _hard_sphere_sf(q, 50.0, 0.25)
    dI = 0.01 * I

    # Start well away from the truth so a no-op parameter cannot pass by luck.
    pop = _hard_sphere_pop(radius=40.0, volume_fraction=0.15,
                           mean_size=40.0, sdeviation=0.08, scale=0.12)
    res = ModelingEngine().fit(_config(pop), q, I, dI)
    got = res.config.populations[0]
    np.testing.assert_allclose(got.sf_params['radius'], 50.0, rtol=0.10)
    np.testing.assert_allclose(got.sf_params['volume_fraction'], 0.25, rtol=0.10)
    np.testing.assert_allclose(got.dist_params['mean_size'], 50.0, rtol=0.05)
    np.testing.assert_allclose(got.scale, 0.20, rtol=0.10)
