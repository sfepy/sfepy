"""
Test computing volumes by volume or surface integrals.
"""
import numpy as nm
import pytest

from sfepy import data_dir
import sfepy.base.testing as tst

filename_mesh = data_dir + '/meshes/3d/elbow.mesh'

fields = {
    'scalar' : ('real', 'scalar', 'Omega', 1),
    'vector' : ('real', 'vector', 'Omega', 1),
}

integrals = {
    'i' : 2,
}

regions = {
    'Omega' : 'all',
    'Gamma' : ('vertices of surface', 'facet'),
}

expressions = {
    'volume_p' : 'ev_volume.i.Omega(p)',
    'volume_u' : 'ev_volume.i.Omega(u)',
    'surface_p' : 'ev_volume_surface.i.Gamma(p)',
    'surface_u' : 'ev_volume_surface.i.Gamma(u)',
}

@pytest.fixture(scope='module')
def problem():
    import sys
    from sfepy.discrete import Problem
    from sfepy.base.conf import ProblemConf

    conf = ProblemConf.from_dict(globals(), sys.modules[__name__])

    problem = Problem.from_conf(conf, init_equations=False)
    return problem

def test_volume(problem):
    from sfepy.discrete import FieldVariable

    ok = True

    field_map = {'u' : 'vector', 'p' : 'scalar'}

    volumes = {}
    avg = 0.0
    for key, term in expressions.items():
        var_name = key[-1]
        field = problem.fields[field_map[var_name]]
        var = FieldVariable(var_name, 'parameter', field,
                            primary_var_name='(set-to-None)')

        val = problem.evaluate(term, **{var_name : var})

        volumes[key] = val
        avg += val

    avg /= len(volumes)

    for key, val in volumes.items():
        err = nm.abs(avg - val) / nm.abs(avg)
        _ok = err < 1e-12
        tst.report('"'"%s"'" - volume: %e' % (key, val))
        tst.report('"'"%s"'" - relative volume difference: %e -> %s'
                    % (key, err, _ok))
        ok = ok and _ok

    assert ok

def test_volume_tl(problem):
    from sfepy.discrete import FieldVariable

    fu = problem.fields['vector']
    fq = problem.fields['scalar']

    var_u = FieldVariable('u', 'parameter', fu,
                          primary_var_name='(set-to-None)')
    var_q = FieldVariable('q', 'test', fq,
                          primary_var_name='(set-to-None)')

    var_u.set_data(nm.linspace(0, 0.004, var_u.n_dof))

    vval = problem.evaluate('dw_tl_volume.i.Omega( q, u )',
                            term_mode='volume', q=var_q, u=var_u)

    sval = problem.evaluate('ev_tl_volume_surface.i.Gamma( u )',
                            u=var_u)

    ok = abs(vval - sval) < 1e-14

    tst.report('TL: by volume: %e == by surface: %e -> %s' %
                (vval, sval, ok))

    assert ok
