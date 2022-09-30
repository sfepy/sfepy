"""
Test computing volumes by volume or surface integrals.
"""
import numpy as nm
import pytest

from sfepy import data_dir
import sfepy.base.testing as tst

filename_mesh = data_dir + '/meshes/elements/3_8_1.mesh'

fields = {
    'scalar': ('real', 'scalar', 'Omega', 1),
    'vector': ('real', 'vector', 'Omega', 1),
}

integrals = {
    'i': 2,
}

regions = {
    'Omega': 'all',
}

expressions = {
    'grad_p': 'ev_grad.i.Omega(p)',
    'grad_p_de': 'de_grad.i.Omega(p)',
    'grad_u': 'ev_grad.i.Omega(u)',
    'grad_u_de': 'de_grad.i.Omega(u)',
    'def_grad_u': 'ev_def_grad.i.Omega(u)',
    'volume': 'ev_volume.i.Omega(p)',
}


@pytest.fixture(scope='module')
def problem():
    import sys
    from sfepy.discrete import Problem
    from sfepy.base.conf import ProblemConf

    conf = ProblemConf.from_dict(globals(), sys.modules[__name__])

    problem = Problem.from_conf(conf, init_equations=False)
    return problem


def test_gradients(problem):
    from sfepy.discrete import FieldVariable

    expected_grad = {
        'p': nm.array([[2., -3., 4.]]).T,
        'u': nm.array([[2., -3., 0.], [4., 5., -1.], [1., 0., 1.]]),
    }

    ok = True

    field_map = {'u': 'vector', 'p': 'scalar'}

    coors = problem.domain.mesh.coors
    ecoors = nm.average(coors, axis=0).reshape((1, -1))
    dim = coors.shape[1]

    values = {
        'u': nm.dot(coors, expected_grad['u'].T),
        'p': nm.dot(coors, expected_grad['p']),
    }

    volume = None
    outs = {}
    fields = {vn: problem.fields[field_map[vn]] for vn in values.keys()}
    vars = {vn: FieldVariable(vn, 'parameter', fields[vn],
                              primary_var_name='(set-to-None)')
            for vn in values.keys()}
    for key, term in expressions.items():
        vn = term[-2]
        var = vars[vn]

        if key == 'def_grad_u':
            var.set_data(values[vn] - coors)
        else:
            var.set_data(values[vn])

        val = problem.evaluate(term, **{vn: var})

        if key == 'volume':
            volume = val
        else:
            outs[f"{term.split('.')[0]}({vn})"] = val

    outs['evaluate_at(p)'] =\
        vars['p'].field.evaluate_at(ecoors, values['p'],
                                    mode='grad')[0, ...] * volume
    outs['evaluate_at(u)'] =\
        vars['u'].field.evaluate_at(ecoors, values['u'],
                                    mode='grad')[0, ...] * volume

    ok = True
    for k, v in outs.items():
        vn = k[-2]
        exp_v = expected_grad[vn]

        _ok = nm.allclose(exp_v, v, atol=1e-12) and exp_v.shape == v.shape

        tst.report(f'{k}:')
        tst.report(f'  expected values: {exp_v}')
        tst.report(f'  received values: {v}')
        tst.report(f'  difference: {nm.abs(exp_v - v)}')
        tst.report(f'  expected shape: {exp_v.shape}')
        tst.report(f'  received shape: {v.shape}')
        ok = ok and _ok

    assert ok
