import pytest

from sfepy import data_dir
import sfepy.base.testing as tst

filename_mesh = data_dir + '/meshes/2d/special/circle_in_square.mesh'

dim = 2

fields = {
    'scalar_field' : ('real', 'scalar', 'Omega', 1),
    'vector_field' : ('real', 'vector', 'Omega', 1),
}

variables = {
    'us'  : ('unknown field',   'scalar_field', 0),
    'ts'  : ('test field',      'scalar_field', 'us'),
    'ps1' : ('parameter field', 'scalar_field', 'us'),
    'ps2' : ('parameter field', 'scalar_field', 'us'),
    'uv'  : ('unknown field',   'vector_field', 1),
    'tv'  : ('test field',      'vector_field', 'uv'),
    'pv1' : ('parameter field', 'vector_field', 'uv'),
    'pv2' : ('parameter field', 'vector_field', 'uv'),
}

regions = {
    'Omega' : 'all',
    'Left' : ('vertices in (x < -0.499)', 'facet'),
}

integrals = {
    'i' : 2,
}

materials = {
    'm' : 'get_pars',
    'm2' : ({'K' : [[3.0, 0.1], [0.3, 1.0]]},),
}

equations = {
    'eq' : """dw_diffusion.i.Omega(m2.K, ts, us) = 0"""
}

def get_pars(ts, coor, mode=None, term=None, **kwargs):
    if mode == 'qp':
        n_nod, dim = coor.shape
        sym = (dim + 1) * dim // 2

        if 'biot' in term.name:
            val = nm.zeros((sym, 1), dtype=nm.float64)
            val[:dim] = 0.132
            val[dim:sym] = 0.092
        elif '_dot' in term.name:
            val = 1.0 / nm.array([3.8], dtype=nm.float64)
        elif 'diffusion' in term.name:
            val = nm.eye(dim, dtype=nm.float64)
        else:
            raise ValueError

        return {'val' : nm.tile(val, (coor.shape[0], 1, 1))}

functions = {
    'get_pars' : (get_pars,),
}

# (eval term prefix, parameter corresponding to test variable, 'd' variables,
# 'dw' variables (test must be paired with unknown, which should be at
# index 2!), mat mode)
test_terms = [
    ('%s_biot.i.Omega(m.val, %s, %s)',
     ('dw', 'ps1', ('pv1', 'ps1'), ('pv1', 'ts', 'us', 'uv', 'tv'))),
    ('%s_biot.i.Omega(m.val, %s, %s)',
     ('dw', 'pv1', ('pv1', 'ps1'), ('tv', 'ps1', 'uv', 'us', 'ts'))),
    ('%s_diffusion.i.Omega(m.val, %s, %s)',
     ('dw', 'ps1', ('ps1', 'ps2'), ('ts', 'ps1', 'us'))),
    ('%s_dot.i.Omega(m.val, %s, %s)',
     ('dw', 'ps1', ('ps1', 'ps2'), ('ts', 'ps1', 'us'))),
]

import numpy as nm

def _integrate(var, val_qp):
    from sfepy.discrete import Integral
    from sfepy.discrete.common.mappings import get_jacobian

    integral = Integral('i', 2)
    det = get_jacobian(var.field, integral)
    val = (val_qp * det).sum(axis=1) / det.sum(axis=1)

    return val

@pytest.fixture(scope='module')
def problem():
    import sys
    from sfepy.discrete import Problem
    from sfepy.base.conf import ProblemConf

    conf = ProblemConf.from_dict(globals(), sys.modules[__name__])

    problem = Problem.from_conf(conf, init_equations=False)
    return problem

def test_consistency_d_dw(problem):
    from sfepy.discrete import Variables

    ok = True
    for aux in test_terms:
        term_template, (prefix, par_name, d_vars, dw_vars) = aux
        tst.report(term_template, prefix, par_name, d_vars, dw_vars)

        term1 = term_template % ((prefix,) + d_vars)

        variables = Variables.from_conf(problem.conf.variables, problem.fields)

        for var_name in d_vars:
            var = variables[var_name]
            n_dof = var.field.n_nod * var.field.shape[0]
            aux = nm.arange(n_dof, dtype=nm.float64)
            var.set_data(aux)

        if prefix == 'd':
            val1 = problem.evaluate(term1, var_dict=variables.as_dict())

        else:
            val1 = problem.evaluate(term1, call_mode='d_eval',
                                    var_dict=variables.as_dict())

        tst.report('%s: %s' % (term1, val1))

        term2 = term_template % (('dw',) + dw_vars[:2])

        vec, vv = problem.evaluate(term2, mode='weak',
                                   var_dict=variables.as_dict(),
                                   ret_variables=True)

        pvec = vv.get_vec_part(vec, dw_vars[2])
        val2 = nm.dot(variables[par_name](), pvec)
        tst.report('%s: %s' % (term2, val2))

        err = nm.abs(val1 - val2) / nm.abs(val1)
        _ok = err < 1e-12
        tst.report('relative difference: %e -> %s' % (err, _ok))

        ok = ok and _ok

    assert ok

def test_eval_matrix(problem):
    problem.set_equations()
    problem.time_update(ebcs={}, epbcs={})

    var = problem.get_variables()['us']

    vec = nm.arange(var.n_dof, dtype=var.dtype)
    var.set_data(vec)

    val1 = problem.evaluate('dw_diffusion.i.Omega(m2.K, us, us)',
                            mode='eval', us=var)

    mtx = problem.evaluate('dw_diffusion.i.Omega(m2.K, ts, us)',
                           mode='weak', dw_mode='matrix')

    val2 = nm.dot(vec, mtx * vec)

    ok = (nm.abs(val1 - val2) / nm.abs(val1)) < 1e-14
    tst.report('eval: %s, weak: %s, ok: %s' % (val1, val2, ok))

    assert ok

def test_vector_matrix(problem):
    problem.set_equations()
    problem.time_update()

    state = problem.create_state()
    state.apply_ebc()

    aux1 = problem.evaluate("dw_diffusion.i.Omega(m2.K, ts, us)",
                            mode='weak', dw_mode='vector')

    problem.time_update(ebcs={}, epbcs={})

    mtx = problem.evaluate("dw_diffusion.i.Omega(m2.K, ts, us)",
                           mode='weak', dw_mode='matrix')
    aux2g = mtx * state()
    problem.time_update(ebcs=problem.conf.ebcs,
                        epbcs=problem.conf.epbcs)
    aux2 = problem.equations.reduce_vec(aux2g, follow_epbc=True)

    ok = tst.compare_vectors(aux1, aux2,
                             label1='vector mode',
                             label2='matrix mode')
    if not ok:
        tst.report('failed')

    assert ok

def test_surface_evaluate(problem):
    from sfepy.discrete import FieldVariable

    us = problem.create_variables(['us'])['us']
    vec = nm.ones(us.n_dof, dtype=us.dtype)
    us.set_data(vec)

    expr = 'ev_integrate.i.Left(us)'
    val = problem.evaluate(expr, us=us)
    ok1 = nm.abs(val - 1.0) < 1e-15
    tst.report('with unknown: %s, value: %s, ok: %s'
               % (expr, val, ok1))

    ps1 = FieldVariable('ps1', 'parameter', us.get_field(),
                        primary_var_name='(set-to-None)')
    ps1.set_data(vec)

    expr = 'ev_integrate.i.Left(ps1)'
    val = problem.evaluate(expr, ps1=ps1)
    ok2 = nm.abs(val - 1.0) < 1e-15
    tst.report('with parameter: %s, value: %s, ok: %s'
               % (expr, val, ok2))

    assert ok1 and ok2

def test_ev_grad(problem):
    var = problem.create_variables(['us'])['us']
    val = nm.arange(var.n_dof, dtype=var.dtype)
    var.set_data(val)

    val1 = problem.evaluate('ev_grad.i.Omega(us)', us=var, mode='el_avg')
    tst.report('ev_grad(el_avg): min, max:', val1.min(), val1.max())

    aux = problem.evaluate('ev_grad.i.Omega(us)', us=var, mode='qp')
    val2 = _integrate(var, aux)
    val2.shape = val1.shape
    tst.report('ev_grad(qp): min, max:', val2.min(), val2.max())

    ok = tst.compare_vectors(val1, val2,
                             label1='de mode',
                             label2='dq mode')
    if not ok:
        tst.report('failed')

    assert ok

def test_ev_div(problem):
    var = problem.create_variables(['uv'])['uv']
    val = nm.arange(var.n_dof, dtype=var.dtype)
    var.set_data(val)

    val1 = problem.evaluate('ev_div.i.Omega(uv)', uv=var, mode='el_avg')
    tst.report('ev_div(el_avg): min, max:', val1.min(), val1.max())

    aux = problem.evaluate('ev_div.i.Omega(uv)', uv=var, mode='qp')
    val2 = _integrate(var, aux)
    val2.shape = val1.shape
    tst.report('ev_div(qp): min, max:', val2.min(), val2.max())

    ok = tst.compare_vectors(val1, val2,
                             label1='de mode',
                             label2='dq mode')
    if not ok:
        tst.report('failed')

    assert ok
