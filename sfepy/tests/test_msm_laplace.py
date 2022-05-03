import numpy as nm
import pytest

from sfepy.base.base import assert_
from sfepy import data_dir
import sfepy.base.testing as tst

filename_mesh = data_dir + '/meshes/2d/special/circle_in_square.mesh'

dim = 2

field_1 = {
    'name' : 'a_harmonic_field',
    'dtype' : 'real',
    'shape' : 'scalar',
    'region' : 'Omega',
    'approx_order' : 1,
}

variables = {
    't': ('unknown field', 'a_harmonic_field', 0),
    's': ('test field',    'a_harmonic_field', 't'),
}

regions = {
    'Omega' : 'all',
    'Gamma' : ('vertices of surface', 'facet'),
}

ebcs = {
    't_left' : ('Gamma', {'t.0' : 'ebc'}),
}

integral_1 = {
    'name' : 'i',
    'order' : 2,
}

coef = 2.0
materials = {
    'coef' : ({'val' : coef},),
    'rhs' : 'rhs',
}

equations = {
    'Temperature' :
    """dw_laplace.i.Omega(coef.val, s, t)
       = - dw_volume_lvf.i.Omega(rhs.val, s)""",
}

solutions = {
    'sincos' : ('t', 'sin(3.0 * x) * cos(4.0 * y)',
                '-25.0 * %s * sin(3.0 * x) * cos(4.0 * y)' % coef),
    'poly' : ('t', '(x**2) + (y**2)', '4.0 * %s' % coef),
    'polysin' : ('t', '((x - 0.5)**3) * sin(5.0 * y)',
                 '%s * (6.0 * (x - 0.5) * sin(5.0 * y)'
                 ' - 25.0 * ((x - 0.5)**3) * sin(5.0 * y))' % coef),
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 1,
    'eps_a'      : 1e-10,
}

solution = ['']
def ebc(ts, coor, **kwargs):
    expression = solution[0]
    val = tst.eval_coor_expression(expression, coor)
    return nm.atleast_1d(val)

def rhs(ts, coor, mode=None, expression=None, **kwargs):
    if mode == 'qp':
        if expression is None:
            expression = '0.0 * x'

        val = tst.eval_coor_expression(expression, coor)
        val.shape = (val.shape[0], 1, 1)

        return {'val' : val}

functions = {
    'ebc' : (ebc,),
    'rhs' : (rhs,),
}

@pytest.fixture(scope='module')
def problem():
    import sys
    from sfepy.discrete import Problem
    from sfepy.base.conf import ProblemConf

    conf = ProblemConf.from_dict(globals(), sys.modules[__name__])

    problem = Problem.from_conf(conf)
    return problem

def _build_rhs(sols):
    for sol in sols.values():
        assert_(len(sol) == 3)
    return sols

def test_msm_laplace(problem, output_dir):
    import os.path as op

    variables = problem.get_variables()
    materials = problem.get_materials()

    sols = _build_rhs(solutions)

    ok = True
    for sol_name, sol in sols.items():
        tst.report('testing', sol_name)
        var_name, sol_expr, rhs_expr=sol

        tst.report('sol:', sol_expr)
        tst.report('rhs:', rhs_expr)
        globals()['solution'][0] = sol_expr
        materials['rhs'].function.set_extra_args(expression=rhs_expr)
        problem.time_update()
        state = problem.solve(save_results=False)
        coor = variables[var_name].field.get_coor()
        ana_sol = tst.eval_coor_expression(sol_expr, coor)
        num_sol = state(var_name)

        ana_norm = nm.linalg.norm(ana_sol, nm.inf)
        ret = tst.compare_vectors(ana_sol, num_sol,
                                  allowed_error=ana_norm * 1e-2,
                                  label1='analytical %s' % var_name,
                                  label2='numerical %s' % var_name,
                                  norm=nm.inf)
        if not ret:
            tst.report('variable %s: failed' % var_name)

        fname = op.join(output_dir, 'test_msm_laplace_%s.vtk')
        out = {}
        astate = state.copy()
        astate.set_state(ana_sol)
        aux = astate.create_output()
        out['ana_t'] = aux['t']
        aux = state.create_output()
        out['num_t'] = aux['t']

        problem.domain.mesh.write(fname % sol_name, io='auto', out=out)

        ok = ok and ret

    assert_(ok)
