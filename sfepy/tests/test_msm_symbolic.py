import numpy as nm
import pytest

try:
    import sfepy.linalg.sympy_operators as sops
except ImportError:
    sops = None

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

material_1 = {
    'name' : 'coef',
    'values' : {
        'val' : 12.0,
        'K' : [[1.0, 0.3], [0.3, 2.0]],
    }
}

material_2 = {
    'name' : 'rhs',
    'function' : 'rhs',
}

equations = {
    'Laplace' :
    """2 * dw_laplace.i.Omega(coef.val, s, t)
    """,
    'Diffusion' :
    """3 * dw_diffusion.i.Omega(coef.K, s, t)
    """,
}
equations_rhs = {
    'Laplace' :
    """= - dw_volume_lvf.i.Omega(rhs.val, s)""",
    'Diffusion' :
    """= - dw_volume_lvf.i.Omega(rhs.val, s)""",
}

solutions = {
    'sincos' : ('t', 'sin(3.0 * x) * cos(4.0 * y)'),
    'poly' : ('t', '(x**2) + (y**2)'),
    'polysin' : ('t', '((x - 0.5)**3) * sin(5.0 * y)'),
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
def ebc(ts, coor, solution=None):
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
    'ebc' : (lambda ts, coor, **kwargs:
             ebc(ts, coor, solution=solution),),
    'rhs' : (rhs,),
}

def _build_rhs(equation, sols):
    rhss = {}
    tst.report('%s:' % equation.name)
    tst.report('evaluating terms, "<=" is solution, "=>" is the rhs:')
    for term in equation.terms:
        if not hasattr(term, 'symbolic'):
            tst.report('term %s has no symbolic description!' % term.name)
            raise ValueError
        expr = term.symbolic['expression']
        arg_map = term.symbolic['map']
        tst.report('%s(%s)' % (term.name, ', '.join(term.ats)))
        tst.report('multiplicator: %f' % term.sign)
        tst.report('  symbolic:', expr)
        tst.report('  using argument map:', arg_map)

        for sol_name, sol in sols.items():
            rhs = _eval_term(sol[1], term, sops)
            srhs = "(%s * (%s))" % (term.sign, rhs)
            rhss.setdefault(sol_name, []).append(srhs)

    for key, val in rhss.items():
        rhss[key] = '+'.join(val)

    return rhss

def _eval_term(sol, term, sops):
    """Works for scalar, single unknown terms only!"""
    expr = term.symbolic['expression']
    arg_map = term.symbolic['map']
    env = {'x' : sops.Symbol('x'),
           'y' : sops.Symbol('y'),
           'z' : sops.Symbol('z'),
           'dim' : dim}
    for key, val in arg_map.items():
        if val == 'state':
            env[key] = sol
        else:
            env[key] = term.get_args([val])[0]

        if 'material' in val:
            # Take the first value - constant in all QPs.
            aux = env[key][0,0]
            if nm.prod(aux.shape) == 1:
                env[key] = aux.squeeze()
            else:
                import sympy
                env[key] = sympy.Matrix(aux)

    tst.report('  <= ', sol)
    sops.set_dim(dim)
    val = str(eval(expr, sops.__dict__, env))
    tst.report('   =>', val)
    return val

def _test_msm_symbolic(problem, equations, output_dir):
    import os.path as op

    if sops is None:
        tst.report('cannot import sympy, skipping')
        return True

    ok = True
    for eq_name, equation in equations.items():
        problem.set_equations({eq_name : equation})
        problem.update_materials()

        rhss = _build_rhs(problem.equations[eq_name], solutions)
        erhs = problem.conf.equations_rhs[eq_name]

        problem.set_equations({eq_name : equation + erhs})
        variables = problem.get_variables()
        materials = problem.get_materials()
        rhs_mat = materials['rhs']

        for sol_name, sol in problem.conf.solutions.items():
            tst.report('testing', sol_name)
            var_name, sol_expr = sol
            rhs_expr = rhss[sol_name]

            tst.report('sol:', sol_expr)
            tst.report('rhs:', rhs_expr)
            globals()['solution'][0] = sol_expr
            rhs_mat.function.set_extra_args(expression=rhs_expr)
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

            fname = op.join(output_dir, 'test_msm_symbolic_%s.vtk')
            out = {}
            astate = state.copy()
            astate.set_state(ana_sol)
            aux = astate.create_output()
            out['ana_t'] = aux['t']
            aux = state.create_output()
            out['num_t'] = aux['t']

            problem.domain.mesh.write(fname % '_'.join((sol_name, eq_name)),
                                      io='auto', out=out)

            ok = ok and ret

    return ok

@pytest.fixture(scope='module')
def problem():
    import sys
    from sfepy.discrete import Problem
    from sfepy.base.conf import ProblemConf

    conf = ProblemConf.from_dict(globals(), sys.modules[__name__])

    problem = Problem.from_conf(conf)
    return problem

def _get_equations(problem, name):
    """
    Choose a sub-problem from all equations.
    """
    return {name : problem.conf.equations[name]}

def test_msm_symbolic_laplace(problem, output_dir):
    equations = _get_equations(problem, 'Laplace')
    assert _test_msm_symbolic(problem, equations, output_dir)

def test_msm_symbolic_diffusion(problem, output_dir):
    equations = _get_equations(problem, 'Diffusion')
    assert _test_msm_symbolic(problem, equations, output_dir)
