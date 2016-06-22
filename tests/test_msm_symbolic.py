from __future__ import absolute_import
from sfepy import data_dir
import six

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

import numpy as nm
try:
    import sfepy.linalg.sympy_operators as sops
except ImportError as exc:
    sops = None

from sfepy.base.testing import TestCommon

output_name = 'test_msm_symbolic_%s.vtk'

solution = ['']
def ebc(ts, coor, solution=None):
    expression = solution[0]
    val = TestCommon.eval_coor_expression(expression, coor)
    return nm.atleast_1d(val)

def rhs(ts, coor, mode=None, expression=None, **kwargs):
    if mode == 'qp':
        if expression is None:
            expression = '0.0 * x'

        val = TestCommon.eval_coor_expression(expression, coor)
        val.shape = (val.shape[0], 1, 1)
        return {'val' : val}

functions = {
    'ebc' : (lambda ts, coor, **kwargs:
             ebc(ts, coor, solution=solution),),
    'rhs' : (rhs,),
}

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        from sfepy.discrete import Problem

        problem = Problem.from_conf(conf, init_equations=False)
        test = Test(problem=problem, conf=conf, options=options)
        return test

    def _build_rhs(self, equation, sols):
        rhss = {}
        self.report('%s:' % equation.name)
        self.report('evaluating terms, "<=" is solution, "=>" is the rhs:')
        for term in equation.terms:
            if not hasattr(term, 'symbolic'):
                self.report('term %s has no symbolic description!' % term.name)
                raise ValueError
            expr = term.symbolic['expression']
            arg_map = term.symbolic['map']
            self.report('%s(%s)' %\
                         (term.name, ', '.join(term.ats)))
            self.report('multiplicator: %f' % term.sign)
            self.report('  symbolic:', expr)
            self.report('  using argument map:', arg_map)

            for sol_name, sol in six.iteritems(sols):
                rhs = self._eval_term(sol[1], term, sops)
                srhs = "(%s * (%s))" % (term.sign, rhs)
                rhss.setdefault(sol_name, []).append(srhs)

        for key, val in six.iteritems(rhss):
            rhss[key] = '+'.join(val)

        return rhss

    def _eval_term(self, sol, term, sops):
        """Works for scalar, single unknown terms only!"""
        expr = term.symbolic['expression']
        arg_map = term.symbolic['map']
        env = {'x' : sops.Symbol('x'),
               'y' : sops.Symbol('y'),
               'z' : sops.Symbol('z'),
               'dim' : dim}
        for key, val in six.iteritems(arg_map):
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

        self.report('  <= ', sol)
        sops.set_dim(dim)
        val = str(eval(expr, sops.__dict__, env))
        self.report('   =>', val)
        return val

    def _test_msm_symbolic(self, equations):
        import os.path as op

        if sops is None:
            self.report('cannot import sympy, skipping')
            return True

        problem  = self.problem
        ok = True
        for eq_name, equation in six.iteritems(equations):
            problem.set_equations({eq_name : equation})
            problem.update_materials()

            rhss = self._build_rhs(problem.equations[eq_name],
                                   self.conf.solutions)
            erhs = problem.conf.equations_rhs[eq_name]

            problem.set_equations({eq_name : equation + erhs})
            variables = problem.get_variables()
            materials = problem.get_materials()
            rhs_mat = materials['rhs']

            for sol_name, sol in six.iteritems(problem.conf.solutions):
                self.report('testing', sol_name)
                var_name, sol_expr = sol
                rhs_expr = rhss[sol_name]

                self.report('sol:', sol_expr)
                self.report('rhs:', rhs_expr)
                globals()['solution'][0] = sol_expr
                rhs_mat.function.set_extra_args(expression=rhs_expr)
                problem.time_update()
                state = problem.solve()
                coor = variables[var_name].field.get_coor()
                ana_sol = self.eval_coor_expression(sol_expr, coor)
                num_sol = state(var_name)

                ana_norm = nm.linalg.norm(ana_sol, nm.inf)
                ret = self.compare_vectors(ana_sol, num_sol,
                                           allowed_error=ana_norm * 1e-2,
                                           label1='analytical %s' % var_name,
                                           label2='numerical %s' % var_name,
                                           norm=nm.inf)
                if not ret:
                    self.report('variable %s: failed' % var_name)

                fname = op.join(self.options.out_dir, self.conf.output_name)
                out = {}
                astate = state.copy()
                astate.set_full(ana_sol)
                aux = astate.create_output_dict()
                out['ana_t'] = aux['t']
                aux = state.create_output_dict()
                out['num_t'] = aux['t']

                problem.domain.mesh.write(fname % '_'.join((sol_name, eq_name)),
                                          io='auto', out=out)

                ok = ok and ret

        return ok

    def _get_equations(self, name):
        """Choose a sub-problem from all equations."""
        return {name : self.problem.conf.equations[name]}

    def test_msm_symbolic_laplace(self):
        return self._test_msm_symbolic(self._get_equations('Laplace'))

    def test_msm_symbolic_diffusion(self):
        return self._test_msm_symbolic(self._get_equations('Diffusion'))
