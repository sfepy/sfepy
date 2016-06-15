# c: 07.05.2007, r: 25.06.2008
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

coef = 2.0
materials = {
    'coef' : ({'val' : coef},),
    'rhs' : 'rhs',
}

equations = {
    'Temperature' :
    """dw_laplace.i.Omega( coef.val, s, t )
       = - dw_volume_lvf.i.Omega( rhs.val, s )""",
}

solutions = {
    'sincos' : ('t', 'sin( 3.0 * x ) * cos( 4.0 * y )',
                '-25.0 * %s * sin( 3.0 * x ) * cos( 4.0 * y )' % coef),
    'poly' : ('t', '(x**2) + (y**2)', '4.0 * %s' % coef),
    'polysin' : ('t', '((x - 0.5)**3) * sin( 5.0 * y )',
                 '%s * (6.0 * (x - 0.5) * sin( 5.0 * y ) - 25.0 * ((x - 0.5)**3) * sin( 5.0 * y ))' % coef),
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
from sfepy.base.testing import TestCommon
from sfepy.base.base import debug, pause, assert_
output_name = 'test_msm_laplace_%s.vtk'

##
# c: 07.05.2007, r: 09.05.2008
solution = ['']
def ebc(ts, coor, **kwargs):
    expression = solution[0]
    val = TestCommon.eval_coor_expression( expression, coor )
    return nm.atleast_1d( val )

def rhs(ts, coor, mode=None, expression=None, **kwargs):
    if mode == 'qp':
        if expression is None:
            expression = '0.0 * x'

        val = TestCommon.eval_coor_expression( expression, coor )
        val.shape = (val.shape[0], 1, 1)

        return {'val' : val}

functions = {
    'ebc' : (ebc,),
    'rhs' : (rhs,),
}

##
# c: 07.05.2008
class Test( TestCommon ):

    ##
    # c: 07.05.2007, r: 07.05.2008
    def from_conf( conf, options ):
        from sfepy.discrete import Problem

        problem = Problem.from_conf(conf)
        test = Test( problem = problem,
                     conf = conf, options = options )
        return test
    from_conf = staticmethod( from_conf )

    ##
    # c: 09.05.2007, r: 25.06.2008
    def _build_rhs( self, sols ):
        for sol in six.itervalues(sols):
            assert_( len( sol ) == 3 )
        return sols

    ##
    # c: 07.05.2007, r: 09.05.2008
    def test_msm_laplace( self ):
        import os.path as op
        problem  = self.problem

        variables = problem.get_variables()
        materials = problem.get_materials()

        sols = self._build_rhs( self.conf.solutions )

        ok = True
        for sol_name, sol in six.iteritems(sols):
            self.report( 'testing', sol_name )
            var_name, sol_expr, rhs_expr = sol

            self.report( 'sol:', sol_expr )
            self.report( 'rhs:', rhs_expr )
            globals()['solution'][0] = sol_expr
            materials['rhs'].function.set_extra_args(expression=rhs_expr)
            problem.time_update()
            state = problem.solve()
            coor = variables[var_name].field.get_coor()
            ana_sol = self.eval_coor_expression( sol_expr, coor )
            num_sol = state(var_name)

            ana_norm = nm.linalg.norm( ana_sol, nm.inf )
            ret = self.compare_vectors( ana_sol, num_sol,
                                       allowed_error = ana_norm * 1e-2,
                                       label1 = 'analytical %s' % var_name,
                                       label2 = 'numerical %s' % var_name,
                                       norm = nm.inf )
            if not ret:
                self.report( 'variable %s: failed' % var_name )

            fname = op.join( self.options.out_dir, self.conf.output_name )
            out = {}
            astate = state.copy()
            astate.set_full(ana_sol)
            aux = astate.create_output_dict()
            out['ana_t'] = aux['t']
            aux = state.create_output_dict()
            out['num_t'] = aux['t']

            problem.domain.mesh.write( fname % sol_name, io = 'auto', out = out )

            ok = ok and ret

        return ok
