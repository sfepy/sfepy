# c: 07.05.2007, r: 25.06.2008
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/special/circle_in_square.mesh'

dim = 2

field_1 = {
    'name' : 'a_harmonic_field',
    'dim' : (1,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '2_3_P1'}
}

variables = {
    't': ('unknown field', 'a_harmonic_field', 0),
    's': ('test field',    'a_harmonic_field', 't'),
}

regions = {
    'Omega' : ('all', {}),
    'Left' : ('nodes in (x < 0.001) & (y < 0.001)', {}),
    'Right' : ('nodes in (x > 0.999)', {}),
    'Gamma' : ('nodes of surface', {}),
}

ebcs = {
    't_left' : ('Gamma', {'t.0' : 'ebc'}),
#    't_right' : ('Right', {'t.0' : 'ebc'}),
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d2',
}

coef = 2.0
materials = {
    'coef' : ('Omega', {'val' : coef}),
    'rhs' : ('Omega', None, 'rhs'),
}

equations = {
    'Temperature' :
    """dw_laplace.i1.Omega( coef.val, s, t )
       = - dw_volume_lvf.i1.Omega( rhs.val, s )""",
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
    'eps_r'      : 1.0,
    'macheps'   : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'ls_red'     : 0.1,
    'ls_red_warp' : 0.001,
    'ls_on'      : 1.1,
    'ls_min'     : 1e-5,
    'check'     : 0,
    'delta'     : 1e-6,
    'is_plot'    : False,
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
}

fe = {
    'chunk_size' : 100000
}

import numpy as nm
from sfepy.base.testing import TestCommon
from sfepy.base.base import debug, pause, assert_
output_name = 'test_msm_laplace_%s.vtk'

##
# c: 07.05.2007, r: 09.05.2008
solution = ['']
def ebc(ts, coor, bc):
    expression = solution[0]
    val = TestCommon.eval_coor_expression( expression, coor )
    return nm.atleast_1d( val )

##
# c: 07.05.2007, r: 09.05.2008
def rhs(ts, coor, mode=None, region=None, ig=None, expression=None):
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
        from sfepy.fem import ProblemDefinition

        problem = ProblemDefinition.from_conf( conf )
        test = Test( problem = problem,
                     conf = conf, options = options )
        return test
    from_conf = staticmethod( from_conf )

    ##
    # c: 09.05.2007, r: 25.06.2008
    def _build_rhs( self, sols ):
        for sol in sols.itervalues():
            assert_( len( sol ) == 3 )
        return sols

    ##
    # c: 07.05.2007, r: 09.05.2008
    def test_msm_laplace( self ):
        import os.path as op
        problem  = self.problem

        sols = self._build_rhs( self.conf.solutions )

        ok = True
        for sol_name, sol in sols.iteritems():
            self.report( 'testing', sol_name )
            var_name, sol_expr, rhs_expr = sol

            self.report( 'sol:', sol_expr )
            self.report( 'rhs:', rhs_expr )
            globals()['solution'][0] = sol_expr
            problem.materials['rhs'].function.set_extra_args(expression=rhs_expr)
            problem.time_update()
            problem.equations.reset_term_caches()
            vec = problem.solve()
            coor = problem.variables[var_name].field.get_coor()
            ana_sol = self.eval_coor_expression( sol_expr, coor )
            num_sol = problem.variables.get_state_part_view( vec, var_name )

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
            aux = problem.state_to_output( ana_sol )
            out['ana_t'] = aux['t']
            aux = problem.state_to_output( num_sol )
            out['num_t'] = aux['t']

            problem.domain.mesh.write( fname % sol_name, io = 'auto', out = out )

            ok = ok and ret

        return ok
