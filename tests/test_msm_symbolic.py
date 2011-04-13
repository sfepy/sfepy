# c: 07.05.2007, r: 08.07.2008
from sfepy import data_dir

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
    """2 * dw_laplace.i1.Omega( coef.val, s, t )
    """,
    'Diffusion' :
    """3 * dw_diffusion.i1.Omega( coef.K, s, t )
    """,
}
equations_rhs = {
    'Laplace' :
    """= - dw_volume_lvf.i1.Omega( rhs.val, s )""",
    'Diffusion' :
    """= - dw_volume_lvf.i1.Omega( rhs.val, s )""",
}

solutions = {
    'sincos' : ('t', 'sin( 3.0 * x ) * cos( 4.0 * y )'),
    'poly' : ('t', '(x**2) + (y**2)'),
    'polysin' : ('t', '((x - 0.5)**3) * sin( 5.0 * y )'),
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
try:
    import sympy_operators as sops
except ImportError, exc:
    sops = None
    
from sfepy.base.testing import TestCommon
from sfepy.base.base import debug, pause
output_name = 'test_msm_symbolic_%s.vtk'

##
# c: 07.05.2007, r: 09.05.2008
solution = ['']
def ebc(ts, coor, solution=None):
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
    'ebc' : (lambda ts, coor, **kwargs:
             ebc(ts, coor, solution=solution),),
    'rhs' : (rhs,),
}

##
# c: 07.05.2008
class Test( TestCommon ):

    ##
    # c: 07.05.2007, r: 25.06.2008
    def from_conf( conf, options ):
        from sfepy.fem import ProblemDefinition

        problem = ProblemDefinition.from_conf( conf, init_equations = False )
        test = Test( problem = problem,
                     conf = conf, options = options )
        return test
    from_conf = staticmethod( from_conf )

    ##
    # c: 09.05.2007, r: 08.07.2008
    def _build_rhs( self, equation, sols ):

        problem  = self.problem
        rhss = {}
        self.report( '%s:' % equation.name )
        self.report( 'evaluating terms, "<=" is solution, "=>" is the rhs:' )
        for term in equation.terms:
            if not hasattr( term, 'symbolic' ):
                self.report( 'term %s has no symbolic description!' % term.name )
                raise ValueError
            expr = term.symbolic['expression']
            arg_map = term.symbolic['map']
            self.report( '%s( %s )' %\
                         (term.name, ', '.join( term.ats )) )
            self.report( 'multiplicator: %f' % term.sign )
            self.report( '  symbolic:', expr )
            self.report( '  using argument map:', arg_map )

            for sol_name, sol in sols.iteritems():
                rhs = self._eval_term( sol[1], term, sops )
                srhs = "(%s * (%s))" % (term.sign, rhs)
                rhss.setdefault( sol_name, [] ).append( srhs )

        for key, val in rhss.iteritems():
            rhss[key] = '+'.join( val )

        return rhss

    ##
    # c: 09.05.2007, r: 25.06.2008
    def _eval_term( self, sol, term, sops ):
        """Works for scalar, single unknown terms only!"""
        expr = term.symbolic['expression']
        arg_map = term.symbolic['map']
        env = {'x' : sops.Symbol( 'x' ),
               'y' : sops.Symbol( 'y' ),
               'z' : sops.Symbol( 'z' ),
               'dim' : dim}
        for key, val in arg_map.iteritems():
            if val == 'state':
                env[key] = sol
            else:
                term.set_current_group(0)
                env[key] = term.get_args( [val] )[0]

            if val[:8] == 'material':
                # Take the first value - constant in all QPs.
                aux = env[key][0,0]
                if nm.prod( aux.shape ) == 1:
                    env[key] = aux.squeeze()
                else:
                    import sympy
                    env[key] = sympy.Matrix( aux )

#        print env

        self.report( '  <= ', sol )
        sops.set_dim( dim )
        val = str( eval( expr, sops.__dict__, env ) )
        self.report( '   =>', val )
        return val

    ##
    # c: 07.05.2007, r: 30.06.2008
    def _test_msm_symbolic( self, equations ):
        import os.path as op

        if sops is None:
            self.report( 'cannot import sympy, skipping' )
            return True

        problem  = self.problem
        ok = True
        for eq_name, equation in equations.iteritems():
            problem.set_equations( {eq_name : equation} )
            problem.update_materials()

            rhss = self._build_rhs( problem.equations[eq_name],
                                   self.conf.solutions )
            erhs = problem.conf.equations_rhs[eq_name]  

            problem.set_equations( {eq_name : equation + erhs} )
            variables = problem.get_variables()
            materials = problem.get_materials()
            rhs_mat = materials['rhs']

            for sol_name, sol in problem.conf.solutions.iteritems():
                self.report( 'testing', sol_name )
                var_name, sol_expr = sol
                rhs_expr = rhss[sol_name]

                self.report( 'sol:', sol_expr )
                self.report( 'rhs:', rhs_expr )
                globals()['solution'][0] = sol_expr
                rhs_mat.function.set_extra_args(expression=rhs_expr)
                problem.time_update()
                problem.equations.reset_term_caches()
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

                problem.domain.mesh.write( fname % '_'.join( (sol_name, eq_name) ),
                                           io = 'auto', out = out )

                ok = ok and ret

        return ok

    ##
    # c: 30.06.2008, r: 30.06.2008
    def _get_equations( self, name ):
        """Choose a sub-problem from all equations."""
        return {name : self.problem.conf.equations[name]}
        
    ##
    # c: 30.06.2008, r: 30.06.2008
    def test_msm_symbolic_laplace( self ):
        return self._test_msm_symbolic( self._get_equations( 'Laplace' ) )

    ##
    # c: 30.06.2008, r: 30.06.2008
    def test_msm_symbolic_diffusion( self ):
        return self._test_msm_symbolic( self._get_equations( 'Diffusion' ) )
