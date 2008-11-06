from sfepy.base.base import *
from sfepy.fem import eval_term_op
from utils import iter_sym

class MiniAppBase( Struct ):
    def __init__( self, problem, kwargs ):
        Struct.__init__( self, problem = problem, **kwargs )

class CorrectorsRS( MiniAppBase ):

    def __call__( self, pis, problem = None, save_hook = None ):
        problem = get_default( problem, self.problem )

        var_names = self.variables
        pi_name = var_names[2]
        problem.select_variables( var_names )
        problem.set_equations( self.equations )

        problem.select_bcs( ebc_names = self.ebcs, epbc_names = self.epbcs )

        dim = problem.domain.mesh.dim
        states_rs = nm.zeros( (dim, dim), dtype = nm.object )
        for ir in range( dim ):
            for ic in range( dim ):
                pi = pis[ir,ic]
                problem.variables[pi_name].data_from_data( pi )

                state = problem.create_state_vector()
                problem.apply_ebc( state )
                state = problem.solve()
                assert_( problem.variables.has_ebc( state ) )
                states_rs[ir,ic] = state

                if save_hook is not None:
                    save_hook( state, problem, ir, ic )

        return Struct( name = 'Steady RS correctors',
                       states_rs = states_rs,
                       di = problem.variables.di )

class CoefE( MiniAppBase ):

    def __call__( self, pis, corrs_rs, volume, problem = None ):
        problem = get_default( problem, self.problem )
        var_names = self.variables

        problem.select_variables( var_names )

        dim = problem.domain.mesh.dim
        sym = (dim + 1) * dim / 2
        coef = nm.zeros( (sym, sym), dtype = nm.float64 )

        u_name = problem.variables[var_names[0]].primary_var_name
        indx = corrs_rs.di.indx[u_name]

        for ir, (irr, icr) in enumerate( iter_sym( dim ) ):
            omega1 = corrs_rs.states_rs[irr,icr][indx]
            pi1 = pis[irr,icr] + omega1
            problem.variables[var_names[0]].data_from_data( pi1 )

            for ic, (irc, icc) in enumerate( iter_sym( dim ) ):
                omega2 = corrs_rs.states_rs[irc,icc][indx]
                pi2 = pis[irc,icc] + omega2
                problem.variables[var_names[1]].data_from_data( pi2 )

                val = eval_term_op( None, self.expression,
                                    problem, call_mode = 'd_eval' )

                coef[ir,ic] = val

        coef /= volume

        return coef
