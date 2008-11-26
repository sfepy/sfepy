from sfepy.base.base import *
from sfepy.fem import eval_term_op
from utils import iter_sym

class MiniAppBase( Struct ):
    def __init__( self, name, problem, kwargs ):
        Struct.__init__( self, name = name, problem = problem, **kwargs )

class CorrDimDim( MiniAppBase ):
    """ __init__() kwargs:
        {
             'variables' : [],
             'ebcs' : [],
             'epbcs' : [],
             'equations' : {},
        },
    """

    def get_variables( self, ir, ic, data ):
            raise StopIteration

    def __call__( self, problem = None, data = None, save_hook = None ):
        problem = get_default( problem, self.problem )

        problem.select_variables( self.variables )
        problem.set_equations( self.equations )

        problem.select_bcs( ebc_names = self.ebcs, epbc_names = self.epbcs )

        dim = problem.domain.mesh.dim
        states = nm.zeros( (dim, dim), dtype = nm.object )
        for ir in range( dim ):
            for ic in range( dim ):
                for name, val in self.get_variables( ir, ic, data ):
                    problem.variables[name].data_from_data( val )

                state = problem.create_state_vector()
                problem.apply_ebc( state )
                state = problem.solve()
                assert_( problem.variables.has_ebc( state ) )
                states[ir,ic] = state

                if save_hook is not None:
                    save_hook( state, problem, ir, ic )

        return Struct( name = self.name,
                       states = states,
                       di = problem.variables.di )

class CorrectorsRS( CorrDimDim ):
    def get_variables( self, ir, ic, data ):
        """data: pis"""
        yield (self.variables[2], data[ir,ic])

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
            omega1 = corrs_rs.states[irr,icr][indx]
            pi1 = pis[irr,icr] + omega1
            problem.variables[var_names[0]].data_from_data( pi1 )

            for ic, (irc, icc) in enumerate( iter_sym( dim ) ):
                omega2 = corrs_rs.states[irc,icc][indx]
                pi2 = pis[irc,icc] + omega2
                problem.variables[var_names[1]].data_from_data( pi2 )

                val = eval_term_op( None, self.expression,
                                    problem, call_mode = 'd_eval' )

                coef[ir,ic] = val

        coef /= volume

        return coef
