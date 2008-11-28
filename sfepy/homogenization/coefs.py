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
    """Steady state correctors $\bar{\omega}^{rs}$."""

    def get_variables( self, ir, ic, data ):
        """data: pis"""
        yield (self.variables[2], data[ir,ic])

class CoefSymSym( MiniAppBase ):
    
    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )

        dim, sym = problem.get_dim( get_sym = True )
        coef = nm.zeros( (sym, sym), dtype = nm.float64 )

        for ir, (irr, icr) in enumerate( iter_sym( dim ) ):
            for name, val in self.get_variables( problem, irr, icr, data,
                                                 'row' ):
                problem.variables[name].data_from_data( val )

            for ic, (irc, icc) in enumerate( iter_sym( dim ) ):
                for name, val in self.get_variables( problem, irc, icc, data,
                                                     'col' ):
                    problem.variables[name].data_from_data( val )

                val = eval_term_op( None, self.expression,
                                    problem, call_mode = 'd_eval' )

                coef[ir,ic] = val

        coef /= volume

        return coef

class ElasticCoef( CoefSymSym ):
    """Homogenized elastic tensor $E_{ijkl}$."""

    mode2var = {'row' : 0, 'col' : 1}

    def get_variables( self, problem, ir, ic, data, mode ):
        var_name = self.variables[self.mode2var[mode]]
        u_name = problem.variables[var_name].primary_var_name

        corrs, pis = data['corrs'], data['pis']
        indx = corrs.di.indx[u_name]

        omega = corrs.states[ir,ic][indx]
        pi = pis[ir,ic] + omega

        yield (var_name, pi)
