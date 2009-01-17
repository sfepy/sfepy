from sfepy.base.base import *
from sfepy.fem import eval_term_op
from utils import iter_sym, create_pis, create_scalar_pis

class MiniAppBase( Struct ):
    def any_from_conf( name, problem, kwargs ):
        cls = kwargs.get( 'class' )
        obj = cls( name, problem, kwargs )
        return obj
    any_from_conf = staticmethod( any_from_conf )

    def __init__( self, name, problem, kwargs ):
        Struct.__init__( self, name = name, problem = problem, **kwargs )
        self.set_default_attr( 'requires', [] )

    def make_save_hook( self, base_name, format,
                        post_process_hook = None, file_per_var = None ):
        return None

class ShapeDimDim( MiniAppBase ):
    
    def __call__( self, problem = None, data = None, save_hook = None ):
        problem = get_default( problem, self.problem )

        return create_pis( problem, self.variables[0] )

class ShapeDim( MiniAppBase ):
    
    def __call__( self, problem = None, data = None, save_hook = None ):
        problem = get_default( problem, self.problem )

        return create_scalar_pis( problem, self.variables[0] )

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

    def make_save_hook( self, base_name, format,
                        post_process_hook = None, file_per_var = None ):
        def save_correctors( state, problem, ir, ic ):
            problem.save_state( (base_name + '_%d%d.' % (ir, ic)) + format,
                                state,
                                post_process_hook = post_process_hook,
                                file_per_var = file_per_var )
        return save_correctors

class CorrDim( MiniAppBase ):
    def get_variables( self, ir, data ):
            raise StopIteration

    def __call__( self, problem = None, data = None, save_hook = None ):
        problem = get_default( problem, self.problem )

        problem.select_variables( self.variables )
        problem.set_equations( self.equations )

        problem.select_bcs( ebc_names = self.ebcs, epbc_names = self.epbcs )

        dim = problem.domain.mesh.dim
        states = nm.zeros( (dim,), dtype = nm.object )
        for ir in range( dim ):
            for name, val in self.get_variables( ir, data ):
                problem.variables[name].data_from_data( val )

            state = problem.create_state_vector()
            problem.apply_ebc( state )
            state = problem.solve()
            assert_( problem.variables.has_ebc( state ) )
            states[ir] = state

            if save_hook is not None:
                save_hook( state, problem, ir )

        return Struct( name = self.name,
                       states = states,
                       di = problem.variables.di )

    def make_save_hook( self, base_name, format,
                        post_process_hook = None, file_per_var = None ):
        def save_correctors( state, problem, ir ):
            problem.save_state( (base_name + '_%d.' % (ir,)) + format,
                                state,
                                post_process_hook = post_process_hook,
                                file_per_var = file_per_var )
        return save_correctors

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

class CoefDimSym( MiniAppBase ):

    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )

        dim, sym = problem.get_dim( get_sym = True )
        coef = nm.zeros( (dim, sym), dtype = nm.float64 )

        for ir in range( dim ):
            for name, val in self.get_variables( problem, ir, None, data,
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
    
class CoefDimDim( MiniAppBase ):
    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )

        dim = problem.get_dim()
        coef = nm.zeros( (dim, dim), dtype = nm.float64 )

        for ir in range( dim ):
            for name, val in self.get_variables( problem, ir, None, data,
                                                 'row' ):
                problem.variables[name].data_from_data( val )

            for ic in range( dim ):
                for name, val in self.get_variables( problem, None, ic, data,
                                                     'col' ):
                    problem.variables[name].data_from_data( val )

                val = eval_term_op( None, self.expression,
                                    problem, call_mode = 'd_eval' )

                coef[ir,ic] = val

        coef /= volume

        return coef
