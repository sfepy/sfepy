from sfepy.base.base import *
from sfepy.fem import eval_term_op
from utils import iter_sym, create_pis, create_scalar_pis

class MiniAppBase( Struct ):
    def any_from_conf( name, problem, kwargs ):
        try:
            cls = kwargs['class']
        except KeyError:
            raise KeyError("set 'class' for MiniApp %s!" % name)
        obj = cls( name, problem, kwargs )
        return obj
    any_from_conf = staticmethod( any_from_conf )

    def __init__( self, name, problem, kwargs ):
        Struct.__init__( self, name = name, problem = problem, **kwargs )
        self.set_default_attr( 'requires', [] )
        self.set_default_attr( 'is_linear', False )

    def init_solvers(self, problem):
        """For linear problems, assemble the matrix and try to presolve the
        linear system."""
        if self.is_linear:
            output('linear problem, trying to presolve...')
            tt = time.clock()

            ev = problem.get_evaluator( mtx = problem.mtx_a )

            state = problem.create_state_vector()
            try:
                mtx_a = ev.eval_tangent_matrix( state, is_full = True )
            except ValueError:
                raise ValueError('matrix evaluation failed, giving up...')

            problem.set_linear(True)
            problem.init_solvers(mtx=mtx_a, presolve=True)

            output( '...done in %.2f s' % (time.clock() - tt) )

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

class CorrNN( MiniAppBase ):
    """ __init__() kwargs:
        {
             'variables' : [],
             'ebcs' : [],
             'epbcs' : [],
             'equations' : {},
        },
    """

    def __init__( self, Ndim, name, problem, kwargs ):
        MiniAppBase.__init__( self, name, problem, kwargs )
        self.set_default_attr( 'Ndim', Ndim )

    def get_variables( self, ir, ic, data ):
            raise StopIteration

    def __call__( self, problem = None, data = None, save_hook = None ):
        problem = get_default( problem, self.problem )

        problem.select_variables( self.variables )
        problem.set_equations( self.equations )

        problem.select_bcs( ebc_names = self.ebcs, epbc_names = self.epbcs )

        self.init_solvers(problem)

        states = nm.zeros( (self.Ndim, self.Ndim), dtype = nm.object )
        for ir in range( self.Ndim ):
            for ic in range( self.Ndim ):
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

class CorrN( MiniAppBase ):
    def get_variables( self, ir, data ):
            raise StopIteration

    def __init__( self, Ndim, name, problem, kwargs ):
        MiniAppBase.__init__( self, name, problem, kwargs )
        self.set_default_attr( 'Ndim', Ndim )

    def __call__( self, problem = None, data = None, save_hook = None ):
        problem = get_default( problem, self.problem )

        problem.select_variables( self.variables )
        problem.set_equations( self.equations )

        problem.select_bcs( ebc_names = self.ebcs, epbc_names = self.epbcs )

        self.init_solvers(problem)

        states = nm.zeros( (self.Ndim,), dtype = nm.object )
        for ir in range( self.Ndim ):
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

class CorrDimDim( CorrNN ):

    def __init__( self, name, problem, kwargs ):
        CorrNN.__init__( self, problem.domain.mesh.dim,
                         name, problem, kwargs )

class CorrDim( CorrN ):

    def __init__( self, name, problem, kwargs ):
        CorrN.__init__( self, problem.domain.mesh.dim,
                        name, problem, kwargs )

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
    
class CoefNN( MiniAppBase ):

    def __init__( self, Ndim, name, problem, kwargs ):
        MiniAppBase.__init__( self, name, problem, kwargs )
        self.set_default_attr( 'Ndim', Ndim )

    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )

        coef = nm.zeros( (self.Ndim, self.Ndim), dtype = nm.float64 )

        for ir in range( self.Ndim ):
            for name, val in self.get_variables( problem, ir, None, data,
                                                 'row' ):
                problem.variables[name].data_from_data( val )

            for ic in range( self.Ndim ):
                for name, val in self.get_variables( problem, None, ic, data,
                                                     'col' ):
                    problem.variables[name].data_from_data( val )

                val = eval_term_op( None, self.expression,
                                    problem, call_mode = 'd_eval' )

                coef[ir,ic] = val

        coef /= volume

        return coef

class CoefN( MiniAppBase ):

    def __init__( self, Ndim, name, problem, kwargs ):
        MiniAppBase.__init__( self, name, problem, kwargs )
        self.set_default_attr( 'Ndim', Ndim )

    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )

        coef = nm.zeros( (self.Ndim, ), dtype = nm.float64 )

        for ir in range( self.Ndim ):
            for name, val in self.get_variables( problem, ir, data ):
                problem.variables[name].data_from_data( val )

            val = eval_term_op( None, self.expression,
                                problem, call_mode = 'd_eval' )
            coef[ir] = val

        coef /= volume

        return coef

class CoefDimDim( CoefNN ):

    def __init__( self, name, problem, kwargs ):
        CoefNN.__init__( self, problem.domain.mesh.dim,
                         name, problem, kwargs )

class CoefDim( CoefN ):

    def __init__( self, name, problem, kwargs ):
        CoefN.__init__( self, problem.domain.mesh.dim,
                        name, problem, kwargs )
