from sfepy.base.base import *
from sfepy.fem import eval_term_op, assemble_by_blocks
from sfepy.solvers import Solver
from sfepy.base.la import MatrixAction, eig
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

        self.init_solvers(problem)

        dim = problem.domain.mesh.dim
        states = nm.zeros( (dim, dim), dtype = nm.object )
        for ir in range( dim ):
            for ic in range( dim ):
                for name, val in self.get_variables( ir, ic, data ):
                    problem.variables[name].data_from_data( val )

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

        self.init_solvers(problem)

        dim = problem.domain.mesh.dim
        states = nm.zeros( (dim,), dtype = nm.object )
        for ir in range( dim ):
            for name, val in self.get_variables( ir, data ):
                problem.variables[name].data_from_data( val )

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

class CorrOne( MiniAppBase ):
    def get_variables( self, data ):
            raise StopIteration

    def __call__( self, problem = None, data = None, save_hook = None ):
        problem = get_default( problem, self.problem )

        problem.select_variables( self.variables )
        problem.set_equations( self.equations )

        problem.select_bcs( ebc_names = self.ebcs, epbc_names = self.epbcs )

        for name, val in self.get_variables( data ):
            problem.variables[name].data_from_data( val )

        state = problem.solve()
        assert_( problem.variables.has_ebc( state ) )

        if save_hook is not None:
            save_hook( state, problem )

        return Struct( name = self.name,
                       state = state,
                       di = problem.variables.di )

    def make_save_hook( self, base_name, format,
                        post_process_hook = None, file_per_var = None ):
        def save_correctors( state, problem ):
            problem.save_state( '.'.join( (base_name, format) ),
                                state,
                                post_process_hook = post_process_hook,
                                file_per_var = file_per_var )
        return save_correctors

class PressureEigenvalueProblem( MiniAppBase ):
    """Pressure eigenvalue problem solver for time-dependent correctors."""

    def presolve( self, mtx ):
        """Prepare A^{-1} B^T for the Schur complement."""

        mtx_a = mtx['A']
        mtx_bt = mtx['BT']
        output( 'full A size: %.3f MB' % (8.0 * nm.prod( mtx_a.shape ) / 1e6) )
        output( 'full B size: %.3f MB' % (8.0 * nm.prod( mtx_bt.shape ) / 1e6) )

        ls = Solver.any_from_conf( self.problem.ls_conf,
                                   presolve = True, mtx = mtx_a )
        if self.mode == 'explicit':
            tt = time.clock()
            mtx_aibt = nm.zeros( mtx_bt.shape, dtype = mtx_bt.dtype )
            for ic in xrange( mtx_bt.shape[1] ):
                mtx_aibt[:,ic] = ls( mtx_bt[:,ic].toarray().squeeze() )
            output( 'mtx_aibt: %.2f s' % (time.clock() - tt) )
            action_aibt = MatrixAction.from_array( mtx_aibt )
        else:
            ##
            # c: 30.08.2007, r: 13.02.2008
            def fun_aibt( vec ):
                # Fix me for sparse mtx_bt...
                rhs = sc.dot( mtx_bt, vec )
                out = ls( rhs )
                return out
            action_aibt = MatrixAction.from_function( fun_aibt,
                                                    (mtx_a.shape[0],
                                                     mtx_bt.shape[1]),
                                                    nm.float64 )
        mtx['action_aibt'] = action_aibt

    def solve_pressure_eigenproblem( self, mtx, eig_problem = None,
                                     n_eigs = 0, check = False ):
        """G = B*AI*BT or B*AI*BT+D"""

        def get_slice( n_eigs, nn ):
            if n_eigs > 0:
                ii = slice( 0, n_eigs )
            elif n_eigs < 0:
                ii = slice( nn + n_eigs, nn )
            else:
                ii = slice( 0, 0 )
            return ii

        eig_problem = get_default( eig_problem, self.eig_problem )
        n_eigs = get_default( n_eigs, self.n_eigs )
        check = get_default( check, self.check )

        mtx_c, mtx_b, action_aibt = mtx['C'], mtx['B'], mtx['action_aibt']
        mtx_g = mtx_b * action_aibt.to_array() # mtx_b must be sparse!
        if eig_problem == 'B*AI*BT+D':
            mtx_g += mtx['D'].toarray()

        mtx['G'] = mtx_g
        output( mtx_c.shape, mtx_g.shape )

        eigs, mtx_q = eig( mtx_c.toarray(), mtx_g, method = 'eig.sgscipy' )

        if check:
            ee = nm.diag( sc.dot( mtx_q.T * mtx_c, mtx_q ) ).squeeze()
            oo = nm.diag( sc.dot( sc.dot( mtx_q.T,  mtx_g ), mtx_q ) ).squeeze()
            try:
                assert_( nm.allclose( ee, eigs ) )
                assert_( nm.allclose( oo, nm.ones_like( eigs ) ) )
            except ValueError:
                debug()

        nn = mtx_c.shape[0]
        if isinstance( n_eigs, tuple ):
            output( 'required number of eigenvalues: (%d, %d)' % n_eigs )
            if sum( n_eigs ) < nn:
                ii0 = get_slice( n_eigs[0], nn )
                ii1 = get_slice( -n_eigs[1], nn )
                eigs = nm.concatenate( (eigs[ii0], eigs[ii1] ) )
                mtx_q = nm.concatenate( (mtx_q[:,ii0], mtx_q[:,ii1]), 1 ) 
        else:
            output( 'required number of eigenvalues: %d' % n_eigs )
            if (n_eigs != 0) and (abs( n_eigs ) < nn):
                ii = get_slice( n_eigs, nn )
                eigs = eigs[ii]
                mtx_q = mtx_q[:,ii]

##         from sfepy.base.plotutils import pylab, iplot
##         pylab.semilogy( eigs )
##         pylab.figure( 2 )
##         iplot( eigs )
##         pylab.show()
##         debug()

        out = Struct( eigs = eigs, mtx_q = mtx_q )
        return out

    def __call__( self, problem = None, data = None, save_hook = None ):
        problem = get_default( problem, self.problem )

        mtx = assemble_by_blocks( self.equations, problem,
                                  ebcs = self.ebcs, epbcs = self.epbcs,
                                  restore_variables = False )
        self.presolve( mtx )

        evp = self.solve_pressure_eigenproblem( mtx )
        return Struct( name = self.name,
                       ebcs = self.ebcs, epbcs = self.epbcs,
                       mtx = mtx, evp = evp )

    def make_save_hook( self, base_name, format,
                        post_process_hook = None, file_per_var = None ):
        return None

class TSTimes( MiniAppBase ):
    """Coefficient-like class, returns times of the time stepper."""
    def __call__( self, volume = None, problem = None, data = None ):
        problem = get_default( problem, self.problem )
        return problem.get_time_solver().ts.times

class VolumeFractions( MiniAppBase ):
    """Coefficient-like class, returns volume fractions of given regions within
    the whole domain."""
    def __call__( self, volume = None, problem = None, data = None ):
        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )

        vf = {}
        for region_name in self.regions:
            vkey = 'volume_%s' % region_name
            key = 'fraction_%s' % region_name

            val = eval_term_op( None, self.expression % region_name, problem )
            vf[vkey] = nm.asarray( val, dtype = nm.float64 )
            vf[key] = vf[vkey] / volume

        return vf

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

class CoefSym( MiniAppBase ):
    
    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )

        dim, sym = problem.get_dim( get_sym = True )
        coef = nm.zeros( (sym,), dtype = nm.float64 )

        for name, val in self.get_variables( problem, None, None, data, 'col' ):
            problem.variables[name].data_from_data( val )

        if isinstance( self.expression, tuple ):
            coef[:] = eval_term_op( None, self.expression[0],
                                    problem, shape = (sym,), mode = 'const' )
            expression = self.expression[1]
        else:
            expression = self.expression

        for ii, (ir, ic) in enumerate( iter_sym( dim ) ):
            for name, val in self.get_variables( problem, ir, ic, data, 'row' ):
                problem.variables[name].data_from_data( val )

            val = eval_term_op( None, expression,
                                problem, call_mode = 'd_eval' )
            coef[ii] += val

        coef /= volume

        return coef

class CoefOne( MiniAppBase ):
    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )

        coef = nm.zeros( (1,), dtype = nm.float64 )

        for name, val in self.get_variables( problem, data ):
            problem.variables[name].data_from_data( val )

        if isinstance( self.expression, tuple ):
            coef[:] = eval_term_op( None, self.expression[0],
                                    problem, shape = (1,), mode = 'const' )
            expression = self.expression[1]
        else:
            expression = self.expression

        val = eval_term_op( None, expression,
                            problem, call_mode = 'd_eval' )
        coef += val

        coef /= volume

        return coef
