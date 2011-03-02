import os
import time

import numpy as nm
import numpy.linalg as nla
import scipy as sc

from sfepy.base.base import output, assert_, get_default, debug, Struct
from sfepy.fem import assemble_by_blocks
from sfepy.fem.evaluate import eval_equations
from sfepy.solvers.ts import TimeStepper
from sfepy.fem.meshio import HDF5MeshIO
from sfepy.solvers import Solver, eig
from sfepy.linalg import MatrixAction
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

        self.problem.clear_equations()
        self.set_default_attr( 'requires', [] )
        self.set_default_attr( 'is_linear', False )

    def init_solvers(self, problem):
        """For linear problems, assemble the matrix and try to presolve the
        linear system."""
        if self.is_linear:
            output('linear problem, trying to presolve...')
            tt = time.clock()

            ev = problem.get_evaluator()

            state = problem.create_state()
            try:
                mtx_a = ev.eval_tangent_matrix(state(), is_full=True)
            except ValueError:
                raise ValueError('matrix evaluation failed, giving up...')

            problem.set_linear(True)
            problem.init_solvers(mtx=mtx_a, presolve=True)

            output( '...done in %.2f s' % (time.clock() - tt) )

class CorrSolution(Struct):
    """
    Class for holding solutions of corrector problems.
    """

    def iter_solutions(self):
        if hasattr(self, 'components'):
            for indx in self.components:
                key = ('%d' * len(indx)) % indx
                yield key, self.states[indx]

        else:
            yield '', self.state

class CorrMiniApp( MiniAppBase ):

    def __init__( self, name, problem, kwargs ):
        MiniAppBase.__init__( self, name, problem, kwargs )
        self.output_dir = self.problem.output_dir
        self.set_default_attr( 'save_name', '(not_set)' )
        self.set_default_attr( 'dump_name', self.save_name )
        self.set_default_attr( 'dump_variables', [] )
        self.set_default_attr( 'save_variables', self.dump_variables )

        self.save_name = os.path.normpath( os.path.join( self.output_dir,
                                                         self.save_name ) )
        self.dump_name = os.path.normpath( os.path.join( self.output_dir,
                                                         self.dump_name ) )
    def setup_output( self, save_format = None, dump_format = None,
                      post_process_hook = None, file_per_var = None ):
        """Instance attributes have precedence!"""
        self.set_default_attr( 'dump_format', dump_format )
        self.set_default_attr( 'save_format', save_format )
        self.set_default_attr( 'post_process_hook', post_process_hook )
        self.set_default_attr( 'file_per_var', file_per_var )

    def get_save_name_base( self ):
        return self.save_name

    def get_dump_name_base( self ):
        return self.get_save_name_base()

    def get_save_name( self ):
        return '.'.join( (self.get_save_name_base(), self.save_format))

    def get_dump_name( self ):
        if self.dump_format is not None:
            return '.'.join( (self.get_dump_name_base(), self.dump_format))
        else:
            return None

    def get_output(self, corr_sol, is_dump=False, extend=True):
        variables = self.problem.get_variables()
        to_output = variables.state_to_output

        if is_dump:
            var_names = self.dump_variables
            extend = False

        else:
            var_names = self.save_variables

        out = {}
        for key, sol in corr_sol.iter_solutions():
            for var_name in var_names:
                skey = var_name + '_' + key

                dof_vector = sol[var_name]

                if is_dump:
                        var = variables[var_name]
                        shape = (var.n_dof / var.n_components,
                                 var.n_components)
                        out[skey] = Struct(name = 'dump', mode = 'nodes',
                                           data = dof_vector,
                                           dofs = var.dofs,
                                           shape = shape,
                                           var_name = var_name)

                else:
                    aux = to_output(dof_vector,
                                    var_info={var_name: (True, var_name)},
                                    extend=extend)
                    if self.post_process_hook is not None:
                        aux = self.post_process_hook(aux, self.problem,
                                                     None,
                                                     extend=extend)

                    for _key, val in aux.iteritems():
                        new_key = _key + '_' + key
                        out[new_key] = val

        return out

    def save(self, state, problem):
        save_name = self.get_save_name()
        if save_name is not None:
            extend = not self.file_per_var
            out = self.get_output(state, extend=extend)

            problem.save_state(save_name, out=out,
                               file_per_var=self.file_per_var)

        dump_name = self.get_dump_name()
        if dump_name is not None:
            problem.save_state(dump_name,
                               out=self.get_output(state, is_dump=True),
                               file_per_var=False)

class ShapeDimDim( CorrMiniApp ):

    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        clist, pis = create_pis( problem, self.variables[0] )

        corr_sol = CorrSolution(name = self.name,
                                states = pis,
                                components = clist)
        return corr_sol

class ShapeDim( CorrMiniApp ):

    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        clist, pis = create_scalar_pis( problem, self.variables[0] )

        corr_sol = CorrSolution(name = self.name,
                                states = pis,
                                components = clist)
        return corr_sol

class CorrNN( CorrMiniApp ):
    """ __init__() kwargs:
        {
             'ebcs' : [],
             'epbcs' : [],
             'equations' : {},
             'set_variables' : None,
        },
    """

    def set_variables_default(variables, ir, ic, set_var, data):
        for (var, req, comp) in set_var:
            variables[var].data_from_any(data[req].states[ir,ic][comp])

    set_variables_default = staticmethod(set_variables_default)

    def __init__( self, name, problem, kwargs ):
        """When dim is not in kwargs, problem dimension is used."""
        CorrMiniApp.__init__( self, name, problem, kwargs )
        self.set_default_attr( 'dim', problem.get_dim() )

    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        problem.set_equations( self.equations )

        problem.select_bcs(ebc_names = self.ebcs, epbc_names = self.epbcs,
                           lcbc_names = self.get_default_attr('lcbcs', []))

        problem.update_materials(problem.ts)

        self.init_solvers(problem)

        variables = problem.get_variables()

        states = nm.zeros( (self.dim, self.dim), dtype = nm.object )
        clist = []
        for ir in range( self.dim ):
            for ic in range( self.dim ):
                if isinstance(self.set_variables, list):
                    self.set_variables_default(variables, ir, ic,
                                               self.set_variables, data)
                else:
                    self.set_variables(variables, ir, ic, **data)

                state = problem.solve()
                assert_(state.has_ebc())
                states[ir,ic] = state.get_parts()

                clist.append( (ir, ic) )

        corr_sol = CorrSolution(name = self.name,
                                states = states,
                                components = clist)

        self.save(corr_sol, problem)

        return corr_sol

class CorrN( CorrMiniApp ):

    def set_variables_default(variables, ir, set_var, data):
        for (var, req, comp) in set_var:
            variables[var].data_from_any(data[req].states[ir][comp])

    set_variables_default = staticmethod(set_variables_default)

    def __init__( self, name, problem, kwargs ):
        """When dim is not in kwargs, problem dimension is used."""
        CorrMiniApp.__init__( self, name, problem, kwargs )
        self.set_default_attr( 'dim', problem.get_dim() )

    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        problem.set_equations( self.equations )

        problem.select_bcs(ebc_names = self.ebcs, epbc_names = self.epbcs,
                           lcbc_names = self.get_default_attr('lcbcs', []))

        problem.update_materials(problem.ts)

        self.init_solvers(problem)

        variables = problem.get_variables()

        states = nm.zeros( (self.dim,), dtype = nm.object )
        clist = []
        for ir in range( self.dim ):
            if isinstance(self.set_variables, list):
                self.set_variables_default(variables, ir,
                                           self.set_variables, data)
            else:
                self.set_variables(variables, ir, **data)
            state = problem.solve()
            assert_(state.has_ebc())
            states[ir] = state.get_parts()

            clist.append((ir,))

        corr_sol = CorrSolution(name = self.name,
                                states = states,
                                components = clist)

        self.save(corr_sol, problem)

        return corr_sol

class CorrDimDim( CorrNN ):
    pass

class CorrDim( CorrN ):
    pass

class CorrOne( CorrMiniApp ):

    def set_variables_default(variables, set_var, data):
        for (var, req, comp) in set_var:
            variables[var].data_from_any(data[req].state[comp])

    set_variables_default = staticmethod(set_variables_default)

    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        problem.set_equations( self.equations )

        problem.select_bcs(ebc_names = self.ebcs, epbc_names = self.epbcs,
                           lcbc_names = self.get_default_attr('lcbcs', []))

        problem.update_materials(problem.ts)

        self.init_solvers(problem)

        variables = problem.get_variables()

        if hasattr(self, 'set_variables'):
            if isinstance(self.set_variables, list):
                self.set_variables_default(variables, self.set_variables,
                                           data)
            else:
                self.set_variables(variables, **data)

        state = problem.solve()
        assert_(state.has_ebc())

        corr_sol = CorrSolution(name = self.name,
                                state = state.get_parts())

        self.save(corr_sol, problem)

        return corr_sol

class CorrEqPar(CorrOne):
    """
    The corrector which equation can be parametrized via 'eq_pars',
    the dimension is given by the number of parameters.

    Example:

        'equations': 'dw_diffusion.5.Y(mat.k, q, p) =
                      dw_surface_integrate.5.%s(q)',
        'eq_pars': ('bYMp', 'bYMm'),
        'class': cb.CorrEqPar,

    """

    def __init__( self, name, problem, kwargs ):
        """When dim is not in kwargs, problem dimension is used."""
        CorrMiniApp.__init__(self, name, problem, kwargs)
        self.set_default_attr('dim', len(self.eq_pars))

    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        states = nm.zeros( (self.dim,), dtype = nm.object )
        clist = []

        eqns ={}
        for ir in range( self.dim ):
            for key_eq, val_eq in self.equations.iteritems():
                eqns[key_eq] = val_eq % self.eq_pars[ir]

            problem.set_equations(eqns)

            problem.select_bcs(ebc_names = self.ebcs, epbc_names = self.epbcs,
                               lcbc_names = self.get_default_attr('lcbcs', []))

            problem.update_materials(problem.ts)

            self.init_solvers(problem)

            variables = problem.get_variables()

            if hasattr(self, 'set_variables'):
                if isinstance(self.set_variables, list):
                    self.set_variables_default(variables, self.set_variables,
                                               data)
                else:
                    self.set_variables(variables, **data)

            state = problem.solve()
            assert_(state.has_ebc())

            states[ir] = state.get_parts()
            clist.append((ir,))

        corr_sol = CorrSolution(name = self.name,
                                states = states,
                                components = clist)

        self.save(corr_sol, problem)

        return corr_sol

class PressureEigenvalueProblem( CorrMiniApp ):
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

    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        mtx = assemble_by_blocks( self.equations, problem,
                                  ebcs = self.ebcs, epbcs = self.epbcs )
        self.presolve( mtx )

        evp = self.solve_pressure_eigenproblem( mtx )
        return Struct( name = self.name,
                       ebcs = self.ebcs, epbcs = self.epbcs,
                       mtx = mtx, evp = evp )

class TCorrectorsViaPressureEVP( CorrMiniApp ):
    """
    Time correctors via the pressure eigenvalue problem.
    """

    def compute_correctors( self, evp, state0, ts, dump_name, save_name,
                            problem = None, vec_g = None ):
        problem = get_default( problem, self.problem )

        problem.set_equations( self.equations )

        problem.select_bcs(ebc_names = self.ebcs, epbc_names = self.epbcs,
                           lcbc_names = self.get_default_attr('lcbcs', []))

        variables = problem.get_variables()
        get_state = variables.get_state_part_view

        eigs = evp.evp.eigs
        mtx_q = evp.evp.mtx_q
        mtx = evp.mtx

        nr, nc = mtx_q.shape

        if vec_g is not None:
            output( 'nonzero pressure EBC: max = %e, min = %e' \
                    % (vec_g.max(), vec_g.min()) )
            one = nm.ones( (nc,), dtype = nm.float64 )

        ##
        # follow_epbc = False -> R1 = - R2 as required. ? for other correctors?
        sstate0 = variables.strip_state_vector(state0, follow_epbc=False)
        vu, vp = self.dump_variables
        vec_p0 = get_state( sstate0, vp, True )
##         print state0
##         print vec_p0
##         print vec_p0.min(), vec_p0.max(), nla.norm( vec_p0 )
##         debug()

        # xi0 = Q^{-1} p(0) = Q^T G p(0)
        vec_xi0 = sc.dot( mtx_q.T, sc.dot( mtx['G'],
                                           vec_p0[:,nm.newaxis] ) ).squeeze()
        action_aibt = mtx['action_aibt']

        e_e_qg = 0.0
        iee_e_qg = 0.0
        format = '====== time %%e (step %%%dd of %%%dd) ====='\
                 % ((ts.n_digit,) * 2)
        for step, time in ts:
            output( format % (time, step + 1, ts.n_step) )

            e_e = nm.exp( - eigs * time )
            e_e_qp = e_e * vec_xi0 # exp(-Et) Q^{-1} p(0)

            if vec_g is not None:
                Qg = sc.dot( mtx_q.T, vec_g )
                e_e_qg = e_e * Qg
                iee_e_qg = ((one - e_e) / eigs) * Qg

            vec_p = sc.dot( mtx_q, e_e_qp + iee_e_qg )
            vec_dp = - sc.dot( mtx_q, (eigs * e_e_qp - e_e_qg) )
            vec_u = action_aibt( vec_dp )
##             bbb = sc.dot( vec_dp.T, - mtx['C'] * vec_p0 )

            vec_u = variables[vu].get_full(vec_u)
            vec_p = variables[vp].get_full(vec_p)
            # BC nodes - time derivative of constant is zero!
            vec_dp = variables[vp].get_full(vec_dp, force_value=0.0)
##             aaa = sc.dot( vec_xi0.T, eigs * (eigs * e_e_qp) )
##             print aaa
##             print bbb

            self.save( dump_name, save_name, vec_u, vec_p, vec_dp, ts, problem )

    def save( self, dump_name, save_name, vec_u, vec_p, vec_dp, ts, problem ):
        """
        1. saves raw correctors into hdf5 files (filename)
        2. saves correctors transformed to output for visualization
        """
        vu, vp = self.dump_variables
        out = {vu : Struct( name = 'dump', mode = 'nodes', data = vec_u,
                            dofs = None, var_name = vu ),
               vp : Struct( name = 'dump', mode = 'nodes', data = vec_p,
                            dofs = None, var_name = vp ),
               'd'+vp : Struct( name = 'dump', mode = 'nodes', data = vec_dp,
                                dofs = None, var_name = vp )}

        problem.save_state( dump_name, out = out, file_per_var = False, ts = ts )

        # For visualization...
        out = {}
        extend = not self.file_per_var

        variables = self.problem.get_variables()
        to_output = variables.state_to_output

        out.update( to_output( vec_u, var_info = {vu : (True, vu)},
                               extend = extend ) )
        out.update( to_output( vec_p, var_info = {vp : (True, vp)},
                               extend = extend ) )
        out.update( to_output( vec_dp, var_info = {vp : (True, 'd'+vp)},
                               extend = extend ) )
        if self.post_process_hook is not None:
            out = self.post_process_hook( out, problem,
                                          {vu : vec_u,
                                           vp : vec_p, 'd'+vp : vec_dp},
                                          extend = extend )
        problem.save_state( save_name, out = out,
                            file_per_var = self.file_per_var, ts = ts )

    def verify_correctors( self, initial_state, filename, problem = None ):

        problem = get_default( problem, self.problem )

        problem.set_equations( self.verify_equations )

        io = HDF5MeshIO( filename )
        ts = TimeStepper( *io.read_time_stepper() )

        variables = self.problem.get_variables()
        get_state = variables.get_state_part_view

        vu, vp = self.dump_variables
        vdp = self.verify_variables[-1]

        p0 = get_state( initial_state, vp )

        format = '====== time %%e (step %%%dd of %%%dd) ====='\
                 % ((ts.n_digit,) * 2)
        vv = variables
        ok = True
        for step, time in ts:
            output( format % (time, step + 1, ts.n_step) )

            data = io.read_data( step )
            if step == 0:
                assert_( nm.allclose( data[vp].data, p0 ) )

            problem.select_bcs(ebc_names = self.ebcs, epbc_names = self.epbcs,
                               lcbc_names = self.get_default_attr('lcbcs', []))

            state0 = problem.create_state()
            state0.set_full(data[vu].data, vu)
            state0.set_full(data[vp].data, vp)
            vv[vdp].data_from_data( data['d'+vp].data )

            state = problem.solve( state0 = state0, ts = ts )
            state, state0 = state(), state0()
            err = nla.norm( state - state0 )
            output( state.min(), state.max() )
            output( state0.min(), state0.max() )
            output( '>>>>>', err )

            ok = ok and (err < 1e-15)
            problem.advance( ts )

        return ok

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

        vf = {}
        for region_name in self.regions:
            vkey = 'volume_%s' % region_name
            key = 'fraction_%s' % region_name

            equations, variables = problem.create_evaluable(self.expression % region_name)
            val = eval_equations(equations, variables)

            vf[vkey] = nm.asarray( val, dtype = nm.float64 )
            vf[key] = vf[vkey] / volume

        return vf

class CoefSymSym( MiniAppBase ):

    def set_variables_default(variables, ir, ic, mode, set_var, data):
        mode2var = {'row' : 0, 'col' : 1}
        idx = mode2var[mode]

        val = data[set_var[idx][1]].states[ir, ic][set_var[idx][2]]
        variables[set_var[idx][0]].data_from_any(val)

    set_variables_default = staticmethod(set_variables_default)

    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        dim, sym = problem.get_dim( get_sym = True )
        coef = nm.zeros( (sym, sym), dtype = nm.float64 )

        equations, variables = problem.create_evaluable(self.expression)

        for ir, (irr, icr) in enumerate( iter_sym( dim ) ):
            if isinstance(self.set_variables, list):
                self.set_variables_default(variables, irr, icr, 'row',
                                           self.set_variables, data)
            else:
                self.set_variables(variables, irr, icr, 'row', **data)

            for ic, (irc, icc) in enumerate( iter_sym( dim ) ):
                if isinstance(self.set_variables, list):
                    self.set_variables_default(variables, irr, icr, 'col',
                                               self.set_variables, data)
                else:
                    self.set_variables(variables, irc, icc, 'col', **data)

                val = eval_equations(equations, variables)

                coef[ir,ic] = val

        coef /= volume

        return coef

class CoefFMSymSym( MiniAppBase ):
    """
    Base class for fading memory dim x dim coefficients.
    """
    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        dim, sym = problem.get_dim( get_sym = True )

        aux = self.get_filename( data, 0, 0 )
        ts = TimeStepper( *HDF5MeshIO( aux ).read_time_stepper() )

        coef = nm.zeros( (ts.n_step, sym, sym), dtype = nm.float64 )

        equations, variables = problem.create_evaluable(self.expression)

        for ir, (irr, icr) in enumerate( iter_sym( dim ) ):
            io = HDF5MeshIO( self.get_filename( data, irr, icr ) )

            for step, time in ts:
                self.set_variables(variables, io, step, None, None,
                                   'row', **data)

                for ic, (irc, icc) in enumerate( iter_sym( dim ) ):
                    self.set_variables(variables, None, None, irc, icc,
                                       'col', **data)

                    val = eval_equations(equations, variables)

                    coef[step,ir,ic] = val

        coef /= volume

        return coef

class CoefDimSym( MiniAppBase ):

    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        dim, sym = problem.get_dim( get_sym = True )
        coef = nm.zeros( (dim, sym), dtype = nm.float64 )

        equations, variables = problem.create_evaluable(self.expression)

        for ir in range( dim ):
            self.set_variables(variables, ir, None, 'row', **data)

            for ic, (irc, icc) in enumerate( iter_sym( dim ) ):
                self.set_variables(variables, irc, icc, 'col', **data)

                val = eval_equations(equations, variables)

                coef[ir,ic] = val

        coef /= volume

        return coef

class CoefNN( MiniAppBase ):

    def set_variables_default(variables, ir, ic, mode, set_var, data):
        mode2var = {'row' : 0, 'col' : 1}

        if mode == 'col':
            ir = ic
        idx = mode2var[mode]

        val = data[set_var[idx][1]].states[ir][set_var[idx][2]]
        variables[set_var[idx][0]].data_from_any(val)

    set_variables_default = staticmethod(set_variables_default)

    def __init__( self, name, problem, kwargs ):
        """When dim is not in kwargs, problem dimension is used."""
        MiniAppBase.__init__( self, name, problem, kwargs )
        self.set_default_attr( 'dim', problem.get_dim() )

    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        coef = nm.zeros((self.dim, self.dim), dtype = nm.float64)
        equations, variables = problem.create_evaluable(self.expression)

        if isinstance(self.set_variables, list):
            for ir in range(self.dim):
                self.set_variables_default(variables, ir, None, 'row',
                                           self.set_variables, data)
                for ic in range(self.dim):
                    self.set_variables_default(variables, None, ic, 'col',
                                               self.set_variables, data)
                    val = eval_equations(equations, variables)
                    coef[ir,ic] = val
        else:
            for ir in range(self.dim):
                self.set_variables(variables, ir, None, 'row', **data)
                for ic in range(self.dim):
                    self.set_variables(variables, None, ic, 'col', **data)
                    val = eval_equations(equations, variables)
                    coef[ir,ic] = val

        coef /= volume

        return coef

class CoefN( MiniAppBase ):

    def set_variables_default(variables, ir, set_var, data):
        for (var, req, comp) in set_var:
            variables[var].data_from_any(data[req].states[ir][comp])

    set_variables_default = staticmethod(set_variables_default)

    def __init__( self, name, problem, kwargs ):
        """When dim is not in kwargs, problem dimension is used."""
        MiniAppBase.__init__( self, name, problem, kwargs )
        self.set_default_attr( 'dim', problem.get_dim() )

    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        coef = nm.zeros((self.dim,), dtype = nm.float64)
        equations, variables = problem.create_evaluable(self.expression)

        for ir in range(self.dim):
            if isinstance(self.set_variables, list):
                self.set_variables_default(variables, ir, self.set_variables,
                                           data)
            else:
                self.set_variables(variables, ir, **data)

            val = eval_equations(equations, variables)

            coef[ir] = val

        coef /= volume

        return coef

class CoefDimDim( CoefNN ):
    pass

class CoefDim( CoefN ):
    pass

class CoefSym( MiniAppBase ):

    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        dim, sym = problem.get_dim( get_sym = True )
        coef = nm.zeros( (sym,), dtype = nm.float64 )

        equations, variables = problem.create_evaluable(self.expression)

        self.set_variables(variables, None, None, 'col', **data)

        for ii, (ir, ic) in enumerate( iter_sym( dim ) ):
            self.set_variables(variables, ir, ic, 'row', **data)

            val = eval_equations(equations, variables)
            coef[ii] = val

        coef /= volume

        return coef

class CoefFMSym( MiniAppBase ):

    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        dim, sym = problem.get_dim( get_sym = True )

        aux = self.get_filename( data, 0, 0 )
        ts = TimeStepper( *HDF5MeshIO( aux ).read_time_stepper() )

        coef = nm.zeros( (ts.n_step, sym), dtype = nm.float64 )

        equations, variables = problem.create_evaluable(self.expression)

        self.set_variables(variables, None, None, 'col', **data)

        for ii, (ir, ic) in enumerate( iter_sym( dim ) ):
            io = HDF5MeshIO( self.get_filename( data, ir, ic ) )
            for step, time in ts:
                self.set_variables(variables, io, step, 'row', **data)

                val = eval_equations(equations, variables)

                coef[step,ii] = val

        coef /= volume

        return coef

class CoefOne( MiniAppBase ):

    def set_variables_default(variables, set_var, data):
        for (var, req, comp) in set_var:
            variables[var].data_from_any(data[req].state[comp])

    set_variables_default = staticmethod(set_variables_default)

    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )
        equations, variables = problem.create_evaluable(self.expression)

        if isinstance(self.set_variables, list):
            self.set_variables_default(variables, self.set_variables,
                                       data)
        else:
            self.set_variables(variables, **data)
        val = eval_equations(equations, variables)

        coef = val / volume

        return coef

class CoefFMOne( MiniAppBase ):
    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        aux = self.get_filename( data )
        io = HDF5MeshIO( self.get_filename( data ) )
        ts = TimeStepper( *io.read_time_stepper() )

        coef = nm.zeros( (ts.n_step, 1), dtype = nm.float64 )

        equations, variables = problem.create_evaluable(self.expression)

        self.set_variables(variables, None, None, 'col', **data)

        for step, time in ts:
            self.set_variables(variables, io, step, 'row', **data)

            val = eval_equations(equations, variables)

            coef[step] = val

        coef /= volume

        return coef

class CoefSum( MiniAppBase ):

    def __call__( self, volume, problem = None, data = None ):

        coef = nm.zeros_like(data[self.requires[0]])
        for i in range(len(self.requires)):
            coef += data[self.requires[i]]

        return coef

class CoefEval( MiniAppBase ):
    """
    Evaluate expression.
    """
    def __call__( self, volume, problem = None, data = None ):

        expr = self.expression
        for i in range( len(self.requires) ):
            expr = expr.replace(self.requires[i],
                                "data['%s']" % self.requires[i])

        coef = eval(expr)

        return coef
