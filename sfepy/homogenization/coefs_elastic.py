import numpy as nm

from sfepy.base.base import output, assert_, get_default, Struct
from sfepy.homogenization.coefs_base import CoefOne, CorrDim, \
     TCorrectorsViaPressureEVP, \
     CoefFMSymSym, CoefFMSym, CoefFMOne, \
     CorrMiniApp, CorrSolution
from sfepy.fem.meshio import HDF5MeshIO
from sfepy.solvers.ts import TimeStepper

class CorrectorsPermeability( CorrDim ):

    def __call__( self, problem = None, data = None, save_hook = None ):
        problem = get_default( problem, self.problem )

        equations = {}
        for key, eq in self.equations.iteritems():
            equations[key] = eq % tuple( self.regions )

        index = [0]
        problem.set_equations( equations, user={self.index_name : index} )

        problem.select_bcs( ebc_names = self.ebcs, epbc_names = self.epbcs )
        problem.update_materials(problem.ts)

        self.init_solvers(problem)

        variables = problem.get_variables()

        dim = problem.get_dim()
        states = nm.zeros( (dim,), dtype = nm.object )
        clist = []
        for ir in range( dim ):
            index[0] = ir # Set the index - this is visible in the term.

            state = problem.solve()
            assert_(state.has_ebc())
            states[ir] = variables.get_state_parts()

            clist.append( (ir,) )

        corr_sol = CorrSolution(name = self.name,
                                states = states,
                                components = clist)

        self.save(corr_sol, problem)

        return corr_sol

class PressureRHSVector( CorrMiniApp ):
    
    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        problem.select_variables( self.variables )
        problem.set_equations( self.equations )
        problem.select_bcs( ebc_names = self.ebcs, epbc_names = self.epbcs )

        state = problem.create_state()
        state.apply_ebc()

        vec = eval_term_op( state, self.equations.values()[0],
                            problem, dw_mode = 'vector' )
##         print vec.max(), vec.min()

        return vec

class TCorrectorsRSViaPressureEVP( TCorrectorsViaPressureEVP ):

    def get_save_name_base( self ):
        return self.save_name + '_%d%d'

    def get_dump_name_base( self ):
        return self.dump_name + '_%d%d'

    def __call__( self, problem = None, data = None ):
        """data: corrs_rs, evp"""
        problem = get_default( problem, self.problem )
        ts = problem.get_time_solver().ts

        corrs, evp = [data[ii] for ii in self.requires]

        assert_( evp.ebcs == self.ebcs )
        assert_( evp.epbcs == self.epbcs )

        dim = problem.get_dim()
        
        filenames = nm.zeros( (dim, dim), dtype = nm.object )

        solve = self.compute_correctors
        for ir in range( dim ):
            for ic in range( dim ):
                filenames[ir,ic] = self.get_dump_name() % (ir,ic)
                savename = self.get_save_name() % (ir,ic)
                solve( evp, -corrs.states[ir,ic], ts,
                       filenames[ir,ic], savename )


        if self.check:
            output( 'verifying correctors %s...' % self.name )
            verify = self.verify_correctors
            ok = True
            for ir in range( dim ):
                for ic in range( dim ):
                    oo = verify( -corrs.states[ir,ic], filenames[ir,ic] )
                    ok = ok and oo
            output( '...done, ok: %s' % ok )

        return Struct( name = self.name,
                       filenames = filenames )

class TCorrectorsPressureViaPressureEVP( TCorrectorsViaPressureEVP ):

    def get_save_name_base( self ):
        return self.save_name

    def get_dump_name_base( self ):
        return self.dump_name

    def __call__( self, problem = None, data = None, save_hook = None ):
        """data: corrs_pressure, evp, optionally vec_g"""
        problem = get_default( problem, self.problem )
        ts = problem.get_time_solver().ts

        corrs, evp = [data[ii] for ii in self.requires[:2]]
        if len(self.requires) == 3:
            vec_g = data[self.requires[2]]
        else:
            vec_g = None

        assert_( evp.ebcs == self.ebcs )
        assert_( evp.epbcs == self.epbcs )

        filename = self.get_dump_name()
        savename = self.get_save_name()

        solve = self.compute_correctors
        solve( evp, corrs.state, ts, filename, savename, vec_g = vec_g )

        if self.check:
            output( 'verifying correctors %s...' % self.name )
            verify = self.verify_correctors
            ok = verify( corrs.state, filename )
            output( '...done, ok: %s' % ok )

        return Struct( name = self.name,
                       filename = filename )

class ViscousFMCoef( CoefFMSymSym ):
    """Homogenized viscous fading memory tensor $H_{ijkl}$."""
 
    def get_filename( self, data, ir, ic ):
        tcorrs = data[self.requires[1]]
        return tcorrs.filenames[ir,ic]

    def get_variables( self, problem, io, step, ir, ic, data, mode ):

        corrs = data[self.requires[0]]

        if mode == 'row':
            var_name = self.variables[0]
            c_name = problem.variables[var_name].primary_var_name
            dpc = io.read_data( step )['d'+c_name].data
            yield var_name, dpc

        else:
            var_name = self.variables[1]
            c_name = problem.variables[var_name].primary_var_name
            indx = corrs.di.indx[c_name]
            pc = corrs.states[ir,ic][indx]
            yield var_name, pc

class RBiotCoef( CoefFMSym ):
    """Homogenized fading memory Biot-like coefficient."""

    def get_filename( self, data, ir, ic ):
        tcorrs = data[self.requires[1]]
        self.ir, self.ic = ir, ic
        return tcorrs.filename

    def get_variables( self, problem, io, step, data, mode ):

        if mode == 'col':
            return
        
        else:
            pis = data[self.requires[0]]
            step_data = io.read_data( step )

            # omega.
            var_name = self.variables[0]
            c_name = problem.variables[var_name].primary_var_name
            yield var_name, step_data[c_name].data

            var_name = self.variables[1]
            yield var_name, pis.states[self.ir,self.ic]

            var_name = self.variables[2]
            c_name = problem.variables[var_name].primary_var_name
            pc = step_data['d'+c_name].data

            if self.ir == self.ic:
                yield (var_name, pc)
            else:
                yield (var_name, nm.zeros_like(pc))

class BiotFMCoef( CoefFMSym ):
    """Fading memory Biot coefficient."""

    def get_filename( self, data, ir, ic ):
        tcorrs = data[self.requires[0]]
        return tcorrs.filenames[ir,ic]

    def get_variables( self, problem, io, step, data, mode ):

        if mode == 'col':
            for var_name, val in  generate_ones( self.problem,
                                                 self.variables[0:2] ):
                yield var_name, val

        else:
            step_data = io.read_data( step )

            var_name = self.variables[2]
            c_name = problem.variables[var_name].primary_var_name
            yield var_name, step_data[c_name].data

            var_name = self.variables[3]
            c_name = problem.variables[var_name].primary_var_name
            yield var_name, step_data['d'+c_name].data

class BiotFM2Coef( CoefFMSym ):
    """Fading memory Biot coefficient, alternative form."""

    def get_filename( self, data, ir, ic ):
        tcorrs = data[self.requires[1]]
        return tcorrs.filenames[ir,ic]

    def get_variables( self, problem, io, step, data, mode ):

        if mode == 'col':
            var_name = self.variables[0]
            c_name = problem.variables[var_name].primary_var_name

            corrs = data[self.requires[0]]
            indx = corrs.di.indx[c_name]
            p0 = corrs.state[indx]
            yield var_name, p0

        else:
            step_data = io.read_data( step )

            var_name = self.variables[1]
            c_name = problem.variables[var_name].primary_var_name
            yield var_name, step_data[c_name].data

class GBarCoef( CoefOne ):
    """
    Asymptotic Barenblatt coefficient.

    data = [p^{\infty}]

    Note:

    solving "dw_diffusion.i1.Y3( m.K, qc, pc ) = 0" solve, in fact
    "C p^{\infty} = \hat{C} \hat{\pi}" with the result "\hat{p^{\infty}}",
    where the rhs comes from E(P)BC.
    - it is preferable to computing directly by
    "\hat{p^{\infty}} = \hat{C^-1 \strip(\hat{C} \hat{\pi})}", as it checks
    explicitly the rezidual.
    """

    def __call__( self, volume, problem = None, data = None ):
        expression, region_name = self.expression

        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )
        problem.set_equations( {'eq' : expression} )
        problem.time_update( conf_ebc = {}, conf_epbc = {}, conf_lcbc = {} )

        pi_inf = data[self.requires[0]]

        coef = nm.zeros( (1,), dtype = nm.float64 )

        vec = eval_term_op( pi_inf.state, expression,
                            problem, new_geometries = False, dw_mode = 'vector' )

        reg = problem.domain.regions[region_name]
        field = problem.variables[self.variables[1]].field
        nods = field.get_dofs_in_region(reg, merge=True)
        coef[0] = vec[nods].sum()

        coef /= volume

        return coef

def eval_boundary_diff_vel_grad( problem, uc, pc, equation, region_name,
                                 pi = None ):
    get_state = problem.variables.get_state_part_view

    problem.set_equations( {'eq_2' : equation} )

    state = problem.create_state()
    state.set_full(uc, 'uc')
    state.set_full(pc, 'pc')
    if pi is not None:
        problem.variables['Pi'].data_from_data( pi )

    problem.time_update( conf_ebc = {}, conf_epbc = {} )

    state.apply_ebc()
    aux = problem.get_evaluator().eval_residual(state())

    pc = get_state( aux, 'pc', True )
    pc = problem.variables.make_full_vec( pc, 'pc', 0 )

    field = problem.variables['pc'].field

    reg = problem.domain.regions[region_name]
    nods = field.get_dofs_in_region(reg, merge=True)
    val = pc[nods].sum()

#    assert nm.all( problem.variables.di.ptr == problem.variables.adi.ptr )
#    problem.time_update() # Restore EBC.

    return val

class GPlusCoef( CoefFMOne ):

    def get_filename( self, data ):
        tcorrs = data[self.requires[0]]
        return tcorrs.filename

    def __call__( self, volume, problem = None, data = None ):
        expression, region_name, aux_eq = self.expression

        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )

        aux = self.get_filename( data )
        io = HDF5MeshIO( self.get_filename( data ) )
        ts = TimeStepper( *io.read_time_stepper() )

        coef = nm.zeros( (ts.n_step, 1), dtype = nm.float64 )
        
        for step, time in ts:
            step_data = io.read_data( step )

            var_name = self.variables[0]
            c_name = problem.variables[var_name].primary_var_name
            omega = step_data[c_name].data
            problem.variables[var_name].data_from_data( omega )

            val1 = eval_term_op( None, expression,
                                problem, call_mode = 'd_eval' )

            pc = step_data[self.variables[3]].data
            val2 = eval_boundary_diff_vel_grad( problem, omega, pc, aux_eq,
                                                region_name )
            coef[step] = val1 + val2

        coef /= volume

        return coef

class FMRBiotModulus( CoefFMOne ):
    """Fading memory reciprocal Biot modulus."""

    def get_filename( self, data ):
        tcorrs = data[self.requires[0]]
        return tcorrs.filename
    
    def get_variables( self, problem, io, step, data, mode ):

        if mode == 'col':
            for var_name, val in  generate_ones( self.problem,
                                                 self.variables[0:2] ):
                yield var_name, val

        else:
            step_data = io.read_data( step )

            var_name = self.variables[2]
            c_name = problem.variables[var_name].primary_var_name
            yield var_name, step_data[c_name].data

            var_name = self.variables[3]
            c_name = problem.variables[var_name].primary_var_name
            yield var_name, step_data['d'+c_name].data
