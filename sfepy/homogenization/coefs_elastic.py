from sfepy.base.base import *
from sfepy.homogenization.coefs_base import CoefSymSym, CoefSym, CorrDimDim,\
     CoefOne, CorrOne, CorrDim, CoefDimDim, ShapeDimDim,\
     PressureEigenvalueProblem, TCorrectorsViaPressureEVP,\
     CoefFMSymSym, CoefFMSym, TSTimes, VolumeFractions

class CorrectorsElasticRS( CorrDimDim ):

    def get_variables( self, ir, ic, data ):
        """data: pis"""
        pis = data[self.requires[0]]
        yield (self.variables[-1], pis[ir,ic])

class CorrectorsPressure( CorrOne ):

    def get_variables( self, data ):
        """data: None"""
        vv = self.problem.variables

        name_y = self.variables[-2]
        name_ym = self.variables[-1]

        one = nm.ones( (vv[name_y].field.n_nod,), dtype = nm.float64 )
        yield name_y, one

        one_m = nm.ones( (vv[name_ym].field.n_nod,), dtype = nm.float64 )
        yield name_ym, one_m

class CorrectorsPermeability( CorrDim ):

    def __call__( self, problem = None, data = None, save_hook = None ):
        problem = get_default( problem, self.problem )

        problem.select_variables( self.variables )
        equations = {}
        for key, eq in self.equations.iteritems():
            equations[key] = eq % tuple( self.regions )
        problem.set_equations( equations )

        problem.select_bcs( ebc_names = self.ebcs, epbc_names = self.epbcs )

        dim = problem.get_dim()
        states = nm.zeros( (dim,), dtype = nm.object )
        for ir in range( dim ):
            state = problem.solve( ir = ir )
            assert_( problem.variables.has_ebc( state ) )
            states[ir] = state

            if save_hook is not None:
                save_hook( state, problem, ir )

        return Struct( name = self.name,
                       states = states,
                       di = problem.variables.di )

    def make_save_hook( self, base_name, format,
                        post_process_hook = None, file_per_var = None ):
        return CorrDim.make_save_hook( self, base_name, format,
                                       post_process_hook = None,
                                       file_per_var = None )

class TCorrectorsRSViaPressureEVP( TCorrectorsViaPressureEVP ):

    def __call__( self, problem = None, data = None, save_hook = None ):
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
                filename = self.save_name + ("_%d%d" % (ir,ic)) + '.h5'
                solve( evp, -corrs.states[ir,ic], ts, filename, save_hook )
                filenames[ir,ic] = os.path.join( problem.output_dir, filename )

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


class ElasticCoef( CoefSymSym ):
    """Homogenized elastic tensor $E_{ijkl}$."""

    mode2var = {'row' : 0, 'col' : 1}

    def get_variables( self, problem, ir, ic, data, mode ):

        pis, corrs = [data[ii] for ii in self.requires]

        var_name = self.variables[self.mode2var[mode]]
        c_name = problem.variables[var_name].primary_var_name

        indx = corrs.di.indx[c_name]
        omega = corrs.states[ir,ic][indx]
        pi = pis[ir,ic] + omega
        yield (var_name, pi)

class ViscousFMCoef( CoefFMSymSym ):
    """Homogenized viscous fading memory tensor $H_{ijkl}$."""
 
    def get_filename( self, data, ir, ic ):
        tcorrs = data[self.requires[1]]
        return tcorrs.filenames[ir,ic]

    def get_variables( self, problem, io, step, ir, ic, data, mode ):

        corrs = data[self.requires[0]]

        if mode == 'row':
            dpc = io.read_data( step )['dp'].data
            yield self.variables[0], dpc

        else:
            var_name = self.variables[1]
            c_name = problem.variables[var_name].primary_var_name
            indx = corrs.di.indx[c_name]
            pc = corrs.states[ir,ic][indx]
            yield var_name, pc

class ElasticBiotCoef( CoefSym ):
    """Homogenized elastic Biot coefficient."""

    def get_variables( self, problem, ir, ic, data, mode ):

        if mode == 'col':
            var_name = self.variables[0]
            one_var = problem.variables[var_name]
            one = nm.ones( (one_var.field.n_nod,), dtype = nm.float64 )
            yield var_name, one

        else:
            var_name = self.variables[1]
            c_name = problem.variables[var_name].primary_var_name

            corrs = [data[ii] for ii in self.requires][0]
            indx = corrs.di.indx[c_name]
            omega = corrs.states[ir,ic][indx]
            yield var_name, omega

class BiotFMCoef( CoefFMSym ):
    """Fading memory Biot coefficient."""

    def get_filename( self, data, ir, ic ):
        tcorrs = data[self.requires[0]]
        return tcorrs.filenames[ir,ic]

    def get_variables( self, problem, io, step, data, mode ):

        if mode == 'col':
            var_name = self.variables[0]
            one_var = problem.variables[var_name]
            one = nm.ones( (one_var.field.n_nod,), dtype = nm.float64 )
            yield var_name, one

            var_name = self.variables[1]
            one_var = problem.variables[var_name]
            one_m = nm.ones( (one_var.field.n_nod,), dtype = nm.float64 )
            yield var_name, one_m

        else:
            step_data = io.read_data( step )

            yield self.variables[2], step_data['u'].data
            yield self.variables[3], step_data['dp'].data

class IRBiotModulus( CoefOne ):
    """Homogenized instantaneous reciprocal Biot modulus."""

    def get_variables( self, problem, data ):

        var_name = self.variables[0]
        one_var = problem.variables[var_name]
        one = nm.ones( (one_var.field.n_nod,), dtype = nm.float64 )
        yield var_name, one

        var_name = self.variables[1]
        one_var = problem.variables[var_name]
        one_m = nm.ones( (one_var.field.n_nod,), dtype = nm.float64 )
        yield var_name, one_m
        
        corrs = [data[ii] for ii in self.requires][0]

        for ii in [2, 3]:
            # omega and pp.
            var_name = self.variables[ii]
            c_name = problem.variables[var_name].primary_var_name
            indx = corrs.di.indx[c_name]
            val = corrs.state[indx]
            yield var_name, val

class DiffusionCoef( CoefDimDim ):
    """Homogenized diffusion coefficient."""

    mode2var = {'row' : 0, 'col' : 1}

    def get_variables( self, problem, ir, ic, data, mode ):

        pressure = problem.variables[self.variables[0]]
        coor = pressure.field.get_coor()
        
        corrs = [data[ii] for ii in self.requires][0]

        var_name = self.variables[self.mode2var[mode]]
        c_name = problem.variables[var_name].primary_var_name
        indx = corrs.di.indx[c_name]

        if mode == 'row':
            ii = ir
        else:
            ii = ic

        val = coor[:,ii] + corrs.states[ii][indx]

        yield var_name, val
