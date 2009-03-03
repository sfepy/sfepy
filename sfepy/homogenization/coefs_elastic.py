from sfepy.base.base import *
from sfepy.homogenization.coefs_base import CoefSymSym, CoefSym, CorrDimDim,\
     CoefOne, CorrOne, CorrDim, CoefDimDim, ShapeDimDim,\
     PressureEigenvalueProblem, TCorrectorsViaPressureEVP,\
     CoefFMSymSym, CoefFMSym, CoefFMOne, TSTimes, VolumeFractions

def generate_ones( problem, var_names ):
    for var_name in var_names:
        one_var = problem.variables[var_name]
        one = nm.ones( (one_var.field.n_nod,), dtype = nm.float64 )
        yield var_name, one

class CorrectorsElasticRS( CorrDimDim ):

    def get_variables( self, ir, ic, data ):
        """data: pis"""
        pis = data[self.requires[0]]
        yield (self.variables[-1], pis[ir,ic])

class CorrectorsPressure( CorrOne ):

    def get_variables( self, data ):
        """data: None"""
        return generate_ones( self.problem, self.variables[-2:] )

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

            self.save( state, problem, ir )

        return Struct( name = self.name,
                       states = states,
                       di = problem.variables.di )


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
        """data: corrs_pressure, evp"""
        problem = get_default( problem, self.problem )
        ts = problem.get_time_solver().ts

        corrs, evp = [data[ii] for ii in self.requires]

        assert_( evp.ebcs == self.ebcs )
        assert_( evp.epbcs == self.epbcs )

        filename = self.get_dump_name()
        savename = self.get_save_name()

        solve = self.compute_correctors
        solve( evp, corrs.state, ts, filename, savename )

        if self.check:
            output( 'verifying correctors %s...' % self.name )
            verify = self.verify_correctors
            ok = verify( corrs.state, filename )
            output( '...done, ok: %s' % ok )

        return Struct( name = self.name,
                       filename = filename )

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

class ElasticBiotCoef( CoefSym ):
    """Homogenized elastic Biot coefficient."""

    def get_variables( self, problem, ir, ic, data, mode ):

        if mode == 'col':
            for var_name, val in  generate_ones( self.problem,
                                                 self.variables[0:1] ):
                yield var_name, val

        else:
            var_name = self.variables[1]
            c_name = problem.variables[var_name].primary_var_name

            corrs = data[self.requires[0]]
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

class IRBiotModulus( CoefOne ):
    """Homogenized instantaneous reciprocal Biot modulus."""

    def get_variables( self, problem, data ):

        for var_name, val in  generate_ones( self.problem,
                                             self.variables[0:2] ):
            yield var_name, val
        
        corrs = data[self.requires[0]]

        for ii in [2, 3]:
            # omega and pp.
            var_name = self.variables[ii]
            c_name = problem.variables[var_name].primary_var_name
            indx = corrs.di.indx[c_name]
            val = corrs.state[indx]
            yield var_name, val

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

class DiffusionCoef( CoefDimDim ):
    """Homogenized diffusion coefficient."""

    mode2var = {'row' : 0, 'col' : 1}

    def get_variables( self, problem, ir, ic, data, mode ):

        pressure = problem.variables[self.variables[0]]
        coor = pressure.field.get_coor()
        
        corrs = data[self.requires[0]]

        var_name = self.variables[self.mode2var[mode]]
        c_name = problem.variables[var_name].primary_var_name
        indx = corrs.di.indx[c_name]

        if mode == 'row':
            ii = ir
        else:
            ii = ic

        val = coor[:,ii] + corrs.states[ii][indx]

        yield var_name, val
