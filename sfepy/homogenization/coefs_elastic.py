from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import output, assert_, get_default, Struct
from sfepy.homogenization.coefs_base import CorrSolution, \
     TCorrectorsViaPressureEVP, CorrMiniApp
from sfepy.solvers.ts import TimeStepper
from six.moves import range

class PressureRHSVector( CorrMiniApp ):
    
    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        problem.select_variables( self.variables )
        problem.set_equations( self.equations )
        problem.select_bcs(ebc_names = self.ebcs, epbc_names = self.epbcs,
                           lcbc_names=self.get('lcbcs', []))

        state = problem.create_state()
        state.apply_ebc()

        eqs = problem.equations
        eqs.set_variables_from_state(state.vec)
        vec = eqs.create_stripped_state_vector()
        eqs.time_update_materials(problem.get_timestepper())
        eqs.evaluate(mode='weak', dw_mode='vector', asm_obj=vec)

        return vec

class TCorrectorsRSViaPressureEVP( TCorrectorsViaPressureEVP ):

    def get_save_name_base( self ):
        return self.save_name + '_%d%d'

    def get_dump_name_base( self ):
        return self.dump_name + '_%d%d'

    def __call__( self, problem = None, data = None ):
        """data: corrs_rs, evp"""
        problem = get_default( problem, self.problem )
        self.init_solvers(problem)
        ts = problem.get_timestepper()

        corrs, evp = [data[ii] for ii in self.requires]

        assert_( evp.ebcs == self.ebcs )
        assert_( evp.epbcs == self.epbcs )

        dim = problem.get_dim()

        filenames = nm.zeros( (dim, dim), dtype = nm.object )

        self.setup_equations(self.equations)

        solve = self.compute_correctors
        states = nm.zeros((dim, dim), dtype=nm.object)
        clist = []
        for ir in range( dim ):
            for ic in range( dim ):
                filenames[ir,ic] = self.get_dump_name() % (ir,ic)
                savename = self.get_save_name() % (ir,ic)
                states[ir,ic] = solve(evp, -1.0, corrs.states[ir,ic], ts,
                                filenames[ir,ic], savename)
                clist.append((ir, ic))

        corr_sol = CorrSolution(name=self.name,
                                states=states,
                                n_tstep=ts.n_step,
                                components=clist)

        if self.check:
            self.setup_equations(self.verify_equations)
            self.init_solvers(problem)

            output( 'verifying correctors %s...' % self.name )
            verify = self.verify_correctors
            ok = True
            for ir in range( dim ):
                for ic in range( dim ):
                    oo = verify(-1.0, corrs.states[ir,ic], filenames[ir,ic])
                    ok = ok and oo
            output( '...done, ok: %s' % ok )

        return corr_sol


class TCorrectorsPressureViaPressureEVP( TCorrectorsViaPressureEVP ):

    def get_save_name_base( self ):
        return self.save_name

    def get_dump_name_base( self ):
        return self.dump_name

    def __call__( self, problem = None, data = None, save_hook = None ):
        """data: corrs_pressure, evp, optionally vec_g"""
        problem = get_default( problem, self.problem )
        self.init_solvers(problem)
        ts = problem.get_timestepper()

        corrs, evp = [data[ii] for ii in self.requires[:2]]
        if len(self.requires) == 3:
            vec_g = data[self.requires[2]]
        else:
            vec_g = None

        assert_( evp.ebcs == self.ebcs )
        assert_( evp.epbcs == self.epbcs )

        filename = self.get_dump_name()
        savename = self.get_save_name()

        self.setup_equations(self.equations)

        solve = self.compute_correctors
        state = solve(evp, 1.0, corrs.state, ts, filename, savename,
                      vec_g=vec_g)

        corr_sol = CorrSolution(name=self.name,
                                state=state,
                                n_tstep=ts.n_step)

        if self.check:
            self.setup_equations(self.verify_equations)
            self.init_solvers(problem)

            output( 'verifying correctors %s...' % self.name )
            verify = self.verify_correctors
            ok = verify(1.0, corrs.state, filename)
            output( '...done, ok: %s' % ok )

        return corr_sol
