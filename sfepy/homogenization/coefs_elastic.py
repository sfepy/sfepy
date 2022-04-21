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
        vec = eqs.create_reduced_vec()
        eqs.time_update_materials(problem.get_timestepper())
        eqs.evaluate(mode='weak', dw_mode='vector', asm_obj=vec)

        return vec

class TCorrectorsRSViaPressureEVP( TCorrectorsViaPressureEVP ):

    def __call__( self, problem = None, data = None ):
        """data: corrs_rs, evp"""
        problem = get_default( problem, self.problem )
        self.init_solvers(problem)
        ts = problem.get_timestepper()

        corrs, evp = [data[ii] for ii in self.requires]

        assert_( evp.ebcs == self.ebcs )
        assert_( evp.epbcs == self.epbcs )

        dim = problem.get_dim()

        self.setup_equations(self.equations)

        solve = self.compute_correctors
        states = nm.zeros((dim, dim), dtype=object)
        clist = []
        for ir in range( dim ):
            for ic in range( dim ):
                states[ir,ic] = solve(evp, -1.0, corrs.states[ir,ic], ts)
                clist.append((ir, ic))

        corr_sol = CorrSolution(name=self.name,
                                states=states,
                                n_step=ts.n_step,
                                components=clist)

        self.save(corr_sol, problem, ts)

        return corr_sol


class TCorrectorsPressureViaPressureEVP( TCorrectorsViaPressureEVP ):

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

        self.setup_equations(self.equations)

        solve = self.compute_correctors
        state = solve(evp, 1.0, corrs.state, ts, vec_g=vec_g)

        corr_sol = CorrSolution(name=self.name,
                                state=state,
                                n_step=ts.n_step)

        self.save(corr_sol, problem, ts)

        return corr_sol
