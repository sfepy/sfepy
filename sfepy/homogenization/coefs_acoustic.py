import numpy as nm

from sfepy.base.base import assert_, get_default, Struct
from sfepy.fem import eval_term_op
from sfepy.homogenization.coefs_base import CorrMiniApp, MiniAppBase

class CorrVector( CorrMiniApp ):

    def get_variables( self, problem, i, data ):
        yield ( self.variables[i], data[self.requires[i]].states )

    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        problem.select_variables( self.variables )
        problem.set_equations( self.equations )

        problem.select_bcs( ebc_names = self.ebcs, epbc_names = self.epbcs )

        self.init_solvers(problem)

        for ival in range( len( self.requires ) ):
            for name, val in self.get_variables( problem, ival, data ):
                problem.variables[name].data_from_data( val )

        state = problem.create_state()
        state.apply_ebc()
        state = problem.solve()
        assert_(state.has_ebc())

        self.save( state, problem )

        ostate = problem.state_to_output(state)
        ostate = ostate[ostate.keys()[0]].data
        ndim = ostate.shape[1]

        states = nm.zeros( (ndim,), dtype = nm.object )

        clist = []
        for ir in range( ndim ):
            states[ir] = ostate[:,ir]
            clist.append( (ir,) )

        return Struct( name = self.name,
                       states = states )

class CoefGammaT( MiniAppBase ):

    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )

        coef = nm.zeros( (1,), dtype = nm.float64 )

        ldata = data[self.requires[0]]
        var0 = problem.variables[0]

        ndata = nm.zeros( (var0.n_nod,), dtype=nm.float64 )
        for ivar in ldata.di.vnames:
            ndata += ldata.states[ldata.di.indx[ivar]] 

        problem.variables[0].data_from_data( ndata )

        val = eval_term_op( None, self.expression,
                            problem, call_mode = 'd_eval' )

        coef = val / volume

        return coef
