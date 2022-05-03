from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import assert_, get_default, Struct
from sfepy.homogenization.coefs_base import CorrMiniApp, CoefN
import six
from six.moves import range

class CorrRegion( CorrMiniApp ):

    def __init__( self, name, problem, kwargs ):
        CorrMiniApp.__init__( self, name, problem, kwargs )
        self.set_default('Nreg', len(list(self.regions.values())[0]))
        self.set_default('ebcs_list', False)

    def get_variables( self, ir, data ):
        return iter([])
        
    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        states = nm.zeros((self.Nreg,), dtype=object)
        clist = []
        for ir in range( self.Nreg ):

            problem.select_variables( self.variables )
            
            for name, val in self.get_variables( ir, data ):
                problem.variables[name].set_data(val)

            equations = {}
        
            for keye, vale in six.iteritems(self.equations):
                for keyr, valr in six.iteritems(self.regions):
                    vale = vale.replace( keyr, valr[ir] )
            
                equations[keye] = vale    

            problem.set_equations( equations )

            if ( self.ebcs_list ):
                problem.select_bcs( ebc_names = self.ebcs[ir], epbc_names = self.epbcs )
            else:
                problem.select_bcs( ebc_names = self.ebcs, epbc_names = self.epbcs )

            self.init_solvers(problem)

            state = problem.solve()
            assert_(state.has_ebc())
            states[ir] = state()
            clist.append( (ir,) )

        self.save( states, problem, clist )

        return Struct( name = self.name,
                       states = states,
                       di = problem.variables.di )

class CoefRegion( CoefN ):

    def __init__( self, name, problem, kwargs ):
        CoefN.__init__( self, name, problem, kwargs )
        self.corr_dim = len(list(kwargs['regions'].values())[0])
        
    def get_variables( self, problem, ir, data ):

        corr = data[self.requires[-1]]
        yield (self.variables[0], corr.states[ir])

