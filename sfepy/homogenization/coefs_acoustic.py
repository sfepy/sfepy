from sfepy.base.base import *
from sfepy.fem import eval_term_op
from sfepy.homogenization.coefs_base import VolumeFractions, \
     CorrMiniApp, MiniAppBase, CorrN, CoefN, CoefNN, \
     create_scalar_pis

class ShapeDimM1( CorrMiniApp ):
    
    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        pis = create_scalar_pis( problem, self.variables[0] )

        return Struct( name = self.name,
                       states = pis[0:-1].copy() )

class CorrDimM1( CorrN ):

    def __init__( self, name, problem, kwargs ):
        CorrN.__init__( self, name, problem, kwargs )
        self.corr_dim -= 1
        
    def get_variables( self, ir, data ):
        """data: pis"""
        pis = data[self.requires[0]]

        yield (self.variables[0], pis.states[ir])

class CorrScalar( CorrMiniApp ):

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
            
        state = problem.create_state_vector()
        problem.apply_ebc( state )
        state = problem.solve()
        assert_( problem.variables.has_ebc( state ) )

        self.save( state, problem )

        return Struct( name = self.name,
                       states = state,
                       di = problem.variables.di )

class CorrVector( CorrScalar ):

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
            
        state = problem.create_state_vector()
        problem.apply_ebc( state )
        state = problem.solve()
        assert_( problem.variables.has_ebc( state ) )

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

class CoefDimM1PDimM1P( CoefNN ):

    mode2var = {'row' : 0, 'col' : 1}
    
    def __init__( self, name, problem, kwargs ):
        CoefNN.__init__( self, name, problem, kwargs )
        self.corr_dim -= 1
    
    def get_variables( self, problem, ir, ic, data, mode ):

        corr = data[self.requires[-1]]
        pis = data[self.requires[0]]
        
        var_name = self.variables[self.mode2var[mode]]
        c_name = problem.variables[var_name].primary_var_name
        indx = corr.di.indx[c_name]

        if mode =='col':
            ir = ic

        ret = pis.states[ir] + corr.states[ir][indx]

        yield (var_name, ret)

class CoefDimM1DimM1( CoefNN ):

    mode2var = {'row' : 0, 'col' : 1}

    def __init__( self, name, problem, kwargs ):
        CoefNN.__init__( self, name, problem, kwargs )
        self.corr_dim -= 1
        
    def get_variables( self, problem, ir, ic, data, mode ):

        nreq = len(self.requires)
        for i in range(nreq-1):
            if i == (nreq-2):
                idx = nreq - 2 + self.mode2var[mode]
                corr = data[self.requires[idx]]
                var_name = self.variables[idx]
                if mode =='col':
                    ir = ic

                out = corr.states[ir].copy()
            else:
                var_name = self.variables[i]
                out = data[self.requires[i]].states.copy()

            yield (var_name, out)

class CoefDimM1( CoefN ):

    def __init__( self, name, problem, kwargs ):
        CoefN.__init__( self, name, problem, kwargs )
        self.corr_dim -= 1
        
    def get_variables( self, problem, ir, data ):

        nreq = len(self.requires) 
        for i in range(nreq):
            if i == (nreq-1):
                out = data[self.requires[i]].states[ir].copy()
            else:
                out = data[self.requires[i]].states.copy()

            yield (self.variables[i], out)

class CoefScalar( MiniAppBase ):

    def get_variables( self, problem, data ):
        for i in range(len(self.requires)):
            yield ( self.variables[i], data[self.requires[i]].states )
        
    def __call__( self, volume, problem = None, data = None ):
        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )

        coef = nm.zeros( (1,), dtype = nm.float64 )

        for name, val in self.get_variables( problem, data ):
            problem.variables[name].data_from_data( val )

        val = eval_term_op( None, self.expression,
                            problem, call_mode = 'd_eval' )

        coef = val / volume

        return coef

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
