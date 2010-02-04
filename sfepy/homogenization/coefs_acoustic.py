from sfepy.base.base import *
from sfepy.fem import eval_term_op
from sfepy.homogenization.coefs_base import VolumeFractions, \
     CorrMiniApp, MiniAppBase, CorrN, CoefN, CoefNN, \
     create_scalar_pis

class ShapeDimM1( CorrMiniApp ):
    
    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        pis = create_scalar_pis( problem, self.variables[0] )

        return pis[0:-1].copy()

class CorrDimM1( CorrN ):

    def __init__( self, name, problem, kwargs ):
        CorrN.__init__( self, name, problem, kwargs )
        self.corr_dim -= 1
        
    def get_variables( self, ir, data ):
        """data: pis"""
        pis = data[self.requires[0]]

        yield (self.variables[0], pis[ir])

class CorrScalar( CorrMiniApp ):

    def get_variables( self, problem, i, data ):
        yield ( self.variables[i], data[self.requires[i]].states )

    def __call__( self, problem = None, data = None ):
        problem = get_default( problem, self.problem )

        problem.select_variables( self.variables )
        problem.set_equations( self.equations )

        problem.select_bcs( ebc_names = self.ebcs, epbc_names = self.epbcs )

        self.init_solvers(problem)
        states = nm.zeros( (1,), dtype = nm.object )

        for ival in range( len( self.requires ) ):
            for name, val in self.get_variables( problem, ival, data ):
                problem.variables[name].data_from_data( val )
            
        state = problem.create_state_vector()
        problem.apply_ebc( state )
        state = problem.solve()
        assert_( problem.variables.has_ebc( state ) )

        self.save( state, problem,
                   self.get_save_name(), self.get_dump_name() )

        return Struct( name = self.name,
                       states = state,
                       di = problem.variables.di )

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

        ret = pis[ir] + corr.states[ir][indx]

        yield (var_name, ret)

class CoefDimM1DimM1( CoefNN ):

    mode2var = {'row' : 0, 'col' : 1}

    def __init__( self, name, problem, kwargs ):
        CoefNN.__init__( self, name, problem, kwargs )
        self.corr_dim -= 1
        
    def get_variables( self, problem, ir, ic, data, mode ):

        corr = data[self.requires[-1]]

        var_name = self.variables[self.mode2var[mode]]
        c_name = problem.variables[var_name].primary_var_name
        indx = corr.di.indx[c_name]

        if mode =='col':
            ir = ic

        yield (var_name, corr.states[ir].copy())

class CoefDimM1( CoefN ):

    def __init__( self, name, problem, kwargs ):
        CoefN.__init__( self, name, problem, kwargs )
        self.corr_dim -= 1
        
    def get_variables( self, problem, ir, data ):

        corr = data[self.requires[-1]]
        yield (self.variables[0], corr.states[ir])

class CoefScalar( MiniAppBase ):

    def get_variables( self, problem, data ):
        yield ( self.variables[0], data[self.requires[0]].states )
        
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
