from terms import *
from utils import chooseScalarOrInEl

##
# 21.11.2006, c
class MassTerm( Term ):
    r""":description: Inertial forces term (constant density).
    :definition: $\int_{\Omega} \rho \ul{v} \cdot \frac{\ul{u} -
    \ul{u}_0}{\dt}$
    :arguments: material.rho : $\rho$, ts.dt : $\dt$, parameter : $\ul{u}_0$"""
    name = 'dw_mass'
    argTypes = ('ts', 'material', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual')]

    ##
    # created:       21.11.2006
    # last revision: 19.12.2007
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )

    ##
    # c: 01.04.2008, r: 01.04.2008
    def getShape( self, diffVar, chunkSize, apr, apc = None ):
        self.dataShape = apr.getVDataShape( self.integralName )
        nEl, nQP, dim, nEP = self.dataShape

        if diffVar is None:
            return (chunkSize, 1, dim * nEP, 1), 0
        elif diffVar == self.getArgName( 'state' ):
            return (chunkSize, 1, dim * nEP, dim * nEP), 1
        else:
            raise StopIteration
        
    ##
    # c: 01.04.2008, r: 01.04.2008
    def buildCFunArgs( self, mat, state, ap, vg ):
        # terms.dw_mass_rho_in_el is missing
        mat, self.function = chooseScalarOrInEl( mat, nm.float64,
                                                 terms.dw_mass,
                                                 NotImplemented )
        ts, state0 = self.getArgs( ['ts', 'parameter'], **kwargs )

        vec, indx = state()
        vec0, indx0 = state0()
        dvec = vec[indx] - vec0[indx0]
        rhodt = mat / ts.dt
        bf = ap.getBase( 'v', 0, self.integralName )
        return rhodt, dvec, 0, bf, vg, ap.econn

    ##
    # c: 21.11.2006, r: 01.04.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        mat, virtual, state = self.getArgs( ['material', 'virtual', 'state'],
                                            **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )

        shape, mode = self.getShape( diffVar, chunkSize, ap )
        fargs = self.buildCFunArgs( mat, state, ap, vg )

        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, *fargs + (chunk, mode) )
            yield out, chunk, status

##
# c: 25.09.2007
class MassVectorTerm( MassTerm ):
    r""":description: Vector field mass matrix/rezidual.
    :definition: $\int_{\Omega} \rho\ \ul{v} \cdot \ul{u}$
    """
    name = 'dw_mass_vector'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # c: 01.04.2008, r: 01.04.2008
    def buildCFunArgs( self, mat, state, ap, vg ):
        mat, self.function = chooseScalarOrInEl( mat, nm.float64,
                                                 terms.dw_mass,
                                                 NotImplemented )
        vec, indx = state()
        bf = ap.getBase( 'v', 0, self.integralName )
        return mat, vec, indx.start, bf, vg, ap.econn

##
# 04.09.2007, c
class MassScalarTerm( Term ):
    r""":description: Scalar field mass matrix/rezidual.
    :definition: $\int_{\Omega} q p$
    """
    name = 'dw_mass_scalar'
    argTypes = ('virtual', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # 04.09.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_mass_scalar )
        
    ##
    # c: 01.04.2008, r: 01.04.2008
    def getShape( self, diffVar, chunkSize, apr, apc = None ):
        self.dataShape = apr.getVDataShape( self.integralName )
        nEl, nQP, dim, nEP = self.dataShape

        if diffVar is None:
            return (chunkSize, 1, nEP, 1), 0
        elif diffVar == self.getArgName( 'state' ):
            return (chunkSize, 1, nEP, nEP), 1
        else:
            raise StopIteration
        
    ##
    # c: 01.04.2008, r: 01.04.2008
    def buildCFunArgs( self, state, ap, vg, **kwargs ):
        vec, indx = state()
        bf = ap.getBase( 'v', 0, self.integralName )
        return vec, indx.start, bf, vg, ap.econn

    ##
    # c: 04.09.2007, r: 01.04.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        if self.name.endswith( '_r' ):
            virtual, state = self.getArgs( ['virtual', 'parameter'], **kwargs )
        else:
            virtual, state = self.getArgs( ['virtual', 'state'], **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )

        shape, mode = self.getShape( diffVar, chunkSize, ap )
        fargs = self.buildCFunArgs( state, ap, vg, **kwargs )

        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, *fargs + (chunk, mode) )
            yield out, chunk, status

##
# 06.02.2008, c
class MassScalarRTerm( MassScalarTerm ):
    r""":description: Scalar field mass rezidual --- $r$ is assumed to be known.
    :definition: $\int_{\Omega} q r$
    """
    name = 'dw_mass_scalar_r'
    argTypes = ('virtual', 'parameter')
    geometry = [(Volume, 'virtual')]

##
# c: 01.02.2008
class MassScalarVariableTerm( MassScalarTerm ):
    r""":description: Scalar field mass matrix/rezidual with coefficient $c$
    defined in nodes.
    :definition: $\int_{\Omega} c q p$
    """
    name = 'dw_mass_scalar_variable'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]
    useCaches = {'mat_in_qp' : [['material']]}

    ##
    # c: 01.02.2008, r: 01.02.2008
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_mass_scalar_variable )
        
    ##
    # c: 01.04.2008, r: 01.04.2008
    def buildCFunArgs( self, state, ap, vg, **kwargs ):
        nEl, nQP = self.dataShape[:2]
        
        mat, = self.getArgs( ['material'], **kwargs )
        cache = self.getCache( 'mat_in_qp', 0 )
        matQP = cache( 'matqp', self.getCurrentGroup(), 0,
                       mat = mat, ap = ap,
                       assumedShapes = [(nEl, nQP, 1, 1)],
                       modeIn = 'vertex' )

        vec, indx = state()
        bf = ap.getBase( 'v', 0, self.integralName )
        return matQP, vec, indx.start, bf, vg, ap.econn

##
# 05.09.2007, c
class MassScalarFineCoarseTerm( Term ):
    r""":description: Scalar field mass matrix/rezidual for coarse to fine grid
    interpolation. Field $p_H$ belong to the coarse grid, test field $q_h$ to
    the fine grid.
    :definition: $\int_{\Omega} q_h p_H$
    """
    name = 'dw_mass_scalar_fine_coarse'
    argTypes = ('virtual', 'state', 'iemaps', 'pbase' )
    geometry = [(Volume, 'virtual')]

    ##
    # 05.09.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign,
                       terms.dw_mass_scalar_fine_coarse )
        
    ##
    # c: 05.09.2007, r: 01.04.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        virtual, state, iemaps, pbase = self.getArgs( **kwargs )
        apr, vgr = virtual.getCurrentApproximation()
        apc, vgc = virtual.getCurrentApproximation()
        nEl, nQP, dim, nEPR = apr.getVDataShape()
        
        if diffVar is None:
            shape = (chunkSize, 1, nEPR, 1)
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            nEPC = apc.getVDataShape()[3]
            shape = (chunkSize, 1, nEPR, nEPC)
            mode = 1
        else:
            raise StopIteration

        vec, indx = state()

        cbfs = pbase[self.charFun.ig]
        iemap = iemaps[self.charFun.ig]
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec, indx.start, apr.bf['v'], cbfs,
                                    vgr, apc.econn, iemap, chunk, mode )
            
            yield out, chunk, status
