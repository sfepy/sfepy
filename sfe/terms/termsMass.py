from terms import *

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

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_mass )
        
    ##
    # 21.11.2006, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        ts, mat, virtual, state, state0 = self.getArgs( **kwargs )
        ap, vg = virtual.getCurrentApproximation()
        nEl, nQP, dim, nEP = ap.getVDataShape()
        
        if diffVar is None:
            shape = (chunkSize, 1, dim * nEP, 1)
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            shape = (chunkSize, 1, dim * nEP, dim * nEP)
            mode = 1
        else:
            raise StopIteration

        vec, indx = state()
        vec0, indx0 = state0()
        dvec = vec[indx] - vec0[indx0]
        rhodt = mat / ts.dt

        for out, chunk in vectorChunkGenerator( nEl, chunkSize, shape ):
            status = self.function( out, rhodt, dvec, 0, ap.bf['v'],
                                    vg, ap.econn, chunk, mode )
##             from pylab import spy, show
##             if mode == 1:
##                 print out[0].squeeze()
##                 spy( nm.abs( out[0].squeeze() ) )
##                 show()
            
            yield out, chunk, status

##
# created:       25.09.2007
class MassVectorTerm( MassTerm ):
    r""":description: Vector field mass matrix/rezidual.
    :definition: $\int_{\Omega} \rho\ \ul{v} \cdot \ul{u}$
    """
    name = 'dw_mass_vector'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # created:       25.09.2007
    # last revision: 11.12.2007
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        mat, virtual, state = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )
#        debug()
        if diffVar is None:
            shape = (chunkSize, 1, dim * nEP, 1)
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            shape = (chunkSize, 1, dim * nEP, dim * nEP)
            mode = 1
        else:
            raise StopIteration

        vec, indx = state()
        bf = ap.getBase( 'v', 0, self.integralName )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, float( mat ), vec, indx.start,
                                    bf, vg, ap.econn, chunk, mode )
            
            yield out, chunk, status

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
    # 04.09.2007, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        virtual, state = self.getArgs( **kwargs )
        ap, vg = virtual.getCurrentApproximation()
        nEl, nQP, dim, nEP = ap.getVDataShape()
        
        if diffVar is None:
            shape = (chunkSize, 1, nEP, 1)
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            shape = (chunkSize, 1, nEP, nEP)
            mode = 1
        else:
            raise StopIteration

        vec, indx = state()

        for out, chunk in vectorChunkGenerator( nEl, chunkSize, shape ):
            status = self.function( out, vec, indx.start, ap.bf['v'],
                                    vg, ap.econn, chunk, mode )
            
            yield out, chunk, status

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
    # 05.09.2007, c
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
        for out, chunk in vectorChunkGenerator( nEl, chunkSize, shape ):
            status = self.function( out, vec, indx.start, apr.bf['v'], cbfs,
                                    vgr, apc.econn, iemap, chunk, mode )
            
            yield out, chunk, status
