from terms import *

##
# 10.07.2007, c
class LinearPointSpringTerm( Term ):
    r""":description: Linear springs constraining movement of FE nodes in a
    reagion; use as a relaxed Dirichlet boundary conditions.
    :definition: $\ul{f}^i = -k \ul{u}^i \quad \forall \mbox{ FE node } i
    \mbox{ in region }$
    """
    name = 'dw_point_lspring'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Point, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        self.dofConnType = 'point'
        
    ##
    # 10.07.2007, c
    # 18.07.2007
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        """TODO: projection to direction"""
        mat, virtual, state = self.getArgs( **kwargs )
        ap, pg = virtual.getApproximation( self.getCurrentGroup(), 'Point' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )

        if self.charFun.iCurrent > 0:
            raise StopIteration
            
        vecU = state.getStateInRegion( self.region )
        nNod = vecU.shape[0]
##         print vecU.shape
##         pause()

        if diffVar is None:
            shape = (chunkSize, 1, dim, 1)
            for out, chunk in vectorChunkGenerator( nNod, chunkSize, shape ):
                out[:,0,:,0] = - mat.stiffness * vecU[chunk]
                yield out, chunk, 0

        elif diffVar == self.getArgName( 'state' ):
            shape = (chunkSize, 1, dim, dim)
            eye = nm.eye( dim, dim, dtype = nm.float64 )
            eye.shape = (1, 1) + eye.shape
            for out, chunk in vectorChunkGenerator( nNod, chunkSize, shape ):
                out[...] = - mat.stiffness * eye
                yield out, chunk, 0

        else:
            raise StopIteration
