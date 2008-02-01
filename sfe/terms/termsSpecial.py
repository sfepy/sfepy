from terms import *

##
# c: 01.02.2008, r: 01.02.2008
class MassScalarVariableTerm( Term ):
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
    # c: 01.02.2008, r: 01.02.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        mat, virtual, state = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )
        
        if diffVar is None:
            shape = (chunkSize, 1, nEP, 1)
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            shape = (chunkSize, 1, nEP, nEP)
            mode = 1
        else:
            raise StopIteration

        cache = self.getCache( 'mat_in_qp', 0 )
        matQP = cache( 'matqp', self.getCurrentGroup(), 0,
                       mat = mat, ap = ap,
                       assumedShapes = [(nEl, nQP, 1, 1)],
                       modeIn = 'vertex' )
##         print matQP
##         pause()

        vec, indx = state()
        bf = ap.getBase( 'v', 0, self.integralName )
        for out, chunk in vectorChunkGenerator( nEl, chunkSize, shape ):
            status = self.function( out, matQP, vec, indx.start, bf,
                                    vg, ap.econn, chunk, mode )
            
            yield out, chunk, status
