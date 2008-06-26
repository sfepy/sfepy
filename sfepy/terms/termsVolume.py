from terms import *

##
# 18.09.2006, c
class LinearVolumeForceTerm( Term ):
    r""":description: Vector or scalar linear volume forces (weak form) --- a
    right-hand side source term.
    :definition: $\int_{\Omega} \ul{f} \cdot \ul{v}$ or $\int_{\Omega} f q$
    """
    name = 'dw_volume_lvf'
    argTypes = ('material', 'virtual')
    geometry = [(Volume, 'virtual')]
    useCaches = {'mat_in_qp' : [['material']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_volume_lvf )
        
    ##
    # c: 18.09.2006, r: 07.05.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        force, virtual = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )
        vdim = virtual.field.dim[0]
        
        if diffVar is None:
            shape = (chunkSize, 1, vdim * nEP, 1)
            mode = 0
        else:
            raise StopIteration

        cache = self.getCache( 'mat_in_qp', 0 )
        mat = nm.asarray( force, dtype = nm.float64 )
        if mat.ndim == 1:
            mat = nm.ascontiguousarray( mat[...,nm.newaxis] )
        matQP = cache( 'matqp', self.getCurrentGroup(), 0,
                       mat = mat, ap = ap,
                       assumedShapes = [(nEl, nQP, vdim, 1)],
                       modeIn = None )
#        print matQP
        bf = ap.getBase( 'v', 0, self.integralName )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, bf, matQP, vg, chunk )
            yield out, chunk, status
