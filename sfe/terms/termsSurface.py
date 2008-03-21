from terms import *

##
# c: 13.11.2007, r: 15.01.2008
def fixTractionShape( tr, nEl ):
    tr = nm.array( tr, ndmin = 1 )
    if tr.ndim < 2:
        tr = tr[:,nm.newaxis]
    if tr.shape[0] == 1:
        tr = nm.tile( tr, (nEl,) + tr.shape[1:] )
    return nm.ascontiguousarray( tr )

##
# 22.08.2006, c
class LinearTractionTerm( Term ):
    r""":description: Linear traction forces (weak form).
    :definition: $\int_{\Gamma} \ul{v} \cdot \ull{\sigma} \cdot \ul{n}$, where,
    depending on dimension of 'material' argument, $\ull{\sigma} \cdot
    \ul{n}$ is $\bar{p} \ull{I} \cdot \ul{n}$ for given scalar pressure,
    $\ul{f}$ for traction vector, and itself for a stress tensor
    """
    name = 'dw_surface_ltr'
    argTypes = ('material', 'virtual')
    geometry = [(Surface, 'virtual')]

    # 11.10.2006
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_surface_ltr )
        self.dofConnType = 'surface'

    ##
    # c: 05.09.2006, r: 15.01.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        """
        Works in scalar, vector and tensor mode.
        Tractions defined in vertices -> using 'vertex' subset of leconn
        """
        traction, virtual = self.getArgs( **kwargs )
        ap, sg = virtual.getApproximation( self.getCurrentGroup(), 'Surface' )
        if diffVar is None:
            shape = (chunkSize, 1, sg.dim * sg.nFP, 1)
        else:
            raise StopIteration

        sd = ap.surfaceData[self.region.name]
        bf = ap.getBase( sd.faceType, 0, self.integralName )
        gbf = ap.getBase( sd.faceType, 0, self.integralName,
                          fromGeometry = True )

        traction = fixTractionShape( traction, sd.nodes.shape[0] )
##        sg.str( sys.stdout, 0 )
##         print ap.bf[sd.faceType]
##         pause()

##         if ap.bf[sd.faceType].shape != gbf.shape:
##             raise NotImplementedError, 'tractions on P1 edges only!'
        nEP = gbf.shape[2]
        leconn = sd.leconn[:,:nEP].copy()
        for out, chunk in self.charFun( chunkSize, shape ):
            lchunk = self.charFun.getLocalChunk()
#            print out.shape, lchunk.shape
            status = self.function( out, bf, gbf,
                                    traction, sg, leconn, lchunk )
##             print out
##             print nm.sum( out )
##             pause()
            yield out, lchunk, status
