from terms import *

##
# 18.09.2006, c
class LinearVolumeForceTerm( Term ):
    r""":description: Linear volume forces (weak form).
    :definition: $ \int_{\Omega} \ul{v} \cdot \ul{f}$
    """
    name = 'dw_volume_lvf'
    argTypes = ('material', 'virtual')
    geometry = [(Volume, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
        self.function = terms.dw_volume_lvf
        
    ##
    # 18.09.2006, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        """Full domain volume only! (no ap.l(e)conn...)"""
        force, virtual = self.getArgs( **kwargs )
        ap, vg = virtual.getCurrentApproximation()
        nEl, nQP, dim, nEP = ap.getVDataShape()
        
        if diffVar is None:
            shape = (chunkSize, 1, dim * nEP, 1)
            mode = 0
        else:
            raise StopIteration

        gbf = ap.gbf['v']
        for out, chunk in vectorChunkGenerator( nEl, chunkSize, shape ):
            status = self.function( out, ap.bf['v'], gbf,
                                    force, vg, ap.sub.conn, chunk )
            yield out, chunk, status
