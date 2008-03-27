from terms import *
from utils import fixScalarConstant, fixScalarInEl

##
# 28.11.2005, c
class LaplaceTerm( Term ):
    r""":description: Laplace term with $c$ constant or constant per element.
    :definition: $c \int_{\Omega}\nabla s \cdot \nabla r$
    or $\sum_{K \in \Tcal_h}\int_{T_K} c_K\ \nabla s \cdot \nabla r$
    """
    name = 'dw_laplace'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # c: 28.11.2005, r: 27.03.2008
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )

    ##
    # c: 27.03.2008, r: 27.03.2008
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
    # c: 27.03.2008, r: 27.03.2008
    def buildCFunArgs( self, mat, state, ap, vg ):
        vec, indx = state()
        matArg = fixScalarConstant( mat, nm.float64 )
        if matArg is None:
            matArg = fixScalarInEl( mat, self.dataShape[0], nm.float64 )
            self.function = terms.dw_st_pspg_p
        else:
            self.function = terms.term_laplace_asm
        fargs = vec, indx.start, matArg, vg, ap.econn

        return self.function, fargs
        
    ##
    # c: 28.11.2005, r: 27.03.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        material, virtual, state = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )

        shape, mode = self.getShape( diffVar, chunkSize, ap )
        function, fargs = self.buildCFunArgs( material, state, ap, vg )
        
        for out, chunk in self.charFun( chunkSize, shape ):
            status = function( out, *fargs + (chunk, mode) )
            yield out, chunk, status
