from terms import *

##
# 28.11.2005, c
class LaplaceTerm( Term ):
    r""":description: Laplace term (constant parameter).
    :definition: $c \int_{\Omega}\nabla s : \nabla r $
    """
    name = 'dw_laplace'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
        self.function = terms.term_laplace_asm
        
    ##
    # c: 28.11.2005, r: 20.02.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        material, virtual, state = self.getArgs( **kwargs )
##         print self.__class__.name, material.name, virtual.name, state.name

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

##         print nEl, nQP, dim, nEP
##         print material
##         pause()

        vec, indx = state()
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec, indx.start, nm.float64( material ),
                                    vg, ap.econn, chunk, mode )
            yield out, chunk, status
