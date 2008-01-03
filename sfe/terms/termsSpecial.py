from terms import *

class VTerm( Term ):
    name = 'dq_mass_v'
    argTypes = ('material', 'parameter')
    geometry = [(Volume, 'parameter')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        """chunkSize = number of elements to process in one call"""
        material, psi = self.getArgs( **kwargs )
        ap, vg = psi.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, nQP, dim, dim)
        else:
            raise StopIteration

        # vec[indx] is the state vector.
        vec, indx = psi()

        # Get volume base function
        bf = ap.getBase( 'v', 0, self.integralName )
        # Get volume base function gradient
        bfg = ap.getBase( 'v', 1, self.integralName )

        # Element connectivity
        conn = ap.econn
        for out, chunk in self.charFun( chunkSize, shape ):
            # Fill-in 'out', its shape is 'shape'
#            out[...] = ...

            # 'chunk' are the indices of elements (into ap.econn)
            # 'status' == 0 means success
            status = 0
            yield out, chunk, status
