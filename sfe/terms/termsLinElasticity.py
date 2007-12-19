from terms import *

##
# 07.03.2006, c
class SDCCTerm( Term ):
    r""":description: Homogeneous isotropic linear elasticity term.
    :definition: $\int_{\Omega}  D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})$
    with $D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
    \lambda \ \delta_{ij} \delta_{kl}$ 
    """
    name = 'dw_sdcc'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # 07.03.2006, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
        self.function = terms.dw_lin_elasticity
        
    ##
    # 07.03.2006, c
    # 26.07.2006
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        material, virtual, state = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )
        
        if diffVar is None:
            shape = (chunkSize, 1, dim * nEP, 1)
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            shape = (chunkSize, 1, dim * nEP, dim * nEP)
            mode = 1
        else:
            raise StopIteration
        
        vec, indx = state()
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec, indx.start,
                                    material.lam, material.mu,
                                    vg, ap.econn, chunk, mode )
            yield out, chunk, status

##
# 21.09.2006, c
class SDCCStrainTerm( Term ):
    r""":description: Cauchy strain tensor averaged in elements.
    :definition: vector of $\forall K \in \Tcal_h: \int_{T_K} \ull{e}(\ul{w})$
    """
    name = 'de_sdcc_strain'
    argTypes = ('parameter',)
    geometry = [(Volume, 'parameter')]

    ##
    # 05.10.2007
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )

        self.function = terms.de_cauchy_strain
        
    ##
    # 21.09.2006, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        parameter, = self.getArgs( **kwargs )

        ap, vg = parameter.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, dim * (dim + 1) / 2, 1)
        else:
            raise StopIteration

        vec, indx = parameter()
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec, indx.start,
                                    vg, ap.econn, chunk )
            yield out, chunk, status
