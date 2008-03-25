from terms import *

##
# 07.03.2006, c
class LinearElasticIsotropicTerm( Term ):
    r""":description: isotropic linear elasticity term.
    :definition: $\int_{\Omega}  D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})$
    with $D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
    \lambda \ \delta_{ij} \delta_{kl}$ 
    """
    name = 'dw_lin_elastic_iso'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # c: 07.03.2006, r: 21.03.2008
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_lin_elastic_iso )
        
    ##
    # c: 21.03.2008, r: 25.03.2008
    def getShape( self, diffVar, chunkSize, apr, apc = None ):
        self.dataShape = apr.getVDataShape( self.integralName )
        nEl, nQP, dim, nEP = self.dataShape
        
        if diffVar is None:
            return (chunkSize, 1, dim * nEP, 1), 0
        elif diffVar == self.getArgName( 'state' ):
            return (chunkSize, 1, dim * nEP, dim * nEP), 1
        else:
            raise StopIteration

    ##
    # c: 21.03.2008, r: 21.03.2008
    def buildCFunArgs( self, mat, state, ap, vg ):
        vec, indx = state()
        lam, mu = map( nm.float64, [mat[ii] for ii in ['lambda', 'mu']] )
        return vec, indx.start, lam, mu, vg, ap.econn

    ##
    # c: 07.03.2006, r: 21.03.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        material, virtual, state = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )

        shape, mode = self.getShape( diffVar, chunkSize, ap )
        fargs = self.buildCFunArgs( material, state, ap, vg )
        
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, *fargs + (chunk, mode) )
            yield out, chunk, status

##
# c: 02.08.2006
class LinearElasticTerm( LinearElasticIsotropicTerm ):
    r""":description: general linear elasticity term, with $D_{ijkl}$ given in
    the usual matrix form exploiting symmetry: in 3D it is $6\times6$ with the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it is $3\times3$ with
    the indices ordered as $[11, 22, 12]$.
    :definition: $\int_{\Omega}  D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})$
     
    """
    name = 'dw_lin_elastic'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]
    useCaches = {'cauchy_strain' : [['state', {'strain' : (1,1)}]]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_lin_elastic )

    ##
    # c: 25.03.2008, r: 25.03.2008
    def buildCFunArgs( self, mat, state, ap, vg ):
        matQP = mat[nm.newaxis,:,:].repeat( self.dataShape[1], 0 )
        cache = self.getCache( 'cauchy_strain', 0 )
        strain = cache( 'strain', self.getCurrentGroup(), 0, state = state )
        return 1.0, strain, matQP, vg

##
# 21.09.2006, c
class CauchyStrainTerm( Term ):
    r""":description: Cauchy strain tensor averaged in elements.
    :definition: vector of $\forall K \in \Tcal_h: \int_{T_K} \ull{e}(\ul{w})$
    """
    name = 'de_cauchy_strain'
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
