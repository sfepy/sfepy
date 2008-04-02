from terms import *
from termsNavierStokes import GradTerm, DivTerm

##
# 01.08.2006, c
class BiotGradTerm( GradTerm ):
    r""":description: Biot gradient-like term (weak form) with $\alpha_{ij}$
    given in vector form exploiting symmetry: in 3D it has the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it has
    the indices ordered as $[11, 22, 12]$.
    :definition: $\int_{\Omega}  p\ \alpha_{ij} e_{ij}(\ul{v})$
    """
    name = 'dw_biot_grad'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    useCaches = {'state_in_volume_qp' : [['state']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_biot_grad )

    ##
    # c: 31.03.2008, r: 31.03.2008
    def buildCFunArgs( self, state, apc, vgr, **kwargs ):
        cache = self.getCache( 'state_in_volume_qp', 0 )
        vec_qp = cache( 'state', self.getCurrentGroup(), 0, state = state )

        mat, = self.getArgs( ['material'], **kwargs )
        matQP = mat[nm.newaxis,:,nm.newaxis].repeat( self.dataShape[1], 0 )
#        print matQP
        bf = apc.getBase( 'v', 0, self.integralName )
        return 1.0, vec_qp, bf, matQP, vgr

##
# 20.09.2006, c
class BiotGradRTerm( BiotGradTerm ):
    r""":description: Biot gradient-like term (weak form) with $\alpha_{ij}$
    given in vector form exploiting symmetry: in 3D it has the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it has
    the indices ordered as $[11, 22, 12]$. The argument $r$ is a known field
    (to use on a right-hand side).
    :definition: $\int_{\Omega}  r\ \alpha_{ij} e_{ij}(\ul{v})$
    """
    name = 'dw_biot_grad_r'
    argTypes = ('material', 'virtual', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'parameter')]
    useCaches = {'state_in_volume_qp' : [['parameter']]}

##
# 01.08.2006, c
class BiotDivDtTerm( DivTerm ):
    r""":description: Biot divergence-like rate term (weak form) with
    $\alpha_{ij}$ given in vector form exploiting symmetry: in 3D it has the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it has
    the indices ordered as $[11, 22, 12]$.
    :definition: $\int_{\Omega} q\ \alpha_{ij} \frac{e_{ij}(\ul{u}) -
    e_{ij}(\ul{u_0})}{\dt}$
    """
    name = 'dw_biot_div_dt'
    argTypes = ('ts', 'material', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    useCaches = {'cauchy_strain' : [['state', {'strain' : (2,2),
                                               'dstrain' : (1,1)}]]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_biot_div )

    ##
    # c: 31.03.2008, r: 02.04.2008
    def buildCFunArgs( self, state, apr, apc, vgc, **kwargs ):
        cache = self.getCache( 'cauchy_strain', 0 )
        dstrain = cache( 'dstrain', self.getCurrentGroup(), 0, state = state )

        ts, mat = self.getArgs( ['ts', 'material'], **kwargs )
        matQP = mat[nm.newaxis,:,nm.newaxis].repeat( self.dataShape[1], 0 )
#        print matQP
        bf = apr.getBase( 'v', 0, self.integralName )
        return 1.0 / ts.dt, dstrain, bf, matQP, vgc


##
# 20.09.2006, c
class BiotDivTerm( DivTerm ):
    r""":description: Biot divergence-like term (weak form) with
    $\alpha_{ij}$ given in vector form exploiting symmetry: in 3D it has the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it has
    the indices ordered as $[11, 22, 12]$.
    :definition: $\int_{\Omega} q\ \alpha_{ij} e_{ij}(\ul{u})$
    """
    name = 'dw_biot_div'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    useCaches = {'cauchy_strain' : [['state']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_biot_div )

    ##
    # c: 31.03.2008, r: 31.03.2008
    def buildCFunArgs( self, state, apr, apc, vgc, **kwargs ):
        cache = self.getCache( 'cauchy_strain', 0 )
        strain = cache( 'strain', self.getCurrentGroup(), 0, state = state )
        mat, = self.getArgs( ['material'], **kwargs )
        matQP = mat[nm.newaxis,:,nm.newaxis].repeat( self.dataShape[1], 0 )
        bf = apr.getBase( 'v', 0, self.integralName )
        return 1.0, strain, bf, matQP, vgc

##
# c: 05.03.2008
class BiotDivRTerm( BiotDivTerm ):
    r""":description: Biot divergence-like term (weak form) with
    $\alpha_{ij}$ given in vector form exploiting symmetry: in 3D it has the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it has
    the indices ordered as $[11, 22, 12]$. The argument $\ul{w}$ is a known field
    (to use on a right-hand side).
    :definition: $\int_{\Omega} q\ \alpha_{ij} e_{ij}(\ul{w})$
    """
    name = 'dw_biot_div_r'
    argTypes = ('material', 'virtual', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'parameter')]
    useCaches = {'cauchy_strain' : [['parameter']]}

##
# c: 05.03.2008
class BiotDivRIntegratedTerm( Term ):
    r""":description: Integrated Biot divergence-like term (weak form) with
    $\alpha_{ij}$ given in vector form exploiting symmetry: in 3D it has the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it has
    the indices ordered as $[11, 22, 12]$.
    :definition: $\int_{\Omega} r\ \alpha_{ij} e_{ij}(\ul{w})$
    """
    name = 'd_biot_div'
    argTypes = ('material', 'parameter_1', 'parameter_2')
    geometry = [(Volume, 'parameter_1'), (Volume, 'parameter_2')]
    useCaches = {'state_in_volume_qp' : [['parameter_1']],
                 'cauchy_strain' : [['parameter_2']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_biot_div )

    ##
    # c: 05.03.2008, r: 05.03.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        mat, parP, parU = self.getArgs( **kwargs )
        apr, vgr = parP.getApproximation( self.getCurrentGroup(), 'Volume' )
        apc, vgc = parU.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEPR = apr.getVDataShape( self.integralName )
        
        shape = (chunkSize, 1, 1, 1)
        
        cache = self.getCache( 'state_in_volume_qp', 0 )
        vecQP = cache( 'state', self.getCurrentGroup(), 0, state = parP )

        cache = self.getCache( 'cauchy_strain', 0 )
        strain = cache( 'strain', self.getCurrentGroup(), 0, state = parU )

        matQP = mat[nm.newaxis,:,nm.newaxis].repeat( nQP, 0 )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, 1.0, vecQP, strain, matQP, vgc, chunk )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status
