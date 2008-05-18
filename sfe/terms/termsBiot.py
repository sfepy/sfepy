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
# c: 18.05.2008
class BiotGradDtTerm( GradTerm ):
    r""":description: Biot gradient-like term (weak form) with time-discretized
    $\dot{p}$ and $\alpha_{ij}$
    given in vector form exploiting symmetry: in 3D it has the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it has
    the indices ordered as $[11, 22, 12]$.
    :definition: $\int_{\Omega}  \frac{p - p_0}{\dt}\ \alpha_{ij} e_{ij}(\ul{v})$
    :arguments: ts.dt : $\dt$, parameter : $p_0$
    """
    name = 'dw_biot_grad_dt'
    argTypes = ('ts', 'material', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    useCaches = {'state_in_volume_qp' : [['state', {'state' : (2, 2)}]]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_biot_grad )

    ##
    # c: 18.05.2008, r: 18.05.2008
    def buildCFunArgs( self, state, apc, vgr, **kwargs ):
        cache = self.getCache( 'state_in_volume_qp', 0 )
        vec_qp = cache( 'state', self.getCurrentGroup(), 0, state = state )
        vec0_qp = cache( 'state', self.getCurrentGroup(), 1,
                         state = parameter )

        ts, mat = self.getArgs( ['ts', 'material'], **kwargs )
        matQP = mat[nm.newaxis,:,nm.newaxis].repeat( self.dataShape[1], 0 )
#        print matQP
        bf = apc.getBase( 'v', 0, self.integralName )
        return 1.0 / ts.dt, vec_qp - vec0_qp, bf, matQP, vgr

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
# c: 03.04.2008
class BiotGradTHTerm( BiotGradTerm ):
    r""":definition: $\int_{\Omega} \left [\int_0^t
    \alpha_{ij}(t-\tau)\,p(\tau)) \difd{\tau} \right]\,e_{ij}(\ul{v})$
    """
    name = 'dw_biot_grad_th'
    argTypes = ('ts', 'material', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    useCaches = {'state_in_volume_qp' : [['state', {'state' : (-1,-1)}]]}

    ##
    # c: 03.04.2008, r: 04.04.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        """history for now is just state_0, it is not used anyway, as the
        history is held in the dstrain cache"""
        ts, mats, virtual, state, history = self.getArgs( **kwargs )
        apr, vgr = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        apc, vgc = state.getApproximation( self.getCurrentGroup(), 'Volume' )

        shape, mode = self.getShape( diffVar, chunkSize, apr, apc )
        nEl, nQP, dim, nEP = self.dataShape

        if (ts.step == 0) and (mode == 0):
            raise StopIteration

        bf = apc.getBase( 'v', 0, self.integralName )

        if mode == 1:
            matQP = mats[0][nm.newaxis,:,nm.newaxis].repeat( nQP, 0 )
            for out, chunk in self.charFun( chunkSize, shape ):
                status = self.function( out, 1.0, nm.empty( 0 ), bf,
                                        matQP, vgr, chunk, 1 )
                yield out, chunk, status
        else:
            cache = self.getCache( 'state_in_volume_qp', 0 )
            for out, chunk in self.charFun( chunkSize, shape, zero = True ):
                out1 = nm.empty_like( out )
                for ii, mat in enumerate( mats ):
                    matQP = mat[nm.newaxis,:,nm.newaxis].repeat( nQP, 0 )
                    vec_qp = cache( 'state', self.getCurrentGroup(), ii,
                                    state = state, history = history )
                    status = self.function( out1, 1.0, vec_qp, bf,
                                            matQP, vgr, chunk, 0 )
                    out += out1
                yield out, chunk, status

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

##
# c: 03.04.2008
class BiotDivTHTerm( BiotDivTerm ):
    r""":definition: $\int_{\Omega} \left [\int_0^t
    \alpha_{ij}(t-\tau) \tdiff{e_{kl}(\ul{u}(\tau))}{\tau}
    \difd{\tau} \right] q$
    """
    name = 'dw_biot_div_th'
    argTypes = ('ts', 'material', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    useCaches = {'cauchy_strain' : [['state', {'strain' : (2,2),
                                               'dstrain' : (-1,-1)}]]}

    ##
    # c: 03.04.2008, r: 04.04.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        """history for now is just state_0, it is not used anyway, as the
        history is held in the dstrain cache"""
        ts, mats, virtual, state, history = self.getArgs( **kwargs )
        apr, vgr = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        apc, vgc = state.getApproximation( self.getCurrentGroup(), 'Volume' )

        shape, mode = self.getShape( diffVar, chunkSize, apr, apc )
        nEl, nQP, dim, nEP = self.dataShape

        if (ts.step == 0) and (mode == 0):
            raise StopIteration

        bf = apr.getBase( 'v', 0, self.integralName )

        if mode == 1:
            matQP = mats[0][nm.newaxis,:,nm.newaxis].repeat( nQP, 0 )
            for out, chunk in self.charFun( chunkSize, shape ):
                status = self.function( out, 1.0 / ts.dt, nm.empty( 0 ), bf,
                                        matQP, vgc, chunk, 1 )
                yield out, chunk, status
        else:
            cache = self.getCache( 'cauchy_strain', 0 )
            for out, chunk in self.charFun( chunkSize, shape, zero = True ):
                out1 = nm.empty_like( out )
                for ii, mat in enumerate( mats ):
                    matQP = mat[nm.newaxis,:,nm.newaxis].repeat( nQP, 0 )
                    dstrain = cache( 'dstrain', self.getCurrentGroup(), ii,
                                     state = state, history = history )
                    status = self.function( out1, 1.0 / ts.dt, dstrain,
                                            bf, matQP, vgc, chunk, 0 )
                    out += out1
                yield out, chunk, status
