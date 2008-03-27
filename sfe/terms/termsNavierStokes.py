from terms import *
from utils import fixScalarInEl

##
# 24.10.2005, c
class DivGradTerm( Term ):
    r""":description: Diffusion term.
    :definition: $\int_{\Omega} \nu\ \nabla \ul{v} : \nabla \ul{u}$
    """
    name = 'dw_div_grad'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # 24.10.2005, c
    # 15.11.2005
    # 16.11.2005
    # 10.01.2006
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.term_ns_asmDivGrad )
        
    ##
    # c: 26.10.2005, r: 20.02.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        material, virtual, state = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, dim * nEP, 1 )
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            shape = (chunkSize, 1, dim * nEP, dim * nEP )
            mode = 1
        else:
            raise StopIteration

        vec, indx = state()
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec, indx.start, nm.float64( material ),
                                    vg, ap.econn, chunk, mode )
            yield out, chunk, status

##
# 20.12.2005, c
class ConvectTerm( Term ):
    r""":description: Nonlinear convective term.
    :definition: $\int_{\Omega} ((\ul{u} \cdot \nabla) \ul{u}) \cdot \ul{v}$
    """
    name = 'dw_convect'
    argTypes = ('virtual', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # 20.12.2005, c
    # 10.01.2006
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.term_ns_asmConvect )
        
    ##
    # 20.12.2005, c
    # 25.07.2006
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        virtual, state = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, dim * nEP, 1 )
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            shape = (chunkSize, 1, dim * nEP, dim * nEP )
            mode = 1
        else:
            raise StopIteration

        vec, indx = state()
        bf = ap.getBase( 'v', 0, self.integralName )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec, indx.start, bf,
                                    vg, ap.econn, chunk, mode )
            yield out, chunk, status

##
# 25.07.2007, c
class LinearConvectTerm( Term ):
    r""":description: Linearized convective term.
    :definition: $\int_{\Omega} ((\ul{b} \cdot \nabla) \ul{u}) \cdot \ul{v}$
    """
    name = 'dw_lin_convect'
    argTypes = ('virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # 25.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_lin_convect )
        
    ##
    # 25.07.2007, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        virtual, par, state = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, dim * nEP, 1 )
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            shape = (chunkSize, 1, dim * nEP, dim * nEP )
            mode = 1
        else:
            raise StopIteration

        vec1, i1 = par()
        vec2, i2 = state()
        bf = ap.getBase( 'v', 0, self.integralName )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec1, i1.start, vec2, i2.start,
                                    bf, vg, ap.econn, chunk, mode )
            yield out, chunk, status

##
# 30.07.2007, c
class LinearConvectQTerm( Term ):
    r""":description: Linearized convective term evaluated in quadrature points.
    :definition: $((\ul{b} \cdot \nabla) \ul{u})|_{qp}$
    """
    name = 'dq_lin_convect'
    argTypes = ('parameter', 'state')
    geometry = [(Volume, 'state')]

    ##
    # 30.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_lin_convect )
        
    ##
    # 30.07.2007, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        par, state = self.getArgs( **kwargs )
        ap, vg = state.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, nQP, dim, 1 )
            mode = 2
        else:
            raise StopIteration

        vec1, i1 = par()
        vec2, i2 = state()
        bf = ap.getBase( 'v', 0, self.integralName )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec1, i1.start, vec2, i2.start,
                                    bf, vg, ap.econn, chunk, mode )
            yield out, chunk, status

##
# 15.12.2005, c
class GradTerm( Term ):
    r""":description: Gradient term (weak form).
    :definition: $\int_{\Omega}  p\ \nabla \cdot \ul{v}$
    """
    name = 'dw_grad'
    argTypes = ('virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # 15.12.2005, c
    # 10.01.2006
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.term_ns_asmGrad )
        
    ##
    # 15.12.2005, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        virtual, state = self.getArgs( **kwargs )
        apr, vgr = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        apc, vgc = state.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEPR = apr.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, dim * nEPR, 1 )
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            nEPC = apc.getVDataShape( self.integralName )[3]
            shape = (chunkSize, 1, dim * nEPR, nEPC )
            mode = 1
        else:
            raise StopIteration

        vec, indx = state()
        bf = apc.getBase( 'v', 0, self.integralName )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec, indx.start, bf,
                                    vgr, apc.econn, chunk, mode )
            yield out, chunk, status

##
# 30.07.2007, c
class GradQTerm( Term ):
    r""":description: Gradient term (weak form) in quadrature points.
    :definition: $(\nabla p)|_{qp}$
    """
    name = 'dq_grad'
    argTypes = ('state',)
    geometry = [(Volume, 'state')]

    ##
    # 30.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dq_grad_scalar )

    ##
    # 30.07.2007, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        state = self.getArgs( **kwargs )
        ap, vg = state.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, nQP, dim, 1 )
            mode = 0
        else:
            raise StopIteration

        vec, indx = state()
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec, indx.start,
                                    vg, apc.econn, chunk )
            yield out, chunk, status

##
# 09.03.2007, c
class GradDtTerm( GradTerm ):
    r""":description: Gradient term (weak form) with time-discretized $\dot{p}$.
    :definition: $ \int_{\Omega}  \frac{p - p_0}{\dt} \nabla \cdot \ul{v}$
    :arguments: ts.dt : $\dt$, parameter : $p_0$
    """
    name = 'dw_gradDt'
    argTypes = ('ts', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # 09.03.2007, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        ts, virtual, state, state_0 = self.getArgs( **kwargs )
        apr, vgr = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        apc, vgc = state.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEPR = apr.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, dim * nEPR, 1 )
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            nEPC = apc.getVDataShape( self.integralName )[3]
            shape = (chunkSize, 1, dim * nEPR, nEPC )
            mode = 1
        else:
            raise StopIteration

        vec, indx = state()
        vec0, indx0 = state_0()
        dvec = vec[indx] - vec0[indx0]
        idt = 1.0/ts.dt 
        bf = apc.getBase( 'v', 0, self.integralName )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, dvec, 0, bf,
                                    vgr, apc.econn, chunk, mode )
            out *= idt
            yield out, chunk, status

##
# 14.12.2005, c
class DivTerm( Term ):
    r""":description: Divergence term (weak form).
    :definition: $\int_{\Omega} q\ \nabla \cdot \ul{u}$
    """
    name = 'dw_div'
    argTypes = ('virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # 14.12.2005, c
    # 10.01.2006
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.term_ns_asmDiv )

    ##
    # 14.12.2005, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        virtual, state = self.getArgs( **kwargs )
        apr, vgr = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        apc, vgc = state.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEPR = apr.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, nEPR, 1 )
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            nEPC = apc.getVDataShape( self.integralName )[3]
            shape = (chunkSize, 1, nEPR, dim * nEPC )
            mode = 1
        else:
            raise StopIteration

        vec, indx = state()
        bf = apr.getBase( 'v', 0, self.integralName )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec, indx.start, bf,
                                    vgc, apc.econn, chunk, mode )
            yield out, chunk, status

##
# 27.02.2007, c
class DivRTerm( DivTerm ):
    r""":description: Divergence term (weak form) with a known field (to use
    on a right-hand side).
    :definition: $\int_{\Omega} q\ \nabla \cdot \ul{w}$
    """
    name = 'dw_div_r'
    argTypes = ('virtual', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'parameter')]

##
# 13.03.2007, c
class DivIntegratedTerm( Term ):
    r""":description: Integrated divergence term (weak form).
    :definition: $\int_{\Omega} \bar{p}\ \nabla \cdot \ul{w}$
    """
    name = 'd_div'
    argTypes = ('parameter_1', 'parameter_2')
    geometry = [(Volume, 'parameter_1'), (Volume, 'parameter_2')]
    useCaches = {'state_in_volume_qp' : [['parameter_1']],
                 'div_vector' : [['parameter_2']]}

    ##
    # 13.03.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )

    ##
    # c: 13.03.2007, r: 17.01.2008
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        par1, par2 = self.getArgs( **kwargs )
        apc, vgc = par2.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEPC = apc.getVDataShape( self.integralName )
        shape = (0,)

        cache = self.getCache( 'state_in_volume_qp', 0 )
        vec1 = cache( 'state', self.getCurrentGroup(), 0, state = par1 )
        cache = self.getCache( 'div_vector', 0 )
        div2 = cache( 'div', self.getCurrentGroup(), 0, state = par2 )

        for out, chunk in self.charFun( chunkSize, shape ):
            out = nm.sum( vec1[chunk] * div2[chunk] * vgc.variable( 1 ) )
            yield out, chunk, 0

##
# 26.07.2007, c
class GradDivStabilizationTerm( Term ):
    r""":description: Grad-div stabilization term ($\gamma$ is a global
    stabilization parameter).
    :definition: $\gamma \int_{\Omega} (\nabla\cdot\ul{u}) \cdot
    (\nabla\cdot\ul{v})$
    """
    name = 'dw_st_grad_div'
    argTypes = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # 26.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_grad_div )
        
    ##
    # 26.07.2007, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        gamma, virtual, state = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, dim * nEP, 1 )
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            shape = (chunkSize, 1, dim * nEP, dim * nEP )
            mode = 1
        else:
            raise StopIteration

        vec, indx = state()
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec, indx.start, float( gamma ),
                                    vg, ap.econn, chunk, mode )
            yield out, chunk, status

##
# 31.07.2007, c
from termsLaplace import LaplaceTerm
class PSPGPStabilizationTerm( LaplaceTerm ):
    r""":description: PSPG stabilization term, pressure part ($\tau$ is a local
    stabilization parameter), alias to Laplace term dw_laplace.
    :definition: $\sum_{K \in \Tcal_h}\int_{T_K} \tau_K\ \nabla p \cdot \nabla q$
    """
    name = 'dw_st_pspg_p'

##
# 31.07.2007, c
class PSPGCStabilizationTerm( Term ):
    r""":description: PSPG stabilization term, convective part ($\tau$ is a local
    stabilization parameter).
    :definition: $\sum_{K \in \Tcal_h}\int_{T_K} \tau_K\ ((\ul{b}
    \cdot \nabla) \ul{u}) \cdot \nabla q $
    """
    name = 'dw_st_pspg_c'
    argTypes = ('material', 'virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # 31.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_pspg_c )
        
    ##
    # 31.07.2007, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        tau, virtual, par, state = self.getArgs( **kwargs )
        apr, vgr = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        apc, vgc = state.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEPR = apr.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, nEPR, 1 )
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            nEPC = apc.getVDataShape( self.integralName )[3]
            shape = (chunkSize, 1, nEPR, dim * nEPC )
            mode = 1
        else:
            raise StopIteration

        tauInEl = fixScalarInEl( tau, nEl, nm.float64 )
            
        vec1, i1 = par()
        vec2, i2 = state()
        bf = apc.getBase( 'v', 0, self.integralName )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec1, i1.start, vec2, i2.start,
                                    tauInEl, bf, vgr, vgc,
                                    apc.econn, chunk, mode )
            yield out, chunk, status

##
# 31.07.2007, c
class SUPGPStabilizationTerm( Term ):
    r""":description: SUPG stabilization term, pressure part ($\delta$ is a local
    stabilization parameter).
    :definition: $\sum_{K \in \Tcal_h}\int_{T_K} \delta_K\  \nabla p\cdot
    ((\ul{b} \cdot \nabla) \ul{v})$
    """
    name = 'dw_st_supg_p'
    argTypes = ('material', 'virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # 31.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_supg_p )
        
    ##
    # 31.07.2007, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        delta, virtual, par, state = self.getArgs( **kwargs )
        apr, vgr = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        apc, vgc = state.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEPR = apr.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, dim * nEPR, 1 )
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            nEPC = apc.getVDataShape( self.integralName )[3]
            shape = (chunkSize, 1, dim * nEPR, nEPC )
            mode = 1
        else:
            raise StopIteration

        deltaInEl = fixScalarInEl( delta, nEl, nm.float64 )
            
        vec1, i1 = par()
        vec2, i2 = state()
        bf = apr.getBase( 'v', 0, self.integralName )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec1, i1.start, vec2, i2.start,
                                    deltaInEl, bf, vgr, vgc,
                                    apr.econn, apc.econn, chunk, mode )
            yield out, chunk, status

##
# 31.07.2007, c
class SUPGCStabilizationTerm( Term ):
    r""":description: SUPG stabilization term, convective part ($\delta$ is a
    local stabilization parameter).
    :definition: $\sum_{K \in \Tcal_h}\int_{T_K} \delta_K\ ((\ul{b}
    \cdot \nabla) \ul{u})\cdot ((\ul{b} \cdot \nabla) \ul{v})$
    """
    name = 'dw_st_supg_c'
    argTypes = ('material', 'virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # 31.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_supg_c )
        
    ##
    # 31.07.2007, c
    def __call__( self, diffVar = None, chunkSize = None, **kwargs ):
        delta, virtual, par, state = self.getArgs( **kwargs )
        ap, vg = virtual.getApproximation( self.getCurrentGroup(), 'Volume' )
        nEl, nQP, dim, nEP = ap.getVDataShape( self.integralName )

        if diffVar is None:
            shape = (chunkSize, 1, dim * nEP, 1 )
            mode = 0
        elif diffVar == self.getArgName( 'state' ):
            shape = (chunkSize, 1, dim * nEP, dim * nEP )
            mode = 1
        else:
            raise StopIteration

        deltaInEl = fixScalarInEl( delta, nEl, nm.float64 )
            
        vec1, i1 = par()
        vec2, i2 = state()
        bf = ap.getBase( 'v', 0, self.integralName )
        for out, chunk in self.charFun( chunkSize, shape ):
            status = self.function( out, vec1, i1.start, vec2, i2.start,
                                    deltaInEl, bf, vg,
                                    ap.econn, chunk, mode )
            yield out, chunk, status
