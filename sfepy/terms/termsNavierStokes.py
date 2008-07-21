from terms import *
from utils import fix_scalar_in_el

##
# 24.10.2005, c
class DivGradTerm( Term ):
    r""":description: Diffusion term.
    :definition: $\int_{\Omega} \nu\ \nabla \ul{v} : \nabla \ul{u}$
    """
    name = 'dw_div_grad'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # 24.10.2005, c
    # 15.11.2005
    # 16.11.2005
    # 10.01.2006
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.term_ns_asm_div_grad )
        
    ##
    # c: 26.10.2005, r: 20.02.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        material, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_ep, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            shape = (chunk_size, 1, dim * n_ep, dim * n_ep )
            mode = 1
        else:
            raise StopIteration

        vec, indx = state()
        for out, chunk in self.char_fun( chunk_size, shape ):
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
    arg_types = ('virtual', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # 20.12.2005, c
    # 10.01.2006
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.term_ns_asm_convect )
        
    ##
    # 20.12.2005, c
    # 25.07.2006
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_ep, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            shape = (chunk_size, 1, dim * n_ep, dim * n_ep )
            mode = 1
        else:
            raise StopIteration

        vec, indx = state()
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
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
    arg_types = ('virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # 25.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_lin_convect )
        
    ##
    # 25.07.2007, c
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, par, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_ep, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            shape = (chunk_size, 1, dim * n_ep, dim * n_ep )
            mode = 1
        else:
            raise StopIteration

        vec1, i1 = par()
        vec2, i2 = state()
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
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
    arg_types = ('parameter', 'state')
    geometry = [(Volume, 'state')]

    ##
    # 30.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_lin_convect )
        
    ##
    # 30.07.2007, c
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par, state = self.get_args( **kwargs )
        ap, vg = state.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, n_qp, dim, 1 )
            mode = 2
        else:
            raise StopIteration

        vec1, i1 = par()
        vec2, i2 = state()
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
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
    arg_types = ('virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # c: 15.12.2005, 31.03.2008
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_grad )
        
    ##
    # c: 31.03.2008, r: 31.03.2008
    def get_shape( self, diff_var, chunk_size, apr, apc = None ):
        self.data_shape = apr.get_v_data_shape( self.integral_name ) 
        n_el, n_qp, dim, n_epr = self.data_shape

        if diff_var is None:
            return (chunk_size, 1, dim * n_epr, 1 ), 0
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = apc.get_v_data_shape( self.integral_name )[3]
            return (chunk_size, 1, dim * n_epr, n_epc ), 1
        else:
            raise StopIteration

    ##
    # c: 31.03.2008, r: 31.03.2008
    def build_c_fun_args( self, state, apc, vgr, **kwargs ):
        vec, indx = state()
        bf = apc.get_base( 'v', 0, self.integral_name )
        return 1.0, vec, indx.start, bf, vgr, apc.econn

    ##
    # c: 15.12.2005, r: 04.07.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, state = self.get_args( ['virtual', 'state'], **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(), 'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(), 'Volume' )

        shape, mode = self.get_shape( diff_var, chunk_size, apr, apc )
        fargs = self.build_c_fun_args( state, apc, vgr, **kwargs )
        
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, *fargs + (chunk, mode) )
            yield out, chunk, status

##
# 30.07.2007, c
class GradQTerm( Term ):
    r""":description: Gradient term (weak form) in quadrature points.
    :definition: $(\nabla p)|_{qp}$
    """
    name = 'dq_grad'
    arg_types = ('state',)
    geometry = [(Volume, 'state')]

    ##
    # 30.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dq_grad_scalar )

    ##
    # 30.07.2007, c
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        state = self.get_args( **kwargs )
        ap, vg = state.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, n_qp, dim, 1 )
            mode = 0
        else:
            raise StopIteration

        vec, indx = state()
        for out, chunk in self.char_fun( chunk_size, shape ):
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
    name = 'dw_grad_dt'
    arg_types = ('ts', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # c: 31.03.2008, r: 31.03.2008
    def build_c_fun_args( self, state, apc, vgr, **kwargs ):
        ts, state_0 = self.get_args( ['ts', 'parameter'], **kwargs )

        vec, indx = state()
        vec0, indx0 = state_0()
        dvec = vec[indx] - vec0[indx0]
        idt = 1.0/ts.dt 

        bf = apc.get_base( 'v', 0, self.integral_name )
        return idt, dvec, 0, bf, vgr, apc.econn

##
# 14.12.2005, c
class DivTerm( Term ):
    r""":description: Divergence term (weak form).
    :definition: $\int_{\Omega} q\ \nabla \cdot \ul{u}$
    """
    name = 'dw_div'
    arg_types = ('virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # 14.12.2005, c
    # 10.01.2006
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_div )

    ##
    # c: 31.03.2008, r: 31.03.2008
    def get_shape( self, diff_var, chunk_size, apr, apc = None ):
        self.data_shape = apr.get_v_data_shape( self.integral_name ) 
        n_el, n_qp, dim, n_epr = self.data_shape

        if diff_var is None:
            return (chunk_size, 1, n_epr, 1 ), 0
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = apc.get_v_data_shape( self.integral_name )[3]
            return (chunk_size, 1, n_epr, dim * n_epc ), 1
        else:
            raise StopIteration

    ##
    # c: 31.03.2008, r: 31.03.2008
    def build_c_fun_args( self, state, apr, apc, vgc, **kwargs ):
        vec, indx = state()
        bf = apr.get_base( 'v', 0, self.integral_name )
        return vec, indx.start, bf, vgc, apc.econn

    ##
    # c: 14.12.2005, r: 04.07.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, state = self.get_args( ['virtual', 'state'], **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(), 'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(), 'Volume' )

        shape, mode = self.get_shape( diff_var, chunk_size, apr, apc )
        fargs = self.build_c_fun_args( state, apr, apc, vgc, **kwargs )
        
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, *fargs + (chunk, mode) )
            yield out, chunk, status

##
# 13.03.2007, c
class DivIntegratedTerm( Term ):
    r""":description: Integrated divergence term (weak form).
    :definition: $\int_{\Omega} \bar{p}\ \nabla \cdot \ul{w}$
    """
    name = 'd_div'
    arg_types = ('parameter_1', 'parameter_2')
    geometry = [(Volume, 'parameter_1'), (Volume, 'parameter_2')]
    use_caches = {'state_in_volume_qp' : [['parameter_1']],
                 'div_vector' : [['parameter_2']]}

    ##
    # 13.03.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )

    ##
    # c: 13.03.2007, r: 17.01.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par1, par2 = self.get_args( **kwargs )
        apc, vgc = par2.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_epc = apc.get_v_data_shape( self.integral_name )
        shape = (0,)

        cache = self.get_cache( 'state_in_volume_qp', 0 )
        vec1 = cache( 'state', self.get_current_group(), 0, state = par1 )
        cache = self.get_cache( 'div_vector', 0 )
        div2 = cache( 'div', self.get_current_group(), 0, state = par2 )

        for out, chunk in self.char_fun( chunk_size, shape ):
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
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # 26.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_grad_div )
        
    ##
    # 26.07.2007, c
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        gamma, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_ep, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            shape = (chunk_size, 1, dim * n_ep, dim * n_ep )
            mode = 1
        else:
            raise StopIteration

        vec, indx = state()
        for out, chunk in self.char_fun( chunk_size, shape ):
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
    arg_types = ('material', 'virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # 31.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_pspg_c )
        
    ##
    # 31.07.2007, c
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        tau, virtual, par, state = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(), 'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_epr = apr.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, n_epr, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = apc.get_v_data_shape( self.integral_name )[3]
            shape = (chunk_size, 1, n_epr, dim * n_epc )
            mode = 1
        else:
            raise StopIteration

        tau_in_el = fix_scalar_in_el( tau, n_el, nm.float64 )
            
        vec1, i1 = par()
        vec2, i2 = state()
        bf = apc.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec1, i1.start, vec2, i2.start,
                                    tau_in_el, bf, vgr, vgc,
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
    arg_types = ('material', 'virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # 31.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_supg_p )
        
    ##
    # 31.07.2007, c
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        delta, virtual, par, state = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(), 'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_epr = apr.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_epr, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = apc.get_v_data_shape( self.integral_name )[3]
            shape = (chunk_size, 1, dim * n_epr, n_epc )
            mode = 1
        else:
            raise StopIteration

        delta_in_el = fix_scalar_in_el( delta, n_el, nm.float64 )
            
        vec1, i1 = par()
        vec2, i2 = state()
        bf = apr.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec1, i1.start, vec2, i2.start,
                                    delta_in_el, bf, vgr, vgc,
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
    arg_types = ('material', 'virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # 31.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_supg_c )
        
    ##
    # 31.07.2007, c
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        delta, virtual, par, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_ep, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            shape = (chunk_size, 1, dim * n_ep, dim * n_ep )
            mode = 1
        else:
            raise StopIteration

        delta_in_el = fix_scalar_in_el( delta, n_el, nm.float64 )
            
        vec1, i1 = par()
        vec2, i2 = state()
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec1, i1.start, vec2, i2.start,
                                    delta_in_el, bf, vg,
                                    ap.econn, chunk, mode )
            yield out, chunk, status
