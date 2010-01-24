from sfepy.terms.terms import *
from sfepy.terms.terms_base import CouplingVectorScalar

class DivGradTerm( Term ):
    r"""
    :Description:
    Diffusion term.

    :Definition:
    .. math::
        \int_{\Omega} \nu\ \nabla \ul{v} : \nabla \ul{u}
    """
    name = 'dw_div_grad'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.term_ns_asm_div_grad )
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( **kwargs )
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

        vec = state()
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec, 0, mat,
                                    vg, ap.econn, chunk, mode )
            yield out, chunk, status

class ConvectTerm( Term ):
    r"""
    :Description:
    Nonlinear convective term.

    :Definition:
    .. math::
        \int_{\Omega} ((\ul{u} \cdot \nabla) \ul{u}) \cdot \ul{v}
    """
    name = 'dw_convect'
    arg_types = ('virtual', 'state')
    geometry = [(Volume, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.term_ns_asm_convect )
        
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

        vec = state()
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec, 0, bf,
                                    vg, ap.econn, chunk, mode )
            yield out, chunk, status

class LinearConvectTerm( Term ):
    r"""
    :Description:
    Linearized convective term.

    :Definition:
    .. math::
        \int_{\Omega} ((\ul{b} \cdot \nabla) \ul{u}) \cdot \ul{v}
    """
    name = 'dw_lin_convect'
    arg_types = ('virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_lin_convect )
        
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

        vec1 = par()
        vec2 = state()
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec1, 0, vec2, 0,
                                    bf, vg, ap.econn, chunk, mode )
            yield out, chunk, status

class LinearConvectQTerm( Term ):
    r"""
    :Description:
    Linearized convective term evaluated in quadrature points.

    :Definition:
    .. math::
        ((\ul{b} \cdot \nabla) \ul{u})|_{qp}
    """
    name = 'dq_lin_convect'
    arg_types = ('parameter', 'state')
    geometry = [(Volume, 'state')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_lin_convect )
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par, state = self.get_args( **kwargs )
        ap, vg = state.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, n_qp, dim, 1 )
            mode = 2
        else:
            raise StopIteration

        vec1 = par()
        vec2 = state()
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec1, 0, vec2, 0,
                                    bf, vg, ap.econn, chunk, mode )
            yield out, chunk, status

class StokesGrad( CouplingVectorScalar ):

    def get_fargs_grad( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, state = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(),
                                              'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(),
                                            'Volume' )

        self.set_data_shape( apr, apc )
        shape, mode = self.get_shape_grad( diff_var, chunk_size )

        vec = self.get_vector( state )
        bf = apc.get_base( 'v', 0, self.integral_name )
        return (1.0, vec, 0, bf, vgr, apc.econn), shape, mode

class StokesDiv( CouplingVectorScalar ):

    def get_fargs_div( self, diff_var = None, chunk_size = None, **kwargs ):
        state, virtual = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(),
                                              'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(),
                                            'Volume' )

        self.set_data_shape( apr, apc )
        shape, mode = self.get_shape_div( diff_var, chunk_size )

        vec = self.get_vector( state )
        bf = apr.get_base( 'v', 0, self.integral_name )
        return (vec, 0, bf, vgc, apc.econn), shape, mode

class StokesEval( CouplingVectorScalar ):

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        par_v, par_s = self.get_args( **kwargs )
        aps, vgs = par_s.get_approximation( self.get_current_group(),
                                            'Volume' )
        apv, vgv = par_v.get_approximation( self.get_current_group(),
                                            'Volume' )
        self.set_data_shape( aps, apv )
        return (par_v, par_s, vgv), (chunk_size, 1, 1, 1), 0

    def d_eval( self, out, par_v, par_s, vgv, chunk ):
        cache = self.get_cache( 'state_in_volume_qp', 0 )
        vec = cache( 'state', self.get_current_group(), 0,
                     state = par_s, get_vector = self.get_vector )
        cache = self.get_cache( 'div_vector', 0 )
        div = cache( 'div', self.get_current_group(), 0, state = par_v )

        out_qp = vec[chunk] * div[chunk]

        status = vgv.integrate_chunk( out, out_qp, chunk )
        
        return status

class StokesTerm( StokesDiv, StokesGrad, StokesEval, Term ):
    r"""
    :Description:
    Stokes problem coupling term. Corresponds to weak forms of gradient and
    divergence terms. Can be evaluated.

    :Definition:
    .. math::
        \int_{\Omega} p\ \nabla \cdot \ul{v} \mbox{ , } \int_{\Omega} q\ \nabla
        \cdot \ul{u}
    """
    name = 'dw_stokes'
    arg_types = (('virtual', 'state'),
                 ('state', 'virtual'),
                 ('parameter_v', 'parameter_s'))
    geometry = ([(Volume, 'virtual'), (Volume, 'state')],
                [(Volume, 'virtual'), (Volume, 'state')],
                [(Volume, 'parameter_v'), (Volume, 'parameter_s')])
    modes = ('grad', 'div', 'eval')

    def set_arg_types( self ):
        """Dynamically inherits from either StokesGrad or
        StokesDiv."""
        if self.mode == 'grad':
            self.function = terms.dw_grad
            use_method_with_name( self, self.get_fargs_grad, 'get_fargs' )
        elif self.mode == 'div':
            self.function = terms.dw_div
            use_method_with_name( self, self.get_fargs_div, 'get_fargs' )
        else:
            self.function = self.d_eval
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )
            self.use_caches = {'state_in_volume_qp' : [['parameter_s']],
                               'div_vector' : [['parameter_v']]}

class GradQTerm( Term ):
    r"""
    :Description:
    Gradient term (weak form) in quadrature points.

    :Definition:
    .. math::
        (\nabla p)|_{qp}
    """
    name = 'dq_grad'
    arg_types = ('state',)
    geometry = [(Volume, 'state')]

    ##
    # 30.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dq_grad )

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        state, = self.get_args( **kwargs )
        ap, vg = state.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, n_qp, dim, 1 )
            mode = 0
        else:
            raise StopIteration

        vec = state()
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec, 0, vg, ap.econn[chunk] )
            yield out, chunk, status

class GradETerm( Term ):
    r"""
    :Description:
    Gradient term (weak form) in averaged in elements.
    
    :Definition:
    .. math::
        \mbox{vector of } \forall K \in \Tcal_h: \int_{T_K} \nabla p /
        \int_{T_K} 1 \mbox{ or } \int_{T_K} \nabla \ul{w} /
        \int_{T_K} 1
    """
    name = 'de_grad'
    arg_types = ('state',)
    geometry = [(Volume, 'state')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.de_grad )

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        state, = self.get_args( **kwargs )
        ap, vg = state.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        vdim = ap.dim[0]
        
        if diff_var is None:
            shape = (chunk_size, 1, dim, vdim )
            mode = 0
        else:
            raise StopIteration

        ac = nm.ascontiguousarray

        vec = state()
        for out, chunk in self.char_fun( chunk_size, shape ):
            if state.is_real():
                status = self.function( out, vec, 0, vg, ap.econn, chunk )
            else:
                status_r = self.function(out, ac(vec.real), 0,
                                         vg, ap.econn, chunk)
                out_imag = nm.zeros_like(out)
                status_i = self.function(out_imag, ac(vec.imag), 0,
                                         vg, ap.econn, chunk)

                status = status_r or status_i
                out = out + 1j * out_imag

            out1 = out / vg.variable( 2 )[chunk]
            yield out1, chunk, status
                
##
# 26.07.2007, c
class GradDivStabilizationTerm( Term ):
    r"""
    :Description:
    Grad-div stabilization term ( :math:`\gamma` is a global stabilization
    parameter).

    :Definition:
    .. math::
        \gamma \int_{\Omega} (\nabla\cdot\ul{u}) \cdot (\nabla\cdot\ul{v})
    """
    name = 'dw_st_grad_div'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # 26.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_grad_div )
        
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

        vec = state()
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec, 0, gamma,
                                    vg, ap.econn, chunk, mode )
            yield out, chunk, status

##
# 31.07.2007, c
from sfepy.terms.termsLaplace import LaplaceTerm
class PSPGPStabilizationTerm( LaplaceTerm ):
    r"""
    :Description:
    PSPG stabilization term, pressure part ( :math:`\tau` is a local
    stabilization parameter), alias to Laplace term dw_laplace.

    :Definition:
    .. math::
        \sum_{K \in \Tcal_h}\int_{T_K} \tau_K\ \nabla p \cdot \nabla q
    """
    name = 'dw_st_pspg_p'

##
# 31.07.2007, c
class PSPGCStabilizationTerm( Term ):
    r"""
    :Description:
    PSPG stabilization term, convective part ( :math:`\tau` is a local
    stabilization parameter).

    :Definition:
    .. math::
        \sum_{K \in \Tcal_h}\int_{T_K} \tau_K\ ((\ul{b} \cdot \nabla) \ul{u})
        \cdot \nabla q
    """
    name = 'dw_st_pspg_c'
    arg_types = ('material', 'virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # 31.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_pspg_c )
        
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

        vec1 = par()
        vec2 = state()
        bf = apc.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec1, 0, vec2, 0,
                                    tau, bf, vgr, vgc,
                                    apc.econn, chunk, mode )
            yield out, chunk, status

##
# 31.07.2007, c
class SUPGPStabilizationTerm( Term ):
    r"""
    :Description:
    SUPG stabilization term, pressure part ( :math:`\delta` is a local
    stabilization parameter).

    :Definition:
    .. math::
        \sum_{K \in \Tcal_h}\int_{T_K} \delta_K\ \nabla p\cdot ((\ul{b} \cdot
        \nabla) \ul{v})
    """
    name = 'dw_st_supg_p'
    arg_types = ('material', 'virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    ##
    # 31.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_supg_p )
        
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

        vec1 = par()
        vec2 = state()
        bf = apr.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec1, 0, vec2, 0,
                                    delta, bf, vgr, vgc,
                                    apr.econn, apc.econn, chunk, mode )
            yield out, chunk, status

##
# 31.07.2007, c
class SUPGCStabilizationTerm( Term ):
    r"""
    :Description:
    SUPG stabilization term, convective part ( :math:`\delta` is a local
    stabilization parameter).

    :Definition:
    .. math::
        \sum_{K \in \Tcal_h}\int_{T_K} \delta_K\ ((\ul{b} \cdot \nabla)
        \ul{u})\cdot ((\ul{b} \cdot \nabla) \ul{v})
    """
    name = 'dw_st_supg_c'
    arg_types = ('material', 'virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual')]

    ##
    # 31.07.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_supg_c )
        
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

        vec1 = par()
        vec2 = state()
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec1, 0, vec2, 0,
                                    delta, bf, vg,
                                    ap.econn, chunk, mode )
            yield out, chunk, status
