import numpy as nm

from sfepy.base.base import use_method_with_name, assert_
from sfepy.terms.terms import Term, terms
from sfepy.terms.terms_base import CouplingVectorScalar

class DivGradTerm( Term ):
    r"""
    :Description:
    Diffusion term.

    :Definition:
    .. math::
        \int_{\Omega} \nu\ \nabla \ul{v} : \nabla \ul{u}

    :Arguments:
        material : :math:`\nu` (viscosity),
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_div_grad'
    arg_types = ('material', 'virtual', 'state')

    function = staticmethod(terms.term_ns_asm_div_grad)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

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

    :Arguments:
        virtual : :math:`\ul{v}`,
        state   : :math:`\ul{u}`
    """
    name = 'dw_convect'
    arg_types = ('virtual', 'state')

    function = staticmethod(terms.term_ns_asm_convect)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, state = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_ep, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            shape = (chunk_size, 1, dim * n_ep, dim * n_ep )
            mode = 1
        else:
            raise StopIteration

        vec = state()
        bf = ap.get_base('v', 0, self.integral)
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

    :Arguments:
        virtual   : :math:`\ul{v}`,
        parameter : :math:`\ul{b}`,
        state     : :math:`\ul{u}`
    """
    name = 'dw_lin_convect'
    arg_types = ('virtual', 'parameter', 'state')

    function = staticmethod(terms.dw_lin_convect)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, par, state = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

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
        bf = ap.get_base('v', 0, self.integral)
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

    :Arguments:
        parameter : :math:`\ul{b}`,
        state     : :math:`\ul{u}`
    """
    name = 'dq_lin_convect'
    arg_types = ('parameter', 'state')

    function = staticmethod(terms.dw_lin_convect)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par, state = self.get_args( **kwargs )
        ap, vg = self.get_approximation(state)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

        if diff_var is None:
            shape = (chunk_size, n_qp, dim, 1 )
            mode = 2
        else:
            raise StopIteration

        vec1 = par()
        vec2 = state()
        bf = ap.get_base('v', 0, self.integral)
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec1, 0, vec2, 0,
                                    bf, vg, ap.econn, chunk, mode )
            yield out, chunk, status

class StokesGrad( CouplingVectorScalar ):

    def get_fargs_grad( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, state = self.get_args(['virtual', 'state'], **kwargs)
        apr, vgr = self.get_approximation(virtual)
        apc, vgc = self.get_approximation(state)

        self.set_data_shape( apr, apc )
        shape, mode = self.get_shape_grad( diff_var, chunk_size )

        vec = self.get_vector( state )
        bf = apc.get_base('v', 0, self.integral)

        if 'material' in self.arg_types:
            coef, = self.get_args(['material'], **kwargs)

        else:
            coef = nm.ones((1, self.data_shape_r[1], 1, 1), dtype=nm.float64)

        return (coef, vec, 0, bf, vgr, apc.econn), shape, mode

class StokesDiv( CouplingVectorScalar ):

    def get_fargs_div( self, diff_var = None, chunk_size = None, **kwargs ):
        state, virtual = self.get_args(['state', 'virtual'], **kwargs)
        apr, vgr = self.get_approximation(virtual)
        apc, vgc = self.get_approximation(state)

        self.set_data_shape( apr, apc )
        shape, mode = self.get_shape_div( diff_var, chunk_size )

        vec = self.get_vector( state )
        bf = apr.get_base('v', 0, self.integral)

        if 'material' in self.arg_types:
            coef, = self.get_args(['material'], **kwargs)

        else:
            coef = nm.ones((1, self.data_shape_r[1], 1, 1), dtype=nm.float64)

        return (coef, vec, 0, bf, vgc, apc.econn), shape, mode

class StokesEval( CouplingVectorScalar ):

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        par_v, par_s = self.get_args(['parameter_v', 'parameter_s'], **kwargs)
        aps, vgs = self.get_approximation(par_s)
        apv, vgv = self.get_approximation(par_v)

        self.set_data_shape( aps, apv )

        return (par_v, par_s, vgv), (chunk_size, 1, 1, 1), 0

    def d_eval( self, out, par_v, par_s, vgv, chunk ):
        cache = self.get_cache( 'state_in_volume_qp', 0 )
        vec = cache('state', self, 0,
                    state=par_s, get_vector=self.get_vector)
        cache = self.get_cache( 'div_vector', 0 )
        div = cache('div', self, 0, state=par_v)

        out_qp = vec[chunk] * div[chunk]

        if 'material' in self.arg_types:
            coef, = self.get_args(['material'], **kwargs)
            out_qp *= coef[chunk]

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

    :Arguments 1:
        virtual : :math:`\ul{v}`,
        state   : :math:`\ul{p}`

    :Arguments 2:
        state   : :math:`\ul{u}`,
        virtual : :math:`\ul{q}`

    :Arguments 3:
        parameter_v : :math:`\ul{u}`,
        parameter_s : :math:`p`
    """
    name = 'dw_stokes'
    arg_types = (('virtual', 'state'),
                 ('state', 'virtual'),
                 ('parameter_v', 'parameter_s'))
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

class StokesWTerm(StokesTerm):
    r"""
    :Description:
    Stokes problem coupling term weighted by a scalar function. Corresponds to
    weighted weak forms of gradient and divergence terms. Can be evaluated.

    :Definition:
    .. math::
        \int_{\Omega} c\ p\ \nabla \cdot \ul{v} \mbox{ , }
        \int_{\Omega} c\ q\ \nabla \cdot \ul{u}

    :Arguments 1:
        material : :math:`c`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`p`

    :Arguments 2:
        material : :math:`c`,
        state    : :math:`\ul{u}`,
        virtual  : :math:`q`

    :Arguments 3:
        material    : :math:`c`,
        parameter_v : :math:`\ul{u}`,
        parameter_s : :math:`p`
    """
    name = 'dw_stokes_w'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_v', 'parameter_s'))
    modes = ('grad', 'div', 'eval')

class GradQTerm( Term ):
    r"""
    :Description:
    Gradient term (weak form) in quadrature points.

    :Definition:
    .. math::
        (\nabla p)|_{qp}

    :Arguments:
        state : :math:`p`
    """
    name = 'dq_grad'
    arg_types = ('state',)

    function = staticmethod(terms.dq_grad)

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        state, = self.get_args( **kwargs )
        ap, vg = self.get_approximation(state)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

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
        \mbox{vector for } K \from \Ical_h: \int_{T_K} \nabla p /
        \int_{T_K} 1 \mbox{ or } \int_{T_K} \nabla \ul{w} /
        \int_{T_K} 1

    :Arguments:
        state : :math:`p` or :math:`\ul{w}`
    """
    name = 'de_grad'
    arg_types = ('state',)

    function = staticmethod(terms.de_grad)

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        state, = self.get_args( **kwargs )
        ap, vg = self.get_approximation(state)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

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

class DivQTerm(Term):
    r"""
    :Description:
    Divergence term (weak form) in quadrature points.

    :Definition:
    .. math::
        (\nabla \cdot \ul{u})|_{qp}

    :Arguments:
        state : :math:`\ul{u}`
    """
    name = 'dq_div'
    arg_types = ('state',)

    function = staticmethod(terms.dq_div_vector)

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        state, = self.get_args(**kwargs)
        ap, vg = self.get_approximation(state)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

        assert_(state.n_components == dim)

        if diff_var is None:
            shape = (chunk_size, n_qp, 1, 1)
            mode = 0
        else:
            raise StopIteration

        vec = state()
        for out, chunk in self.char_fun(chunk_size, shape):
            status = self.function(out, vec, 0, vg, ap.econn[chunk])
            yield out, chunk, status

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

    :Arguments:
        material : :math:`\gamma`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_st_grad_div'
    arg_types = ('material', 'virtual', 'state')

    function = staticmethod(terms.dw_st_grad_div)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        gamma, virtual, state = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

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
        \sum_{K \in \Ical_h}\int_{T_K} \tau_K\ \nabla p \cdot \nabla q

    :Arguments:
        material : :math:`\tau_K`,
        virtual  : :math:`q`,
        state    : :math:`p`
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
        \sum_{K \in \Ical_h}\int_{T_K} \tau_K\ ((\ul{b} \cdot \nabla) \ul{u})
        \cdot \nabla q

    :Arguments:
        material  : :math:`\tau_K`,
        virtual   : :math:`q`,
        parameter : :math:`\ul{b}`,
        state     : :math:`\ul{u}`
    """
    name = 'dw_st_pspg_c'
    arg_types = ('material', 'virtual', 'parameter', 'state')

    function = staticmethod(terms.dw_st_pspg_c)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        tau, virtual, par, state = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)
        apr, vgr = self.get_approximation(virtual)
        apc, vgc = self.get_approximation(state)
        n_el, n_qp, dim, n_epr = apr.get_v_data_shape(self.integral)

        if diff_var is None:
            shape = (chunk_size, 1, n_epr, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = apc.get_v_data_shape(self.integral)[3]
            shape = (chunk_size, 1, n_epr, dim * n_epc )
            mode = 1
        else:
            raise StopIteration

        vec1 = par()
        vec2 = state()
        bf = apc.get_base('v', 0, self.integral)
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
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ \nabla p\cdot ((\ul{b} \cdot
        \nabla) \ul{v})

    :Arguments:
        material  : :math:`\delta_K`,
        virtual   : :math:`\ul{v}`,
        parameter : :math:`\ul{b}`,
        state     : :math:`p`
    """
    name = 'dw_st_supg_p'
    arg_types = ('material', 'virtual', 'parameter', 'state')

    function = staticmethod(terms.dw_st_supg_p)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        delta, virtual, par, state = self.get_args( **kwargs )
        apr, vgr = self.get_approximation(virtual)
        apc, vgc = self.get_approximation(state)
        n_el, n_qp, dim, n_epr = apr.get_v_data_shape(self.integral)

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_epr, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = apc.get_v_data_shape(self.integral)[3]
            shape = (chunk_size, 1, dim * n_epr, n_epc )
            mode = 1
        else:
            raise StopIteration

        vec1 = par()
        vec2 = state()
        bf = apr.get_base('v', 0, self.integral)
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
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ ((\ul{b} \cdot \nabla)
        \ul{u})\cdot ((\ul{b} \cdot \nabla) \ul{v})

    :Arguments:
        material  : :math:`\delta_K`,
        virtual   : :math:`\ul{v}`,
        parameter : :math:`\ul{b}`,
        state     : :math:`\ul{u}`
    """
    name = 'dw_st_supg_c'
    arg_types = ('material', 'virtual', 'parameter', 'state')

    function = staticmethod(terms.dw_st_supg_c)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        delta, virtual, par, state = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

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
        bf = ap.get_base('v', 0, self.integral)
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec1, 0, vec2, 0,
                                    delta, bf, vg,
                                    ap.econn, chunk, mode )
            yield out, chunk, status
