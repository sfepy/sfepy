import numpy as nm

from sfepy.base.base import use_method_with_name
from sfepy.terms.terms import Term, terms
from sfepy.terms.terms_base import ScalarScalar

class DiffusionTerm(Term):
    r"""
    :Description:
    General diffusion term with permeability :math:`K_{ij}`. Can be
    evaluated. Can use derivatives.

    :Definition:
    .. math::
        \int_{\Omega} K_{ij} \nabla_i q \nabla_j p \mbox{ , } \int_{\Omega}
        K_{ij} \nabla_i \bar{p} \nabla_j r

    :Arguments 1:
        material : :math:`K_{ij}`,
        virtual  : :math:`q`,
        state    : :math:`p`

    :Arguments 2:
        material    : :math:`K_{ij}`,
        parameter_1 : :math:`\bar{p}`,
        parameter_2 : :math:`r`
    """
    name = 'dw_diffusion'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    modes = ('weak', 'eval')
    symbolic = {'expression': 'div( K * grad( u ) )',
                'map' : {'u' : 'state', 'K' : 'material'}}

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        if mode == 'weak':
            if diff_var is None:
                grad = self.get(state, 'grad')
                fmode = 0

            else:
                grad = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return grad, mat, vg, fmode

        elif mode == 'eval':
            grad1 = self.get(virtual, 'grad')
            grad2 = self.get(state, 'grad')

            return grad1, grad2, mat, vg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, mat, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, 1, 1), state.dtype

    def set_arg_types(self):
        if self.mode == 'weak':
            self.function = terms.dw_diffusion

        else:
            self.function = terms.d_diffusion

class LaplaceTerm(DiffusionTerm):
    r"""
    :Description:
    Laplace term with :math:`c` coefficient. Can be
    evaluated. Can use derivatives.

    :Definition:
    .. math::
        \int_{\Omega} c \nabla q \cdot \nabla p \mbox{ , } \int_{\Omega}
        c \nabla \bar{p} \cdot \nabla r

    :Arguments 1:
        material : :math:`c`,
        virtual  : :math:`q`,
        state    : :math:`p`

    :Arguments 2:
        material    : :math:`c`,
        parameter_1 : :math:`\bar{p}`,
        parameter_2 : :math:`r`
    """
    name = 'dw_laplace'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    modes = ('weak', 'eval')
    symbolic = {'expression': 'c * div( grad( u ) )',
                'map' : {'u' : 'state', 'c' : 'material'}}

    def set_arg_types(self):
        if self.mode == 'weak':
            self.function = terms.dw_laplace

        else:
            self.function = terms.d_laplace

class PermeabilityRTerm( Term ):
    r"""
    :Description:
    Special-purpose diffusion-like term with permeability :math:`K_{ij}` (to
    use on the right-hand side).

    :Definition:
    .. math::
        \int_{\Omega} K_{ij} \nabla_j q

    :Arguments:
        material : :math:`K_{ij}`,
        virtual  : :math:`q`,
        index    : :math:`i`
    """
    name = 'dw_permeability_r'
    arg_types = ('material', 'virtual', 'index')

    function = staticmethod(terms.dw_permeability_r)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, index = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

        if diff_var is None:
            shape = (chunk_size, 1, n_ep, 1)
        else:
            raise StopIteration

        if isinstance(index, list):
            index = index[0]

        mat = nm.ascontiguousarray(mat[...,index:index+1])
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, mat, vg, ap.econn, chunk )
            yield out, chunk, status

class DiffusionRTerm( PermeabilityRTerm ):
    r"""
    :Description:
    Diffusion-like term with material parameter :math:`K_{j}` (to
    use on the right-hand side).

    :Definition:
    .. math::
        \int_{\Omega} K_{j} \nabla_j q

    :Arguments:
        material : :math:`K_j`,
        virtual  : :math:`q`
    """
    name = 'dw_diffusion_r'
    arg_types = ('material', 'virtual')

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

        if diff_var is None:
            shape = (chunk_size, 1, n_ep, 1)
        else:
            raise StopIteration

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, mat, vg, ap.econn, chunk )
            yield out, chunk, status

class DiffusionCoupling(ScalarScalar, Term):
    r"""
    :Description:
    Diffusion copupling term with material parameter :math:`K_{j}`.

    :Definition:
    .. math::
        \int_{\Omega}  p K_{j} \nabla_j q

    :Arguments:
        material : :math:`K_{j}`,
        virtual  : :math:`q`,
        state    : :math:`p`
    """
    name = 'dw_diffusion_coupling'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_1', 'parameter_2'))
    modes = ('weak0', 'weak1', 'eval')
    functions = {'weak': terms.dw_diffusion_coupling,
                 'eval': terms.d_diffusion_coupling}

    def get_fargs_weak( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        term_mode = self.mode[-1]
        vec = self.get_vector( state )
        bf = ap.get_base('v', 0, self.integral)

        n_el, n_qp, dim, n_ep = self.data_shape
        if state.is_real():
            return (vec, 0, mat, vg, bf, ap.econn), shape, mode, term_mode
        else:
            ac = nm.ascontiguousarray
            mode += 1j
            return [(ac(vec.real), 0, mat, bf, vg, ap.econn),
                    (ac(vec.imag), 0, mat, bf, vg, ap.econn)], shape, mode, term_mode

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, par1, par2 = self.get_args( **kwargs )
        ap, vg = self.get_approximation(par1)

        self.set_data_shape( ap )
        n_el, n_qp, dim, n_ep = self.data_shape
        bf = ap.get_base('v', 0, self.integral)

        return (par1(), par2(), mat, bf, vg), (chunk_size, 1, 1, 1), 0

    def set_arg_types( self ):
        if self.mode[:-1] == 'weak':
            self.function = self.functions['weak']
            use_method_with_name( self, self.get_fargs_weak, 'get_fargs' )
        else:
            self.function = self.functions['eval']
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )

class DiffusionVelocityTerm( Term ):
    r"""
    :Description:
    Diffusion velocity averaged in elements.

    :Definition:
    .. math::
        \mbox{vector for } K \from \Ical_h: \int_{T_K} -K_{ij} \nabla_j r
        / \int_{T_K} 1

    :Arguments:
        material  : :math:`K_{ij}`,
        parameter : :math:`r`
    """
    name = 'de_diffusion_velocity'
    arg_types = ('material', 'parameter')

    function = staticmethod(terms.de_diffusion_velocity)

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        grad = self.get(parameter, 'grad')

        fmode = {'eval' : 0, 'el_avg' : 1}.get(mode, 1)

        return grad, mat, vg, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_el, 1, dim, 1), parameter.dtype

class DiffusionIntegrateTerm( Term ):
    r"""
    :Description:
    Diffusion integrate term.

    :Definition:
    .. math::
        \int_{\Omega} K_{ij} \nabla_j \bar{p}

    :Arguments:
        material: :math:`\uv{K}`,
        parameter:  :math:`\bar{p}`,
    """
    name = 'di_diffusion_integrate'
    arg_types = ('material', 'parameter')

    @staticmethod
    def function(out, grad, mat, vg):

        aux = nm.sum(mat * grad, axis=3)[:,:,:,nm.newaxis]
        status = vg.integrate(out, nm.ascontiguousarray(aux), 0)

        return status

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)
        grad = self.get(parameter, 'grad')

        return grad, mat, vg

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, dim, 1), parameter.dtype

    # def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
    #     mat, parameter = self.get_args( **kwargs )
    #     ap, vg = self.get_approximation(parameter)
    #     n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

    #     if diff_var is None:
    #         shape = (chunk_size, 1, dim, 1)
    #     else:
    #         raise StopIteration

    #     vec = parameter()
    #     for out, chunk in self.char_fun( chunk_size, shape ):
    #         status = self.function(out, vec, mat, vg, ap.econn, chunk)
    #         out1 = nm.sum(out, 0)
    #         out1.shape = (dim,)
    #         yield out1, chunk, status

class SurfaceFluxTerm(Term):
    r"""
    :Description:
    Surface flux term.

    :Definition:
    .. math::
        \int_{\Gamma} \ul{n} \cdot K_{ij} \nabla_j \bar{p}

    :Arguments:
        material: :math:`\ul{K}`,
        parameter:  :math:`\bar{p}`,
    """
    name = 'd_surface_flux'
    arg_types = ('material', 'parameter')
    integration = 'surface_extra'

    function = staticmethod(terms.d_surface_flux)

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(parameter)

        grad = self.get(parameter, 'grad')

        fmode = {'eval' : 0, 'el_avg' : 1}.get(mode, 1)

        return grad, mat, sg, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_fa, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_fa, 1, 1, 1), parameter.dtype
