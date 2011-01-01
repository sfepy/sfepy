import numpy as nm

from sfepy.base.base import use_method_with_name
from sfepy.terms.terms import Term, terms
from sfepy.terms.terms_base import ScalarScalar

class DiffusionTerm( ScalarScalar, Term ):
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

    def get_fargs_weak( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        vec = self.get_vector( state )

        n_el, n_qp, dim, n_ep = self.data_shape
        if state.is_real():
            return (vec, 0, mat, vg, ap.econn), shape, mode
        else:
            ac = nm.ascontiguousarray
            mode += 1j
            return [(ac(vec.real), 0, mat, vg, ap.econn),
                    (ac(vec.imag), 0, mat, vg, ap.econn)], shape, mode
    
    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, par1, par2 = self.get_args( **kwargs )
        ap, vg = self.get_approximation(par1)

        self.set_data_shape( ap )
        n_el, n_qp, dim, n_ep = self.data_shape

        cache = self.get_cache( 'grad_scalar', 0 )
        gp1 = cache('grad', self, 0, state=par1, get_vector=self.get_vector)
        cache = self.get_cache( 'grad_scalar', 1 )
        gp2 = cache('grad', self, 0, state=par2, get_vector=self.get_vector)

        return (gp1, gp2, mat, vg), (chunk_size, 1, 1, 1), 0

    def set_arg_types( self ):
        if self.mode == 'weak':
            self.function = terms.dw_diffusion
            use_method_with_name( self, self.get_fargs_weak, 'get_fargs' )
        else:
            self.function = terms.d_diffusion
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )
            self.use_caches = {'grad_scalar' : [['parameter_1'],
                                                ['parameter_2']]}

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
            use_method_with_name(self, self.get_fargs_weak, 'get_fargs')
        else:
            self.function = terms.d_laplace
            use_method_with_name(self, self.get_fargs_eval, 'get_fargs')
            self.use_caches = {'grad_scalar' : [['parameter_1'],
                                                ['parameter_2']]}

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
    arg_types = ('material','parameter')

    function = staticmethod(terms.de_diffusion_velocity)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, parameter = self.get_args( **kwargs )
        ap, vg = self.get_approximation(parameter)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

        if diff_var is None:
            shape = (chunk_size, 1, dim, 1)
        else:
            raise StopIteration

        vec = parameter()
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec, 0,
                                    mat, vg, ap.econn, chunk )
            out1 = out / vg.variable( 2 )[chunk]
            yield out1, chunk, status
