from sfepy.terms.terms import *
from sfepy.terms.terms_base import ScalarScalar

class LaplaceLayerPSA1Term(ScalarScalar, Term):
    r"""
    :Description:
    Sensitivity analysis term -- in-plane directions.

    :Definition:
    .. math::
        \int_{\Omega} \partial_\alpha w_k\, \partial_k \ul{u}\, \partial_\alpha
        \ul{v}, \alpha = 1,\dots,N-1

    :Arguments:
        parametr_1: :math:`\ul{u}`,
        parametr_2: :math:`\ul{v}`,
        parametr_3: :math:`w`
    """
    name = 'd_llaplace_p_sa1'
    arg_types = ('parameter_1', 'parameter_2', 'parameter_3')

    function = staticmethod(terms.d_llaplace_p_sa)

    sa_mode = 0

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        par1, par2, par3 = self.get_args( ['parameter_1', 'parameter_2',
                                           'parameter_3'], **kwargs )
        ap, vg = self.get_approximation(par1)

        self.set_data_shape( ap )

        fargs = ( par1(), par2(), par3(), vg, ap.econn, self.sa_mode )
        return fargs, (chunk_size, 1, 1, 1), 0
    
class LaplaceLayerPSA2Term(LaplaceLayerPSA1Term):
    r"""
    :Description:
    Sensitivity analysis term -- in-plane directions.

    :Definition:
    .. math::
        \int_{\Omega} \dvg w \partial_\alpha \ul{u}\, \partial_\alpha \ul{v},
        \alpha = 1,\dots,N-1

    :Arguments:
        parametr_1: :math:`\ul{u}`,
        parametr_2: :math:`\ul{v}`,
        parametr_3: :math:`w`
    """

    name = 'd_llaplace_p_sa2'

    sa_mode = 1

class LaplaceLayerTSA1Term(LaplaceLayerPSA1Term):
    r"""
    :Description:
    Sensitivity analysis term -- transversal direction.

    :Definition:
    .. math::
        \int_{\Omega} \partial_N w_k\, \partial_k \ul{v}\,\partial_N \ul{u}

    :Arguments:
        parametr_1: :math:`\ul{u}`,
        parametr_2: :math:`\ul{v}`,
        parametr_3: :math:`w`
    """

    name = 'd_llaplace_t_sa1'

    function = staticmethod(terms.d_llaplace_t_sa)

    sa_mode = 0

class LaplaceLayerTSA2Term(LaplaceLayerTSA1Term):
    r"""
    :Description:
    Sensitivity analysis acoustic term (transversal direction).

    :Definition:
    .. math::
        \int_{\Omega} \dvg w \partial_N \ul{v}\,\partial_N \ul{u}

    :Arguments:
        parametr_1: :math:`\ul{u}`,
        parametr_2: :math:`\ul{v}`,
        parametr_3: :math:`w`
    """

    name = 'd_llaplace_t_sa2'

    sa_mode = 1

class LaplaceLayerTerm( ScalarScalar, Term ):
    r"""
    :Description:
    Acoustic term.

    :Definition:
    .. math::
        \int_{\Omega} (c_1 \partial_\alpha \ul{v}\,\partial_\alpha \ul{u} + c_2
        \partial_N \ul{v}\,\partial_N \ul{u} ), \alpha = 1,\dots,N-1

    :Arguments:
        material_1: :math:`c_1`,
        material_2: :math:`c_2`,
        virtual:    :math:`\ul{v}`,
        state:      :math:`\ul{u}`
    """
    name = 'dw_llaplace'
    arg_types = (('material', 'material', 'virtual', 'state'),
                 ('material', 'material', 'parameter_1', 'parameter_2'))
    modes = ('weak', 'eval')
    functions = {'weak': terms.dw_llaplace,
                 'eval': terms.d_llaplace }
        
    def get_fargs_weak( self, diff_var = None, chunk_size = None, **kwargs ):
        mat1, mat2, virtual, state = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)
        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        return (state(), mat1, mat2, vg, ap.econn), shape, mode

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        mat1, mat2, par1, par2 = self.get_args( **kwargs )
        ap, vg = self.get_approximation(par1)
        self.set_data_shape( ap )

        return ((par1(), par2(), mat1, mat2, vg, ap.econn),
                (chunk_size, 1, 1, 1), 0)

    def set_arg_types( self ):
        if self.mode == 'weak':
            self.function = self.functions['weak']
            use_method_with_name( self, self.get_fargs_weak, 'get_fargs' )
        else:
            self.function = self.functions['eval']
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )

class AcousticSurfaceTerm( ScalarScalar, Term ):
    r"""
    :Description:
    Acoustic surface term (in-plane directions).

    :Definition:
    .. math::
        \int_{\Gamma} \ul{n}\cdot \partial_{\alpha, 1/h z} \ul{y},
        \alpha = 1,\dots,N-1

    :Arguments:
        material_1: :math:`p_1`,
        material_2: :math:`p_2`,
        parameter:  :math:`z`
    """
    name = 'd_acoustic_surface'
    arg_types = ('material', 'material', 'parameter')
    integration = 'surface_extra'

    function = staticmethod(terms.d_acoustic_surface)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        mat1, mat2, par = self.get_args( **kwargs )
        ap, sg = self.get_approximation(par)

        self.set_data_shape( ap )

        fargs = (par(), mat1, mat2, sg, ap.econn)

        return fargs, (chunk_size, 1, 1, 1), 0

class AcousticIntegrateTerm( ScalarScalar, Term ):
    r"""
    :Description:
    Integration of acoustic term (in-plane directions).

    :Definition:
    .. math::
        \int_{\Omega} m  \partial_\alpha \ul{u},
        \alpha = 1,\dots,N-1

    :Arguments:
        material: :math:`m`,
        virtual:  :math:`\ul{v}`
    """
    name = 'dw_acoustic_integrate'
    arg_types = ('material', 'virtual')

    function = staticmethod(terms.dw_acoustic_integrate)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)
        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        return (mat, vg, ap.econn), shape, mode

class AcousticEvalAlphaTerm( Term ):
    r"""
    :Description:
    Evaluation of acoustic term (in-plane directions).

    :Definition:
    .. math::
        \int_{\Omega} \partial_{\alpha} \ul{y},
        \alpha = 1,\dots,N-1

    :Arguments:
        parameter: :math:`y`
    """
    name = 'd_acoustic_alpha'
    arg_types = ('parameter',)

    function = staticmethod(terms.d_acoustic_alpha)

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par, = self.get_args( **kwargs )
        ap, vg = self.get_approximation()
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

        if diff_var is None:
            shape = (1, 1, dim-1, 1)
        else:
            raise StopIteration

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, par(), vg, ap.econn, chunk )
            out1 = nm.sum( out, 0 )
            yield out1, chunk, status
