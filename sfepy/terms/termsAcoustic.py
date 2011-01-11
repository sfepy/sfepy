import numpy as nm

from sfepy.base.base import use_method_with_name
from sfepy.terms.terms import Term, terms
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

class SurfaceLaplaceLayerTerm(ScalarScalar, Term):
    r"""
    :Description:
    Acoustic 'layer' term - derivatives in surface directions.

    :Definition:
    .. math::
        \int_{\Gamma} c \partial_\alpha \ul{v}\,\partial_\alpha \ul{u}, \alpha = 1,\dots,N-1

    :Arguments:
        material: :math:`c`,
        virtual:  :math:`\ul{v}`,
        state:    :math:`\ul{u}`
    """
    name = 'dw_surface_llaplace'
    arg_types = ('material', 'virtual', 'state')
    integration = 'surface'

    function = staticmethod(terms.dw_surf_llaplace)

    def get_fargs(self, diff_var = None, chunk_size = None, **kwargs):
        coef, virtual, state = self.get_args(**kwargs)
        ap, sg = self.get_approximation(virtual)
        aps, sgs = self.get_approximation(state)

        self.set_data_shape(ap)
        shape, mode = self.get_shape(diff_var, chunk_size)

        vec = self.get_vector(state)
        sd = aps.surface_data[self.region.name]

        bfg = ap.get_base(sd.face_type, 1, self.integral)
        econn = sd.get_connectivity(state.is_surface)

        if state.is_real():
            fargs = vec, coef, bfg, sgs, econn

        else:
            ac = nm.ascontiguousarray
            fargs = [(ac(vec.real), coef, bfg, sgs, econn),
                     (ac(vec.imag), coef, bfg, sgs, econn)]
            mode += 1j

        return fargs, shape, mode

class SurfaceCoupleLayerTerm(ScalarScalar, Term):
    r"""
    :Description:
    Acoustic 'layer' term - derivatives in surface directions.

    :Definition:
    .. math::
        \int_{\Gamma} c \ul{v}\,\partial_\alpha \ul{u},
        \int_{\Gamma} c \partial_\alpha \ul{u}\,\ul{v}, \alpha = 1,\dots,N-1

    :Arguments 1:
        material: :math:`c`,
        virtual:  :math:`\ul{v}`,
        state:    :math:`\ul{u}`

    :Arguments 2:
        material: :math:`c`,
        state:    :math:`\ul{u}`
        virtual:  :math:`\ul{v}`,

    """
    name = 'dw_surface_lcouple'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'))
    modes = ('bv_ns', 'nv_bs')
    integration = 'surface'

    function = staticmethod(terms.dw_surf_lcouple)

    def get_fargs(self, diff_var = None, chunk_size = None, **kwargs):
        if self.mode == 'nv_bs':
            coef, state, virtual = self.get_args(**kwargs)
        else:
            coef, virtual, state = self.get_args(**kwargs)

        ap, sg = self.get_approximation(virtual)
        aps, sgs = self.get_approximation(state)

        self.set_data_shape(ap)
        shape, mode = self.get_shape(diff_var, chunk_size)

        vec = self.get_vector(state)
        sd = aps.surface_data[self.region.name]

        bf = ap.get_base(sd.face_type, 0, self.integral)
        bfg = ap.get_base(sd.face_type, 1, self.integral)
        econn = sd.get_connectivity(state.is_surface)

        aux = coef.shape

        if self.mode == 'nv_bs':
            bf, bfg = bfg, bf
            coef.reshape((aux[0], aux[1], 1, nm.max(aux[2:])))

        else:
            coef.reshape((aux[0], aux[1], nm.max(aux[2:]), 1))

        if state.is_real():
            fargs = vec, coef, bf, bfg, sgs, econn

        else:
            ac = nm.ascontiguousarray
            fargs = [(ac(vec.real), coef, bf, bfg, sgs, econn),
                     (ac(vec.imag), coef, bf, bfg, sgs, econn)]
            mode += 1j

        return fargs, shape, mode
