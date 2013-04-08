import numpy as nm

from sfepy.linalg import dot_sequences
from sfepy.terms.terms import Term, terms

class DiffusionTerm(Term):
    r"""
    General diffusion term with permeability :math:`K_{ij}`. Can be
    evaluated. Can use derivatives.

    :Definition:

    .. math::
        \int_{\Omega} K_{ij} \nabla_i q \nabla_j p \mbox{ , } \int_{\Omega}
        K_{ij} \nabla_i \bar{p} \nabla_j r

    :Arguments 1:
        - material : :math:`K_{ij}`
        - virtual  : :math:`q`
        - state    : :math:`p`

    :Arguments 2:
        - material    : :math:`K_{ij}`
        - parameter_1 : :math:`\bar{p}`
        - parameter_2 : :math:`r`
    """
    name = 'dw_diffusion'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = {'material' : 'D, D', 'virtual' : (1, 'state'),
                  'state' : 1, 'parameter_1' : 1, 'parameter_2' : 1}
    modes = ('weak', 'eval')
    symbolic = {'expression': 'div( K * grad( u ) )',
                'map' : {'u' : 'state', 'K' : 'material'}}

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        if mat is None:
            if self.name == 'dw_laplace':
                n_el, n_qp, _, _, _ = self.get_data_shape(state)
                mat = nm.ones((1, n_qp, 1, 1), dtype=nm.float64)

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
    Laplace term with :math:`c` coefficient. Can be
    evaluated. Can use derivatives.

    :Definition:

    .. math::
        \int_{\Omega} c \nabla q \cdot \nabla p \mbox{ , } \int_{\Omega}
        c \nabla \bar{p} \cdot \nabla r

    :Arguments 1:
        - material : :math:`c`
        - virtual  : :math:`q`
        - state    : :math:`p`

    :Arguments 2:
        - material    : :math:`c`
        - parameter_1 : :math:`\bar{p}`
        - parameter_2 : :math:`r`
    """
    name = 'dw_laplace'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'opt_material' : 'D, D', 'virtual' : (1, 'state'),
                   'state' : 1, 'parameter_1' : 1, 'parameter_2' : 1},
                  {'opt_material' : None}]
    modes = ('weak', 'eval')
    symbolic = {'expression': 'c * div( grad( u ) )',
                'map' : {'u' : 'state', 'c' : 'opt_material'}}

    def set_arg_types(self):
        if self.mode == 'weak':
            self.function = terms.dw_laplace

        else:
            self.function = terms.d_laplace

class PermeabilityRTerm(Term):
    r"""
    Special-purpose diffusion-like term with permeability :math:`K_{ij}` (to
    use on the right-hand side).

    :Definition:

    .. math::
        \int_{\Omega} K_{ij} \nabla_j q

    :Arguments:
        - material : :math:`K_{ij}`
        - virtual  : :math:`q`
        - index    : :math:`i`
    """
    name = 'dw_permeability_r'
    arg_types = ('material', 'virtual', 'index')

    function = staticmethod(terms.dw_permeability_r)

    def get_fargs(self, mat, virtual, index,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(virtual)

        if isinstance(index, list):
            index = index[0]

        mat = nm.ascontiguousarray(mat[...,index:index+1])

        return mat, vg

class DiffusionRTerm(Term):
    r"""
    Diffusion-like term with material parameter :math:`K_{j}` (to
    use on the right-hand side).

    :Definition:

    .. math::
        \int_{\Omega} K_{j} \nabla_j q

    :Arguments:
        - material : :math:`K_j`
        - virtual  : :math:`q`
    """
    name = 'dw_diffusion_r'
    arg_types = ('material', 'virtual')
    arg_shapes = {'material' : 'D, 1', 'virtual' : (1, None)}
    function = staticmethod(terms.dw_permeability_r)

    def get_fargs(self, mat, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):

        vg, _ = self.get_mapping(virtual)
        return mat, vg

class DiffusionCoupling(Term):
    r"""
    Diffusion copupling term with material parameter :math:`K_{j}`.

    :Definition:

    .. math::
        \int_{\Omega}  p K_{j} \nabla_j q

    :Arguments:
        - material : :math:`K_{j}`
        - virtual  : :math:`q`
        - state    : :math:`p`
    """
    name = 'dw_diffusion_coupling'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = {'material' : 'D, 1', 'virtual' : (1, 'state'),
                  'state' : 1, 'parameter_1' : 1, 'parameter_2' : 1}
    modes = ('weak0', 'weak1', 'eval')

    @staticmethod
    def d_fun(out, mat, val, grad, vg):
        out_qp = grad * mat * val

        status = vg.integrate(out, out_qp)

        return status

    @staticmethod
    def dw_fun(out, val, mat, bf, vg, fmode):

        if fmode == 0:
            status = terms.mulATB_integrate(out, vg.bfg, mat * val, vg)

        elif fmode == 1:
            status = terms.mulATB_integrate(out, bf * mat, val, vg)

        elif fmode == 2:
            status = terms.mulATB_integrate(out, vg.bfg, mat * bf, vg)

        elif fmode == 3:
            status = terms.mulATB_integrate(out, mat * bf, vg.bfg, vg)

        return status

    def get_fargs( self, mat, virtual, state,
                   mode=None, term_mode=None, diff_var=None, **kwargs):
        ap, vg = self.get_approximation(virtual)

        if mode == 'weak':

            aps, vgs = self.get_approximation(state)
            bf = aps.get_base('v', 0, self.integral)

            if diff_var is None:
                if self.mode == 'weak0':
                    val = self.get(state, 'val')
                    fmode = 0

                else:
                    val = self.get(virtual, 'grad')
                    fmode = 1

            else:
                val = nm.array([0], ndmin=4, dtype=nm.float64)
                if self.mode == 'weak0':
                    fmode = 2

                else:
                    fmode = 3

            return val, mat, bf, vg, fmode

        elif mode == 'eval':

            grad = self.get(virtual, 'grad')
            val = self.get(state, 'val')

            return mat, val, grad, vg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, mat, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, 1, 1), state.dtype

    def set_arg_types( self ):
        if self.mode[:-1] == 'weak':
            self.function = self.dw_fun

        else:
            self.function = self.d_fun

class DiffusionVelocityTerm( Term ):
    r"""
    Evaluate diffusion velocity.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:

    .. math::
        - \int_{\Omega} K_{ij} \nabla_j \bar{p}

    .. math::
        \mbox{vector for } K \from \Ical_h: - \int_{T_K} K_{ij} \nabla_j \bar{p}
        / \int_{T_K} 1

    .. math::
        - K_{ij} \nabla_j \bar{p}

    :Arguments:
        - material  : :math:`K_{ij}`
        - parameter : :math:`\bar{p}`
    """
    name = 'ev_diffusion_velocity'
    arg_types = ('material', 'parameter')
    arg_shapes = {'material' : 'D, D', 'parameter' : 1}

    @staticmethod
    def function(out, grad, mat, vg, fmode):
        dvel = dot_sequences(mat, grad)

        if fmode == 2:
            out[:] = dvel
            status = 0

        else:
            status = vg.integrate(out, dvel, fmode)

        out *= -1.0

        return status

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)
        grad = self.get(parameter, 'grad')
        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)

        return grad, mat, vg, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, dim, 1), parameter.dtype

class SurfaceFluxTerm(Term):
    r"""
    Surface flux term.

    Supports 'eval', 'el_avg' and 'el' evaluation modes.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{n} \cdot K_{ij} \nabla_j \bar{p}

    .. math::
        \mbox{vector for } K \from \Ical_h: \int_{T_K} \ul{n}
        \cdot K_{ij} \nabla_j \bar{p}\ / \int_{T_K} 1

    .. math::
        \mbox{vector for } K \from \Ical_h: \int_{T_K} \ul{n}
        \cdot K_{ij} \nabla_j \bar{p}

    :Arguments:
        - material: :math:`\ul{K}`
        - parameter:  :math:`\bar{p}`,
    """
    name = 'd_surface_flux'
    arg_types = ('material', 'parameter')
    arg_shapes = {'material' : 'D, D', 'parameter' : 1}
    integration = 'surface_extra'

    function = staticmethod(terms.d_surface_flux)

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(parameter)

        grad = self.get(parameter, 'grad')

        fmode = {'eval' : 0, 'el_avg' : 1, 'el' : 0}.get(mode, 0)

        return grad, mat, sg, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_fa, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_fa, 1, 1, 1), parameter.dtype
