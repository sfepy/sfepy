import numpy as nm

from sfepy.base.base import assert_
from sfepy.linalg import dot_sequences
from sfepy.terms.terms import Term, terms
from sfepy.terms.terms_dot import ScalarDotMGradScalarTerm

class DiffusionTerm(Term):
    r"""
    General diffusion term with permeability :math:`K_{ij}`. Can be
    evaluated. Can use derivatives.

    :Definition:

    .. math::
        \int_{\Omega} K_{ij} \nabla_i q \nabla_j p

    :Arguments:
        - material: :math:`K_{ij}`
        - virtual/parameter_1: :math:`q`
        - state/parameter_2: :math:`p`
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
            if self.name in ('dw_laplace', 'dw_st_pspg_p'):
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

class SDDiffusionTerm(Term):
    r"""
    Diffusion sensitivity analysis term.

    :Definition:

    .. math::
        \int_{\Omega} \hat{K}_{ij} \nabla_i q\, \nabla_j p

    .. math::
        \hat{K}_{ij} = K_{ij}\left(
            \delta_{ik}\delta_{jl} \nabla \cdot \ul{\Vcal}
          - \delta_{ik}{\partial \Vcal_j \over \partial x_l}
          - \delta_{jl}{\partial \Vcal_i \over \partial x_k}\right)

    :Arguments:
        - material:    :math:`K_{ij}`
        - parameter_q: :math:`q`
        - parameter_p: :math:`p`
        - parameter_mv: :math:`\ul{\Vcal}`
    """
    name = 'ev_sd_diffusion'
    arg_types = ('material', 'parameter_q', 'parameter_p',
                 'parameter_mv')
    arg_shapes = {'material' : 'D, D',
                  'parameter_q' : 1, 'parameter_p' : 1,
                  'parameter_mv' : 'D'}

    function = staticmethod(terms.d_sd_diffusion)

    def get_fargs(self, mat, parameter_q, parameter_p, parameter_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter_p)

        grad_q = self.get(parameter_q, 'grad')
        grad_p = self.get(parameter_p, 'grad')
        grad_mv = self.get(parameter_mv, 'grad')
        div_mv = self.get(parameter_mv, 'div')

        return grad_q, grad_p, grad_mv, div_mv, mat, vg

    def get_eval_shape(self, mat, parameter_q, parameter_p, parameter_mv,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter_q)

        return (n_el, 1, 1, 1), parameter_q.dtype

class LaplaceTerm(DiffusionTerm):
    r"""
    Laplace term with :math:`c` coefficient. Can be
    evaluated. Can use derivatives.

    :Definition:

    .. math::
        \int_{\Omega} c \nabla q \cdot \nabla p

    :Arguments 1:
        - material: :math:`c`
        - virtual/parameter_1: :math:`q`
        - state/parameter_2: :math:`p`
    """
    name = 'dw_laplace'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : (1, 'state'),
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
    function = staticmethod(terms.dw_diffusion_r)

    def get_fargs(self, mat, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):

        vg, _ = self.get_mapping(virtual)
        return mat, vg

class DiffusionCoupling(Term):
    r"""
    Diffusion copupling term with material parameter :math:`K_{j}`.

    :Definition:

    .. math::
        \int_{\Omega}  p K_{j} \nabla_j q \mbox{ , }
        \int_{\Omega}  q K_{j} \nabla_j p

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
        out_qp = val * dot_sequences(mat, grad, 'ATB')

        status = vg.integrate(out, out_qp)

        return status

    @staticmethod
    def dw_fun(out, val, mat, bf, vg, fmode):

        if fmode == 0:
            status = terms.mulAB_integrate(out, vg.bfg, mat * val, vg.cmap,
                                           mode='ATB')

        elif fmode == 1:
            status = terms.mulAB_integrate(out, bf * mat, val, vg.cmap,
                                           mode='ATB')

        elif fmode == 2:
            status = terms.mulAB_integrate(out, vg.bfg, mat * bf, vg.cmap,
                                           mode='ATB')

        elif fmode == 3:
            status = terms.mulAB_integrate(out, mat * bf, vg.bfg, vg.cmap,
                                           mode='ATB')

        return status

    def get_fargs( self, mat, virtual, state,
                   mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(virtual)

        if mode == 'weak':

            vgs, _ = self.get_mapping(state)

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

            return val, mat, vgs.bf, vg, fmode

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
        - \int_{\cal{D}} K_{ij} \nabla_j p

    :Arguments:
        - material  : :math:`K_{ij}`
        - parameter : :math:`p`
    """
    name = 'ev_diffusion_velocity'
    arg_types = ('material', 'parameter')
    arg_shapes = {'material' : 'D, D', 'parameter' : 1}
    integration = ('cell', 'facet_extra')

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

    Supports 'eval', 'el_eval' and 'el_avg' evaluation modes.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{n} \cdot K_{ij} \nabla_j p

    :Arguments:
        - material: :math:`\ul{K}`
        - parameter:  :math:`p`,
    """
    name = 'ev_surface_flux'
    arg_types = ('material', 'parameter')
    arg_shapes = {'material' : 'D, D', 'parameter' : 1}
    integration = 'facet_extra'

    function = staticmethod(terms.d_surface_flux)

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(parameter)

        grad = self.get(parameter, 'grad')

        fmode = {'eval' : 0, 'el_avg' : 1}.get(mode, 0)

        return grad, mat, sg, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_fa, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_fa, 1, 1, 1), parameter.dtype

class SurfaceFluxOperatorTerm(Term):
    r"""
    Surface flux operator term.

    :Definition:

    .. math::
        \int_{\Gamma} q \ul{n} \cdot \ull{K} \cdot \nabla p

    :Arguments:
        - material : :math:`\ull{K}`
        - virtual  : :math:`q`
        - state    : :math:`p`
    """
    name = 'dw_surface_flux'
    arg_types = ('opt_material', 'virtual', 'state')
    arg_shapes = [{'opt_material' : 'D, D', 'virtual' : (1, 'state'),
                   'state' : 1},
                  {'opt_material' : None}]
    integration = 'facet_extra'
    function = staticmethod(terms.dw_surface_flux)

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(state)
        sd = state.field.extra_data[f'sd_{self.region.name}']
        self.get_mapping(virtual) # Creates BQP for the basis.
        bf = virtual.field.eval_basis(sd.bkey, 0, self.integral)

        if mat is None:
            _, n_qp, dim, _, _ = self.get_data_shape(state)
            mat = nm.empty((1, n_qp, dim, dim), dtype=nm.float64)
            mat[..., :, :] = nm.eye(dim, dtype=nm.float64)

        if diff_var is None:
            grad = self.get(state, 'grad', integration='facet_extra')
            fmode = 0

        else:
            grad = nm.array([0], ndmin=4, dtype=nm.float64)
            fmode = 1

        return grad, mat,  bf, sg, sd.fis, fmode

class ConvectVGradSTerm(Term):
    r"""
    Scalar gradient term with convective velocity.

    :Definition:

    .. math::
        \int_{\Omega} q (\ul{u} \cdot \nabla p)

    :Arguments:
        - virtual  : :math:`q`
        - state_v  : :math:`\ul{u}`
        - state_s  : :math:`p`
    """
    name = 'dw_convect_v_grad_s'
    arg_types = ('virtual', 'state_v', 'state_s')
    arg_shapes = [{'virtual' : (1, 'state_s'), 'state_v' : 'D', 'state_s' : 1}]
    function = staticmethod(terms.dw_convect_v_grad_s)

    def get_fargs(self, virtual, state_v, state_s,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vvg, _ = self.get_mapping(state_v)
        svg, _ = self.get_mapping(state_s)

        if diff_var is None:
            grad_s = self.get(state_s, 'grad')
            val_v = self.get(state_v, 'val')
            fmode = 0

        elif diff_var == state_s.name:
            grad_s = nm.array([0], ndmin=4, dtype=nm.float64)
            val_v = self.get(state_v, 'val')
            fmode = 1

        else:
            assert_(diff_var == state_v.name)
            grad_s = self.get(state_s, 'grad')
            val_v = nm.array([0], ndmin=4, dtype=nm.float64)
            fmode = 2

        return val_v, grad_s, vvg, svg, fmode

class AdvectDivFreeTerm(ScalarDotMGradScalarTerm):
    r"""
    Advection of a scalar quantity :math:`p` with the advection velocity
    :math:`\ul{y}` given as a material parameter (a known function of space and
    time).

    The advection velocity has to be divergence-free!

    :Definition:

    .. math::
        \int_{\Omega} \nabla \cdot (\ul{y} p) q
        = \int_{\Omega} (\underbrace{(\nabla \cdot \ul{y})}_{\equiv 0}
        + \ul{y} \cdot \nabla) p) q

    :Arguments:
        - material : :math:`\ul{y}`
        - virtual  : :math:`q`
        - state    : :math:`p`
    """
    name = 'dw_advect_div_free'
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material' : 'D, 1', 'virtual' : ('1', 'state'),
                  'state' : '1'}
    mode = 'grad_state'


class NonlinearDiffusionTerm(Term):
    r"""
    The diffusion term with a scalar coefficient given by a user
    supplied function of the state variable.

    :Definition:

    .. math::
        \int_{\Omega} \nabla q \cdot \nabla p f(p)

    :Arguments:
        - fun : :math:`f(p)`
        - dfun : :math:`\partial f(p) / \partial p`
        - virtual : :math:`q`
        - state : :math:`p`
    """
    name = 'dw_nl_diffusion'
    arg_types = ('fun', 'dfun', 'virtual', 'state')
    arg_shapes = {'fun'     : lambda x: x,
                  'dfun'    : lambda x: x,
                  'virtual' : (1, 'state'),
                  'state'   : 1}

    @staticmethod
    def function(out, out_qp, geo):
        status = geo.integrate(out, out_qp)
        return status

    def get_fargs(self, fun, dfun, var1, var2,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg1, _ = self.get_mapping(var1)
        vg2, _ = self.get_mapping(var2)

        if diff_var is None:
            geo = vg1
            val_grad_qp = self.get(var2, 'grad')
            val_qp = fun(self.get(var2, 'val'))
            out_qp = dot_sequences(vg1.bfg, val_grad_qp*val_qp,'ATB')

        else:
            geo = vg1
            val_grad_qp = self.get(var2, 'grad')
            val_d_qp = dfun(self.get(var2, 'val'))
            val_qp = fun(self.get(var2, 'val'))
            out_qp = (dot_sequences(vg1.bfg, vg2.bfg*val_qp,'ATB') +
                      dot_sequences(vg1.bfg, val_grad_qp*val_d_qp,'ATB')*vg2.bf)

        return out_qp, geo
