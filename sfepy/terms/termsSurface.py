import numpy as nm

from sfepy.terms.terms import Term, terms
from sfepy.linalg import dot_sequences

##
# 22.08.2006, c
class LinearTractionTerm( Term ):
    r"""
    Linear traction forces (weak form), where, depending on dimension of
    'material' argument, :math:`\ull{\sigma} \cdot \ul{n}` is
    :math:`\bar{p} \ull{I} \cdot \ul{n}` for a given scalar pressure,
    :math:`\ul{f}` for a traction vector, and itself for a stress tensor.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{v} \cdot \ull{\sigma} \cdot \ul{n}

    :Arguments:
        - material : :math:`\ull{\sigma}`
        - virtual  : :math:`\ul{v}`
    """
    name = 'dw_surface_ltr'
    arg_types = ('material', 'virtual')
    integration = 'surface'

    function = staticmethod(terms.dw_surface_ltr)

    def get_fargs(self, traction, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)

        return traction, sg

class SufaceNormalDotTerm(Term):
    r"""
    "Scalar traction" term, (weak form).

    :Definition:

    .. math::
        \int_{\Gamma} q \ul{c} \cdot \ul{n}

    :Arguments:
        - material : :math:`\ul{c}`
        - virtual  : :math:`q`
    """
    name = 'dw_surface_ndot'
    arg_types = (('material', 'virtual'),
                 ('material', 'parameter'))
    modes = ('weak', 'eval')
    integration = 'surface'

    @staticmethod
    def dw_fun(out, material, bf, sg):
        bf_t = nm.tile(bf.transpose((0, 1, 3, 2)), (out.shape[0], 1, 1, 1))
        bf_t = nm.ascontiguousarray(bf_t)
        aux = dot_sequences(material, sg.normal)
        status = sg.integrate(out, bf_t * aux)
        return status

    @staticmethod
    def d_fun(out, material, val, sg):
        aux = dot_sequences(material, sg.normal)
        status = sg.integrate(out, val * aux)
        return status

    def get_fargs(self, mat, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)

        if mode == 'weak':
            return mat, sg.bf, sg

        elif mode == 'eval':
            val = self.get(virtual, 'val')
            return mat, val, sg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, mat, virtual,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(virtual)

        return (n_el, 1, 1, 1), virtual.dtype

    def set_arg_types( self ):
        if self.mode == 'weak':
            self.function = self.dw_fun

        else:
            self.function = self.d_fun

class SDSufaceNormalDotTerm(Term):
    r"""
    Sensitivity of scalar traction.

    :Definition:

    .. math::
        \int_{\Gamma} p \ul{c} \cdot \ul{n} \nabla \cdot \ul{\Vcal}

    :Arguments:
        - material : :math:`\ul{c}`
        - parameter : :math:`p`
        - parameter_mesh_velocity : :math:`\ul{\Vcal}`
    """
    name = 'd_sd_surface_ndot'
    arg_types = ('material', 'parameter', 'parameter_mesh_velocity')
    integration = 'surface'

    @staticmethod
    def d_fun(out, material, val_p, div_mv, sg):
        aux = dot_sequences(material, sg.normal)
        aux2 = dot_sequences(aux, div_mv)
        status = sg.integrate(out, val_p * aux2)
        return status

    def get_fargs(self, mat, par, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(par)

        val_p = self.get(par, 'val')
        div_mv = self.get(par_mv, 'div')

        return mat, val_p, div_mv, sg

    def get_eval_shape(self, mat, virtual,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(virtual)

        return (n_el, 1, 1, 1), virtual.dtype

class SurfaceJumpTerm(Term):
    r"""
    Interface jump condition.

    :Definition:

    .. math::
        \int_{\Gamma} q (p_1 - p_2 - c)

    :Arguments:
        - material : :math:`c`
        - virtual  : :math:`q`
        - state_1  : :math:`p_1`
        - state_2  : :math:`p_2`
    """
    name = 'dw_jump'
    arg_types = ('material', 'virtual', 'state_1', 'state_2')
    integration = 'surface'

    @staticmethod
    def function(out, jump, bf1, bf2, sg, fmode):
        bf_t = nm.tile(sg.bf.transpose((0, 1, 3, 2)), (out.shape[0], 1, 1, 1))

        if fmode == 0:
            vec = bf_t * jump

        elif fmode == 1:
            vec = bf_t * bf1

        else:
            vec = - bf_t * bf2

        status = sg.integrate(out, vec)

        return status

    def get_fargs(self, coef, virtual, state1, state2,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)
        sg1, _ = self.get_mapping(state1)
        sg2, _ = self.get_mapping(state2)

        if diff_var is None:
            val1 = self.get(state1, 'val')
            val2 = self.get(state2, 'val')
            jump = val1 - val2 - coef
            fmode = 0

        else:
            jump = None

            if diff_var == state1.name:
                fmode = 1

            else:
                fmode = 2

        return jump, sg1.bf, sg2.bf, sg, fmode
