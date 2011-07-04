import numpy as nm

from sfepy.terms.terms import Term, terms

##
# 22.08.2006, c
class LinearTractionTerm( Term ):
    r"""
    :Description:
    Linear traction forces (weak form), where, depending on dimension of
    'material' argument, :math:`\ull{\sigma} \cdot \ul{n}` is
    :math:`\bar{p} \ull{I} \cdot \ul{n}` for a given scalar pressure,
    :math:`\ul{f}` for a traction vector, and itself for a stress tensor.

    :Definition:
    .. math::
        \int_{\Gamma} \ul{v} \cdot \ull{\sigma} \cdot \ul{n}

    :Arguments:
        material : :math:`\ull{\sigma}`,
        virtual  : :math:`\ul{v}`
    """
    name = 'dw_surface_ltr'
    arg_types = ('material', 'virtual')
    integration = 'surface'

    function = staticmethod(terms.dw_surface_ltr)

    def get_fargs(self, traction, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)

        return sg.bf, traction, sg

class SurfaceJumpTerm(Term):
    r"""
    :Description:
    Interface jump condition.

    :Definition:
    .. math::
        \int_{\Gamma} q (p_1 - p_2 - c)

    :Arguments:
        material : :math:`c`,
        virtual  : :math:`q`,
        state_1  : :math:`p_1`,
        state_2  : :math:`p_2`
    """
    name = 'dw_jump'
    arg_types = ('material', 'virtual', 'state_1', 'state_2')
    integration = 'surface'

    @staticmethod
    def function(out, jump, bf1, bf2, sg, fmode):
        bf_t = nm.tile(sg.bf.transpose((0, 2, 1)), (out.shape[0], 1, 1, 1))

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

            if diff_var == self.get_arg_name('state_1'):
                fmode = 1

            else:
                fmode = 2

        return jump, sg1.bf, sg2.bf, sg, fmode
