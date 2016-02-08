import numpy as nm

from sfepy.terms.terms import Term, terms

class DiffusionSATerm(Term):
    r"""
    Diffusion sensitivity analysis term.

    :Definition:

    .. math::
        \int_{\Omega} \left[ (\dvg \ul{\Vcal}) K_{ij} \nabla_i q\, \nabla_j p -
        K_{ij} (\nabla_j \ul{\Vcal} \nabla q) \nabla_i p - K_{ij} \nabla_j q
        (\nabla_i \ul{\Vcal} \nabla p)\right]

    :Arguments:
        - material:    :math:`K_{ij}`
        - parameter_q: :math:`q`
        - parameter_p: :math:`p`
        - parameter_v: :math:`\ul{\Vcal}`
    """
    name = 'd_diffusion_sa'
    arg_types = ('material', 'parameter_q', 'parameter_p', 'parameter_v')
    arg_shapes = {'material' : 'D, D',
                  'parameter_q' : 1, 'parameter_p' : 1, 'parameter_v' : 'D'}

    function = staticmethod(terms.d_diffusion_sa)

    def get_fargs(self, mat, parameter_q, parameter_p, parameter_v,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter_p)

        grad_q = self.get(parameter_q, 'grad')
        grad_p = self.get(parameter_p, 'grad')
        grad_v = self.get(parameter_v, 'grad')
        div_v = self.get(parameter_v, 'div')

        return grad_q, grad_p, grad_v, div_v, mat, vg

    def get_eval_shape(self, mat, parameter_q, parameter_p, parameter_v,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter_q)

        return (n_el, 1, 1, 1), parameter_q.dtype
