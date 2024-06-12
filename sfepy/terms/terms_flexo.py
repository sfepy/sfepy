"""
Flexoelectricity related terms.
"""
import numpy as nm

from sfepy.terms.terms_multilinear import ETermBase

def make_grad2strain():
    g2s = nm.array([
        [1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 1, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 1, 0],
    ], dtype=nm.float64)
    g2s = g2s.reshape((1, 1, 6, 9))

    return g2s

class MixedStrainGradElasticTerm(ETermBase):
    r"""
    Flexoelectric strain gradient elasticity term, mixed formulation.

    Additional evaluation modes:

      - `'strain'` - compute strain from the displacement gradient (state)
        variable.

    :Definition:

    .. math::
        \int_{\Omega} a_{ijklmn}\ e_{ij,k}(\ull{\delta w}) \ e_{lm,n}(\ull{w})

    :Arguments:
        - material: :math:`a_{ijklmn}`
        - virtual/parameter_1: :math:`\ull{\delta w}`
        - state/parameter_2: :math:`\ull{w}`
    """
    name = 'de_m_sg_elastic'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = {'material' : 'SD, SD', 'virtual' : ('D2', 'state'),
                  'state' : 'D2', 'parameter_1' : 'D2', 'parameter_2' : 'D2'}
    modes = ('weak', 'eval')

    def get_function(self, mat, virtual, state, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        aux = make_grad2strain()
        n_qp = mat.shape[1]
        mat = mat.reshape((-1, n_qp, 3, 6, 3, 6))

        if term_mode is None:
            return self.make_function(
                'kIlJ,Ii,Jj,i.k,j.l',
                (mat, 'mat'), (aux, 'aux1'), (aux, 'aux2'), virtual, state,
                mode=mode, diff_var=diff_var,
            )

        elif term_mode == 'strain':
            return self.make_function(
                'Ii,i', (aux, 'aux1'), state, mode=mode, diff_var=None,
            )

class MixedFlexoCouplingTerm(ETermBase):
    r"""
    Flexoelectric coupling term, mixed formulation.

    :Definition:

    .. math::
        \int_{\Omega} f_{ijkl}\ e_{jk,l}(\ull{\delta w}) \nabla_i p \\
        \int_{\Omega} f_{ijkl}\ e_{jk,l}(\ull{w}) \nabla_i q

    :Arguments 1:
        - material: :math:`f_{ijkl}`
        - virtual/parameter_t: :math:`\ull{\delta w}`
        - state/parameter_s: :math:`p`

    :Arguments 2:
        - material: :math:`f_{ijkl}`
        - state    : :math:`\ull{w}`
        - virtual  : :math:`q`
    """
    name = 'de_m_flexo_coupling'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_t', 'parameter_s'))
    arg_shapes = [{'material' : 'D, SD',
                   'virtual/dw-p' : ('D2', None), 'state/dw-p' : 1,
                   'virtual/dp-w' : (1, None), 'state/dp-w' : 'D2',
                   'parameter_t' : 'D2', 'parameter_s' : 1}]
    modes = ('dw-p', 'dp-w', 'eval')

    def get_function(self, mat, tvar, svar, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        aux = make_grad2strain()
        n_qp = mat.shape[1]
        mat = mat.reshape((-1, n_qp, 3, 3, 6))

        fun = self.make_function(
            'jkI,Ii,i.k,0.j',
            (mat, 'mat'), (aux, 'aux'), tvar, svar, diff_var=diff_var,
        )

        return fun
