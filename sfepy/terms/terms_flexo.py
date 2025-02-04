"""
Flexoelectricity related terms.
"""
import numpy as nm

from sfepy.mechanics.tensors import dim2sym
from sfepy.terms.terms_multilinear import ETermBase

def make_grad2strain(dim):
    if dim == 3:
        g2s = nm.array([
           # 11 12 13 21 22 23 31 32 33
            [1, 0, 0, 0, 0, 0, 0, 0, 0], # 11
            [0, 0, 0, 0, 1, 0, 0, 0, 0], # 22
            [0, 0, 0, 0, 0, 0, 0, 0, 1], # 33
            [0, 1, 0, 1, 0, 0, 0, 0, 0], # 12
            [0, 0, 1, 0, 0, 0, 1, 0, 0], # 13
            [0, 0, 0, 0, 0, 1, 0, 1, 0], # 23
        ], dtype=nm.float64)
        g2s = g2s.reshape((1, 1, 6, 9))

    elif dim == 2:
        g2s = nm.array([
           # 11 12 21 22
            [1, 0, 0, 0], # 11
            [0, 0, 0, 1], # 22
            [0, 1, 1, 0], # 12
        ], dtype=nm.float64)
        g2s = g2s.reshape((1, 1, 3, 4))

    elif dim == 1:
        g2s = nm.array([
            [1],
        ], dtype=nm.float64)
        g2s = g2s.reshape((1, 1, 1, 1))

    else:
        raise ValueError(f'space dimension must be 1, 2, or 3! (is {dim})')

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
        dim = self.region.dim
        sym = dim2sym(dim)
        aux = make_grad2strain(dim)
        n_qp = mat.shape[1]
        mat = mat.reshape((-1, n_qp, dim, sym, dim, sym))

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

        elif term_mode == 'grad_strain':
            return self.make_function(
                'Jj,j.l', (aux, 'aux2'), state, mode=mode, diff_var=None,
            )

        elif term_mode == 'double_stress':
            return self.make_function(
                'kIlJ,Ii,Jj,j.l',
                (mat, 'mat'), (aux, 'aux1'), (aux, 'aux2'), state,
                mode=mode, diff_var=None,
            )

        else:
            raise ValueError('unsupported term mode in %s! (%s)'
                             % (self.name, term_mode))

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
        dim = self.region.dim
        sym = dim2sym(dim)
        aux = make_grad2strain(dim)
        n_qp = mat.shape[1]
        mat = mat.reshape((-1, n_qp, dim, dim, sym))

        if term_mode is None:
            fun = self.make_function(
                'jkI,Ii,i.k,0.j',
                (mat, 'mat'), (aux, 'aux'), tvar, svar,
                mode=mode, diff_var=diff_var,
            )

        elif term_mode == 'electric_displacement':
            fun = self.make_function(
                'jkI,Ii,i.k',
                (mat, 'mat'), (aux, 'aux'), tvar,
                mode=mode, diff_var=diff_var,
            )

        elif term_mode == 'double_stress':
            fun = self.make_function(
                'jkI,Ii,0.j',
                (mat, 'mat'), (aux, 'aux'), svar,
                mode=mode, diff_var=diff_var,
            )

        else:
            raise ValueError('unsupported term mode in %s! (%s)'
                             % (self.name, term_mode))

        return fun

class MixedFlexoTerm(ETermBase):
    r"""
    Mixed formulation displacement gradient consistency term.

    :Definition:

    .. math::
        \int_{\Omega} v_{i,j} a_{ij} \\
        \int_{\Omega} u_{i,j} \delta a_{ij}

    :Arguments 1:
        - virtual/parameter_v: :math:`\ul{v}`
        - state/parameter_t: :math:`\ull{a}`

    :Arguments 2:
        - state    : :math:`\ul{u}`
        - virtual  : :math:`\ull{\delta a}`
    """
    name = 'de_m_flexo'
    arg_types = (('virtual', 'state'),
                 ('state', 'virtual'),
                 ('parameter_v', 'parameter_t'))
    arg_shapes = [{'virtual/du-a' : ('D', None), 'state/du-a' : 'D2',
                   'virtual/da-u' : ('D2', None), 'state/da-u' : 'D',
                   'parameter_v' : 'D', 'parameter_t' : 'D2'},
                  {'opt_material' : None}]
    modes = ('du-a', 'da-u', 'eval')

    def get_function(self, vvar, tvar, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        fun = self.make_function(
            'v(i.j)->I,I', vvar, tvar, mode=mode, diff_var=diff_var,
        )

        return fun
