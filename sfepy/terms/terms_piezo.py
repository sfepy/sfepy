import numpy as nm

from sfepy.linalg import dot_sequences
from sfepy.terms.terms import Term, terms

class PiezoCouplingTerm(Term):
    r"""
    Piezoelectric coupling term. Can be evaluated.

    :Definition:

    .. math::
        \int_{\Omega} g_{kij}\ e_{ij}(\ul{v}) \nabla_k p \mbox{ , }
        \int_{\Omega} g_{kij}\ e_{ij}(\ul{u}) \nabla_k q

    :Arguments 1:
        - material : :math:`g_{kij}`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`p`

    :Arguments 2:
        - material : :math:`g_{kij}`
        - state    : :math:`\ul{u}`
        - virtual  : :math:`q`

    :Arguments 3:
        - material    : :math:`g_{kij}`
        - parameter_v : :math:`\ul{u}`
        - parameter_s : :math:`p`
    """
    name = 'dw_piezo_coupling'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_v', 'parameter_s'))
    arg_shapes = {'material' : 'D, S',
                  'virtual/grad' : ('D', None), 'state/grad' : 1,
                  'virtual/div' : (1, None), 'state/div' : 'D',
                  'parameter_v' : 'D', 'parameter_s' : 1}
    modes = ('grad', 'div', 'eval')

    def get_fargs(self, mat, vvar, svar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        if self.mode == 'grad':
            qp_var, qp_name = svar, 'grad'

        else:
            qp_var, qp_name = vvar, 'cauchy_strain'

        vvg, _ = self.get_mapping(vvar)

        if mode == 'weak':
            aux = nm.array([0], ndmin=4, dtype=nm.float64)
            if diff_var is None:
                # grad or strain according to mode.
                val_qp = self.get(qp_var, qp_name)
                fmode = 0

            else:
                val_qp = aux
                fmode = 1

            if self.mode == 'grad':
                strain, grad = aux, val_qp

            else:
                strain, grad = val_qp, aux
                fmode += 2

            return strain, grad, mat, vvg, fmode

        elif mode == 'eval':
            strain = self.get(vvar, 'cauchy_strain')
            grad = self.get(svar, 'grad')

            return strain, grad, mat, vvg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, mat, vvar, svar,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(vvar)

        return (n_el, 1, 1, 1), vvar.dtype

    def set_arg_types( self ):
        self.function = {
            'grad' : terms.dw_piezo_coupling,
            'div' : terms.dw_piezo_coupling,
            'eval' : terms.d_piezo_coupling,
        }[self.mode]

class PiezoStressTerm(Term):
    r"""
    Evaluate piezoelectric stress tensor.

    It is given in the usual vector form exploiting symmetry: in 3D it has 6
    components with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in
    2D it has 3 components with the indices ordered as :math:`[11, 22, 12]`.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:

    .. math::
        \int_{\Omega} g_{kij} \nabla_k p

    :Arguments:
        - material  : :math:`g_{kij}`
        - parameter : :math:`p`
    """
    name = 'ev_piezo_stress'
    arg_types = ('material', 'parameter')
    arg_shapes = {'material' : 'D, S', 'parameter' : '1'}

    @staticmethod
    def function(out, val_qp, vg, fmode):
        if fmode == 2:
            out[:] = val_qp
            status = 0

        else:
            status = vg.integrate(out, val_qp, fmode)

        return status

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        grad = self.get(parameter, 'grad')
        val_qp = dot_sequences(mat, grad, mode='ATB')

        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)

        return val_qp, vg, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, dim * (dim + 1) // 2, 1), parameter.dtype


class PiezoStrainTerm(PiezoStressTerm):
    r"""
    Evaluate piezoelectric strain tensor.

    It is given in the usual vector form exploiting symmetry: in 3D it has 6
    components with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in
    2D it has 3 components with the indices ordered as :math:`[11, 22, 12]`.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:

    .. math::
        \int_{\Omega} g_{kij} e_{ij}(\ul{u})

    :Arguments:
        - material  : :math:`g_{kij}`
        - parameter : :math:`\ul{u}`
    """
    name = 'ev_piezo_strain'
    arg_shapes = {'material' : 'D, S', 'parameter' : 'D'}

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        strain = self.get(parameter, 'cauchy_strain')
        val_qp = dot_sequences(mat, strain, mode='AB')

        fmode = {'eval': 0, 'el_avg': 1, 'qp': 2}.get(mode, 1)

        return val_qp, vg, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, dim, 1), parameter.dtype
