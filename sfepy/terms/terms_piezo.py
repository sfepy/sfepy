import numpy as nm

from sfepy.linalg import dot_sequences
from sfepy.terms.terms import Term, terms
from sfepy.terms.terms_multilinear import ETermBase

class PiezoCouplingTerm(Term):
    r"""
    Piezoelectric coupling term. Can be evaluated.

    :Definition:

    .. math::
        \int_{\Omega} g_{kij}\ e_{ij}(\ul{v}) \nabla_k p\\
        \int_{\Omega} g_{kij}\ e_{ij}(\ul{u}) \nabla_k q

    :Arguments 1:
        - material: :math:`g_{kij}`
        - virtual/parameter_v: :math:`\ul{v}`
        - state/parameter_s: :math:`p`

    :Arguments 2:
        - material : :math:`g_{kij}`
        - state    : :math:`\ul{u}`
        - virtual  : :math:`q`
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

class SDPiezoCouplingTerm(ETermBase):
    r"""
    Sensitivity (shape derivative) of the piezoelectric coupling term.

    :Definition:

    .. math::
        \int_{\Omega} \hat{g}_{kij}\ e_{ij}(\ul{u}) \nabla_k p

    .. math::
        \hat{g}_{kij} = g_{kij}(\nabla \cdot \ul{\Vcal})
        - g_{kil}{\partial \Vcal_j \over \partial x_l}
        - g_{lij}{\partial \Vcal_k \over \partial x_l}

    :Arguments:
        - material    : :math:`g_{kij}`
        - parameter_u : :math:`\ul{u}`
        - parameter_p : :math:`p`
        - parameter_mv : :math:`\ul{\Vcal}`
    """
    name = 'ev_sd_piezo_coupling'
    arg_types = ('material', 'parameter_u', 'parameter_p', 'parameter_mv')
    arg_shapes = {'material': 'D, S', 'parameter_u': 'D', 'parameter_p': 1,
                  'parameter_mv': 'D'}
    geometries = ['2_3', '2_4', '3_4', '3_8']

    def get_function(self, mat, par_u, par_p, par_mv,
                     mode=None, term_mode=None, diff_var=None, **kwargs):
        grad_mv = self.get(par_mv, 'grad')
        div_mv = self.get(par_mv, 'div')

        nel, nqp, dim, _ = mat.shape
        remap = nm.array([0, 2, 2, 1] if dim == 2 else [0, 3, 4, 3, 1, 5, 4, 5, 2])
        mat_f = mat[:, :, :, remap].reshape((nel, nqp, dim, dim, dim))

        ## grad_mv = [v11, v21; v12, v22]
        sa_mat_f = mat_f * div_mv[..., nm.newaxis]
        sa_mat_f -= nm.einsum('qpkil,qplj->qpkij', mat_f, grad_mv,
                              optimize='greedy')
        sa_mat_f -= nm.einsum('qplij,qplk->qpkij', mat_f, grad_mv,
                              optimize='greedy')

        return self.make_function(
            'kij,i.j,0.k', (sa_mat_f, 'material'), par_u, par_p,
            mode=mode, diff_var=diff_var,
        )

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
