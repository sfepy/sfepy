import numpy as nm

from sfepy.linalg import dot_sequences
from sfepy.terms.terms import Term, terms
from sfepy.terms.terms_th import THTerm, ETHTerm
from sfepy.terms.terms_elastic import CauchyStressTerm

class BiotTerm(Term):
    r"""
    Biot coupling term with :math:`\alpha_{ij}` given in:

    * vector form exploiting symmetry - in 3D it has the
      indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in 2D it has
      the indices ordered as :math:`[11, 22, 12]`,

    * matrix form - non-symmetric coupling parameter.

    Corresponds to weak forms of Biot gradient and divergence terms.
    Can be evaluated. Can use derivatives.

    :Definition:

    .. math::
        \int_{\Omega}  p\ \alpha_{ij} e_{ij}(\ul{v}) \mbox{ , } \int_{\Omega}
        q\ \alpha_{ij} e_{ij}(\ul{u})

    :Arguments 1:
        - material : :math:`\alpha_{ij}`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`p`

    :Arguments 2:
        - material : :math:`\alpha_{ij}`
        - state    : :math:`\ul{u}`
        - virtual  : :math:`q`

    :Arguments 3:
        - material    : :math:`\alpha_{ij}`
        - parameter_v : :math:`\ul{u}`
        - parameter_s : :math:`p`
    """
    name = 'dw_biot'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_v', 'parameter_s'))
    arg_shapes = [{'material' : 'S, 1',
                  'virtual/grad' : ('D', None), 'state/grad' : 1,
                  'virtual/div' : (1, None), 'state/div' : 'D',
                  'parameter_v' : 'D', 'parameter_s' : 1},
                  {'material' : 'D, D'}]

    modes = ('grad', 'div', 'eval')

    def get_fargs(self, mat, vvar, svar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):

        sym_mode = False if mat.shape[-2] == mat.shape[-1] > 1 else True
        if not sym_mode:
            sh = mat.shape
            # the gradient given by 'self.get' is transposed
            mat = nm.swapaxes(mat, 2, 3)
            mat = mat.reshape(sh[:2] + (sh[2]**2, 1))

        if self.mode == 'grad':
            qp_var, qp_name = svar, 'val'

        else:
            if sym_mode:
                qp_var, qp_name = vvar, 'cauchy_strain'
            else:
                qp_var, qp_name = vvar, 'grad'

        if mode == 'weak':
            vvg, _ = self.get_mapping(vvar)
            svg, _ = self.get_mapping(svar)

            if diff_var is None:
                val_qp = self.get(qp_var, qp_name)
                if qp_name == 'grad':
                    sh = val_qp.shape
                    val_qp = val_qp.reshape(sh[:2] + (sh[2]**2, 1))

                fmode = 0

            else:
                val_qp = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return 1.0, val_qp, mat, svg, vvg, fmode

        elif mode == 'eval':
            vvg, _ = self.get_mapping(vvar)

            if sym_mode:
                strain = self.get(vvar, 'cauchy_strain')
            else:
                strain = self.get(vvar, 'grad')
                sh = strain.shape
                strain = strain.reshape(sh[:2] + (sh[2]**2, 1))

            pval = self.get(svar, 'val')

            return 1.0, pval, strain, mat, vvg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, mat, vvar, svar,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(vvar)

        return (n_el, 1, 1, 1), vvar.dtype

    def set_arg_types(self):
        self.function = {
            'grad' : terms.dw_biot_grad,
            'div' : terms.dw_biot_div,
            'eval' : terms.d_biot_div,
        }[self.mode]

class BiotStressTerm(CauchyStressTerm):
    r"""
    Evaluate Biot stress tensor.

    It is given in the usual vector form exploiting symmetry: in 3D it has 6
    components with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in
    2D it has 3 components with the indices ordered as :math:`[11, 22, 12]`.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:

    .. math::
        - \int_{\Omega} \alpha_{ij} p

    :Arguments:
        - material  : :math:`\alpha_{ij}`
        - parameter : :math:`p`
    """
    name = 'ev_biot_stress'
    arg_types = ('material', 'parameter')
    arg_shapes = {'material' : 'S, 1', 'parameter' : 1}
    integration = 'cell'

    @staticmethod
    def function(out, val_qp, mat, vg, fmode):
        if fmode == 2:
            out[:] = dot_sequences(mat, val_qp)
            status = 0

        else:

            status = terms.de_cauchy_stress(out, val_qp, mat, vg.cmap, fmode)

        out *= -1.0

        return status

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        val_qp = self.get(parameter, 'val')

        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)

        return val_qp, mat, vg, fmode

class BiotTHTerm(BiotTerm, THTerm):
    r"""
    Fading memory Biot term. Can use derivatives.

    :Definition:

    .. math::
        \begin{array}{l}
        \int_{\Omega} \left [\int_0^t \alpha_{ij}(t-\tau)\,p(\tau)) \difd{\tau}
        \right]\,e_{ij}(\ul{v}) \mbox{ ,} \\
        \int_{\Omega} \left [\int_0^t
        \alpha_{ij}(t-\tau) e_{kl}(\ul{u}(\tau)) \difd{\tau} \right] q
        \end{array}

    :Arguments 1:
        - ts       : :class:`TimeStepper` instance
        - material : :math:`\alpha_{ij}(\tau)`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`p`

    :Arguments 2:
        - ts       : :class:`TimeStepper` instance
        - material : :math:`\alpha_{ij}(\tau)`
        - state    : :math:`\ul{u}`
        - virtual  : :math:`q`
    """
    name = 'dw_biot_th'
    arg_types = (('ts', 'material', 'virtual', 'state'),
                 ('ts', 'material', 'state', 'virtual'))
    arg_shapes = {'material' : '.: N, S, 1',
                  'virtual/grad' : ('D', None), 'state/grad' : 1,
                  'virtual/div' : (1, None), 'state/div' : 'D'}
    modes = ('grad', 'div')

    def get_fargs(self, ts, mats, vvar, svar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        if self.mode == 'grad':
            qp_var, qp_name = svar, 'val'

        else:
            qp_var, qp_name = vvar, 'cauchy_strain'

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(svar)

        if mode == 'weak':
            vvg, _ = self.get_mapping(vvar)
            svg, _ = self.get_mapping(svar)

            if diff_var is None:
                def iter_kernel():
                    for ii, mat in enumerate(mats):
                        val_qp = self.get(qp_var, qp_name, step=-ii)
                        mat = nm.tile(mat, (n_el, n_qp, 1, 1))
                        yield ii, (ts.dt, val_qp, mat, svg, vvg, 0)
                fargs = iter_kernel

            else:
                val_qp = nm.array([0], ndmin=4, dtype=nm.float64)
                mat = nm.tile(mats[0], (n_el, n_qp, 1, 1))
                fargs = ts.dt, val_qp, mat, svg, vvg, 1

            return fargs

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

class BiotETHTerm(BiotTerm, ETHTerm):
    r"""
    This term has the same definition as dw_biot_th, but assumes an
    exponential approximation of the convolution kernel resulting in much
    higher efficiency. Can use derivatives.

    :Definition:

    .. math::
        \begin{array}{l}
        \int_{\Omega} \left [\int_0^t \alpha_{ij}(t-\tau)\,p(\tau)) \difd{\tau}
        \right]\,e_{ij}(\ul{v}) \mbox{ ,} \\
        \int_{\Omega} \left [\int_0^t
        \alpha_{ij}(t-\tau) e_{kl}(\ul{u}(\tau)) \difd{\tau} \right] q
        \end{array}

    :Arguments 1:
        - ts         : :class:`TimeStepper` instance
        - material_0 : :math:`\alpha_{ij}(0)`
        - material_1 : :math:`\exp(-\lambda \Delta t)` (decay at :math:`t_1`)
        - virtual    : :math:`\ul{v}`
        - state      : :math:`p`

    :Arguments 2:
        - ts         : :class:`TimeStepper` instance
        - material_0 : :math:`\alpha_{ij}(0)`
        - material_1 : :math:`\exp(-\lambda \Delta t)` (decay at :math:`t_1`)
        - state      : :math:`\ul{u}`
        - virtual    : :math:`q`
    """
    name = 'dw_biot_eth'
    arg_types = (('ts', 'material_0', 'material_1', 'virtual', 'state'),
                 ('ts', 'material_0', 'material_1', 'state', 'virtual'))
    arg_shapes = {'material_0' : 'S, 1', 'material_1' : '1, 1',
                  'virtual/grad' : ('D', None), 'state/grad' : 1,
                  'virtual/div' : (1, None), 'state/div' : 'D'}
    modes = ('grad', 'div')

    def get_fargs(self, ts, mat0, mat1, vvar, svar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        if self.mode == 'grad':
            qp_var, qp_name, iv = svar, 'val', 4

        else:
            qp_var, qp_name, iv = vvar, 'cauchy_strain', 3

        if mode == 'weak':
            vvg, _, key = self.get_mapping(vvar, return_key=True)
            svg, _ = self.get_mapping(svar)

            if diff_var is None:
                val_qp = self.get(qp_var, qp_name)

                key += tuple(self.arg_names[ii] for ii in [1, 2, iv])
                data = self.get_eth_data(key, qp_var, mat1, val_qp)

                val = data.history + data.values
                fargs = (ts.dt, val, mat0, svg, vvg, 0)

            else:
                val_qp = nm.array([0], ndmin=4, dtype=nm.float64)
                fargs = (ts.dt, val_qp, mat0, svg, vvg, 1)

            return fargs

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))
