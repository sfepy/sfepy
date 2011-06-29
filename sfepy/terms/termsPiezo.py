import numpy as nm

from sfepy.terms.terms import Term, terms

class PiezoCouplingTerm(Term):
    r"""
    :Description:
    Piezoelectric coupling term. Can be evaluated.

    :Definition:
    .. math::
        \int_{\Omega} g_{kij}\ e_{ij}(\ul{v}) \nabla_k p \mbox{ , }
        \int_{\Omega} g_{kij}\ e_{ij}(\ul{u}) \nabla_k q

    :Arguments 1:
        material : :math:`g_{kij}`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`p`

    :Arguments 2:
        material : :math:`g_{kij}`,
        state    : :math:`\ul{u}`,
        virtual  : :math:`q`

    :Arguments 3:
        material    : :math:`g_{kij}`,
        parameter_v : :math:`\ul{u}`,
        parameter_s : :math:`p`
    """
    name = 'dw_piezo_coupling'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_v', 'parameter_s'))
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

    def get_eval_shape(self, vvar, svar,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(vvar)

        return (n_el, 1, 1, 1), vvar.dtype

    def set_arg_types( self ):
        self.function = {
            'grad' : terms.dw_piezo_coupling,
            'div' : terms.dw_piezo_coupling,
            'eval' : terms.d_piezo_coupling,
        }[self.mode]
