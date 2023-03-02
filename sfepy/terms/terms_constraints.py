from __future__ import absolute_import
import numpy as nm

from sfepy.terms.terms import Term
from sfepy.linalg import dot_sequences

class NonPenetrationTerm(Term):
    r"""
    Non-penetration condition in the weak sense.

    :Definition:

    .. math::
        \int_{\Gamma} c \lambda \ul{n} \cdot \ul{v} \mbox{ , }
        \int_{\Gamma} c \hat\lambda \ul{n} \cdot \ul{u} \\
        \int_{\Gamma} \lambda \ul{n} \cdot \ul{v} \mbox{ , }
        \int_{\Gamma} \hat\lambda \ul{n} \cdot \ul{u}

    :Arguments 1:
        - material : :math:`c` (optional)
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\lambda`

    :Arguments 2:
        - material : :math:`c` (optional)
        - state    : :math:`\ul{u}`
        - virtual  : :math:`\hat\lambda`
    """
    name = 'dw_non_penetration'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'state', 'virtual'))
    arg_shapes = [{'opt_material' : '1, 1',
                   'virtual/grad' : ('D', None), 'state/grad' : 1,
                   'virtual/div' : (1, None), 'state/div' : 'D'},
                  {'opt_material' : None}]
    modes = ('grad', 'div')
    integration = 'facet'

    @staticmethod
    def function(out, val_qp, ebf, bf, mat, sg, diff_var, mode):
        """
        `ebf` belongs to vector variable, `bf` to scalar variable.
        """
        normals = sg.normal
        n_fa = out.shape[0]

        if diff_var is None:
            if mode == 'grad':
                ebf_t = nm.tile(ebf.transpose((0, 1, 3, 2)), (n_fa, 1, 1, 1))

                nl = normals * val_qp
                eftnl = mat * dot_sequences(ebf_t, nl)
                status = sg.integrate(out, eftnl, 0)

            else:
                bf_t = nm.tile(bf.transpose((0, 1, 3, 2)), (n_fa, 1, 1, 1))

                ntu = (normals * val_qp).sum(axis=-2)[...,None]
                ftntu = mat * (bf_t * ntu)

                status = sg.integrate(out, ftntu, 0)

        else:
            ebf_t = nm.tile(ebf.transpose((0, 1, 3, 2)), (n_fa, 1, 1, 1))
            bf_ = nm.tile(bf, (n_fa, 1, 1, 1))

            eftn = mat * dot_sequences(ebf_t, normals)
            eftnf = eftn * bf_

            if mode == 'grad':
                status = sg.integrate(out, eftnf, 0)

            else:
                ftntef = nm.ascontiguousarray(eftnf.transpose((0, 1, 3, 2)))
                status = sg.integrate(out, ftntef, 0)

        return status

    def get_fargs(self, mat, vvar, svar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        from sfepy.discrete.variables import expand_basis

        if self.mode == 'grad':
            qp_var = svar

        else:
            qp_var = vvar

        val_qp = self.get(qp_var, 'val')

        vsg, _ = self.get_mapping(vvar)
        ssg, _ = self.get_mapping(svar)
        n_fa, n_qp, dim, n_fn, n_c = self.get_data_shape(vvar)

        if mat is None:
            mat = nm.ones((1, n_qp, 1, 1), dtype=nm.float64)

        ebf = expand_basis(vsg.bf, dim)

        return val_qp, ebf, ssg.bf, mat, vsg, diff_var, self.mode

class NonPenetrationPenaltyTerm(Term):
    r"""
    Non-penetration condition in the weak sense using a penalty.

    :Definition:

    .. math::
        \int_{\Gamma} c (\ul{n} \cdot \ul{v}) (\ul{n} \cdot \ul{u})

    :Arguments:
        - material : :math:`c`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_non_penetration_p'
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material' : '1, 1',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    integration = 'facet'

    @staticmethod
    def function(out, val_qp, ebf, mat, sg, diff_var):
        normals = sg.normal
        n_fa = out.shape[0]

        ebf_t = nm.tile(ebf.transpose((0, 1, 3, 2)), (n_fa, 1, 1, 1))

        if diff_var is None:
            nu = dot_sequences(normals, val_qp, 'ATB')
            nt = dot_sequences(ebf_t, normals)
            entnu = mat * nt * nu
            status = sg.integrate(out, entnu, 0)

        else:
            nt = dot_sequences(ebf_t, normals)
            entn = mat * dot_sequences(nt, nt, 'ABT')
            status = sg.integrate(out, entn, 0)

        return status

    def get_fargs(self, mat, vvar, svar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        from sfepy.discrete.variables import expand_basis

        if diff_var is None:
            val_qp = self.get(svar, 'val')

        else:
            val_qp = nm.array([0], ndmin=4, dtype=nm.float64)

        sg, _ = self.get_mapping(vvar)
        n_fa, n_qp, dim, n_fn, n_c = self.get_data_shape(vvar)

        ebf = expand_basis(sg.bf, dim)

        return val_qp, ebf, mat, sg, diff_var
