import numpy as nm

from sfepy.terms.terms import Term
from sfepy.linalg import dot_sequences

class NonPenetrationTerm(Term):
    r"""
    Non-penetration condition in the weak sense.

    :Definition:

    .. math::
        \int_{\Gamma} \lambda \ul{n} \cdot \ul{v} \mbox{ , }
        \int_{\Gamma} \hat\lambda \ul{n} \cdot \ul{u}

    :Arguments 1:
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\lambda`

    :Arguments 2:
        - state    : :math:`\ul{u}`
        - virtual  : :math:`\hat\lambda`
    """
    name = 'dw_non_penetration'
    arg_types = (('virtual', 'state'),
                 ('state', 'virtual'))
    modes = ('grad', 'div')
    integration = 'surface'

    @staticmethod
    def function(out, val_qp, ebf, bf, sg, diff_var, mode):
        """
        `ebf` belongs to vector variable, `bf` to scalar variable.
        """
        normals = sg.normal
        n_fa = out.shape[0]

        if diff_var is None:
            if mode == 'grad':
                ebf_t = nm.tile(ebf.transpose((0, 2, 1)), (n_fa, 1, 1, 1))

                nl = normals * val_qp
                eftnl = dot_sequences(ebf_t, nl)
                status = sg.integrate(out, eftnl, 0)

            else:
                bf_t = nm.tile(bf.transpose((0, 2, 1)), (n_fa, 1, 1, 1))

                ntu = (normals * val_qp).sum(axis=-2)[...,None]
                ftntu = (bf_t * ntu)

                status = sg.integrate(out, ftntu, 0)

        else:
            ebf_t = nm.tile(ebf.transpose((0, 2, 1)), (n_fa, 1, 1, 1))
            bf_ = nm.tile(bf, (n_fa, 1, 1, 1))

            eftn = dot_sequences(ebf_t, normals)
            eftnf = eftn * bf_

            if mode == 'grad':
                status = sg.integrate(out, eftnf, 0)

            else:
                ftntef = nm.ascontiguousarray(eftnf.transpose((0, 1, 3, 2)))
                status = sg.integrate(out, ftntef, 0)

        return status

    def get_fargs(self, vvar, svar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        if self.mode == 'grad':
            qp_var = svar

        else:
            qp_var = vvar

        val_qp = self.get(qp_var, 'val')

        vsg, _ = self.get_mapping(vvar)
        ssg, _ = self.get_mapping(svar)
        n_fa, n_qp, dim, n_fn, n_c = self.get_data_shape(vvar)

        # Expand base corresponding to \ul{u} for all dofs.
        bf = vsg.bf
        ebf = nm.zeros((bf.shape[0], dim, n_fn * dim), dtype=nm.float64)
        for ir in xrange(dim):
            ebf[:, ir, ir*n_fn:(ir+1)*n_fn] = bf[:, 0, :]

        return val_qp, ebf, ssg.bf, vsg, diff_var, self.mode
