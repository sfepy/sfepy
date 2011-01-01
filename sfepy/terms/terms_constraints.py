import numpy as nm

from sfepy.terms.terms import Term
from sfepy.linalg import dot_sequences

class NonPenetrationTerm(Term):
    r"""
    :Description:
    Non-penetration condition in the weak sense.

    :Definition:
    .. math::
        \int_{\Gamma} \lambda \ul{n} \cdot \ul{v} \mbox{ , }
        \int_{\Gamma} \hat\lambda \ul{n} \cdot \ul{u}

    :Arguments 1:
        virtual  : :math:`\ul{v}`,
        state    : :math:`\lambda`

    :Arguments 2:
        state    : :math:`\ul{u}`,
        virtual  : :math:`\hat\lambda`
    """
    name = 'dw_non_penetration'
    arg_types = (('virtual', 'state'),
                 ('state', 'virtual'))
    modes = ('grad', 'div')
    integration = 'surface'
    use_caches = {'state_in_surface_qp' : [['state']]}

    def __call__(self, diff_var=None, **kwargs):
        if self.mode == 'grad':
            virtual, state = self.get_args(**kwargs)
            ap, sg = self.get_approximation(virtual)

            cache = self.get_cache('state_in_surface_qp', 0)

        else:
            state, virtual = self.get_args(**kwargs)
            ap, sg = self.get_approximation(state)

            cache = self.get_cache('state_in_surface_qp', 0)

        n_fa, n_qp, dim, n_fp = ap.get_s_data_shape(self.integral,
                                                    self.region.name)
        rdim, cdim = virtual.n_components, state.n_components

        if diff_var is None:
            shape = (n_fa, 1, rdim * n_fp, 1)

        elif diff_var == self.get_arg_name('state'):
            shape = (n_fa, 1, rdim * n_fp, cdim * n_fp)

        else:
            raise StopIteration

        sd = ap.surface_data[self.region.name]

        # ap corresponds to \ul{u} field.
        bf = ap.get_base(sd.face_type, 0, integral=self.integral)
        ebf = nm.zeros((bf.shape[0], dim, n_fp * dim),
                       dtype=nm.float64)
        for ir in xrange(dim):
            ebf[:, ir, ir*n_fp:(ir+1)*n_fp] = bf[:,0,:]

        normals = sg.variable(0)

        out = nm.zeros(shape, dtype=nm.float64)

        lchunk = nm.arange(n_fa, dtype=nm.int32)

        if diff_var is None:
            vec_qp = cache('state', self, 0,
                           state=state, get_vector=self.get_vector)

            if self.mode == 'grad':
                ebf_t = nm.tile(ebf.transpose((0, 2, 1)), (n_fa, 1, 1, 1))

                nl = normals * vec_qp
                eftnl = dot_sequences(ebf_t, nl, use_rows=True)
                sg.integrate_chunk(out, eftnl, lchunk, 1)

            else:
                bf_t = nm.tile(bf.transpose((0, 2, 1)), (n_fa, 1, 1, 1))

                ntu = (normals * vec_qp).sum(axis=-2)[...,None]
                ftntu = (bf_t * ntu)

                sg.integrate_chunk(out, ftntu, lchunk, 1)

        else:
            ebf_t = nm.tile(ebf.transpose((0, 2, 1)), (n_fa, 1, 1, 1))
            bf_ = nm.tile(bf, (n_fa, 1, 1, 1))

            eftn = dot_sequences(ebf_t, normals, use_rows=True)
            eftnf = eftn * bf_

            if self.mode == 'grad':
                sg.integrate_chunk(out, eftnf, lchunk, 1)

            else:
                ftntef = nm.ascontiguousarray(eftnf.transpose((0, 1, 3, 2)))
                sg.integrate_chunk(out, ftntef, lchunk, 1)

        yield out, lchunk, 0
