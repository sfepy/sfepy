import numpy as nm
from numpy.lib.stride_tricks import as_strided

from sfepy.base.base import assert_
from sfepy.linalg import norm_l2_along_axis as norm
from sfepy.linalg import dot_sequences
from sfepy.fem.mappings import VolumeMapping
from sfepy.terms.terms import Term

class TLMembraneTerm(Term):
    r"""
    :Arguments:
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_tl_membrane'
    arg_types = ('material', 'virtual', 'state')
    integration = 'surface'

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        """
        Membrane term evaluation function.
        """
        mat, vv, vu = self.get_args(**kwargs)
        ap, sg = self.get_approximation(vv)
        sd = ap.surface_data[self.region.name]

        # Coordinates of element nodes.
        coors = vu.field.coors[sd.econn]

        # Local coordinate system.
        t1 = coors[:, 1, :] - coors[:, 0, :]
        t2 = coors[:, -1, :] - coors[:, 0, :]
        n = nm.cross(t1, t2)
        t2 = nm.cross(n, t1)

        t1 = t1 / norm(t1)[:, None]
        t2 = t2 / norm(t2)[:, None]
        n = n / norm(n)[:, None]

        # Coordinate transformation matrix (transposed!).
        mtx_t = nm.concatenate((t1[:, :, None],
                                t2[:, :, None],
                                n[:, :, None]), axis=2)

        # Displacements of element nodes.
        vec_u = vu.get_state_in_region(self.region, igs=[self.char_fun.ig])
        el_u = vec_u[sd.leconn]

        # Transform displacements to the local coordinate system.
        el_u_loc = dot_sequences(mtx_t, el_u)

        # Transform coordinates to the local coordinate system.
        coors_loc = dot_sequences((coors - coors[:, 0:1, :]), mtx_t)

        # Strip 'z' component (should be 0 now...).
        assert_(nm.allclose(coors_loc[:, :, -1], 0.0, rtol=1e-12, atol=1e-12))
        coors_loc = coors_loc[:, :, :-1].copy()

        # Mapping from transformed element to reference element.
        sh = coors_loc.shape
        seq_coors = coors_loc.reshape((sh[0] * sh[1], sh[2]))
        seq_conn = nm.arange(seq_coors.shape[0], dtype=nm.int32)
        seq_conn.shape = sh[:2]

        gel = vu.field.gel.surface_facet
        vm = VolumeMapping(seq_coors, seq_conn, gel=gel, order=1)

        qp = self.integral.get_qp(gel.name)
        geo = vm.get_mapping(*qp)

        # Transformed base function gradient w.r.t. material coordinates
        # in quadrature points.
        bfg = geo.variable(0)

        # Repeat el_u_loc by number of quadrature points.
        sh = list(el_u_loc.shape)
        sh.insert(1, bfg.shape[1])

        strides = list(el_u_loc.strides)
        strides.insert(1, 0)
        el_u_loc_qp = as_strided(el_u_loc, shape=sh, strides=strides)

        # Transformed (in-plane) displacement gradient with
        # shape (n_el, n_qp, dim-1 (-> j), dim (-> i)), du_i/dx_j.
        du = dot_sequences(bfg, el_u_loc_qp)

        from debug import debug; debug()

        for out, chunk in self.char_fun( chunk_size, shape ):
            lchunk = self.char_fun.get_local_chunk()

            out = None

            yield out, lchunk, 0
