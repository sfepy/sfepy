import numpy as nm

from sfepy.base.base import assert_, get_default, Struct

class FESurface(Struct):
    """Description of a surface of a finite element domain."""

    def __init__(self, name, region, efaces, volume_econn, ig):
        """nodes[leconn] == econn"""
        """nodes are sorted by node number -> same order as region.vertices"""
        self.name = get_default(name, 'surface_data_%s' % region.name)

        face_indices = region.fis[ig]
        
        faces = efaces[face_indices[:,1]]
        if faces.size == 0:
            raise ValueError('region with group with no faces! (%s)'
                             % region.name)

        try:
            ee = volume_econn[face_indices[:,0]]

        except:
            raise ValueError('missing region face indices! (%s)'
                             % region.name)

        econn = nm.empty(faces.shape, dtype=nm.int32)
        for ir, face in enumerate( faces ):
#            print ir, face, ee[ir,face]
            econn[ir] = ee[ir,face]

        ef = econn.flat
        nodes = nm.unique(ef)

        aux = -nm.ones((nm.max( ef ) + 1,), dtype=nm.int32)
        aux[nodes] = nm.arange(len(nodes), dtype=nm.int32)
        leconn = aux[econn].copy()
        assert_(nm.alltrue(nodes[leconn] == econn))
        
        n_fa, n_fp = face_indices.shape[0], faces.shape[1]
        face_type = 's%d' % n_fp

        # Store bkey in SurfaceData, so that base function can be
        # queried later.
        bkey = 'b%s' % face_type[1:]

        self.ig = ig
        self.econn = econn
        self.fis = face_indices
        self.n_fa, self.n_fp = n_fa, n_fp
        self.nodes = nodes
        self.leconn = leconn
        self.face_type = face_type
        self.bkey = bkey
        self.meconn = self.mleconn = None

    def setup_mirror_connectivity(self, region):
        """
        Setup mirror surface connectivity required to integrate over a
        mirror region.

        Notes
        -----
        Works only for linear 2D problems so far.
        """
        dim = region.domain.shape.dim

        facets = region.domain.get_facets(force_faces=True)[1]
        mregion, ig_map, ig_map_i = region.get_mirror_region()

        mig = ig_map_i[self.ig]

        off = nm.cumsum(nm.r_[0, facets.n_obj])
        ii = region.get_surface_entities(self.ig) - off[self.ig]
        mii = mregion.get_surface_entities(mig) - off[mig]

        ori = facets.oris[self.ig][ii]
        mori = facets.oris[mig][mii]

        n_fp = facets.n_fps[self.ig]

        if dim == 2:
            assert_(((ori + mori) == 1).all())

            if self.n_fp > n_fp:
                raise NotImplementedError

            self.meconn = self.econn[:, ::-1].copy()
            self.mleconn = self.leconn[:, ::-1].copy()

        else:
            raise NotImplementedError

    def get_connectivity(self, local=False, is_trace=False):
        """
        Return the surface element connectivity.

        Parameters
        ----------
        local : bool
            If True, return local connectivity w.r.t. surface nodes,
            otherwise return global connectivity w.r.t. all mesh nodes.
        is_trace : bool
            If True, return mirror connectivity according to `local`.
        """
        if not is_trace:
            if local:
                return self.leconn

            else:
                return self.econn

        else:
            if local:
                return self.mleconn

            else:
                return self.meconn
