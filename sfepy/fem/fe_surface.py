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
        """

        facets = region.domain.get_facets(force_faces=True)[1]
        mregion, ig_map, ig_map_i = region.get_mirror_region()
        mig = ig_map_i[self.ig]

        grp_int = facets.group_interfaces

        ofaces = region.get_surface_entities(self.ig)
        mfaces = mregion.get_surface_entities(mig)

        ooris = facets.oris[self.ig]
        moris = facets.oris[mig]

        oori_maps = facets.ori_maps[self.ig]
        mori_maps = facets.ori_maps[mig]

        self.meconn = nm.zeros_like(self.econn)
        self.mleconn = nm.zeros_like(self.leconn)

        for idx, ifc in enumerate(ofaces):
            aux = nm.where(grp_int[:,0] == ifc)[0]
            if len(aux) == 0:
                aux = nm.where(grp_int[:,1] == ifc)[0]
                mifc = grp_int[aux[0],0]

            else:
                mifc = grp_int[aux[0],1]

            midx = nm.where(mfaces == mifc)[0][0]

            ogslice = facets.indx[self.ig]
            mgslice = facets.indx[mig]
            omap = oori_maps[ooris[ifc - ogslice.start]][1]
            mmap = mori_maps[moris[mifc - mgslice.start]][1]
            self.meconn[idx,omap] = self.econn[midx,mmap]
            self.mleconn[idx,omap] = self.leconn[midx,mmap]

        nods = region.domain.mesh.coors

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
