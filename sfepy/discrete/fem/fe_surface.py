import numpy as nm

from sfepy.base.base import get_default, Struct
from sfepy.discrete.fem.facets import build_orientation_map
from sfepy.discrete.fem.utils import prepare_remap

class FESurface(Struct):
    """Description of a surface of a finite element domain."""

    def __init__(self, name, region, efaces, volume_econn, volume_region=None):
        """nodes[leconn] == econn"""
        """nodes are sorted by node number -> same order as region.vertices"""
        self.name = get_default(name, 'surface_data_%s' % region.name)

        face_indices = region.get_facet_indices()

        faces = efaces[face_indices[:,1]]
        if faces.size == 0 and not region.is_empty:
            raise ValueError('region with no faces! (%s)' % region.name)

        if volume_region is None:
            ii = face_indices[:, 0]

        else:
            ii = volume_region.get_cell_indices(face_indices[:, 0])

        try:
            ee = volume_econn[ii]

        except:
            raise ValueError('missing region face indices! (%s)'
                             % region.name)

        econn = nm.empty(faces.shape, dtype=nm.int32)
        for ir, face in enumerate( faces ):
            econn[ir] = ee[ir,face]

        nodes = nm.unique(econn)
        if len(nodes):
            remap = prepare_remap(nodes, nodes.max() + 1)
            leconn = remap[econn].copy()

        else:
            leconn = econn.copy()

        n_fa, n_fp = face_indices.shape[0], faces.shape[1]
        face_type = 's%d' % n_fp

        # Store bkey in SurfaceData, so that base function can be
        # queried later.
        bkey = 'b%s' % face_type[1:]

        self.econn = econn
        self.fis = nm.ascontiguousarray(face_indices.astype(nm.int32))
        self.n_fa, self.n_fp = n_fa, n_fp
        self.nodes = nodes
        self.leconn = leconn
        self.face_type = face_type
        self.bkey = bkey
        self.meconn = self.mleconn = None

        if self.n_fp <= 4:
            self.ori_map, _, _ = build_orientation_map(self.n_fp)

        else:
            self.ori_map = None

    def setup_mirror_connectivity(self, region, mirror_name):
        """
        Setup mirror surface connectivity required to integrate over a
        mirror region.

        1. Get orientation of the faces:
           a) for own elements -> ooris
           b) for mirror elements -> moris

        2. orientation -> permutation.
        """
        mregion = region.get_mirror_region(mirror_name)

        oo = self.ori_map
        ori_map = nm.zeros((nm.max(list(oo.keys())) + 1, self.n_fp), dtype=nm.int32)
        ori_map[list(oo.keys())] = nm.array([ii[1] for ii in oo.values()])

        conn = region.domain.cmesh.get_conn_as_graph(region.dim, region.dim - 1)
        oris = region.domain.cmesh.facet_oris

        econn = self.econn
        ofis = region.get_facet_indices()
        if mirror_name in region.mirror_maps\
            and region.mirror_maps[mirror_name] is not None:
            mirror_map = region.mirror_maps[mirror_name]
            ofis = ofis[mirror_map]
            econn = econn[mirror_map]
        ooris = oris[conn.indptr[ofis[:, 0]] + ofis[:, 1]]
        mfis = mregion.get_facet_indices()
        moris = oris[conn.indptr[mfis[:, 0]] + mfis[:, 1]]

        omap = ori_map[ooris]
        mmap = ori_map[moris]

        n_el, n_ep = self.econn.shape
        ii = nm.repeat(nm.arange(n_el)[:, None], n_ep, 1)
        self.meconn = nm.empty_like(self.econn)
        self.meconn[ii, omap] = econn[ii, mmap]

        nodes = nm.unique(self.meconn)
        remap = prepare_remap(nodes, nodes.max() + 1)
        self.mleconn = remap[self.meconn].copy()

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
