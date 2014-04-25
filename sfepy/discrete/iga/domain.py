"""
Computational domain for isogeometric analysis.
"""
import os.path as op

import numpy as nm

from sfepy.base.base import Struct
from sfepy.discrete.common.domain import Domain
import sfepy.discrete.iga as iga
import sfepy.discrete.iga.io as io

class IGDomain(Domain):
    """
    Bezier extraction based NURBS domain for isogeometric analysis.
    """

    @staticmethod
    def from_file(filename):
        """
        filename : str
            The name of the IGA domain file.
        """
        (knots, degrees, cps, weights, cs, conn,
         bcps, bweights, bconn, regions) = io.read_iga_data(filename)

        nurbs = Struct(name='nurbs', knots=knots, degrees=degrees, cps=cps,
                       weights=weights, cs=cs, conn=conn)
        bmesh = Struct(name='bmesh', cps=bcps, weights=bweights, conn=bconn)

        name = op.splitext(filename)[0]
        domain = IGDomain(name, nurbs=nurbs, bmesh=bmesh, regions=regions)
        return domain

    def __init__(self, name, nurbs, bmesh, regions=None, **kwargs):
        """
        Create an IGA domain.

        Parameters
        ----------
        name : str
            The domain name.
        """
        Domain.__init__(self, name, nurbs=nurbs, bmesh=bmesh, regions=regions,
                        **kwargs)
        from sfepy.discrete.fem.geometry_element import create_geometry_elements
        from sfepy.discrete.fem import Mesh
        from sfepy.discrete.fem.extmods.cmesh import CMesh
        from sfepy.discrete.fem.utils import prepare_remap

        self.facets = iga.get_bezier_element_entities(nurbs.degrees)

        tconn = iga.get_bezier_topology(bmesh.conn, nurbs.degrees)
        itc = nm.unique(tconn)

        remap = prepare_remap(itc, bmesh.conn.max() + 1)

        ltcoors = bmesh.cps[itc]
        ltconn = remap[tconn]

        n_nod, dim = ltcoors.shape
        n_el = ltconn.shape[0]
        self.shape = Struct(n_nod=n_nod, dim=dim, tdim=0, n_el=n_el, n_gr=1)

        desc = '%d_%d' % (dim, 2**dim)
        mat_id = nm.zeros(ltconn.shape[0], dtype=nm.int32)
        self.mesh = Mesh.from_data(self.name + '_topo', ltcoors, None, [ltconn],
                                   [mat_id], [desc])

        self.cmesh = CMesh.from_mesh(self.mesh)
        gels = create_geometry_elements()
        self.cmesh.set_local_entities(gels)
        self.cmesh.setup_entities()

        self.shape.tdim = self.cmesh.tdim

        if regions is not None:
            self.vertex_set_bcs = {}
            for key, val in self.regions.iteritems():
                self.vertex_set_bcs[key] = remap[val]

        self.cell_offsets = {0 : 0}
