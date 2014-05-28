"""
Computational domain for isogeometric analysis.
"""
import os.path as op

import numpy as nm

from sfepy.base.base import Struct
from sfepy.linalg import cycle
from sfepy.discrete.common.domain import Domain
import sfepy.discrete.iga as iga
import sfepy.discrete.iga.io as io

class NurbsPatch(Struct):
    """
    Single NURBS patch data.
    """

    def __init__(self, knots, degrees, cps,
                 weights, cs, conn):
        Struct.__init__(self, name='nurbs', knots=knots, degrees=degrees,
                        cps=cps, weights=weights, cs=cs, conn=conn)
        self.n_els = [len(ii) for ii in cs]
        self.dim = len(self.n_els)

    def _get_ref_coors_1d(self, pars, axis):
        uk = nm.unique(self.knots[axis])
        indices = nm.searchsorted(uk[1:], pars)
        ref_coors = nm.empty_like(pars)
        for ii in xrange(len(uk) - 1):
            ispan = nm.where(indices == ii)[0]
            pp = pars[ispan]
            ref_coors[ispan] = (pp - uk[ii]) / (uk[ii+1] - uk[ii])

        return uk, indices, ref_coors

    def __call__(self, u=None, v=None, w=None, field=None):
        """
        Igakit-like interface for NURBS evaluation.
        """
        pars = [u]
        if v is not None: pars += [v]
        if w is not None: pars += [w]

        indices = []
        uks = []
        rcs = []
        for ia, par in enumerate(pars):
            uk, indx, rc = self._get_ref_coors_1d(par, ia)
            indices.append(indx)
            uks.append(uk)
            rcs.append(rc)

        shape = [len(ii) for ii in pars]
        n_vals = nm.prod(shape)

        if field is None:
            out = nm.zeros((n_vals, self.dim), dtype=nm.float64)

        else:
            out = nm.zeros((n_vals, field.shape[1]), dtype=nm.float64)

        for ip, igrid in enumerate(cycle(shape)):
            iis = [indices[ii][igrid[ii]] for ii in xrange(self.dim)]
            ie = iga.get_raveled_index(iis, self.n_els)

            rc = [rcs[ii][igrid[ii]] for ii in xrange(self.dim)]

            bf, bfg, det = iga.eval_nurbs_basis_tp(rc, ie,
                                                   self.cps, self.weights,
                                                   self.degrees, self.cs,
                                                   self.conn)
            ec = self.conn[ie]

            if field is None:
                out[ip, :] = nm.dot(bf, self.cps[ec])

            else:
                out[ip, :] = nm.dot(bf, field[ec])

        return out

    def evaluate(self, field, u=None, v=None, w=None):
        """
        Igakit-like interface for NURBS evaluation.
        """
        return self(u, v, w, field)

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

        nurbs = NurbsPatch(knots, degrees, cps, weights, cs, conn)
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

        self.gel = gels[desc]

        if regions is not None:
            self.vertex_set_bcs = {}
            for key, val in self.regions.iteritems():
                self.vertex_set_bcs[key] = remap[val]

        self.cell_offsets = {0 : 0}

        self.reset_regions()
