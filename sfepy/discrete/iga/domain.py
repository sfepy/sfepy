"""
Computational domain for isogeometric analysis.
"""
from __future__ import absolute_import
import os.path as op

import numpy as nm

from sfepy.base.base import assert_, Struct
from sfepy.discrete.common.domain import Domain
from sfepy.discrete.iga import iga
from sfepy.discrete.iga import io
from sfepy.discrete.iga.extmods.igac import eval_in_tp_coors
import six
from six.moves import range

class NurbsPatch(Struct):
    """
    Single NURBS patch data.
    """

    def __init__(self, knots, degrees, cps,
                 weights, cs, conn):
        degrees = nm.asarray(degrees, dtype=nm.int32)
        cs = [nm.asarray(cc, dtype=nm.float64) for cc in cs]
        if cs[0].ndim == 3:
            cs = [nm.ascontiguousarray(cc[:, None, ...]) for cc in cs]

        Struct.__init__(self, name='nurbs', knots=knots, degrees=degrees,
                        cps=cps, weights=weights, cs=cs, conn=conn)
        self.n_els = [len(ii) for ii in cs]
        self.dim = len(self.n_els)

    def _get_ref_coors_1d(self, pars, axis):
        uk = nm.unique(self.knots[axis])
        indices = nm.searchsorted(uk[1:], pars)
        ref_coors = nm.empty_like(pars)
        for ii in range(len(uk) - 1):
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
        rcs = []
        for ia, par in enumerate(pars):
            uk, indx, rc = self._get_ref_coors_1d(par, ia)
            indices.append(indx.astype(nm.uint32))
            rcs.append(rc)

        out = eval_in_tp_coors(field, indices,
                               rcs, self.cps, self.weights,
                               self.degrees,
                               self.cs, self.conn)

        return out

    def evaluate(self, field, u=None, v=None, w=None):
        """
        Igakit-like interface for NURBS evaluation.
        """
        return self(u, v, w, field)

    def _to_igakit(self):
        import igakit.cad as cad

        n_efuns = self.degrees + 1
        nks = nm.array([len(ii) for ii in self.knots])
        shape = tuple(nks - n_efuns)

        cps = self.cps.reshape(shape + (-1,))
        weights = self.weights.reshape(shape)

        return cad.NURBS(self.knots, cps, weights=weights)

    def _from_igakit(self, inurbs):
        cs = iga.compute_bezier_extraction(inurbs.knots, inurbs.degree)
        n_els = [len(ii) for ii in cs]
        conn, bconn = iga.create_connectivity(n_els, inurbs.knots,
                                              inurbs.degree)

        cps = inurbs.points[..., :self.dim].copy()
        cps = cps.reshape((-1, self.dim))

        return NurbsPatch(inurbs.knots, inurbs.degree, cps,
                          inurbs.weights.ravel(), cs, conn)

    def elevate(self, times=0):
        """
        Elevate the patch degrees several `times` by one.

        Returns
        -------
        nurbs : NurbsPatch instance
            Either `self` if `times` is zero, or a new instance.
        """
        if times == 0: return self

        aux = self._to_igakit()
        for ia in range(self.dim):
            aux.elevate(ia, times)
            assert_(nm.isfinite(aux.points).all(),
                    'igakit degree elevation failed for axis %d!' % ia)

        return self._from_igakit(aux)

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
        data = io.read_iga_data(filename)
        name = op.splitext(filename)[0]
        return IGDomain.from_data(*(data + (name,)))

    @staticmethod
    def read_domain_from_hdf5(fd, group):
        """
        Create a domain from the given hdf5 data group.

        fd: tables.File
            HDF5 file handle to read the mesh from.
        group: tables.group.Group
            HDF5 data group (of file fd) to read the mesh from.
        """
        data = io.read_iga_data(fd, group)
        return IGDomain.from_data(*data)

    def write_domain_to_hdf5(self, fd, group):
        """
        Save the domain to a hdf5 file.

        fd: tables.File
            HDF5 file handle to write the mesh to.
        group: tables.group.Group
            HDF5 data group (of file fd) to write the mesh to.
        """
        io.write_iga_data(fd, group, *(self._get_io_data() + (self.name,)))

    def _get_io_data(self):
        """
        Return the data describing the domain for storing the domain in a hdf5
        file.

        TODO - data for regions recreating
        """
        knots = self.nurbs.knots
        degrees = self.nurbs.degrees
        cps = self.nurbs.cps
        weights = self.nurbs.weights
        cs = self.nurbs.cs
        conn = self.nurbs.conn
        bcps = self.bmesh.cps
        bweights = self.bmesh.weights
        bconn = self.bmesh.conn

        return (knots, degrees, cps, weights, cs, conn,
                bcps, bweights, bconn, self.vertex_set_bcs)

    @staticmethod
    def from_data(knots, degrees, cps, weights, cs, conn,
                  bcps, bweights, bconn, regions, name='iga_domain_from_data'):
        """
        Create the IGA domain from the given data.
        """

        nurbs = NurbsPatch(knots, degrees, cps, weights, cs, conn)
        bmesh = Struct(name='bmesh', cps=bcps, weights=bweights, conn=bconn)

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
        from sfepy.discrete.fem.utils import prepare_remap

        tconn = iga.get_bezier_topology(bmesh.conn, nurbs.degrees)
        itc = nm.unique(tconn)

        remap = prepare_remap(itc, bmesh.conn.max() + 1)

        ltcoors = bmesh.cps[itc]
        ltconn = remap[tconn]

        n_nod, dim = ltcoors.shape
        n_el = ltconn.shape[0]
        self.shape = Struct(n_nod=n_nod, dim=dim, tdim=0, n_el=n_el)

        desc = '%d_%d' % (dim, bmesh.conn.shape[1])
        mat_id = nm.zeros(bmesh.conn.shape[0], dtype=nm.int32)
        eval_mesh = Mesh.from_data(self.name + '_eval', nurbs.cps, None,
                                   [nurbs.conn], [mat_id], [desc])
        self.eval_mesh = eval_mesh

        desc = '%d_%d' % (dim, 2**dim)
        mat_id = nm.zeros(ltconn.shape[0], dtype=nm.int32)
        self.mesh = Mesh.from_data(self.name + '_topo', ltcoors, None, [ltconn],
                                   [mat_id], [desc])
        self.cmesh = self.mesh.cmesh
        self.cmesh_tdim = self.mesh.cmesh_tdim
        gels = create_geometry_elements()
        self.cmesh.set_local_entities(gels)
        self.cmesh.setup_entities()

        self.shape.tdim = self.cmesh.tdim

        self.gel = gels[desc]

        if regions is not None:
            self.vertex_set_bcs = {}
            for key, val in six.iteritems(self.regions):
                self.vertex_set_bcs[key] = remap[val]

        self.reset_regions()
