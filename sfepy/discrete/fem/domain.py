"""
Computational domain, consisting of the mesh and regions.
"""
import time

import numpy as nm

from sfepy.base.base import output, Struct
from geometry_element import GeometryElement
from sfepy.discrete.common.domain import Domain
from sfepy.discrete.fem.refine import refine_2_3, refine_2_4, refine_3_4, refine_3_8
from sfepy.discrete.fem.fe_surface import FESurface
from sfepy.discrete.fem.mesh import make_inverse_connectivity
import fea

class FEDomain(Domain):
    """
    Domain is divided into groups, whose purpose is to have homogeneous
    data shapes.
    """

    def __init__(self, name, mesh, verbose=False, **kwargs):
        """Create a Domain.

        Parameters
        ----------
        name : str
            Object name.
        mesh : Mesh
            A mesh defining the domain.
        """
        Domain.__init__(self, name, mesh=mesh, verbose=verbose, **kwargs)

        self.geom_els = geom_els = {}
        for ig, desc in enumerate(mesh.descs):
            gel = GeometryElement(desc)
            # Create geometry elements of dimension - 1.
            gel.create_surface_facet()

            geom_els[desc] = gel

        self.geom_interps = interps = {}
        for gel in geom_els.itervalues():
            key = gel.get_interpolation_name()

            gel.interp = interps.setdefault(key, fea.Interpolant(key, gel))
            gel = gel.surface_facet
            if gel is not None:
                key = gel.get_interpolation_name()
                gel.interp = interps.setdefault(key, fea.Interpolant(key, gel))

        self.vertex_set_bcs = self.mesh.nodal_bcs

        self.mat_ids_to_i_gs = {}
        for ig, mat_id in enumerate(mesh.mat_ids):
            self.mat_ids_to_i_gs[mat_id[0]] = ig

        n_nod, dim = self.mesh.coors.shape
        self.shape = Struct(n_nod=n_nod, dim=dim, tdim=0,
                            n_el=0,
                            n_gr=len(self.mesh.conns))

        self.setup_groups()
        self.fix_element_orientation()
        self.reset_regions()
        self.clear_surface_groups()

        from sfepy.discrete.fem.geometry_element import create_geometry_elements
        from sfepy.discrete.fem.extmods.cmesh import CMesh
        self.cmesh = CMesh.from_mesh(mesh)
        gels = create_geometry_elements()
        self.cmesh.set_local_entities(gels)
        self.cmesh.setup_entities()

        self.shape.tdim = self.cmesh.tdim
        self.cell_offsets = self.mesh.el_offsets

    def setup_groups(self):
        self.groups = {}
        for ii in range(self.shape.n_gr):
            gel = self.geom_els[self.mesh.descs[ii]] # Shortcut.
            conn = self.mesh.conns[ii]
            vertices = nm.unique(conn)

            n_vertex = vertices.shape[0]
            n_el, n_ep = conn.shape
            n_edge = gel.n_edge
            n_edge_total = n_edge * n_el

            if gel.dim == 3:
                n_face = gel.n_face
                n_face_total = n_face * n_el
            else:
                n_face = n_face_total = 0

            shape = Struct(n_vertex=n_vertex, n_el=n_el, n_ep=n_ep,
                           n_edge=n_edge, n_edge_total=n_edge_total,
                           n_face=n_face, n_face_total=n_face_total,
                           dim=self.mesh.dims[ii])

            self.groups[ii] = Struct(ig=ii, vertices=vertices, conn=conn,
                                      gel=gel, shape=shape)
            self.shape.n_el += n_el

    def iter_groups(self, igs=None):
        if igs is None:
            for ig in xrange(self.shape.n_gr): # sorted by ig.
                yield self.groups[ig]
        else:
            for ig in igs:
                yield ig, self.groups[ig]

    def get_cell_offsets(self):
        offs = {}
        for ig in range(self.shape.n_gr):
            offs[ig] = self.cell_offsets[ig]

        return offs

    def get_mesh_coors(self, actual=False):
        """
        Return the coordinates of the underlying mesh vertices.
        """
        if actual and hasattr(self.mesh, 'coors_act'):
            return self.mesh.coors_act
        else:
            return self.mesh.coors

    def get_mesh_bounding_box(self):
        """
        Return the bounding box of the underlying mesh.

        Returns
        -------
        bbox : ndarray (2, dim)
            The bounding box with min. values in the first row and max. values
            in the second row.
        """
        return self.mesh.get_bounding_box()

    def get_diameter(self):
        """
        Return the diameter of the domain.

        Notes
        -----
        The diameter corresponds to the Friedrichs constant.
        """
        bbox = self.get_mesh_bounding_box()
        return (bbox[1,:] - bbox[0,:]).max()

    def fix_element_orientation(self):
        """
        Ensure element nodes ordering giving positive element volume.

        The groups with elements of lower dimension than the space dimension
        are skipped.
        """
        from extmods.cmesh import orient_elements

        coors = self.mesh.coors
        for ii, group in self.groups.iteritems():
            if group.shape.dim < self.shape.dim: continue

            ori, conn = group.gel.orientation, group.conn

            itry = 0
            while itry < 2:
                flag = -nm.ones(conn.shape[0], dtype=nm.int32)

                # Changes orientation if it is wrong according to swap*!
                # Changes are indicated by positive flag.
                orient_elements(flag, conn, coors,
                                ori.roots, ori.vecs,
                                ori.swap_from, ori.swap_to)

                if nm.alltrue(flag == 0):
                    if itry > 0: output('...corrected')
                    itry = -1
                    break

                output('warning: bad element orientation, trying to correct...')
                itry += 1

            if itry == 2 and flag[0] != -1:
                raise RuntimeError('elements cannot be oriented! (%d, %s)'
                                   % (ii, self.mesh.descs[ii]))
            elif flag[0] == -1:
                output('warning: element orienation not checked')

    def get_element_diameters(self, ig, cells, vg, mode, square=True):
        group = self.groups[ig]
        diameters = nm.empty((len(cells), 1, 1, 1), dtype=nm.float64)
        if vg is None:
            diameters.fill(1.0)
        else:
            vg.get_element_diameters(diameters, group.gel.edges,
                                     self.get_mesh_coors().copy(), group.conn,
                                     cells.astype(nm.int32), mode)
        if square:
            out = diameters.squeeze()
        else:
            out = nm.sqrt(diameters.squeeze())

        return out

    def get_evaluate_cache(self, cache=None, share_geometry=False):
        """
        Get the evaluate cache for :func:`Variable.evaluate_at()
        <sfepy.discrete.variables.Variable.evaluate_at()>`.

        Parameters
        ----------
        cache : Struct instance, optional
            Optionally, use the provided instance to store the cache data.
        share_geometry : bool
            Set to True to indicate that all the probes will work on the same
            domain. Certain data are then computed only for the first probe and
            cached.

        Returns
        -------
        cache : Struct instance
            The evaluate cache.
        """
        try:
            from scipy.spatial import cKDTree as KDTree
        except ImportError:
            from scipy.spatial import KDTree

        if cache is None:
            cache = Struct(name='evaluate_cache')

        tt = time.clock()
        if (cache.get('iconn', None) is None) or not share_geometry:
            mesh = self.mesh
            offsets, iconn = make_inverse_connectivity(mesh.conns, mesh.n_nod,
                                                       ret_offsets=True)
            ii = nm.where(offsets[1:] == offsets[:-1])[0]
            if len(ii):
                raise ValueError('some vertices not in any element! (%s)' % ii)

            cache.offsets = offsets
            cache.iconn = iconn
        output('iconn: %f s' % (time.clock()-tt))

        tt = time.clock()
        if (cache.get('kdtree', None) is None) or not share_geometry:
            cache.kdtree = KDTree(mesh.coors)
        output('kdtree: %f s' % (time.clock()-tt))

        return cache

    def clear_surface_groups(self):
        """
        Remove surface group data.
        """
        self.surface_groups = {}

    def create_surface_group(self, region):
        """
        Create a new surface group corresponding to `region` if it does
        not exist yet.

        Notes
        -----
        Surface groups define surface facet connectivity that is needed
        for :class:`sfepy.discrete.fem.mappings.SurfaceMapping`.
        """
        for ig in region.igs:
            groups = self.surface_groups.setdefault(ig, {})
            if region.name not in groups:
                group = self.groups[ig]
                gel_faces = group.gel.get_surface_entities()

                name = 'surface_group_%s_%d' % (region.name, ig)
                surface_group = FESurface(name, region, gel_faces,
                                          group.conn, ig)

                groups[region.name] = surface_group

    def refine(self):
        """
        Uniformly refine the domain mesh.

        Returns
        -------
        domain : FEDomain instance
            The new domain with the refined mesh.

        Notes
        -----
        Works only for meshes with single element type! Does not
        preserve node groups!
        """

        names = set()
        for group in self.groups.itervalues():
            names.add(group.gel.name)

        if len(names) != 1:
            msg = 'refine() works only for meshes with single element type!'
            raise NotImplementedError(msg)

        el_type = names.pop()
        if el_type == '2_3':
            mesh = refine_2_3(self.mesh, self.cmesh)

        elif el_type == '2_4':
            mesh = refine_2_4(self.mesh, self.cmesh)

        elif el_type == '3_4':
            mesh = refine_3_4(self.mesh, self.cmesh)

        elif el_type == '3_8':
            mesh = refine_3_8(self.mesh, self.cmesh)

        else:
            msg = 'unsupported element type! (%s)' % el_type
            raise NotImplementedError(msg)

        domain = FEDomain(self.name + '_r', mesh)

        return domain
