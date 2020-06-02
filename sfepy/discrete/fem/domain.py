"""
Computational domain, consisting of the mesh and regions.
"""
from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import output, Struct
from .geometry_element import GeometryElement
from sfepy.discrete import Domain, PolySpace
from sfepy.discrete.fem.refine import refine_2_3, refine_2_4, refine_3_4, \
    refine_3_8, refine_1_2
from sfepy.discrete.fem.fe_surface import FESurface
import six

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

        if len(mesh.descs) > 1:
            msg = 'meshes with several cell kinds are not supported!'
            raise NotImplementedError(msg)

        self.geom_els = geom_els = {}
        for ig, desc in enumerate(mesh.descs):
            gel = GeometryElement(desc)

            if gel.dim > 0:
                # Create geometry elements of dimension - 1.
                gel.create_surface_facet()

            geom_els[desc] = gel

        for gel in six.itervalues(geom_els):
            key = gel.get_interpolation_name()

            gel.poly_space = PolySpace.any_from_args(key, gel, 1)
            gel = gel.surface_facet
            if gel is not None:
                key = gel.get_interpolation_name()
                gel.poly_space = PolySpace.any_from_args(key, gel, 1)

        self.vertex_set_bcs = self.mesh.nodal_bcs

        self.cmesh = self.mesh.cmesh

        # Must be before creating derived connectivities.
        self.fix_element_orientation()

        from sfepy.discrete.fem.geometry_element import create_geometry_elements
        gels = create_geometry_elements()
        self.cmesh.set_local_entities(gels)
        self.cmesh.setup_entities()

        n_nod, dim = self.mesh.coors.shape
        self.shape = Struct(n_nod=n_nod, dim=dim, tdim=self.cmesh.tdim,
                            n_el=self.cmesh.n_el,
                            n_gr=len(self.geom_els))

        self.reset_regions()
        self.clear_surface_groups()

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
        Ensure element vertices ordering giving positive cell volumes.
        """
        from sfepy.discrete.common.extmods.cmesh import orient_elements

        if self.cmesh.tdim != self.cmesh.dim:
            output('warning: mesh with topological dimension %d lower than'
                   ' space dimension %d' % (self.cmesh.tdim, self.cmesh.dim))
            output('- element orientation not checked!')
            return

        cmesh = self.cmesh
        for key, gel in six.iteritems(self.geom_els):
            ori = gel.orientation

            cells = nm.where(cmesh.cell_types == cmesh.key_to_index[gel.name])
            cells = cells[0].astype(nm.uint32)

            itry = 0
            while itry < 2:
                flag = -nm.ones(self.cmesh.n_el, dtype=nm.int32)

                # Changes orientation if it is wrong according to swap*!
                # Changes are indicated by positive flag.
                orient_elements(flag, self.cmesh, cells, gel.dim,
                                ori.roots, ori.vecs,
                                ori.swap_from, ori.swap_to)

                if nm.alltrue(flag == 0):
                    if itry > 0: output('...corrected')
                    itry = -1
                    break

                output('warning: bad element orientation, trying to correct...')
                itry += 1

            if itry == 2 and flag[0] != -1:
                raise RuntimeError('elements cannot be oriented! (%s)' % key)

            elif flag[0] == -1:
                output('warning: element orienation not checked')

    def get_conn(self, ret_gel=False):
        """
        Get the cell-vertex connectivity and, if `ret_gel` is True, also the
        corresponding reference geometry element.
        """
        conn = self.cmesh.get_conn(self.cmesh.tdim, 0).indices
        conn = conn.reshape((self.cmesh.n_el, -1)).astype(nm.int32)

        if ret_gel:
            gel = list(self.geom_els.values())[0]

            return conn, gel

        else:
            return conn

    def get_element_diameters(self, cells, vg, mode, square=True):
        diameters = nm.empty((len(cells), 1, 1, 1), dtype=nm.float64)
        if vg is None:
            diameters.fill(1.0)
        else:
            conn, gel = self.get_conn(ret_gel=True)
            vg.get_element_diameters(diameters, gel.edges,
                                     self.get_mesh_coors().copy(), conn,
                                     cells.astype(nm.int32), mode)
        if square:
            out = diameters.squeeze()
        else:
            out = nm.sqrt(diameters.squeeze())

        return out

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
        groups = self.surface_groups
        if region.name not in groups:
            conn, gel = self.get_conn(ret_gel=True)
            gel_faces = gel.get_surface_entities()

            name = 'surface_group_%s' % (region.name)
            surface_group = FESurface(name, region, gel_faces, conn)

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
        if len(self.geom_els) != 1:
            msg = 'refine() works only for meshes with single element type!'
            raise NotImplementedError(msg)

        el_type = list(self.geom_els.values())[0].name
        if el_type == '1_2':
            mesh = refine_1_2(self.mesh)

        elif el_type == '2_3':
            mesh = refine_2_3(self.mesh)

        elif el_type == '2_4':
            mesh = refine_2_4(self.mesh)

        elif el_type == '3_4':
            mesh = refine_3_4(self.mesh)

        elif el_type == '3_8':
            mesh = refine_3_8(self.mesh)

        else:
            msg = 'unsupported element type! (%s)' % el_type
            raise NotImplementedError(msg)

        domain = FEDomain(self.name + '_r', mesh)

        return domain
