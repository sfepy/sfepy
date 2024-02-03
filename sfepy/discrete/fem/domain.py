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

            if gel.dim > 0:
                # Create geometry elements of dimension - 1.
                gel.create_surface_facet()

            geom_els[desc] = gel

        for gel in geom_els.values():
            key = gel.get_interpolation_name()

            gel.poly_space = PolySpace.any_from_args(key, gel, 1)
            gel = gel.surface_facet
            if isinstance(gel, dict):
                for v in gel.values():
                    key = v.get_interpolation_name()
                    v.poly_space = PolySpace.any_from_args(key, v, 1)
            elif gel is not None:
                key = gel.get_interpolation_name()
                gel.poly_space = PolySpace.any_from_args(key, gel, 1)

        self.vertex_set_bcs = self.mesh.nodal_bcs
        self.cmesh_tdim = self.mesh.cmesh_tdim

        # Must be before creating derived connectivities.
        self.fix_element_orientation()

        from sfepy.discrete.fem.geometry_element import create_geometry_elements
        gels = create_geometry_elements()
        max_tdim = 0
        for cmesh in self.mesh.cmesh_tdim:
            if cmesh is not None:
                cmesh.set_local_entities(gels)
                cmesh.setup_entities()
                max_tdim = max(max_tdim, cmesh.tdim)

        self.cmesh = self.cmesh_tdim[max_tdim]
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

    def fix_element_orientation(self, geom_els=None, force_check=False):
        """
        Ensure element vertices ordering giving positive cell volumes.
        """
        from sfepy.discrete.common.extmods.cmesh import orient_elements

        if geom_els is None:
            geom_els = self.geom_els

        for key, gel in geom_els.items():
            ori = gel.orientation
            cmesh = self.cmesh_tdim[gel.dim]
            if (cmesh.tdim != cmesh.dim) and (not force_check):
                output('warning: mesh with topological dimension %d lower than'
                       ' space dimension %d' % (cmesh.tdim, cmesh.dim))
                output('- element orientation not checked!')
                return

            cells = nm.where(cmesh.cell_types == cmesh.key_to_index[gel.name])
            cells = cells[0].astype(nm.uint32)

            itry = 0
            while itry < 2:
                flag = -nm.ones(cmesh.n_el, dtype=nm.int32)

                # Changes orientation if it is wrong according to swap*!
                # Changes are indicated by positive flag.
                orient_elements(flag, cmesh, cells, gel.dim,
                                ori.roots, ori.vecs,
                                ori.swap_from, ori.swap_to)

                flag = flag[cells]

                if nm.all(flag == 0):
                    if itry > 0: output('...corrected')
                    itry = -1
                    break

                output('warning: bad element orientation, trying to correct...')
                itry += 1

            if itry == 2 and flag[0] != -1:
                raise RuntimeError('elements cannot be oriented! (%s)' % key)

            elif flag[0] == -1:
                output('warning: element orienation not checked')

    def get_conn(self, ret_gel=False, tdim=None, cells=None):
        """
        Get the cell-vertex connectivity and, if `ret_gel` is True, also the
        corresponding reference geometry element. If `tdim` is not None get
        the connectivity of the cells with topological dimension `tdim`.
        """
        cmesh = self.cmesh if tdim is None else self.cmesh_tdim[tdim]
        if cells is None:
            conn = cmesh.get_conn(cmesh.tdim, 0).indices
            conn = conn.reshape((cmesh.n_el, -1)).astype(nm.int32)
        else:
            conn = cmesh.get_conn(cmesh.tdim, 0)
            aux = conn.offsets[cells]
            nv = set(conn.offsets[cells + 1] - aux)
            if len(nv) != 1:
                raise ValueError('different cell types in a region!')

            nc = cells.shape[0]
            nv = nv.pop()
            idxs = (nm.repeat(conn.offsets[cells], nv) +
                    nm.tile(nm.arange(nv), nc))
            conn = conn.indices[idxs].reshape((nc, nv))

        name = f'{cmesh.tdim}_{conn.shape[1]}'
        if ret_gel:
            gel = None
            for gv in self.geom_els.values():
                if gv.name == name:
                    gel = gv
                    break

            return conn, gel

        else:
            return conn

    def get_element_diameters(self, cells, volume, mode, square=True):
        conn, gel = self.get_conn(ret_gel=True)
        coors = self.get_mesh_coors()
        lconn = conn[cells]
        ed = gel.edges

        if mode == 0 or mode == 2:
            vv = coors[lconn[:, ed[:, 1]]] - coors[lconn[:, ed[:, 0]]]
            v0 = nm.max(nm.sum(vv**2, axis=2), axis=1)
            out = v0

        if mode == 1 or mode == 2:
            v1 = nm.power(0.16 * volume, 1.0 / gel.dim).squeeze()
            out = v1

        if mode == 2:
            out = nm.max(nm.stack([v0, v1]), axis=0)

        return out if square else nm.sqrt(out)

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
        for :class:`sfepy.discrete.fem.mappings.FEMapping`.
        """
        groups = self.surface_groups
        if region.name not in groups:
            name = 'surface_group_%s' % (region.name)
            surface_group = FESurface.from_region(name, region)

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
