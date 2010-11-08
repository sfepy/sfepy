"""
Friction-slip model formulated as the implicit complementarity problem.

To integrate over a (dual) mesh, one needs:

* coordinates of element vertices
* element connectivity
* local base for each element
  * constant in each sub-triangle of the dual mesh

Data for each dual element:

* connectivity of its sub-triangles
* base directions t_1, t_2

Normal stresses:

* Assemble the rezidual and apply the LCBC operator described below.

Solution in \hat{V}_h^c:

* construct a restriction operator via LCBC just like in the no-penetration case
* use the substitution:
  u_1 = n_1 * w
  u_2 = n_2 * w
  u_3 = n_3 * w
  The new DOF is `w`.
* for the record, no-penetration does:
  w_1 = - (1 / n_1) * (u_2 * n_2 + u_3 * n_3)
  w_2 = u_2
  w_3 = u_3
"""
import numpy as nm
from copy import copy

from sfepy.base.base import Struct, assert_
from sfepy.base.compat import unique
import sfepy.linalg as la
from sfepy.fem import Mesh, Field, Region
from sfepy.fem.mappings import VolumeMapping, SurfaceMapping
from sfepy.fem.conditions import EssentialBC
from sfepy.fem.fe_surface import FESurface
from sfepy.fem.utils import compute_nodal_normals

def edge_data_to_output(coors, conn, e_sort, data):
    out = nm.zeros_like(coors)
    out[conn[e_sort,0]] = data
    return Struct(name='output_data',
                  mode='vertex', data=out,
                  dofs=None)

class DualMesh(Struct):
    """Dual mesh corresponding to a (surface) region."""

    def __init__(self, region):
        """
        Assume a single GeometryElement type in all groups, linear
        approximation.

        Works for one group only for the moment.
        """
        domain = region.domain

        self.dim = domain.shape.dim

        self.region = copy(region)
        self.region.setup_face_indices()

        self.mesh_coors = domain.mesh.coors

        # add_to_regions=True due to Field implementation shortcomings.
        omega = domain.create_region('Omega', 'all', add_to_regions=True)
        self.field = Field('displacements', nm.float64, (3,), omega, 1)

        self.gel = domain.geom_els.values()[0]
        self.sgel = self.gel.surface_facet

        face_key = 's%d' % self.sgel.n_vertex

        # Coordinate interpolation to face centres.
        self.ps = self.gel.interp.poly_spaces[face_key]
        centre = self.ps.node_coors.sum(axis=0) / self.ps.n_nod
        self.bf = self.ps.eval_base(centre[None,:])

        self.surfaces = surfaces = {}
        self.dual_surfaces = dual_surfaces = {}

        n_nod = n_fa = 0
        for ig, conn in enumerate(domain.mesh.conns):
            surface = FESurface(None, self.region, self.gel.faces, conn, ig)
            surfaces[ig] = surface

            dual_surface = self.describe_dual_surface(surface)
            dual_surfaces[ig] = dual_surface

            n_nod += dual_surface.n_nod
            n_fa += dual_surface.n_fa

        # Domain-like shape.
        self.shape = Struct(n_nod=n_nod, n_el=n_fa, dim=self.dim)

    def describe_dual_surface(self, surface):
        n_fa, n_edge = surface.n_fa, self.sgel.n_edge

        mesh_coors = self.mesh_coors

        # Face centres.
        fcoors = mesh_coors[surface.econn]
        centre_coors = nm.dot(self.bf.squeeze(), fcoors)

        surface_coors = mesh_coors[surface.nodes]

        dual_coors = nm.r_[surface_coors, centre_coors]
        coor_offset = surface.nodes.shape[0]

        # Normals in primary mesh nodes.
        nodal_normals = compute_nodal_normals(surface.nodes, self.region,
                                              self.field)

        ee = surface.leconn[:,self.sgel.edges].copy()
        edges_per_face = ee.copy()
        sh = edges_per_face.shape
        ee.shape = edges_per_face.shape = (sh[0] * sh[1], sh[2])
        edges_per_face.sort(axis=1)

        eo = nm.empty((sh[0] * sh[1],), dtype=nm.object)
        eo[:] = [tuple(ii) for ii in edges_per_face]

        ueo, e_sort, e_id = unique(eo, return_index=True, return_inverse=True)
        ueo = edges_per_face[e_sort]

        # edge centre, edge point 1, face centre, edge point 2
        conn = nm.empty((n_edge * n_fa, 4), dtype=nm.int32)
        conn[:,0] = e_id
        conn[:,1] = ee[:,0]
        conn[:,2] = nm.repeat(nm.arange(n_fa, dtype=nm.int32), n_edge) \
                    + coor_offset
        conn[:,3] = ee[:,1]

        # face centre, edge point 2, edge point 1
        tri_conn = nm.ascontiguousarray(conn[:,[2,1,3]])

        # Ensure orientation - outward normal.
        cc = dual_coors[tri_conn]
        v1 = cc[:,1] - cc[:,0]
        v2 = cc[:,2] - cc[:,0]

        normals = nm.cross(v1, v2)
        nn = nodal_normals[surface.leconn].sum(axis=1).repeat(n_edge, 0)
        centre_normals = (1.0 / surface.n_fp) * nn
        centre_normals /= la.norm_l2_along_axis(centre_normals)[:,None]
        dot = nm.sum(normals * centre_normals, axis=1)

        assert_((dot > 0.0).all())

        # Prepare mapping from reference triangle e_R to a
        # triangle within reference face e_D.
        gel = self.gel.surface_facet
        ref_coors = gel.coors
        ref_centre = nm.dot(self.bf.squeeze(), ref_coors)
        cc = nm.r_[ref_coors, ref_centre[None,:]]
        rconn = nm.empty((n_edge, 3), dtype=nm.int32)
        rconn[:,0] = gel.n_vertex
        rconn[:,1] = gel.edges[:,0]
        rconn[:,2] = gel.edges[:,1]

        map_er_ed = VolumeMapping(cc, rconn, gel=gel)

        # Prepare mapping from reference triangle e_R to a
        # physical triangle e.
        map_er_e = SurfaceMapping(dual_coors, tri_conn, gel=gel)

        # Compute triangle basis (edge) vectors.
        nn = surface.nodes[ueo]
        edge_coors = mesh_coors[nn]

        edge_centre_coors = 0.5 * edge_coors.sum(axis=1)

        edge_normals = 0.5 * nodal_normals[ueo].sum(axis=1)

        edge_normals /= la.norm_l2_along_axis(edge_normals)[:,None]

        nn = surface.nodes[ueo]
        edge_dirs = edge_coors[:,1] - edge_coors[:,0]
        edge_dirs /= la.norm_l2_along_axis(edge_dirs)[:,None]

        edge_ortho = nm.cross(edge_normals, edge_dirs)
        edge_ortho /= la.norm_l2_along_axis(edge_ortho)[:,None]

        # Primary face - dual sub-faces map.
        # i-th row: indices to conn corresponding to sub-faces of i-th face.
        face_map = nm.arange(n_fa * n_edge, dtype=nm.int32)
        face_map.shape = (n_fa, n_edge)

        # The actual connectivity for assembling (unique nodes per master
        # faces).
        asm_conn = e_id[face_map]

        n_nod = ueo.shape[0] # One node per unique edge.
        n_components = self.dim - 1

        dual_surface = Struct(name = 'dual_surface_description',
                              dim = self.dim,
                              n_dual_fa = conn.shape[0],
                              n_dual_fp = self.dim,
                              n_fa = n_fa,
                              n_edge = n_edge,
                              n_nod = n_nod,
                              n_components = n_components,
                              n_dof = n_nod * n_components,
                              dual_coors = dual_coors,
                              coor_offset = coor_offset,
                              e_sort = e_sort,
                              conn = conn,
                              tri_conn = tri_conn,
                              map_er_e = map_er_e,
                              map_er_ed = map_er_ed,
                              face_map = face_map,
                              asm_conn = asm_conn,
                              nodal_normals = nodal_normals,
                              edge_centre_coors = edge_centre_coors,
                              edge_normals = edge_normals,
                              edge_dirs = edge_dirs,
                              edge_ortho = edge_ortho)

        return dual_surface

    def iter_groups(self, igs=None):
        """
        Domain-like functionality.
        """
        if igs is None:
            for ig, dual_surface in self.dual_surfaces.iteritems():
                vertices = nm.arange(dual_surface.n_nod, dtype=nm.int32)
                group = Struct(ig=ig, vertices=vertices,
                               conn=dual_surface.asm_conn)
                yield group

        else:
            for ig in igs:
                dual_surface = self.dual_surfaces[ig]
                vertices = nm.arange(dual_surface.n_nod, dtype=nm.int32)
                group = Struct(ig=ig, vertices=vertices,
                               conn=dual_surface.asm_conn)
                yield ig, group

    def create_friction_bcs(self, dof_name):
        """
        Fix friction DOFs on surface boundary edges, i.e. that are not
        shared by two friction surface faces.
        """
        bcs = []
        for ig, dual_surface in self.dual_surfaces.iteritems():
            e_id = dual_surface.conn[:, 0]
            ii = nm.where(nm.bincount(e_id) == 1)

            region = Region('__friction_%d' % ig, '', self, '')
            region.set_vertices(ii)
            region.is_complete = True

            dofs = {'%s.all' % dof_name : 0.0}

            bc = EssentialBC('__friction_%d' % ig, region, dofs)
            bcs.append(bc)

        return bcs

    def save(self, filename):
        coors = []
        conns = []
        mat_ids = []
        offset = 0
        for ig, dual_surface in self.dual_surfaces.iteritems():
            cc = dual_surface.dual_coors
            coors.append(cc)

            conn = dual_surface.conn[:,1:].copy() + offset
            conns.append(conn)
            
            mat_id = nm.empty((conn.shape[0],), dtype=nm.int32)
            mat_id[:] = ig
            mat_ids.append(mat_id)

            offset += cc.shape[0]

        coors = nm.concatenate(coors, axis=0)

        dual_mesh = Mesh.from_data('dual_mesh', coors, None, conns,
                                   mat_ids, ['2_3'] * len(conns))
        dual_mesh.write(filename, io='auto')

    def save_axes(self, filename):
        coors = []
        conns = []
        mat_ids = []
        offset = 0
        for ig, dual_surface in self.dual_surfaces.iteritems():
            cc = nm.r_[dual_surface.edge_centre_coors,
                       dual_surface.dual_coors]
            coors.append(cc)

            conn = dual_surface.conn.copy() + offset
            conn[:,1:] += dual_surface.edge_centre_coors.shape[0]
            conns.append(conn)
            
            mat_id = nm.empty((conn.shape[0],), dtype=nm.int32)
            mat_id[:] = ig
            mat_ids.append(mat_id)

            offset += cc.shape[0]

        coors = nm.concatenate(coors, axis=0)

        out = {}
        for ig, dual_surface in self.dual_surfaces.iteritems():
            eto = edge_data_to_output
            out['en_%d' % ig] = eto(coors, conns[ig], dual_surface.e_sort,
                                    dual_surface.edge_normals)
            out['ed_%d' % ig] = eto(coors, conns[ig], dual_surface.e_sort,
                                    dual_surface.edge_dirs)
            out['eo_%d' % ig] = eto(coors, conns[ig], dual_surface.e_sort,
                                    dual_surface.edge_ortho)

        dual_mesh = Mesh.from_data('dual_mesh_vectors', coors, None, conns,
                                   mat_ids, ['2_4'] * len(conns))
        dual_mesh.write(filename, io='auto', out=out)

