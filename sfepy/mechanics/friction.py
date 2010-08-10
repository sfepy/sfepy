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

Surface integral ??

Normal stresses:

* Just assemble the rezidual as in r.h.s. of (91) and put nodal normals into
  contact nodes into \hat{v}.

Solution in \hat{V}_h^c (not needed!):

* construct a restriction operator via LCBC just like in the no-penetration case
* use the substitution:
  w_1 = u_1
  w_2 = (n_2 / n_1) * u_1
  w_3 = (n_3 / n_1) * u_1
* for the record, no-penetration does:
  w_1 = - (1 / n_1) * (u_2 * n_2 + u_3 * n_3)
  w_2 = u_2
  w_3 = u_3
"""
from sfepy.base.base import *
import sfepy.base.la as la
from sfepy.fem import Mesh, Domain, Field, Variables
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

        self.region = copy(region)
        self.region.setup_face_indices(domain.fa)

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
        
        el_map = {}

        for ig, conn in enumerate(domain.mesh.conns):
            surface = FESurface(None, self.region, self.gel.faces, conn, ig)
            surfaces[ig] = surface

            dual_surface = self.describe_dual_surface(surface)
            dual_surfaces[ig] = dual_surface

            print dual_surface

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

        ueo, e_sort, e_id = nm.unique1d(eo, return_index=True,
                                        return_inverse=True)
        ueo = edges_per_face[e_sort]

        # edge centre, edge point 1, face centre, edge point 2.
        conn = nm.empty((n_edge * n_fa, 4), dtype=nm.int32)
        conn[:,0] = e_id
        conn[:,1] = ee[:,0]
        conn[:,2] = nm.repeat(nm.arange(n_fa, dtype=nm.int32), n_edge) \
                    + coor_offset
        conn[:,3] = ee[:,1]

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

        dual_surface = Struct(name = 'dual_surface_description',
                              dual_coors = dual_coors,
                              coor_offset = coor_offset,
                              e_sort = e_sort,
                              conn = conn,
                              nodal_normals = nodal_normals,
                              edge_centre_coors = edge_centre_coors,
                              edge_normals = edge_normals,
                              edge_dirs = edge_dirs,
                              edge_ortho = edge_ortho)

        return dual_surface

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


