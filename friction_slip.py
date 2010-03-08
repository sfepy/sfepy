#!/usr/bin/env python
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
import sys
import os.path as op

from sfepy.base.base import *
import sfepy.base.la as la
from sfepy.fem import Mesh, Domain, Field, Fields, Variables
from sfepy.fem.fe_surface import FESurface
from sfepy.fem.variables import compute_nodal_normals

def edge_data_to_output(coors, conn, e_sort, data):
    out = nm.zeros_like(coors)
    out[conn[e_sort,0]] = data
    return Struct(name='output_data',
                  mode='vertex', data=out,
                  dofs=None)

def define_dual_mesh(region, region_omega):
    """
    Assume a single GeometryElement type in all groups, linear
    approximation.
    """
    domain = region.domain

    ed = domain.ed
    ## edge_map = {}
    ## for ig, edges in region.edges.iteritems():
        

    ##     data = ed.data[edges]
    ##     for row in data:
    ##         key = tuple(row[3:])
    ##         els = edge_map.setdefault(key, [])
    ##         els.append(row[:3])

    ##     print edge_map
    ##     pause()

    field = Field('displacements', nm.float64, (3,), region_omega, 1)

    gel = domain.geom_els.values()[0]
    sgel = gel.surface_facet
    
    face_key = 's%d' % sgel.n_vertex

    # Coordinate interpolation to face centres.
    ps = gel.interp.poly_spaces[face_key]
    centre = 0.5 * ps.bbox.sum(axis=0)
    bf = ps.eval_base(centre[None,:])

    dual_coors = {}
    dual_conns = {}
    surfaces = {}
    centre_coors = {}
    coor_offsets = {}
    normals = {}
    el_map = {}

    region.setup_face_indices(domain.fa)
    mesh_coors = domain.mesh.coors
    for ig, conn in enumerate(domain.mesh.conns):
        surface = FESurface(None, region, gel.faces, conn, ig)

        print surface

        n_fa, n_edge = surface.n_fa, sgel.n_edge

        surfaces[ig] = surface

        # Face centres.
        fcoors = mesh_coors[surface.econn]
        vals = nm.dot(bf.squeeze(), fcoors)
        centre_coors[ig] = vals

        surface_coors = mesh_coors[surface.nodes]

        dual_coors[ig] = nm.r_[surface_coors, vals]
        coor_offsets[ig] = surface.nodes.shape[0]

        # Normals in primary mesh nodes.
        nodal_normals, imap = compute_nodal_normals(surface.nodes, region,
                                                    field, return_imap=True)
        normals[ig] = nodal_normals

        ee = surface.leconn[:,sgel.edges].copy()
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
        conn = nm.empty((n_edge * n_fa, n_edge), dtype=nm.int32)
        conn[:,0] = e_id
        conn[:,1] = ee[:,0]
        conn[:,2] = nm.repeat(nm.arange(n_fa, dtype=nm.int32), n_edge) \
                    + coor_offsets[ig]
        conn[:,3] = ee[:,1]

        nn = surface.nodes[ueo]
        edge_coors = mesh_coors[nn]

        centre_coors = 0.5 * edge_coors.sum(axis=1)

        edge_normals = 0.5 * nodal_normals[ueo].sum(axis=1)

        edge_normals /= la.norm_l2_along_axis(edge_normals)[:,None]

        nn = surface.nodes[ueo]
        edge_dirs = edge_coors[:,1] - edge_coors[:,0]
        edge_dirs /= la.norm_l2_along_axis(edge_dirs)[:,None]

        edge_ortho = nm.cross(edge_normals, edge_dirs)
        edge_ortho /= la.norm_l2_along_axis(edge_ortho)[:,None]

        print edge_normals
        print edge_dirs
        print edge_ortho

        print centre_coors

        dm_coors = dual_coors[ig]
        dm_conn = conn[:,1:].copy()
        mat_id = nm.zeros((dm_conn.shape[0],), dtype=nm.int32)
        dual_mesh = Mesh.from_data('dual_vis', dm_coors, None, [dm_conn],
                                   [mat_id], ['2_3'])
        dual_mesh.write('aux3.mesh', io='auto')

        dm_coors = nm.r_[centre_coors, dual_coors[ig]]
        dm_conn = conn.copy()
        dm_conn[:,1:] += centre_coors.shape[0]
        mat_id = nm.zeros((dm_conn.shape[0],), dtype=nm.int32)

        out = {}
        out['en'] = edge_data_to_output(dm_coors, dm_conn, e_sort, edge_normals)
        out['ed'] = edge_data_to_output(dm_coors, dm_conn, e_sort, edge_dirs)
        out['eo'] = edge_data_to_output(dm_coors, dm_conn, e_sort, edge_ortho)

        dual_mesh = Mesh.from_data('dual_vis', dm_coors, None, [dm_conn],
                                   [mat_id], ['2_4'])
        dual_mesh.write('aux4.vtk', io='auto', out=out)

        aux = Mesh.from_surface([surface.econn], domain.mesh)
        aux.write('surface.mesh', io='auto')

        debug()

        ## cc = centre_coors[0]
        ## nn = edge_normals[0]
        ## dd = edge_dirs[0]
        ## oo = edge_ortho[0]
        
    print centre_coors
    print dual_coors
    print coor_offsets
    print normals
    debug()


def main():

    mesh = Mesh('mesh', 'meshes/3d/block.mesh')
    print mesh

    domain = Domain('domain', mesh)
    print domain

    reg_omega = domain.create_region('Omega', 'all')
    reg = domain.create_region('Surface',
                               'nodes of surface',
                               {'can_cells' : True})

    dual_mesh = define_dual_mesh(reg, reg_omega)

if __name__ == '__main__':
    main()
