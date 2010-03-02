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

        edges_per_face = surface.leconn[:,sgel.edges].copy()
        sh = edges_per_face.shape
        edges_per_face.shape = (sh[0] * sh[1], sh[2])
        edges_per_face.sort(axis=1)

        edge_normals = 0.5 * nodal_normals[edges_per_face].sum(axis=1)

        edge_normals /= la.norm_l2_along_axis(edge_normals)[:,None]

        nn = surface.nodes[edges_per_face]
        edge_dirs = mesh_coors[nn[:,1]] = mesh_coors[nn[:,0]]
        edge_dirs /= la.norm_l2_along_axis(edge_dirs)[:,None]

        edge_ortho = nm.cross(edge_normals, edge_dirs)

        print edge_normals
        print edge_dirs
        print edge_ortho

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
    reg = domain.create_region('Screw',
                               ' +n '.join(['nodes of group %d'
                                            % ii for ii in [1,2,3,4]]),
                               {'can_cells' : True})

    dual_mesh = define_dual_mesh(reg, reg_omega)

if __name__ == '__main__':
    main()
