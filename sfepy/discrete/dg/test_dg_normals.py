"""
11.12.2018
"""
import numpy as nm
import matplotlib.pyplot as plt

from sfepy.discrete.fem import FEDomain, Mesh
from sfepy.mesh.mesh_generators import gen_block_mesh

import sfepy.postprocess.plot_cmesh as pc


mesh = gen_block_mesh([1, 1], [4, 4], [0, 0], verbose=False)
#mesh = gen_block_mesh([1, 1, 1], [4, 4, 3], [0, 0, 0], verbose=False)
#mesh = Mesh.from_file('meshes/2d/its2D.mesh')

domain = FEDomain('domain', mesh)
omega = domain.create_region('Omega', 'all')

cmesh = mesh.cmesh
dim = cmesh.dim
cmesh.setup_connectivity(dim, dim)
cmesh.setup_connectivity(dim, dim - 1)
cmesh.setup_connectivity(dim - 1, dim)
c2c = cmesh.get_conn(dim, dim)

ccs, coffs = cmesh.get_incident(dim, omega.cells, dim, ret_offsets=True)

scale = 0.05 * nm.diff(mesh.get_bounding_box(), axis=0).min()

normals = cmesh.get_facet_normals()

def get_facet_neighbours(cmesh, cells):
    dim = cmesh.dim
    c2fi, c2fo = cmesh.get_incident(dim-1, cells, dim, ret_offsets=True)

    neighbours = []
    for ic, o1 in enumerate(c2fo[:-1]):
        o2 = c2fo[ic+1]

        print ic, c2fi[o1:o2]

        c2ci, c2co = cmesh.get_incident(dim, c2fi[o1:o2], dim-1,
                                        ret_offsets=True)
        nbrs = []
        for ifa, of1 in enumerate(c2co[:-1]):
            of2 = c2co[ifa+1]
            if of2 == (of1 + 1):
                print ifa, of1, of2, c2ci[of1]
                # Surface facet.
                nbrs.append(c2ci[of1])

            else:
                print ifa, of1, of2, c2ci[of1], c2ci[of2-1]
                if c2ci[of1] == cells[ic]:
                    nbrs.append(c2ci[of2-1])

                else:
                    nbrs.append(c2ci[of1])

        neighbours.append(nbrs)

    return neighbours

neighbours = get_facet_neighbours(cmesh, omega.cells)
print neighbours

def plot_facet_normals(ax, cmesh, normals, scale, neighbours):
    dim = cmesh.dim
    ax = pc._get_axes(ax, dim)

    edim = dim - 1
    coors = cmesh.get_centroids(edim)
    coors = pc._to2d(coors)

    cmesh.setup_connectivity(dim, edim)
    c2f = cmesh.get_conn(dim, edim)
    for ic, o1 in enumerate(c2f.offsets[:-1]):
        o2 = c2f.offsets[ic+1]
        for ifal, ifa in enumerate(c2f.indices[o1:o2]):
            print ic, ifal, ifa
            cc = nm.array([coors[ifa], coors[ifa] + scale * normals[o1 + ifal]])
            print cc
            ax.plot(*cc.T, color='m')
            ax.text(*cc[1], s=neighbours[ic][ifal],
                     horizontalalignment='center', verticalalignment='center')

ax = pc.plot_cmesh(
    None, cmesh,
    wireframe_opts = {'color' : 'k', 'linewidth' : 2},
    entities_opts=[
        {'color' : 'k', 'label_global' : 12, 'label_local' : 8, 'size' : 20},
        {'color' : 'b', 'label_global' : 12, 'label_local' : 8, 'size' : 10},
        {'color' : 'r', 'label_global' : 20, 'size' : 20},
    ])

ax = plot_facet_normals(ax, cmesh, normals, scale, neighbours)

plt.show()
