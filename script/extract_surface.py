#!/usr/bin/env python
# 05.10.2005, c
"""
Given a mesh file, this script extracts its surface and prints it to stdout in
form of a list where each row is [element, face, component]. A component
corresponds to a contiguous surface region - for example, a cubical mesh with a
spherical hole has two surface components. Two surface faces sharing a single
node belong to one component.

With '-m' option, a mesh of the surface is created and saved in
'<original path>/surf_<original mesh file name>.mesh'.
"""
import sys
from optparse import OptionParser

import numpy as nm
import scipy.sparse as sp

import sfepy
from sfepy.base.base import output
from sfepy.base.ioutils import edit_filename
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.discrete.fem.extmods.cmesh import create_mesh_graph, graph_components

def _get_facets(vertices, offsets, ii, n_fp):
    facets = []
    for ic in range(n_fp):
        facets.append(vertices[offsets[ii] + ic][:, None])

    facets = nm.concatenate(facets, axis=1)

    return nm.ascontiguousarray(facets.astype(nm.int32))

def get_surface_faces(domain):
    cmesh = domain.cmesh
    faces = cmesh.get_surface_facets()
    vertices_f, offs_f = cmesh.get_incident(0, faces,
                                            cmesh.dim - 1, ret_offsets=True)

    n_fp = nm.diff(offs_f)
    surf_faces = []

    itri = nm.where(n_fp == 3)[0]
    if itri.size:
        surf_faces.append(_get_facets(vertices_f, offs_f, itri, 3))

    itet = nm.where(n_fp == 4)[0]
    if itet.size:
        surf_faces.append(_get_facets(vertices_f, offs_f, itet, 4))

    cells_c, offs_c = cmesh.get_incident(cmesh.dim, faces, cmesh.dim - 1,
                                         ret_offsets=True)
    ids = cmesh.get_local_ids(faces, cmesh.dim - 1, cells_c, offs_c,
                              cmesh.dim)
    lst = nm.c_[cells_c, ids]

    return lst, surf_faces

def surface_graph(surf_faces, n_nod):
    nnz, prow, icol = create_mesh_graph(n_nod, n_nod, len(surf_faces),
                                        surf_faces, surf_faces)
    data = nm.empty((nnz,), dtype=nm.int32)
    data.fill(2)
    return sp.csr_matrix((data, icol, prow), (n_nod, n_nod))

def surface_components(gr_s, surf_faces):
    """
    Determine surface components given surface mesh connectivity graph.
    """
    n_nod = gr_s.shape[0]
    n_comp, flag = graph_components(n_nod, gr_s.indptr, gr_s.indices)

    comps = []
    for ii, face in enumerate(surf_faces):
        comp = flag[face[:,0]]
        comps.append(comp)

    return n_comp, comps

usage = """%prog [options] filename_in|- filename_out|-

'-' is for stdin, stdout
""" + __doc__.rstrip()

def main():
    parser = OptionParser(usage=usage, version="%prog " + sfepy.__version__)
    parser.add_option("-m", "--mesh",
                      action="store_true", dest="save_mesh",
                      default=True,
                      help="save surface mesh [default: %default]")
    parser.add_option("-n", "--no-surface",
                      action="store_true", dest="no_surface",
                      default=False,
                      help="do not output surface [default: %default]")
    (options, args) = parser.parse_args()

    if (len(args) == 2):
        filename_in = args[0];
        filename_out = args[1];
    else:
        parser.print_help(),
        return

    if (filename_in == '-'):
        file_in = sys.stdin
    else:
        file_in = open(filename_in, "r");

    mesh = Mesh.from_file(filename_in)

    if (filename_in != '-'):
        file_in.close()

    domain = FEDomain('domain', mesh)
    domain.setup_groups()

    if domain.has_faces():
        domain.fix_element_orientation()

        lst, surf_faces = get_surface_faces(domain)

        surf_mesh = Mesh.from_surface(surf_faces, mesh)

        if options.save_mesh:
            aux = edit_filename(filename_in, prefix='surf_', new_ext='.mesh')
            surf_mesh.write(aux, io='auto')

        if options.no_surface:
            return

        gr_s = surface_graph(surf_faces, mesh.n_nod)

        n_comp, comps = surface_components(gr_s, surf_faces)
        output('number of surface components:', n_comp)

        ccs, comps = comps, nm.zeros((0,1), nm.int32)
        for cc in ccs:
            comps = nm.concatenate((comps, cc[:,nm.newaxis]), 0)

        out = nm.concatenate((lst, comps), 1)

        if (filename_out == '-'):
            file_out = sys.stdout
        else:
            file_out = open(filename_out, "w");
        for row in out:
            file_out.write('%d %d %d\n' % (row[0], row[1], row[2]))
        if (filename_out != '-'):
            file_out.close()

if __name__=='__main__':
    main()
