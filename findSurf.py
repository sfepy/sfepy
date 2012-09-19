#!/usr/bin/env python
# 05.10.2005, c
"""
Given a mesh file, this script extracts its surface and prints it to stdout in
form of a list where each row is [group, element, face, component]. A component
corresponds to a contiguous surface region - for example, a cubical mesh with a
spherical hole has two surface components. Two surface faces sharing a single
node belong to one component.

With '-m' option, a mesh of the surface is created and saved in
'surf_<original mesh file name>.mesh'.

Try ./find_surf.py --help to see more options.
"""
import sys
from optparse import OptionParser

import numpy as nm
import scipy.sparse as sp

import sfepy
from sfepy.base.ioutils import edit_filename
from sfepy.fem import Mesh, Domain
from sfepy.fem.extmods.mesh import create_mesh_graph, graph_components

##
# 29.08.2007, c
def surface_graph( surf_faces, n_nod ):
    nnz, prow, icol = create_mesh_graph(n_nod, n_nod, len(surf_faces),
                                        surf_faces, surf_faces)
    data = nm.empty( (nnz,), dtype = nm.int32 )
    data.fill( 2 )
    return sp.csr_matrix( (data, icol, prow), (n_nod, n_nod) )

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

'-' is for stdin, stdout"""

##
# c: 05.10.2005, r: 09.07.2008
def main():
    parser = OptionParser(usage = usage, version = "%prog " + sfepy.__version__)
    parser.add_option( "-m", "--mesh",
                       action = "store_true", dest = "save_mesh",
                       default = True,
                       help = "save surface mesh [default: %default]" )
    parser.add_option( "-n", "--no-surface",
                       action = "store_true", dest = "no_surface",
                       default = False,
                       help = "do not output surface [default: %default]" )
    (options, args) = parser.parse_args()

    if (len( args ) == 2):
        filename_in = args[0];
        filename_out = args[1];
    else:
        parser.print_help(),
        return

    if (filename_in == '-'):
        file_in = sys.stdin
    else:
        file_in = open( filename_in, "r" ); 

    mesh = Mesh.from_file( filename_in )

    if (filename_in != '-'):
        file_in.close()

    domain = Domain('domain', mesh)
    domain.setup_groups()

    if domain.has_faces():
        domain.fix_element_orientation()
        domain.setup_facets(create_edges=False)

        lst, surf_faces = domain.surface_faces()

        surf_mesh = Mesh.from_surface( surf_faces, mesh )

        if options.save_mesh:
            aux = edit_filename(filename_in, prefix='surf_', new_ext='.mesh')
            surf_mesh.write(aux, io='auto')

        if options.no_surface:
            return

        gr_s = surface_graph( surf_faces, mesh.n_nod )
##         import sfepy.base.plotutils as plu
##         plu.spy( gr_s )
##         plu.pylab.show()

        n_comp, comps = surface_components( gr_s, surf_faces )
#        print 'components:', n_comp

        ccs, comps = comps, nm.zeros( (0,1), nm.int32 )
        for cc in ccs:
            comps = nm.concatenate( (comps, cc[:,nm.newaxis]), 0 )

        out = nm.concatenate( (lst, comps), 1 )

        if (filename_out == '-'):
            file_out = sys.stdout
        else:
            file_out = open( filename_out, "w" ); 
        for row in out:
            file_out.write( '%d %d %d %d\n' % (row[0], row[1], row[2], row[3]) )
        if (filename_out != '-'):
            file_out.close()
    

if __name__=='__main__':
    main()
