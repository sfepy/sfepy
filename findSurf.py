#!/usr/bin/env python
# 05.10.2005, c
import sys
import os.path as op
from optparse import OptionParser

import init_sfe

from sfe.base.base import *
from sfe.fem.mesh import Mesh
from sfe.fem.domain import Domain
from sfe.fem.extmods.fem import rawGraph
from sfe.fem.extmods.meshutils import graphComponents

##
# 29.08.2007, c
def surfaceGraph( surfFaces, nNod ):
    ret, prow, icol = rawGraph( nNod, nNod,
                                len( surfFaces), surfFaces, surfFaces )
    nnz = prow[-1]
    data = nm.empty( (nnz,), dtype = nm.int32 )
    data.fill( 2 )
    return sp.csr_matrix( (data, icol, prow), (nNod, nNod) )

##
# 29.08.2007, c
def surfaceComponents( grS, surfFaces ):
    nNod = grS.shape[0]
    flag = nm.empty( (nNod,), dtype = nm.int32 )
    pos = nm.empty_like( flag )
    ret, nComp = graphComponents( flag, grS.indptr, grS.indices, pos )

    comps = []
    for ii, face in enumerate( surfFaces ):
        comp = flag[face[:,0]]
        comps.append( comp )
    print nComp
    return nComp, comps

usage = """%prog [options] fileNameIn|- fileNameOut|-

'-' is for stdin, stdout"""

# 17.10.2005
version = open( op.join( init_sfe.install_dir, 'VERSION' ) ).readlines()[0][:-1]

# 05.10.2005, c
# 06.10.2005
# 11.11.2005
# 29.08.2007
def main():
    parser = OptionParser( usage = usage, version = "%prog " + version )
    parser.add_option( "-m", "--mesh",
                       action = "store_true", dest = "saveMesh",
                       default = True,
                       help = "save surface mesh [default: %default]" )
    parser.add_option( "-n", "--no-surface",
                       action = "store_true", dest = "noSurface",
                       default = False,
                       help = "do not output surface [default: %default]" )
    (options, args) = parser.parse_args()

    if (len( args ) == 2):
        fileNameIn = args[0];
        fileNameOut = args[1];
    else:
        parser.print_help(),
        return

    if (fileNameIn == '-'):
        fileIn = sys.stdin
    else:
        fileIn = open( fileNameIn, "r" ); 

    mesh = Mesh.fromFile( fileNameIn )

    if (fileNameIn != '-'):
        fileIn.close()

    domain = Domain.fromMesh( mesh, 'eldesc' )
    domain.setupSubDomains()

    if domain.subs[0].nFace:
        domain.setupNeighbourLists( createEdgeList = False )

        lst, surfFaces = domain.surfaceFaces()

        surfMesh = Mesh.fromSurface( surfFaces, mesh )

        if options.saveMesh:
            surfMesh.write( "surf_" + op.basename( fileNameIn ) )

        if options.noSurface:
            return

        nNod = mesh.nod0.shape[0]
        grS = surfaceGraph( surfFaces, nNod )
##         import sfe.base.plotutils as plu
##         plu.spy( grS )
##         plu.pylab.show()

        nComp, comps = surfaceComponents( grS, surfFaces )
#        print 'components:', nComp

        ccs, comps = comps, nm.zeros( (0,1), nm.int32 )
        for cc in ccs:
            comps = nm.concatenate( (comps, cc[:,nm.newaxis]), 0 )

        out = nm.concatenate( (lst, comps), 1 )

        if (fileNameOut == '-'):
            fileOut = sys.stdout
        else:
            fileOut = open( fileNameOut, "w" ); 
        for row in out:
            fileOut.write( '%d %d %d %d\n' % (row[0], row[1], row[2], row[3]) )
        if (fileNameOut != '-'):
            fileOut.close()
    

if __name__=='__main__':
    main()
