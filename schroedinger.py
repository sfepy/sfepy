#!/usr/bin/env python
# 12.01.2007, c 
import os.path as op
from optparse import OptionParser

import init_sfe
from sfe.base.base import *
from sfe.base.conf import ProblemConf
from sfe.base.la import eig
from sfe.fem.evaluate import evalTermOP
import sfe.base.ioutils as io
from sfe.fem.problemDef import ProblemDefinition
from sfe.homogenization.phono import processOptions

##
# c: 01.02.2008, r: 08.02.2008
def solveEigenProblem( conf, options ):

    if options.outputFileNameTrunk:
        ofnTrunk = options.outputFileNameTrunk
    else:
        ofnTrunk = io.getTrunk( conf.fileName_mesh )

    pb = ProblemDefinition.fromConf( conf )
    dim = pb.domain.mesh.dim

    pb.timeUpdate()

    dummy = pb.createStateVector()
    mtxA = evalTermOP( dummy, conf.equations['lhs'], pb,
                       dwMode = 'matrix', tangentMatrix = pb.mtxA )
    mtxB = evalTermOP( dummy, conf.equations['rhs'], pb,
                       dwMode = 'matrix', tangentMatrix = pb.mtxA.copy() )

##     mtxA.save( 'a.txt', format='%d %d %.12f\n' )
##     mtxB.save( 'b.txt', format='%d %d %.12f\n' )
    print 'computing resonance frequencies...'
    tt = [0]
    eigs, mtxSPhi = eig( mtxA.toarray(), mtxB.toarray(), returnTime = tt )
    print 'done in %.2f s' % tt[0]
    print eigs
##     import sfe.base.plotutils as plu
##     plu.spy( mtxB, eps = 1e-12 )
##     plu.pylab.show()
##     pause()
    nEigs = eigs.shape[0]
    opts = processOptions( conf.options, nEigs )

    mtxPhi = nm.empty( (pb.variables.di.ptr[-1], mtxSPhi.shape[1]),
                       dtype = nm.float64 )
    for ii in xrange( nEigs ):
        mtxPhi[:,ii] = pb.variables.makeFullVec( mtxSPhi[:,ii] )

    out = {}
    for ii in xrange( nEigs ):
        if (ii > opts.save[0]) and (ii < (nEigs - opts.save[1])): continue
        aux = pb.stateToOutput( mtxPhi[:,ii] )
        key = aux.keys()[0]
        out[key+'%03d' % ii] = aux[key]

    pb.domain.mesh.write( ofnTrunk + '.vtk', io = 'auto', out = out )

    fd = open( ofnTrunk + '_eigs.txt', 'w' )
    eigs.tofile( fd, ' ' )
    fd.close()

    return Struct( pb = pb, eigs = eigs, mtxPhi = mtxPhi )


usage = """%prog [options] fileNameIn"""

help = {
    'fileName' :
    'basename of output file(s) [default: <basename of input file>]',
}

##
# c: 01.02.2008, r: 01.02.2008
def main():
    version = open( op.join( init_sfe.install_dir,
                             'VERSION' ) ).readlines()[0][:-1]

    parser = OptionParser( usage = usage, version = "%prog " + version )
    parser.add_option( "-o", "", metavar = 'fileName',
                       action = "store", dest = "outputFileNameTrunk",
                       default = None, help = help['fileName'] )

    options, args = parser.parse_args()

    if (len( args ) == 1):
        fileNameIn = args[0];
    else:
        parser.print_help(),
        return
    
    required = ['fileName_mesh', 'field_[0-9]+', 'ebc|nbc', 'fe', 'equations',
                'region_[0-9]+|regions', 'integral_[0-9]+', 'variables',
                'material_[0-9]+']
    other = ['functions', 'modules', 'epbc', 'lcbc']

    conf = ProblemConf.fromFile( fileNameIn, required, other )
##     print conf
##     pause()

    evp = solveEigenProblem( conf, options )

if __name__ == '__main__':
    main()
