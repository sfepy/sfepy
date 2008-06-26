#!/usr/bin/env python
# 26.09.2006, c 
import os.path as op
from optparse import OptionParser

import init_sfepy
from sfepy.base.base import *
from sfepy.fem.mesh import Mesh
from sfepy.fem.meshio import HDF5MeshIO
from sfepy.solvers.ts import TimeStepper
from sfepy.base.ioutils import getTrunk, writeDictHDF5

##
# c: 26.09.2006, r: 23.06.2008
def dumpToVTK( fileName, options, steps = None ):
    output( 'dumping to VTK...' )
    
    mesh = Mesh.fromFile( fileName )

    io = HDF5MeshIO( fileName )
    ts = TimeStepper( *io.readTimeStepper() )
    nDigit = int( nm.log10( ts.nStep - 1 ) + 1 )
    format = '%%%dd of %%%dd' % (nDigit, nDigit)
    suffix = '.%%0%dd.vtk' % nDigit

    if options.outputFileNameTrunk:
        ofnTrunk = options.outputFileNameTrunk
    else:
        ofnTrunk = getTrunk( fileName )

    if steps is None:
        iterator = ts.iterFrom( options.step0 )
    else:
        iterator = [(step, ts.times[step]) for step in steps]

    for step, time in iterator:
        output( format % (step, ts.nStep - 1) )
        out = io.readData( step )
        if out is None: break
        mesh.write( ofnTrunk + suffix % step, io = 'auto', out = out )

    output( '...done' )
    return suffix

##
# c: 26.09.2006, r: 23.06.2008
def extractTimeHistory( fileName, options ):
    output( 'extracting selected data...' )

    el = options.extractList
    output( 'extraction list:', el )

    ##
    # Parse extractions.
    pes = OneTypeList( Struct )
    for chunk in el.split( ',' ):
        aux =  chunk.strip().split()
        pes.append( Struct( var = aux[0],
                            mode = aux[1],
                            indx = map( int, aux[2:] ),
                            igs = None ) )

    ##
    # Verify array limits, set igs for element data, shift indx.
    mesh = Mesh.fromFile( fileName )
    nEl, nEls, offs = mesh.nEl, mesh.nEls, mesh.elOffsets
    for pe in pes:
        if pe.mode == 'n':
            for ii in pe.indx:
                if (ii < 0) or (ii >= mesh.nNod):
                    raise IndexError, 'node index 0 <= %d < %d'\
                          % (ii, mesh.nNod)

        if pe.mode == 'e':
            pe.igs = []
            for ii, ie in enumerate( pe.indx[:] ):
                if (ie < 0) or (ie >= nEl):
                    raise IndexError, 'element index 0 <= %d < %d'\
                          % (ie, nEl)
                ig = (ie < nEls).argmax()
                pe.igs.append( ig )
                pe.indx[ii] = ie - offs[ig]

##     print pes

    ##
    # Extract data.
    # Assumes only one element group (ignores igs)!
    io = HDF5MeshIO( fileName )
    ths = {}
    for pe in pes:
        mode, nname = io.readDataHeader( pe.var )
        print mode, nname
        if ((pe.mode == 'n' and mode == 'vertex') or
            (pe.mode == 'e' and mode == 'cell')):
            th = io.readTimeHistory( nname, pe.indx )

        elif pe.mode == 'e' and mode == 'vertex':
            conn = mesh.conns[0]
            th = {}
            for iel in pe.indx:
                ips = conn[iel]
                th[iel] = io.readTimeHistory( nname, ips )
        else:
            raise RuntimeError, 'cannot extract cell data %s in nodes' % pe.var
            
        ths[pe.var] = th
    return ths

##
# 27.09.2006, c
def averageVertexVarInCells( thsIn ):
    ths = dict.fromkeys( thsIn.keys() )
    for var, th in thsIn.iteritems():
        aux = dict.fromkeys( th.keys() )
        for ir, data in th.iteritems():
            if isinstance( data, dict ):
                for ic, ndata in data.iteritems():
                    if aux[ir] is None:
                        aux[ir] = ndata
                    else:
                        aux[ir] += ndata
                aux[ir] /= float( len( data ) )
            else:
                aux[ir] = data
        ths[var] = aux

    return ths

usage = """%prog [options] fileNameIn"""

help = {
    'fileName' :
    'basename of output file(s) [default: <basename of input file>]',
    'dump' :
    'dump to sequence of VTK files',
    'from' :
    'start dumping from step ii [default: %default]',
    'extract' :
    'extract variables according to extraction list',
    'average' :
    'average vertex variable into cells ("e" extraction mode)'
}

##
# c: 26.09.2006, r: 23.06.2008
def main():
    version = open( op.join( init_sfepy.install_dir,
                             'VERSION' ) ).readlines()[0][:-1]

    parser = OptionParser( usage = usage, version = "%prog " + version )
    parser.add_option( "-o", "", metavar = 'fileName',
                       action = "store", dest = "outputFileNameTrunk",
                       default = None, help = help['fileName'] )
    parser.add_option( "-d", "--dump",
                       action = "store_true", dest = "dump",
                       default = False, help = help['dump'] )
    parser.add_option( "-f", "--from", type = int, metavar = 'ii',
                       action = "store", dest = "step0",
                       default = 0, help = help['from'] )
    parser.add_option( "-e", "--extract", metavar = 'list',
                       action = "store", dest = "extractList",
                       default = None, help = help['extract'] )
    parser.add_option( "-a", "--average",
                       action = "store_true", dest = "average",
                       default = False, help = help['average'] )

    (options, args) = parser.parse_args()
    print options
    print args
    if (len( args ) == 1):
        fileNameIn = args[0];
    else:
        parser.print_help(),
        return

    if options.dump:
        dumpToVTK( fileNameIn, options )

    if options.extractList:
        ths = extractTimeHistory( fileNameIn, options )
##         print ths

        if options.average:
            ths = averageVertexVarInCells( ths )
##             print ths

        if options.outputFileNameTrunk:
            ts = TimeStepper( *HDF5MeshIO( fileNameIn ).readTimeStepper() )
            ths.update( {'times' : ts.times, 'dt' : ts.dt} )
            writeDictHDF5( options.outputFileNameTrunk + '.h5', ths )

if __name__ == '__main__':
    main()
