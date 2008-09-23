#!/usr/bin/env python
# 26.09.2006, c 
import os.path as op
from optparse import OptionParser

import init_sfepy
from sfepy.base.base import *
from sfepy.fem.mesh import Mesh
from sfepy.fem.meshio import HDF5MeshIO
from sfepy.solvers.ts import TimeStepper
from sfepy.base.ioutils import get_trunk, write_dict_hdf5

##
# c: 26.09.2006, r: 09.07.2008
def dump_to_vtk( filename, options, steps = None ):
    output( 'dumping to VTK...' )
    
    mesh = Mesh.from_file( filename )

    io = HDF5MeshIO( filename )
    ts = TimeStepper( *io.read_time_stepper() )

    if options.output_filename_trunk:
        ofn_trunk = options.output_filename_trunk
    else:
        ofn_trunk = get_trunk( filename )

    if steps is None:
        iterator = ts.iter_from( options.step0 )
    else:
        iterator = [(step, ts.times[step]) for step in steps]

    for step, time in iterator:
        output( ts.format % (step, ts.n_step - 1) )
        out = io.read_data( step )
        if out is None: break
        mesh.write( ofn_trunk + ts.suffix % step + '.vtk',
                    io = 'auto', out = out )

    output( '...done' )
    return ts.suffix

##
# c: 26.09.2006, r: 23.06.2008
def extract_time_history( filename, options ):
    output( 'extracting selected data...' )

    el = options.extract_list
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
    mesh = Mesh.from_file( filename )
    n_el, n_els, offs = mesh.n_el, mesh.n_els, mesh.el_offsets
    for pe in pes:
        if pe.mode == 'n':
            for ii in pe.indx:
                if (ii < 0) or (ii >= mesh.n_nod):
                    raise IndexError, 'node index 0 <= %d < %d'\
                          % (ii, mesh.n_nod)

        if pe.mode == 'e':
            pe.igs = []
            for ii, ie in enumerate( pe.indx[:] ):
                if (ie < 0) or (ie >= n_el):
                    raise IndexError, 'element index 0 <= %d < %d'\
                          % (ie, n_el)
                ig = (ie < n_els).argmax()
                pe.igs.append( ig )
                pe.indx[ii] = ie - offs[ig]

##     print pes

    ##
    # Extract data.
    # Assumes only one element group (ignores igs)!
    io = HDF5MeshIO( filename )
    ths = {}
    for pe in pes:
        mode, nname = io.read_data_header( pe.var )
        print mode, nname
        if ((pe.mode == 'n' and mode == 'vertex') or
            (pe.mode == 'e' and mode == 'cell')):
            th = io.read_time_history( nname, pe.indx )

        elif pe.mode == 'e' and mode == 'vertex':
            conn = mesh.conns[0]
            th = {}
            for iel in pe.indx:
                ips = conn[iel]
                th[iel] = io.read_time_history( nname, ips )
        else:
            raise RuntimeError, 'cannot extract cell data %s in nodes' % pe.var
            
        ths[pe.var] = th
    return ths

##
# 27.09.2006, c
def average_vertex_var_in_cells( ths_in ):
    ths = dict.fromkeys( ths_in.keys() )
    for var, th in ths_in.iteritems():
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

usage = """%prog [options] filename_in"""

help = {
    'filename' :
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
    parser.add_option( "-o", "", metavar = 'filename',
                       action = "store", dest = "output_filename_trunk",
                       default = None, help = help['filename'] )
    parser.add_option( "-d", "--dump",
                       action = "store_true", dest = "dump",
                       default = False, help = help['dump'] )
    parser.add_option( "-f", "--from", type = int, metavar = 'ii',
                       action = "store", dest = "step0",
                       default = 0, help = help['from'] )
    parser.add_option( "-e", "--extract", metavar = 'list',
                       action = "store", dest = "extract_list",
                       default = None, help = help['extract'] )
    parser.add_option( "-a", "--average",
                       action = "store_true", dest = "average",
                       default = False, help = help['average'] )

    (options, args) = parser.parse_args()
    print options
    print args
    if (len( args ) == 1):
        filename_in = args[0];
    else:
        parser.print_help(),
        return

    if options.dump:
        dump_to_vtk( filename_in, options )

    if options.extract_list:
        ths = extract_time_history( filename_in, options )
##         print ths

        if options.average:
            ths = average_vertex_var_in_cells( ths )
##             print ths

        if options.output_filename_trunk:
            ts = TimeStepper( *HDF5MeshIO( filename_in ).read_time_stepper() )
            ths.update( {'times' : ts.times, 'dt' : ts.dt} )
            write_dict_hdf5( options.output_filename_trunk + '.h5', ths )

if __name__ == '__main__':
    main()
