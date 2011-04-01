#!/usr/bin/env python
# 26.09.2006, c 
"""
Examples
--------

$ ./extractor.py -e "p e 0 1999" bone.h5
$ ./extractor.py -e "p e 0 1999" bone.h5 -a
$ ./extractor.py -e "p e 0 1999" bone.h5 -o extracted.h5
$ ./extractor.py -e "p e 0 1999" bone.h5 -o extracted.h5 -a
"""
import os
from optparse import OptionParser

import sfepy
from sfepy.base.base import nm, dict_to_struct, get_default
from sfepy.base.ioutils import get_trunk
import sfepy.postprocess.time_history as th

usage = """%prog [options] filename_in

Extract information from a SfePy multi-time-step results file (HDF5 format).
"""

help = {
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'dump' :
    'dump to sequence of VTK files',
    'same_dir' :
    'store the dumped VTK files in the directory of filename_in',
    'from' :
    'start dumping from time step ii [default: %default]',
    'to' :
    'stop dumping at time step ii [default: <last step>]',
    'step' :
    'use every ii-th step for dumping [default: %default]',
    'extract' :
    'extract variables according to extraction list',
    'average' :
    'average vertex variable into cells ("e" extraction mode)'
}

##
# c: 26.09.2006, r: 23.06.2008
def main():
    parser = OptionParser(usage = usage, version = "%prog " + sfepy.__version__)
    parser.add_option( "-o", "", metavar = 'filename',
                       action = "store", dest = "output_filename_trunk",
                       default = None, help = help['filename'] )
    parser.add_option( "-d", "--dump",
                       action = "store_true", dest = "dump",
                       default = False, help = help['dump'] )
    parser.add_option( "", "--same-dir",
                       action = "store_true", dest = "same_dir",
                       default = False, help = help['same_dir'] )
    parser.add_option( "-f", "--from", type = int, metavar = 'ii',
                       action = "store", dest = "step_from",
                       default = 0, help = help['from'] )
    parser.add_option( "-t", "--to", type = int, metavar = 'ii',
                       action = "store", dest = "step_to",
                       default = None, help = help['to'] )
    parser.add_option( "-s", "--step", type = int, metavar = 'ii',
                       action = "store", dest = "step_by",
                       default = 1, help = help['step'] )
    parser.add_option( "-e", "--extract", metavar = 'list',
                       action = "store", dest = "extract",
                       default = None, help = help['extract'] )
    parser.add_option( "-a", "--average",
                       action = "store_true", dest = "average",
                       default = False, help = help['average'] )

    (options, args) = parser.parse_args()

    if (len( args ) == 1):
        filename_in = args[0];
    else:
        parser.print_help(),
        return

    if options.dump:
        trunk = get_default(options.output_filename_trunk,
                            get_trunk(filename_in))
        if options.same_dir:
            trunk = os.path.join(os.path.dirname(filename_in),
                                 os.path.basename(trunk))
        
        if options.step_to is None:
            th.dump_to_vtk(filename_in,
                           output_filename_trunk=trunk,
                           step0=options.step_from)

        else:
            th.dump_to_vtk(filename_in,
                           output_filename_trunk=trunk,
                           steps=nm.arange(options.step_from,
                                           options.step_to + 1,
                                           options.step_by, dtype=nm.int))

    if options.extract:
        ths, ts = th.extract_time_history(filename_in, options.extract)
##         print ths

        if options.average:
            ths = th.average_vertex_var_in_cells( ths )
##             print ths

        if options.output_filename_trunk:
            th.save_time_history(ths, ts, options.output_filename_trunk + '.h5')

        else:
            print dict_to_struct(ths, flag=(1, 1, 1)).str_all()

if __name__ == '__main__':
    main()
