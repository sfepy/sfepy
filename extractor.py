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
from optparse import OptionParser

import sfepy
from sfepy.base.base import dict_to_struct
from sfepy.postprocess.time_history \
     import dump_to_vtk, extract_time_history, average_vertex_var_in_cells, \
            save_time_history

usage = """%prog [options] filename_in

Extract information from a SfePy multi-time-step results file (HDF5 format).
"""

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
    parser = OptionParser(usage = usage, version = "%prog " + sfepy.__version__)
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
        dump_to_vtk(filename_in, ofn_trunk=options.output_filename_trunk,
                    step0=options.step0)

    if options.extract:
        ths, ts = extract_time_history(filename_in, options.extract)
##         print ths

        if options.average:
            ths = average_vertex_var_in_cells( ths )
##             print ths

        if options.output_filename_trunk:
            save_time_history(ths, ts, options.output_filename_trunk + '.h5')

        else:
            print dict_to_struct(ths, flag=(1, 1, 1)).str_all()

if __name__ == '__main__':
    main()
