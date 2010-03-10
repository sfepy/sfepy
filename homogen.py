#!/usr/bin/env python
# 11.07.2006, c 
from optparse import OptionParser

import sfepy
from sfepy.base.base import *
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.homogenization.homogen_app import HomogenizationApp

usage = """%prog [options] filename_in"""

help = {
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
}

def main():
    parser = OptionParser(usage = usage, version = "%prog " + sfepy.__version__)
    parser.add_option( "-o", "", metavar = 'filename',
                       action = "store", dest = "output_filename_trunk",
                       default = None, help = help['filename'] )
    (options, args) = parser.parse_args()

    if (len( args ) == 1):
        filename_in = args[0];
    else:
        parser.print_help(),
        return

    required, other = get_standard_keywords()
    required.remove( 'equations' )

    conf = ProblemConf.from_file( filename_in, required, other )

    app = HomogenizationApp( conf, options, 'homogen:' )
    opts = conf.options
    if hasattr( opts, 'parametric_hook' ): # Parametric study.
        parametric_hook = getattr( conf, opts.parametric_hook )
        app.parametrize( parametric_hook )
    app()

if __name__ == '__main__':
    main()
