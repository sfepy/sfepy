#!/usr/bin/env python
# 12.01.2007, c 
import os.path as op
import shutil
from optparse import OptionParser

import sfepy
from sfepy.base.base import *
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.applications import SimpleApp

##
# 26.03.2007, c
def print_terms():
    import sfepy.terms as t
    tt = t.term_table
    ct = t.cache_table
    print 'Terms: %d available:' % len( tt )
    print sorted( tt.keys() )
    print 'Term caches: %d available:' % len( ct )
    print sorted( ct.keys() )

usage = """%prog [options] filename_in"""

help = {
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'output_format' :
    'output file format, one of: {vtk, h5, mesh} [default: %default]',
    'save_ebc' :
    "save problem state showing EBC (Dirichlet conditions)",
    'save_regions' :
    "save problem regions as meshes",
    'save_field_meshes' :
    "save meshes of problem fields (with extra DOF nodes)",
    'save_region_field_meshes' :
    "save meshes of regions of problem fields (with extra DOF nodes)",
    'solve_not' :
    "do not solve (use in connection with --save-*)",
    'list' :
    "list data, what can be one of: {terms}",
}

def main():
    parser = OptionParser(usage = usage, version = "%prog " + sfepy.__version__)
    parser.add_option( "-o", "", metavar = 'filename',
                       action = "store", dest = "output_filename_trunk",
                       default = None, help = help['filename'] )
    parser.add_option( "", "--format", metavar = 'format',
                       action = "store", dest = "output_format",
                       default = "vtk", help = help['output_format'] )
    parser.add_option( "", "--save-ebc",
                       action = "store_true", dest = "save_ebc",
                       default = False, help = help['save_ebc'] )
    parser.add_option( "", "--save-regions",
                       action = "store_true", dest = "save_regions",
                       default = False, help = help['save_regions'] )
    parser.add_option( "", "--save-field-meshes",
                       action = "store_true", dest = "save_field_meshes",
                       default = False, help = help['save_field_meshes'] )
    parser.add_option( "", "--save-region-field-meshes",
                       action = "store_true", dest = "save_region_field_meshes",
                       default = False, help = help['save_region_field_meshes'] )
    parser.add_option( "", "--solve-not",
                       action = "store_true", dest = "solve_not",
                       default = False, help = help['solve_not'] )
    parser.add_option( "", "--list", metavar = 'what',
                       action = "store", dest = "_list",
                       default = None, help = help['list'] )

    options, args = parser.parse_args()
#    print options; pause()

    if (len( args ) == 1):
        filename_in = args[0];
    else:
        if options._list == 'terms':
            print_terms()
        else:
            parser.print_help(),
        return
    
    required, other = get_standard_keywords()
    if options.solve_not:
        required.remove( 'equations' )
        required.remove( 'solver_[0-9]+|solvers' )
        other.extend( ['equations'] )

    conf = ProblemConf.from_file( filename_in, required, other )
    opts = conf.options
    output_prefix = get_default_attr( opts, 'output_prefix', 'sfepy:' )

    app = SimpleApp( conf, options, output_prefix )
    if hasattr( opts, 'parametric_hook' ): # Parametric study.
        parametric_hook = getattr( conf, opts.parametric_hook )
        app.parametrize( parametric_hook )
    app()

if __name__ == '__main__':
    main()
