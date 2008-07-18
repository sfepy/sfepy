#!/usr/bin/env python
# 12.01.2007, c 
import os.path as op
from optparse import OptionParser

import init_sfepy
from sfepy.base.base import *
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.solvers.generic import solve_direct

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

usage = """%prog [options] file_name_in"""

help = {
    'file_name' :
    'basename of output file(s) [default: <basename of input file>]',
    'dump' :
    "dump problem state [default: %default]",
    'save_ebc' :
    "save problem state showing EBC (Dirichlet conditions) [default: %default]",
    'save_regions' :
    "save problem regions as meshes [default: %default]",
    'save_field_meshes' :
    "save meshes of problem fields (with extra DOF nodes) [default: %default]",
    'save_region_field_meshes' :
    "save meshes of regions of problem fields (with extra DOF nodes)"
    "[default: %default]",
    'solve_not' :
    "do not solve (use in connection with --save-*) [default: %default]",
    'list' :
    "list data according to what, what can be one of: terms",
}

##
# c: 12.01.2007, r: 02.04.2008
def main():
    version = open( op.join( init_sfepy.install_dir,
                             'VERSION' ) ).readlines()[0][:-1]

    parser = OptionParser( usage = usage, version = "%prog " + version )
    parser.add_option( "-o", "", metavar = 'file_name',
                       action = "store", dest = "output_file_name_trunk",
                       default = None, help = help['file_name'] )
    parser.add_option( "", "--dump",
                       action = "store_true", dest = "dump",
                       default = False, help = help['dump'] )
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
        file_name_in = args[0];
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

    conf = ProblemConf.from_file( file_name_in, required, other )
##     print conf
##     pause()

    opts = conf.options
    if hasattr( opts, 'output_prefix' ):
        set_output_prefix( opts.output_prefix )

    dpb, vec_dp, data = solve_direct( conf, options )

    if hasattr( opts, 'post_process_hook_final' ): # User postprocessing.
        hook = getattr( conf, opts.post_process_hook_final )
        hook( dpb, vec_dp, data )

if __name__ == '__main__':
    main()
