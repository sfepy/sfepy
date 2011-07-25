#!/usr/bin/env python
# 12.01.2007, c 
from optparse import OptionParser

import sfepy
from sfepy.base.base import output, get_default_attr
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
    'conf' :
    'override problem description file items, written as python'
    ' dictionary without surrouding braces',
    'options' : 'override options item of problem description,'
    ' written as python dictionary without surrouding braces',
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'output_format' :
    'output file format, one of: {vtk, h5, mesh} [default: vtk]',
    'log' :
    "log all messages to specified file (existing file will be overwritten!)",
    'quiet' :
    "do not print any messages to screen",
    'save_ebc' :
    "save problem state showing EBC (Dirichlet conditions)",
    'save_regions' :
    "save problem regions as meshes",
    'save_regions_as_groups' :
    "save problem regions in a single mesh but mark them by using different"
    " element/node group numbers",
    'save_field_meshes' :
    "save meshes of problem fields (with extra DOF nodes)",
    'solve_not' :
    "do not solve (use in connection with --save-*)",
    'list' :
    "list data, what can be one of: {terms}",
}

def main():
    parser = OptionParser(usage = usage, version = "%prog " + sfepy.__version__)
    parser.add_option('-c', '--conf', metavar='"key : value, ..."',
                      action='store', dest='conf', type='string',
                      default=None, help= help['conf'])
    parser.add_option('-O', '--options', metavar='"key : value, ..."',
                      action='store', dest='app_options', type='string',
                      default=None, help=help['options'])
    parser.add_option( "-o", "", metavar = 'filename',
                       action = "store", dest = "output_filename_trunk",
                       default = None, help = help['filename'] )
    parser.add_option( "", "--format", metavar = 'format',
                       action = "store", dest = "output_format",
                       default = None, help = help['output_format'] )
    parser.add_option( "", "--log", metavar = 'file',
                       action = "store", dest = "log",
                       default = None, help = help['log'] )
    parser.add_option( "-q", "--quiet",
                       action = "store_true", dest = "quiet",
                       default = False, help = help['quiet'] )
    parser.add_option( "", "--save-ebc",
                       action = "store_true", dest = "save_ebc",
                       default = False, help = help['save_ebc'] )
    parser.add_option( "", "--save-regions",
                       action = "store_true", dest = "save_regions",
                       default = False, help = help['save_regions'] )
    parser.add_option( "", "--save-regions-as-groups",
                       action = "store_true", dest = "save_regions_as_groups",
                       default = False, help = help['save_regions_as_groups'] )
    parser.add_option( "", "--save-field-meshes",
                       action = "store_true", dest = "save_field_meshes",
                       default = False, help = help['save_field_meshes'] )
    parser.add_option( "", "--solve-not",
                       action = "store_true", dest = "solve_not",
                       default = False, help = help['solve_not'] )
    parser.add_option( "", "--list", metavar = 'what',
                       action = "store", dest = "_list",
                       default = None, help = help['list'] )

    options, args = parser.parse_args()

    if (len( args ) == 1):
        filename_in = args[0];
    else:
        if options._list == 'terms':
            print_terms()
        else:
            parser.print_help(),
        return

    output.set_output(filename=options.log,
                      quiet=options.quiet,
                      combined=options.log is not None)

    required, other = get_standard_keywords()
    if options.solve_not:
        required.remove( 'equations' )
        required.remove( 'solver_[0-9]+|solvers' )
        other.extend( ['equations'] )

    conf = ProblemConf.from_file_and_options(filename_in, options,
                                             required, other)

    opts = conf.options
    output_prefix = get_default_attr( opts, 'output_prefix', 'sfepy:' )

    app = SimpleApp( conf, options, output_prefix )
    if hasattr( opts, 'parametric_hook' ): # Parametric study.
        parametric_hook = getattr( conf, opts.parametric_hook )
        app.parametrize( parametric_hook )
    app()

if __name__ == '__main__':
    main()
