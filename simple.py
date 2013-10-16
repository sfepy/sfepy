#!/usr/bin/env python
# 12.01.2007, c
"""
Solve partial differential equations given in a SfePy problem definition file.

Example problem definition files can be found in ``examples/`` directory of the
SfePy top-level directory. This script works with all the examples except those
in ``examples/standalone/``.

Both normal and parametric study runs are supported. A parametric study allows
repeated runs for varying some of the simulation parameters - see
``examples/diffusion/poisson_parametric_study.py`` file.
"""
from optparse import OptionParser

import sfepy
from sfepy.base.base import output
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.applications import PDESolverApp

def print_terms():
    import sfepy.terms as t
    tt = t.term_table
    print 'Terms: %d available:' % len(tt)
    print sorted(tt.keys())

def print_solvers():
    from sfepy.solvers import solver_table
    print 'Solvers: %d available:' % len(solver_table)
    print sorted(solver_table.keys())

usage = """%prog [options] filename_in\n""" + __doc__.rstrip()

help = {
    'conf' :
    'override problem description file items, written as python'
    ' dictionary without surrounding braces',
    'options' : 'override options item of problem description,'
    ' written as python dictionary without surrounding braces',
    'define' : 'pass given arguments written as python dictionary'
    ' without surrounding braces to define() function of problem description'
    ' file',
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'output_format' :
    'output file format, one of: {vtk, h5} [default: vtk]',
    'log' :
    'log all messages to specified file (existing file will be overwritten!)',
    'quiet' :
    'do not print any messages to screen',
    'save_ebc' :
    'save a zero solution with applied EBCs (Dirichlet boundary conditions)',
    'save_ebc_nodes' :
    'save a zero solution with added non-zeros in EBC (Dirichlet boundary'
    ' conditions) nodes - scalar variables are shown using colors,'
    ' vector variables using arrows with non-zero components corresponding'
    ' to constrained components',
    'save_regions' :
    'save problem regions as meshes',
    'save_regions_as_groups' :
    'save problem regions in a single mesh but mark them by using different'
    ' element/node group numbers',
    'save_field_meshes' :
    'save meshes of problem fields (with extra DOF nodes)',
    'solve_not' :
    'do not solve (use in connection with --save-*)',
    'list' :
    'list data, what can be one of: {terms, solvers}',
}

def main():
    parser = OptionParser(usage=usage, version='%prog ' + sfepy.__version__)
    parser.add_option('-c', '--conf', metavar='"key : value, ..."',
                      action='store', dest='conf', type='string',
                      default=None, help= help['conf'])
    parser.add_option('-O', '--options', metavar='"key : value, ..."',
                      action='store', dest='app_options', type='string',
                      default=None, help=help['options'])
    parser.add_option('-d', '--define', metavar='"key : value, ..."',
                      action='store', dest='define_args', type='string',
                      default=None, help=help['define'])
    parser.add_option('-o', '', metavar='filename',
                      action='store', dest='output_filename_trunk',
                      default=None, help=help['filename'])
    parser.add_option('', '--format', metavar='format',
                      action='store', dest='output_format',
                      default=None, help=help['output_format'])
    parser.add_option('', '--log', metavar='file',
                      action='store', dest='log',
                      default=None, help=help['log'])
    parser.add_option('-q', '--quiet',
                      action='store_true', dest='quiet',
                      default=False, help=help['quiet'])
    parser.add_option('', '--save-ebc',
                      action='store_true', dest='save_ebc',
                      default=False, help=help['save_ebc'])
    parser.add_option('', '--save-ebc-nodes',
                      action='store_true', dest='save_ebc_nodes',
                      default=False, help=help['save_ebc_nodes'])
    parser.add_option('', '--save-regions',
                      action='store_true', dest='save_regions',
                      default=False, help=help['save_regions'])
    parser.add_option('', '--save-regions-as-groups',
                      action='store_true', dest='save_regions_as_groups',
                      default=False, help=help['save_regions_as_groups'])
    parser.add_option('', '--save-field-meshes',
                      action='store_true', dest='save_field_meshes',
                      default=False, help=help['save_field_meshes'])
    parser.add_option('', '--solve-not',
                      action='store_true', dest='solve_not',
                      default=False, help=help['solve_not'])
    parser.add_option('', '--list', metavar='what',
                      action='store', dest='_list',
                      default=None, help=help['list'])

    options, args = parser.parse_args()

    if (len(args) == 1):
        filename_in = args[0];
    else:
        if options._list == 'terms':
            print_terms()

        elif options._list == 'solvers':
            print_solvers()

        else:
            parser.print_help(),
        return

    output.set_output(filename=options.log,
                      quiet=options.quiet,
                      combined=options.log is not None)

    required, other = get_standard_keywords()
    if options.solve_not:
        required.remove('equations')
        required.remove('solver_[0-9]+|solvers')
        other.extend(['equations'])

    conf = ProblemConf.from_file_and_options(filename_in, options,
                                             required, other,
                                             define_args=options.define_args)

    opts = conf.options
    output_prefix = opts.get('output_prefix', 'sfepy:')

    app = PDESolverApp(conf, options, output_prefix)
    if hasattr(opts, 'parametric_hook'): # Parametric study.
        parametric_hook = conf.get_function(opts.parametric_hook)
        app.parametrize(parametric_hook)
    app()

if __name__ == '__main__':
    main()
