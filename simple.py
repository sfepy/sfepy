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
from __future__ import print_function
from __future__ import absolute_import
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import sfepy
from sfepy.base.base import output
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.applications import PDESolverApp, EVPSolverApp

def print_terms():
    import sfepy.terms as t
    tt = t.term_table
    print('Terms: %d available:' % len(tt))
    print(sorted(tt.keys()))

def print_solvers():
    from sfepy.solvers import solver_table
    print('Solvers: %d available:' % len(solver_table))
    print(sorted(solver_table.keys()))

helps = {
    'debug':
    'automatically start debugger when an exception is raised',
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
    'save_restart' :
    'if given, save restart files according to the given mode.',
    'load_restart' :
    'if given, load the given restart file',
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
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + sfepy.__version__)
    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        default=False, help=helps['debug'])
    parser.add_argument('-c', '--conf', metavar='"key : value, ..."',
                        action='store', dest='conf', type=str,
                        default=None, help= helps['conf'])
    parser.add_argument('-O', '--options', metavar='"key : value, ..."',
                        action='store', dest='app_options', type=str,
                        default=None, help=helps['options'])
    parser.add_argument('-d', '--define', metavar='"key : value, ..."',
                        action='store', dest='define_args', type=str,
                        default=None, help=helps['define'])
    parser.add_argument('-o', metavar='filename',
                        action='store', dest='output_filename_trunk',
                        default=None, help=helps['filename'])
    parser.add_argument('--format', metavar='format',
                        action='store', dest='output_format',
                        default=None, help=helps['output_format'])
    parser.add_argument('--save-restart', metavar='mode', type=int,
                        action='store', dest='save_restart',
                        default=None, help=helps['save_restart'])
    parser.add_argument('--load-restart', metavar='filename',
                        action='store', dest='load_restart',
                        default=None, help=helps['load_restart'])
    parser.add_argument('--log', metavar='file',
                        action='store', dest='log',
                        default=None, help=helps['log'])
    parser.add_argument('-q', '--quiet',
                        action='store_true', dest='quiet',
                        default=False, help=helps['quiet'])
    parser.add_argument('--save-ebc',
                        action='store_true', dest='save_ebc',
                        default=False, help=helps['save_ebc'])
    parser.add_argument('--save-ebc-nodes',
                        action='store_true', dest='save_ebc_nodes',
                        default=False, help=helps['save_ebc_nodes'])
    parser.add_argument('--save-regions',
                        action='store_true', dest='save_regions',
                        default=False, help=helps['save_regions'])
    parser.add_argument('--save-regions-as-groups',
                        action='store_true', dest='save_regions_as_groups',
                        default=False, help=helps['save_regions_as_groups'])
    parser.add_argument('--save-field-meshes',
                        action='store_true', dest='save_field_meshes',
                        default=False, help=helps['save_field_meshes'])
    parser.add_argument('--solve-not',
                        action='store_true', dest='solve_not',
                        default=False, help=helps['solve_not'])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--list', metavar='what',
                        action='store', dest='_list',
                        default=None, help=helps['list'])
    group.add_argument('filename_in', nargs='?')
    options, petsc_opts = parser.parse_known_args()

    if options._list is not None:
        if options._list == 'terms':
            print_terms()

        elif options._list == 'solvers':
            print_solvers()

        return

    if options.debug:
        from sfepy.base.base import debug_on_error; debug_on_error()

    filename_in = options.filename_in
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

    opts.save_restart = options.save_restart
    opts.load_restart = options.load_restart

    if conf.options.get('evps') is None:
        app = PDESolverApp(conf, options, output_prefix)

    else:
        app = EVPSolverApp(conf, options, output_prefix)

    if hasattr(opts, 'parametric_hook'): # Parametric study.
        parametric_hook = conf.get_function(opts.parametric_hook)
        app.parametrize(parametric_hook)
    app()

if __name__ == '__main__':
    main()
