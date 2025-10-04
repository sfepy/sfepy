#!/usr/bin/env python
# 12.01.2007, c
"""
Solve partial differential equations given in a SfePy problem definition file.

Example problem definition files can be found in ``sfepy/examples/`` directory
of the SfePy top-level directory.

In the examples below it is supposed that sfepy is installed. When using the
in-place build, replace ``sfepy-run`` by ``python3 sfepy/scripts/simple.py``.

The supported application kinds (--app option) are:

- bvp - boundary value problem. Example::

    sfepy-run sfepy/examples/diffusion/poisson.py

- homogen - calculation of local microscopic problems (correctors) and
  homogenized coefficients. Example::

    sfepy-run sfepy/examples/homogenization/perfusion_micro.py

- bvp-mM - micro-macro boundary value problem. Solve a coupled two-scale
  problem in parallel using MPI. One computational node is solving a
  macroscopic equation while the others are solving local microscopic problems
  and homogenized coefficients. The --app option is required in this case.
  Example::

    mpiexec -n 4 sfepy-run --app=bvp-mM --debug-mpi sfepy/examples/homogenization/nonlinear_hyperelastic_mM.py

- evp - eigenvalue problem. Example::

    sfepy-run sfepy/examples/quantum/well.py

- phonon - phononic band gaps. Example::

    sfepy-run sfepy/examples/phononic/band_gaps.py --phonon-plot

Both normal and parametric study runs are supported. A parametric study allows
repeated runs for varying some of the simulation parameters - see
``sfepy/examples/diffusion/poisson_parametric_study.py`` file.
"""
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import sfepy
from sfepy.base.base import output, Struct
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
    'app' :
    'override application kind, normally determined automatically.'
    ' The supported kinds are:'
    ' bvp (boundary value problem),'
    ' homogen (correctors, homogenized coefficients),'
    ' bvp-mM (micro-macro boundary value problem,'
    '         homogenized coefficients computed in parallel using MPI),'
    ' evp (eigenvalue problem),'
    ' phonon (phononic band gaps)',
    'debug':
    'automatically start debugger when an exception is raised',
    'debug_mpi': 'log MPI communication (mM mode only)',
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
    'solve_not' :
    'do not solve (use in connection with --save-*)',
    'detect_band_gaps' :
    'detect frequency band gaps',
    'analyze_dispersion' :
    'analyze dispersion properties (low frequency domain)',
    'plot' :
    'plot frequency band gaps, assumes -b',
    'phase_velocity' :
    'compute phase velocity (frequency-independent mass only)',
    'list' :
    'list data, what can be one of: {terms, solvers}',
}

def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + sfepy.__version__)
    parser.add_argument('-a', '--app', action='store', dest='app',
                        choices=['bvp', 'homogen', 'bvp-mM', 'evp', 'phonon'],
                        default=None, help= helps['app'])
    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        default=False, help=helps['debug'])
    parser.add_argument('--debug-mpi',
                        action='store_true', dest='debug_mpi',
                        default=False, help=helps['debug_mpi'])
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
    parser.add_argument('--solve-not',
                        action='store_true', dest='solve_not',
                        default=False, help=helps['solve_not'])
    parser.add_argument('--phonon-band-gaps',
                        action='store_true', dest='detect_band_gaps',
                        default=False, help=helps['detect_band_gaps'])
    parser.add_argument('--phonon-dispersion',
                        action='store_true', dest='analyze_dispersion',
                        default=False, help=helps['analyze_dispersion'])
    parser.add_argument('--phonon-plot',
                        action='store_true', dest='plot',
                        default=False, help=helps['plot'])
    parser.add_argument('--phonon-phase-velocity',
                        action='store_true', dest='phase_velocity',
                        default=False, help=helps['phase_velocity'])
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

    if not (options.analyze_dispersion or options.detect_band_gaps):
            options.plot = False

    if options.debug:
        from sfepy.base.base import debug_on_error; debug_on_error()

    filename_in = options.filename_in
    output.set_output(filename=options.log,
                      quiet=options.quiet,
                      combined=options.log is not None)

    required, other = get_standard_keywords()
    required.remove('equations')
    if options.solve_not:
        required.remove('solver_[0-9]+|solvers')
        other.extend(['equations'])
    if options.detect_band_gaps and not options.analyze_dispersion:
        if 'solver_[0-9]+|solvers' in required:
            required.remove('solver_[0-9]+|solvers')
    if options.phase_velocity:
        required = [ii for ii in required if 'ebc' not in ii]

    conf = ProblemConf.from_file_and_options(filename_in, options,
                                             required, other,
                                             define_args=options.define_args)
    if conf.options.get('coefs') is None:
        if conf.get('equations') is None:
            ValueError('required missing: equations')

    app_mode = options.app
    if app_mode is None:
        if conf.options.get('coefs') is not None:
            if 'band_gaps' in conf.get(conf.options.coefs).keys():
                app_mode = 'phonon'

            else:
                app_mode = 'homogen'

        elif conf.options.get('evps') is not None:
            app_mode = 'evp'

        else:
            app_mode = 'bvp'

    opts = conf.options

    opts.save_restart = options.save_restart
    opts.load_restart = options.load_restart

    if app_mode == 'bvp':
        output_prefix = opts.get('output_prefix', 'sfepy:')
        app = PDESolverApp(conf, options, output_prefix)

    elif app_mode == 'homogen':
        from sfepy.homogenization.homogen_app import HomogenizationApp

        output_prefix = opts.get('output_prefix', 'homogen:')
        app = HomogenizationApp(conf, options, output_prefix)

    elif app_mode == 'bvp-mM':
        import sfepy.base.multiproc_mpi as multi_mpi

        if options.debug_mpi:
            multi_mpi.set_logging_level('debug')

        if multi_mpi.mpi_rank == multi_mpi.mpi_master:
            nslaves = multi_mpi.cpu_count() - 1
            opts.n_mpi_homog_slaves = nslaves
            output_prefix = opts.get('output_prefix', 'sfepy:')

            app = PDESolverApp(conf, options, output_prefix)
            if hasattr(opts, 'parametric_hook'):  # Parametric study.
                parametric_hook = conf.get_function(opts.parametric_hook)
                app.parametrize(parametric_hook)
            app()

            multi_mpi.master_send_task('finalize', None)
            return

        else:
            # MPI slave mode - calculate homogenized coefficients
            homogen_app = None
            done = False
            rank = multi_mpi.mpi_rank
            while not done:
                task, data = multi_mpi.slave_get_task('main slave loop')

                if task == 'init':  # data: micro_file, n_micro
                    output.set_output(filename='homog_app_mpi_%d.log' % rank,
                                      quiet=True)
                    micro_file, n_micro = data[:2]
                    required, other = get_standard_keywords()
                    required.remove('equations')
                    conf = ProblemConf.from_file(micro_file, required, other,
                                                 verbose=False)
                    options = Struct(output_filename_trunk=None)
                    homogen_app = HomogenizationApp(conf, options, 'micro:',
                                                    n_micro=n_micro)
                elif task == 'calculate':  # data: rel_def_grad, ts, iteration
                    macro_data, ts, iteration = data[:3]
                    homogen_app.setup_macro_data(macro_data)
                    homogen_app(ret_all=True, itime=ts.step, iiter=iteration)
                elif task == 'finalize':
                    done = True
            return

    elif app_mode == 'evp':
        output_prefix = opts.get('output_prefix', 'sfepy:')
        app = EVPSolverApp(conf, options, output_prefix)

    elif app_mode == 'phonon':
        from sfepy.homogenization.band_gaps_app import AcousticBandGapsApp

        output_prefix = opts.get('output_prefix', 'phonon:')
        app = AcousticBandGapsApp(conf, options, output_prefix)

    if hasattr(opts, 'parametric_hook'): # Parametric study.
        parametric_hook = conf.get_function(opts.parametric_hook)
        app.parametrize(parametric_hook)
    app()

if __name__ == '__main__':
    main()
