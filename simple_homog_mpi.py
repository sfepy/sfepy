#!/usr/bin/env python
# 06.06.2017, c
"""
Solve a coupled two-scale problem in parallel. One computational node is solving
a macroscopic equation while the others are solving local microscopic problems
and homogenized coefficients.

Run this script as::

    mpiexec -n 4 simple_homog_mpi.py examples/homogenization/nonlinear_hyperelastic_mM.py

"""
from __future__ import print_function
from __future__ import absolute_import
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from sfepy.base.base import output, Struct
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.applications import PDESolverApp
from sfepy.homogenization.homogen_app import HomogenizationApp
import simple

import sfepy.base.multiproc_mpi as multi_mpi


helps = {k: simple.helps[k] for k in ['debug', 'conf', 'options', 'define',
                                      'filename', 'output_format', 'log',
                                      'quiet']}
helps.update({
    'debug_mpi': 'log MPI communication',
})


def main():
    # if multi_mpi.cpu_count() < 2:
    #     raise ValueError('MPI mode - the number of nodes is less then 2!')

    if multi_mpi.mpi_rank == multi_mpi.mpi_master:
        # MPI master node - solve problem at macro scale
        parser = ArgumentParser(description=__doc__,
                                formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument('--debug',
                            action='store_true', dest='debug',
                            default=False, help=helps['debug'])
        parser.add_argument('--debug_mpi',
                            action='store_true', dest='debug_mpi',
                            default=False, help=helps['debug_mpi'])
        parser.add_argument('-c', '--conf', metavar='"key : value, ..."',
                            action='store', dest='conf', type=str,
                            default=None, help=helps['conf'])
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
        parser.add_argument('--log', metavar='file',
                            action='store', dest='log',
                            default=None, help=helps['log'])
        parser.add_argument('-q', '--quiet',
                            action='store_true', dest='quiet',
                            default=False, help=helps['quiet'])
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('filename_in', nargs='?')
        options = parser.parse_args()

        for k in ['save_ebc', 'save_ebc_nodes', 'save_regions',
                  'save_regions_as_groups', 'save_field_meshes', 'solve_not']:
            setattr(options, k, False)

        if options.debug:
            from sfepy.base.base import debug_on_error; debug_on_error()

        if options.debug_mpi:
            multi_mpi.set_logging_level('debug')

        filename_in = options.filename_in
        output.set_output(filename=options.log,
                          quiet=options.quiet,
                          combined=options.log is not None)

        required, other = get_standard_keywords()
        conf = ProblemConf.from_file_and_options(filename_in, options,
            required, other, define_args=options.define_args)

        opts = conf.options
        nslaves = multi_mpi.cpu_count() - 1
        opts.n_mpi_homog_slaves = nslaves
        output_prefix = opts.get('output_prefix', 'sfepy:')

        app = PDESolverApp(conf, options, output_prefix)
        if hasattr(opts, 'parametric_hook'):  # Parametric study.
            parametric_hook = conf.get_function(opts.parametric_hook)
            app.parametrize(parametric_hook)
        app()

        multi_mpi.master_send_task('finalize', None)
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


if __name__ == '__main__':
    main()
