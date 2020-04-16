#!/usr/bin/env python
# -*- coding:utf-8 -*-
""""
Script for running DG conf files.
"""
import sys
import os
import time
from os.path import join as pjoin

import argparse
import numpy as nm

from sfepy.applications.pde_solver_app import PDESolverApp
from sfepy.base.conf import ProblemConf
from sfepy.base.ioutils import ensure_path
from sfepy.base.base import (get_default, output, assert_,
                             Struct, basestr, IndexedStruct)

from sfepy.discrete.dg.my_utils.plot_1D_dg import load_and_plot_fun
from run_dg_utils import clear_folder

from examples.dg.run_dg_utils import  calculate_num_order, outputs_folder, \
    plot_conv_results, build_attrs_string, output, compute_erros, configure_output

from examples.dg.example_dg_common import get_1Dmesh_hook, diffusion_schemes_explicit

def create_argument_parser():
    parser = argparse.ArgumentParser(
        description='Run SfePy DG example conf python files',
        epilog='(c) 2019  Man-machine Interaction at NTC UWB, ' +
               '\nauthor: Tomas Zitka, email: zitkat@ntc.zcu.cz')

    parser.add_argument("problem_file",
                        help="File containing problem configuration, "+
                             "can include define method and dim variable",
                        metavar='path')

    parser.add_argument("--output",
                        help="Output directory, can contain {}, " +
                             "exampe name is than plugged in.",
                        dest="output_dir",
                        metavar='path', type=str, default=None)

    parser.add_argument("-m", "--mesh", help="File with starting mesh to use",
                        default=None, metavar='path', action='store',
                        dest='mesh_file',)

    parser.add_argument('-dp', '--doplot', help="To plot 1D case",
                        action="store_true")

    parser.add_argument("-v", "--verbose", help="To be verbose or",
                        default=False, action='store_true', dest='verbose',)

    parser.add_argument("--noscreenlog", help="Do not print log to screen",
                        default=False, action='store_true',
                        dest='no_output_screen', )

    parser.add_argument('--logfile', type=str,
                        action='store', dest='output_log_name',
                        default=None, help="Path to log output file.")

    parser.add_argument('--advelo', metavar='float', type=float,
                        action='store', dest='advelo',
                        default=1, help="Advection velocity")

    parser.add_argument('--adflux', metavar='float', type=float,
                        action='store', dest='adflux',
                        default=0, help="Advection flux parameter, " +
                                        "\n0 - upwind, " +
                                        "\n1 - central")

    parser.add_argument("--limit", help="Use 1D or 2D moment limiter",
                        default=False, action='store_true', dest='limit',)

    parser.add_argument('--cw', metavar='float', type=float,
                        action='store', dest='cw',
                        default=1, help="Diffusion penalty coefficient")

    parser.add_argument('--diffcoef', metavar='float', type=float,
                        action='store', dest='diffcoef',
                        default=1, help="Diffusion coeffcient")

    parser.add_argument('--diffscheme', type=str,
                        choices=diffusion_schemes_explicit.keys(),
                        action='store', dest='diffscheme',
                        default="symmetric", help="Scheme to use for diffusion")

    parser.add_argument('--cfl', metavar='float', type=float,
                        action='store', dest='cfl',
                        default=0.314, help="CFL coefficient")

    parser.add_argument('--dt', metavar='float', type=float,
                        action='store', dest='dt',
                        default=None, help="Time step size, overrides CFL "+
                                           "coefficient")

    parser.add_argument('--order', metavar="int" , default=1,
                        help='Approximation order', type=int)

    return parser


def get_parametrized_conf(filename, args):
    import importlib
    problem_module_name = "examples.dg." + filename.replace(".py", "").strip("\\.") \
        .replace("\\", ".").replace("/", ".")
    print(problem_module_name)
    problem_module = importlib.import_module(problem_module_name)

    if hasattr(problem_module, "define"):
        mod = sys.modules[__name__]
        # noinspection PyCallingNonCallable
        problem_conf = ProblemConf.from_dict(
            problem_module.define(
                filename_mesh=args.mesh_file,
                approx_order=args.order,

                flux=args.adflux,
                limit=args.limit,

                Cw=args.cw,
                diffusion_coef=args.diffcoef,
                diff_scheme_name=args.diffscheme,

                CFL=args.cfl,
                dt=args.dt,
            ), mod, verbose=args.verbose)
    else:
        output("Problem file is not parametrized, arguments ignored")
        problem_conf = ProblemConf.from_file(filename, verbose=False)
        problem_conf.verbose = args.verbose

    return problem_conf


def main(argv):
    if argv is None:
        argv = sys.argv[1:]

    parser = create_argument_parser()
    args = parser.parse_args(argv)

    # configure_output({'output_screen': not args.no_output_screen,
    #                   'output_log_name': "last_run.txt"})

    conf_file_name = args.problem_file

    output("Processing conf file {}".format(conf_file_name))
    pc = get_parametrized_conf(conf_file_name, args)

    output("----------------Running--------------------------")
    output("{}: {}".format(pc.example_name, time.asctime()))

    if args.output_dir is None:
        output_folder = pjoin(outputs_folder, "output")
    elif "{}" in args.output_dir:
        output_folder = args.output_dir.format(pc.example_name)
    else:
        output_folder = args.output_dir

    output_name_trunk_folder = pjoin(output_folder, pc.example_name,
                                     str(pc.approx_order) + "/")

    output_name_trunk_name = pc.example_name + str(pc.approx_order)
    output_name_trunk = pjoin(output_name_trunk_folder, output_name_trunk_name)
    ensure_path(output_name_trunk_folder)
    output_format = "{}.*.{}".format(output_name_trunk,
                                      pc.options.output_format
                                      if hasattr(pc.options, "output_format")
                                      else "vtk")

    output("Output set to {}, clearing.".format(output_format))
    clear_folder(output_format, confirm=False)

    sa = PDESolverApp(pc, Struct(output_filename_trunk=output_name_trunk,
                                 save_ebc=False,
                                 save_ebc_nodes=False,
                                 save_region=False,
                                 save_regions=False,
                                 save_regions_as_groups=False,
                                 save_field_meshes=False,
                                 solve_not=False), "sfepy")
    tt = time.clock()
    sa()
    elapsed = time.clock() - tt
    output("{}: {}".format(pc.example_name, time.asctime()))
    output("------------------Finished------------------\n\n")


    if pc.dim == 1 and args.doplot:
        if pc.transient:
            load_times = min(pc.options.save_times, sa.problem.ts.n_step)
            load_and_plot_fun(output_name_trunk_folder, output_name_trunk_name,
                              pc.t0, pc.t1, load_times, pc.approx_order,
                              pc.get_ic, exact=getattr(pc, "analytic_sol", None),
                              polar=False, compare=True)
        else:
            load_times = 1
            load_and_plot_fun(output_name_trunk_folder, output_name_trunk_name,
                              pc.t0, pc.t1, load_times, pc.approx_order)


if __name__ == '__main__':
    main(sys.argv[1:])
