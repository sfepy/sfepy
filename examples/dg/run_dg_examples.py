#!python
# -*- coding:utf-8 -*-
""""
Script for running DG conf files.
"""
import sys
import os
from os.path import join as pjoin

import glob
import time
import subprocess as sub
import logging
import argparse
import numpy as nm

from sfepy.applications.pde_solver_app import PDESolverApp
from sfepy.base.conf import ProblemConf
from sfepy.base.ioutils import ensure_path
from sfepy.base.base import (get_default, output, assert_,
                             Struct, basestr, IndexedStruct)

from sfepy.discrete.dg.my_utils.plot_1D_dg import load_and_plot_fun
from sfepy.discrete.dg.my_utils.plot_1D_dg import clear_folder

parser = argparse.ArgumentParser(
    description='Run SfePy DG example conf python files',
    epilog='(c) 2019 by T. Zitka , Man-machine Interaction at NTC UWB')
parser.add_argument("conf_file", help="File containing problem configuration")
parser.add_argument('-p', '--plot', help="To plot 1D case", action="store_true")


def main(argv):
    if argv is None:
        argv = sys.argv[1:]

    args = parser.parse_args(argv)
    conf_file_name = args.conf_file

    output("Processing conf file {}".format(conf_file_name))
    pc = ProblemConf.from_file(conf_file_name)
    pc.verbose = False

    output("Running {}".format(pc.example_name))

    output_folder = "output"
    output_name_trunk_folder = pjoin(output_folder, pc.example_name,
                                     str(pc.approx_order) + "/")
    output_name_trunk_name = pc.example_name + str(pc.approx_order)
    output_name_trunk = pjoin(output_name_trunk_folder, output_name_trunk_name)
    ensure_path(output_name_trunk_folder)
    output_format = "{}.*.{}".format(output_name_trunk,
                                      pc.options.output_format
                                      if hasattr(pc.options, "output_format")
                                      else "vtk")

    output("Output set to {}, clearing ...".format(output_format))
    clear_folder(output_format, confirm=False)

    sa = PDESolverApp(pc, Struct(output_filename_trunk=output_name_trunk,
                                 save_ebc=False,
                                 save_ebc_nodes=False,
                                 save_region=False,
                                 save_regions=False,
                                 save_regions_as_groups=False,
                                 save_field_meshes=False,
                                 solve_not=False), "sfepy")

    sa()

    if pc.dim == 1 and args.plot:
        if pc.transient:
            load_times = min(pc.options.save_times, sa.problem.ts.n_step)
            load_and_plot_fun(output_name_trunk_folder, output_name_trunk_name,
                              pc.t0, pc.t1, load_times, pc.approx_order,
                              pc.get_ic, exact=getattr(pc, "analytic_sol", None))
        else:
            load_times = 1
            load_and_plot_fun(output_name_trunk_folder, output_name_trunk_name,
                              pc.t0, pc.t1, load_times, pc.approx_order)


if __name__ == '__main__':
    main(sys.argv[1:])
