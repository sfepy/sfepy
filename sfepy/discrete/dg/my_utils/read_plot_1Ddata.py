#!python
# -*- coding:utf-8 -*-
""""
Script for plotting 1D DG FEM data stored in VTK files
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

import numpy as nm
import matplotlib.pyplot as plt
import argparse

from os.path import join as pjoin
import os
import sys
from glob import glob

from sfepy.base.base import (get_default, output, assert_,
                             Struct, basestr, IndexedStruct)
from sfepy.discrete.dg.my_utils.visualizer import load_1D_vtks, plot1D_DG_sol
from sfepy.discrete.dg.dg_field import get_unraveler, get_raveler
from sfepy.discrete.dg.my_utils.visualizer import \
    load_state_1D_vtk, plot_1D_legendre_dofs, reconstruct_legendre_dofs



def clear_folder(clear_format, confirm=False):
    """
    Deletes files matching the format
    :param clear_format:
    :param confirm:
    :return:
    """
    files = glob(clear_format)
    doit = True
    if confirm:
        for file in files:
            output("Will delete file {}".format(file))
        doit = input("--------------\nDelete files [Y/n]? ").strip() == "Y"

    if doit:
        for file in files:
            os.remove(file)



def load_and_plot_fun(folder, filename, t0, t1, tn, approx_order, ic_fun=None, compare=False, polar=False):

    # load time data
    lmesh, u = load_1D_vtks(folder, filename, order=approx_order)
    plot1D_DG_sol(lmesh, t0, t1, u, tn=tn, ic=ic_fun,
                  delay=100, polar=polar)

    if compare:
        # TODO get name of first and last time step file
        coors, u_end = load_state_1D_vtk(pjoin(folder, filename + "_594.vtk"), order=approx_order)
        coors, u_start = load_state_1D_vtk(pjoin(folder, filename + "_000.vtk"), order=approx_order)

        plot_1D_legendre_dofs(coors, [u_start.swapaxes(0, 1)[:, :, 0], u_end.swapaxes(0, 1)[:, :, 0]])

        plt.figure("reconstructed")
        ww_s, xx = reconstruct_legendre_dofs(coors, None, u_end)
        ww_e, _ = reconstruct_legendre_dofs(coors, None, u_start)
        plt.plot(xx, ww_s[:, 0])
        plt.plot(xx, ww_e[:, 0])
        plt.show()


def main(argv):
    parser = argparse.ArgumentParser(description='Plotting of 1D DG data in VTK files',
                                     epilog='(c) 2019 by T. Zitka , Man-machine Interaction at NTC UWB')
    parser.add_argument("input_name", help="""Folder or name of the example in output folder with VTK data, file names 
        in format <name>.[0-9]*.vtk. If not provided ask for the name of the example""", nargs="?")
    parser.add_argument("-t0", "--start_time", type=float, default=0, help="Start time of the simulation")
    parser.add_argument("-t1", "--end_time", type=float, default=1., help="End time of the simulation")
    parser.add_argument("-o", "--order", type=int, default=None, help="""Order of the approximation, when example folder 
        cantains more orders this chooses the one""")
    parser.add_argument("-cf", "--compare-final", action="store_true", help="To compare starting and final time - " +
                                                                            "usefull for periodic boundary problems")
    parser.add_argument("-p", "--polar", help="Plot in polar projection", action="store_true")

    if argv is None:
        argv = sys.argv[1:]
    args = parser.parse_args(argv)

    t0 = args.start_time
    t1 = args.end_time
    cf = args.compare_final
    pol = args.polar
    if args.input_name is None:
        input_name = str(input("Please provide name of the example or path to data: "))
    else:
        input_name = args.input_name

    if os.path.isdir(input_name):
        full_infolder_path = os.path.abspath(input_name)
    else:
        input_name = pjoin("output", input_name)
        if os.path.isdir(input_name):
            full_infolder_path = os.path.abspath(input_name)
        else:
            print("Example {} not found in {}".format(input_name, os.path.abspath(output)))
            return

    print("Input folder found: {}".format(full_infolder_path))
    print("Input folder contains results for orders {}"
          .format([os.path.basename(fol) for fol in glob(pjoin(full_infolder_path, "[0-9]*"))]))

    if args.order is None:
        order = str(input("Please provide order of approximation, default is 1: "))
        if len(order) == 0:
            order = 1
        else:
            try:
                order = int(order)
            except ValueError:
                print("Value {} for order not understood!".format(order))
                return
    else:
        order = args.order
    print("Looking for results of order {}".format(order))
    full_infolder_path = pjoin(full_infolder_path, str(order))
    if not os.path.isdir(full_infolder_path):
        print("Input folder with order {} not found".format(full_infolder_path))
        return

    contents = glob(pjoin(full_infolder_path, "*.vtk"))
    tn = len(contents)  # we assume ll th contents are time step data files
    if tn == 0:
        print("Input folder {} is empty!".format(full_infolder_path))
        return
    base_name = os.path.basename(contents[0]).split(".")[0]
    print("Found {} filse, basename is {}".format(tn, base_name))
    print("Plotting ...")
    load_and_plot_fun(full_infolder_path, base_name, t0, t1, tn, order, compare=cf, polar=pol)










if __name__ == '__main__':

    main(sys.argv[1:])