import numpy as nm
import matplotlib.pyplot as plt

from os.path import join as pjoin
import os
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



def load_and_plot_fun(folder, filename, t0, t1, approx_order, ic_fun=None):

    # load time data
    lmesh, u = load_1D_vtks(folder, filename, order=approx_order)
    plot1D_DG_sol(lmesh, t0, t1, u, tn=30, ic=ic_fun,
                  delay=100, polar=False)


    coors, u_end = load_state_1D_vtk(pjoin(folder, filename + "_end.vtk"), order=approx_order)
    coors, u_start = load_state_1D_vtk(pjoin(folder, filename + "_start.vtk"), order=approx_order)


    plot_1D_legendre_dofs(coors, [u_start.swapaxes(0, 1)[:, :, 0], u_end.swapaxes(0, 1)[:, :, 0]])

    plt.figure("reconstructed")
    ww_s, xx = reconstruct_legendre_dofs(coors, None, u_end)
    ww_e, _ = reconstruct_legendre_dofs(coors, None, u_start)

    plt.plot(xx, ww_s[:, 0])
    plt.plot(xx, ww_e[:, 0])
    plt.show()


if __name__ == '__main__':
    # TODO get neccesary data for plotting, maybe from args?
    # load_1D_vtks(folder, filename, t0, t1, approx_order, ic_fun)
    pass