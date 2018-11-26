# -*- coding: utf-8 -*-
"""
Module for animating solutions in 1D.
Can also save them but requieres ffmpeg package
see save_animation method.
"""

import matplotlib.animation as animation
from matplotlib import pyplot as plt
import numpy as nm
from numpy import newaxis as nax
from matplotlib import pylab as plt
from matplotlib import colors
from os.path import join as pjoin
from toolz import accumulate


__author__ = 'tomas_zitka'

ffmpeg_path = 'C:\\Users\\tomas\\bin\\ffmpeg\\bin\\ffmpeg.exe'  # for saving animations


def animate1d(Y, X, T, ax=None, fig=None, ylims=None, labs=None, plott=None):
    """
    Animates solution of 1D problem into current figure.
    Keep reference to returned animation object otherwise
    it is discarded

    :param Y: solution, array |T| x |X| x n, where n is dimension of the solution
    :param X: space interval discetization
    :param T: time interval discretization
    :param ax: specify axes to plot to
    :param fig: specifiy figure to plot to
    :param ylims: limits for y axis, default are 10% offsets of Y extremes
    :param labs: labels to use for parts of the solution
    :param plott: plot type - how to plot data: tested plot, step

    :return: the animation object, keep it to see the animation, used for savig too
    """

    ax, fig, time_text = setup_axis(X, Y, ax, fig, ylims)

    if not isinstance(Y, nm.ndarray):
        Y = nm.stack(Y, axis=2)

    lines = setup_lines(ax, Y.shape, labs, plott)

    def animate(i):
        ax.legend()
        time_text.set_text("t= {0:3.2f} / {1:3.3}".format(T[i], T[-1]))
        # from sfepy.base.base import debug;
        # debug()
        if len(Y.shape) > 2:
            for ln, l in enumerate(lines):
                l.set_data(X, Y[i].swapaxes(0, 1)[ln])
            return tuple(lines) + (time_text,)
        # https://stackoverflow.com/questions/20624408/matplotlib-animating-multiple-lines-and-text
        else:
            lines.set_data(X, Y[i])
            return lines, time_text

    delay = int(nm.round(2000 * (T[-1] - T[0]) / len(T)))
    delay = 1000
    anim = animation.FuncAnimation(fig, animate, frames=len(T), interval=delay,
                                   blit=True, repeat=True, repeat_delay=250)

    return anim


def setup_axis(X, Y, ax=None, fig=None, ylims=None):
    """
    Setup axis, including timer for animation or snaps
    :param X: space disctretization to get limits
    :param Y: solution to get limits
    :param ax: ax where to put everything, if None current axes are used
    :param fig: fig where to put everything, if None current figure is used
    :param ylims: custom ylims, if None y axis limits are calculated from Y
    :return: ax, fig, time_text object to fill in text
    """
    if ax is None:
        fig = plt.gcf()
        ax = plt.gca()
    if ylims is None:
        lowery = nm.min(Y) - nm.min(Y) / 10
        uppery = nm.max(Y) + nm.max(Y) / 10
    else:
        lowery = ylims[0]
        uppery = ylims[1]
    ax.set_ylim(lowery, uppery)
    ax.set_xlim(X[0], X[-1])
    time_text = ax.text(X[0] + nm.sign(X[0]) * X[0] / 10, uppery - uppery / 10, 'empty', fontsize=15)
    return ax, fig, time_text


def setup_lines(ax, Yshape, labs, plott):
    """
    Sets up artist for animation or solution snaps
    :param ax: axes to use for artist
    :param Yshape: shape of the solution array
    :param labs: labels for the solution
    :param plott: type of plot to use i.e. steps or plot
    :return:
    """
    if plott is None:
        plott = ax.plot
    else:
        plott = ax.__getattribute__(plott)

    if len(Yshape) > 2:
        lines = [plott([], [], lw=2)[0] for foo in range(Yshape[2])]
        for i, l in enumerate(lines):
            if labs is None:
                l.set_label("q" + str(i + 1) + "(x, t)")
            else:
                l.set_label(labs[i])
    else:
        lines, = plott([], [], lw=2)
        if labs is None:
            lines.set_label("q(x, t)")
        else:
            lines.set_label(labs)
    return lines


def save_animation(anim, filename):
    """
    Saves animation as .mp4, requires ffmeg package
    :param anim: animation object
    :param filename: name of the file, without the .mp4 ending
    :return: None
    """
    plt.rcParams['animation.ffmpeg_path'] = ffmpeg_path
    writer = animation.FFMpegWriter(fps=24)
    anim.save(filename + ".mp4", writer=writer)


def sol_frame(Y, X, T, t0=.5, ax=None, fig=None, ylims=None, labs=None, plott=None):
    """
    Creates snap of solution at specified time frame t0, basically gets one frame from animate1d,
    but colors wont be the same :-(
    :param Y: solution, array |T| x |X| x n, where n is dimension of the solution
    :param X: space interval discetization
    :param T: time interval discretization
    :param t0: time to take snap at
    :param ax: specify axes to plot to
    :param fig: specifiy figure to plot to
    :param ylims: limits for y axis, default are 10% offsets of Y extremes
    :param labs: labels to use for parts of the solution
    :param plott: plot type - how to plot data: tested plot, step
    """

    ax, fig, time_text = setup_axis(X, Y, ax, fig, ylims)

    if not isinstance(Y, nm.ndarray):
        Y = nm.stack(Y, axis=2)

    lines = setup_lines(ax, Y.shape, labs, plott)

    nt0 = nm.abs(T - t0).argmin()

    ax.legend()
    time_text.set_text("t= {0:3.2f} / {1:3.3}".format(T[nt0], T[-1]))
    if len(Y.shape) > 2:
        for ln, l in enumerate(lines):
            l.set_data(X, Y[nt0].swapaxes(0, 1)[ln])
    else:
        lines.set_data(X, Y[nt0])
    return fig


def save_sol_snap(Y, X, T, t0=.5, filename=None, name=None, ylims=None, labs=None, plott=None):
    """
    Wrapper for sol_frame, saves the frame to file specified.
    :param name: name of the solution e.g. name of the solver used
    :param filename: name of the file, overrides automatic generation
    :param Y: solution, array |T| x |X| x n, where n is dimension of the solution
    :param X: space interval discetization
    :param T: time interval discretization
    :param t0: time to take snap at
    :param ylims: limits for y axis, default are 10% offsets of Y extremes
    :param labs: labels to use for parts of the solution
    :param plott: plot type - how to plot data: tested plot, step
    """

    if filename is None:
        filename = "{0}_solsnap{1:3.2f}-{2:3.3}".format(name, t0, T[-1]).replace(".", "_")
        if name is None:
            name = "unknown_solver"
        filename = "{0}_solsnap{1:3.2f}-{2:3.3}".format(name, t0, T[-1]).replace(".", "_")
        filename = pjoin("semestralka", "figs", filename)

    fig = plt.figure(filename)

    snap1 = sol_frame(Y, X, T, t0=t0, ylims=ylims, labs=labs, plott=None)
    if not isinstance(Y, nm.ndarray):
        plt.plot(X, Y[0][0], label="q(x, 0)")
    else:
        if len(Y.shape) > 2:
            plt.plot(X, Y[0, :, 0], label="q(x, 0)")
        else:
            plt.plot(X, Y[0, :], label="q(x, 0)")
    plt.legend()
    snap1.savefig(filename)
    return fig


def plotsXT(Y1, Y2, YE, extent, lab1=None, lab2=None, lab3=None):
    """
    Plots Y1 and Y2 to one axes and YE to the second axes,
    Y1 and Y2 are presumed to two solution and YE their error - hence the names
    :param Y1:
    :param Y2:
    :param YE:
    :param extent:
    :return:
    """

    # >> Plot contours
    cmap1 = plt.cm.get_cmap("bwr")
    cmap1.set_bad('white')
    # cmap2 = plt.cm.get_cmap("BrBG")
    # cmap2.set_bad('white')
    bounds = nm.arange(-1, 1, .05)
    norm1 = colors.BoundaryNorm(bounds, cmap1.N)
    # norm2 = colors.BoundaryNorm(bounds, cmap2.N)

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharey=True)
    fig.suptitle("X-T plane plot")
    if lab1 is not None:
        ax1.set(title=lab1)
    # FIXME X-T plane plot of exact solution is not ok
    c1 = ax1.imshow(Y1, extent=extent,
                    cmap=cmap1, norm=norm1,
                    interpolation='none',
                    origin='lower')
    ax1.grid()
    if lab2 is not None:
        ax2.set(title=lab2)
    c2 = ax2.imshow(Y2, extent=extent,
                    cmap=cmap1, norm=norm1,
                    interpolation='none',
                    origin='lower')
    ax2.grid()

    if lab3 is not None:
        ax3.set(title=lab3)
    c3 = ax3.imshow(YE, extent=extent,
                    cmap="bwr", norm=norm1,
                    interpolation='none',
                    origin='lower')
    ax3.grid()
    fig.colorbar(c3, ax=[ax1, ax2, ax3])


def load_vtks(fold, name, tn, order, tns=None):
    """
    Reads series of .vtk files and crunches them into form
    suitable for plot10_DG_sol.

    In fact reads cell data of shape mesh.coors - 1 i.e. mesh.n_el

    :param fold: folder where to look for files
    :param name: used in {name}.i.vtk, i = 0,1, ... tn - 1
    :param tn: number of time steps, i.e. number of files
    :param order: order of approximation used
    :return: space coors, solution data
    """

    from sfepy.discrete.fem.meshio import VTKMeshIO

    if tns is None:
        tns = tn

    digs = int(nm.ceil(nm.log10(tn)))   # number of digits in filename
    full_name = ".".join((name, ("{:0" + str(digs) + "d}"), "vtk"))
    io = VTKMeshIO(pjoin(fold, full_name.format(0)))
    coors = io.read_coors()[:, 0, None]
    u = nm.zeros((order + 1, coors.shape[0] - 1, tns, 1))

    for i in range(tns):
        io = VTKMeshIO(pjoin(fold, full_name.format(i)))
        data = io.read_data(step=0)  # parameter "step" does nothing for VTKMeshIO, but is obligatory
        for ii in range(order + 1):
            u[ii, :, i, 0] = data['u{}'.format(ii)].data

    return coors, u

def plot1D_DG_sol(coors, t0, t1, u, ic=lambda x: 0.0, tn=None, dt=None):
    """
    Plots solution produced by DG to 1D problem, handles discontinuities,
    u are vectors of coefficients, for each order one

    :param coors: coordinates of the mesh
    :param t0: starting time
    :param t1: final time
    :param tn: number of time steps, must correspond to dimention of u
    :param u: shape(u) = (order, space_steps, t_steps, 1)
    :param ic: analytical initial condition, optional
    :return: nothing
    """
    XN1 = coors[-1]
    X1 = coors[0]
    n_nod = len(coors)

    figs, axs = plt.subplots()
    X = (coors[1:] + coors[:-1]) / 2
    if tn is not None and dt is not None:
        T = nm.array(nm.cumsum(nm.ones(tn) * dt))
    elif tn is not None:
        T = nm.linspace(t0, t1, tn)
    elif dt is not None:
        tn = float(t1 - t0) / dt
        T = nm.linspace(t0, t1, tn)
    # sic = TSSolver.initial_cond

    # Plot mesh
    plt.vlines(coors[:, 0], ymin=0, ymax=.5, colors="grey")
    plt.vlines((coors[0], coors[-1]), ymin=0, ymax=.5, colors="k")
    plt.vlines(X, ymin=0, ymax=.3, colors="grey", linestyles="--")

    # Plot IC and its sampling
    # TODO get IC sampling, from where?
    c0 = plt.plot(X, u[0, :, 0, 0], label="IC-0", marker=".", ls="")[0].get_color()
    c1 = plt.plot(X, u[1, :, 0, 0], label="IC-1", marker=".", ls="")[0].get_color()
    # # plt.plot(coors, .1*alones(n_nod), marker=".", ls="")
    plt.step(coors[1:], u[0, :, 0,  0], color=c0)
    plt.step(coors[1:], u[1, :, 0,  0], color=c1)
    # plt.plot(coors[1:], sic[1, :], label="IC-1", color=c1)
    xs = nm.linspace(X1, XN1, 500)[:, None]
    plt.plot(xs, ic(xs), label="IC-ex")

    # Animate sampled solution
    anim = animate1d(u[:, :, :, 0].T, coors[1:], T, axs, figs, ylims=[-1, 2], plott="step")
    plt.xlim(coors[0] - .1, coors[-1] + .1)
    plt.legend(loc="upper left")
    plt.title("Sampled solution")

    figr, axr = plt.subplots()

    # Plot mesh
    plt.vlines(coors[:, 0], ymin=0, ymax=.5, colors="grey")
    plt.vlines((coors[0], coors[-1]), ymin=0, ymax=.5, colors="k")
    plt.vlines(X, ymin=0, ymax=.3, colors="grey", linestyles="--")

    # Plot discontinuously!
    # (order, space_steps, t_steps, 1)
    ww = nm.zeros((3 * n_nod - 1, tn, 1))
    ww[0, :] = u[0, 0, :] - u[1, 0, :]  # left bc
    ww[-1, :] = u[0, -1, :] + u[1, -1, :]  # right bc

    ww[0:-2:3] = u[0, :, :] - u[1, :, :]  # left edges of elements
    ww[1:-1:3] = u[0, :, :] + u[1, :, :]  # right edges of elements
    ww[2::3, :] = nm.NaN  # NaNs ensure plotting of discontinuities at element borders

    # nodes for plotting reconstructed solution
    xx = nm.zeros((3 * n_nod - 1, 1))
    xx[0] = coors[0]
    xx[-1] = coors[-1]
    # the ending ones are still a bit odd, but hey, it works!
    xx[1:-1] = nm.repeat(coors[1:], 3)[:, None]
    # plt.vlines(xx, ymin=0, ymax=.3, colors="green")

    # plot reconstructed IC
    plt.plot(xx, ww[:, 0], label="IC")

    # Animate reconstructed
    anim_disc = animate1d(ww[:, :, 0].T, xx, T, axr, figr, ylims=[-1, 2])
    plt.xlim(coors[0] - .1, coors[-1] + .1)
    plt.legend(loc="upper left")
    plt.title("Reconstructed solution")

    # sol_frame(u[:, :, :, 0].T, nm.append(coors, coors[-1]), T, t0=0., ylims=[-1, 1], plott="step")
    plt.show()