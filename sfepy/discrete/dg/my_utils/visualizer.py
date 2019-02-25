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


def animate1d(Y, X, T, ax=None, fig=None, ylims=None, labs=None, plott=None, delay=None):
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

    if delay is None:
        delay = int(nm.round(2000 * (T[-1] - T[0]) / len(T)))
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


def load_state_1D_vtk(name, order):
    """
    Load one VTK file containing state in time
    :param name:
    :param order:
    :return:
    """

    from sfepy.discrete.fem.meshio import VTKMeshIO
    from glob import glob
    from os.path import join as pjoin
    io = VTKMeshIO(name)
    data = io.read_data(step=0)
    coors = io.read_coors()[:, 0, None]
    u = nm.zeros((order + 1, coors.shape[0] - 1, 1, 1))
    for ii in range(order + 1):
        u[ii, :, 0, 0] = data['u_modal{}'.format(ii)].data

    return coors, u


def load_1D_vtks(fold, name, order, tns=None):
    """
    Reads series of .vtk files and crunches them into form
    suitable for plot10_DG_sol.

    Attempts to read cell data for variables u0, u1 ...

    Resulting solution data have shape:
    (order, nspace_steps, ntime_steps, 1)

    :param fold: folder where to look for files
    :param name: used in {name}.i.vtk, i = 0,1, ... tns - 1
    :param tns: number of time steps, i.e. number of files to read
    :param order: order of approximation used in u1, u1 ...u{order}
    :return: space coors, solution data
    """

    from sfepy.discrete.fem.meshio import VTKMeshIO
    from glob import glob
    from os.path import join as pjoin
    files = glob(pjoin(fold, name) + ".[0-9]*")

    if len(files) == 0:
        print("No files for name {} found in {}".format(name, fold))
        return
    io = VTKMeshIO(files[0])
    coors = io.read_coors()[:, 0, None]

    tn = len(files)
    nts = sorted([int(f.split(".")[-2]) for f in files])

    digs = len(files[0].split(".")[-2])
    full_name_form = ".".join((pjoin(fold, name), ("{:0" + str(digs) + "d}"), "vtk"))

    u = nm.zeros((order + 1, coors.shape[0] - 1, tn, 1))
    for i, nt in enumerate(nts):
        io = VTKMeshIO(full_name_form.format(nt))
        data = io.read_data(step=0)  # parameter "step" does nothing for VTKMeshIO, but is obligatory
        for ii in range(order + 1):
            u[ii, :, i, 0] = data['u_modal{}'.format(ii)].data

    return coors, u


def plot1D_DG_sol(coors, t0, t1, u,
                  tn=None, dt=None, ic=lambda x: 0.0, delay=None, polar=False):
    """
    Animates solution to 1D problem produced by DG:
        1. animates DOF values in elements as steps
        2. animates reconstructed solution with discontinuities

    :param coors: coordinates of the mesh
    :param t0: starting time
    :param t1: final time
    :param u: vectors of DOFs, for each order one, shape(u) = (order, nspace_steps, ntime_steps, 1)
    :param ic: analytical initial condition, optional
    :param tn: number of time steps to plot, starting at 0, if None and dt is not None run animation through
        all time steps, spaced dt within [t0, tn]
    :param dt: time step size, if None and tn is not None computed as (t1- t0) / tn otherwise set to 1
        if dt and tn are both None, t0 and t1 are ignored and solution is animated as if in time 0 ... ntime_steps
    :return: anim object of DOFs, anim object of reconstruction
    """


    # Setup space coordinates
    XN = coors[-1]
    X1 = coors[0]
    Xvol = XN - X1
    X = (coors[1:] + coors[:-1]) / 2
    XS = nm.linspace(X1, XN, 500)[:, None]
    ics = ic(XS)

    if polar:
        coors *= 2*nm.pi
        X *= 2*nm.pi
        XS *= 2*nm.pi

    # Setup times
    if tn is not None and dt is not None:
        T = nm.array(nm.cumsum(nm.ones(tn) * dt))
    elif tn is not None:
        T, dt = nm.linspace(t0, t1, tn, retstep=True)
    elif dt is not None:
        tn = int(nm.ceil(float(t1 - t0) / dt))
        T = nm.linspace(t0, t1, tn)
    else:
        T = nm.arange(nm.shape(u)[2])

    n_nod = len(coors)
    n_el_nod = nm.shape(u)[0]
    # prepend u[:, 0, ...] to all time frames for plotting steps
    u_step = nm.append(u[:, 0:1, :, 0], u[:, :, :, 0], axis=1)

    # Plot DOFs directly
    figs = plt.figure()
    if polar:
        axs = plt.subplot(111, projection='polar')
        axs.set_theta_direction('clockwise')

    else:
        axs = plt.subplot(111)

    # Plot mesh
    axs.vlines(coors[:, 0], ymin=0, ymax=.5, colors="grey")
    axs.vlines((X1, XN), ymin=0, ymax=.5, colors="k")
    axs.vlines(X, ymin=0, ymax=.3, colors="grey", linestyles="--")

    axs.plot([X1, XN], [1, 1], 'k')

    # Plot IC and its sampling
    for i in range(n_el_nod):
        c0 = axs.plot(X, u[i, :, 0, 0], label="IC-{}".format(i), marker=".", ls="")[0].get_color()
        # c1 = plt.plot(X, u[1, :, 0, 0], label="IC-1", marker=".", ls="")[0].get_color()
        # # plt.plot(coors, .1*alones(n_nod), marker=".", ls="")
        axs.step(coors[1:], u[i, :, 0,  0], color=c0)
        # plt.step(coors[1:], u[1, :, 0,  0], color=c1)
        # plt.plot(coors[1:], sic[1, :], label="IC-1", color=c1)
    axs.plot(XS, ics, label="IC-ex")

    # Animate sampled solution DOFs directly
    anim_dofs = animate1d(u_step.T, coors, T, axs, figs, ylims=[-1, 2], plott="step", delay=delay)
    if not polar:
        axs.set_xlim(coors[0] - .1 * Xvol, coors[-1] + .1 * Xvol)
    axs.legend(loc="upper left")
    axs.set_title("Sampled solution")


    # Plot reconstructed solution
    figr = plt.figure()
    if polar:
        axr = plt.subplot(111, projection='polar')
        axr.set_theta_direction('clockwise')
    else:
        axr = plt.subplot(111)

    # Plot mesh
    axr.vlines(coors[:, 0], ymin=0, ymax=.5, colors="grey")
    axr.vlines((X1, XN), ymin=0, ymax=.5, colors="k")
    axr.vlines(X, ymin=0, ymax=.3, colors="grey", linestyles="--")

    axr.plot([X1, XN], [1, 1], 'k')
    # Plot discontinuously!
    # (order, space_steps, t_steps, 1)
    ww, xx = reconstruct_legendre_dofs(coors, tn, u)
    # plt.vlines(xx, ymin=0, ymax=.3, colors="green")

    # plot reconstructed IC
    axr.plot(xx, ww[:, 0], label="IC")

    # Animate reconstructed solution
    anim_recon = animate1d(ww[:, :, 0].T, xx, T, axr, figr, ylims=[-1, 2], delay=delay)
    if not polar:
        axr.set_xlim(coors[0] - .1 * Xvol, coors[-1] + .1 * Xvol)
    axr.legend(loc="upper left")
    axr.set_title("Reconstructed solution")

    # sol_frame(u[:, :, :, 0].T, nm.append(coors, coors[-1]), T, t0=0., ylims=[-1, 1], plott="step")
    plt.show()
    return anim_dofs, anim_recon


def plot_1D_legendre_dofs(coors, dofss, fun=None):
    """
    Plots values of DOFs as steps
    :param coors: coordinates of nodes of the mesh
    :param dofss: iterable of different projections' DOFs into legendre space
    :param fun: analytical function to plot
    :return:
    """
    X = (coors[1:] + coors[:-1]) / 2
    plt.figure("DOFs for function fun")
    plt.gcf().clf()
    for ii, dofs in enumerate(dofss):
        for i in range(dofs.shape[1]):
            c0 = plt.plot(X, dofs[:, i], label="fun-{}dof-{}".format(ii, i), marker=".", ls="")[0].get_color()
            # # plt.plot(coors, .1*alones(n_nod), marker=".", ls="")
            plt.step(coors[1:], dofs[:, i], color=c0)
            # plt.plot(coors[1:], sic[1, :], label="IC-1", color=c1)

    if fun is not None:
        xs = nm.linspace(nm.min(coors), nm.max(coors), 500)[:, None]
        plt.plot(xs, fun(xs), label="fun-ex")
    plt.legend()
    # plt.show()


def reconstruct_legendre_dofs(coors, tn, u):
    """
    Creates solution and coordinates vector which when plotted by as

        plot(xx, ww)

    represent solution reconstructed from DOFs in Legendre poly space at
    cell borders

    So far work only for order 1
    # TODO reconstruct solution on finer mesh to display curvature in higher order -> different function
    :param coors: coors of nodes of the mesh
    :param u: vectors of DOFs, for each order one, shape(u) = (order, nspace_steps, ntime_steps, 1)
    :param tn: number of time steps to reconstruct, if None all steps are reconstructed
    :return: ww - solution values vector, shape is (3 * nspace_steps - 1, ntime_steps, 1),
             xx - corresponding coordinates vector, shape is (3 * nspace_steps - 1, 1)
    """


    XN = coors[-1]
    X1 = coors[0]
    n_nod = len(coors) - 1
    if tn is None:
        tn = nm.shape(u)[2]
    n_el_nod = nm.shape(u)[0]

    ww = nm.zeros((3 * n_nod - 1, tn, 1))

    for i in range(n_el_nod):
        ww[0:-1:3] = ww[0:-1:3] + (-1)**i * u[i, :, :]  # left edges of elements
        ww[1::3] = ww[1::3] + u[i, :, :]  # right edges of elements
    ww[2::3, :] = nm.NaN  # NaNs ensure plotting of discontinuities at element borders

    # nodes for plotting reconstructed solution
    xx = nm.zeros((3 * n_nod - 1, 1))
    xx[0] = X1
    xx[-1] = XN
    # the ending is still a bit odd, but hey, it works!
    xx[1:-1] = nm.repeat(coors[1:-1], 3)[:, None]
    return ww, xx
