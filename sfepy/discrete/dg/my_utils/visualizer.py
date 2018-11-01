# -*- coding: utf-8 -*-
"""
Module for animating solutions in 1D.
Can also save them but requieres ffmpeg package
see save_animation method.
"""

import matplotlib.animation as animation
from matplotlib import pyplot as plt
import numpy as np
from numpy import newaxis as nax
from matplotlib import pylab as plt
from matplotlib import colors
from os.path import join as pjoin


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

    if not isinstance(Y, np.ndarray):
        Y = np.stack(Y, axis=2)

    lines = setup_lines(ax, Y.shape, labs, plott)

    def animate(i):
        ax.legend()
        time_text.set_text("t= {0:3.2f} / {1:3.3}".format(T[i], T[-1]))
        if len(Y.shape) > 2:
            for ln, l in enumerate(lines):
                l.set_data(X, Y[i].swapaxes(0, 1)[ln])
            return tuple(lines) + (time_text,)
        # https://stackoverflow.com/questions/20624408/matplotlib-animating-multiple-lines-and-text
        else:
            lines.set_data(X, Y[i])
            return lines, time_text

    delay = int(np.round(2000 * (T[-1] - T[0])/len(T)))
    # delay = 1000
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
        lowery = np.min(Y) - np.min(Y) / 10
        uppery = np.max(Y) + np.max(Y) / 10
    else:
        lowery = ylims[0]
        uppery = ylims[1]
    ax.set_ylim(lowery, uppery)
    ax.set_xlim(X[0], X[-1])
    time_text = ax.text(X[0] + np.sign(X[0]) * X[0] / 10, uppery - uppery / 10, 'empty', fontsize=15)
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

    if not isinstance(Y, np.ndarray):
        Y = np.stack(Y, axis=2)

    lines = setup_lines(ax, Y.shape, labs, plott)

    nt0 = np.abs(T - t0).argmin()

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
    if not isinstance(Y, np.ndarray):
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
    bounds = np.arange(-1, 1, .05)
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


