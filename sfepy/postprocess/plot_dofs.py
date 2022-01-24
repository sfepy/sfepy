"""
Functions to visualize the mesh connectivity with global and local DOF
numberings.
"""
import numpy as nm

import matplotlib.pyplot as plt

def _get_axes(ax, dim):
    if ax is None:
        fig = plt.figure()
        if dim == 3:
            from mpl_toolkits.mplot3d import axes3d
            axes3d # Make pyflakes happy...

            ax = fig.add_subplot(111, projection='3d')

        else:
            ax = fig.add_subplot(111)

    return ax

def _to2d(coors):
    if coors.shape[1] == 1:
        coors = nm.c_[coors, nm.zeros_like(coors)]

    return coors

def plot_points(ax, coors, vals=None, point_size=20,
                show_colorbar=False):
    """
    Plot points with given coordinates, optionally colored using `vals` values.
    """
    dim = coors.shape[1]
    ax = _get_axes(ax, dim)

    colors = 'b' if vals is None else vals

    coors = _to2d(coors)
    sc = ax.scatter(*coors.T, s=point_size, c=colors, alpha=1)

    if show_colorbar and (vals is not None):
        plt.colorbar(sc)

    return ax

def plot_mesh(ax, coors, conn, edges, color='k', **plot_kwargs):
    """
    Plot a finite element mesh as a wireframe.
    """
    dim = coors.shape[1]
    ax = _get_axes(ax, dim)
    coors = _to2d(coors)

    for el in conn:
        eds = el[edges]

        for ed in eds:
            cc = coors[ed]

            ax.plot(*cc.T, color=color, **plot_kwargs)

    return ax

def plot_global_dofs(ax, coors, econn):
    """
    Plot global DOF numbers given in an extended connectivity.

    The DOF numbers are plotted for each element, so on common facets they are
    plotted several times - this can be used to check the consistency of the
    global DOF connectivity.
    """
    dim = coors.shape[1]
    ax = _get_axes(ax, dim)
    coors = _to2d(coors)

    for el in econn:
        for gdof in el:
            ax.text(*coors[gdof], s='%d' % gdof,
                    color='g', fontsize=12, weight='bold')

    return ax

def plot_local_dofs(ax, coors, econn):
    """
    Plot local DOF numbers corresponding to an extended connectivity.
    """
    dim = coors.shape[1]
    ax = _get_axes(ax, dim)
    coors = _to2d(coors)

    eps = 0.1
    oeps = 1.0 - eps
    for el in econn:
        # Element centre.
        centre = coors[el].sum(0) / el.shape[0]

        for ldof, gdof in enumerate(el):
            # Shift labels towards the centre.
            cc = oeps * coors[gdof] + eps * centre

            ax.text(*cc, s='%d' % ldof,
                    color='b', fontsize=10, weight='light')

    return ax

def plot_nodes(ax, coors, econn, ref_nodes, dofs):
    """
    Plot Lagrange reference element nodes corresponding to global DOF numbers
    given in an extended connectivity.
    """
    dim = coors.shape[1]
    ax = _get_axes(ax, dim)
    coors = _to2d(coors)

    eps = 0.2
    oeps = 1.0 - eps
    for el in econn:
        # Element centre.
        centre = coors[el].sum(0) / el.shape[0]

        for gdof in dofs:
            if not gdof in el:
                continue
            ldof = nm.where(el == gdof)[0]
            node = ref_nodes[ldof]

            # Shift labels towards the centre.
            cc = oeps * coors[gdof] + eps * centre

            ax.text(*cc, s='%s' % node,
                    color='r', fontsize=8, weight='light')

    return ax
