"""
Functions to visualize the CMesh geometry and topology.
"""
import matplotlib.pyplot as plt

from sfepy.postprocess.plot_dofs import _get_axes

def plot_wireframe(ax, cmesh, color='k', show=False):
    """
    Plot a finite element mesh as a wireframe using edges connectivity.
    """
    coors = cmesh.coors
    dim = cmesh.dim

    edges = cmesh.get_conn(1, 0)

    ax = _get_axes(ax, dim)

    for edge_vertices in edges.indices.reshape((edges.num, 2)):
        cc = coors[edge_vertices]
        if dim == 3:
            ax.plot(cc[:, 0], cc[:, 1], cc[:, 2], color=color)

        else:
            ax.plot(cc[:, 0], cc[:, 1], color=color)

    if show:
        plt.show()

    return ax

def plot_entities(ax, cmesh, edim, color='b', size=10, show=False):
    """
    Plot mesh topology entities using scatter plot.
    """
    coors = cmesh.get_centroids(edim)
    dim = cmesh.dim

    ax = _get_axes(ax, dim)

    if dim == 3:
        ax.scatter(coors[:, 0], coors[:, 1], coors[:, 2], s=size, c=color)

    else:
        ax.scatter(coors[:, 0], coors[:, 1], s=size, c=color)

    if show:
        plt.show()

    return ax

def label_global_entities(ax, cmesh, edim, color='b', fontsize=10, show=False):
    """
    Label mesh topology entities using global ids.
    """
    coors = cmesh.get_centroids(edim)
    dim = cmesh.dim

    ax = _get_axes(ax, dim)

    for ii, cc in enumerate(coors):
        if dim == 3:
            ax.text(cc[0], cc[1], cc[2], ii,
                    color=color, fontsize=fontsize)

        else:
            ax.text(cc[0], cc[1], ii,
                    color=color, fontsize=fontsize)

    if show:
        plt.show()

    return ax

def label_local_entities(ax, cmesh, edim, color='b', fontsize=10, show=False):
    """
    Label mesh topology entities using cell-local ids.
    """
    coors = cmesh.get_centroids(edim)
    dim = cmesh.dim
    centres = cmesh.get_centroids(dim)

    conn = cmesh.get_conn(dim, edim)
    off = conn.offsets

    ax = _get_axes(ax, dim)

    eps = 0.1
    oeps = 1.0 - eps
    for ii in xrange(conn.num):
        for ic, ie in enumerate(conn.indices[off[ii]:off[ii+1]]):
            # Shift labels towards the cell centre.
            cc = oeps * coors[ie] + eps * centres[ii]
            if dim == 3:
                ax.text(cc[0], cc[1], cc[2], ic,
                        color=color, fontsize=fontsize)

            else:
                ax.text(cc[0], cc[1], ic,
                        color=color, fontsize=fontsize)

    if show:
        plt.show()

    return ax
