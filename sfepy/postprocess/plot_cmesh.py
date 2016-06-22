"""
Functions to visualize the CMesh geometry and topology.
"""
from __future__ import absolute_import
from sfepy.postprocess.plot_dofs import _get_axes, _to2d
from six.moves import range

def plot_wireframe(ax, cmesh, color='k'):
    """
    Plot a finite element mesh as a wireframe using edges connectivity.
    """
    coors = cmesh.coors
    coors = _to2d(coors)
    dim = cmesh.dim

    ax = _get_axes(ax, dim)

    edges = cmesh.get_conn(1, 0)
    for edge_vertices in edges.indices.reshape((edges.num, 2)):
        cc = coors[edge_vertices]
        ax.plot(*cc.T, color=color)

    return ax

def plot_entities(ax, cmesh, edim, color='b', size=10):
    """
    Plot mesh topology entities using scatter plot.
    """
    coors = cmesh.get_centroids(edim)
    coors = _to2d(coors)
    dim = cmesh.dim

    ax = _get_axes(ax, dim)
    ax.scatter(*coors.T, s=size, c=color)

    return ax

def label_global_entities(ax, cmesh, edim, color='b', fontsize=10):
    """
    Label mesh topology entities using global ids.
    """
    coors = cmesh.get_centroids(edim)
    coors = _to2d(coors)
    dim = cmesh.dim

    ax = _get_axes(ax, dim)

    for ii, cc in enumerate(coors):
        ax.text(*cc.T, s=ii, color=color, fontsize=fontsize)

    return ax

def label_local_entities(ax, cmesh, edim, color='b', fontsize=10):
    """
    Label mesh topology entities using cell-local ids.
    """
    coors = cmesh.get_centroids(edim)
    coors = _to2d(coors)
    dim = cmesh.dim
    centres = cmesh.get_centroids(dim)

    cmesh.setup_connectivity(dim, edim)
    conn = cmesh.get_conn(dim, edim)
    off = conn.offsets

    ax = _get_axes(ax, dim)

    eps = 0.1
    oeps = 1.0 - eps
    for ii in range(conn.num):
        for ic, ie in enumerate(conn.indices[off[ii]:off[ii+1]]):
            # Shift labels towards the cell centre.
            cc = oeps * coors[ie] + eps * centres[ii]
            ax.text(*cc.T, s=ic, color=color, fontsize=fontsize)

    return ax
