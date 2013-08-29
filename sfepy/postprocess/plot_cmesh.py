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
            ax.plot(cc[:, 0], cc[:, 1], cc[:, 2], color)

        else:
            ax.plot(cc[:, 0], cc[:, 1], color)

    if show:
        plt.show()

    return ax
