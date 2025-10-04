"""
Functions to visualize the CMesh geometry and topology.
"""
from sfepy.postprocess.plot_dofs import _get_axes, _to2d

def plot_wireframe(ax, cmesh, color='k', **kwargs):
    """
    Plot a finite element mesh as a wireframe using edges connectivity.
    """
    coors = cmesh.coors
    coors = _to2d(coors)
    dim = cmesh.dim

    ax = _get_axes(ax, dim)

    cmesh.setup_connectivity(1, 0)
    edges = cmesh.get_conn(1, 0)
    for edge_vertices in edges.indices.reshape((edges.num, 2)):
        cc = coors[edge_vertices]
        ax.plot(*cc.T, color=color, **kwargs)

    return ax

def plot_entities(ax, cmesh, edim, color='b', size=10, **kwargs):
    """
    Plot mesh topology entities using scatter plot.
    """
    coors = cmesh.get_centroids(edim)
    coors = _to2d(coors)
    dim = cmesh.dim

    ax = _get_axes(ax, dim)
    ax.scatter(*coors.T, s=size, c=color, **kwargs)

    return ax

def label_global_entities(ax, cmesh, edim, color='b', fontsize=10, **kwargs):
    """
    Label mesh topology entities using global ids.
    """
    coors = cmesh.get_centroids(edim)
    coors = _to2d(coors)
    dim = cmesh.dim

    ax = _get_axes(ax, dim)

    for ii, cc in enumerate(coors):
        ax.text(*cc.T, s=ii, color=color, fontsize=fontsize,
                horizontalalignment='center', verticalalignment='center',
                **kwargs)

    return ax

def label_local_entities(ax, cmesh, edim, color='b', fontsize=10, **kwargs):
    """
    Label mesh topology entities using cell-local ids.
    """
    coors = cmesh.get_centroids(edim)
    coors = _to2d(coors)
    dim = cmesh.dim
    centres = cmesh.get_centroids(dim)
    centres = _to2d(centres)

    cmesh.setup_connectivity(dim, edim)
    conn = cmesh.get_conn(dim, edim)
    off = conn.offsets

    ax = _get_axes(ax, dim)

    eps = 0.015 * fontsize
    oeps = 1.0 - eps
    if dim == 1:
        centres[:, 1] -= eps
    for ii in range(conn.num):
        for ic, ie in enumerate(conn.indices[off[ii]:off[ii+1]]):
            # Shift labels towards the cell centre.
            cc = oeps * coors[ie] + eps * centres[ii]
            ax.text(*cc.T, s=ic, color=color, fontsize=fontsize,
                    horizontalalignment='center', verticalalignment='center',
                    **kwargs)

    return ax

def plot_cmesh(ax, cmesh, wireframe_opts=None, entities_opts=None):
    """
    Convenience function for plotting all entities of a finite element mesh.

    Pass `plot()` arguments to `wireframe_opts` dict.

    Pass `'color'`, `'label_global'`, `'label_global'` for `text()` color and
    font sizes arguments and `'size'` for `scatter()` to each dict for
    topological entities in `entities_opts` list.

    Examples
    --------
    >>> # 2D mesh.
    >>> plot_cmesh(None, cmesh,
                   wireframe_opts = {'color' : 'k', 'linewidth' : 2},
                   entities_opts=[
          {'color' : 'k', 'label_local' : 8, 'size' : 20},
          {'color' : 'b', 'label_global' : 12, 'label_local' : 8, 'size' : 10},
          {'color' : 'r', 'label_global' : 12, 'size' : 20},
          ])
    """
    dim = cmesh.dim
    ax = _get_axes(ax, dim)

    if wireframe_opts is None: wireframe_opts = {}
    if wireframe_opts.get('color') is not None:
        ax = plot_wireframe(ax, cmesh, **wireframe_opts)

    if entities_opts is None: entities_opts = [{}] * dim
    for ii, entity_opts in enumerate(entities_opts):
        if entity_opts.get('color') is not None:
            fsg = entity_opts.pop('label_global', 0)
            fsl = entity_opts.pop('label_local', 0)
            size = entity_opts.pop('size', 10)
            ax = plot_entities(ax, cmesh, ii, size=size, **entity_opts)
            if fsg:
                ax = label_global_entities(ax, cmesh, ii, fontsize=fsg,
                                           **entity_opts)
            if fsl:
                ax = label_local_entities(ax, cmesh, ii, fontsize=fsl,
                                          **entity_opts)

    return ax
