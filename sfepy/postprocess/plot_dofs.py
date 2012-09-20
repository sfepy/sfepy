"""
Functions to visualize the mesh connectivity with global and local DOF
numberings.
"""
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

def plot_mesh(ax, coors, conn, edges, show=False):
    """
    Plot a finite element mesh as a wireframe.
    """
    dim = coors.shape[1]
    ax = _get_axes(ax, dim)

    for el in conn:
        eds = el[edges]

        for ed in eds:
            cc = coors[ed]

            if dim == 3:
                ax.plot(cc[:, 0], cc[:, 1], cc[:, 2], 'k')

            else:
                ax.plot(cc[:, 0], cc[:, 1], 'k')

    if show:
        plt.show()

    return ax

def plot_global_dofs(ax, coors, econn, show=False):
    """
    Plot global DOF numbers given in an extended connectivity.

    The DOF numbers are plotted for each element, so on common facets they are
    plotted several times - this can be used to check the consistency of the
    global DOF connectivity.
    """
    dim = coors.shape[1]
    ax = _get_axes(ax, dim)

    for el in econn:
        for gdof in el:
            if dim == 3:
                cx, cy, cz = coors[gdof]
                ax.text(cx, cy, cz, '%d' % gdof,
                        color='g', fontsize=12, weight='bold')

            else:
                cx, cy = coors[gdof]
                ax.text(cx, cy, '%d' % gdof,
                        color='g', fontsize=12, weight='bold')

    if show:
        plt.show()

    return ax

def plot_local_dofs(ax, coors, econn, show=False):
    """
    Plot local DOF numbers corresponding to an extended connectivity.
    """
    dim = coors.shape[1]
    ax = _get_axes(ax, dim)

    eps = 0.1
    oeps = 1.0 - eps
    for el in econn:
        # Element centre.
        centre = coors[el].sum(0) / el.shape[0]

        for ldof, gdof in enumerate(el):
            # Shift labels towards the centre.
            cc = oeps * coors[gdof] + eps * centre

            if dim == 3:
                cx, cy, cz = cc
                ax.text(cx, cy, cz, '%d' % ldof,
                        color='b', fontsize=10, weight='light')

            else:
                cx, cy = cc
                ax.text(cx, cy, '%d' % ldof,
                        color='b', fontsize=10, weight='light')

    if show:
        plt.show()

    return ax
