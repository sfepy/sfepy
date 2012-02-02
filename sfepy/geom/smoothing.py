from sfepy.fem import Domain
import scipy.sparse as sps
import numpy as nm

def smooth_mesh(mesh, n_iter=4, lam=0.6307, mu=-0.6347,
                weights=None, bconstr=True):
    """
    FE mesh smoothing.

    Based on:

    [1] Steven K. Boyd, Ralph Muller, Smooth surface meshing for automated
    finite element model generation from 3D image data, Journal of
    Biomechanics, Volume 39, Issue 7, 2006, Pages 1287-1295,
    ISSN 0021-9290, 10.1016/j.jbiomech.2005.03.006.
    (http://www.sciencedirect.com/science/article/pii/S0021929005001442)

    Parameters
    ----------
    mesh : mesh
        FE mesh.
    n_iter : integer, optional
        Number of iteration steps.
    lam : float, optional
        Smoothing factor, see [1].
    mu : float, optional
        Unshrinking factor, see [1].
    weights : array, optional
        Edge weights, see [1].
    bconstr: logical, optional
        Boundary constraints, if True only surface smoothing performed.

    Returns
    -------
    coors : array
        Coordinates of mesh nodes.
    """

    domain = Domain('mesh', mesh)

    n_nod = mesh.n_nod
    edges = domain.ed

    if weights is None:
        # initialize all vertices as inner - hierarchy = 2
        node_group = nm.ones((n_nod,), dtype=nm.int16) * 2
        # boundary vertices - set hierarchy = 4
        if bconstr:
            # get "nodes of surface"
            if domain.fa: # 3D.
                fa = domain.fa
            else:
                fa = domain.ed

            flag = fa.mark_surface_facets()
            ii = nm.where(flag > 0)[0]
            aux = nm.unique(fa.facets[ii])
            if aux[0] == -1: # Triangular faces have -1 as 4. point.
                aux = aux[1:]

            node_group[aux] = 4

        # generate costs matrix
        costs = sps.lil_matrix((n_nod, n_nod), dtype=nm.double)
        for ied in range(edges.mtx.shape[0]):
            cc = edges.mtx.getrow(ied).indices[0]
            n1, n2 = edges.facets[cc,:]
            if node_group[n2] >= node_group[n1]:
                costs[n1, n2] = 1.0

            if node_group[n1] >= node_group[n2]:
                costs[n2, n1] = 1.0

        # generate weights matrix
        aux = sps.lil_matrix((n_nod, n_nod), dtype=nm.double)
        aux.setdiag(1.0 / costs.sum(1))
        weights = (aux * costs).tocsr()

    coors = tsmooth(mesh.coors, weights, lam, mu, n_iter)

    return coors

def laplacian(coors, weights):

    n_nod = coors.shape[0]
    displ = (weights - sps.identity(n_nod)) * coors

    return displ

def tsmooth(coors0, weights, lam, mu, n_iter):

    coors = coors0.copy()

    for ii in range(n_iter):
        displ = laplacian(coors, weights)
        if nm.mod(ii, 2) == 0:
            coors += lam * displ;
        else:
            coors += mu * displ;

    return coors
