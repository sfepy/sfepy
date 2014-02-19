from sfepy.discrete.fem import Domain
import scipy.sparse as sps
import numpy as nm
from sfepy.base.compat import factorial
from sfepy.base.base import output
from numpy.core import intc
from numpy.linalg import lapack_lite

def elems_q2t(el):

    nel, nnd = el.shape
    if nnd > 4:
        q2t = nm.array([[0, 2, 3, 6],
                        [0, 3, 7, 6],
                        [0, 7, 4, 6],
                        [0, 5, 6, 4],
                        [1, 5, 6, 0],
                        [1, 6, 2, 0]])

    else:
        q2t = nm.array([[0, 1, 2],
                        [0, 2, 3]])

    ns, nn = q2t.shape
    nel *= ns

    out = nm.zeros((nel, nn), dtype=nm.int32);

    for ii in range(ns):
        idxs = nm.arange(ii, nel, ns)

        out[idxs,:] = el[:, q2t[ii,:]]

    return nm.ascontiguousarray(out)

def smooth_mesh(mesh, n_iter=4, lam=0.6307, mu=-0.6347,
                weights=None, bconstr=True,
                volume_corr=False):
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
    volume_corr: logical, optional
        Correct volume after smoothing process.

    Returns
    -------
    coors : array
        Coordinates of mesh nodes.
    """

    def laplacian(coors, weights):

        n_nod = coors.shape[0]
        displ = (weights - sps.identity(n_nod)) * coors

        return displ

    def taubin(coors0, weights, lam, mu, n_iter):

        coors = coors0.copy()

        for ii in range(n_iter):
            displ = laplacian(coors, weights)
            if nm.mod(ii, 2) == 0:
                coors += lam * displ
            else:
                coors += mu * displ

        return coors

    def get_volume(el, nd):
        from sfepy.linalg.utils import dets_fast

        dim = nd.shape[1]
        nnd = el.shape[1]

        etype = '%d_%d' % (dim, nnd)
        if etype == '2_4' or etype == '3_8':
            el = elems_q2t(el)

        nel = el.shape[0]

        #bc = nm.zeros((dim, ), dtype=nm.double)
        mul = 1.0 / factorial(dim)
        if dim == 3:
            mul *= -1.0

        mtx = nm.ones((nel, dim + 1, dim + 1), dtype=nm.double)
        mtx[:,:,:-1] = nd[el,:]
        vols = mul * dets_fast(mtx.copy())
        vol = vols.sum()
        bc = nm.dot(vols, mtx.sum(1)[:,:-1] / nnd)

        bc /= vol

        return vol, bc

    import time

    output('smoothing...')
    tt = time.clock()

    if weights is None:
        n_nod = mesh.n_nod
        domain = Domain('mesh', mesh)
        cmesh = domain.cmesh

        # initiate all vertices as inner - hierarchy = 2
        node_group = nm.ones((n_nod,), dtype=nm.int16) * 2
        # boundary vertices - set hierarchy = 4
        if bconstr:
            # get "vertices of surface"
            facets = cmesh.get_surface_facets()
            f_verts = cmesh.get_incident(0, facets, cmesh.dim - 1)
            node_group[f_verts] = 4

        # generate costs matrix
        e_verts = cmesh.get_conn(1, 0).indices
        fc1, fc2 = e_verts[0::2], e_verts[1::2]
        idxs = nm.where(node_group[fc2] >= node_group[fc1])
        rows1 = fc1[idxs]
        cols1 = fc2[idxs]
        idxs = nm.where(node_group[fc1] >= node_group[fc2])
        rows2 = fc2[idxs]
        cols2 = fc1[idxs]
        crows = nm.concatenate((rows1, rows2))
        ccols = nm.concatenate((cols1, cols2))
        costs = sps.coo_matrix((nm.ones_like(crows), (crows, ccols)),
                               shape=(n_nod, n_nod),
                               dtype=nm.double)

        # generate weights matrix
        idxs = range(n_nod)
        aux = sps.coo_matrix((1.0 / nm.asarray(costs.sum(1)).squeeze(),
                              (idxs, idxs)),
                             shape=(n_nod, n_nod),
                             dtype=nm.double)

        #aux.setdiag(1.0 / costs.sum(1))
        weights = (aux.tocsc() * costs.tocsc()).tocsr()

    coors = taubin(mesh.coors, weights, lam, mu, n_iter)

    output('...done in %.2f s' % (time.clock() - tt))

    if volume_corr:
        output('rescaling...')
        volume0, bc = get_volume(mesh.conns[0], mesh.coors)
        volume, _ = get_volume(mesh.conns[0], coors)

        scale = volume0 / volume
        output('scale factor: %.2f' % scale)

        coors = (coors - bc) * scale + bc

        output('...done in %.2f s' % (time.clock() - tt))

    return coors
