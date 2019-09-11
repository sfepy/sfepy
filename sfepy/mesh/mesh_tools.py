from __future__ import absolute_import
from sfepy.discrete.fem import Mesh, FEDomain
import scipy.sparse as sps
import numpy as nm
from sfepy.base.compat import factorial
from sfepy.base.base import output
from six.moves import range

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

def triangulate(mesh, verbose=False):
    """
    Triangulate a 2D or 3D tensor product mesh: quadrilaterals->triangles,
    hexahedrons->tetrahedrons.

    Parameters
    ----------
    mesh : Mesh
        The input mesh.

    Returns
    -------
    mesh : Mesh
        The triangulated mesh.
    """
    conns = None
    for k, new_desc in [('3_8', '3_4'), ('2_4', '2_3')]:
        if k in mesh.descs:
            conns = mesh.get_conn(k)
            break

    if conns is not None:
        nelo = conns.shape[0]
        output('initial mesh: %d elements' % nelo, verbose=verbose)

        new_conns = elems_q2t(conns)
        nn = new_conns.shape[0] // nelo
        new_cgroups = nm.repeat(mesh.cmesh.cell_groups, nn)

        output('new mesh: %d elements' % new_conns.shape[0], verbose=verbose)
        mesh = Mesh.from_data(mesh.name, mesh.coors,
                              mesh.cmesh.vertex_groups,
                              [new_conns], [new_cgroups], [new_desc])

    return mesh

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
        vols = mul * dets_fast(mtx)
        vol = vols.sum()
        bc = nm.dot(vols, mtx.sum(1)[:,:-1] / nnd)

        bc /= vol

        return vol, bc

    from sfepy.base.timing import Timer

    output('smoothing...')
    timer = Timer(start=True)

    if weights is None:
        n_nod = mesh.n_nod
        domain = FEDomain('mesh', mesh)
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
        idxs = list(range(n_nod))
        aux = sps.coo_matrix((1.0 / nm.asarray(costs.sum(1)).squeeze(),
                              (idxs, idxs)),
                             shape=(n_nod, n_nod),
                             dtype=nm.double)

        #aux.setdiag(1.0 / costs.sum(1))
        weights = (aux.tocsc() * costs.tocsc()).tocsr()

    coors = taubin(mesh.coors, weights, lam, mu, n_iter)

    output('...done in %.2f s' % timer.stop())

    if volume_corr:
        output('rescaling...')
        timer.start()
        volume0, bc = get_volume(mesh.conns[0], mesh.coors)
        volume, _ = get_volume(mesh.conns[0], coors)

        scale = volume0 / volume
        output('scale factor: %.2f' % scale)

        coors = (coors - bc) * scale + bc

        output('...done in %.2f s' % timer.stop())

    return coors

def expand2d(mesh2d, dist, rep):
    """
    Expand 2D planar mesh into 3D volume,
    convert triangular/quad mesh to tetrahedrons/hexahedrons.

    Parameters
    ----------
    mesh2d : Mesh
        The 2D mesh.
    dist : float
        The elements size in the 3rd direction.
    rep : int
        The number of elements in the 3rd direction.

    Returns
    -------
    mesh3d : Mesh
        The 3D mesh.
    """
    if len(mesh2d.descs) > 1:
        raise ValueError('More than one cell type (%s). Not supported!'
                         % ', '.join(mesh2d.descs))

    nel = mesh2d.n_el
    nnd = mesh2d.n_nod
    et = mesh2d.descs[0]
    coors = mesh2d.coors
    conn = mesh2d.get_conn(et)

    zcoor = nm.arange(rep + 1) * dist
    coors3d = nm.hstack([nm.tile(coors, (rep + 1, 1)),
                         nm.tile(zcoor, (nnd,1)).T.flatten()[:,nm.newaxis]])
    ngroups = nm.tile(mesh2d.cmesh.vertex_groups, (rep + 1,))

    if et == '2_4':
        descs3d = '3_8'
        conn3d = nm.zeros((nel * rep, 8), dtype=nm.int32)
        mats3d = nm.tile(mesh2d.cmesh.cell_groups, (1, rep)).squeeze()

    elif et == '2_3':
        descs3d = '3_4'
        conn3d = nm.zeros((3 * nel * rep, 4), dtype=nm.int32)
        mats3d = nm.tile(mesh2d.cmesh.cell_groups, (1, 3 * rep)).squeeze()

    for ii in range(rep):
        bgn0 = nnd * ii
        bgn1 = bgn0 + nnd
        if et == '2_4':
            bge0 = nel * ii
            bge1 = bge0 + nel
            conn3d[bge0:bge1,:4] = conn + bgn0
            conn3d[bge0:bge1,4:] = conn + bgn1

        elif et == '2_3':
            # 0 1 2 5
            bge0 = 3 * nel * ii
            bge1 = bge0 + nel
            conn3d[bge0:bge1,:] = nm.array([conn[:,0] + bgn0,
                                            conn[:,1] + bgn0,
                                            conn[:,2] + bgn0,
                                            conn[:,2] + bgn1]).T
            # 0 1 5 4
            bge0 += nel
            bge1 += nel
            conn3d[bge0:bge1,:] = nm.array([conn[:,0] + bgn0,
                                            conn[:,1] + bgn0,
                                            conn[:,2] + bgn1,
                                            conn[:,1] + bgn1]).T
            # 0 4 5 3
            bge0 += nel
            bge1 += nel
            conn3d[bge0:bge1,:] = nm.array([conn[:,0] + bgn0,
                                            conn[:,1] + bgn1,
                                            conn[:,2] + bgn1,
                                            conn[:,0] + bgn1]).T

    mesh3d = Mesh.from_data('mesh', coors3d, ngroups, [conn3d],
                            [mats3d], [descs3d])

    return mesh3d
