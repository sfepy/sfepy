import numpy as nm
import scipy.sparse as sps
from scipy.sparse.csgraph import connected_components

from sfepy.base.compat import factorial
from sfepy.base.base import output
from sfepy.discrete.equations import create_dof_graph
from sfepy.discrete.fem import Mesh, FEDomain

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

def merge_lines(mesh, eps=1e-18):
    """
    Merge edges of an edge-only mesh that are in the same direction w.r.t. the
    tolerance `eps`.
    """
    coors, ngroups, conns, mat_ids, ctype = mesh
    conns = conns[0]

    # vertices to edges map
    n_v = coors.shape[0]
    n_e = conns.shape[0]
    row = nm.repeat(nm.arange(n_e), 2)
    aux = sps.coo_matrix((nm.ones((n_e * 2,), dtype=bool),
                          (row, conns.flatten())), shape=(n_e, n_v))
    v2e = aux.tocsc()
    n_epv = nm.diff(v2e.indptr)

    # directional vectors of edges
    de = coors[conns[:, 1], :] - coors[conns[:, 0], :]
    de = de / nm.linalg.norm(de, axis=1)[:, nm.newaxis]

    eflag = nm.ones((n_e, ), dtype=bool)
    valid_e = nm.where(eflag)[0]
    e_remove = []

    while len(valid_e) > 0:
        ie = valid_e[0]
        d = de[ie]
        buff = [(ie, conns[ie, 0]), (ie, conns[ie, 1])]
        eflag[ie] = False  # invalidate edge
        while len(buff) > 0:
            e, v = buff.pop(-1)
            if n_epv[v] == 2:
                idx = v2e.indptr[v]
                aux = v2e.indices[idx]
                next_e = v2e.indices[idx + 1] if aux == e else aux

                if not eflag[next_e]:  # valid edge?
                    continue
                if nm.linalg.norm(de[next_e] - d) < eps\
                        or nm.linalg.norm(de[next_e] + d) < eps:
                    next_ec = conns[next_e, :]
                    new_v = next_ec[0] if next_ec[1] == v else next_ec[1]
                    idx = 0 if conns[e, 0] == v else 1
                    conns[e, idx] = new_v  # reconnect edge

                    idx = v2e.indptr[new_v]
                    aux = v2e.indices[idx]
                    idx += 0 if aux == next_e else 1
                    v2e.indices[idx] = e  # update v2e map

                    buff.append((e, new_v))  # continue in searching
                    eflag[next_e] = False  # invalidate edge
                    e_remove.append(next_e)

        valid_e = nm.where(eflag)[0]

    if len(e_remove) > 0:
        # remove unused edges and vertices
        eflag.fill(True)
        eflag[nm.asarray(e_remove)] = False
        remap = -nm.ones((n_v, ), dtype=nm.int64)
        remap[conns[eflag, :]] = 1
        vidx = nm.where(remap > 0)[0]
        remap[vidx] = nm.arange(len(vidx))
        conns_new = remap[conns[eflag, :]]

        return coors[vidx, :], ngroups[vidx],\
            [conns_new], [mat_ids[0][eflag]], ctype
    else:
        return mesh

def extract_edges(mesh, eps=1e-16):
    """
    Extract outline edges of a given mesh.
    The outline edge is an edge for which norm(nvec_1 - nvec_2) < eps,
    where nvec_1 and nvec_2 are the normal vectors of the incident facets.

    Parameters
    ----------
    mesh : Mesh
        The 3D or 2D mesh.
    eps : float
        The tolerance parameter of the outline edge searching algorithm.

    Returns
    -------
    mesh_out : tuple
        The data of the outline mesh, Mesh.from_data() format, i.e.
        (coors, ngroups, ed_conns, mat_ids, descs).
    """
    domain = FEDomain('domain', mesh)
    cmesh = domain.cmesh

    output('Mesh - dimension: %d, vertices: %d, elements: %d'
           % (mesh.dim, mesh.n_nod, mesh.n_el))

    if mesh.dim == 2:
        oedges = cmesh.get_surface_facets()
        mesh_coors = nm.hstack([cmesh.coors,
                                nm.zeros((cmesh.coors.shape[0], 1))])

    elif mesh.dim == 3:
        cmesh.setup_connectivity(1, 2)
        cmesh.setup_connectivity(3, 2)

        sfaces = cmesh.get_surface_facets()
        _, idxs = nm.unique(cmesh.get_conn(3, 2).indices, return_index=True)

        normals = cmesh.get_facet_normals()[idxs, :]
        se_map, se_off = cmesh.get_incident(1, sfaces, 2, ret_offsets=True)
        sedges = nm.unique(se_map)
        n_se = sedges.shape[0]

        # remap surface edges to continuous range
        se_remap = -nm.ones(sedges.max() + 1)
        se_remap[sedges] = nm.arange(n_se)
        se_map0 = se_remap[se_map]

        # surface face/edge connectivity matrix (n_surf x n_edge)
        n_ef = nm.diff(se_off)[0]  # = 2
        n_sf = se_map.shape[0] // n_ef
        row = nm.repeat(nm.arange(n_sf), n_ef)
        sf2e = sps.coo_matrix((nm.ones((n_sf * n_ef,), dtype=bool),
                               (row , se_map0)), shape=(n_sf, n_se))
        # edge to face map (n_edge x 2)
        se2f = sf2e.tocsc().indices.reshape((sedges.shape[0], 2))

        snormals = normals[sfaces]
        err = nm.linalg.norm(snormals[se2f[:, 0]] - snormals[se2f[:, 1]],
                             axis=1)
        oedges = sedges[nm.where(err > eps)[0]]
        mesh_coors = cmesh.coors

    else:
        raise NotImplementedError

    # save outline mesh
    if oedges.shape[0] > 0:
        ec_idxs = nm.unique(cmesh.get_incident(0, oedges, 1))
        ed_coors = mesh_coors[ec_idxs, :]
        ngroups = nm.zeros((ed_coors.shape[0],), dtype=nm.int16)

        aux = cmesh.get_conn(1, 0).indices

        ed_conns = aux.reshape((aux.shape[0] // 2, 2))[oedges, :]
        ec_remap = -nm.ones((ec_idxs.max() + 1, ), dtype=nm.int64)
        ec_remap[ec_idxs] = nm.arange(ec_idxs.shape[0])
        ed_conns = ec_remap[ed_conns]
        mat_ids = nm.ones((ed_conns.shape[0],), dtype=nm.int16)

        mesh_out = ed_coors, ngroups, [ed_conns], [mat_ids], ['3_2']

        return mesh_out

    else:
        raise ValueError('no outline edges found (eps=%e)!' % eps)

def _get_facets(vertices, offsets, ii, n_fp):
    facets = []
    for ic in range(n_fp):
        facets.append(vertices[offsets[ii] + ic][:, None])

    facets = nm.concatenate(facets, axis=1)

    return nm.ascontiguousarray(facets.astype(nm.int32))

def get_surface_faces(domain):
    cmesh = domain.cmesh
    faces = cmesh.get_surface_facets()
    vertices_f, offs_f = cmesh.get_incident(0, faces,
                                            cmesh.dim - 1, ret_offsets=True)

    n_fp = nm.diff(offs_f)
    surf_faces = []

    itri = nm.where(n_fp == 3)[0]
    if itri.size:
        surf_faces.append(_get_facets(vertices_f, offs_f, itri, 3))

    itet = nm.where(n_fp == 4)[0]
    if itet.size:
        surf_faces.append(_get_facets(vertices_f, offs_f, itet, 4))

    cells_c, offs_c = cmesh.get_incident(cmesh.dim, faces, cmesh.dim - 1,
                                         ret_offsets=True)
    ids = cmesh.get_local_ids(faces, cmesh.dim - 1, cells_c, offs_c,
                              cmesh.dim)
    lst = nm.c_[cells_c, ids]

    return lst, surf_faces

def surface_graph(surf_faces, n_nod):
    graph = 0
    for conn in surf_faces:
        graph += create_dof_graph(conn, conn, (n_nod, n_nod)).tocsr()

    return graph

def surface_components(gr_s, surf_faces):
    """
    Determine surface components given surface mesh connectivity graph.
    """
    n_comp, flag = connected_components(gr_s)

    comps = []
    for ii, face in enumerate(surf_faces):
        comp = flag[face[:,0]]
        comps.append(comp)

    return n_comp, comps
