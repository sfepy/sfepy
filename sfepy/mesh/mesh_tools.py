import numpy as nm
import scipy.sparse as sps
from scipy.sparse.csgraph import connected_components
from scipy.spatial.transform import Rotation
from scipy.spatial import cKDTree

from sfepy.base.compat import factorial
from sfepy.base.base import output
from sfepy.discrete.equations import create_dof_graph
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.discrete.fem.utils import prepare_translate

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

    out = nm.zeros((nel, nn), dtype=nm.int32)

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
        displ = (weights - sps.eye_array(n_nod)) @ coors

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
        costs = sps.coo_array((nm.ones_like(crows), (crows, ccols)),
                              shape=(n_nod, n_nod),
                              dtype=nm.double)

        # generate weights matrix
        idxs = list(range(n_nod))
        aux = sps.coo_array((1.0 / nm.asarray(costs.sum(1)).squeeze(),
                             (idxs, idxs)),
                            shape=(n_nod, n_nod),
                            dtype=nm.double)

        #aux.setdiag(1.0 / costs.sum(1))
        weights = (aux.tocsc() @ costs.tocsc()).tocsr()

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
    aux = sps.coo_array((nm.ones((n_e * 2,), dtype=bool),
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
        sf2e = sps.coo_array((nm.ones((n_sf * n_ef,), dtype=bool),
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


def get_cell_vertices_only(mesh):
    """
    Remove vertices not used in any cell.

    Parameters
    ----------
    mesh: Mesh
        FE mesh

    Returns
    -------
    out: Mesh
        FE mesh
    """
    data = list(mesh._get_io_data())
    vertices = nm.unique([nm.unique(conn.ravel()) for conn in data[2]])

    remap = prepare_translate(vertices, nm.arange(len(vertices)))

    coors = data[0][vertices]
    ngroups = data[1][vertices]
    conns = [remap[conn] for conn in data[2]]

    return Mesh.from_data(mesh.name + '_cleaned', coors, ngroups,
                          conns, data[3], data[4])


def get_mesh_by_cgroup(mesh, value, cell_vertices_only=True):
    """
    Extract mesh cells using cell_group value(s).

    Parameters
    ----------
    mesh: Mesh
        FE mesh
    value: int, list, or tuple
        Select cells with a given value, list, or range: (value[0], value[1])
    cell_vertices_only: bool
        If True, remove free vertices

    Returns
    -------
    out: Mesh
        FE mesh
    """
    new_conns, new_mat_ids, new_descs = [], [], []

    for desc in mesh.descs:
        conns, cidxs = mesh.get_conn(desc, ret_cells=True)
        mat_ids = mesh.cmesh.cell_groups[cidxs]

        if isinstance(value, tuple):
            ncidxs = nm.logical_and(mat_ids >= value[0], mat_ids <= value[1])
        elif isinstance(value, list):
            ncidxs = nm.logical_or.reduce([mat_ids == k for k in value])
        else:
            ncidxs = mat_ids == value

        if ncidxs.sum() > 0:
            new_conns.append(conns[ncidxs])
            new_mat_ids.append(mat_ids[ncidxs])
            new_descs.append(desc)

    out = Mesh.from_data(mesh.name, mesh.coors,
                         mesh.cmesh.vertex_groups,
                         new_conns, new_mat_ids, new_descs)

    if cell_vertices_only:
        out = get_cell_vertices_only(out)

    return out


def get_mesh_by_ngroup(mesh, value, cell_vertices_only=True):
    """
    Extract mesh cells using cell_group value(s).

    Parameters
    ----------
    mesh: Mesh
        FE mesh
    value: int, list, or tuple
        Select nodes with a given value, list, or range: (value[0], value[1])
    cell_vertices_only: bool
        If True, remove free vertices

    Returns
    -------
    out: Mesh
        FE mesh
    """
    vgroups = mesh.cmesh.vertex_groups
    if isinstance(value, tuple):
        vidxs = nm.logical_and(vgroups >= value[0], vgroups <= value[1])
    elif isinstance(value, list):
        vidxs = nm.logical_or.reduce([vgroups == k for k in value])
    else:
        vidxs = vgroups == value

    new_vgroups = vgroups[vidxs]
    new_coors = mesh.coors[vidxs]

    remap = -nm.ones((len(vidxs),), dtype=nm.int64)
    remap[vidxs] = nm.arange(vidxs.sum())

    new_conns, new_mat_ids, new_descs = [], [], []
    for desc in mesh.descs:
        conns, cidxs = mesh.get_conn(desc, ret_cells=True)
        mat_ids = mesh.cmesh.cell_groups[cidxs]
        ncidxs = (remap[conns] >= 0).sum(axis=1) == conns.shape[1]

        if ncidxs.sum() > 0:
            new_conns.append(remap[conns[ncidxs]])
            new_mat_ids.append(mat_ids[ncidxs])
            new_descs.append(desc)

    out = Mesh.from_data(mesh.name, new_coors, new_vgroups,
                         new_conns, new_mat_ids, new_descs)

    if cell_vertices_only:
        out = get_cell_vertices_only(out)

    return out


def stack_descs(conns, mat_ids, descs):
    sdescs = set(descs)
    sconns = {k: [] for k in sdescs}
    smat_ids = {k: [] for k in sdescs}

    for k, desc in enumerate(descs):
        sconns[desc].append(conns[k])
        smat_ids[desc].append(mat_ids[k])

    return ([nm.vstack(sconns[k]) for k in sdescs],
            [nm.hstack(smat_ids[k]) for k in sdescs],
            list(sdescs))


def extrude(mesh, cline, twist=None, scale=None, nvec=None,
            wedge_to_tetra=True):
    """
    Create a solid 3D mesh from a given planar 2D mesh by extruding it.
    The new points in each layer lie in a plane perpendicular to one of the
    cline segments. Each layer can be twisted and scaled.

    Parameters
    ----------
    mesh: Mesh
        2D planar FE mesh (tri or quad elements)
    cline: list of coordinates
        Points of central line.
    twist: float, or array
        Angle of twist in each layer.
    scale: float, or array
        Scale factor in each layer.
    nvec: array
        Normal vectors of layers.
    wedge_to_tetra: bool
        If True, convert wedge elements to tetrahedrons.

    Returns
    -------
    out: Mesh
        3D FE mesh
    """
    def get_rotation(v0, v1):
        cross = nm.cross(v0, v1)
        dot = nm.dot(v0, v1)

        if nm.linalg.norm(cross) < 1e-18:
            # Already aligned or opposite
            if dot > 0:
                rotation = Rotation.identity()
            else:
                # pi/2 rotation around any perpendicular axis
                axis = nm.array([0, 1, 0])
                if nm.allclose(v0, axis):
                    axis = nm.array([1, 0, 0])
                axis = nm.cross(v0, axis)
                axis /= nm.linalg.norm(axis)
                rotation = Rotation.from_rotvec(nm.pi * axis)
        else:
            axis = cross / nm.linalg.norm(cross)
            angle = nm.arccos(nm.clip(dot, -1.0, 1.0))
            rotation = Rotation.from_rotvec(angle * axis)

        return rotation

    coors0 = mesh.coors
    if coors0.shape[1] == 2:
        coors0 = nm.hstack([coors0, coors0[:, [0]] * 0])

    ccoor = nm.mean(coors0, axis=0)
    coors0 = coors0 - ccoor[None, :]

    nvec0 = nm.array([0, 0, 1.])

    cline = nm.array(cline)
    if len(cline.shape) == 1:
        dz, nz = cline[0], int(cline[1])
        cline = nm.tile(ccoor, (nz + 1, 1))
        cline[:, 2] = nm.arange(nz + 1) * dz

    if nvec is not None:
        nvec = nm.array(nvec)
        if len(nvec.shape) == 0:
            nvec = nm.repeat(nvec, len(cline) - 1)
    else:
        nvec = cline[1:] - cline[:-1]

    nvec = nvec / nm.linalg.norm(nvec, axis=1)[:, None]

    if twist is not None:
        twist = nm.array(twist)
        if len(twist.shape) == 0:
            twist = nm.cumsum(nm.repeat(twist, len(cline)))

    if scale is not None:
        scale = nm.array(scale)
        if len(scale.shape) == 0:
            scale = nm.cumsum(nm.repeat(scale, len(cline)))

    new_coors = [coors0 + cline[[0], :]]
    new_conns, new_mat_ids, new_descs = [], [], []
    mat_ids = mesh.cmesh.cell_groups
    nnd = coors0.shape[0]
    new_vgroups = nm.hstack([mesh.cmesh.vertex_groups] * len(cline))

    for k in range(len(cline) - 1):
        rotation = get_rotation(nvec0, nvec[k])
        coors = rotation.apply(coors0)

        if twist is not None:
            rotation = Rotation.from_rotvec(nvec[k] * twist[k])
            coors = rotation.apply(coors)

        if scale is not None:
            coors *= scale[k] 

        new_coors.append(coors + cline[[k + 1], :])

        for desc in mesh.descs:
            conns, cidxs = mesh.get_conn(desc, ret_cells=True)
            econns = nm.hstack([conns + nnd * k, conns + nnd * (k + 1)])
            emat_ids = mat_ids[cidxs]
            if econns.shape[1] == 6:
                if wedge_to_tetra:  # wedges -> 3 x tetrahedron
                    econns1 = nm.empty((3* econns.shape[0], 4),
                                       dtype=econns.dtype)
                    econns1[0::3] = econns[:, [0, 1, 2, 5]]
                    econns1[1::3] = econns[:, [0, 1, 5, 4]]
                    econns1[2::3] = econns[:, [0, 4, 5, 3]]
                    econns = econns1
                    emat_ids = nm.repeat(mat_ids, 3)

            new_conns.append(econns)
            new_mat_ids.append(emat_ids)
            new_descs.append(f'3_{econns.shape[1]}')

    out = Mesh.from_data(mesh.name,
                         nm.vstack(new_coors), new_vgroups,
                         *stack_descs(new_conns, new_mat_ids, new_descs))

    return out


def expand2d(mesh2d, dist, rep):
    """
    Expand a 2D planar mesh into a 3D volume,
    convert triangular/quad elements to tetrahedrons/hexahedrons.

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
    return extrude(mesh2d, (dist, rep))


def get_unique_coor_map(coors, eps=1e-18):
    tree = cKDTree(coors)
    pairs = nm.array(list(tree.query_pairs(r=eps)))

    if len(pairs) > 0:
        n, row, col = len(coors), pairs[:, 0], pairs[:, 1]
        graph = sps.csr_array((nm.ones(len(row)), (row, col)), shape=(n, n))
        _, remap = connected_components(graph + graph.T)

        return remap


def merge_nodes(mesh, eps=1e-12):
    """
    Merge duplicate mesh nodes.

    Parameters
    ----------
    mesh: Mesh
        FE mesh.
    eps: float
        Tolerance for duplicity search.

    Returns
    -------
    mesh: Mesh
        FE mesh.
    """
    remap = get_unique_coor_map(mesh.coors, eps=eps)

    if remap is not None:
        _, nidxs = nm.unique(remap, return_index=True)
        new_conns, new_mat_ids = [], []
        for desc in mesh.descs:
            conns, cidxs = mesh.get_conn(desc, ret_cells=True)
            new_mat_ids.append(mesh.cmesh.cell_groups[cidxs])
            new_conns.append(remap[conns])

        return Mesh.from_data(mesh.name, mesh.coors[nidxs],
                              mesh.cmesh.vertex_groups[nidxs],
                              new_conns, new_mat_ids, mesh.descs)
    else:
        return mesh


def revolve(mesh, nphi, p, v=(1, 0, 0), phi_max=None, wedge_to_tetra=True):
    """
    Create a solid 3D mesh from a given planar 2D mesh by revolving it
    around the axis defined by a vector.

    Parameters
    ----------
    mesh: Mesh
        2D planar FE mesh.
    nphi: int
        Number of elements in the circumferential direction
    p: list, tuple, or numpy.ndarray
        Coordinates of a point on the axis of rotation
    v: list, tuple, or numpy.ndarray
        Directional vector of the axis of rotation
    phi_max: float or None
        Angle of the revolution, if None, 360 degrees is used (closed loop)
    wedge_to_tetra: bool
        If True, convert wedge elements to tetrahedrons.

    Returns
    -------
    out: Mesh
        3D FE mesh.
    """
    if phi_max is None:
        phi_max = 360
        phi = nm.linspace(0, 2*nm.pi, nphi + 1)
    else:
        phi = nm.linspace(0, nm.deg2rad(phi_max), nphi + 1)

    p = nm.array(p)[None, :]
    v = nm.array(v)

    coors = mesh.coors
    ccoor = nm.mean(coors, axis=0)
    if coors.shape[1] == 2:
        ccoor = nm.array([ccoor[0], ccoor[1], 0])

    cline = nm.vstack([Rotation.from_rotvec(v * ph).apply(ccoor - p)
                       for ph in phi]) + p
 
    nvec0 = nm.array([0, 0, 1])
    nvec = nm.vstack([Rotation.from_rotvec(v * ph).apply(nvec0)
                      for ph in phi[1:]])

    out = extrude(mesh, cline, nvec=nvec, wedge_to_tetra=wedge_to_tetra)

    if phi_max == 360:
        out = merge_nodes(out)

    return out


def mirror(mesh, p, v):
    """
    Duplicate the mesh by mirroring it. The mirror plane is defined by a point
    in the plane and by a normal vector to that plane.

    Parameters
    ----------
    mesh: Mesh
        FE mesh.
    p: list, tuple, or numpy.ndarray
        Coordinates of a point at the mirror plane.
    v: list, tuple, or numpy.ndarray
        Normal vector of the mirror plane.

    Returns
    -------
    out: Mesh
        FE mesh.
    """
    mirror_map = {
        '2_3': nm.array([1, 0, 2]),
        '2_4': nm.array([1, 0, 3, 2]),
        '3_4': nm.array([1, 0, 2, 3]),
        '3_8': nm.array([1, 0, 3, 2, 5, 4, 7, 6]),
    }
    coors = mesh.coors
    if coors.shape[1] == 2:
        coors = nm.hstack([coors, coors[:, [0]] * 0])

    p = nm.asarray(p)[None, :]
    v = nm.asarray(v)

    v = v / nm.linalg.norm(v)

    vc = coors - p
    dist = nm.dot(vc, v)

    new_coors = nm.vstack([
        coors,
        coors - 2 * dist[:, None] * v
    ])

    new_vgroups = nm.hstack([mesh.cmesh.vertex_groups] * 2)
    new_conns, new_mat_ids, new_descs = [], [], []
    i1 = coors.shape[0]

    for desc in mesh.descs:
        conns, cidxs = mesh.get_conn(desc, ret_cells=True)
        new_mat_ids.append(nm.hstack([mesh.cmesh.cell_groups[cidxs]] * 2))
        new_descs.append(desc)
        new_conns.append(nm.vstack([conns, conns[:, mirror_map[desc]] + i1]))

    out = Mesh.from_data(mesh.name, new_coors, new_vgroups,
                         new_conns, new_mat_ids, mesh.descs)

    return merge_nodes(out)
