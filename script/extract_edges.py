#!/usr/bin/env python
"""
Extract outline edges of a given mesh and save them into
'<original path>/edge_<original mesh file name>.vtk'
or into a user defined output file.
The outline edge is an edge for which norm(nvec1 - nvec2) < eps,
where nvec1 and nvec2 are the normal vectors of the incident facets.
"""
from __future__ import absolute_import
import numpy as nm
from scipy.sparse import coo_matrix
import sys
sys.path.append('.')
from argparse import ArgumentParser
from sfepy.base.base import output, Struct
from sfepy.base.ioutils import edit_filename
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.discrete.fem.meshio import MeshioLibIO


def merge_lines(mesh, eps=1e-18):
    coors, ngroups, conns, mat_ids, ctype = mesh
    conns = conns[0]

    # vertices to edges map     
    n_v = coors.shape[0]
    n_e = conns.shape[0]
    row = nm.repeat(nm.arange(n_e), 2)
    aux = coo_matrix((nm.ones((n_e * 2,), dtype=nm.bool),
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
        sf2e = coo_matrix((nm.ones((n_sf * n_ef,), dtype=nm.bool),
                           (row, se_map0)), shape=(n_sf, n_se))
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

helps = {
    'eps': 'tolerance parameter of the edge search algorithm (default: 1e-12)',
    'filename-out': 'name of output file',
}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.add_argument('--eps', action='store', dest='eps',
                        default=1e-12, help=helps['eps'])
    parser.add_argument('-o', '--filename-out',
                        action='store', dest='filename_out',
                        default=None, help=helps['filename-out'])
    parser.add_argument('filename')
    options = parser.parse_args()

    filename = options.filename

    mesh = Mesh.from_file(filename)
    mesh_out = extract_edges(mesh, eps=float(options.eps))
    mesh_out = merge_lines(mesh_out)

    filename_out = options.filename_out
    if filename_out is None:
        filename_out = edit_filename(filename, prefix='edge_', new_ext='.vtk')

    output('Outline mesh - vertices: %d, edges: %d, output filename: %s'
           % (mesh_out[0].shape[0], mesh_out[2][0].shape[0], filename_out))

    # hack to write '3_2' elements - edges
    io = MeshioLibIO()
    aux_mesh = Struct()
    aux_mesh._get_io_data = lambda: mesh_out
    aux_mesh.n_el = mesh_out[2][0].shape[0]
    io.write(filename_out, aux_mesh)

if __name__ == '__main__':
    main()
