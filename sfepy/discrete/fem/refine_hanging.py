"""
Functions for a mesh refinement with hanging nodes.

Notes
-----
Using LCBCs with hanging nodes is not supported.
"""

import numpy as nm

from sfepy.base.base import assert_
from sfepy.discrete import Functions, Function
from sfepy.discrete.fem import Mesh, FEDomain

# Rows = facets of reference cell, columns = [sub_cell_i, local facet_i]
refine_edges_2_4 = nm.array([[[0, 0], [1, 3]],
                             [[1, 0], [2, 3]],
                             [[2, 0], [3, 3]],
                             [[3, 0], [0, 3]]])

refine_faces_3_8 = nm.array([[[0, 0], [1, 0], [2, 0], [3, 0]],
                             [[0, 1], [3, 2], [4, 2], [7, 1]],
                             [[0, 2], [1, 1], [4, 1], [5, 2]],
                             [[4, 0], [5, 0], [6, 0], [7, 0]],
                             [[1, 2], [2, 1], [5, 1], [6, 2]],
                             [[2, 2], [3, 1], [6, 1], [7, 2]]])

refine_edges_3_8 = nm.array([[[0, 0], [1, 3]],
                             [[1, 0], [2, 3]],
                             [[2, 0], [3, 3]],
                             [[3, 0], [0, 3]],
                             [[4, 3], [5, 0]],
                             [[5, 3], [6, 0]],
                             [[6, 3], [7, 0]],
                             [[7, 3], [4, 0]],
                             [[0, 8], [4, 8]],
                             [[1, 8], [5, 8]],
                             [[2, 8], [6, 8]],
                             [[3, 8], [7, 8]]])

def find_level_interface(domain, refine_flag):
    """
    Find facets of the coarse mesh that are on the coarse-refined cell
    boundary.

    ids w.r.t. current mesh:
    - facets: global, local w.r.t. cells[:, 0], local w.r.t. cells[:, 1]

    - interface cells:
      - cells[:, 0] - cells to refine
      - cells[:, 1] - their facet sharing neighbors (w.r.t. both meshes)
      - cells[:, 2] - facet kind: 0 = face, 1 = edge
    """
    if not refine_flag.any():
        facets = nm.zeros((0, 3), dtype=nm.uint32)
        cells = nm.zeros((0, 3), dtype=nm.uint32)
        return facets, cells, 0, None, None

    def _get_refine(coors, domain=None):
        return nm.nonzero(refine_flag)[0]

    def _get_coarse(coors, domain=None):
        return nm.nonzero(1 - refine_flag)[0]

    get_refine = Function('get_refine', _get_refine)
    get_coarse = Function('get_coarse', _get_coarse)
    functions = Functions([get_refine, get_coarse])
    region0 = domain.create_region('coarse', 'cells by get_coarse',
                                   functions=functions, add_to_regions=False,
                                   allow_empty=True)
    region1 = domain.create_region('refine', 'cells by get_refine',
                                   functions=functions, add_to_regions=False)

    cmesh = domain.mesh.cmesh
    dim = cmesh.dim

    if dim == 2:
        oe = 0

        facets = nm.intersect1d(region0.facets, region1.facets)

        cmesh.setup_connectivity(dim - 1, dim)
        cells, offs = cmesh.get_incident(dim, facets, dim - 1,
                                         ret_offsets=True)
        assert_((nm.diff(offs) == 2).all())

        ii = cmesh.get_local_ids(facets, dim - 1, cells, offs, dim)

        ii = ii.reshape((-1, 2))
        cells = cells.reshape((-1, 2))

        ii = nm.where(refine_flag[cells], ii[:, :1], ii[:, 1:])
        cells = nm.where(refine_flag[cells], cells[:, :1], cells[:, 1:])

        facets = nm.c_[facets, ii]

        cells = nm.c_[cells, nm.zeros_like(cells[:, 1])]

    else: # if dim == 3:
        gel = domain.geom_els['3_8']
        epf = gel.get_edges_per_face()

        cmesh.setup_connectivity(dim, dim)
        fc, coffs = cmesh.get_incident(dim, region1.cells, dim,
                                       ret_offsets=True)
        cc = nm.repeat(region1.cells, nm.diff(coffs))
        aux = nm.c_[cc, fc]
        """
        nnn[:, 0] cells to refine, nnn[:, 1] non-refined neighbours, nnn[:, 2]
        neighbour kind : 0 face, 1 edge.
        """
        nn = aux[refine_flag[fc] == 0]

        cf = nn[:, 0].copy().astype(nm.uint32)
        cc = nn[:, 1].copy().astype(nm.uint32)

        vc, vco = cmesh.get_incident(0, cc, dim, ret_offsets=True)
        vf, vfo = cmesh.get_incident(0, cf, dim, ret_offsets=True)
        vc = vc.reshape((-1, 8))
        vf = vf.reshape((-1, 8))

        nnn = []
        oe = 0
        ov = nn.shape[0]
        for ii in range(vc.shape[0]):
            aux = set(vc[ii]).intersection(vf[ii])
            nc = len(aux)
            if nc == 1:
                nnn.append((0, 0, 2))
                ov -= 1

            elif nc == 4:
                nnn.append((nn[ii, 0], nn[ii, 1], 0))
                oe += 1

            else:
                nnn.append((nn[ii, 0], nn[ii, 1], 1))
        nnn = nm.array(nnn)

        if nnn.shape[0] == 0:
            facets = nm.zeros((0, 3), dtype=nm.uint32)
            cells = nm.zeros((0, 4), dtype=nm.uint32)

            return facets, cells, 0, region0, region1

        # Sort by neighbour kind, skip vertex-only neighbours.
        ii = nm.argsort(nnn[:, 2])
        nnn = nnn[ii][:ov]
        cf = cf[ii][:ov]
        cc = cc[ii][:ov]

        ec, eco = cmesh.get_incident(1, cc, dim, ret_offsets=True)
        ef, efo = cmesh.get_incident(1, cf, dim, ret_offsets=True)
        ec = ec.reshape((-1, 12))
        ef = ef.reshape((-1, 12))

        fc, fco = cmesh.get_incident(2, cc, dim, ret_offsets=True)
        ff, ffo = cmesh.get_incident(2, cf, dim, ret_offsets=True)
        fc = fc.reshape((-1, 6))
        ff = ff.reshape((-1, 6))

        emask = nm.zeros((domain.shape.n_el, 12), dtype=bool)

        ffs = []
        for ii in range(oe):
            facet = nm.intersect1d(fc[ii], ff[ii])[0]
            i1 = nm.where(ff[ii] == facet)[0][0]
            i0 = nm.where(fc[ii] == facet)[0][0]
            ffs.append((facet, i1, i0))

            emask[nnn[ii, 0], epf[i1]] = True

        for ii in range(oe, nnn.shape[0]):
            facet = nm.intersect1d(ec[ii], ef[ii])[0]
            i1 = nm.where(ef[ii] == facet)[0][0]
            i0 = nm.where(ec[ii] == facet)[0][0]
            ffs.append((facet, i1, i0))

        ffs = nm.array(ffs)

        ie = nm.where(nnn[:, 2] == 1)[0]
        ennn = nnn[ie]
        effs = ffs[ie]
        omit = ie[emask[ennn[:, 0], effs[:, 1]]]

        valid = nm.ones(nnn.shape[0], dtype=bool)
        valid[omit] = False

        cells = nnn[valid]
        facets = ffs[valid]

    return facets, cells, oe, region0, region1

def refine_region(domain0, region0, region1):
    """
    Coarse cell sub_cells[ii, 0] in mesh0 is split into sub_cells[ii, 1:] in
    mesh1.

    The new fine cells are interleaved among the original coarse cells so that
    the indices of the coarse cells do not change.

    The cell groups are preserved. The vertex groups are preserved only in the
    coarse (non-refined) cells.
    """
    if region1 is None:
        return domain0, None

    mesh0 = domain0.mesh
    mesh1 = Mesh.from_region(region1, mesh0)
    domain1 = FEDomain('d', mesh1)
    domain1r = domain1.refine()
    mesh1r = domain1r.mesh

    n_cell = region1.shape.n_cell
    n_sub = 4 if mesh0.cmesh.tdim == 2 else 8

    sub_cells = nm.empty((n_cell, n_sub + 1), dtype=nm.uint32)
    sub_cells[:, 0] = region1.cells
    sub_cells[:, 1] = region1.cells
    aux = nm.arange((n_sub - 1) * n_cell, dtype=nm.uint32)
    sub_cells[:, 2:] = mesh0.n_el + aux.reshape((n_cell, -1))

    coors0, vgs0, conns0, mat_ids0, descs0 = mesh0._get_io_data()
    coors, vgs, _conns, _mat_ids, descs = mesh1r._get_io_data()

    # Preserve vertex groups of non-refined cells.
    vgs[:len(vgs0)] = vgs0

    def _interleave_refined(c0, c1):
        if c1.ndim == 1:
            c0 = c0[:, None]
            c1 = c1[:, None]

        n_row, n_col = c1.shape
        n_new = region0.shape.n_cell + n_row

        out = nm.empty((n_new, n_col), dtype=c0.dtype)
        out[region0.cells] = c0[region0.cells]
        out[region1.cells] = c1[::n_sub]
        aux = c1.reshape((-1, n_col * n_sub))
        out[mesh0.n_el:] = aux[:, n_col:].reshape((-1, n_col))

        return out

    conn = _interleave_refined(conns0[0], _conns[0])
    mat_id = _interleave_refined(mat_ids0[0], _mat_ids[0]).squeeze()

    mesh = Mesh.from_data('a', coors, vgs, [conn], [mat_id], descs)
    domain = FEDomain('d', mesh)

    return domain, sub_cells

def find_facet_substitutions(facets, cells, sub_cells, refine_facets):
    """
    Find facet substitutions in connectivity.

    sub = [coarse cell, coarse facet, fine1 cell, fine1 facet, fine2 cell,
           fine2 facet]
    """
    subs = []
    for ii, fac in enumerate(facets):
        fine = cells[ii, 0]
        coarse = cells[ii, 1]

        isub = nm.searchsorted(sub_cells[:, 0], fine)

        refined = sub_cells[isub, 1:]
        rf = refine_facets[fac[1]]
        used = refined[rf[:, 0]]
        fused = rf[:, 1]

        master = [coarse, fac[2]]
        slave = list(zip(used, fused))
        sub = nm.r_[[master], slave].ravel()

        subs.append(sub)

    subs = nm.array(subs)
    return subs

def refine(domain0, refine, subs=None, ret_sub_cells=False):
    desc = domain0.mesh.descs[0]
    assert_(desc in ['2_4', '3_8'])

    facets, cells, oe, region0, region1 = find_level_interface(domain0, refine)

    if region1 is None:
        return domain0, None

    domain, sub_cells = refine_region(domain0, region0, region1)

    if facets.shape[0] > 0:
        desc = domain0.mesh.descs[0]
        conn0 = domain0.mesh.get_conn(desc)
        conn1 = domain.mesh.get_conn(desc)

        assert_((conn0[cells[:, 1]] == conn1[cells[:, 1]]).all())

    desc = domain0.mesh.descs[0]
    if desc == '2_4':
        subs1 = find_facet_substitutions(facets, cells,
                                         sub_cells, refine_edges_2_4)
        if subs is None:
            subs = subs1 if len(subs1) else None

        elif len(subs1):
            subs = nm.r_[subs, subs1]

    else:
        subs1f = find_facet_substitutions(facets[:oe], cells[:oe],
                                          sub_cells, refine_faces_3_8)

        subs1e = find_facet_substitutions(facets[oe:], cells[oe:],
                                          sub_cells, refine_edges_3_8)

        if subs is None:
            subs = (subs1f if len(subs1f) else None,
                    subs1e if len(subs1f) else None) # !!!

        elif len(subs1f):
            subsf, subse = subs

            subsf = nm.r_[subsf, subs1f]

            if len(subse):
                if len(subs1e):
                    subse = nm.r_[subse, subs1e]

            subs = (subsf, subse)

        if (isinstance(subs, tuple)
            and (subs[0] is None) and (subs[1] is None)): subs = None

    out = (domain, subs)
    if ret_sub_cells:
        out += (sub_cells,)

    return out
