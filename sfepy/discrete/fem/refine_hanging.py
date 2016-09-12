"""
Functions for a mesh refinement with hanging nodes.

Notes
-----
Using LCBCs with hanging nodes is not supported.
"""
from six.moves import range

import numpy as nm

from sfepy.base.base import assert_
from sfepy.discrete import Integral, Functions, Function
from sfepy.discrete.fem import Mesh, FEDomain, Field

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
      - cells[:, 1] - their facet sharing neighbors
      - cells[:, 2] - facet kind: 0 = face, 1 = edge
      - cells[:, 3] - their facet sharing neighbors w.r.t. the locally refined
        mesh.
    """
    if not refine_flag.any():
        facets = nm.zeros((0, 3), dtype=nm.uint32)
        cells = nm.zeros((0, 4), dtype=nm.uint32)
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

        # Indices of non-refined cells in the level-1 mesh.
        ii = nm.searchsorted(region0.cells, cells[:, 1])
        cells = nm.c_[cells, nm.zeros_like(ii), ii]

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

        emask = nm.zeros((domain.shape.n_el, 12), dtype=nm.bool)

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

        valid = nm.ones(nnn.shape[0], dtype=nm.bool)
        valid[omit] = False

        cells = nnn[valid]
        facets = ffs[valid]

        # Indices of non-refined cells in the level-1 mesh.
        ii = nm.searchsorted(region0.cells, cells[:, 1])
        cells = nm.c_[cells, ii]

    return facets, cells, oe, region0, region1

def refine_region(domain0, region0, region1):
    """
    Coarse cell sub_cells[ii, 0] in mesh0 is split into sub_cells[ii, 1:] in
    mesh1.
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
    aux = nm.arange(n_sub * n_cell, dtype=nm.uint32).reshape((-1, n_sub))
    sub_cells[:, 1:] = region0.shape.n_cell + aux

    coors0, vgs0, conns0, mat_ids0, descs0 = mesh0._get_io_data()

    coors, vgs, _conns, _mat_ids, descs = mesh1r._get_io_data()

    conn0 = mesh0.get_conn(domain0.mesh.descs[0])
    conns = [nm.r_[conn0[region0.cells], _conns[0]]]
    mat_ids = [nm.r_[mat_ids0[0][region0.cells], _mat_ids[0]]]
    mesh = Mesh.from_data('a', coors, vgs, conns, mat_ids, descs)
    domain = FEDomain('d', mesh)

    # Preserve orientations of coarse cell edges and faces.
    if domain0.cmesh.edge_oris is not None:
        n_cell = domain0.shape.n_el
        oris0 = domain0.cmesh.edge_oris.reshape((n_cell, -1))[region0.cells]

        n_cell = domain.shape.n_el
        domain.cmesh.edge_oris.reshape((n_cell, -1))[:oris0.shape[0]] = oris0

    if domain0.cmesh.face_oris is not None:
        n_cell = domain0.shape.n_el
        oris0 = domain0.cmesh.face_oris.reshape((n_cell, -1))[region0.cells]

        n_cell = domain.shape.n_el
        domain.cmesh.face_oris.reshape((n_cell, -1))[:oris0.shape[0]] = oris0

    return domain, sub_cells

def find_facet_substitutions(facets, cells, sub_cells, oris, refine_facets):
    """
    Find facet substitutions in connectivity as well as coarse cell facets
    orientations.

    sub = [coarse cell, coarse facet, fine1 cell, fine1 facet, fine2 cell,
           fine2 facet]
    """
    subs = []
    coris = nm.zeros(facets.shape[0], dtype=nm.uint32)
    for ii, fac in enumerate(facets):
        fine = cells[ii, 0]
        coarse = cells[ii, 3]

        isub = nm.searchsorted(sub_cells[:, 0], fine)

        coris[ii] = oris[isub, fac[1]]

        refined = sub_cells[isub, 1:]
        rf = refine_facets[fac[1]]
        used = refined[rf[:, 0]]
        fused = rf[:, 1]

        master = [coarse, fac[2]]
        slave = zip(used, fused)
        sub = nm.r_[[master], slave].ravel()

        # !!!!!
        print ii, fac, fine, coarse, isub, refined, used
        subs.append(sub)

    subs = nm.array(subs)
    return subs, coris

def force_facet_orientations(domain, gsubs, coris, dim):
    """
    Force coarse cell facets orientations to fine cell facets.
    """
    if gsubs is None: return

    cmesh = domain.mesh.cmesh
    foris = cmesh.get_orientations(dim)
    num = foris.shape[0] / cmesh.n_el

    for ii, subs in enumerate(gsubs):
        ic = subs[2::2] * num + subs[3::2]

        foris[ic] = coris[ii]

def refine(domain0, refine, gsubs=None):
    desc = domain0.mesh.descs[0]
    assert_(desc in ['2_4', '3_8'])

    facets, cells, oe, region0, region1 = find_level_interface(domain0, refine)

    print nm.c_[facets, cells]

    domain, sub_cells = refine_region(domain0, region0, region1)

    #_plot(domain.mesh.cmesh)

    if facets.shape[0] > 0:
        desc = domain0.mesh.descs[0]
        conn0 = domain0.mesh.get_conn(desc)
        conn1 = domain.mesh.get_conn(desc)

        print conn1

        print conn1[sub_cells[:, 1:]]
        print conn0[sub_cells[:, 0]]

        print cells[:, 2]
        print conn0[cells[:, 1]]
        print conn1[cells[:, 3]]
        assert_((conn0[cells[:, 1]] == conn1[cells[:, 3]]).all())

    desc = domain0.mesh.descs[0]
    if desc == '2_4':
        cmesh0 = domain0.mesh.cmesh
        oris = cmesh0.edge_oris.reshape((-1, 4))[region1.cells]

        gsubs1, coris = find_facet_substitutions(facets, cells,
                                                 sub_cells, oris,
                                                 refine_edges_2_4)
        force_facet_orientations(domain, gsubs1, coris, 1)

        if gsubs is None:
            gsubs = gsubs1 if len(gsubs1) else None

        elif len(gsubs1):
            mods = nm.zeros(domain.shape.n_el + 1, dtype=nm.int32)
            mods[refine > 0] = -1
            mods = nm.cumsum(mods)

            gsubs[:, [0, 2, 4]] += mods[gsubs[:, [0, 2, 4]]]

            gsubs = nm.r_[gsubs, gsubs1]

    else:
        cmesh0 = domain0.mesh.cmesh
        foris = cmesh0.face_oris.reshape((-1, 6))[region1.cells]
        eoris = cmesh0.edge_oris.reshape((-1, 12))[region1.cells]

        gsubs1f, fcoris = find_facet_substitutions(facets[:oe], cells[:oe],
                                                   sub_cells, foris,
                                                   refine_faces_3_8)
        force_facet_orientations(domain, gsubs1f, fcoris, 2)

        gsubs1e, ecoris = find_facet_substitutions(facets[oe:], cells[oe:],
                                                   sub_cells, eoris,
                                                   refine_edges_3_8)
        force_facet_orientations(domain, gsubs1e, ecoris, 1)

        if gsubs is None:
            gsubs = (gsubs1f if len(gsubs1f) else None,
                     gsubs1e if len(gsubs1f) else None) # !!!

        elif len(gsubs1):
            gsubsf, gsubse = gsubs

            mods = nm.zeros(domain.shape.n_el + 1, dtype=nm.int32)
            mods[refine > 0] = -1
            mods = nm.cumsum(mods)

            gsubsf[:, [0, 2, 4, 6, 8]] += mods[gsubsf[:, [0, 2, 4, 6, 8]]]
            gsubsf = nm.r_[gsubsf, gsubs1f]

            gsubse[:, [0, 2, 4]] += mods[gsubse[:, [0, 2, 4]]]
            gsubse = nm.r_[gsubse, gsubs1e]

            gsubs = (gsubs1f, gsubs1e)

        if gsubs == (None, None): gsubs = None

    return domain, gsubs

def eval_basis_transform(field, gsubs):
    """
    """
    transform = nm.tile(nm.eye(field.econn.shape[1]),
                        (field.econn.shape[0], 1, 1))
    if gsubs is None:
        return transform

    gel = field.gel
    ao = field.approx_order

    conn = [gel.conn]
    mesh = Mesh.from_data('a', gel.coors, None, [conn], [nm.array([0])],
                          [gel.name])
    cdomain = FEDomain('d', mesh)
    comega = cdomain.create_region('Omega', 'all')
    rcfield = Field.from_args('rc', field.dtype, 1, comega, approx_order=ao)

    fdomain = cdomain.refine()
    fomega = fdomain.create_region('Omega', 'all')
    rffield = Field.from_args('rf', field.dtype, 1, fomega, approx_order=ao)

    def assign_transform(transform, bf, gsubs, ef):
        if not len(gsubs): return

        n_sub = (gsubs.shape[1] - 2) // 2

        for ii, sub in enumerate(gsubs):
            print ii, sub
            for ij in range(n_sub):
                ik = 2 * (ij + 1)
                print ij, ik

                fface = ef[sub[ik+1]]

                mtx = transform[sub[ik]]
                print sub[ik], mtx
                ix, iy = nm.meshgrid(fface, fface)

                cbf = bf[iy, 0, ix]

                mtx[ix, iy] = cbf
                print mtx

    fcoors = rffield.get_coor()

    coors = fcoors[rffield.econn[0]]
    integral = Integral('i', coors=coors, weights=nm.ones_like(coors[:, 0]))

    rcfield.clear_qp_base()
    bf = rcfield.get_base('v', False, integral)

    if gel.name == '2_4':
        fgsubs = gsubs
        egsubs = None

        assign_transform(transform, bf, fgsubs, rffield.efaces)

    else:
        fgsubs = gsubs[0]
        egsubs = gsubs[1]

        assign_transform(transform, bf, fgsubs, rffield.efaces)
        if egsubs is not None:
            assign_transform(transform, bf, egsubs, rffield.eedges)

    assert_((nm.abs(transform.sum(1) - 1.0) < 1e-15).all())

    return transform
