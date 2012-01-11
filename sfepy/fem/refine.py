"""
Basic uniform mesh refinement functions.
"""
import numpy as nm

from sfepy.linalg import cycle
from sfepy.fem import Mesh

def refine_2_3(mesh_in, ed):
    """
    Refines mesh out of triangles by cutting cutting each edge in half
    and making 4 new finer triangles out of one coarser one.
    """
    # Unique edge centres.
    e_coors, e_uid = ed.get_coors()
    e_centres = 0.5 * nm.sum(e_coors, axis=1)

    # New coordinates after the original ones.
    coors = nm.r_[mesh_in.coors, e_centres]

    conns = []
    mat_ids = []
    for ig, conn in enumerate(mesh_in.conns):
        indx = ed.indx[ig]
        n_el  = conn.shape[0]

        e_nodes = ed.uid_i[indx].reshape((n_el, 3)) + mesh_in.n_nod

        c = nm.c_[conn, e_nodes].T

        new_conn = nm.vstack([c[0], c[3], c[5],
                              c[3], c[4], c[5],
                              c[1], c[4], c[3],
                              c[2], c[5], c[4]]).T
        new_conn = new_conn.reshape((4 * n_el, 3))
        conns.append(new_conn)

        new_mat_id = mesh_in.mat_ids[ig].repeat(4)
        mat_ids.append(new_mat_id)

    mesh = Mesh.from_data(mesh_in.name + '_r', coors, None, conns,
                          mat_ids, mesh_in.descs )

    return mesh

def refine_2_4(mesh_in, ed):
    """
    Refines mesh out of quadrilaterals by cutting cutting each edge in
    half and making 4 new finer quadrilaterals out of one coarser one.
    """
    # Unique edge centres.
    e_coors, e_uid = ed.get_coors()
    e_centres = 0.5 * nm.sum(e_coors, axis=1)

    # Unique element centres.
    coors = mesh_in.get_element_coors()
    centres = 0.25 * nm.sum(coors, axis=1)

    # New coordinates after the original ones.
    coors = nm.r_[mesh_in.coors, e_centres, centres]

    o1 = mesh_in.n_nod
    o2 = o1 + e_centres.shape[0]

    conns = []
    mat_ids = []
    for ig, conn in enumerate(mesh_in.conns):
        e_indx = ed.indx[ig]
        off = mesh_in.el_offsets[ig]
        n_el  = conn.shape[0]

        e_nodes = ed.uid_i[e_indx].reshape((n_el, 4)) + o1
        nodes = nm.arange(n_el) + off + o2

        c = nm.c_[conn, e_nodes, nodes].T

        new_conn = nm.vstack([c[0], c[4], c[8], c[7],
                              c[1], c[5], c[8], c[4],
                              c[2], c[6], c[8], c[5],
                              c[3], c[7], c[8], c[6]]).T
        new_conn = new_conn.reshape((4 * n_el, 4))
        conns.append(new_conn)

        new_mat_id = mesh_in.mat_ids[ig].repeat(4)
        mat_ids.append(new_mat_id)

    mesh = Mesh.from_data(mesh_in.name + '_r', coors, None, conns,
                          mat_ids, mesh_in.descs )

    return mesh

def refine_3_4(mesh_in, ed):
    """
    Refines tetrahedra by cutting each edge in half and making 8 new
    finer tetrahedra out of one coarser one. Old nodal coordinates come
    first in `coors`, then the new ones. The new tetrahedra are similar
    to the old one, no degeneration is supposed to occur as at most 3
    congruence classes of tetrahedra appear, even when re-applied
    iteratively (provided that `conns` are not modified between two
    applications - ordering of vertices in tetrahedra matters not only
    for positivity of volumes).

    References:

    - Juergen Bey: Simplicial grid refinement: on Freudenthal s algorithm and 
      the optimal number of congruence classes, Numer.Math. 85 (2000), 
      no. 1, 1--29, or
    - Juergen Bey: Tetrahedral grid refinement, Computing 55 (1995), 
      no. 4, 355--378, or
      http://citeseer.ist.psu.edu/bey95tetrahedral.html
    """
    # Unique edge centres.
    e_coors, e_uid = ed.get_coors()
    e_centres = 0.5 * nm.sum(e_coors, axis=1)

    # New coordinates after the original ones.
    coors = nm.r_[mesh_in.coors, e_centres]

    conns = []
    mat_ids = []
    for ig, conn in enumerate(mesh_in.conns):
        indx = ed.indx[ig]
        n_el  = conn.shape[0]

        e_nodes = ed.uid_i[indx].reshape((n_el, 6)) + mesh_in.n_nod

        c = nm.c_[conn, e_nodes].T

        new_conn = nm.vstack([c[0], c[4], c[6], c[7],
                              c[4], c[1], c[5], c[8],
                              c[6], c[5], c[2], c[9],
                              c[7], c[8], c[9], c[3],
                              c[4], c[6], c[7], c[8],
                              c[4], c[6], c[8], c[5],
                              c[6], c[7], c[8], c[9],
                              c[6], c[5], c[9], c[8]]).T
        new_conn = new_conn.reshape((8 * n_el, 4))
        conns.append(new_conn)

        new_mat_id = mesh_in.mat_ids[ig].repeat(8)
        mat_ids.append(new_mat_id)

    mesh = Mesh.from_data(mesh_in.name + '_r', coors, None, conns,
                          mat_ids, mesh_in.descs )

    return mesh

def refine_3_8(mesh_in, ed, fa):
    """
    Refines hexahedral mesh by cutting cutting each edge in half and
    making 8 new finer hexahedrons out of one coarser one.
    """
    # Unique edge centres.
    e_coors, e_uid = ed.get_coors()
    e_centres = 0.5 * nm.sum(e_coors, axis=1)

    # Unique face centres.
    f_coors, f_uid = fa.get_coors()
    f_centres = 0.25 * nm.sum(f_coors, axis=1)

    # Unique element centres.
    coors = mesh_in.get_element_coors()
    centres = 0.125 * nm.sum(coors, axis=1)

    # New coordinates after the original ones.
    coors = nm.r_[mesh_in.coors, e_centres, f_centres, centres]

    o1 = mesh_in.n_nod
    o2 = o1 + e_centres.shape[0]
    o3 = o2 + f_centres.shape[0]

    st = nm.vstack

    conns = []
    mat_ids = []
    for ig, conn in enumerate(mesh_in.conns):
        e_indx = ed.indx[ig]
        f_indx = fa.indx[ig]
        off = mesh_in.el_offsets[ig]
        n_el  = conn.shape[0]

        e_nodes = ed.uid_i[e_indx].reshape((n_el, 12)) + o1
        f_nodes = fa.uid_i[f_indx].reshape((n_el, 6)) + o2
        nodes = nm.arange(n_el) + off + o3

        c = nm.c_[conn, e_nodes, f_nodes, nodes].T

        new_conn = st([c[0], c[8], c[20], c[11], c[16], c[22], c[26], c[21],
                       c[1], c[9], c[20], c[8], c[17], c[24], c[26], c[22],
                       c[2], c[10], c[20], c[9], c[18], c[25], c[26], c[24],
                       c[3], c[11], c[20], c[10], c[19], c[21], c[26], c[25],
                       c[4], c[15], c[23], c[12], c[16], c[21], c[26], c[22],
                       c[5], c[12], c[23], c[13], c[17], c[22], c[26], c[24],
                       c[6], c[13], c[23], c[14], c[18], c[24], c[26], c[25],
                       c[7], c[14], c[23], c[15], c[19], c[25], c[26], c[21]]).T
        new_conn = new_conn.reshape((8 * n_el, 8))
        conns.append(new_conn)

        new_mat_id = mesh_in.mat_ids[ig].repeat(8)
        mat_ids.append(new_mat_id)

    mesh = Mesh.from_data(mesh_in.name + '_r', coors, None, conns,
                          mat_ids, mesh_in.descs )

    return mesh

def refine_reference(geometry, level):
    """
    Refine reference element given by `geometry`.

    Notes
    -----
    The error edges must be generated in the order of the connectivity
    of the previous (lower) level.
    """
    gcoors, gconn = geometry.coors, geometry.conn
    if level == 0:
        return gcoors, gconn

    c1d = gcoors[geometry.edges[0], 0]
    n1d = 2**level + 1

    ip = nm.linspace(c1d[0], c1d[1], n1d)

    if geometry.name == '2_3':
        n_edge = 3

        coors = nm.zeros((n1d * (n1d + 1) / 2, 2), dtype=nm.float64)
        g = nm.zeros((n1d, n1d), dtype=nm.int32)
        ii = 0
        for y in range(n1d):
            for x in range(n1d - y):
                g[x, y] = ii
                coors[ii] = ip[[x, y]]
                ii += 1

        conn = []
        for y in range(n1d - 1):
            for x in range(n1d - y - 1):
                conn.append([g[x, y], g[x+1, y], g[x, y+1]])
            for x in range(n1d - y - 2):
                conn.append([g[x+1, y], g[x+1, y+1], g[x, y+1]])

        error_edges = []
        ap = error_edges.append
        for y0 in range(0, n1d - 1, 2):
            y1 = y0 + 1
            y2 = y0 + 2
            for x0 in range(0, n1d - y0 - 1, 2):
                x1 = x0 + 1
                x2 = x0 + 2
                ap([g[x0, y0], g[x1, y0], g[x2, y0]])
                ap([g[x0, y0], g[x0, y1], g[x0, y2]])
                ap([g[x2, y0], g[x1, y1], g[x0, y2]])

            for x0 in range(0, n1d - y0 - 3, 2):
                x1 = x0 + 1
                x2 = x0 + 2
                ap([g[x2, y0], g[x2, y1], g[x2, y2]])
                ap([g[x2, y0], g[x1, y1], g[x0, y2]])
                ap([g[x0, y2], g[x1, y2], g[x2, y2]])

    elif geometry.name == '3_4':
        n_edge = 6

        coors = nm.zeros(((n1d * (n1d + 1) * (n1d + 2)) / 6, 3),
                         dtype=nm.float64)
        g = nm.zeros((n1d, n1d, n1d), dtype=nm.int32)
        ii = 0
        for z in range(n1d):
            for y in range(n1d - z):
                for x in range(n1d - y - z):
                    g[x, y, z] = ii
                    coors[ii] = ip[[x, y, z]]
                    ii += 1

        conn = []
        ap = conn.append
        for z in range(n1d - 1):
            for y in range(n1d - z - 1):
                for x in range(n1d - y - z - 1):
                    # 'abde'
                    ap([g[x, y, z], g[x+1, y, z],
                        g[x, y+1, z], g[x, y, z+1]])

                for x in range(n1d - y - z - 2):
                    # 'bfde'
                    ap([g[x+1, y, z], g[x+1, y, z+1],
                        g[x, y+1, z], g[x, y, z+1]])
                    # 'dfhe'
                    ap([g[x, y+1, z], g[x+1, y, z+1],
                        g[x, y+1, z+1], g[x, y, z+1]])
                    # 'bcdf'
                    ap([g[x+1, y, z], g[x+1, y+1, z],
                        g[x, y+1, z], g[x+1, y, z+1]])
                    # 'dchf'
                    ap([g[x, y+1, z], g[x+1, y+1, z],
                        g[x, y+1, z+1], g[x+1, y, z+1]])

                for x in range(n1d - y - z - 3):
                    # 'fgch'
                    ap([g[x+1, y, z+1], g[x+1, y+1, z+1],
                        g[x+1, y+1, z], g[x, y+1, z+1]])

        error_edges = []
        ap = error_edges.append
        for z0 in range(0, n1d - 1, 2):
            z1 = z0 + 1
            z2 = z0 + 2
            for y0 in range(0, n1d - z0 - 1, 2):
                y1 = y0 + 1
                y2 = y0 + 2
                for x0 in range(0, n1d - y0 - z0 - 1, 2):
                    x1 = x0 + 1
                    x2 = x0 + 2
                    # 'abde'
                    ap([g[x0, y0, z0], g[x1, y0, z0], g[x2, y0, z0]])
                    ap([g[x0, y0, z0], g[x0, y1, z0], g[x0, y2, z0]])
                    ap([g[x2, y0, z0], g[x1, y1, z0], g[x0, y2, z0]])
                    ap([g[x0, y0, z0], g[x0, y0, z1], g[x0, y0, z2]])
                    ap([g[x2, y0, z0], g[x1, y0, z1], g[x0, y0, z2]])
                    ap([g[x0, y2, z0], g[x0, y1, z1], g[x0, y0, z2]])
                for x0 in range(0, n1d - y0 - z0 - 3, 2):
                    x1 = x0 + 1
                    x2 = x0 + 2
                    # 'bfde'
                    ap([g[x2, y0, z0], g[x2, y0, z1], g[x2, y0, z2]])
                    ap([g[x2, y0, z0], g[x1, y1, z0], g[x0, y2, z0]])
                    ap([g[x2, y0, z0], g[x1, y0, z1], g[x0, y0, z2]])
                    ap([g[x0, y2, z0], g[x1, y1, z1], g[x2, y0, z2]])
                    ap([g[x0, y2, z0], g[x0, y1, z1], g[x0, y0, z2]])
                    ap([g[x0, y0, z2], g[x1, y0, z2], g[x2, y0, z2]])
                    # 'dfhe'
                    ap([g[x0, y2, z0], g[x1, y1, z1], g[x2, y0, z2]])
                    ap([g[x0, y2, z0], g[x0, y1, z1], g[x0, y0, z2]])
                    ap([g[x0, y2, z0], g[x0, y2, z1], g[x0, y2, z2]])
                    ap([g[x0, y0, z2], g[x1, y0, z2], g[x2, y0, z2]])
                    ap([g[x0, y0, z2], g[x0, y1, z2], g[x0, y2, z2]])
                    ap([g[x0, y2, z2], g[x1, y1, z2], g[x2, y0, z2]])
                    # 'bcdf'
                    ap([g[x2, y0, z0], g[x2, y1, z0], g[x2, y2, z0]])
                    ap([g[x2, y0, z0], g[x2, y0, z1], g[x2, y0, z2]])
                    ap([g[x2, y0, z0], g[x1, y1, z0], g[x0, y2, z0]])
                    ap([g[x0, y2, z0], g[x1, y2, z0], g[x2, y2, z0]])
                    ap([g[x0, y2, z0], g[x1, y1, z1], g[x2, y0, z2]])
                    ap([g[x2, y0, z2], g[x2, y1, z1], g[x2, y2, z0]])
                    # 'dchf'
                    ap([g[x0, y2, z0], g[x1, y2, z0], g[x2, y2, z0]])
                    ap([g[x0, y2, z0], g[x1, y1, z1], g[x2, y0, z2]])
                    ap([g[x0, y2, z0], g[x0, y2, z1], g[x0, y2, z2]])
                    ap([g[x0, y2, z2], g[x1, y1, z2], g[x2, y0, z2]])
                    ap([g[x0, y2, z2], g[x1, y2, z1], g[x2, y2, z0]])
                    ap([g[x2, y2, z0], g[x2, y1, z1], g[x2, y0, z2]])

                for x0 in range(0, n1d - y0 - z0 - 5, 2):
                    x1 = x0 + 1
                    x2 = x0 + 2
                    # 'fgch'
                    ap([g[x2, y2, z2], g[x1, y2, z2], g[x0, y2, z2]])
                    ap([g[x2, y2, z2], g[x2, y1, z2], g[x2, y0, z2]])
                    ap([g[x2, y2, z2], g[x2, y2, z1], g[x2, y2, z0]])
                    ap([g[x2, y0, z2], g[x2, y1, z1], g[x2, y2, z0]])
                    ap([g[x2, y0, z2], g[x1, y1, z2], g[x0, y2, z2]])
                    ap([g[x0, y2, z2], g[x1, y2, z1], g[x2, y2, z0]])

    elif geometry.name == '2_4':
        n_edge = 6

        shape = (n1d, n1d)
        coors = nm.zeros((nm.prod(shape), 2), dtype=nm.float64)
        g = nm.zeros(shape, dtype=nm.int32)
        for ii, (y, x) in enumerate(cycle(shape)):
            g[x, y] = ii
            coors[ii] = ip[[x, y]]

        conn = []
        gc = nm.zeros(nm.array(shape) - 1, dtype=nm.int32)
        for ii, (y, x) in enumerate(cycle(nm.array(shape) - 1)):
            gc[x, y] = ii
            conn.append([g[x  , y  ], g[x+1, y  ],
                         g[x+1, y+1], g[x  , y+1]])

        nesting = []
        apn = nesting.append
        error_edges = []
        ap = error_edges.append
        for y0 in range(0, n1d - 1, 2):
            y1 = y0 + 1
            y2 = y0 + 2
            for x0 in range(0, n1d - 1, 2):
                x1 = x0 + 1
                x2 = x0 + 2
                ap([g[x0, y0], g[x1, y0], g[x2, y0]])
                ap([g[x0, y1], g[x1, y1], g[x2, y1]])
                ap([g[x0, y2], g[x1, y2], g[x2, y2]])
                ap([g[x0, y0], g[x0, y1], g[x0, y2]])
                ap([g[x1, y0], g[x1, y1], g[x1, y2]])
                ap([g[x2, y0], g[x2, y1], g[x2, y2]])

                apn([gc[x0, y0], gc[x1, y0], gc[x0, y1], gc[x1, y1]])

    elif geometry.name == '3_8':
        n_edge = 27

        shape = (n1d, n1d, n1d)
        coors = nm.zeros((nm.prod(shape), 3), dtype=nm.float64)
        g = nm.zeros(shape, dtype=nm.int32)
        for ii, (z, y, x) in enumerate(cycle(shape)):
            g[x, y, z] = ii
            coors[ii] = ip[[x, y, z]]

        conn = []
        gc = nm.zeros(nm.array(shape) - 1, dtype=nm.int32)
        for ii, (z, y, x) in enumerate(cycle(nm.array(shape) - 1)):
            gc[x, y, z] = ii
            conn.append([g[x  , y  , z  ], g[x+1, y  , z  ],
                         g[x+1, y+1, z  ], g[x  , y+1, z  ],
                         g[x  , y  , z+1], g[x+1, y  , z+1],
                         g[x+1, y+1, z+1], g[x  , y+1, z+1]])

        nesting = []
        apn = nesting.append
        error_edges = []
        ap = error_edges.append
        for z0 in range(0, n1d - 1, 2):
            z1 = z0 + 1
            z2 = z0 + 2
            for y0 in range(0, n1d - 1, 2):
                y1 = y0 + 1
                y2 = y0 + 2
                for x0 in range(0, n1d - 1, 2):
                    x1 = x0 + 1
                    x2 = x0 + 2
                    ap([g[x0, y0, z0], g[x1, y0, z0], g[x2, y0, z0]])
                    ap([g[x0, y1, z0], g[x1, y1, z0], g[x2, y1, z0]])
                    ap([g[x0, y2, z0], g[x1, y2, z0], g[x2, y2, z0]])
                    ap([g[x0, y0, z0], g[x0, y1, z0], g[x0, y2, z0]])
                    ap([g[x1, y0, z0], g[x1, y1, z0], g[x1, y2, z0]])
                    ap([g[x2, y0, z0], g[x2, y1, z0], g[x2, y2, z0]])

                    ap([g[x0, y0, z1], g[x1, y0, z1], g[x2, y0, z1]])
                    ap([g[x0, y1, z1], g[x1, y1, z1], g[x2, y1, z1]])
                    ap([g[x0, y2, z1], g[x1, y2, z1], g[x2, y2, z1]])
                    ap([g[x0, y0, z1], g[x0, y1, z1], g[x0, y2, z1]])
                    ap([g[x1, y0, z1], g[x1, y1, z1], g[x1, y2, z1]])
                    ap([g[x2, y0, z1], g[x2, y1, z1], g[x2, y2, z1]])

                    ap([g[x0, y0, z2], g[x1, y0, z2], g[x2, y0, z2]])
                    ap([g[x0, y1, z2], g[x1, y1, z2], g[x2, y1, z2]])
                    ap([g[x0, y2, z2], g[x1, y2, z2], g[x2, y2, z2]])
                    ap([g[x0, y0, z2], g[x0, y1, z2], g[x0, y2, z2]])
                    ap([g[x1, y0, z2], g[x1, y1, z2], g[x1, y2, z2]])
                    ap([g[x2, y0, z2], g[x2, y1, z2], g[x2, y2, z2]])

                    ap([g[x0, y0, z0], g[x0, y0, z1], g[x0, y0, z2]])
                    ap([g[x1, y0, z0], g[x1, y0, z1], g[x1, y0, z2]])
                    ap([g[x2, y0, z0], g[x2, y0, z1], g[x2, y0, z2]])

                    ap([g[x0, y1, z0], g[x0, y1, z1], g[x0, y1, z2]])
                    ap([g[x1, y1, z0], g[x1, y1, z1], g[x1, y1, z2]])
                    ap([g[x2, y1, z0], g[x2, y1, z1], g[x2, y1, z2]])

                    ap([g[x0, y2, z0], g[x0, y2, z1], g[x0, y2, z2]])
                    ap([g[x1, y2, z0], g[x1, y2, z1], g[x1, y2, z2]])
                    ap([g[x2, y2, z0], g[x2, y2, z1], g[x2, y2, z2]])

                    apn([gc[x0, y0, z0], gc[x1, y0, z0],
                         gc[x0, y1, z0], gc[x1, y1, z0],
                         gc[x0, y0, z1], gc[x1, y0, z1],
                         gc[x0, y1, z1], gc[x1, y1, z1]])

    else:
        raise ValueError('unsupported geometry! (%s)' % geometry.name)

    conn = nm.array(conn, dtype=nm.int32)
    nesting = nm.array(nesting, dtype=nm.int32)
    error_edges = nm.array(error_edges, dtype=nm.int32)

    sh = error_edges.shape
    error_edges.shape = (sh[0] / n_edge, n_edge, sh[1])

    return coors, conn, nesting, error_edges
