import numpy as nm

from sfepy.fem.geometry_element import geometry_data
from sfepy.fem import Mesh

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
