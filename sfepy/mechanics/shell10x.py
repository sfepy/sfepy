"""
Functions implementing the shell10x element.
"""
import numpy as nm

from sfepy.linalg import norm_l2_along_axis as norm
from sfepy.linalg import dot_sequences as ddot

def create_elastic_tensor(young, poisson, shear_correction=True):
    """
    Create the elastic tensor with the applied shear correction (the default)
    for the shell10x element.
    """
    from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson

    # Static condensation so that \sigma_{33} = 0.
    mtx = stiffness_from_youngpoisson(3, young, poisson, plane='stress')
    mtx[2, :3] = mtx[:2, 2] = 0.0

    if shear_correction:
        mtx[4, 4] *= 5.0 / 6.0
        mtx[5, 5] *= 5.0 / 6.0

    return mtx

def create_transformation_matrix(coors):
    r"""
    Create a transposed coordinate transformation matrix, that transforms 3D
    coordinates of quadrilateral cell vertices so that the transformed vertices
    of a plane cell are in the :math:`x-y` plane. The rotation is performed
    w.r.t. the centres of quadrilaterals.

    Parameters
    ----------
    coors : array
        The coordinates of cell vertices, shape `(n_el, 4, 3)`.

    Returns
    -------
    mtx_t : array
        The transposed transformation matrix :math:`T`, i.e.
        :math:`X_{inplane} = T^T X_{3D}`.

    Notes
    -----
    :math:`T = [t_1, t_2, n]`, where :math:`t_1`, :math:`t_2`, are unit
    in-plane (column) vectors and :math:`n` is the unit normal vector,
    all mutually orthonormal.
    """
    # Local coordinate system.
    t1 = coors[:, 1, :] + coors[:, 2, :] - coors[:, 0, :] - coors[:, 3, :]
    t2 = coors[:, 2, :] + coors[:, 3, :] - coors[:, 0, :] - coors[:, 1, :]
    n = nm.cross(t1, t2)
    t2 = nm.cross(n, t1)

    t1 = t1 / norm(t1)[:, None]
    t2 = t2 / norm(t2)[:, None]
    n = n / norm(n)[:, None]

    # Coordinate transformation matrix (transposed!).
    mtx_t = nm.concatenate((t1[:, :, None],
                            t2[:, :, None],
                            n[:, :, None]), axis=2)

    return mtx_t

def transform_asm_matrices(out, mtx_t, blocks):
    """
    Transform matrix assembling contributions to global coordinate system, one
    node at a time.

    Parameters
    ----------
    out : array
        The array of matrices, transformed in-place.
    mtx_t : array
        The array of transposed transformation matrices :math:`T`, see
        :func:`create_transformation_matrix`.
    blocks : array
        The DOF blocks that are
    """
    for ir in blocks:
        ir = ir[:, None]
        for ic in blocks:
            fn = out[:, 0, ir, ic]
            out[:, 0, ir, ic] = ddot(ddot(mtx_t, fn, 'AB'), mtx_t, 'ABT')

def create_local_bases(coors):
    """
    Create local orthonormal bases in each vertex of quadrilateral cells.

    Parameters
    ----------
    coors : array
        The coordinates of cell vertices, shape `(n_el, 4, 3)`.

    Returns
    -------
    ebs : array
        The local bases, shape `(n_el, 4, 3, 3)`. The basis vectors are rows of
        the (..., 3, 3) blocks.
    """
    ebs = nm.zeros((coors.shape[0], 4, 3, 3), dtype=nm.float64)
    ebs[:, 0, 0, :] = coors[:, 1, :] - coors[:, 0, :]
    ebs[:, 0, 1, :] = coors[:, 3, :] - coors[:, 0, :]
    ebs[:, 0, 2, :] = nm.cross(ebs[:, 0, 0, :], ebs[:, 0, 1, :])
    ebs[:, 1, 0, :] = coors[:, 1, :] - coors[:, 0, :]
    ebs[:, 1, 1, :] = coors[:, 2, :] - coors[:, 1, :]
    ebs[:, 1, 2, :] = nm.cross(ebs[:, 1, 0, :], ebs[:, 1, 1, :])
    ebs[:, 2, 0, :] = coors[:, 2, :] - coors[:, 3, :]
    ebs[:, 2, 1, :] = coors[:, 2, :] - coors[:, 1, :]
    ebs[:, 2, 2, :] = nm.cross(ebs[:, 2, 0, :], ebs[:, 2, 1, :])
    ebs[:, 3, 0, :] = coors[:, 2, :] - coors[:, 3, :]
    ebs[:, 3, 1, :] = coors[:, 3, :] - coors[:, 0, :]
    ebs[:, 3, 2, :] = nm.cross(ebs[:, 3, 0, :], ebs[:, 3, 1, :])

    e2 = nm.array([0, 1, 0], dtype=nm.float64)
    for ii in range(4):
        ebs[:, ii, 0, :] = nm.cross(e2, ebs[:, ii, 2, :])
        ebs[:, ii, 1, :] = nm.cross(ebs[:, ii, 2, :], ebs[:, ii, 0, :])
        ebs[:, ii] = ebs[:, ii] / nm.linalg.norm(ebs[:, ii], axis=2)[..., None]

    return ebs

def create_rotation_ops(ebs):
    """
    Create operators associated to rotation DOFs.

    Parameters
    ----------
    ebs : array
        The local bases, shape `(n_el, 4, 3, 3)`.

    Returns
    -------
    rops : array
        The rotation operators, shape `(n_el, 4, 3, 3)`.
    """
    rops = nm.empty((ebs.shape[0], 4, 3, 3), dtype=nm.float64)
    for ii in range(4):
        aux = nm.einsum('cj,ck->cjk', ebs[:, ii, 0], ebs[:, ii, 1])
        rops[:, ii] = aux - aux.transpose((0, 2, 1))

    return rops

def create_strain_transform(mtx_ts):
    r"""
    Create strain tensor transformation matrices, given coordinate
    transformation matrices.

    Notes
    -----
    Expresses :math:`T E T^T` in terms of symmetrix storage as :math:`Q e`,
    with the ordering of components:
    :math:`e = [e_{11}, e_{22}, e_{33}, 2 e_{12}, 2 e_{13}, 2 e_{23}]`.
    """
    mtx_qs = nm.empty((mtx_ts.shape[0], mtx_ts.shape[1], 6, 6),
                      dtype=nm.float64)
    for ie in range(mtx_ts.shape[0]):
        for iq in range(mtx_ts.shape[1]):
            l = mtx_ts[ie, iq, :, 0]
            m = mtx_ts[ie, iq, :, 1]
            n = mtx_ts[ie, iq, :, 2]

            mtx_q = [
                [l[0]**2,     m[0]**2,     n[0]**2,     l[0]*m[0],           n[0]*l[0],           m[0]*n[0]],
                [l[1]**2,     m[1]**2,     n[1]**2,     l[1]*m[1],           n[1]*l[1],           m[1]*n[1]],
                [l[2]**2,     m[2]**2,     n[2]**2,     l[2]*m[2],           n[2]*l[2],           m[2]*n[2]],
                [2*l[0]*l[1], 2*m[0]*m[1], 2*n[0]*n[1], l[0]*m[1]+l[1]*m[0], n[0]*l[1]+n[1]*l[0], m[0]*n[1]+m[1]*n[0]],
                [2*l[2]*l[0], 2*m[2]*m[0], 2*n[2]*n[0], l[2]*m[0]+l[0]*m[2], n[2]*l[0]+n[0]*l[2], m[2]*n[0]+m[0]*n[2]],
                [2*l[1]*l[2], 2*m[1]*m[2], 2*n[1]*n[2], l[1]*m[2]+l[2]*m[1], n[1]*l[2]+n[2]*l[1], m[1]*n[2]+m[2]*n[1]],
                ]
            mtx_qs[ie, iq] = mtx_q

    return mtx_qs

def get_mapping_data(ebs, rops, ps, coors_loc, qp_coors, qp_weights,
                     special_dx3=False):
    r"""
    Compute reference element mapping data for shell10x elements.

    Notes
    -----
    The code assumes that the quadrature points are w.r.t. (:math:`t` =
    thickness of the shell) :math:`[0, 1] \times [0, 1] \times [-t/2, t/2]`
    reference cell and the quadrature weights are multiplied by :math:`t`.
    """
    n_el = coors_loc.shape[0]
    n_qp = qp_weights.shape[0]

    bfu = ps.eval_basis(qp_coors[:, :2].copy())
    bfgu = ps.eval_basis(qp_coors[:, :2].copy(), diff=True)

    nh = ebs[..., -1, :]

    # h = tilde dzeta = t * (dzeta - 0.5), d = dzeta in [0, 1]
    # dh/dd = t
    h = qp_coors[:, -1, None, None]

    dxdxi = nm.empty((n_el, n_qp, 3, 3), dtype=nm.float64)
    coors_loc_3d = coors_loc[:, None, ...] + (nh[:, None, ...] * h)
    dxdxi[..., :2, :] = nm.einsum('qij,cqjk->cqik', bfgu, coors_loc_3d)

    # The factor of t from dh/dd is not here because quadrature weights were
    # multiplied by t. This seems more accurate than having t here and the
    # usual weights.
    # dx_col/dxi_row
    dxdxi[..., 2:, :] = nm.einsum('qij,cjk->cqik', bfu, nh) # * 2
    if special_dx3:
        dxdxi[..., 2:, :] = [0.0, 0.0, 1.0]

    det = nm.linalg.det(dxdxi)

    # *2 to have [-1, 1] reference cell numbers...
    # dxi_col/dx_row
    dxidx = nm.linalg.inv(dxdxi)

    if special_dx3:
        return dxidx, det

    bfg = nm.zeros((n_el, n_qp, 3, 3, 24), dtype=nm.float64)

    bfg[..., 0, :2, :4] = bfgu
    bfg[..., 1, :2, 4:8] = bfgu
    bfg[..., 2, :2, 8:12] = bfgu

    n = None

    # bfg[:, 0, :2, 12:16] = (bfgu * rops[..., 0, 0][:, n, n, :]) * h = 0
    bfg[..., 0, :2, 16:20] = (bfgu * rops[..., 0, 1][:, n, n, :]) * h
    bfg[..., 0, :2, 20:24] = (bfgu * rops[..., 0, 2][:, n, n, :]) * h
    # bfg[:, 0, 2:, 12:16] = bfu * rops[..., 0, 0][:, n, n, :] = 0
    bfg[..., 0, 2:, 16:20] = bfu * rops[..., 0, 1][:, n, n, :]
    bfg[..., 0, 2:, 20:24] = bfu * rops[..., 0, 2][:, n, n, :]

    bfg[..., 1, :2, 12:16] = (bfgu * rops[..., 1, 0][:, n, n, :]) * h
    # bfg[:, 1, :2, 16:20] = (bfgu * rops[..., 1, 1][:, n, n, :]) * h = 0
    bfg[..., 1, :2, 20:24] = (bfgu * rops[..., 1, 2][:, n, n, :]) * h
    bfg[..., 1, 2:, 12:16] = bfu * rops[..., 1, 0][:, n, n, :]
    # bfg[:, 1, 2:, 16:20] = bfu * rops[..., 1, 1][:, n, n, :] = 0
    bfg[..., 1, 2:, 20:24] = bfu * rops[..., 1, 2][:, n, n, :]

    bfg[..., 2, :2, 12:16] = (bfgu * rops[..., 2, 0][:, n, n, :]) * h
    bfg[..., 2, :2, 16:20] = (bfgu * rops[..., 2, 1][:, n, n, :]) * h
    # bfg[:, 2, :2, 20:24] = (bfgu * rops[..., 2, 2][:, n, n, :]) * h = 0
    bfg[..., 2, 2:, 12:16] = bfu * rops[..., 2, 0][:, n, n, :]
    bfg[..., 2, 2:, 16:20] = bfu * rops[..., 2, 1][:, n, n, :]
    # bfg[:, 2, 2:, 20:24] = bfu * rops[..., 2, 2][:, n, n, :] = 0

    # ... related to [0, 1] vs. [-1, 1]:
    # bfg[:, :, 2, 12:] /= 0.5 when *2 in dxidx above

    bfgm = nm.einsum('cqij,cqkjl->cqikl', dxidx, bfg)

    return coors_loc_3d, bfu, bfgm, dxidx, det

def get_dsg_strain(coors_loc, qp_coors):
    r"""
    Compute DSG strain components.

    Returns
    -------
    dsg : array
        The strain matrix components corresponding to :math:`e_{13}`,
        :math:`e_{23}`, shape `(n_el, n_qp, 2, 24)`.

    Notes
    -----
    Involves :math:`w`, :math:`\alpha`, :math:`\beta` DOFs.
    """
    dsg = nm.zeros((coors_loc.shape[0], qp_coors.shape[0], 2, 24),
                   dtype=nm.float64)

    x1, x2, x3, x4 = coors_loc.transpose((1, 0, 2))

    o = nm.ones((coors_loc.shape[0], 1), dtype=nm.float64)

    v13_1 = nm.c_[-o, o, -0.5 * (x2[:, 1] - x1[:, 1]),
                   -0.5 * (x2[:, 1] - x1[:, 1]),
                   0.5 * (x2[:, 0] - x1[:, 0]),
                   0.5 * (x2[:, 0] - x1[:, 0])]
    i13_1 = nm.r_[8, 9, 12, 13, 16, 17]
    v13_2 = nm.c_[o, -o, -0.5 * (x3[:, 1] - x4[:, 1]),
                  -0.5 * (x3[:, 1] - x4[:, 1]),
                  0.5 * (x3[:, 0] - x4[:, 0]),
                  0.5 * (x3[:, 0] - x4[:, 0])]
    i13_2 = nm.r_[10, 11, 14, 15, 18, 19]

    v23_1 = nm.c_[-o, o, -0.5 * (x4[:, 1] - x1[:, 1]),
                   -0.5 * (x4[:, 1] - x1[:, 1]),
                   0.5 * (x4[:, 0] - x1[:, 0]),
                   0.5 * (x4[:, 0] - x1[:, 0])]
    i23_1 = nm.r_[8, 11, 12, 15, 16, 19]
    v23_2 = nm.c_[-o, o, -0.5 * (x3[:, 1] - x2[:, 1]),
                   -0.5 * (x3[:, 1] - x2[:, 1]),
                   0.5 * (x3[:, 0] - x2[:, 0]),
                   0.5 * (x3[:, 0] - x2[:, 0])]
    i23_2 = nm.r_[9, 10, 13, 14, 17, 18]

    c1, c2 = qp_coors[:, 0], qp_coors[:, 1]

    # Assuming qp in [-1, 1]
    # dsg[:, 0, i13_1] = 0.25 * (1.0 - c2)[:, None] * v13_1
    # dsg[:, 0, i13_2] = 0.25 * (1.0 + c2)[:, None] * v13_2
    # dsg[:, 1, i23_1] = 0.25 * (1.0 - c1)[:, None] * v23_1
    # dsg[:, 1, i23_2] = 0.25 * (1.0 + c1)[:, None] * v23_2

    # Assuming qp in [0, 1]
    n = None
    dsg[..., 0, i13_1] = (1.0 - c2)[n, :, n] * v13_1[:, n, :]
    dsg[..., 0, i13_2] = (c2)[n, :, n] * v13_2[:, n, :]
    dsg[..., 1, i23_1] = (1.0 - c1)[n, :, n] * v23_1[:, n, :]
    dsg[..., 1, i23_2] = (c1)[n, :, n] * v23_2[:, n, :]

    return dsg

def create_strain_matrix(bfgm, dxidx, dsg):
    """
    Create the strain operator matrix.
    """
    tg = create_strain_transform(dxidx)

    mtx_b = nm.concatenate((bfgm[..., 0:1, 0, :],
                            bfgm[..., 1:2, 1, :],
                            bfgm[..., 2:3, 2, :],
                            bfgm[..., 0:1, 1, :] + bfgm[..., 1:2, 0, :],
                            nm.einsum('ciaj,cijk->ciak', tg[..., 4:5, 4:], dsg),
                            nm.einsum('ciaj,cijk->ciak', tg[..., 5:6, 4:], dsg)),
                           2)

    return mtx_b

def add_eas_dofs(mtx_b, qp_coors, det, det0, dxidx0):
    """
    Add additional strain components [Andelfinger and Ramm] (7 parameters to be
    condensed out).
    """
    tg0 = create_strain_transform(dxidx0)

    # mtx_b is the same as Bk -> mtx_ee must be the same as Ba
    # -> qp coors transformation below is needed.
    qp_coors = 8 * (qp_coors - 0.5)

    c1, c2 = qp_coors[:, 0], qp_coors[:, 1]
    cc = c1 * c2

    mtx_e = nm.zeros((qp_coors.shape[0], 6, 7), dtype=nm.float64)

    mtx_e[:, 0, 0] = mtx_e[:, 3, 2] = c1
    mtx_e[:, 1, 1] = mtx_e[:, 3, 3] = c2
    mtx_e[:, 0, 4] = mtx_e[:, 1, 5] = mtx_e[:, 3, 6] = cc

    dd = (det0 / det)[..., None, None]
    mtx_ee = dd * nm.einsum('cij,qjk->cqik', tg0[:, 0], mtx_e)

    mtx_be = nm.concatenate((mtx_b, mtx_ee), 3)

    return mtx_be

def rotate_elastic_tensor(mtx_d, bfu, ebs):
    """
    Rotate the elastic tensor into the local coordinate system of each cell.
    The local coordinate system results from interpolation of `ebs` with the
    bilinear basis.
    """
    # ! mtx_ts is not orthonormal for non-planar cells!
    mtx_ts = nm.einsum('qi,cijk->cqjk', bfu[:, 0, :], ebs)
    # print nm.einsum('cqji,cqjk->cqik', mtx_ts, mtx_ts)

    # ? double entries here?
    tt = create_strain_transform(mtx_ts)

    mtx_dr = ddot(ddot(tt, mtx_d, 'ATB'), tt, 'AB')
    return mtx_dr

def create_drl_transform(ebs):
    """
    Create the transformation matrix for locking of the drilling rotations.
    """
    n_el = ebs.shape[0]
    mtx_drl = nm.tile(nm.eye(24, dtype=nm.float64), (n_el, 1, 1))
    for ii in range(4):
        nh = ebs[:, ii, -1, :]
        mtx = nm.c_[nm.c_[1.0 - nh[:, 0]**2,
                          -nh[:, 0]*nh[:, 1],
                          -nh[:, 0]*nh[:, 2]],
                    nm.c_[-nh[:, 0]*nh[:, 1],
                           1.0 - nh[:, 1]**2,
                           -nh[:, 1]*nh[:, 2]],
                    nh].reshape(n_el, 3, 3)
        mtx_drl[:, 12+ii:12+ii+9:4, 12+ii:12+ii+9:4] = mtx

    return mtx_drl

def lock_drilling_rotations(mtx, ebs, coefs):
    """
    Lock the drilling rotations in the stiffness matrix.
    """
    mtx_drl = create_drl_transform(ebs)
    mtx_idrl = nm.linalg.inv(mtx_drl)
    mtx_tr = ddot(ddot(mtx_drl, mtx), mtx_idrl)

    idrl = nm.arange(20, 24)
    mtx_tr[:, idrl, idrl] = coefs[:, None]

    mtx2 = ddot(ddot(mtx_idrl, mtx_tr), mtx_drl)

    return mtx2
