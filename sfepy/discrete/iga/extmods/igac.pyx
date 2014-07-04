# -*- Mode: Python -*-
cimport cython

import numpy as np
cimport numpy as np

from sfepy.discrete.fem.extmods.types cimport int32, uint32, float64

from sfepy.discrete.fem.extmods._fmfield cimport (FMField,
                                                  array2fmfield4,
                                                  array2fmfield2,
                                                  array2fmfield1,
                                                  array2pint2,
                                                  array2pint1,
                                                  array2puint1,
                                                  fmf_alloc,
                                                  fmf_free)

cdef extern from 'nurbs.h':
    cdef void _ravel_multi_index \
         'ravel_multi_index'(uint32 *index, uint32 *indices,
                             uint32 *shape, uint32 num)
    cdef void _unravel_index \
         'unravel_index'(uint32 *indices, uint32 index,
                         uint32 *shape, uint32 num)

    cdef int32 _eval_bernstein_basis \
         'eval_bernstein_basis'(FMField *funs, FMField *ders,
                                float64 x, uint32 degree)
    cdef int32 _eval_nurbs_basis_tp \
         'eval_nurbs_basis_tp'(FMField *R, FMField *dR_dx, FMField *det,
                               FMField *dR_dxi,
                               FMField *dx_dxi, FMField *dxi_dx,
                               FMField *B, FMField *dB_dxi,
                               FMField *N, FMField *dN_dxi,
                               FMField *qp, uint32 ie, FMField *control_points,
                               FMField *weights, int32 *degrees, int32 dim,
                               FMField *cs,
                               int32 *conn, int32 n_el, int32 n_ep)

def eval_bernstein_basis(np.ndarray funs not None,
                         np.ndarray ders not None,
                         float64 x,
                         uint32 degree):
    cdef int32 ret
    cdef FMField _funs[1], _ders[1]

    array2fmfield1(_funs, funs)
    array2fmfield1(_ders, ders)

    ret = _eval_bernstein_basis(_funs, _ders, x, degree)
    return ret

def eval_mapping_data_in_qp(np.ndarray[float64, mode='c', ndim=2] qps not None,
                            np.ndarray[float64, mode='c', ndim=2]
                            control_points not None,
                            np.ndarray[float64, mode='c', ndim=1]
                            weights not None,
                            np.ndarray[int32, mode='c', ndim=1]
                            degrees not None,
                            cs not None,
                            np.ndarray[int32, mode='c', ndim=2] conn not None,
                            np.ndarray[uint32, mode='c', ndim=1] cells=None):
    """
    Evaluate data required for the isogeometric domain reference mapping in the
    given quadrature points. The quadrature points are the same for all Bezier
    elements and should correspond to the Bernstein basis degree.

    Parameters
    ----------
    qps : array
        The quadrature points coordinates with components in [0, 1] reference
        element domain.
    control_points : array
        The NURBS control points.
    weights : array
        The NURBS weights.
    degrees : sequence of ints or int
        The basis degrees in each parametric dimension.
    cs : list of lists of 2D arrays
        The element extraction operators in each parametric dimension.
    conn : array
        The connectivity of the global NURBS basis.
    cells : array, optional
        If given, use only the given Bezier elements.

    Returns
    -------
    bfs : array
        The NURBS shape functions in the physical quadrature points of all
        elements.
    bfgs : array
        The NURBS shape functions derivatives w.r.t. the physical coordinates
        in the physical quadrature points of all elements.
    dets : array
        The Jacobians of the mapping to the unit reference element in the
        physical quadrature points of all elements.
    """
    cdef uint32 ii, ie, n_qp, n_efun
    cdef int32 n_el, n_ep, dim, aux
    cdef uint32 *_cells
    cdef int32 *_degrees, *_conn
    cdef FMField _bf[1], _bfg[1], _det[1]
    cdef FMField _bfg_dxi[1], _dx_dxi[1], _dxi_dx[1]
    cdef FMField _qp[1], _control_points[1], _weights[1]
    cdef FMField _cs[3]
    cdef FMField _B[3], _dB_dxi[3], _N[3], _dN_dxi[3]
    cdef np.ndarray[float64, mode='c', ndim=4] bfs, bfgs, dets

    if cells is None:
        cells = np.arange(conn.shape[0], dtype=np.uint32)

    degrees = np.asarray(degrees, dtype=np.int32)

    n_el = len(cells)
    n_qp = qps.shape[0]
    dim = control_points.shape[1]
    n_efuns = degrees + 1
    n_efun = np.prod(n_efuns)

    # Output Jacobians.
    dets = np.empty((n_el, n_qp, 1, 1), dtype=np.float64)

    # Output shape functions.
    bfs = np.empty((n_el, n_qp, 1, n_efun), dtype=np.float64)

    # Output gradients of shape functions.
    bfgs = np.empty((n_el, n_qp, dim, n_efun), dtype=np.float64)

    # Setup C termporary arrays.
    bfg_dxi = np.empty((1, 1, dim, n_efun), dtype=np.float64)
    dx_dxi = np.empty((1, 1, dim, dim), dtype=np.float64)
    dxi_dx = np.empty((1, 1, dim, dim), dtype=np.float64)

    for ii in range(0, dim):
        fmf_alloc(_B + ii, 1, 1, n_efuns[ii], 1)
        fmf_alloc(_dB_dxi + ii, 1, 1, n_efuns[ii], 1)
        fmf_alloc(_N + ii, 1, 1, n_efuns[ii], 1)
        fmf_alloc(_dN_dxi + ii, 1, 1, n_efuns[ii], 1)

    # Assign to C structures.
    array2fmfield4(_bfg_dxi, bfg_dxi)
    array2fmfield4(_dx_dxi, dx_dxi)
    array2fmfield4(_dxi_dx, dxi_dx)
    array2fmfield2(_control_points, control_points)
    array2fmfield1(_weights, weights)
    for ii in range(dim):
        array2fmfield4(_cs + ii, cs[ii])

    array2pint1(&_degrees, &dim, degrees)
    array2pint2(&_conn, &aux, &n_ep, conn)

    _bf.offset = _bfg.offset = _det.offset = _qp.offset = -1
    _bf.nAlloc = _bfg.nAlloc = _det.nAlloc = _qp.nAlloc = -1
    _bf.nCell = _bfg.nCell = _det.nCell = _qp.nCell = 1
    _bf.nLev = _bfg.nLev = _det.nLev = _qp.nLev = 1
    _bf.nRow = _det.nRow = _qp.nRow = 1
    _bfg.nRow = dim
    _bf.nCol = _bfg.nCol = n_efun
    _det.nCol = 1
    _qp.nCol = dim

    _bf.val = _bf.val0 = &bfs[0, 0, 0, 0]
    _bfg.val = _bfg.val0 = &bfgs[0, 0, 0, 0]
    _det.val = _det.val0 = &dets[0, 0, 0, 0]
    _qp.val = _qp.val0 = &qps[0, 0]

    # Loop over elements.
    _cells = &cells[0]
    for iseq in range(0, n_el):
        ie = _cells[iseq]

        _qp.val = _qp.val0

        # Loop over quadrature points.
        for iqp in range(0, n_qp):
            _eval_nurbs_basis_tp(_bf, _bfg, _det,
                                 _bfg_dxi,
                                 _dx_dxi, _dxi_dx,
                                 _B, _dB_dxi, _N, _dN_dxi,
                                 _qp, ie,
                                 _control_points, _weights,
                                 _degrees, dim, _cs, _conn, n_el, n_ep)
            _bf.val += n_efun
            _bfg.val += dim * n_efun
            _det.val += 1
            _qp.val += dim

    for ii in range(0, dim):
        fmf_free(_B + ii)
        fmf_free(_dB_dxi + ii)
        fmf_free(_N + ii)
        fmf_free(_dN_dxi + ii)

    return bfs, bfgs, dets

def eval_variable_in_qp(np.ndarray[float64, mode='c', ndim=2] variable not None,
                        np.ndarray[float64, mode='c', ndim=2] qps not None,
                        np.ndarray[float64, mode='c', ndim=2]
                        control_points not None,
                        np.ndarray[float64, mode='c', ndim=1] weights not None,
                        np.ndarray[int32, mode='c', ndim=1] degrees not None,
                        cs not None,
                        np.ndarray[int32, mode='c', ndim=2] conn not None,
                        np.ndarray[uint32, mode='c', ndim=1] cells=None):
    """
    Evaluate a field variable in the given quadrature points. The quadrature
    points are the same for all Bezier elements and should correspond to the
    Bernstein basis degree. The field variable is defined by its DOFs - the
    coefficients of the NURBS basis.

    Parameters
    ----------
    variable : array
        The DOF values of the variable with n_c components, shape (:, n_c).
    qps : array
        The quadrature points coordinates with components in [0, 1] reference
        element domain.
    control_points : array
        The NURBS control points.
    weights : array
        The NURBS weights.
    degrees : sequence of ints or int
        The basis degrees in each parametric dimension.
    cs : list of lists of 2D arrays
        The element extraction operators in each parametric dimension.
    conn : array
        The connectivity of the global NURBS basis.
    cells : array, optional
        If given, use only the given Bezier elements.

    Returns
    -------
    coors : array
        The physical coordinates of the quadrature points of all elements.
    vals : array
        The field variable values in the physical quadrature points.
    dets : array
        The Jacobians of the mapping to the unit reference element in the
        physical quadrature points.
    """
    cdef uint32 ii, ie, n_qp, n_efun, nc, ir, ic
    cdef int32 n_el, n_ep, dim, aux
    cdef uint32 *_cells
    cdef int32 *_degrees, *_conn, *ec
    cdef float64 val
    cdef FMField _bf[1], _bfg[1], _det[1], _vals[1], _coors[1]
    cdef FMField _bfg_dxi[1], _dx_dxi[1], _dxi_dx[1]
    cdef FMField _qp[1], _variable[1], _control_points[1], _weights[1]
    cdef FMField _cs[3]
    cdef FMField _B[3], _dB_dxi[3], _N[3], _dN_dxi[3]
    cdef np.ndarray[float64, mode='c', ndim=2] coors, vals, dets

    if cells is None:
        cells = np.arange(conn.shape[0], dtype=np.uint32)

    n_el = len(cells)
    n_qp = qps.shape[0]
    dim = control_points.shape[1]
    n_efuns = degrees + 1
    n_efun = np.prod(n_efuns)
    nc = variable.shape[1]

    # Output values of the variable.
    vals = np.empty((n_el * n_qp, nc), dtype=np.float64)

    # Output physical coordinates of QPs.
    coors = np.empty((n_el * n_qp, dim), dtype=np.float64)

    # Output Jacobians.
    dets = np.empty((n_el * n_qp, 1), dtype=np.float64)

    # Setup C termporary arrays.
    bf = np.empty((n_efun,), dtype=np.float64)
    bfg = np.empty((dim, n_efun), dtype=np.float64)
    bfg_dxi = np.empty((1, 1, dim, n_efun), dtype=np.float64)
    dx_dxi = np.empty((1, 1, dim, dim), dtype=np.float64)
    dxi_dx = np.empty((1, 1, dim, dim), dtype=np.float64)

    for ii in range(0, dim):
        fmf_alloc(_B + ii, 1, 1, n_efuns[ii], 1)
        fmf_alloc(_dB_dxi + ii, 1, 1, n_efuns[ii], 1)
        fmf_alloc(_N + ii, 1, 1, n_efuns[ii], 1)
        fmf_alloc(_dN_dxi + ii, 1, 1, n_efuns[ii], 1)

    # Assign to C structures.
    array2fmfield1(_bf, bf)
    array2fmfield2(_bfg, bfg)
    array2fmfield4(_bfg_dxi, bfg_dxi)
    array2fmfield4(_dx_dxi, dx_dxi)
    array2fmfield4(_dxi_dx, dxi_dx)
    array2fmfield2(_variable, variable)
    array2fmfield2(_control_points, control_points)
    array2fmfield1(_weights, weights)
    for ii in range(dim):
        array2fmfield4(_cs + ii, cs[ii])

    array2pint1(&_degrees, &dim, degrees)
    array2pint2(&_conn, &aux, &n_ep, conn)

    _vals.offset = _coors.offset = _det.offset = _qp.offset = -1
    _vals.nAlloc = _coors.nAlloc = _det.nAlloc = _qp.nAlloc = -1
    _vals.nCell = _coors.nCell = _det.nCell = _qp.nCell = 1
    _vals.nLev = _coors.nLev = _det.nLev = _qp.nLev = 1
    _vals.nRow = _coors.nRow = _det.nRow = _qp.nRow = 1
    _vals.nCol = nc
    _coors.nCol = dim
    _det.nCol = 1
    _qp.nCol = dim

    _vals.val = _vals.val0 = &vals[0, 0]
    _coors.val = _coors.val0 = &coors[0, 0]
    _det.val = _det.val0 = &dets[0, 0]
    _qp.val = _qp.val0 = &qps[0, 0]

    # Loop over elements.
    _cells = &cells[0]
    for iseq in range(0, n_el):
        ie = _cells[iseq]

        ec = _conn + n_ep * ie

        _qp.val = _qp.val0

        # Loop over quadrature points.
        for iqp in range(0, n_qp):
            _eval_nurbs_basis_tp(_bf, _bfg, _det,
                                 _bfg_dxi,
                                 _dx_dxi, _dxi_dx,
                                 _B, _dB_dxi, _N, _dN_dxi,
                                 _qp, ie,
                                 _control_points, _weights,
                                 _degrees, dim, _cs, _conn, n_el, n_ep)

            # vals[ii, :] = np.dot(bf, variable[ec])
            for ir in range(0, nc):
                _vals.val[ir] = 0.0

                for ic in range(0, n_efun):
                    val = _variable.val[ec[ic] * nc + ir]
                    _vals.val[ir] += _bf.val[ic] * val

            # coors[ii, :] = np.dot(bf, control_points[ec])
            for ir in range(0, dim):
                _coors.val[ir] = 0.0

                for ic in range(0, n_efun):
                    val = _control_points.val[ec[ic] * dim + ir]
                    _coors.val[ir] += _bf.val[ic] * val

            _vals.val += nc
            _coors.val += dim
            _det.val += 1
            _qp.val += dim

    return coors, vals, dets

def eval_in_tp_coors(np.ndarray[float64, mode='c', ndim=2] variable,
                     indices not None,
                     ref_coors not None,
                     np.ndarray[float64, mode='c', ndim=2]
                     control_points not None,
                     np.ndarray[float64, mode='c', ndim=1] weights not None,
                     np.ndarray[int32, mode='c', ndim=1] degrees not None,
                     cs not None,
                     np.ndarray[int32, mode='c', ndim=2] conn not None):
    """
    Evaluate a field variable (if given) or the NURBS geometry in the given
    tensor-product reference coordinates. The field variable is defined by its
    DOFs - the coefficients of the NURBS basis.

    Parameters
    ----------
    variable : array
        The DOF values of the variable with n_c components, shape (:, n_c).
    indices : list of arrays
        The indices of knot spans for each axis, defining the Bezier element
        numbers.
    ref_coors : list of arrays
        The reference coordinates in [0, 1] for each knot span for each axis,
        defining the reference coordinates in the Bezier elements given by
        `indices`.
    control_points : array
        The NURBS control points.
    weights : array
        The NURBS weights.
    degrees : sequence of ints or int
        The basis degrees in each parametric dimension.
    cs : list of lists of 2D arrays
        The element extraction operators in each parametric dimension.
    conn : array
        The connectivity of the global NURBS basis.

    Returns
    -------
    out : array
        The field variable values or NURBS geometry coordinates for the given
        reference coordinates.
    """
    cdef uint32 ii, ip, ie, n_efun, nc, ir, ic, n_vals, uaux
    cdef int32 n_el, n_ep, dim, aux
    cdef int32 *_degrees, *_conn, *ec
    cdef uint32 igrid[3], shape[3], iis[3], n_els[3]
    cdef uint32 *_indices[3], **puaux
    cdef FMField _bf[1], _bfg[1], _det[1], _vals[1]
    cdef FMField _bfg_dxi[1], _dx_dxi[1], _dxi_dx[1]
    cdef FMField _rc[1], _control_points[1], _weights[1]
    cdef FMField _cs[3], _ref_coors[3]
    cdef np.ndarray[float64, mode='c', ndim=2] _evals
    cdef FMField _B[3], _dB_dxi[3], _N[3], _dN_dxi[3]

    dim = control_points.shape[1]
    n_efuns = degrees + 1
    n_efun = np.prod(n_efuns)

    n_vals = 1
    for ii in range(0, dim):
        shape[ii] = len(ref_coors[ii])
        n_vals *= shape[ii]

    # Output values.
    if variable is None:
        nc = dim
        out = np.zeros((n_vals, dim), dtype=np.float64)
        _evals = control_points

    else:
        nc = variable.shape[1]
        out = np.zeros((n_vals, nc), dtype=np.float64)
        _evals = variable

    # Setup C termporary arrays.
    rc = np.empty((dim,), dtype=np.float64)
    det = np.empty((1,), dtype=np.float64)
    bf = np.empty((n_efun,), dtype=np.float64)
    bfg = np.empty((dim, n_efun), dtype=np.float64)
    bfg_dxi = np.empty((1, 1, dim, n_efun), dtype=np.float64)
    dx_dxi = np.empty((1, 1, dim, dim), dtype=np.float64)
    dxi_dx = np.empty((1, 1, dim, dim), dtype=np.float64)

    for ii in range(0, dim):
        fmf_alloc(_B + ii, 1, 1, n_efuns[ii], 1)
        fmf_alloc(_dB_dxi + ii, 1, 1, n_efuns[ii], 1)
        fmf_alloc(_N + ii, 1, 1, n_efuns[ii], 1)
        fmf_alloc(_dN_dxi + ii, 1, 1, n_efuns[ii], 1)

    # Assign to C structures.
    array2fmfield1(_rc, rc)
    array2fmfield1(_det, det)
    array2fmfield1(_bf, bf)
    array2fmfield2(_bfg, bfg)
    array2fmfield4(_bfg_dxi, bfg_dxi)
    array2fmfield4(_dx_dxi, dx_dxi)
    array2fmfield4(_dxi_dx, dxi_dx)
    array2fmfield2(_control_points, control_points)
    array2fmfield1(_weights, weights)
    for ii in range(dim):
        array2fmfield4(_cs + ii, cs[ii])
        n_els[ii] = (_cs + ii).nCell;

    array2pint1(&_degrees, &dim, degrees)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    for ii in range(0, dim):
        puaux = _indices + ii
        array2puint1(puaux, &uaux, indices[ii])
        array2fmfield1(_ref_coors + ii, ref_coors[ii])

    for ip in range(0, n_vals):
        _unravel_index(igrid, ip, shape, dim)

        for ii in range(0, dim):
            iis[ii] = _indices[ii][igrid[ii]]
            _rc.val[ii] = (_ref_coors + ii).val[igrid[ii]]
        _ravel_multi_index(&ie, iis, n_els, dim)

        _eval_nurbs_basis_tp(_bf, _bfg, _det,
                             _bfg_dxi,
                             _dx_dxi, _dxi_dx,
                             _B, _dB_dxi, _N, _dN_dxi,
                             _rc, ie,
                             _control_points, _weights,
                             _degrees, dim, _cs, _conn, n_el, n_ep)

        # vals[ip, :] = np.dot(bf, variable[ec])
        ec = _conn + n_ep * ie;
        array2fmfield1(_vals, out[ip, :])
        for ir in range(0, nc):
            _vals.val[ir] = 0.0

            for ic in range(0, n_efun):
                _vals.val[ir] += _bf.val[ic] * _evals[ec[ic], ir]

    return out
