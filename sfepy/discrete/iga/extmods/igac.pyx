# -*- Mode: Python -*-
# cython: language_level=3
cimport cython

import numpy as np
cimport numpy as np

from sfepy.discrete.common.extmods.types cimport int32, uint32, float64

from sfepy.discrete.common.extmods._fmfield cimport (FMField,
                                                     array2fmfield4,
                                                     array2fmfield3,
                                                     array2fmfield2,
                                                     array2fmfield1,
                                                     array2pint2,
                                                     array2pint1,
                                                     array2puint1,
                                                     fmf_pretend_nc,
                                                     fmf_alloc,
                                                     fmf_free)

cdef extern from 'common.h':
    void *pyalloc(size_t size)
    void pyfree(void *pp)

cdef extern from 'nurbs.h':
    ctypedef struct NURBSContext:
        int32 (*get_xi_dist)(float64 *pdist, FMField *xi,
                             FMField *point, FMField *e_coors,
                             void *_ctx)
        int32 (*eval_basis)(FMField *out, FMField *coors, int32 diff,
                            void *_ctx)
        int32 iel # >= 0.
        int32 is_dx # 1 => apply reference mapping to gradient.
        FMField e_coors_max[1] # Buffer for coordinates of element nodes.

        FMField control_points[1]
        FMField weights[1]
        int32 *degrees
        int32 dim
        FMField cs[3]
        int32 *conn
        int32 n_cell
        int32 n_efun

        FMField bf[1]
        FMField bfg[1]

        FMField R[1]
        FMField dR_dxi[1]
        FMField dR_dx[1]

        FMField B[3]
        FMField dB_dxi[3]
        FMField N[3]
        FMField dN_dxi[3]

        int32 reuse

        int32 has_bernstein
        int32 is_nurbs

        int32 i_max
        float64 newton_eps

    cdef void _ravel_multi_index \
         'ravel_multi_index'(uint32 *index, uint32 *indices,
                             uint32 *shape, uint32 num)
    cdef void _unravel_index \
         'unravel_index'(uint32 *indices, uint32 index,
                         uint32 *shape, uint32 num)

    void _print_context_nurbs \
         'print_context_nurbs'(NURBSContext *ctx)

    int32 _get_xi_dist \
          'get_xi_dist'(float64 *pdist, FMField *xi,
                        FMField *point, FMField *e_coors,
                        void *_ctx)

    int32 _eval_basis_nurbs \
          'eval_basis_nurbs'(FMField *out, FMField *coors, int32 diff,
                             void *_ctx)

    cdef int32 _eval_bernstein_basis \
         'eval_bernstein_basis'(FMField *funs, FMField *ders,
                                float64 x, uint32 degree)
    cdef int32 _eval_bspline_basis_tp \
         'eval_bspline_basis_tp'(FMField *R, FMField *dR_dx, FMField *det,
                                 FMField *dR_dxi,
                                 FMField *dx_dxi, FMField *dxi_dx,
                                 FMField *B, FMField *dB_dxi,
                                 FMField *N, FMField *dN_dxi,
                                 FMField *qp, uint32 ie,
                                 FMField *control_points,
                                 int32 *degrees, int32 dim,
                                 FMField *cs,
                                 int32 *conn, int32 n_el, int32 n_ep,
                                 int32 has_bernstein, int32 is_dx)
    cdef int32 _eval_nurbs_basis_tp \
         'eval_nurbs_basis_tp'(FMField *R, FMField *dR_dx, FMField *det,
                               FMField *dR_dxi,
                               FMField *dx_dxi, FMField *dxi_dx,
                               FMField *B, FMField *dB_dxi,
                               FMField *N, FMField *dN_dxi,
                               FMField *qp, uint32 ie, FMField *control_points,
                               FMField *weights, int32 *degrees, int32 dim,
                               FMField *cs,
                               int32 *conn, int32 n_el, int32 n_ep,
                               int32 has_bernstein, int32 is_dx)

cdef class CNURBSContext:

    cdef NURBSContext *ctx

    # Store arrays to prevent their deallocation in Python.
    cdef readonly np.ndarray e_coors_max # Auxiliary buffer.

    cdef readonly np.ndarray bf # Auxiliary buffer.
    cdef readonly np.ndarray bfg # Auxiliary buffer.

    cdef readonly np.ndarray R # Auxiliary buffer.
    cdef readonly np.ndarray dR_dxi # Auxiliary buffer.
    cdef readonly np.ndarray dR_dx # Auxiliary buffer.

    cdef readonly np.ndarray bufBN # Auxiliary buffer.

    property iel:

        def __get__(self):
            return self.ctx.iel

        def __set__(self, int32 iel):
            assert iel < self.ctx.n_cell
            self.ctx.iel = iel

    def __cinit__(self,
                  np.ndarray[float64, mode='c', ndim=2] control_points not None,
                  np.ndarray[float64, mode='c', ndim=1] weights not None,
                  np.ndarray[int32, mode='c', ndim=1] degrees not None,
                  cs not None,
                  np.ndarray[int32, mode='c', ndim=2] conn not None,
                  int32 i_max=100,
                  float64 newton_eps=1e-8):
        cdef NURBSContext *ctx
        cdef int32 ii, num
        cdef float64 *buf
        cdef np.ndarray[float64, mode='c', ndim=2] _e_coors_max
        cdef np.ndarray[float64, mode='c', ndim=1] _bf
        cdef np.ndarray[float64, mode='c', ndim=2] _bfg
        cdef np.ndarray[float64, mode='c', ndim=1] _R
        cdef np.ndarray[float64, mode='c', ndim=2] _dR_dxi
        cdef np.ndarray[float64, mode='c', ndim=2] _dR_dx
        cdef np.ndarray[float64, mode='c', ndim=1] _bufBN

        ctx = self.ctx = <NURBSContext *> pyalloc(sizeof(NURBSContext))

        if ctx is NULL:
            raise MemoryError()

        ctx.get_xi_dist = &_get_xi_dist
        ctx.eval_basis = &_eval_basis_nurbs
        ctx.iel = 0
        ctx.is_dx = 1

        array2fmfield2(ctx.control_points, control_points)
        array2fmfield1(ctx.weights, weights)

        array2pint1(&ctx.degrees, &ctx.dim, degrees)

        for ii in range(ctx.dim):
            array2fmfield4(ctx.cs + ii, cs[ii])

        array2pint2(&ctx.conn, &ctx.n_cell, &ctx.n_efun, conn)

        _e_coors_max = self.e_coors_max = np.zeros((ctx.n_efun, ctx.dim),
                                                   dtype=np.float64)
        array2fmfield2(ctx.e_coors_max, _e_coors_max)

        _bf = self.bf = np.zeros((ctx.n_efun,), dtype=np.float64)
        array2fmfield1(ctx.bf, _bf)

        _bfg = self.bfg = np.zeros((ctx.dim, ctx.n_efun),
                                   dtype=np.float64)
        array2fmfield2(ctx.bfg, _bfg)

        _R = self.R = np.zeros((ctx.n_efun,), dtype=np.float64)
        array2fmfield1(ctx.R, _R)

        _dR_dxi = self.dR_dxi = np.zeros((ctx.dim, ctx.n_efun),
                                         dtype=np.float64)
        array2fmfield2(ctx.dR_dxi, _dR_dxi)

        _dR_dx = self.dR_dx = np.zeros((ctx.dim, ctx.n_efun),
                                       dtype=np.float64)
        array2fmfield2(ctx.dR_dx, _dR_dx)

        n_efuns = degrees + 1

        _bufBN = self.bufBN = np.zeros((4 * n_efuns.sum(),), dtype=np.float64)
        buf = &_bufBN[0]
        for ii in range(0, ctx.dim):
            num = n_efuns[ii]
            fmf_pretend_nc(ctx.B + ii, 1, 1, num, 1, buf)
            buf += num
            fmf_pretend_nc(ctx.dB_dxi + ii, 1, 1, num, 1, buf)
            buf += num
            fmf_pretend_nc(ctx.N + ii, 1, 1, num, 1, buf)
            buf += num
            fmf_pretend_nc(ctx.dN_dxi + ii, 1, 1, num, 1, buf)
            buf += num

        ctx.reuse = 0

        ctx.is_nurbs = is_nurbs(weights)

        ctx.i_max = i_max
        ctx.newton_eps = newton_eps

    def __dealloc__(self):
        pyfree(self.ctx)

    def __str__(self):
        return 'CNURBSContext'

    def cprint(self):
        _print_context_nurbs(self.ctx)

    def evaluate(self, np.ndarray[float64, mode='c', ndim=2] coors not None,
                 int32 diff=False, **kwargs):
        cdef int32 n_coor = coors.shape[0]
        cdef int32 n_efun = self.ctx.n_efun
        cdef int32 dim = coors.shape[1]
        cdef int32 bdim
        cdef FMField[1] _out, _coors

        ctx = self.ctx

        if diff:
            bdim = dim

        else:
            bdim = 1

        cdef np.ndarray[float64, ndim=3] out = np.zeros((n_coor, bdim, n_efun),
                                                        dtype=np.float64)

        array2fmfield3(_out, out)
        array2fmfield2(_coors, coors)

        self.ctx.eval_basis(_out, _coors, diff, ctx)

        return out

def is_nurbs(np.ndarray[float64, mode='c', ndim=1] weights not None):
    """
    Return True if some weights are not one.
    """
    return np.any(weights != 1.0)

def eval_bernstein_basis(np.ndarray funs not None,
                         np.ndarray ders not None,
                         float64 x,
                         uint32 degree):
    cdef int32 ret
    cdef FMField[1] _funs, _ders

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
    cdef uint32 ii, ie, n_qp, n_efun, nf
    cdef int32 n_el, n_ep, dim, aux
    cdef uint32 *_cells
    cdef (int32 *) _degrees, _conn
    cdef FMField[1] _bf, _bfg, _det
    cdef FMField[1] _bfg_dxi, _dx_dxi, _dxi_dx
    cdef FMField[1] _qp, _control_points, _weights
    cdef FMField _cs[3]
    cdef FMField[3] _B, _dB_dxi, _N, _dN_dxi
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

    for ii in range(0, <uint32>dim):
        fmf_alloc(_B + ii, n_qp, 1, n_efuns[ii], 1)
        fmf_alloc(_dB_dxi + ii, n_qp, 1, n_efuns[ii], 1)
        fmf_alloc(_N + ii, 1, 1, n_efuns[ii], 1)
        fmf_alloc(_dN_dxi + ii, 1, 1, n_efuns[ii], 1)

    # Assign to C structures.
    array2fmfield4(_bfg_dxi, bfg_dxi)
    array2fmfield4(_dx_dxi, dx_dxi)
    array2fmfield4(_dxi_dx, dxi_dx)
    array2fmfield2(_control_points, control_points)
    array2fmfield1(_weights, weights)
    for ii in range(<uint32>dim):
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

    # Pre-compute 1D Bernstein basis B, dB/dxi.
    for iqp in range(0, n_qp):
        for ii in range(0, <uint32>dim):
            nf = n_efuns[ii]
            _eval_bernstein_basis(_B + ii, _dB_dxi + ii,
                                  _qp.val[ii], _degrees[ii])
            (_B + ii).val += nf
            (_dB_dxi + ii).val += nf
        _qp.val += dim

    if is_nurbs(weights):
        # Loop over elements.
        _cells = &cells[0]
        for iseq in range(0, n_el):
            ie = _cells[iseq]

            _qp.val = _qp.val0
            for ii in range(0, <uint32>dim):
                (_B + ii).val = (_B + ii).val0
                (_dB_dxi + ii).val = (_dB_dxi + ii).val0

            # Loop over quadrature points.
            for iqp in range(0, n_qp):
                _eval_nurbs_basis_tp(_bf, _bfg, _det,
                                     _bfg_dxi,
                                     _dx_dxi, _dxi_dx,
                                     _B, _dB_dxi, _N, _dN_dxi,
                                     _qp, ie,
                                     _control_points, _weights,
                                     _degrees, dim, _cs, _conn, n_el, n_ep,
                                     1, 1)
                _bf.val += n_efun
                _bfg.val += dim * n_efun
                _det.val += 1
                _qp.val += dim
                for ii in range(0, <uint32>dim):
                    nf = n_efuns[ii]
                    (_B + ii).val += nf
                    (_dB_dxi + ii).val += nf

    else:
        # Loop over elements.
        _cells = &cells[0]
        for iseq in range(0, n_el):
            ie = _cells[iseq]

            _qp.val = _qp.val0
            for ii in range(0, <uint32>dim):
                (_B + ii).val = (_B + ii).val0
                (_dB_dxi + ii).val = (_dB_dxi + ii).val0

            # Loop over quadrature points.
            for iqp in range(0, n_qp):
                _eval_bspline_basis_tp(_bf, _bfg, _det,
                                       _bfg_dxi,
                                       _dx_dxi, _dxi_dx,
                                       _B, _dB_dxi, _N, _dN_dxi,
                                       _qp, ie,
                                       _control_points,
                                       _degrees, dim, _cs, _conn, n_el, n_ep,
                                       1, 1)
                _bf.val += n_efun
                _bfg.val += dim * n_efun
                _det.val += 1
                _qp.val += dim
                for ii in range(0, <uint32>dim):
                    nf = n_efuns[ii]
                    (_B + ii).val += nf
                    (_dB_dxi + ii).val += nf

    for ii in range(0, <uint32>dim):
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
    cdef uint32 ii, ie, n_qp, n_efun, nf, nc, ir, ic
    cdef int32 n_el, n_ep, dim, aux
    cdef uint32 *_cells
    cdef (int32 *) _degrees, _conn, ec
    cdef float64 val
    cdef FMField[1] _bf, _bfg, _det, _vals, _coors
    cdef FMField[1] _bfg_dxi, _dx_dxi, _dxi_dx
    cdef FMField[1] _qp, _variable, _control_points, _weights
    cdef FMField _cs[3]
    cdef FMField[3] _B, _dB_dxi, _N, _dN_dxi
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

    for ii in range(0, <uint32>dim):
        fmf_alloc(_B + ii, n_qp, 1, n_efuns[ii], 1)
        fmf_alloc(_dB_dxi + ii, n_qp, 1, n_efuns[ii], 1)
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
    for ii in range(<uint32>dim):
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

    # Pre-compute 1D Bernstein basis B, dB/dxi.
    for iqp in range(0, n_qp):
        for ii in range(0, <uint32>dim):
            nf = n_efuns[ii]
            _eval_bernstein_basis(_B + ii, _dB_dxi + ii,
                                  _qp.val[ii], _degrees[ii])
            (_B + ii).val += nf
            (_dB_dxi + ii).val += nf
        _qp.val += dim

    if is_nurbs(weights):
        # Loop over elements.
        _cells = &cells[0]
        for iseq in range(0, n_el):
            ie = _cells[iseq]

            ec = _conn + n_ep * ie

            _qp.val = _qp.val0
            for ii in range(0, <uint32>dim):
                (_B + ii).val = (_B + ii).val0
                (_dB_dxi + ii).val = (_dB_dxi + ii).val0

            # Loop over quadrature points.
            for iqp in range(0, n_qp):
                _eval_nurbs_basis_tp(_bf, _bfg, _det,
                                     _bfg_dxi,
                                     _dx_dxi, _dxi_dx,
                                     _B, _dB_dxi, _N, _dN_dxi,
                                     _qp, ie,
                                     _control_points, _weights,
                                     _degrees, dim, _cs, _conn, n_el, n_ep,
                                     1, 1)

                # vals[ii, :] = np.dot(bf, variable[ec])
                for ir in range(0, nc):
                    _vals.val[ir] = 0.0

                    for ic in range(0, n_efun):
                        val = _variable.val[ec[ic] * nc + ir]
                        _vals.val[ir] += _bf.val[ic] * val

                # coors[ii, :] = np.dot(bf, control_points[ec])
                for ir in range(0, <uint32>dim):
                    _coors.val[ir] = 0.0

                    for ic in range(0, n_efun):
                        val = _control_points.val[ec[ic] * dim + ir]
                        _coors.val[ir] += _bf.val[ic] * val

                _vals.val += nc
                _coors.val += dim
                _det.val += 1
                _qp.val += dim
                for ii in range(0, <uint32>dim):
                    nf = n_efuns[ii]
                    (_B + ii).val += nf
                    (_dB_dxi + ii).val += nf

    else:
        # Loop over elements.
        _cells = &cells[0]
        for iseq in range(0, n_el):
            ie = _cells[iseq]

            ec = _conn + n_ep * ie

            _qp.val = _qp.val0
            for ii in range(0, <uint32>dim):
                (_B + ii).val = (_B + ii).val0
                (_dB_dxi + ii).val = (_dB_dxi + ii).val0

            # Loop over quadrature points.
            for iqp in range(0, n_qp):
                _eval_bspline_basis_tp(_bf, _bfg, _det,
                                       _bfg_dxi,
                                       _dx_dxi, _dxi_dx,
                                       _B, _dB_dxi, _N, _dN_dxi,
                                       _qp, ie,
                                       _control_points,
                                       _degrees, dim, _cs, _conn, n_el, n_ep,
                                       1, 1)

                # vals[ii, :] = np.dot(bf, variable[ec])
                for ir in range(0, nc):
                    _vals.val[ir] = 0.0

                    for ic in range(0, n_efun):
                        val = _variable.val[ec[ic] * nc + ir]
                        _vals.val[ir] += _bf.val[ic] * val

                # coors[ii, :] = np.dot(bf, control_points[ec])
                for ir in range(0, <uint32>dim):
                    _coors.val[ir] = 0.0

                    for ic in range(0, n_efun):
                        val = _control_points.val[ec[ic] * dim + ir]
                        _coors.val[ir] += _bf.val[ic] * val

                _vals.val += nc
                _coors.val += dim
                _det.val += 1
                _qp.val += dim
                for ii in range(0, <uint32>dim):
                    nf = n_efuns[ii]
                    (_B + ii).val += nf
                    (_dB_dxi + ii).val += nf

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
    cdef (int32 *) _degrees, _conn, ec
    cdef uint32[3] igrid, shape, iis, n_els
    cdef uint32 *_indices[3]
    cdef uint32 **puaux
    cdef FMField[1] _bf, _bfg, _det, _vals, _out
    cdef FMField[1] _bfg_dxi, _dx_dxi, _dxi_dx
    cdef FMField[1] _rc, _control_points, _weights
    cdef FMField[3] _cs, _ref_coors
    cdef np.ndarray[float64, mode='c', ndim=2] _evals, out
    cdef FMField[3] _B, _dB_dxi, _N, _dN_dxi

    dim = control_points.shape[1]
    n_efuns = degrees + 1
    n_efun = np.prod(n_efuns)

    n_vals = 1
    for ii in range(0, <uint32>dim):
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

    for ii in range(0, <uint32>dim):
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
    array2fmfield2(_vals, _evals)
    for ii in range(<uint32>dim):
        array2fmfield4(_cs + ii, cs[ii])
        n_els[ii] = (_cs + ii).nCell;

    array2pint1(&_degrees, &dim, degrees)
    array2pint2(&_conn, &n_el, &n_ep, conn)

    _out.offset = _out.nAlloc = -1
    _out.nCell = _out.nLev = _out.nRow = 1
    _out.nCol = nc

    _out.val = _out.val0 = &out[0, 0]

    for ii in range(0, <uint32>dim):
        puaux = _indices + ii
        array2puint1(puaux, &uaux, indices[ii])
        array2fmfield1(_ref_coors + ii, ref_coors[ii])

    if is_nurbs(weights):
        for ip in range(0, n_vals):
            _unravel_index(igrid, ip, shape, dim)

            for ii in range(0, <uint32>dim):
                iis[ii] = _indices[ii][igrid[ii]]
                _rc.val[ii] = (_ref_coors + ii).val[igrid[ii]]
            _ravel_multi_index(&ie, iis, n_els, dim)

            _eval_nurbs_basis_tp(_bf, _bfg, _det,
                                 _bfg_dxi,
                                 _dx_dxi, _dxi_dx,
                                 _B, _dB_dxi, _N, _dN_dxi,
                                 _rc, ie,
                                 _control_points, _weights,
                                 _degrees, dim, _cs, _conn, n_el, n_ep,
                                 0, 1)

            # vals[ip, :] = np.dot(bf, variable[ec])
            ec = _conn + n_ep * ie;
            for ir in range(0, nc):
                _out.val[ir] = 0.0

                for ic in range(0, n_efun):
                    _out.val[ir] += _bf.val[ic] * _vals.val[ec[ic] * nc + ir]

            _out.val += nc

    else:
        for ip in range(0, n_vals):
            _unravel_index(igrid, ip, shape, dim)

            for ii in range(0, <uint32>dim):
                iis[ii] = _indices[ii][igrid[ii]]
                _rc.val[ii] = (_ref_coors + ii).val[igrid[ii]]
            _ravel_multi_index(&ie, iis, n_els, dim)

            _eval_bspline_basis_tp(_bf, _bfg, _det,
                                   _bfg_dxi,
                                   _dx_dxi, _dxi_dx,
                                   _B, _dB_dxi, _N, _dN_dxi,
                                   _rc, ie,
                                   _control_points,
                                   _degrees, dim, _cs, _conn, n_el, n_ep,
                                   0, 1)

            # vals[ip, :] = np.dot(bf, variable[ec])
            ec = _conn + n_ep * ie;
            for ir in range(0, nc):
                _out.val[ir] = 0.0

                for ic in range(0, n_efun):
                    _out.val[ir] += _bf.val[ic] * _vals.val[ec[ic] * nc + ir]

            _out.val += nc

    return out
