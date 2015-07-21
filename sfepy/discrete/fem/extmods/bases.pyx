# -*- Mode: Python -*-
"""
Polynomial base functions and related utilities.
"""
cimport cython

cimport numpy as np
import numpy as np

cimport sfepy.discrete.common.extmods._fmfield as _f
from sfepy.discrete.common.extmods._fmfield cimport FMField

from sfepy.discrete.common.extmods.types cimport int32, float64, complex128

cdef extern from 'math.h':
    cdef float64 sqrt(float x)

cdef extern from 'common.h':
    void *pyalloc(size_t size)
    void pyfree(void *pp)
    int Max_i 'Max'(int a, int b)
    double Max_f 'Max'(double a, double b)
    double Min_f 'Min'(double a, double b)

cdef extern from 'lagrange.h':
    ctypedef struct LagrangeContext:
        int32 (*get_xi_dist)(float64 *pdist, FMField *xi,
                             FMField *point, FMField *e_coors,
                             void *_ctx)

        FMField bc[1]
        FMField mtx_i[1]
        FMField base1d[1]
        FMField ref_coors[1]
        int32 *nodes
        int32 n_nod
        int32 n_col
        int32 tdim
        float64 eps
        int32 check_errors
        int32 i_max
        float64 newton_eps
        float64 vmin
        float64 vmax

    void _print_context_lagrange \
         'print_context_lagrange'(LagrangeContext *ctx)

    int32 _get_barycentric_coors \
          'get_barycentric_coors'(FMField *bc, FMField *coors,
                                  LagrangeContext *ctx)

    int32 _get_xi_dist \
          'get_xi_dist'(float64 *pdist, FMField *xi,
                        FMField *point, FMField *e_coors,
                        void *_ctx)

    int32 _get_xi_simplex \
          'get_xi_simplex'(FMField *xi, FMField *dest_point, FMField *e_coors,
                           LagrangeContext *ctx)

    int32 _get_xi_tensor \
          'get_xi_tensor'(FMField *xi, FMField *dest_point, FMField *e_coors,
                          LagrangeContext *ctx)

    int32 _eval_lagrange_simplex \
          'eval_lagrange_simplex'(FMField *out, int32 order, int32 diff,
                                  LagrangeContext *ctx)

    int32 _eval_lagrange_tensor_product \
          'eval_lagrange_tensor_product'(FMField *out, int32 order, int32 diff,
                                         LagrangeContext *ctx)

cdef class CLagrangeContext:

    cdef LagrangeContext *ctx

    def __cinit__(self,
                  np.ndarray[float64, mode='c', ndim=2] bc=None,
                  np.ndarray[float64, mode='c', ndim=2] mtx_i=None,
                  np.ndarray[float64, mode='c', ndim=1] base1d=None,
                  np.ndarray[float64, mode='c', ndim=2] ref_coors=None,
                  np.ndarray[int32, mode='c', ndim=2] nodes=None,
                  int32 tdim=0,
                  float64 eps=1e-15,
                  int32 check_errors=0,
                  int32 i_max=100,
                  float64 newton_eps=1e-8,
                  float64 vmin=0.0,
                  float64 vmax=1.0):
        cdef LagrangeContext *ctx

        ctx = self.ctx = <LagrangeContext *> pyalloc(sizeof(LagrangeContext))

        if ctx is NULL:
            raise MemoryError()

        ctx.get_xi_dist = &_get_xi_dist

        if bc is not None:
            _f.array2fmfield2(ctx.bc, bc)

        if mtx_i is not None:
            _f.array2fmfield2(ctx.mtx_i, mtx_i)

        if base1d is not None:
            _f.array2fmfield3(ctx.base1d, base1d)

        if ref_coors is not None:
            _f.array2fmfield2(ctx.ref_coors, ref_coors)

        if nodes is not None:
            ctx.nodes = &nodes[0, 0]
            ctx.n_nod = nodes.shape[0]
            ctx.n_col = nodes.shape[1]

        ctx.tdim = tdim if tdim > 0 else ref_coors.shape[1]

        ctx.eps = eps
        ctx.check_errors = check_errors

        ctx.i_max = i_max
        ctx.newton_eps = newton_eps

        ctx.vmin = vmin
        ctx.vmax = vmax

    def __dealloc__(self):
        pyfree(self.ctx)

    def __str__(self):
        return 'CLagrangeContext'

    def cprint(self):
        _print_context_lagrange(self.ctx)

@cython.boundscheck(False)
def get_barycentric_coors(np.ndarray[float64, mode='c', ndim=2] coors not None,
                          np.ndarray[float64, mode='c', ndim=2] mtx_i not None,
                          float64 eps=1e-8,
                          int check_errors=False):
    """
    Get barycentric (area in 2D, volume in 3D) coordinates of points.

    Parameters
    ----------
    coors : array
        The coordinates of the points, shape `(n_coor, dim)`.
    mtx_i : array
        The inverse of simplex coordinates matrix, shape `(dim + 1, dim + 1)`.
    eps : float
        The tolerance for snapping out-of-simplex point back to the simplex.
    check_errors : bool
        If True, raise ValueError if a barycentric coordinate is outside
        the snap interval `[-eps, 1 + eps]`.

    Returns
    -------
    bc : array
        The barycentric coordinates, shape `(n_coor, dim + 1)`. Then
        reference element coordinates `xi = dot(bc, ref_coors)`.
    """
    cdef int n_coor = coors.shape[0]
    cdef int n_v = mtx_i.shape[0]
    cdef LagrangeContext ctx[1]
    cdef FMField _coors[1]
    cdef np.ndarray[float64, ndim=2] bc = np.zeros((n_coor, n_v),
                                                   dtype=np.float64)

    ctx.eps = eps
    ctx.check_errors = check_errors
    _f.array2fmfield2(ctx.bc, bc)
    _f.array2fmfield2(ctx.mtx_i, mtx_i)
    _f.array2fmfield2(_coors, coors)

    _get_barycentric_coors(ctx.bc, _coors, ctx)
    return bc

@cython.boundscheck(False)
def eval_lagrange_simplex(np.ndarray[float64, mode='c', ndim=2] coors not None,
                          np.ndarray[float64, mode='c', ndim=2] mtx_i not None,
                          np.ndarray[int32, mode='c', ndim=2] nodes not None,
                          int order, int diff=False,
                          float64 eps=1e-15,
                          int check_errors=True):
    """
    Evaluate Lagrange base polynomials in given points on simplex domain.

    Parameters
    ----------
    coors : array
        The coordinates of the points, shape `(n_coor, dim)`.
    mtx_i : array
        The inverse of simplex coordinates matrix, shape `(dim + 1, dim + 1)`.
    nodes : array
        The description of finite element nodes, shape `(n_nod, dim + 1)`.
    order : int
        The polynomial order.
    diff : bool
        If True, return base function derivatives.
    eps : float
        The tolerance for snapping out-of-simplex point back to the simplex.
    check_errors : bool
        If True, raise ValueError if a barycentric coordinate is outside
        the snap interval `[-eps, 1 + eps]`.

    Returns
    -------
    out : array
        The evaluated base functions, shape `(n_coor, 1 or dim, n_nod)`.
    """
    cdef int bdim
    cdef int n_coor = coors.shape[0]
    cdef int dim = mtx_i.shape[0] - 1
    cdef int n_nod = nodes.shape[0]
    cdef LagrangeContext ctx[1]
    cdef FMField _out[1], _coors[1]
    cdef np.ndarray[float64, ndim=2] bc = np.zeros((n_coor, dim + 1),
                                                   dtype=np.float64)

    assert mtx_i.shape[0] == nodes.shape[1]

    if diff:
        bdim = dim

    else:
        bdim = 1

    cdef np.ndarray[float64, ndim=3] out = np.zeros((n_coor, bdim, n_nod),
                                                    dtype=np.float64)

    ctx.eps = eps
    ctx.check_errors = check_errors
    ctx.nodes = &nodes[0, 0]
    ctx.n_col = nodes.shape[1]
    _f.array2fmfield2(ctx.bc, bc)
    _f.array2fmfield2(ctx.mtx_i, mtx_i)
    _f.array2fmfield3(_out, out)
    _f.array2fmfield2(_coors, coors)

    _get_barycentric_coors(ctx.bc, _coors, ctx)
    _eval_lagrange_simplex(_out, order, diff, ctx)

    return out

@cython.boundscheck(False)
def eval_lagrange_tensor_product(np.ndarray[float64, mode='c', ndim=2]
                                 coors not None,
                                 np.ndarray[float64, mode='c', ndim=2]
                                 mtx_i not None,
                                 np.ndarray[int32, mode='c', ndim=2]
                                 nodes not None,
                                 int order, int diff=False,
                                 float64 eps=1e-15,
                                 int check_errors=True):
    """
    Evaluate Lagrange base polynomials in given points on tensor product
    domain.

    Parameters
    ----------
    coors : array
        The coordinates of the points, shape `(n_coor, dim)`.
    mtx_i : array
        The inverse of 1D simplex coordinates matrix, shape `(2, 2)`.
    nodes : array
        The description of finite element nodes, shape `(n_nod, 2 * dim)`.
    order : int
        The polynomial order.
    diff : bool
        If True, return base function derivatives.
    eps : float
        The tolerance for snapping out-of-simplex point back to the simplex.
    check_errors : bool
        If True, raise ValueError if a barycentric coordinate is outside
        the snap interval `[-eps, 1 + eps]`.

    Returns
    -------
    out : array
        The evaluated base functions, shape `(n_coor, 1 or dim, n_nod)`.
    """
    cdef int ii, idim, im, ic
    cdef int n_coor = coors.shape[0]
    cdef int n_nod = nodes.shape[0]
    cdef int dim = coors.shape[1]
    cdef LagrangeContext ctx[1]
    cdef FMField _out[1], _coors[1]
    cdef int32 *_nodes = &nodes[0, 0]
    cdef np.ndarray[float64, ndim=3] bc = np.zeros((dim, n_coor, 2),
                                                   dtype=np.float64)
    cdef np.ndarray[float64, ndim=3] base1d = np.zeros((n_coor, 1, n_nod),
                                                       dtype=np.float64)
    if diff:
        bdim = dim

    else:
        bdim = 1

    cdef np.ndarray[float64, ndim=3] out = np.zeros((n_coor, bdim, n_nod),
                                                    dtype=np.float64)

    ctx.eps = eps
    ctx.check_errors = check_errors
    ctx.nodes = &nodes[0, 0]
    ctx.n_col = nodes.shape[1]
    _f.array2fmfield2(ctx.mtx_i, mtx_i)
    _f.array2fmfield3(ctx.base1d, base1d)
    _f.fmf_pretend_nc(ctx.bc, dim, 1, n_coor, 2, &bc[0, 0, 0])
    _f.array2fmfield3(_out, out)

    for ii in range(0, dim):
        _f.FMF_SetCell(ctx.bc, ii)
         # slice [:,ii:ii+1]
        _f.fmf_pretend_nc(_coors, 1, 1, coors.shape[0], coors.shape[1],
                          &coors[0, ii])
        _get_barycentric_coors(ctx.bc, _coors, ctx)

    _eval_lagrange_tensor_product(_out, order, diff, ctx)

    return out

from libc.stdio cimport FILE, stdout
@cython.boundscheck(False)
@cython.cdivision(True)
cpdef evaluate_in_rc(np.ndarray[float64, mode='c', ndim=3] out,
                     np.ndarray[float64, mode='c', ndim=2] ref_coors,
                     np.ndarray[int32, mode='c', ndim=1] cells,
                     np.ndarray[int32, mode='c', ndim=1] status,
                     np.ndarray[float64, mode='c', ndim=2] source_vals,
                     np.ndarray[int32, mode='c', ndim=2] conn,
                     np.ndarray[float64, mode='c', ndim=2] eref_coors,
                     np.ndarray[int32, mode='c', ndim=2] nodes,
                     int32 order,
                     np.ndarray[float64, mode='c', ndim=2] mesh_coors,
                     np.ndarray[int32, mode='c', ndim=2] mesh_conn,
                     np.ndarray[int32, mode='c', ndim=2] geo_nodes,
                     int32 diff,
                     np.ndarray[float64, mode='c', ndim=2] mtx_i,
                     float64 qp_eps):
    """
    Evaluate source field DOF values or gradients in the given reference
    element coordinates using the given interpolation.

    1. Evaluate basis functions or gradients of basis functions in the
    reference coordinates. For gradients, tranform the values to the material
    coordinates.
    2. Interpolate source values using the basis functions/gradients.

    Interpolation uses field approximation connectivity.
    """
    cdef int32 ip, ib, ic, iel, n_v, n_ep
    cdef int32 n_cp = 0
    cdef int32 ii, ik
    cdef int32 n_point = ref_coors.shape[0]
    cdef int32 dim = ref_coors.shape[1]
    cdef int32 dpn = out.shape[1]
    cdef int32 bdim = out.shape[2]
    cdef int32 *_cells = &cells[0]
    cdef int32 *_status = &status[0]
    cdef int32 *_conn, *_mesh_conn
    cdef float64 vmin, vmax, aux
    cdef LagrangeContext ctx[1], geo_ctx[1]
    cdef FMField _ref_coors[1], _out[1], bf[1], src[1], cell_coors[1]
    cdef FMField _source_vals[1], _mesh_coors[1]
    cdef FMField mtxMR[1], mtxMRI[1], bfg[1], gbfg[1], gbase1d[1]
    cdef float64 *buf, *buf_b1d_max, *buf_bf_max, *buf_src_max, *buf_bfg
    cdef float64 buf6[6]
    cdef float64 buf9_1[9]
    cdef float64 buf9_2[9]
    cdef float64 buf24_1[24]
    cdef float64 buf24_2[24]

    # Prepare buffers.
    n_ep = conn.shape[1]

    buf = <float64 *> pyalloc(n_ep * (bdim + 1 + dpn + dim) * sizeof(float64))
    buf_b1d_max = buf
    buf_bf_max = buf_b1d_max + n_ep
    buf_src_max = buf_bf_max + n_ep * bdim
    buf_bfg = buf_src_max + n_ep * dpn

    _f.array2fmfield2(_source_vals, source_vals)
    _f.fmf_pretend_nc(_out, n_point, 1, dpn, bdim, &out[0, 0, 0])
    _f.fmf_pretend_nc(_ref_coors, n_point, 1, 1, dim, &ref_coors[0, 0])

    _conn = &conn[0, 0]

    ctx.nodes = &nodes[0, 0]
    ctx.n_col = nodes.shape[1]
    ctx.eps = qp_eps
    ctx.check_errors = 0
    _f.array2fmfield2(ctx.mtx_i, mtx_i)
    _f.array2fmfield2(ctx.ref_coors, eref_coors)

    n_v = ctx.ref_coors.nRow
    vmin = ctx.ref_coors.val[0]
    vmax = ctx.ref_coors.val[dim]

    if n_v == (dim + 1):
        _f.fmf_pretend_nc(ctx.bc, 1, 1, 1, dim + 1, buf6)

    else:
        _f.fmf_pretend_nc(ctx.bc, dim, 1, 1, 2, buf6)
        _f.fmf_pretend_nc(ctx.base1d, 1, 1, 1, n_ep, buf_b1d_max)

    _f.fmf_pretend_nc(bf, 1, 1, bdim, n_ep, buf_bf_max)

    if diff:
        geo_ctx.nodes = &geo_nodes[0, 0]
        geo_ctx.n_col = geo_nodes.shape[1]
        geo_ctx.eps = qp_eps
        geo_ctx.check_errors = 0
        _f.array2fmfield2(geo_ctx.mtx_i, mtx_i)
        _f.array2fmfield2(geo_ctx.ref_coors, eref_coors)

        _mesh_conn = &mesh_conn[0, 0]
        n_cp = mesh_conn.shape[1]

        _f.array2fmfield2(_mesh_coors, mesh_coors)
        _f.fmf_pretend_nc(cell_coors, 1, 1, n_cp, dim, buf24_1)
        _f.fmf_pretend_nc(mtxMR, 1, 1, dim, dim, buf9_1)
        _f.fmf_pretend_nc(mtxMRI, 1, 1, dim, dim, buf9_2)
        _f.fmf_pretend_nc(gbfg, 1, 1, dim, n_cp, buf24_2)
        _f.fmf_pretend_nc(bfg, 1, 1, dim, n_ep, buf_bfg)

        if n_v == (dim + 1):
            # Shares buffer with ctx!
            _f.fmf_pretend_nc(geo_ctx.bc, 1, 1, 1, dim + 1, buf6)

        else:
            # Shares buffer with ctx!
            _f.fmf_pretend_nc(geo_ctx.bc, dim, 1, 1, 2, buf6)
            # Shares buffer with ctx!
            _f.fmf_pretend_nc(geo_ctx.base1d, 1, 1, 1, n_cp, buf_b1d_max)

    else:
        _mesh_conn = NULL

    _f.fmf_pretend_nc(src, 1, 1, dpn, n_ep, buf_src_max)

    # Point (destination coordinate) loop.
    for ip in range(0, n_point):
        _f.FMF_SetCell(_out, ip)
        _f.FMF_SetCell(_ref_coors, ip)

        if _status[ip] <= 1:
            iel = _cells[ip]

            if n_v == (dim + 1):
                _get_barycentric_coors(ctx.bc, _ref_coors, ctx)
                _eval_lagrange_simplex(bf, order, diff, ctx)

            else:
                for ii in range(0, dim):
                    _f.FMF_SetCell(ctx.bc, ii)
                    # slice [:,ii:ii+1]
                    ctx.bc.val[1] = (_ref_coors.val[ii] - vmin) / (vmax - vmin)
                    ctx.bc.val[0] = 1.0 - ctx.bc.val[1]
                _eval_lagrange_tensor_product(bf, order, diff, ctx)

            _f.ele_extractNodalValuesDBD(src, _source_vals,
                                         _conn + n_ep * iel)

            if diff == 0:
                for ic in range(0, dpn):
                    aux = 0.0
                    for ik in range(0, n_ep):
                        aux += bf.val[ik] * src.val[n_ep*ic+ik]

                    _out.val[ic] = aux

            else:
                if n_v == (dim + 1):
                    _eval_lagrange_simplex(gbfg, 1, diff, geo_ctx)

                else:
                    _eval_lagrange_tensor_product(gbfg, 1, diff, geo_ctx)

                _f.ele_extractNodalValuesNBN(cell_coors, _mesh_coors,
                                             _mesh_conn + n_cp * iel)

                # Jacobi matrix from reference to material system.
                _f.fmf_mulATBT_1n(mtxMR, cell_coors, gbfg)
                # Inverse of Jacobi matrix reference to material system.
                _f.geme_invert3x3(mtxMRI, mtxMR)
                # Base function gradient w.r.t. material system.
                _f.fmf_mulATB_nn(bfg, mtxMRI, bf)

                for ib in range(0, bdim):
                    for ic in range(0, dpn):
                        aux = 0.0
                        for ik in range(0, n_ep):
                            aux += bfg.val[n_ep*ib+ik] * src.val[n_ep*ic+ik]

                        _out.val[dpn*ib+ic] = aux

        else:
            _f.fmf_fillC(_out, 0.0)

    pyfree(buf)
