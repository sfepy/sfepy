# -*- Mode: Python -*-
"""
Polynomial base functions and related utilities.
"""
cimport cython

cimport numpy as np
import numpy as np

cimport _fmfield as _f

from types cimport int32, float64, complex128

cdef extern from 'math.h':
    cdef float64 sqrt(float x)

@cython.boundscheck(False)
cdef _get_barycentric_coors(_f.FMField *bc, _f.FMField *coors,
                            _f.FMField *mtx_i,
                            float64 eps=1e-8, int check_errors=False):
    """
    Get barycentric (area in 2D, volume in 3D) coordinates of points.

    Parameters
    ----------
    bc : FMField
        The barycentric coordinates, shape `(n_coor, dim + 1)`. Then
        reference element coordinates `xi = dot(bc, ref_coors)`.
    coors : FMField
        The coordinates of the points, shape `(n_coor, dim)`.
    mtx_i : FMField
        The inverse of simplex coordinates matrix, shape `(dim + 1, dim + 1)`.
    eps : float
        The tolerance for snapping out-of-simplex point back to the simplex.
    check_errors : bool
        If True, raise ValueError if a barycentric coordinate is outside
        the snap interval `[-eps, 1 + eps]`.
    """
    cdef int32 error
    cdef int32 ir, ic, ii
    cdef int32 n_coor = coors.nRow
    cdef int32 nc = coors.nCol
    cdef int32 n_v = mtx_i.nRow
    cdef int32 dim = n_v - 1
    cdef float64 val

    for ir in range(0, n_coor):
        for ic in range(0, n_v):
            val = 0.0
            for ii in range(0, dim):
                val += mtx_i.val[n_v*ic+ii] * coors.val[nc*ir+ii]
            val += mtx_i.val[n_v*ic+dim]

            error = False
            if val < 0.0:
                if val > -eps:
                    val = 0.0

                else:
                    error = True

            if val > 1.0:
                if val < (1.0 + eps):
                    val = 1.0

                else:
                    error = True

            if check_errors and error:
                msg = 'point %d outside of element!' % ir
                raise ValueError(msg)

            bc.val[n_v*ir+ic] = val

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
    cdef _f.FMField _bc[1], _coors[1], _mtx_i[1]
    cdef np.ndarray[float64, ndim=2] bc = np.zeros((n_coor, n_v),
                                                   dtype=np.float64)

    _f.array2fmfield2(_bc, bc)
    _f.array2fmfield2(_coors, coors)
    _f.array2fmfield2(_mtx_i, mtx_i)
    _get_barycentric_coors(_bc, _coors, _mtx_i, eps, check_errors)
    return bc

@cython.boundscheck(False)
@cython.cdivision(True)
cdef _eval_lagrange_simplex(_f.FMField *out, _f.FMField *bc, _f.FMField *mtx_i,
                            int32 *nodes, int32 n_col,
                            int order, int diff=False):
    """
    Evaluate Lagrange base polynomials in given points on simplex domain.
    """
    cdef int32 n_coor = bc.nRow
    cdef int32 n_v = bc.nCol
    cdef int32 dim = n_v - 1
    cdef int32 n_nod = out.nCol
    cdef int32 ii, ir, ic, i1, i2, inod, n_i1, n_ii
    cdef float64 dval, dd, vv, bci1, bcii
    cdef float64 *pout

    assert n_coor == out.nLev

    if not diff:
        for ic in range(0, n_coor):
            pout = _f.FMF_PtrLevel(out, ic)

            for inod in range(0, n_nod):
                pout[inod] = 1.0
                for i1 in range(0, n_v):
                    n_i1 = nodes[n_col*inod+i1]
                    bci1 = bc.val[n_v*ic+i1]

                    for i2 in range(0, n_i1):
                        pout[inod] *= (order * bci1 - i2) / (i2 + 1.0)

    else:
        for ic in range(0, n_coor):
            pout = _f.FMF_PtrLevel(out, ic)

            for inod in range(0, n_nod):
                pout[inod] = 0.0

                for ii in range(0, n_v):
                    vv = 1.0
                    bcii = bc.val[n_v*ic+ii]

                    for i1 in range(0, n_v):
                        if i1 == ii: continue
                        n_i1 = nodes[n_col*inod+i1]
                        bci1 = bc.val[n_v*ic+i1]

                        for i2 in range(0, n_i1):
                            vv *= (order * bci1 - i2) / (i2 + 1.0)

                    dval = 0.0
                    n_ii = nodes[n_col*inod+ii]
                    for i1 in range(0, n_ii):
                        dd = 1.0

                        for i2 in range(0, n_ii):
                            if i1 == i2: continue

                            dd *= (order * bcii - i2) / (i2 + 1.0)

                        dval += dd * order / (i1 + 1.0)

                    for ir in range(0, dim):
                        pout[n_nod*ir+inod] += vv * dval * mtx_i.val[n_v*ii+ir]

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
    cdef _f.FMField _out[1], _bc[1], _coors[1], _mtx_i[1]
    cdef int32 *_nodes = &nodes[0, 0]
    cdef np.ndarray[float64, ndim=2] bc = np.zeros((n_coor, dim + 1),
                                                   dtype=np.float64)

    assert mtx_i.shape[0] == nodes.shape[1]

    if diff:
        bdim = dim

    else:
        bdim = 1

    cdef np.ndarray[float64, ndim=3] out = np.zeros((n_coor, bdim, n_nod),
                                                    dtype=np.float64)

    _f.array2fmfield3(_out, out)
    _f.array2fmfield2(_bc, bc)
    _f.array2fmfield2(_coors, coors)
    _f.array2fmfield2(_mtx_i, mtx_i)

    _get_barycentric_coors(_bc, _coors, _mtx_i, eps, check_errors)
    _eval_lagrange_simplex(_out, _bc, _mtx_i, _nodes, nodes.shape[1],
                           order, diff)

    return out

@cython.boundscheck(False)
cdef _eval_lagrange_tensor_product(_f.FMField *out, _f.FMField *bc,
                                   _f.FMField *mtx_i, _f.FMField *base1d,
                                   int32 *nodes, int32 n_col,
                                   int order, int diff=False):
    """
    Evaluate Lagrange base polynomials in given points on tensor product
    domain.
    """
    cdef int32 ii, idim, im, ic
    cdef int32 n_coor = out.nLev
    cdef int32 nr = out.nRow
    cdef int32 n_nod = out.nCol
    cdef int32 dim = bc.nCell
    cdef int32 out_size = n_coor * nr * n_nod
    cdef int32 *pnodes

    _f.fmf_fillC(out, 1.0)

    if not diff:
        for ii in range(0, dim):
            # slice [:,2*ii:2*ii+2]
            pnodes = nodes + 2 * ii
            # slice [:,ii:ii+1]
            _f.FMF_SetCell(bc, ii)

            _eval_lagrange_simplex(base1d, bc, mtx_i, pnodes, n_col,
                                   order, diff)

            for im in range(0, out.cellSize):
                out.val[im] *= base1d.val[im]

    else:
        for ii in range(0, dim):
            # slice [:,2*ii:2*ii+2]
            pnodes = nodes + 2 * ii
            # slice [:,ii:ii+1]
            _f.FMF_SetCell(bc, ii)

            for idim in range(0, dim):
                if ii == idim:
                    _eval_lagrange_simplex(base1d, bc, mtx_i, pnodes, n_col,
                                           order, diff)

                else:
                    _eval_lagrange_simplex(base1d, bc, mtx_i, pnodes, n_col,
                                           order, False)

                # slice [:,idim:idim+1,:]
                for im in range(0, n_coor):
                    for ic in range(0, n_nod):
                        out.val[nr*n_nod*im+n_nod*idim+ic] \
                             *= base1d.val[n_nod*im+ic]

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
    cdef _f.FMField _out[1], _bc[1], _coors[1], _mtx_i[1], _base1d[1]
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

    _f.array2fmfield3(_out, out)
    _bc.nAlloc = -1
    _f.fmf_pretend(_bc, dim, 1, n_coor, 2, &bc[0, 0, 0])
    _f.array2fmfield2(_mtx_i, mtx_i)
    _f.array2fmfield3(_base1d, base1d)
    _coors.nAlloc = -1

    for ii in range(0, dim):
        _f.FMF_SetCell(_bc, ii)
        # slice [:,ii:ii+1]
        _f.fmf_pretend(_coors, 1, 1, coors.shape[0], coors.shape[1],
                       &coors[0, ii])
        _get_barycentric_coors(_bc, _coors, _mtx_i, eps, check_errors)

    _eval_lagrange_tensor_product(_out, _bc, _mtx_i, _base1d,
                                  _nodes, nodes.shape[1], order, diff)

    return out

cdef _get_xi_simplex(_f.FMField *xi, _f.FMField *bc, _f.FMField *dest_point,
                     _f.FMField *ref_coors, _f.FMField *e_coors):
    """
    Get reference simplex coordinates `xi` of `dest_point` given spatial
    element coordinates `e_coors` and coordinates of reference simplex
    vertices `ref_coors`. Output also the corresponding barycentric
    coordinates `bc`.
    """
    cdef int32 idim, ii
    cdef int32 n_v = e_coors.nRow
    cdef int32 dim = e_coors.nCol
    cdef _f.FMField mtx[1], mtx_i[1], rhs[1], bc[1]
    cdef float64 buf16[16], buf16_2[16], buf4[4]

    mtx.nAlloc = -1
    mtx_i.nAlloc = -1
    rhs.nAlloc = -1
    bc.nAlloc = -1

    _f.fmf_pretend(mtx, 1, 1, n_v, n_v, buf16)
    _f.fmf_pretend(mtx_i, 1, 1, n_v, n_v, buf16_2)
    _f.fmf_pretend(rhs, 1, 1, n_v, 1, buf4)

    for idim in range(dim):
        for ii in range(n_v):
            mtx.val[n_v*idim+ii] = e_coors.val[dim*ii+idim]
        rhs.val[idim] = dest_point.val[idim]

    for ii in range(n_v):
        mtx.val[n_v*dim+ii] = 1.0
    rhs.val[dim] = 1.0

    if dim == 3:
        _f.geme_invert4x4(mtx_i, mtx)

    else:
        _f.geme_invert3x3(mtx_i, mtx)

    _f.fmf_mulAB_nn(bc, mtx_i, rhs)
    _f.fmf_mulATB_nn(xi, bc, ref_coors)

cdef int _get_xi_tensor(_f.FMField *xi,
                        _f.FMField *dest_point, _f.FMField *e_coors,
                        _f.FMField *mtx_i,
                        _f.FMField *base1d, int32 *nodes, int32 n_col,
                        float64 vmin, float64 vmax,
                        int i_max, float64 newton_eps):
    """
    Get reference tensor product element coordinates using Newton method.

    Uses linear 1D base functions.
    """
    cdef int32 idim, ok = 0
    cdef int32 dim = e_coors.nCol
    cdef _f.FMField bc[1], bf[1], bfg[1], xint[1], res[1], mtx[1], imtx[1]
    cdef float64 buf6[6], buf8[8], buf24[24]
    cdef float64 buf3_1[3], buf3_2[3], buf9_1[9], buf9_2[9]

    _f.fmf_pretend(bc, dim, 1, 1, 2, buf6)
    _f.fmf_pretend(bf, 1, 1, 1, dim**2, buf8)
    _f.fmf_pretend(bfg, 1, 1, dim, dim**2, buf24)
    _f.fmf_pretend(res, 1, 1, 1, dim, buf3_1)
    _f.fmf_pretend(xint, 1, 1, 1, dim, buf3_2)
    _f.fmf_pretend(mtx, 1, 1, dim, dim, buf9_1)
    _f.fmf_pretend(imtx, 1, 1, dim, dim, buf9_2)

    ii = 0
    _f.fmf_fillC(xi, 0.5 * (vmin + vmax))
    while ii < i_max:
        # Base(xi).
        for ii in range(0, dim):
            _f.FMF_SetCell(bc, ii)
            # slice [:,ii:ii+1]
            bc.val[0] = (xi.val[ii] - vmin) / (vmax - vmin)
            bc.val[1] = 1.0 - bc.val[0]

        _eval_lagrange_tensor_product(bf, bc, mtx_i, base1d,
                                      nodes, n_col, 1, 0)

        # X(xi).
        _f.fmf_mulAB_n1(xint, bf, e_coors)
        # Rezidual.
        _f.fmf_subAB_nn(res, dest_point, xint)

        err = 0.0
        for idim in range(0, dim):
            err += res.val[idim] * res.val[idim]
        err = sqrt(err)

        if err < newton_eps:
            break

        # grad Base(xi).
        _eval_lagrange_tensor_product(bfg, bc, mtx_i, base1d,
                                      nodes, n_col, 1, 1)
        # - Matrix.
        _f.fmf_mulAB_n1(mtx, bfg, e_coors)

        _f.geme_invert3x3(imtx, mtx)

        _f.fmf_mulAB_nn(xint, res, imtx)
        _f.fmf_addAB_nn(xi, xi, xint)
        ii += 1

    if ii == i_max:
        ok = 2
        # Newton did not converge.
        # Use centre if nothing better is available, but do not spoil
        # possible d_min (see below).
        _f.fmf_fillC(xi, 0.5 * (vmin + vmax))

    return ok
