# -*- Mode: Python -*-
# cython: language_level=3
cimport cython

cimport numpy as np
import numpy as np

cimport sfepy.discrete.common.extmods._fmfield as _f
from sfepy.discrete.common.extmods._fmfield cimport FMField

cimport sfepy.discrete.common.extmods.cmesh as cm

from sfepy.discrete.common.extmods.types cimport int32, float64, complex128

cdef extern from 'common.h':
    void *pyalloc(size_t size)
    void pyfree(void *pp)

cdef extern from 'mesh.h':
    ctypedef struct Mesh:
        pass

cdef extern from 'refcoors.h':
    int32 _refc_find_ref_coors_convex \
        'refc_find_ref_coors_convex'(FMField *ref_coors,
                                     int32 *cells, int32 n_cells,
                                     int32 *status, int32 n_status,
                                     FMField *coors,
                                     Mesh *mesh,
                                     FMField *centroids,
                                     FMField *normals0,
                                     FMField *normals1,
                                     int32 *ics, int32 n_ics,
                                     int32 allow_extrapolation,
                                     float64 qp_eps,
                                     float64 close_limit,
                                     void *_ctx)

    int32 _refc_find_ref_coors \
        'refc_find_ref_coors'(FMField *ref_coors,
                              int32 *cells, int32 n_cells,
                              int32 *status, int32 n_status,
                              FMField *coors,
                              Mesh *mesh,
                              int32 *candidates, int32 n_candidates,
                              int32 *offsets, int32 n_offsets,
                              int32 allow_extrapolation,
                              float64 qp_eps,
                              float64 close_limit,
                              void *_ctx)

    ctypedef struct BasisContext:
        int32 (*get_xi_dist)(float64 *pdist, FMField *xi,
                             FMField *point, FMField *e_coors,
                             void *_ctx)
        int32 (*eval_basis)(FMField *out, FMField *coors, int32 diff,
                            void *_ctx)
        int32 iel # >= 0.
        int32 is_dx # 1 => apply reference mapping to gradient.

from libc.stdio cimport FILE, stdout

cdef class CBasisContext:

    cdef void *ctx

@cython.boundscheck(False)
def find_ref_coors_convex(
        np.ndarray[float64, mode='c', ndim=2] ref_coors not None,
        np.ndarray[int32, mode='c', ndim=1] cells not None,
        np.ndarray[int32, mode='c', ndim=1] status not None,
        np.ndarray[float64, mode='c', ndim=2] coors not None,
        cm.CMesh cmesh not None,
        np.ndarray[float64, mode='c', ndim=2] centroids not None,
        np.ndarray[float64, mode='c', ndim=2] normals0 not None,
        np.ndarray[float64, mode='c', ndim=2] normals1,
        np.ndarray[int32, mode='c', ndim=1] ics,
        int allow_extrapolation,
        float64 qp_eps,
        float64 close_limit,
        _ctx
    ):
    cdef int32 n_cells, n_status, n_ics, n_nodes
    cdef (int32 *) _cells, _status, _ics
    cdef CBasisContext ctx = <CBasisContext> _ctx
    cdef FMField[1] _ref_coors, _coors, _centroids, _normals1, _normals0

    _f.array2fmfield2(_ref_coors, ref_coors)
    _f.array2fmfield2(_coors, coors)
    _f.array2fmfield2(_centroids, centroids)
    _f.array2fmfield2(_normals0, normals0)
    if normals1 is not None:
        _f.array2fmfield2(_normals1, normals1)

    _f.array2pint1(&_cells, &n_cells, cells)
    _f.array2pint1(&_status, &n_status, status)
    _f.array2pint1(&_ics, &n_ics, ics)

    _refc_find_ref_coors_convex(_ref_coors,
                                _cells, n_cells,
                                _status, n_status,
                                _coors,
                                cmesh.mesh,
                                _centroids,
                                _normals0,
                                _normals1,
                                _ics, n_ics,
                                allow_extrapolation, qp_eps, close_limit,
                                ctx.ctx)

@cython.boundscheck(False)
def find_ref_coors(
        np.ndarray[float64, mode='c', ndim=2] ref_coors not None,
        np.ndarray[int32, mode='c', ndim=1] cells not None,
        np.ndarray[int32, mode='c', ndim=1] status not None,
        np.ndarray[float64, mode='c', ndim=2] coors not None,
        cm.CMesh cmesh not None,
        np.ndarray[int32, mode='c', ndim=1] candidates,
        np.ndarray[int32, mode='c', ndim=1] offsets,
        int allow_extrapolation,
        float64 qp_eps,
        float64 close_limit,
        _ctx
    ):
    cdef int32 n_cells, n_status, n_candidates, n_offsets, n_nodes
    cdef (int32 *) _cells, _status, _candidates, _offsets
    cdef CBasisContext ctx = <CBasisContext> _ctx
    cdef FMField[1] _ref_coors, _coors

    _f.array2fmfield2(_ref_coors, ref_coors)
    _f.array2fmfield2(_coors, coors)

    _f.array2pint1(&_cells, &n_cells, cells)
    _f.array2pint1(&_status, &n_status, status)
    _f.array2pint1(&_candidates, &n_candidates, candidates)
    _f.array2pint1(&_offsets, &n_offsets, offsets)

    _refc_find_ref_coors(_ref_coors,
                         _cells, n_cells,
                         _status, n_status,
                         _coors,
                         cmesh.mesh,
                         _candidates, n_candidates,
                         _offsets, n_offsets,
                         allow_extrapolation, qp_eps, close_limit,
                         ctx.ctx)

@cython.boundscheck(False)
cpdef evaluate_in_rc(np.ndarray[float64, mode='c', ndim=3] out,
                     np.ndarray[float64, mode='c', ndim=2] ref_coors,
                     np.ndarray[int32, mode='c', ndim=1] cells,
                     np.ndarray[int32, mode='c', ndim=1] status,
                     np.ndarray[float64, mode='c', ndim=2] source_vals,
                     np.ndarray[int32, mode='c', ndim=2] conn,
                     int32 diff, _ctx):
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
    cdef (int32 *) _conn, _mesh_conn
    cdef float64 aux
    cdef CBasisContext __ctx = <CBasisContext> _ctx
    cdef BasisContext *ctx = <BasisContext *> __ctx.ctx
    cdef FMField[1] _ref_coors, _out, bf, src
    cdef FMField _source_vals[1]
    cdef (float64 *) buf, buf_bf_max, buf_src_max

    # Prepare buffers.
    n_ep = conn.shape[1]

    if diff:
        assert bdim == dim

    buf = <float64 *> pyalloc(n_ep * (bdim + dpn) * sizeof(float64))
    buf_bf_max = buf
    buf_src_max = buf_bf_max + n_ep * bdim

    _f.array2fmfield2(_source_vals, source_vals)
    _f.fmf_pretend_nc(_out, n_point, 1, dpn, bdim, &out[0, 0, 0])
    _f.fmf_pretend_nc(_ref_coors, n_point, 1, 1, dim, &ref_coors[0, 0])

    _conn = &conn[0, 0]

    _f.fmf_pretend_nc(bf, 1, 1, bdim, n_ep, buf_bf_max)
    _f.fmf_pretend_nc(src, 1, 1, dpn, n_ep, buf_src_max)

    ctx.is_dx = 1

    # Point (destination coordinate) loop.
    for ip in range(0, n_point):
        _f.FMF_SetCell(_out, ip)
        _f.FMF_SetCell(_ref_coors, ip)

        if _status[ip] <= 1:
            iel = _cells[ip]

            ctx.iel = iel

            ctx.eval_basis(bf, _ref_coors, diff, <void *> ctx)

            _f.ele_extractNodalValuesDBD(src, _source_vals,
                                         _conn + n_ep * iel)

            if diff == 0:
                for ic in range(0, dpn):
                    aux = 0.0
                    for ik in range(0, n_ep):
                        aux += bf.val[ik] * src.val[n_ep*ic+ik]

                    _out.val[ic] = aux

            else:
                for ic in range(0, dpn):
                    for ib in range(0, bdim):
                        aux = 0.0
                        for ik in range(0, n_ep):
                            aux += bf.val[n_ep*ib+ik] * src.val[n_ep*ic+ik]

                        _out.val[bdim*ic+ib] = aux

        else:
            _f.fmf_fillC(_out, 0.0)

    pyfree(buf)
