# -*- Mode: Python -*-
cimport cython

cimport numpy as np
import numpy as np

cimport _fmfield as _f
from _fmfield cimport FMField

cimport cmesh as cm

from types cimport int32, float64, complex128

cdef extern from 'mesh.h':
    ctypedef struct Mesh:
        pass

cdef extern from 'lagrange.h':
    ctypedef struct LagrangeContext:
        FMField *bc
        FMField *mtx_i
        FMField *base1d
        FMField *ref_coors
        int32 *nodes
        int32 n_col
        float64 eps
        int32 check_errors
        int32 i_max
        float64 newton_eps
        float64 vmin
        float64 vmax

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
                                     float64 close_limit, LagrangeContext *ctx)

    int32 _refc_find_ref_coors \
        'refc_find_ref_coors'(FMField *ref_coors,
                              int32 *cells, int32 n_cells,
                              int32 *status, int32 n_status,
                              FMField *coors,
                              Mesh *mesh,
                              int32 *candidates, int32 n_candidates,
                              int32 *offsets, int32 n_offsets,
                              int32 allow_extrapolation,
                              float64 close_limit, LagrangeContext *ctx)

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
        np.ndarray[float64, mode='c', ndim=2] eref_coors,
        np.ndarray[int32, mode='c', ndim=2] nodes,
        np.ndarray[float64, mode='c', ndim=2] mtx_i,
        int allow_extrapolation,
        float64 close_limit, float64 qp_eps,
        int i_max, float64 newton_eps
    ):
    cdef int32 n_cells, n_status, n_ics, n_nodes
    cdef int32 *_cells, *_status, *_ics
    cdef LagrangeContext ctx[1]
    cdef FMField _ref_coors[1], _coors[1], _centroids[1], _normals0[1], \
        _normals1[1]

    ctx.eps = qp_eps
    ctx.i_max = i_max
    ctx.newton_eps = newton_eps

    _f.array2fmfield2(_ref_coors, ref_coors)
    _f.array2fmfield2(_coors, coors)
    _f.array2fmfield2(_centroids, centroids)
    _f.array2fmfield2(_normals0, normals0)
    if normals1 is not None:
        _f.array2fmfield2(_normals1, normals1)
    _f.array2fmfield2(ctx.ref_coors, eref_coors)
    _f.array2fmfield2(ctx.mtx_i, mtx_i)

    _f.array2pint1(&_cells, &n_cells, cells)
    _f.array2pint1(&_status, &n_status, status)
    _f.array2pint1(&_ics, &n_ics, ics)
    _f.array2pint2(&ctx.nodes, &n_nodes, &ctx.n_col, nodes)

    _refc_find_ref_coors_convex(_ref_coors,
                                _cells, n_cells,
                                _status, n_status,
                                _coors,
                                cmesh.mesh,
                                _centroids,
                                _normals0,
                                _normals1,
                                _ics, n_ics,
                                allow_extrapolation, close_limit,
                                ctx)

@cython.boundscheck(False)
def find_ref_coors(np.ndarray[float64, mode='c', ndim=2] ref_coors not None,
                   np.ndarray[int32, mode='c', ndim=1] cells not None,
                   np.ndarray[int32, mode='c', ndim=1] status not None,
                   np.ndarray[float64, mode='c', ndim=2] coors not None,
                   cm.CMesh cmesh not None,
                   np.ndarray[int32, mode='c', ndim=1] candidates,
                   np.ndarray[int32, mode='c', ndim=1] offsets,
                   np.ndarray[float64, mode='c', ndim=2] eref_coors,
                   np.ndarray[int32, mode='c', ndim=2] nodes,
                   np.ndarray[float64, mode='c', ndim=2] mtx_i,
                   int allow_extrapolation,
                   float64 close_limit, float64 qp_eps,
                   int i_max, float64 newton_eps):
    cdef int32 n_cells, n_status, n_candidates, n_offsets, n_nodes
    cdef int32 *_cells, *_status, *_candidates, *_offsets
    cdef LagrangeContext ctx[1]
    cdef FMField _ref_coors[1], _coors[1]

    ctx.eps = qp_eps
    ctx.i_max = i_max
    ctx.newton_eps = newton_eps

    _f.array2fmfield2(_ref_coors, ref_coors)
    _f.array2fmfield2(_coors, coors)
    _f.array2fmfield2(ctx.ref_coors, eref_coors)
    _f.array2fmfield2(ctx.mtx_i, mtx_i)

    ctx.vmin = ctx.ref_coors.val[0]
    ctx.vmax = ctx.ref_coors.val[ref_coors.shape[1]]

    _f.array2pint1(&_cells, &n_cells, cells)
    _f.array2pint1(&_status, &n_status, status)
    _f.array2pint1(&_candidates, &n_candidates, candidates)
    _f.array2pint1(&_offsets, &n_offsets, offsets)
    _f.array2pint2(&ctx.nodes, &n_nodes, &ctx.n_col, nodes)

    _refc_find_ref_coors(_ref_coors,
                         _cells, n_cells,
                         _status, n_status,
                         _coors,
                         cmesh.mesh,
                         _candidates, n_candidates,
                         _offsets, n_offsets,
                         allow_extrapolation, close_limit,
                         ctx)
