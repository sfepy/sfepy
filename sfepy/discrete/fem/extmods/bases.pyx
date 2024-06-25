# -*- Mode: Python -*-
# cython: language_level=3
"""
Polynomial base functions and related utilities.
"""
cimport cython

cimport numpy as np
import numpy as np

cimport sfepy.discrete.common.extmods._fmfield as _f
from sfepy.discrete.common.extmods._fmfield cimport FMField

from sfepy.discrete.common.extmods.types cimport int32, float64, complex128

cdef extern from 'common.h':
    void *pyalloc(size_t size)
    void pyfree(void *pp)

cdef extern from 'lagrange.h':
    ctypedef struct LagrangeContext:
        int32 (*get_xi_dist)(float64 *pdist, FMField *xi,
                             FMField *point, FMField *e_coors,
                             void *_ctx)
        int32 (*eval_basis)(FMField *out, FMField *coors, int32 diff,
                            void *_ctx)
        int32 iel # >= 0.
        int32 is_dx # 1 => apply reference mapping to gradient.
        FMField e_coors_max[1] # Buffer for coordinates of element nodes.

        LagrangeContext *geo_ctx

        int32 order
        int32 is_bubble
        int32 tdim
        int32 *nodes
        int32 n_nod
        int32 n_col

        FMField ref_coors[1]
        float64 vmin
        float64 vmax

        FMField mesh_coors[1]
        int32 *mesh_conn
        int32 n_cell
        int32 n_cp

        FMField mtx_i[1]

        FMField *bc
        FMField base1d[1]
        FMField mbfg[1]

        float64 eps
        int32 check_errors
        int32 i_max
        float64 newton_eps

    void _print_context_lagrange \
         'print_context_lagrange'(LagrangeContext *ctx)

    int32 _get_xi_dist \
          'get_xi_dist'(float64 *pdist, FMField *xi,
                        FMField *point, FMField *e_coors,
                        void *_ctx)

    int32 _eval_basis_lagrange \
          'eval_basis_lagrange'(FMField *out, FMField *coors, int32 diff,
                                void *_ctx)

cdef class CLagrangeContext:

    cdef LagrangeContext *ctx

    # Store arrays to prevent their deallocation in Python.
    cdef readonly CLagrangeContext _geo_ctx
    cdef readonly np.ndarray mesh_coors
    cdef readonly np.ndarray mesh_conn
    cdef readonly np.ndarray e_coors_max # Auxiliary buffer.
    cdef readonly np.ndarray base1d # Auxiliary buffer.
    cdef readonly np.ndarray mbfg # Auxiliary buffer.

    property is_bubble:

        def __get__(self):
            return self.ctx.is_bubble

        def __set__(self, int32 is_bubble):
            self.ctx.is_bubble = is_bubble

    property iel:

        def __get__(self):
            return self.ctx.iel

        def __set__(self, int32 iel):
            assert iel < self.ctx.n_cell
            self.ctx.iel = iel

    property geo_ctx:

        def __set__(self, _ctx):
            cdef CLagrangeContext __ctx = <CLagrangeContext> _ctx
            cdef LagrangeContext *ctx = <LagrangeContext *> __ctx.ctx

            self._geo_ctx = __ctx
            self.ctx.geo_ctx = ctx

    def __cinit__(self,
                  int32 order=1,
                  int32 is_bubble=0,
                  int32 tdim=0,
                  np.ndarray[int32, mode='c', ndim=2] nodes=None,
                  np.ndarray[float64, mode='c', ndim=2] ref_coors=None,
                  np.ndarray mesh_coors=None,
                  np.ndarray mesh_conn=None,
                  np.ndarray[float64, mode='c', ndim=2] mtx_i=None,
                  float64 eps=1e-15,
                  int32 check_errors=0,
                  int32 i_max=100,
                  float64 newton_eps=1e-8):
        cdef LagrangeContext *ctx
        cdef np.ndarray[float64, mode='c', ndim=2] _mesh_coors
        cdef np.ndarray[int32, mode='c', ndim=2] _mesh_conn
        cdef np.ndarray[float64, mode='c', ndim=2] _e_coors_max
        cdef np.ndarray[float64, mode='c', ndim=1] _base1d
        cdef np.ndarray[float64, mode='c', ndim=2] _mbfg

        ctx = self.ctx = <LagrangeContext *> pyalloc(sizeof(LagrangeContext))

        if ctx is NULL:
            raise MemoryError()

        ctx.get_xi_dist = &_get_xi_dist
        ctx.eval_basis = &_eval_basis_lagrange
        ctx.iel = 0
        ctx.is_dx = 0

        ctx.order = order
        ctx.is_bubble = is_bubble

        ctx.tdim = tdim

        _e_coors_max = self.e_coors_max = np.zeros((8, 3), dtype=np.float64)
        _f.array2fmfield2(ctx.e_coors_max, _e_coors_max)

        if nodes is not None:
            ctx.nodes = &nodes[0, 0]
            ctx.n_nod = nodes.shape[0]
            ctx.n_col = nodes.shape[1]

            _base1d = self.base1d = np.zeros((ctx.n_nod,), dtype=np.float64)
            _f.array2fmfield1(ctx.base1d, _base1d)

            _mbfg = self.mbfg = np.zeros((ref_coors.shape[1], ctx.n_nod),
                                         dtype=np.float64)
            _f.array2fmfield2(ctx.mbfg, _mbfg)

        else:
            raise ValueError('nodes argument is required!')

        if ref_coors is not None:
            _f.array2fmfield2(ctx.ref_coors, ref_coors)

            ctx.vmin = ref_coors[0, 0]
            ctx.vmax = ref_coors[1, 0] if tdim > 0 else ctx.vmin

        else:
            raise ValueError('ref_coors argument is required!')

        if mesh_coors is not None:
            _mesh_coors = self.mesh_coors = mesh_coors
            _f.array2fmfield2(ctx.mesh_coors, _mesh_coors)

        else:
            _f.fmf_pretend_nc(ctx.mesh_coors, 0, 0, 0, 0, NULL)

        if mesh_conn is not None:
            _mesh_conn = self.mesh_conn = mesh_conn

            ctx.mesh_conn = &_mesh_conn[0, 0]
            ctx.n_cell = mesh_conn.shape[0]
            ctx.n_cp = mesh_conn.shape[1]

        else:
            ctx.mesh_conn = NULL
            ctx.n_cell = ctx.n_cp = 0

        if mtx_i is not None:
            _f.array2fmfield2(ctx.mtx_i, mtx_i)

        else:
            raise ValueError('mtx_i argument is required!')

        ctx.eps = eps
        ctx.check_errors = check_errors

        ctx.i_max = i_max
        ctx.newton_eps = newton_eps

    def __dealloc__(self):
        pyfree(self.ctx)

    def __str__(self):
        return 'CLagrangeContext'

    def cprint(self):
        _print_context_lagrange(self.ctx)

    def evaluate(self, np.ndarray[float64, mode='c', ndim=2] coors not None,
                 int32 diff=False,
                 float64 eps=1e-15,
                 int32 check_errors=True):
        cdef int32 n_coor = coors.shape[0]
        cdef int32 n_nod = self.ctx.n_nod
        cdef int32 dim = coors.shape[1]
        cdef int32 bdim
        cdef FMField[1] _out, _coors

        ctx = self.ctx

        ctx.check_errors = check_errors
        ctx.eps = eps

        if diff:
            bdim = dim

        else:
            bdim = 1

        cdef np.ndarray[float64, ndim=3] out = np.zeros((n_coor, bdim, n_nod),
                                                        dtype=np.float64)

        _f.array2fmfield3(_out, out)
        _f.array2fmfield2(_coors, coors)

        self.ctx.eval_basis(_out, _coors, diff, ctx)

        return out
