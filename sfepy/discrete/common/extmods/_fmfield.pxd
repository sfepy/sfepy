# -*- Mode: Python -*-
# cython: language_level=3
from libc.stdio cimport FILE, stdout

cimport numpy as np
import numpy as np

from sfepy.discrete.common.extmods.types cimport int32, uint32, float64

cdef extern from "fmfield.h":
    ctypedef struct FMField:
        int32 nCell
        int32 nLev
        int32 nRow
        int32 nCol
        float64 *val0
        float64 *val
        int32 nAlloc
        int32 cellSize

        int32 offset
        int32 nColFull

    cdef int32 fmf_alloc(FMField *obj, int32 nCell, int32 nLev,
                         int32 nRow, int32 nCol)
    cdef int32 fmf_createAlloc(FMField **p_obj, int32 nCell, int32 nLev,
                               int32 nRow, int32 nCol)
    cdef int32 fmf_createAllocInit(FMField **p_obj, int32 nCell, int32 nLev,
                                   int32 nRow, int32 nCol, float64 *val)
    cdef int32 fmf_createAllocCopy(FMField **p_obj, FMField *obj)
    cdef int32 fmf_free(FMField *obj)
    cdef int32 fmf_freeDestroy(FMField **p_obj)

    cdef int32 fmf_pretend(FMField *obj,
                           int32 nCell, int32 nLev,
                           int32 nRow, int32 nCol,
                           float64 *data)
    cdef int32 fmf_pretend_nc(FMField *obj,
                              int32 nCell, int32 nLev, int32 nRow, int32 nCol,
                              float64 *data)
    cdef int32 fmfr_pretend(FMField *obj,
                            int32 nLev, int32 nRow, int32 nCol,
                            float64 *data, int32 offset, int32 nColFull)
    cdef int32 fmf_set_qp(FMField *qp_obj, int32 iqp, FMField *obj)

    cdef int32 fmf_print(FMField *obj, FILE *file, int32 mode)

    cdef int32 fmf_fillC(FMField *obj, float64 val)
    cdef int32 fmfr_fillC(FMField *obj, float64 val)
    cdef int32 fmfc_fillC(FMField *obj, float64 val)
    cdef int32 fmfc_fill(FMField *obj, float64 *val)

    cdef int32 fmf_mulC(FMField *obj, float64 val)
    cdef int32 fmfc_mulC(FMField *obj, float64 val)
    cdef int32 fmf_mul(FMField *obj, float64 *val)

    cdef int32 fmf_mulAC(FMField *objR, FMField *objA, float64 val)
    cdef int32 fmf_mulAF(FMField *objR, FMField *objA, float64 *val)
    cdef int32 fmf_mulATF(FMField *objR, FMField *objA, float64 *val)
    cdef int32 fmf_mulAB_nn(FMField *objR, FMField *objA, FMField *objB)
    cdef int32 fmf_mulAB_n1(FMField *objR, FMField *objA, FMField *objB)
    cdef int32 fmf_mulAB_1n(FMField *objR, FMField *objA, FMField *objB)
    cdef int32 fmf_mulATB_nn(FMField *objR, FMField *objA, FMField *objB)
    cdef int32 fmf_mulATB_1n(FMField *objR, FMField *objA, FMField *objB)
    cdef int32 fmf_mulABT_nn(FMField *objR, FMField *objA, FMField *objB)
    cdef int32 fmf_mulATBT_nn(FMField *objR, FMField *objA, FMField *objB)
    cdef int32 fmf_mulATBT_1n(FMField *objR, FMField *objA, FMField *objB)

    cdef int32 fmf_addAB_nn(FMField *objR, FMField *objA, FMField *objB)
    cdef int32 fmf_subAB_nn(FMField *objR, FMField *objA, FMField *objB)
    cdef int32 fmfc_addAB_nn(FMField *objR, FMField *objA, FMField *objB)
    cdef int32 fmf_averageCACB(FMField *objR, float64 c1, FMField *objA,
                               float64 c2, FMField *objB)
    cdef int32 fmfc_averageCACB(FMField *objR, float64 c1, FMField *objA,
                                float64 c2, FMField *objB)
    cdef int32 fmfc_normalize(FMField *objR, FMField *objA)

    cdef int32 fmf_addAmulF(FMField *objR, FMField *objA, float64 *val)
    cdef int32 fmfc_addAmulF(FMField *objR, FMField *objA, float64 *val)
    cdef int32 fmf_copyAmulC(FMField *objR, FMField *objA, float64 val)
    cdef int32 fmfc_copyAmulF(FMField *objR, FMField *objA, float64 *val)

    cdef int32 fmfr_addA_blockNC(FMField *objR, FMField *objA, int32 row, int32 col)
    cdef int32 fmfr_addAT_blockNC(FMField *objR, FMField *objA, int32 row, int32 col)

    cdef int32 fmf_sumLevelsMulF(FMField *objR, FMField *objA, float64 *val)
    cdef int32 fmf_sumLevelsTMulF(FMField *objR, FMField *objA, float64 *val)
    cdef int32 fmfr_sumLevelsMulF(FMField *objR, FMField *objA, float64 *val)
    cdef int32 fmfr_sumLevelsTMulF(FMField *objR, FMField *objA, float64 *val)

    cdef int32 fmf_copy(FMField *objR, FMField *objA)
    cdef int32 fmfr_copy(FMField *objR, FMField *objA)
    cdef int32 fmfc_copy(FMField *objR, FMField *objA)

    cdef int32 fmf_print(FMField *obj, FILE *file, int32 mode)
    cdef int32 fmfr_print(FMField *obj, FILE *file, int32 mode)
    cdef int32 fmf_save(FMField *obj, char *fileName, int32 mode)
    cdef int32 fmfr_save(FMField *obj, char *fileName, int32 mode)
    cdef int32 fmfc_save(FMField *obj, char *fileName, int32 mode)

    cdef int32 fmf_gMtx2VecDUL3x3(FMField *objR, FMField *objA)
    cdef int32 fmf_gMtx2VecDLU3x3(FMField *objR, FMField *objA)

    cdef void FMF_SetCell(FMField *obj, int ii)
    cdef float64 *FMF_PtrLevel(FMField *obj, int ii)

cdef extern from "geommech.h":
    cdef int32 geme_invert3x3(FMField *mtxI, FMField *mtx)
    cdef int32 geme_invert4x4(FMField *mtxI, FMField *mtx)
    cdef int32 geme_tensor2vectorS3(FMField *vec, FMField *mtx)
    cdef int32 geme_det3x3(float64 *det, FMField *mtx)
    cdef int32 geme_trace3x3(float64 *tr, FMField *mtx)
    cdef int32 geme_norm3(float64 *out, FMField *mtx)
    cdef int32 geme_eig3x3(float64 *out, FMField *mtx)
    cdef int32 geme_mulAVSB3(FMField *out, FMField *vs, FMField *_in)

    cdef int32 geme_mulT2ST2S_T4S_ikjl(FMField *t4, FMField *t21, FMField *t22)
    cdef int32 geme_mulT2ST2S_T4S_iljk(FMField *t4, FMField *t21, FMField *t22)
    cdef int32 geme_mulT2S_AA(FMField *R, FMField *A)

    cdef int32 geme_elementVolume(float64 *volume, float64 *jacobian, int32 nQP)

    cdef int32 geme_buildOpOmega_VS3(float64 *pomega, float64 *pdir,
                                     int32 nItem, int32 dim, int32 sym)
    cdef int32 geme_projectToDir(float64 *pdef, float64 *pomega,
                                 float64 *pstrain, int32 nItem, int32 size)

    cdef int32 bf_act(FMField *out, FMField *bf, FMField *_in)
    cdef int32 bf_ract(FMField *out, FMField *bf, FMField *_in)
    cdef int32 bf_actt(FMField *out, FMField *bf, FMField *_in)
    cdef int32 bf_actt_c1(FMField *out, FMField *bf, FMField *_in)
    cdef int32 bf_buildFTF(FMField *ftf, FMField *ftf1)

    cdef int32 geme_invar1(float64 *invar, FMField *mtx)
    cdef int32 geme_invar2(float64 *invar, FMField *mtx)

    cdef void debug_printConn(int32 *conn, int32 nEP)

    cdef int32 ele_extractNodalValuesNBN(FMField *out, FMField *_in,
                                         int32 *conn)

    cdef int32 ele_extractNodalValuesDBD(FMField *out, FMField *_in,
                                         int32 *conn)

cdef int array2fmfield4(FMField *out,
                        np.ndarray[float64, mode='c', ndim=4] arr) except -1
cdef int array2fmfield3(FMField *out,
                        np.ndarray[float64, mode='c', ndim=3] arr) except -1
cdef int array2fmfield2(FMField *out,
                        np.ndarray[float64, mode='c', ndim=2] arr) except -1
cdef int array2fmfield1(FMField *out,
                        np.ndarray[float64, mode='c', ndim=1] arr) except -1
cdef int array2pint2(int32 **out, int32 *n_row, int32 *n_col,
                     np.ndarray[int32, mode='c', ndim=2] arr) except -1
cdef int array2pint1(int32 **out, int32 *n_row,
                     np.ndarray[int32, mode='c', ndim=1] arr) except -1
cdef int array2puint2(uint32 **out, uint32 *n_row, uint32 *n_col,
                      np.ndarray[uint32, mode='c', ndim=2] arr) except -1
cdef int array2puint1(uint32 **out, uint32 *n_row,
                      np.ndarray[uint32, mode='c', ndim=1] arr) except -1
