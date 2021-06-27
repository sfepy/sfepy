/*!
  @par Revision history:
  - 06.02.2001, c
  - 05.03.2003, adopted from rcfem2
*/
#ifndef _FMFIELD_H_
#define _FMFIELD_H_

#include "common.h"
BEGIN_C_DECLS

#define MachEps 1e-16

/*!
  @par Revision history:
  - 06.02.2001, c
  - 27.04.2001
  - 31.03.2003
*/
typedef struct FMField {
  int32 nCell;
  int32 nLev;
  int32 nRow;
  int32 nCol;
  float64 *val0;
  float64 *val;
  int32 nAlloc;
  int32 cellSize;

  int32 offset;
  int32 nColFull;
} FMField;

/*!
  FMField cell pointer manipulation.

  @par Revision history:
  - 27.04.2001, c
  - 30.04.2001
  - 11.07.2002
*/
#define FMF_PtrFirst( obj ) ((obj)->val0)
#define FMF_PtrCurrent( obj ) ((obj)->val)
#define FMF_PtrCell( obj, n ) ((obj)->val0 + (n) * (obj)->cellSize)
#define FMF_PtrCellX1( obj, n ) ((obj)->nCell > 1 \
  ? (obj)->val0 + (n) * (obj)->cellSize : (obj)->val0)
#define FMF_SetFirst( obj ) ((obj)->val = (obj)->val0)
#define FMF_SetCell( obj, n ) ((obj)->val = (obj)->val0 + (n) * (obj)->cellSize)
#define FMF_SetCellX1( obj, n ) do {\
    if ((obj)->nCell > 1) ((obj)->val = (obj)->val0 + (n) * (obj)->cellSize); \
  } while (0)
#define FMF_SetCellNext( obj ) ((obj)->val += (obj)->cellSize)

/*!
  Access to row @ir of level @a il of FMField @a obj.

  @par Revision history:
  - 22.04.2001, c
*/
#define FMF_PtrRowOfLevel( obj, il, ir ) ((obj)->val + \
(obj)->nCol * ((obj)->nRow * (il) + (ir)))

/*!
  Access to level @a il of FMField @a obj.

  @par Revision history:
  - 22.04.2001, c
*/
#define FMF_PtrLevel( obj, il ) ((obj)->val + \
(obj)->nCol * (obj)->nRow * (il))

int32 fmf_alloc( FMField *obj, int32 nCell, int32 nLev,
		 int32 nRow, int32 nCol );
int32 fmf_createAlloc( FMField **p_obj, int32 nCell, int32 nLev,
		       int32 nRow, int32 nCol );
int32 fmf_createAllocInit( FMField **p_obj, int32 nCell, int32 nLev,
			   int32 nRow, int32 nCol, float64 *val );
int32 fmf_createAllocCopy( FMField **p_obj, FMField *obj );
int32 fmf_free( FMField *obj );
int32 fmf_freeDestroy( FMField **p_obj );

int32 fmf_pretend( FMField *obj,
		   int32 nCell, int32 nLev, int32 nRow, int32 nCol,
		   float64 *data );
int32 fmf_pretend_nc( FMField *obj,
                      int32 nCell, int32 nLev, int32 nRow, int32 nCol,
                      float64 *data );
int32 fmfr_pretend( FMField *obj,
		    int32 nLev, int32 nRow, int32 nCol,
		    float64 *data, int32 offset, int32 nColFull );
int32 fmf_set_qp(FMField *qp_obj, int32 iqp, FMField *obj);
int32 fmf_getDim( FMField *obj, int32 *p_nCell, int32 *p_nLev,
		  int32 *p_nRow, int32 *p_nCol );

int32 fmf_fillC( FMField *obj, float64 val );
int32 fmfr_fillC( FMField *obj, float64 val );
int32 fmfc_fillC( FMField *obj, float64 val );
int32 fmfc_fill( FMField *obj, float64 *val );

int32 fmf_mulC( FMField *obj, float64 val );
int32 fmfc_mulC( FMField *obj, float64 val );
int32 fmf_mul( FMField *obj, float64 *val );

int32 fmf_mulAC( FMField *objR, FMField *objA, float64 val );
int32 fmf_mulATC( FMField *objR, FMField *objA, float64 val );
int32 fmf_mulAF( FMField *objR, FMField *objA, float64 *val );
int32 fmf_mulATF( FMField *objR, FMField *objA, float64 *val );
int32 fmf_mulAB_nn( FMField *objR, FMField *objA, FMField *objB );
int32 fmf_mulAB_n1( FMField *objR, FMField *objA, FMField *objB );
int32 fmf_mulAB_1n( FMField *objR, FMField *objA, FMField *objB );
int32 fmf_mulATB_nn( FMField *objR, FMField *objA, FMField *objB );
int32 fmf_mulATB_1n( FMField *objR, FMField *objA, FMField *objB );
int32 fmf_mulABT_nn( FMField *objR, FMField *objA, FMField *objB );
int32 fmf_mulATBT_nn( FMField *objR, FMField *objA, FMField *objB );
int32 fmf_mulATBT_1n( FMField *objR, FMField *objA, FMField *objB );

int32 fmf_addAB_nn( FMField *objR, FMField *objA, FMField *objB );
int32 fmf_subAB_nn( FMField *objR, FMField *objA, FMField *objB );
int32 fmfc_addAB_nn( FMField *objR, FMField *objA, FMField *objB );
int32 fmf_averageCACB( FMField *objR, float64 c1, FMField *objA,
		       float64 c2, FMField *objB );
int32 fmfc_averageCACB( FMField *objR, float64 c1, FMField *objA,
			float64 c2, FMField *objB );
int32 fmfc_normalize( FMField *objR, FMField *objA );

int32 fmf_addAmulF( FMField *objR, FMField *objA, float64 *val );
int32 fmfc_addAmulF( FMField *objR, FMField *objA, float64 *val );
int32 fmf_copyAmulC( FMField *objR, FMField *objA, float64 val );
int32 fmfc_copyAmulF( FMField *objR, FMField *objA, float64 *val );

int32 fmfr_addA_blockNC( FMField *objR, FMField *objA, int32 row, int32 col );
int32 fmfr_addAT_blockNC( FMField *objR, FMField *objA, int32 row, int32 col );

int32 fmf_sumLevelsMulF( FMField *objR, FMField *objA, float64 *val );
int32 fmf_sumLevelsTMulF( FMField *objR, FMField *objA, float64 *val );
int32 fmfr_sumLevelsMulF( FMField *objR, FMField *objA, float64 *val );
int32 fmfr_sumLevelsTMulF( FMField *objR, FMField *objA, float64 *val );

int32 fmf_copy( FMField *objR, FMField *objA );
int32 fmfr_copy( FMField *objR, FMField *objA );
int32 fmfc_copy( FMField *objR, FMField *objA );

int32 fmf_print( FMField *obj, FILE *file, int32 mode );
int32 fmfr_print( FMField *obj, FILE *file, int32 mode );
int32 fmf_save( FMField *obj, const char *fileName, int32 mode );
int32 fmfr_save( FMField *obj, const char *fileName, int32 mode );
int32 fmfc_save( FMField *obj, const char *fileName, int32 mode );

int32 fmf_gMtx2VecDUL3x3( FMField *objR, FMField *objA );
int32 fmf_gMtx2VecDLU3x3( FMField *objR, FMField *objA );

END_C_DECLS

#endif /* Header */
