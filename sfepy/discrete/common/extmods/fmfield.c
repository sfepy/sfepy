#include "fmfield.h"

#undef __FUNC__
#define __FUNC__ "fmf_alloc"
/*!
  @par Revision history:
  - 27.04.2001, c
  - 03.05.2001
  - 06.03.2003, adopted from rcfem2
  - 31.03.2003
*/
int32 fmf_alloc( FMField *obj, int32 nCell, int32 nLev,
		 int32 nRow, int32 nCol )
{
  obj->nCell = nCell;
  obj->nLev = nLev;
  obj->nRow = nRow;
  obj->nCol = nCol;
  obj->cellSize = obj->nLev * obj->nRow * obj->nCol;
  obj->nAlloc = obj->nCell * obj->cellSize;
  obj->val0 = obj->val = alloc_mem( float64, obj->nAlloc );

  obj->nColFull = obj->nCol;
  obj->offset = 0;

  return( RET_OK );
}


/*!
  @par Revision history:
  - 06.02.2001, c
  - 27.04.2001
  - 06.03.2003, adopted from rcfem2
*/
int32 fmf_createAlloc( FMField **p_obj, int32 nCell, int32 nLev,
		       int32 nRow, int32 nCol )
{
  *p_obj = alloc_mem( FMField, 1 );
  fmf_alloc( *p_obj, nCell, nLev, nRow, nCol );

  return( RET_OK );
}

/*!
  @par Revision history:
  - 13.03.2003, c
*/
int32 fmf_createAllocInit( FMField **p_obj, int32 nCell, int32 nLev,
			   int32 nRow, int32 nCol, float64 *val )
{
  fmf_createAlloc( p_obj, nCell, nLev, nRow, nCol );
  memcpy( (*p_obj)->val0, val, (*p_obj)->nAlloc * sizeof( float64 ) );

  return( RET_OK );
}

/*!
  @par Revision history:
  - 14.03.2003, c
*/
int32 fmf_createAllocCopy( FMField **p_obj, FMField *obj )
{
  fmf_createAllocInit( p_obj, obj->nCell, obj->nLev, obj->nRow, obj->nCol,
		       obj->val0 );

  return( RET_OK );
}

/*!
  @par Revision history:
  - 06.02.2001, c
  - 30.04.2001
  - 27.05.2001
  - 06.03.2003, adopted from rcfem2
*/
int32 fmf_free( FMField *obj )
{
  if (obj == 0) return( RET_OK );

  if (obj->nAlloc >= 0) {
    free_mem( obj->val0 );
  } else {
    errput( ErrHead "FMField was pretended\n" );
  }

  return( RET_OK );
}

/*!
  @par Revision history:
  - 06.02.2001, c
  - 06.03.2003, adopted from rcfem2
  - 20.03.2003
*/
int32 fmf_freeDestroy( FMField **p_obj )
{
  if ((*p_obj) == 0) return( RET_OK );

  fmf_free( *p_obj );
  free_mem( *p_obj );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_pretend"
/*!
  @par Revision history:
  - 06.02.2001, c
  - 05.03.2003, adopted from rcfem2
  - 31.03.2003
*/
int32 fmf_pretend( FMField *obj,
		   int32 nCell, int32 nLev, int32 nRow, int32 nCol,
		   float64 *data )
{
  if (obj->nAlloc >= 0) {
    errput( ErrHead "ERR_Overwrite\n" );
  }
  
  obj->nCell = nCell;
  obj->nLev = nLev;
  obj->nRow = nRow;
  obj->nCol = nCol;
  obj->val =  obj->val0 = data;
  obj->nAlloc = -1;
  obj->offset = 0;
  obj->nColFull = obj->nCol;
  obj->cellSize = obj->nLev * obj->nRow * obj->nCol;

  return( RET_OK );
}

/*
  No `(obj->nAlloc >= 0)` check - assumes `obj` to be a FMField with
  unallocated data meant to point to some buffer.
*/
int32 fmf_pretend_nc( FMField *obj,
                      int32 nCell, int32 nLev, int32 nRow, int32 nCol,
                      float64 *data )
{
  obj->nCell = nCell;
  obj->nLev = nLev;
  obj->nRow = nRow;
  obj->nCol = nCol;
  obj->val =  obj->val0 = data;
  obj->nAlloc = -1;
  obj->offset = 0;
  obj->nColFull = obj->nCol;
  obj->cellSize = obj->nLev * obj->nRow * obj->nCol;

  return( RET_OK );
}


#undef __FUNC__
#define __FUNC__ "fmfr_pretend"
/*!
  @par Revision history:
  - 31.03.2003, c
  - 04.04.2003
*/
int32 fmfr_pretend( FMField *obj,
		    int32 nLev, int32 nRow, int32 nCol,
		    float64 *data, int32 offset, int32 nColFull )
{
  if (obj->nAlloc >= 0) {
    errput( ErrHead "ERR_Overwrite\n" );
  }
  
  obj->nCell = 1;
  obj->nLev = nLev;
  obj->nRow = nRow;
  obj->nCol = nCol;
  obj->val = obj->val0 = data;
  obj->nAlloc = -1;
  obj->offset = offset;
  obj->nColFull = nColFull;
  obj->cellSize = obj->nLev * obj->nRow * obj->nCol;

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_set_qp"
int32 fmf_set_qp(FMField *qp_obj, int32 iqp, FMField *obj)
{
#ifdef DEBUG_FMF
  if ((qp_obj->nRow != obj->nRow) || (qp_obj->nCol != obj->nCol)
      || (qp_obj->nLev != 1)) {
    errput(ErrHead "ERR_BadMatch: (1 %d %d) == (%d %d %d)|iqp\n",
           qp_obj->nLev, qp_obj->nRow, qp_obj->nCol,
           obj->nLev, obj->nRow, obj->nCol);
  }
#endif

  qp_obj->val = obj->val + obj->nRow * obj->nCol * iqp;

  return(RET_OK);
}

#undef __FUNC__
#define __FUNC__ "fmf_getDim"
/*!
  @par Revision history:
  - 16.04.2003, c
*/
int32 fmf_getDim( FMField *obj, int32 *p_nCell, int32 *p_nLev,
		  int32 *p_nRow, int32 *p_nCol )
{
  *p_nCell = obj->nCell;
  *p_nLev = obj->nLev;
  *p_nRow = obj->nRow;
  *p_nCol = obj->nCol;

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_fillC"
/*!
  @par Revision history:
  - 06.02.2001, c
  - 05.03.2003, adopted from rcfem2
*/
int32 fmf_fillC( FMField *obj, float64 val )
{
  int32 i;

  for (i = 0; i < (obj->nLev * obj->nRow * obj->nCol); i++) {
    obj->val[i] = val;
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfr_fillC"
/*!
  @par Revision history:
  - 07.04.2003, c
  - 08.04.2003
*/
int32 fmfr_fillC( FMField *obj, float64 val )
{
  int32 il, ir, ic, wr, hr;
  float64 *pr;

  wr = obj->nColFull;
  hr = obj->nRow;
  for (il = 0; il < obj->nLev; il++) {
    pr = obj->val + obj->offset + wr * hr * il;
    for (ir = 0; ir < hr; ir++) {
      for (ic = 0; ic < obj->nCol; ic++) {
/*  	output( "%d %d %d %d (%d) %d %d\n", */
/*  		il, ir, ic, wr, obj->nCol, hr, pr - obj->val ); */
	pr[wr*ir+ic] = val;
      }
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfc_fillC"
/*!
  @par Revision history:
  - 08.10.2001, c
  - 05.03.2003, adopted from rcfem2
*/
int32 fmfc_fillC( FMField *obj, float64 val )
{
  int32 i;

  for (i = 0; i < (obj->nCell * obj->nLev * obj->nRow * obj->nCol); i++) {
    obj->val0[i] = val;
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfc_fill"
/*!
  @par Revision history:
  - 21.09.2003, c
*/
int32 fmfc_fill( FMField *obj, float64 *val )
{
  memcpy( obj->val0, val,
	  obj->nCell * obj->nLev * obj->nRow * obj->nCol
	  * sizeof( float64 ) );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_mulC"
/*!
  @par Revision history:
  - 06.02.2001, c
  - 05.03.2003, adopted from rcfem2
*/
int32 fmf_mulC( FMField *obj, float64 val )
{
  int32 i;

  for (i = 0; i < (obj->nLev * obj->nRow * obj->nCol); i++) {
    obj->val[i] *= val;
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfc_mulC"
/*!
  @par Revision history:
  - 03.08.2006, c
*/
int32 fmfc_mulC( FMField *obj, float64 val )
{
  int32 i;

  for (i = 0; i < (obj->nCell * obj->nLev * obj->nRow * obj->nCol); i++) {
    obj->val0[i] *= val;
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_mul"
/*!
  @par Revision history:
  - 06.02.2001, c
  - 05.03.2003, adopted from rcfem2
*/
int32 fmf_mul( FMField *obj, float64 *val )
{
  int32 i, il;
  float64 *pr;

  for (il = 0; il < (obj->nLev); il++) {
    pr = obj->val + obj->nCol * obj->nRow * il;
    for (i = 0; i < (obj->nRow * obj->nCol); i++) {
      pr[i] *= val[il];
    }
  }

  return( RET_OK );
}


#undef __FUNC__
#define __FUNC__ "fmf_mulAC"
/*!
  objR = const * objA

  @par Revision history:
  - 12.05.2005, c
*/
int32 fmf_mulAC( FMField *objR, FMField *objA, float64 val )
{
  int32 i, il;
  float64 *pr, *pa;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objA->nCol)
      || (objR->nLev != objA->nLev)) {
    errput( ErrHead "ERR_BadMatch\n" );
  }
#endif

  for (il = 0; il < (objR->nLev); il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pa = objA->val + objA->nCol * objA->nRow * il;
    for (i = 0; i < (objR->nRow * objR->nCol); i++) {
      pr[i] = pa[i] * val;
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_mulATC"
/*!
  objR = const * objA^T
*/
int32 fmf_mulATC( FMField *objR, FMField *objA, float64 val )
{
  int32 ir, ic, il;
  int32 wa;
  float64 *pr, *pa;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nCol) || (objR->nCol != objA->nRow)
      || (objR->nLev != objA->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) == (%d %d %d)^T * (1)\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol );
  }
#endif

  wa = objA->nCol;
  for (il = 0; il < (objR->nLev); il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pa = objA->val + objA->nCol * objA->nRow * il;
    for (ir = 0; ir < objR->nRow; ir++) {
      for (ic = 0; ic < objR->nCol; ic++) {
	pr[ic] = pa[wa*ic+ir] * val;
      }
      pr += objR->nCol;
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_mulAF"
/*!
  objR = val * objA

  @par Revision history:
  - 23.02.2001, c
  - 05.03.2003, adopted from rcfem2
*/
int32 fmf_mulAF( FMField *objR, FMField *objA, float64 *val )
{
  int32 i, il;
  float64 *pr, *pa;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objA->nCol)
      || (objR->nLev != objA->nLev)) {
    errput( ErrHead "ERR_BadMatch\n" );
  }
#endif

  for (il = 0; il < (objR->nLev); il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pa = objA->val + objA->nCol * objA->nRow * il;
    for (i = 0; i < (objR->nRow * objR->nCol); i++) {
      pr[i] = pa[i] * val[il];
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_mulATF"
/*!
  objR = val * objA^T.
  Can be used instead of transposition.

  @par Revision history:
  - 24.04.2001, c
  - 05.03.2003, adopted from rcfem2
  - 07.04.2003
*/
int32 fmf_mulATF( FMField *objR, FMField *objA, float64 *val )
{
  int32 ir, ic, il;
  int32 wa;
  float64 *pr, *pa;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nCol) || (objR->nCol != objA->nRow)
      || (objR->nLev != objA->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) == (%d %d %d)^T * (?, 1)\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol );
  }
#endif

  wa = objA->nCol;
  for (il = 0; il < (objR->nLev); il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pa = objA->val + objA->nCol * objA->nRow * il;
    for (ir = 0; ir < objR->nRow; ir++) {
      for (ic = 0; ic < objR->nCol; ic++) {
	pr[ic] = pa[wa*ic+ir] * val[il];
      }
      pr += objR->nCol;
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_mulAB_nn"
/*!
  objR = objA * objB

  @par Revision history:
  - 06.02.2001, c
  - 23.09.2002
  - 05.03.2003, adopted from rcfem2
*/
int32 fmf_mulAB_nn( FMField *objR, FMField *objA, FMField *objB )
{
  int32 i, j, k, il, wr, wa, wb;
  float64 *pr, *pa, *pb;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objB->nCol)
      || (objA->nCol != objB->nRow) || (objR->nLev != objA->nLev)
      || (objR->nLev != objB->nLev) ) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) == (%d %d %d) * (%d %d %d)\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol,
	    objB->nLev, objB->nRow, objB->nCol );
  }
#endif

  wr = objR->nCol;
  wa = objA->nCol;
  wb = objB->nCol;
  for (il = 0; il < objR->nLev; il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pa = objA->val + objA->nCol * objA->nRow * il;
    pb = objB->val + objB->nCol * objB->nRow * il;
    for (i = 0; i < objR->nRow; i++) {
      for (j = 0; j < objR->nCol; j++) {
	pr[wr*i+j] = 0.0;
	for (k = 0; k < objA->nCol; k++) {
	  pr[wr*i+j] += pa[wa*i+k] * pb[wb*k+j];
	}
      }
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_mulAB_n1"
/*!
  objR = objA * objB, @a objB has exactly 1 level, which is used to
  multiply all levels of @a objA.

  @par Revision history:
  - 20.04.2001, c
  - 23.09.2002
  - 05.03.2003, adopted from rcfem2
*/
int32 fmf_mulAB_n1( FMField *objR, FMField *objA, FMField *objB )
{
  int32 i, j, k, il, wr, wa, wb;
  float64 *pr, *pa, *pb;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objB->nCol)
      || (objA->nCol != objB->nRow) || (objR->nLev != objA->nLev)
      || (objB->nLev != 1) ) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) = (%d %d %d) * (%d(1) %d %d)\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol,
	    objB->nLev, objB->nRow, objB->nCol );
  }
#endif

  wr = objR->nCol;
  wa = objA->nCol;
  wb = objB->nCol;
  pb = objB->val;
  for (il = 0; il < objR->nLev; il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pa = objA->val + objA->nCol * objA->nRow * il;
    for (i = 0; i < objR->nRow; i++) {
      for (j = 0; j < objR->nCol; j++) {
	pr[wr*i+j] = 0.0;
	for (k = 0; k < objA->nCol; k++) {
	  pr[wr*i+j] += pa[wa*i+k] * pb[wb*k+j];
	}
      }
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_mulAB_1n"
/*!
  objR = objA * objB, @a objA has exactly 1 level, which is used to
  multiply all levels of @a objB.
*/
int32 fmf_mulAB_1n( FMField *objR, FMField *objA, FMField *objB )
{
  int32 i, j, k, il, wr, wa, wb;
  float64 *pr, *pa, *pb;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objB->nCol)
      || (objA->nCol != objB->nRow) || (objR->nLev != objB->nLev)
      || (objA->nLev != 1) ) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) = (%d %d %d) * (%d(1) %d %d)\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol,
	    objB->nLev, objB->nRow, objB->nCol );
  }
#endif

  wr = objR->nCol;
  wa = objA->nCol;
  wb = objB->nCol;
  pa = objA->val;
  for (il = 0; il < objR->nLev; il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pb = objB->val + objB->nCol * objB->nRow * il;
    for (i = 0; i < objR->nRow; i++) {
      for (j = 0; j < objR->nCol; j++) {
	pr[wr*i+j] = 0.0;
	for (k = 0; k < objA->nCol; k++) {
	  pr[wr*i+j] += pa[wa*i+k] * pb[wb*k+j];
	}
      }
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_mulATB_nn"
/*!
  objR = objA^T * objB

  @par Revision history:
  - 06.02.2001, c
  - 23.09.2002
  - 05.03.2003, adopted from rcfem2
*/
int32 fmf_mulATB_nn( FMField *objR, FMField *objA, FMField *objB )
{
  int32 i, j, k, il, wr, wa, wb;
  float64 *pr, *pa, *pb;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nCol) || (objR->nCol != objB->nCol)
      || (objA->nRow != objB->nRow)
      || (objR->nLev != objA->nLev) || (objR->nLev != objB->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) = (%d %d %d)^T * (%d %d %d)\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol,
	    objB->nLev, objB->nRow, objB->nCol );
  }
#endif

  wr = objR->nCol;
  wa = objA->nCol;
  wb = objB->nCol;
  for (il = 0; il < objR->nLev; il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pa = objA->val + objA->nCol * objA->nRow * il;
    pb = objB->val + objB->nCol * objB->nRow * il;
    for (i = 0; i < objR->nRow; i++) {
      for (j = 0; j < objR->nCol; j++) {
	pr[wr*i+j] = 0.0;
	for (k = 0; k < objA->nRow; k++) {
	  pr[wr*i+j] += pa[wa*k+i] * pb[wb*k+j];
	}
      }
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_mulATB_1n"
/*!
  objR = objA^T * objB

  @par Revision history:
  - 20.12.2005, c
*/
int32 fmf_mulATB_1n( FMField *objR, FMField *objA, FMField *objB )
{
  int32 i, j, k, il, wr, wa, wb;
  float64 *pr, *pa, *pb;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nCol) || (objR->nCol != objB->nCol)
      || (objA->nRow != objB->nRow)
      || (1 != objA->nLev) || (objR->nLev != objB->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) = (%d %d %d)^T * (%d %d %d)\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol,
	    objB->nLev, objB->nRow, objB->nCol );
  }
#endif

  wr = objR->nCol;
  wa = objA->nCol;
  wb = objB->nCol;
  pa = objA->val;
  for (il = 0; il < objR->nLev; il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pb = objB->val + objB->nCol * objB->nRow * il;
    for (i = 0; i < objR->nRow; i++) {
      for (j = 0; j < objR->nCol; j++) {
	pr[wr*i+j] = 0.0;
	for (k = 0; k < objA->nRow; k++) {
	  pr[wr*i+j] += pa[wa*k+i] * pb[wb*k+j];
	}
      }
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_mulABT_nn"
/*!
  objR = objA * objB^T

  @par Revision history:
  - 06.02.2001, c
  - 05.03.2003, adopted from rcfem2
*/
int32 fmf_mulABT_nn( FMField *objR, FMField *objA, FMField *objB )
{
  int32 i, j, k, il, wr, wa, wb;
  float64 *pr, *pa, *pb;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objB->nRow)
      || (objA->nCol != objB->nCol)
      || (objR->nLev != objA->nLev) || (objR->nLev != objB->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) = (%d %d %d) * (%d %d %d)^T\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol,
	    objB->nLev, objB->nRow, objB->nCol );
  }
#endif

  wr = objR->nCol;
  wa = objA->nCol;
  wb = objB->nCol;
  for (il = 0; il < objR->nLev; il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pa = objA->val + objA->nCol * objA->nRow * il;
    pb = objB->val + objB->nCol * objB->nRow * il;
    for (i = 0; i < objR->nRow; i++) {
      for (j = 0; j < objR->nCol; j++) {
	pr[wr*i+j] = 0.0;
	for (k = 0; k < objA->nCol; k++) {
	  pr[objR->nCol*i+j] += pa[wa*i+k] * pb[wb*j+k];
	}
      }
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_mulATBT_nn"
/*!
  objR = objA^T * objB^T

  @par Revision history:
  - 06.02.2001, c
  - 05.03.2003, adopted from rcfem2
*/
int32 fmf_mulATBT_nn( FMField *objR, FMField *objA, FMField *objB )
{
  int32 i, j, k, il, wr, wa, wb;
  float64 *pr, *pa, *pb;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nCol) || (objR->nCol != objB->nRow)
      || (objA->nRow != objB->nCol)
      || (objR->nLev != objA->nLev) || (objR->nLev != objB->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) = (%d %d %d)^T * (%d %d %d)^T\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol,
	    objB->nLev, objB->nRow, objB->nCol );
  }
#endif

  wr = objR->nCol;
  wa = objA->nCol;
  wb = objB->nCol;
  for (il = 0; il < objR->nLev; il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pa = objA->val + objA->nCol * objA->nRow * il;
    pb = objB->val + objB->nCol * objB->nRow * il;
    for (i = 0; i < objR->nRow; i++) {
      for (j = 0; j < objR->nCol; j++) {
	pr[wr*i+j] = 0.0;
	for (k = 0; k < objA->nRow; k++) {
	  pr[wr*i+j] += pa[wa*k+i] * pb[wb*j+k];
	}
      }
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_mulATBT_1n"
/*!
  objR = objA^T * objB^T, @a objA has exactly 1 level, which is used to
  multiply all levels of @a objB.

  @par Revision history:
  - 09.02.2001, c
  - 22.02.2001
  - 05.03.2003, adopted from rcfem2
*/
int32 fmf_mulATBT_1n( FMField *objR, FMField *objA, FMField *objB )
{
  int32 i, j, k, il, wr, wa, wb;
  float64 *pr, *pa, *pb;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nCol) || (objR->nCol != objB->nRow)
      || (objA->nRow != objB->nCol)
      || (1 != objA->nLev) || (objR->nLev != objB->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) = (%d(1) %d %d)^T * (%d %d %d)^T\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol,
	    objB->nLev, objB->nRow, objB->nCol );
    errput( ErrHead "ERR_BadMatch\n" );
  }
#endif

  wr = objR->nCol;
  wa = objA->nCol;
  wb = objB->nCol;
  pa = objA->val;
  for (il = 0; il < objR->nLev; il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pb = objB->val + objB->nCol * objB->nRow * il;
    for (i = 0; i < objR->nRow; i++) {
      for (j = 0; j < objR->nCol; j++) {
	pr[wr*i+j] = 0.0;
	for (k = 0; k < objA->nRow; k++) {
	  pr[wr*i+j] += pa[wa*k+i] * pb[wb*j+k];
	}
      }
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_addAB_nn"
/*!
  objR = objA + objB

  @par Revision history:
  - 06.02.2001, c
  - 05.03.2003, adopted from rcfem2
*/
int32 fmf_addAB_nn( FMField *objR, FMField *objA, FMField *objB )
{
  int32 i;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objB->nCol)
      || (objA->nRow != objB->nRow)
      || (objR->nLev != objA->nLev) || (objR->nLev != objB->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) = (%d %d %d) + (%d %d %d)\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol,
	    objB->nLev, objB->nRow, objB->nCol );
  }
#endif

  for (i = 0; i < objR->nLev * objR->nRow * objR->nCol; i++) {
//      if (debug == 1) {
//        printf( "%d %p %p %p %d %f %f\n",
//  	      objR->nLev * objR->nRow * objR->nCol,
//  	      this, a, b, i, objA->val[i], objB->val[i] );
//        debug = 0;
//      }
    
    objR->val[i] = objA->val[i] + objB->val[i];
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_subAB_nn"
/*!
  objR = objA - objB

  @par Revision history:
  - 04.04.2003, c
*/
int32 fmf_subAB_nn( FMField *objR, FMField *objA, FMField *objB )
{
  int32 ii;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objB->nCol)
      || (objA->nRow != objB->nRow)
      || (objR->nLev != objA->nLev) || (objR->nLev != objB->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) = (%d %d %d) + (%d %d %d)\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol,
	    objB->nLev, objB->nRow, objB->nCol );
  }
#endif

  for (ii = 0; ii < objR->nLev * objR->nRow * objR->nCol; ii++) {
    objR->val[ii] = objA->val[ii] - objB->val[ii];
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfc_addAB_nn"
/*!
  objR = objA + objB

  @par Revision history:
  - 15.07.2002, c
  - 05.03.2003, adopted from rcfem2
*/
int32 fmfc_addAB_nn( FMField *objR, FMField *objA, FMField *objB )
{
  int32 i;

#ifdef DEBUG_FMF
  if ((objR->nCell != objA->nCell) || (objR->nCell != objB->nCell)
      || (objR->nRow != objA->nRow) || (objR->nCol != objB->nCol)
      || (objA->nRow != objB->nRow)
      || (objR->nLev != objA->nLev) || (objR->nLev != objB->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d %d) = (%d %d %d %d) + (%d %d %d %d)\n",
	    objR->nCell, objR->nLev, objR->nRow, objR->nCol,
	    objA->nCell, objA->nLev, objA->nRow, objA->nCol,
	    objB->nCell, objB->nLev, objB->nRow, objB->nCol );
  }
#endif

  for (i = 0; i < (objR->nCell * objR->nLev * objR->nRow * objR->nCol); i++) {
    objR->val0[i] = objA->val0[i] + objB->val0[i];
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_averageCACB"
/*!
  objR = c1 * objA + c2 * objB.

  @par Revision history:
  - 09.11.2001, c
  - 05.03.2003, adopted from rcfem2
*/
int32 fmf_averageCACB( FMField *objR, float64 c1, FMField *objA,
		       float64 c2, FMField *objB )
{
  int32 ii;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objB->nCol)
      || (objA->nRow != objB->nRow)
      || (objR->nLev != objA->nLev) || (objR->nLev != objB->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) = c1 * (%d %d %d) + c2 * (%d %d %d)\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol,
	    objB->nLev, objB->nRow, objB->nCol );
  }
#endif

  for (ii = 0; ii < objR->nLev * objR->nRow * objR->nCol; ii++) {
    objR->val[ii] = c1 * objA->val[ii] + c2 * objB->val[ii];
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfc_averageCACB"
/*!
  objR = c1 * objA + c2 * objB.

  @par Revision history:
*/
int32 fmfc_averageCACB( FMField *objR, float64 c1, FMField *objA,
			float64 c2, FMField *objB )
{
  int32 ii;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objB->nCol)
      || (objA->nRow != objB->nRow)
      || (objR->nLev != objA->nLev) || (objR->nLev != objB->nLev)
      || (objR->nCell != objA->nCell) || (objR->nCell != objB->nCell)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d %d) = c1 * (%d %d %d %d)"
	    " + c2 * (%d %d %d %d)\n",
	    objR->nCell, objR->nLev, objR->nRow, objR->nCol,
	    objA->nCell, objA->nLev, objA->nRow, objA->nCol,
	    objB->nCell, objB->nLev, objB->nRow, objB->nCol );
  }
#endif

  for (ii = 0; ii < objR->nAlloc; ii++) {
    objR->val0[ii] = c1 * objA->val0[ii] + c2 * objB->val0[ii];
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfc_normalize"
/*!
  Vectors of zero norms are copied without warning...

  @par Revision history:
  - 21.03.2003, c
*/
int32 fmfc_normalize( FMField *objR, FMField *objA )
{
  int32 ii, il, ic, dim;
  float64 norm;
  float64 *pr, *pa;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objA->nCol)
      || (objR->nLev != objA->nLev) || (objR->nCell != objA->nCell)) {
    errput( ErrHead "ERR_BadMatch\n" );
  }
#endif

  pr = objR->val0;
  pa = objA->val0;
  dim = objR->nRow * objR->nCol;
  for (ic = 0; ic < objR->nCell; ic++) {
    for (il = 0; il < (objR->nLev); il++) {
      norm = 0.0;
      for (ii = 0; ii < dim; ii++) {
	norm += pa[ii] * pa[ii];
      }
      if (norm > MachEps) {
	for (ii = 0; ii < dim; ii++) {
	  pr[ii] = pa[ii] / norm;
	}
      } else {
	pr[ii] = pa[ii];
      }
      pr += dim;
      pa += dim;
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_addAmulF"
/*!
  @par Revision history:
  - 23.05.2005, c
  - 25.05.2005
*/
int32 fmf_addAmulF( FMField *objR, FMField *objA, float64 *val )
{
  int32 ii, il;
  float64 *pr, *pa;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objA->nCol)
      || (objR->nLev != objA->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) = (%d %d %d) * C\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol );
  }
#endif

  pr = objR->val;
  pa = objA->val;
  for (il = 0; il < (objR->nLev); il++) {
    for (ii = 0; ii < (objR->nRow * objR->nCol); ii++) {
      pr[ii] += pa[ii] * val[il];
    }
    pr += objR->nCol * objR->nRow;
    pa += objA->nCol * objA->nRow;
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfc_addAmulF"
/*!
  @par Revision history:
  - 15.10.2001, c
  - 17.10.2001
  - 05.03.2003, adopted from rcfem2
*/
int32 fmfc_addAmulF( FMField *objR, FMField *objA, float64 *val )
{
  int32 i, il, ic;
  float64 *pr, *pa;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objA->nCol)
      || (objR->nLev != objA->nLev) || (objR->nCell != objA->nCell)) {
    errput( ErrHead "ERR_BadMatch\n" );
  }
#endif

  pr = objR->val0;
  pa = objA->val0;
  for (ic = 0; ic < objR->nCell; ic++) {
    for (il = 0; il < (objR->nLev); il++) {
      for (i = 0; i < (objR->nRow * objR->nCol); i++) {
	pr[i] += pa[i] * val[il];
      }
      pr += objR->nCol * objR->nRow;
      pa += objA->nCol * objA->nRow;
    }
    val += objR->nLev;
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_copyAmulC"
/*!
  @par Revision history:
  - 02.04.2003, c
*/
int32 fmf_copyAmulC( FMField *objR, FMField *objA, float64 val )
{
  int32 ii;
  float64 *pr, *pa;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objA->nCol)
      || (objR->nLev != objA->nLev)) {
    errput( ErrHead "ERR_BadMatch\n" );
  }
#endif

  pr = objR->val;
  pa = objA->val;
  for (ii = 0; ii < (objR->nLev * objR->nRow * objR->nCol); ii++) {
    pr[ii] = pa[ii] * val;
  }

  return( RET_OK );
}
#undef __FUNC__
#define __FUNC__ "fmfc_copyAmulF"
/*!
  @par Revision history:
  - 18.10.2001, c
  - 05.03.2003, adopted from rcfem2
*/
int32 fmfc_copyAmulF( FMField *objR, FMField *objA, float64 *val )
{
  int32 i, il, ic;
  float64 *pr, *pa;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objA->nCol)
      || (objR->nLev != objA->nLev) || (objR->nCell != objA->nCell)) {
    errput( ErrHead "ERR_BadMatch\n" );
  }
#endif

  pr = objR->val0;
  pa = objA->val0;
  for (ic = 0; ic < objR->nCell; ic++) {
    for (il = 0; il < (objR->nLev); il++) {
      for (i = 0; i < (objR->nRow * objR->nCol); i++) {
	pr[i] = pa[i] * val[il];
      }
      pr += objR->nCol * objR->nRow;
      pa += objA->nCol * objA->nRow;
    }
    val += objR->nLev;
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfr_addA_blockNC"
/*!
  Adds @a objA as a block into position (@a row, @a col) of @a objR.
  Rectangular submatrix support for @a objR.

  No size checks!

  @par Revision history:
  - 24.04.2001, c
  - 05.03.2003, adopted from rcfem2
  - 31.03.2003
*/
int32 fmfr_addA_blockNC( FMField *objR, FMField *objA, int32 row, int32 col )
{
  int32 il, ir, ic;
  int32 wa, wr, ha, hr;
  float64 *pr, *pa;

  wa = objA->nCol;
  wr = objR->nColFull;
  ha = objA->nRow;
  hr = objR->nRow;
  for (il = 0; il < objR->nLev; il++) {
    pr = objR->val + objR->offset + wr * hr * il + wr * row + col;
    pa = objA->val + wa * ha * il;
    for (ir = 0; ir < ha; ir++) {
      for (ic = 0; ic < wa; ic++) {
	pr[ic] += pa[ic];
      }
      pa += wa;
      pr += wr;
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfr_addAT_blockNC"
/*!
  Adds @a objA^T as a block into position (@a row, @a col) of @a objR.
  Can be used instead of transposition.
  Rectangular submatrix support for @a objR.

  No size checks!

  @par Revision history:
  - 24.04.2001, c
  - 05.03.2003, adopted from rcfem2
  - 31.03.2003
  - 21.01.2004
*/
int32 fmfr_addAT_blockNC( FMField *objR, FMField *objA, int32 row, int32 col )
{
  int32 il, ir, ic;
  int32 wa, wr, ha, hr;
  float64 *pr, *pa;

  wa = objA->nCol;
  wr = objR->nColFull;
  ha = objA->nRow;
  hr = objR->nRow;
/*   output( "%d %d %d %d %d %d\n", hr, wr, ha, wa, row, col ); */
  for (il = 0; il < objR->nLev; il++) {
    pr = objR->val + objR->offset + wr * hr * il + wr * row + col;
    pa = objA->val + wa * ha * il;
    for (ir = 0; ir < wa; ir++) {
      for (ic = 0; ic < ha; ic++) {
/* 	output( "ir: %d ic: %d -> %d\n", ir, ic, pr + ic - objR->val ); */
	pr[ic] += pa[wa*ic+ir];
      }
      pr += wr;
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_sumLevelsMulF"
/*!
  @par Revision history:
  - 26.04.2001, c
  - 19.01.2004, adopted from rcfem2
*/
int32 fmf_sumLevelsMulF( FMField *objR, FMField *objA, float64 *val )
{
  int32 il, i;
  float64 *pa;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objA->nCol)
      || (objR->nLev != 1)) {
    errput( ErrHead "ERR_BadMatch (%d == %d, %d == %d, %d == 1)\n",
	    objR->nRow, objA->nRow, objR->nCol, objA->nCol, objR->nLev );
  }
#endif

  fmf_fillC( objR, 0.0 );
  for (il = 0; il < objA->nLev; il++) {
    pa = objA->val + objA->nCol * objA->nRow * il;
    for (i = 0; i < (objR->nRow * objR->nCol); i++) {
      objR->val[i] += pa[i] * val[il];
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_sumLevelsTMulF"
/*!
  @a objA^T.

  @par Revision history:
  - 23.05.2005, c
*/
int32 fmf_sumLevelsTMulF( FMField *objR, FMField *objA, float64 *val )
{
  int32 il, ir, ic, wr, wc;
  float64 *pa, *pr;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nCol) || (objR->nCol != objA->nRow)
      || (objR->nLev != 1)) {
    errput( ErrHead "ERR_BadMatch (%d == %d, %d == %d, %d == 1)\n",
	    objR->nRow, objA->nCol, objR->nCol, objA->nRow, objR->nLev );
  }
#endif

  wr = objR->nCol;
  pr = objR->val;
  wc = objA->nCol;

  fmf_fillC( objR, 0.0 );
  for (il = 0; il < objA->nLev; il++) {
    pa = objA->val + objA->nCol * objA->nRow * il;
    for (ir = 0; ir < objR->nRow; ir++) {
      for (ic = 0; ic < objR->nCol; ic++) {
	pr[wr*ir+ic] += pa[wc*ic+ir] * val[il];
      }
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfr_sumLevelsMulF"
/*!
  Rectangular submatrix support for @a objR.

  @par Revision history:
  - 26.04.2001, c
  - 05.03.2003, adopted from rcfem2
  - 31.03.2003
*/
int32 fmfr_sumLevelsMulF( FMField *objR, FMField *objA, float64 *val )
{
  int32 il, ir, ic, ii;
  float64 *pa, *pr;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objA->nCol)
      || (objR->nLev != 1)) {
    errput( ErrHead "ERR_BadMatch (%d == %d, %d == %d, %d == 1)\n",
	    objR->nRow, objA->nRow, objR->nCol, objA->nCol, objR->nLev );
  }
#endif

/*    output( "%d %d\n", objR->offset, objR->nColFull ); */
  pr = objR->val + objR->offset;
  for (ir = 0; ir < objR->nRow; ir++) {
    for (ic = 0; ic < objR->nCol; ic++) {
      pr[ic] = 0.0;
    }
    pr += objR->nColFull;
  }
  for (il = 0; il < objA->nLev; il++) {
    pr = objR->val + objR->offset;
    pa = objA->val + objA->nCol * objA->nRow * il;
    ii = 0;
    for (ir = 0; ir < objR->nRow; ir++) {
      for (ic = 0; ic < objR->nCol; ic++, ii++) {
	pr[ic] += pa[ii] * val[il];
      }
      pr += objR->nColFull;
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfr_sumLevelsTMulF"
/*!
  Rectangular submatrix support for @a objR. @a objA^T.

  @par Revision history:
  - 23.05.2005, c
*/
int32 fmfr_sumLevelsTMulF( FMField *objR, FMField *objA, float64 *val )
{
  int32 il, ir, ic, wr, wa;
  float64 *pa, *pr;

#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nCol) || (objR->nCol != objA->nRow)
      || (objR->nLev != 1)) {
    errput( ErrHead "ERR_BadMatch (%d == %d, %d == %d, %d == 1)\n",
	    objR->nRow, objA->nCol, objR->nCol, objA->nRow, objR->nLev );
  }
#endif

  wr  = objR->nColFull;
  wa  = objA->nCol;

  pr = objR->val + objR->offset;
  for (ir = 0; ir < objR->nRow; ir++) {
    for (ic = 0; ic < objR->nCol; ic++) {
      pr[ic] = 0.0;
    }
    pr += wr;
  }
  for (il = 0; il < objA->nLev; il++) {
    pr = objR->val + objR->offset;
    pa = objA->val + objA->nCol * objA->nRow * il;
    for (ir = 0; ir < objR->nRow; ir++) {
      for (ic = 0; ic < objR->nCol; ic++) {
	pr[ic] += pa[wa*ic+ir] * val[il];
      }
      pr += objR->nColFull;
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_copy"
/*!
  @par Revision history:
  - 06.02.2001, c
  - 17.10.2001
  - 05.03.2003, adopted from rcfem2
*/
int32 fmf_copy( FMField *objR, FMField *objA )
{
  if (objR->cellSize != objA->cellSize) {
    errput(ErrHead "ERR_BadMatch: (%d %d %d) = (%d %d %d)\n",
           objR->nLev, objR->nRow, objR->nCol,
           objA->nLev, objA->nRow, objA->nCol);
  }
  memcpy( objR->val, objA->val, objA->cellSize * sizeof( float64 ) );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfr_copy"
/*!
  @par Revision history:
  - 04.04.2003, c
  - 07.04.2003
  - 20.01.2004
*/
int32 fmfr_copy( FMField *objR, FMField *objA )
{
  int32 il, ir, ic;
  int32 wa, wr, ha, hr;
  float64 *pa, *pr;
#ifdef DEBUG_FMF
  if ((objR->nRow != objA->nRow) || (objR->nCol != objA->nCol)
      || (objR->nLev != objA->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) = (%d %d %d)\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol );
  }
#endif
  wa = objA->nCol;
  wr = objR->nColFull;
  ha = objA->nRow;
  hr = objR->nRow;
  for (il = 0; il < objR->nLev; il++) {
    pr = objR->val + objR->offset + wr * hr * il;
    pa = objA->val + wa * ha * il;
    for (ir = 0; ir < ha; ir++) {
      for (ic = 0; ic < wa; ic++) {
	pr[wr*ir+ic] = pa[wa*ir+ic];
      }
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfc_copy"
/*!
  @par Revision history:
  - 17.10.2001, c
  - 05.03.2003, adopted from rcfem2
*/
int32 fmfc_copy( FMField *objR, FMField *objA )
{
  // Cannot use nAlloc because of pretended FMFields...
  if ((objR->nCell * objR->nLev * objR->nRow * objR->nCol)
      != (objA->nCell * objA->nLev * objA->nRow * objA->nCol)) {
    errput( ErrHead "ERR_BadMatch\n" );
  }
  memcpy( objR->val0, objA->val0,
	  objA->nCell * objA->nLev * objA->nRow * objA->nCol
	  * sizeof( float64 ) );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_print"
/*!
  @par Revision history:
  - 06.02.2001, c
  - 14.11.2001
  - 05.03.2003, adopted from rcfem2
  - 18.03.2003
  - 08.04.2003
*/
int32 fmf_print( FMField *obj, FILE *file, int32 mode )
{
  int32 i, j, il;

  if (mode == 0) {
    fprintf( file, FI32" "FI32" "FI32"\n", obj->nLev, obj->nRow, obj->nCol );
    for (il = 0; il < obj->nLev; il++) {
      fprintf( file, FI32"\n", il );
      for (i = 0; i < obj->nRow; i++) {
	for (j = 0; j < obj->nCol; j++) {
/* 	  fprintf( file, " %.12e", obj->val[obj->nCol*(obj->nRow*il+i)+j] ); */
	  fprintf( file, " %.6e", obj->val[obj->nCol*(obj->nRow*il+i)+j] );
	}
	fprintf( file, "\n" );
      }
    }
  } else if (mode == 1) {
    fprintf( file, "nCell: "FI32" nLev: "FI32" nRow: "FI32" nCol: "FI32"\n",
	     obj->nCell, obj->nLev, obj->nRow, obj->nCol );
    fprintf( file, "offset: "FI32" nColFull: "FI32" nAlloc: "FI32
	     " cellSize "FI32"\n",
	     obj->offset, obj->nColFull, obj->nAlloc, obj->cellSize );
  } else {
    errput( ErrHead "ERR_Switch!\n" );
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfr_print"
/*!
  @par Revision history:
  - 04.12.2003, c
*/
int32 fmfr_print( FMField *obj, FILE *file, int32 mode )
{
  int32 i, j, il;

  if (mode == 0) {
    fprintf( file, FI32" "FI32" "FI32" "FI32" "FI32"\n",
	     obj->nLev, obj->nRow, obj->nCol, obj->offset, obj->nColFull );
    for (il = 0; il < obj->nLev; il++) {
      fprintf( file, FI32"\n", il );
      for (i = 0; i < obj->nRow; i++) {
	for (j = 0; j < obj->nCol; j++) {
	  fprintf( file, " %.12e",
		   obj->val[obj->offset+obj->nColFull*(obj->nRow*il+i)+j] );
	}
	fprintf( file, "\n" );
      }
    }
  } else if (mode == 1) {
    fmf_print( obj, file, mode );
  } else {
    errput( ErrHead "ERR_Switch!\n" );
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_save"
/*!
  @par Revision history:
  - 06.02.2001, c
  - 05.03.2003, adopted from rcfem2
  - 18.03.2003
*/
int32 fmf_save( FMField *obj, const char *fileName, int32 mode )
{
  FILE *file;

  if ((file = fopen( fileName, "w" )) == 0) {
    errput( ErrHead "ERR_FileOpen\n" );
  }

  fmf_print( obj, file, mode );

  fclose( file );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfr_save"
/*!
  @par Revision history:
  - 15.04.2004, c
*/
int32 fmfr_save( FMField *obj, const char *fileName, int32 mode )
{
  FILE *file;

  if ((file = fopen( fileName, "w" )) == 0) {
    errput( ErrHead "ERR_FileOpen\n" );
  }

  fmfr_print( obj, file, mode );

  fclose( file );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmfc_save"
/*!
  @par Revision history:
  - 15.07.2002, c
  - 26.07.2002
  - 05.03.2003, adopted from rcfem2
  - 18.03.2003
*/
int32 fmfc_save( FMField *obj, const char *fileName, int32 mode )
{
  int32 ii;
  FILE *file;

  if ((file = fopen( fileName, "w" )) == 0) {
    errput( ErrHead "ERR_FileOpen\n" );
  }

  if (mode == 0) {
    FMF_SetFirst( obj );
    for (ii = 0; ii < obj->nCell; ii++) {
      fmf_print( obj, file, 0 );
      FMF_SetCellNext( obj );
    }
  } else if (mode == 1) {
    fprintf( file, FI32"\n", obj->nAlloc );
    for (ii = 0; ii < obj->nAlloc; ii++) {
      fprintf( file, FI32" %.12e\n", ii, obj->val0[ii] );
    }
  }

  fclose( file );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_gMtx2VecDUL3x3"
/*!
  @par Revision history:
  - 20.04.2001, c
  - 22.04.2001
  - 17.09.2001
  - 20.01.2004, adopted from rcfem2
*/
int32 fmf_gMtx2VecDUL3x3( FMField *objR, FMField *objA )
{
  int32 i, il;
  static int32 order[][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
			     {0, 3, 1, 2, 0, 0, 0, 0, 0},
			     {0, 4, 8, 1, 2, 5, 3, 6, 7}};
  int32 *pord;
  float64 *pr, *pa;

#ifdef DEBUG_FMF
  if ((objR->nRow != (objA->nRow * objA->nCol)) || (objR->nLev != objA->nLev)
    || (objR->nCol != 1) || (objA->nRow > 3) || (objA->nCol != objA->nRow)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) <-> (%d %d %d)\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol );
  }
#endif

  pord = order[objA->nRow-1];

  for (il = 0; il < objR->nLev; il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pa = objA->val + objA->nCol * objA->nRow * il;
    for (i = 0; i < (objR->nRow); i++) {
      pr[i] = pa[pord[i]];
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "fmf_gMtx2VecDLU3x3"
/*!
  @par Revision history:
  - 22.04.2001, c
  - 17.09.2001
  - 20.01.2004, adopted from rcfem2
*/
int32 fmf_gMtx2VecDLU3x3( FMField *objR, FMField *objA )
{
  int32 i, il;
  static int32 order[][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
			     {0, 3, 2, 1, 0, 0, 0, 0, 0},
			     {0, 4, 8, 3, 6, 7, 1, 2, 5}};
  int32 *pord;
  float64 *pr, *pa;

#ifdef DEBUG_FMF
  if ((objR->nRow != (objA->nRow * objA->nCol)) || (objR->nLev != objA->nLev)
    || (objR->nCol != 1) || (objA->nRow > 3) || (objA->nCol != objA->nRow)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) <-> (%d %d %d)\n",
	    objR->nLev, objR->nRow, objR->nCol,
	    objA->nLev, objA->nRow, objA->nCol );
  }
#endif

  pord = order[objA->nRow-1];

  for (il = 0; il < objR->nLev; il++) {
    pr = objR->val + objR->nCol * objR->nRow * il;
    pa = objA->val + objA->nCol * objA->nRow * il;
    for (i = 0; i < (objR->nRow); i++) {
      pr[i] = pa[pord[i]];
    }
  }

  return( RET_OK );
}
