/*
  These sorting functions are copied almost directly from numarray 
   with a few modifications (complex comparisons compare the imaginary 
   part if the real parts are equal, for example), and the names
   are changed. 

   The original sorting code is due to Charles R. Harris who wrote
   it for numarray.
*/

/* Quick sort is usually the fastest, but the worst case scenario can
   be slower than the merge and heap sorts.  The merge sort requires
   extra memory and so for large arrays may not be useful. 
 
   Specialized to sort 2D int32 matrix by rows by Robert Cimrman.
*/

#include "sort.h"

#define PYA_QS_STACK 100
#define MAX_NCOL 100
#define SMALL_QUICKSORT 15
#define STDC_LT(a,b) ((a) < (b))
#define STDC_LE(a,b) ((a) <= (b))
#define STDC_EQ(a,b) ((a) == (b))
#define SWAP(a,b) {SWAP_temp = (b); (b)=(a); (a) = SWAP_temp;}

/*     output( "%d %d %d\n", ii, rowA[ic[ii]], rowB[ic[ii]] );\ */
#define MTXR_LT( ret, rowA, rowB ) do {\
  ret = 0;\
  for (ii = 0; ii < nC; ii++) {\
    if (rowA[ic[ii]] < rowB[ic[ii]]) {\
      ret = 1;\
      break;\
    } else if (rowA[ic[ii]] > rowB[ic[ii]]) {\
      break;\
    }\
  }\
} while (0)
#define MTXR_Push( row ) do {\
  for (jj = 0; jj < nCol; jj++) {\
    SWAP_temp[jj] = row[jj];\
  }\
} while (0)
#define MTXR_Pop( row ) do {\
  for (jj = 0; jj < nCol; jj++) {\
    row[jj] = SWAP_temp[jj];\
  }\
} while (0)
#define MTXR_Copy( rowTo, rowFrom ) do {\
  for (jj = 0; jj < nCol; jj++) {\
    rowTo[jj] = rowFrom[jj];\
  }\
} while (0)

#undef __FUNC__
#define __FUNC__ "int32_sortRows"
/*!
  @par Revision history:
  - 07.02.2006, c
*/
int32 int32_sortRows( int32 *array, int32 nRow, int32 nCol,
		      int32 *iSortCol, int32 nSortCol )
{
  int32 *SWAP_temp;
  int32 ii, jj, inew, ichain;
  int32 *perm, *permI;

  perm = allocMem( int32, nRow );
  permI = allocMem( int32, nRow );
  SWAP_temp = allocMem( int32, nCol );

  /*
    Arg-sort rows: sorted = original[perm], sorted[permI] = original.
  */
  for (ii = 0; ii < nRow; ii++) {
    perm[ii] = ii;
  }

  int32_mtx_aquicksort( array, nRow, nCol, iSortCol, nSortCol, perm );

  for (ii = 0; ii < nRow; ii++) {
    permI[perm[ii]] = ii;
  }
/*   output( "argsort done.\n" ); */
  /*
    Swap rows in place.
  */
  for (ii = 0; ii < nRow; ii++) {
    if (ii == perm[ii]) continue;

    MTXR_Push( (array + nCol * ii) );
    inew = perm[ii];
    MTXR_Copy( (array + nCol * ii), (array + nCol * inew) );
    perm[ii] = ii; /* Mark done (useless). */
    ichain = permI[ii];

/*     output( "%d %d %d\n", ii, inew, ichain ); */
    while (inew != ichain) {
      MTXR_Pop( (array + nCol * inew) );
      MTXR_Push( (array + nCol * ichain) );
      MTXR_Copy( (array + nCol * ichain), (array + nCol * inew) );

      perm[ichain] = ichain; /* Mark done. */
      ichain = permI[ichain];
/*       output( "   %d\n", ichain ); */
/*       getchar(); */
    }
    MTXR_Pop( (array + nCol * inew) );
    perm[inew] = inew; /* Mark done. */
  }
/*   output( "swap done.\n" ); */

  freeMem( perm );
  freeMem( permI );
  freeMem( SWAP_temp );

  return( RET_OK );
}

/*!
  @par Revision history:
  - 07.02.2006, from numpy
*/
int32 int32_quicksort(int32 *start, intp num, void *unused)
{
  int32 *pl = start;
  int32 *pr = start + num - 1;
  int32 vp, SWAP_temp;
  int32 *stack[PYA_QS_STACK], **sptr = stack, *pm, *pi, *pj, *pt;

  for(;;) {
    while ((pr - pl) > SMALL_QUICKSORT) {
      /* quicksort partition */
      pm = pl + ((pr - pl) >> 1);
      if (STDC_LT(*pm,*pl)) SWAP(*pm,*pl);
      if (STDC_LT(*pr,*pm)) SWAP(*pr,*pm);
      if (STDC_LT(*pm,*pl)) SWAP(*pm,*pl);
      vp = *pm;
      pi = pl;
      pj = pr - 1;
      SWAP(*pm,*pj);
      for(;;) {
	do ++pi; while (STDC_LT(*pi,vp));
	do --pj; while (STDC_LT(vp,*pj));
	if (pi >= pj)  break;
	SWAP(*pi,*pj);
      }
      SWAP(*pi,*(pr-1));
      /* push largest partition on stack */
      if (pi - pl < pr - pi) {
	*sptr++ = pi + 1;
	*sptr++ = pr;
	pr = pi - 1;
      }else{
	*sptr++ = pl;
	*sptr++ = pi - 1;
	pl = pi + 1;
      }
    }
    /* insertion sort */
    for(pi = pl + 1; pi <= pr; ++pi) {
      vp = *pi;
      for(pj = pi, pt = pi - 1; \
	    pj > pl && STDC_LT(vp, *pt);) {
	*pj-- = *pt--;
      }
      *pj = vp;
    }
    if (sptr == stack) break;
    pr = *(--sptr);
    pl = *(--sptr);
  }
  return 0;
}

/*!
  @par Revision history:
  - 07.02.2006, from numpy
*/
int32 int32_aquicksort(int32 *v, intp* tosort, intp num, void *unused)
{
  int32 vp;
  intp *pl, *pr, SWAP_temp;
  intp *stack[PYA_QS_STACK], **sptr=stack, *pm, *pi, *pj, *pt, vi;

  pl = tosort;
  pr = tosort + num - 1;

  for(;;) {
    while ((pr - pl) > SMALL_QUICKSORT) {
      /* quicksort partition */
      pm = pl + ((pr - pl) >> 1);
      if (STDC_LT(v[*pm],v[*pl])) SWAP(*pm,*pl);
      if (STDC_LT(v[*pr],v[*pm])) SWAP(*pr,*pm);
      if (STDC_LT(v[*pm],v[*pl])) SWAP(*pm,*pl);
      vp = v[*pm];
      pi = pl;
      pj = pr - 1;
      SWAP(*pm,*pj);
      for(;;) {
	do ++pi; while (STDC_LT(v[*pi],vp));
	do --pj; while (STDC_LT(vp,v[*pj]));
	if (pi >= pj)  break;
	SWAP(*pi,*pj);
      }
      SWAP(*pi,*(pr-1));
      /* push largest partition on stack */
      if (pi - pl < pr - pi) {
	*sptr++ = pi + 1;
	*sptr++ = pr;
	pr = pi - 1;
      }else{
	*sptr++ = pl;
	*sptr++ = pi - 1;
	pl = pi + 1;
      }
    }
    /* insertion sort */
    for(pi = pl + 1; pi <= pr; ++pi) {
      vi = *pi;
      vp = v[vi];
      for(pj = pi, pt = pi - 1; \
	    pj > pl && STDC_LT(vp, v[*pt]);)
	{
	  *pj-- = *pt--;
	}
      *pj = vi;
    }
    if (sptr == stack) break;
    pr = *(--sptr);
    pl = *(--sptr);
  }
  return 0;
}


/*!
  @par Revision history:
  - 07.02.2006, c
*/
int32 int32_mtx_aquicksort( int32 *v, int32 nRow, int32 nCol,
			    int32 *ic, int32 nC,
			    intp* tosort )
{
  int32 ret, ii;
  int32 *vp;
  intp *pl, *pr, SWAP_temp;
  intp *stack[PYA_QS_STACK], **sptr=stack, *pm, *pi, *pj, *pt, vi;

  pl = tosort;
  pr = tosort + nRow - 1;

  for(;;) {
    while ((pr - pl) > SMALL_QUICKSORT) {
      /* quicksort partition */
      pm = pl + ((pr - pl) >> 1);

      MTXR_LT( ret, (v + nCol * (*pm)), (v + nCol * (*pl)) );
      if (ret) SWAP(*pm,*pl);

      MTXR_LT( ret, (v + nCol * (*pr)), (v + nCol * (*pm)) );
      if (ret) SWAP(*pr,*pm);

      MTXR_LT( ret, (v + nCol * (*pm)), (v + nCol * (*pl)) );
      if (ret) SWAP(*pm,*pl);
      vp = v + nCol * (*pm);
      pi = pl;
      pj = pr - 1;
      SWAP(*pm,*pj);
      for(;;) {
	do {
	  ++pi;
	  MTXR_LT( ret, (v + nCol * (*pi)), vp );
	} while (ret);
	do {
	  --pj;
	  MTXR_LT( ret, vp, (v + nCol * (*pj)) );
	} while (ret);
	if (pi >= pj)  break;
	SWAP(*pi,*pj);
      }
      SWAP(*pi,*(pr-1));
      /* push largest partition on stack */
      if (pi - pl < pr - pi) {
	*sptr++ = pi + 1;
	*sptr++ = pr;
	pr = pi - 1;
      }else{
	*sptr++ = pl;
	*sptr++ = pi - 1;
	pl = pi + 1;
      }
    }
    /* insertion sort */
    for(pi = pl + 1; pi <= pr; ++pi) {
      vi = (*pi);
      vp = v + nCol * vi;
      
      pj = pi;
      pt = pi - 1;
      while (pj > pl) {
	MTXR_LT( ret, vp, (v + nCol * (*pt)) );
	if (!ret) break;
	*pj-- = *pt--;
      }
      *pj = vi;
    }
    if (sptr == stack) break;
    pr = *(--sptr);
    pl = *(--sptr);
  }
  return 0;
}

