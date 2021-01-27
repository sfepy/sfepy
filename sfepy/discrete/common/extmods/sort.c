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

/*     output( "%d %d %d\n", ii, row_a[ic[ii]], row_b[ic[ii]] );\ */
#define MTXR_LT( ret, row_a, row_b ) do {\
  ret = 0;\
  for (ii = 0; ii < n_c; ii++) {\
    if (row_a[ic[ii]] < row_b[ic[ii]]) {\
      ret = 1;\
      break;\
    } else if (row_a[ic[ii]] > row_b[ic[ii]]) {\
      break;\
    }\
  }\
} while (0)
#define MTXR_Push( row ) do {\
  for (jj = 0; jj < n_col; jj++) {\
    SWAP_temp[jj] = row[jj];\
  }\
} while (0)
#define MTXR_Pop( row ) do {\
  for (jj = 0; jj < n_col; jj++) {\
    row[jj] = SWAP_temp[jj];\
  }\
} while (0)
#define MTXR_Copy( row_to, row_from ) do {\
  for (jj = 0; jj < n_col; jj++) {\
    row_to[jj] = row_from[jj];\
  }\
} while (0)

#undef __FUNC__
#define __FUNC__ "int32_sort_rows"
/*!
  @par Revision history:
  - 07.02.2006, c
*/
int32 int32_sort_rows( int32 *array, int32 n_row, int32 n_col,
		      int32 *i_sort_col, int32 n_sort_col )
{
  int32 *SWAP_temp;
  int32 ii, jj, inew, ichain;
  int32 *perm, *perm_i;

  perm = alloc_mem( int32, n_row );
  perm_i = alloc_mem( int32, n_row );
  SWAP_temp = alloc_mem( int32, n_col );

  /*
    Arg-sort rows: sorted = original[perm], sorted[perm_i] = original.
  */
  for (ii = 0; ii < n_row; ii++) {
    perm[ii] = ii;
  }

  int32_mtx_aquicksort( array, n_row, n_col, i_sort_col, n_sort_col, (intp *) perm );

  for (ii = 0; ii < n_row; ii++) {
    perm_i[perm[ii]] = ii;
  }
/*   output( "argsort done.\n" ); */
  /*
    Swap rows in place.
  */
  for (ii = 0; ii < n_row; ii++) {
    if (ii == perm[ii]) continue;

    MTXR_Push( (array + n_col * ii) );
    inew = perm[ii];
    MTXR_Copy( (array + n_col * ii), (array + n_col * inew) );
    perm[ii] = ii; /* Mark done (useless). */
    ichain = perm_i[ii];

/*     output( "%d %d %d\n", ii, inew, ichain ); */
    while (inew != ichain) {
      MTXR_Pop( (array + n_col * inew) );
      MTXR_Push( (array + n_col * ichain) );
      MTXR_Copy( (array + n_col * ichain), (array + n_col * inew) );

      perm[ichain] = ichain; /* Mark done. */
      ichain = perm_i[ichain];
/*       output( "   %d\n", ichain ); */
/*       getchar(); */
    }
    MTXR_Pop( (array + n_col * inew) );
    perm[inew] = inew; /* Mark done. */
  }
/*   output( "swap done.\n" ); */

  free_mem( perm );
  free_mem( perm_i );
  free_mem( SWAP_temp );

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
int32 int32_mtx_aquicksort( int32 *v, int32 n_row, int32 n_col,
			    int32 *ic, int32 n_c,
			    intp* tosort )
{
  int32 ret, ii;
  int32 *vp;
  intp *pl, *pr, SWAP_temp;
  intp *stack[PYA_QS_STACK], **sptr=stack, *pm, *pi, *pj, *pt, vi;

  pl = tosort;
  pr = tosort + n_row - 1;

  for(;;) {
    while ((pr - pl) > SMALL_QUICKSORT) {
      /* quicksort partition */
      pm = pl + ((pr - pl) >> 1);

      MTXR_LT( ret, (v + n_col * (*pm)), (v + n_col * (*pl)) );
      if (ret) SWAP(*pm,*pl);

      MTXR_LT( ret, (v + n_col * (*pr)), (v + n_col * (*pm)) );
      if (ret) SWAP(*pr,*pm);

      MTXR_LT( ret, (v + n_col * (*pm)), (v + n_col * (*pl)) );
      if (ret) SWAP(*pm,*pl);
      vp = v + n_col * (*pm);
      pi = pl;
      pj = pr - 1;
      SWAP(*pm,*pj);
      for(;;) {
	do {
	  ++pi;
	  MTXR_LT( ret, (v + n_col * (*pi)), vp );
	} while (ret);
	do {
	  --pj;
	  MTXR_LT( ret, vp, (v + n_col * (*pj)) );
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
      vp = v + n_col * vi;
      
      pj = pi;
      pt = pi - 1;
      while (pj > pl) {
	MTXR_LT( ret, vp, (v + n_col * (*pt)) );
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

