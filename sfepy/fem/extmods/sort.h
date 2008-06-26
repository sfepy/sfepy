#ifndef _SORT_H_
#define _SORT_H_

#include "common.h"
BEGIN_C_DECLS

typedef unsigned int intp;

int32 int32_sortRows( int32 *array, int32 nRow, int32 nCol,
		      int32 *iSortCol, int32 nSortCol );
int32 int32_quicksort(int32 *start, intp num, void *unused);
int32 int32_aquicksort(int32 *v, intp* tosort, intp num, void *unused);
int32 int32_mtx_aquicksort( int32 *v, int32 nRow, int32 nCol,
			    int32 *ic, int32 nC,
			    intp* tosort );

END_C_DECLS

#endif /* Header */
