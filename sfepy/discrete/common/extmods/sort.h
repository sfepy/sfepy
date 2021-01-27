#ifndef _SORT_H_
#define _SORT_H_

#include "common.h"
BEGIN_C_DECLS

typedef uint32 intp;

int32 int32_sort_rows( int32 *array, int32 n_row, int32 n_col,
		      int32 *i_sort_col, int32 n_sort_col );
int32 int32_quicksort(int32 *start, intp num, void *unused);
int32 int32_aquicksort(int32 *v, intp* tosort, intp num, void *unused);
int32 int32_mtx_aquicksort( int32 *v, int32 n_row, int32 n_col,
			    int32 *ic, int32 n_c,
			    intp* tosort );

END_C_DECLS

#endif /* Header */
