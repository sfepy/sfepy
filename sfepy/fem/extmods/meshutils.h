#ifndef _MESHUTILS_H_
#define _MESHUTILS_H_

#include "common.h"

int32 orient_elements( int32 *flag, int32 flag_n_row,
		      int32 *conn, int32 conn_n_row, int32 conn_n_col,
		      float64 *coors, int32 coors_n_row, int32 coors_n_col,
		      int32 *v_roots, int32 v_roots_n_row,
		      int32 *v_vecs, int32 v_vecs_n_row, int32 v_vecs_n_col,
		      int32 *swap_from, int32 swap_from_n_row, int32 swap_from_n_col,
		      int32 *swap_to, int32 swap_to_n_row, int32 swap_to_n_col );

int32 graph_components( int32 *p_n_comp,
		       int32 *flag, int32 flag_len,
		       int32 *row, int32 row_len,
		       int32 *col, int32 col_len,
		       int32 *pos, int32 pos_len );

#endif /* Header */
