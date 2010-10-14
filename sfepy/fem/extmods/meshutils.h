#ifndef _MESHUTILS_H_
#define _MESHUTILS_H_

#include "common.h"

int32 assign_edge_nodes( int32 *p_iseq,  
		       int32 *econn, int32 econn_n_row, int32 econn_n_col,
		       int32 *cnt_en, int32 cnt_en_n_row, int32 cnt_en_n_col,
		       int32 *ori, int32 ori_n_row, int32 ori_n_col,
		       int32 *ntt, int32 ntt_n_row, int32 ntt_n_col,
		       int32 *uid, int32 uid_n_row,
		       int32 **edge_desc, int32 *edge_desc_n_row, int32 n_edge,
		       int32 cptr0 );

int32 interp_vertex_data( float64 *nod_out, int32 nod_out_n_row, int32 nod_out_n_col,
			int32 *econn, int32 econn_n_row, int32 econn_n_col,
			float64 *nod_in, int32 nod_in_n_row, int32 nod_in_n_col,
			int32 *conn, int32 conn_n_row, int32 conn_n_col,
			float64 *bf, int32 bf_n_row, int32 bf_n_col,
			int32 omit_cols );

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
