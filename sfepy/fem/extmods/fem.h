/*!
  @par Revision history:
  - 21.11.2005, c
*/
#ifndef _FEM_H_
#define _FEM_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"

int32 assemble_vector( FMField *vec, FMField *vecInEls,
		       int32 *iels, int32 iels_len,
		       float64 sign, int32 *conn, int32 nEl, int32 nEP );

int32 assemble_vector_complex( FMField *vec_r, FMField *vec_i,
			       FMField *vecInEls_r, FMField *vecInEls_i,
			       int32 *iels, int32 iels_len,
			       float64 sign_r, float64 sign_i,
			       int32 *conn, int32 nEl, int32 nEP );

int32 assemble_matrix( FMField *mtx,
		      int32 *prows, int32 prows_len,
		      int32 *cols, int32 cols_len,
		      FMField *mtxInEls,
		      int32 *iels, int32 iels_len, float64 sign,
		      int32 *connR, int32 nElR, int32 nEPR,
		      int32 *connC, int32 nElC, int32 nEPC );

int32 assemble_matrix_complex( FMField *mtx_r, FMField *mtx_i,
			       int32 *prows, int32 prows_len,
			       int32 *cols, int32 cols_len,
			       FMField *mtxInEls_r, FMField *mtxInEls_i,
			       int32 *iels, int32 iels_len,
			       float64 sign_r, float64 sign_i,
			       int32 *connR, int32 nElR, int32 nEPR,
			       int32 *connC, int32 nElC, int32 nEPC );

int32 raw_graph( int32 *p_nRow, int32 **p_prow,
		int32 *p_nnz, int32 **p_icol,
		int32 nRow, int32 nCol, int32 nGr,
		int32 *nElR, int32 *nEPR, int32 **connR,
		int32 *nElC, int32 *nEPC, int32 **connC );

int32 eval_lagrange_simplex( FMField *out, FMField *coors,
			     int32 *nodes, int32 nNod, int32 nCol,
			     int32 order, int32 diff,
			     FMField *mtx_i, FMField *bc,
			     int32 suppress_errors, float64 eps );

int32 eval_lagrange_tensor_product( FMField *out, FMField *coors,
				    int32 *nodes, int32 nNod, int32 nCol,
				    int32 order, int32 diff,
				    FMField *mtx_i, FMField *bc, FMField *base1d,
				    int32 suppress_errors, float64 eps );

int32 inverse_element_mapping( FMField *out,
			       FMField *coors, FMField *e_coors,
			       FMField *ref_coors, int32 i_max, float64 eps );

int32 evaluate_at( FMField *out,
		   int32 *cells, int32 n_cells, int32 n_cells_col,
		   int32 *status, int32 n_status,
		   FMField *dest_coors, FMField *source_vals,
		   int32 *ics, int32 n_ics,
		   int32 *offsets, int32 n_offsets,
		   int32 *iconn0, int32 n_iconn0,
		   FMField *mesh_coors,
		   int32 *nEls0, int32 *nEPs0, int32 **conns0,
		   int32 *nEls, int32 *nEPs, int32 **conns,
		   int32 n_ref_coorss, FMField *ref_coorss,
		   int32 *nNod, int32 *nCol, int32 **nodess,
		   int32 *orders, int32 n_orders,
		   int32 n_mtx_is, FMField *mtx_is,
		   int32 allow_extrapolation,
		   float64 close_limit, float64 qp_eps,
		   int32 i_max, float64 newton_eps );

END_C_DECLS

#endif /* Header */
