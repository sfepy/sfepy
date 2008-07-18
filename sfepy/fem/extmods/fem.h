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

int32 assemble_matrix( FMField *mtx,
		      int32 *prows, int32 prows_len,
		      int32 *cols, int32 cols_len,
		      FMField *mtxInEls,
		      int32 *iels, int32 iels_len, float64 sign,
		      int32 *connR, int32 nElR, int32 nEPR,
		      int32 *connC, int32 nElC, int32 nEPC );

int32 raw_graph( int32 *p_nRow, int32 **p_prow,
		int32 *p_nnz, int32 **p_icol,
		int32 nRow, int32 nCol, int32 nGr,
		int32 *nElR, int32 *nEPR, int32 **connR,
		int32 *nElC, int32 *nEPC, int32 **connC );

END_C_DECLS

#endif /* Header */
