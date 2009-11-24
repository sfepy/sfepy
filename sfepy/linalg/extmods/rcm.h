/*!
  @par Revision history:
  - 04.06.2001, adopted for rcfem2
*/
#ifndef _RCM_H_
#define _RCM_H_

#include "common.h"
BEGIN_C_DECLS

/*!
  @file
  @version v.1.0
 
  @short Reversed Cuthill-McKee permutation of a sparse matrix.

  "these subroutines appear in the book "computer solution of large sparse
  positive definite systems" by george and liu, prentice hall 19
  they form the basic building blocks for part of the sparse matrix
  package known as sparspak, currently distributed under license 
  by the university of waterloo. a few of the routines in the
  package are revised/improved versions of those that are published
  in the book cited above."

  careful! anything free comes with no guarantee
  from netlib, sat oct 28 12:12:27 edt 1989

  Adopted from code by Ales Janka, 17.4.1998.

  @par Revision history:
  - 11.10.2000, 0.27.0, c
  - 16.10.2000, 0.27.1
*/

void rcm_genrcm(int32 *perm, int32 neqns, int32 *xadj, int32 n_ptr,
		int32 *adjncy, int32 n_indx);

int32 gr_permuteInPlace(int32 *row, int32 n_row,
			int32 *col, int32 n_col,
			int32 *perm, int32 n_perm,
			int32 *permI, int32 n_permI);

END_C_DECLS

#endif /* Header */
