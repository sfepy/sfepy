#include "terms.h"

/*!
  @par Revision history:
  - 04.08.2006, c
  - 24.04.2007
*/
void debug_printConn( int32 *conn, int32 nEP )
{
  int32 ii;

  for (ii = 0; ii < nEP; ii++) {
    printf( "%d ", conn[ii] );
  }
  printf( "\n" );
}

#undef __FUNC__
#define __FUNC__ "ele_extractNodalValuesNBN"
/*!
  Extract values from element nodes, for example coordinates of nodes.
  Each node must have @ dim values!

  Input vector order is node-by-node.
  The extraction order is node-by-node.

  Use when local shape is (nEP, dpn).

  @par Revision history:
  - 26.10.2005, adapted from mafest1
*/
int32 ele_extractNodalValuesNBN( FMField *out, FMField *in,
				 int32 *conn )
{
  int32 inod, idof;

  for (inod = 0; inod < out->nRow; inod++) {
    for (idof = 0; idof < out->nCol; idof++ ) {
      out->val[out->nCol*inod+idof] = in->val[out->nCol*conn[inod]+idof];
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "ele_extractNodalValuesDBD"
/*!
  Extract values from element nodes, for example coordinates of nodes.
  Each node must have @ dim values!

  Input vector order is node-by-node.
  The extraction order is dof-by-dof.

  Use when local shape is (dpn, nEP).

  @par Revision history:
  - 14.12.2005, c
*/
int32 ele_extractNodalValuesDBD( FMField *out, FMField *in,
				 int32 *conn )
{
  int32 inod, idof;

  for (idof = 0; idof < out->nRow; idof++ ) {
    for (inod = 0; inod < out->nCol; inod++) {
      out->val[out->nCol*idof+inod] = in->val[out->nRow*conn[inod]+idof];
    }
  }

  return( RET_OK );
}
