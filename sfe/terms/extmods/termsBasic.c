#include "termsBasic.h"
#include "terms.h"

#undef __FUNC__
#define __FUNC__ "dq_state_in_qp"
/*!
  @par Revision history:
  - 24.04.2007, c
*/
int32 dq_state_in_qp( FMField *out, FMField *state, int32 offset,
		      FMField *bf,
		      int32 *conn, int32 nEl, int32 nEP )
{
  int32 ii, ret = RET_OK;
  FMField *st = 0;

  if (nEP !=  bf->nCol) {
    errput( "nEP mismatch: %d == %d!", nEP, bf->nCol );
  }
/*   if (out->nRow != 1) { */
/*     errput( "state_in_qp works only for scalars!" ); */
/*   } */

  state->val = FMF_PtrFirst( state ) + offset;

  fmf_createAlloc( &st, 1, 1, out->nRow, nEP );

  for (ii = 0; ii < nEl; ii++) {
    FMF_SetCell( out, ii );

    ele_extractNodalValuesDBD( st, state, conn + nEP * ii );
    bf_act( out, bf, st );
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &st ); 

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_grad_scalar"
/*!
  @par Revision history:
  - 30.07.2007, from dw_hdpm_cache()
*/
int32 dq_grad_scalar( FMField *out, FMField *state, int32 offset,
		      VolumeGeometry *vg,
		      int32 *conn, int32 nEl, int32 nEP )
{
  int32 ii, nQP, ret = RET_OK;
  FMField *st = 0;

  state->val = FMF_PtrFirst( state ) + offset;

  nQP = vg->bfGM->nLev;

  fmf_createAlloc( &st, 1, 1, nEP, 1 );

  for (ii = 0; ii < nEl; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, ii );

    ele_extractNodalValuesNBN( st, state, conn + nEP * ii );
    fmf_mulAB_n1( out, vg->bfGM, st );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &st ); 

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_div_vector"
/*!
  @par Revision history:
  - 30.07.2007, from dw_hdpm_cache()
*/
int32 dq_div_vector( FMField *out, FMField *state, int32 offset,
		     VolumeGeometry *vg,
		     int32 *conn, int32 nEl, int32 nEP )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *st = 0;
  FMField gcl[1], stv[1];

  state->val = FMF_PtrFirst( state ) + offset;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  fmf_createAlloc( &st, 1, 1, dim, nEP );
  stv->nAlloc = -1;
  fmf_pretend( stv, 1, 1, nEP * dim, 1, st->val );

  gcl->nAlloc = -1;
  fmf_pretend( gcl, 1, nQP, 1, nEP * dim, vg->bfGM->val0 );

  for (ii = 0; ii < nEl; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( gcl, ii );

    ele_extractNodalValuesDBD( st, state, conn + nEP * ii );
    fmf_mulAB_n1( out, gcl, stv );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &st ); 

  return( ret );
}
