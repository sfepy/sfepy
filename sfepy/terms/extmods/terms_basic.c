#include "terms_basic.h"
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
    FMF_SetCellX1( bf, ii );

    ele_extractNodalValuesDBD( st, state, conn + nEP * ii );
    bf_act( out, bf, st );
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &st );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_grad"
int32 dq_grad( FMField *out, FMField *state, int32 offset,
	       Mapping *vg, int32 *conn, int32 nEl, int32 nEP )
{
  int32 ii, ret = RET_OK;
  FMField *st = 0;

  state->val = FMF_PtrFirst( state ) + offset;

  fmf_createAlloc( &st, 1, 1, nEP, out->nCol );

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
		     Mapping *vg,
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

#undef __FUNC__
#define __FUNC__ "d_volume_surface"
int32 d_volume_surface( FMField *out, FMField *in,
			Mapping *sg,
			int32 *conn, int32 nEl, int32 nEP )
{
  int32 ii, dim, nQP, nFP, ret = RET_OK;
  FMField *lcoor, *aux, *aux2;
  float64 val;

  nFP = sg->bf->nCol;
  nQP = sg->det->nLev;
  dim = sg->normal->nRow;
  val = 1.0/dim;

  fmf_createAlloc( &lcoor, 1, 1, nFP, dim );
  fmf_createAlloc( &aux, 1, nQP, 1, dim );
  fmf_createAlloc( &aux2, 1, nQP, 1, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( sg->normal, ii );
    FMF_SetCell( sg->det, ii );
    FMF_SetCellX1( sg->bf, ii );

    ele_extractNodalValuesNBN( lcoor, in, conn + nEP * ii );
    fmf_mulAB_n1( aux, sg->bf, lcoor );
    fmf_mulAB_nn( aux2, aux, sg->normal );
    fmf_sumLevelsMulF( out, aux2, sg->det->val );
    fmf_mulC( out, val );
    ERR_CheckGo( ret );
  }  /* for (ii) */

 end_label:
  fmf_freeDestroy( &aux );
  fmf_freeDestroy( &aux2 );
  fmf_freeDestroy( &lcoor );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "di_surface_moment"

int32 di_surface_moment( FMField *out, FMField *in,
			 Mapping *sg,
			 int32 *conn, int32 nEl, int32 nEP )
{
  int32 ii, dim, nQP, nFP, ret = RET_OK;
  FMField *lcoor, *aux, *aux2;

  nFP = sg->bf->nCol;
  nQP = sg->det->nLev;
  dim = sg->normal->nRow;

  fmf_createAlloc( &lcoor, 1, 1, nFP, dim );
  fmf_createAlloc( &aux, 1, nQP, 1, dim );
  fmf_createAlloc( &aux2, 1, nQP, dim, dim );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( sg->normal, ii );
    FMF_SetCell( sg->det, ii );
    FMF_SetCellX1( sg->bf, ii );

    ele_extractNodalValuesNBN( lcoor, in, conn + nEP * ii );
    fmf_mulAB_n1( aux, sg->bf, lcoor );
    fmf_mulAB_nn( aux2, sg->normal, aux );
    fmf_sumLevelsMulF( out, aux2, sg->det->val );

    ERR_CheckGo( ret );
  }  /* for (ii) */

 end_label:
  fmf_freeDestroy( &aux );
  fmf_freeDestroy( &aux2 );
  fmf_freeDestroy( &lcoor );

  return( ret );
}
