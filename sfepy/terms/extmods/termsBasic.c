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
#define __FUNC__ "dq_grad"
int32 dq_grad( FMField *out, FMField *state, int32 offset,
	       VolumeGeometry *vg, int32 *conn, int32 nEl, int32 nEP )
{
  int32 ii, nQP, ret = RET_OK;
  FMField *st = 0;

  state->val = FMF_PtrFirst( state ) + offset;

  nQP = vg->bfGM->nLev;

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
#define __FUNC__ "de_grad"
int32 de_grad( FMField *out, FMField *state, int32 offset,
	       VolumeGeometry *vg, int32 *conn, int32 nEl, int32 nEP,
	       int32 *elList, int32 elList_nRow )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *st = 0, *out_qp = 0;

  state->val = FMF_PtrFirst( state ) + offset;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  fmf_createAlloc( &st, 1, 1, nEP, out->nCol );
  fmf_createAlloc( &out_qp, 1, nQP, dim, out->nCol );

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];
    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( vg->det, iel );

    ele_extractNodalValuesNBN( st, state, conn + nEP * iel );
    fmf_mulAB_n1( out_qp, vg->bfGM, st );
    fmf_sumLevelsMulF( out, out_qp, vg->det->val );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &st );
  fmf_freeDestroy( &out_qp );

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

#undef __FUNC__
#define __FUNC__ "dw_volume_wdot_scalar"
/*!
  @par Revision history:
  - 04.08.2006, c
  - 01.12.2006
*/
int32 dw_volume_wdot_scalar( FMField *out, float64 coef, FMField *state_qp,
			     FMField *bf, FMField *mtxD, VolumeGeometry *vg,
			     int32 *elList, int32 elList_nRow,
			     int32 isDiff )
{
  int32 ii, iel, nQP, nEP, ret = RET_OK;
  FMField *ftd = 0, *ftdf = 0, *dfp = 0, *ftdfp = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;

  if (isDiff) {
    fmf_createAlloc( &ftd, 1, nQP, nEP, 1 );
    fmf_createAlloc( &ftdf, 1, nQP, nEP, nEP );
  } else {
    fmf_createAlloc( &dfp, 1, nQP, 1, 1 );
    fmf_createAlloc( &ftdfp, 1, nQP, nEP, 1 );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg->det, iel );

    if (isDiff) {
      fmf_mulATB_nn( ftd, bf, mtxD );
      fmf_mulAB_nn( ftdf, ftd, bf );
      fmf_sumLevelsMulF( out, ftdf, vg->det->val );
    } else {
      FMF_SetCell( state_qp, iel );
      fmf_mulAB_nn( dfp, mtxD, state_qp );
      fmf_mulATB_nn( ftdfp, bf, dfp );
      fmf_sumLevelsMulF( out, ftdfp, vg->det->val );
    }
    ERR_CheckGo( ret );
  }

  // E.g. 1/dt.
  if (isDiff == 2) {
    coef = -coef;
  }
  fmfc_mulC( out, coef );

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftd ); 
    fmf_freeDestroy( &ftdf ); 
  } else {
    fmf_freeDestroy( &dfp ); 
    fmf_freeDestroy( &ftdfp ); 
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_volume_surface"

int32 d_volume_surface( FMField *out, FMField *in,
			FMField *bf, SurfaceGeometry *sg,
			int32 *conn, int32 nEl, int32 nEP,
			int32 *elList, int32 elList_nRow )
{
  int32 ii, iel, dim, nQP, nFP, ret = RET_OK;
  FMField *lcoor, *aux, *aux2;
  float64 val;

  nFP = bf->nCol;
  nQP = sg->det->nLev;
  dim = sg->normal->nRow;
  val = 1.0/dim;

  fmf_createAlloc( &lcoor, 1, 1, nFP, dim );
  fmf_createAlloc( &aux, 1, nQP, 1, dim );
  fmf_createAlloc( &aux2, 1, nQP, 1, 1 );

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( sg->normal, iel );
    FMF_SetCell( sg->det, iel );

    ele_extractNodalValuesNBN( lcoor, in, conn + nEP * iel );
    fmf_mulAB_n1( aux, bf, lcoor );
    fmf_mulAB_nn( aux2, aux, sg->normal );    
    fmf_sumLevelsMulF( out, aux2, sg->det->val );
    fmf_mulC( out, val );
  }  /* for (ii) */

 end_label:
  fmf_freeDestroy( &aux );
  fmf_freeDestroy( &aux2 );
  fmf_freeDestroy( &lcoor );

  return( ret );
}
