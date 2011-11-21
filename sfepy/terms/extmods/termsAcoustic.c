#include "termsAcoustic.h"
#include "terms.h"

#undef __FUNC__
#define __FUNC__ "d_diffusion_sa"
int32 d_diffusion_sa( FMField *out,
		      FMField *grad_q, FMField *grad_p,
		      FMField *grad_w, FMField *div_w,
		      FMField *mtxD, VolumeGeometry *vg )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *aux2 = 0, *aux3 = 0, *aux4 = 0, *out0 = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  FMF_SetFirst( out );

  fmf_createAlloc( &aux2, 1, nQP, dim, 1 );
  fmf_createAlloc( &aux3, 1, nQP, 1, 1 );
  fmf_createAlloc( &aux4, 1, nQP, dim, 1 );
  fmf_createAlloc( &out0, 1, nQP, 1, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCell( mtxD, ii );
    FMF_SetCell( grad_q, ii );
    FMF_SetCell( grad_p, ii );
    FMF_SetCell( grad_w, ii );
    FMF_SetCell( div_w, ii );

    /* div w K_ij grad_j q grad_i p */
    fmf_mulAB_nn( aux2, mtxD, grad_p );
    fmf_mulATB_nn( aux3, grad_q, aux2 );
    fmf_mulAB_nn( out0, div_w, aux3 );

    /* grad_k q K_ij grad_j w_k grad_i p */
    fmf_mulATB_nn( aux4, grad_w, aux2 );
    fmf_mulATB_nn( aux3, grad_q, aux4 );
    fmf_subAB_nn( out0, out0, aux3 );

    /* grad_k q K_ij grad_j w_k grad_i p */
    fmf_mulAB_nn( aux2, grad_w, grad_p );
    fmf_mulAB_nn( aux4, mtxD, aux2 );
    fmf_mulATB_nn( aux3, grad_q, aux4 );
    fmf_subAB_nn( out0, out0, aux3 );

    fmf_sumLevelsMulF( out, out0, vg->det->val );

    FMF_SetCellNext( out );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &out0 );
  fmf_freeDestroy( &aux2 );
  fmf_freeDestroy( &aux3 );
  fmf_freeDestroy( &aux4 );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_surf_laplace"
int32 dw_surf_laplace( FMField *out, FMField *state, FMField *coef,
		       FMField *gbf, SurfaceGeometry *sg,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow,
		       int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *bb = 0, *bb0 = 0, *aux = 0, *det0 = 0, *gbf0 = 0;
  FMField *st = 0;

  nQP = gbf->nLev;
  dim = gbf->nRow;

  fmf_createAlloc( &bb, 1, nQP, nEP, nEP );
  fmf_createAlloc( &aux, 1, nQP, dim, nEP );
  fmf_createAlloc( &gbf0, 1, nQP, dim, nEP );
  fmf_createAlloc( &det0, 1, nQP, 1, 1 );

  fmf_fillC(det0, 1.0 / nQP);
  for (ii = 0; ii < det0->nLev; ii++)
    det0->val[ii] /= sg->det->val[ii];
  fmf_mulAF(gbf0, gbf, det0->val);

  if (!isDiff) {
    state->val = FMF_PtrFirst( state );
    fmf_createAlloc( &bb0, 1, 1, nEP, nEP );
    fmf_createAlloc( &st, 1, 1, nEP, 1 );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    if (coef->nCell > 1) {
      FMF_SetCell( coef, ii );
    }

    FMF_SetCell( out, ii );
    FMF_SetCell( sg->det, iel );

    fmf_mulAB_nn( aux, coef, gbf0 );
    fmf_mulATB_nn( bb, gbf0, aux );

    if (isDiff) {
      fmf_sumLevelsMulF( out, bb, sg->det->val );
    }
    else {
      ele_extractNodalValuesNBN( st, state, conn + nEP * iel );
      fmf_sumLevelsMulF( bb0, bb, sg->det->val );
      fmf_mulAB_nn( out, bb0, st );
    }

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &bb );
  fmf_freeDestroy( &aux );
  fmf_freeDestroy(&gbf0);
  fmf_freeDestroy(&det0);
  if (!isDiff) {
    fmf_freeDestroy( &bb0 );
    fmf_freeDestroy( &st );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_surf_laplace"
int32 d_surf_laplace( FMField *out, FMField *stateP, FMField *stateQ, FMField *coef,
		      FMField *gbf, SurfaceGeometry *sg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow)
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *aux1 = 0, *aux2 = 0, *aux3 = 0, *det0 = 0, *gbf0;
  FMField *st = 0;

  nQP = gbf->nLev;
  dim = gbf->nRow;

  fmf_createAlloc( &aux1, 1, nQP, dim, 1 );
  fmf_createAlloc( &aux2, 1, nQP, dim, 1 );
  fmf_createAlloc( &aux3, 1, nQP, 1, 1 );
  fmf_createAlloc( &st, 1, 1, nEP, 1 );
  fmf_createAlloc( &gbf0, 1, nQP, dim, nEP );
  fmf_createAlloc( &det0, 1, nQP, 1, 1 );

  fmf_fillC(det0, 1.0 / nQP);
  for (ii = 0; ii < det0->nLev; ii++)
    det0->val[ii] /= sg->det->val[ii];
  fmf_mulAF(gbf0, gbf, det0->val);

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    if (coef->nCell > 1) {
      FMF_SetCell( coef, ii );
    }

    FMF_SetCell( out, ii );
    FMF_SetCell( sg->det, iel );

    ele_extractNodalValuesNBN(st, stateP, conn + nEP * iel);
    fmf_mulAB_n1(aux1, gbf0, st);
    fmf_mulAB_nn(aux2, coef, aux1);
    ele_extractNodalValuesNBN(st, stateQ, conn + nEP * iel);
    fmf_mulAB_n1(aux1, gbf0, st);
    fmf_mulATB_nn(aux3, aux1, aux2);
    fmf_sumLevelsMulF(out, aux3, sg->det->val);

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &aux1 );
  fmf_freeDestroy( &aux2 );
  fmf_freeDestroy( &aux3 );
  fmf_freeDestroy( &st );
  fmf_freeDestroy(&gbf0);
  fmf_freeDestroy(&det0);

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_surf_lcouple"
int32 dw_surf_lcouple( FMField *out, FMField *state, FMField *coef,
		       FMField *bf, FMField *gbf, SurfaceGeometry *sg,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow,
		       int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *bb = 0, *bb0 = 0, *aux = 0, *st = 0, *det0 = 0, *gbf0 = 0;

  nQP = gbf->nLev;
  dim = gbf->nRow;

  fmf_createAlloc( &bb, 1, nQP, nEP, nEP );
  fmf_createAlloc( &aux, 1, nQP, dim, nEP );
  fmf_createAlloc( &gbf0, 1, nQP, dim, nEP );
  fmf_createAlloc( &det0, 1, nQP, 1, 1 );

  fmf_fillC(det0, 1.0 / nQP);
  for (ii = 0; ii < det0->nLev; ii++)
    det0->val[ii] /= sg->det->val[ii];
  fmf_mulAF(gbf0, gbf, det0->val);

  if (!isDiff) {
    state->val = FMF_PtrFirst( state );
    fmf_createAlloc( &bb0, 1, 1, nEP, nEP );
    fmf_createAlloc( &st, 1, 1, nEP, 1 );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    if (coef->nCell > 1) {
      FMF_SetCell( coef, ii );
    }

    FMF_SetCell( out, ii );
    FMF_SetCell( sg->det, iel );

    fmf_mulAB_nn( aux, coef, bf );
    fmf_mulATB_nn( bb, gbf0, aux );

    if (isDiff) {
      fmf_sumLevelsMulF( out, bb, sg->det->val );
    }
    else {
      ele_extractNodalValuesNBN( st, state, conn + nEP * iel );
      fmf_sumLevelsMulF( bb0, bb, sg->det->val );
      fmf_mulAB_nn( out, bb0, st );
    }

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &bb );
  fmf_freeDestroy( &aux );
  fmf_freeDestroy(&gbf0);
  fmf_freeDestroy(&det0);
  if (!isDiff) {
    fmf_freeDestroy( &bb0 );
    fmf_freeDestroy( &st );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_surf_lcouple"
int32 d_surf_lcouple( FMField *out, FMField *stateP, FMField *stateQ, FMField *coef,
		       FMField *bf, FMField *gbf, SurfaceGeometry *sg,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow)
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *aux1 = 0, *aux2 = 0, *aux3 = 0, *aux4 = 0, *st = 0, *det0 = 0, *gbf0 = 0;

  nQP = gbf->nLev;
  dim = gbf->nRow;

  fmf_createAlloc(&aux1, 1, nQP, dim, 1);
  fmf_createAlloc(&aux2, 1, nQP, 1, 1);
  fmf_createAlloc(&aux3, 1, nQP, 1, 1);
  fmf_createAlloc(&aux4, 1, nQP, 1, 1);
  fmf_createAlloc(&st, 1, 1, nEP, 1);
  fmf_createAlloc( &gbf0, 1, nQP, dim, nEP );
  fmf_createAlloc( &det0, 1, nQP, 1, 1 );

  fmf_fillC(det0, 1.0 / nQP);
  for (ii = 0; ii < det0->nLev; ii++)
    det0->val[ii] /= sg->det->val[ii];
  fmf_mulAF(gbf0, gbf, det0->val);

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    if (coef->nCell > 1) {
      FMF_SetCell( coef, ii );
    }

    FMF_SetCell( out, ii );
    FMF_SetCell( sg->det, iel );

    ele_extractNodalValuesNBN(st, stateP, conn + nEP * iel);
    fmf_mulAB_n1(aux1, gbf0, st);
    fmf_mulAB_nn(aux2, coef, aux1);
    ele_extractNodalValuesNBN(st, stateQ, conn + nEP * iel);
    fmf_mulAB_n1(aux3, bf, st);
    fmf_mulATB_nn(aux4, aux3, aux2);
    fmf_sumLevelsMulF(out, aux4, sg->det->val);

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy(&aux1);
  fmf_freeDestroy(&aux2);
  fmf_freeDestroy(&aux3);
  fmf_freeDestroy(&aux4);
  fmf_freeDestroy(&st);
  fmf_freeDestroy(&gbf0);
  fmf_freeDestroy(&det0);

  return( ret );
}
