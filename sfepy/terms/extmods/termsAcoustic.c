#include "termsAcoustic.h"
#include "terms.h"

#undef __FUNC__
#define __FUNC__ "d_diffusion_sa"
int32 d_diffusion_sa( FMField *out,
		      FMField *stateQ, FMField *stateP, FMField *stateW,
		      FMField *mtxD,
		      VolumeGeometry *vg, VolumeGeometry *vg_w,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *conn_w, int32 nEl_w, int32 nEP_w,
		      int32 *elList, int32 elList_nRow )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *nodval = 0, *nodval_w = 0, *nodval_w0 = 0;
  FMField *gd_q = 0, *gd_p = 0, *gd_w = 0, *div_w = 0;
  FMField *aux2 = 0, *aux3 = 0, *aux4 = 0, *out0 = 0;
  FMField divop[1], nodval_wv[1];

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  FMF_SetFirst( out );

  fmf_createAlloc( &div_w, 1, nQP, 1, 1 );
  fmf_createAlloc( &gd_w, 1, nQP, dim, dim );
  fmf_createAlloc( &gd_q, 1, nQP, dim, 1 );
  fmf_createAlloc( &gd_p, 1, nQP, dim, 1 );
  fmf_createAlloc( &aux2, 1, nQP, dim, 1 );
  fmf_createAlloc( &aux3, 1, nQP, 1, 1 );
  fmf_createAlloc( &aux4, 1, nQP, dim, 1 );
  fmf_createAlloc( &out0, 1, nQP, 1, 1 );
  fmf_createAlloc( &nodval, 1, 1, nEP, 1 );
  fmf_createAlloc( &nodval_w, 1, 1, nEP_w, dim );
  fmf_createAlloc( &nodval_w0, 1, 1, dim, nEP_w );

  divop->nAlloc = -1;
  fmf_pretend( divop, 1, nQP, 1, nEP_w * dim, vg_w->bfGM->val0 );
  nodval_wv->nAlloc = -1;
  fmf_pretend( nodval_wv, 1, 1, nEP_w * dim, 1, nodval_w0->val );

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( divop, iel );
    FMF_SetCell( vg->det, iel );

    /* grad p */
    ele_extractNodalValuesNBN( nodval, stateP, conn + nEP * iel );
    fmf_mulAB_n1( gd_p, vg->bfGM, nodval );
    /* grad q */
    ele_extractNodalValuesNBN( nodval, stateQ, conn + nEP * iel );
    fmf_mulAB_n1( gd_q, vg->bfGM, nodval );
    /* w */
    ele_extractNodalValuesNBN( nodval_w, stateW, conn_w + nEP_w * iel );
    ele_extractNodalValuesDBD( nodval_w0, stateW, conn_w + nEP_w * iel );
    /* grad w */
    fmf_mulAB_n1( gd_w, vg->bfGM, nodval_w );
    /* div w */
    fmf_mulAB_n1( div_w, divop, nodval_wv );

    /* div w K_ij grad_j q grad_i p */
    fmf_mulAB_nn( aux2, mtxD, gd_p );
    fmf_mulATB_nn( aux3, gd_q, aux2 );
    fmf_mulAB_nn( out0, div_w, aux3 );

    /* grad_k q K_ij grad_j w_k grad_i p */
    fmf_mulATB_nn( aux4, gd_w, aux2 );
    fmf_mulATB_nn( aux3, gd_q, aux4 );
    fmf_subAB_nn( out0, out0, aux3 );

    /* grad_k q K_ij grad_j w_k grad_i p */
    fmf_mulAB_nn( aux2, gd_w, gd_p );
    fmf_mulAB_nn( aux4, mtxD, aux2 );
    fmf_mulATB_nn( aux3, gd_q, aux4 );
    fmf_subAB_nn( out0, out0, aux3 );

    fmf_sumLevelsMulF( out, out0, vg->det->val );

    FMF_SetCellNext( out );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &out0 );
  fmf_freeDestroy( &gd_q );
  fmf_freeDestroy( &gd_p );
  fmf_freeDestroy( &gd_w );
  fmf_freeDestroy( &div_w );
  fmf_freeDestroy( &aux2 );
  fmf_freeDestroy( &aux3 );
  fmf_freeDestroy( &aux4 );
  fmf_freeDestroy( &nodval );
  fmf_freeDestroy( &nodval_w );
  fmf_freeDestroy( &nodval_w0 );

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
