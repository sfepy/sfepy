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
int32 dw_surf_laplace( FMField *out, FMField *grad, FMField *coef,
		       FMField *gbf, SurfaceGeometry *sg,
		       int32 isDiff )
{
  int32 ii, nEP, dim, nQP, ret = RET_OK;
  FMField *aux1 = 0, *aux2 = 0;

  nQP = gbf->nLev;
  dim = gbf->nRow;
  nEP = gbf->nCol;

  fmf_createAlloc( &aux1, 1, nQP, nEP, dim );

  if (isDiff)
    fmf_createAlloc( &aux2, 1, nQP, nEP, nEP );
  else
    fmf_createAlloc( &aux2, 1, nQP, nEP, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    if (coef->nCell > 1)
      FMF_SetCell( coef, ii );

    FMF_SetCell( out, ii );
    FMF_SetCell( sg->det, ii );

    fmf_mulATB_nn(aux1, gbf, coef);

    if (isDiff) {
      fmf_mulAB_nn(aux2, aux1, gbf);
    }
    else {
      FMF_SetCell(grad, ii);
      fmf_mulAB_nn(aux2, aux1, grad);
    }

    fmf_sumLevelsMulF(out, aux2, sg->det->val);

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy(&aux1);
  fmf_freeDestroy(&aux2);

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_surf_laplace"
int32 d_surf_laplace( FMField *out, FMField *gradP, FMField *gradQ, FMField *coef,
		      SurfaceGeometry *sg )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *aux1 = 0, *aux2 = 0;

  nQP = gradP->nLev;
  dim = coef->nRow;

  fmf_createAlloc( &aux1, 1, nQP, dim, 1 );
  fmf_createAlloc( &aux2, 1, nQP, 1, 1 );

  for (ii = 0; ii < out->nCell; ii++) {

    if (coef->nCell > 1)
      FMF_SetCell( coef, ii );

    FMF_SetCell( out, ii );
    FMF_SetCell( sg->det, ii );
    FMF_SetCell( gradP, ii );
    FMF_SetCell( gradQ, ii );

    fmf_mulAB_nn(aux1, coef, gradP);
    fmf_mulATB_nn(aux2, gradQ, aux1);
    fmf_sumLevelsMulF(out, aux2, sg->det->val);

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &aux1 );
  fmf_freeDestroy( &aux2 );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_surf_lcouple"
int32 dw_surf_lcouple( FMField *out, FMField *state, FMField *coef,
		       FMField *bf, FMField *gbf, SurfaceGeometry *sg,
		       int32 isDiff )
{
  int32 ii, nEP, dim, nQP, ret = RET_OK;
  FMField *aux1 = 0, *aux2 = 0;

  nQP = gbf->nLev;
  nEP = gbf->nCol;
  dim = coef->nCol;

  fmf_createAlloc(&aux1, 1, nQP, nEP, dim);

  if (isDiff)
    fmf_createAlloc(&aux2, 1, nQP, nEP, nEP);
  else
    fmf_createAlloc(&aux2, 1, nQP, nEP, 1);

  for (ii = 0; ii < out->nCell; ii++) {
    if (coef->nCell > 1)
      FMF_SetCell( coef, ii );

    FMF_SetCell( out, ii );
    FMF_SetCell( sg->det, ii );

    fmf_mulATB_nn(aux1, gbf, coef);

    if (isDiff) {
      fmf_mulAB_nn(aux2, aux1, bf);
    }
    else {
      FMF_SetCell(state, ii);
      fmf_mulAB_nn(aux2, aux1, state);
    }

    fmf_sumLevelsMulF(out, aux2, sg->det->val);

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy(&aux1);
  fmf_freeDestroy(&aux2);

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_surf_lcouple"
int32 d_surf_lcouple( FMField *out, FMField *stateP, FMField *gradQ, FMField *coef,
		      SurfaceGeometry *sg )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *aux1 = 0, *aux2 = 0;

  nQP = stateP->nLev;
  dim = coef->nRow;

  fmf_createAlloc(&aux1, 1, nQP, dim, 1);
  fmf_createAlloc(&aux2, 1, nQP, 1, 1);

  for (ii = 0; ii < out->nCell; ii++) {
    if (coef->nCell > 1)
      FMF_SetCell( coef, ii );

    FMF_SetCell( out, ii );
    FMF_SetCell( sg->det, ii );
    FMF_SetCell(stateP, ii);
    FMF_SetCell(gradQ, ii);

    fmf_mulAB_nn(aux1, coef, stateP);
    fmf_mulATB_nn(aux2, gradQ, aux1);
    fmf_sumLevelsMulF(out, aux2, sg->det->val);

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy(&aux1);
  fmf_freeDestroy(&aux2);

  return( ret );
}
