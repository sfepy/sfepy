#include "terms_piezo.h"
#include "terms.h"
#include "form_sdcc.h"

#undef __FUNC__
#define __FUNC__ "dw_piezo_coupling"
int32 dw_piezo_coupling( FMField *out, FMField *strain, FMField *charge_grad,
			 FMField *mtxG, Mapping *vg,
			 int32 mode )
{
  int32 ii, nEPU, nEPP, dim, sym, nQP, ret = RET_OK;
  FMField *gtgp = 0, *gtgtgp = 0, *ge = 0, *gctge = 0, *gg = 0, *gctgg = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;
  nEPU = vg->bfGM->nCol;
  sym = (dim + 1) * dim / 2;

  if ((mode == 0) || (mode == 1)) {
    nEPP = out->nCol;
  } else {
    nEPP = out->nRow;
  }

  if (mode == 0) { 
    fmf_createAlloc( &gtgp, 1, nQP, sym, 1 );
    fmf_createAlloc( &gtgtgp, 1, nQP, dim * nEPU, 1 );
  } else if (mode == 2) {
    fmf_createAlloc( &ge, 1, nQP, dim, 1 );
    fmf_createAlloc( &gctge, 1, nQP, nEPP, 1 );
  } else {
    fmf_createAlloc( &gg, 1, nQP, dim, dim * nEPU );
    fmf_createAlloc( &gctgg, 1, nQP, nEPP, dim * nEPU );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCellX1( mtxG, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );

    if (mode == 0) { // vector - grad
      FMF_SetCell( charge_grad, ii );

      fmf_mulATB_nn( gtgp, mtxG, charge_grad );
      form_sdcc_actOpGT_VS3( gtgtgp, vg->bfGM, gtgp );
      fmf_sumLevelsMulF( out, gtgtgp, vg->det->val );

    } else if (mode == 2) { // vector - div
      FMF_SetCell( strain, ii );

      fmf_mulAB_nn( ge, mtxG, strain );
      fmf_mulATB_nn( gctge, vg->bfGM, ge );
      fmf_sumLevelsMulF( out, gctge, vg->det->val );

    } else { // matrix - div, grad
      form_sdcc_actOpG_RM3( gg, mtxG, vg->bfGM );
      fmf_mulATB_nn( gctgg, vg->bfGM, gg );
      if (mode == 1) { // matrix - grad
	fmf_sumLevelsTMulF( out, gctgg, vg->det->val );
      } else { // matrix - div
	fmf_sumLevelsMulF( out, gctgg, vg->det->val );
      }
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (mode == 0) {
    fmf_freeDestroy( &gtgp );
    fmf_freeDestroy( &gtgtgp );
  } else if (mode == 2) {
    fmf_freeDestroy( &ge );
    fmf_freeDestroy( &gctge );
  } else {
    fmf_freeDestroy( &gg );
    fmf_freeDestroy( &gctgg );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_piezo_coupling"
int32 d_piezo_coupling( FMField *out, FMField *strain, FMField *charge_grad,
			FMField *mtxG, Mapping *vg )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *ge = 0, *gptge = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  fmf_createAlloc( &ge, 1, nQP, dim, 1 );
  fmf_createAlloc( &gptge, 1, nQP, 1, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCellX1( mtxG, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCell( strain, ii );
    FMF_SetCell( charge_grad, ii );

    fmf_mulAB_nn( ge, mtxG, strain );
    fmf_mulATB_nn( gptge, charge_grad, ge );
    fmf_sumLevelsMulF( out, gptge, vg->det->val );
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &ge );
  fmf_freeDestroy( &gptge );

  return( ret );
}
