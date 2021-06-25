#include "terms_electric.h"
#include "terms.h"
#include "geommech.h"

#undef __FUNC__
#define __FUNC__ "dw_electric_source"
int32 dw_electric_source( FMField *out, FMField *grad, FMField *coef,
			  Mapping *vg )
{
  int32 ii, nQP, nEP, ret = RET_OK;
  FMField *gp2 = 0, *bftgp2 = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bf->nCol;

  fmf_createAlloc( &gp2, 1, nQP, 1, 1 );
  fmf_createAlloc( &bftgp2, 1, nQP, nEP, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCellX1( coef, ii );
    FMF_SetCell( grad, ii );
    FMF_SetCellX1( vg->bf, ii );

    fmf_mulATB_nn( gp2, grad, grad );
    fmf_mulATB_nn( bftgp2, vg->bf, gp2 );
    fmf_sumLevelsMulF( out, bftgp2, vg->det->val );
    fmf_mulC( out, coef->val[0] );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &gp2 );
  fmf_freeDestroy( &bftgp2 );

  return( ret );
}
