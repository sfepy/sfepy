#include "termsElectric.h"
#include "terms.h"
#include "geommech.h"

#undef __FUNC__
#define __FUNC__ "dw_electric_source"
int32 dw_electric_source( FMField *out,
			  FMField *state,
			  FMField *coef, FMField *bf,
			  VolumeGeometry *vgc,
			  int32 *conn, int32 nEl, int32 nEP,
			  int32 *elList, int32 elList_nRow,
			  int32 mode )
{
  int32 ii, iel, dim, nQP, nEPR, nEPC, ret = RET_OK;
  FMField *p = 0, *gp = 0, *gp2 = 0, *bftgp2 = 0;

  nQP = vgc->bfGM->nLev;
  dim = vgc->bfGM->nRow;
  nEPC = vgc->bfGM->nCol;
  nEPR = bf->nCol;

  state->val = FMF_PtrFirst( state );

  fmf_createAlloc( &p, 1, 1, nEPC, 1 );
  fmf_createAlloc( &gp, 1, nQP, dim, 1 );
  fmf_createAlloc( &gp2, 1, nQP, 1, 1 );
  fmf_createAlloc( &bftgp2, 1, nQP, nEPR, 1 );

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vgc->bfGM, iel );
    FMF_SetCell( vgc->det, iel );
    FMF_SetCell( coef, ii );

    ele_extractNodalValuesNBN( p, state, conn + nEPC * iel );
    fmf_mulAB_n1( gp, vgc->bfGM, p );
    fmf_mulATB_nn( gp2, gp, gp );
    fmf_mulATB_nn( bftgp2, bf, gp2 );
    fmf_sumLevelsMulF( out, bftgp2, vgc->det->val );
    fmf_mulC( out, coef->val[0] );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &p );
  fmf_freeDestroy( &gp );
  fmf_freeDestroy( &gp2 ); 
  fmf_freeDestroy( &bftgp2 );

  return( ret );
}
