#include "terms_biot.h"
#include "terms.h"
#include "form_sdcc.h"

#undef __FUNC__
#define __FUNC__ "dw_biot_grad"
/*!
  @par Revision history:
  - 03.08.2006, c
  - 01.12.2006
*/
int32 dw_biot_grad( FMField *out, float64 coef, FMField *pressure_qp,
		    FMField *mtxD, Mapping *svg, Mapping *vvg,
		    int32 isDiff )
{
  int32 ii, nEPU, nEP, dim, nQP, ret = RET_OK;
  FMField *dfp = 0, *gtdfp = 0, *gtd = 0, *gtdf = 0;

  nQP = vvg->bfGM->nLev;
  dim = vvg->bfGM->nRow;
  nEPU = vvg->bfGM->nCol;
  nEP = svg->bf->nCol;

/*   fmf_print( mtxD, stdout, 0 ); */

  if (isDiff == 1) {
    fmf_createAlloc( &gtd, 1, nQP, dim * nEPU, 1 );
    fmf_createAlloc( &gtdf, 1, nQP, dim * nEPU, nEP );
  } else {
    int32 sym = (dim + 1) * dim / 2;
    fmf_createAlloc( &dfp, 1, nQP, sym, 1 );
    fmf_createAlloc( &gtdfp, 1, nQP, dim * nEPU, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( mtxD, ii );
    FMF_SetCell( vvg->bfGM, ii );
    FMF_SetCell( vvg->det, ii );

    if (isDiff == 1) {
      FMF_SetCellX1( svg->bf, ii );
      form_sdcc_actOpGT_M3( gtd, vvg->bfGM, mtxD );
      fmf_mulAB_nn( gtdf, gtd, svg->bf );
      fmf_sumLevelsMulF( out, gtdf, vvg->det->val );
    } else {
      FMF_SetCell( pressure_qp, ii );
      fmf_mulAB_nn( dfp, mtxD, pressure_qp );
      form_sdcc_actOpGT_VS3( gtdfp, vvg->bfGM, dfp );
      fmf_sumLevelsMulF( out, gtdfp, vvg->det->val );
    }
    ERR_CheckGo( ret );
  }

  // E.g. 1/dt.
  fmfc_mulC( out, coef );

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &gtd );
    fmf_freeDestroy( &gtdf );
  } else {
    fmf_freeDestroy( &dfp );
    fmf_freeDestroy( &gtdfp );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_biot_div"
/*!
  @par Revision history:
  - 04.08.2006, c
  - 01.12.2006
*/
int32 dw_biot_div( FMField *out, float64 coef, FMField *strain,
		   FMField *mtxD, Mapping *svg, Mapping *vvg,
                   int32 isDiff )
{
  int32 ii, nEPP, nEP, dim, sym, nQP, ret = RET_OK;
  FMField *dtg = 0, *ftdtg = 0, *dtgu = 0, *ftdtgu = 0;
  FMField drow[1];

  nQP = vvg->bfGM->nLev;
  nEP = vvg->bfGM->nCol;
  dim = vvg->bfGM->nRow;
  sym = (dim + 1) * dim / 2;
  nEPP = svg->bf->nCol;

/*   fmf_print( mtxD, stdout, 0 ); */

  if (isDiff == 1) {
    fmf_createAlloc( &dtg, 1, nQP, 1, dim * nEP );
    fmf_createAlloc( &ftdtg, 1, nQP, nEPP, dim * nEP );

    drow->nAlloc = -1;
    fmf_pretend( drow, 1, nQP, 1, sym, mtxD->val );
  } else {
    fmf_createAlloc( &dtgu, 1, nQP, 1, 1 );
    fmf_createAlloc( &ftdtgu, 1, nQP, nEPP, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( mtxD, ii );
    FMF_SetCell( vvg->bfGM, ii );
    FMF_SetCell( vvg->det, ii );
    FMF_SetCellX1( svg->bf, ii );

    if (isDiff == 1) {
      drow->val = mtxD->val;
      form_sdcc_actOpG_RM3( dtg, drow, vvg->bfGM );
      fmf_mulATB_nn( ftdtg, svg->bf, dtg );
      fmf_sumLevelsMulF( out, ftdtg, vvg->det->val );
    } else {
      FMF_SetCell( strain, ii );
      fmf_mulATB_nn( dtgu, mtxD, strain );
      fmf_mulATB_nn( ftdtgu, svg->bf, dtgu );
      fmf_sumLevelsMulF( out, ftdtgu, vvg->det->val );
    }
    ERR_CheckGo( ret );
  }

  // E.g. 1/dt.
  fmfc_mulC( out, coef );

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &dtg );
    fmf_freeDestroy( &ftdtg );
  } else {
    fmf_freeDestroy( &dtgu );
    fmf_freeDestroy( &ftdtgu );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_biot_div"
/*!
  @par Revision history:
  - c: 05.03.2008, r: 05.03.2008
*/
int32 d_biot_div( FMField *out, float64 coef, FMField *state, FMField *strain,
		  FMField *mtxD, Mapping *vg )
{
  int32 ii, nQP, ret = RET_OK;
  FMField *dtgu = 0, *pftdtgu = 0;

  nQP = vg->bfGM->nLev;

/*   fmf_print( mtxD, stdout, 0 ); */

  fmf_createAlloc( &dtgu, 1, nQP, 1, 1 );
  fmf_createAlloc( &pftdtgu, 1, nQP, 1, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( mtxD, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCell( state, ii );
    FMF_SetCell( strain, ii );

    fmf_mulATB_nn( dtgu, mtxD, strain );
    fmf_mulATB_nn( pftdtgu, state, dtgu );
    fmf_sumLevelsMulF( out, pftdtgu, vg->det->val );
    ERR_CheckGo( ret );
  }

  // E.g. 1/dt.
  fmfc_mulC( out, coef );

 end_label:
  fmf_freeDestroy( &dtgu );
  fmf_freeDestroy( &pftdtgu );

  return( ret );
}
