#include "termsBiot.h"
#include "terms.h"
#include "formSDCC.h"

#undef __FUNC__
#define __FUNC__ "dw_biot_grad"
/*!
  @par Revision history:
  - 03.08.2006, c
  - 01.12.2006
*/
int32 dw_biot_grad( FMField *out, float64 coef, FMField *pressure_qp,
		    FMField *bf, FMField *mtxD, VolumeGeometry *vg,
		    int32 *elList, int32 elList_nRow,
		    int32 isDiff )
{
  int32 ii, iel, nEPU, nEP, dim, nQP, ret = RET_OK;
  FMField *dfp = 0, *gtdfp = 0, *gtd = 0, *gtdf = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;
  nEPU = vg->bfGM->nCol;
  nEP = bf->nCol;

/*   fmf_print( mtxD, stdout, 0 ); */

  if (isDiff == 1) { 
    fmf_createAlloc( &gtd, 1, nQP, dim * nEPU, 1 );
    fmf_createAlloc( &gtdf, 1, nQP, dim * nEPU, nEP );
  } else {
    int32 sym = (dim + 1) * dim / 2;
    fmf_createAlloc( &dfp, 1, nQP, sym, 1 );
    fmf_createAlloc( &gtdfp, 1, nQP, dim * nEPU, 1 );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( vg->det, iel );
      
    if (isDiff == 1) { 
      form_sdcc_actOpGT_M3( gtd, vg->bfGM, mtxD );
      fmf_mulAB_nn( gtdf, gtd, bf );
      fmf_sumLevelsMulF( out, gtdf, vg->det->val );
    } else {
      FMF_SetCell( pressure_qp, iel );
      fmf_mulAB_nn( dfp, mtxD, pressure_qp );
      form_sdcc_actOpGT_VS3( gtdfp, vg->bfGM, dfp );
      fmf_sumLevelsMulF( out, gtdfp, vg->det->val );
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
		   FMField *bf, FMField *mtxD, VolumeGeometry *vg,
		   int32 *elList, int32 elList_nRow,
		   int32 isDiff )
{
  int32 ii, iel, nEPP, nEP, dim, sym, nQP, ret = RET_OK;
  FMField *dtg = 0, *ftdtg = 0, *dtgu = 0, *ftdtgu = 0;
  FMField drow[1];

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;
  dim = vg->bfGM->nRow;
  sym = (dim + 1) * dim / 2;
  nEPP = bf->nCol;

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

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( vg->det, iel );
      
    if (isDiff == 1) {
      form_sdcc_actOpG_RM3( dtg, drow, vg->bfGM );
      fmf_mulATB_nn( ftdtg, bf, dtg );
      fmf_sumLevelsMulF( out, ftdtg, vg->det->val );
    } else {
      FMF_SetCell( strain, iel );
      fmf_mulATB_nn( dtgu, mtxD, strain );
      fmf_mulATB_nn( ftdtgu, bf, dtgu );
      fmf_sumLevelsMulF( out, ftdtgu, vg->det->val );
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
		  FMField *mtxD, VolumeGeometry *vg,
		  int32 *elList, int32 elList_nRow )
{
  int32 ii, iel, nQP, ret = RET_OK;
  FMField *dtgu = 0, *pftdtgu = 0;

  nQP = vg->bfGM->nLev;

/*   fmf_print( mtxD, stdout, 0 ); */

  fmf_createAlloc( &dtgu, 1, nQP, 1, 1 );
  fmf_createAlloc( &pftdtgu, 1, nQP, 1, 1 );

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg->det, iel );
    FMF_SetCell( state, iel );
    FMF_SetCell( strain, iel );
      
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
