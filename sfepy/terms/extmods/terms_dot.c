#include "terms.h"
#include "geommech.h"

#undef __FUNC__
#define __FUNC__ "dw_volume_dot_vector"
/*!
  @par Revision history:
  - 21.11.2006, c
*/
int32 dw_volume_dot_vector( FMField *out, FMField *coef, FMField *val_qp,
                            FMField *rbf, FMField *cbf, VolumeGeometry *vg,
                            int32 isDiff )
{
  int32 ii, dim, nQP, nEPR, nEPC, ret = RET_OK;
  FMField *ftfu = 0, *ftf1 = 0, *ftf = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;
  nEPR = rbf->nCol;
  nEPC = cbf->nCol;

  if (isDiff) {
    fmf_createAlloc( &ftf, 1, nQP, nEPR * dim, nEPC * dim );
    fmf_createAlloc( &ftf1, 1, nQP, nEPR, nEPC );

    fmf_mulATB_nn( ftf1, rbf, cbf );
  } else {
    fmf_createAlloc( &ftfu, 1, nQP, dim * nEPR, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    if (coef->nCell > 1) {
      FMF_SetCell( coef, ii );
    }
    FMF_SetCell( vg->det, ii );

    if (isDiff) {
      bf_buildFTF( ftf, ftf1 );
      fmf_mul( ftf, coef->val );
      fmf_sumLevelsMulF( out, ftf, vg->det->val );
    } else {
      FMF_SetCell( val_qp, ii );

      bf_actt( ftfu, rbf, val_qp );
      fmf_mul( ftfu, coef->val );
      fmf_sumLevelsMulF( out, ftfu, vg->det->val );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftf1 );
    fmf_freeDestroy( &ftf );
  } else {
    fmf_freeDestroy( &ftfu );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_surface_dot_vector"
int32 dw_surface_dot_vector( FMField *out, FMField *coef, FMField *val_qp,
                             FMField *rbf, FMField *cbf, SurfaceGeometry *sg,
                             int32 isDiff )
{
  int32 ii, dim, nQP, nEPR, nEPC, ret = RET_OK;
  FMField *ftfu = 0, *ftf1 = 0, *ftf = 0;

  nQP = sg->normal->nLev;
  dim = sg->normal->nRow;
  nEPR = rbf->nCol;
  nEPC = cbf->nCol;

  if (isDiff) {
    fmf_createAlloc( &ftf, 1, nQP, nEPR * dim, nEPC * dim );
    fmf_createAlloc( &ftf1, 1, nQP, nEPR, nEPC );

    fmf_mulATB_nn( ftf1, rbf, cbf );
  } else {
    fmf_createAlloc( &ftfu, 1, nQP, dim * nEPR, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    if (coef->nCell > 1) {
      FMF_SetCell( coef, ii );
    }
    FMF_SetCell( sg->det, ii );

    if (isDiff) {
      bf_buildFTF( ftf, ftf1 );
      fmf_mul( ftf, coef->val );
      fmf_sumLevelsMulF( out, ftf, sg->det->val );
    } else {
      FMF_SetCell( val_qp, ii );

      bf_actt( ftfu, rbf, val_qp );
      fmf_mul( ftfu, coef->val );
      fmf_sumLevelsMulF( out, ftfu, sg->det->val );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftf1 );
    fmf_freeDestroy( &ftf );
  } else {
    fmf_freeDestroy( &ftfu );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_volume_dot_scalar"
/*!
  @par Revision history:
  - 01.02.2008, c
*/
int32 dw_volume_dot_scalar( FMField *out, FMField *coef, FMField *val_qp,
                            FMField *rbf, FMField *cbf, VolumeGeometry *vg,
                            int32 isDiff )
{
  int32 ii, nQP, nEPR, nEPC, ret = RET_OK;
  FMField *ftfp = 0, *ftf = 0, *cftf = 0;

  nQP = vg->bfGM->nLev;
  nEPR = rbf->nCol;
  nEPC = cbf->nCol;

  if (isDiff) {
    fmf_createAlloc( &ftf, 1, nQP, nEPR, nEPC );
    fmf_createAlloc( &cftf, 1, nQP, nEPR, nEPC );

    fmf_mulATB_nn( ftf, rbf, cbf );
  } else {
    fmf_createAlloc( &ftfp, 1, nQP, nEPR, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( vg->det, ii );
    if (coef->nCell > 1) {
      FMF_SetCell( coef, ii );
    }

    if (isDiff) {
      fmf_mulAF( cftf, ftf, coef->val );
      fmf_sumLevelsMulF( out, cftf, vg->det->val );
    } else {
      FMF_SetCell( val_qp, ii );

      bf_actt( ftfp, rbf, val_qp );
      fmf_mul( ftfp, coef->val );
      fmf_sumLevelsMulF( out, ftfp, vg->det->val );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftf );
    fmf_freeDestroy( &cftf );
  } else {
    fmf_freeDestroy( &ftfp );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_surface_dot_scalar"
/*!
  @par Revision history:
  - 09.03.2009, c
*/
int32 dw_surface_dot_scalar( FMField *out, FMField *coef, FMField *val_qp,
                             FMField *rbf, FMField *cbf, SurfaceGeometry *sg,
                             int32 isDiff )
{
  int32 ii, nQP, nEPR, nEPC, ret = RET_OK;
  FMField *ftfp = 0, *ftf = 0, *cftf = 0;

  nQP = sg->normal->nLev;
  nEPR = rbf->nCol;
  nEPC = cbf->nCol;

  if (isDiff) {
    fmf_createAlloc( &ftf, 1, nQP, nEPR, nEPC );
    fmf_createAlloc( &cftf, 1, nQP, nEPR, nEPC );

    fmf_mulATB_nn( ftf, rbf, cbf );
  } else {
    fmf_createAlloc( &ftfp, 1, nQP, nEPR, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( sg->det, ii );
    if (coef->nCell > 1) {
      FMF_SetCell( coef, ii );
    }

    if (isDiff) {
      fmf_mulAF( cftf, ftf, coef->val );
      fmf_sumLevelsMulF( out, cftf, sg->det->val );
    } else {
      FMF_SetCell( val_qp, ii );

      bf_actt( ftfp, rbf, val_qp );
      fmf_mul( ftfp, coef->val );
      fmf_sumLevelsMulF( out, ftfp, sg->det->val );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftf );
    fmf_freeDestroy( &cftf );
  } else {
    fmf_freeDestroy( &ftfp );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_v_dot_grad_s_vw"
int32 dw_v_dot_grad_s_vw( FMField *out, FMField *coef, FMField *grad,
                          FMField *vbf, VolumeGeometry *vvg,
                          VolumeGeometry *svg, int32 isDiff )
{
  int32 ii, nEPV, nEPS, dim, nQP, ret = RET_OK;
  FMField *ftg = 0;

  nQP = vvg->bfGM->nLev;
  dim = vvg->bfGM->nRow;
  nEPS = svg->bfGM->nCol;
  nEPV = vbf->nCol;

  if (isDiff == 1) {
    fmf_createAlloc( &ftg, 1, nQP, dim * nEPV, nEPS );
  } else {
    fmf_createAlloc( &ftg, 1, nQP, dim * nEPV, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    if (coef->nCell > 1) {
      FMF_SetCell( coef, ii );
    }
    FMF_SetCell( vvg->det, ii );

    if (isDiff == 1) {
      FMF_SetCell( svg->bfGM, ii );

      bf_actt( ftg, vbf, svg->bfGM );
    } else {
      FMF_SetCell( grad, ii );

      bf_actt_c1( ftg, vbf, grad );
    }
    fmf_mul( ftg, coef->val );
    fmf_sumLevelsMulF( out, ftg, vvg->det->val );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &ftg );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_v_dot_grad_s_sw"
int32 dw_v_dot_grad_s_sw( FMField *out, FMField *coef, FMField *val_qp,
                          FMField *vbf, VolumeGeometry *vvg,
                          VolumeGeometry *svg, int32 isDiff )
{
  int32 ii, nEPV, nEPS, dim, nQP, ret = RET_OK;
  FMField *ftg = 0, *gtf = 0;

  nQP = vvg->bfGM->nLev;
  dim = vvg->bfGM->nRow;
  nEPS = svg->bfGM->nCol;
  nEPV = vbf->nCol;

  if (isDiff == 1) {
    fmf_createAlloc( &ftg, 1, nQP, dim * nEPV, nEPS );
  } else {
    fmf_createAlloc( &gtf, 1, nQP, nEPS, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    if (coef->nCell > 1) {
      FMF_SetCell( coef, ii );
    }
    FMF_SetCell( svg->bfGM, ii );
    FMF_SetCell( vvg->det, ii );

    if (isDiff == 1) {
      bf_actt( ftg, vbf, svg->bfGM );
      fmf_mul( ftg, coef->val );

      fmf_sumLevelsTMulF( out, ftg, vvg->det->val );
    } else {
      FMF_SetCell( val_qp, ii );

      fmf_mulATB_nn( gtf, svg->bfGM, val_qp );
      fmf_mul( gtf, coef->val );
      fmf_sumLevelsMulF( out, gtf, vvg->det->val );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftg );
  } else {
    fmf_freeDestroy( &gtf );
  }

  return( ret );
}
