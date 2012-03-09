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
