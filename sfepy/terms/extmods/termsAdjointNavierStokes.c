#include "termsNavierStokes.h"
#include "termsAdjointNavierStokes.h"
#include "terms.h"
#include "geommech.h"


#undef __FUNC__
#define __FUNC__ "sub_mul_gradddgrad_scalar"
/*!
  @par Revision history:
  - 24.10.2007, c
*/
int32 sub_mul_gradddgrad_scalar( FMField *out,
				 FMField *grad1, FMField *grad2,
				 FMField *scalar )
{
  int32 ir, ic, iqp;
  int32 nQP = scalar->nLev;
  int32 d2 = grad1->nRow;
  int32 dim = sqrt( d2 );

  for (iqp = 0; iqp < nQP; iqp++) {
    for (ir = 0; ir < dim; ir++) {
      for (ic = 0; ic < dim; ic++) {
	out->val[iqp] -= scalar->val[iqp] * grad1->val[d2*iqp+dim*ir+ic]
	  * grad2->val[d2*iqp+dim*ic+ir];
      }
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "dw_adj_convect1"
/*!
  @par Revision history:
  - 12.12.2005, c
  - 14.12.2005
  - 22.03.2006
  - 27.07.2006
*/
int32 dw_adj_convect1( FMField *out, FMField *stateW, FMField *gradU,
                       FMField *bf, VolumeGeometry *vg, int32 isDiff )
{
  int32 ii, dim, nQP, nEP, ret = RET_OK;
  FMField *gf = 0, *ftgf = 0, *ftgfu = 0;
  FMField *gfu = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;
  dim = vg->bfGM->nRow;

  if (isDiff) {
    fmf_createAlloc( &gf, 1, nQP, dim, dim * nEP );
    fmf_createAlloc( &ftgf, 1, nQP, dim * nEP, dim * nEP );
  } else {
    fmf_createAlloc( &gfu, 1, nQP, dim, 1 );
    fmf_createAlloc( &ftgfu, 1, nQP, dim * nEP, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( gradU, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );

    if (isDiff) {
      bf_ract( gf, bf, gradU );
      bf_actt( ftgf, bf, gf );
      fmf_sumLevelsMulF( out, ftgf, vg->det->val );
    } else {
      FMF_SetCell( stateW, ii );
      fmf_mulAB_nn( gfu, gradU, stateW );
      bf_actt( ftgfu, bf, gfu );
      fmf_sumLevelsMulF( out, ftgfu, vg->det->val );
    }
    ERR_CheckGo( ret );
  }


 end_label:
  if (isDiff) {
    fmf_freeDestroy( &gf );
    fmf_freeDestroy( &ftgf );
  } else {
    fmf_freeDestroy( &gfu );
    fmf_freeDestroy( &ftgfu );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_adj_convect2"
/*!
  @par Revision history:
  - 20.12.2005, c
  - 22.03.2006
  - 27.07.2006
*/
int32 dw_adj_convect2( FMField *out, FMField *stateW, FMField *stateU,
                       FMField *bf, VolumeGeometry *vg, int32 isDiff )
{
  int32 ii, dim, nQP, nEP, ret = RET_OK;
  FMField *vtg = 0, *ftvtg = 0;
  FMField *futvtg = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;
  dim = vg->bfGM->nRow;

  fmf_createAlloc( &vtg, 1, nQP, dim, nEP * dim );
  if (isDiff) {
    fmf_createAlloc( &ftvtg, 1, nQP, dim * nEP, dim * nEP );
  } else {
    fmf_createAlloc( &futvtg, 1, nQP, 1, dim * nEP );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( stateU, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );

    convect_build_vtg( vtg, vg->bfGM, stateU );

    if (isDiff) {
      bf_actt( ftvtg, bf, vtg );
      fmf_sumLevelsTMulF( out, ftvtg, vg->det->val );
    } else {
      FMF_SetCell( stateW, ii );
      fmf_mulATB_nn( futvtg, stateW, vtg );

      fmf_sumLevelsTMulF( out, futvtg, vg->det->val );
    }
    ERR_CheckGo( ret );
  }


 end_label:
  fmf_freeDestroy( &vtg );

  if (isDiff) {
    fmf_freeDestroy( &ftvtg );
  } else {
    fmf_freeDestroy( &futvtg );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_st_adj_supg_c"
/*!
  @par Revision history:
  - 30.10.2007, c
*/
int32 dw_st_adj_supg_c( FMField *out, FMField *stateW,
			FMField *stateU, FMField *gradU,
			FMField *coef, FMField *bf, VolumeGeometry *vg,
			int32 *conn, int32 nEl, int32 nEP,
			int32 isDiff )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *stW = 0, *gUfU = 0, *fUTg = 0;
  FMField *gUfUTg = 0, *fTgUfUTg = 0;
  FMField *gUfUTgT = 0, *fTgUfUTgT = 0;
  FMField *outdqp = 0, *outqp = 0, *out1qp = 0, *out2qp = 0;
  FMField stWv[1];

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  stateW->val = FMF_PtrFirst( stateW );

  fmf_createAlloc( &gUfU, 1, nQP, dim, 1 );
  fmf_createAlloc( &gUfUTgT, 1, nQP, dim, nEP * dim );
  fmf_createAlloc( &fTgUfUTgT, 1, nQP, nEP * dim, nEP * dim );

  fmf_createAlloc( &fUTg, 1, nQP, dim, nEP * dim );
  fmf_createAlloc( &gUfUTg, 1, nQP, dim, nEP * dim );
  fmf_createAlloc( &fTgUfUTg, 1, nQP, nEP * dim, nEP * dim );

  if (isDiff == 1) {
    fmf_createAlloc( &outdqp, 1, nQP, dim * nEP, dim * nEP );
  } else {
    fmf_createAlloc( &stW, 1, 1, dim, nEP );
    stWv->nAlloc = -1;
    fmf_pretend( stWv, 1, 1, nEP * dim, 1, stW->val );

    fmf_createAlloc( &out1qp, 1, nQP, dim * nEP, 1 );
    fmf_createAlloc( &out2qp, 1, nQP, dim * nEP, 1 );
    fmf_createAlloc( &outqp, 1, nQP, dim * nEP, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( stateU, ii );
    FMF_SetCell( gradU, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCell( coef, ii );

    // u grad u.
    fmf_mulAB_nn( gUfU, gradU, stateU );
    // (u grad u, grad).
    convect_build_vtbg( gUfUTgT, vg->bfGM, gUfU );
    // (u grad u, v grad).
    bf_actt( fTgUfUTgT, bf, gUfUTgT );

    // u grad.
    convect_build_vtg( fUTg, vg->bfGM, stateU );
    // (grad u, u^T grad).
    fmf_mulAB_nn( gUfUTg, gradU, fUTg );
    // (v grad u, u grad).
    bf_actt( fTgUfUTg, bf, gUfUTg );

    if (isDiff == 1) {
      fmf_addAB_nn( outdqp, fTgUfUTgT, fTgUfUTg );
      fmf_sumLevelsMulF( out, outdqp, vg->det->val );
    } else {
      ele_extractNodalValuesDBD( stW, stateW, conn + nEP * ii );

      // (u grad u, v grad w).
      fmf_mulAB_n1( out1qp, fTgUfUTgT, stWv );

      // (v grad u, u grad w).
      fmf_mulAB_n1( out2qp, fTgUfUTg, stWv );

      fmf_addAB_nn( outqp, out1qp, out2qp );
      fmf_sumLevelsMulF( out, outqp, vg->det->val );

    }
    fmf_mulC( out, coef->val[0] );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &gUfU );
  fmf_freeDestroy( &gUfUTgT );
  fmf_freeDestroy( &fTgUfUTgT );
  fmf_freeDestroy( &fUTg );
  fmf_freeDestroy( &gUfUTg );
  fmf_freeDestroy( &fTgUfUTg );
  if (isDiff) {
    fmf_freeDestroy( &outdqp );
  } else {
    fmf_freeDestroy( &stW );
    fmf_freeDestroy( &out1qp );
    fmf_freeDestroy( &out2qp );
    fmf_freeDestroy( &outqp );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_st_adj1_supg_p"
/*!
  @par Revision history:
  - 30.10.2007, c
*/
int32 dw_st_adj1_supg_p( FMField *out, FMField *stateW, FMField *gradP,
			 FMField *coef, FMField *bf_w, VolumeGeometry *vg_w,
			 int32 *conn_w, int32 nEl_w, int32 nEP_w,
			 int32 isDiff )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *stW = 0, *gPTgT = 0, *fTgPTgT = 0;
  FMField *outqp = 0;
  FMField stWv[1];

  nQP = vg_w->bfGM->nLev;
  dim = vg_w->bfGM->nRow;

  stateW->val = FMF_PtrFirst( stateW );

  fmf_createAlloc( &gPTgT, 1, nQP, dim, nEP_w * dim );
  fmf_createAlloc( &fTgPTgT, 1, nQP, nEP_w * dim, nEP_w * dim );

  if (isDiff == 0) {
    fmf_createAlloc( &outqp, 1, nQP, nEP_w * dim, 1 );

    fmf_createAlloc( &stW, 1, 1, dim, nEP_w );
    stWv->nAlloc = -1;
    fmf_pretend( stWv, 1, 1, nEP_w * dim, 1, stW->val );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( gradP, ii );
    FMF_SetCell( vg_w->bfGM, ii );
    FMF_SetCell( vg_w->det, ii );
    FMF_SetCell( coef, ii );

    // (grad p, grad).
    convect_build_vtbg( gPTgT, vg_w->bfGM, gradP );
    // (grad p, v grad).
    bf_actt( fTgPTgT, bf_w, gPTgT );

    if (isDiff == 1) {
      fmf_sumLevelsMulF( out, fTgPTgT, vg_w->det->val );
    } else {
      ele_extractNodalValuesDBD( stW, stateW, conn_w + nEP_w * ii );

      // (grad p, v grad w).
      fmf_mulAB_n1( outqp, fTgPTgT, stWv );
      fmf_sumLevelsMulF( out, outqp, vg_w->det->val );
    }
    fmf_mulC( out, coef->val[0] );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &gPTgT );
  fmf_freeDestroy( &fTgPTgT );
  if (isDiff == 0) {
    fmf_freeDestroy( &stW );
    fmf_freeDestroy( &outqp );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_st_adj2_supg_p"
/*!
  @par Revision history:
  - 30.10.2007, c
*/
int32 dw_st_adj2_supg_p( FMField *out, FMField *gradU, FMField *stateR,
			 FMField *coef, FMField *bf_u,
			 VolumeGeometry *vg_u, VolumeGeometry *vg_r,
			 int32 *conn_r, int32 nEl_r, int32 nEP_r,
			 int32 isDiff )
{
  int32 ii, dim, nQP, nEP_u, ret = RET_OK;
  FMField *stR = 0, *gUTg = 0, *fTgUTg = 0;
  FMField *outqp = 0;
  FMField stRv[1];

  nQP = vg_u->bfGM->nLev;
  nEP_u = vg_u->bfGM->nCol;
  dim = vg_u->bfGM->nRow;

  stateR->val = FMF_PtrFirst( stateR );

  fmf_createAlloc( &gUTg, 1, nQP, dim, nEP_r );
  fmf_createAlloc( &fTgUTg, 1, nQP, nEP_u * dim, nEP_r );

  if (isDiff == 0) {
    fmf_createAlloc( &outqp, 1, nQP, nEP_u * dim, 1 );

    fmf_createAlloc( &stR, 1, 1, 1, nEP_r );
    stRv->nAlloc = -1;
    fmf_pretend( stRv, 1, 1, nEP_r, 1, stR->val );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( gradU, ii );
    FMF_SetCell( vg_u->bfGM, ii );
    FMF_SetCell( vg_r->bfGM, ii );
    FMF_SetCell( vg_u->det, ii );
    FMF_SetCell( coef, ii );

    // (grad u, grad).
    fmf_mulATB_nn( gUTg, gradU, vg_r->bfGM );

    // (v grad u, grad).
    bf_actt( fTgUTg, bf_u, gUTg );

    if (isDiff == 1) {
      fmf_sumLevelsMulF( out, fTgUTg, vg_u->det->val );
    } else {
      ele_extractNodalValuesDBD( stR, stateR, conn_r + nEP_r * ii );

      // (v grad u, grad r).
      fmf_mulAB_n1( outqp, fTgUTg, stRv );
      fmf_sumLevelsMulF( out, outqp, vg_u->det->val );
    }
    fmf_mulC( out, coef->val[0] );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &gUTg );
  fmf_freeDestroy( &fTgUTg );
  if (isDiff == 0) {
    fmf_freeDestroy( &stR );
    fmf_freeDestroy( &outqp );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_of_nsMinGrad"
/*!
  @par Revision history:
  - 23.01.2006, c
  - 02.03.2006
  - 27.07.2006
*/
int32 d_of_nsMinGrad( FMField *out, FMField *grad,
		      FMField *viscosity, VolumeGeometry *vg )
{
  int32 ii, nQP, ret = RET_OK;
  float64 aux;
  FMField *out1 = 0, *gvel2 = 0;

  nQP = vg->bfGM->nLev;

  fmf_createAlloc( &out1, 1, 1, 1, 1 );
  fmf_createAlloc( &gvel2, 1, nQP, 1, 1 );

  FMF_SetFirst( out );
  aux = 0.0;
  for (ii = 0; ii < grad->nCell; ii++) {
    FMF_SetCell( grad, ii );
    FMF_SetCell( viscosity, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );

    fmf_mulATB_nn( gvel2, grad, grad );
    fmf_mul( gvel2, viscosity->val );
    fmf_sumLevelsMulF( out1, gvel2, vg->det->val );
    aux += out1->val[0];

    ERR_CheckGo( ret );
  }

  out->val[0] = aux * 0.5;

 end_label:
  fmf_freeDestroy( &out1 );
  fmf_freeDestroy( &gvel2 );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_of_nsSurfMinDPress"
/*!
  @par Revision history:
  - 23.03.2007, c
*/
int32 d_of_nsSurfMinDPress( FMField *out, FMField *pressure,
                            float64 weight, float64 bpress,
			    FMField *bf, SurfaceGeometry *sg, int32 isDiff )
{
  int32 ii, iqp, nQP, ret = RET_OK;
  float64 aux;
  FMField *out1 = 0, *pressQP = 0;

  nQP = sg->det->nLev;

  if (isDiff == 0) {
    fmf_createAlloc( &out1, 1, 1, 1, 1 );

    fmf_createAlloc( &pressQP, 1, nQP, 1, 1 );

    aux = 0.0;
    for (ii = 0; ii < pressure->nCell; ii++) {
      FMF_SetCell( pressure, ii );
      FMF_SetCell( sg->det, ii );

      for (iqp = 0; iqp < nQP; iqp++) {
	pressQP->val[iqp] -= pressure->val[iqp] - bpress;
      }
      fmf_sumLevelsMulF( out1, pressQP, sg->det->val );
      aux += out1->val[0];

      ERR_CheckGo( ret );
    }
    out->val[0] = aux * weight;
  } else {
    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCell( sg->det, ii );
      fmf_sumLevelsTMulF( out, bf, sg->det->val );

      ERR_CheckGo( ret );
    }
    fmfc_mulC( out, weight );
  }

 end_label:
  if (isDiff == 0) {
    fmf_freeDestroy( &out1 );
    fmf_freeDestroy( &pressQP );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_sd_div"
/*!
  @par Revision history:
  - 06.03.2006, c
  - 24.10.2007
*/
int32 d_sd_div( FMField *out,
		FMField *stateU, int32 offsetU,
		FMField *stateP, int32 offsetP,
		FMField *vecMV, int32 offsetMV,
		FMField *bf_p,
		VolumeGeometry *vg_u,
		VolumeGeometry *vg_p,
		VolumeGeometry *vg_mv,
		int32 *conn_u, int32 nEl_u, int32 nEP_u,
		int32 *conn_p, int32 nEl_p, int32 nEP_p,
		int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
		int32 *elList, int32 elList_nRow,
		int32 mode )
{
  int32 ii, iel, dim, nQP, d2, ret = RET_OK;
  FMField *stp = 0, *stu = 0, *fp = 0, *aux11 = 0, *divU = 0;
  FMField *mv = 0, *divMV = 0, *gu = 0, *gmv = 0;
  FMField stpv[1], stuv[1], mvv[1], gclu[1], gclmv[1];

  nQP = vg_u->bfGM->nLev;
  dim = vg_u->bfGM->nRow;
  d2 = dim * dim;

  stateU->val = FMF_PtrFirst( stateU ) + offsetU;
  stateP->val = FMF_PtrFirst( stateP ) + offsetP;

  fmf_createAlloc( &stp, 1, 1, 1, nEP_p );
  stpv->nAlloc = -1;
  fmf_pretend( stpv, 1, 1, nEP_p, 1, stp->val );

  fmf_createAlloc( &stu, 1, 1, dim, nEP_u );
  stuv->nAlloc = -1;
  fmf_pretend( stuv, 1, 1, nEP_u * dim, 1, stu->val );

  gclu->nAlloc = -1;
  fmf_pretend( gclu, 1, nQP, 1, nEP_u * dim, vg_u->bfGM->val0 );

  fmf_createAlloc( &fp, 1, nQP, 1, 1 );
  fmf_createAlloc( &divU, 1, nQP, 1, 1 );
  fmf_createAlloc( &aux11, 1, nQP, 1, 1 );

  if (mode == 1) { 
    vecMV->val = FMF_PtrFirst( vecMV ) + offsetMV;

    gclmv->nAlloc = -1;
    fmf_pretend( gclmv, 1, nQP, 1, nEP_mv * dim, vg_mv->bfGM->val0 );

    fmf_createAlloc( &mv, 1, 1, dim, nEP_mv );
    mvv->nAlloc = -1;
    fmf_pretend( mvv, 1, 1, nEP_mv * dim, 1, mv->val );

    fmf_createAlloc( &divMV, 1, nQP, 1, 1 );
    fmf_createAlloc( &gu, 1, nQP, d2, 1 );
    fmf_createAlloc( &gmv, 1, nQP, d2, 1 );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg_u->det, iel );
    FMF_SetCell( gclu, iel );

    ele_extractNodalValuesDBD( stp, stateP, conn_p + nEP_p * iel );
    ele_extractNodalValuesDBD( stu, stateU, conn_u + nEP_u * iel );

    fmf_mulAB_n1( fp, bf_p, stpv );
    fmf_mulAB_n1( divU, gclu, stuv );
    fmf_mulAB_nn( aux11, fp, divU );
      
    if (mode == 1) { 
      FMF_SetCell( gclmv, iel );
      FMF_SetCell( vg_u->bfGM, iel );
      FMF_SetCell( vg_mv->bfGM, iel );

      ele_extractNodalValuesDBD( mv, vecMV, conn_mv + nEP_mv * iel );
      fmf_mulAB_n1( divMV, gclmv, mvv );
      fmf_mul( aux11, divMV->val );

      divgrad_act_g_m( gu, vg_u->bfGM, stuv );
      divgrad_act_g_m( gmv, vg_mv->bfGM, mvv );
      sub_mul_gradddgrad_scalar( aux11, gu, gmv, fp );
    }
    fmf_sumLevelsMulF( out, aux11, vg_u->det->val );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &stp );
  fmf_freeDestroy( &stu );
  fmf_freeDestroy( &fp );
  fmf_freeDestroy( &aux11 );
  fmf_freeDestroy( &divU );
  if (mode == 1) {
    fmf_freeDestroy( &mv );
    fmf_freeDestroy( &divMV );
    fmf_freeDestroy( &gu );
    fmf_freeDestroy( &gmv );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_sd_div_grad"
/*!
  @par Revision history:
  - 08.01.2006, c
  - 09.01.2006
  - 27.02.2006
  - 06.03.2006
  - 07.03.2006
*/
int32 d_sd_div_grad( FMField *out,
		     FMField *stateU, int32 offsetU,
		     FMField *stateW, int32 offsetW,
		     FMField *vecMV, int32 offsetMV,
		     float64 viscosity,
		     VolumeGeometry *vg_u,
		     VolumeGeometry *vg_mv,
		     int32 *conn_u, int32 nEl_u, int32 nEP_u,
		     int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
		     int32 *elList, int32 elList_nRow,
		     int32 mode )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *st = 0, *vel = 0, *mv = 0, *gw = 0, *gvel = 0, *gmv = 0;
  FMField *divMV = 0, *uvel = 0, *aux = 0, *aux1 = 0, *aux2 = 0, *aux3 = 0;
  FMField gcl[1], stv[1], velv[1], mvv[1], gwm[1], gvelm[1], gmvm[1], aux3m[1];

  nQP = vg_u->bfGM->nLev;
  dim = vg_u->bfGM->nRow;

  stateU->val = FMF_PtrFirst( stateU ) + offsetU;
  stateW->val = FMF_PtrFirst( stateW ) + offsetW;

  fmf_createAlloc( &uvel, 1, nQP, 1, 1 );

  fmf_createAlloc( &st, 1, 1, dim, nEP_u );
  stv->nAlloc = -1;
  fmf_pretend( stv, 1, 1, nEP_u * dim, 1, st->val );

  fmf_createAlloc( &vel, 1, 1, dim, nEP_u );
  velv->nAlloc = -1;
  fmf_pretend( velv, 1, 1, nEP_u * dim, 1, vel->val );

  fmf_createAlloc( &gw, 1, nQP, dim * dim, 1 );
  gwm->nAlloc = -1;
  fmf_pretend( gwm, 1, nQP, dim, dim, gw->val );

  fmf_createAlloc( &gvel, 1, nQP, dim * dim, 1 );
  gvelm->nAlloc = -1;
  fmf_pretend( gvelm, 1, nQP, dim, dim, gvel->val );

  if (mode == 1) {
    vecMV->val = FMF_PtrFirst( vecMV ) + offsetMV;

    gcl->nAlloc = -1;
    fmf_pretend( gcl, 1, nQP, 1, nEP_mv * dim, vg_mv->bfGM->val0 );

    fmf_createAlloc( &divMV, 1, nQP, 1, 1 );

    fmf_createAlloc( &mv, 1, 1, dim, nEP_mv );
    mvv->nAlloc = -1;
    fmf_pretend( mvv, 1, 1, nEP_mv * dim, 1, mv->val );

    fmf_createAlloc( &gmv, 1, nQP, dim * dim, 1 );
    gmvm->nAlloc = -1;
    fmf_pretend( gmvm, 1, nQP, dim, dim, gmv->val );

    fmf_createAlloc( &aux, 1, 1, 1, 1 );
    fmf_createAlloc( &aux1, 1, nQP, 1, 1 );
    fmf_createAlloc( &aux2, 1, nQP, 1, 1 );
    fmf_createAlloc( &aux3, 1, nQP, dim * dim, 1 );
    aux3m->nAlloc = -1;
    fmf_pretend( aux3m, 1, nQP, dim, dim, aux3->val );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg_u->bfGM, iel );
    FMF_SetCell( vg_u->det, iel );

    // w.      
    ele_extractNodalValuesDBD( st, stateW, conn_u + nEP_u * iel );
    // u.
    ele_extractNodalValuesDBD( vel, stateU, conn_u + nEP_u * iel );

    divgrad_act_g_m( gw, vg_u->bfGM, stv );
    divgrad_act_g_m( gvel, vg_u->bfGM, velv );

    if (mode == 0) {
      // (gu:gw)
      fmf_mulATB_nn( uvel, gw, gvel );
      fmf_sumLevelsMulF( out, uvel, vg_u->det->val );
      fmf_mulC( out, viscosity );

    } else if (mode == 1) {
      FMF_SetCell( vg_mv->bfGM, iel );
      FMF_SetCell( gcl, iel );
      // nu.
      ele_extractNodalValuesDBD( mv, vecMV, conn_mv + nEP_mv * iel );
      
      divgrad_act_g_m( gmv, vg_mv->bfGM, mvv );

      // (gu:gw) div nu
      fmf_mulAB_n1( divMV, gcl, mvv );
      fmf_mulATB_nn( uvel, gw, gvel );
      fmf_mulAB_nn( aux1, uvel, divMV );
      fmf_sumLevelsMulF( out, aux1, vg_u->det->val );

      // (gu * gnu) : gw
      fmf_mulAB_nn( aux3m, gvelm, gmvm );
      fmf_mulATB_nn( aux1, aux3, gw );

      // (gw * gnu) : gu
      fmf_mulAB_nn( aux3m, gwm, gmvm );
      fmf_mulATB_nn( aux2, aux3, gvel );

      fmf_addAB_nn( aux1, aux1, aux2 );
      fmf_sumLevelsMulF( aux, aux1, vg_u->det->val );
      fmf_subAB_nn( out, out, aux );

      fmf_mulC( out, viscosity );
    }
/*     fmfc_save( out, "000", 0 ); */
/*     sys_pause(); */

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &uvel );
  fmf_freeDestroy( &st );
  fmf_freeDestroy( &vel );
  fmf_freeDestroy( &gw );
  fmf_freeDestroy( &gvel );
  if (mode == 1) {
    fmf_freeDestroy( &divMV );
    fmf_freeDestroy( &mv );
    fmf_freeDestroy( &gmv );

    fmf_freeDestroy( &aux );
    fmf_freeDestroy( &aux1 );
    fmf_freeDestroy( &aux2 );
    fmf_freeDestroy( &aux3 );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_sd_convect"
/*!
  @par Revision history:
  - 22.03.2006, c
*/
int32 d_sd_convect( FMField *out,
		    FMField *stateU, int32 offsetU,
		    FMField *stateW, int32 offsetW,
		    FMField *vecMV, int32 offsetMV,
		    FMField *bf_u, FMField *bf_w,
		    VolumeGeometry *vg_u,
		    VolumeGeometry *vg_w,
		    VolumeGeometry *vg_mv,
		    int32 *conn_u, int32 nEl_u, int32 nEP_u,
		    int32 *conn_w, int32 nEl_w, int32 nEP_w,
		    int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
		    int32 *elList, int32 elList_nRow,
		    int32 mode )
{
  int32 ii, iel, dim, nQP, d2, ret = RET_OK;
  FMField *stw = 0, *stu = 0, *fw = 0, *fu = 0, *aux11 = 0, *aux = 0;
  FMField *mv = 0, *divMV = 0, *gu = 0, *gmv = 0;
  FMField *fwgu = 0, *fwgugmv = 0;
  FMField gum[1], stuv[1], mvv[1], gmvm[1], gclmv[1];

  nQP = vg_u->bfGM->nLev;
  dim = vg_u->bfGM->nRow;
  d2 = dim * dim;

  stateU->val = FMF_PtrFirst( stateU ) + offsetU;
  stateW->val = FMF_PtrFirst( stateW ) + offsetW;

  fmf_createAlloc( &stw, 1, 1, dim, nEP_w );

  fmf_createAlloc( &stu, 1, 1, dim, nEP_u );
  stuv->nAlloc = -1;
  fmf_pretend( stuv, 1, 1, nEP_u * dim, 1, stu->val );

  fmf_createAlloc( &fu, 1, nQP, dim, 1 );
  fmf_createAlloc( &fw, 1, nQP, dim, 1 );

  fmf_createAlloc( &gu, 1, nQP, d2, 1 );
  gum->nAlloc = -1;
  fmf_pretend( gum, 1, nQP, dim, dim, gu->val );

  fmf_createAlloc( &fwgu, 1, nQP, 1, dim );

  fmf_createAlloc( &aux11, 1, nQP, 1, 1 );

  if (mode == 1) { 
    vecMV->val = FMF_PtrFirst( vecMV ) + offsetMV;

    gclmv->nAlloc = -1;
    fmf_pretend( gclmv, 1, nQP, 1, nEP_mv * dim, vg_mv->bfGM->val0 );

    fmf_createAlloc( &mv, 1, 1, dim, nEP_mv );
    mvv->nAlloc = -1;
    fmf_pretend( mvv, 1, 1, nEP_mv * dim, 1, mv->val );

    fmf_createAlloc( &divMV, 1, nQP, 1, 1 );

    fmf_createAlloc( &gmv, 1, nQP, d2, 1 );
    gmvm->nAlloc = -1;
    fmf_pretend( gmvm, 1, nQP, dim, dim, gmv->val );

    fmf_createAlloc( &fwgugmv, 1, nQP, 1, dim );

    fmf_createAlloc( &aux, 1, nQP, 1, 1 );

  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg_u->det, iel );
    FMF_SetCell( vg_u->bfGM, iel );

    ele_extractNodalValuesDBD( stw, stateW, conn_w + nEP_w * iel );
    ele_extractNodalValuesDBD( stu, stateU, conn_u + nEP_u * iel );

    divgrad_act_g_m( gu, vg_u->bfGM, stuv );
    bf_act( fu, bf_u, stu );
    bf_act( fw, bf_w, stw );

    fmf_mulATB_nn( fwgu, fw, gum );
    fmf_mulAB_nn( aux11, fwgu, fu );
    
    if (mode == 1) { 
      FMF_SetCell( gclmv, iel );
      FMF_SetCell( vg_mv->bfGM, iel );

      ele_extractNodalValuesDBD( mv, vecMV, conn_mv + nEP_mv * iel );
      fmf_mulAB_n1( divMV, gclmv, mvv );
      fmf_mul( aux11, divMV->val );

      divgrad_act_g_m( gmv, vg_mv->bfGM, mvv );
      fmf_mulAB_nn( fwgugmv, fwgu, gmvm );
      fmf_mulAB_nn( aux, fwgugmv, fu );
      fmf_subAB_nn( aux11, aux11, aux );
    }
    fmf_sumLevelsMulF( out, aux11, vg_u->det->val );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &stw );
  fmf_freeDestroy( &stu );
  fmf_freeDestroy( &fu );
  fmf_freeDestroy( &fw );
  fmf_freeDestroy( &aux11 );
  fmf_freeDestroy( &gu );
  fmf_freeDestroy( &fwgu );
  if (mode == 1) {
    fmf_freeDestroy( &mv );
    fmf_freeDestroy( &divMV );
    fmf_freeDestroy( &gmv );
    fmf_freeDestroy( &fwgugmv );
    fmf_freeDestroy( &aux );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_sd_testPQ"
/*!
  mode == 0: \int pq
  mode == 1: \int pq div V

  vg, conn .. pressure-like, assuming the same nodes for p, q, mv

  @par Revision history:
  - 24.02.2006, c
  - 27.02.2006
*/
int32 d_sd_testPQ( FMField *out,
		   FMField *stateP, int32 offsetP,
		   FMField *stateQ, int32 offsetQ,
		   FMField *vecMV, int32 offsetMV,
		   FMField *bf, VolumeGeometry *vg,
		   int32 *conn, int32 nEl, int32 nEP,
		   int32 *elList, int32 elList_nRow,
		   int32 mode )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *stp = 0, *stq = 0, *fp = 0, *fq = 0, *pq = 0, *mv = 0, *divMV = 0;
  FMField stpv[1], stqv[1], mvv[1], gcl[1];

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  stateP->val = FMF_PtrFirst( stateP ) + offsetP;
  stateQ->val = FMF_PtrFirst( stateQ ) + offsetQ;

  fmf_createAlloc( &stp, 1, 1, 1, nEP );
  stpv->nAlloc = -1;
  fmf_pretend( stpv, 1, 1, nEP, 1, stp->val );

  fmf_createAlloc( &stq, 1, 1, 1, nEP );
  stqv->nAlloc = -1;
  fmf_pretend( stqv, 1, 1, nEP, 1, stq->val );

  fmf_createAlloc( &fp, 1, nQP, 1, 1 );
  fmf_createAlloc( &fq, 1, nQP, 1, 1 );
  fmf_createAlloc( &pq, 1, nQP, 1, 1 );

  if (mode == 1) { 
    vecMV->val = FMF_PtrFirst( vecMV ) + offsetMV;

    gcl->nAlloc = -1;
    fmf_pretend( gcl, 1, nQP, 1, nEP * dim, vg->bfGM->val0 );

    fmf_createAlloc( &mv, 1, 1, dim, nEP );
    mvv->nAlloc = -1;
    fmf_pretend( mvv, 1, 1, nEP * dim, 1, mv->val );

    fmf_createAlloc( &divMV, 1, nQP, 1, 1 );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg->det, iel );

    ele_extractNodalValuesDBD( stp, stateP, conn + nEP * iel );
    ele_extractNodalValuesDBD( stq, stateQ, conn + nEP * iel );

    fmf_mulAB_n1( fp, bf, stpv );
    fmf_mulAB_n1( fq, bf, stqv );
    fmf_mulAB_nn( pq, fp, fq );
      
    if (mode == 1) { 
      FMF_SetCell( gcl, iel );

      ele_extractNodalValuesDBD( mv, vecMV, conn + nEP * iel );
      fmf_mulAB_n1( divMV, gcl, mvv );
      fmf_mul( pq, divMV->val );

/*       fmf_print( stp, stdout, 0 ); */
/*       fmf_print( stq, stdout, 0 ); */
/*       fmf_print( mv, stdout, 0 ); */
/*       fmf_print( divMV, stdout, 0 ); */
/*       fmf_print( vecMV, stdout, 0 ); */
/*       sys_pause(); */
    }
    fmf_sumLevelsMulF( out, pq, vg->det->val );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &stp );
  fmf_freeDestroy( &stq );
  fmf_freeDestroy( &fp );
  fmf_freeDestroy( &fq );
  fmf_freeDestroy( &pq );
  if (mode == 1) {
    fmf_freeDestroy( &mv );
    fmf_freeDestroy( &divMV );
  }

  return( ret );
}


#undef __FUNC__
#define __FUNC__ "d_sd_st_grad_div"
/*!
  @par Revision history:
  - 24.10.2007, c
*/
int32 d_sd_st_grad_div( FMField *out,
			FMField *stateU, int32 offsetU,
			FMField *stateW, int32 offsetW,
			FMField *vecMV, int32 offsetMV,
			float64 coef,
			VolumeGeometry *vg_u,
			VolumeGeometry *vg_mv,
			int32 *conn_u, int32 nEl_u, int32 nEP_u,
			int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
			int32 *elList, int32 elList_nRow,
			int32 mode )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *stU = 0, *stW = 0, *stMV = 0, *gU = 0, *gW = 0, *gMV = 0;
  FMField *divU = 0, *divW = 0, *divMV = 0;
  FMField *scalar1 = 0, *scalar2 = 0;
  FMField gcl_u[1], gcl_mv[1], stUv[1], stWv[1], stMVv[1];

  nQP = vg_u->bfGM->nLev;
  dim = vg_u->bfGM->nRow;

  stateU->val = FMF_PtrFirst( stateU ) + offsetU;
  stateW->val = FMF_PtrFirst( stateW ) + offsetW;

  fmf_createAlloc( &stU, 1, 1, dim, nEP_u );
  stUv->nAlloc = -1;
  fmf_pretend( stUv, 1, 1, nEP_u * dim, 1, stU->val );

  fmf_createAlloc( &stW, 1, 1, dim, nEP_u );
  stWv->nAlloc = -1;
  fmf_pretend( stWv, 1, 1, nEP_u * dim, 1, stW->val );

  gcl_u->nAlloc = -1;
  fmf_pretend( gcl_u, 1, nQP, 1, nEP_u * dim, vg_u->bfGM->val0 );

  fmf_createAlloc( &divU, 1, nQP, 1, 1 );
  fmf_createAlloc( &divW, 1, nQP, 1, 1 );
  fmf_createAlloc( &scalar1, 1, nQP, 1, 1 );

  if (mode == 1) {
    vecMV->val = FMF_PtrFirst( vecMV ) + offsetMV;

    fmf_createAlloc( &stMV, 1, 1, dim, nEP_mv );
    stMVv->nAlloc = -1;
    fmf_pretend( stMVv, 1, 1, nEP_mv * dim, 1, stMV->val );

    gcl_mv->nAlloc = -1;
    fmf_pretend( gcl_mv, 1, nQP, 1, nEP_mv * dim, vg_mv->bfGM->val0 );

    fmf_createAlloc( &divMV, 1, nQP, 1, 1 );
    fmf_createAlloc( &scalar2, 1, nQP, 1, 1 );

    fmf_createAlloc( &gU, 1, nQP, dim * dim, 1 );
    fmf_createAlloc( &gW, 1, nQP, dim * dim, 1 );
    fmf_createAlloc( &gMV, 1, nQP, dim * dim, 1 );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

/*     output( "\n%d %d %d\n", mode, ii, iel ); */

    FMF_SetCell( out, ii );
    FMF_SetCell( vg_u->bfGM, iel );
    FMF_SetCell( gcl_u, iel );
    FMF_SetCell( vg_u->det, iel );

    // u.
    ele_extractNodalValuesDBD( stU, stateU, conn_u + nEP_u * iel );
    // div u.
    fmf_mulAB_n1( divU, gcl_u, stUv );

    // w.      
    ele_extractNodalValuesDBD( stW, stateW, conn_u + nEP_u * iel );
    // div w.
    fmf_mulAB_n1( divW, gcl_u, stWv );

    if (mode == 0) {
      // (div u, div w)
      fmf_mulATB_nn( scalar1, divU, divW );
      fmf_sumLevelsMulF( out, scalar1, vg_u->det->val );

    } else if (mode == 1) {
      FMF_SetCell( vg_mv->bfGM, iel );
      FMF_SetCell( gcl_mv, iel );

      // nu.
      ele_extractNodalValuesDBD( stMV, vecMV, conn_mv + nEP_mv * iel );
      // div nu.
      fmf_mulAB_n1( divMV, gcl_mv, stMVv );
      
      // grad u.
      divgrad_act_g_m( gU, vg_u->bfGM, stUv );
      // grad w.
      divgrad_act_g_m( gW, vg_u->bfGM, stWv );
      // grad mv
      divgrad_act_g_m( gMV, vg_mv->bfGM, stMVv );

      // div u div w div nu ...
      fmf_mulATB_nn( scalar1, divU, divW );
      fmf_mulATB_nn( scalar2, divMV, scalar1 );

      // ... - (gu : gnu^T) div w ...
      sub_mul_gradddgrad_scalar( scalar2, gMV, gU, divW );

      // ... - (gw : gnu^T) div u
      sub_mul_gradddgrad_scalar( scalar2, gMV, gW, divU );

      fmf_sumLevelsMulF( out, scalar2, vg_u->det->val );
    }
/*     fmfc_save( out, "000", 0 ); */
/*     sys_pause(); */

    ERR_CheckGo( ret );
  }

  fmfc_mulC( out, coef );

 end_label:
  fmf_freeDestroy( &stU );
  fmf_freeDestroy( &stW );
  fmf_freeDestroy( &divU );
  fmf_freeDestroy( &divW );
  fmf_freeDestroy( &scalar1 );
  if (mode == 1) { 
    fmf_freeDestroy( &stMV );
    fmf_freeDestroy( &divMV );
    fmf_freeDestroy( &gU );
    fmf_freeDestroy( &gW );
    fmf_freeDestroy( &gMV );
    fmf_freeDestroy( &scalar2 );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_sd_st_supg_c"
/*!
  @par Revision history:
  - 25.10.2007, c
  - 29.10.2007
*/
int32 d_sd_st_supg_c( FMField *out,
		      FMField *stateU, int32 offsetU,
		      FMField *stateB, int32 offsetB,
		      FMField *stateW, int32 offsetW,
		      FMField *vecMV, int32 offsetMV,
		      FMField *bf_u,
		      FMField *coef,
		      VolumeGeometry *vg_u,
		      VolumeGeometry *vg_mv,
		      int32 *conn_u, int32 nEl_u, int32 nEP_u,
		      int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
		      int32 *elList, int32 elList_nRow,
		      int32 mode )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *stU = 0, *stB = 0, *stW = 0, *stMV = 0, *gU = 0, *gW = 0, *gMV = 0;
  FMField *divMV = 0, *fB = 0;
  FMField *fBgU = 0, *fBgW = 0, *fBgMV = 0, *fBgMVgU = 0, *fBgMVgW = 0;
  FMField *scalar1 = 0, *scalar2 = 0;
  FMField gcl_mv[1], stUv[1], stBv[1], stWv[1], stMVv[1];
  FMField gUm[1], gWm[1], gMVm[1];

  nQP = vg_u->bfGM->nLev;
  dim = vg_u->bfGM->nRow;

  stateU->val = FMF_PtrFirst( stateU ) + offsetU;
  stateB->val = FMF_PtrFirst( stateB ) + offsetB;
  stateW->val = FMF_PtrFirst( stateW ) + offsetW;

  fmf_createAlloc( &stU, 1, 1, dim, nEP_u );
  stUv->nAlloc = -1;
  fmf_pretend( stUv, 1, 1, nEP_u * dim, 1, stU->val );

  fmf_createAlloc( &stB, 1, 1, dim, nEP_u );
  stBv->nAlloc = -1;
  fmf_pretend( stBv, 1, 1, nEP_u * dim, 1, stB->val );

  fmf_createAlloc( &stW, 1, 1, dim, nEP_u );
  stWv->nAlloc = -1;
  fmf_pretend( stWv, 1, 1, nEP_u * dim, 1, stW->val );

  fmf_createAlloc( &gU, 1, nQP, dim * dim, 1 );
  gUm->nAlloc = -1;
  fmf_pretend( gUm, 1, nQP, dim, dim, gU->val );

  fmf_createAlloc( &gW, 1, nQP, dim * dim, 1 );
  gWm->nAlloc = -1;
  fmf_pretend( gWm, 1, nQP, dim, dim, gW->val );

  fmf_createAlloc( &fB, 1, nQP, dim, 1 );
  fmf_createAlloc( &fBgU, 1, nQP, 1, dim );
  fmf_createAlloc( &fBgW, 1, nQP, 1, dim );

  fmf_createAlloc( &scalar1, 1, nQP, 1, 1 );

  if (mode == 1) {
    vecMV->val = FMF_PtrFirst( vecMV ) + offsetMV;

    fmf_createAlloc( &stMV, 1, 1, dim, nEP_mv );
    stMVv->nAlloc = -1;
    fmf_pretend( stMVv, 1, 1, nEP_mv * dim, 1, stMV->val );

    gcl_mv->nAlloc = -1;
    fmf_pretend( gcl_mv, 1, nQP, 1, nEP_mv * dim, vg_mv->bfGM->val0 );

    fmf_createAlloc( &divMV, 1, nQP, 1, 1 );
    fmf_createAlloc( &scalar2, 1, nQP, 1, 1 );

    fmf_createAlloc( &gMV, 1, nQP, dim * dim, 1 );
    gMVm->nAlloc = -1;
    fmf_pretend( gMVm, 1, nQP, dim, dim, gMV->val );

    fmf_createAlloc( &fBgMV, 1, nQP, 1, dim );
    fmf_createAlloc( &fBgMVgU, 1, nQP, 1, dim );
    fmf_createAlloc( &fBgMVgW, 1, nQP, 1, dim );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

/*     output( "\n%d %d %d\n", mode, ii, iel ); */

    FMF_SetCell( out, ii );
    FMF_SetCell( vg_u->bfGM, iel );
    FMF_SetCell( vg_u->det, iel );
    FMF_SetCell( coef, ii );

    // u.
    ele_extractNodalValuesDBD( stU, stateU, conn_u + nEP_u * iel );
    // grad u.
    divgrad_act_g_m( gU, vg_u->bfGM, stUv );

    // b.
    ele_extractNodalValuesDBD( stB, stateB, conn_u + nEP_u * iel );
    bf_act( fB, bf_u, stB );

    // w.
    ele_extractNodalValuesDBD( stW, stateW, conn_u + nEP_u * iel );
    // grad w.
    divgrad_act_g_m( gW, vg_u->bfGM, stWv );

    // b grad u.
    fmf_mulATBT_nn( fBgU, fB, gUm );

    // b grad w.
    fmf_mulATBT_nn( fBgW, fB, gWm );

    // (b grad u, b grad w).
    fmf_mulABT_nn( scalar1, fBgU, fBgW );

    if (mode == 0) {
      fmf_sumLevelsMulF( out, scalar1, vg_u->det->val );

    } else if (mode == 1) {
      FMF_SetCell( vg_mv->bfGM, iel );
      FMF_SetCell( gcl_mv, iel );

      // nu.
      ele_extractNodalValuesDBD( stMV, vecMV, conn_mv + nEP_mv * iel );
      // div nu.
      fmf_mulAB_n1( divMV, gcl_mv, stMVv );
      // grad nu.
      divgrad_act_g_m( gMV, vg_mv->bfGM, stMVv );
      
      // b grad nu.
      fmf_mulATBT_nn( fBgMV, fB, gMVm );

      // (b grad u, b grad w) div nu ...
      fmf_mulATB_nn( scalar2, divMV, scalar1 );

      // ... - (b grad nu grad u, b grad w) ...
      fmf_mulABT_nn( fBgMVgU, fBgMV, gUm );
      fmf_mulABT_nn( scalar1, fBgMVgU, fBgW );
      fmf_subAB_nn( scalar2, scalar2, scalar1 );

      // ... - (b grad nu grad w, b grad u) ...
      fmf_mulABT_nn( fBgMVgW, fBgMV, gWm );
      fmf_mulABT_nn( scalar1, fBgMVgW, fBgU );
      fmf_subAB_nn( scalar2, scalar2, scalar1 );

      fmf_sumLevelsMulF( out, scalar2, vg_u->det->val );
    }
/*     fmfc_save( out, "000", 0 ); */
/*     sys_pause(); */

    fmf_mulC( out, coef->val[0] );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &stU );
  fmf_freeDestroy( &stB );
  fmf_freeDestroy( &stW );
  fmf_freeDestroy( &gU );
  fmf_freeDestroy( &gW );
  fmf_freeDestroy( &scalar1 );
  fmf_freeDestroy( &fB );
  fmf_freeDestroy( &fBgU );
  fmf_freeDestroy( &fBgW );
  if (mode == 1) { 
    fmf_freeDestroy( &stMV );
    fmf_freeDestroy( &divMV );
    fmf_freeDestroy( &gMV );
    fmf_freeDestroy( &scalar2 );
    fmf_freeDestroy( &fBgMV );
    fmf_freeDestroy( &fBgMVgU );
    fmf_freeDestroy( &fBgMVgW );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_sd_st_pspg_c"
/*!
  @par Revision history:
  - 25.10.2007, c
  - 29.10.2007
*/
int32 d_sd_st_pspg_c( FMField *out,
		      FMField *stateU, int32 offsetU,
		      FMField *stateB, int32 offsetB,
		      FMField *stateR, int32 offsetR,
		      FMField *vecMV, int32 offsetMV,
		      FMField *bf_u,
		      FMField *coef,
		      VolumeGeometry *vg_u,
		      VolumeGeometry *vg_r,
		      VolumeGeometry *vg_mv,
		      int32 *conn_u, int32 nEl_u, int32 nEP_u,
		      int32 *conn_r, int32 nEl_r, int32 nEP_r,
		      int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
		      int32 *elList, int32 elList_nRow,
		      int32 mode )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *stU = 0, *stB = 0, *stR = 0, *stMV = 0, *gU = 0, *gR = 0, *gMV = 0;
  FMField *divMV = 0, *fB = 0;
  FMField *fBgU = 0, *fBgMV = 0, *fBgMVgU = 0, *gMVgR = 0;
  FMField *scalar1 = 0, *scalar2 = 0;
  FMField gcl_mv[1], stUv[1], stBv[1], stRv[1], stMVv[1];
  FMField gUm[1], gMVm[1];

  nQP = vg_u->bfGM->nLev;
  dim = vg_u->bfGM->nRow;

  stateU->val = FMF_PtrFirst( stateU ) + offsetU;
  stateB->val = FMF_PtrFirst( stateB ) + offsetB;
  stateR->val = FMF_PtrFirst( stateR ) + offsetR;

  fmf_createAlloc( &stU, 1, 1, dim, nEP_u );
  stUv->nAlloc = -1;
  fmf_pretend( stUv, 1, 1, nEP_u * dim, 1, stU->val );

  fmf_createAlloc( &stB, 1, 1, dim, nEP_u );
  stBv->nAlloc = -1;
  fmf_pretend( stBv, 1, 1, nEP_u * dim, 1, stB->val );

  fmf_createAlloc( &stR, 1, 1, 1, nEP_r );
  stRv->nAlloc = -1;
  fmf_pretend( stRv, 1, 1, nEP_r, 1, stR->val );

  fmf_createAlloc( &gU, 1, nQP, dim * dim, 1 );
  gUm->nAlloc = -1;
  fmf_pretend( gUm, 1, nQP, dim, dim, gU->val );

  fmf_createAlloc( &gR, 1, nQP, dim, 1 );

  fmf_createAlloc( &fB, 1, nQP, dim, 1 );
  fmf_createAlloc( &fBgU, 1, nQP, 1, dim );

  fmf_createAlloc( &scalar1, 1, nQP, 1, 1 );

  if (mode == 1) {
    vecMV->val = FMF_PtrFirst( vecMV ) + offsetMV;

    fmf_createAlloc( &stMV, 1, 1, dim, nEP_mv );
    stMVv->nAlloc = -1;
    fmf_pretend( stMVv, 1, 1, nEP_mv * dim, 1, stMV->val );

    gcl_mv->nAlloc = -1;
    fmf_pretend( gcl_mv, 1, nQP, 1, nEP_mv * dim, vg_mv->bfGM->val0 );

    fmf_createAlloc( &divMV, 1, nQP, 1, 1 );
    fmf_createAlloc( &scalar2, 1, nQP, 1, 1 );

    fmf_createAlloc( &gMV, 1, nQP, dim * dim, 1 );
    gMVm->nAlloc = -1;
    fmf_pretend( gMVm, 1, nQP, dim, dim, gMV->val );

    fmf_createAlloc( &fBgMV, 1, nQP, 1, dim );
    fmf_createAlloc( &fBgMVgU, 1, nQP, 1, dim );
    fmf_createAlloc( &gMVgR, 1, nQP, dim, 1 );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

/*     output( "\n%d %d %d\n", mode, ii, iel ); */

    FMF_SetCell( out, ii );
    FMF_SetCell( vg_u->bfGM, iel );
    FMF_SetCell( vg_r->bfGM, iel );
    FMF_SetCell( vg_u->det, iel );
    FMF_SetCell( coef, ii );

    // u.
    ele_extractNodalValuesDBD( stU, stateU, conn_u + nEP_u * iel );
    // grad u.
    divgrad_act_g_m( gU, vg_u->bfGM, stUv );

    // b.
    ele_extractNodalValuesDBD( stB, stateB, conn_u + nEP_u * iel );
    bf_act( fB, bf_u, stB );

    // r.
    ele_extractNodalValuesDBD( stR, stateR, conn_r + nEP_r * iel );
    // grad r.
    fmf_mulAB_n1( gR, vg_r->bfGM, stRv );

    // b grad u.
    fmf_mulATBT_nn( fBgU, fB, gUm );

    // (grad r, b grad u).
    fmf_mulATBT_nn( scalar1, gR, fBgU );

    if (mode == 0) {
      fmf_sumLevelsMulF( out, scalar1, vg_u->det->val );

    } else if (mode == 1) {
      FMF_SetCell( vg_mv->bfGM, iel );
      FMF_SetCell( gcl_mv, iel );

      // nu.
      ele_extractNodalValuesDBD( stMV, vecMV, conn_mv + nEP_mv * iel );
      // div nu.
      fmf_mulAB_n1( divMV, gcl_mv, stMVv );
      // grad nu.
      divgrad_act_g_m( gMV, vg_mv->bfGM, stMVv );
      
      // b grad nu.
      fmf_mulATBT_nn( fBgMV, fB, gMVm );

      // (grad r, b grad u) div nu ...
      fmf_mulATB_nn( scalar2, divMV, scalar1 );

      // ... - (grad nu grad r, b grad u) ...
      fmf_mulATB_nn( gMVgR, gMVm, gR );
      fmf_mulATBT_nn( scalar1, gMVgR, fBgU );
      fmf_subAB_nn( scalar2, scalar2, scalar1 );

      // ... - (grad r, b grad nu grad u) ...
      fmf_mulABT_nn( fBgMVgU, fBgMV, gUm );
      fmf_mulATBT_nn( scalar1, gR, fBgMVgU );
      fmf_subAB_nn( scalar2, scalar2, scalar1 );

      fmf_sumLevelsMulF( out, scalar2, vg_u->det->val );
    }
/*     fmfc_save( out, "000", 0 ); */
/*     sys_pause(); */

    fmf_mulC( out, coef->val[0] );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &stU );
  fmf_freeDestroy( &stB );
  fmf_freeDestroy( &stR );
  fmf_freeDestroy( &gU );
  fmf_freeDestroy( &gR );
  fmf_freeDestroy( &scalar1 );
  fmf_freeDestroy( &fB );
  fmf_freeDestroy( &fBgU );
  if (mode == 1) { 
    fmf_freeDestroy( &stMV );
    fmf_freeDestroy( &divMV );
    fmf_freeDestroy( &gMV );
    fmf_freeDestroy( &scalar2 );
    fmf_freeDestroy( &fBgMV );
    fmf_freeDestroy( &fBgMVgU );
    fmf_freeDestroy( &gMVgR );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_sd_st_pspg_p"
/*!
  @par Revision history:
  - 25.10.2007, c
  - 29.10.2007
*/
int32 d_sd_st_pspg_p( FMField *out,
		      FMField *stateP, int32 offsetP,
		      FMField *stateR, int32 offsetR,
		      FMField *vecMV, int32 offsetMV,
		      FMField *coef,
		      VolumeGeometry *vg_p,
		      VolumeGeometry *vg_mv,
		      int32 *conn_p, int32 nEl_p, int32 nEP_p,
		      int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
		      int32 *elList, int32 elList_nRow,
		      int32 mode )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *stP = 0, *stR = 0, *stMV = 0, *gP = 0, *gR = 0, *gMV = 0;
  FMField *divMV = 0;
  FMField *gMVgR = 0, *gMVgP = 0;
  FMField *scalar1 = 0, *scalar2 = 0;
  FMField gcl_mv[1], stPv[1], stRv[1], stMVv[1];
  FMField gMVm[1];

  nQP = vg_p->bfGM->nLev;
  dim = vg_mv->bfGM->nRow;

  stateP->val = FMF_PtrFirst( stateP ) + offsetP;
  stateR->val = FMF_PtrFirst( stateR ) + offsetR;

  fmf_createAlloc( &stP, 1, 1, 1, nEP_p );
  stPv->nAlloc = -1;
  fmf_pretend( stPv, 1, 1, nEP_p, 1, stP->val );

  fmf_createAlloc( &stR, 1, 1, 1, nEP_p );
  stRv->nAlloc = -1;
  fmf_pretend( stRv, 1, 1, nEP_p, 1, stR->val );

  fmf_createAlloc( &gP, 1, nQP, dim, 1 );
  fmf_createAlloc( &gR, 1, nQP, dim, 1 );

  fmf_createAlloc( &scalar1, 1, nQP, 1, 1 );

  if (mode == 1) {
    vecMV->val = FMF_PtrFirst( vecMV ) + offsetMV;

    fmf_createAlloc( &stMV, 1, 1, dim, nEP_mv );
    stMVv->nAlloc = -1;
    fmf_pretend( stMVv, 1, 1, nEP_mv * dim, 1, stMV->val );

    gcl_mv->nAlloc = -1;
    fmf_pretend( gcl_mv, 1, nQP, 1, nEP_mv * dim, vg_mv->bfGM->val0 );

    fmf_createAlloc( &divMV, 1, nQP, 1, 1 );
    fmf_createAlloc( &scalar2, 1, nQP, 1, 1 );

    fmf_createAlloc( &gMV, 1, nQP, dim * dim, 1 );
    gMVm->nAlloc = -1;
    fmf_pretend( gMVm, 1, nQP, dim, dim, gMV->val );

    fmf_createAlloc( &gMVgP, 1, nQP, dim, 1 );
    fmf_createAlloc( &gMVgR, 1, nQP, dim, 1 );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

/*     output( "\n%d %d %d\n", mode, ii, iel ); */

    FMF_SetCell( out, ii );
    FMF_SetCell( vg_p->bfGM, iel );
    FMF_SetCell( vg_p->det, iel );
    FMF_SetCell( coef, ii );

    // p.
    ele_extractNodalValuesDBD( stP, stateP, conn_p + nEP_p * iel );
    // grad p.
    fmf_mulAB_n1( gP, vg_p->bfGM, stPv );

    // r.
    ele_extractNodalValuesDBD( stR, stateR, conn_p + nEP_p * iel );
    // grad r.
    fmf_mulAB_n1( gR, vg_p->bfGM, stRv );

    // (grad r, grad p).
    fmf_mulATB_nn( scalar1, gR, gP );

    if (mode == 0) {
      fmf_sumLevelsMulF( out, scalar1, vg_p->det->val );

    } else if (mode == 1) {
      FMF_SetCell( vg_mv->bfGM, iel );
      FMF_SetCell( gcl_mv, iel );

      // nu.
      ele_extractNodalValuesDBD( stMV, vecMV, conn_mv + nEP_mv * iel );
      // div nu.
      fmf_mulAB_n1( divMV, gcl_mv, stMVv );
      // grad nu.
      divgrad_act_g_m( gMV, vg_mv->bfGM, stMVv );

      // grad nu grad r.
      fmf_mulATB_nn( gMVgR, gMVm, gR );

      // grad nu grad p.
      fmf_mulATB_nn( gMVgP, gMVm, gP );
      
      // (grad r, grad p) div nu ...
      fmf_mulATB_nn( scalar2, divMV, scalar1 );

      // ... - (grad nu grad r, grad p) ...
      fmf_mulATB_nn( scalar1, gMVgR, gP );
      fmf_subAB_nn( scalar2, scalar2, scalar1 );

      // ... - (grad nu grad p, grad r) ...
      fmf_mulATB_nn( scalar1, gMVgP, gR );
      fmf_subAB_nn( scalar2, scalar2, scalar1 );

      fmf_sumLevelsMulF( out, scalar2, vg_p->det->val );
    }
/*     fmfc_save( out, "000", 0 ); */
/*     sys_pause(); */

    fmf_mulC( out, coef->val[0] );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &stP );
  fmf_freeDestroy( &stR );
  fmf_freeDestroy( &gP );
  fmf_freeDestroy( &gR );
  fmf_freeDestroy( &scalar1 );
  if (mode == 1) { 
    fmf_freeDestroy( &stMV );
    fmf_freeDestroy( &divMV );
    fmf_freeDestroy( &gMV );
    fmf_freeDestroy( &scalar2 );
    fmf_freeDestroy( &gMVgP );
    fmf_freeDestroy( &gMVgR );
  }

  return( ret );
}
