#include "terms_navier_stokes.h"
#include "terms_adj_navier_stokes.h"
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
                       Mapping *vg, int32 isDiff )
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
    FMF_SetCell( vg->det, ii );
    FMF_SetCellX1( vg->bf, ii );

    if (isDiff) {
      bf_ract( gf, vg->bf, gradU );
      bf_actt( ftgf, vg->bf, gf );
      fmf_sumLevelsMulF( out, ftgf, vg->det->val );
    } else {
      FMF_SetCell( stateW, ii );
      fmf_mulAB_nn( gfu, gradU, stateW );
      bf_actt( ftgfu, vg->bf, gfu );
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
                       Mapping *vg, int32 isDiff )
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
      FMF_SetCellX1( vg->bf, ii );
      bf_actt( ftvtg, vg->bf, vtg );

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
			FMField *coef, Mapping *vg,
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
    FMF_SetCellX1( coef, ii );
    FMF_SetCellX1( vg->bf, ii );

    // u grad u.
    fmf_mulAB_nn( gUfU, gradU, stateU );
    // (u grad u, grad).
    convect_build_vtbg( gUfUTgT, vg->bfGM, gUfU );
    // (u grad u, v grad).
    bf_actt( fTgUfUTgT, vg->bf, gUfUTgT );

    // u grad.
    convect_build_vtg( fUTg, vg->bfGM, stateU );
    // (grad u, u^T grad).
    fmf_mulAB_nn( gUfUTg, gradU, fUTg );
    // (v grad u, u grad).
    bf_actt( fTgUfUTg, vg->bf, gUfUTg );

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
			 FMField *coef, Mapping *vg_w,
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
    FMF_SetCellX1( coef, ii );
    FMF_SetCellX1( vg_w->bf, ii );

    // (grad p, grad).
    convect_build_vtbg( gPTgT, vg_w->bfGM, gradP );
    // (grad p, v grad).
    bf_actt( fTgPTgT, vg_w->bf, gPTgT );

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
			 FMField *coef,
			 Mapping *vg_u, Mapping *vg_r,
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
    FMF_SetCell( vg_r->bfGM, ii );
    FMF_SetCell( vg_u->det, ii );
    FMF_SetCellX1( coef, ii );
    FMF_SetCellX1( vg_u->bf, ii );

    // (grad u, grad).
    fmf_mulATB_nn( gUTg, gradU, vg_r->bfGM );

    // (v grad u, grad).
    bf_actt( fTgUTg, vg_u->bf, gUTg );

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
		      FMField *viscosity, Mapping *vg )
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
    FMF_SetCellX1( viscosity, ii );
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
			    Mapping *sg, int32 isDiff )
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
      FMF_SetCellX1( sg->bf, ii );
      fmf_sumLevelsTMulF( out, sg->bf, sg->det->val );

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
int32 d_sd_div( FMField *out, FMField *divU, FMField *gradU,
                FMField *stateP, FMField *divMV, FMField *gradMV,
                Mapping *vg_u, int32 mode )
{
  int32 ii, nQP, ret = RET_OK;
  FMField *aux11 = 0;

  nQP = vg_u->bfGM->nLev;

  fmf_createAlloc( &aux11, 1, nQP, 1, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( stateP, ii );
    FMF_SetCell( divU, ii );
    FMF_SetCell( vg_u->det, ii );

    fmf_mulAB_nn( aux11, stateP, divU );

    if (mode == 1) {
      FMF_SetCell( gradU, ii );
      FMF_SetCell( divMV, ii );
      FMF_SetCell( gradMV, ii );

      fmf_mul( aux11, divMV->val );
      sub_mul_gradddgrad_scalar( aux11, gradU, gradMV, stateP );
    }
    fmf_sumLevelsMulF( out, aux11, vg_u->det->val );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &aux11 );

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
int32 d_sd_div_grad( FMField *out, FMField *gradU, FMField *gradW,
                     FMField *divMV, FMField *gradMV, FMField *viscosity,
		     Mapping *vg_u, int32 mode )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *uvel = 0, *aux = 0, *aux1 = 0, *aux2 = 0, *aux3 = 0;
  FMField gum[1], gwm[1], gmvm[1], aux3m[1];

  nQP = vg_u->bfGM->nLev;
  dim = vg_u->bfGM->nRow;

  fmf_createAlloc( &uvel, 1, nQP, 1, 1 );

  if (mode == 1) {
    fmf_createAlloc( &aux, 1, 1, 1, 1 );
    fmf_createAlloc( &aux1, 1, nQP, 1, 1 );
    fmf_createAlloc( &aux2, 1, nQP, 1, 1 );
    fmf_createAlloc( &aux3, 1, nQP, dim * dim, 1 );
    aux3m->nAlloc = -1;
    fmf_pretend( aux3m, 1, nQP, dim, dim, aux3->val );

    gum->nAlloc = -1;
    fmf_pretend( gum, gradU->nCell, nQP, dim, dim, gradU->val0 );

    gwm->nAlloc = -1;
    fmf_pretend( gwm, gradW->nCell, nQP, dim, dim, gradW->val0 );

    gmvm->nAlloc = -1;
    fmf_pretend( gmvm, gradMV->nCell, nQP, dim, dim, gradMV->val0 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( gradU, ii );
    FMF_SetCell( gradW, ii );
    FMF_SetCellX1( viscosity, ii );
    FMF_SetCell( vg_u->det, ii );

    // (gu:gw)
    fmf_mulATB_nn( uvel, gradW, gradU );

    if (mode == 0) {
      fmf_mul( uvel, viscosity->val );
      fmf_sumLevelsMulF( out, uvel, vg_u->det->val );

    } else if (mode == 1) {
      FMF_SetCell( divMV, ii );
      FMF_SetCell( gum, ii );
      FMF_SetCell( gwm, ii );
      FMF_SetCell( gmvm, ii );

      // (gu:gw) div nu
      fmf_mulAB_nn( aux1, uvel, divMV );
      fmf_mul( aux1, viscosity->val );
      fmf_sumLevelsMulF( out, aux1, vg_u->det->val );

      // (gu * gnu) : gw
      fmf_mulAB_nn( aux3m, gum, gmvm );
      fmf_mulATB_nn( aux1, aux3, gradW );

      // (gw * gnu) : gu
      fmf_mulAB_nn( aux3m, gwm, gmvm );
      fmf_mulATB_nn( aux2, aux3, gradU );

      fmf_addAB_nn( aux1, aux1, aux2 );
      fmf_mul( aux1, viscosity->val );
      fmf_sumLevelsMulF( aux, aux1, vg_u->det->val );
      fmf_subAB_nn( out, out, aux );
    }

    ERR_CheckGo( ret );
  }

 end_label:
  if (mode == 1) {
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
int32 d_sd_convect( FMField *out, FMField *stateU, FMField *gradU,
		    FMField *stateW, FMField *divMV, FMField *gradMV,
		    Mapping *vg_u, int32 mode )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *aux11 = 0, *aux = 0, *fwgu = 0, *fwgugmv = 0;
  FMField gum[1], gmvm[1];

  nQP = vg_u->bfGM->nLev;
  dim = vg_u->bfGM->nRow;

  gum->nAlloc = -1;
  fmf_pretend( gum, gradU->nCell, nQP, dim, dim, gradU->val0 );

  fmf_createAlloc( &fwgu, 1, nQP, 1, dim );
  fmf_createAlloc( &aux11, 1, nQP, 1, 1 );

  if (mode == 1) {
    gmvm->nAlloc = -1;
    fmf_pretend( gmvm, gradMV->nCell, nQP, dim, dim, gradMV->val0 );

    fmf_createAlloc( &fwgugmv, 1, nQP, 1, dim );
    fmf_createAlloc( &aux, 1, nQP, 1, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( stateU, ii );
    FMF_SetCell( stateW, ii );
    FMF_SetCell( gum, ii );
    FMF_SetCell( vg_u->det, ii );

    fmf_mulATB_nn( fwgu, stateW, gum );
    fmf_mulAB_nn( aux11, fwgu, stateU );

    if (mode == 1) {
      FMF_SetCell( divMV, ii );
      FMF_SetCell( gmvm, ii );

      fmf_mul( aux11, divMV->val );

      fmf_mulAB_nn( fwgugmv, fwgu, gmvm );
      fmf_mulAB_nn( aux, fwgugmv, stateU );
      fmf_subAB_nn( aux11, aux11, aux );
    }
    fmf_sumLevelsMulF( out, aux11, vg_u->det->val );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &aux11 );
  fmf_freeDestroy( &fwgu );
  if (mode == 1) {
    fmf_freeDestroy( &fwgugmv );
    fmf_freeDestroy( &aux );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_sd_volume_dot"
/*!
  mode == 0: \int pq
  mode == 1: \int pq div V

  Works for both scalars and vectors.

  @par Revision history:
  - 24.02.2006, c
  - 27.02.2006
*/
int32 d_sd_volume_dot( FMField *out, FMField *stateP, FMField *stateQ,
                       FMField *divMV, Mapping *vg, int32 mode )
{
  int32 ii, nQP, ret = RET_OK;
  FMField *pq = 0;

  nQP = vg->bfGM->nLev;

  fmf_createAlloc( &pq, 1, nQP, 1, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( stateP, ii );
    FMF_SetCell( stateQ, ii );
    FMF_SetCell( vg->det, ii );

    fmf_mulATB_nn( pq, stateP, stateQ );

    if (mode == 1) {
      FMF_SetCell( divMV, ii );

      fmf_mul( pq, divMV->val );
    }
    fmf_sumLevelsMulF( out, pq, vg->det->val );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &pq );

  return( ret );
}


#undef __FUNC__
#define __FUNC__ "d_sd_st_grad_div"
/*!
  @par Revision history:
  - 24.10.2007, c
*/
int32 d_sd_st_grad_div( FMField *out, FMField *divU, FMField *gradU,
			FMField *divW, FMField *gradW, FMField *divMV,
			FMField *gradMV, FMField *coef,
			Mapping *vg_u, int32 mode )
{
  int32 ii, nQP, ret = RET_OK;
  FMField *scalar1 = 0, *scalar2 = 0;

  nQP = vg_u->bfGM->nLev;

  fmf_createAlloc( &scalar1, 1, nQP, 1, 1 );

  if (mode == 1) {
    fmf_createAlloc( &scalar2, 1, nQP, 1, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCellX1( coef, ii );
    FMF_SetCell( divU, ii );
    FMF_SetCell( divW, ii );
    FMF_SetCell( vg_u->det, ii );

    if (mode == 0) {
      // (div u, div w)
      fmf_mulATB_nn( scalar1, divU, divW );
      fmf_mul( scalar1, coef->val );
      fmf_sumLevelsMulF( out, scalar1, vg_u->det->val );

    } else if (mode == 1) {
      FMF_SetCell( divMV, ii );
      FMF_SetCell( gradU, ii );
      FMF_SetCell( gradW, ii );
      FMF_SetCell( gradMV, ii );

      // div u div w div nu ...
      fmf_mulATB_nn( scalar1, divU, divW );
      fmf_mulATB_nn( scalar2, divMV, scalar1 );

      // ... - (gu : gnu^T) div w ...
      sub_mul_gradddgrad_scalar( scalar2, gradMV, gradU, divW );

      // ... - (gw : gnu^T) div u
      sub_mul_gradddgrad_scalar( scalar2, gradMV, gradW, divU );

      fmf_mul( scalar2, coef->val );
      fmf_sumLevelsMulF( out, scalar2, vg_u->det->val );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &scalar1 );
  if (mode == 1) {
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
int32 d_sd_st_supg_c( FMField *out, FMField *stateB, FMField *gradU,
                      FMField *gradW, FMField *divMV, FMField *gradMV,
		      FMField *coef, Mapping *vg_u, int32 mode )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *fBgU = 0, *fBgW = 0, *fBgMV = 0, *fBgMVgU = 0, *fBgMVgW = 0;
  FMField *scalar1 = 0, *scalar2 = 0;

  nQP = vg_u->bfGM->nLev;
  dim = vg_u->bfGM->nRow;

  fmf_createAlloc( &fBgU, 1, nQP, 1, dim );
  fmf_createAlloc( &fBgW, 1, nQP, 1, dim );

  fmf_createAlloc( &scalar1, 1, nQP, 1, 1 );

  if (mode == 1) {
    fmf_createAlloc( &scalar2, 1, nQP, 1, 1 );

    fmf_createAlloc( &fBgMV, 1, nQP, 1, dim );
    fmf_createAlloc( &fBgMVgU, 1, nQP, 1, dim );
    fmf_createAlloc( &fBgMVgW, 1, nQP, 1, dim );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( gradU, ii );
    FMF_SetCell( stateB, ii );
    FMF_SetCell( gradW, ii );
    FMF_SetCell( vg_u->det, ii );
    FMF_SetCellX1( coef, ii );

    // b grad u.
    fmf_mulATBT_nn( fBgU, stateB, gradU );

    // b grad w.
    fmf_mulATBT_nn( fBgW, stateB, gradW );

    // (b grad u, b grad w).
    fmf_mulABT_nn( scalar1, fBgU, fBgW );

    if (mode == 0) {
      fmf_mul( scalar1, coef->val );
      fmf_sumLevelsMulF( out, scalar1, vg_u->det->val );

    } else if (mode == 1) {
      FMF_SetCell( divMV, ii );
      FMF_SetCell( gradMV, ii );

      // b grad nu.
      fmf_mulATBT_nn( fBgMV, stateB, gradMV );

      // (b grad u, b grad w) div nu ...
      fmf_mulATB_nn( scalar2, divMV, scalar1 );

      // ... - (b grad nu grad u, b grad w) ...
      fmf_mulABT_nn( fBgMVgU, fBgMV, gradU );
      fmf_mulABT_nn( scalar1, fBgMVgU, fBgW );
      fmf_subAB_nn( scalar2, scalar2, scalar1 );

      // ... - (b grad nu grad w, b grad u) ...
      fmf_mulABT_nn( fBgMVgW, fBgMV, gradW );
      fmf_mulABT_nn( scalar1, fBgMVgW, fBgU );
      fmf_subAB_nn( scalar2, scalar2, scalar1 );

      fmf_mul( scalar2, coef->val );
      fmf_sumLevelsMulF( out, scalar2, vg_u->det->val );
    }

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &scalar1 );
  fmf_freeDestroy( &fBgU );
  fmf_freeDestroy( &fBgW );
  if (mode == 1) {
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
int32 d_sd_st_pspg_c( FMField *out, FMField *stateB, FMField *gradU,
		      FMField *gradR, FMField *divMV, FMField *gradMV,
		      FMField *coef, Mapping *vg_u, int32 mode )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *fBgU = 0, *fBgMV = 0, *fBgMVgU = 0, *gMVgR = 0;
  FMField *scalar1 = 0, *scalar2 = 0;

  nQP = vg_u->bfGM->nLev;
  dim = vg_u->bfGM->nRow;

  fmf_createAlloc( &fBgU, 1, nQP, 1, dim );
  fmf_createAlloc( &scalar1, 1, nQP, 1, 1 );

  if (mode == 1) {
    fmf_createAlloc( &scalar2, 1, nQP, 1, 1 );

    fmf_createAlloc( &fBgMV, 1, nQP, 1, dim );
    fmf_createAlloc( &fBgMVgU, 1, nQP, 1, dim );
    fmf_createAlloc( &gMVgR, 1, nQP, dim, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( stateB, ii );
    FMF_SetCell( gradU, ii );
    FMF_SetCell( gradR, ii );
    FMF_SetCell( vg_u->det, ii );
    FMF_SetCellX1( coef, ii );

    // b grad u.
    fmf_mulATBT_nn( fBgU, stateB, gradU );

    // (grad r, b grad u).
    fmf_mulATBT_nn( scalar1, gradR, fBgU );

    if (mode == 0) {
      fmf_mul( scalar1, coef->val );
      fmf_sumLevelsMulF( out, scalar1, vg_u->det->val );

    } else if (mode == 1) {
      FMF_SetCell( divMV, ii );
      FMF_SetCell( gradMV, ii );

      // b grad nu.
      fmf_mulATBT_nn( fBgMV, stateB, gradMV );

      // (grad r, b grad u) div nu ...
      fmf_mulATB_nn( scalar2, divMV, scalar1 );

      // ... - (grad nu grad r, b grad u) ...
      fmf_mulATB_nn( gMVgR, gradMV, gradR );
      fmf_mulATBT_nn( scalar1, gMVgR, fBgU );
      fmf_subAB_nn( scalar2, scalar2, scalar1 );

      // ... - (grad r, b grad nu grad u) ...
      fmf_mulABT_nn( fBgMVgU, fBgMV, gradU );
      fmf_mulATBT_nn( scalar1, gradR, fBgMVgU );
      fmf_subAB_nn( scalar2, scalar2, scalar1 );

      fmf_mul( scalar2, coef->val );
      fmf_sumLevelsMulF( out, scalar2, vg_u->det->val );
    }

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &scalar1 );
  fmf_freeDestroy( &fBgU );
  if (mode == 1) {
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
int32 d_sd_st_pspg_p( FMField *out, FMField *gradR, FMField *gradP,
		      FMField *divMV, FMField *gradMV, FMField *coef,
		      Mapping *vg_p, int32 mode )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *gMVgR = 0, *gMVgP = 0;
  FMField *scalar1 = 0, *scalar2 = 0;

  nQP = vg_p->bfGM->nLev;
  dim = gradR->nRow;

  fmf_createAlloc( &scalar1, 1, nQP, 1, 1 );

  if (mode == 1) {
    fmf_createAlloc( &scalar2, 1, nQP, 1, 1 );

    fmf_createAlloc( &gMVgP, 1, nQP, dim, 1 );
    fmf_createAlloc( &gMVgR, 1, nQP, dim, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( gradR, ii );
    FMF_SetCell( gradP, ii );
    FMF_SetCell( vg_p->det, ii );
    FMF_SetCellX1( coef, ii );

    // (grad r, grad p).
    fmf_mulATB_nn( scalar1, gradR, gradP );

    if (mode == 0) {
      fmf_mul( scalar1, coef->val );
      fmf_sumLevelsMulF( out, scalar1, vg_p->det->val );

    } else if (mode == 1) {
      FMF_SetCell( divMV, ii );
      FMF_SetCell( gradMV, ii );

      // grad nu grad r.
      fmf_mulATB_nn( gMVgR, gradMV, gradR );

      // grad nu grad p.
      fmf_mulATB_nn( gMVgP, gradMV, gradP );

      // (grad r, grad p) div nu ...
      fmf_mulATB_nn( scalar2, divMV, scalar1 );

      // ... - (grad nu grad r, grad p) ...
      fmf_mulATB_nn( scalar1, gMVgR, gradP );
      fmf_subAB_nn( scalar2, scalar2, scalar1 );

      // ... - (grad nu grad p, grad r) ...
      fmf_mulATB_nn( scalar1, gMVgP, gradR );
      fmf_subAB_nn( scalar2, scalar2, scalar1 );

      fmf_mul( scalar2, coef->val );
      fmf_sumLevelsMulF( out, scalar2, vg_p->det->val );
    }

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &scalar1 );
  if (mode == 1) {
    fmf_freeDestroy( &scalar2 );
    fmf_freeDestroy( &gMVgP );
    fmf_freeDestroy( &gMVgR );
  }

  return( ret );
}
