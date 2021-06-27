#include "terms_diffusion.h"
#include "terms.h"

#undef __FUNC__
#define __FUNC__ "laplace_build_gtg"
/*!
  @par Revision history:
  - 28.11.2005, c
  - 30.05.2007
*/
int32 laplace_build_gtg( FMField *out, FMField *gc )
{
  int32 iqp, ir, ic, nEP, nQP, nCol;
  float64 *pout, *pg1, *pg2, *pg3;

  nEP = gc->nCol;
  nQP = gc->nLev;
  nCol = out->nCol;

  fmf_fillC( out, 0.0 );
  switch (gc->nRow) {
  case 3:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      pg3 = pg2 + nEP;

      pout = FMF_PtrLevel( out, iqp );

      for (ir = 0; ir < nEP; ir++) {
	for (ic = 0; ic < nEP; ic++) {
	  pout[ic] = pg1[ir] * pg1[ic] + pg2[ir] * pg2[ic] + pg3[ir] * pg3[ic];
	}
	pout += nCol;
      }
    }
    break;

  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;

      pout = FMF_PtrLevel( out, iqp );

      for (ir = 0; ir < nEP; ir++) {
	for (ic = 0; ic < nEP; ic++) {
	  pout[ic] = pg1[ir] * pg1[ic] + pg2[ir] * pg2[ic];
	}
	pout += nCol;
      }
    }
    break;

  case 1:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel(gc, iqp);

      pout = FMF_PtrLevel(out, iqp);

      for (ir = 0; ir < nEP; ir++) {
        for (ic = 0; ic < nEP; ic++) {
          pout[ic] = pg1[ir] * pg1[ic];
        }
        pout += nCol;
      }
    }
    break;

  default:
    errput( ErrHead "ERR_Switch\n" );
    return( RET_Fail );
  }
  return( RET_OK );
}


#undef __FUNC__
#define __FUNC__ "laplace_act_g_m"
/*!
  @par Revision history:
  - 28.11.2005, c
  - 30.05.2007
*/
int32 laplace_act_g_m( FMField *out, FMField *gc, FMField *mtx )
{
  int32 iqp, ic, ik, nEP, nQP, nCol;
  float64 val1, val2, val3;
  float64 *pout, *pmtx, *pg1, *pg2, *pg3;

  nEP = gc->nCol;
  nQP = gc->nLev;
  nCol = mtx->nCol;

/*   output( "%d %d %d %d\n", nEP, nQP, nCol, dim ); */

  switch (gc->nRow) {
  case 3:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      pg3 = pg2 + nEP;
      pout = FMF_PtrLevel( out, iqp );

      if (mtx->nLev == nQP) {
	pmtx = FMF_PtrLevel( mtx, iqp );
      } else {
	pmtx = FMF_PtrCurrent( mtx );
      }
      for (ic = 0; ic < nCol; ic++) {
	val1 = val2 = val3 = 0.0;
	for (ik = 0; ik < nEP; ik++) {
/* 	    output( "%d %d %d %d\n", iqp, ic, ik ); */
	  val1 += pg1[ik] * pmtx[ic+nCol*ik];
	  val2 += pg2[ik] * pmtx[ic+nCol*ik];
	  val3 += pg3[ik] * pmtx[ic+nCol*ik];
	}
	pout[ic+0] = val1;
	pout[ic+1] = val2;
	pout[ic+2] = val3;
      }
    }
    break;

  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      pout = FMF_PtrLevel( out, iqp );

      if (mtx->nLev == nQP) {
	pmtx = FMF_PtrLevel( mtx, iqp );
      } else {
	pmtx = FMF_PtrCurrent( mtx );
      }
      for (ic = 0; ic < nCol; ic++) {
	val1 = val2 = 0.0;
	for (ik = 0; ik < nEP; ik++) {
/* 	    output( "%d %d %d %d\n", iqp, ic, ik ); */
	  val1 += pg1[ik] * pmtx[ic+nCol*ik];
	  val2 += pg2[ik] * pmtx[ic+nCol*ik];
	}
	pout[ic+0] = val1;
	pout[ic+1] = val2;
      }
    }
    break;

  case 1:
      for (iqp = 0; iqp < nQP; iqp++) {
          pg1 = FMF_PtrLevel(gc, iqp);
          pout = FMF_PtrLevel(out, iqp);

          if (mtx->nLev == nQP)
              pmtx = FMF_PtrLevel(mtx, iqp);
          else
              pmtx = FMF_PtrCurrent(mtx);

          for (ic = 0; ic < nCol; ic++) {
              val1 = 0.0;
              for (ik = 0; ik < nEP; ik++) {
                  val1 += pg1[ik] * pmtx[ic + nCol*ik];
              }
              pout[ic + 0] = val1;
          }
      }
      break;

  default:
    errput( ErrHead "ERR_Switch\n" );
    return( RET_Fail );
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "laplace_act_gt_m"
/*!
  @par Revision history:
  - 28.11.2005, c
  - 30.05.2007
*/
int32 laplace_act_gt_m( FMField *out, FMField *gc, FMField *mtx )
{
  int32 iqp, iep, ii, nEP, nQP, nCol;
  float64 *pout, *pmtx, *pg1, *pg2, *pg3;

  nEP = gc->nCol;
  nQP = gc->nLev;
  nCol = mtx->nCol;

  switch (gc->nRow) {
  case 3:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      pg3 = pg2 + nEP;

      pmtx = FMF_PtrLevel( mtx, iqp );
      for (iep = 0; iep < nEP; iep++) {
	pout = FMF_PtrLevel( out, iqp ) + nCol * iep;
	for (ii = 0; ii < nCol; ii++) {
	  pout[ii]
	    = pg1[iep] * pmtx[0*nCol+ii]
	    + pg2[iep] * pmtx[1*nCol+ii]
	    + pg3[iep] * pmtx[2*nCol+ii];
	}
      }
    }
    break;

  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;

      pmtx = FMF_PtrLevel( mtx, iqp );
      for (iep = 0; iep < nEP; iep++) {
	pout = FMF_PtrLevel( out, iqp ) + nCol * iep;
	for (ii = 0; ii < nCol; ii++) {
	  pout[ii]
	    = pg1[iep] * pmtx[0*nCol+ii]
	    + pg2[iep] * pmtx[1*nCol+ii];
	}
      }
    }
    break;

  case 1:
      for (iqp = 0; iqp < nQP; iqp++) {
          pg1 = FMF_PtrLevel(gc, iqp);

          pmtx = FMF_PtrLevel(mtx, iqp);

          for (iep = 0; iep < nEP; iep++) {
              pout = FMF_PtrLevel(out, iqp) + nCol * iep;

              for (ii = 0; ii < nCol; ii++) {
                  pout[ii] = pg1[iep] * pmtx[0 * nCol + ii];
              }
          }
      }
      break;

  default:
    errput( ErrHead "ERR_Switch\n" );
    return( RET_Fail );
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "dw_laplace"
/*!
  @par Revision history:
  - 28.11.2005, c
  - 09.12.2005
*/
int32 dw_laplace( FMField *out, FMField *grad,
		  FMField *coef, Mapping *vg,
		  int32 isDiff )
{
  int32 ii, nQP, nEP, ret = RET_OK;
  FMField *gtg = 0, *gtgu = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;

  if (isDiff) {
    fmf_createAlloc( &gtg, 1, nQP, nEP, nEP );
  } else {
    fmf_createAlloc( &gtgu, 1, nQP, nEP, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCellX1( coef, ii );

    if (isDiff) {
      laplace_build_gtg( gtg, vg->bfGM );
      fmf_mulAF( gtg, gtg, coef->val );
      fmf_sumLevelsMulF( out, gtg, vg->det->val );
    } else {
      FMF_SetCell( grad, ii );
      laplace_act_gt_m( gtgu, vg->bfGM, grad );
      fmf_mulAF( gtgu, gtgu, coef->val );
      fmf_sumLevelsMulF( out, gtgu, vg->det->val );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &gtg );
  } else {
    fmf_freeDestroy( &gtgu );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_laplace"
int32 d_laplace( FMField *out, FMField *gradP1, FMField *gradP2,
		 FMField *coef, Mapping *vg )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *dgp2 = 0, *gp1tdgp2 = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  fmf_createAlloc( &dgp2, 1, nQP, dim, 1 );
  fmf_createAlloc( &gp1tdgp2, 1, nQP, 1, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCell( gradP1, ii );
    FMF_SetCell( gradP2, ii );
    FMF_SetCellX1( coef, ii );

    fmf_mulAF( dgp2, gradP2, coef->val );
    fmf_mulATB_nn( gp1tdgp2, gradP1, dgp2 );
    fmf_sumLevelsMulF( out, gp1tdgp2, vg->det->val );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &dgp2 );
  fmf_freeDestroy( &gp1tdgp2 );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_diffusion"
/*!
  @par Revision history:
  - c: 03.08.2006, r: 23.01.2008
*/
int32 dw_diffusion( FMField *out, FMField *grad,
		    FMField *mtxD, Mapping *vg,
		    int32 isDiff )
{
  int32 ii, dim, nQP, nEP, ret = RET_OK;
  FMField *gtd = 0, *gtdg = 0, *dgp = 0, *gtdgp = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;
  dim = vg->bfGM->nRow;

  if (isDiff) {
    fmf_createAlloc( &gtd, 1, nQP, nEP, dim );
    fmf_createAlloc( &gtdg, 1, nQP, nEP, nEP );
  } else {
    fmf_createAlloc( &dgp, 1, nQP, dim, 1 );
    fmf_createAlloc( &gtdgp, 1, nQP, nEP, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCellX1( mtxD, ii );

    if (isDiff) {
      fmf_mulATB_nn( gtd, vg->bfGM, mtxD );
      fmf_mulAB_nn( gtdg, gtd, vg->bfGM );
      fmf_sumLevelsMulF( out, gtdg, vg->det->val );
    } else {
      FMF_SetCell( grad, ii );
      fmf_mulAB_nn( dgp, mtxD, grad );
      fmf_mulATB_nn( gtdgp, vg->bfGM, dgp );
      fmf_sumLevelsMulF( out, gtdgp, vg->det->val );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &gtd );
    fmf_freeDestroy( &gtdg );
  } else {
    fmf_freeDestroy( &dgp );
    fmf_freeDestroy( &gtdgp );
  }

  return( ret );
}


#undef __FUNC__
#define __FUNC__ "d_diffusion"
/*!
  @par Revision history:
  - c: 12.03.2007, r: 23.01.2008
*/
int32 d_diffusion( FMField *out, FMField *gradP1, FMField *gradP2,
		   FMField *mtxD, Mapping *vg )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *dgp2 = 0, *gp1tdgp2 = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  fmf_createAlloc( &dgp2, 1, nQP, dim, 1 );
  fmf_createAlloc( &gp1tdgp2, 1, nQP, 1, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCell( gradP1, ii );
    FMF_SetCell( gradP2, ii );
    FMF_SetCellX1( mtxD, ii );

    fmf_mulAB_nn( dgp2, mtxD, gradP2 );
    fmf_mulATB_nn( gp1tdgp2, gradP1, dgp2 );
    fmf_sumLevelsMulF( out, gp1tdgp2, vg->det->val );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &dgp2 );
  fmf_freeDestroy( &gp1tdgp2 );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_sd_diffusion"
int32 d_sd_diffusion(FMField *out,
                     FMField *grad_q, FMField *grad_p,
                     FMField *grad_w, FMField *div_w,
                     FMField *mtxD, Mapping *vg)
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
    FMF_SetCellX1( mtxD, ii );
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
#define __FUNC__ "dw_diffusion_r"
/*!
  @par Revision history:
  - c: 23.04.2007, r: 23.01.2008
*/
int32 dw_diffusion_r( FMField *out, FMField *mtxD, Mapping *vg )
{
  int32 ii, nQP, nEP, ret = RET_OK;
  FMField *gtd = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;

  fmf_createAlloc( &gtd, 1, nQP, nEP, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCellX1( mtxD, ii );

    fmf_mulATB_nn( gtd, vg->bfGM, mtxD );
    fmf_sumLevelsMulF( out, gtd, vg->det->val );
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &gtd );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_surface_flux"
int32 d_surface_flux( FMField *out, FMField *grad,
                      FMField *mtxD, Mapping *sg, int32 mode )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *dgp = 0, *ntdgp = 0;

  nQP = sg->normal->nLev;
  dim = sg->normal->nRow;

  fmf_createAlloc( &dgp, 1, nQP, dim, 1 );
  fmf_createAlloc( &ntdgp, 1, nQP, 1, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( grad, ii );
    FMF_SetCellX1( mtxD, ii );
    FMF_SetCell( sg->normal, ii );
    FMF_SetCell( sg->det, ii );

    fmf_mulAB_nn( dgp, mtxD, grad );
    fmf_mulATB_nn( ntdgp, sg->normal, dgp );

    fmf_sumLevelsMulF( out, ntdgp, sg->det->val );
    if (mode == 1) {
      FMF_SetCell( sg->volume, ii );
      fmf_mulC( out, 1.0 / sg->volume->val[0] );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &dgp );
  fmf_freeDestroy( &ntdgp );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_surface_flux"
int32 dw_surface_flux(FMField *out, FMField *grad,
                      FMField *mat, FMField *bf, Mapping *sg,
                      int32 *fis, int32 nFa, int32 nFP, int32 mode)
{
  int32 ii, ifa, dim, nQP, nEP, ret = RET_OK;
  FMField *ntk = 0, *ntkg = 0, *out_qp = 0;

  nQP = sg->normal->nLev;
  dim = sg->normal->nRow;
  nEP = sg->bfGM->nCol;

  fmf_createAlloc(&ntk, 1, nQP, 1, dim);
  if (mode) {
    fmf_createAlloc(&ntkg, 1, nQP, 1, nEP);
    fmf_createAlloc(&out_qp, 1, nQP, nEP, nEP);
  } else {
    fmf_createAlloc(&ntkg, 1, nQP, 1, 1);
    fmf_createAlloc(&out_qp, 1, nQP, nEP, 1);
  }

  for (ii = 0; ii < out->nCell; ii++) {
    ifa = fis[ii*nFP+1];

    FMF_SetCell(out, ii);
    FMF_SetCellX1(mat, ii);
    FMF_SetCell(sg->det, ii);
    FMF_SetCell(sg->normal, ii);
    FMF_SetCell(bf, ifa);

    fmf_mulATB_nn(ntk, sg->normal, mat);

    if (mode) {
      FMF_SetCell(sg->bfGM, ii);

      fmf_mulAB_nn(ntkg, ntk, sg->bfGM);
    } else {
      FMF_SetCell(grad, ii);

      fmf_mulAB_nn(ntkg, ntk, grad);
    }
    fmf_mulATB_nn(out_qp, bf, ntkg);
    fmf_sumLevelsMulF(out, out_qp, sg->det->val);

    ERR_CheckGo(ret);
  }

 end_label:
  fmf_freeDestroy(&ntk);
  fmf_freeDestroy(&ntkg);
  fmf_freeDestroy(&out_qp);

  return(ret);
}

#undef __FUNC__
#define __FUNC__ "dw_convect_v_grad_s"
int32 dw_convect_v_grad_s( FMField *out, FMField *val_v, FMField *grad_s,
                           Mapping *vvg, Mapping *svg,
                           int32 isDiff )
{
  int32 ii, nEPV, nEPS, dim, nQP, ret = RET_OK;
  FMField *aux = 0, *out_qp = 0, *gst = 0;

  nQP = vvg->bfGM->nLev;
  dim = vvg->bfGM->nRow;
  nEPS = svg->bfGM->nCol;
  nEPV = vvg->bf->nCol;

  if (isDiff == 0) {
    fmf_createAlloc( &aux, 1, nQP, 1, 1 );
    fmf_createAlloc( &out_qp, 1, nQP, nEPS, 1 );

  } else if (isDiff == 1) { // d/ds.
    fmf_createAlloc( &aux, 1, nQP, 1, nEPS );
    fmf_createAlloc( &out_qp, 1, nQP, nEPS, nEPS );

  } else { // d/dv.
    fmf_createAlloc( &aux, 1, nQP, 1, dim * nEPV );
    fmf_createAlloc( &out_qp, 1, nQP, nEPS, dim * nEPV );
    fmf_createAlloc( &gst, 1, nQP, 1, dim );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCellX1( svg->bf, ii );

    if (isDiff == 0) {
      FMF_SetCell( val_v, ii );
      FMF_SetCell( grad_s, ii );

      // v^T Psi^T G_c s.
      fmf_mulATB_nn( aux, val_v, grad_s );
      // Phi^T v^T Psi^T G_c s.
      fmf_mulATB_nn( out_qp, svg->bf, aux );

    } else if (isDiff == 1) { // d/ds.
      FMF_SetCell( val_v, ii );
      FMF_SetCell( svg->bfGM, ii );

      // v^T Psi^T G_c.
      fmf_mulATB_nn( aux, val_v, svg->bfGM );
      // Phi^T v^T Psi^T G_c.
      fmf_mulATB_nn( out_qp, svg->bf, aux );

    } else { // d/dv.
      FMF_SetCell( grad_s, ii );
      FMF_SetCellX1( vvg->bf, ii );

      // s^T G_c^T - transpose grad_s.
      fmf_mulATC( gst, grad_s, 1.0 );
      // s^T G_c^T Psi.
      bf_ract( aux, vvg->bf, gst );
      // Phi^T s^T G_c^T Psi.
      fmf_mulATB_nn( out_qp, svg->bf, aux );
    }

    fmf_sumLevelsMulF( out, out_qp, vvg->det->val );
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &aux );
  fmf_freeDestroy( &out_qp );
  fmf_freeDestroy( &gst );

  return( ret );
}
