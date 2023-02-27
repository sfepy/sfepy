#include "terms_navier_stokes.h"
#include "terms.h"
#include "geommech.h"

/*
  Below the DBD order is used.
*/

#undef __FUNC__
#define __FUNC__ "divgrad_build_gtg"
/*!
  @par Revision history:
  - 26.10.2005, c
  - 24.05.2007
*/
int32 divgrad_build_gtg( FMField *out, FMField *gcv, FMField *gcs )
{
  int32 iqp, ir, ic, dim, nEPv, nEPs, nQP, nCol;
  float64 *pout1, *pout2, *pout3, *pg1v, *pg2v, *pg3v, *pg1s, *pg2s, *pg3s;

  nEPv = gcv->nCol;
  nEPs = gcs->nCol;
  nQP = gcv->nLev;
  nCol = out->nCol;
  dim = gcv->nRow;

#ifdef DEBUG_FMF
  if ((out->nCol != (dim * nEPs))
       || (out->nRow != (dim * nEPv)) || (out->nLev != gcv->nLev)
       || (out->nLev != gcs->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d), (%d %d %d), (%d %d %d)\n",
	    out->nLev, out->nRow, out->nCol,
	    gcv->nLev, gcv->nRow, gcv->nCol,
	    gcs->nLev, gcs->nRow, gcs->nCol );
  }
#endif

  fmf_fillC( out, 0.0 );
  switch (dim) {
  case 3:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1v = FMF_PtrLevel( gcv, iqp );
      pg2v = pg1v + nEPv;
      pg3v = pg2v + nEPv;

      pg1s = FMF_PtrLevel( gcs, iqp );
      pg2s = pg1s + nEPs;
      pg3s = pg2s + nEPs;

      pout1 = FMF_PtrLevel( out, iqp );
      pout2 = pout1 + (nCol + 1) * nEPv;
      pout3 = pout2 + (nCol + 1) * nEPv;

      for (ir = 0; ir < nEPv; ir++) {
        for (ic = 0; ic < nEPs; ic++) {
          pout1[ic] = pout2[ic] = pout3[ic]
            = pg1v[ir] * pg1s[ic] + pg2v[ir] * pg2s[ic] + pg3v[ir] * pg3s[ic];
        }
        pout1 += nCol;
        pout2 += nCol;
        pout3 += nCol;
      }
    }
    break;
    
  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1v = FMF_PtrLevel( gcv, iqp );
      pg2v = pg1v + nEPv;

      pg1s = FMF_PtrLevel( gcs, iqp );
      pg2s = pg1s + nEPs;

      pout1 = FMF_PtrLevel( out, iqp );
      pout2 = pout1 + (nCol + 1) * nEPv;

      for (ir = 0; ir < nEPv; ir++) {
        for (ic = 0; ic < nEPs; ic++) {
          pout1[ic] = pout2[ic]
            = pg1v[ir] * pg1s[ic] + pg2v[ir] * pg2s[ic];
        }
        pout1 += nCol;
        pout2 += nCol;
      }
    }
    break;

  case 1:
    for (iqp = 0; iqp < nQP; iqp++){
      pg1v = FMF_PtrLevel(gcv, iqp);
      pg1s = FMF_PtrLevel(gcs, iqp);
      pout1 = FMF_PtrLevel(out, iqp);

      for (ir = 0; ir < nEPv; ir++){
        for (ic = 0; ic < nEPs; ic++){
          pout1[ic] = pg1v[ir] * pg1s[ic];
        }
        pout1 += nCol;
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
#define __FUNC__ "divgrad_act_g_m"
/*!
  @par Revision history:
  - 26.10.2005, c
  - 28.11.2005
  - 14.12.2005
  - 24.05.2007
*/
int32 divgrad_act_g_m( FMField *out, FMField *gc, FMField *mtx )
{
  int32 iqp, ir, ic, ik, dim, nEP, nQP, nCol;
  float64 val1, val2, val3;
  float64 *pout, *pmtx, *pg1, *pg2, *pg3;

  nEP = gc->nCol;
  nQP = gc->nLev;
  nCol = mtx->nCol;
  dim = gc->nRow;

#ifdef DEBUG_FMF
  if ((out->nCol != mtx->nCol) || (out->nRow != (dim * dim))
       || (mtx->nRow != (dim * nEP)) || (out->nLev != gc->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d), (%d %d %d), (%d %d %d)\n",
	    out->nLev, out->nRow, out->nCol,
	    gc->nLev, gc->nRow, gc->nCol,
	    mtx->nLev, mtx->nRow, mtx->nCol );
  }
#endif

/*   output( "%d %d %d %d\n", nEP, nQP, nCol, dim ); */

  switch (dim) {
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
      for (ir = 0; ir < dim; ir++) {
	for (ic = 0; ic < nCol; ic++) {
	  val1 = val2 = val3 = 0.0;
	  for (ik = 0; ik < nEP; ik++) {
/* 	    output( "%d %d %d %d\n", iqp, ir, ic, ik ); */
	    val1 += pg1[ik] * pmtx[ic+nCol*ik];
	    val2 += pg2[ik] * pmtx[ic+nCol*ik];
	    val3 += pg3[ik] * pmtx[ic+nCol*ik];
	  }
	  pout[dim*nCol*ir+ic+0] = val1;
	  pout[dim*nCol*ir+ic+nCol] = val2;
	  pout[dim*nCol*ir+ic+2*nCol] = val3;
	}
	pmtx += nCol * nEP; 
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
      for (ir = 0; ir < dim; ir++) {
	for (ic = 0; ic < nCol; ic++) {
	  val1 = val2 = 0.0;
	  for (ik = 0; ik < nEP; ik++) {
/* 	    output( "%d %d %d %d\n", iqp, ir, ic, ik ); */
	    val1 += pg1[ik] * pmtx[ic+nCol*ik];
	    val2 += pg2[ik] * pmtx[ic+nCol*ik];
	  }
	  pout[dim*nCol*ir+ic+0] = val1;
	  pout[dim*nCol*ir+ic+nCol] = val2;
	}
	pmtx += nCol * nEP; 
      }
    }
    break;

  case 1:
    for (iqp = 0; iqp < nQP; iqp++){
      pg1 = FMF_PtrLevel(gc, iqp);
      pout = FMF_PtrLevel(out, iqp);

      if (mtx->nLev == nQP) {
        pmtx = FMF_PtrLevel(mtx, iqp);
      }
      else {
        pmtx = FMF_PtrCurrent(mtx);
      }

      for (ir = 0; ir < dim; ir++){
        for (ic = 0; ic < nCol; ic++){
          val1 = 0.0;
          for (ik = 0; ik < nEP; ik++){
            val1 += pg1[ik] * pmtx[ic + nCol*ik];
          }
          pout[dim*nCol*ir + ic + 0] = val1;
        }
        pmtx += nCol * nEP;
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
#define __FUNC__ "divgrad_act_gt_m"
/*!
  @par Revision history:
  - 26.10.2005, c
*/
int32 divgrad_act_gt_m( FMField *out, FMField *gc, FMField *mtx )
{
  int32 iqp, iep, ii, dim, nEP, nQP, nCol;
  float64 *pout1, *pout2, *pout3, *pmtx, *pg1, *pg2, *pg3;

  nEP = gc->nCol;
  nQP = gc->nLev;
  nCol = mtx->nCol;
  dim = gc->nRow;

#ifdef DEBUG_FMF
  if ((mtx->nRow != (dim * dim)) || (out->nRow != (dim * nEP))
      || (out->nCol != mtx->nCol) || (out->nLev != gc->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d), (%d %d %d), (%d %d %d)\n",
	    out->nLev, out->nRow, out->nCol,
	    gc->nLev, gc->nRow, gc->nCol,
	    mtx->nLev, mtx->nRow, mtx->nCol );
  }
#endif

  switch (dim) {
  case 3:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      pg3 = pg2 + nEP;
      
      pmtx = FMF_PtrLevel( mtx, iqp );
      for (iep = 0; iep < nEP; iep++) {
	pout1 = FMF_PtrLevel( out, iqp ) + nCol * iep;
	pout2 = pout1 + nCol * nEP;
	pout3 = pout2 + nCol * nEP;
	for (ii = 0; ii < nCol; ii++) {
	  pout1[ii]
	    = pg1[iep] * pmtx[0*nCol+ii]
	    + pg2[iep] * pmtx[1*nCol+ii]
	    + pg3[iep] * pmtx[2*nCol+ii];
	  pout2[ii]
	    = pg1[iep] * pmtx[3*nCol+ii]
	    + pg2[iep] * pmtx[4*nCol+ii]
	    + pg3[iep] * pmtx[5*nCol+ii];
	  pout3[ii]
	    = pg1[iep] * pmtx[6*nCol+ii]
	    + pg2[iep] * pmtx[7*nCol+ii]
	    + pg3[iep] * pmtx[8*nCol+ii];
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
	pout1 = FMF_PtrLevel( out, iqp ) + nCol * iep;
	pout2 = pout1 + nCol * nEP;
	for (ii = 0; ii < nCol; ii++) {
	  pout1[ii]
	    = pg1[iep] * pmtx[0*nCol+ii]
	    + pg2[iep] * pmtx[1*nCol+ii];
	  pout2[ii]
	    = pg1[iep] * pmtx[2*nCol+ii]
	    + pg2[iep] * pmtx[3*nCol+ii];
	}
      }
    }
    break;
  case 1:
    for (iqp = 0; iqp < nQP; iqp++){
      pg1 = FMF_PtrLevel(gc, iqp);
      pmtx = FMF_PtrLevel(mtx, iqp);

      for (iep = 0; iep < nEP; iep++){
        pout1 = FMF_PtrLevel(out, iqp) + nCol * iep;
        for (ii = 0; ii < nCol; ii++){
          pout1[ii] = pg1[iep] * pmtx[0 * nCol + ii];
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
#define __FUNC__ "divgrad_act_bg_m"
/*!
  @par Revision history:
  - 21.03.2006, c
*/
int32 divgrad_act_bg_m( FMField *out, FMField *gc, FMField *mtx )
{
  int32 iqp, ir, ic, ik, dim, nEP, nQP, nCol;
  float64 val1, val2, val3;
  float64 *pout, *pmtx, *pg1, *pg2, *pg3;

  nEP = gc->nCol;
  nQP = gc->nLev;
  nCol = mtx->nCol;
  dim = gc->nRow;

#ifdef DEBUG_FMF
  if ((out->nCol != mtx->nCol) || (out->nRow != (dim * dim))
       || (mtx->nRow != (dim * nEP)) || (out->nLev != gc->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d), (%d %d %d), (%d %d %d)\n",
	    out->nLev, out->nRow, out->nCol,
	    gc->nLev, gc->nRow, gc->nCol,
	    mtx->nLev, mtx->nRow, mtx->nCol );
  }
#endif

/*   output( "%d %d %d %d\n", nEP, nQP, nCol, dim ); */

  switch (dim) {
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
      for (ir = 0; ir < dim; ir++) {
	for (ic = 0; ic < nCol; ic++) {
	  val1 = val2 = val3 = 0.0;
	  for (ik = 0; ik < nEP; ik++) {
/* 	    output( "%d %d %d %d\n", iqp, ir, ic, ik ); */
	    val1 += pg1[ik] * pmtx[ic+nCol*ik];
	    val2 += pg2[ik] * pmtx[ic+nCol*ik];
	    val3 += pg3[ik] * pmtx[ic+nCol*ik];
	  }
	  pout[0*dim*nCol+nCol*ir+ic] = val1;
	  pout[1*dim*nCol+nCol*ir+ic] = val2;
	  pout[2*dim*nCol+nCol*ir+ic] = val3;
	}
	pmtx += nCol * nEP; 
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
#define __FUNC__ "convect_build_vtbg"
/*!
  U^T \bar{G}

  @par Revision history:
  - 26.10.2005, c
*/
int32 convect_build_vtbg( FMField *out, FMField *gc, FMField *fv )
{
  int32 iqp, iep, ii, dim, nEP, nQP;
  float64 *pout1, *pout2, *pout3, *pfv, *pg1, *pg2, *pg3;

  nEP = gc->nCol;
  nQP = gc->nLev;
  dim = gc->nRow;

#ifdef DEBUG_FMF
  if ((fv->nRow != dim) || (fv->nCol != 1)|| (out->nRow != dim) 
      || (out->nCol != (dim * nEP)) || (out->nLev != gc->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d), (%d %d %d), (%d %d %d)\n",
	    out->nLev, out->nRow, out->nCol,
	    gc->nLev, gc->nRow, gc->nCol,
	    fv->nLev, fv->nRow, fv->nCol );
  }
#endif

  switch (dim) {
  case 3:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      pg3 = pg2 + nEP;

      pout1 = FMF_PtrLevel( out, iqp );
      pout2 = pout1 + dim * nEP;
      pout3 = pout2 + dim * nEP;

      pfv = FMF_PtrLevel( fv, iqp );
      for (ii = 0; ii < dim; ii++) {
	for (iep = 0; iep < nEP; iep++) {
	  pout1[iep] = pg1[iep] * pfv[ii];
	  pout2[iep] = pg2[iep] * pfv[ii];
	  pout3[iep] = pg3[iep] * pfv[ii];
	}
	pout1 += nEP;
	pout2 += nEP;
	pout3 += nEP;
      }
    }
    break;
  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;

      pout1 = FMF_PtrLevel( out, iqp );
      pout2 = pout1 + dim * nEP;

      pfv = FMF_PtrLevel( fv, iqp );
      for (ii = 0; ii < dim; ii++) {
	for (iep = 0; iep < nEP; iep++) {
	  pout1[iep] = pg1[iep] * pfv[ii];
	  pout2[iep] = pg2[iep] * pfv[ii];
	}
	pout1 += nEP;
	pout2 += nEP;
      }
    }
    break;
  case 1:
    for (iqp = 0; iqp < nQP; iqp++){
      pg1 = FMF_PtrLevel(gc, iqp);
      pout1 = FMF_PtrLevel(out, iqp);
      pfv = FMF_PtrLevel(fv, iqp);

      for (ii = 0; ii < dim; ii++){
        for (iep = 0; iep < nEP; iep++){
          pout1[iep] = pg1[iep] * pfv[ii];
        }
        pout1 += nEP;
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
#define __FUNC__ "convect_build_vtg"
/*!
  U^T G

  @par Revision history:
  - 21.03.2006, c
*/
int32 convect_build_vtg( FMField *out, FMField *gc, FMField *fv )
{
  int32 iqp, iep, dim, nEP, nQP;
  float64 *pout1, *pout2, *pout3, *pfv, *pg1, *pg2, *pg3;

  nEP = gc->nCol;
  nQP = gc->nLev;
  dim = gc->nRow;

#ifdef DEBUG_FMF
  if ((fv->nRow != dim) || (fv->nCol != 1)|| (out->nRow != dim) 
      || (out->nCol != (dim * nEP)) || (out->nLev != gc->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d), (%d %d %d), (%d %d %d)\n",
	    out->nLev, out->nRow, out->nCol,
	    gc->nLev, gc->nRow, gc->nCol,
	    fv->nLev, fv->nRow, fv->nCol );
  }
#endif

  switch (dim) {
  case 3:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      pg3 = pg2 + nEP;

      pout1 = FMF_PtrLevel( out, iqp );
      pout2 = pout1 + dim * nEP + nEP;
      pout3 = pout2 + dim * nEP + nEP;

      pfv = FMF_PtrLevel( fv, iqp );
      for (iep = 0; iep < nEP; iep++) {
	pout1[iep] = pout2[iep] = pout3[iep]
	  = pg1[iep] * pfv[0] + pg2[iep] * pfv[1] + pg3[iep] * pfv[2];
      }
    }
    break;
  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;

      pout1 = FMF_PtrLevel( out, iqp );
      pout2 = pout1 + dim * nEP + nEP;

      pfv = FMF_PtrLevel( fv, iqp );
      for (iep = 0; iep < nEP; iep++) {
	pout1[iep] = pout2[iep]
	  = pg1[iep] * pfv[0] + pg2[iep] * pfv[1];
      }
    }
    break;
  case 1:
    for (iqp = 0; iqp < nQP; iqp++){
      pg1 = FMF_PtrLevel(gc, iqp);
      pout1 = FMF_PtrLevel(out, iqp);
      pfv = FMF_PtrLevel(fv, iqp);

      for (iep = 0; iep < nEP; iep++){
        pout1[iep] = pg1[iep] * pfv[0];
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
#define __FUNC__ "term_ns_asm_div_grad"
/*!
  @par Revision history:
  - 25.10.2005, c
  - 26.10.2005
  - 09.12.2005
  - 14.12.2005
*/
int32 term_ns_asm_div_grad( FMField *out, FMField *grad,
			    FMField *viscosity, Mapping *vgv, Mapping *vgs,
			    int32 isDiff )
{
  int32 ii, dim, nQP, nEPv, nEPs, ret = RET_OK;
  FMField *gtg = 0, *gtgu = 0;

  nQP = vgv->bfGM->nLev;
  nEPv = vgv->bfGM->nCol;
  nEPs = vgs->bfGM->nCol;
  dim = vgv->bfGM->nRow;

  if (isDiff) {
    fmf_createAlloc( &gtg, 1, nQP, dim * nEPv, dim * nEPs );
  } else {
    fmf_createAlloc( &gtgu, 1, nQP, dim * nEPv, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCellX1( viscosity, ii );
    FMF_SetCell( vgv->bfGM, ii );
    FMF_SetCell( vgv->det, ii );

    if (isDiff) {
      FMF_SetCell( vgs->bfGM, ii );

      divgrad_build_gtg( gtg, vgv->bfGM, vgs->bfGM );
      fmf_mul( gtg, viscosity->val );
      fmf_sumLevelsMulF( out, gtg, vgv->det->val );
    } else {
      FMF_SetCell( grad, ii );
      divgrad_act_gt_m( gtgu, vgv->bfGM, grad );
      fmf_mul( gtgu, viscosity->val );
      fmf_sumLevelsMulF( out, gtgu, vgv->det->val );
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
#define __FUNC__ "term_ns_asm_convect"
/*!
  @par Revision history:
  - 20.12.2005, c
  - 30.07.2007
*/
int32 term_ns_asm_convect( FMField *out, FMField *grad, FMField *state,
                           Mapping *vg, int32 isDiff )
{
  int32 ii, dim, nQP, nEP, ret = RET_OK;
  FMField *guf = 0, *ftguf = 0, *utg = 0, *ftutg = 0, *gufu = 0;
  FMField *ftgufu = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;
  dim = vg->bfGM->nRow;

  if (isDiff) {
    fmf_createAlloc( &guf, 1, nQP, dim, dim * nEP );
    fmf_createAlloc( &ftguf, 1, nQP, dim * nEP, dim * nEP );
    fmf_createAlloc( &utg, 1, nQP, dim, nEP * dim );
    fmf_createAlloc( &ftutg, 1, nQP, dim * nEP, dim * nEP );
  } else {
    fmf_createAlloc( &gufu, 1, nQP, dim, 1 );
    fmf_createAlloc( &ftgufu, 1, nQP, dim * nEP, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( state, ii );
    FMF_SetCell( grad, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCellX1( vg->bf, ii );

    if (isDiff) {
      bf_ract( guf, vg->bf, grad );
      bf_actt( ftguf, vg->bf, guf );

      convect_build_vtg( utg, vg->bfGM, state );
      bf_actt( ftutg, vg->bf, utg );

      fmf_addAB_nn( ftguf, ftguf, ftutg );
      fmf_sumLevelsMulF( out, ftguf, vg->det->val );
    } else {
      fmf_mulAB_nn( gufu, grad, state );
      bf_actt( ftgufu, vg->bf, gufu );
      fmf_sumLevelsMulF( out, ftgufu, vg->det->val );
    }
    ERR_CheckGo( ret );
  }


 end_label:
  if (isDiff) {
    fmf_freeDestroy( &guf );
    fmf_freeDestroy( &ftguf );
    fmf_freeDestroy( &utg );
    fmf_freeDestroy( &ftutg );
  } else {
    fmf_freeDestroy( &gufu );
    fmf_freeDestroy( &ftgufu );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "term_dw_lin_convect"
/*!
  isDiff == 2 for qp mode.

  @par Revision history:
  - 25.07.2007, c
  - 30.07.2007
*/
int32 dw_lin_convect( FMField *out, FMField *grad, FMField *stateB,
		      Mapping *vg, int32 isDiff )
{
  int32 ii, dim, nQP, nEP, ret = RET_OK;
  FMField *btg = 0, *ftbtg = 0, *gufb = 0, *ftgufb = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;
  dim = vg->bfGM->nRow;

  if (isDiff == 1) {
    fmf_createAlloc( &btg, 1, nQP, dim, nEP * dim );
    fmf_createAlloc( &ftbtg, 1, nQP, dim * nEP, dim * nEP );
  } else {
    if (isDiff == 0) {
      fmf_createAlloc( &gufb, 1, nQP, dim, 1 );
      fmf_createAlloc( &ftgufb, 1, nQP, dim * nEP, 1 );
    }
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCellX1( stateB, ii );
    FMF_SetCell( grad, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCellX1( vg->bf, ii );

    if (isDiff == 1) {
      convect_build_vtg( btg, vg->bfGM, stateB );
      bf_actt( ftbtg, vg->bf, btg );

      fmf_sumLevelsMulF( out, ftbtg, vg->det->val );

    } else {
      if (isDiff == 0) {
	fmf_mulAB_nn( gufb, grad, stateB );
	bf_actt( ftgufb, vg->bf, gufb );
	fmf_sumLevelsMulF( out, ftgufb, vg->det->val );
      } else {
	fmf_mulAB_nn( out, grad, stateB );
      }
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &btg );
  fmf_freeDestroy( &ftbtg );
  fmf_freeDestroy( &gufb );
  fmf_freeDestroy( &ftgufb );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_div"
/*!
  @par Revision history:
  - 14.12.2005, c
*/
int32 dw_div( FMField *out, FMField *coef, FMField *div,
	      Mapping *svg, Mapping *vvg, int32 isDiff )
{
  int32 ii, nEPP, dim, nQP, nEP, ret = RET_OK;
  FMField *ftgu = 0, *ftg = 0;
  FMField gcl[1];

  nQP = vvg->bfGM->nLev;
  nEP = vvg->bfGM->nCol;
  dim = vvg->bfGM->nRow;
  nEPP = svg->bf->nCol;

  gcl->nAlloc = -1;
  fmf_pretend( gcl, vvg->bfGM->nCell, nQP, 1, nEP * dim, vvg->bfGM->val0 );

  if (isDiff == 1) {
    fmf_createAlloc( &ftg, 1, nQP, nEPP, dim * nEP );
  } else {
    fmf_createAlloc( &ftgu, 1, nQP, nEPP, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( gcl, ii );
    FMF_SetCell( vvg->det, ii );
    FMF_SetCellX1( coef, ii );
    FMF_SetCellX1( svg->bf, ii );

    if (isDiff == 1) {
      fmf_mulATB_nn( ftg, svg->bf, gcl );
      fmf_mulAF( ftg, ftg, coef->val );
      fmf_sumLevelsMulF( out, ftg, vvg->det->val );
    } else {
      FMF_SetCell( div, ii );
      fmf_mulATB_nn( ftgu, svg->bf, div );
      fmf_mulAF( ftgu, ftgu, coef->val );
      fmf_sumLevelsMulF( out, ftgu, vvg->det->val );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftg );
  } else {
    fmf_freeDestroy( &ftgu );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_grad"
/*!
  @par Revision history:
  - c: 15.12.2005, r: 31.03.2008
*/
int32 dw_grad( FMField *out, FMField *coef, FMField *state,
	       Mapping *svg, Mapping *vvg, int32 isDiff )
{
  int32 ii, nEPU, dim, nQP, nEP, ret = RET_OK;
  FMField *gtfp = 0, *gtf = 0;
  FMField gcl[1];

  nQP = vvg->bfGM->nLev;
  dim = vvg->bfGM->nRow;
  nEPU = vvg->bfGM->nCol;
  nEP = svg->bf->nCol;

  gcl->nAlloc = -1;
  fmf_pretend( gcl, vvg->bfGM->nCell, nQP, 1, nEPU * dim, vvg->bfGM->val0 );

  if (isDiff == 1) {
    fmf_createAlloc( &gtf, 1, nQP, dim * nEPU, nEP );
  } else {
    fmf_createAlloc( &gtfp, 1, nQP, dim * nEPU, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( gcl, ii );
    FMF_SetCell( vvg->det, ii );
    FMF_SetCellX1( coef, ii );

    if (isDiff == 1) {
      FMF_SetCellX1( svg->bf, ii );
      fmf_mulATB_nn( gtf, gcl, svg->bf );
      fmf_mulAF( gtf, gtf, coef->val );
      fmf_sumLevelsMulF( out, gtf, vvg->det->val );
    } else {
      FMF_SetCell( state, ii );
      fmf_mulATB_nn( gtfp, gcl, state );
      fmf_mulAF( gtfp, gtfp, coef->val );
      fmf_sumLevelsMulF( out, gtfp, vvg->det->val );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &gtf );
  } else {
    fmf_freeDestroy( &gtfp );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_st_pspg_c"
/*!
  @par Revision history:
  - 31.07.2007, c
*/
int32 dw_st_pspg_c( FMField *out,
		    FMField *stateB, FMField *stateU,
		    FMField *coef,
		    Mapping *vg_p, Mapping *vg_u,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 isDiff )
{
  int32 ii, dim, nQP, nEPP, ret = RET_OK;
  FMField *stu = 0, *btg = 0, *gtbtg = 0, *btgu = 0;
  FMField *gtbtgu = 0;
  FMField stuv[1];

  nQP = vg_u->bfGM->nLev;
  dim = vg_u->bfGM->nRow;
  nEPP = vg_p->bfGM->nCol;

  stateU->val = FMF_PtrFirst( stateU );

  fmf_createAlloc( &btg, 1, nQP, dim, nEP * dim );

  if (isDiff == 1) {
    fmf_createAlloc( &gtbtg, 1, nQP, nEPP, dim * nEP );
  } else {
    fmf_createAlloc( &stu, 1, 1, dim, nEP );
    fmf_createAlloc( &btgu, 1, nQP, dim, 1 );
    fmf_createAlloc( &gtbtgu, 1, nQP, nEPP, 1 );

    stuv->nAlloc = -1;
    fmf_pretend( stuv, 1, 1, nEP * dim, 1, stu->val );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( stateB, ii );
    FMF_SetCell( vg_u->bfGM, ii );
    FMF_SetCell( vg_p->bfGM, ii );
    FMF_SetCell( vg_u->det, ii );
    FMF_SetCellX1( coef, ii );

    convect_build_vtg( btg, vg_u->bfGM, stateB );

    if (isDiff == 1) {
      fmf_mulATB_nn( gtbtg, vg_p->bfGM, btg );
      fmf_mul( gtbtg, coef->val );
      fmf_sumLevelsMulF( out, gtbtg, vg_u->det->val );
    } else {
      ele_extractNodalValuesDBD( stu, stateU, conn + nEP * ii );

      fmf_mulAB_n1( btgu, btg, stuv );
      fmf_mulATB_nn( gtbtgu, vg_p->bfGM, btgu );
      fmf_mul( gtbtgu, coef->val );
      fmf_sumLevelsMulF( out, gtbtgu, vg_u->det->val );
    }

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &btg );
  if (isDiff) {
    fmf_freeDestroy( &gtbtg );
  } else {
    fmf_freeDestroy( &stu );
    fmf_freeDestroy( &btgu );
    fmf_freeDestroy( &gtbtgu );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_st_supg_p"
/*!
  @par Revision history:
  - 31.07.2007, c
*/
int32 dw_st_supg_p( FMField *out,
		    FMField *stateB, FMField *gradP,
		    FMField *coef,
		    Mapping *vg_u, Mapping *vg_p,
		    int32 isDiff )
{
  int32 ii, dim, nQP, nEP_u, nEP_p, ret = RET_OK;
  FMField *btg = 0, *gtbg = 0, *gtbgp = 0;

  nQP = vg_u->bfGM->nLev;
  nEP_u = vg_u->bfGM->nCol;
  dim = vg_u->bfGM->nRow;
  nEP_p = vg_p->bfGM->nCol;

  fmf_createAlloc( &btg, 1, nQP, dim, nEP_u * dim );

  if (isDiff == 1) { 
    fmf_createAlloc( &gtbg, 1, nQP, dim * nEP_u, nEP_p );
  } else {
    fmf_createAlloc( &gtbgp, 1, nQP, dim * nEP_u, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( vg_u->bfGM, ii );
    FMF_SetCell( vg_p->bfGM, ii );
    FMF_SetCell( vg_u->det, ii );
    FMF_SetCellX1( coef, ii );
    FMF_SetCell( stateB, ii );

    convect_build_vtg( btg, vg_u->bfGM, stateB );

    if (isDiff == 1) {
      fmf_mulATB_nn( gtbg, btg, vg_p->bfGM );
      fmf_mul( gtbg, coef->val );
      fmf_sumLevelsMulF( out, gtbg, vg_u->det->val );
    } else {
      FMF_SetCell( gradP, ii );
      fmf_mulATB_nn( gtbgp, btg, gradP );
      fmf_mul( gtbgp, coef->val );
      fmf_sumLevelsMulF( out, gtbgp, vg_u->det->val );
    }

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &btg );
  if (isDiff) {
    fmf_freeDestroy( &gtbg );
  } else {
    fmf_freeDestroy( &gtbgp );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_st_supg_c"
/*!
  @par Revision history:
  - 31.07.2007, c
*/
int32 dw_st_supg_c( FMField *out,
		    FMField *stateB, FMField *stateU,
		    FMField *coef, Mapping *vg,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 isDiff )
{
  int32 ii, dim, nQP, ret = RET_OK;
  FMField *stu = 0, *btg = 0, *gtbbtg = 0, *btgu = 0, *gtbbtgu = 0;
  FMField stuv[1];

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  stateU->val = FMF_PtrFirst( stateU );

  fmf_createAlloc( &btg, 1, nQP, dim, nEP * dim );

  if (isDiff == 1) {
    fmf_createAlloc( &gtbbtg, 1, nQP, dim * nEP, dim * nEP );
  } else {
    fmf_createAlloc( &stu, 1, 1, dim, nEP );
    fmf_createAlloc( &btgu, 1, nQP, dim, 1 );
    fmf_createAlloc( &gtbbtgu, 1, nQP, dim * nEP, 1 );

    stuv->nAlloc = -1;
    fmf_pretend( stuv, 1, 1, nEP * dim, 1, stu->val );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCellX1( coef, ii );
    FMF_SetCell( stateB, ii );

    convect_build_vtg( btg, vg->bfGM, stateB );

    if (isDiff == 1) {
      fmf_mulATB_nn( gtbbtg, btg, btg );
      fmf_mul( gtbbtg, coef->val );
      fmf_sumLevelsMulF( out, gtbbtg, vg->det->val );
    } else {
      ele_extractNodalValuesDBD( stu, stateU, conn + nEP * ii );

      fmf_mulAB_n1( btgu, btg, stuv );
      fmf_mulATB_nn( gtbbtgu, btg, btgu );
      fmf_mul( gtbbtgu, coef->val );
      fmf_sumLevelsMulF( out, gtbbtgu, vg->det->val );
    }

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &btg );
  if (isDiff) {
    fmf_freeDestroy( &gtbbtg );
  } else {
    fmf_freeDestroy( &stu );
    fmf_freeDestroy( &btgu );
    fmf_freeDestroy( &gtbbtgu );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_st_grad_div"
/*!
  @par Revision history:
  - 26.07.2007, c
*/
int32 dw_st_grad_div( FMField *out, FMField *div,
		      FMField *coef, Mapping *vg,
		      int32 isDiff )
{
  int32 ii, dim, nQP, nEP, ret = RET_OK;
  FMField *gtgu = 0, *gtg = 0;
  FMField gcl[1];

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;
  dim = vg->bfGM->nRow;

  gcl->nAlloc = -1;
  fmf_pretend( gcl, vg->bfGM->nCell, nQP, 1, nEP * dim, vg->bfGM->val0 );

  if (isDiff == 1) {
    fmf_createAlloc( &gtg, 1, nQP, dim * nEP, dim * nEP );
  } else {
    fmf_createAlloc( &gtgu, 1, nQP, dim * nEP, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCellX1( coef, ii );
    FMF_SetCell( gcl, ii );
    FMF_SetCell( vg->det, ii );

    if (isDiff == 1) {
      fmf_mulATB_nn( gtg, gcl, gcl );
      fmf_mul( gtg, coef->val );
      fmf_sumLevelsMulF( out, gtg, vg->det->val );
    } else {
      FMF_SetCell( div, ii );
      fmf_mulATB_nn( gtgu, gcl, div );
      fmf_mul( gtgu, coef->val );
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
