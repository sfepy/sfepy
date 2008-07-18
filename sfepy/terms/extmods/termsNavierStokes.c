#include "termsNavierStokes.h"
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
int32 divgrad_build_gtg( FMField *out, FMField *gc )
{
  int32 iqp, ir, ic, dim, nEP, nQP, nCol;
  float64 *pout1, *pout2, *pout3, *pg1, *pg2, *pg3;

  nEP = gc->nCol;
  nQP = gc->nLev;
  nCol = out->nCol;
  dim = gc->nRow;

#ifdef DEBUG_FMF
  if ((out->nCol != (dim * nEP))
       || (out->nRow != (dim * nEP)) || (out->nLev != gc->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d), (%d %d %d)\n",
	    out->nLev, out->nRow, out->nCol,
	    gc->nLev, gc->nRow, gc->nCol );
  }
#endif

  fmf_fillC( out, 0.0 );
  switch (dim) {
  case 3:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      pg3 = pg2 + nEP;

      pout1 = FMF_PtrLevel( out, iqp );
      pout2 = pout1 + (nCol + 1) * nEP;
      pout3 = pout2 + (nCol + 1) * nEP;

      for (ir = 0; ir < nEP; ir++) {
	for (ic = 0; ic < nEP; ic++) {
	  pout1[ic] = pout2[ic] = pout3[ic]
	    = pg1[ir] * pg1[ic] + pg2[ir] * pg2[ic] + pg3[ir] * pg3[ic];
	}
	pout1 += nCol;
	pout2 += nCol;
	pout3 += nCol;
      }
    }
    break;
    
  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;

      pout1 = FMF_PtrLevel( out, iqp );
      pout2 = pout1 + (nCol + 1) * nEP;

      for (ir = 0; ir < nEP; ir++) {
	for (ic = 0; ic < nEP; ic++) {
	  pout1[ic] = pout2[ic]
	    = pg1[ir] * pg1[ic] + pg2[ir] * pg2[ic];
	}
	pout1 += nCol;
	pout2 += nCol;
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
int32 term_ns_asm_div_grad( FMField *out, FMField *state, int32 offset,
			  float64 viscosity, VolumeGeometry *vg,
			  int32 *conn, int32 nEl, int32 nEP,
			  int32 *elList, int32 elList_nRow,
			  int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *st = 0, *gtg = 0, *gu = 0, *gtgu = 0;
  FMField stv[1];

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

/*   output( "%d %d %d %d %d %d\n", offset, nEl, nEP, nQP, dim, elList_nRow ); */

  state->val = FMF_PtrFirst( state ) + offset;

/*     fmf_createAlloc( &gtg, 1, nQP, dim * nEP, dim * nEP ); */
  if (isDiff) {
    fmf_createAlloc( &gtg, 1, nQP, dim * nEP, dim * nEP );
  } else {
    fmf_createAlloc( &st, 1, 1, dim, nEP );
    fmf_createAlloc( &gu, 1, nQP, dim * dim, 1 );
    fmf_createAlloc( &gtgu, 1, nQP, dim * nEP, 1 );
    stv->nAlloc = -1;
    fmf_pretend( stv, 1, 1, nEP * dim, 1, st->val );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];
/*     output( "%d\n", iel ); */

    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( vg->det, iel );

    if (isDiff) {
      divgrad_build_gtg( gtg, vg->bfGM );
      fmf_sumLevelsMulF( out, gtg, vg->det->val );
      fmf_mulC( out, viscosity );
/*       fmf_print( gtg, stdout, 0 ); */
/*       sys_pause(); */
    } else {
      ele_extractNodalValuesDBD( st, state, conn + nEP * iel );

/*       divgrad_build_gtg( gtg, vg->bfGM ); ERR_CheckGo( ret ); */
/*       fmf_mulAB_n1( gtgu, gtg, stv ); */

      divgrad_act_g_m( gu, vg->bfGM, stv );
      divgrad_act_gt_m( gtgu, vg->bfGM, gu );
      fmf_sumLevelsMulF( out, gtgu, vg->det->val );
      fmf_mulC( out, viscosity );

      if (ii == 0) {
/* 	fmf_print( st, stdout, 0 ); */
/* 	fmf_print( vg->det, stdout, 0 ); */
/* 	fmf_print( out, stdout, 0 ); */
/* 	output( "%f\n", viscosity ); */
      }
/*       sys_pause(); */
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &gtg ); 
  } else {
    fmf_freeDestroy( &st ); 
    fmf_freeDestroy( &gu ); 
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
int32 term_ns_asm_convect( FMField *out, FMField *state, int32 offset,
			  FMField *bf, VolumeGeometry *vg,
			  int32 *conn, int32 nEl, int32 nEP,
			  int32 *elList, int32 elList_nRow,
			  int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *st = 0, *guf = 0, *ftguf = 0, *utg = 0, *ftutg = 0, *gufu = 0;
  FMField *fu = 0, *gu = 0, *ftgufu = 0;
  FMField stv[1], gum[1];

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

/*   output( "%d %d %d %d %d %d\n", offset, nEl, nEP, nQP, dim, elList_nRow ); */
  state->val = FMF_PtrFirst( state ) + offset;

  fmf_createAlloc( &st, 1, 1, dim, nEP );
  fmf_createAlloc( &fu, 1, nQP, dim, 1 );
  fmf_createAlloc( &gu, 1, nQP, dim * dim, 1 );

  stv->nAlloc = -1;
  fmf_pretend( stv, 1, 1, nEP * dim, 1, st->val );
  gum->nAlloc = -1;
  fmf_pretend( gum, 1, nQP, dim, dim, gu->val );
  
  if (isDiff) {
    fmf_createAlloc( &guf, 1, nQP, dim, dim * nEP );
    fmf_createAlloc( &ftguf, 1, nQP, dim * nEP, dim * nEP );
    fmf_createAlloc( &utg, 1, nQP, dim, nEP * dim );
    fmf_createAlloc( &ftutg, 1, nQP, dim * nEP, dim * nEP );
  } else {
    fmf_createAlloc( &gufu, 1, nQP, dim, 1 );
    fmf_createAlloc( &ftgufu, 1, nQP, dim * nEP, 1 );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

/*     output( "%d\n", iel ); */

    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( vg->det, iel );

    ele_extractNodalValuesDBD( st, state, conn + nEP * iel );
    bf_act( fu, bf, st );
    divgrad_act_g_m( gu, vg->bfGM, stv );

/*     fmf_print( vel, stdout, 0 ); */
/*     fmf_print( gu, stdout, 0 ); */
/*     sys_pause(); */

    if (isDiff) {
      bf_ract( guf, bf, gum );
      bf_actt( ftguf, bf, guf );

      convect_build_vtg( utg, vg->bfGM, fu );
      bf_actt( ftutg, bf, utg );

      fmf_addAB_nn( ftguf, ftguf, ftutg );
      fmf_sumLevelsMulF( out, ftguf, vg->det->val );
/*       fmf_print( out, stdout, 0 ); */
/*       sys_pause(); */
    } else {
      fmf_mulAB_nn( gufu, gum, fu );
      bf_actt( ftgufu, bf, gufu );
      fmf_sumLevelsMulF( out, ftgufu, vg->det->val );
    }
    ERR_CheckGo( ret );
  }


 end_label:
  fmf_freeDestroy( &st ); 
  fmf_freeDestroy( &fu ); 
  fmf_freeDestroy( &gu ); 

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
int32 dw_lin_convect( FMField *out,
		      FMField *stateB, int32 offsetB,
		      FMField *stateU, int32 offsetU,
		      FMField *bf, VolumeGeometry *vg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow,
		      int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *stu = 0, *stb = 0, *fb = 0, *btg = 0, *ftbtg = 0, *gu = 0, *gufb = 0;
  FMField *ftgufb = 0;
  FMField stuv[1], gum[1];

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

/*   output( "%d %d %d %d %d %d\n", offset, nEl, nEP, nQP, dim, elList_nRow ); */
  stateB->val = FMF_PtrFirst( stateB ) + offsetB;
  stateU->val = FMF_PtrFirst( stateU ) + offsetU;

  fmf_createAlloc( &stb, 1, 1, dim, nEP );
  fmf_createAlloc( &fb, 1, nQP, dim, 1 );

  if (isDiff == 1) {
    fmf_createAlloc( &btg, 1, nQP, dim, nEP * dim );
    fmf_createAlloc( &ftbtg, 1, nQP, dim * nEP, dim * nEP );
  } else {
    fmf_createAlloc( &stu, 1, 1, dim, nEP );
    fmf_createAlloc( &gu, 1, nQP, dim * dim, 1 );

    if (isDiff == 0) {
      fmf_createAlloc( &gufb, 1, nQP, dim, 1 );
      fmf_createAlloc( &ftgufb, 1, nQP, dim * nEP, 1 );
    }

    stuv->nAlloc = -1;
    fmf_pretend( stuv, 1, 1, nEP * dim, 1, stu->val );
    gum->nAlloc = -1;
    fmf_pretend( gum, 1, nQP, dim, dim, gu->val );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

/*     output( "%d\n", iel ); */

    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( vg->det, iel );

    ele_extractNodalValuesDBD( stb, stateB, conn + nEP * iel );
    bf_act( fb, bf, stb );

    if (isDiff == 1) {
      convect_build_vtg( btg, vg->bfGM, fb );
      bf_actt( ftbtg, bf, btg );

      fmf_sumLevelsMulF( out, ftbtg, vg->det->val );
/*       fmf_print( out, stdout, 0 ); */
/*       sys_pause(); */

    } else {
      ele_extractNodalValuesDBD( stu, stateU, conn + nEP * iel );
      divgrad_act_g_m( gu, vg->bfGM, stuv );
    
      if (isDiff == 0) {
	fmf_mulAB_nn( gufb, gum, fb );
	bf_actt( ftgufb, bf, gufb );
	fmf_sumLevelsMulF( out, ftgufb, vg->det->val );
      } else {
	fmf_mulAB_nn( out, gum, fb );
      }
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &stb );
  fmf_freeDestroy( &fb );
  fmf_freeDestroy( &btg ); 
  fmf_freeDestroy( &ftbtg ); 
  fmf_freeDestroy( &stu ); 
  fmf_freeDestroy( &gu ); 
  fmf_freeDestroy( &gufb ); 
  fmf_freeDestroy( &ftgufb ); 

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_div"
/*!
  @a state, @a conn are velocity-like.

  @par Revision history:
  - 14.12.2005, c
*/
int32 dw_div( FMField *out, FMField *state, int32 offset,
	      FMField *bf, VolumeGeometry *vg,
	      int32 *conn, int32 nEl, int32 nEP,
	      int32 *elList, int32 elList_nRow,
	      int32 isDiff )
{
  int32 ii, iel, nEPP, dim, nQP, ret = RET_OK;
  FMField *gu = 0, *ftgu = 0, *ftg = 0, *st = 0;
  FMField gcl[1], stv[1];

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;
  nEPP = bf->nCol;

  state->val = FMF_PtrFirst( state ) + offset;

  gcl->nAlloc = -1;
  fmf_pretend( gcl, 1, nQP, 1, nEP * dim, vg->bfGM->val0 );

  if (isDiff == 1) { 
    fmf_createAlloc( &ftg, 1, nQP, nEPP, dim * nEP );
  } else {
    fmf_createAlloc( &gu, 1, nQP, 1, 1 );
    fmf_createAlloc( &ftgu, 1, nQP, nEPP, 1 );
    fmf_createAlloc( &st, 1, 1, dim, nEP );
    stv->nAlloc = -1;
    fmf_pretend( stv, 1, 1, nEP * dim, 1, st->val );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( gcl, iel );
    FMF_SetCell( vg->det, iel );
      
    if (isDiff == 1) { 
      fmf_mulATB_nn( ftg, bf, gcl );
      fmf_sumLevelsMulF( out, ftg, vg->det->val );
    } else {
      ele_extractNodalValuesDBD( st, state, conn + nEP * iel );

      fmf_mulAB_n1( gu, gcl, stv );
      fmf_mulATB_nn( ftgu, bf, gu );
      fmf_sumLevelsMulF( out, ftgu, vg->det->val );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftg );
  } else {
    fmf_freeDestroy( &st );
    fmf_freeDestroy( &gu );
    fmf_freeDestroy( &ftgu );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_grad"
/*!
  @a state, @a conn are pressure-like.

  @par Revision history:
  - c: 15.12.2005, r: 31.03.2008
*/
int32 dw_grad( FMField *out, float64 coef, FMField *state, int32 offset,
	       FMField *bf, VolumeGeometry *vg,
	       int32 *conn, int32 nEl, int32 nEP,
	       int32 *elList, int32 elList_nRow,
	       int32 isDiff )
{
  int32 ii, iel, nEPU, dim, nQP, ret = RET_OK;
  FMField *fp = 0, *gtfp = 0, *gtf = 0, *st = 0;
  FMField gcl[1];

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;
  nEPU = vg->bfGM->nCol;

  state->val = FMF_PtrFirst( state ) + offset;

  gcl->nAlloc = -1;
  fmf_pretend( gcl, 1, nQP, 1, nEPU * dim, vg->bfGM->val0 );

  if (isDiff == 1) { 
    fmf_createAlloc( &gtf, 1, nQP, dim * nEPU, nEP );
  } else {
    fmf_createAlloc( &fp, 1, nQP, 1, 1 );
    fmf_createAlloc( &gtfp, 1, nQP, dim * nEPU, 1 );
    fmf_createAlloc( &st, 1, 1, nEP, 1 );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( gcl, iel );
    FMF_SetCell( vg->det, iel );
      
    if (isDiff == 1) { 
      fmf_mulATB_nn( gtf, gcl, bf );
      fmf_sumLevelsMulF( out, gtf, vg->det->val );
    } else {
      ele_extractNodalValuesNBN( st, state, conn + nEP * iel );

      fmf_mulAB_n1( fp, bf, st );
      fmf_mulATB_nn( gtfp, gcl, fp );
      fmf_sumLevelsMulF( out, gtfp, vg->det->val );
    }
    ERR_CheckGo( ret );
  }
  fmfc_mulC( out, coef );

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &gtf );
  } else {
    fmf_freeDestroy( &st );
    fmf_freeDestroy( &fp );
    fmf_freeDestroy( &gtfp );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_st_pspg_p"
/*!
  @par Revision history:
  - 09.01.2006, c
  - 10.01.2006
  - 31.07.2007
*/
int32 dw_st_pspg_p( FMField *out, FMField *state, int32 offset,
		    FMField *coef, VolumeGeometry *vg,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 *elList, int32 elList_nRow,
		    int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *gtg = 0, *gp = 0, *gtgp = 0, *st = 0;
  FMField stv[1];

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  state->val = FMF_PtrFirst( state ) + offset;
  if (isDiff == 1) { 
    fmf_createAlloc( &gtg, 1, nQP, nEP, nEP );
  } else {
    fmf_createAlloc( &gp, 1, nQP, dim, 1 );
    fmf_createAlloc( &gtgp, 1, nQP, nEP, 1 );
    fmf_createAlloc( &st, 1, 1, 1, nEP );
    stv->nAlloc = -1;
    fmf_pretend( stv, 1, 1, nEP, 1, st->val );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( vg->det, iel );
    FMF_SetCell( coef, iel );
      
    if (isDiff == 1) { 
      fmf_mulATB_nn( gtg, vg->bfGM, vg->bfGM );
      fmf_sumLevelsMulF( out, gtg, vg->det->val );
    } else {
      ele_extractNodalValuesDBD( st, state, conn + nEP * iel );

      fmf_mulAB_n1( gp, vg->bfGM, stv );
      fmf_mulATB_nn( gtgp, vg->bfGM, gp );
      fmf_sumLevelsMulF( out, gtgp, vg->det->val );
    }
    fmf_mulC( out, coef->val[0] );

    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &gtg );
  } else {
    fmf_freeDestroy( &st );
    fmf_freeDestroy( &gp );
    fmf_freeDestroy( &gtgp );
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
		    FMField *stateB, int32 offsetB,
		    FMField *stateU, int32 offsetU,
		    FMField *coef, FMField *bf_u,
		    VolumeGeometry *vg_p, VolumeGeometry *vg_u,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 *elList, int32 elList_nRow,
		    int32 isDiff )
{
  int32 ii, iel, dim, nQP, nEPP, ret = RET_OK;
  FMField *stu = 0, *stb = 0, *btg = 0, *gtbtg = 0, *btgu = 0;
  FMField *gtbtgu = 0, *fb = 0;
  FMField stuv[1];

  nQP = vg_u->bfGM->nLev;
  dim = vg_u->bfGM->nRow;
  nEPP = vg_p->bfGM->nCol;

  stateB->val = FMF_PtrFirst( stateB ) + offsetB;
  stateU->val = FMF_PtrFirst( stateU ) + offsetU;

  fmf_createAlloc( &stb, 1, 1, dim, nEP );
  fmf_createAlloc( &fb, 1, nQP, dim, 1 );
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

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg_u->bfGM, iel );
    FMF_SetCell( vg_p->bfGM, iel );
    FMF_SetCell( vg_u->det, iel );
    FMF_SetCell( coef, iel );

    ele_extractNodalValuesDBD( stb, stateB, conn + nEP * iel );
    bf_act( fb, bf_u, stb );
    convect_build_vtg( btg, vg_u->bfGM, fb );
      
    if (isDiff == 1) { 
      fmf_mulATB_nn( gtbtg, vg_p->bfGM, btg );
      fmf_sumLevelsMulF( out, gtbtg, vg_u->det->val );
    } else {
      ele_extractNodalValuesDBD( stu, stateU, conn + nEP * iel );

      fmf_mulAB_n1( btgu, btg, stuv );
      fmf_mulATB_nn( gtbtgu, vg_p->bfGM, btgu );
      fmf_sumLevelsMulF( out, gtbtgu, vg_u->det->val );
    }
    fmf_mulC( out, coef->val[0] );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &stb );
  fmf_freeDestroy( &fb );
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
		    FMField *stateB, int32 offsetB,
		    FMField *stateP, int32 offsetP,
		    FMField *coef, FMField *bf_u,
		    VolumeGeometry *vg_u, VolumeGeometry *vg_p,
		    int32 *conn_u, int32 nEl_u, int32 nEP_u,
		    int32 *conn_p, int32 nEl_p, int32 nEP_p,
		    int32 *elList, int32 elList_nRow,
		    int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *stp = 0, *stb = 0, *btg = 0, *gtbg = 0, *gp = 0;
  FMField *gtbgp = 0, *fb = 0;
  FMField stpv[1];

  nQP = vg_u->bfGM->nLev;
  dim = vg_u->bfGM->nRow;

  stateB->val = FMF_PtrFirst( stateB ) + offsetB;
  stateP->val = FMF_PtrFirst( stateP ) + offsetP;

  fmf_createAlloc( &stb, 1, 1, dim, nEP_u );
  fmf_createAlloc( &fb, 1, nQP, dim, 1 );
  fmf_createAlloc( &btg, 1, nQP, dim, nEP_u * dim );

  if (isDiff == 1) { 
    fmf_createAlloc( &gtbg, 1, nQP, dim * nEP_u, nEP_p );
  } else {
    fmf_createAlloc( &stp, 1, 1, 1, nEP_p );
    fmf_createAlloc( &gp, 1, nQP, dim, 1 );
    fmf_createAlloc( &gtbgp, 1, nQP, dim * nEP_u, 1 );

    stpv->nAlloc = -1;
    fmf_pretend( stpv, 1, 1, nEP_p, 1, stp->val );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg_u->bfGM, iel );
    FMF_SetCell( vg_p->bfGM, iel );
    FMF_SetCell( vg_u->det, iel );
    FMF_SetCell( coef, iel );

    ele_extractNodalValuesDBD( stb, stateB, conn_u + nEP_u * iel );
    bf_act( fb, bf_u, stb );
    convect_build_vtg( btg, vg_u->bfGM, fb );
      
    if (isDiff == 1) { 
      fmf_mulATB_nn( gtbg, btg, vg_p->bfGM );
      fmf_sumLevelsMulF( out, gtbg, vg_u->det->val );
    } else {
      ele_extractNodalValuesDBD( stp, stateP, conn_p + nEP_p * iel );

      fmf_mulAB_n1( gp, vg_p->bfGM, stpv );
      fmf_mulATB_nn( gtbgp, btg, gp );
      fmf_sumLevelsMulF( out, gtbgp, vg_u->det->val );
    }
    fmf_mulC( out, coef->val[0] );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &stb );
  fmf_freeDestroy( &fb );
  fmf_freeDestroy( &btg ); 
  if (isDiff) {
    fmf_freeDestroy( &gtbg );
  } else {
    fmf_freeDestroy( &stp );
    fmf_freeDestroy( &gp );
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
		    FMField *stateB, int32 offsetB,
		    FMField *stateU, int32 offsetU,
		    FMField *coef, FMField *bf, VolumeGeometry *vg,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 *elList, int32 elList_nRow,
		    int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *stu = 0, *stb = 0, *btg = 0, *gtbbtg = 0, *btgu = 0;
  FMField *gtbbtgu = 0, *fb = 0;
  FMField stuv[1];

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  stateB->val = FMF_PtrFirst( stateB ) + offsetB;
  stateU->val = FMF_PtrFirst( stateU ) + offsetU;

  fmf_createAlloc( &stb, 1, 1, dim, nEP );
  fmf_createAlloc( &fb, 1, nQP, dim, 1 );
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

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( vg->det, iel );
    FMF_SetCell( coef, iel );

    ele_extractNodalValuesDBD( stb, stateB, conn + nEP * iel );
    bf_act( fb, bf, stb );
    convect_build_vtg( btg, vg->bfGM, fb );
      
    if (isDiff == 1) { 
      fmf_mulATB_nn( gtbbtg, btg, btg );
      fmf_sumLevelsMulF( out, gtbbtg, vg->det->val );
    } else {
      ele_extractNodalValuesDBD( stu, stateU, conn + nEP * iel );

      fmf_mulAB_n1( btgu, btg, stuv );
      fmf_mulATB_nn( gtbbtgu, btg, btgu );
      fmf_sumLevelsMulF( out, gtbbtgu, vg->det->val );
    }
    fmf_mulC( out, coef->val[0] );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &stb );
  fmf_freeDestroy( &fb );
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
int32 dw_st_grad_div( FMField *out, FMField *state, int32 offset,
		      float64 gamma, VolumeGeometry *vg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow,
		      int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *gu = 0, *gtgu = 0, *gtg = 0, *st = 0;
  FMField gcl[1], stv[1];

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  state->val = FMF_PtrFirst( state ) + offset;

  gcl->nAlloc = -1;
  fmf_pretend( gcl, 1, nQP, 1, nEP * dim, vg->bfGM->val0 );

  if (isDiff == 1) { 
    fmf_createAlloc( &gtg, 1, nQP, dim * nEP, dim * nEP );
  } else {
    fmf_createAlloc( &gu, 1, nQP, 1, 1 );
    fmf_createAlloc( &gtgu, 1, nQP, dim * nEP, 1 );
    fmf_createAlloc( &st, 1, 1, dim, nEP );
    stv->nAlloc = -1;
    fmf_pretend( stv, 1, 1, nEP * dim, 1, st->val );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( gcl, iel );
    FMF_SetCell( vg->det, iel );
      
    if (isDiff == 1) { 
      fmf_mulATB_nn( gtg, gcl, gcl );
      fmf_sumLevelsMulF( out, gtg, vg->det->val );
    } else {
      ele_extractNodalValuesDBD( st, state, conn + nEP * iel );

      fmf_mulAB_n1( gu, gcl, stv );
      fmf_mulATB_nn( gtgu, gcl, gu );
      fmf_sumLevelsMulF( out, gtgu, vg->det->val );
    }
    ERR_CheckGo( ret );
  }
  fmfc_mulC( out, gamma );

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &gtg );
  } else {
    fmf_freeDestroy( &st );
    fmf_freeDestroy( &gu );
    fmf_freeDestroy( &gtgu );
  }

  return( ret );
}
