#include "form_sdcc.h"

#undef __FUNC__
#define __FUNC__ "form_sdcc_strainCauchy_VS"
/*!
  Cauchy strain tensor stored as a vector, symmetric storage,
  doubled non-diagonal entries.
  
  @par Revision history:
  - 30.04.2001, c
  - 17.03.2003, adopted from rcfem2
  - 07.08.2006
*/
int32 form_sdcc_strainCauchy_VS( FMField *strain, FMField *dg )
{
  int32 iqp, nQP;
  float64 *pstrain, *pdg;

#ifdef DEBUG_FMF
  if ((strain->nLev != dg->nLev) || (strain->nCol != 1)
      || (dg->nCol != dg->nCol)
      || (strain->nRow != ((dg->nRow + 1) * dg->nRow / 2))) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d) <- (%d %d %d)\n",
	    strain->nLev, strain->nRow, strain->nCol,
	    dg->nLev, dg->nRow, dg->nCol );
  }
#endif

  nQP = dg->nLev;

  switch (dg->nRow) {
  case 1:
    for (iqp = 0; iqp < nQP; iqp++) {
      pstrain = FMF_PtrLevel( strain, iqp );
      pdg = FMF_PtrLevel( dg, iqp );
      pstrain[0] = pdg[0];
    }
    break;
  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pstrain = FMF_PtrLevel( strain, iqp );
      pdg = FMF_PtrLevel( dg, iqp );
      pstrain[0] = pdg[0];
      pstrain[1] = pdg[3];
      pstrain[2] = pdg[1] + pdg[2];
    }
    break;
  case 3:
    for (iqp = 0; iqp < nQP; iqp++) {
      pstrain = FMF_PtrLevel( strain, iqp );
      pdg = FMF_PtrLevel( dg, iqp );
      pstrain[0] = pdg[0];
      pstrain[1] = pdg[4];
      pstrain[2] = pdg[8];
      pstrain[3] = pdg[1] + pdg[3];
      pstrain[4] = pdg[2] + pdg[6];
      pstrain[5] = pdg[5] + pdg[7];
    }
    break;
  default:
    errput( ErrHead "ERR_Switch\n" );
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "form_sdcc_actOpGT_VS3"
/*!
  @fn int32 form_sdcc_actOpGT_VS3( FMField *diff, FMField *vec, FMField *gc )

  @f$(G_{SD})^T@f$ operation (i.e. divergence of columns of a tensor)
  used on a symmetric tensor (symmetric vector storage - e.g. stress vector).

  @par Revision history:
  - 30.04.2001, c
  - 03.10.2001
  - 06.06.2002
  - 17.03.2003, adopted from rcfem2
*/
int32 form_sdcc_actOpGT_VS3( FMField *diff, FMField *gc, FMField *vec )
{
  int32 iqp, iep, nEP, nQP;
  float64 *pdiff1, *pdiff2, *pdiff3, *pvec, *pg1, *pg2, *pg3;

  nEP = gc->nCol;
  nQP = gc->nLev;

  switch (gc->nRow) {
  case 3:
    for (iqp = 0; iqp < nQP; iqp++) {
      pdiff1 = FMF_PtrLevel( diff, iqp );
      pdiff2 = pdiff1 + nEP;
      pdiff3 = pdiff2 + nEP;
      pvec = FMF_PtrLevel( vec, iqp );
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      pg3 = pg2 + nEP;
      for (iep = 0; iep < nEP; iep++) {
	pdiff1[iep]
	  = pg1[iep] * pvec[0]
	  + pg2[iep] * pvec[3]
	  + pg3[iep] * pvec[4];
	pdiff2[iep]
	  = pg1[iep] * pvec[3]
	  + pg2[iep] * pvec[1]
	  + pg3[iep] * pvec[5];
	pdiff3[iep]
	  = pg1[iep] * pvec[4]
	  + pg2[iep] * pvec[5]
	  + pg3[iep] * pvec[2];
      }
    }
    break;
  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pdiff1 = FMF_PtrLevel( diff, iqp );
      pdiff2 = pdiff1 + nEP;
      pvec = FMF_PtrLevel( vec, iqp );
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      for (iep = 0; iep < nEP; iep++) {
	pdiff1[iep]
	  = pg1[iep] * pvec[0]
	  + pg2[iep] * pvec[2];
	pdiff2[iep]
	  = pg1[iep] * pvec[2]
	  + pg2[iep] * pvec[1];
      }
    }
    break;
  case 1:
    for (iqp = 0; iqp < nQP; iqp++){
      pdiff1 = FMF_PtrLevel(diff, iqp);
      pvec = FMF_PtrLevel(vec, iqp);
      pg1 = FMF_PtrLevel(gc, iqp);

      for (iep = 0; iep < nEP; iep++){
        pdiff1[iep] = pg1[iep] * pvec[0];
      }
    }
    break;
  default:
    errput( ErrHead "ERR_Switch\n" );
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "form_sdcc_actOpGT_M3"
/*!
  @fn int32 form_sdcc_actOpGT_M3( FMField *diff, FMField *vec, FMField *gc )

  @f$(G_{SD})^T@f$ operation (i.e. divergence of columns of a tensor)
  used on a matrix (e.g. tangent modulus) from left.

  @par Revision history:
  - 03.10.2001, c
  - 06.06.2002
  - 17.03.2003, adopted from rcfem2
*/
int32 form_sdcc_actOpGT_M3( FMField *diff, FMField *gc, FMField *mtx )
{
  int32 iqp, iep, ii, nEP, nQP, nCol;
  float64 *pdiff1, *pdiff2, *pdiff3, *pvec, *pg1, *pg2, *pg3;

  nEP = gc->nCol;
  nQP = gc->nLev;
  nCol = mtx->nCol;

  switch (gc->nRow) {
  case 3:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      pg3 = pg2 + nEP;
      
      pvec = FMF_PtrLevel( mtx, iqp );
      for (iep = 0; iep < nEP; iep++) {
	pdiff1 = FMF_PtrLevel( diff, iqp ) + nCol * iep;
	pdiff2 = pdiff1 + nCol * nEP;
	pdiff3 = pdiff2 + nCol * nEP;
	for (ii = 0; ii < nCol; ii++) {
	  pdiff1[ii]
	    = pg1[iep] * pvec[0*nCol+ii]
	    + pg2[iep] * pvec[3*nCol+ii]
	    + pg3[iep] * pvec[4*nCol+ii];
	  pdiff2[ii]
	    = pg1[iep] * pvec[3*nCol+ii]
	    + pg2[iep] * pvec[1*nCol+ii]
	    + pg3[iep] * pvec[5*nCol+ii];
	  pdiff3[ii]
	    = pg1[iep] * pvec[4*nCol+ii]
	    + pg2[iep] * pvec[5*nCol+ii]
	    + pg3[iep] * pvec[2*nCol+ii];
	}
      }
    }
    break;
  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      
      pvec = FMF_PtrLevel( mtx, iqp );
      for (iep = 0; iep < nEP; iep++) {
	pdiff1 = FMF_PtrLevel( diff, iqp ) + nCol * iep;
	pdiff2 = pdiff1 + nCol * nEP;
	for (ii = 0; ii < nCol; ii++) {
	  pdiff1[ii]
	    = pg1[iep] * pvec[0*nCol+ii]
	    + pg2[iep] * pvec[2*nCol+ii];
	  pdiff2[ii]
	    = pg1[iep] * pvec[2*nCol+ii]
	    + pg2[iep] * pvec[1*nCol+ii];
	}
      }
    }
    break;
  case 1:
    for (iqp = 0; iqp < nQP; iqp++){
      pg1 = FMF_PtrLevel(gc, iqp);
      pvec = FMF_PtrLevel(mtx, iqp);

      for (iep = 0; iep < nEP; iep++){
        pdiff1 = FMF_PtrLevel(diff, iqp) + nCol * iep;
        for (ii = 0; ii < nCol; ii++){
          pdiff1[ii] = pg1[iep] * pvec[0 * nCol + ii];
        }
      }
    }
    break;
  default:
    errput( ErrHead "ERR_Switch\n" );
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "form_sdcc_actOpG_RM3"
/*!
  @fn int32 form_sdcc_actOpG_RM3( FMField *diff, FMField *vec, FMField *gc )

  @f$(G_{SD})@f$ operation (i.e. divergence of columns of a tensor)
  used on a matrix (e.g. tangent modulus) from right.

  @par Revision history:
  - 03.10.2001, c
  - 06.06.2002
  - 17.03.2003, adopted from rcfem2
*/
int32 form_sdcc_actOpG_RM3( FMField *diff, FMField *mtx, FMField *gc )
{
  int32 iqp, iep, ii, nEP, nQP, nRow;
  float64 *pdiff1, *pdiff2, *pdiff3, *pvec, *pg1, *pg2, *pg3;

  nEP = gc->nCol;
  nQP = gc->nLev;
  nRow = mtx->nRow;

  switch (gc->nRow) {
  case 3:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      pg3 = pg2 + nEP;
      
      for (ii = 0; ii < nRow; ii++) {
	pvec = FMF_PtrLevel( mtx, iqp ) + mtx->nCol * ii;
	pdiff1 = FMF_PtrLevel( diff, iqp ) + diff->nCol * ii;
	pdiff2 = pdiff1 + nEP;
	pdiff3 = pdiff2 + nEP;
	for (iep = 0; iep < nEP; iep++) {
	  pdiff1[iep]
	    = pg1[iep] * pvec[0]
	    + pg2[iep] * pvec[3]
	    + pg3[iep] * pvec[4];
	  pdiff2[iep]	     
	    = pg1[iep] * pvec[3]
	    + pg2[iep] * pvec[1]
	    + pg3[iep] * pvec[5];
	  pdiff3[iep]	     
	    = pg1[iep] * pvec[4]
	    + pg2[iep] * pvec[5]
	    + pg3[iep] * pvec[2];
	}
      }
    }
    break;
  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg1 = FMF_PtrLevel( gc, iqp );
      pg2 = pg1 + nEP;
      
      for (ii = 0; ii < nRow; ii++) {
	pvec = FMF_PtrLevel( mtx, iqp ) + mtx->nCol * ii;
	pdiff1 = FMF_PtrLevel( diff, iqp ) + diff->nCol * ii;
	pdiff2 = pdiff1 + nEP;
	for (iep = 0; iep < nEP; iep++) {
	  pdiff1[iep]
	    = pg1[iep] * pvec[0]
	    + pg2[iep] * pvec[2];
	  pdiff2[iep]	     
	    = pg1[iep] * pvec[2]
	    + pg2[iep] * pvec[1];
	}
      }
    }
    break;
  case 1:
    for (iqp = 0; iqp < nQP; iqp++){
      pg1 = FMF_PtrLevel(gc, iqp);

      for (ii = 0; ii < nRow; ii++){
        pvec = FMF_PtrLevel(mtx, iqp) + mtx->nCol * ii;
        pdiff1 = FMF_PtrLevel(diff, iqp) + diff->nCol * ii;
        for (iep = 0; iep < nEP; iep++){
          pdiff1[iep] = pg1[iep] * pvec[0];
        }
      }
    }
    break;
  default:
    errput( ErrHead "ERR_Switch\n" );
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "build_nonsym_grad"
int32 build_nonsym_grad(FMField *out, FMField *gc)
{
  int32 iqp, ic, dim, idim, nEP, nQP;
  float64 *pout1, *pout2, *pout3, *pg;

  nEP = gc->nCol;
  nQP = gc->nLev;
  dim = gc->nRow;

#ifdef DEBUG_FMF
  if ((out->nCol != (dim * nEP))
      || (out->nRow != (dim * dim)) || (out->nLev != gc->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d), (%d %d %d)\n",
	    out->nLev, out->nRow, out->nCol,
	    gc->nLev, gc->nRow, gc->nCol );
  }
#endif

  fmf_fillC(out, 0.0);
  switch (dim) {
  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg = FMF_PtrLevel(gc, iqp);
      pout1 = FMF_PtrLevel(out, iqp);
      pout2 = pout1 + 5 * nEP;
      for (idim = 0; idim < dim; idim++) {
        for (ic = 0; ic < nEP; ic++) {
          pout1[ic] = pg[ic];
	  pout2[ic] = pg[ic];
        } /* for (ic) */
	pout1 += 2 * nEP;
	pout2 += 2 * nEP;
	pg += nEP;
      } /* for (idim) */
    } /* for (iqp) */
    break;
  case 3:
    for (iqp = 0; iqp < nQP; iqp++) {
      pg = FMF_PtrLevel(gc, iqp);
      pout1 = FMF_PtrLevel(out, iqp);
      pout2 = pout1 + 10 * nEP;
      pout3 = pout2 + 10 * nEP;
      for (idim = 0; idim < dim; idim++) {
        for (ic = 0; ic < nEP; ic++) {
          pout1[ic] = pg[ic];
          pout2[ic] = pg[ic];
          pout3[ic] = pg[ic];
	} /* for (ic) */
	pout1 += 3 * nEP;
	pout2 += 3 * nEP;
	pout3 += 3 * nEP;
	pg += nEP;
      } /* for (idim) */
    } /* for (iqp) */
    break;
  default:
    errput( ErrHead "ERR_Switch\n" );
    return( RET_Fail );
  }

  return(RET_OK);
}