#include "termsHyperElasticity.h"
#include "terms.h"

/*
  notes: disG -> mtxF => removed '1.0 +' in form_tlcc_buildOpB_VS3()
*/

static float64 trace[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};

#undef __FUNC__
#define __FUNC__ "form_tlcc_strainGreen_VS"
/*!
  Cauchy strain tensor stored as a vector, symmetric storage.
  
  @par Revision history:
  - 13.11.2001, c
  - 14.11.2001
  - 29.11.2001
  - 20.03.2003, adopted from rcfem2
*/
int32 form_tlcc_strainGreen_VS( FMField *strain, FMField *mtxF )
{
  int32 iqp, ii, ir, ic, ik;
  int32 sym = strain->nRow;
  int32 dim = mtxF->nRow;
  float64 *pstrain, *pf;

  for (iqp = 0; iqp < strain->nLev; iqp++) {
    pstrain = FMF_PtrLevel( strain, iqp );
    pf = FMF_PtrLevel( mtxF, iqp );
    for (ii = 0; ii < sym; ii++) {
      ir = t2i[ii];
      ic = t2j[ii];
      pstrain[ii] = 0.0;
      for (ik = 0; ik < dim; ik++) {
	pstrain[ii] += pf[dim*ik+ir] * pf[dim*ik+ic];
      }
    }
    for (ik = 0; ik < dim; ik++) {
      pstrain[ik] -= 1.0;
      pstrain[ik] *= 0.5;
    }
  }

  return( RET_OK );
}
#undef __FUNC__
#define __FUNC__ "form_tlcc_buildOpB_VS3"
/*!
  @par Revision history:
  - 13.11.2001, c
  - 16.11.2001
  - 20.03.2003, adopted from rcfem2
  - 26.08.2008: updated to work with mtxF instead of disG:
    - removed '1+' at diagonal items
    - indices to pd are transposed w.r.t. mafest1
*/
int32 form_tlcc_buildOpB_VS3( FMField *out, FMField *mtxF, FMField *gc )
{
  int32 iqp, iep, nQP, nEP;
  float64 *pgc, *pout, *pd, *pg[3];

  nQP = gc->nLev;
  nEP = gc->nCol;

  fmf_fillC( out, 0.0 );
  for (iqp = 0; iqp < nQP; iqp++) {
    pgc = FMF_PtrLevel( gc, iqp );
    pg[0] = pgc;
    pg[1] = pgc + nEP;
    pg[2] = pgc + 2 * nEP;

    pd = FMF_PtrLevel( mtxF, iqp );

    // Row 1.
    pout = FMF_PtrLevel( out, iqp );
    for (iep = 0; iep < nEP; iep++) {
      pout[iep] = (pd[0]) * pg[0][iep];
      pout[iep+nEP] = (pd[3]) * pg[0][iep];
      pout[iep+2*nEP] = (pd[6]) * pg[0][iep];
    }
    // Row 2.
    pout += 3 * nEP;
    for (iep = 0; iep < nEP; iep++) {
      pout[iep] = (pd[1]) * pg[1][iep];
      pout[iep+nEP] = (pd[4]) * pg[1][iep];
      pout[iep+2*nEP] = (pd[7]) * pg[1][iep];
    }
    // Row 3.
    pout += 3 * nEP;
    for (iep = 0; iep < nEP; iep++) {
      pout[iep] = (pd[2]) * pg[2][iep];
      pout[iep+nEP] = (pd[5]) * pg[2][iep];
      pout[iep+2*nEP] = (pd[8]) * pg[2][iep];
    }
    // Row 4.
    pout += 3 * nEP;
    for (iep = 0; iep < nEP; iep++) {
      pout[iep] = (pd[1]) * pg[0][iep]
	+ (pd[0]) * pg[1][iep];
      pout[iep+nEP] = (pd[4]) * pg[0][iep]
	+ (pd[3]) * pg[1][iep];
      pout[iep+2*nEP] = (pd[7]) * pg[0][iep]
	+ (pd[6]) * pg[1][iep];
    }
    // Row 5.
    pout += 3 * nEP;
    for (iep = 0; iep < nEP; iep++) {
      pout[iep] = (pd[2]) * pg[0][iep]
	+ (pd[0]) * pg[2][iep];
      pout[iep+nEP] = (pd[5]) * pg[0][iep]
	+ (pd[3]) * pg[2][iep];
      pout[iep+2*nEP] = (pd[8]) * pg[0][iep]
	+ (pd[6]) * pg[2][iep];
    }
    // Row 6.
    pout += 3 * nEP;
    for (iep = 0; iep < nEP; iep++) {
      pout[iep] = (pd[2]) * pg[1][iep]
	+ (pd[1]) * pg[2][iep];
      pout[iep+nEP] = (pd[5]) * pg[1][iep]
	+ (pd[4]) * pg[2][iep];
      pout[iep+2*nEP] = (pd[8]) * pg[1][iep]
	+ (pd[7]) * pg[2][iep];
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "form_tlcc_buildOpKtsC_VS3"
/*!
  @par Revision history:
  - 13.11.2001, c
  - 20.03.2003, adopted from rcfem2
*/
int32 form_tlcc_buildOpKtsC_VS3( FMField *out, FMField *tau, FMField *gc )
{
  int32 iqp, ir, ic, nQP, nEP;
  float64 *pgc, *pout, *ptau, *pg[3];

  nQP = gc->nLev;
  nEP = gc->nCol;

  fmf_fillC( out, 0.0 );
  for (iqp = 0; iqp < nQP; iqp++) {
    pgc = FMF_PtrLevel( gc, iqp );
    pg[0] = pgc;
    pg[1] = pgc + nEP;
    pg[2] = pgc + 2 * nEP;

    ptau = FMF_PtrLevel( tau, iqp );

    pout = FMF_PtrLevel( out, iqp );
    for (ir = 0; ir < nEP; ir++) {
      for (ic = 0; ic < nEP; ic++) {
	pout[ic]
	  = ptau[0] * pg[0][ir] * pg[0][ic]
	  + ptau[3] * pg[1][ir] * pg[0][ic]
	  + ptau[4] * pg[2][ir] * pg[0][ic]
	  + ptau[3] * pg[0][ir] * pg[1][ic]
	  + ptau[1] * pg[1][ir] * pg[1][ic]
	  + ptau[5] * pg[2][ir] * pg[1][ic]
	  + ptau[4] * pg[0][ir] * pg[2][ic]
	  + ptau[5] * pg[1][ir] * pg[2][ic]
	  + ptau[2] * pg[2][ir] * pg[2][ic];
      }
      pout += nEP;
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "dq_finite_strain_tl"
int32 dq_finite_strain_tl( FMField *mtxF, FMField *detF, FMField *vecCS,
			   FMField *trC, FMField *in2C, FMField *vecInvCS,
			   FMField *vecES,
			   FMField *state, int32 offset, VolumeGeometry *vg,
			   int32 *conn, int32 nEl, int32 nEP )
{
  int32 ii, id, iqp, nQP, dim, ret = RET_OK;
  FMField *st = 0, *mtx1 = 0, *mtx2 = 0;

  state->val = FMF_PtrFirst( state ) + offset;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  fmf_createAlloc( &st, 1, 1, nEP, dim );
  fmf_createAlloc( &mtx1, 1, nQP, dim, dim );
  fmf_createAlloc( &mtx2, 1, nQP, dim, dim );

  for (ii = 0; ii < nEl; ii++) {
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( mtxF, ii );
    FMF_SetCell( detF, ii );
    FMF_SetCell( vecCS, ii );
    FMF_SetCell( trC, ii );
    FMF_SetCell( in2C, ii );
    FMF_SetCell( vecInvCS, ii );
    FMF_SetCell( vecES, ii );

    // Deformation gradient.
    ele_extractNodalValuesNBN( st, state, conn + nEP * ii );
    fmf_mulATBT_1n( mtxF, st, vg->bfGM );
/*     fmf_mulAB_n1( mtxF, vg->bfGM, st ); */
    for (iqp = 0; iqp < nQP; iqp++) {
      for (id = 0; id < dim; id++) {
	mtxF->val[dim*(dim*iqp+id)+id] += 1.0;
      }
    }

    // Determinant of deformation gradient.
    geme_det3x3( detF->val, mtxF );
    for (iqp = 0; iqp < nQP; iqp++) {
      if (detF->val[iqp] <= MachEps) {
	errput( "warp violation %e at (iel: %d, iqp: %d)!\n",
		detF->val[iqp], ii, iqp );
      }
    }
    // Right Cauchy-Green tensor C = F^T F.
    fmf_mulATB_nn( mtx1, mtxF, mtxF );
    geme_tensor2vectorS3( vecCS, mtx1 );
    // trace C.
    geme_trace3x3( trC->val, mtx1 );
    // I_2 of C.
    fmf_mulATB_nn( mtx2, mtx1, mtx1 );
    geme_trace3x3( in2C->val, mtx2 );
    for (iqp = 0; iqp < nQP; iqp++) {
      in2C->val[iqp] = 0.5 * (trC->val[iqp] *  trC->val[iqp]
			      - in2C->val[iqp]);
    }
    // C^{-1} as a vector, symmetric storage.
    geme_invert3x3( mtx2, mtx1 );
    geme_tensor2vectorS3( vecInvCS, mtx2 );

    form_tlcc_strainGreen_VS( vecES, mtxF );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &mtx1 ); 
  fmf_freeDestroy( &mtx2 ); 

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_tl_he_rtm"
int32 dw_tl_he_rtm( FMField *out,
		    FMField *stress, FMField *tan_mod, FMField *mtxF,
		    VolumeGeometry *vg,
		    int32 *elList, int32 elList_nRow, int32 isDiff )
{
  int32 ii, iel, sym, nRow, nQP, ret = RET_OK;
  FMField *mtxB = 0, *out_qp = 0, *btd = 0, *btdb = 0, *ktsc = 0, *iktsc = 0;

  nQP = vg->bfGM->nLev;
  sym = stress->nRow;
  nRow = out->nRow; // dim * nEP.

  fmf_createAlloc( &mtxB, 1, nQP, sym, nRow );

  if (isDiff) {
    int32 nEP = vg->nEP;

    fmf_createAlloc( &btd, 1, nQP, nRow, sym );
    fmf_createAlloc( &btdb, 1, nQP, nRow, nRow );
    fmf_createAlloc( &ktsc, 1, nQP, nEP, nEP );
    fmf_createAlloc( &iktsc, 1, 1, nEP, nEP );

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];

      FMF_SetCell( out, ii );
      FMF_SetCell( stress, iel );
      FMF_SetCell( tan_mod, iel );
      FMF_SetCell( mtxF, iel );
      FMF_SetCell( vg->bfGM, iel );
      FMF_SetCell( vg->det, iel );

      /* B^T D B. */
      form_tlcc_buildOpB_VS3( mtxB, mtxF, vg->bfGM );
      fmf_mulATB_nn( btd, mtxB, tan_mod );
      fmf_mulAB_nn( btdb, btd, mtxB );
      fmf_sumLevelsMulF( out, btdb, vg->det->val );

      /* + K_{t\sigma}. */
      form_tlcc_buildOpKtsC_VS3( ktsc, stress, vg->bfGM );
      fmf_sumLevelsMulF( iktsc, ktsc, vg->det->val );

      fmfr_addA_blockNC( out, iktsc, 0, 0 );
      fmfr_addA_blockNC( out, iktsc, nEP, nEP );
      fmfr_addA_blockNC( out, iktsc, 2 * nEP, 2 * nEP );

      ERR_CheckGo( ret );
    }
  } else {
    fmf_createAlloc( &out_qp, 1, nQP, nRow, 1 );
    
    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];

      FMF_SetCell( out, ii );
      FMF_SetCell( stress, iel );
      FMF_SetCell( mtxF, iel );
      FMF_SetCell( vg->bfGM, iel );
      FMF_SetCell( vg->det, iel );

      form_tlcc_buildOpB_VS3( mtxB, mtxF, vg->bfGM );
      fmf_mulATB_nn( out_qp, mtxB, stress );
      fmf_sumLevelsMulF( out, out_qp, vg->det->val );

      ERR_CheckGo( ret );
    }
  }

 end_label:
  fmf_freeDestroy( &mtxB );

  if (isDiff) {
    fmf_freeDestroy( &btd );
    fmf_freeDestroy( &btdb );
    fmf_freeDestroy( &ktsc );
    fmf_freeDestroy( &iktsc );
  } else {
    fmf_freeDestroy( &out_qp );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_tl_he_stress_neohook"
int32 dq_tl_he_stress_neohook( FMField *out, FMField *mat,
			       FMField *detF, FMField *trC, FMField *vecInvCS )
{
  int32 ii, iqp, ir, sym, nQP, ret = RET_OK;
  float64 detF23;
  float64 *pmu, *pstress, *ptrC, *pinvC, *pdetF;

  nQP = detF->nLev;
  sym = out->nRow;
  
  if (sym != 6) errput( "Hyperelastic materials are 3D only!\n" );
  ERR_CheckGo( ret );

/*   output( "%d\n", sym ); */
  for (ii = 0; ii < out->nCell; ii++) {
    pdetF = FMF_PtrCell( detF, ii );
    ptrC = FMF_PtrCell( trC, ii );
    pinvC = FMF_PtrCell( vecInvCS, ii );

    pstress = FMF_PtrCell( out, ii );
    pmu = FMF_PtrCell( mat, ii );
    for (iqp = 0; iqp < nQP; iqp++) {
      detF23 = exp( -2.0/3.0 * log( pdetF[iqp] ) );
      for (ir = 0; ir < sym; ir++) {
	pstress[ir]
	  = pmu[iqp] * detF23 * (trace[ir] - ptrC[iqp]/3.0 * pinvC[ir]);
/* 	output( "%d %d %d %f %f %f %f\n", ii, iqp, ir, */
/* 		pmu[iqp], ptrC[iqp], pinvC[ir], pstress[ir] ); */
      }
      pstress += sym;
      pinvC += sym;
    }
  }

 end_label:
  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_tl_he_stress_bulk"
int32 dq_tl_he_stress_bulk( FMField *out,FMField *mat,
			    FMField *detF, FMField *vecInvCS )
{
  int32 ii, iqp, ir, sym, nQP, ret = RET_OK;
  float64 *pbulk, *pstress, *pinvC, *pdetF;

  nQP = detF->nLev;
  sym = out->nRow;
  
  if (sym != 6) errput( "Hyperelastic materials are 3D only!\n" );
  ERR_CheckGo( ret );

  for (ii = 0; ii < out->nCell; ii++) {
    pdetF = FMF_PtrCell( detF, ii );
    pinvC = FMF_PtrCell( vecInvCS, ii );

    pstress = FMF_PtrCell( out, ii );
    pbulk = FMF_PtrCell( mat, ii );

    for (iqp = 0; iqp < nQP; iqp++) {
      // Volumetric part.
      for (ir = 0; ir < sym; ir++) {
	pstress[ir]
	  = pbulk[iqp] * pdetF[iqp] * (pdetF[iqp] - 1.0) * pinvC[ir];
      }
      pstress += sym;
      pinvC += sym;
    }
  }

 end_label:
  return( ret );
}


#undef __FUNC__
#define __FUNC__ "dq_tl_he_tan_mod_neohook"
int32 dq_tl_he_tan_mod_neohook( FMField *out, FMField *mat,
				FMField *detF, FMField *trC, FMField *vecInvCS )
{
  int32 ii, nQP, ir, ic, iqp, sym, ret = RET_OK;
  float64 cc, c1, c2, c3, detF23;
  float64 *pd;
  float64 *pmu;
  float64 *pinvC, *ptrC, *pdetF, *pinvC2_ikjl, *pinvC2_iljk;
  FMField *invC2_ikjl = 0, *invC2_iljk = 0;

  sym = out->nRow;
  nQP = out->nLev;

  fmf_createAlloc( &invC2_ikjl, 1, nQP, sym, sym );
  fmf_createAlloc( &invC2_iljk, 1, nQP, sym, sym );

  pinvC2_ikjl = FMF_PtrCurrent( invC2_ikjl );
  pinvC2_iljk = FMF_PtrCurrent( invC2_iljk );

  for (ii = 0; ii < out->nCell; ii++) {
    pdetF = FMF_PtrCell( detF, ii );
    ptrC = FMF_PtrCell( trC, ii );
    pinvC = FMF_PtrCell( vecInvCS, ii );
    FMF_SetCell( vecInvCS, ii );
    pd = FMF_PtrCell( out, ii );
    pmu = FMF_PtrCell( mat, ii );

    geme_mulT2ST2S_T4S_ikjl( invC2_ikjl, vecInvCS, vecInvCS );
    geme_mulT2ST2S_T4S_iljk( invC2_iljk, vecInvCS, vecInvCS );


    for (iqp = 0; iqp < nQP; iqp++) {
      detF23 = exp( -2.0/3.0 * log( pdetF[iqp] ) );
      cc = pmu[iqp] * detF23;
      c1 = 2.0/9.0 * cc * ptrC[iqp];
      c2 = - 2.0/3.0 * cc;
      c3 = cc / 3.0 * ptrC[iqp];
      for (ir = 0; ir < sym; ir++) {
	for (ic = 0; ic < sym; ic++) {
	  pd[sym*ir+ic]
	    = c1 * pinvC[sym*iqp+ir] * pinvC[sym*iqp+ic]
	    + c2 * ((pinvC[sym*iqp+ir] * trace[ic])
		  + (pinvC[sym*iqp+ic] * trace[ir]))
	    + c3 * (pinvC2_ikjl[sym*(sym*iqp+ir)+ic]
		    + pinvC2_iljk[sym*(sym*iqp+ir)+ic]);
	}
      }
      pd += sym * sym;
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &invC2_ikjl );
  fmf_freeDestroy( &invC2_iljk );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_tl_he_tan_mod_bulk"
int32 dq_tl_he_tan_mod_bulk( FMField *out, FMField *mat,
			     FMField *detF, FMField *vecInvCS )
{
  int32 ii, nQP, ir, ic, iqp, sym, ret = RET_OK;
  float64 cbulk21, cbulk22;
  float64 *pd;
  float64 *pbulk;
  float64 *pinvC, *pdetF, *pinvC2_ikjl, *pinvC2_iljk;
  FMField *invC2_ikjl = 0, *invC2_iljk = 0;

  sym = out->nRow;
  nQP = out->nLev;

  fmf_createAlloc( &invC2_ikjl, 1, nQP, sym, sym );
  fmf_createAlloc( &invC2_iljk, 1, nQP, sym, sym );

  pinvC2_ikjl = FMF_PtrCurrent( invC2_ikjl );
  pinvC2_iljk = FMF_PtrCurrent( invC2_iljk );

  for (ii = 0; ii < out->nCell; ii++) {
    pdetF = FMF_PtrCell( detF, ii );
    pinvC = FMF_PtrCell( vecInvCS, ii );
    FMF_SetCell( vecInvCS, ii );
    pd = FMF_PtrCell( out, ii );
    pbulk = FMF_PtrCell( mat, ii );

    geme_mulT2ST2S_T4S_ikjl( invC2_ikjl, vecInvCS, vecInvCS );
    geme_mulT2ST2S_T4S_iljk( invC2_iljk, vecInvCS, vecInvCS );

    for (iqp = 0; iqp < nQP; iqp++) {
      cbulk21 = pbulk[iqp] * (pdetF[iqp] * (pdetF[iqp] - 1.0));
      cbulk22 = pbulk[iqp] * pdetF[iqp] * pdetF[iqp];
      for (ir = 0; ir < sym; ir++) {
	for (ic = 0; ic < sym; ic++) {
	  pd[sym*ir+ic]
	    = (cbulk21 + cbulk22) * pinvC[sym*iqp+ir] * pinvC[sym*iqp+ic]
	    - cbulk21 * (pinvC2_ikjl[sym*(sym*iqp+ir)+ic]
			+ pinvC2_iljk[sym*(sym*iqp+ir)+ic]);
	}
      }
      pd += sym * sym;
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &invC2_ikjl );
  fmf_freeDestroy( &invC2_iljk );

  return( ret );
}
