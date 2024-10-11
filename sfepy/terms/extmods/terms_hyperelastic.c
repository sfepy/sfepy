#include "terms_hyperelastic.h"
#include "terms.h"
#include "form_sdcc.h"

/*
  notes: disG -> mtxF => removed '1.0 +' in form_tlcc_buildOpB_VS3()
*/

static float64 trace3d[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
static float64 trace2d[3] = {1.0, 1.0, 0.0};
static float64 trace1d[1] = {1.0};

#undef __FUNC__
#define __FUNC__ "get_trace"
float64 *get_trace( int32 sym )
{
  if (sym == 1)
    return( trace1d );
  if (sym == 3)
    return( trace2d );
  if (sym == 6)
    return( trace3d );

  return( NULL );
}

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
  int32 *t2i = 0, *t2j = 0;

  switch (dim) {
  case 1:
    t2i = t2i1D;
    t2j = t2j1D;
    break;
  case 2:
    t2i = t2i2D;
    t2j = t2j2D;
    break;
  case 3:
    t2i = t2i3D;
    t2j = t2j3D;
    break;
  default:
    errput( ErrHead "ERR_Switch\n" );
  }

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
    }
    for (ii = 0; ii < sym; ii++) {
      pstrain[ii] *= 0.5;
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
  int32 iqp, iep, nQP, nEP, dim;
  float64 *pgc, *pout, *pd, *pg[3];

  nQP = gc->nLev;
  nEP = gc->nCol;
  dim = gc->nRow;

  fmf_fillC( out, 0.0 );
  switch (dim) {
  case 3:
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
    } /* for (iqp) */
    break;
  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pgc = FMF_PtrLevel( gc, iqp );
      pg[0] = pgc;
      pg[1] = pgc + nEP;

      pd = FMF_PtrLevel( mtxF, iqp );

      // Row 1.
      pout = FMF_PtrLevel( out, iqp );
      for (iep = 0; iep < nEP; iep++) {
        pout[iep] = (pd[0]) * pg[0][iep];
        pout[iep+nEP] = (pd[2]) * pg[0][iep];
      }
      // Row 2.
      pout += 2 * nEP;
      for (iep = 0; iep < nEP; iep++) {
        pout[iep] = (pd[1]) * pg[1][iep];
        pout[iep+nEP] = (pd[3]) * pg[1][iep];
      }
      // Row 3.
      pout += 2 * nEP;
      for (iep = 0; iep < nEP; iep++) {
        pout[iep] = (pd[1]) * pg[0][iep]
          + (pd[0]) * pg[1][iep];
        pout[iep+nEP] = (pd[3]) * pg[0][iep]
          + (pd[2]) * pg[1][iep];
      }
    } /* for (iqp) */
    break;
  case 1:
    for (iqp = 0; iqp < nQP; iqp++) {
      pgc = FMF_PtrLevel( gc, iqp );
      pg[0] = pgc;

      pd = FMF_PtrLevel( mtxF, iqp );

      // Row 1.
      pout = FMF_PtrLevel( out, iqp );
      for (iep = 0; iep < nEP; iep++) {
        pout[iep] = (pd[0]) * pg[0][iep];
      }
    } /* for (iqp) */
    break;
  } /* switch */

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
  int32 iqp, ir, ic, nQP, nEP, dim;
  float64 *pgc, *pout, *ptau, *pg[3];

  nQP = gc->nLev;
  nEP = gc->nCol;
  dim = gc->nRow;

  fmf_fillC( out, 0.0 );
  switch (dim) {
  case 3:
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
    } /* for (iqp) */
    break;
  case 2:
    for (iqp = 0; iqp < nQP; iqp++) {
      pgc = FMF_PtrLevel( gc, iqp );
      pg[0] = pgc;
      pg[1] = pgc + nEP;

      ptau = FMF_PtrLevel( tau, iqp );

      pout = FMF_PtrLevel( out, iqp );
      for (ir = 0; ir < nEP; ir++) {
        for (ic = 0; ic < nEP; ic++) {
          pout[ic]
            = ptau[0] * pg[0][ir] * pg[0][ic]
            + ptau[2] * pg[1][ir] * pg[0][ic]
            + ptau[2] * pg[0][ir] * pg[1][ic]
            + ptau[1] * pg[1][ir] * pg[1][ic];
        }
        pout += nEP;
      }
    } /* for (iqp) */
    break;
  case 1:
    for (iqp = 0; iqp < nQP; iqp++) {
      pgc = FMF_PtrLevel( gc, iqp );
      pg[0] = pgc;

      ptau = FMF_PtrLevel( tau, iqp );

      pout = FMF_PtrLevel( out, iqp );
      for (ir = 0; ir < nEP; ir++) {
        for (ic = 0; ic < nEP; ic++) {
          pout[ic]
            = ptau[0] * pg[0][ir] * pg[0][ic];
        }
        pout += nEP;
      }
    } /* for (iqp) */
    break;
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "dq_finite_strain"
int32 dq_finite_strain( FMField *mtxF, FMField *detF, FMField *vecCS,
                        FMField *trC, FMField *in2C, FMField *vecInvCS,
                        FMField *vecES,
                        FMField *state, int32 offset, Mapping *vg,
                        int32 *conn, int32 nEl, int32 nEP, int32 mode_ul)
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
    FMF_SetCell( vecES, ii );
    if ( vecInvCS > 0 ) {
      FMF_SetCell( vecInvCS, ii );
    }

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
    if (mode_ul) {
      // Left Cauchy-Green tensor b = F F^T.
      fmf_mulABT_nn( mtx1, mtxF, mtxF );
    } else {
      // Right Cauchy-Green tensor C = F^T F.
      fmf_mulATB_nn( mtx1, mtxF, mtxF );
    }
    geme_tensor2vectorS3( vecCS, mtx1 );
    // trace C.
/*     geme_trace3x3( trC->val, mtx1 ); */
    geme_invar1( trC->val, mtx1 );
    // I_2 of C.
    geme_invar2( in2C->val, mtx1 );
    if ( vecInvCS > 0 ) {
      // C^{-1} as a vector, symmetric storage.
      geme_invert3x3( mtx2, mtx1 );
      geme_tensor2vectorS3( vecInvCS, mtx2 );
    }
    form_tlcc_strainGreen_VS( vecES, mtxF );

    ERR_CheckGo( ret );
  }
 end_label:
  errclear(); // Prevent false memory errors in mem_free_mem().

  fmf_freeDestroy( &st );
  fmf_freeDestroy( &mtx1 );
  fmf_freeDestroy( &mtx2 );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_finite_strain_tl"
int32 dq_finite_strain_tl( FMField *mtxF, FMField *detF, FMField *vecCS,
                           FMField *trC, FMField *in2C, FMField *vecInvCS,
                           FMField *vecES,
                           FMField *state, int32 offset, Mapping *vg,
                           int32 *conn, int32 nEl, int32 nEP )
{
  return( dq_finite_strain( mtxF, detF, vecCS, trC, in2C, vecInvCS, vecES,
                            state, offset, vg, conn, nEl, nEP, 0 ) );
}

#undef __FUNC__
#define __FUNC__ "dq_finite_strain_ul"
int32 dq_finite_strain_ul( FMField *mtxF, FMField *detF, FMField *vecBS,
                           FMField *trB, FMField *in2B, FMField *vecES,
                           FMField *state, int32 offset, Mapping *vg,
                           int32 *conn, int32 nEl, int32 nEP )
{
  return( dq_finite_strain( mtxF, detF, vecBS, trB, in2B, 0, vecES,
                            state, offset, vg, conn, nEl, nEP, 1 ) );
}

#undef __FUNC__
#define __FUNC__ "dq_tl_finite_strain_surface"
int32 dq_tl_finite_strain_surface( FMField *mtxF, FMField *detF, FMField *mtxFI,
                                   FMField *state, int32 offset,
                                   Mapping *sg,
                                   int32 *fis, int32 nFa, int32 nFP,
                                   int32 *conn, int32 nEl, int32 nEP)
{
  int32 ii, iel, id, iqp, nQP, dim, ret = RET_OK;
  FMField *st = 0;

  state->val = FMF_PtrFirst( state ) + offset;

  nQP = sg->bfGM->nLev;
  dim = sg->bfGM->nRow;

  fmf_createAlloc( &st, 1, 1, nEP, dim );

  for (ii = 0; ii < nFa; ii++) {
    iel = fis[ii*nFP+0];

    FMF_SetCell( sg->bfGM, ii );
    FMF_SetCell( mtxF, ii );
    FMF_SetCell( mtxFI, ii );
    FMF_SetCell( detF, ii );

    // Deformation gradient.
    ele_extractNodalValuesNBN( st, state, conn + nEP * iel );
    fmf_mulATBT_1n( mtxF, st, sg->bfGM );
    for (iqp = 0; iqp < nQP; iqp++) {
      for (id = 0; id < dim; id++) {
        mtxF->val[dim*(dim*iqp+id)+id] += 1.0;
      }
    }

    // Determinant of deformation gradient.
    geme_det3x3( detF->val, mtxF );
    geme_invert3x3( mtxFI, mtxF );

    ERR_CheckGo( ret );
  }
 end_label:
  errclear(); // Prevent false memory errors in mem_free_mem().

  fmf_freeDestroy( &st );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_he_rtm"
int32 dw_he_rtm( FMField *out,
                 FMField *stress, FMField *tan_mod,
                 FMField *mtxF, FMField *detF,
                 Mapping *vg,
                 int32 isDiff, int32 mode_ul )
{
  int32 ii, j, sym, nRow, nQP, ret = RET_OK, dim;
  FMField *aux = 0, *out_qp = 0, *btd = 0, *btdb = 0, *ktsc = 0, *iktsc = 0;

  nQP = vg->bfGM->nLev;
  sym = stress->nRow;
  nRow = out->nRow; // dim * nEP.
  dim = vg->dim;

  if (mode_ul) {
    fmf_createAlloc( &aux, 1, 1, 1, nQP );
  }
  else {
    fmf_createAlloc( &aux, 1, nQP, sym, nRow );
  }

  if (isDiff) {
    int32 nEP = vg->bfGM->nCol;

    fmf_createAlloc( &btd, 1, nQP, nRow, sym );
    fmf_createAlloc( &btdb, 1, nQP, nRow, nRow );
    fmf_createAlloc( &ktsc, 1, nQP, nEP, nEP );
    fmf_createAlloc( &iktsc, 1, 1, nEP, nEP );

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCell( stress, ii );
      FMF_SetCell( tan_mod, ii );
      FMF_SetCell( vg->bfGM, ii );

      FMF_SetCell( vg->det, ii );

      /* B^T D B. */
      if (mode_ul) {
        /* ULF */
        FMF_SetCell( detF, ii );
        for (j = 0; j < nQP; j++) /* det/J */
          aux->val[j] = vg->det->val[j] / detF->val[j];

        form_sdcc_actOpGT_M3( btd, vg->bfGM, tan_mod );
        form_sdcc_actOpG_RM3( btdb, btd, vg->bfGM );
        fmf_sumLevelsMulF( out, btdb, aux->val );
        /* + K_{t\sigma}. */
        form_tlcc_buildOpKtsC_VS3( ktsc, stress, vg->bfGM );
        fmf_sumLevelsMulF( iktsc, ktsc, aux->val );
      }
      else {
        /* TLF */
        FMF_SetCell( mtxF, ii );
        form_tlcc_buildOpB_VS3( aux, mtxF, vg->bfGM );
        fmf_mulATB_nn( btd, aux, tan_mod );
        fmf_mulAB_nn( btdb, btd, aux );
        fmf_sumLevelsMulF( out, btdb, vg->det->val );

        /* + K_{t\sigma}. */
        form_tlcc_buildOpKtsC_VS3( ktsc, stress, vg->bfGM );
        fmf_sumLevelsMulF( iktsc, ktsc, vg->det->val );
      }

      fmfr_addA_blockNC( out, iktsc, 0, 0 );
      if (dim > 1)
        fmfr_addA_blockNC( out, iktsc, nEP, nEP );
      if (dim > 2)
        fmfr_addA_blockNC( out, iktsc, 2 * nEP, 2 * nEP );

      ERR_CheckGo( ret );
    }
  } else {
    fmf_createAlloc( &out_qp, 1, nQP, nRow, 1 );

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCell( stress, ii );
      FMF_SetCell( vg->bfGM, ii );
      FMF_SetCell( vg->det, ii );

      if (mode_ul) {
        /* ULF */
        FMF_SetCell( detF, ii );
        for (j = 0; j < nQP; j++) /* det/J */
          aux->val[j] = vg->det->val[j] / detF->val[j];

        form_sdcc_actOpGT_M3( out_qp, vg->bfGM, stress );
        fmf_sumLevelsMulF( out, out_qp, aux->val );
      }
      else {
        /* TLF */
        FMF_SetCell( mtxF, ii );
        form_tlcc_buildOpB_VS3( aux, mtxF, vg->bfGM );
        fmf_mulATB_nn( out_qp, aux, stress );
        fmf_sumLevelsMulF( out, out_qp, vg->det->val );
      }

      ERR_CheckGo( ret );
    }
  }

 end_label:

  fmf_freeDestroy( &aux );

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
#define __FUNC__ "de_he_rtm"
int32 de_he_rtm( FMField *out,
                 FMField *stress, FMField *detF,
                 Mapping *vg,
                 int32 *elList, int32 elList_nRow,
                 int32 mode_ul )
{
  int32 ii, iel, nQP, ret = RET_OK, j;
  FMField *aux = 0;

  nQP = vg->det->nLev;

  if (mode_ul) {
    fmf_createAlloc( &aux, 1, 1, 1, nQP );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( stress, iel );
    FMF_SetCell( vg->det, iel );
    FMF_SetCell( detF, iel );

    if (mode_ul) {
      /* ULF */
      for (j = 0; j < nQP; j++) /* det/J */
        aux->val[j] = vg->det->val[j] / detF->val[j];
      fmf_sumLevelsMulF( out, stress, aux->val );
    }
    else {
      /* TLF */
      fmf_sumLevelsMulF( out, stress, vg->det->val );
    }

    ERR_CheckGo( ret );
  }

 end_label:

  fmf_freeDestroy( &aux );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_he_stress_bulk"
int32 dq_he_stress_bulk( FMField *out, FMField *mat,
                         FMField *detF, FMField *vecCG,
                         int32 mode_ul )
{
  int32 ii, iqp, ir, sym, nQP, ret = RET_OK;
  float64 *pbulk, *pstress, *pCG, *pdetF, *ptrace;

  nQP = detF->nLev;
  sym = out->nRow;
  ptrace = get_trace( sym );

  for (ii = 0; ii < out->nCell; ii++) {
    pdetF = FMF_PtrCell( detF, ii );

    pstress = FMF_PtrCell( out, ii );
    pbulk = FMF_PtrCellX1( mat, ii );

    if (mode_ul) {
      /* ULF */
      for (iqp = 0; iqp < nQP; iqp++) {
        // Volumetric part.
        for (ir = 0; ir < sym; ir++) {
          /* Crisfield II., (13.35) */
          pstress[ir]
            = pbulk[iqp] * pdetF[iqp] * (pdetF[iqp] - 1.0) * ptrace[ir];
        }
        pstress += sym;
      }
    }
    else {
      /* TLF */
      pCG = FMF_PtrCell( vecCG, ii );
      for (iqp = 0; iqp < nQP; iqp++) {
        // Volumetric part.
        for (ir = 0; ir < sym; ir++) {
          /* Crisfield II., (13.35) */
          pstress[ir]
            = pbulk[iqp] * pdetF[iqp] * (pdetF[iqp] - 1.0) * pCG[ir];
        }
        pstress += sym;
        pCG += sym;
      }
    } /* if (mode_ul) */
    ERR_CheckGo( ret );
  }

 end_label:
  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_tl_he_stress_bulk"
int32 dq_tl_he_stress_bulk( FMField *out,FMField *mat,
                            FMField *detF, FMField *vecInvCS )
{
  return( dq_he_stress_bulk( out, mat, detF, vecInvCS, 0 ) );
}

#undef __FUNC__
#define __FUNC__ "dq_ul_he_stress_bulk"
int32 dq_ul_he_stress_bulk( FMField *out,FMField *mat,
                            FMField *detF )
{
  return( dq_he_stress_bulk( out, mat, detF, NULL, 1 ) );
}

#undef __FUNC__
#define __FUNC__ "dq_tl_he_stress_bulk_active"
int32 dq_tl_he_stress_bulk_active( FMField *out,FMField *mat,
                                   FMField *detF, FMField *vecCG )
{
  int32 ii, iqp, ir, sym, nQP, ret = RET_OK;
  float64 *pmat, *pstress, *pCG, *pdetF;

  nQP = detF->nLev;
  sym = out->nRow;

  for (ii = 0; ii < out->nCell; ii++) {
    pdetF = FMF_PtrCell( detF, ii );

    pstress = FMF_PtrCell( out, ii );
    pmat = FMF_PtrCellX1( mat, ii );

    pCG = FMF_PtrCell( vecCG, ii );
    for (iqp = 0; iqp < nQP; iqp++) {
      // Volumetric part.
      for (ir = 0; ir < sym; ir++) {
        pstress[ir] = pmat[iqp] * pdetF[iqp] * pCG[ir];
      }
      pstress += sym;
      pCG += sym;
    }
    ERR_CheckGo( ret );
  }

 end_label:
  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_he_stress_neohook"
int32 dq_he_stress_neohook( FMField *out, FMField *mat,
                            FMField *detF, FMField *trC, FMField *vecCG,
                            int32 mode_ul )
{
  int32 ii, iqp, ir, sym, nQP, ret = RET_OK;
  float64 detF23;
  float64 *pmu, *pstress, *ptrC, *pCG, *pdetF, *ptrace;

  nQP = detF->nLev;
  sym = out->nRow;
  ptrace = get_trace( sym );

  /*   output( "%d\n", sym ); */
  for (ii = 0; ii < out->nCell; ii++) {
    pdetF = FMF_PtrCell( detF, ii );
    ptrC = FMF_PtrCell( trC, ii );
    pCG = FMF_PtrCell( vecCG, ii );

    pstress = FMF_PtrCell( out, ii );
    pmu = FMF_PtrCellX1( mat, ii );

    if (mode_ul) {
      /* ULF */
      for (iqp = 0; iqp < nQP; iqp++) {
        detF23 = exp( -2.0/3.0 * log( pdetF[iqp] ) );
        for (ir = 0; ir < sym; ir++) {
          /* Crisfield II., (13.35) */
          pstress[ir]
            = pmu[iqp] * detF23 * (pCG[ir] - ptrC[iqp]/3.0 * ptrace[ir]);
        }
        pstress += sym;
        pCG += sym;
      }
    }
    else {
      /* TLF */
      for (iqp = 0; iqp < nQP; iqp++) {
        detF23 = exp( -2.0/3.0 * log( pdetF[iqp] ) );
        for (ir = 0; ir < sym; ir++) {
          /* Crisfield II., (13.35) */
          pstress[ir]
            = pmu[iqp] * detF23 * (ptrace[ir] - ptrC[iqp]/3.0 * pCG[ir]);
          /*    output( "%d %d %d %f %f %f %f\n", ii, iqp, ir, */
          /*            pmu[iqp], ptrC[iqp], pinvC[ir], pstress[ir] ); */
        }
        pstress += sym;
        pCG += sym;
      }
    } /* if (mode_ul) */
    ERR_CheckGo( ret );
  }

 end_label:
  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_tl_he_stress_neohook"
int32 dq_tl_he_stress_neohook( FMField *out, FMField *mat,
                               FMField *detF, FMField *trC, FMField *vecInvCS )
{
  return ( dq_he_stress_neohook( out, mat, detF, trC, vecInvCS, 0 ) );
}

#undef __FUNC__
#define __FUNC__ "dq_ul_he_stress_neohook"
int32 dq_ul_he_stress_neohook( FMField *out, FMField *mat,
                               FMField *detF, FMField *trB, FMField *vecBS )
{
  return ( dq_he_stress_neohook( out, mat, detF, trB, vecBS, 1 ) );
}

#undef __FUNC__
#define __FUNC__ "dq_tl_he_stress_mooney_rivlin"
int32 dq_tl_he_stress_mooney_rivlin( FMField *out, FMField *mat,
                                     FMField *detF, FMField *trC,
                                     FMField *vecInvCS, FMField *vecCS,
                                     FMField *in2C )
{
  int32 ii, iqp, ir, sym, nQP, ret = RET_OK;
  float64 detF23;
  float64 *pkappa, *pstress, *ptrC, *pinvC, *pdetF, *pC, *pin2C, *ptrace;

  nQP = detF->nLev;
  sym = out->nRow;
  ptrace = get_trace( sym );

/*   output( "%d\n", sym ); */
  for (ii = 0; ii < out->nCell; ii++) {
    pdetF = FMF_PtrCell( detF, ii );
    ptrC = FMF_PtrCell( trC, ii );
    pinvC = FMF_PtrCell( vecInvCS, ii );
    pC = FMF_PtrCell( vecCS, ii );
    pin2C = FMF_PtrCell( in2C, ii );

    pstress = FMF_PtrCell( out, ii );
    pkappa = FMF_PtrCellX1( mat, ii );
    for (iqp = 0; iqp < nQP; iqp++) {
      detF23 = exp( -2.0/3.0 * log( pdetF[iqp] ) );
      for (ir = 0; ir < sym; ir++) {
        pstress[ir]
          = pkappa[iqp] * detF23 * detF23
          * (ptrC[iqp] * ptrace[ir] - pC[ir]
             - (2.0/3.0) * pin2C[iqp] * pinvC[ir]);
      }
      pstress += sym;
      pinvC += sym;
      pC += sym;
    }
    ERR_CheckGo( ret );
  }

 end_label:
  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_ul_he_stress_mooney_rivlin"
/*!
  Mooney-Rivlin ( c1, c2 )
  Kirchhoff stress tensor
  tau =  B1 * b + B2 * b*b + B3 * I
*/

int32 dq_ul_he_stress_mooney_rivlin( FMField *out, FMField *mat,
                                     FMField *detF, FMField *trB,
                                     FMField *vecBS, FMField *in2B )
{
  int32 ii, iqp, ir, sym, nQP, ret = RET_OK;
  float64 detF23;
  float64 *pkappa, *pstress, *ptrB, *pB, *pBB, *pdetF,  *pin2B, *ptrace;
  FMField *vecBB;

  nQP = detF->nLev;
  sym = out->nRow;
  ptrace = get_trace( sym );

  fmf_createAlloc( &vecBB, 1, nQP, sym, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    pdetF = FMF_PtrCell( detF, ii );
    ptrB = FMF_PtrCell( trB, ii );
    pin2B = FMF_PtrCell( in2B, ii );
    pB = FMF_PtrCell( vecBS, ii );
    pBB = vecBB->val0;

    FMF_SetCell( vecBS, ii );
    geme_mulT2S_AA( vecBB, vecBS );

    pstress = FMF_PtrCell( out, ii );
    pkappa = FMF_PtrCellX1( mat, ii );
    for (iqp = 0; iqp < nQP; iqp++) {
      detF23 = exp( -2.0/3.0 * log( pdetF[iqp] ) );
      for (ir = 0; ir < sym; ir++) {
        pstress[ir]
          = pkappa[iqp] * detF23 * detF23
          * (ptrB[iqp] * pB[ir] - pBB[ir]
             - (2.0/3.0) * pin2B[iqp] * ptrace[ir]);
      }
      pstress += sym;
      pBB += sym;
      pB += sym;
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &vecBB );

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
    FMF_SetCell( vecInvCS, ii );
    pdetF = FMF_PtrCell( detF, ii );
    pinvC = FMF_PtrCell( vecInvCS, ii );
    pd = FMF_PtrCell( out, ii );
    pbulk = FMF_PtrCellX1( mat, ii );

    geme_mulT2ST2S_T4S_ikjl( invC2_ikjl, vecInvCS, vecInvCS );
    geme_mulT2ST2S_T4S_iljk( invC2_iljk, vecInvCS, vecInvCS );

    for (iqp = 0; iqp < nQP; iqp++) {
      cbulk21 = pbulk[iqp] * pdetF[iqp] * (pdetF[iqp] - 1.0);
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

#undef __FUNC__
#define __FUNC__ "dq_ul_he_tan_mod_bulk"
int32 dq_ul_he_tan_mod_bulk( FMField *out, FMField *mat, FMField *detF )
{
  int32 ii, nQP, ir, ic, iqp, sym, ret = RET_OK;
  float64 cbulk21, cbulk22;
  float64 *pd;
  float64 *pbulk;
  float64 *pdetF, *p_ikjl, *p_iljk, *ptrace;
  FMField *ikjl = 0, *iljk = 0;
  FMField traceVec[1];

  sym = out->nRow;
  nQP = out->nLev;
  ptrace = get_trace( sym );

  fmf_createAlloc( &ikjl, 1, 1, sym, sym );
  fmf_createAlloc( &iljk, 1, 1, sym, sym );
  traceVec->nAlloc = -1;
  fmf_pretend( traceVec, 1, 1, sym, 1, ptrace );

  p_ikjl = FMF_PtrCurrent( ikjl );
  p_iljk = FMF_PtrCurrent( iljk );

  for (ii = 0; ii < out->nCell; ii++) {
    pdetF = FMF_PtrCell( detF, ii );
    pd = FMF_PtrCell( out, ii );
    pbulk = FMF_PtrCellX1( mat, ii );

    geme_mulT2ST2S_T4S_ikjl( ikjl, traceVec, traceVec );
    geme_mulT2ST2S_T4S_iljk( iljk, traceVec, traceVec );

    for (iqp = 0; iqp < nQP; iqp++) {
      cbulk21 = pbulk[iqp] * pdetF[iqp] * (pdetF[iqp] - 1.0);
      cbulk22 = pbulk[iqp] * pdetF[iqp] * pdetF[iqp];
      for (ir = 0; ir < sym; ir++) {
        for (ic = 0; ic < sym; ic++) {
          pd[sym*ir+ic]
            = (cbulk21 + cbulk22) * ptrace[ir] * ptrace[ic]
            - cbulk21 * (p_ikjl[sym*ir+ic]
                        + p_iljk[sym*ir+ic]);
        }
      }
      pd += sym * sym;
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &ikjl );
  fmf_freeDestroy( &iljk );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_tl_he_tan_mod_bulk_active"
int32 dq_tl_he_tan_mod_bulk_active( FMField *out, FMField *mat,
                                    FMField *detF, FMField *vecInvCS )
{
  int32 ii, nQP, ir, ic, iqp, sym, ret = RET_OK;
  float64 cmat;
  float64 *pd;
  float64 *pmat;
  float64 *pinvC, *pdetF, *pinvC2_ikjl, *pinvC2_iljk;
  FMField *invC2_ikjl = 0, *invC2_iljk = 0;

  sym = out->nRow;
  nQP = out->nLev;

  fmf_createAlloc( &invC2_ikjl, 1, nQP, sym, sym );
  fmf_createAlloc( &invC2_iljk, 1, nQP, sym, sym );

  pinvC2_ikjl = FMF_PtrCurrent( invC2_ikjl );
  pinvC2_iljk = FMF_PtrCurrent( invC2_iljk );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( vecInvCS, ii );
    pdetF = FMF_PtrCell( detF, ii );
    pinvC = FMF_PtrCell( vecInvCS, ii );
    pd = FMF_PtrCell( out, ii );
    pmat = FMF_PtrCellX1( mat, ii );

    geme_mulT2ST2S_T4S_ikjl( invC2_ikjl, vecInvCS, vecInvCS );
    geme_mulT2ST2S_T4S_iljk( invC2_iljk, vecInvCS, vecInvCS );

    for (iqp = 0; iqp < nQP; iqp++) {
      cmat = pmat[iqp] * pdetF[iqp];
      for (ir = 0; ir < sym; ir++) {
        for (ic = 0; ic < sym; ic++) {
          pd[sym*ir+ic]
            = cmat * pinvC[sym*iqp+ir] * pinvC[sym*iqp+ic]
            - cmat * (pinvC2_ikjl[sym*(sym*iqp+ir)+ic]
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
#define __FUNC__ "dq_tl_he_tan_mod_neohook"
int32 dq_tl_he_tan_mod_neohook( FMField *out, FMField *mat,
                                FMField *detF, FMField *trC, FMField *vecInvCS )
{
  int32 ii, nQP, ir, ic, iqp, sym, ret = RET_OK;
  float64 cc, c1, c2, c3, detF23;
  float64 *pd;
  float64 *pmu;
  float64 *pinvC, *ptrC, *pdetF, *pinvC2_ikjl, *pinvC2_iljk, *ptrace;
  FMField *invC2_ikjl = 0, *invC2_iljk = 0;

  sym = out->nRow;
  nQP = out->nLev;
  ptrace = get_trace( sym );

  fmf_createAlloc( &invC2_ikjl, 1, nQP, sym, sym );
  fmf_createAlloc( &invC2_iljk, 1, nQP, sym, sym );

  pinvC2_ikjl = FMF_PtrCurrent( invC2_ikjl );
  pinvC2_iljk = FMF_PtrCurrent( invC2_iljk );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( vecInvCS, ii );
    pdetF = FMF_PtrCell( detF, ii );
    ptrC = FMF_PtrCell( trC, ii );
    pinvC = FMF_PtrCell( vecInvCS, ii );
    pd = FMF_PtrCell( out, ii );
    pmu = FMF_PtrCellX1( mat, ii );

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
            + c2 * ((pinvC[sym*iqp+ir] * ptrace[ic])
                  + (pinvC[sym*iqp+ic] * ptrace[ir]))
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
#define __FUNC__ "dq_ul_he_tan_mod_neohook"
int32 dq_ul_he_tan_mod_neohook( FMField *out, FMField *mat,
                                FMField *detF, FMField *trB, FMField *vecBS )
{
  int32 ii, nQP, ir, ic, iqp, sym, ret = RET_OK;
  float64 cc, c1, c2, c3, detF23;
  float64 *pd;
  float64 *pmu;
  float64 *pBG, *ptrB, *pdetF, *pd_ikjl, *pd_iljk, *ptrace;
  FMField *d_ikjl = 0, *d_iljk = 0;
  FMField traceVec[1];

  sym = out->nRow;
  nQP = out->nLev;
  ptrace = get_trace( sym );

  fmf_createAlloc( &d_ikjl, 1, 1, sym, sym );
  fmf_createAlloc( &d_iljk, 1, 1, sym, sym );
  traceVec->nAlloc = -1;
  fmf_pretend( traceVec, 1, 1, sym, 1, ptrace );

  pd_ikjl = FMF_PtrCurrent( d_ikjl );
  pd_iljk = FMF_PtrCurrent( d_iljk );

  geme_mulT2ST2S_T4S_ikjl( d_ikjl, traceVec, traceVec );
  geme_mulT2ST2S_T4S_iljk( d_iljk, traceVec, traceVec );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( vecBS, ii );
    pdetF = FMF_PtrCell( detF, ii );
    ptrB = FMF_PtrCell( trB, ii );
    pBG = FMF_PtrCell( vecBS, ii );
    pd = FMF_PtrCell( out, ii );
    pmu = FMF_PtrCellX1( mat, ii );

    for (iqp = 0; iqp < nQP; iqp++) {
      detF23 = exp( -2.0/3.0 * log( pdetF[iqp] ) );
      cc = pmu[iqp] * detF23;
      c1 = 2.0/9.0 * cc * ptrB[iqp];
      c2 = - 2.0/3.0 * cc;
      c3 = cc / 3.0 * ptrB[iqp];
      for (ir = 0; ir < sym; ir++) {
        for (ic = 0; ic < sym; ic++) {
          pd[sym*ir+ic]
            = c1 * ptrace[ir] * ptrace[ic]
            + c2 * ((pBG[sym*iqp+ir] * ptrace[ic])
                  + (pBG[sym*iqp+ic] * ptrace[ir]))
            + c3 * (pd_ikjl[sym*ir+ic]
                    + pd_iljk[sym*ir+ic]);
        }
      }
      pd += sym * sym;
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &d_ikjl );
  fmf_freeDestroy( &d_iljk );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_tl_he_tan_mod_mooney_rivlin"
int32 dq_tl_he_tan_mod_mooney_rivlin( FMField *out, FMField *mat,
                                      FMField *detF, FMField *trC,
                                      FMField *vecInvCS, FMField *vecCS,
                                      FMField *in2C )
{
  int32 ii, nQP, ir, ic, iqp, sym, dim, ret = RET_OK;
  float64 c1, c2, c3, c4, c5, c6, detF23, detF43;
  float64 *pd;
  float64 *pkappa, *ptrace;
  float64 *pinvC, *ptrC, *pdetF, *pC, *pin2C, *pinvC2_ikjl, *pinvC2_iljk;
  FMField *invC2_ikjl = 0, *invC2_iljk = 0;

  sym = out->nRow;
  nQP = out->nLev;
  dim = sym2dim( sym );
  ptrace = get_trace( sym );

  fmf_createAlloc( &invC2_ikjl, 1, nQP, sym, sym );
  fmf_createAlloc( &invC2_iljk, 1, nQP, sym, sym );

  pinvC2_ikjl = FMF_PtrCurrent( invC2_ikjl );
  pinvC2_iljk = FMF_PtrCurrent( invC2_iljk );

  for (ii = 0; ii < out->nCell; ii++) {
    pdetF = FMF_PtrCell( detF, ii );
    ptrC = FMF_PtrCell( trC, ii );
    pinvC = FMF_PtrCell( vecInvCS, ii );
    pin2C = FMF_PtrCell( in2C, ii );
    pC = FMF_PtrCell( vecCS, ii );
    FMF_SetCell( vecInvCS, ii );
    pd = FMF_PtrCell( out, ii );
    pkappa = FMF_PtrCellX1( mat, ii );

    geme_mulT2ST2S_T4S_ikjl( invC2_ikjl, vecInvCS, vecInvCS );
    geme_mulT2ST2S_T4S_iljk( invC2_iljk, vecInvCS, vecInvCS );

    for (iqp = 0; iqp < nQP; iqp++) {
      detF23 = exp( -2.0/3.0 * log( pdetF[iqp] ) );
      detF43 = detF23 * detF23;
      c1 = 8.0/9.0 * pkappa[iqp] * detF43 * pin2C[iqp];
      c2 = - 4.0/3.0 * pkappa[iqp] * detF43 * ptrC[iqp];
      c3 = 2.0/3.0 * pkappa[iqp] * detF43 * pin2C[iqp];
      c4 = 2.0 * pkappa[iqp] * detF43;
      c5 = - pkappa[iqp] * detF43;
      c6 = 4.0/3.0 * pkappa[iqp] * detF43;
      for (ir = 0; ir < sym; ir++) {
        for (ic = 0; ic < sym; ic++) {
          pd[sym*ir+ic]
            = c1 * pinvC[sym*iqp+ir] * pinvC[sym*iqp+ic]
            + c2 * ((pinvC[sym*iqp+ir] * ptrace[ic])
                  + (pinvC[sym*iqp+ic] * ptrace[ir]))
            + c3 * (pinvC2_ikjl[sym*(sym*iqp+ir)+ic]
                  + pinvC2_iljk[sym*(sym*iqp+ir)+ic])
            + c4 * (ptrace[ir] * ptrace[ic])
            + c6 * ((pC[sym*iqp+ir] * pinvC[sym*iqp+ic])
                  + (pC[sym*iqp+ic] * pinvC[sym*iqp+ir]))
          ;
        }
      }
      for (ir = 0; ir < dim; ir++) {
        pd[sym*ir+ir] += 2.0 * c5;
      }
      for (ir = dim; ir < sym; ir++) {
        pd[sym*ir+ir] += c5;
      }
      pd += sym*sym;
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &invC2_ikjl );
  fmf_freeDestroy( &invC2_iljk );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_ul_he_tan_mod_mooney_rivlin"
int32 dq_ul_he_tan_mod_mooney_rivlin( FMField *out, FMField *mat,
                                      FMField *detF, FMField *trB,
                                      FMField *vecBS, FMField *in2B )
{
  int32 ii, nQP, ir, ic, iqp, sym, ret = RET_OK;
  float64 F1, F2, F3, F4, F5, F6, detF23, detF43;
  float64 *pd, *pkappa, *ptrace;
  float64 *ptrB, *pdetF, *pB, *pin2B, *pB_ikjl, *pB_iljk, *pd_ikjl, *pd_iljk, *pBB;
  FMField *B_ikjl = 0, *B_iljk = 0, *d_ikjl = 0, *d_iljk = 0, *vecBB = 0;
  FMField traceVec[1];

  sym = out->nRow;
  nQP = out->nLev;
  ptrace = get_trace( sym );

  fmf_createAlloc( &B_ikjl, 1, nQP, sym, sym );
  fmf_createAlloc( &B_iljk, 1, nQP, sym, sym );
  fmf_createAlloc( &d_ikjl, 1, 1, sym, sym );
  fmf_createAlloc( &d_iljk, 1, 1, sym, sym );
  traceVec->nAlloc = -1;
  fmf_pretend( traceVec, 1, 1, sym, 1, ptrace );
  fmf_createAlloc( &vecBB, 1, nQP, sym, 1 );

  pB_ikjl = FMF_PtrCurrent( B_ikjl );
  pB_iljk = FMF_PtrCurrent( B_iljk );
  pd_ikjl = FMF_PtrCurrent( d_ikjl );
  pd_iljk = FMF_PtrCurrent( d_iljk );

  geme_mulT2ST2S_T4S_ikjl( d_ikjl, traceVec, traceVec );
  geme_mulT2ST2S_T4S_iljk( d_iljk, traceVec, traceVec );

  for (ii = 0; ii < out->nCell; ii++) {
    pdetF = FMF_PtrCell( detF, ii );
    ptrB = FMF_PtrCell( trB, ii );
    pin2B = FMF_PtrCell( in2B, ii );
    pB = FMF_PtrCell( vecBS, ii );
    FMF_SetCell( vecBS, ii );
    pd = FMF_PtrCell( out, ii );
    pkappa = FMF_PtrCellX1( mat, ii );
    pBB = vecBB->val0;

    geme_mulT2ST2S_T4S_ikjl( B_ikjl, vecBS, vecBS );
    geme_mulT2ST2S_T4S_iljk( B_iljk, vecBS, vecBS );
    geme_mulT2S_AA( vecBB, vecBS );

    /* Crisfield II., (13.40), (13.41), (13.100) */
    for (iqp = 0; iqp < nQP; iqp++) {
      detF23 = exp( -2.0/3.0 * log( pdetF[iqp] ) );
      detF43 = detF23 * detF23;
      F1 = 16.0/9.0 * pkappa[iqp] * detF43 * pin2B[iqp];
      F2 = -8.0/3.0 * pkappa[iqp] * detF43 * ptrB[iqp];
      F3 = 4.0/3.0 * pkappa[iqp] * detF43 * pin2B[iqp];
      F4 = 4.0 * pkappa[iqp] * detF43;
      F5 = -2.0 * pkappa[iqp] * detF43;
      F6 = 8.0/3.0 * pkappa[iqp] * detF43;
      for (ir = 0; ir < sym; ir++) {
        for (ic = 0; ic < sym; ic++) {
          pd[sym*ir+ic]
            = F1 * ptrace[ir] * ptrace[ic]
            + F2 * (pB[ir] * ptrace[ic] + ptrace[ir] * pB[ic])
            + F3 * (pd_ikjl[sym*ir+ic] + pd_iljk[sym*ir+ic])
            + F4 * pB[ir] * pB[ic]
            + F5 * (pB_ikjl[sym*(sym*iqp+ir)+ic] + pB_iljk[sym*(sym*iqp+ir)+ic])
            + F6 * (pBB[ir] * ptrace[ic] + pBB[ic] * ptrace[ir]);
        }
      }
      pd += sym*sym;
      pB += sym;
      pBB += sym;
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &B_ikjl );
  fmf_freeDestroy( &B_iljk );
  fmf_freeDestroy( &d_ikjl );
  fmf_freeDestroy( &d_iljk );
  fmf_freeDestroy( &vecBB );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_he_stress_bulk_pressure"
int32 dq_he_stress_bulk_pressure( FMField *out, FMField *pressure_qp,
                                  FMField *detF, FMField *vecInvCS,
                                  int32 mode_ul )
{
  int32 ii, iqp, ir, sym, nQP, ret = RET_OK;
  float64 aux;
  float64 *pp = 0, *pstress = 0, *pCI = 0, *pdetF = 0, *ptrace;

  nQP = detF->nLev;
  sym = out->nRow;
  ptrace = get_trace( sym );

  for (ii = 0; ii < out->nCell; ii++) {
    pstress = FMF_PtrCell( out, ii );
    pp = FMF_PtrCell( pressure_qp, ii );
    pdetF = FMF_PtrCell( detF, ii );

    if (mode_ul) {
      for (iqp = 0; iqp < nQP; iqp++) {
        aux = - pp[iqp] * pdetF[iqp];
        for (ir = 0; ir < sym; ir++) {
          pstress[ir] = aux * ptrace[ir];
        }
        pstress += sym;
      }
    }
    else {
      pCI = FMF_PtrCell( vecInvCS, ii );
      for (iqp = 0; iqp < nQP; iqp++) {
        aux = - pp[iqp] * pdetF[iqp];
        for (ir = 0; ir < sym; ir++) {
          pstress[ir] = aux * pCI[ir];
        }
        pstress += sym;
        pCI += sym;
      }
    } /* if mode_ul */
    ERR_CheckGo( ret );
  }

 end_label:
  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_tl_stress_bulk_pressure"
int32 dq_tl_stress_bulk_pressure( FMField *out, FMField *pressure_qp,
                                  FMField *detF, FMField *vecInvCS )
{
  return( dq_he_stress_bulk_pressure( out, pressure_qp, detF, vecInvCS, 0 ) );
}

#undef __FUNC__
#define __FUNC__ "dq_ul_stress_bulk_pressure"
int32 dq_ul_stress_bulk_pressure( FMField *out, FMField *pressure_qp,
                                  FMField *detF )
{
  return( dq_he_stress_bulk_pressure( out, pressure_qp, detF, NULL, 1 ) );
}

#undef __FUNC__
#define __FUNC__ "dq_tl_tan_mod_bulk_pressure_u"
int32 dq_tl_tan_mod_bulk_pressure_u( FMField *out, FMField *pressure_qp,
                                     FMField *detF, FMField *vecInvCS )
{
  int32 ii, nQP, ir, ic, iqp, sym, ret = RET_OK;
  float64 cc;
  float64 *pd, *ppress;
  float64 *pinvC, *pdetF, *pinvC2_ikjl, *pinvC2_iljk;
  FMField *invC2_ikjl = 0, *invC2_iljk = 0;

  sym = out->nRow;
  nQP = out->nLev;

  fmf_createAlloc( &invC2_ikjl, 1, nQP, sym, sym );
  fmf_createAlloc( &invC2_iljk, 1, nQP, sym, sym );

  pinvC2_ikjl = FMF_PtrCurrent( invC2_ikjl );
  pinvC2_iljk = FMF_PtrCurrent( invC2_iljk );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( vecInvCS, ii );
    pdetF = FMF_PtrCell( detF, ii );
    pinvC = FMF_PtrCell( vecInvCS, ii );
    pd = FMF_PtrCell( out, ii );
    ppress = FMF_PtrCell( pressure_qp, ii );

    geme_mulT2ST2S_T4S_ikjl( invC2_ikjl, vecInvCS, vecInvCS );
    geme_mulT2ST2S_T4S_iljk( invC2_iljk, vecInvCS, vecInvCS );

    for (iqp = 0; iqp < nQP; iqp++) {
      cc = pdetF[iqp] * ppress[iqp];
      for (ir = 0; ir < sym; ir++) {
        for (ic = 0; ic < sym; ic++) {
          pd[sym*ir+ic] =
            - cc * pinvC[sym*iqp+ir] * pinvC[sym*iqp+ic]
            + cc * (pinvC2_ikjl[sym*(sym*iqp+ir)+ic]
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
#define __FUNC__ "dq_ul_tan_mod_bulk_pressure_u"
int32 dq_ul_tan_mod_bulk_pressure_u( FMField *out, FMField *pressure_qp,
                                     FMField *detF )
{
  int32 ii, nQP, ir, ic, iqp, sym, ret = RET_OK;
  float64 cc;
  float64 *pd, *ppress;
  float64 *pdetF, *p_ikjl, *p_iljk, *ptrace;
  FMField *ikjl = 0, *iljk = 0;
  FMField traceVec[1];

  sym = out->nRow;
  nQP = out->nLev;
  ptrace = get_trace( sym );

  fmf_createAlloc( &ikjl, 1, 1, sym, sym );
  fmf_createAlloc( &iljk, 1, 1, sym, sym );
  traceVec->nAlloc = -1;
  fmf_pretend( traceVec, 1, 1, sym, 1, ptrace );

  p_ikjl = FMF_PtrCurrent( ikjl );
  p_iljk = FMF_PtrCurrent( iljk );

  for (ii = 0; ii < out->nCell; ii++) {
    pdetF = FMF_PtrCell( detF, ii );
    pd = FMF_PtrCell( out, ii );
    ppress = FMF_PtrCell( pressure_qp, ii );

    geme_mulT2ST2S_T4S_ikjl( ikjl, traceVec, traceVec );
    geme_mulT2ST2S_T4S_iljk( iljk, traceVec, traceVec );

    for (iqp = 0; iqp < nQP; iqp++) {
      cc = pdetF[iqp] * ppress[iqp];
      for (ir = 0; ir < sym; ir++) {
        for (ic = 0; ic < sym; ic++) {
          pd[sym*ir+ic] =
            - cc * ptrace[ir] * ptrace[ic]
            + cc * (p_ikjl[sym*ir+ic]
                  + p_iljk[sym*ir+ic]);
        }
      }
      pd += sym * sym;
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &ikjl );
  fmf_freeDestroy( &iljk );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_tl_volume"
int32 dw_tl_volume( FMField *out, FMField *mtxF,
                    FMField *vecInvCS, FMField *detF,
                    Mapping *vgs, Mapping *vgv,
                    int32 transpose, int32 mode )
{
  int32 ii, nQP, nEP, nRow, sym, ret = RET_OK;
  FMField *aux = 0, *jcitb = 0, *fjcitb = 0;

  if (mode == 0) {
    fmf_createAlloc( &aux, 1, vgs->bf->nLev, vgs->bf->nRow, vgs->bf->nCol );

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( vgv->det, ii );
      FMF_SetCellX1( vgs->bf, ii );

      FMF_SetCell( out, ii );
      FMF_SetCell( detF, ii );

      fmf_copy( aux, vgs->bf );
      fmf_mul( aux, detF->val );

      fmf_sumLevelsTMulF( out, aux, vgv->det->val );
      ERR_CheckGo( ret );
    }

  } else if ((mode == 1) || (mode == -1)) {
    nQP = vgv->bfGM->nLev;
    sym = vecInvCS->nRow;
    nEP = vgs->bf->nCol;
    nRow = vgv->bfGM->nRow * vgv->bfGM->nCol; /* dim * nEP */

    fmf_createAlloc( &aux, 1, nQP, sym, nRow );
    fmf_createAlloc( &jcitb, 1, nQP, 1, nRow );
    fmf_createAlloc( &fjcitb, 1, nQP, nEP, nRow );

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( vgv->bfGM, ii );
      FMF_SetCell( vgv->det, ii );
      FMF_SetCellX1( vgs->bf, ii );

      FMF_SetCell( out, ii );
      FMF_SetCell( mtxF, ii );
      FMF_SetCell( detF, ii );
      FMF_SetCell( vecInvCS, ii );

      form_tlcc_buildOpB_VS3( aux, mtxF, vgv->bfGM );
      fmf_mulATB_nn( jcitb, vecInvCS, aux );
      fmf_mul( jcitb, detF->val );
      fmf_mulATB_nn( fjcitb, vgs->bf, jcitb );

      if (transpose) {
        fmf_sumLevelsTMulF( out, fjcitb, vgv->det->val );
      } else {
        fmf_sumLevelsMulF( out, fjcitb, vgv->det->val );
      }
      fmf_mulC( out, mode );
      ERR_CheckGo( ret );
    }
  } else if (mode == 2){ // de_volume

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( vgv->det, ii );

      FMF_SetCell( out, ii );
      FMF_SetCell( detF, ii );

      fmf_sumLevelsMulF( out, detF, vgv->det->val );
      ERR_CheckGo( ret );
    }
  } else { // mode == 3, de_rel_volume

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( vgv->det, ii );
      FMF_SetCell( vgv->volume, ii );

      FMF_SetCell( out, ii );
      FMF_SetCell( detF, ii );

      fmf_sumLevelsMulF( out, detF, vgv->det->val );
      fmf_mulC( out, 1.0 / vgv->volume->val[0] );

      ERR_CheckGo( ret );
    }
  }

 end_label:
  fmf_freeDestroy( &aux );
  fmf_freeDestroy( &jcitb );
  fmf_freeDestroy( &fjcitb );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_ul_volume"
int32 dw_ul_volume( FMField *out, FMField *detF,
                    Mapping *vgs, Mapping *vgv,
                    int32 transpose, int32 mode )
{
  int32 ii, iqp, nQP, nEPu, nEPp, dim, ret = RET_OK;
  FMField *aux = 0, aux2[1];

  nQP = vgv->bfGM->nLev;
  nEPu = vgv->bfGM->nCol;
  dim = vgv->bfGM->nRow;
  nEPp = vgs->bf->nCol;

  if (mode == 0) {
    fmf_createAlloc( &aux, 1, nQP, 1, 1 );

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( detF, ii );
      FMF_SetCell( vgv->det, ii );
      FMF_SetCellX1( vgs->bf, ii );
      FMF_SetCell( out, ii );

      for (iqp = 0; iqp < nQP; iqp++) {
        aux->val[iqp] = (1.0 - 1.0 / detF->val[iqp]) * vgv->det->val[iqp];
      }
      /* \int N^T (1 - 1/J) */
      fmf_sumLevelsTMulF( out, vgs->bf, aux->val );
      ERR_CheckGo( ret );
    }

  } else if ((mode == 1) || (mode == -1)) {
    fmf_createAlloc( &aux, 1, nQP, nEPp, dim * nEPu );
    aux2->nAlloc = -1;
    fmf_pretend( aux2, 1, nQP, 1, dim * nEPu, NULL );

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( vgv->bfGM, ii );
      FMF_SetCell( vgv->det, ii );
      FMF_SetCellX1( vgs->bf, ii );
      FMF_SetCell( out, ii );
      aux2->val = vgv->bfGM->val;

      /* \int N^T \delta B*/
      fmf_mulATB_nn( aux, vgs->bf, aux2 );

      if (transpose) {
        fmf_sumLevelsTMulF( out, aux, vgv->det->val );
      } else {
        fmf_sumLevelsMulF( out, aux, vgv->det->val );
      }
      if (mode == -1) {
        fmf_mulC(out, -1.0);
      }

      ERR_CheckGo( ret );
    }
  } else if (mode == 2){ // de_volume

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( vgv->det, ii );
      FMF_SetCell( detF, ii );
      FMF_SetCell( out, ii );

      fmf_sumLevelsMulF( out, detF, vgv->det->val );
      ERR_CheckGo( ret );
    }
  } else { // mode == 3, de_rel_volume

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( vgv->det, ii );
      FMF_SetCell( vgv->volume, ii );
      FMF_SetCell( detF, ii );
      FMF_SetCell( out, ii );

      fmf_sumLevelsMulF( out, detF, vgv->det->val );
      fmf_mulC( out, 1.0 / vgv->volume->val[0] );

      ERR_CheckGo( ret );
    }
  }

 end_label:
  fmf_freeDestroy( &aux );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_tl_diffusion"
int32 dw_tl_diffusion( FMField *out, FMField *pressure_grad,
                       FMField *mtxD, FMField *ref_porosity,
                       FMField *mtxF, FMField *detF,
                       Mapping *vg, int32 mode )
{
  int32 ii, iqp, dim, nEP, nQP, ret = RET_OK;
  float64 val;
  FMField *gtd = 0, *gtdg = 0, *dgp = 0, *gtdgp = 0;
  FMField *coef = 0, *perm = 0, *mtxFI = 0, *aux = 0, *mtxK = 0, *w_qp = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;
  nEP = vg->bfGM->nCol;

  fmf_createAlloc( &coef, 1, nQP, 1, 1 );
  fmf_createAlloc( &perm, 1, nQP, dim, dim );
  fmf_createAlloc( &mtxFI, 1, nQP, dim, dim );
  fmf_createAlloc( &aux, 1, nQP, dim, dim );

  if (mode < 2) {
    fmf_createAlloc( &mtxK, 1, nQP, dim, dim );

    if (mode == 1) {
      fmf_createAlloc( &gtd, 1, nQP, nEP, dim );
      fmf_createAlloc( &gtdg, 1, nQP, nEP, nEP );
    } else {
      fmf_createAlloc( &dgp, 1, nQP, dim, 1 );
      fmf_createAlloc( &gtdgp, 1, nQP, nEP, 1 );
    }
  } else {
    fmf_createAlloc( &w_qp, 1, nQP, dim, 1 );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCell( mtxF, ii );
    FMF_SetCell( detF, ii );
    FMF_SetCellX1( mtxD, ii );
    FMF_SetCellX1( ref_porosity, ii );

    // max(0, (1 + (J - 1) / N_f))^2
    for (iqp = 0; iqp < nQP; iqp++) {
      val = 1.0 + ((detF->val[iqp] - 1.0) / ref_porosity->val[iqp]);
      if (val <= 0.0) val = 0.0;
      coef->val[iqp] = val * val;
    }
    // Actual permeability.
    fmf_mulAF( perm, mtxD, coef->val );

    geme_invert3x3( mtxFI, mtxF );

    if (mode < 2) {
      // Regular assembling.

      // Transformed permeability.
      fmf_mulAB_nn( aux, mtxFI, perm );
      fmf_mulABT_nn( mtxK, aux, mtxFI );
      fmf_mul( mtxK, detF->val );

      if (mode == 1) {
        fmf_mulATB_nn( gtd, vg->bfGM, mtxK );
        fmf_mulAB_nn( gtdg, gtd, vg->bfGM );
        fmf_sumLevelsMulF( out, gtdg, vg->det->val );
      } else {
        FMF_SetCell( pressure_grad, ii );
        fmf_mulAB_nn( dgp, mtxK, pressure_grad );
        fmf_mulATB_nn( gtdgp, vg->bfGM, dgp );
        fmf_sumLevelsMulF( out, gtdgp, vg->det->val );
      }
    } else {
      // Diffusion velocity averaged in elements.
      FMF_SetCell( vg->volume, ii );
      FMF_SetCell( pressure_grad, ii );

      fmf_mulABT_nn( aux, perm, mtxFI );
      fmf_mulAB_nn( w_qp, aux, pressure_grad );
      fmf_sumLevelsMulF( out, w_qp, vg->det->val );
      fmf_mulC( out, -1.0 / vg->volume->val[0] );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &coef );
  fmf_freeDestroy( &perm );
  fmf_freeDestroy( &mtxFI );
  fmf_freeDestroy( &aux );

  if (mode < 2) {
    fmf_freeDestroy( &mtxK );

    if (mode == 1) {
      fmf_freeDestroy( &gtd );
      fmf_freeDestroy( &gtdg );
    } else {
      fmf_freeDestroy( &dgp );
      fmf_freeDestroy( &gtdgp );
    }
  } else {
      fmf_freeDestroy( &w_qp );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_tl_surface_flux"
int32 d_tl_surface_flux( FMField *out, FMField *pressure_grad,
                         FMField *mtxD, FMField *ref_porosity,
                         FMField *mtxFI, FMField *detF,
                         Mapping *sg, int32 mode )
{
  int32 ii, iqp, dim, nQP, ret = RET_OK;
  float64 val;
  FMField *coef = 0, *perm = 0, *aux = 0, *mtxK = 0, *kgp = 0, *ntkgp = 0;

  nQP = sg->normal->nLev;
  dim = sg->normal->nRow;

  fmf_createAlloc( &coef, 1, nQP, 1, 1 );
  fmf_createAlloc( &perm, 1, nQP, dim, dim );
  fmf_createAlloc( &aux, 1, nQP, dim, dim );
  fmf_createAlloc( &mtxK, 1, nQP, dim, dim );
  fmf_createAlloc( &kgp, 1, nQP, dim, 1 );
  fmf_createAlloc( &ntkgp, 1, nQP, 1, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( pressure_grad, ii );
    FMF_SetCellX1( mtxD, ii );
    FMF_SetCell( mtxFI, ii );
    FMF_SetCell( detF, ii );
    FMF_SetCellX1( ref_porosity, ii );
    FMF_SetCell( sg->normal, ii );
    FMF_SetCell( sg->det, ii );

    // max(0, (1 + (J - 1) / N_f))^2
    for (iqp = 0; iqp < nQP; iqp++) {
      val = 1.0 + ((detF->val[iqp] - 1.0) / ref_porosity->val[iqp]);
      if (val <= 0.0) val = 0.0;
      coef->val[iqp] = val * val;
    }
    // Actual permeability.
    fmf_mulAF( perm, mtxD, coef->val );

    // Transformed permeability.
    fmf_mulAB_nn( aux, mtxFI, perm );
    fmf_mulABT_nn( mtxK, aux, mtxFI );
    fmf_mul( mtxK, detF->val );

    // Flux.
    fmf_mulAB_nn( kgp, mtxK, pressure_grad );
    fmf_mulATB_nn( ntkgp, sg->normal, kgp );

    fmf_sumLevelsMulF( out, ntkgp, sg->det->val );
    if (mode == 1) {
      FMF_SetCell( sg->volume, ii );
      fmf_mulC( out, 1.0 / sg->volume->val[0] );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &coef );
  fmf_freeDestroy( &perm );
  fmf_freeDestroy( &aux );
  fmf_freeDestroy( &mtxK );
  fmf_freeDestroy( &kgp );
  fmf_freeDestroy( &ntkgp );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_tl_surface_traction"
int32 dw_tl_surface_traction( FMField *out, FMField *traction,
                              FMField *detF, FMField *mtxFI,
                              FMField *bf, Mapping *sg,
                              int32 *fis, int32 nFa, int32 nFP,
                              int32 mode )
{
  int32 ii, iqp, idr, idc, iep, ifa, nEP, nQP, dim, ret = RET_OK;
  float64 *pn2, *pbfBGS, *paux;
  FMField *n2 = 0, *stn2 = 0, *trq = 0;
  FMField *trdq = 0, *aux = 0, *staux = 0, *bfBGS = 0;

  dim = mtxFI->nRow;
  nQP = mtxFI->nLev;
  nEP = sg->bfGM->nCol;

/*    output( "%d %d %d\n", dim, nQP, nEP ); */

  fmf_createAlloc( &n2, 1, nQP, dim, 1 );
  if (mode == 0) {
    fmf_createAlloc( &stn2, 1, nQP, dim, 1 );
    fmf_createAlloc( &trq, 1, nQP, dim * nEP, 1 );
  } else {
    fmf_createAlloc( &bfBGS, 1, nQP, dim, nEP );
    fmf_createAlloc( &aux, 1, nQP, dim, dim * nEP );
    fmf_createAlloc( &staux, 1, nQP, dim, dim * nEP );
    fmf_createAlloc( &trdq, 1, nQP, dim * nEP, dim * nEP );
  }

  for (ii = 0; ii < out->nCell; ii++) {
    ifa = fis[ii*nFP+1];

    FMF_SetCell( out, ii );
    FMF_SetCellX1( traction, ii );
    FMF_SetCell( detF, ii );
    FMF_SetCell( mtxFI, ii );
    FMF_SetCell( sg->normal, ii );
    FMF_SetCell( sg->det, ii );
    FMF_SetCell( bf, ifa );

    fmf_mulATB_nn( n2, mtxFI, sg->normal );

    if (mode == 0) {
      fmf_mulATB_nn( stn2, traction, n2 );
      fmf_mul( stn2, detF->val );
      bf_actt( trq, bf, stn2 );
      fmf_sumLevelsMulF( out, trq, sg->det->val );

    } else {
      FMF_SetCell( sg->bfGM, ii );

      fmf_mulATB_nn( bfBGS, mtxFI, sg->bfGM );

      for (iqp = 0; iqp < nQP; iqp++) {
        pn2 = FMF_PtrLevel( n2, iqp );
        pbfBGS = FMF_PtrLevel( bfBGS, iqp );

        for (idr = 0; idr < dim; idr++) {
          paux = FMF_PtrRowOfLevel( aux, iqp, idr );

          for (idc = 0; idc < dim; idc++) {
            for (iep = 0; iep < nEP; iep++) {
              paux[nEP*idc+iep]
                = (pn2[idr] * pbfBGS[nEP*idc+iep]
                   - pn2[idc] * pbfBGS[nEP*idr+iep]) * detF->val[iqp];
            }
          }
        }
      }
      fmf_mulATB_nn( staux, traction, aux );
      bf_actt( trdq, bf, staux );
      fmf_sumLevelsMulF( out, trdq, sg->det->val );
    }

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &n2 );
  if (mode == 0) {
    fmf_freeDestroy( &stn2 );
    fmf_freeDestroy( &trq );
  } else {
    fmf_freeDestroy( &bfBGS );
    fmf_freeDestroy( &aux );
    fmf_freeDestroy( &staux );
    fmf_freeDestroy( &trdq );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_tl_volume_surface"
int32 d_tl_volume_surface( FMField *out, FMField *coors,
                           FMField *detF, FMField *mtxFI,
                           FMField *bf, Mapping *sg,
                           int32 *conn, int32 nFa, int32 nFP )
{
  int32 ii, nQP, dim, ret = RET_OK;
  float64 val;
  FMField *aux = 0, *coors_qp = 0, *n2 = 0, *aux2 = 0;

  dim = mtxFI->nRow;
  nQP = mtxFI->nLev;

  val = 1.0 / dim;

  fmf_createAlloc( &aux, 1, 1, nFP, dim );
  fmf_createAlloc( &coors_qp, 1, nQP, 1, dim );
  fmf_createAlloc( &n2, 1, nQP, dim, 1 );
  fmf_createAlloc( &aux2, 1, nQP, 1, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( detF, ii );
    FMF_SetCell( mtxFI, ii );
    FMF_SetCell( sg->normal, ii );
    FMF_SetCell( sg->det, ii );
    FMF_SetCellX1( sg->bf, ii );

    ele_extractNodalValuesNBN( aux, coors, conn + nFP * ii );
    fmf_mulAB_n1( coors_qp, sg->bf, aux );

    fmf_mulATB_nn( n2, mtxFI, sg->normal );
    fmf_mulAB_nn( aux2, coors_qp, n2 );
    fmf_mul( aux2, detF->val );
    fmf_sumLevelsMulF( out, aux2, sg->det->val );
    fmf_mulC( out, val );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &aux );
  fmf_freeDestroy( &coors_qp );
  fmf_freeDestroy( &n2 );
  fmf_freeDestroy( &aux2 );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_def_grad"
int32 dq_def_grad( FMField *out, FMField *state, Mapping *vg,
                   int32 *conn, int32 nEl, int32 nEP, int32 mode )
{
  int32 ii, id, iqp, nQP, dim, ret = RET_OK;
  FMField *st = 0, *mtxF = 0;

  state->val = FMF_PtrFirst( state );

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  fmf_createAlloc( &st, 1, 1, nEP, dim );
  if (mode == 1) {
    fmf_createAlloc( &mtxF, 1, nQP, dim, dim );
  }

  for (ii = 0; ii < nEl; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( vg->bfGM, ii );

    // Deformation gradient.
    ele_extractNodalValuesNBN( st, state, conn + nEP * ii );
    if (mode == 1) {
      fmf_mulATBT_1n( mtxF, st, vg->bfGM );

      for (iqp = 0; iqp < nQP; iqp++) {
        for (id = 0; id < dim; id++) {
          mtxF->val[dim*(dim*iqp+id)+id] += 1.0;
        }
      }

      // Determinant of deformation gradient.
      geme_det3x3( out->val, mtxF );

    } else {
      fmf_mulATBT_1n( out, st, vg->bfGM );

      for (iqp = 0; iqp < nQP; iqp++) {
        for (id = 0; id < dim; id++) {
          out->val[dim*(dim*iqp+id)+id] += 1.0;
        }
      }
    }

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &st );
  fmf_freeDestroy( &mtxF );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "he_residuum_from_mtx"
int32 he_residuum_from_mtx(FMField *out, FMField *mtxD,
                           FMField *state,
                           int32 *conn, int32 nEl, int32 nEP,
                           int32 *elList, int32 elList_nRow)
{
  int32 ii, iel, ret = RET_OK, dim;
  FMField *st = 0;
  FMField pst[1];

  dim = mtxD->nRow / nEP;

  fmf_createAlloc( &st, 1, 1, dim, nEP );
  pst->nAlloc = -1;
  fmf_pretend( pst, 1, 1, nEP * dim, 1, st->val );

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCellX1( mtxD, ii );
    ele_extractNodalValuesDBD( st, state, conn + nEP * iel );

    fmf_mulAB_nn( out, mtxD, pst );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &st );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "he_eval_from_mtx"
int32 he_eval_from_mtx(FMField *out, FMField *mtxD,
                       FMField *stateV, FMField *stateU,
                       int32 *conn, int32 nEl, int32 nEP,
                       int32 *elList, int32 elList_nRow)
{
  int32 ii, iel, ret = RET_OK, dim;
  FMField *st = 0, *aux = 0;
  FMField pst[1];

  dim = mtxD->nRow / nEP;

  fmf_createAlloc( &st, 1, 1, dim, nEP );
  pst->nAlloc = -1;
  fmf_pretend( pst, 1, 1, nEP * dim, 1, st->val );
  fmf_createAlloc( &aux, 1, 1, nEP * dim, 1 );

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCellX1( mtxD, ii );

    ele_extractNodalValuesDBD( st, stateU, conn + nEP * iel );
    fmf_mulAB_nn( aux, mtxD, pst );
    ele_extractNodalValuesDBD( st, stateV, conn + nEP * iel );
    fmf_mulATB_nn( out, pst, aux );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &st );
  fmf_freeDestroy( &aux );

  return( ret );
}
