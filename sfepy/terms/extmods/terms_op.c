#include "terms_op.h"
#include "terms.h"

#undef __FUNC__
#define __FUNC__ "mulAB_integrate"
int32 mulAB_integrate(FMField *out, FMField *A, FMField *B,
                      Mapping *vg, int32 mode)
{
  int32 ii, ret = RET_OK;
  FMField *aux = 0;
  int32 (*fmul)(FMField*, FMField*, FMField*) = 0;

  if (mode == 0) {
    fmul = &fmf_mulATB_nn;
    fmf_createAlloc(&aux, 1, A->nLev, A->nCol, B->nCol);
  } else if (mode == 1) {
    fmul = &fmf_mulAB_nn;
    fmf_createAlloc(&aux, 1, A->nLev, A->nRow, B->nCol);
  } else if (mode == 2) {
    fmul = &fmf_mulABT_nn;
    fmf_createAlloc(&aux, 1, A->nLev, A->nRow, B->nRow);
  } else if (mode == 3) {
    fmul = &fmf_mulATBT_nn;
    fmf_createAlloc(&aux, 1, A->nLev, A->nCol, B->nRow);
  } else {
    errput("unknown multiplication mode!\n");
  }

  FMF_SetFirst(A);
  FMF_SetFirst(B);
  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell(out, ii);
    FMF_SetCell(vg->det, ii);
    FMF_SetCellX1(A, ii);
    FMF_SetCellX1(B, ii);
    (*fmul)(aux, A, B);
    fmf_sumLevelsMulF(out, aux, vg->det->val);

    ERR_CheckGo(ret);
  }

 end_label:
  fmf_freeDestroy(&aux);

  return(ret);
}

#undef __FUNC__
#define __FUNC__ "actBfT"
int32 actBfT(FMField *out, FMField *bf, FMField *A)
{
  int32 iqp, iel, ir, ic, ii, nEP, dim;
  float64 *pout, *pbf, *pA;

  nEP = bf->nCol;
  dim = A->nRow;

#ifdef DEBUG_FMF
  if ((out->nCol != A->nCol) || (out->nRow != (dim * nEP))
      || (out->nLev != bf->nLev) || (out->nCell != A->nCell)
      || (A->nLev != bf->nLev)) {
    errput(ErrHead "ERR_BadMatch: (%d, %d %d %d), (%d %d %d), (%d, %d %d %d)\n",
           out->nCell, out->nLev, out->nRow, out->nCol,
           bf->nLev, bf->nRow, bf->nCol,
           A->nCell, A->nLev, A->nRow, A->nCol);
  }
#endif

  fmf_fillC(out, 0.0);
  for (iel = 0; iel < out->nCell; iel++) {
    FMF_SetCell(out, iel);
    FMF_SetCellX1(A, iel);
    for (iqp = 0; iqp < bf->nLev; iqp++) {
      pbf = FMF_PtrLevel(bf, iqp);
      pout = FMF_PtrLevel(out, iqp);
      for (ii = 0; ii < nEP; ii++) {
        pA = FMF_PtrLevel(A, iqp);
        for (ir = 0; ir < dim; ir++) {
          for (ic = 0; ic < A->nCol; ic++) {
            pout[ic] = pbf[ii] * pA[ic];
          }
          pout += out->nCol;
          pA += A->nCol;
        }
      }
    } // for (iqp)
  } // for (iel)
  return(RET_OK);
}

#undef __FUNC__
#define __FUNC__ "sym2nonsym"
int32 sym2nonsym(FMField *out, FMField *A)
{
  int32 iqp, iel, ir, ic, sym, nonsym, dim, ii;
  float64 *pout, *pA, *pA0;
  int32 *map, *pmap;
  int32 nonsym_tab_2d[4] = {0, 2, 2, 1};
  int32 nonsym_tab_3d[9] = {0, 3, 4, 3, 1, 5, 4, 5, 2};

  sym = A->nRow;
  switch (sym) {
  case 3:
    dim = 2;
    nonsym = 4;
    map = nonsym_tab_2d;
    break;
  case 6:
    dim = 3;
    nonsym = 9;
    map = nonsym_tab_3d;
    break;
  default:
    errput( ErrHead "ERR_Switch\n" );
    return( RET_Fail );
  }

#ifdef DEBUG_FMF
  if ((out->nCell != A->nCell) || (out->nLev != A->nLev)
      || (out->nCol != nonsym) || (out->nRow != nonsym)) {
    errput(ErrHead "ERR_BadMatch: (%d, %d %d %d), (%d %d %d %d)\n",
           out->nCell, out->nLev, out->nRow, out->nCol,
           A->nCell, A->nLev, A->nRow, A->nCol);
  }
#endif

  if (A->nCol == 1) { /* stress mode */
    for (iel = 0; iel < out->nCell; iel++) {
      FMF_SetCell(out, iel);
      FMF_SetCellX1(A, iel);
      for (iqp = 0; iqp < out->nLev; iqp++) {
        pout = FMF_PtrLevel(out, iqp);
        pA = FMF_PtrLevel(A, iqp);
        for (ii = 0; ii < dim; ii++) {
          pmap = map;
          for (ir = 0; ir < dim; ir++) {
            for (ic = 0; ic < dim; ic++)
              pout[ic] = pA[pmap[ic]];
            pmap += dim;
            pout += dim * dim;
          } // for (ir)
          pout += dim;
        } // for (ii)
      } // for (iqp)
    } // for (iel)
  }
  else { /* tan mode */
    for (iel = 0; iel < out->nCell; iel++) {
      FMF_SetCell(out, iel);
      FMF_SetCellX1(A, iel);
      for (iqp = 0; iqp < out->nLev; iqp++) {
        pout = FMF_PtrLevel(out, iqp);
        pA0 = FMF_PtrLevel(A, iqp);
        for (ir = 0; ir < nonsym; ir++) {
          pA = pA0 + map[ir] * sym;
          for (ic = 0; ic < nonsym; ic++)
              pout[ic] = pA[map[ic]];
          pout += nonsym;
        } // for (ir)
      } // for (iqp)
    } // for (iel)
  }

  return(RET_OK);
}
