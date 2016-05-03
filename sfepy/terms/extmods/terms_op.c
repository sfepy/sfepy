#include "terms_op.h"
#include "terms.h"

#undef __FUNC__
#define __FUNC__ "mulAB_integrate"
int32 mulAB_integrate(FMField *out, FMField *A, FMField *B,
                      Mapping *vg, int32 mode)
{
  int32 ii, ret = RET_OK;
  FMField *aux = 0;
  int32 (*fmul)(FMField*, FMField*, FMField*);

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
    FMF_SetCell(A, iel);
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
