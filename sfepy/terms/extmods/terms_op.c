#include "terms_op.h"
#include "terms.h"

#undef __FUNC__
#define __FUNC__ "mulATB_integrate"
int32 mulATB_integrate(FMField *out, FMField *A, FMField *B,
		       VolumeGeometry *vg)
{
  int32 ii, ret = RET_OK;
  FMField *aux = 0;

  fmf_createAlloc(&aux, 1, A->nLev, A->nCol, B->nCol);

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell(out, ii);
    FMF_SetCell(vg->det, ii);
    FMF_SetCell(A, ii);
    FMF_SetCell(B, ii);
    fmf_mulATB_nn(aux, A, B);
    fmf_sumLevelsMulF(out, aux, vg->det->val);

    ERR_CheckGo(ret);
  }

 end_label:
  fmf_freeDestroy(&aux);

  return( ret );
}
