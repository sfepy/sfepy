#include "terms_volume.h"
#include "terms.h"
#include "geommech.h"

#undef __FUNC__
#define __FUNC__ "dw_volume_lvf"
/*!
  nFEP .. field.

  @par Revision history:
  - c: 18.09.2006, r: 12.05.2008
*/
int32 dw_volume_lvf( FMField *out, FMField *forceQP, Mapping *vg )
{
  int32 ii, dim, nQP, nFEP, ret = RET_OK;
  FMField *outQP = 0;

  nFEP = vg->bf->nCol;
  nQP = vg->det->nLev;
  dim = forceQP->nRow;

  fmf_createAlloc( &outQP, 1, nQP, dim * nFEP, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCellX1( forceQP, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCellX1( vg->bf, ii );

    bf_actt_c1( outQP, vg->bf, forceQP );
    fmf_sumLevelsMulF( out, outQP, vg->det->val );
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &outQP );

  return( ret );
}
