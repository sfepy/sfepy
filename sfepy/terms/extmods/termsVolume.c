#include "termsVolume.h"
#include "terms.h"
#include "geommech.h"

#undef __FUNC__
#define __FUNC__ "dw_volume_lvf"
/*!
  nFEP .. field.

  @par Revision history:
  - c: 18.09.2006, r: 12.05.2008
*/
int32 dw_volume_lvf( FMField *out, FMField *bf, FMField *forceQP,
		     VolumeGeometry *vg, int32 *elList, int32 elList_nRow )
{
  int32 ii, iel, dim, nQP, nFEP, ret = RET_OK;
  FMField *outQP = 0;

  nFEP = bf->nCol;
  nQP = vg->det->nLev;
  dim = forceQP->nRow;

  fmf_createAlloc( &outQP, 1, nQP, dim * nFEP, 1 );

/*   output( "nFEP: %d, nQP: %d, dim: %d\n", nFEP, nQP, dim ); */
/*   fmf_print( bf, stdout, 1 ); */
/*   fmf_print( forceQP, stdout, 1 ); */
/*   fmf_print( outQP, stdout, 1 ); */

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( forceQP, ii );
    FMF_SetCell( vg->det, iel );

    bf_actt_c1( outQP, bf, forceQP );
    fmf_sumLevelsMulF( out, outQP, vg->det->val );
/*     fmf_print( forceQP, stdout, 0 ); */
/*     fmf_print( outQP, stdout, 0 ); */
/*     sys_pause(); */
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &outQP ); 

  return( ret );
}
