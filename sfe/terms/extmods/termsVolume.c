#include "termsVolume.h"
#include "terms.h"
#include "geommech.h"

#undef __FUNC__
#define __FUNC__ "dw_volume_lvf"
/*!
  nEP .. geometry, nFEP .. field.

  @par Revision history:
  - 18.09.2006, c
*/
int32 dw_volume_lvf( FMField *out, FMField *bf, FMField *gbf,
		     FMField *force, VolumeGeometry *vg,
		     int32 *conn, int32 nEl, int32 nEP,
		     int32 *elList, int32 elList_nRow )
{
  int32 ii, iel, dim, nQP, nFEP, ret = RET_OK;
  FMField *vf = 0, *vfQP = 0, *outQP = 0;

  nFEP = bf->nCol;
  nQP = vg->det->nLev;
  dim = vg->bfGM->nRow;

  FMF_SetFirst( force );
/*   fmf_print( force, stdout, 0 ); */

  fmf_createAlloc( &vf, 1, 1, dim, nEP );
  fmf_createAlloc( &vfQP, 1, nQP, dim, 1 );
  fmf_createAlloc( &outQP, 1, nQP, dim * nFEP, 1 );

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg->det, iel );

    ele_extractNodalValuesDBD( vf, force, conn + nEP * iel );
    bf_act( vfQP, gbf, vf );
    bf_actt_c1( outQP, bf, vfQP );
    fmf_sumLevelsMulF( out, outQP, vg->det->val );
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &vf ); 
  fmf_freeDestroy( &vfQP ); 
  fmf_freeDestroy( &outQP ); 

  return( ret );
}
/*     fmf_print( vf, stdout, 0 ); */
/*     fmf_print( gbf, stdout, 0 ); */
/*     fmf_print( bf, stdout, 0 ); */
/*     fmf_print( vfQP, stdout, 0 ); */
/*     sys_pause(); */
