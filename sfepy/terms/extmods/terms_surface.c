#include "terms_surface.h"
#include "terms.h"
#include "geommech.h"

#undef __FUNC__
#define __FUNC__ "dw_surface_ltr"
/*!
  @par Revision history:
  - 06.09.2006, c
  - 11.10.2006
*/
int32 dw_surface_ltr( FMField *out, FMField *traction, Mapping *sg )
{
  int32 ii, dim, sym, nQP, nFP, ret = RET_OK;
  FMField *outQP = 0, *pn = 0, *stn = 0;

  nFP = sg->bf->nCol;
  nQP = sg->det->nLev;
  dim = sg->normal->nRow;
  sym = (dim + 1) * dim / 2;

  fmf_createAlloc( &outQP, 1, nQP, dim * nFP, 1 );

  if (traction->nRow == 0) {
      for (ii = 0; ii < out->nCell; ii++) {
        FMF_SetCell( out, ii );
        FMF_SetCell( sg->normal, ii );
        FMF_SetCell( sg->det, ii );
        FMF_SetCellX1( sg->bf, ii );

        bf_actt( outQP, sg->bf, sg->normal );

        fmf_sumLevelsMulF( out, outQP, sg->det->val );
        ERR_CheckGo( ret );
      }
  } else if (traction->nRow == 1) { // Pressure.
    fmf_createAlloc( &pn, 1, nQP, dim, 1 );

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCellX1( traction, ii );
      FMF_SetCell( sg->normal, ii );
      FMF_SetCell( sg->det, ii );
      FMF_SetCellX1( sg->bf, ii );

      fmf_mulAB_nn( pn, sg->normal, traction );
      bf_actt( outQP, sg->bf, pn );

      fmf_sumLevelsMulF( out, outQP, sg->det->val );
      ERR_CheckGo( ret );
    }

  } else if (traction->nRow == dim) { // Traction vector.

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCellX1( traction, ii );
      FMF_SetCell( sg->normal, ii );
      FMF_SetCell( sg->det, ii );
      FMF_SetCellX1( sg->bf, ii );

      bf_actt( outQP, sg->bf, traction );
      fmf_sumLevelsMulF( out, outQP, sg->det->val );
      ERR_CheckGo( ret );
    }

  } else if (traction->nRow == sym) { // Traction tensor.
    fmf_createAlloc( &stn, 1, nQP, dim, 1 );

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCellX1( traction, ii );
      FMF_SetCell( sg->normal, ii );
      FMF_SetCell( sg->det, ii );
      FMF_SetCellX1( sg->bf, ii );

      geme_mulAVSB3( stn, traction, sg->normal );
      bf_actt( outQP, sg->bf, stn );

      fmf_sumLevelsMulF( out, outQP, sg->det->val );
      ERR_CheckGo( ret );
    }

  } else {
    errput( ErrHead "ERR_Switch\n" );
  }

 end_label:
  fmf_freeDestroy( &outQP );
  if (traction->nCol == 1) {
    fmf_freeDestroy( &pn );
  } else if (traction->nCol == sym) {
    fmf_freeDestroy( &stn );
  }

  return( ret );
}

/*       fmf_print( trac, stdout, 0 ); */
/*       fmf_print( tracQP, stdout, 0 ); */
/*       fmf_print( pn, stdout, 0 ); */
/*       fmf_print( stn, stdout, 0 ); */
/*       fmf_print( outQP, stdout, 0 ); */
/*       fmf_print( sg->normal, stdout, 0 ); */
/*       fmf_print( sg->det, stdout, 0 ); */
/*       sys_pause(); */
