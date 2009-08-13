#include "termsHDPM.h"
#include "formSDCC.h"
#include "geommech.h"
#include "terms.h"

#undef __FUNC__
#define __FUNC__ "dw_hdpm_cache"
/*!
  @par Revision history:
  - 28.11.2006, c
  - 29.11.2006
  - 12.03.2007
  - 13.03.2007
  - 24.04.2007
*/
int32 dw_hdpm_cache( FMField *out, FMField *state, int32 offset,
		     FMField *bf, VolumeGeometry *vg,
		     int32 *conn, int32 nEl, int32 nEP,
		     int32 mode )
{
  int32 ii, dim, sym, nQP, ret = RET_OK;
  FMField *st = 0, *disG = 0;

  state->val = FMF_PtrFirst( state ) + offset;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  if (mode == 0) {
    sym = (dim + 1) * dim / 2;
    fmf_createAlloc( &st, 1, 1, nEP, dim );
    fmf_createAlloc( &disG, 1, nQP, dim, dim );

    for (ii = 0; ii < nEl; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCell( vg->bfGM, ii );

      ele_extractNodalValuesNBN( st, state, conn + nEP * ii );
      fmf_mulAB_n1( disG, vg->bfGM, st );
      form_sdcc_strainCauchy_VS( out, disG );
/*       fmf_print( out, stdout, 0 ); */
/*       sys_pause(); */
      ERR_CheckGo( ret );
    }

  } else if (mode == 2) {
    fmf_createAlloc( &st, 1, 1, nEP, 1 );

    for (ii = 0; ii < nEl; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCell( vg->bfGM, ii );

      ele_extractNodalValuesNBN( st, state, conn + nEP * ii );
      fmf_mulAB_n1( out, vg->bfGM, st );

      ERR_CheckGo( ret );
    }
  } else if (mode == 3) {
    FMField gcl[1], stv[1];

    fmf_createAlloc( &st, 1, 1, dim, nEP );
    stv->nAlloc = -1;
    fmf_pretend( stv, 1, 1, nEP * dim, 1, st->val );

    gcl->nAlloc = -1;
    fmf_pretend( gcl, 1, nQP, 1, nEP * dim, vg->bfGM->val0 );

    for (ii = 0; ii < nEl; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCell( gcl, ii );

      ele_extractNodalValuesDBD( st, state, conn + nEP * ii );
      fmf_mulAB_n1( out, gcl, stv );

      ERR_CheckGo( ret );
    }
  }

 end_label:
  fmf_freeDestroy( &st ); 
  fmf_freeDestroy( &disG ); 

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_hdpm_surfdvel"
/*!
  @par Revision history:
  - 04.05.2007, c
  - 09.05.2007
*/
int32 d_hdpm_surfdvel( FMField *out, FMField *state, int32 offset,
		       FMField *mtxD, SurfaceGeometry *sg,
		       int32 *fis, int32 nFa, int32 nFP,
		       int32 *faceList, int32 faceList_nRow,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow )
{
  // remove fis!
  int32 ii, iel, ifa, dim, nQP, ret = RET_OK;
  FMField *st = 0, *gp = 0, *dgp = 0, *ntdgp = 0;

  nQP = sg->normal->nLev;
  dim = sg->normal->nRow;

  if (elList_nRow != faceList_nRow) {
    errput( "dimensions mismatch! (%d == %d)\n", elList_nRow, faceList_nRow );
  }

/*   output( "%d %d %d %d %d %d\n", offset, nEl, nEP, nQP, dim, elList_nRow ); */
  state->val = FMF_PtrFirst( state ) + offset;

  fmf_createAlloc( &st, 1, 1, nEP, 1 );
  fmf_createAlloc( &gp, 1, nQP, dim, 1 );
  fmf_createAlloc( &dgp, 1, nQP, dim, 1 );
  fmf_createAlloc( &ntdgp, 1, nQP, 1, 1 );

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];
    ifa = faceList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( sg->normal, ifa );
    FMF_SetCell( sg->bfBGM, ifa );
    FMF_SetCell( sg->det, ifa );

#ifdef DEBUG_FMF
    if (fis[ifa*nFP+0] != iel) {
      errput( "%d == %d\n", fis[ifa*nFP+0], iel );
    }
#endif

    ele_extractNodalValuesNBN( st, state, conn + nEP * iel );
    fmf_mulAB_n1( gp, sg->bfBGM, st );
    fmf_mulAB_nn( dgp, mtxD, gp );
    fmf_mulATB_nn( ntdgp, sg->normal, dgp );

/*     fmf_print( mtxD, stdout, 0 ); */
/*     fmf_print( sg->normal, stdout, 0 ); */
/*     fmf_print( sg->bfBGM, stdout, 0 ); */
/*     fmf_print( st, stdout, 0 ); */
/*     fmf_print( gp, stdout, 0 ); */
/*     fmf_print( dgp, stdout, 0 ); */
/*     output( "%d %d %d\n", ii, iel, ifa ); */
/*     sys_pause(); */
    
    fmf_sumLevelsMulF( out, ntdgp, sg->det->val );
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &st ); 
  fmf_freeDestroy( &gp ); 
  fmf_freeDestroy( &dgp ); 
  fmf_freeDestroy( &ntdgp ); 

  return( ret );
}
