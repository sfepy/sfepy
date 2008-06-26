#include "termsVolume.h"
#include "terms.h"
#include "geommech.h"

#undef __FUNC__
#define __FUNC__ "dw_mass"
/*!
  @par Revision history:
  - 21.11.2006, c
*/
int32 dw_mass( FMField *out, float64 coef, FMField *state, int32 offset,
	       FMField *bf, VolumeGeometry *vg,
	       int32 *conn, int32 nEl, int32 nEP,
	       int32 *elList, int32 elList_nRow,
	       int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *st = 0, *fu = 0, *ftfu = 0, *ftf1 = 0, *ftf = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

/*   output( "%d %d %d %d %d %d\n", offset, nEl, nEP, nQP, dim, elList_nRow ); */
  if (isDiff) {
    fmf_createAlloc( &ftf, 1, nQP, nEP * dim, nEP * dim );
    fmf_createAlloc( &ftf1, 1, nQP, nEP, nEP );

    fmf_mulATB_nn( ftf1, bf, bf );
    bf_buildFTF( ftf, ftf1 );

/*     fmf_print( bf, stdout, 0 ); */
/*     fmf_print( ftf1, stdout, 0 ); */
/*     fmf_print( ftf, stdout, 0 ); */
/*     sys_pause(); */

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];

      FMF_SetCell( out, ii );
      FMF_SetCell( vg->det, iel );

      fmf_sumLevelsMulF( out, ftf, vg->det->val );
/*       fmf_print( out, stdout, 0 ); */
/*       sys_pause(); */

      ERR_CheckGo( ret );
    }
  } else {
    state->val = FMF_PtrFirst( state ) + offset;

    fmf_createAlloc( &st, 1, 1, dim, nEP );
    fmf_createAlloc( &fu, 1, nQP, dim, 1 );
    fmf_createAlloc( &ftfu, 1, nQP, dim * nEP, 1 );

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];

      FMF_SetCell( out, ii );
      FMF_SetCell( vg->det, iel );

      ele_extractNodalValuesDBD( st, state, conn + nEP * iel );

      bf_act( fu, bf, st );
      bf_actt( ftfu, bf, fu );
      fmf_sumLevelsMulF( out, ftfu, vg->det->val );

      ERR_CheckGo( ret );
    }
  }

  // E.g. 1/dt.
  fmfc_mulC( out, coef );

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftf1 );
    fmf_freeDestroy( &ftf );
  } else {
    fmf_freeDestroy( &st ); 
    fmf_freeDestroy( &fu ); 
    fmf_freeDestroy( &ftfu ); 
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_mass_scalar"
/*!
  @par Revision history:
  - 04.09.2007, c
*/
int32 dw_mass_scalar( FMField *out, FMField *state, int32 offset,
		      FMField *bf, VolumeGeometry *vg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow,
		      int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *st = 0, *fp = 0, *ftfp = 0, *ftf = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

/*   output( "%d %d %d %d %d %d\n", offset, nEl, nEP, nQP, dim, elList_nRow ); */
  if (isDiff) {
    fmf_createAlloc( &ftf, 1, nQP, nEP, nEP );

    fmf_mulATB_nn( ftf, bf, bf );

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];

      FMF_SetCell( out, ii );
      FMF_SetCell( vg->det, iel );

      fmf_sumLevelsMulF( out, ftf, vg->det->val );

      ERR_CheckGo( ret );
    }
  } else {
    state->val = FMF_PtrFirst( state ) + offset;

    fmf_createAlloc( &st, 1, 1, 1, nEP );
    fmf_createAlloc( &fp, 1, nQP, 1, 1 );
    fmf_createAlloc( &ftfp, 1, nQP, nEP, 1 );

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];

      FMF_SetCell( out, ii );
      FMF_SetCell( vg->det, iel );

      ele_extractNodalValuesDBD( st, state, conn + nEP * iel );

      bf_act( fp, bf, st );
      bf_actt( ftfp, bf, fp );
      fmf_sumLevelsMulF( out, ftfp, vg->det->val );

      ERR_CheckGo( ret );
    }
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftf );
  } else {
    fmf_freeDestroy( &st ); 
    fmf_freeDestroy( &fp ); 
    fmf_freeDestroy( &ftfp ); 
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_mass_scalar_variable"
/*!
  @par Revision history:
  - 01.02.2008, c
*/
int32 dw_mass_scalar_variable( FMField *out, FMField *coef,
			       FMField *state, int32 offset,
			       FMField *bf, VolumeGeometry *vg,
			       int32 *conn, int32 nEl, int32 nEP,
			       int32 *elList, int32 elList_nRow,
			       int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *st = 0, *fp = 0, *ftfp = 0, *ftf = 0, *cftf = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

/*   output( "%d %d %d %d %d %d\n", offset, nEl, nEP, nQP, dim, elList_nRow ); */
  if (isDiff) {
    fmf_createAlloc( &ftf, 1, nQP, nEP, nEP );
    fmf_createAlloc( &cftf, 1, nQP, nEP, nEP );

    fmf_mulATB_nn( ftf, bf, bf );

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];

      FMF_SetCell( out, ii );
      FMF_SetCell( coef, iel );
      FMF_SetCell( vg->det, iel );

      fmf_mulAF( cftf, ftf, coef->val );
      fmf_sumLevelsMulF( out, cftf, vg->det->val );

      ERR_CheckGo( ret );
    }
  } else {
    state->val = FMF_PtrFirst( state ) + offset;

    fmf_createAlloc( &st, 1, 1, 1, nEP );
    fmf_createAlloc( &fp, 1, nQP, 1, 1 );
    fmf_createAlloc( &ftfp, 1, nQP, nEP, 1 );

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];

      FMF_SetCell( out, ii );
      FMF_SetCell( coef, iel );
      FMF_SetCell( vg->det, iel );

      ele_extractNodalValuesDBD( st, state, conn + nEP * iel );

      bf_act( fp, bf, st );
      bf_actt( ftfp, bf, fp );
      fmf_mul( ftfp, coef->val );
      fmf_sumLevelsMulF( out, ftfp, vg->det->val );

      ERR_CheckGo( ret );
    }
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftf );
    fmf_freeDestroy( &cftf );
  } else {
    fmf_freeDestroy( &st ); 
    fmf_freeDestroy( &fp ); 
    fmf_freeDestroy( &ftfp ); 
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_mass_scalar_fine_coarse"
/*!
  @par Revision history:
  - 05.09.2007, c
*/
int32 dw_mass_scalar_fine_coarse( FMField *out, FMField *state, int32 offset,
				  FMField *bf, FMField *cbfs,
				  VolumeGeometry *vg,
				  int32 *conn, int32 nEl, int32 nEP,
				  int32 *iemap, int32 iemap_nRow,
				  int32 *elList, int32 elList_nRow,
				  int32 isDiff )
{
  int32 ii, iel, dim, nQP, nEPR, ret = RET_OK;
  FMField *st = 0, *fp = 0, *ftfp = 0, *ftf = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;
  nEPR = vg->bfGM->nCol;

/*   output( "%d %d %d %d %d %d\n", offset, nEl, nEP, nQP, dim, elList_nRow ); */
  if (isDiff) {
    fmf_createAlloc( &ftf, 1, nQP, nEPR, nEP );

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];

      FMF_SetCell( out, ii );
      FMF_SetCell( vg->det, iel );
      FMF_SetCell( cbfs, iel );

      fmf_mulATB_nn( ftf, bf, cbfs );
      fmf_sumLevelsMulF( out, ftf, vg->det->val );

      ERR_CheckGo( ret );
    }
  } else {
    state->val = FMF_PtrFirst( state ) + offset;

    fmf_createAlloc( &st, 1, 1, 1, nEP );
    fmf_createAlloc( &fp, 1, nQP, 1, 1 );
    fmf_createAlloc( &ftfp, 1, nQP, nEPR, 1 );

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];

      FMF_SetCell( out, ii );
      FMF_SetCell( vg->det, iel );
      FMF_SetCell( cbfs, iel );

      ele_extractNodalValuesDBD( st, state, conn + nEP * iemap[iel] );

      bf_act( fp, cbfs, st );
      bf_actt( ftfp, bf, fp );
      fmf_sumLevelsMulF( out, ftfp, vg->det->val );

      ERR_CheckGo( ret );
    }
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &ftf );
  } else {
    fmf_freeDestroy( &st ); 
    fmf_freeDestroy( &fp ); 
    fmf_freeDestroy( &ftfp ); 
  }

  return( ret );
}
