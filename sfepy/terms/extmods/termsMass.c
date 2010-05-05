#include "termsVolume.h"
#include "terms.h"
#include "geommech.h"

#undef __FUNC__
#define __FUNC__ "dw_mass"
/*!
  @par Revision history:
  - 21.11.2006, c
*/
int32 dw_mass( FMField *out, FMField *coef, FMField *state, int32 offset,
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

/*     fmf_print( bf, stdout, 0 ); */
/*     fmf_print( ftf1, stdout, 0 ); */
/*     fmf_print( ftf, stdout, 0 ); */
/*     sys_pause(); */

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];

      FMF_SetCell( out, ii );
      FMF_SetCell( coef, ii );
      FMF_SetCell( vg->det, iel );

      bf_buildFTF( ftf, ftf1 );
      fmf_mul( ftf, coef->val );
      fmf_sumLevelsMulF( out, ftf, vg->det->val );
/*       printf("%d %d\n", ii, iel); */
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
      FMF_SetCell( coef, ii );
      FMF_SetCell( vg->det, iel );

      ele_extractNodalValuesDBD( st, state, conn + nEP * iel );

      bf_act( fu, bf, st );
      bf_actt( ftfu, bf, fu );

      fmf_mul( ftfu, coef->val );
      fmf_sumLevelsMulF( out, ftfu, vg->det->val );

      ERR_CheckGo( ret );
    }
  }


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
  - 01.02.2008, c
*/
int32 dw_mass_scalar( FMField *out, FMField *coef,
		      FMField *state, FMField *bf, VolumeGeometry *vg,
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
      FMF_SetCell( vg->det, iel );
      if (coef->nCell > 1) {
	FMF_SetCell( coef, ii );
      }

      fmf_mulAF( cftf, ftf, coef->val );
      fmf_sumLevelsMulF( out, cftf, vg->det->val );

      ERR_CheckGo( ret );
    }
  } else {
    fmf_createAlloc( &st, 1, 1, 1, nEP );
    fmf_createAlloc( &fp, 1, nQP, 1, 1 );
    fmf_createAlloc( &ftfp, 1, nQP, nEP, 1 );

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];

      FMF_SetCell( out, ii );
      FMF_SetCell( vg->det, iel );
      if (coef->nCell > 1) {
	FMF_SetCell( coef, ii );
      }

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
#define __FUNC__ "d_mass_scalar"
/*!
  @par Revision history:
  - 04.09.2007, c
*/
int32 d_mass_scalar( FMField *out, FMField *coef,
		     FMField *stateP, FMField *stateQ,
		     FMField *bf, VolumeGeometry *vg,
		     int32 *conn, int32 nEl, int32 nEP,
		     int32 *elList, int32 elList_nRow )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  FMField *st = 0, *fp = 0, *ftfp = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  fmf_createAlloc( &st, 1, 1, 1, nEP );
  fmf_createAlloc( &fp, 1, nQP, 1, 1 );
  fmf_createAlloc( &ftfp, 1, nQP, nEP, 1 );

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];
    
    FMF_SetCell( out, ii );
    FMF_SetCell( vg->det, iel );
    if (coef->nCell > 1) {
      FMF_SetCell( coef, ii );
    }
    
    ele_extractNodalValuesDBD( st, stateP, conn + nEP * iel );
    
    bf_act( fp, bf, st );
    bf_actt( ftfp, bf, fp );

    ele_extractNodalValuesDBD( st, stateQ, conn + nEP * iel );
    fmf_mulATB_1n( fp, st, ftfp );
    fmf_mul( fp, coef->val );

    fmf_sumLevelsMulF( out, fp, vg->det->val );
    
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &st ); 
  fmf_freeDestroy( &fp ); 
  fmf_freeDestroy( &ftfp ); 

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_surf_mass_scalar"
/*!
  @par Revision history:
  - 09.03.2009, c
*/
int32 dw_surf_mass_scalar( FMField *out, FMField *coef,
			   FMField *state, int32 offset,
			   FMField *bf, SurfaceGeometry *sg,
			   int32 *conn, int32 nEl, int32 nEP,
			   int32 *elList, int32 elList_nRow,
			   int32 isDiff )
{
  int32 ii, iel, nFP, ret = RET_OK;
  FMField *st = 0, *fp = 0, *ftfp = 0, *ftf = 0;

  nFP = bf->nCol;

  /* output( "%d %d %d %d %d\n", offset, nEl, nEP, bf->nLev, elList_nRow ); */

  if (isDiff) {
    fmf_createAlloc( &ftf, 1, sg->nQP, nFP, nFP );

    fmf_mulATB_nn( ftf, bf, bf );

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];
      
      FMF_SetCell( out, ii );
      FMF_SetCell( sg->det, iel );
      if (coef->nCell > 1) {
	FMF_SetCell( coef, ii );
      }

      fmf_mulAF( ftf, ftf, coef->val );
      fmf_sumLevelsMulF( out, ftf, sg->det->val );

      ERR_CheckGo( ret );
    }
  } else {
    state->val = FMF_PtrFirst( state ) + offset;

    fmf_createAlloc( &st, 1, 1, 1, nFP );
    fmf_createAlloc( &fp, 1, sg->nQP, 1, 1 );
    fmf_createAlloc( &ftfp, 1, sg->nQP, nFP, 1 );

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];

      FMF_SetCell( out, ii );
      FMF_SetCell( sg->det, iel );
      if (coef->nCell > 1) {
	FMF_SetCell( coef, ii );
      }

      ele_extractNodalValuesDBD( st, state, conn + nEP * iel );
      bf_act( fp, bf, st );
      bf_actt( ftfp, bf, fp );
      fmf_mulAF( ftfp, ftfp, coef->val );
      fmf_sumLevelsMulF( out, ftfp, sg->det->val );

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
