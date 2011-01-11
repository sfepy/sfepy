#include "termsAcoustic.h"
#include "terms.h"

#undef __FUNC__
#define __FUNC__ "d_llaplace_p_sa"
/* 09.04.2009, c */
int32 d_llaplace_p_sa( FMField *out,
		       FMField *stateU, FMField *stateV, FMField *stateW,
		       VolumeGeometry *vg,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 mode,
		       int32 *elList, int32 elList_nRow )
{
  int32 ii, iel, dim, nQP, ret = RET_OK, iqp, iep, jj;
  FMField *nodval = 0, *nodval2 = 0;
  FMField *aux = 0, *aux2 = 0, *aux3 = 0, *aux4 = 0, *aux5 = 0;
  float64 *pg1, *pg2, *pg3, *pdiff, *pnodval;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  FMF_SetFirst( out );

  fmf_createAlloc( &aux, 1, nQP, nEP, nEP );
  fmf_createAlloc( &aux2, 1, 1, nEP, nEP );
  fmf_createAlloc( &aux3, 1, 1, nEP, 1 );
  fmf_createAlloc( &nodval, 1, 1, nEP, 1 );
  fmf_createAlloc( &nodval2, 1, 1, nEP, dim );
  if ( mode == 1 ) {
      fmf_createAlloc( &aux4, 1, 1, 1, nQP );
  }
  else {
    fmf_createAlloc( &aux4, 1, nQP, nEP, nEP );
    fmf_createAlloc( &aux5, 1, nQP, nEP, nEP );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( vg->det, iel );
    ele_extractNodalValuesNBN( nodval2, stateW, conn + nEP * iel );

    if ( mode == 1 )
      fmf_fillC( aux4, 0.0 );

    switch ( dim ) {
    case 3:
      for (iqp = 0; iqp < nQP; iqp++) {
	pg1 = FMF_PtrLevel( vg->bfGM, iqp );
	pg2 = pg1 + nEP;
	pg3 = pg2 + nEP;

	for (iep = 0; iep < nEP; iep++) {
	  /* Ba^T * Ba */
	  pdiff = FMF_PtrLevel( aux, iqp ) + nEP * iep;
	  for (jj = 0; jj < nEP; jj++)
	    pdiff[jj] = pg1[iep] * pg1[jj] + pg2[iep] * pg2[jj];
	  if ( mode == 1 ) {
	    /* div(V) */
	    pnodval = nodval2->val + iep*dim;
	    aux4->val[iqp] += vg->det->val[iqp] *
	      ( pnodval[0]*pg1[iep] +
		pnodval[1]*pg2[iep] +
		pnodval[2]*pg3[iep] );
	  }
	} /* for (iep) */
      } /* for (iqp) */
      break;
    case 2:
      for (iqp = 0; iqp < nQP; iqp++) {
	pg1 = FMF_PtrLevel( vg->bfGM, iqp );
	pg2 = pg1 + nEP;

	for (iep = 0; iep < nEP; iep++) {
	  /* Ba^T * Ba */
	  pdiff = FMF_PtrLevel( aux, iqp ) + nEP * iep;
	  for (jj = 0; jj < nEP; jj++)
	    pdiff[jj] = pg1[iep] * pg1[jj];
	  if ( mode == 1 ) {
	    /* div(V) */
	    pnodval = nodval2->val + iep*dim;
	    aux4->val[iqp] += vg->det->val[iqp] *
	      ( pnodval[0]*pg1[iep] + pnodval[1]*pg2[iep] );
	  }
	} /* for (iel) */
      } /* for (iqp) */
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }

    if ( mode == 1 ) {
      fmf_sumLevelsMulF( aux2, aux, aux4->val );
    }
    else {
      fmf_mulAB_1n( aux4, nodval2, vg->bfGM );
      fmf_mulAB_nn( aux5, aux, aux4);
      fmf_sumLevelsMulF( aux2, aux5, vg->det->val );
    }

    ele_extractNodalValuesNBN( nodval, stateU, conn + nEP * iel );
    fmf_mulAB_nn( aux3, aux2, nodval );
    ele_extractNodalValuesNBN( nodval, stateV, conn + nEP * iel );
    fmf_mulATB_nn( out, nodval, aux3 );

    FMF_SetCellNext( out );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &aux );
  fmf_freeDestroy( &aux2 );
  fmf_freeDestroy( &aux3 );
  fmf_freeDestroy( &aux4 );
  fmf_freeDestroy( &aux5 );
  fmf_freeDestroy( &nodval );
  fmf_freeDestroy( &nodval2 );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_llaplace_t_sa"
/* 09.04.2009, c */
int32 d_llaplace_t_sa( FMField *out,
		       FMField *stateU, FMField *stateV, FMField *stateW,
		       VolumeGeometry *vg,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 mode,
		       int32 *elList, int32 elList_nRow )
{
  int32 ii, iel, dim, nQP, ret = RET_OK, iqp, iep, jj;
  FMField *nodval = 0, *nodval2 = 0;
  FMField *aux = 0, *aux2 = 0, *aux3 = 0, *aux4 = 0, *aux5 = 0;
  float64 *pg1, *pg2, *pg3, *pdiff, *pnodval;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  FMF_SetFirst( out );

  fmf_createAlloc( &aux, 1, nQP, nEP, nEP );
  fmf_createAlloc( &aux2, 1, 1, nEP, nEP );
  fmf_createAlloc( &aux3, 1, 1, nEP, 1 );
  fmf_createAlloc( &nodval, 1, 1, nEP, 1 );
  fmf_createAlloc( &nodval2, 1, 1, nEP, dim );
  if ( mode == 1 ) {
    fmf_createAlloc( &aux4, 1, 1, 1, nQP );
  }
  else {
    fmf_createAlloc( &aux4, 1, nQP, nEP, nEP );
    fmf_createAlloc( &aux5, 1, nQP, nEP, nEP );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( vg->det, iel );
    ele_extractNodalValuesNBN( nodval2, stateW, conn + nEP * iel );

    if ( mode == 1 )
      fmf_fillC( aux4, 0.0 );

    switch ( dim ) {
    case 3:
      for (iqp = 0; iqp < nQP; iqp++) {
	pg1 = FMF_PtrLevel( vg->bfGM, iqp );
	pg2 = pg1 + nEP;
	pg3 = pg2 + nEP;

	for (iep = 0; iep < nEP; iep++) {
	  /* Bz^T * Bz */
	  pdiff = FMF_PtrLevel( aux, iqp ) + nEP * iep;
	  for (jj = 0; jj < nEP; jj++)
	    pdiff[jj] = pg3[iep] * pg3[jj];
	  if ( mode == 1 ) {
	    /* div(V) */
	    pnodval = nodval2->val + iep*dim;
	    aux4->val[iqp] += vg->det->val[iqp] *
	      ( pnodval[0]*pg1[iep] +
		pnodval[1]*pg2[iep] +
		pnodval[2]*pg3[iep] );
	  }
	} /* for (iep) */
      } /* for (iqp) */
      break;
    case 2:
      for (iqp = 0; iqp < nQP; iqp++) {
	pg1 = FMF_PtrLevel( vg->bfGM, iqp );
	pg2 = pg1 + nEP;

	for (iep = 0; iep < nEP; iep++) {
	  /* Bz^T * Bz */
	  pdiff = FMF_PtrLevel( aux, iqp ) + nEP * iep;
	  for (jj = 0; jj < nEP; jj++)
	    pdiff[jj] = pg2[iep] * pg2[jj];
	  if ( mode == 1 ) {
	    /* div(V) */
	    pnodval = nodval2->val + iep*dim;
	    aux4->val[iqp] += vg->det->val[iqp] *
	      ( pnodval[0]*pg1[iep] + pnodval[1]*pg2[iep] );
	  }
	} /* for (iep) */
      } /* for (iqp) */
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }

    if ( mode == 1 ) {
      fmf_sumLevelsMulF( aux2, aux, aux4->val );
    }
    else {
      fmf_mulAB_1n( aux4, nodval2, vg->bfGM );
      fmf_mulAB_nn( aux5, aux, aux4);
      fmf_sumLevelsMulF( aux2, aux5, vg->det->val );
    }

    ele_extractNodalValuesNBN( nodval, stateU, conn + nEP * iel );
    fmf_mulAB_nn( aux3, aux2, nodval );
    ele_extractNodalValuesNBN( nodval, stateV, conn + nEP * iel );
    fmf_mulATB_nn( out, nodval, aux3 );
    FMF_SetCellNext( out );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &aux );
  fmf_freeDestroy( &aux2 );
  fmf_freeDestroy( &aux3 );
  fmf_freeDestroy( &aux4 );
  fmf_freeDestroy( &aux5 );
  fmf_freeDestroy( &nodval );
  fmf_freeDestroy( &nodval2 );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_surf_llaplace"
int32 dw_surf_llaplace( FMField *out, FMField *state, FMField *coef,
			FMField *gbf, SurfaceGeometry *sg,
		        int32 *conn, int32 nEl, int32 nEP,
		        int32 *elList, int32 elList_nRow,
		        int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK;
  //int32 iqp, iep, jj;
  FMField *bb = 0, *bb0 = 0, *aux = 0;
  FMField *st = 0;
  //float64 *pg1, *pg2, *pdiff;

  nQP = gbf->nLev;
  //  dim = gbf->nRow + 1;
  dim = gbf->nRow;

  /*  output( "%d %d %d %d %d %d\n", nEl, nEP, nQP, dim, elList_nRow, isDiff); */

  fmf_createAlloc( &bb, 1, nQP, nEP, nEP );
  fmf_createAlloc( &aux, 1, nQP, dim, nEP );

  if (!isDiff) {
    state->val = FMF_PtrFirst( state );
    fmf_createAlloc( &bb0, 1, 1, nEP, nEP );
    fmf_createAlloc( &st, 1, 1, nEP, 1 );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    if (coef->nCell > 1) {
      FMF_SetCell( coef, ii );
    }

    FMF_SetCell( out, ii );
    FMF_SetCell( sg->det, iel );

    fmf_mulAB_1n( aux, coef, gbf );
    fmf_mulATB_nn( bb, gbf, aux );

    if (isDiff) {
      fmf_sumLevelsMulF( out, bb, sg->det->val );
    }
    else {
      ele_extractNodalValuesNBN( st, state, conn + nEP * iel );
      fmf_sumLevelsMulF( bb0, bb, sg->det->val );
      fmf_mulAB_nn( out, bb0, st );
    }

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &bb );
  fmf_freeDestroy( &aux );
  if (!isDiff) {
    fmf_freeDestroy( &bb0 );
    fmf_freeDestroy( &st );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_surf_lcouple"
int32 dw_surf_lcouple( FMField *out, FMField *state, FMField *coef,
		       FMField *bf, FMField *gbf, SurfaceGeometry *sg,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow,
		       int32 isDiff )
{
  int32 ii, iel, nn, nQP, ret = RET_OK;
  FMField *bb = 0, *bb0 = 0, *aux = 0;
  FMField *st = 0;

  nQP = gbf->nLev;
  nn = gbf->nRow;

  /* output( "%d %d %d %d %d %d\n", nEl, nEP, nQP, dim, elList_nRow, isDiff); */

  fmf_createAlloc( &bb, 1, nQP, nEP, nEP );
  fmf_createAlloc( &aux, 1, nQP, nn, nEP );

  if (!isDiff) {
    state->val = FMF_PtrFirst( state );
    fmf_createAlloc( &bb0, 1, 1, nEP, nEP );
    fmf_createAlloc( &st, 1, 1, nEP, 1 );
  }

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    if (coef->nCell > 1) {
      FMF_SetCell( coef, ii );
    }

    FMF_SetCell( out, ii );
    FMF_SetCell( sg->det, iel );

    fmf_mulAB_1n( aux, coef, bf );
    fmf_mulATB_nn( bb, gbf, aux );

    if (isDiff) {
      fmf_sumLevelsMulF( out, bb, sg->det->val );
    }
    else {
      ele_extractNodalValuesNBN( st, state, conn + nEP * iel );
      fmf_sumLevelsMulF( bb0, bb, sg->det->val );
      fmf_mulAB_nn( out, bb0, st );
    }

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &bb );
  fmf_freeDestroy( &aux );
  if (!isDiff) {
    fmf_freeDestroy( &bb0 );
    fmf_freeDestroy( &st );
  }

  return( ret );
}
