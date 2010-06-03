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
#define __FUNC__ "dw_llaplace"
int32 dw_llaplace( FMField *out, FMField *state,
		   FMField *coef, FMField *coef2, VolumeGeometry *vg,
		   int32 *conn, int32 nEl, int32 nEP,
		   int32 *elList, int32 elList_nRow,
		   int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK, iqp, iep, jj;
  FMField *gtd11 = 0, *gtd11a = 0;
  FMField *st = 0;
  float64 *pg1, *pg2, *pg3, *pdiff;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  // output( "%d %d %d %d %d %d %d\n", offset, nEl, nEP, nQP, dim, elList_nRow, isDiff);
  if (isDiff) {
    fmf_createAlloc( &gtd11, 1, nQP, nEP, nEP );
  } else {
    state->val = FMF_PtrFirst( state );
    
    fmf_createAlloc( &gtd11, 1, nQP, nEP, nEP );
    fmf_createAlloc( &gtd11a, 1, 1, nEP, nEP );
    fmf_createAlloc( &st, 1, 1, nEP, 1 );
  }
  
  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];
    
    FMF_SetCell( out, ii );
    FMF_SetCell( coef, ii );
    FMF_SetCell( coef2, ii );
    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( vg->det, iel );
    
    switch ( dim ) {
    case 3:
      for (iqp = 0; iqp < nQP; iqp++) {
	pg1 = FMF_PtrLevel( vg->bfGM, iqp );
	pg2 = pg1 + nEP;
	pg3 = pg2 + nEP;

	for (iep = 0; iep < nEP; iep++) {
	  pdiff = FMF_PtrLevel( gtd11, iqp ) + nEP * iep;
	  for (jj = 0; jj < nEP; jj++) {
	    pdiff[jj] = coef->val[iqp] * (pg1[iep] * pg1[jj] + pg2[iep] * pg2[jj]);
	    pdiff[jj] += coef2->val[iqp] * pg3[iep] * pg3[jj];
	  }
	} /* for (iep) */
      } /* for (iqp) */
      break;
    case 2:
      for (iqp = 0; iqp < nQP; iqp++) {
	pg1 = FMF_PtrLevel( vg->bfGM, iqp );
	pg2 = pg1 + nEP;

	for (iep = 0; iep < nEP; iep++) {
	  pdiff = FMF_PtrLevel( gtd11, iqp ) + nEP * iep;
	  for (jj = 0; jj < nEP; jj++) {
	    pdiff[jj] = coef->val[iqp] * pg1[iep] * pg1[jj];
	    pdiff[jj] += coef2->val[iqp] * pg2[iep] * pg2[jj];
	  }
	} /* for (iep) */
      } /* for (iqp) */
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }

    if (isDiff) {
      fmf_sumLevelsMulF( out, gtd11, vg->det->val );
    }
    else {
      ele_extractNodalValuesNBN( st, state, conn + nEP * iel );
      fmf_sumLevelsMulF( gtd11a, gtd11, vg->det->val );
      fmf_mulAB_nn( out, gtd11a, st );
    }

    ERR_CheckGo( ret );
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &gtd11 );
  } else {
    fmf_freeDestroy( &gtd11 );
    fmf_freeDestroy( &gtd11a );
    fmf_freeDestroy( &st );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_llaplace"
int32 d_llaplace( FMField *out, FMField *stateU, FMField *stateV,
		  FMField *coef, FMField *coef2, VolumeGeometry *vg,
		  int32 *conn, int32 nEl, int32 nEP,
		  int32 *elList, int32 elList_nRow )
{
  int32 ii, iel, dim, nQP, ret = RET_OK, iqp, iep, jj;
  FMField *nodval = 0;
  FMField *aux = 0, *aux2 = 0, *aux3 = 0;
  float64 *pg1, *pg2, *pg3, *pdiff;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  FMF_SetFirst( out );

  fmf_createAlloc( &aux, 1, nQP, nEP, nEP );
  fmf_createAlloc( &aux2, 1, 1, nEP, nEP );
  fmf_createAlloc( &aux3, 1, 1, nEP, 1 );
  fmf_createAlloc( &nodval, 1, 1, nEP, 1 );

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];
    
    FMF_SetCell( coef, ii );
    FMF_SetCell( coef2, ii );
    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( vg->det, iel );

    switch ( dim ) {
    case 3:
      for (iqp = 0; iqp < nQP; iqp++) {
	pg1 = FMF_PtrLevel( vg->bfGM, iqp );
	pg2 = pg1 + nEP;
	pg3 = pg2 + nEP;
	
	for (iep = 0; iep < nEP; iep++) {
	  pdiff = FMF_PtrLevel( aux, iqp ) + nEP * iep;
	  for (jj = 0; jj < nEP; jj++) {
	    pdiff[jj] = coef->val[iqp] * (pg1[iep] * pg1[jj] + pg2[iep] * pg2[jj]);
	    pdiff[jj] += coef2->val[iqp] * pg3[iep] * pg3[jj];
	  }
	} /* for (iep)*/
      } /* for (iqp)*/
      break;
    case 2:
      for (iqp = 0; iqp < nQP; iqp++) {
	pg1 = FMF_PtrLevel( vg->bfGM, iqp );
	pg2 = pg1 + nEP;

	for (iep = 0; iep < nEP; iep++) {
	  pdiff = FMF_PtrLevel( aux, iqp ) + nEP * iep;
	  for (jj = 0; jj < nEP; jj++) {
	    pdiff[jj] = coef->val[iqp] * pg1[iep] * pg1[jj];
	    pdiff[jj] += coef2->val[iqp] * pg2[iep] * pg2[jj];
	  }
	} /* for (iep)*/
      } /* for (iqp)*/
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }
    
    fmf_sumLevelsMulF( aux2, aux, vg->det->val );

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
  fmf_freeDestroy( &nodval );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_acoustic_surface"
int32 d_acoustic_surface( FMField *out, FMField *in,
			  FMField *coef, FMField *coef2, SurfaceGeometry *sg,
			  int32 *conn, int32 nEl, int32 nEP,
			  int32 *elList, int32 elList_nRow )
{
  int32 ii, iqp, iep, iel, dim, nQP, ret = RET_OK;
  FMField *nodval = 0, *aux = 0, *aux2 = 0, *aux3 = 0;
  float64 *pdiff;
  
  nQP = sg->det->nLev;
  dim = sg->normal->nRow;

  fmf_createAlloc( &aux, 1, nQP, dim, nEP );
  fmf_createAlloc( &aux2, 1, nQP, dim, 1 );
  fmf_createAlloc( &aux3, 1, nQP, 1, 1 );
  fmf_createAlloc( &nodval, 1, 1, nEP, 1 );

  //  output( "(%d %d %d) (%d %d %d)\n", sg->normal->nLev, sg->normal->nRow, sg->normal->nCol, sg->bfBGM->nLev, sg->bfBGM->nRow, sg->bfBGM->nCol );

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( coef, ii );
    FMF_SetCell( coef2, ii );
    FMF_SetCell( sg->normal, ii );
    FMF_SetCell( sg->det, ii );
    FMF_SetCell( sg->bfBGM, ii );

    fmf_copy( aux, sg->bfBGM );

    for (iqp = 0; iqp < nQP; iqp++) {
      /* x */
      pdiff = FMF_PtrLevel( aux, iqp );
      for (iep = 0; iep < nEP; iep++)
	pdiff[iep] *= coef->val[iqp];
      /* y */
      if ( dim > 2 ) {
	pdiff += nEP;
	for (iep = 0; iep < nEP; iep++)
	  pdiff[iep] *= coef->val[iqp];
      }
      /* z */
      pdiff += nEP;
      for (iep = 0; iep < nEP; iep++)
	pdiff[iep] *= coef2->val[iqp];
    } /* for (iqp)*/

    ele_extractNodalValuesNBN( nodval, in, conn + nEP * iel );
    fmf_mulAB_n1( aux2, aux, nodval );
    fmf_mulATB_nn( aux3, sg->normal, aux2 );
    fmf_sumLevelsMulF( out, aux3, sg->det->val );

    ERR_CheckGo( ret );
  } /* for (ii) */

 end_label:
  fmf_freeDestroy( &aux );
  fmf_freeDestroy( &aux2 );
  fmf_freeDestroy( &aux3 );
  fmf_freeDestroy( &nodval );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_acoustic_integrate"
int32 dw_acoustic_integrate( FMField *out, FMField *coef, VolumeGeometry *vg,
			     int32 *conn, int32 nEl, int32 nEP,
			     int32 *elList, int32 elList_nRow,
			     int32 isDiff )
{
  int32 ii, iel, dim, nQP, ret = RET_OK, iqp, iep;
  FMField *aux = 0;
  float64 *pg1, *pg2, *pdiff, *pcoef;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  fmf_createAlloc( &aux, 1, nQP, nEP, 1 );
  
  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];
    
    FMF_SetCell( out, ii );
    FMF_SetCell( coef, ii );
    FMF_SetCell( vg->bfGM, iel );
    FMF_SetCell( vg->det, iel );
    
    switch ( dim ) {
    case 3:
      for (iqp = 0; iqp < nQP; iqp++) {
	pg1 = FMF_PtrLevel( vg->bfGM, iqp );
	pg2 = pg1 + nEP;
	pdiff = FMF_PtrLevel( aux, iqp );
	pcoef = FMF_PtrLevel( coef, iqp );
	for (iep = 0; iep < nEP; iep++) {
	  pdiff[iep] = pcoef[0] * pg1[iep] +  pcoef[1] * pg2[iep];
	} /* for (iep) */
      } /* for (iqp) */
      break;
    case 2:
      for (iqp = 0; iqp < nQP; iqp++) {
	pg1 = FMF_PtrLevel( vg->bfGM, iqp );
	pdiff = FMF_PtrLevel( aux, iqp );
	for (iep = 0; iep < nEP; iep++) {
	  pdiff[iep] = coef->val[iqp] * pg1[iep];
	} /* for (iep) */
      } /* for (iqp) */
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }
    
    fmf_sumLevelsMulF( out, aux, vg->det->val );

    ERR_CheckGo( ret ); 
  }

 end_label:
  fmf_freeDestroy( &aux );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_acoustic_integrate"
int32 d_acoustic_alpha( FMField *out, FMField *in,
			VolumeGeometry *vg,
			int32 *conn, int32 nEl, int32 nEP,
			int32 *elList, int32 elList_nRow )
{
  int32 ii, j, iel, dim, nQP, ret = RET_OK;
  FMField *nodval = 0, *aux = 0, *aux2 = 0;
  
  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

  FMF_SetFirst( out );

  fmf_createAlloc( &aux, 1, nQP, dim, 1 );
  fmf_createAlloc( &aux2, 1, 1, dim, 1 );
  fmf_createAlloc( &nodval, 1, 1, nEP, 1 );

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];

    FMF_SetCell( out, ii );
    FMF_SetCell( vg->det, iel );
    FMF_SetCell( vg->bfGM, iel );

    ele_extractNodalValuesNBN( nodval, in, conn + nEP * iel );
    fmf_mulAB_n1( aux, vg->bfGM, nodval );
    fmf_sumLevelsMulF( aux2, aux, vg->det->val );

    for (j = 0; j < (dim-1); j++)
      out->val[j] = aux2->val[j];

    FMF_SetCellNext( out );

    ERR_CheckGo( ret );
  } /* for (ii) */

 end_label:
  fmf_freeDestroy( &aux );
  fmf_freeDestroy( &aux2 );
  fmf_freeDestroy( &nodval );

  return( ret );
}
