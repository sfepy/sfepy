#include "formSDCC.h"
#include "terms_elastic.h"
#include "terms.h"

#undef __FUNC__
#define __FUNC__ "mat_le_tanModuli11"
/*!
  @par Revision history:
  - 17.03.2003, c
  - 31.01.2006
  - 06.02.2006
  - 07.03.2006, adopted from mafest1
*/
int32 mat_le_tanModuli11( FMField *mtx, FMField *lam, FMField *mu, int32 mode  )
#define MAT_LE_AuxMacro1_3D \
    do { for (iqp = 0; iqp < nQP; iqp++) { \
      _lam = lam->val[iqp]; \
      _mu = mu->val[iqp]; \
      pd[0] = _lam + 2.0 * _mu; \
      pd[1] = _lam; \
      pd[2] = _lam; \
      pd[6] = _lam; \
      pd[7] = _lam + 2.0 * _mu; \
      pd[8] = _lam; \
      pd[12] = _lam; \
      pd[13] = _lam; \
      pd[14] = _lam + 2.0 * _mu; \
      pd[21] = _mu; \
      pd[28] = _mu; \
      pd[35] = _mu; \
      pd += 36; \
    } } while (0)
#define MAT_LE_AuxMacro2_3D \
    do { for (iqp = 0; iqp < nQP; iqp++) { \
      _lam = lam->val[iqp]; \
      _mu = mu->val[iqp]; \
      mu23 = _mu * (2.0/3.0); \
      mu43 = 2.0 * mu23; \
      pd[0] = mu43; \
      pd[1] = -mu23; \
      pd[2] = -mu23; \
      pd[6] = -mu23; \
      pd[7] = mu43; \
      pd[8] = -mu23; \
      pd[12] = -mu23; \
      pd[13] = -mu23; \
      pd[14] = mu43; \
      pd[21] = _mu; \
      pd[28] = _mu; \
      pd[35] = _mu; \
      pd += 36; \
    } } while (0)
#define MAT_LE_AuxMacro1_2D \
    do { for (iqp = 0; iqp < nQP; iqp++) { \
      _lam = lam->val[iqp]; \
      _mu = mu->val[iqp]; \
      pd[0] = _lam + 2.0 * _mu; \
      pd[1] = _lam; \
      pd[3] = _lam; \
      pd[4] = _lam + 2.0 * _mu; \
      pd[8] = _mu; \
      pd += 9; \
    } } while (0)
#define MAT_LE_AuxMacro2_2D \
    do { for (iqp = 0; iqp < nQP; iqp++) { \
      _lam = lam->val[iqp]; \
      _mu = mu->val[iqp]; \
      mu23 = _mu * (2.0/3.0); \
      mu43 = 2.0 * mu23; \
      pd[0] = mu43; \
      pd[1] = -mu23; \
      pd[3] = -mu23; \
      pd[4] = mu43; \
      pd[8] = _mu; \
      pd += 9; \
    } } while (0)
{
  float64 *pd;
  float64 _lam, _mu;
  int32 nQP, iqp, sym;

  nQP = mtx->nLev;
  sym = mtx->nRow;

  pd = FMF_PtrCurrent( mtx );

  if (sym == 6) {
    if (1) {
      MAT_LE_AuxMacro1_3D;
    } else {
      float64 mu23, mu43;
      MAT_LE_AuxMacro2_3D;
    }
  } else if (sym == 3) {
    if (1) {
      MAT_LE_AuxMacro1_2D;
    } else {
      float64 mu23, mu43;
      MAT_LE_AuxMacro2_2D;
    }
  }

  return( RET_OK );
}
#undef MAT_LE_AuxMacro1_3D
#undef MAT_LE_AuxMacro2_3D
#undef MAT_LE_AuxMacro1_2D
#undef MAT_LE_AuxMacro2_2D

#undef __FUNC__
#define __FUNC__ "mat_le_stress"
/*!
  In mixed case, builds only the deviatoric part of stress.

  @par Revision history:
  - 17.03.2003, c
  - 31.01.2006
  - 06.02.2006
  - 07.03.2006, adopted from mafest1
*/
int32 mat_le_stress( FMField *stress, FMField *strain,
		     FMField *lam, FMField *mu )
{
  int32 iell, iqp;
  int32 sym, nQP;
  float64 *pstress, *pstrain;
  float64 _lam, _mu, mu23, mu43, l2m;

  nQP = stress->nLev;
  sym = stress->nRow;

  if (sym == 6) {
    for (iell = 0; iell < stress->nCell; iell++) {
      FMF_SetCell( lam, iell );
      FMF_SetCell( mu, iell );
      pstress = FMF_PtrCell( stress, iell );
      pstrain = FMF_PtrCell( strain, iell );
      if (1) {
	for (iqp = 0; iqp < nQP; iqp++) {
	  _lam = lam->val[iqp];
	  _mu = mu->val[iqp];
	  l2m = 2.0 * _mu + _lam;
	  pstress[0] = l2m * pstrain[0] + _lam * (pstrain[1] + pstrain[2]);
	  pstress[1] = l2m * pstrain[1] + _lam * (pstrain[0] + pstrain[2]);
	  pstress[2] = l2m * pstrain[2] + _lam * (pstrain[0] + pstrain[1]);
	  pstress[3] = _mu * pstrain[3];
	  pstress[4] = _mu * pstrain[4];
	  pstress[5] = _mu * pstrain[5];
	  pstress += sym;
	  pstrain += sym;
	}
      } else {
	for (iqp = 0; iqp < nQP; iqp++) {
	  _lam = lam->val[iqp];
	  _mu = mu->val[iqp];
	  mu23 = _mu * (2.0/3.0);
	  mu43 = 2.0 * mu23;
	  pstress[0] = mu43 * pstrain[0] - mu23 * (pstrain[1] + pstrain[2]);
	  pstress[1] = mu43 * pstrain[1] - mu23 * (pstrain[0] + pstrain[2]);
	  pstress[2] = mu43 * pstrain[2] - mu23 * (pstrain[0] + pstrain[1]);
	  pstress[3] = _mu * pstrain[3];
	  pstress[4] = _mu * pstrain[4];
	  pstress[5] = _mu * pstrain[5];
	  pstress += sym;
	  pstrain += sym;
	}
      }
    }
  } else if (sym == 3) {
    for (iell = 0; iell < stress->nCell; iell++) {
      pstress = FMF_PtrCell( stress, iell );
      pstrain = FMF_PtrCell( strain, iell );
      if (1) {
	for (iqp = 0; iqp < nQP; iqp++) {
	  _lam = lam->val[iqp];
	  _mu = mu->val[iqp];
	  l2m = 2.0 * _mu + _lam;
	  pstress[0] = l2m * pstrain[0] + _lam * (pstrain[1]);
	  pstress[1] = l2m * pstrain[1] + _lam * (pstrain[0]);
	  pstress[2] = _mu * pstrain[2];
	  pstress += sym;
	  pstrain += sym;
	}
      } else {
	for (iqp = 0; iqp < nQP; iqp++) {
	  _lam = lam->val[iqp];
	  _mu = mu->val[iqp];
	  mu23 = _mu * (2.0/3.0);
	  mu43 = 2.0 * mu23;
	  pstress[0] = mu43 * pstrain[0] - mu23 * pstrain[1];
	  pstress[1] = mu43 * pstrain[1] - mu23 * pstrain[0];
	  pstress[2] = _mu * pstrain[2];
	  pstress += sym;
	  pstrain += sym;
	}
      }
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "dw_lin_elastic_iso"
/*!
  @par Revision history:
  - 07.03.2006, c
*/
int32 dw_lin_elastic_iso( FMField *out, FMField *strain,
			  FMField *lam, FMField *mu, Mapping *vg,
			  int32 isDiff )
{
  int32 ii, dim, sym, nQP, nEP, ret = RET_OK;
  FMField *stress = 0;
  FMField *res = 0, *d11 = 0, *gtd11 = 0, *gtd11g = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;
  dim = vg->bfGM->nRow;
  sym = (dim + 1) * dim / 2;

  if (isDiff) {
    fmf_createAlloc( &d11, 1, nQP, sym, sym );
    fmf_createAlloc( &gtd11, 1, nQP, nEP * dim, sym );
    fmf_createAlloc( &gtd11g, 1, nQP, nEP * dim, nEP * dim );

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCell( lam, ii );
      FMF_SetCell( mu, ii );
      FMF_SetCell( vg->bfGM, ii );
      FMF_SetCell( vg->det, ii );

      mat_le_tanModuli11( d11, lam, mu, 0 );
/*       fmf_print( d11, stdout, 0 ); */
/*       sys_pause(); */

      form_sdcc_actOpGT_M3( gtd11, vg->bfGM, d11 );
      form_sdcc_actOpG_RM3( gtd11g, gtd11, vg->bfGM );
      fmf_sumLevelsMulF( out, gtd11g, vg->det->val );

      ERR_CheckGo( ret );
    }
  } else {
    fmf_createAlloc( &stress, strain->nCell, nQP, sym, 1 );
    fmf_createAlloc( &res, 1, nQP, dim * nEP, 1 );

    mat_le_stress( stress, strain, lam, mu );

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCell( stress, ii );
      FMF_SetCell( vg->bfGM, ii );
      FMF_SetCell( vg->det, ii );

      form_sdcc_actOpGT_VS3( res, vg->bfGM, stress );
      fmf_sumLevelsMulF( out, res, vg->det->val );
      ERR_CheckGo( ret );
    }
  }

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &d11 );
    fmf_freeDestroy( &gtd11 );
    fmf_freeDestroy( &gtd11g );
  } else {
    fmf_freeDestroy( &res );
    fmf_freeDestroy( &stress );
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_lin_elastic"
/*!
  @par Revision history:
  - 03.08.2006, c
  - 29.11.2006
*/
int32 dw_lin_elastic( FMField *out, float64 coef, FMField *strain,
		      FMField *mtxD, Mapping *vg,
		      int32 isDiff )
{
  int32 ii, dim, sym, nQP, nEP, ret = RET_OK;
  FMField *stress = 0;
  FMField *res = 0, *gtd = 0, *gtdg = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;
  dim = vg->bfGM->nRow;
  sym = (dim + 1) * dim / 2;

/*       fmf_print( mtxD, stdout, 0 ); */
/*   output( "%d %d %d %d %d %d\n", offset, nEl, nEP, nQP, dim, elList_nRow ); */
  if (isDiff) {
    fmf_createAlloc( &gtd, 1, nQP, nEP * dim, sym );
    fmf_createAlloc( &gtdg, 1, nQP, nEP * dim, nEP * dim );

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCell( mtxD, ii );
      FMF_SetCell( vg->bfGM, ii );
      FMF_SetCell( vg->det, ii );

      form_sdcc_actOpGT_M3( gtd, vg->bfGM, mtxD );
      form_sdcc_actOpG_RM3( gtdg, gtd, vg->bfGM );
      fmf_sumLevelsMulF( out, gtdg, vg->det->val );

      ERR_CheckGo( ret );
    }
  } else {
    fmf_createAlloc( &stress, 1, nQP, sym, 1 );
    fmf_createAlloc( &res, 1, nQP, dim * nEP, 1 );

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCell( mtxD, ii );
      FMF_SetCell( vg->bfGM, ii );
      FMF_SetCell( vg->det, ii );
      FMF_SetCell( strain, ii );

      fmf_mulAB_nn( stress, mtxD, strain );
      form_sdcc_actOpGT_VS3( res, vg->bfGM, stress );
      fmf_sumLevelsMulF( out, res, vg->det->val );
      ERR_CheckGo( ret );
    }
  }

  // E.g. 1/dt.
  fmfc_mulC( out, coef );

 end_label:
  if (isDiff) {
    fmf_freeDestroy( &gtd );
    fmf_freeDestroy( &gtdg );
  } else {
    fmf_freeDestroy( &res ); 
    fmf_freeDestroy( &stress ); 
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "d_lin_elastic"
/*!
  @par Revision history:
  - 02.03.2007, c
*/
int32 d_lin_elastic( FMField *out, float64 coef, FMField *strainV,
		     FMField *strainU, FMField *mtxD, Mapping *vg )
{
  int32 ii, dim, sym, nQP, nEP, ret = RET_OK;
  FMField *std = 0, *stds = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;
  dim = vg->bfGM->nRow;
  sym = mtxD->nRow;

  fmf_createAlloc( &std, 1, nQP, 1, sym );
  fmf_createAlloc( &stds, 1, nQP, 1, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( mtxD, ii );
    FMF_SetCell( vg->det, ii );
    FMF_SetCell( strainV, ii );
    FMF_SetCell( strainU, ii );

    fmf_mulATB_nn( std, strainV, mtxD );
    fmf_mulAB_nn( stds, std, strainU );
    fmf_sumLevelsMulF( out, stds, vg->det->val );

    ERR_CheckGo( ret );
  }

  // E.g. 1/dt.
  fmfc_mulC( out, coef );

 end_label:
  fmf_freeDestroy( &std );
  fmf_freeDestroy( &stds );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_lin_prestress"
int32 dw_lin_prestress( FMField *out, FMField *stress, Mapping *vg )
{
  int32 ii, dim, nQP, nEP, ret = RET_OK;
  FMField *res = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;
  dim = vg->bfGM->nRow;

  fmf_createAlloc( &res, 1, nQP, dim * nEP, 1 );

    for (ii = 0; ii < out->nCell; ii++) {
      FMF_SetCell( out, ii );
      FMF_SetCell( vg->bfGM, ii );
      FMF_SetCell( vg->det, ii );
      FMF_SetCell( stress, ii );

      form_sdcc_actOpGT_VS3( res, vg->bfGM, stress );
      fmf_sumLevelsMulF( out, res, vg->det->val );
      ERR_CheckGo( ret );
    }

 end_label:
    fmf_freeDestroy( &res );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dw_lin_strain_fib"
int32 dw_lin_strain_fib( FMField *out, FMField *mtxD, FMField *mat,
			 Mapping *vg )
{
  int32 ii, dim, sym, nQP, nEP, ret = RET_OK;
  FMField *aux1 = 0, *aux2 = 0;

  nQP = vg->bfGM->nLev;
  nEP = vg->bfGM->nCol;
  dim = vg->bfGM->nRow;
  sym = (dim + 1) * dim / 2;

  fmf_createAlloc( &aux1, 1, nQP, nEP * dim, sym );
  fmf_createAlloc( &aux2, 1, nQP, nEP * dim, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( mtxD, ii );
    FMF_SetCell( mat, ii );
    FMF_SetCell( vg->bfGM, ii );
    FMF_SetCell( vg->det, ii );

    form_sdcc_actOpGT_M3( aux1, vg->bfGM, mtxD );
    fmf_mulAB_nn( aux2, aux1, mat );
    fmf_sumLevelsMulF( out, aux2, vg->det->val );
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &aux1 );
  fmf_freeDestroy( &aux2 );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "de_cauchy_strain"
/*!
  @par Revision history:
  - 21.09.2006, c
*/
int32 de_cauchy_strain( FMField *out, FMField *strain,
			Mapping *vg, int32 mode )
{
  int32 ii, dim, sym, nQP, ret = RET_OK;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;
  sym = (dim + 1) * dim / 2;

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( strain, ii );
    FMF_SetCell( vg->det, ii );

    fmf_sumLevelsMulF( out, strain, vg->det->val );
    if (mode == 1) {
      FMF_SetCell( vg->volume, ii );
      fmf_mulC( out, 1.0 / vg->volume->val[0] );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  return( ret );
}

#undef __FUNC__
#define __FUNC__ "de_cauchy_stress"
/*!
  @par Revision history:
  - c: 25.03.2008
*/
int32 de_cauchy_stress( FMField *out, FMField *strain,
			FMField *mtxD,  Mapping *vg,
			int32 mode )
{
  int32 ii, dim, sym, nQP, ret = RET_OK;
  FMField *stress = 0;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;
  sym = (dim + 1) * dim / 2;

  fmf_createAlloc( &stress, 1, nQP, sym, 1 );

  for (ii = 0; ii < out->nCell; ii++) {
    FMF_SetCell( out, ii );
    FMF_SetCell( mtxD, ii );
    FMF_SetCell( strain, ii );
    FMF_SetCell( vg->det, ii );

    fmf_mulAB_nn( stress, mtxD, strain );
    fmf_sumLevelsMulF( out, stress, vg->det->val );
    if (mode == 1) {
      FMF_SetCell( vg->volume, ii );
      fmf_mulC( out, 1.0 / vg->volume->val[0] );
    }
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &stress );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "dq_cauchy_strain"
/*!
  @par Revision history:
  - 30.07.2007, from dw_hdpm_cache()
*/
int32 dq_cauchy_strain( FMField *out, FMField *state, int32 offset,
			Mapping *vg,
			int32 *conn, int32 nEl, int32 nEP )
{
  int32 ii, dim, sym, nQP, ret = RET_OK;
  FMField *st = 0, *disG = 0;

  state->val = FMF_PtrFirst( state ) + offset;

  nQP = vg->bfGM->nLev;
  dim = vg->bfGM->nRow;

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

 end_label:
  fmf_freeDestroy( &st ); 
  fmf_freeDestroy( &disG ); 

  return( ret );
}
