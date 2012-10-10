/*!
  @par Revision history:
  - 28.11.2005, c
*/
#ifndef _TERMSLINELASTICITY_H_
#define _TERMSLINELASTICITY_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "refmaps.h"

int32 dw_lin_elastic_iso( FMField *out, FMField *strain,
			  FMField *lam, FMField *mu, Mapping *vg,
			  int32 isDiff );
int32 dw_lin_elastic( FMField *out, float64 coef, FMField *strain,
		      FMField *mtxD, Mapping *vg,
		      int32 isDiff );
int32 d_lin_elastic( FMField *out, float64 coef, FMField *strainV,
		     FMField *strainU, FMField *mtxD, Mapping *vg );

int32 dw_lin_prestress( FMField *out, FMField *stress, Mapping *vg );

int32 dw_lin_strain_fib( FMField *out, FMField *mtxD, FMField *mat,
			 Mapping *vg );

int32 de_cauchy_strain( FMField *out, FMField *strain,
			Mapping *vg, int32 mode );
int32 de_cauchy_stress( FMField *out, FMField *strain,
			FMField *mtxD,  Mapping *vg,
			int32 mode );
int32 dq_cauchy_strain( FMField *out, FMField *state, int32 offset,
			Mapping *vg,
			int32 *conn, int32 nEl, int32 nEP );


END_C_DECLS

#endif /* Header */
