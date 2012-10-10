/*!
  @par Revision history:
  - 26.02.2007, c
*/
#ifndef _TERMSADJOINTNAVIERSTOKES_H_
#define _TERMSADJOINTNAVIERSTOKES_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "refmaps.h"


int32 dw_adj_convect1( FMField *out, FMField *stateW, FMField *gradU,
                       Mapping *vg, int32 isDiff );

int32 dw_adj_convect2( FMField *out, FMField *stateW, FMField *stateU,
                       Mapping *vg, int32 isDiff );

int32 dw_st_adj_supg_c( FMField *out, FMField *stateW,
			FMField *stateU, FMField *gradU,
			FMField *coef, Mapping *vg,
			int32 *conn, int32 nEl, int32 nEP,
			int32 isDiff );

int32 dw_st_adj1_supg_p( FMField *out, FMField *stateW, FMField *gradP,
			 FMField *coef, Mapping *vg_w,
			 int32 *conn_w, int32 nEl_w, int32 nEP_w,
			 int32 isDiff );

int32 dw_st_adj2_supg_p( FMField *out, FMField *gradU, FMField *stateR,
			 FMField *coef,
			 Mapping *vg_u, Mapping *vg_r,
			 int32 *conn_r, int32 nEl_r, int32 nEP_r,
			 int32 isDiff );

int32 d_of_nsMinGrad( FMField *out, FMField *grad,
		      FMField *viscosity, Mapping *vg );

int32 d_of_nsSurfMinDPress( FMField *out, FMField *pressure,
                            float64 weight, float64 bpress,
			    Mapping *sg, int32 isDiff );

int32 d_sd_div( FMField *out, FMField *divU, FMField *gradU,
                FMField *stateP, FMField *divMV, FMField *gradMV,
                Mapping *vg_u, int32 mode );

int32 d_sd_div_grad( FMField *out, FMField *gradU, FMField *gradW,
                     FMField *divMV, FMField *gradMV, FMField *viscosity,
		     Mapping *vg_u, int32 mode );

int32 d_sd_convect( FMField *out, FMField *stateU, FMField *gradU,
		    FMField *stateW, FMField *divMV, FMField *gradMV,
		    Mapping *vg_u, int32 mode );

int32 d_sd_volume_dot( FMField *out, FMField *stateP, FMField *stateQ,
                       FMField *divMV, Mapping *vg, int32 mode );

int32 d_sd_st_grad_div( FMField *out, FMField *divU, FMField *gradU,
			FMField *divW, FMField *gradW, FMField *divMV,
			FMField *gradMV, FMField *coef,
			Mapping *vg_u, int32 mode );

int32 d_sd_st_supg_c( FMField *out, FMField *stateB, FMField *gradU,
                      FMField *gradW, FMField *divMV, FMField *gradMV,
		      FMField *coef, Mapping *vg_u, int32 mode );

int32 d_sd_st_pspg_c( FMField *out, FMField *stateB, FMField *gradU,
		      FMField *gradR, FMField *divMV, FMField *gradMV,
		      FMField *coef, Mapping *vg_u, int32 mode );

int32 d_sd_st_pspg_p( FMField *out, FMField *gradR, FMField *gradP,
		      FMField *divMV, FMField *gradMV, FMField *coef,
		      Mapping *vg_p, int32 mode );

END_C_DECLS

#endif /* Header */
