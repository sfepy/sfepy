/*!
  @par Revision history:
  - 26.10.2005, c
*/
#ifndef _TERMSNAVIERSTOKES_H_
#define _TERMSNAVIERSTOKES_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"

int32 divgrad_build_gtg( FMField *out, FMField *gc );
int32 divgrad_act_g_m( FMField *out, FMField *gc, FMField *mtx );
int32 divgrad_act_gt_m( FMField *out, FMField *gc, FMField *mtx );
int32 divgrad_act_bg_m( FMField *out, FMField *gc, FMField *mtx );
int32 convect_build_vtbg( FMField *out, FMField *gc, FMField *fv );
int32 convect_build_vtg( FMField *out, FMField *gc, FMField *fv );

int32 term_ns_asm_div_grad( FMField *out, FMField *grad,
			    FMField *viscosity, VolumeGeometry *vg,
			    int32 isDiff );

int32 term_ns_asm_convect( FMField *out, FMField *grad, FMField *state,
                           VolumeGeometry *vg, int32 isDiff );

int32 dw_lin_convect( FMField *out, FMField *grad, FMField *stateB,
		      VolumeGeometry *vg, int32 isDiff );

int32 dw_div( FMField *out, FMField *coef, FMField *div,
	      VolumeGeometry *svg, VolumeGeometry *vvg, int32 isDiff );

int32 dw_grad( FMField *out, FMField *coef, FMField *state,
	       VolumeGeometry *svg, VolumeGeometry *vvg, int32 isDiff );

int32 dw_st_pspg_c( FMField *out,
		    FMField *stateB, FMField *stateU,
		    FMField *coef,
		    VolumeGeometry *vg_p, VolumeGeometry *vg_u,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 isDiff );

int32 dw_st_supg_p( FMField *out,
		    FMField *stateB, FMField *gradP,
		    FMField *coef,
		    VolumeGeometry *vg_u, VolumeGeometry *vg_p,
		    int32 isDiff );

int32 dw_st_supg_c( FMField *out,
		    FMField *stateB, FMField *stateU,
		    FMField *coef, VolumeGeometry *vg,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 isDiff );

int32 dw_st_grad_div( FMField *out, FMField *div,
		      FMField *coef, VolumeGeometry *vg,
		      int32 isDiff );

END_C_DECLS

#endif /* Header */
