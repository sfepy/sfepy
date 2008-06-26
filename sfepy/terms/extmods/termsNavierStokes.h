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

int32 term_ns_asmDivGrad( FMField *out, FMField *state, int32 offset,
			  float64 viscosity, VolumeGeometry *vg,
			  int32 *conn, int32 nEl, int32 nEP,
			  int32 *elList, int32 elList_nRow,
			  int32 isDiff );


int32 term_ns_asmConvect( FMField *out, FMField *state, int32 offset,
			  FMField *bf, VolumeGeometry *vg,
			  int32 *conn, int32 nEl, int32 nEP,
			  int32 *elList, int32 elList_nRow,
			  int32 isDiff );

int32 dw_lin_convect( FMField *out,
		      FMField *stateB, int32 offsetB,
		      FMField *stateU, int32 offsetU,
		      FMField *bf, VolumeGeometry *vg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow,
		      int32 isDiff );

int32 dw_div( FMField *out, FMField *state, int32 offset,
	      FMField *bf, VolumeGeometry *vg,
	      int32 *conn, int32 nEl, int32 nEP,
	      int32 *elList, int32 elList_nRow,
	      int32 isDiff );

int32 dw_grad( FMField *out, float64 coef, FMField *state, int32 offset,
	       FMField *bf, VolumeGeometry *vg,
	       int32 *conn, int32 nEl, int32 nEP,
	       int32 *elList, int32 elList_nRow,
	       int32 isDiff );

int32 dw_st_pspg_p( FMField *out, FMField *state, int32 offset,
		    FMField *coef, VolumeGeometry *vg,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 *elList, int32 elList_nRow,
		    int32 isDiff );

int32 dw_st_pspg_c( FMField *out,
		    FMField *stateB, int32 offsetB,
		    FMField *stateU, int32 offsetU,
		    FMField *coef, FMField *bf_u,
		    VolumeGeometry *vg_p, VolumeGeometry *vg_u,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 *elList, int32 elList_nRow,
		    int32 isDiff );

int32 dw_st_supg_p( FMField *out,
		    FMField *stateB, int32 offsetB,
		    FMField *stateP, int32 offsetP,
		    FMField *coef, FMField *bf_u,
		    VolumeGeometry *vg_u, VolumeGeometry *vg_p,
		    int32 *conn_u, int32 nEl_u, int32 nEP_u,
		    int32 *conn_p, int32 nEl_p, int32 nEP_p,
		    int32 *elList, int32 elList_nRow,
		    int32 isDiff );

int32 dw_st_supg_c( FMField *out,
		    FMField *stateB, int32 offsetB,
		    FMField *stateU, int32 offsetU,
		    FMField *coef, FMField *bf, VolumeGeometry *vg,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 *elList, int32 elList_nRow,
		    int32 isDiff );

int32 dw_st_grad_div( FMField *out, FMField *state, int32 offset,
		      float64 gamma, VolumeGeometry *vg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow,
		      int32 isDiff );

END_C_DECLS

#endif /* Header */
