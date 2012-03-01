/*!
  @par Revision history:
  - 26.02.2007, c
*/
#ifndef _TERMSADJOINTNAVIERSTOKES_H_
#define _TERMSADJOINTNAVIERSTOKES_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"


int32 dw_adj_convect1( FMField *out, FMField *stateW, FMField *gradU,
                       FMField *bf, VolumeGeometry *vg, int32 isDiff );

int32 dw_adj_convect2( FMField *out, FMField *stateW, FMField *stateU,
                       FMField *bf, VolumeGeometry *vg, int32 isDiff );

int32 dw_st_adj_supg_c( FMField *out, FMField *stateW,
			FMField *stateU, FMField *gradU,
			FMField *coef, FMField *bf, VolumeGeometry *vg,
			int32 *conn, int32 nEl, int32 nEP,
			int32 isDiff );

int32 dw_st_adj1_supg_p( FMField *out, FMField *stateW, FMField *gradP,
			 FMField *coef, FMField *bf_w, VolumeGeometry *vg_w,
			 int32 *conn_w, int32 nEl_w, int32 nEP_w,
			 int32 isDiff );

int32 dw_st_adj2_supg_p( FMField *out, FMField *gradU, FMField *stateR,
			 FMField *coef, FMField *bf_u,
			 VolumeGeometry *vg_u, VolumeGeometry *vg_r,
			 int32 *conn_r, int32 nEl_r, int32 nEP_r,
			 int32 isDiff );

int32 d_of_nsMinGrad( FMField *out, FMField *grad,
		      FMField *viscosity, VolumeGeometry *vg );

int32 d_of_nsSurfMinDPress( FMField *out, FMField *pressure,
                            float64 weight, float64 bpress,
			    FMField *bf, SurfaceGeometry *sg, int32 isDiff );

int32 d_sd_div( FMField *out,
		FMField *stateU, int32 offsetU,
		FMField *stateP, int32 offsetP,
		FMField *vecMV, int32 offsetMV,
		FMField *bf_p,
		VolumeGeometry *vg_u,
		VolumeGeometry *vg_p,
		VolumeGeometry *vg_mv,
		int32 *conn_u, int32 nEl_u, int32 nEP_u,
		int32 *conn_p, int32 nEl_p, int32 nEP_p,
		int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
		int32 *elList, int32 elList_nRow,
		int32 mode );

int32 d_sd_div_grad( FMField *out,
		     FMField *stateU, int32 offsetU,
		     FMField *stateW, int32 offsetW,
		     FMField *vecMV, int32 offsetMV,
		     float64 viscosity,
		     VolumeGeometry *vg_u,
		     VolumeGeometry *vg_mv,
		     int32 *conn_u, int32 nEl_u, int32 nEP_u,
		     int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
		     int32 *elList, int32 elList_nRow,
		     int32 mode );

int32 d_sd_convect( FMField *out,
		    FMField *stateU, int32 offsetU,
		    FMField *stateW, int32 offsetW,
		    FMField *vecMV, int32 offsetMV,
		    FMField *bf_u, FMField *bf_w,
		    VolumeGeometry *vg_u,
		    VolumeGeometry *vg_w,
		    VolumeGeometry *vg_mv,
		    int32 *conn_u, int32 nEl_u, int32 nEP_u,
		    int32 *conn_w, int32 nEl_w, int32 nEP_w,
		    int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
		    int32 *elList, int32 elList_nRow,
		    int32 mode );

int32 d_sd_testPQ( FMField *out,
		   FMField *stateP, int32 offsetP,
		   FMField *stateQ, int32 offsetQ,
		   FMField *vecMV, int32 offsetMV,
		   FMField *bf, VolumeGeometry *vg,
		   int32 *conn, int32 nEl, int32 nEP,
		   int32 *elList, int32 elList_nRow,
		   int32 mode );

int32 d_sd_st_grad_div( FMField *out,
			FMField *stateU, int32 offsetU,
			FMField *stateW, int32 offsetW,
			FMField *vecMV, int32 offsetMV,
			float64 gamma,
			VolumeGeometry *vg_u,
			VolumeGeometry *vg_mv,
			int32 *conn_u, int32 nEl_u, int32 nEP_u,
			int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
			int32 *elList, int32 elList_nRow,
			int32 mode );

int32 d_sd_st_supg_c( FMField *out,
		      FMField *stateU, int32 offsetU,
		      FMField *stateB, int32 offsetB,
		      FMField *stateW, int32 offsetW,
		      FMField *vecMV, int32 offsetMV,
		      FMField *bf_u,
		      FMField *coef,
		      VolumeGeometry *vg_u,
		      VolumeGeometry *vg_mv,
		      int32 *conn_u, int32 nEl_u, int32 nEP_u,
		      int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
		      int32 *elList, int32 elList_nRow,
		      int32 mode );

int32 d_sd_st_pspg_c( FMField *out,
		      FMField *stateU, int32 offsetU,
		      FMField *stateB, int32 offsetB,
		      FMField *stateR, int32 offsetR,
		      FMField *vecMV, int32 offsetMV,
		      FMField *bf_u,
		      FMField *coef,
		      VolumeGeometry *vg_u,
		      VolumeGeometry *vg_r,
		      VolumeGeometry *vg_mv,
		      int32 *conn_u, int32 nEl_u, int32 nEP_u,
		      int32 *conn_r, int32 nEl_r, int32 nEP_r,
		      int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
		      int32 *elList, int32 elList_nRow,
		      int32 mode );

int32 d_sd_st_pspg_p( FMField *out,
		      FMField *stateP, int32 offsetP,
		      FMField *stateR, int32 offsetR,
		      FMField *vecMV, int32 offsetMV,
		      FMField *coef,
		      VolumeGeometry *vg_p,
		      VolumeGeometry *vg_mv,
		      int32 *conn_p, int32 nEl_p, int32 nEP_p,
		      int32 *conn_mv, int32 nEl_mv, int32 nEP_mv,
		      int32 *elList, int32 elList_nRow,
		      int32 mode );

END_C_DECLS

#endif /* Header */
