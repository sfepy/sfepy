/* -*- C -*- */
%module terms

%{
#include "termsBasic.h"
#include "termsLaplace.h"
#include "termsLinElasticity.h"
#include "termsNavierStokes.h"
#include "termsBiot.h"
#include "termsSurface.h"
#include "termsMass.h"
#include "termsVolume.h"
#ifndef ISRELEASE
  #include "termsAdjointNavierStokes.h"
  #include "termsHDPM.h"
#endif
%}

%include "types.h"

%include common.i
%include array.i
%include fmfield.i

%apply (FMField *in) {
    (FMField *out),
    (FMField *state),
    (FMField *velocity),
    (FMField *gbf),
    (FMField *cbfs),
    (FMField *bf),
    (FMField *bf_p),
    (FMField *bf_u),
    (FMField *bf_w),
    (FMField *meshVelocity),
    (FMField *diameter2),
    (FMField *stateU),
    (FMField *stateB),
    (FMField *stateW),
    (FMField *stateP),
    (FMField *stateQ),
    (FMField *stateR),
    (FMField *vecMV),
    (FMField *mtxD),
    (FMField *traction),
    (FMField *force),
    (FMField *coef),
    (FMField *strain),
    (FMField *strainV),
    (FMField *strainU),
    (FMField *gradP1),
    (FMField *gradP2),
    (FMField *pressure),
    (FMField *pressure_qp)
};

%apply (int32 *array, int32 nRow, int32 nCol) {
    (int32 *conn, int32 nEl, int32 nEP),
    (int32 *conn_u, int32 nEl_u, int32 nEP_u),
    (int32 *conn_w, int32 nEl_w, int32 nEP_w),
    (int32 *conn_p, int32 nEl_p, int32 nEP_p),
    (int32 *conn_r, int32 nEl_r, int32 nEP_r),
    (int32 *conn_mv, int32 nEl_mv, int32 nEP_mv),
    (int32 *fis, int32 nFa, int32 nFP)
};
%apply (int32 *array, int32 len) {
    (int32 *elList, int32 elList_nRow),
    (int32 *iemap, int32 iemap_nRow),
    (int32 *faceList, int32 faceList_nRow)
};

int32 dq_state_in_qp( FMField *out, FMField *state, int32 offset,
		      FMField *bf,
		      int32 *conn, int32 nEl, int32 nEP );

int32 dq_grad_scalar( FMField *out, FMField *state, int32 offset,
		      VolumeGeometry *vg,
		      int32 *conn, int32 nEl, int32 nEP );

int32 dq_div_vector( FMField *out, FMField *state, int32 offset,
		     VolumeGeometry *vg,
		     int32 *conn, int32 nEl, int32 nEP );

int32 dw_laplace( FMField *out, FMField *state, int32 offset,
		  float64 coef, VolumeGeometry *vg,
		  int32 *conn, int32 nEl, int32 nEP,
		  int32 *elList, int32 elList_nRow,
		  int32 isDiff );
int32 dw_diffusion( FMField *out, float64 coef, FMField *state, int32 offset,
		    FMField *mtxD, VolumeGeometry *vg,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 *elList, int32 elList_nRow,
		    int32 isDiff );
int32 d_diffusion( FMField *out, float64 coef, FMField *gradP1, FMField *gradP2,
		   FMField *mtxD, VolumeGeometry *vg,
		   int32 *elList, int32 elList_nRow );
int32 dw_permeability_r( FMField *out, FMField *mtxD, VolumeGeometry *vg,
			 int32 *conn, int32 nEl, int32 nEP,
			 int32 *elList, int32 elList_nRow );
int32 de_diffusion_velocity( FMField *out, FMField *state, int32 offset,
			     FMField *mtxD, VolumeGeometry *vg,
			     int32 *conn, int32 nEl, int32 nEP,
			     int32 *elList, int32 elList_nRow );

int32 dw_lin_elastic_iso( FMField *out, FMField *state, int32 offset,
			  float64 lam, float64 mu, VolumeGeometry *vg,
			  int32 *conn, int32 nEl, int32 nEP,
			  int32 *elList, int32 elList_nRow,
			  int32 isDiff );
int32 dw_lin_elastic( FMField *out, float64 coef, FMField *strain,
		      FMField *mtxD, VolumeGeometry *vg,
		      int32 *elList, int32 elList_nRow,
		      int32 isDiff );
int32 d_lin_elastic( FMField *out, float64 coef, FMField *strainV,
		     FMField *strainU, FMField *mtxD, VolumeGeometry *vg,
		     int32 *elList, int32 elList_nRow );

int32 de_cauchy_strain( FMField *out, FMField *state, int32 offset,
			VolumeGeometry *vg,
			int32 *conn, int32 nEl, int32 nEP,
			int32 *elList, int32 elList_nRow );
int32 de_cauchy_stress( FMField *out, FMField *strain,
			FMField *mtxD,  VolumeGeometry *vg,
			int32 *elList, int32 elList_nRow );
int32 dq_cauchy_strain( FMField *out, FMField *state, int32 offset,
			VolumeGeometry *vg,
			int32 *conn, int32 nEl, int32 nEP );

int32 dw_surface_ltr( FMField *out, FMField *bf, FMField *gbf,
		      FMField *traction, SurfaceGeometry *sg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow );

int32 dw_volume_lvf( FMField *out, FMField *bf, FMField *gbf,
		     FMField *force, VolumeGeometry *vg,
		     int32 *conn, int32 nEl, int32 nEP,
		     int32 *elList, int32 elList_nRow );

int32 dw_mass( FMField *out, float64 coef, FMField *state, int32 offset,
	       FMField *bf, VolumeGeometry *vg,
	       int32 *conn, int32 nEl, int32 nEP,
	       int32 *elList, int32 elList_nRow,
	       int32 isDiff );

int32 dw_mass_scalar( FMField *out, FMField *state, int32 offset,
		      FMField *bf, VolumeGeometry *vg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow,
		      int32 isDiff );

int32 dw_mass_scalar_variable( FMField *out, FMField *coef,
			       FMField *state, int32 offset,
			       FMField *bf, VolumeGeometry *vg,
			       int32 *conn, int32 nEl, int32 nEP,
			       int32 *elList, int32 elList_nRow,
			       int32 isDiff );

int32 dw_mass_scalar_fine_coarse( FMField *out, FMField *state, int32 offset,
				  FMField *bf, FMField *cbfs,
				  VolumeGeometry *vg,
				  int32 *conn, int32 nEl, int32 nEP,
				  int32 *iemap, int32 iemap_nRow,
				  int32 *elList, int32 elList_nRow,
				  int32 isDiff );

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

int32 dw_biot_grad( FMField *out, float64 coef, FMField *pressure_qp,
		    FMField *bf, FMField *mtxD, VolumeGeometry *vg,
		    int32 *elList, int32 elList_nRow,
		    int32 isDiff );

int32 dw_biot_div( FMField *out, float64 coef, FMField *strain,
		   FMField *bf, FMField *mtxD, VolumeGeometry *vg,
		   int32 *elList, int32 elList_nRow,
		   int32 isDiff );

int32 d_biot_div( FMField *out, float64 coef, FMField *state, FMField *strain,
		  FMField *mtxD, VolumeGeometry *vg,
		  int32 *elList, int32 elList_nRow );

#ifndef ISRELEASE

int32 dw_adj_convect1( FMField *out, FMField *state, int32 offset,
		       FMField *velocity, int32 voffset, FMField *bf,
		       VolumeGeometry *vg,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow,
		       int32 isDiff );

int32 dw_adj_convect2( FMField *out, FMField *state, int32 offset,
		       FMField *velocity, int32 voffset, FMField *bf,
		       VolumeGeometry *vg,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow,
		       int32 isDiff );

int32 dw_st_adj_supg_c( FMField *out,
			FMField *stateU, int32 offsetU,
			FMField *stateW, int32 offsetW,
			FMField *coef, FMField *bf, VolumeGeometry *vg,
			int32 *conn, int32 nEl, int32 nEP,
			int32 *elList, int32 elList_nRow,
			int32 isDiff );

int32 dw_st_adj1_supg_p( FMField *out,
			 FMField *stateP, int32 offsetP,
			 FMField *stateW, int32 offsetW,
			 FMField *coef, FMField *bf_w,
			 VolumeGeometry *vg_w, VolumeGeometry *vg_p,
			 int32 *conn_w, int32 nEl_w, int32 nEP_w,
			 int32 *conn_p, int32 nEl_p, int32 nEP_p,
			 int32 *elList, int32 elList_nRow,
			 int32 isDiff );

int32 dw_st_adj2_supg_p( FMField *out,
			 FMField *stateU, int32 offsetU,
			 FMField *stateR, int32 offsetR,
			 FMField *coef, FMField *bf_u,
			 VolumeGeometry *vg_u, VolumeGeometry *vg_r,
			 int32 *conn_u, int32 nEl_u, int32 nEP_u,
			 int32 *conn_r, int32 nEl_r, int32 nEP_r,
			 int32 *elList, int32 elList_nRow,
			 int32 isDiff );

int32 d_of_nsMinGrad( FMField *out, FMField *velocity, int32 offset,
		      float64 viscosity, VolumeGeometry *vg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow );

int32 d_of_nsSurfMinDPress( FMField *out, FMField *pressure, int32 offset,
			    float64 weight, float64 bpress,
			    FMField *bf, SurfaceGeometry *sg,
			    int32 *conn, int32 nEl, int32 nEP,
			    int32 *elList, int32 elList_nRow, int32 isDiff );

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

int32 dw_hdpm_g( FMField *out, float64 coef, FMField *pressure_qp,
		 FMField *bf, FMField *mtxD, VolumeGeometry *vg,
		 int32 *elList, int32 elList_nRow,
		 int32 isDiff );

int32 de_hdpm_dvel( FMField *out, FMField *state, int32 offset,
		    FMField *mtxD, VolumeGeometry *vg,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 *elList, int32 elList_nRow );

int32 d_hdpm_surfdvel( FMField *out, FMField *state, int32 offset,
		       FMField *mtxD, SurfaceGeometry *sg,
		       int32 *fis, int32 nFa, int32 nFP,
		       int32 *faceList, int32 faceList_nRow,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow );

#endif // ISRELEASE
