/* -*- C -*- */
%module terms

%{
#include "termsAdjointNavierStokes.h"
#include "termsBasic.h"
#include "termsLaplace.h"
#include "termsLinElasticity.h"
#include "termsHyperElasticity.h"
#include "termsNavierStokes.h"
#include "termsBiot.h"
#include "termsPiezo.h"
#include "termsElectric.h"
#include "termsSurface.h"
#include "termsMass.h"
#include "termsVolume.h"
#include "termsAcoustic.h"
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
    (FMField *stateV),
    (FMField *stateB),
    (FMField *stateW),
    (FMField *stateP),
    (FMField *stateQ),
    (FMField *stateR),
    (FMField *state0),
    (FMField *state1),
    (FMField *state2),
    (FMField *vecMV),
    (FMField *mtxD),
    (FMField *ref_porosity),
    (FMField *mtxG),
    (FMField *traction),
    (FMField *forceQP),
    (FMField *coef),
    (FMField *coef2),
    (FMField *strain),
    (FMField *strainV),
    (FMField *strainU),
    (FMField *gradP1),
    (FMField *gradP2),
    (FMField *pressure),
    (FMField *pressure_grad),
    (FMField *charge_grad),
    (FMField *pressure_qp),
    (FMField *state_qp),
    (FMField *lam),
    (FMField *mu),
    (FMField *viscosity),
    (FMField *stress),
    (FMField *tan_mod),
    (FMField *mtxF),
    (FMField *mtxFI),
    (FMField *detF),
    (FMField *vecCS),
    (FMField *trC),
    (FMField *in2C),
    (FMField *vecInvCS),
    (FMField *vecES),
    (FMField *vecBS),
    (FMField *trB),
    (FMField *in2B),
    (FMField *mat)
};

%apply (int32 *array, int32 n_row, int32 n_col) {
    (int32 *conn, int32 nEl, int32 nEP),
    (int32 *conn_u, int32 nEl_u, int32 nEP_u),
    (int32 *conn_w, int32 nEl_w, int32 nEP_w),
    (int32 *conn_p, int32 nEl_p, int32 nEP_p),
    (int32 *conn_r, int32 nEl_r, int32 nEP_r),
    (int32 *conn_mv, int32 nEl_mv, int32 nEP_mv),
    (int32 *conn1, int32 nEl1, int32 nEP1),
    (int32 *conn2, int32 nEl2, int32 nEP2),
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

int32 dq_grad( FMField *out, FMField *state, int32 offset,
	       VolumeGeometry *vg, int32 *conn, int32 nEl, int32 nEP );

int32 de_grad( FMField *out, FMField *state, int32 offset,
	       VolumeGeometry *vg, int32 *conn, int32 nEl, int32 nEP,
	       int32 *elList, int32 elList_nRow );

int32 dq_div_vector( FMField *out, FMField *state, int32 offset,
		     VolumeGeometry *vg,
		     int32 *conn, int32 nEl, int32 nEP );

int32 d_volume_surface( FMField *out, FMField *in,
			FMField *bf, SurfaceGeometry *sg,
			int32 *conn, int32 nEl, int32 nEP,
			int32 *elList, int32 elList_nRow );

int32 di_surface_moment( FMField *out, FMField *in,
			 FMField *bf, SurfaceGeometry *sg,
			 int32 *conn, int32 nEl, int32 nEP,
			 int32 *elList, int32 elList_nRow );

int32 dq_finite_strain_tl( FMField *mtxF, FMField *detF, FMField *vecCS,
			   FMField *trC, FMField *in2C, FMField *vecInvCS,
			   FMField *vecES,
			   FMField *state, int32 offset, VolumeGeometry *vg,
			   int32 *conn, int32 nEl, int32 nEP );

int32 dq_finite_strain_ul( FMField *mtxF, FMField *detF, FMField *vecBS,
			   FMField *trB, FMField *in2B, FMField *vecES,
			   FMField *state, FMField *state0,
			   int32 offset, VolumeGeometry *vg,
			   int32 *conn, int32 nEl, int32 nEP );

int32 dq_tl_finite_strain_surface( FMField *mtxF, FMField *detF, FMField *mtxFI,
				   FMField *state, int32 offset,
				   SurfaceGeometry *sg,
				   int32 *fis, int32 nFa, int32 nFP,
				   int32 *conn, int32 nEl, int32 nEP);

int32 dq_tl_he_stress_bulk( FMField *out,FMField *mat,
			    FMField *detF, FMField *vecInvCS );

int32 dq_ul_he_stress_bulk( FMField *out,FMField *mat,
			    FMField *detF );

int32 dq_tl_he_stress_neohook( FMField *out, FMField *mat,
			       FMField *detF, FMField *trC, FMField *vecInvCS );

int32 dq_ul_he_stress_neohook( FMField *out, FMField *mat,
			       FMField *detF, FMField *trB, FMField *vecBS );

int32 dq_tl_he_stress_mooney_rivlin( FMField *out, FMField *mat,
				     FMField *detF, FMField *trC,
				     FMField *vecInvCS, FMField *vecCS,
				     FMField *in2C );

int32 dq_ul_he_stress_mooney_rivlin( FMField *out, FMField *mat,
				     FMField *detF, FMField *trB,
				     FMField *vecBS, FMField *in2B );

int32 dq_tl_he_tan_mod_bulk( FMField *out, FMField *mat,
			     FMField *detF, FMField *vecInvCS );

int32 dq_ul_he_tan_mod_bulk( FMField *out, FMField *mat, FMField *detF );

int32 dq_tl_he_tan_mod_neohook( FMField *out, FMField *mat,
				FMField *detF, FMField *trC, FMField *vecInvCS );

int32 dq_ul_he_tan_mod_neohook( FMField *out, FMField *mat,
				FMField *detF, FMField *trB, FMField *vecBS );

int32 dq_tl_he_tan_mod_mooney_rivlin( FMField *out, FMField *mat,
				      FMField *detF, FMField *trC,
				      FMField *vecInvCS, FMField *vecCS,
				      FMField *in2C );

int32 dq_ul_he_tan_mod_mooney_rivlin( FMField *out, FMField *mat,
				      FMField *detF, FMField *trB,
				      FMField *vecBS, FMField *in2B );

int32 dw_he_rtm( FMField *out,
		 FMField *stress, FMField *tan_mod,
		 FMField *mtxF, FMField *detF,
		 VolumeGeometry *vg,
		 int32 *elList, int32 elList_nRow, int32 isDiff, int32 mode_ul );

int32 dq_tl_stress_bulk_pressure( FMField *out, FMField *pressure_qp,
				  FMField *detF, FMField *vecInvCS );
int32 dq_tl_tan_mod_bulk_pressure_u( FMField *out, FMField *pressure_qp,
				     FMField *detF, FMField *vecInvCS );

int32 dw_tl_volume( FMField *out, FMField *bf, FMField *mtxF,
		    FMField *vecInvCS, FMField *detF,
		    VolumeGeometry *vg, int32 transpose,
		    int32 *elList, int32 elList_nRow,
		    int32 mode );

int32 dw_tl_diffusion( FMField *out, FMField *pressure_grad,
		       FMField *mtxD, FMField *ref_porosity,
		       FMField *mtxF, FMField *detF,
		       VolumeGeometry *vg,
		       int32 *elList, int32 elList_nRow,
		       int32 mode );

int32 dw_tl_surface_traction( FMField *out, FMField *traction,
			      FMField *detF, FMField *mtxFI,
			      FMField *bf, SurfaceGeometry *sg,
			      int32 *fis, int32 nFa, int32 nFP,
			      int32 *elList, int32 elList_nRow,
			      int32 mode );

int32 dq_def_grad( FMField *out, FMField *state, VolumeGeometry *vg,
		   int32 *conn, int32 nEl, int32 nEP,
		   int32 *elList, int32 elList_nRow, int32 mode );

int32 dw_volume_wdot_scalar( FMField *out, float64 coef, FMField *state_qp,
			     FMField *bf, FMField *mtxD, VolumeGeometry *vg,
			     int32 *elList, int32 elList_nRow,
			     int32 isDiff );

int32 dw_laplace( FMField *out, FMField *state, int32 offset,
		  FMField *coef, VolumeGeometry *vg,
		  int32 *conn, int32 nEl, int32 nEP,
		  int32 *elList, int32 elList_nRow,
		  int32 isDiff );
int32 d_laplace( FMField *out, FMField *gradP1, FMField *gradP2,
		 FMField *coef, VolumeGeometry *vg,
		 int32 *elList, int32 elList_nRow );
int32 dw_diffusion( FMField *out, FMField *state, int32 offset,
		    FMField *mtxD, VolumeGeometry *vg,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 *elList, int32 elList_nRow,
		    int32 isDiff );
int32 d_diffusion( FMField *out, FMField *gradP1, FMField *gradP2,
		   FMField *mtxD, VolumeGeometry *vg,
		   int32 *elList, int32 elList_nRow );
int32 dw_permeability_r( FMField *out, FMField *mtxD, VolumeGeometry *vg,
			 int32 *conn, int32 nEl, int32 nEP,
			 int32 *elList, int32 elList_nRow );
int32 dw_diffusion_coupling( FMField *out, FMField *state, int32 offset,
			     FMField *mtxD, FMField *bf, VolumeGeometry *vg,
			     int32 *conn, int32 nEl, int32 nEP,
			     int32 *elList, int32 elList_nRow,
			     int32 isDiff, int32 mode);
int32 d_diffusion_coupling( FMField *out, FMField *stateP, FMField *stateQ,
			    FMField *mtxD, FMField *bf, VolumeGeometry *vg,
			    int32 *conn, int32 nEl, int32 nEP,
			    int32 *elList, int32 elList_nRow,
			    int32 isDiff, int32 mode);
int32 de_diffusion_velocity( FMField *out, FMField *state, int32 offset,
			     FMField *mtxD, VolumeGeometry *vg,
			     int32 *conn, int32 nEl, int32 nEP,
			     int32 *elList, int32 elList_nRow );
int32 d_diffusion_integrate( FMField *out, FMField *in,
			     FMField *mtxD, VolumeGeometry *vg,
			     int32 *conn, int32 nEl, int32 nEP,
			     int32 *elList, int32 elList_nRow );
int32 d_surf_diffusion_integrate( FMField *out, FMField *in,
				  FMField *mtxD, SurfaceGeometry *sg,
				  int32 *conn, int32 nEl, int32 nEP,
				  int32 *elList, int32 elList_nRow );

int32 dw_lin_elastic_iso( FMField *out, FMField *state, int32 offset,
			  FMField *lam, FMField *mu, VolumeGeometry *vg,
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

int32 dw_lin_prestress( FMField *out, FMField *stress, VolumeGeometry *vg,
			int32 *elList, int32 elList_nRow, int32 isDiff );

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
		      int32 *elList, int32 elList_nRow );

int32 dw_jump( FMField *out, FMField *coef, FMField *state1, FMField *state2,
	       FMField *bf, SurfaceGeometry *sg,
	       int32 *conn1, int32 nEl1, int32 nEP1,
	       int32 *conn2, int32 nEl2, int32 nEP2,
	       int32 *elList, int32 elList_nRow, int32 mode );

int32 dw_volume_lvf( FMField *out, FMField *bf, FMField *forceQP,
		     VolumeGeometry *vg, int32 *elList, int32 elList_nRow );

int32 dw_mass( FMField *out, FMField *coef, FMField *state, int32 offset,
	       FMField *bf, VolumeGeometry *vg,
	       int32 *conn, int32 nEl, int32 nEP,
	       int32 *elList, int32 elList_nRow,
	       int32 isDiff );

int32 dw_mass_scalar( FMField *out, FMField *coef,
		      FMField *state, FMField *bf, VolumeGeometry *vg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow,
		      int32 isDiff );

int32 d_mass_scalar( FMField *out, FMField *coef,
		     FMField *stateP, FMField *stateQ,
		     FMField *bf, VolumeGeometry *vg,
		     int32 *conn, int32 nEl, int32 nEP,
		     int32 *elList, int32 elList_nRow );

int32 dw_surf_mass_scalar( FMField *out, FMField *coef,
			   FMField *state, int32 offset,
			   FMField *bf, SurfaceGeometry *sg,
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

int32 term_ns_asm_div_grad( FMField *out, FMField *state, int32 offset,
			    FMField *viscosity, VolumeGeometry *vg,
			    int32 *conn, int32 nEl, int32 nEP,
			    int32 *elList, int32 elList_nRow,
			    int32 isDiff );


int32 term_ns_asm_convect( FMField *out, FMField *state, int32 offset,
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

int32 dw_div( FMField *out, FMField *coef, FMField *state, int32 offset,
	      FMField *bf, VolumeGeometry *vg,
	      int32 *conn, int32 nEl, int32 nEP,
	      int32 *elList, int32 elList_nRow,
	      int32 isDiff );

int32 dw_grad( FMField *out, FMField *coef, FMField *state, int32 offset,
	       FMField *bf, VolumeGeometry *vg,
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
		      FMField *coef, VolumeGeometry *vg,
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

int32 dw_piezo_coupling( FMField *out, FMField *strain, FMField *charge_grad,
			 FMField *mtxG, VolumeGeometry *vg,
			 int32 *elList, int32 elList_nRow,
			 int32 mode );

int32 d_piezo_coupling( FMField *out, FMField *strain, FMField *charge_grad,
			FMField *mtxG, VolumeGeometry *vg,
			int32 *elList, int32 elList_nRow );

int32 dw_electric_source( FMField *out,
			  FMField *state,
			  FMField *coef, FMField *bf,
			  VolumeGeometry *vgc,
			  int32 *conn, int32 nEl, int32 nEP,
			  int32 *elList, int32 elList_nRow,
			  int32 mode );

int32 d_diffusion_sa( FMField *out,
		      FMField *stateQ, FMField *stateP, FMField *stateW,
		      FMField *mtxD,
		      VolumeGeometry *vg, VolumeGeometry *vg_w,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *conn_w, int32 nEl_w, int32 nEP_w,
		      int32 *elList, int32 elList_nRow );

int32 dw_surf_laplace( FMField *out, FMField *state, FMField *coef,
		       FMField *gbf, SurfaceGeometry *sg,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow,
		       int32 isDiff );

int32 d_surf_laplace( FMField *out, FMField *stateP, FMField *stateQ, FMField *coef,
		      FMField *gbf, SurfaceGeometry *sg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow);

int32 dw_surf_lcouple( FMField *out, FMField *state, FMField *coef,
		       FMField *bf, FMField *gbf, SurfaceGeometry *sg,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow,
		       int32 isDiff );

int32 d_surf_lcouple(FMField *out, FMField *stateP, FMField *stateQ, FMField *coef,
		     FMField *bf, FMField *gbf, SurfaceGeometry *sg,
		     int32 *conn, int32 nEl, int32 nEP,
		     int32 *elList, int32 elList_nRow);

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

int32 d_hdpm_surfdvel( FMField *out, FMField *state, int32 offset,
		       FMField *mtxD, SurfaceGeometry *sg,
		       int32 *fis, int32 nFa, int32 nFP,
		       int32 *faceList, int32 faceList_nRow,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow );
