#ifndef _TERMSHYPERELASTICITY_H_
#define _TERMSHYPERELASTICITY_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "refmaps.h"

int32 dq_finite_strain_tl( FMField *mtxF, FMField *detF, FMField *vecCS,
			   FMField *trC, FMField *in2C, FMField *vecInvCS,
			   FMField *vecES,
			   FMField *state, int32 offset, Mapping *vg,
			   int32 *conn, int32 nEl, int32 nEP );

int32 dq_finite_strain_ul( FMField *mtxF, FMField *detF, FMField *vecBS,
			   FMField *trB, FMField *in2B, FMField *vecES,
			   FMField *state, int32 offset, Mapping *vg,
			   int32 *conn, int32 nEl, int32 nEP );

int32 dq_tl_finite_strain_surface( FMField *mtxF, FMField *detF, FMField *mtxFI,
				   FMField *state, int32 offset,
				   Mapping *sg,
				   int32 *fis, int32 nFa, int32 nFP,
				   int32 *conn, int32 nEl, int32 nEP);

int32 dq_tl_he_stress_bulk( FMField *out,FMField *mat,
			    FMField *detF, FMField *vecInvCS );

int32 dq_ul_he_stress_bulk( FMField *out,FMField *mat,
			    FMField *detF );

int32 dq_tl_he_stress_bulk_active( FMField *out,FMField *mat,
                                   FMField *detF, FMField *vecCG );

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

int32 dq_tl_he_tan_mod_bulk_active( FMField *out, FMField *mat,
                                    FMField *detF, FMField *vecInvCS );

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
		 Mapping *vg,
		 int32 isDiff, int32 mode_ul );

int32 de_he_rtm( FMField *out,
		 FMField *stress, FMField *detF,
		 Mapping *vg,
		 int32 *elList, int32 elList_nRow,
		 int32 mode_ul );

int32 dq_tl_stress_bulk_pressure( FMField *out, FMField *pressure_qp,
				  FMField *detF, FMField *vecInvCS );
int32 dq_ul_stress_bulk_pressure( FMField *out, FMField *pressure_qp,
				  FMField *detF);
int32 dq_tl_tan_mod_bulk_pressure_u( FMField *out, FMField *pressure_qp,
				     FMField *detF, FMField *vecInvCS );
int32 dq_ul_tan_mod_bulk_pressure_u( FMField *out, FMField *pressure_qp,
				     FMField *detF );

int32 dw_tl_volume( FMField *out, FMField *mtxF,
		    FMField *vecInvCS, FMField *detF,
		    Mapping *vgs, Mapping *vgv,
                    int32 transpose, int32 mode );
int32 dw_ul_volume( FMField *out, FMField *detF,
		    Mapping *vgs, Mapping *vgv,
                    int32 transpose, int32 mode );

int32 dw_tl_diffusion( FMField *out, FMField *pressure_grad,
		       FMField *mtxD, FMField *ref_porosity,
		       FMField *mtxF, FMField *detF,
		       Mapping *vg, int32 mode );

int32 d_tl_surface_flux( FMField *out, FMField *pressure_grad,
                         FMField *mtxD, FMField *ref_porosity,
                         FMField *mtxFI, FMField *detF,
                         Mapping *sg, int32 mode );

int32 dw_tl_surface_traction( FMField *out, FMField *traction,
			      FMField *detF, FMField *mtxFI,
			      FMField *bf, Mapping *sg,
			      int32 *fis, int32 nFa, int32 nFP,
			      int32 mode );

int32 d_tl_volume_surface( FMField *out, FMField *coors,
                           FMField *detF, FMField *mtxFI,
                           FMField *bf, Mapping *sg,
                           int32 *conn, int32 nFa, int32 nFP );

int32 dq_def_grad( FMField *out, FMField *state, Mapping *vg,
		   int32 *conn, int32 nEl, int32 nEP, int32 mode );

int32 he_residuum_from_mtx(FMField *out, FMField *mtxD,
			   FMField *state,
			   int32 *conn, int32 nEl, int32 nEP,
			   int32 *elList, int32 elList_nRow);
int32 he_eval_from_mtx(FMField *out, FMField *mtxD,
		       FMField *stateV, FMField *stateU,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow);

END_C_DECLS

#endif /* Header */
