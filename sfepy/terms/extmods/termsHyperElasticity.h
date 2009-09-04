#ifndef _TERMSHYPERELASTICITY_H_
#define _TERMSHYPERELASTICITY_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"

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

END_C_DECLS

#endif /* Header */
