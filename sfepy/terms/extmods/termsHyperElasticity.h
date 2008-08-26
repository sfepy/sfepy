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

int32 dq_tl_he_stress_neohook( FMField *out, FMField *mat,
			       FMField *detF, FMField *trC, FMField *vecInvCS );

int32 dq_tl_he_stress_bulk( FMField *out,FMField *mat,
			    FMField *detF, FMField *vecInvCS );

int32 dq_tl_he_tan_mod_neohook( FMField *out, FMField *mat,
				FMField *detF, FMField *trC, FMField *vecInvCS );

int32 dq_tl_he_tan_mod_bulk( FMField *out, FMField *mat,
			     FMField *detF, FMField *vecInvCS );

int32 dw_tl_he_rtm( FMField *out,
		    FMField *stress, FMField *tan_mod, FMField *mtxF,
		    VolumeGeometry *vg,
		    int32 *elList, int32 elList_nRow, int32 isDiff );

END_C_DECLS

#endif /* Header */
