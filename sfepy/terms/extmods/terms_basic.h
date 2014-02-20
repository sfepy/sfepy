/*!
  @par Revision history:
  - 18.09.2006, c
*/
#ifndef _TERMSBASIC_H_
#define _TERMSBASIC_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "refmaps.h"

int32 dq_state_in_qp( FMField *out, FMField *state, int32 offset,
		      FMField *bf,
		      int32 *conn, int32 nEl, int32 nEP );
int32 dq_grad( FMField *out, FMField *state, int32 offset,
	       Mapping *vg, int32 *conn, int32 nEl, int32 nEP );

int32 dq_div_vector( FMField *out, FMField *state, int32 offset,
		     Mapping *vg,
		     int32 *conn, int32 nEl, int32 nEP );

int32 d_volume_surface( FMField *out, FMField *in,
			Mapping *sg,
			int32 *conn, int32 nEl, int32 nEP );

int32 di_surface_moment( FMField *out, FMField *in,
			 Mapping *sg,
			 int32 *conn, int32 nEl, int32 nEP );

END_C_DECLS

#endif /* Header */
