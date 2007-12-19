/*!
  @par Revision history:
  - 18.09.2006, c
*/
#ifndef _TERMSBASIC_H_
#define _TERMSBASIC_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"

int32 dq_state_in_qp( FMField *out, FMField *state, int32 offset,
		      FMField *bf,
		      int32 *conn, int32 nEl, int32 nEP );
int32 dq_grad_scalar( FMField *out, FMField *state, int32 offset,
		      VolumeGeometry *vg,
		      int32 *conn, int32 nEl, int32 nEP );
int32 dq_div_vector( FMField *out, FMField *state, int32 offset,
		     VolumeGeometry *vg,
		     int32 *conn, int32 nEl, int32 nEP );

END_C_DECLS

#endif /* Header */
