/*!
  @par Revision history:
  - 28.11.2005, c
*/
#ifndef _TERMSLAPLACE_H_
#define _TERMSLAPLACE_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"

int32 term_laplace_asm( FMField *out, FMField *state, int32 offset,
			float64 coef, VolumeGeometry *vg,
			int32 *conn, int32 nEl, int32 nEP,
			int32 *elList, int32 elList_nRow,
			int32 isDiff );

END_C_DECLS

#endif /* Header */
