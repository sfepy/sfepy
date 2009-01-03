/*!
  @par Revision history:
  - 03.01.2009, c
*/
#ifndef _TERMSELECTRIC_H_
#define _TERMSELECTRIC_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"

int32 dw_electric_source( FMField *out,
			  FMField *state,
			  FMField *coef, FMField *bf,
			  VolumeGeometry *vgc,
			  int32 *conn, int32 nEl, int32 nEP,
			  int32 *elList, int32 elList_nRow,
			  int32 mode );

END_C_DECLS

#endif /* Header */

