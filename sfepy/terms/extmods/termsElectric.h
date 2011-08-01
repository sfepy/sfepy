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

int32 dw_electric_source( FMField *out, FMField *grad, FMField *coef,
			  FMField *bf, VolumeGeometry *vg );

END_C_DECLS

#endif /* Header */

