/*!
  @par Revision history:
  - 06.09.2006, c
*/
#ifndef _TERMSSURFACE_H_
#define _TERMSSURFACE_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"

int32 dw_surface_ltr( FMField *out, FMField *bf,
		      FMField *traction, SurfaceGeometry *sg );

END_C_DECLS

#endif /* Header */
