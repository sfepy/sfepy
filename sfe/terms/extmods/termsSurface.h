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

int32 dw_surface_ltr( FMField *out, FMField *bf, FMField *gbf,
		      FMField *traction, SurfaceGeometry *sg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow );

END_C_DECLS

#endif /* Header */
