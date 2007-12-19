/*!
  @par Revision history:
  - 18.09.2006, c
*/
#ifndef _TERMSVOLUME_H_
#define _TERMSVOLUME_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"

int32 dw_volume_lvf( FMField *out, FMField *bf, FMField *gbf,
		     FMField *force, VolumeGeometry *vg,
		     int32 *conn, int32 nEl, int32 nEP,
		     int32 *elList, int32 elList_nRow );

END_C_DECLS

#endif /* Header */
