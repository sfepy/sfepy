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

int32 dw_volume_lvf( FMField *out, FMField *bf, FMField *forceQP,
		     VolumeGeometry *vg, int32 *elList, int32 elList_nRow );

END_C_DECLS

#endif /* Header */
