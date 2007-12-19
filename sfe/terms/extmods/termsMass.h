/*!
  @par Revision history:
  - 21.11.2006, c
*/
#ifndef _TERMSMASS_H_
#define _TERMSMASS_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"

int32 dw_mass( FMField *out, float64 coef, FMField *state, int32 offset,
	       FMField *bf, VolumeGeometry *vg,
	       int32 *conn, int32 nEl, int32 nEP,
	       int32 *elList, int32 elList_nRow,
	       int32 isDiff );

int32 dw_mass_scalar( FMField *out, FMField *state, int32 offset,
		      FMField *bf, VolumeGeometry *vg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow,
		      int32 isDiff );

int32 dw_mass_scalar_fine_coarse( FMField *out, FMField *state, int32 offset,
				  FMField *bf, FMField *cbfs,
				  VolumeGeometry *vg,
				  int32 *conn, int32 nEl, int32 nEP,
				  int32 *iemap, int32 iemap_nRow,
				  int32 *elList, int32 elList_nRow,
				  int32 isDiff );

END_C_DECLS

#endif /* Header */
