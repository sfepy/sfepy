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

int32 dw_laplace( FMField *out, FMField *state, int32 offset,
		  FMField *coef, VolumeGeometry *vg,
		  int32 *conn, int32 nEl, int32 nEP,
		  int32 *elList, int32 elList_nRow,
		  int32 isDiff );
int32 d_laplace( FMField *out, FMField *gradP1, FMField *gradP2,
		 FMField *coef, VolumeGeometry *vg,
		 int32 *elList, int32 elList_nRow );

int32 dw_diffusion( FMField *out, FMField *state, int32 offset,
		    FMField *mtxD, VolumeGeometry *vg,
		    int32 *conn, int32 nEl, int32 nEP,
		    int32 *elList, int32 elList_nRow,
		    int32 isDiff );
int32 d_diffusion( FMField *out, FMField *gradP1, FMField *gradP2,
		   FMField *mtxD, VolumeGeometry *vg,
		   int32 *elList, int32 elList_nRow );
int32 dw_permeability_r( FMField *out, FMField *mtxD, VolumeGeometry *vg,
			 int32 *conn, int32 nEl, int32 nEP,
			 int32 *elList, int32 elList_nRow );
int32 de_diffusion_velocity( FMField *out, FMField *state, int32 offset,
			     FMField *mtxD, VolumeGeometry *vg,
			     int32 *conn, int32 nEl, int32 nEP,
			     int32 *elList, int32 elList_nRow );

END_C_DECLS

#endif /* Header */
