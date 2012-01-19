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

int32 dw_laplace( FMField *out, FMField *grad,
		  FMField *coef, VolumeGeometry *vg,
                  int32 isDiff );
int32 d_laplace( FMField *out, FMField *gradP1, FMField *gradP2,
		 FMField *coef, VolumeGeometry *vg );

int32 dw_diffusion( FMField *out, FMField *grad,
		    FMField *mtxD, VolumeGeometry *vg,
		    int32 isDiff );
int32 d_diffusion( FMField *out, FMField *gradP1, FMField *gradP2,
		   FMField *mtxD, VolumeGeometry *vg );
int32 dw_permeability_r( FMField *out, FMField *mtxD, VolumeGeometry *vg );
int32 d_surface_flux( FMField *out, FMField *grad,
                      FMField *mtxD, SurfaceGeometry *sg, int32 mode );

END_C_DECLS

#endif /* Header */
