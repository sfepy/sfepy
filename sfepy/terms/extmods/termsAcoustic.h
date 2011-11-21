/*!
  @par Revision history:
  - 15.07.2008, c
*/
#ifndef _TERMSACOUSTIC_H_
#define _TERMSACOUSTIC_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"

int32 d_diffusion_sa( FMField *out,
		      FMField *grad_q, FMField *grad_p,
		      FMField *grad_w, FMField *div_w,
		      FMField *mtxD, VolumeGeometry *vg );

int32 dw_surf_laplace( FMField *out, FMField *state, FMField *coef,
			FMField *gbf, SurfaceGeometry *sg,
		        int32 *conn, int32 nEl, int32 nEP,
		        int32 *elList, int32 elList_nRow,
		        int32 isDiff );

int32 d_surf_laplace( FMField *out, FMField *stateP, FMField *stateQ, FMField *coef,
		      FMField *gbf, SurfaceGeometry *sg,
		      int32 *conn, int32 nEl, int32 nEP,
		      int32 *elList, int32 elList_nRow);

int32 dw_surf_lcouple( FMField *out, FMField *state, FMField *coef,
		       FMField *bf, FMField *gbf, SurfaceGeometry *sg,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow,
		       int32 isDiff );

int32 d_surf_lcouple(FMField *out, FMField *stateP, FMField *stateQ, FMField *coef,
		     FMField *bf, FMField *gbf, SurfaceGeometry *sg,
		     int32 *conn, int32 nEl, int32 nEP,
		     int32 *elList, int32 elList_nRow);

END_C_DECLS

#endif /* Header */
