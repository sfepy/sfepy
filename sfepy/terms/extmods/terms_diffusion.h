/*!
  @par Revision history:
  - 28.11.2005, c
*/
#ifndef _TERMSLAPLACE_H_
#define _TERMSLAPLACE_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "refmaps.h"

int32 dw_laplace( FMField *out, FMField *grad,
		  FMField *coef, Mapping *vg,
                  int32 isDiff );
int32 d_laplace( FMField *out, FMField *gradP1, FMField *gradP2,
		 FMField *coef, Mapping *vg );

int32 dw_diffusion( FMField *out, FMField *grad,
		    FMField *mtxD, Mapping *vg,
		    int32 isDiff );
int32 d_diffusion( FMField *out, FMField *gradP1, FMField *gradP2,
		   FMField *mtxD, Mapping *vg );
int32 d_sd_diffusion(FMField *out,
                     FMField *grad_q, FMField *grad_p,
                     FMField *grad_w, FMField *div_w,
                     FMField *mtxD, Mapping *vg);

int32 dw_diffusion_r( FMField *out, FMField *mtxD, Mapping *vg );
int32 d_surface_flux( FMField *out, FMField *grad,
                      FMField *mtxD, Mapping *sg, int32 mode );
int32 dw_surface_flux(FMField *out, FMField *grad,
                      FMField *mat, FMField *bf, Mapping *sg,
                      int32 *fis, int32 nFa, int32 nFP, int32 mode);

int32 dw_convect_v_grad_s( FMField *out, FMField *val_v, FMField *grad_s,
                           Mapping *vvg, Mapping *svg,
                           int32 isDiff );

END_C_DECLS

#endif /* Header */
