#ifndef _TERMS_DOT_H_
#define _TERMS_DOT_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"

int32 dw_volume_dot_vector( FMField *out, FMField *coef, FMField *val_qp,
                            VolumeGeometry *rvg, VolumeGeometry *cvg,
                            int32 isDiff );
int32 dw_surface_dot_vector( FMField *out, FMField *coef, FMField *val_qp,
                             SurfaceGeometry *rsg, SurfaceGeometry *csg,
                             int32 isDiff );
int32 dw_volume_dot_scalar( FMField *out, FMField *coef, FMField *val_qp,
                            VolumeGeometry *rvg, VolumeGeometry *cvg,
                            int32 isDiff );
int32 dw_surface_dot_scalar( FMField *out, FMField *coef, FMField *val_qp,
                             SurfaceGeometry *rsg, SurfaceGeometry *csg,
                             int32 isDiff );

int32 dw_v_dot_grad_s_vw( FMField *out, FMField *coef, FMField *grad,
                          VolumeGeometry *vvg, VolumeGeometry *svg,
                          int32 isDiff );
int32 dw_v_dot_grad_s_sw( FMField *out, FMField *coef, FMField *val_qp,
                          VolumeGeometry *vvg, VolumeGeometry *svg,
                          int32 isDiff );

END_C_DECLS

#endif /* Header */
