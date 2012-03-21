#ifndef _TERMS_DOT_H_
#define _TERMS_DOT_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"

int32 dw_volume_dot_vector( FMField *out, FMField *coef, FMField *val_qp,
                            FMField *rbf, FMField *cbf, VolumeGeometry *vg,
                            int32 isDiff );
int32 dw_surface_dot_vector( FMField *out, FMField *coef, FMField *val_qp,
                             FMField *rbf, FMField *cbf, SurfaceGeometry *sg,
                             int32 isDiff );
int32 dw_volume_dot_scalar( FMField *out, FMField *coef, FMField *val_qp,
                            FMField *rbf, FMField *cbf, VolumeGeometry *vg,
                            int32 isDiff );
int32 dw_surface_dot_scalar( FMField *out, FMField *coef, FMField *val_qp,
                             FMField *rbf, FMField *cbf, SurfaceGeometry *sg,
                             int32 isDiff );

int32 dw_v_dot_grad_s_vw( FMField *out, FMField *coef, FMField *grad,
                          FMField *vbf, VolumeGeometry *vvg,
                          VolumeGeometry *svg, int32 isDiff );
int32 dw_v_dot_grad_s_sw( FMField *out, FMField *coef, FMField *val_qp,
                          FMField *vbf, VolumeGeometry *vvg,
                          VolumeGeometry *svg, int32 isDiff );

END_C_DECLS

#endif /* Header */
