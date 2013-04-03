#ifndef _TERMS_DOT_H_
#define _TERMS_DOT_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "refmaps.h"

int32 dw_surface_v_dot_n_s(FMField *out,
                           FMField *coef, FMField *val_qp,
                           Mapping *rsg,
                           Mapping *csg,
                           int32 isDiff);
int32 dw_surface_s_v_dot_n(FMField *out,
                           FMField *coef, FMField *val_qp,
                           Mapping *rsg,
                           Mapping *csg,
                           int32 isDiff);
int32 dw_volume_dot_vector( FMField *out, FMField *coef, FMField *val_qp,
                            Mapping *rvg, Mapping *cvg,
                            int32 isDiff );
int32 dw_volume_dot_scalar( FMField *out, FMField *coef, FMField *val_qp,
                            Mapping *rvg, Mapping *cvg,
                            int32 isDiff );

int32 dw_v_dot_grad_s_vw( FMField *out, FMField *coef, FMField *grad,
                          Mapping *vvg, Mapping *svg,
                          int32 isDiff );
int32 dw_v_dot_grad_s_sw( FMField *out, FMField *coef, FMField *val_qp,
                          Mapping *vvg, Mapping *svg,
                          int32 isDiff );

END_C_DECLS

#endif /* Header */
