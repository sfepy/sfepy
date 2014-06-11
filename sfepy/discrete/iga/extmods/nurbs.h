#ifndef _NURBS_H_
#define _NURBS_H_ 1

#ifdef __cplusplus
#  define BEGIN_C_DECLS         extern "C" {
#  define END_C_DECLS           }
#else
#  define BEGIN_C_DECLS
#  define END_C_DECLS
#endif

#include "types.h"
#include "version.h"

#include "fmfield.h"
#include "geommech.h"

void ravel_multi_index(uint32 *index, uint32 *indices,
                       uint32 *shape, uint32 num);
void unravel_index(uint32 *indices, uint32 index,
                   uint32 *shape, uint32 num);

int32 eval_bernstein_basis(FMField *funs, FMField *ders,
                           float64 x, uint32 degree);
int32 eval_nurbs_basis_tp(FMField *R, FMField *dR_dx, FMField *det,
                          FMField *dR_dxi,
                          FMField *dx_dxi, FMField *dxi_dx,
                          FMField *B, FMField *dB_dxi,
                          FMField *N, FMField *dN_dxi,
                          FMField *qp, uint32 ie, FMField *control_points,
                          FMField *weights, int32 *degrees, int32 dim,
                          FMField *cs,
                          int32 *conn, int32 n_el, int32 n_ep);

#endif /* !NURBS_H */
