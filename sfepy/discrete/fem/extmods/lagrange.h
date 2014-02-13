#ifndef _LAGRANGE_H_
#define _LAGRANGE_H_ 1

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

int32 get_barycentric_coors(FMField *bc, FMField *coors, FMField *mtx_i,
                            float64 eps, int32 check_errors);

int32 get_xi_simplex(FMField *xi, FMField *bc, FMField *dest_point,
                     FMField *ref_coors, FMField *e_coors);

int32 get_xi_tensor(FMField *xi,
                    FMField *dest_point, FMField *e_coors,
                    FMField *mtx_i,
                    FMField *base1d, int32 *nodes, int32 n_col,
                    float64 vmin, float64 vmax,
                    int32 i_max, float64 newton_eps);

int32 eval_lagrange_simplex(FMField *out, FMField *bc, FMField *mtx_i,
                            int32 *nodes, int32 n_col,
                            int32 order, int32 diff);

int32 eval_lagrange_tensor_product(FMField *out, FMField *bc,
                                   FMField *mtx_i, FMField *base1d,
                                   int32 *nodes, int32 n_col,
                                   int32 order, int32 diff);

#endif /* !LAGRANGE_H */
