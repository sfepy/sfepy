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

typedef struct LagrangeContext {
  int32 (*get_xi_dist)(float64 *pdist, FMField *xi,
                       FMField *point, FMField *e_coors,
                       void *_ctx);

  FMField bc[1];
  FMField mtx_i[1];
  FMField base1d[1];
  FMField ref_coors[1];
  int32 *nodes;
  int32 n_nod;
  int32 n_col;
  int32 tdim;
  float64 eps;
  int32 check_errors;
  int32 i_max;
  float64 newton_eps;
  float64 vmin;
  float64 vmax;
} LagrangeContext;

void print_context_lagrange(void *_ctx);

int32 get_barycentric_coors(FMField *bc, FMField *coors, void *_ctx);

int32 get_xi_dist(float64 *pdist, FMField *xi,
                  FMField *point, FMField *e_coors,
                  void *_ctx);

int32 get_xi_simplex(FMField *xi, FMField *dest_point, FMField *e_coors,
                     void *_ctx);

int32 get_xi_tensor(FMField *xi, FMField *dest_point, FMField *e_coors,
                    void *_ctx);

int32 eval_lagrange_simplex(FMField *out, int32 order, int32 diff,
                            void *_ctx);

int32 eval_lagrange_tensor_product(FMField *out, int32 order, int32 diff,
                                   void *_ctx);

#endif /* !LAGRANGE_H */
