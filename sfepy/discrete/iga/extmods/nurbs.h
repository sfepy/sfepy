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

typedef struct NURBSContext {
  int32 (*get_xi_dist)(float64 *pdist, FMField *xi,
                       FMField *point, FMField *e_coors,
                       void *_ctx);
  int32 (*eval_basis)(FMField *out, FMField *coors, int32 diff,
                      void *_ctx);
  int32 iel; // >= 0.
  int32 is_dx; // 1 => apply reference mapping to gradient.
  FMField e_coors_max[1]; // Buffer for coordinates of element nodes.

  FMField control_points[1];
  FMField weights[1];
  int32 *degrees;
  int32 dim;
  FMField cs[3];
  int32 *conn;
  int32 n_cell;
  int32 n_efun;

  FMField bf[1];
  FMField bfg[1];

  FMField R[1];
  FMField dR_dxi[1];
  FMField dR_dx[1];

  FMField B[3];
  FMField dB_dxi[3];
  FMField N[3];
  FMField dN_dxi[3];

  int32 reuse;

  int32 has_bernstein;
  int32 is_nurbs;

  int32 i_max;
  float64 newton_eps;
} NURBSContext;

void ravel_multi_index(uint32 *index, uint32 *indices,
                       uint32 *shape, uint32 num);
void unravel_index(uint32 *indices, uint32 index,
                   uint32 *shape, uint32 num);

void print_context_nurbs(void *_ctx);

int32 get_xi_dist(float64 *pdist, FMField *xi,
                  FMField *point, FMField *e_coors,
                  void *_ctx);

int32 eval_basis_nurbs(FMField *out, FMField *coors, int32 diff,
                       void *_ctx);

int32 eval_bernstein_basis(FMField *funs, FMField *ders,
                           float64 x, uint32 degree);
int32 eval_bspline_basis_tp(FMField *R, FMField *dR_dx, FMField *det,
                            FMField *dR_dxi,
                            FMField *dx_dxi, FMField *dxi_dx,
                            FMField *B, FMField *dB_dxi,
                            FMField *N, FMField *dN_dxi,
                            FMField *qp, uint32 ie,
                            FMField *control_points,
                            int32 *degrees, int32 dim,
                            FMField *cs,
                            int32 *conn, int32 n_el, int32 n_ep,
                            int32 has_bernstein, int32 is_dx);
int32 eval_nurbs_basis_tp(FMField *R, FMField *dR_dx, FMField *det,
                          FMField *dR_dxi,
                          FMField *dx_dxi, FMField *dxi_dx,
                          FMField *B, FMField *dB_dxi,
                          FMField *N, FMField *dN_dxi,
                          FMField *qp, uint32 ie, FMField *control_points,
                          FMField *weights, int32 *degrees, int32 dim,
                          FMField *cs,
                          int32 *conn, int32 n_el, int32 n_ep,
                          int32 has_bernstein, int32 is_dx);

#endif /* !NURBS_H */
