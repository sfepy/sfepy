#ifndef _REFCOORS_H_
#define _REFCOORS_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "mesh.h"

typedef struct BasisContext {
  int32 (*get_xi_dist)(float64 *pdist, FMField *xi,
                       FMField *dest_point, FMField *e_coors,
                       void *_ctx);
  int32 (*eval_basis)(FMField *out, FMField *coors, int32 diff,
                      void *_ctx);
  int32 iel; // >= 0.
  int32 is_dx; // 1 => apply reference mapping to gradient.
  FMField e_coors_max[1]; // Buffer for coordinates of element nodes.
} BasisContext;

int32 refc_find_ref_coors_convex(FMField *ref_coors,
                                 int32 *cells, int32 n_cells,
                                 int32 *status, int32 n_status,
                                 FMField *coors,
                                 Mesh *mesh,
                                 FMField *centroids,
                                 FMField *normals0,
                                 FMField *normals1,
                                 int32 *ics, int32 n_ics,
                                 int32 allow_extrapolation,
                                 float64 qp_eps,
                                 float64 close_limit,
                                 void *_ctx);

int32 refc_find_ref_coors(FMField *ref_coors,
                          int32 *cells, int32 n_cells,
                          int32 *status, int32 n_status,
                          FMField *coors,
                          Mesh *mesh,
                          int32 *candidates, int32 n_candidates,
                          int32 *offsets, int32 n_offsets,
                          int32 allow_extrapolation,
                          float64 qp_eps,
                          float64 close_limit,
                          void *_ctx);

END_C_DECLS

#endif /* Header */
