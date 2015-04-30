#ifndef _REFCOORS_H_
#define _REFCOORS_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "mesh.h"

int32 refc_find_ref_coors(FMField *ref_coors,
                          int32 *cells, int32 n_cells,
                          int32 *status, int32 n_status,
                          FMField *coors,
                          Mesh *mesh,
                          FMField *centroids,
                          FMField *normals0,
                          FMField *normals1,
                          int32 *ics, int32 n_ics,
                          FMField *eref_coors,
                          int32 *nodes, int32 n_nodes, int32 n_nodes_col,
                          FMField *mtx_i,
                          int32 allow_extrapolation,
                          float64 close_limit, float64 qp_eps,
                          int32 i_max, float64 newton_eps);

END_C_DECLS

#endif /* Header */
