#include "meshutils.h"
#include "geomtrans.h"
#include "sort.h"

int32 *ic;
int32 n_c;

#define SwapValues(a, b, work) do {\
  (work) = (a); (a) = (b); (b) = (work);\
} while (0)


/*!
  For 3D, 2D and 1D.

  @par Revision history:
  - 08.06.2006, c
  - 02.08.2006
*/
int32 orient_elements(int32 *flag, int32 flag_n_row,
                      Mesh *mesh, Indices *cells, int32 dcells,
                      int32 *v_roots, int32 v_roots_n_row,
                      int32 *v_vecs, int32 v_vecs_n_row, int32 v_vecs_n_col,
                      int32 *swap_from,
                      int32 swap_from_n_row, int32 swap_from_n_col,
                      int32 *swap_to,
                      int32 swap_to_n_row, int32 swap_to_n_col)
{
#define SWF(ir, is) (swap_from[swap_from_n_col*ir+is])
#define SWT(ir, is) (swap_to[swap_to_n_col*ir+is])
#define IR(ir) (cell_vertices->indices[v_roots[ir]])
#define IV(ir, iv) (cell_vertices->indices[v_vecs[v_vecs_n_col*ir+iv]])
#define CONN(ip) (cell_vertices->indices[ip])

  int32 D = mesh->topology->max_dim;
  float64 *coors = mesh->geometry->coors;
  Indices cell_vertices[1];
  MeshEntityIterator it0[1];
  MeshConnectivity *cD0 = 0; // D -> 0
  int32 ir, iel, ii, ip0, ip1, ip2, ip3, tmp, nc, tdim;
  float64 v0[3], v1[3], v2[3], v3[3], cross[3], dot[1], val;

  cD0 = mesh->topology->conn[IJ(D, D, 0)];

  nc = mesh->geometry->dim;
  tdim = mesh->topology->max_dim;

  for (mei_init_sub(it0, mesh, cells, dcells); mei_go(it0); mei_next(it0)) {
    iel = it0->entity->ii;
    flag[iel] = 0;

    me_get_incident2(it0->entity, cell_vertices, cD0);

    for (ir = 0; ir < v_roots_n_row; ir++) {
      ip0 = IR(ir);
      ip1 = IV(ir, 0);

      if (tdim == 1) {
        v0[0] = coors[ip0];
        v1[0] = coors[ip1] - v0[0];
        val = v1[0];
      } else if (tdim >= 2) {
        ip2 = IV(ir, 1);
        v1[2] = v2[2] = 0.0;
        for (ii = 0; ii < nc; ii++) {
          v0[ii] = coors[nc*ip0+ii];
          v1[ii] = coors[nc*ip1+ii] - v0[ii];
          v2[ii] = coors[nc*ip2+ii] - v0[ii];
        }
        gtr_cross_product(cross, v1, v2);
        val = cross[2];
      }

      if (tdim == 3) {
        ip3 = IV(ir, 2);
        for (ii = 0; ii < nc; ii++) {
          v3[ii] = coors[nc*ip3+ii] - v0[ii];
        }
        gtr_dot_v3(dot, v3, cross, 3);
        val = dot[0];
      }

      if (val < 0.0) {
        flag[iel]++;
        for (ii = 0; ii < swap_from_n_col; ii++) {
          SwapValues(CONN(SWF(ir, ii)), CONN(SWT(ir, ii)), tmp);
        }
      }
    }
  }

  return(RET_OK);

#undef IR
#undef IV
#undef CONN
#undef SWF
#undef SWT
}
