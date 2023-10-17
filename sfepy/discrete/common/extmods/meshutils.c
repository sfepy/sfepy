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

int32 compare_i32(const void *a, const void *b)
{
  int32 i1, i2;

  i1 = *((int32 *) a);
  i2 = *((int32 *) b);

  return(i1 - i2);
}

#undef __FUNC__
#define __FUNC__ "mesh_nod_in_el_count"
/*!
  @par Revision history:
  - 21.11.2003, c
  - 23.11.2003
*/
int32 mesh_nod_in_el_count(int32 *p_niec_max, int32 *niec,
                           int32 n_nod, int32 n_gr, int32 *n_el,
                           int32 *n_ep, int32 **conn)
{
  int32 ig, iel, iep, in, niec_max;
  int32 *pconn;

  memset(niec, 0, (n_nod + 1) * sizeof(int32));
  for (ig = 0; ig < n_gr; ig++) {
    for (iel = 0; iel < n_el[ig]; iel++) {
      pconn = conn[ig] + n_ep[ig] * iel;
      for (iep = 0; iep < n_ep[ig]; iep++) {
        niec[1+pconn[iep]]++;
        /* output("%d %d %d\n", iep, niec[1+pconn[iep]], pconn[iep]); */
      }
    }
  }

  niec[0] = 0;
  niec_max = 0;
  for (in = 0; in <= n_nod; in++) {
    /* output("%d %d\n", in, niec[in]); */
    niec_max = Max(niec_max, niec[in]);
  }
  *p_niec_max = niec_max;

  return(RET_OK);
}

#undef __FUNC__
#define __FUNC__ "mesh_graph"
/*!
  @par Revision history:
  - 23.05.2003, c
  - 26.05.2003
  - 27.05.2003
  - 28.05.2003
  - 21.11.2003 former mesh_mesh_graph()
  - 23.11.2003
  - 01.03.2004
  - 03.03.2005
  - 07.02.2006
*/
int32 mesh_graph(int32 *p_nnz, int32 **p_prow, int32 **p_icol,
                 int32 n_row, int32 n_col, int32 n_gr, int32 *n_el,
                 int32 *n_epr, int32 **conn_r, int32 *n_epc, int32 **conn_c)
{
  int32 ret = RET_OK, in, ii, ip, ig, iel, iep, ir, ic, nn, np, pr,
    niec_max_r, n_ep_max_c, n_unique, iir, iic, found;
  int32 *niec, *pconn_r, *pconn_c, *eonlist, *nir, *nods, *icol;

/*   output("%d %d %d %d %d %d\n", n_row, n_col, n_gr, n_el[0], n_epr[0], n_epc[0]); */

  /* Get niec (= nodes in elements count) for rows. */
  niec = alloc_mem(int32, n_row + 1);
  ERR_CheckGo(ret);
  mesh_nod_in_el_count(&niec_max_r, niec, n_row, n_gr, n_el, n_epr, conn_r);
/*   output("%d\n", niec_max_r); */

  /* Cummulative sum. */
  for (in = 0; in < n_row; in++) {
    niec[in+1] += niec[in];
  }

/*    output("00\n"); */

  /* eon = elements of nodes */
  nn = 0;
  n_ep_max_c = 0;
  for (ig = 0; ig < n_gr; ig++) {
    nn += n_epr[ig] * n_el[ig];
    n_ep_max_c = Max(n_ep_max_c, n_epc[ig]);
  }
  eonlist = alloc_mem(int32, 2 * nn);

  /* nir is just a buffer here. */
  nir = alloc_mem(int32, n_row + 1);
  ERR_CheckGo(ret);
  memset(nir, 0, (n_row + 1) * sizeof(int32));

/*    output("1\n"); */

  /* Get list of elements each row node is in. */
  for (ig = 0; ig < n_gr; ig++) {
    for (iel = 0; iel < n_el[ig]; iel++) {
      pconn_r = conn_r[ig] + n_epr[ig] * iel;
      for (iep = 0; iep < n_epr[ig]; iep++) {
        np = pconn_r[iep];
        if (np >= 0) {
          eonlist[2*(niec[np]+nir[np])+0] = iel;
          eonlist[2*(niec[np]+nir[np])+1] = ig;
/*      output("  %d %d %d %d\n", np, eonlist[2*(niec[np]+nir[np])+0], */
/*                 eonlist[2*(niec[np]+nir[np])+1], nir[np]); */
          nir[np]++;
        }
      }
    }
  }

/*    output("2\n"); */

  /* nir = number in row. */
  memset(nir, 0, (n_row + 1) * sizeof(int32));

  /* List of column nodes for each row node. */
/*   output("%d, %d\n", n_ep_max_c, niec_max_r * n_ep_max_c); */
  nods = alloc_mem(int32, niec_max_r * n_ep_max_c);

  nn = 0;
  for (in = 0; in < n_row; in++) {
    ii = 0;
/*      output("%d\n", in); */
    for (ip = niec[in]; ip < niec[in+1]; ip++) {
      iel = eonlist[2*(ip)+0];
      ig = eonlist[2*(ip)+1];
/*        output(" %d %d %d\n", ip, ig, iel); */
      for (iep = 0; iep < n_epc[ig]; iep++) {
        np = conn_c[ig][n_epc[ig]*iel+iep];
        if (np >= 0) {
          nods[ii] = np;
/*      output("  %d %d\n", ii, nods[ii]); */
          ii++;
        }
      }
    }
/*     output("%d\n", ii); */

    if (ii > 0) {
/*       qsort(nods, ii, sizeof(int32), &compareI32); */
      int32_quicksort(nods, ii, 0);
      n_unique = 1;
      for (ir = 0; ir < (ii - 1); ir++) {
        if (nods[ir] != nods[ir+1]) {
          n_unique++;
        }
      }
    } else {
      n_unique = 0;
    }
    nn += n_unique;
/*      output(" -> %d\n", n_unique); */

    nir[in] = n_unique;
  }

/*    output("3\n"); */

  *p_nnz = nn;
  *p_prow = niec;
  icol = *p_icol = alloc_mem(int32, nn);
  ERR_CheckGo(ret);

  /* Fill in *p_prow. */
  niec[0] = 0;
  for (in = 0; in < n_row; in++) {
    niec[in+1] = niec[in] + nir[in];
/*      output(" %d\n", niec[in+1]); */
  }

/*   output("4\n"); */
  /* Fill in *p_icol (sorted). */
  memset(nir, 0, (n_row + 1) * sizeof(int32));
  for (ig = 0; ig < n_gr; ig++) {
/*     output("ig %d\n", ig); */
    for (iel = 0; iel < n_el[ig]; iel++) {
      pconn_r = conn_r[ig] + n_epr[ig] * iel;
      pconn_c = conn_c[ig] + n_epc[ig] * iel;
      for (ir = 0; ir < n_epr[ig]; ir++) {
        iir = pconn_r[ir];
        if (iir < 0) continue;
        pr = niec[iir];
/*      output(" %d %d %d\n", iir, pr, niec[iir+1] - pr); */
        for (ic = 0; ic < n_epc[ig]; ic++) {
          iic = pconn_c[ic];
          if (iic < 0) continue;
/*        output("   %d %d\n", iic, nir[iir]); */
          /* This is a bottle-neck! */
          found = 0;
          for (ii = pr; ii < (pr + nir[iir]); ii++) {
            if (icol[ii] == iic) {
              found = 1;
              break;
            }
          }
/*        output("  ? %d\n", found); */
          if (!found) {
            if (nir[iir] < (niec[iir+1] - pr)) {
              icol[pr+nir[iir]] = iic;
              nir[iir]++;
/*            output("  + %d %d\n", nir[iir], niec[iir+1] - pr); */
            } else {
              output("  %d %d\n", nir[iir], niec[iir+1] - pr);
              errput("ERR_VerificationFail\n");
              ERR_CheckGo(ret);
            }
          }
        }
/*      qsort(icol + pr, nir[iir], sizeof(int32), &compareI32); */
        int32_quicksort(icol + pr, nir[iir], 0);
      }
    }
  }

/*   output("5\n"); */
 end_label:

  free_mem(nods);
  free_mem(nir);
  free_mem(eonlist);

  return(ret);
}

/*!
  @par Revision history:
  - 06.03.2005, c
  - 07.03.2005
  - 29.08.2007, from gr_components() (00.01.18)
*/
int32 graph_components(int32 *p_n_comp,
                       int32 *flag, int32 flag_len,
                       int32 *row, int32 row_len,
                       int32 *col, int32 col_len,
                       int32 *pos, int32 pos_len)
{
  // pos is a work array: list of nodes (rows) to process.
  int32 ret = RET_OK, n_tot, n_pos, n_pos_new, n_pos0, n_new, n_stop, n_nod;
  int32 icomp, ii, ir, ic;

  n_nod = row_len - 1;

  n_stop = n_nod;
  for (ir = 0; ir < n_nod; ir++) {
    flag[ir] = -1;
    if ((row[ir+1] - row[ir]) == 0) {
      n_stop--;
      flag[ir] = -2;
    }
  }

  n_tot = 0;
  for (icomp = 0; icomp < n_nod; icomp++) {
    // Find seed.
    ii = 0;
    while ((flag[ii] >= 0) || (flag[ii] == -2)) {
      ii++;
      if (ii >= n_nod) {
        errput("error in graph_components()!\n");
        ERR_CheckGo(ret);
      }
    }
    flag[ii] = icomp;
    pos[0] = ii;
    n_pos0 = 0;
    n_pos_new = n_pos = 1;

    for (ii = 0; ii < n_nod; ii++) {
      n_new = 0;
      for (ir = n_pos0; ir < n_pos; ir++) {
        for (ic = row[pos[ir]]; ic < row[pos[ir]+1]; ic++) {
          if (flag[col[ic]] == -1) {
            flag[col[ic]] = icomp;
            pos[n_pos_new] = col[ic];
            n_pos_new++;
            n_new++;
          }
        }
      }
/*       output("%d %d %d\n", ii, n_pos, n_new); */
      n_pos0 = n_pos;
      n_pos = n_pos_new;
      if (n_new == 0) break;
    }
    n_tot += n_pos;
/*     output("  %d %d %d %d\n", icomp, n_tot, n_stop, n_nod); */
/*     getchar(); */
    if (n_tot == n_stop) {
      *p_n_comp = icomp + 1;
      break;
    }
  }

 end_label:
  return(ret);
}
