#include "geomtrans.h"
#include "refcoors.h"

#undef __FUNC__
#define __FUNC__ "_mul_c_add_v3"
void _mul_c_add_v3(float64 out[3], float64 v0[3],
                   float64 coef, float64 v1[3], int32 dim)
{
  out[0] = v0[0] + coef * v1[0];
  out[1] = v0[1] + coef * v1[1];
  if (dim == 3) {
    out[2] = v0[2] + coef * v1[2];
  }
}

#undef __FUNC__
#define __FUNC__ "_intersect_line_plane"
int32 _intersect_line_plane(float64 *pt, float64 p00[3], float64 p01[3],
                            float64 pp[3], float64 pn[3], int32 dim)
{
  int32 ii;
  float64 dot1, dot2;
  float64 v0[3], r[3];

  _mul_c_add_v3(v0, p01, -1.0, p00, dim);
  dot1 = 0.0;
  for (ii = 0; ii < dim; ii++) {
    dot1 += fabs(v0[ii]);
  }
  if (dot1 < 1e-10) {
    *pt = 0.0;
    goto end_label;
  }

  _mul_c_add_v3(r, p00, -1.0, pp, dim);

  gtr_dot_v3(&dot1, pn, v0, dim);
  if (fabs(dot1) <= 1e-10) {
    *pt = 1e10;
    goto end_label;
  }

  gtr_dot_v3(&dot2, pn, r, dim);
  *pt = - dot2 / dot1;

 end_label:
  return(RET_OK);
}

#undef __FUNC__
#define __FUNC__ "_intersect_line_triangle"
int32 _intersect_line_triangle(float64 *pt, float64 p00[3], float64 p01[3],
                               float64 pps[9], float64 pn[3])
{
  float64 duv, dwv, dwu, duu, dvv, aux, t1, t2;
  float64 u[3], v[3], w[3], tmp[3];

  _intersect_line_plane(pt, p00, p01, pps, pn, 3);

  // u = pps[1] - pps[0].
  _mul_c_add_v3(u, pps + 3, -1.0, pps, 3);
  // v = pps[2] - pps[0].
  _mul_c_add_v3(v, pps + 6, -1.0, pps, 3);
  // w = p00 + t * (p01 - p00) - pps[0].
  _mul_c_add_v3(tmp, p01, -1.0, p00, 3);
  _mul_c_add_v3(w, p00, *pt, tmp, 3);
  _mul_c_add_v3(w, w, -1.0, pps, 3);

  gtr_dot_v3(&duv, u, v, 3);
  gtr_dot_v3(&dwv, w, v, 3);
  gtr_dot_v3(&dwu, w, u, 3);
  gtr_dot_v3(&duu, u, u, 3);
  gtr_dot_v3(&dvv, v, v, 3);

  aux = duv * duv - duu * dvv;
  if (fabs(aux/duu) < 1e-14) {
    *pt = 1e10;
    goto end_label;
  }

  t1 = (duv * dwv - dvv * dwu) / aux;
  t2 = (duv * dwu - duu * dwv) / aux;

  if ((t1 < -1e-10) || (t2 < -1e-10) || ((t1 + t2) > (1.0 + 1e-10))) {
    *pt = 1e10;
  }

 end_label:
  return(RET_OK);
}

#undef __FUNC__
#define __FUNC__ "_get_cell_coors"
void _get_cell_coors(FMField *e_coors, Indices *cell_vertices,
                     float64 *mesh_coors, int32 dim, float64 *buf)
{
  int32 ir, ii, n_v;

  n_v = cell_vertices->num;
  fmf_pretend_nc(e_coors, 1, 1, n_v, dim, buf);
  for (ir = 0; ir < n_v; ir++) {
    for (ii = 0; ii < dim; ii++) {
      e_coors->val[dim*ir+ii] = mesh_coors[dim*cell_vertices->indices[ir]+ii];
    }
  }
}

#undef __FUNC__
#define __FUNC__ "_get_tri_coors"
void _get_tri_coors(float64 *buf9, uint32 *loc_indices, uint32 off,
                    uint32 *tri, float64 *mesh_coors, uint32 *cv_indices)
{
  int32 ir, ik;

  for (ir = 0; ir < 3; ir++) {
    ik = loc_indices[off + tri[ir]];
    buf9[3*ir+0] = mesh_coors[3*cv_indices[ik]+0];
    buf9[3*ir+1] = mesh_coors[3*cv_indices[ik]+1];
    buf9[3*ir+2] = mesh_coors[3*cv_indices[ik]+2];
  }
}

#undef __FUNC__
#define __FUNC__ "refc_find_ref_coors_convex"
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
                                 void *_ctx)
{
  BasisContext *ctx = (BasisContext *) _ctx;
  int32 ip, ic, icell, icell_max = 0, imin, ik, ok, ret = RET_OK;
  uint32 ii;
  int32 xi_ok, hexa_reverse;
  int32 D = mesh->topology->max_dim;
  int32 dim = D - 1;
  uint32 nc = mesh->geometry->dim;
  uint32 tri0[] = {0, 1, 3};
  uint32 tri1[] = {2, 3, 1};
  uint32 cell, cell0, cell00, facet;
  uint32 *noffs, *foffs, aux[2];
  uint32 *cell_types = mesh->topology->cell_types;
  float64 d_min, tmin, tt, dist;
  float64 *mesh_coors = mesh->geometry->coors;
  float64 buf3[3];
  float64 buf9[9];
  FMField point[1], centroid[1], _normals0[1], _normals1[1], e_coors[1], xi[1];
  Indices cell_vertices[1];
  MeshEntity cell_ent[1];
  MeshConnectivity *cD0 = 0; // D -> 0
  MeshConnectivity *c0D = 0; // 0 -> D
  MeshConnectivity *cDd = 0; // cell -> facet
  MeshConnectivity *cdD = 0; // facet -> cell
  MeshConnectivity *loc = 0;
  MeshConnectivity **locs = 0;

  mesh_setup_connectivity(mesh, D, 0);
  cD0 = mesh->topology->conn[IJ(D, D, 0)];

  mesh_setup_connectivity(mesh, 0, D);
  c0D = mesh->topology->conn[IJ(D, 0, D)];

  mesh_setup_connectivity(mesh, D, dim);
  cDd = mesh->topology->conn[IJ(D, D, dim)];
  noffs = cDd->offsets;

  mesh_setup_connectivity(mesh, dim, D);
  cdD = mesh->topology->conn[IJ(D, dim, D)];

  // Local entities - reference cell edges or faces.
  locs = (dim == 1) ? mesh->entities->edges : mesh->entities->faces;

  fmf_pretend_nc(point, coors->nRow, 1, 1, nc, coors->val);
  fmf_pretend_nc(centroid, centroids->nRow, 1, 1, nc, centroids->val);

  fmf_pretend_nc(xi, 1, 1, 1, nc, buf3);
  fmf_fillC(xi, 0.0);

  ctx->is_dx = 0;

  for (ip = 0; ip < coors->nRow; ip++) {
    ic = ics[ip];
    /* output("***** %d %d\n", ip, ic); */

    FMF_SetCell(point, ip);
    /* fmf_print(point, stdout, 0); */

    cell = cell0 = cell00 = c0D->indices[c0D->offsets[ic]];

    ok = icell = hexa_reverse = imin = 0;
    d_min = 1e10;
    while (1) {
      /* output("*** %d %d %d\n", icell, cell, hexa_reverse); */
      FMF_SetCell(centroid, cell);
      /* fmf_print(centroid, stdout, 0); */

      ctx->iel = cell;
      cell_ent->ii = cell;
      me_get_incident2(cell_ent, cell_vertices, cD0);

      loc = locs[cell_types[cell]];
      foffs = loc->offsets;

      if (cell_types[cell] != 4) { // No hexahedron -> planar facet.
        fmf_pretend_nc(_normals0, noffs[cell+1] - noffs[cell], 1, 1, nc,
                       normals0->val + nc * noffs[cell]);

        tmin = 1e10;
        for (ii = 0; ii < loc->num; ii++) {
          FMF_SetCell(_normals0, ii);
          ik = loc->indices[foffs[ii]];

          _intersect_line_plane(&tt, centroid->val, point->val,
                                mesh_coors + nc * cell_vertices->indices[ik],
                                _normals0->val, nc);
          if ((tt >= -qp_eps) && (tt < (tmin + qp_eps))) {
            imin = ii;
            tmin = tt;
          }
        }

        if (tmin >= (1.0 - qp_eps)) {
          _get_cell_coors(e_coors, cell_vertices, mesh_coors, nc,
                          ctx->e_coors_max->val);
          /* fmf_print(e_coors, stdout, 0); */

          xi_ok = ctx->get_xi_dist(&dist, xi, point, e_coors, ctx);

          d_min = Min(dist, d_min);
          if (xi_ok && (dist < qp_eps)) {
            ok = 1;
          }
          break;
        }

      } else { // Hexahedron -> non-planar facet in general.
        fmf_pretend_nc(_normals0, noffs[cell+1] - noffs[cell], 1, 1, nc,
                       normals0->val + nc * noffs[cell]);
        fmf_pretend_nc(_normals1, noffs[cell+1] - noffs[cell], 1, 1, nc,
                       normals1->val + nc * noffs[cell]);
        for (ii = 0; ii < loc->num; ii++) {
          FMF_SetCell(_normals0, ii);
          _get_tri_coors(buf9, loc->indices, foffs[ii],
                         tri0, mesh_coors, cell_vertices->indices);
          _intersect_line_triangle(&tt, centroid->val, point->val,
                                   buf9, _normals0->val);
          if ((tt >= -qp_eps) && (tt < 1e10)) {
            ok = 2;
            imin = ii;
            if ((tt >= (1.0 - qp_eps)) || hexa_reverse) {
              _get_cell_coors(e_coors, cell_vertices, mesh_coors, nc,
                              ctx->e_coors_max->val);

              xi_ok = ctx->get_xi_dist(&dist, xi, point, e_coors, ctx);

              d_min = Min(dist, d_min);
              if (xi_ok && (dist < qp_eps)) {
                ok = 1;
              } else {
                hexa_reverse = 1;
              }
            }
            break;
          }

          FMF_SetCell(_normals1, ii);
          _get_tri_coors(buf9, loc->indices, foffs[ii],
                         tri1, mesh_coors, cell_vertices->indices);
          _intersect_line_triangle(&tt, centroid->val, point->val,
                                   buf9, _normals1->val);
          if ((tt >= -qp_eps) && (tt < 1e10)) {
            ok = 2;
            imin = ii;
            if ((tt >= (1.0 - qp_eps)) || hexa_reverse) {
              _get_cell_coors(e_coors, cell_vertices, mesh_coors, nc,
                              ctx->e_coors_max->val);

              xi_ok = ctx->get_xi_dist(&dist, xi, point, e_coors, ctx);

              d_min = Min(dist, d_min);
              if (xi_ok && (dist < qp_eps)) {
                ok = 1;
              } else {
                hexa_reverse = 1;
              }
            }
            break;
          }
        }
        if (ok == 1) {
          break;
        }
        if (ok == 0) {
          errput("cannot intersect bilinear faces!\n");
          ERR_CheckGo(ret);
        }
      }

      facet = cDd->indices[cDd->offsets[cell] + imin];
      if ((cdD->offsets[facet+1] - cdD->offsets[facet]) == 2) {
        aux[0] = cdD->indices[cdD->offsets[facet]];
        aux[1] = cdD->indices[cdD->offsets[facet]+1];
        cell00 = cell0;
        cell0 = cell;
        cell = (aux[0] == cell) ? aux[1] : aux[0];

        if (cell == cell00) { // Ping-pong between two cells.
          hexa_reverse = 1;
        }

      } else { // Boundary facet.
        ctx->iel = cell;
        cell_ent->ii = cell;
        me_get_incident2(cell_ent, cell_vertices, cD0);

        _get_cell_coors(e_coors, cell_vertices, mesh_coors, nc,
                        ctx->e_coors_max->val);
        xi_ok = ctx->get_xi_dist(&dist, xi, point, e_coors, ctx);

        d_min = Min(dist, d_min);
        if (xi_ok && (dist < qp_eps)) {
          ok = 1;
        } else {
          ok = 0;
        }
        break;
      }

      icell++;
      icell_max = Max(icell, icell_max);
      if (icell > 10000) {
        errput("cannot find containing cell!\n");
        ERR_CheckGo(ret);
      }
    }

    /* fmf_print(xi, stdout, 0); */
    /* output("-> %d %d %d %.3e\n", cell, xi_ok, ok, d_min); */

    cells[ip] = cell;

    if (ok != 1) {
      if (!xi_ok) {
        status[ip] = 4;
      } else if (allow_extrapolation) {
        // Try using minimum distance xi.
        if (sqrt(d_min) < close_limit) {
          status[ip] = 1;
        } else {
          status[ip] = 2;
        }
      } else {
        status[ip] = 3;
      }
    } else {
      status[ip] = 0;
    }

    for (ii = 0; ii < nc; ii++) {
      ref_coors->val[nc*ip+ii] = xi->val[ii];
    }
  }
  /* output("%d\n", icell_max); */

 end_label:
  return(ret);
}

#undef __FUNC__
#define __FUNC__ "refc_find_ref_coors"
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
                          void *_ctx)
{
  BasisContext *ctx = (BasisContext *) _ctx;

  int32 ip, ic, ii, imin, ok, xi_ok, ret = RET_OK;
  int32 D = mesh->topology->max_dim;
  int32 nc = mesh->geometry->dim;
  float64 d_min, dist;
  float64 *mesh_coors = mesh->geometry->coors;
  float64 buf3[3];
  FMField point[1], e_coors[1], xi[1];
  Indices cell_vertices[1];
  MeshEntity cell_ent[1];
  MeshConnectivity *cD0 = 0; // D -> 0

  mesh_setup_connectivity(mesh, D, 0);
  cD0 = mesh->topology->conn[IJ(D, D, 0)];

  fmf_pretend_nc(point, coors->nRow, 1, 1, nc, coors->val);

  fmf_pretend_nc(xi, 1, 1, 1, nc, buf3);
  fmf_fillC(xi, 0.0);

  ctx->is_dx = 0;

  for (ip = 0; ip < coors->nRow; ip++) {
    FMF_SetCell(point, ip);

    if (offsets[ip] == offsets[ip+1]) {
      status[ip] = 5;
      cells[ip] = 0;
      for (ii = 0; ii < nc; ii++) {
        ref_coors->val[nc*ip+ii] = 0.0;
      }
      continue;
    }

    ok = xi_ok = 0;
    d_min = 1e10;
    imin = candidates[offsets[ip]];

    for (ic = offsets[ip]; ic < offsets[ip+1]; ic++) {
      /* output("***** %d %d %d\n", ip, ic, candidates[ic]); */

      ctx->iel = candidates[ic];
      cell_ent->ii = candidates[ic];
      me_get_incident2(cell_ent, cell_vertices, cD0);

      _get_cell_coors(e_coors, cell_vertices, mesh_coors, nc,
                      ctx->e_coors_max->val);
      xi_ok = ctx->get_xi_dist(&dist, xi, point, e_coors, ctx);

      if (xi_ok) {
        if (dist < qp_eps) {
          imin = cell_ent->ii;
          ok = 1;
          break;
        } else if (dist < d_min) {
          d_min = dist;
          imin = cell_ent->ii;
        }
      } else if (dist < d_min) {
        d_min = dist;
        imin = cell_ent->ii;
      }
    }
    /* output("-> %d %d %d %.3e\n", imin, xi_ok, ok, d_min); */

    cells[ip] = imin;

    if (ok != 1) {
      if (!xi_ok) {
        status[ip] = 4;
      } else if (allow_extrapolation) {
        if (sqrt(d_min) < close_limit) {
          status[ip] = 1;
        } else {
          status[ip] = 2;
        }
      } else {
        status[ip] = 3;
      }
    } else {
      status[ip] = 0;
    }

    for (ii = 0; ii < nc; ii++) {
      ref_coors->val[nc*ip+ii] = xi->val[ii];
    }
    ERR_CheckGo(ret);
  }

 end_label:
  return(ret);
}
