#include "nurbs.h"

#ifdef _WIN32
// If inline is used here on Windows, msvc linker does not find this function.
void ravel_multi_index(uint32 *index, uint32 *indices,
                       uint32 *shape, uint32 num)
#else
inline void ravel_multi_index(uint32 *index, uint32 *indices,
                              uint32 *shape, uint32 num)
#endif
{
  uint32 ii, stride = 1;
  uint32 raveled = 0;

  for (ii = num - 1; ii >= 1; ii--) {
    raveled += stride * indices[ii];
    stride *= shape[ii];
  }
  raveled += stride * indices[0];

  *index = raveled;
}

inline void unravel_index(uint32 *indices, uint32 index,
                          uint32 *shape, uint32 num)
{
  uint32 ii;

  for (ii = num - 1; ii > 0; ii--) {
    indices[ii] = index % shape[ii];
    index /= shape[ii];
  }
  indices[0] = index % shape[0];
}

void print_context_nurbs(void *_ctx)
{
  NURBSContext *ctx = (NURBSContext *) _ctx;

  int32 ir, ic;

  output("iel: %d\n", ctx->iel);
  output("is_dx: %d\n", ctx->is_dx);

  output("e_coors_max:\n");
  fmf_print(ctx->e_coors_max, stdout, 1);

  output("control_points:\n");
  fmf_print(ctx->control_points, stdout, 0);
  output("weights:\n");
  fmf_print(ctx->weights, stdout, 0);

  output("degrees:\n");
  for (ir = 0; ir < ctx->dim; ir++) {
    output(" %d", ctx->degrees[ir]);
  }
  output("\n");
  output("dim: %d\n", ctx->dim);

  output("cs:\n");
  for (ir = 0; ir < ctx->dim; ir++) {
    fmf_print(ctx->cs + ir, stdout, 0);
  }

  output("conn:\n");
  for (ir = 0; ir < ctx->n_cell; ir++) {
    for (ic = 0; ic < ctx->n_efun; ic++) {
      output(" %d", ctx->conn[ctx->n_efun*ir+ic]);
    }
    output("\n");
  }
  output("n_cell: %d\n", ctx->n_cell);
  output("n_efun: %d\n", ctx->n_efun);

  output("bf:\n");
  fmf_print(ctx->bf, stdout, 1);
  output("bfg:\n");
  fmf_print(ctx->bfg, stdout, 1);

  output("R:\n");
  fmf_print(ctx->R, stdout, 1);
  output("dR_dxi:\n");
  fmf_print(ctx->dR_dxi, stdout, 1);
  output("dR_dx:\n");
  fmf_print(ctx->dR_dx, stdout, 1);

  output("B:\n");
  for (ir = 0; ir < ctx->dim; ir++) {
    fmf_print(ctx->B + ir, stdout, 1);
  }
  output("dB_dxi:\n");
  for (ir = 0; ir < ctx->dim; ir++) {
    fmf_print(ctx->dB_dxi + ir, stdout, 1);
  }

  output("N:\n");
  for (ir = 0; ir < ctx->dim; ir++) {
    fmf_print(ctx->N + ir, stdout, 1);
  }
  output("dN_dxi:\n");
  for (ir = 0; ir < ctx->dim; ir++) {
    fmf_print(ctx->dN_dxi + ir, stdout, 1);
  }

  output("reuse: %d\n", ctx->reuse);

  output("has_bernstein: %d\n", ctx->has_bernstein);
  output("is_nurbs: %d\n", ctx->is_nurbs);

  output("i_max: %d\n", ctx->i_max);
  output("newton_eps: %.4e\n", ctx->newton_eps);
}


#undef __FUNC__
#define __FUNC__ "get_xi_dist"
// Returns: ok = 1: success, ok = 0: failure.
int32 get_xi_dist(float64 *pdist, FMField *xi,
                  FMField *point, FMField *e_coors,
                  void *_ctx)
{
  NURBSContext *ctx = (NURBSContext *) _ctx;

  int32 i_max = ctx->i_max;
  float64 newton_eps = ctx->newton_eps;

  int32 idim, ii, ok = 0;
  int32 dim = e_coors->nCol;
  FMField *bf = ctx->bf;
  FMField *bfg = ctx->bfg;
  FMField xint[1], res[1], mtx[1], imtx[1], bc[1], base1d[1];
  float64 err, val, dist;
  float64 buf6[6], buf8[8];
  float64 buf3_1[3], buf3_2[3], buf9_1[9], buf9_2[9];

  fmf_pretend_nc(base1d, 1, 1, 1, e_coors->nRow, buf8);

  fmf_pretend_nc(bc, dim, 1, 1, 2, buf6);
  fmf_pretend_nc(res, 1, 1, 1, dim, buf3_1);
  fmf_pretend_nc(xint, 1, 1, 1, dim, buf3_2);
  fmf_pretend_nc(mtx, 1, 1, dim, dim, buf9_1);
  fmf_pretend_nc(imtx, 1, 1, dim, dim, buf9_2);

  ii = 0;
  fmf_fillC(xi, 0.5);
  while (ii < i_max) {
    // Base(xi).
    ctx->reuse = 0;
    ctx->eval_basis(bf, xi, 0, _ctx);

    // X(xi).
    fmf_mulAB_n1(xint, bf, e_coors);
    // Rezidual.
    fmf_subAB_nn(res, point, xint);

    err = 0.0;
    for (idim = 0; idim < dim; idim++) {
      err += res->val[idim] * res->val[idim];
    }
    err = sqrt(err);

    if (err < newton_eps) {
      ok = 1;
      break;
    }

    // grad Base(xi).
    ctx->reuse = 1;
    ctx->eval_basis(bfg, xi, 1, _ctx);

    // - Matrix.
    fmf_mulAB_n1(mtx, bfg, e_coors);

    geme_invert3x3(imtx, mtx);

    fmf_mulAB_nn(xint, res, imtx);
    fmf_addAB_nn(xi, xi, xint);

    ii += 1;
  }
  if (!ok) {
    // Newton did not converge.
    // Use centre if nothing better is available.
    fmf_fillC(xi, 0.5);

    *pdist = 1e10;

  } else {
    dist = 0.0;
    for (ii = 0; ii < dim; ii++) {
      val = Min(Max(xi->val[ii] - 1.0, 0.0), 100.0);
      dist += val * val;
      val = Min(Max(0.0 - xi->val[ii], 0.0), 100.0);
      dist += val * val;
    }
    *pdist = dist;
  }

  return(ok);
}

/*
  Evaluate NURBS basis polynomials in given points.
*/
int32 eval_basis_nurbs(FMField *out, FMField *coors, int32 diff,
                       void *_ctx)
{
  NURBSContext *ctx = (NURBSContext *) _ctx;
  int32 dim = ctx->dim;
  int32 ii, iqp, ret = RET_OK;
  FMField bf[1], _coors[1], qp[1], _det[1], _dx_dxi[1], _dxi_dx[1];
  float64 buf1[1], buf9_1[9], buf9_2[9];

  fmf_pretend_nc(_coors, 1, coors->nRow, 1, coors->nCol, coors->val);
  fmf_pretend_nc(qp, 1, 1, 1, dim, 0);
  fmf_pretend_nc(bf, 1, 1, out->nRow, out->nCol, 0);
  fmf_pretend_nc(_det, 1, 1, 1, 1, buf1);
  fmf_pretend_nc(_dx_dxi, 1, 1, dim, dim, buf9_1);
  fmf_pretend_nc(_dxi_dx, 1, 1, dim, dim, buf9_2);

  if ((ctx->reuse) && (out->nLev > 1)) {
    errput("cannot reuse more than one point! (%d)\n", out->nLev);
    ERR_CheckGo(ret);
  }

  for (iqp = 0; iqp < out->nLev; iqp++) {
    fmf_set_qp(bf, iqp, out);

    if (ctx->reuse) {
      goto copy;
    }

    fmf_set_qp(qp, iqp, _coors);

    for (ii = 0; ii < dim; ii++) {
      eval_bernstein_basis(ctx->B + ii, ctx->dB_dxi + ii, qp->val[ii],
                           ctx->degrees[ii]);
    }

    if (ctx->is_nurbs) {
      eval_nurbs_basis_tp(ctx->R, ctx->dR_dx, _det,
                          ctx->dR_dxi, _dx_dxi, _dxi_dx,
                          ctx->B, ctx->dB_dxi, ctx->N, ctx->dN_dxi,
                          qp, ctx->iel,
                          ctx->control_points, ctx->weights,
                          ctx->degrees, ctx->dim, ctx->cs,
                          ctx->conn, ctx->n_cell, ctx->n_efun,
                          1, ctx->is_dx);
    } else {
      eval_bspline_basis_tp(ctx->R, ctx->dR_dx, _det,
                            ctx->dR_dxi, _dx_dxi, _dxi_dx,
                            ctx->B, ctx->dB_dxi, ctx->N, ctx->dN_dxi,
                            qp, ctx->iel,
                            ctx->control_points,
                            ctx->degrees, ctx->dim, ctx->cs,
                            ctx->conn, ctx->n_cell, ctx->n_efun,
                            1, ctx->is_dx);
    }

  copy:
    if (!diff) {
      fmf_copy(bf, ctx->R);

    } else if (ctx->is_dx) {
      fmf_copy(bf, ctx->dR_dx);

    } else {
      fmf_copy(bf, ctx->dR_dxi);
    }

    ERR_CheckGo(ret);
  }

 end_label:
  return(ret);
}

int32 eval_bernstein_basis(FMField *funs, FMField *ders,
                           float64 x, uint32 degree)
{
  uint32 ip, ifun;
  uint32 n_fun = degree + 1;
  float64 prev, tmp;

  fmf_fillC(funs, 0.0);
  fmf_fillC(ders, 0.0);

  funs->val[0] = 1.0;

  if (degree == 0) {
    return(RET_OK);
  }

  for (ip = 1; ip < n_fun - 1; ip++) {
    prev = 0.0;
    for (ifun = 0; ifun < ip + 1; ifun++) {
      tmp = x * funs->val[ifun];
      funs->val[ifun] = (1.0 - x) * funs->val[ifun] + prev;
      prev = tmp;
    }
  }

  ders->val[0] = degree * (- funs->val[0]);
  for (ifun = 1; ifun < n_fun; ifun++) {
    ders->val[ifun] = degree * (funs->val[ifun - 1] - funs->val[ifun]);
  }

  prev = 0.0;
  for (ifun = 0; ifun < n_fun; ifun++) {
    tmp = x * funs->val[ifun];
    funs->val[ifun] = (1.0 - x) * funs->val[ifun] + prev;
    prev = tmp;
  }

  return(RET_OK);
}

#undef __FUNC__
#define __FUNC__ "eval_bspline_basis_tp"
/*
  dR_dx has shape (dim, n_efun), transposed w.r.t. the Python version!
*/
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
                            int32 has_bernstein, int32 is_dx)
{
  int32 ret = RET_OK;
  uint32 ii, jj, a, i0, i1, i2;
  uint32 n_efuns[3];
  uint32 n_efun = 1;
  uint32 n_els[3];
  uint32 ic[3];
  int32 *ec;
  FMField *C, *N0, *N1, *N2, *dN0_dxi, *dN1_dxi, *dN2_dxi;
  float64 P;

#ifdef DEBUG_FMF
  if (!((dim == qp->nCol) && (dim <= 3))) {
    errput(ErrHead "inconsistent dimension! (%d == $d <= 3)", dim, qp->nCol);
    ERR_CheckGo(ret);
  }
#endif

  for (ii = 0; ii < (uint32)dim; ii++) {
    n_efuns[ii] = degrees[ii] + 1;
    n_efun *= n_efuns[ii];
  }

#ifdef DEBUG_FMF
  if (n_efun != n_ep) {
    errput(ErrHead "inconsistent number of bases! (%d == $d)", n_efun, n_ep);
    ERR_CheckGo(ret);
  }
#endif

  // Element connectivity.
  ec = conn + n_ep * ie;

  // 1D Bernstein basis B, dB/dxi.
  if (!has_bernstein) {
    for (ii = 0; ii < (uint32)dim; ii++) {
      eval_bernstein_basis(B + ii, dB_dxi + ii, qp->val[ii], degrees[ii]);
    }
  }

  // 1D B-spline basis N = CB, dN/dxi = C dB/dxi.
  for (ii = 0; ii < (uint32)dim; ii++) {
    n_els[ii] = (cs + ii)->nCell;
  }

  unravel_index(ic, ie, n_els, dim);

  for (ii = 0; ii < (uint32)dim; ii++) {
    C = cs + ii;
    FMF_SetCell(C, ic[ii]);

    fmf_mulAB_nn(N + ii, C, B + ii);
    fmf_mulAB_nn(dN_dxi + ii, C, dB_dxi + ii);
  }

  ERR_CheckGo(ret);

  // Tensor-product B-spline basis R, dR/dxi.
  a = 0; // Basis function index.
  if (dim == 3) {
    N0 = N + 0;
    N1 = N + 1;
    N2 = N + 2;
    dN0_dxi = dN_dxi + 0;
    dN1_dxi = dN_dxi + 1;
    dN2_dxi = dN_dxi + 2;
    for (i0 = 0; i0 < n_efuns[0]; i0++) {
      for (i1 = 0; i1 < n_efuns[1]; i1++) {
        for (i2 = 0; i2 < n_efuns[2]; i2++) {
          R->val[a] = N0->val[i0] * N1->val[i1] * N2->val[i2];

          dR_dxi->val[a+n_ep*0] = dN0_dxi->val[i0] * N1->val[i1] * N2->val[i2];

          dR_dxi->val[a+n_ep*1] = N0->val[i0] * dN1_dxi->val[i1] * N2->val[i2];

          dR_dxi->val[a+n_ep*2] = N0->val[i0] * N1->val[i1] * dN2_dxi->val[i2];

          a += 1;
        }
      }
    }
  } else if (dim == 2) {
    N0 = N + 0;
    N1 = N + 1;
    dN0_dxi = dN_dxi + 0;
    dN1_dxi = dN_dxi + 1;
    for (i0 = 0; i0 < n_efuns[0]; i0++) {
      for (i1 = 0; i1 < n_efuns[1]; i1++) {
        R->val[a] = N0->val[i0] * N1->val[i1];

        dR_dxi->val[a+n_ep*0] = dN0_dxi->val[i0] * N1->val[i1];

        dR_dxi->val[a+n_ep*1] = N0->val[i0] * dN1_dxi->val[i1];

        a += 1;
      }
    }
  } else {
    // Simple copy here.
    N0 = N + 0;
    dN0_dxi = dN_dxi + 0;
    for (i0 = 0; i0 < n_efuns[0]; i0++) {
        R->val[a] = N0->val[i0];

        dR_dxi->val[a+n_ep*0] = dN0_dxi->val[i0];

        a += 1;
    }
  }

  if (is_dx) {
    // Mapping reference -> physical domain dxi/dx.
    // x = sum P_a R_a, dx/dxi = sum P_a dR_a/dxi, invert.
    for (ii = 0; ii < (uint32)dim; ii++) {
      for (jj = 0; jj < (uint32)dim; jj++) {
        dx_dxi->val[dim*ii+jj] = 0.0;
        for (a = 0; a < (uint32)dR_dxi->nCol; a++) {
          P = control_points->val[dim*ec[a]+ii];

          dx_dxi->val[dim*ii+jj] += P * dR_dxi->val[a+n_ep*jj];
        }
      }
    }
    geme_det3x3(det->val, dx_dxi);

    geme_invert3x3(dxi_dx, dx_dxi);

    // dR/dx.
    fmf_mulATB_nn(dR_dx, dxi_dx, dR_dxi);
  }

 end_label:
  return(ret);
}

#undef __FUNC__
#define __FUNC__ "eval_nurbs_basis_tp"
/*
  dR_dx has shape (dim, n_efun), transposed w.r.t. the Python version!
*/
int32 eval_nurbs_basis_tp(FMField *R, FMField *dR_dx, FMField *det,
                          FMField *dR_dxi,
                          FMField *dx_dxi, FMField *dxi_dx,
                          FMField *B, FMField *dB_dxi,
                          FMField *N, FMField *dN_dxi,
                          FMField *qp, uint32 ie, FMField *control_points,
                          FMField *weights, int32 *degrees, int32 dim,
                          FMField *cs,
                          int32 *conn, int32 n_el, int32 n_ep,
                          int32 has_bernstein, int32 is_dx)
{
  int32 ret = RET_OK;
  uint32 ii, jj, a, i0, i1, i2;
  uint32 n_efuns[3];
  uint32 n_efun = 1;
  uint32 n_els[3];
  uint32 ic[3];
  int32 *ec;
  FMField *C, *N0, *N1, *N2, *dN0_dxi, *dN1_dxi, *dN2_dxi;
  float64 w, W, P;
  float64 dw_dxi[3];

#ifdef DEBUG_FMF
  if (!((dim == qp->nCol) && (dim <= 3))) {
    errput(ErrHead "inconsistent dimension! (%d == $d <= 3)", dim, qp->nCol);
    ERR_CheckGo(ret);
  }
#endif

  for (ii = 0; ii < (uint32)dim; ii++) {
    n_efuns[ii] = degrees[ii] + 1;
    n_efun *= n_efuns[ii];
  }

#ifdef DEBUG_FMF
  if (n_efun != n_ep) {
    errput(ErrHead "inconsistent number of bases! (%d == $d)", n_efun, n_ep);
    ERR_CheckGo(ret);
  }
#endif

  // Element connectivity.
  ec = conn + n_ep * ie;

  // 1D Bernstein basis B, dB/dxi.
  if (!has_bernstein) {
    for (ii = 0; ii < (uint32)dim; ii++) {
      eval_bernstein_basis(B + ii, dB_dxi + ii, qp->val[ii], degrees[ii]);
    }
  }

  // 1D B-spline basis N = CB, dN/dxi = C dB/dxi.
  for (ii = 0; ii < (uint32)dim; ii++) {
    n_els[ii] = (cs + ii)->nCell;
  }

  unravel_index(ic, ie, n_els, dim);

  for (ii = 0; ii < (uint32)dim; ii++) {
    C = cs + ii;
    FMF_SetCell(C, ic[ii]);

    fmf_mulAB_nn(N + ii, C, B + ii);
    fmf_mulAB_nn(dN_dxi + ii, C, dB_dxi + ii);
  }

  ERR_CheckGo(ret);

  // Numerators and denominator for tensor-product NURBS basis R, dR/dxi.
  w = 0; // w_b
  for (ii = 0; ii < (uint32)dim; ii++) {
    dw_dxi[ii] = 0.0; // dw_b/dxi
  }
  a = 0; // Basis function index.
  if (dim == 3) {
    N0 = N + 0;
    N1 = N + 1;
    N2 = N + 2;
    dN0_dxi = dN_dxi + 0;
    dN1_dxi = dN_dxi + 1;
    dN2_dxi = dN_dxi + 2;
    for (i0 = 0; i0 < n_efuns[0]; i0++) {
      for (i1 = 0; i1 < n_efuns[1]; i1++) {
        for (i2 = 0; i2 < n_efuns[2]; i2++) {
          W = weights->val[ec[a]];

          R->val[a] = N0->val[i0] * N1->val[i1] * N2->val[i2] * W;
          w += R->val[a];

          dR_dxi->val[a+n_ep*0] = dN0_dxi->val[i0] * N1->val[i1] * N2->val[i2] * W;
          dw_dxi[0] += dR_dxi->val[a+n_ep*0];

          dR_dxi->val[a+n_ep*1] = N0->val[i0] * dN1_dxi->val[i1] * N2->val[i2] * W;
          dw_dxi[1] += dR_dxi->val[a+n_ep*1];

          dR_dxi->val[a+n_ep*2] = N0->val[i0] * N1->val[i1] * dN2_dxi->val[i2] * W;
          dw_dxi[2] += dR_dxi->val[a+n_ep*2];

          a += 1;
        }
      }
    }
  } else if (dim == 2) {
    N0 = N + 0;
    N1 = N + 1;
    dN0_dxi = dN_dxi + 0;
    dN1_dxi = dN_dxi + 1;
    for (i0 = 0; i0 < n_efuns[0]; i0++) {
      for (i1 = 0; i1 < n_efuns[1]; i1++) {
        W = weights->val[ec[a]];

        R->val[a] = N0->val[i0] * N1->val[i1] * W;
        w += R->val[a];

        dR_dxi->val[a+n_ep*0] = dN0_dxi->val[i0] * N1->val[i1] * W;
        dw_dxi[0] += dR_dxi->val[a+n_ep*0];

        dR_dxi->val[a+n_ep*1] = N0->val[i0] * dN1_dxi->val[i1] * W;
        dw_dxi[1] += dR_dxi->val[a+n_ep*1];

        a += 1;
      }
    }
  } else {
    N0 = N + 0;
    dN0_dxi = dN_dxi + 0;
    for (i0 = 0; i0 < n_efuns[0]; i0++) {
        W = weights->val[ec[a]];

        R->val[a] = N0->val[i0] * W;
        w += R->val[a];

        dR_dxi->val[a+n_ep*0] = dN0_dxi->val[i0] * W;
        dw_dxi[0] += dR_dxi->val[a+n_ep*0];

        a += 1;
    }
  }

  // Finish R <- R / w_b.
  fmf_mulC(R, 1.0 / w);

  // Finish dR/dxi. D == W C dB/dxi, dR/dxi = (D - R dw_b/dxi) / w_b.
  for (a = 0; a < (uint32)dR_dxi->nCol; a++) {
    for (ii = 0; ii < (uint32)dim; ii++) {
      dR_dxi->val[a+n_ep*ii] = (dR_dxi->val[a+n_ep*ii]
                                - R->val[a] * dw_dxi[ii]) / w;
    }
  }

  if (is_dx) {
    // Mapping reference -> physical domain dxi/dx.
    // x = sum P_a R_a, dx/dxi = sum P_a dR_a/dxi, invert.
    for (ii = 0; ii < (uint32)dim; ii++) {
      for (jj = 0; jj < (uint32)dim; jj++) {
        dx_dxi->val[dim*ii+jj] = 0.0;
        for (a = 0; a < (uint32)dR_dxi->nCol; a++) {
          P = control_points->val[dim*ec[a]+ii];

          dx_dxi->val[dim*ii+jj] += P * dR_dxi->val[a+n_ep*jj];
        }
      }
    }
    geme_det3x3(det->val, dx_dxi);

    geme_invert3x3(dxi_dx, dx_dxi);

    // dR/dx.
    fmf_mulATB_nn(dR_dx, dxi_dx, dR_dxi);
  }

 end_label:
  return(ret);
}
